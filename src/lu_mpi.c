/* ===========================================================================
 *  lu_mpi.c — Factorización LU con pivoteo parcial en PARALELO (MPI)
 *
 *  Proyecto 6 · Métodos Numéricos II y Computación Científica · FACET-UNT 2026
 *
 *  ---------------------------------------------------------------------------
 *  IDEA GENERAL
 *
 *  La matriz se reparte por FILAS entre los p procesos. Nadie tiene la matriz
 *  entera: cada proceso aloca solo sus ~n/p filas. Eso es lo que permite
 *  correr con n = 8000 sin que un proceso necesite 512 MB.
 *
 *  Soportamos dos repartos, seleccionables con -d, porque compararlos es uno
 *  de los experimentos del informe:
 *
 *    CÍCLICO (-d cyclic, el bueno)        BLOQUES (-d block, el malo)
 *    dueño(fila i) = i mod p              dueño(fila i) = i / ceil(n/p)
 *
 *    fila 0 -> P0                         fila 0 -> P0
 *    fila 1 -> P1                         fila 1 -> P0
 *    fila 2 -> P2                         fila 2 -> P1
 *    fila 3 -> P3                         fila 3 -> P1
 *    fila 4 -> P0                         fila 4 -> P2
 *    fila 5 -> P1                         fila 5 -> P2
 *    ...                                  ...
 *
 *  ¿Por qué el cíclico es mejor? Porque la eliminación gaussiana va
 *  "consumiendo" la matriz de arriba hacia abajo: en el paso k, las filas
 *  0..k ya están terminadas. Con reparto en bloques, los procesos dueños de
 *  las primeras filas se quedan SIN TRABAJO a mitad de camino y el resto
 *  carga con todo. Con reparto cíclico, la diferencia de carga entre dos
 *  procesos nunca supera una fila, en cualquier paso k.
 *
 *  ---------------------------------------------------------------------------
 *  EL PASO k, EN DETALLE
 *
 *    1. BÚSQUEDA DEL PIVOTE  (colectiva)
 *       Nadie tiene la columna k completa, así que cada proceso busca su
 *       máximo local en |a_ik| entre sus filas con índice global >= k, y
 *       después se hace MPI_Allreduce con MPI_MAXLOC sobre MPI_DOUBLE_INT.
 *       MAXLOC devuelve el par (valor máximo, índice que lo tiene) en una
 *       sola operación: con MPI_MAX habría que hacer una segunda
 *       comunicación para averiguar de quién era.
 *
 *    2. INTERCAMBIO DE FILAS  (punto a punto, o local)
 *       Si las filas k y prow viven en el mismo proceso, es un swap local
 *       y no cuesta nada. Si viven en procesos distintos, un
 *       MPI_Sendrecv_replace entre los dos. Usamos Sendrecv y no Send+Recv
 *       porque Send+Recv simétrico es un deadlock esperando a ocurrir.
 *
 *    3. DIFUSIÓN DE LA FILA PIVOTE  (colectiva)
 *       El dueño de la fila k la manda a todos con MPI_Bcast, porque todos
 *       la necesitan para actualizar sus propias filas. Difundimos solo los
 *       elementos k..n-1: los de la izquierda ya son parte de L y nadie los
 *       usa. Eso ahorra la mitad del volumen total.
 *
 *    4. ACTUALIZACIÓN  (100% local, cero comunicación)
 *       Cada proceso actualiza sus filas con índice global > k. Este es el
 *       trabajo de verdad: O(n³/p) en total.
 *
 *  Repitiendo n veces: O(n) colectivas de O(log p) etapas moviendo O(n²)
 *  datos, contra O(n³/p) de cómputo. El cociente comunicación/cómputo es
 *  del orden de p/n, y de ahí sale la conclusión central del informe:
 *  la eficiencia mejora cuando N crece a p fijo.
 *
 *  ---------------------------------------------------------------------------
 *  Uso:
 *      mpirun -n 8 ./lu_mpi -n 2000 -k diagdom -s 42
 *      mpirun -n 8 ./lu_mpi -n 2000 -d block --csv
 *      mpirun -n 8 ./lu_mpi -n 2000 --split-idle      (mide T_ocioso aparte)
 * ===========================================================================*/

#include "lu.h"
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#define TAG_SWAP 77

/* ===========================================================================
 *  DESCRIPTOR DE LA DISTRIBUCIÓN
 *
 *  Todo el mapeo "índice global <-> (proceso, índice local)" pasa por acá.
 *  Tenerlo en un solo lugar es lo que permite cambiar de reparto con un flag
 *  sin tocar el algoritmo.
 * ===========================================================================*/

typedef struct {
    int n, p, rank;
    int block;        /* 0 = cíclico, 1 = bloques contiguos */
    int rows_per;     /* techo(n/p): filas por proceso en el reparto en bloques */
    int nloc;         /* filas que realmente tiene este proceso */
    int maxloc;       /* filas alocadas (= rows_per, igual para los dos repartos) */
} dist_t;

/* ¿Qué proceso es dueño de la fila global i? */
static inline int owner_of(const dist_t *d, int i)
{
    return d->block ? (i / d->rows_per) : (i % d->p);
}

/* ¿En qué posición local del arreglo de su dueño está la fila global i? */
static inline int local_of(const dist_t *d, int i)
{
    return d->block ? (i - (i / d->rows_per) * d->rows_per) : (i / d->p);
}

/* ¿Qué fila global es la fila local il de ESTE proceso? */
static inline int global_of(const dist_t *d, int il)
{
    return d->block ? (d->rank * d->rows_per + il) : (d->rank + il * d->p);
}

/* Primer índice local cuya fila global es >= k.
 * Evita recorrer desde 0 todas las filas ya terminadas en cada paso:
 * sin esto el lazo de búsqueda del pivote sería O(n²/p) en vez de O(n²/2p). */
static inline int first_local_ge(const dist_t *d, int k)
{
    if (d->block) {
        int il = k - d->rank * d->rows_per;
        return il < 0 ? 0 : il;
    } else {
        int diff = k - d->rank;
        return (diff <= 0) ? 0 : (diff + d->p - 1) / d->p;   /* techo(diff/p) */
    }
}

static void dist_init(dist_t *d, int n, int p, int rank, int block)
{
    d->n        = n;
    d->p        = p;
    d->rank     = rank;
    d->block    = block;
    d->rows_per = (n + p - 1) / p;
    d->maxloc   = d->rows_per;

    if (block) {
        int start = rank * d->rows_per;
        int cnt   = n - start;
        if (cnt < 0)             cnt = 0;
        if (cnt > d->rows_per)   cnt = d->rows_per;
        d->nloc = cnt;
    } else {
        /* En el reparto cíclico, el proceso r tiene las filas r, r+p, r+2p...
         * Son n/p filas, más una si r < n mod p. */
        d->nloc = n / p + (rank < n % p ? 1 : 0);
    }
}

/* ===========================================================================
 *  MAIN
 * ===========================================================================*/

static void usage(int rank, const char *prog)
{
    if (rank == 0) fprintf(stderr,
      "\nUso: mpirun -n <p> %s [opciones]\n"
      "  -n <int>        dimension de la matriz (default 1000)\n"
      "  -k <tipo>       random | diagdom | hilbert | zerodiag (default diagdom)\n"
      "  -s <int>        semilla (default 42)\n"
      "  -r <int>        repeticiones cronometradas (default 1)\n"
      "  -d <reparto>    cyclic | block  (default cyclic)\n"
      "  --split-idle    poner barreras para medir T_ocioso por separado\n"
      "  --csv           una linea CSV en vez del reporte legible\n\n", prog);
}

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    int p, rank;
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    int n = 1000, seed = 42, reps = 1, block = 0, split_idle = 0, csv = 0;
    mat_kind_t kind = MAT_DIAGDOM;

    for (int i = 1; i < argc; i++) {
        if      (!strcmp(argv[i], "-n") && i + 1 < argc) n    = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-s") && i + 1 < argc) seed = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-r") && i + 1 < argc) reps = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-k") && i + 1 < argc) kind = mat_kind_from_string(argv[++i]);
        else if (!strcmp(argv[i], "-d") && i + 1 < argc) block = !strcmp(argv[++i], "block");
        else if (!strcmp(argv[i], "--split-idle")) split_idle = 1;
        else if (!strcmp(argv[i], "--csv"))        csv        = 1;
        else { usage(rank, argv[0]); MPI_Finalize(); return 1; }
    }
    if (n <= 0 || p <= 0) { usage(rank, argv[0]); MPI_Finalize(); return 1; }

    dist_t d;
    dist_init(&d, n, p, rank, block);

    /* ---------------------------------------------------------------------
     *  MEMORIA
     *
     *  A     : las filas locales, se factoriza in-place
     *  A0    : copia de las filas locales ORIGINALES, para el residuo
     *  pivrow: buffer donde llega la fila pivote difundida
     *  b0    : lado derecho original, replicado (n doubles: despreciable)
     *  pb    : copia de b0 que se va permutando durante la factorización
     *  y, x  : vectores de las sustituciones, replicados
     * ------------------------------------------------------------------- */
    size_t locsz = (size_t)d.maxloc * (size_t)n;

    double *A      = (double *)malloc(locsz * sizeof(double));
    double *A0     = (double *)malloc(locsz * sizeof(double));
    double *pivrow = (double *)malloc((size_t)n * sizeof(double));
    double *b0     = (double *)malloc((size_t)n * sizeof(double));
    double *pb     = (double *)malloc((size_t)n * sizeof(double));
    double *yv     = (double *)malloc((size_t)n * sizeof(double));
    double *xv     = (double *)malloc((size_t)n * sizeof(double));

    if (!A || !A0 || !pivrow || !b0 || !pb || !yv || !xv) {
        fprintf(stderr, "[rank %d] ERROR: sin memoria (necesita ~%.1f MB)\n",
                rank, 2.0 * locsz * sizeof(double) / 1048576.0);
        MPI_Abort(MPI_COMM_WORLD, 2);
    }

    /* Cada proceso genera SOLO sus filas. No hay Scatter de la matriz:
     * nunca existe una copia completa en ningún lado. */
    for (int il = 0; il < d.nloc; il++)
        mat_generate_row(&A0[(size_t)il * n], n, global_of(&d, il), kind, (unsigned)seed);
    vec_generate(b0, n, (unsigned)seed);

    double t_tot = 1e300, t_cmp = 0, t_com = 0, t_idl = 0;
    int    singular_at = -1;

    for (int rep = 0; rep < reps; rep++) {

        memcpy(A,  A0, locsz * sizeof(double));
        memcpy(pb, b0, (size_t)n * sizeof(double));

        double c_cmp = 0, c_com = 0, c_idl = 0, ts;

        MPI_Barrier(MPI_COMM_WORLD);          /* largada pareja */
        double t0 = MPI_Wtime();

        /* =================================================================
         *  LAZO PRINCIPAL DE LA FACTORIZACIÓN
         * ===============================================================*/
        for (int k = 0; k < n && singular_at < 0; k++) {

            /* ---- (1a) máximo LOCAL en la columna k --------------------- */
            ts = MPI_Wtime();
            struct { double val; int idx; } loc, glb;
            loc.val = -1.0;
            loc.idx = -1;
            for (int il = first_local_ge(&d, k); il < d.nloc; il++) {
                int g = global_of(&d, il);
                if (g < k) continue;
                double v = fabs(A[(size_t)il * n + k]);
                if (v > loc.val) { loc.val = v; loc.idx = g; }
            }
            c_cmp += MPI_Wtime() - ts;

            /* ---- (1b) máximo GLOBAL: la colectiva del pivoteo ----------
             * MPI_MAXLOC sobre MPI_DOUBLE_INT devuelve (mayor valor, índice
             * del que lo tiene). Los procesos sin candidato mandan val=-1,
             * que nunca gana porque los valores reales son >= 0. */
            if (split_idle) { ts = MPI_Wtime(); MPI_Barrier(MPI_COMM_WORLD); c_idl += MPI_Wtime() - ts; }
            ts = MPI_Wtime();
            MPI_Allreduce(&loc, &glb, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
            c_com += MPI_Wtime() - ts;

            if (glb.val == 0.0 || glb.idx < 0) { singular_at = k; break; }

            int prow = glb.idx;
            int ok   = owner_of(&d, k);
            int op   = owner_of(&d, prow);

            /* ---- (2) intercambio de las filas k y prow ------------------ */
            if (prow != k) {
                if (ok == op) {
                    /* Las dos filas viven acá: swap local, sin comunicación. */
                    if (rank == ok) {
                        double *ra = &A[(size_t)local_of(&d, k)    * n];
                        double *rb = &A[(size_t)local_of(&d, prow) * n];
                        for (int j = 0; j < n; j++) { double t = ra[j]; ra[j] = rb[j]; rb[j] = t; }
                    }
                } else {
                    ts = MPI_Wtime();
                    MPI_Status st;
                    if (rank == ok)
                        MPI_Sendrecv_replace(&A[(size_t)local_of(&d, k) * n], n,
                                             MPI_DOUBLE, op, TAG_SWAP, op, TAG_SWAP,
                                             MPI_COMM_WORLD, &st);
                    else if (rank == op)
                        MPI_Sendrecv_replace(&A[(size_t)local_of(&d, prow) * n], n,
                                             MPI_DOUBLE, ok, TAG_SWAP, ok, TAG_SWAP,
                                             MPI_COMM_WORLD, &st);
                    c_com += MPI_Wtime() - ts;
                }
                /* pb está replicado: todos aplican el mismo swap, sin comunicar. */
                double t = pb[k]; pb[k] = pb[prow]; pb[prow] = t;
            }

            /* ---- (3) difusión de la fila pivote ------------------------- */
            if (rank == ok)
                memcpy(&pivrow[k], &A[(size_t)local_of(&d, k) * n + k],
                       (size_t)(n - k) * sizeof(double));

            if (split_idle) { ts = MPI_Wtime(); MPI_Barrier(MPI_COMM_WORLD); c_idl += MPI_Wtime() - ts; }
            ts = MPI_Wtime();
            /* Solo los elementos k..n-1: los de la izquierda ya son L. */
            MPI_Bcast(&pivrow[k], n - k, MPI_DOUBLE, ok, MPI_COMM_WORLD);
            c_com += MPI_Wtime() - ts;

            /* ---- (4) actualización local: el trabajo de verdad ---------- */
            ts = MPI_Wtime();
            double inv_piv = 1.0 / pivrow[k];
            for (int il = first_local_ge(&d, k); il < d.nloc; il++) {
                int g = global_of(&d, il);
                if (g <= k) continue;                 /* la fila k no se toca */

                double *row = &A[(size_t)il * n];
                double  l   = row[k] * inv_piv;
                row[k] = l;                            /* acá se guarda L */
                if (l == 0.0) continue;

                for (int j = k + 1; j < n; j++)
                    row[j] -= l * pivrow[j];
            }
            c_cmp += MPI_Wtime() - ts;
        }

        double t1 = MPI_Wtime();
        double dt = t1 - t0;

        /* El tiempo del paso es el del proceso MÁS LENTO: un programa
         * paralelo termina cuando termina el último. */
        double dt_max;
        MPI_Allreduce(&dt, &dt_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

        if (dt_max < t_tot) { t_tot = dt_max; t_cmp = c_cmp; t_com = c_com; t_idl = c_idl; }

        int sing_any;
        MPI_Allreduce(&singular_at, &sing_any, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
        if (sing_any >= 0) { singular_at = sing_any; break; }
    }

    if (singular_at >= 0) {
        if (rank == 0)
            fprintf(stderr, "FALLO: matriz singular, pivote nulo en k=%d\n", singular_at);
        MPI_Finalize();
        return 3;
    }

    /* =====================================================================
     *  SUSTITUCIONES TRIANGULARES DISTRIBUIDAS
     *
     *  L·y = Pb  y después  U·x = y.
     *
     *  Cada fila i la tiene un solo proceso, y ese proceso tiene la fila
     *  COMPLETA (los n elementos). Entonces el dueño de la fila i puede
     *  calcular y_i por sí solo, siempre que ya conozca y_0..y_{i-1}.
     *  Por eso después de cada cálculo hace un Bcast de UN double.
     *
     *  Costo: n colectivas de 1 elemento. Es O(n log p) mensajes para
     *  O(n²) de cómputo — despreciable frente a los n³ de la factorización,
     *  pero es el pedazo menos escalable del programa. La forma seria de
     *  mejorarlo es hacerlo por bloques de b filas; lo dejamos así porque
     *  no cambia las conclusiones y el código queda legible.
     * ===================================================================*/
    {
        double ts = MPI_Wtime();

        /* Sustitución hacia adelante: L·y = Pb  (L con unos en la diagonal) */
        for (int i = 0; i < n; i++) {
            int oi = owner_of(&d, i);
            if (rank == oi) {
                const double *row = &A[(size_t)local_of(&d, i) * n];
                double s = pb[i];
                for (int j = 0; j < i; j++) s -= row[j] * yv[j];
                yv[i] = s;
            }
            MPI_Bcast(&yv[i], 1, MPI_DOUBLE, oi, MPI_COMM_WORLD);
        }

        /* Sustitución hacia atrás: U·x = y */
        for (int i = n - 1; i >= 0; i--) {
            int oi = owner_of(&d, i);
            if (rank == oi) {
                const double *row = &A[(size_t)local_of(&d, i) * n];
                double s = yv[i];
                for (int j = i + 1; j < n; j++) s -= row[j] * xv[j];
                xv[i] = s / row[i];
            }
            MPI_Bcast(&xv[i], 1, MPI_DOUBLE, oi, MPI_COMM_WORLD);
        }

        (void)ts;   /* el solve no entra en el tiempo reportado: medimos la
                       factorización, que es lo que pide la consigna */
    }

    /* =====================================================================
     *  RESIDUO, CALCULADO TAMBIÉN EN PARALELO
     *
     *      r = ||A·x − b||∞ / (||A||∞ · ||x||∞)
     *
     *  Usamos A0 (las filas ORIGINALES, sin permutar) y b0 (el lado derecho
     *  original). Cada proceso calcula el máximo sobre sus filas y se hace
     *  un Allreduce con MPI_MAX. Nunca se junta la matriz en un solo lugar.
     * ===================================================================*/
    double loc_rmax = 0.0, loc_amax = 0.0;
    for (int il = 0; il < d.nloc; il++) {
        int g = global_of(&d, il);
        const double *row = &A0[(size_t)il * n];

        double s = 0.0, a = 0.0;
        for (int j = 0; j < n; j++) { s += row[j] * xv[j]; a += fabs(row[j]); }

        double r = fabs(s - b0[g]);
        if (r > loc_rmax) loc_rmax = r;
        if (a > loc_amax) loc_amax = a;
    }

    double rmax, amax;
    MPI_Allreduce(&loc_rmax, &rmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
    MPI_Allreduce(&loc_amax, &amax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

    double xmax = norm_inf_vec(xv, n);
    double res  = (amax * xmax > 0.0) ? rmax / (amax * xmax) : rmax;

    /* Sumas globales de los tiempos por componente, promediadas entre procesos:
     * así el gráfico de barras apiladas del informe representa al conjunto y
     * no a un proceso particular. */
    double s_cmp, s_com, s_idl;
    MPI_Reduce(&t_cmp, &s_cmp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&t_com, &s_com, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    MPI_Reduce(&t_idl, &s_idl, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        double avg_cmp = s_cmp / p, avg_com = s_com / p, avg_idl = s_idl / p;
        double gflops  = lu_flops(n) / t_tot / 1e9;

        if (csv) {
            /* impl,n,p,dist,kind,variant,t_total,t_comp,t_comm,t_idle,gflops,residual,growth
             * (growth queda vacio: calcularlo requiere una reduccion extra
             *  sobre toda la matriz y no aporta al analisis de escalabilidad) */
            printf("mpi,%d,%d,%s,%s,pivot,%.6f,%.6f,%.6f,%.6f,%.4f,%.6e,\n",
                   n, p, block ? "block" : "cyclic", mat_kind_name(kind),
                   t_tot, avg_cmp, avg_com, avg_idl, gflops, res);
        } else {
            printf("\n  Factorizacion LU paralela (MPI)\n");
            printf("  ------------------------------------------------\n");
            printf("  n                     : %d\n", n);
            printf("  procesos              : %d\n", p);
            printf("  reparto               : %s\n", block ? "bloques contiguos" : "ciclico por filas");
            printf("  matriz                : %s\n", mat_kind_name(kind));
            printf("  filas por proceso     : ~%d  (%.1f MB por proceso)\n",
                   d.maxloc, 2.0 * locsz * sizeof(double) / 1048576.0);
            printf("  ------------------------------------------------\n");
            printf("  T total (factorizac.) : %.6f s\n", t_tot);
            printf("    T computo  (prom)   : %.6f s   (%.1f%%)\n", avg_cmp, 100.0 * avg_cmp / t_tot);
            printf("    T comunic. (prom)   : %.6f s   (%.1f%%)\n", avg_com, 100.0 * avg_com / t_tot);
            if (split_idle)
                printf("    T ocioso   (prom)   : %.6f s   (%.1f%%)\n", avg_idl, 100.0 * avg_idl / t_tot);
            printf("  rendimiento           : %.3f GFLOP/s\n", gflops);
            printf("  residuo relativo      : %.6e\n", res);
            printf("  ------------------------------------------------\n");
            printf("  %s\n\n", res < 1e-11
                     ? "OK: residuo en el orden esperado."
                     : "ATENCION: residuo alto. Revisar.");
            if (!split_idle)
                printf("  (Nota: sin --split-idle, el tiempo de espera por\n"
                       "   desbalance queda contado dentro de T comunicacion.)\n\n");
        }
    }

    free(A); free(A0); free(pivrow); free(b0); free(pb); free(yv); free(xv);
    MPI_Finalize();
    return 0;
}
