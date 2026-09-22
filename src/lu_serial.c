/* ===========================================================================
 *  lu_serial.c — Factorización LU con pivoteo parcial, versión SERIAL
 *
 *  Este archivo tiene dos papeles en el proyecto:
 *
 *    1. Es la IMPLEMENTACIÓN DE REFERENCIA. Todo lo que devuelve la versión
 *       MPI se compara contra lo que devuelve esta.
 *    2. Es el T₁ del speedup. Ojo: T₁ es ESTE programa, no la versión MPI
 *       corrida con un proceso. Usar la versión MPI con p=1 como T₁ infla
 *       el speedup porque arrastra el mismo overhead arriba y abajo de la
 *       fracción, y es el primer error que se detecta al revisar un informe.
 *
 *  Uso:
 *      ./lu_serial -n 2000 -k diagdom -s 42
 *      ./lu_serial -n 500  -k zerodiag --nopivot     (falla a propósito)
 * ===========================================================================*/

#include "lu.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* ===========================================================================
 *  EL ALGORITMO
 *
 *  Variante "right-looking" (también llamada kij, o eliminación gaussiana
 *  clásica). Se llama así porque en cada paso k se mira y se modifica TODO
 *  lo que está a la derecha y abajo del pivote.
 *
 *  Estado de la matriz durante el paso k:
 *
 *        columnas:  0 ... k-1 |  k  | k+1 ... n-1
 *        fila 0     U U U U U   U     U U U U U      ← filas ya terminadas
 *        fila 1     L U U U U   U     U U U U U
 *        ...        L L U U U   U     U U U U U
 *        fila k     L L L L L  [p]    u u u u u      ← fila pivote de este paso
 *        fila k+1   L L L L L   a     a a a a a      ┐
 *        ...        L L L L L   a     a a a a a      │ submatriz activa
 *        fila n-1   L L L L L   a     a a a a a      ┘
 *
 *  Los tres movimientos del paso k:
 *
 *    (1) BUSCAR EL PIVOTE: el elemento de mayor valor absoluto en la
 *        columna k, de la fila k para abajo.
 *    (2) INTERCAMBIAR las filas k y p (la fila completa, porque L y U
 *        comparten el mismo arreglo).
 *    (3) ELIMINAR: para cada fila i > k, calcular el multiplicador
 *        l_ik = a_ik / a_kk, guardarlo en el lugar de a_ik (ahí es donde
 *        vive L), y restarle a la fila i la fila k multiplicada por l_ik.
 *
 *  Al terminar los n pasos, el triángulo inferior estricto es L y el
 *  superior con la diagonal es U. La diagonal de L es de unos y no se guarda:
 *  ese lugar lo ocupa la diagonal de U.
 * ===========================================================================*/

int lu_factor(double *A, int n, int *piv)
{
    for (int k = 0; k < n; k++) {

        /* ---- (1) búsqueda del pivote ------------------------------------
         * Recorremos la columna k desde la fila k hacia abajo buscando el
         * máximo en VALOR ABSOLUTO. El valor absoluto importa: lo que
         * queremos es un divisor grande, no un número grande con signo. */
        int    p    = k;
        double maxv = fabs(A[(size_t)k * n + k]);
        for (int i = k + 1; i < n; i++) {
            double v = fabs(A[(size_t)i * n + k]);
            if (v > maxv) { maxv = v; p = i; }
        }
        piv[k] = p;

        /* Si toda la columna es cero, la matriz es singular: no hay LU. */
        if (maxv == 0.0) return k + 1;

        /* ---- (2) intercambio de filas -----------------------------------
         * Intercambiamos la fila ENTERA, no solo la parte de U. La parte
         * izquierda de esas filas ya contiene multiplicadores de L, y esos
         * también tienen que viajar: L corresponde a la matriz permutada. */
        if (p != k) {
            double *rk = &A[(size_t)k * n];
            double *rp = &A[(size_t)p * n];
            for (int j = 0; j < n; j++) {
                double t = rk[j]; rk[j] = rp[j]; rp[j] = t;
            }
        }

        /* ---- (3) eliminación --------------------------------------------
         * Calculamos 1/a_kk una sola vez: una división cuesta del orden de
         * 20 ciclos, una multiplicación 4. Con n−k−1 filas por paso, la
         * diferencia se nota. (Esto cambia infinitesimalmente el redondeo;
         * es la optimización estándar y LAPACK hace lo mismo.) */
        const double *rowk = &A[(size_t)k * n];
        double inv_akk = 1.0 / rowk[k];

        for (int i = k + 1; i < n; i++) {
            double *rowi = &A[(size_t)i * n];

            double l = rowi[k] * inv_akk;
            rowi[k] = l;                 /* acá se guarda L, in-place       */

            if (l == 0.0) continue;      /* fila ya alineada: nada que hacer */

            /* El lazo interno. Acá está el 99% del tiempo de ejecución del
             * programa entero, y acá es donde el compilador vectoriza
             * (SSE/AVX) si le damos -O3. Es también el trabajo que la
             * versión MPI reparte entre procesos. */
            for (int j = k + 1; j < n; j++)
                rowi[j] -= l * rowk[j];
        }
    }
    return 0;
}

/* ---------------------------------------------------------------------------
 *  Misma factorización SIN pivoteo.
 *  No se usa para producir resultados: existe para demostrar en el informe
 *  qué pasa cuando se omite el pivoteo. Con MAT_ZERODIAG divide por cero en
 *  el primer paso; con MAT_HILBERT el residuo se degrada varios órdenes.
 * -------------------------------------------------------------------------*/
int lu_factor_nopivot(double *A, int n)
{
    for (int k = 0; k < n; k++) {
        double akk = A[(size_t)k * n + k];
        if (akk == 0.0) return k + 1;              /* acá muere */

        const double *rowk = &A[(size_t)k * n];
        double inv = 1.0 / akk;

        for (int i = k + 1; i < n; i++) {
            double *rowi = &A[(size_t)i * n];
            double l = rowi[k] * inv;
            rowi[k] = l;
            for (int j = k + 1; j < n; j++)
                rowi[j] -= l * rowk[j];
        }
    }
    return 0;
}

/* ===========================================================================
 *  RESOLUCIÓN DEL SISTEMA
 *
 *  Tenemos P·A = L·U y queremos resolver A·x = b.
 *
 *      A·x = b
 *      P·A·x = P·b          (multiplicamos por P a izquierda)
 *      L·U·x = P·b
 *
 *  Llamando y = U·x, queda en dos etapas:
 *
 *      L·y = P·b     →  sustitución hacia ADELANTE  (L es triangular inferior)
 *      U·x = y       →  sustitución hacia ATRÁS     (U es triangular superior)
 *
 *  Cada etapa cuesta n² flops. La factorización costó (2/3)n³. Por eso
 *  factorizar una vez y resolver muchos lados derechos es tan conveniente:
 *  el segundo b sale prácticamente gratis.
 * ===========================================================================*/

void lu_solve(const double *LU, int n, const int *piv, double *b)
{
    /* ---- aplicar la permutación P a b -----------------------------------
     * Repetimos exactamente los mismos intercambios, en el mismo orden,
     * que hizo la factorización. */
    for (int k = 0; k < n; k++) {
        int p = piv[k];
        if (p != k) { double t = b[k]; b[k] = b[p]; b[p] = t; }
    }

    /* ---- sustitución hacia adelante:  L·y = Pb --------------------------
     * L tiene unos en la diagonal (no están guardados), así que
     *     y_i = b_i − sum_{j<i} l_ij · y_j
     * y no hace falta dividir. Escribimos y sobre b. */
    for (int i = 1; i < n; i++) {
        const double *Li = &LU[(size_t)i * n];
        double s = b[i];
        for (int j = 0; j < i; j++) s -= Li[j] * b[j];
        b[i] = s;
    }

    /* ---- sustitución hacia atrás:  U·x = y ------------------------------
     *     x_i = ( y_i − sum_{j>i} u_ij · x_j ) / u_ii
     * Acá sí dividimos por la diagonal, que es la de U. */
    for (int i = n - 1; i >= 0; i--) {
        const double *Ui = &LU[(size_t)i * n];
        double s = b[i];
        for (int j = i + 1; j < n; j++) s -= Ui[j] * b[j];
        b[i] = s / Ui[i];
    }
}

/* ===========================================================================
 *  PROGRAMA PRINCIPAL
 * ===========================================================================*/

static void usage(const char *prog)
{
    fprintf(stderr,
      "\nUso: %s [opciones]\n"
      "  -n <int>     dimension de la matriz (default 1000)\n"
      "  -k <tipo>    random | diagdom | hilbert | zerodiag  (default diagdom)\n"
      "  -s <int>     semilla (default 42)\n"
      "  -r <int>     repeticiones cronometradas (default 1)\n"
      "  --nopivot    factorizar SIN pivoteo (experimento de estabilidad)\n"
      "  --csv        imprimir una linea CSV en vez del reporte legible\n"
      "  -h           esta ayuda\n\n", prog);
}

int main(int argc, char **argv)
{
    int n = 1000, seed = 42, reps = 1, nopivot = 0, csv = 0;
    mat_kind_t kind = MAT_DIAGDOM;

    for (int i = 1; i < argc; i++) {
        if      (!strcmp(argv[i], "-n") && i + 1 < argc) n    = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-s") && i + 1 < argc) seed = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-r") && i + 1 < argc) reps = atoi(argv[++i]);
        else if (!strcmp(argv[i], "-k") && i + 1 < argc) kind = mat_kind_from_string(argv[++i]);
        else if (!strcmp(argv[i], "--nopivot")) nopivot = 1;
        else if (!strcmp(argv[i], "--csv"))     csv     = 1;
        else { usage(argv[0]); return 1; }
    }
    if (n <= 0) { usage(argv[0]); return 1; }

    size_t nn = (size_t)n * (size_t)n;

    double *A0  = (double *)malloc(nn * sizeof(double));  /* copia original  */
    double *A   = (double *)malloc(nn * sizeof(double));  /* se factoriza    */
    double *b0  = (double *)malloc((size_t)n * sizeof(double));
    double *x   = (double *)malloc((size_t)n * sizeof(double));
    int    *piv = (int    *)malloc((size_t)n * sizeof(int));

    if (!A0 || !A || !b0 || !x || !piv) {
        fprintf(stderr, "ERROR: no alcanza la memoria para n=%d "
                        "(hacen falta ~%.1f MB)\n", n,
                        2.0 * nn * sizeof(double) / 1048576.0);
        return 2;
    }

    mat_generate(A0, n, kind, (unsigned)seed);
    vec_generate(b0, n, (unsigned)seed);

    double t_best = 1e300;
    int    status = 0;

    for (int rep = 0; rep < reps; rep++) {
        /* Restaurar A y b: la factorización es destructiva. Esto queda
         * FUERA del cronómetro a propósito. */
        memcpy(A, A0, nn * sizeof(double));
        memcpy(x, b0, (size_t)n * sizeof(double));

        double t0 = wall_time();
        status = nopivot ? lu_factor_nopivot(A, n) : lu_factor(A, n, piv);
        double t1 = wall_time();

        if (status) break;

        if (nopivot) { for (int i = 0; i < n; i++) piv[i] = i; }
        lu_solve(A, n, piv, x);

        double dt = t1 - t0;                 /* solo la factorización */
        if (dt < t_best) t_best = dt;
    }

    if (status) {
        fprintf(stderr,
            "FALLO: pivote nulo en el paso k=%d.\n"
            "  Con -k zerodiag y --nopivot esto es EL resultado esperado:\n"
            "  es la demostracion de que el pivoteo no es opcional.\n", status - 1);
        return 3;
    }

    double res    = residual_rel(A0, x, b0, n);
    double growth = growth_factor(A, A0, n);
    double gflops = lu_flops(n) / t_best / 1e9;

    if (csv) {
        /* Mismo esquema de columnas que la version MPI, para que las dos
         * escriban en el mismo results/raw.csv:
         * impl,n,p,dist,kind,variant,t_total,t_comp,t_comm,t_idle,gflops,residual,growth */
        printf("serial,%d,1,-,%s,%s,%.6f,%.6f,0.000000,0.000000,%.4f,%.6e,%.4f\n",
               n, mat_kind_name(kind), nopivot ? "nopivot" : "pivot",
               t_best, t_best, gflops, res, growth);
    } else {
        printf("\n  Factorizacion LU serial\n");
        printf("  ---------------------------------------------\n");
        printf("  n                    : %d\n", n);
        printf("  matriz               : %s%s\n", mat_kind_name(kind),
                                                  nopivot ? "  (SIN pivoteo)" : "");
        printf("  memoria de A         : %.1f MB\n", nn * sizeof(double) / 1048576.0);
        printf("  repeticiones         : %d  (se reporta el mejor tiempo)\n", reps);
        printf("  ---------------------------------------------\n");
        printf("  tiempo factorizacion : %.6f s\n", t_best);
        printf("  rendimiento          : %.3f GFLOP/s\n", gflops);
        printf("  residuo relativo     : %.6e\n", res);
        printf("  factor de crecimiento: %.4f\n", growth);
        printf("  ---------------------------------------------\n");
        printf("  %s\n\n", res < 1e-11
                 ? "OK: el residuo esta en el orden esperado."
                 : "ATENCION: residuo alto. Revisar condicionamiento o pivoteo.");
    }

    free(A0); free(A); free(b0); free(x); free(piv);
    return 0;
}
