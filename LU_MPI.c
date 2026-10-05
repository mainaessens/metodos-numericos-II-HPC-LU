/*
 * LU_MPI.c
 * Factorizacion LU paralela con PIVOTEO PARCIAL (P*A = L*U) usando MPI.
 *
 * Idea de la paralelizacion (distribucion CICLICA por filas):
 *   - Las filas de A se reparten entre los P procesos de forma ciclica:
 *       fila 0 -> proceso 0, fila 1 -> proceso 1, ..., fila P -> proceso 0, ...
 *     La fila global i la tiene el proceso (i % P), como su fila local (i / P).
 *   - Se usa reparto ciclico (y no por bloques) para balancear la carga:
 *     a medida que k avanza, las filas de arriba ya no trabajan; con el reparto
 *     ciclico todos los procesos siguen teniendo filas por debajo de k.
 *
 * En cada paso k de la eliminacion gaussiana hacen falta cuatro cosas:
 *   1. PIVOTEO (busqueda local): cada proceso busca, entre SUS filas i >= k,
 *      la de mayor |a[i][k]|.
 *   2. PIVOTEO (eleccion global): el proceso 0 junta los P candidatos con
 *      MPI_Gather, elige el mayor y avisa a todos cual es con MPI_Bcast.
 *   3. INTERCAMBIO: la fila k y la fila del pivote se intercambian. Si estan en
 *      procesos distintos, sus duenios se las mandan con MPI_Send / MPI_Recv.
 *   4. ELIMINACION: el duenio de la fila k la difunde con MPI_Bcast y cada
 *      proceso actualiza sus propias filas que estan debajo de k.
 *
 * El programa mide por separado el tiempo de COMPUTO y el de COMUNICACION,
 * para poder analizar como escala cada uno al aumentar N y la cantidad de procesos.
 *
 * Funciones MPI usadas (todas vistas en la teoria / TP1):
 *   MPI_Init, MPI_Comm_rank, MPI_Comm_size, MPI_Finalize, MPI_Wtime,
 *   MPI_Barrier, MPI_Scatter, MPI_Bcast, MPI_Gather, MPI_Send, MPI_Recv.
 *
 * Compilar:  mpicc -O2 -o LU_MPI LU_MPI.c -lm
 * Ejecutar:  mpirun -np 4 ./LU_MPI 1000   (N debe ser multiplo de la cantidad de procesos)
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

/* Igual que en LU_Serial.c: matriz guardada fila por fila, (i, j) -> A[i*N + j]. */

/* Genera A y b AL AZAR, con valores uniformes en [-1, 1].
   Es exactamente la misma funcion que en LU_Serial.c (misma semilla y mismo
   orden de generacion), asi los dos programas resuelven el mismo sistema y los
   resultados se pueden comparar. */
void generar(double *A, double *b, int N)
{
    int i, j;
    srand(1);
    for (i = 0; i < N; i++)
        for (j = 0; j < N; j++)
            A[i*N + j] = 2.0 * rand() / RAND_MAX - 1.0;
    for (i = 0; i < N; i++)
        b[i] = 2.0 * rand() / RAND_MAX - 1.0;
}

/* Resuelve A x = b usando la factorizacion P*A = L*U (igual que en LU_Serial.c).
   Como se permutaron las filas de A, a b se le aplica la misma permutacion. */
void resolver(double *LU, int *perm, double *b, double *x, int N)
{
    int i, j;
    double suma;
    for (i = 0; i < N; i++) {                /* L y = P*b */
        suma = b[perm[i]];
        for (j = 0; j < i; j++)
            suma -= LU[i*N + j] * x[j];
        x[i] = suma;
    }
    for (i = N - 1; i >= 0; i--) {           /* U x = y */
        suma = x[i];
        for (j = i + 1; j < N; j++)
            suma -= LU[i*N + j] * x[j];
        x[i] = suma / LU[i*N + i];
    }
}

/* Norma 2 del residuo ||A x - b|| con la matriz A original, y en *rel el
   residuo relativo ||A x - b|| / ||b||. */
double residuo(double *A, double *x, double *b, int N, double *rel)
{
    int i, j;
    double r, suma = 0.0, normab = 0.0;
    for (i = 0; i < N; i++) {
        r = -b[i];
        for (j = 0; j < N; j++)
            r += A[i*N + j] * x[j];
        suma += r * r;
        normab += b[i] * b[i];
    }
    suma = sqrt(suma);
    *rel = suma / sqrt(normab);
    return suma;
}

int main(int argc, char *argv[])
{
    int rank, P, N = 1000;
    int nloc;                 /* cantidad de filas que tiene cada proceso = N / P */
    int i, j, k, l, p, duenio, pivote, dk, dp, otro, aux_i, singular = 0;
    int *perm = NULL;                                                 /* solo proceso 0 */
    double *A = NULL, *buf = NULL, *LU = NULL, *b = NULL, *x = NULL;  /* solo proceso 0 */
    double *cand = NULL, *tiempos = NULL;                             /* solo proceso 0 */
    double *Aloc, *filak, *fila_aux, *fila, *mia;                     /* todos los procesos */
    double local[2], mios[3];
    double mult, maximo, aux, res, res_rel;
    double t, t_ini;
    double t_comp = 0.0, t_comm = 0.0, t_total;      /* tiempos de ESTE proceso */
    double t_comp_max, t_comm_max, t_total_max;      /* del proceso mas lento */

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &P);

    if (argc > 1)
        N = atoi(argv[1]);

    /* Todos conocen N y P, asi que todos pueden hacer este control y salir juntos. */
    if (N % P != 0) {
        if (rank == 0)
            printf("Error: N = %d debe ser multiplo de la cantidad de procesos (%d)\n", N, P);
        MPI_Finalize();
        return 1;
    }
    nloc = N / P;

    Aloc     = malloc((size_t) nloc * N * sizeof(double));  /* mis filas */
    filak    = malloc(N * sizeof(double));                  /* para recibir la fila pivote */
    fila_aux = malloc(N * sizeof(double));                  /* para intercambiar filas */

    /* ---------- Proceso 0: genera A y b, y ordena A para repartirla ---------- */
    if (rank == 0) {
        A       = malloc((size_t) N * N * sizeof(double));  /* original (para verificar) */
        buf     = malloc((size_t) N * N * sizeof(double));  /* filas ordenadas por proceso */
        LU      = malloc((size_t) N * N * sizeof(double));  /* resultado final */
        b       = malloc(N * sizeof(double));
        x       = malloc(N * sizeof(double));
        perm    = malloc(N * sizeof(int));
        cand    = malloc(2 * P * sizeof(double));           /* candidatos a pivote */
        tiempos = malloc(3 * P * sizeof(double));           /* tiempos de cada proceso */
        generar(A, b, N);

        for (i = 0; i < N; i++)
            perm[i] = i;

        /* MPI_Scatter manda bloques CONTIGUOS: el 1er bloque al proceso 0, el 2do al 1, ...
           Entonces se copian las filas en buf en ese orden:
           primero las del proceso 0 (filas 0, P, 2P, ...), despues las del 1 (1, P+1, ...), etc. */
        for (p = 0; p < P; p++)
            for (l = 0; l < nloc; l++)
                for (j = 0; j < N; j++)
                    buf[(p*nloc + l)*N + j] = A[(l*P + p)*N + j];
    }

    /* Todos arrancan a medir juntos. */
    MPI_Barrier(MPI_COMM_WORLD);
    t_ini = MPI_Wtime();

    /* ---------- COMUNICACION: reparto inicial de las filas ---------- */
    t = MPI_Wtime();
    MPI_Scatter(buf, nloc * N, MPI_DOUBLE, Aloc, nloc * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    t_comm += MPI_Wtime() - t;

    /* ---------- Factorizacion LU en paralelo, con pivoteo parcial ---------- */
    for (k = 0; k < N - 1; k++) {

        /* --- 1. COMPUTO: cada proceso busca su mejor candidato a pivote ---
           Solo mira SUS filas que estan de la k para abajo. */
        t = MPI_Wtime();
        local[0] = 0.0;      /* el mayor |a[i][k]| que tiene este proceso */
        local[1] = -1.0;     /* en que fila global esta */
        for (l = 0; l < nloc; l++) {
            i = l * P + rank;                  /* numero de fila global */
            if (i < k)
                continue;                      /* esta fila ya esta terminada */
            if (fabs(Aloc[l*N + k]) > local[0]) {
                local[0] = fabs(Aloc[l*N + k]);
                local[1] = (double) i;
            }
        }
        t_comp += MPI_Wtime() - t;

        /* --- 2. COMUNICACION: el proceso 0 junta los P candidatos, elige el
               mayor de todos y le avisa a todos cual es la fila pivote. --- */
        t = MPI_Wtime();
        MPI_Gather(local, 2, MPI_DOUBLE, cand, 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        if (rank == 0) {
            maximo = 0.0;
            pivote = -1;
            for (p = 0; p < P; p++)
                if (cand[2*p] > maximo) {
                    maximo = cand[2*p];
                    pivote = (int) cand[2*p + 1];
                }
        }
        MPI_Bcast(&pivote, 1, MPI_INT, 0, MPI_COMM_WORLD);
        t_comm += MPI_Wtime() - t;

        /* Toda la columna k vale cero: A es singular. Todos reciben pivote = -1,
           asi que todos cortan el bucle juntos. */
        if (pivote < 0) {
            singular = 1;
            break;
        }

        /* --- 3. COMUNICACION: intercambiar la fila k con la fila pivote ---
           Se intercambia la fila ENTERA (tambien la parte de L ya calculada). */
        t = MPI_Wtime();
        if (pivote != k) {
            dk = k % P;             /* duenio de la fila k */
            dp = pivote % P;        /* duenio de la fila pivote */

            if (dk == dp) {
                /* las dos filas estan en el mismo proceso: intercambio local */
                if (rank == dk)
                    for (j = 0; j < N; j++) {
                        aux                      = Aloc[(k/P)*N + j];
                        Aloc[(k/P)*N + j]        = Aloc[(pivote/P)*N + j];
                        Aloc[(pivote/P)*N + j]   = aux;
                    }
            } else if (rank == dk || rank == dp) {
                /* estan en procesos distintos: los dos duenios se mandan su fila.
                   El de rank menor manda primero y despues recibe; el otro al
                   reves. Asi nunca quedan los dos esperando a la vez (deadlock). */
                otro = (rank == dk) ? dp : dk;
                mia  = (rank == dk) ? &Aloc[(k/P)*N] : &Aloc[(pivote/P)*N];

                if (rank < otro) {
                    MPI_Send(mia, N, MPI_DOUBLE, otro, 0, MPI_COMM_WORLD);
                    MPI_Recv(fila_aux, N, MPI_DOUBLE, otro, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                } else {
                    MPI_Recv(fila_aux, N, MPI_DOUBLE, otro, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Send(mia, N, MPI_DOUBLE, otro, 0, MPI_COMM_WORLD);
                }
                for (j = 0; j < N; j++)
                    mia[j] = fila_aux[j];
            }

            /* El proceso 0 anota la permutacion para despues permutar b. */
            if (rank == 0) {
                aux_i = perm[k]; perm[k] = perm[pivote]; perm[pivote] = aux_i;
            }
        }
        t_comm += MPI_Wtime() - t;

        /* --- 4. COMUNICACION: el duenio de la fila k la difunde a todos ---
           Solo se mandan las columnas k..N-1 (las anteriores ya no se usan).
           El duenio la manda directamente desde Aloc; los demas la reciben en filak. */
        duenio = k % P;
        fila = (rank == duenio) ? &Aloc[(k/P)*N] : filak;

        t = MPI_Wtime();
        MPI_Bcast(&fila[k], N - k, MPI_DOUBLE, duenio, MPI_COMM_WORLD);
        t_comm += MPI_Wtime() - t;

        /* --- 5. COMPUTO: cada proceso actualiza SUS filas que estan debajo de k --- */
        t = MPI_Wtime();
        for (l = 0; l < nloc; l++) {
            i = l * P + rank;                 /* numero de fila global */
            if (i <= k)
                continue;                     /* esta fila ya esta terminada */
            mult = Aloc[l*N + k] / fila[k];   /* multiplicador l_ik, |mult| <= 1 */
            Aloc[l*N + k] = mult;             /* se guarda en la parte de L */
            for (j = k + 1; j < N; j++)
                Aloc[l*N + j] -= mult * fila[j];
        }
        t_comp += MPI_Wtime() - t;
    }

    /* ---------- COMUNICACION: juntar el resultado en el proceso 0 ---------- */
    t = MPI_Wtime();
    MPI_Gather(Aloc, nloc * N, MPI_DOUBLE, buf, nloc * N, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    t_comm += MPI_Wtime() - t;

    t_total = MPI_Wtime() - t_ini;

    /* El tiempo del programa es el del proceso mas lento: el proceso 0 junta los
       tres tiempos de cada proceso y se queda con el maximo de cada uno. */
    mios[0] = t_comp; mios[1] = t_comm; mios[2] = t_total;
    MPI_Gather(mios, 3, MPI_DOUBLE, tiempos, 3, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    /* ---------- Proceso 0: reordenar, resolver y verificar ---------- */
    if (rank == 0) {
        t_comp_max = t_comm_max = t_total_max = 0.0;
        for (p = 0; p < P; p++) {
            if (tiempos[3*p]     > t_comp_max)  t_comp_max  = tiempos[3*p];
            if (tiempos[3*p + 1] > t_comm_max)  t_comm_max  = tiempos[3*p + 1];
            if (tiempos[3*p + 2] > t_total_max) t_total_max = tiempos[3*p + 2];
        }

        if (singular) {
            printf("La matriz resulto singular, no se puede factorizar.\n");
        } else {
            /* Deshacer el orden del Scatter: la fila local l del proceso p es la fila global l*P + p. */
            for (p = 0; p < P; p++)
                for (l = 0; l < nloc; l++)
                    for (j = 0; j < N; j++)
                        LU[(l*P + p)*N + j] = buf[(p*nloc + l)*N + j];

            resolver(LU, perm, b, x, N);
            res = residuo(A, x, b, N, &res_rel);

            printf("LU MPI (con pivoteo parcial)  N = %d  procesos = %d\n", N, P);
            printf("Tiempo de computo:        %f s\n", t_comp_max);
            printf("Tiempo de comunicacion:   %f s\n", t_comm_max);
            printf("Tiempo total:             %f s\n", t_total_max);
            printf("Residuo ||Ax - b||:       %e\n", res);
            printf("Residuo relativo:         %e\n", res_rel);
            /* Linea lista para copiar a una planilla (N, P, computo, comunicacion,
               total, residuo, residuo relativo): */
            printf("CSV,%d,%d,%f,%f,%f,%e,%e\n",
                   N, P, t_comp_max, t_comm_max, t_total_max, res, res_rel);
        }

        free(A); free(buf); free(LU); free(b); free(x);
        free(perm); free(cand); free(tiempos);
    }

    free(Aloc); free(filak); free(fila_aux);
    MPI_Finalize();
    return 0;
}
