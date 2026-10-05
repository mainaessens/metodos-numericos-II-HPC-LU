/*
 * LU_MPI.c
 * Factorizacion LU paralela (A = L*U) con MPI, sin pivoteo.
 *
 * Idea de la paralelizacion (distribucion CICLICA por filas):
 *   - Las filas de A se reparten entre los P procesos de forma ciclica:
 *       fila 0 -> proceso 0, fila 1 -> proceso 1, ..., fila P -> proceso 0, ...
 *     La fila global i la tiene el proceso (i % P), como su fila local (i / P).
 *   - En cada paso k de la eliminacion gaussiana:
 *       1. El duenio de la fila k (proceso k % P) la envia a todos (MPI_Bcast).
 *       2. Cada proceso actualiza SUS filas que estan debajo de k.
 *   - Se usa reparto ciclico (y no por bloques) para balancear la carga:
 *     a medida que k avanza, las filas de arriba ya no trabajan; con el reparto
 *     ciclico todos los procesos siguen teniendo filas por debajo de k.
 *
 * El programa mide por separado el tiempo de COMPUTO y el de COMUNICACION,
 * para poder analizar como escala cada uno al aumentar N y la cantidad de procesos.
 *
 * Funciones MPI usadas (todas de la teoria / TP1):
 *   MPI_Init, MPI_Comm_rank, MPI_Comm_size, MPI_Finalize, MPI_Wtime,
 *   MPI_Barrier, MPI_Scatter, MPI_Bcast, MPI_Gather, MPI_Reduce.
 *
 * Compilar:  mpicc -O2 -o LU_MPI LU_MPI.c -lm
 * Ejecutar:  mpirun -np 4 ./LU_MPI 1000     (N debe ser multiplo de la cantidad de procesos)
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

/* Igual que en LU_Serial.c: matriz guardada fila por fila, (i, j) -> A[i*N + j]. */

/* Genera A (diagonalmente dominante, por lo tanto invertible) y b = A*(1,...,1).
   Misma semilla que LU_Serial.c, asi los dos programas factorizan exactamente
   la misma matriz y los resultados se pueden comparar. */
void generar_matriz(double *A, double *b, int N)
{
    int i, j;
    srand(1);
    for (i = 0; i < N; i++) {
        b[i] = 0.0;
        for (j = 0; j < N; j++) {
            A[i*N + j] = (double) rand() / RAND_MAX;
            if (i == j)
                A[i*N + j] += N;
            b[i] += A[i*N + j];
        }
    }
}

/* Resuelve A x = b con la matriz ya factorizada (igual que en LU_Serial.c). */
void resolver(double *LU, double *b, double *x, int N)
{
    int i, j;
    double suma;
    for (i = 0; i < N; i++) {                /* L y = b */
        suma = b[i];
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
    int i, j, k, l, p, duenio;
    double *A = NULL, *buf = NULL, *LU = NULL, *b = NULL, *x = NULL;  /* solo proceso 0 */
    double *Aloc, *filak, *fila;                                      /* todos los procesos */
    double mult, error, res, res_rel;
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

    Aloc  = malloc((size_t) nloc * N * sizeof(double));  /* mis filas */
    filak = malloc(N * sizeof(double));                  /* para recibir la fila pivote */

    /* ---------- Proceso 0: genera A y la ordena para repartir ---------- */
    if (rank == 0) {
        A   = malloc((size_t) N * N * sizeof(double));   /* original (para verificar) */
        buf = malloc((size_t) N * N * sizeof(double));   /* filas ordenadas por proceso */
        LU  = malloc((size_t) N * N * sizeof(double));   /* resultado final */
        b   = malloc(N * sizeof(double));
        x   = malloc(N * sizeof(double));
        generar_matriz(A, b, N);

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

    /* ---------- Factorizacion LU en paralelo ---------- */
    for (k = 0; k < N - 1; k++) {
        duenio = k % P;                       /* proceso que tiene la fila k */

        /* El duenio envia su fila k directamente desde Aloc; los demas la reciben en filak.
           Solo se envian las columnas k..N-1 (las anteriores ya no se usan). */
        if (rank == duenio)
            fila = &Aloc[(k / P) * N];
        else
            fila = filak;

        /* COMUNICACION: difusion de la fila pivote. */
        t = MPI_Wtime();
        MPI_Bcast(&fila[k], N - k, MPI_DOUBLE, duenio, MPI_COMM_WORLD);
        t_comm += MPI_Wtime() - t;

        /* COMPUTO: cada proceso actualiza sus filas que estan DEBAJO de la fila k. */
        t = MPI_Wtime();
        for (l = 0; l < nloc; l++) {
            i = l * P + rank;                 /* numero de fila global */
            if (i <= k)
                continue;                     /* esta fila ya esta terminada */
            mult = Aloc[l*N + k] / fila[k];   /* multiplicador l_ik */
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

    /* El tiempo del programa es el del proceso mas lento: se toma el maximo. */
    MPI_Reduce(&t_comp,  &t_comp_max,  1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&t_comm,  &t_comm_max,  1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
    MPI_Reduce(&t_total, &t_total_max, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);

    /* ---------- Proceso 0: reordenar, resolver y verificar ---------- */
    if (rank == 0) {
        /* Deshacer el orden del Scatter: la fila local l del proceso p es la fila global l*P + p. */
        for (p = 0; p < P; p++)
            for (l = 0; l < nloc; l++)
                for (j = 0; j < N; j++)
                    LU[(l*P + p)*N + j] = buf[(p*nloc + l)*N + j];

        resolver(LU, b, x, N);

        error = 0.0;
        for (i = 0; i < N; i++)
            if (fabs(x[i] - 1.0) > error)
                error = fabs(x[i] - 1.0);
        res = residuo(A, x, b, N, &res_rel);

        printf("LU MPI  N = %d  procesos = %d\n", N, P);
        printf("Tiempo de computo:        %f s\n", t_comp_max);
        printf("Tiempo de comunicacion:   %f s\n", t_comm_max);
        printf("Tiempo total:             %f s\n", t_total_max);
        printf("Residuo ||Ax - b||:       %e\n", res);
        printf("Residuo relativo:         %e\n", res_rel);
        printf("Error maximo |x - 1|:     %e\n", error);
        /* Linea lista para copiar a una planilla (N, P, computo, comunicacion,
           total, residuo, residuo relativo): */
        printf("CSV,%d,%d,%f,%f,%f,%e,%e\n",
               N, P, t_comp_max, t_comm_max, t_total_max, res, res_rel);

        free(A); free(buf); free(LU); free(b); free(x);
    }

    free(Aloc); free(filak);
    MPI_Finalize();
    return 0;
}
