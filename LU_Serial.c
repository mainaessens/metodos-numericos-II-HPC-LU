/*
 * LU_Serial.c
 * Factorizacion LU con PIVOTEO PARCIAL (P*A = L*U), version secuencial.
 * Es la version de referencia: con sus tiempos se calcula el speedup de LU_MPI.
 *
 * Que hace el programa:
 *   1. Genera AL AZAR una matriz A de N x N y un vector b (valores uniformes
 *      en [-1, 1]). Una matriz al azar es invertible con probabilidad 1, pero
 *      NO es diagonalmente dominante: por eso hace falta pivotear.
 *   2. Factoriza P*A = L*U por eliminacion gaussiana con pivoteo parcial.
 *      L y U se guardan en la misma matriz:
 *        - debajo de la diagonal quedan los multiplicadores (L, con 1 en la diagonal)
 *        - en la diagonal y arriba queda U
 *      La permutacion P se guarda como un vector de indices, no como matriz.
 *   3. Resuelve L y = P*b (sustitucion hacia adelante) y U x = y (hacia atras).
 *   4. Mide la precision con el residuo ||A x - b|| y el residuo relativo
 *      ||A x - b|| / ||b||, usando la matriz A ORIGINAL.
 *
 * Compilar:  gcc -O2 -o LU_Serial LU_Serial.c -lm
 * Ejecutar:  ./LU_Serial 1000        (1000 = tamano N de la matriz)
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

/* La matriz se guarda en un vector de N*N doubles, fila por fila:
   el elemento (i, j) esta en A[i*N + j]. */

/* Reloj de pared, en segundos. Se usa el mismo tipo de reloj que MPI_Wtime
   para que los tiempos del serial y del paralelo sean comparables y el
   speedup tenga sentido. */
double tiempo_actual(void)
{
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return (double) t.tv_sec + (double) t.tv_nsec / 1e9;
}

/* Genera A y b AL AZAR, con valores uniformes en [-1, 1].
   La semilla es fija y el orden de generacion es el mismo que en LU_MPI.c,
   asi los dos programas factorizan exactamente la misma matriz y el mismo
   vector, y los resultados se pueden comparar. */
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

/* Factorizacion LU "en el lugar" con PIVOTEO PARCIAL.
   En el paso k:
     1. se busca en la columna k, de la fila k para abajo, el elemento de mayor
        valor absoluto, y esa fila se intercambia con la fila k (pivoteo parcial);
     2. se usa la fila k para hacer cero la columna k de las filas de abajo.
   Dividir siempre por el pivote mas grande posible mantiene todos los
   multiplicadores con |l_ik| <= 1, y eso es lo que evita que los errores de
   redondeo se amplifiquen.
   perm guarda la permutacion: la fila i de la matriz factorizada es la fila
   perm[i] de la A original. Devuelve 0 si A resulta singular, 1 si todo bien.
   Costo: aproximadamente (2/3) N^3 operaciones de punto flotante. */
int factorizar_lu(double *A, int N, int *perm)
{
    int i, j, k, p, aux_i;
    double mult, maximo, aux;

    for (i = 0; i < N; i++)
        perm[i] = i;

    for (k = 0; k < N - 1; k++) {
        /* 1. buscar el pivote: el mayor |A[i][k]| con i >= k */
        p = k;
        maximo = fabs(A[k*N + k]);
        for (i = k + 1; i < N; i++)
            if (fabs(A[i*N + k]) > maximo) {
                maximo = fabs(A[i*N + k]);
                p = i;
            }

        if (maximo == 0.0)
            return 0;                       /* columna nula: A es singular */

        /* 2. intercambiar la fila k con la fila del pivote (la fila entera:
              tambien la parte de L que ya estaba calculada) */
        if (p != k) {
            for (j = 0; j < N; j++) {
                aux        = A[k*N + j];
                A[k*N + j] = A[p*N + j];
                A[p*N + j] = aux;
            }
            aux_i = perm[k]; perm[k] = perm[p]; perm[p] = aux_i;
        }

        /* 3. eliminar: fila_i = fila_i - l_ik * fila_k */
        for (i = k + 1; i < N; i++) {
            mult = A[i*N + k] / A[k*N + k];   /* |mult| <= 1 gracias al pivoteo */
            A[i*N + k] = mult;                /* se guarda en la parte de L */
            for (j = k + 1; j < N; j++)
                A[i*N + j] -= mult * A[k*N + j];
        }
    }

    return A[(N-1)*N + (N-1)] != 0.0;
}

/* Resuelve A x = b usando la factorizacion P*A = L*U.
   Como se permutaron las filas de A, hay que aplicarle a b la misma
   permutacion: la fila i del sistema factorizado corresponde a b[perm[i]]. */
void resolver(double *LU, int *perm, double *b, double *x, int N)
{
    int i, j;
    double suma;

    /* L y = P*b (hacia adelante). L tiene 1 en la diagonal. y se guarda en x. */
    for (i = 0; i < N; i++) {
        suma = b[perm[i]];
        for (j = 0; j < i; j++)
            suma -= LU[i*N + j] * x[j];
        x[i] = suma;
    }

    /* U x = y (hacia atras). */
    for (i = N - 1; i >= 0; i--) {
        suma = x[i];
        for (j = i + 1; j < N; j++)
            suma -= LU[i*N + j] * x[j];
        x[i] = suma / LU[i*N + i];
    }
}

/* Norma 2 del residuo: ||A x - b||, con la matriz A ORIGINAL (sin permutar).
   Tambien devuelve en *rel el residuo relativo ||A x - b|| / ||b||, que es
   el que conviene comparar entre distintos N (el absoluto crece con N). */
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
    int N = 1000, i, *perm;
    double *A, *LU, *b, *x;
    double res, res_rel, t_comp;

    if (argc > 1)
        N = atoi(argv[1]);

    A    = malloc((size_t) N * N * sizeof(double));   /* matriz original */
    LU   = malloc((size_t) N * N * sizeof(double));   /* copia que se factoriza */
    b    = malloc(N * sizeof(double));
    x    = malloc(N * sizeof(double));
    perm = malloc(N * sizeof(int));
    if (A == NULL || LU == NULL || b == NULL || x == NULL || perm == NULL) {
        printf("No hay memoria suficiente para N = %d\n", N);
        return 1;
    }

    generar(A, b, N);
    for (i = 0; i < N * N; i++)
        LU[i] = A[i];

    /* Se mide solo la factorizacion: es la parte que despues se paraleliza. */
    t_comp = tiempo_actual();
    if (!factorizar_lu(LU, N, perm)) {
        printf("La matriz resulto singular, no se puede factorizar.\n");
        return 1;
    }
    t_comp = tiempo_actual() - t_comp;

    resolver(LU, perm, b, x, N);
    res = residuo(A, x, b, N, &res_rel);

    printf("LU serial (con pivoteo parcial)  N = %d\n", N);
    printf("Tiempo de computo:        %f s\n", t_comp);
    printf("Tiempo de comunicacion:   0.000000 s  (no hay: un solo proceso)\n");
    printf("Tiempo total:             %f s\n", t_comp);
    printf("Residuo ||Ax - b||:       %e\n", res);
    printf("Residuo relativo:         %e\n", res_rel);
    /* Linea lista para copiar a una planilla (N, P, computo, comunicacion,
       total, residuo, residuo relativo): */
    printf("CSV,%d,1,%f,0.000000,%f,%e,%e\n", N, t_comp, t_comp, res, res_rel);

    free(A); free(LU); free(b); free(x); free(perm);
    return 0;
}
