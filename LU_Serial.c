/*
 * LU_Serial.c
 * Factorizacion LU secuencial (A = L*U) por eliminacion gaussiana, sin pivoteo.
 *
 * Que hace el programa:
 *   1. Genera una matriz A de N x N (diagonalmente dominante) y un vector b
 *      elegido para que la solucion exacta de A x = b sea x = (1, 1, ..., 1).
 *   2. Factoriza A = L*U. L y U se guardan en la misma matriz:
 *        - debajo de la diagonal quedan los multiplicadores (L, con 1 en la diagonal)
 *        - en la diagonal y arriba queda U
 *   3. Resuelve L y = b (sustitucion hacia adelante) y U x = y (hacia atras).
 *   4. Verifica: residuo ||A x - b|| y error maximo |x_i - 1|.
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

/* Genera A y b. Se suma N a la diagonal para que A sea diagonalmente
   dominante: asi la eliminacion gaussiana SIN pivoteo es estable
   (nunca aparece un pivote cero o muy chico). */
void generar_matriz(double *A, double *b, int N)
{
    int i, j;
    srand(1);                       /* semilla fija: siempre la misma matriz */
    for (i = 0; i < N; i++) {
        b[i] = 0.0;
        for (j = 0; j < N; j++) {
            A[i*N + j] = (double) rand() / RAND_MAX;   /* valor entre 0 y 1 */
            if (i == j)
                A[i*N + j] += N;
            b[i] += A[i*N + j];      /* b = A * (1,1,...,1) */
        }
    }
}

/* Factorizacion LU "en el lugar" (algoritmo kij de eliminacion gaussiana).
   En el paso k se usa la fila k (fila pivote) para hacer cero la columna k
   de todas las filas que estan debajo. */
void factorizar_lu(double *A, int N)
{
    int i, j, k;
    double l;
    for (k = 0; k < N - 1; k++) {
        for (i = k + 1; i < N; i++) {
            l = A[i*N + k] / A[k*N + k];   /* multiplicador l_ik */
            A[i*N + k] = l;                /* se guarda en la parte de L */
            for (j = k + 1; j < N; j++)    /* fila_i = fila_i - l * fila_k */
                A[i*N + j] -= l * A[k*N + j];
        }
    }
}

/* Resuelve A x = b usando la matriz ya factorizada (L y U juntas en LU). */
void resolver(double *LU, double *b, double *x, int N)
{
    int i, j;
    double suma;

    /* L y = b (hacia adelante). L tiene 1 en la diagonal. y se guarda en x. */
    for (i = 0; i < N; i++) {
        suma = b[i];
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

/* Norma 2 del residuo: ||A x - b||, con la matriz A ORIGINAL. */
double residuo(double *A, double *x, double *b, int N)
{
    int i, j;
    double r, suma = 0.0;
    for (i = 0; i < N; i++) {
        r = -b[i];
        for (j = 0; j < N; j++)
            r += A[i*N + j] * x[j];
        suma += r * r;
    }
    return sqrt(suma);
}

int main(int argc, char *argv[])
{
    int N = 1000, i;
    double *A, *LU, *b, *x;
    double error, tiempo;
    clock_t t0, t1;

    if (argc > 1)
        N = atoi(argv[1]);

    A  = malloc((size_t) N * N * sizeof(double));   /* matriz original */
    LU = malloc((size_t) N * N * sizeof(double));   /* copia que se factoriza */
    b  = malloc(N * sizeof(double));
    x  = malloc(N * sizeof(double));
    if (A == NULL || LU == NULL || b == NULL || x == NULL) {
        printf("No hay memoria suficiente para N = %d\n", N);
        return 1;
    }

    generar_matriz(A, b, N);
    for (i = 0; i < N * N; i++)
        LU[i] = A[i];

    t0 = clock();
    factorizar_lu(LU, N);
    t1 = clock();
    tiempo = (double) (t1 - t0) / CLOCKS_PER_SEC;

    resolver(LU, b, x, N);

    error = 0.0;
    for (i = 0; i < N; i++)
        if (fabs(x[i] - 1.0) > error)
            error = fabs(x[i] - 1.0);

    printf("LU serial  N = %d\n", N);
    printf("Tiempo de factorizacion: %f s\n", tiempo);
    printf("Residuo ||Ax - b||:      %e\n", residuo(A, x, b, N));
    printf("Error maximo |x - 1|:    %e\n", error);

    free(A); free(LU); free(b); free(x);
    return 0;
}
