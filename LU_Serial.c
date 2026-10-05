/*
 * LU_Serial.c
 * Factorizacion LU secuencial (A = L*U) por eliminacion gaussiana, sin pivoteo.
 * Es la version de referencia: con sus tiempos se calcula el speedup de LU_MPI.
 *
 * Que hace el programa:
 *   1. Genera una matriz A de N x N (diagonalmente dominante, por lo tanto
 *      invertible) y un vector b elegido para que la solucion exacta de
 *      A x = b sea x = (1, 1, ..., 1).
 *   2. Factoriza A = L*U. L y U se guardan en la misma matriz:
 *        - debajo de la diagonal quedan los multiplicadores (L, con 1 en la diagonal)
 *        - en la diagonal y arriba queda U
 *   3. Resuelve L y = b (sustitucion hacia adelante) y U x = y (hacia atras).
 *   4. Verifica la precision: residuo ||A x - b||, residuo relativo
 *      ||A x - b|| / ||b|| y error maximo |x_i - 1|.
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

/* Genera A y b. Se suma N a la diagonal para que A sea diagonalmente
   dominante: asi A es invertible y la eliminacion gaussiana SIN pivoteo es
   estable (nunca aparece un pivote cero o muy chico). */
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
   de todas las filas que estan debajo.
   Costo: aproximadamente (2/3) N^3 operaciones de punto flotante. */
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

/* Norma 2 del residuo: ||A x - b||, con la matriz A ORIGINAL.
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
    int N = 1000, i;
    double *A, *LU, *b, *x;
    double error, res, res_rel, t_comp;

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

    /* Se mide solo la factorizacion: es la parte que despues se paraleliza. */
    t_comp = tiempo_actual();
    factorizar_lu(LU, N);
    t_comp = tiempo_actual() - t_comp;

    resolver(LU, b, x, N);

    error = 0.0;
    for (i = 0; i < N; i++)
        if (fabs(x[i] - 1.0) > error)
            error = fabs(x[i] - 1.0);
    res = residuo(A, x, b, N, &res_rel);

    printf("LU serial  N = %d\n", N);
    printf("Tiempo de computo:        %f s\n", t_comp);
    printf("Tiempo de comunicacion:   0.000000 s  (no hay: un solo proceso)\n");
    printf("Tiempo total:             %f s\n", t_comp);
    printf("Residuo ||Ax - b||:       %e\n", res);
    printf("Residuo relativo:         %e\n", res_rel);
    printf("Error maximo |x - 1|:     %e\n", error);
    /* Linea lista para copiar a una planilla (N, P, computo, comunicacion,
       total, residuo, residuo relativo): */
    printf("CSV,%d,1,%f,0.000000,%f,%e,%e\n", N, t_comp, t_comp, res, res_rel);

    free(A); free(LU); free(b); free(x);
    return 0;
}
