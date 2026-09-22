/* ===========================================================================
 *  verify.c — Normas, residuo, factor de crecimiento y utilidades de tiempo
 *
 *  Todo lo que responde a la pregunta "¿el resultado está bien?" vive acá.
 * ===========================================================================*/

/* clock_gettime() y CLOCK_MONOTONIC son POSIX, no C estándar. Con -std=c11
 * el compilador los esconde; este define los vuelve a habilitar sin tener
 * que cambiar las banderas del Makefile. Va ANTES de cualquier #include. */
#define _POSIX_C_SOURCE 199309L

#include "lu.h"
#include <math.h>
#include <time.h>
#include <stdlib.h>

/* --------------------------------------------------------------------------
 *  Norma infinito de un vector:  ||v||∞ = max_i |v_i|
 *
 *  Usamos la norma infinito (y no la euclídea) porque es la que aparece en
 *  las cotas clásicas de error del análisis backward, y porque no sufre
 *  overflow al elevar al cuadrado.
 * ------------------------------------------------------------------------*/
double norm_inf_vec(const double *v, int n)
{
    double m = 0.0;
    for (int i = 0; i < n; i++) {
        double a = fabs(v[i]);
        if (a > m) m = a;
    }
    return m;
}

/* --------------------------------------------------------------------------
 *  Norma infinito de una matriz: el máximo de las sumas de filas en módulo.
 *      ||A||∞ = max_i ( sum_j |a_ij| )
 *  Es la norma inducida por la norma infinito de vectores: ||Ax||∞ ≤ ||A||∞·||x||∞
 * ------------------------------------------------------------------------*/
double norm_inf_mat(const double *A, int n)
{
    double m = 0.0;
    for (int i = 0; i < n; i++) {
        double s = 0.0;
        for (int j = 0; j < n; j++) s += fabs(A[(size_t)i * n + j]);
        if (s > m) m = s;
    }
    return m;
}

/* y = A·x  (producto matriz-vector, n² flops) */
void matvec(const double *A, const double *x, double *y, int n)
{
    for (int i = 0; i < n; i++) {
        double s = 0.0;
        const double *Ai = &A[(size_t)i * n];
        for (int j = 0; j < n; j++) s += Ai[j] * x[j];
        y[i] = s;
    }
}

/* --------------------------------------------------------------------------
 *  RESIDUO RELATIVO — la métrica de precisión del proyecto
 *
 *              ||A·x − b||∞
 *      r = ---------------------
 *           ||A||∞ · ||x||∞
 *
 *  Interpretación (esto hay que saberlo decir en la defensa):
 *
 *  r pequeño NO significa "x es la solución exacta". Significa que x es la
 *  solución EXACTA de un sistema perturbado (A+ΔA)x = b con ||ΔA|| ≈ r·||A||.
 *  Eso se llama estabilidad BACKWARD: el algoritmo no introdujo más error
 *  que el que ya tendría el problema si los datos de entrada estuvieran
 *  medidos con precisión de máquina.
 *
 *  El error FORWARD (cuán lejos está x de la solución verdadera) se acota por
 *
 *      ||x − x_exacta|| / ||x_exacta||  ≲  κ(A) · r
 *
 *  donde κ(A) es el número de condición. Por eso con una matriz de Hilbert,
 *  donde κ(A) puede valer 1e10, un residuo de 1e-16 todavía permite un error
 *  relativo de 1e-6 en la solución. Residuo y error NO son lo mismo.
 *
 *  Normalizamos dividiendo por ||A||·||x|| para que el número sea adimensional
 *  y comparable entre distintos n: así lo esperable es siempre O(ε) ≈ 1e-16,
 *  sin importar la escala de la matriz.
 * ------------------------------------------------------------------------*/
double residual_rel(const double *A, const double *x, const double *b, int n)
{
    double *r = (double *)malloc((size_t)n * sizeof(double));
    if (!r) return -1.0;

    matvec(A, x, r, n);                       /* r = A·x            */
    for (int i = 0; i < n; i++) r[i] -= b[i]; /* r = A·x − b        */

    double num = norm_inf_vec(r, n);
    double den = norm_inf_mat(A, n) * norm_inf_vec(x, n);

    free(r);
    return (den > 0.0) ? num / den : num;
}

/* --------------------------------------------------------------------------
 *  FACTOR DE CRECIMIENTO
 *
 *      ρ = max_ij |u_ij| / max_ij |a_ij|
 *
 *  Es la cantidad que controla la estabilidad de la eliminación gaussiana.
 *  La cota clásica del error backward es proporcional a ρ.
 *
 *  - Sin pivoteo, ρ puede ser arbitrariamente grande (el algoritmo es inestable).
 *  - Con pivoteo parcial, la cota teórica es 2^(n−1) — que es horrible — pero
 *    en la práctica ρ casi nunca pasa de unas pocas decenas. Ese hueco entre
 *    la teoría y la práctica es un tema clásico de análisis numérico, y es la
 *    razón por la que el pivoteo parcial se usa en todos lados a pesar de no
 *    tener una garantía teórica buena.
 *
 *  LU es el arreglo in-place devuelto por lu_factor(); A0 es la matriz original.
 * ------------------------------------------------------------------------*/
double growth_factor(const double *LU, const double *A0, int n)
{
    double maxU = 0.0, maxA = 0.0;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double a = fabs(A0[(size_t)i * n + j]);
            if (a > maxA) maxA = a;
            if (j >= i) {                       /* solo el triángulo U */
                double u = fabs(LU[(size_t)i * n + j]);
                if (u > maxU) maxU = u;
            }
        }
    }
    return (maxA > 0.0) ? maxU / maxA : 0.0;
}

/* --------------------------------------------------------------------------
 *  Reloj monótono de pared, en segundos.
 *
 *  Usamos CLOCK_MONOTONIC y no clock() porque clock() mide tiempo de CPU
 *  (que con varios hilos o procesos no es lo que queremos) y porque
 *  CLOCK_MONOTONIC no salta si alguien cambia la hora del sistema.
 *
 *  En la versión MPI usamos MPI_Wtime(), que es el equivalente portable.
 * ------------------------------------------------------------------------*/
double wall_time(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec + 1e-9 * (double)ts.tv_nsec;
}

/* --------------------------------------------------------------------------
 *  Cantidad de operaciones de punto flotante de la factorización LU.
 *
 *  En el paso k se actualiza una submatriz de (n−k−1)×(n−k−1), y cada
 *  elemento cuesta una multiplicación y una resta (2 flops), más (n−k−1)
 *  divisiones para los multiplicadores. Sumando sobre k:
 *
 *      sum_{k=0}^{n-1} 2(n−k−1)²  =  2 · sum_{m=1}^{n-1} m²  ≈  (2/3)n³
 *
 *  El término exacto es (2/3)n³ − (1/2)n² − (1/6)n. Para reportar GFLOP/s
 *  se usa este conteo dividido por el tiempo medido.
 * ------------------------------------------------------------------------*/
double lu_flops(int n)
{
    double N = (double)n;
    return (2.0 / 3.0) * N * N * N - 0.5 * N * N - (1.0 / 6.0) * N;
}
