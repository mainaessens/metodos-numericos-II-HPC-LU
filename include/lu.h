/* ===========================================================================
 *  lu.h — Declaraciones comunes del proyecto de Factorización LU Paralela
 *
 *  Proyecto 6 · Métodos Numéricos II y Computación Científica
 *  FACET — Universidad Nacional de Tucumán — 2026
 *
 *  Este header lo incluyen TANTO la versión serial COMO la versión MPI.
 *  Todo lo que no depende de MPI vive acá.
 *
 *  Convención de almacenamiento que usamos en todo el proyecto:
 *      Las matrices son arreglos 1D contiguos en orden "row-major".
 *      El elemento de la fila i y columna j de una matriz de n columnas es
 *
 *              A[i*n + j]
 *
 *      ¿Por qué 1D y no double** ? Por dos razones:
 *        1. MPI necesita buffers CONTIGUOS en memoria para mandar datos.
 *           Con double** cada fila puede estar en cualquier lado del heap.
 *        2. Localidad espacial: recorrer una fila recorre memoria consecutiva,
 *           que es exactamente lo que la cache line quiere.
 * ===========================================================================*/

#ifndef LU_H
#define LU_H

/* --------------------------------------------------------------------------
 *  1. GENERACIÓN DE MATRICES DE PRUEBA
 * ------------------------------------------------------------------------*/

typedef enum {
    MAT_RANDOM,      /* Aleatoria uniforme en [-1, 1). Condicionamiento medio. */
    MAT_DIAGDOM,     /* Diagonal estrictamente dominante: bien condicionada.
                        Es el caso "fácil": LU es estable incluso sin pivoteo. */
    MAT_HILBERT,     /* Hilbert: a_ij = 1/(i+j+1). Célebremente MAL
                        condicionada: κ(A) crece exponencialmente con n.
                        Sirve para mostrar que residuo chico ≠ solución exacta. */
    MAT_ZERODIAG     /* Aleatoria pero con TODA la diagonal en cero.
                        Sin pivoteo el algoritmo divide por cero en el paso 0.
                        Con pivoteo parcial se resuelve sin problemas.
                        Es el experimento que justifica el pivoteo. */
} mat_kind_t;

mat_kind_t  mat_kind_from_string(const char *s);
const char *mat_kind_name(mat_kind_t k);

/* Llena A (n x n, row-major) con la matriz del tipo pedido. */
void mat_generate(double *A, int n, mat_kind_t kind, unsigned seed);

/* Genera SOLO la fila global i (n elementos) de esa misma matriz.
 *
 * Esta función es la clave de la versión MPI: permite que cada proceso
 * construya únicamente las filas que le tocan, sin que nadie tenga que
 * alocar la matriz completa de n*n doubles. Para n = 8000 eso es la
 * diferencia entre 512 MB por proceso y 512/p MB por proceso.
 *
 * Garantía importante: la fila que devuelve es BIT A BIT la misma que la
 * fila i de mat_generate() con la misma semilla. Por eso podemos comparar
 * el resultado serial contra el paralelo. */
void mat_generate_row(double *row, int n, int i, mat_kind_t kind, unsigned seed);

/* Vector b de prueba. También determinista a partir de la semilla. */
void vec_generate(double *b, int n, unsigned seed);

/* --------------------------------------------------------------------------
 *  2. FACTORIZACIÓN LU SERIAL
 * ------------------------------------------------------------------------*/

/* Factoriza A in-place con pivoteo parcial: P·A = L·U
 *
 *   - Al volver, el triángulo estrictamente inferior de A contiene L
 *     (sin su diagonal, que es de unos y no hace falta guardar) y el
 *     triángulo superior incluyendo la diagonal contiene U.
 *   - piv[] tiene n elementos. piv[k] = índice de la fila con la que se
 *     intercambió la fila k en el paso k. Así se representa P sin guardar
 *     una matriz de permutación de n*n.
 *
 * Devuelve 0 si salió bien, o k+1 si en el paso k toda la columna era cero
 * (matriz singular). */
int lu_factor(double *A, int n, int *piv);

/* Idéntica pero SIN pivoteo. Existe solo para el experimento de estabilidad:
 * con MAT_ZERODIAG falla, y con MAT_HILBERT pierde precisión. */
int lu_factor_nopivot(double *A, int n);

/* Resuelve el sistema usando la factorización ya calculada.
 * b entra como lado derecho y sale sobrescrito con la solución x. */
void lu_solve(const double *LU, int n, const int *piv, double *b);

/* --------------------------------------------------------------------------
 *  3. VERIFICACIÓN NUMÉRICA
 * ------------------------------------------------------------------------*/

double norm_inf_vec(const double *v, int n);          /* max |v_i|            */
double norm_inf_mat(const double *A, int n);          /* max_i sum_j |a_ij|   */

/* Residuo relativo:   ||A·x − b||∞ / (||A||∞ · ||x||∞)
 *
 * Esta es LA métrica de precisión del proyecto. Un valor del orden de
 * 1e-16 (el epsilon de la máquina en doble precisión) significa que el
 * algoritmo se comportó tan bien como es posible en aritmética finita. */
double residual_rel(const double *A, const double *x, const double *b, int n);

/* Factor de crecimiento:  max|u_ij| / max|a_ij|
 * Mide cuánto "crecieron" los elementos durante la eliminación. Es la
 * cantidad que controla la estabilidad de la eliminación gaussiana. */
double growth_factor(const double *LU, const double *A0, int n);

/* --------------------------------------------------------------------------
 *  4. UTILIDADES
 * ------------------------------------------------------------------------*/

double wall_time(void);                                /* segundos, monótono  */
void   matvec(const double *A, const double *x, double *y, int n);  /* y = Ax */

/* Cantidad de operaciones de punto flotante de la factorización.
 * Se usa para reportar GFLOP/s = flops / tiempo. */
double lu_flops(int n);

#endif /* LU_H */
