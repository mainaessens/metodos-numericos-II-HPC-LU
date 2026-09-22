/* ===========================================================================
 *  matgen.c — Generación reproducible de matrices y vectores de prueba
 *
 *  El punto delicado de este archivo es la REPRODUCIBILIDAD posicional:
 *  el valor de a_ij tiene que depender solo de (semilla, i, j) y de nada más.
 *
 *  ¿Por qué importa tanto?
 *
 *  Un generador clásico tipo rand() tiene ESTADO: el valor que devuelve
 *  depende de cuántas veces lo llamaste antes. Si el proceso MPI de rango 2
 *  solo genera sus propias filas, habrá llamado a rand() muchas menos veces
 *  que el programa serial cuando llega a la fila 2, y obtendría números
 *  distintos. Entonces el serial y el paralelo estarían resolviendo sistemas
 *  DIFERENTES, y compararlos no significaría nada.
 *
 *  La solución es un generador SIN ESTADO: una función hash que mezcla
 *  (semilla, i, j) y devuelve un número pseudoaleatorio. Cada proceso puede
 *  calcular cualquier elemento de la matriz sin haber calculado los
 *  anteriores. Usamos splitmix64, que es corto, rápido y de buena calidad
 *  estadística para este uso.
 * ===========================================================================*/

#include "lu.h"
#include <string.h>
#include <stdlib.h>
#include <stdint.h>

/* --------------------------------------------------------------------------
 *  splitmix64: mezcla un entero de 64 bits en otro bien distribuido.
 * ------------------------------------------------------------------------*/
static uint64_t splitmix64(uint64_t x)
{
    x += 0x9E3779B97F4A7C15ULL;
    x  = (x ^ (x >> 30)) * 0xBF58476D1CE4E5B9ULL;
    x  = (x ^ (x >> 27)) * 0x94D049BB133111EBULL;
    return x ^ (x >> 31);
}

/* Número pseudoaleatorio en [0,1) determinado por (seed, i, j). */
static double rnd_ij(unsigned seed, int i, int j)
{
    uint64_t h = (uint64_t)seed * 0x9E3779B97F4A7C15ULL;
    h ^= splitmix64((uint64_t)(i + 1) * 0xD1342543DE82EF95ULL);
    h ^= splitmix64((uint64_t)(j + 1) * 0xA24BAED4963EE407ULL);
    h  = splitmix64(h);
    /* Tomamos los 53 bits altos: es la mantisa exacta de un double. */
    return (double)(h >> 11) * (1.0 / 9007199254740992.0);   /* 2^53 */
}

/* --------------------------------------------------------------------------*/

mat_kind_t mat_kind_from_string(const char *s)
{
    if (!s)                        return MAT_DIAGDOM;
    if (!strcmp(s, "random"))      return MAT_RANDOM;
    if (!strcmp(s, "diagdom"))     return MAT_DIAGDOM;
    if (!strcmp(s, "hilbert"))     return MAT_HILBERT;
    if (!strcmp(s, "zerodiag"))    return MAT_ZERODIAG;
    return MAT_DIAGDOM;
}

const char *mat_kind_name(mat_kind_t k)
{
    switch (k) {
        case MAT_RANDOM:   return "random";
        case MAT_DIAGDOM:  return "diagdom";
        case MAT_HILBERT:  return "hilbert";
        case MAT_ZERODIAG: return "zerodiag";
    }
    return "?";
}

/* --------------------------------------------------------------------------
 *  Genera UNA fila global.
 * ------------------------------------------------------------------------*/
void mat_generate_row(double *row, int n, int i, mat_kind_t kind, unsigned seed)
{
    switch (kind) {

    case MAT_HILBERT:
        /* a_ij = 1/(i+j+1).  Simétrica, definida positiva, y con un número
         * de condición que crece como e^(3.5n): para n = 12 ya perdés todos
         * los dígitos en doble precisión. Es el caso patológico clásico. */
        for (int j = 0; j < n; j++)
            row[j] = 1.0 / (double)(i + j + 1);
        break;

    case MAT_RANDOM:
        for (int j = 0; j < n; j++)
            row[j] = 2.0 * rnd_ij(seed, i, j) - 1.0;      /* en [-1, 1) */
        break;

    case MAT_ZERODIAG:
        /* Aleatoria pero con a_ii = 0 exactamente.
         * lu_factor_nopivot() muere en el paso 0 dividiendo por cero.
         * lu_factor() (con pivoteo) la resuelve sin despeinarse. */
        for (int j = 0; j < n; j++)
            row[j] = 2.0 * rnd_ij(seed, i, j) - 1.0;
        row[i] = 0.0;
        break;

    case MAT_DIAGDOM:
    default:
        /* Diagonal estrictamente dominante: |a_ii| > sum_{j!=i} |a_ij|.
         * Se demuestra que estas matrices no necesitan pivoteo (el pivote
         * natural nunca es cero y el factor de crecimiento queda acotado).
         * Es nuestro caso de referencia "bien portado". */
        {
            double acum = 0.0;
            for (int j = 0; j < n; j++) {
                row[j] = 2.0 * rnd_ij(seed, i, j) - 1.0;
                if (j != i) acum += (row[j] < 0 ? -row[j] : row[j]);
            }
            row[i] = acum + 1.0;
        }
        break;
    }
}

void mat_generate(double *A, int n, mat_kind_t kind, unsigned seed)
{
    for (int i = 0; i < n; i++)
        mat_generate_row(&A[(size_t)i * n], n, i, kind, seed);
}

void vec_generate(double *b, int n, unsigned seed)
{
    /* Usamos la columna "virtual" j = -1 (o sea, la posición n) para que el
     * vector no coincida con ninguna columna de la matriz. */
    for (int i = 0; i < n; i++)
        b[i] = 2.0 * rnd_ij(seed ^ 0xB5297A4Du, i, n) - 1.0;
}
