# 06 · El código, explicado

> Prerrequisitos: [01](01-que-es-la-factorizacion-lu.md) a
> [05](05-lu-en-paralelo.md).
> Los archivos fuente están comentados en detalle; este documento explica la
> **arquitectura** y los puntos que no se entienden solo leyendo.

---

## 1. Por qué el código está partido así

```
include/lu.h        declaraciones comunes
src/matgen.c        generación de matrices        ─┐  no dependen de MPI:
src/verify.c        normas, residuo, reloj        ─┘  se compilan una vez
src/lu_serial.c     LU serial + main
src/lu_mpi.c        LU paralela + main
```

`matgen.c` y `verify.c` **no incluyen `mpi.h`**. Eso permite compilarlos una
sola vez con `gcc` y enlazarlos tanto al binario serial como al paralelo. Si
tuvieran código MPI adentro habría que duplicarlos o compilarlos dos veces.

Que el serial y el paralelo compartan el **mismo generador de matrices** es lo
que hace válida la comparación: los dos resuelven exactamente el mismo sistema.

---

## 2. El truco del generador sin estado

Este es el punto más sutil de todo el proyecto y suele pasar desapercibido.

### El problema

Un generador tipo `rand()` tiene **estado**: el valor que devuelve depende de
cuántas veces lo llamaste antes.

El programa serial genera la matriz fila por fila: cuando llega a la fila 5 ya
llamó a `rand()` `5n` veces. El proceso MPI de rank 2 con reparto cíclico solo
genera las filas 2, 6, 10…: cuando genera la fila 6 llamó a `rand()` apenas
`n` veces.

**Obtendrían matrices distintas.** Y entonces comparar el residuo del serial
contra el del paralelo no significaría nada.

### La solución

Un generador **sin estado**: una función hash que mezcla `(semilla, i, j)` y
devuelve un pseudoaleatorio. En `src/matgen.c`:

```c
static double rnd_ij(unsigned seed, int i, int j)
{
    uint64_t h = (uint64_t)seed * 0x9E3779B97F4A7C15ULL;
    h ^= splitmix64((uint64_t)(i + 1) * 0xD1342543DE82EF95ULL);
    h ^= splitmix64((uint64_t)(j + 1) * 0xA24BAED4963EE407ULL);
    h  = splitmix64(h);
    return (double)(h >> 11) * (1.0 / 9007199254740992.0);   /* [0,1) */
}
```

Cada elemento se puede calcular **sin haber calculado ningún otro**. Cualquier
proceso genera cualquier fila y obtiene bit a bit lo mismo que el serial.

`splitmix64` es un mezclador de 64 bits, corto y de buena calidad estadística.
El `>> 11` toma los 53 bits altos, que son exactamente la mantisa de un
`double`.

Consecuencia práctica: **`mat_generate_row()` es lo que permite que ningún
proceso tenga que alocar la matriz completa.** Sin esto habría que hacer un
`MPI_Scatter` desde el rank 0, y el rank 0 necesitaría `n²` doubles — 512 MB
con `n = 8000`, justo lo que queríamos evitar.

---

## 3. `lu_serial.c` — los tres puntos que importan

### (a) El lazo triple y el `inv_akk`

```c
const double *rowk = &A[(size_t)k * n];
double inv_akk = 1.0 / rowk[k];            /* UNA división por paso */

for (int i = k + 1; i < n; i++) {
    double *rowi = &A[(size_t)i * n];
    double l = rowi[k] * inv_akk;          /* multiplicación, no división */
    rowi[k] = l;
    if (l == 0.0) continue;
    for (int j = k + 1; j < n; j++)
        rowi[j] -= l * rowk[j];
}
```

Calcular `1.0/a_kk` una vez afuera y multiplicar adentro, en vez de dividir
`n−k−1` veces. Una división de punto flotante cuesta del orden de 4 a 5 veces
lo que una multiplicación.

**Honestidad numérica:** `a/b` y `a·(1/b)` no dan bit a bit el mismo resultado
(el `1/b` se redondea una vez más). Es una diferencia de un ulp, sin
consecuencias prácticas, y LAPACK hace exactamente lo mismo. Vale la pena
saberlo por si lo preguntan.

### (b) El `size_t` en los índices

```c
&A[(size_t)i * n]
```

Con `n = 50000`, `i * n` puede pasar de 2³¹ y **desbordar el `int`**. El
resultado es un índice negativo y un segfault, o peor, corrupción silenciosa
de memoria. El cast a `size_t` antes de multiplicar lo evita.

### (c) El cronómetro

```c
memcpy(A, A0, nn * sizeof(double));      /* restaurar: FUERA del cronómetro */
memcpy(x, b0, n  * sizeof(double));

double t0 = wall_time();
status = lu_factor(A, n, piv);           /* solo esto se mide */
double t1 = wall_time();
```

La factorización es **destructiva**: hay que restaurar `A` antes de cada
repetición. Esa copia queda fuera de la medición, igual que la generación de
la matriz y la verificación del residuo.

---

## 4. `lu_mpi.c` — el descriptor de distribución

Todo el mapeo global ↔ local vive en una struct y tres funciones `inline`:

```c
typedef struct {
    int n, p, rank;
    int block;        /* 0 = cíclico, 1 = bloques */
    int rows_per;     /* ceil(n/p) */
    int nloc;         /* filas reales de este proceso */
    int maxloc;       /* filas alocadas */
} dist_t;

static inline int owner_of (const dist_t *d, int i)  { return d->block ? i / d->rows_per : i % d->p; }
static inline int local_of (const dist_t *d, int i)  { ... }
static inline int global_of(const dist_t *d, int il) { return d->block ? d->rank*d->rows_per + il : d->rank + il*d->p; }
```

**Por qué esto importa:** el algoritmo no sabe qué reparto se está usando.
Cambiar de cíclico a bloques es un flag de línea de comandos (`-d block`), no
una rama de código. Eso es lo que hace posible el experimento E5 sin mantener
dos versiones del programa.

`inline` porque estas funciones se llaman millones de veces dentro de lazos
calientes: hay que evitar el costo de la llamada.

### `nloc` vs `maxloc`

- `maxloc = ⌈n/p⌉` es lo que **alocamos**: igual en todos los procesos, para
  simplificar.
- `nloc` es lo que cada proceso **realmente tiene**.

Con `n = 10`, `p = 4`, reparto cíclico: `maxloc = 3`; `nloc` es 3 para P0 y
P1, y 2 para P2 y P3. Nunca se recorre más allá de `nloc`.

Esta separación es exactamente lo que hace que el programa ande con `p = 3, 5,
7`, o sea, con cualquier `p` que no divida a `n`.

---

## 5. Memoria: la cuenta que hay que saber

```c
size_t locsz = (size_t)d.maxloc * (size_t)n;
double *A  = malloc(locsz * sizeof(double));   /* se factoriza    */
double *A0 = malloc(locsz * sizeof(double));   /* original, p/ residuo */
```

Por proceso: `2 · ⌈n/p⌉ · n · 8` bytes.

| n | p | Por proceso | Total (todos los procesos) |
|---|---|---|---|
| 4 000 | 1 | 256 MB | 256 MB |
| 4 000 | 8 | 32 MB | 256 MB |
| 8 000 | 8 | 128 MB | 1 GB |
| 8 000 | 32 | 32 MB | 1 GB |
| 16 000 | 32 | 128 MB | 4 GB |

**Ningún proceso tiene nunca la matriz completa.** Ni al generarla (cada uno
genera sus filas), ni al factorizar, ni al verificar (el residuo se calcula
distribuido con un `Allreduce`).

Los vectores replicados (`b0`, `pb`, `y`, `x`, `pivrow`) son `n` doubles cada
uno: con `n = 8000` eso es 64 KB por vector. Despreciable.

---

## 6. Los dos cronómetros

```c
double c_cmp = 0, c_com = 0, c_idl = 0, ts;

ts = MPI_Wtime();  /* ... búsqueda local ... */   c_cmp += MPI_Wtime() - ts;
ts = MPI_Wtime();  MPI_Allreduce(...);            c_com += MPI_Wtime() - ts;
ts = MPI_Wtime();  /* ... actualización ... */    c_cmp += MPI_Wtime() - ts;
```

Separar `T_cómputo` de `T_comunicación` es lo que permite responder la parte
de la consigna que pregunta cómo escala la comunicación. **Hay que
instrumentarlo desde el principio**: reinstrumentar después de haber hecho las
corridas significa rehacerlas todas.

### La honestidad de `--split-idle`

Hay un problema con medir así: el tiempo de un `MPI_Allreduce` incluye la
espera a que **lleguen todos los procesos**. Si un proceso tiene más trabajo
que otro, los demás esperan dentro de la colectiva, y ese tiempo se contabiliza
como "comunicación" cuando en realidad es **desbalance de carga**.

Con `--split-idle` ponemos un `MPI_Barrier` antes de cada colectiva:

```c
if (split_idle) { ts = MPI_Wtime(); MPI_Barrier(MPI_COMM_WORLD); c_idl += MPI_Wtime() - ts; }
ts = MPI_Wtime();
MPI_Allreduce(...);
c_com += MPI_Wtime() - ts;
```

La espera se va a `c_idl` y lo que queda en `c_com` es comunicación real.

**Pero las barreras cuestan.** Por eso las corridas con `--split-idle` se usan
solo para el gráfico de composición del tiempo (E4), **nunca para medir el
speedup** (E2). Cuando el programa corre sin el flag, el mensaje de salida lo
aclara explícitamente. Decir esto en el informe es exactamente el tipo de rigor
que diferencia un trabajo bueno de uno aprobado.

---

## 7. Por qué el tiempo reportado es el del más lento

```c
double dt = t1 - t0;
double dt_max;
MPI_Allreduce(&dt, &dt_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
```

Un programa paralelo termina cuando termina el **último** proceso. Reportar el
tiempo del rank 0 o el promedio da números optimistas que no corresponden a
nada observable.

En cambio, para los tiempos **por componente** (`c_cmp`, `c_com`, `c_idl`)
usamos el **promedio** entre procesos, porque lo que queremos es la
composición típica del tiempo del conjunto, no la de un proceso particular:

```c
MPI_Reduce(&t_cmp, &s_cmp, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
/* ... y después: avg_cmp = s_cmp / p */
```

---

## 8. El residuo, calculado sin juntar la matriz

```c
for (int il = 0; il < d.nloc; il++) {
    int g = global_of(&d, il);
    const double *row = &A0[(size_t)il * n];        /* fila ORIGINAL */

    double s = 0.0, a = 0.0;
    for (int j = 0; j < n; j++) { s += row[j]*xv[j];  a += fabs(row[j]); }

    double r = fabs(s - b0[g]);
    if (r > loc_rmax) loc_rmax = r;      /* max |(Ax−b)_i| local */
    if (a > loc_amax) loc_amax = a;      /* max suma de fila local */
}
MPI_Allreduce(&loc_rmax, &rmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
MPI_Allreduce(&loc_amax, &amax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
```

Las dos normas infinito (de la matriz y del residuo) son **máximos sobre
filas**, así que se calculan como máximos locales seguidos de un `Allreduce`
con `MPI_MAX`. El vector `x` está replicado, así que `‖x‖∞` es local.

Usamos `A0` (las filas **originales**, sin permutar) y `b0` (el lado derecho
original), porque el residuo se mide contra el sistema que nos dieron, no
contra el permutado.

---

## 9. Validar que el paralelo está bien

```bash
make test
```

Corre el serial y después el MPI con `p = 1, 2, 3, 4, 5, 7`. **El residuo
tiene que ser idéntico en todos los casos.**

Resultado real de una corrida con `n = 300`:

```
serial : 1.232598e-15
p = 1  : 1.232598e-15
p = 2  : 1.232598e-15
p = 3  : 1.232598e-15
p = 4  : 1.232598e-15
p = 5  : 1.232598e-15
p = 7  : 1.232598e-15
```

**Por qué exactamente igual y no "parecido":** la aritmética de punto flotante
no es asociativa, así que en general un algoritmo paralelo puede dar
resultados levemente distintos. Pero en este caso el **orden de las
operaciones aritméticas es idéntico**: cada elemento `a_ij` se actualiza con
la misma secuencia de restas sin importar quién sea su dueño. Lo único que
cambia es quién las ejecuta. Por eso el resultado es bit a bit el mismo.

Si el residuo **cambia con `p`**, hay un bug. Casi siempre es uno de estos
tres:

1. El mapeo global ↔ local está mal para algún `p`.
2. El intercambio de filas no cubre la parte de L.
3. Se difunde la fila pivote antes de haberla intercambiado.

**Probar con `p` primo y con `p` que no divida a `n` no es opcional.** `p = 4`
con `n = 1000` oculta la mitad de los bugs posibles.

---

## 10. Preguntas típicas de defensa sobre el código

1. ¿Por qué el generador de matrices no usa `rand()`?
2. ¿Cuánta memoria usa cada proceso? Deducilo.
3. ¿Por qué separan `nloc` de `maxloc`?
4. ¿Por qué el tiempo reportado es el máximo entre procesos?
5. ¿Qué hace `--split-idle` y por qué no se usa para el speedup?
6. ¿Por qué el residuo da exactamente igual con cualquier `p`?
7. ¿Por qué el cast a `size_t` en los índices?
8. ¿Cómo calculan el residuo sin juntar la matriz en un proceso?

---

**Siguiente:** [07 · Experimentos y métricas](07-experimentos-y-metricas.md)
