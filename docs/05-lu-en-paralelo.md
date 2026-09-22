# 05 · LU en paralelo

> Prerrequisitos: [01](01-que-es-la-factorizacion-lu.md) a
> [04](04-mpi-desde-cero.md).
> **Este es el documento central del proyecto.** Lo que está acá es la sección
> de "Métodos" del informe y la mitad de la defensa.

---

## 1. Qué problema intenta resolver

Tenemos el algoritmo (doc 01), sabemos que hay que pivotear (doc 02), sabemos
dónde está el paralelismo (doc 03) y sabemos comunicar procesos (doc 04).
Falta la decisión de diseño más importante:

> **¿Quién se queda con qué parte de la matriz?**

Esa decisión —y no la elegancia del código— es la que determina si el programa
escala o no.

---

## 2. Idea intuitiva: el problema del reparto

La eliminación gaussiana **consume la matriz de arriba hacia abajo**. En el
paso `k`, las filas `0` a `k` ya están terminadas: nadie las vuelve a tocar.
Solo las filas por debajo de `k` tienen trabajo.

Eso significa que **el trabajo se va concentrando en la parte de abajo de la
matriz**. Si el reparto no tiene eso en cuenta, los procesos que tienen las
primeras filas se quedan sin nada que hacer.

### Reparto A: bloques contiguos

`dueño(fila i) = i / ⌈n/p⌉`

Con `n = 8`, `p = 4`:

```
 fila 0 → P0      ┐
 fila 1 → P0      ┘
 fila 2 → P1      ┐
 fila 3 → P1      ┘
 fila 4 → P2      ┐
 fila 5 → P2      ┘
 fila 6 → P3      ┐
 fila 7 → P3      ┘
```

**En el paso k = 4**, las filas 0–3 ya están terminadas. Entonces:

| Proceso | Filas activas |
|---|---|
| P0 | **0** ← ocioso el resto de la corrida |
| P1 | **0** ← ocioso el resto de la corrida |
| P2 | 2 |
| P3 | 2 |

La mitad de los procesos no hace **nada** durante la segunda mitad del
algoritmo. Y no es que estén un rato de más esperando: quedaron
permanentemente fuera.

### Reparto B: cíclico por filas

`dueño(fila i) = i mod p`

```
 fila 0 → P0
 fila 1 → P1
 fila 2 → P2
 fila 3 → P3
 fila 4 → P0
 fila 5 → P1
 fila 6 → P2
 fila 7 → P3
```

**En el paso k = 4:**

| Proceso | Filas activas |
|---|---|
| P0 | 1 (la fila 4) |
| P1 | 1 (la fila 5) |
| P2 | 1 (la fila 6) |
| P3 | 1 (la fila 7) |

Los cuatro siguen trabajando. Y esto vale **en cualquier paso `k`**: con
reparto cíclico, la diferencia de carga entre dos procesos cualesquiera nunca
supera **una fila**.

### El argumento formal

Con reparto en bloques, en el paso `k` el número de procesos con trabajo es
aproximadamente `p·(n−k)/n`. Promediando sobre todos los pasos, la fracción
media de procesos ocupados es **1/2**: la eficiencia máxima alcanzable es
**50 %**, incluso con comunicación gratis.

Con reparto cíclico, todos los procesos tienen trabajo mientras queden al
menos `p` filas activas, es decir hasta `k = n − p`. La eficiencia no está
limitada por el reparto.

> **Este es el resultado que el experimento E5 mide y que la figura de
> distribución del informe ilustra.** Nuestro código implementa los dos
> repartos (`-d cyclic` y `-d block`) precisamente para poder mostrarlo.

---

## 3. El algoritmo paralelo, paso a paso

Para `k = 0, 1, …, n−1`:

### Paso 1 — Búsqueda local del pivote *(sin comunicación)*

Cada proceso recorre sus filas con índice global `≥ k` y se queda con el mayor
`|a_ik|` y el índice global de esa fila.

```c
struct { double val; int idx; } loc;
loc.val = -1.0; loc.idx = -1;
for (int il = primera_local_ge(k); il < nloc; il++) {
    int g = global_de(il);
    if (g < k) continue;
    double v = fabs(A[il*n + k]);
    if (v > loc.val) { loc.val = v; loc.idx = g; }
}
```

### Paso 2 — Reducción global *(colectiva)*

```c
MPI_Allreduce(&loc, &glb, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);
```

Al volver, **todos** los procesos saben cuál es la fila pivote. Usamos
`Allreduce` y no `Reduce` porque todos necesitan el dato: el dueño de la fila
`k` para saber con quién intercambiar, y todos para saber si la operación les
toca a ellos.

Si `glb.val == 0`, la columna entera es cero: la matriz es singular.

### Paso 3 — Intercambio de filas *(punto a punto, o nada)*

Sea `prow = glb.idx`. Hay tres casos:

| Caso | Qué pasa |
|---|---|
| `prow == k` | Nada. El pivote natural ya era el máximo |
| `dueño(k) == dueño(prow)` | Swap **local**, dentro del proceso. Cero comunicación |
| Distintos dueños | Un `MPI_Sendrecv_replace` de `n` doubles entre los dos |

```c
if (rank == ok)
    MPI_Sendrecv_replace(&A[local_de(k)*n], n, MPI_DOUBLE,
                         op, TAG, op, TAG, MPI_COMM_WORLD, &st);
else if (rank == op)
    MPI_Sendrecv_replace(&A[local_de(prow)*n], n, MPI_DOUBLE,
                         ok, TAG, ok, TAG, MPI_COMM_WORLD, &st);
```

Usamos `Sendrecv_replace` y no `Send` + `Recv` por lo explicado en el
[documento 04](04-mpi-desde-cero.md): dos `MPI_Send` simétricos con `n` grande
son un deadlock garantizado.

**Nota importante:** el vector `b` está **replicado** en todos los procesos
(son `n` doubles, despreciable frente a `n²/p`). Como todos conocen `prow`
después del `Allreduce`, todos aplican el mismo swap a `b` localmente. Cero
comunicación adicional.

### Paso 4 — Difusión de la fila pivote *(colectiva)*

```c
if (rank == ok)
    memcpy(&pivrow[k], &A[local_de(k)*n + k], (n - k) * sizeof(double));

MPI_Bcast(&pivrow[k], n - k, MPI_DOUBLE, ok, MPI_COMM_WORLD);
```

**Difundimos solo los elementos `k..n−1`**, no la fila entera. Los de la
izquierda ya son parte de L y nadie los va a usar. Eso reduce el volumen total
de `n²` a `n²/2`: la mitad de la comunicación, gratis.

### Paso 5 — Actualización local *(sin comunicación)*

```c
double inv_piv = 1.0 / pivrow[k];
for (int il = primera_local_ge(k); il < nloc; il++) {
    int g = global_de(il);
    if (g <= k) continue;
    double *row = &A[il*n];
    double l = row[k] * inv_piv;
    row[k] = l;
    for (int j = k + 1; j < n; j++)
        row[j] -= l * pivrow[j];
}
```

**Acá está el 99 % del tiempo de ejecución.** Es puro cómputo local, sin una
sola llamada a MPI. Es el trabajo que efectivamente se repartió.

---

## 4. El mapeo global ↔ local

Esta es la parte que más bugs genera, así que conviene tenerla en un solo lugar
y bien clara. En `lu_mpi.c` está encapsulada en tres funciones inline:

### Reparto cíclico

```
   dueño(i)        = i mod p
   índice_local(i) = i / p             (división entera)
   global(il)      = rank + il · p
```

Verificación con `p = 4`, `rank = 2`:
`global(0) = 2`, `global(1) = 6`, `global(2) = 10`. ✓
Y `dueño(6) = 6 mod 4 = 2` ✓, `índice_local(6) = 6/4 = 1` ✓.

**Cuántas filas tiene cada proceso:**

```
   nloc = n/p + (rank < n mod p ? 1 : 0)
```

Con `n = 10`, `p = 4`: P0 y P1 tienen 3 filas, P2 y P3 tienen 2. **Ningún
proceso queda vacío y la diferencia máxima es 1.**

### Reparto en bloques

```
   filas_por_proc  = ⌈n/p⌉
   dueño(i)        = i / filas_por_proc
   índice_local(i) = i − dueño(i)·filas_por_proc
   global(il)      = rank · filas_por_proc + il
```

Con `n = 10`, `p = 4`: `filas_por_proc = 3`, así que P0, P1 y P2 tienen 3
filas y P3 tiene **1**. Nótese que el desbalance ya existe antes de empezar.

### La optimización `primera_local_ge(k)`

En cada paso `k` no queremos recorrer desde `il = 0` todas las filas ya
terminadas: eso convertiría la búsqueda del pivote en `O(n²/p)` cuando debería
ser `O(n²/2p)`. La función devuelve el primer índice local cuya fila global es
`≥ k`:

```c
/* cíclico */
int diff = k - rank;
return (diff <= 0) ? 0 : (diff + p - 1) / p;      /* techo(diff/p) */
```

---

## 5. Las sustituciones triangulares en paralelo

Después de factorizar hay que resolver. `L·y = Pb` y después `U·x = y`.

La clave es que **cada fila la tiene un solo proceso, y ese proceso tiene la
fila completa** (los `n` elementos). Entonces el dueño de la fila `i` puede
calcular `y_i` por sí solo, siempre que ya conozca `y_0 … y_{i−1}`:

```c
for (int i = 0; i < n; i++) {
    int oi = dueño(i);
    if (rank == oi) {
        const double *row = &A[local_de(i)*n];
        double s = pb[i];
        for (int j = 0; j < i; j++) s -= row[j] * y[j];
        y[i] = s;
    }
    MPI_Bcast(&y[i], 1, MPI_DOUBLE, oi, MPI_COMM_WORLD);
}
```

Un `Bcast` de **un solo double** por fila. `n` colectivas chiquititas.

**Es el pedazo menos escalable del programa**, y conviene decirlo en el
informe antes de que lo pregunten: son `n` sincronizaciones para `O(n²)` de
cómputo, mientras que la factorización tiene `n` sincronizaciones para
`O(n³)` de cómputo. La proporción es `n` veces peor.

Lo dejamos así por dos razones: (a) `O(n²)` frente a `O(n³)` sigue siendo
despreciable en el tiempo total, y (b) el código queda legible. La forma seria
de mejorarlo es procesar bloques de `b` filas por vez, reduciendo el número de
colectivas de `n` a `n/b`.

---

## 6. Modelo de costos

### Cómputo

```
   T_comp ≈ (2/3)·n³ / (p · γ)
```

donde `γ` es la velocidad de la máquina en flops/segundo.

### Comunicación

Modelamos el costo de un mensaje como `α + β·m`, donde `α` es la **latencia**
(costo fijo por mensaje, del orden de microsegundos) y `β` el inverso del
**ancho de banda** (costo por byte).

Por cada paso `k`:

| Operación | Mensajes | Volumen |
|---|---|---|
| `Allreduce` del pivote | `O(log p)` | 12 bytes |
| `Sendrecv` del intercambio | 1 (a veces 0) | `n` doubles |
| `Bcast` de la fila pivote | `O(log p)` | `n−k` doubles |

Sumando sobre los `n` pasos:

```
   T_com ≈ n·log(p)·α  +  (n²/2)·log(p)·β
           └─ latencia ┘   └─ ancho de banda ┘
```

### La razón que decide todo

```
    T_com        n²·log(p)·β        p·log(p)·β·γ
   ────────  ≈  ───────────────  =  ──────────────
    T_comp       (n³/p)/γ                n
```

O sea, salvo constantes:

```
    T_com       p
   ──────  ~  ─────
    T_comp      n
```

### Qué predice el modelo

| Situación | Predicción |
|---|---|
| `n` crece con `p` fijo | La comunicación se diluye → **la eficiencia mejora** |
| `p` crece con `n` fijo | La comunicación domina → **el speedup se aplana y puede caer** |
| `n/p` chico (pocas filas por proceso) | Régimen dominado por latencia: no escala |
| `n/p` grande | Régimen dominado por cómputo: escala bien |

**El parámetro que gobierna todo es `n/p`: cuántas filas tiene cada proceso.**

Regla práctica: por debajo de unas 50–100 filas por proceso, agregar procesos
deja de ayudar. Con `n = 2000` y `p = 32` son 62 filas por proceso: ya estás
en la zona gris. Con `n = 8000` y `p = 32` son 250: ahí todavía hay margen.

Esto es literalmente la respuesta a la consigna: *"analice cómo escala el
tiempo de cómputo y la comunicación al incrementar N"*.

---

## 7. Lo que hacen las bibliotecas serias

Nuestra implementación es 1D por filas. Las bibliotecas de producción hacen
tres cosas más:

### (a) Distribución 2D block-cyclic

ScaLAPACK organiza los procesos en una **grilla** de `√p × √p` y reparte
bloques de la matriz cíclicamente en las dos dimensiones.

**Por qué es mejor:** con reparto 1D por filas, la difusión del panel involucra
a los `p` procesos. Con una grilla `√p × √p`, solo involucra a los `√p` de la
misma columna de la grilla. El volumen de comunicación baja de `O(n²·p)` a
`O(n²·√p)`.

Con `p = 4` la diferencia es menor; con `p = 1024` es la diferencia entre
escalar y no escalar.

### (b) Algoritmo por bloques

En vez de eliminar columna por columna, se factoriza un **panel** de `b`
columnas y se actualiza la submatriz restante con un producto matriz-matriz
(`dgemm`, BLAS nivel 3). Esto reutiliza los datos en cache y sube el
rendimiento un orden de magnitud. Es la razón por la que LAPACK le gana a
nuestro código serial.

### (c) Pivoteo que evita comunicación (CALU)

Nuestro algoritmo hace **una reducción global por columna**: `n` reducciones.
El pivoteo por torneo (*tournament pivoting*) de CALU elige buenos pivotes para
un panel entero de `b` columnas con una sola reducción en forma de torneo,
bajando el número de reducciones de `n` a `n/b`.

Se demuestra que CALU alcanza las cotas inferiores teóricas de comunicación
para la factorización LU (Grigori, Demmel & Xiang, 2011). Es el estado del
arte, y SLATE —la biblioteca del proyecto exascale del Departamento de Energía
de EE.UU.— lo implementa con la factorización de paneles corriendo en GPU
(Gates et al., 2025).

> **Para el informe:** mencionar estas tres cosas en la introducción (estado
> del arte) y en las conclusiones (trabajo futuro) demuestra que entendieron
> dónde está parado su trabajo. No hace falta implementarlas.

---

## 8. Errores frecuentes

| Error | Síntoma |
|---|---|
| Usar reparto en bloques "porque es más simple" | Eficiencia tope 50 %, sin explicación aparente |
| Alocar la matriz completa en cada proceso | Anda con n=1000, se queda sin memoria con n=8000 |
| Hacer `Scatter` de la matriz desde el rank 0 | El rank 0 necesita `n²` doubles: mismo problema |
| Difundir la fila pivote entera en vez de desde `k` | El doble de comunicación, gratis |
| Intercambiar solo la parte derecha de las filas | `PA ≠ LU`, resultados mal sin error |
| `MPI_Send` + `MPI_Recv` para el intercambio | Deadlock con `n` grande |
| No usar `primera_local_ge(k)` | La búsqueda del pivote pasa a ser `O(n²/p)` |
| Suponer `n` múltiplo de `p` | Se rompe con p = 3, 5, 7. Hay que probarlos |

---

## 9. Preguntas típicas de defensa

1. ¿Por qué eligieron reparto cíclico? Demostralo con un dibujo.
2. ¿Cuál es la eficiencia máxima teórica con reparto en bloques? ¿Por qué?
3. ¿Por qué `MPI_Allreduce` y no `MPI_Reduce` en el pivoteo?
4. ¿Por qué `MPI_MAXLOC` y no `MPI_MAX`?
5. ¿Por qué difunden solo desde la columna `k`?
6. ¿Qué pasa si `n` no es múltiplo de `p`?
7. Escribí el modelo de costos y explicá qué predice al crecer `N`.
8. ¿Qué cambiarían para correr en 1024 procesos?
9. ¿Por qué las sustituciones triangulares son el pedazo menos escalable?
10. ¿Dónde está la comunicación y dónde el cómputo en el paso `k`?

---

**Siguiente:** [06 · El código, explicado](06-el-codigo-explicado.md)
