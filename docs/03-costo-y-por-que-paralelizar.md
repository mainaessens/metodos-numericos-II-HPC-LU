# 03 · Costo computacional y por qué vale la pena paralelizar

> Prerrequisitos: [01](01-que-es-la-factorizacion-lu.md), [02](02-pivoteo-y-estabilidad.md).

---

## 1. Qué problema intenta resolver

Antes de paralelizar cualquier cosa hay que responder tres preguntas:

1. **¿Cuánto cuesta?** Si el programa tarda 0.2 segundos, paralelizarlo es
   perder el tiempo.
2. **¿Dónde está el costo?** No sirve optimizar el 5 % del programa.
3. **¿Ese pedazo se puede repartir?** Si tiene dependencias secuenciales, no.

Este documento contesta las tres para LU. Y las respuestas son exactamente el
contenido de la sección de "Métodos" del informe.

---

## 2. Contando operaciones

### El lazo triple

```c
for (k = 0; k < n; k++)                 //  n pasos
    for (i = k+1; i < n; i++)           //  n−k−1 filas
        for (j = k+1; j < n; j++)       //  n−k−1 columnas
            a[i][j] -= l * a[k][j];     //  1 multiplicación + 1 resta
```

En el paso `k` se actualiza una submatriz de `(n−k−1) × (n−k−1)`, y cada
elemento cuesta **2 flops** (una multiplicación y una resta). Sumando:

```
        n−1                        n−1
   W =   Σ   2·(n−k−1)²    =   2 ·  Σ  m²
        k=0                        m=1
```

Y usando la fórmula de la suma de cuadrados
`Σ_{m=1}^{M} m² = M(M+1)(2M+1)/6`:

```
   W = 2 · (n−1)n(2n−1)/6  ≈  (2/3)·n³
```

**El costo de la factorización LU es `(2/3)n³` operaciones de punto flotante.**

El término exacto que usa nuestro código (`lu_flops()` en `verify.c`) es

```
   W = (2/3)n³ − (1/2)n² − (1/6)n
```

### El costo de las sustituciones

Cada sustitución triangular recorre medio triángulo: `n²/2` elementos, 2 flops
cada uno, o sea `n²` flops por sustitución, `2n²` las dos juntas.

### La comparación que define el proyecto

| n | Factorización `(2/3)n³` | Sustituciones `2n²` | Proporción |
|---|---|---|---|
| 100 | 6.7 × 10⁵ | 2.0 × 10⁴ | 97 % / 3 % |
| 1 000 | 6.7 × 10⁸ | 2.0 × 10⁶ | 99.7 % / 0.3 % |
| 4 000 | 4.3 × 10¹⁰ | 3.2 × 10⁷ | 99.93 % / 0.07 % |
| 10 000 | 6.7 × 10¹¹ | 2.0 × 10⁸ | 99.97 % / 0.03 % |

> **Conclusión 1: todo el esfuerzo de paralelización va en la factorización.**
> Las sustituciones son ruido.

Esto es directamente la Ley de Amdahl aplicada antes de escribir una línea de
código: si el 0.07 % es serial, el speedup máximo teórico es ~1400×. No es el
cuello de botella.

---

## 3. Memoria: el otro límite

Una matriz de `n × n` en `double` (8 bytes) ocupa:

| n | Memoria de A | ¿Entra en una notebook de 8 GB? |
|---|---|---|
| 1 000 | 8 MB | de sobra |
| 2 000 | 32 MB | sí |
| 4 000 | 128 MB | sí |
| 8 000 | 512 MB | sí, ajustado si hay copias |
| 16 000 | 2 GB | no, con copias |
| 50 000 | 20 GB | **no** |

Y el tiempo crece más rápido que la memoria. Duplicar `n`:

- multiplica la memoria por **4**
- multiplica el tiempo por **8**

Con una máquina que haga 5 GFLOP/s en nuestro código escalar:

| n | Flops | Tiempo serial estimado |
|---|---|---|
| 1 000 | 6.7 × 10⁸ | 0.13 s |
| 2 000 | 5.3 × 10⁹ | 1.1 s |
| 4 000 | 4.3 × 10¹⁰ | 8.5 s |
| 8 000 | 3.4 × 10¹¹ | 68 s |
| 16 000 | 2.7 × 10¹² | 9 min |
| 32 000 | 2.2 × 10¹³ | 73 min |

> **Conclusión 2: hay un problema real que justifica HPC.** A partir de
> `n ≈ 8000` la memoria empieza a apretar en una sola máquina, y a partir de
> `n ≈ 16000` el tiempo se vuelve incómodo. Distribuir la matriz entre `p`
> procesos resuelve las dos cosas a la vez.

Guardá esta tabla: en el informe hay que justificar por qué el problema
*merece* HPC, y esto es la justificación.

---

## 4. ¿Dónde está el paralelismo?

Miremos el paso `k` con cuidado.

```
                 columna k
        ┌──────────┬───┬─────────────┐
        │  ya      │   │   ya        │
 fila k │  hecho   │ p │   hecho     │  ← fila pivote
        ├──────────┼───┼─────────────┤
        │          │a₁ │  x x x x x  │  fila k+1   ┐
        │  parte   │a₂ │  x x x x x  │  fila k+2   │ submatriz
        │  de L    │a₃ │  x x x x x  │   ...       │ activa
        │          │a₄ │  x x x x x  │  fila n−1   ┘
        └──────────┴───┴─────────────┘
```

La actualización es:

```
   a_ij  ←  a_ij − (a_ik / a_kk) · a_kj      para todo i,j > k
```

**La actualización de `a_ij` no depende de la de `a_i'j'`.** Cada elemento se
calcula con: su propio valor viejo, el multiplicador de su fila, y el elemento
de la fila pivote de su columna. Nada más.

Son `(n−k−1)²` operaciones **completamente independientes entre sí**. Eso es
paralelismo de datos puro: se pueden repartir entre tantos procesos como
quieras.

### Lo que NO se puede paralelizar

- **El lazo sobre `k`.** El paso `k+1` necesita la matriz ya actualizada por
  el paso `k`. Es estrictamente secuencial: `n` pasos en serie, sin excepción.
- **La elección del pivote dentro de un paso.** Hay que saber cuál es el
  máximo antes de poder dividir.

> **Conclusión 3: el paralelismo está DENTRO de cada paso, no entre pasos.**
> Los `n` pasos son una barrera de sincronización cada uno. Eso significa
> `n` rondas de comunicación, y es lo que va a limitar la escalabilidad.

---

## 5. Cómo se reparte: cuentas de servilleta

Si repartimos las filas entre `p` procesos, en el paso `k` cada proceso
actualiza aproximadamente `(n−k−1)/p` filas. El cómputo total por proceso:

```
   T_comp ≈ (2/3)·n³ / p
```

Y la comunicación: por cada uno de los `n` pasos hay que difundir la fila
pivote, que tiene `n−k` elementos. Sumando:

```
        n−1
   V =   Σ  (n−k)  ≈  n²/2   doubles difundidos en total
        k=0
```

Con `n` mensajes colectivos (cada uno de `O(log p)` etapas en un árbol
binomial).

### La razón que lo explica todo

```
    T_comunicación        n² · (costo por dato)         p
   ─────────────────  ≈  ─────────────────────────  ~  ───
      T_cómputo            (n³/p) · (costo flop)        n
```

**Esta es la fórmula central del proyecto.** Dice:

- Con `n` grande y `p` fijo, la comunicación se vuelve despreciable →
  **la eficiencia mejora con el tamaño del problema**.
- Con `p` grande y `n` fijo, la comunicación domina →
  **el speedup se aplana y puede hasta empeorar**.
- Lo que gobierna la escalabilidad **no es `n` ni `p` por separado, sino la
  relación `n/p`**: cuántas filas le tocan a cada proceso.

Es exactamente lo que pide analizar la consigna: *"analice cómo escala el
tiempo de cómputo y la comunicación al incrementar N"*.

---

## 6. Ley de Amdahl

Si una fracción `s` del programa es irreduciblemente serial:

```
                 1
   S_máx = ────────────────
            s + (1 − s)/p
```

Y cuando `p → ∞`:  `S_máx → 1/s`.

| s | S_máx con p=16 | S_máx con p=∞ |
|---|---|---|
| 0 % | 16.0 | ∞ |
| 1 % | 13.9 | 100 |
| 5 % | 9.1 | 20 |
| 10 % | 6.4 | 10 |
| 25 % | 3.4 | 4 |

En nuestro LU, la parte "serial" no es un bloque de código: es el **overhead
de comunicación y sincronización**, que además **crece con p**. Por eso la
curva real de speedup no solo se aplana: en algún punto empieza a **bajar**.

En la defensa: "Amdahl da una cota superior asumiendo que la parte paralela
escala perfecto. En la práctica el overhead crece con p, así que nuestra curva
se aparta de Amdahl hacia abajo."

---

## 7. Lo que hacemos y lo que hace LAPACK

Nuestro algoritmo es **escalar** y usa **BLAS nivel 2** conceptualmente:
actualizaciones vector-matriz. Por cada dato que trae de memoria hace pocas
operaciones. La razón flops / accesos a memoria es baja, y en una CPU moderna
la memoria es el cuello de botella, no la unidad aritmética.

LAPACK usa un **algoritmo por bloques**:

1. Factoriza un **panel** angosto de `b` columnas (esto es LU escalar, sobre
   una parte chiquita).
2. Actualiza toda la submatriz restante con un **producto matriz-matriz**
   (`dgemm`, BLAS nivel 3).

El `dgemm` reutiliza cada dato cargado en cache `O(b)` veces, lo que cambia la
razón flops/memoria de `O(1)` a `O(b)`. Por eso LAPACK es **10 a 50 veces más
rápido** que nuestro código para el mismo `n`, **con la misma cantidad de
operaciones**.

> **La diferencia no está en el algoritmo, está en el uso de la jerarquía de
> memoria.** Este es un punto que conviene tener listo: muestra que entendieron
> que HPC no es solo "usar más procesadores".

Si sobra tiempo, el experimento E7 compara nuestro serial contra `dgetrf` de
LAPACK. Es un cierre fuerte para el informe.

---

## 8. El experimento que sostiene todo esto

```bash
# Verificá vos mismo que el tiempo crece como n³:
for n in 500 1000 2000; do ./bin/lu_serial -n $n --csv; done
```

Al duplicar `n` el tiempo debería multiplicarse por ~8. Si no lo hace, es
porque para `n` chico el programa no está limitado por el cómputo sino por
otros efectos (cache, overhead de arranque). Ese también es un resultado
digno de comentar en el informe.

---

## 9. Errores frecuentes

| Error | Por qué |
|---|---|
| Decir que LU es `O(n²)` | Es `O(n³)`. `O(n²)` son las sustituciones |
| Reportar GFLOP/s usando el conteo de operaciones *ejecutadas* de una implementación ineficiente | El estándar es reportar contra `(2/3)n³`, el conteo del algoritmo, para que los números sean comparables entre implementaciones |
| Medir tiempos con `n = 200` | A ese tamaño se mide overhead, no el algoritmo |
| Incluir la generación de la matriz en el tiempo | Distorsiona todo, sobre todo con `n` grande |
| Concluir que el programa "no escala" con `n = 500` y `p = 16` | Con `n/p ≈ 31` filas por proceso, la comunicación domina. No es un bug, es la física del problema |

---

## 10. Preguntas típicas de parcial y defensa

1. Deducí el costo `(2/3)n³` a partir del lazo triple.
2. ¿Por qué duplicar `n` multiplica el tiempo por 8 pero la memoria solo por 4?
3. ¿Qué parte del algoritmo no se puede paralelizar y por qué?
4. Si `T_com/T_comp ≈ p/n`, ¿qué pasa con la eficiencia cuando agrandás `n`?
5. ¿Por qué LAPACK es mucho más rápido si hace las mismas operaciones?
6. ¿Cuál es el speedup máximo con 32 procesos si el 3 % del programa es serial?
7. Justificá con números por qué este problema merece HPC.

---

**Siguiente:** [04 · MPI desde cero](04-mpi-desde-cero.md)
