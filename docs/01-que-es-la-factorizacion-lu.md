# 01 · ¿Qué es la factorización LU?

> Prerrequisitos: saber multiplicar matrices. Nada más.

---

## 1. Qué problema intenta resolver

Queremos resolver un sistema de ecuaciones lineales:

```
 2x₁ +  x₂ −  x₃ =   8
−3x₁ −  x₂ + 2x₃ = −11
−2x₁ +  x₂ + 2x₃ =  −3
```

que se escribe en forma matricial como `A·x = b`:

```
A = ⎡  2   1  −1 ⎤        b = ⎡   8 ⎤        x = ⎡ x₁ ⎤
    ⎢ −3  −1   2 ⎥            ⎢ −11 ⎥            ⎢ x₂ ⎥
    ⎣ −2   1   2 ⎦            ⎣  −3 ⎦            ⎣ x₃ ⎦
```

Esto aparece en todos lados: circuitos eléctricos, estructuras, métodos de
elementos finitos, ajuste de curvas por mínimos cuadrados, el paso interno de
casi cualquier método iterativo no lineal. En problemas reales `n` no es 3,
es 10 000 o 10 000 000.

### La forma que NO se usa

`x = A⁻¹·b`. Nunca. Calcular la inversa cuesta unas tres veces más que
resolver el sistema directamente, y además es numéricamente peor: introduce
más error de redondeo. **Invertir una matriz para resolver un sistema es el
error de principiante clásico en cálculo numérico.**

---

## 2. Idea intuitiva

Los sistemas **triangulares** son fáciles de resolver. Mirá este:

```
 2x₁ +  x₂ − x₃ =  8
       3x₂ + x₃ =  5
             2x₃ = 4
```

La última ecuación te da `x₃ = 2` directo. Con ese valor, la segunda te da
`x₂ = (5 − 2)/3 = 1`. Con los dos, la primera te da `x₁`. Vas hacia atrás,
una incógnita por vez, sin resolver nada simultáneo. Eso se llama
**sustitución hacia atrás**, y cuesta apenas `n²` operaciones.

Entonces la estrategia es:

> **Transformar el sistema difícil en dos sistemas triangulares fáciles.**

Eso es exactamente lo que hace la factorización LU:

```
A = L · U
```

donde

- **L** (*lower*) es **triangular inferior**: solo tiene elementos no nulos en
  la diagonal y por debajo. Además, por convención, su diagonal es toda de unos.
- **U** (*upper*) es **triangular superior**: solo diagonal y por encima.

```
L = ⎡  1   0   0 ⎤        U = ⎡ u₁₁  u₁₂  u₁₃ ⎤
    ⎢ l₂₁  1   0 ⎥            ⎢  0   u₂₂  u₂₃ ⎥
    ⎣ l₃₁ l₃₂  1 ⎦            ⎣  0    0   u₃₃ ⎦
```

Y una vez que la tenemos:

```
A·x = b
L·U·x = b            ← sustituimos A por LU
L·(U·x) = b          ← agrupamos
```

Llamando `y = U·x`, el sistema se parte en dos:

```
(1)  L·y = b      →  sustitución hacia ADELANTE   (n² operaciones)
(2)  U·x = y      →  sustitución hacia ATRÁS      (n² operaciones)
```

Dos sistemas triangulares. Fáciles.

---

## 3. Fundamento teórico

### De dónde sale L y U

La factorización LU **no es un método nuevo**: es la eliminación gaussiana de
toda la vida, pero guardando lo que normalmente se tira.

Cuando hacés eliminación gaussiana, en cada paso restás un múltiplo de una
fila a las de abajo para hacer ceros. Ese múltiplo se llama **multiplicador**.
Normalmente lo calculás, lo usás y lo descartás.

La observación clave es:

> **La matriz que queda al final de la eliminación es U.
> Los multiplicadores que fuiste descartando son exactamente las entradas de L.**

Por eso LU no cuesta más que la eliminación gaussiana: es la misma cuenta, con
la contabilidad guardada.

### Existencia

No toda matriz invertible admite factorización LU tal cual. La condición es
que todas las submatrices principales (los "menores líderes") sean no
singulares. Ejemplo mínimo de matriz invertible SIN factorización LU:

```
A = ⎡ 0  1 ⎤     invertible (det = −1), pero a₁₁ = 0 y el primer
    ⎣ 1  0 ⎦     multiplicador sería 1/0.
```

La solución es **permutar filas**, que es de lo que trata el
[documento 02](02-pivoteo-y-estabilidad.md). Con permutaciones, el teorema es
completo:

> **Toda matriz invertible A admite una factorización `P·A = L·U`**,
> donde P es una matriz de permutación.

---

## 4. Fórmulas importantes

### El paso k de la eliminación

Para `k = 0, 1, …, n−1`:

**Multiplicadores** (para cada fila `i > k`):

```
            a_ik
   l_ik = ────────
            a_kk
```

**Actualización** (para cada fila `i > k` y cada columna `j > k`):

```
   a_ij  ←  a_ij  −  l_ik · a_kj
```

Al terminar, `U` es el triángulo superior de `A` (con diagonal) y `L` es el
triángulo inferior estricto, con unos agregados en la diagonal.

### Significado de cada variable

| Símbolo | Qué es |
|---|---|
| `n` | dimensión de la matriz (número de ecuaciones e incógnitas) |
| `k` | paso actual de la eliminación; va de 0 a n−1 |
| `a_kk` | el **pivote** del paso k: el elemento por el que dividimos |
| `l_ik` | multiplicador: cuántas veces le restamos la fila k a la fila i |
| `a_ij` | elemento de la fila i, columna j, tal como está en ese momento |
| `y` | vector intermedio, resultado de la sustitución hacia adelante |
| `x` | la solución que buscamos |

### Las dos sustituciones

**Hacia adelante** (`L·y = b`, con L de diagonal unitaria):

```
   y_i = b_i − Σ_{j<i} l_ij · y_j
```

No hay división porque los elementos de la diagonal de L son unos.

**Hacia atrás** (`U·x = y`):

```
          y_i − Σ_{j>i} u_ij · x_j
   x_i = ──────────────────────────
                    u_ii
```

Acá sí dividimos, por la diagonal de U.

---

## 5. Procedimiento paso a paso

1. Para `k = 0` hasta `n−1`:
   1. Tomar el pivote `a_kk`.
   2. Para cada fila `i` debajo de `k`: calcular `l_ik = a_ik / a_kk`.
   3. Guardar `l_ik` **en el lugar donde estaba `a_ik`** (ese lugar va a ser
      cero de todos modos, así que no se desperdicia memoria — esto se llama
      almacenamiento *in-place*).
   4. Para cada fila `i > k` y cada columna `j > k`: `a_ij ← a_ij − l_ik·a_kj`.
2. Sustitución hacia adelante para obtener `y`.
3. Sustitución hacia atrás para obtener `x`.

---

## 6. Ejemplo resuelto completo

Tomemos el sistema del principio:

```
A = ⎡  2   1  −1 ⎤        b = ⎡   8 ⎤
    ⎢ −3  −1   2 ⎥            ⎢ −11 ⎥
    ⎣ −2   1   2 ⎦            ⎣  −3 ⎦
```

### Paso k = 0

Pivote: `a₀₀ = 2`.

```
l₁₀ = −3 / 2 = −1.5
l₂₀ = −2 / 2 = −1
```

Actualizamos la fila 1:
```
fila1 ← [−3, −1, 2] − (−1.5)·[2, 1, −1]
      = [−3, −1, 2] + [3, 1.5, −1.5]
      = [0, 0.5, 0.5]
```

Y la fila 2:
```
fila2 ← [−2, 1, 2] − (−1)·[2, 1, −1]
      = [−2, 1, 2] + [2, 1, −1]
      = [0, 2, 1]
```

La matriz queda (con los multiplicadores guardados en lugar de los ceros,
que escribimos entre paréntesis):

```
⎡    2      1    −1  ⎤
⎢ (−1.5)   0.5   0.5 ⎥
⎣  (−1)     2     1  ⎦
```

### Paso k = 1

Pivote: `a₁₁ = 0.5`.

```
l₂₁ = 2 / 0.5 = 4
```

```
fila2 ← [0, 2, 1] − 4·[0, 0.5, 0.5] = [0, 0, 1 − 2] = [0, 0, −1]
```

Resultado final:

```
⎡    2      1    −1  ⎤
⎢ (−1.5)   0.5   0.5 ⎥
⎣  (−1)    (4)   −1  ⎦
```

### Leyendo L y U

```
L = ⎡   1    0   0 ⎤        U = ⎡ 2   1   −1  ⎤
    ⎢ −1.5   1   0 ⎥            ⎢ 0  0.5  0.5 ⎥
    ⎣  −1    4   1 ⎦            ⎣ 0   0   −1  ⎦
```

**Verificación** (hacela vos, es rápida y es exactamente lo que te pueden
pedir en un parcial):

```
fila 0 de L·U:  1·[2, 1, −1]                                = [2, 1, −1]      ✓
fila 1 de L·U:  −1.5·[2, 1, −1] + 1·[0, 0.5, 0.5]           = [−3, −1, 2]     ✓
fila 2 de L·U:  −1·[2, 1, −1] + 4·[0, 0.5, 0.5] + 1·[0,0,−1] = [−2, 1, 2]     ✓
```

### Resolviendo el sistema

**Sustitución hacia adelante**, `L·y = b`:

```
y₀ = 8
y₁ = −11 − (−1.5)(8) = −11 + 12 = 1
y₂ = −3 − (−1)(8) − (4)(1) = −3 + 8 − 4 = 1
```

**Sustitución hacia atrás**, `U·x = y`:

```
x₂ = 1 / (−1) = −1
x₁ = (1 − 0.5·(−1)) / 0.5 = 1.5 / 0.5 = 3
x₀ = (8 − 1·3 − (−1)(−1)) / 2 = (8 − 3 − 1) / 2 = 2
```

**Solución: `x = (2, 3, −1)`.** Verificá reemplazando en las tres ecuaciones
originales.

---

## 7. ¿Y por qué no eliminación gaussiana a secas?

Porque **la factorización se reutiliza**.

Si tenés que resolver `A·x = b` para muchos lados derechos distintos —cosa
habitual: el mismo circuito con distintas fuentes, la misma estructura con
distintas cargas, el mismo paso de un método iterativo— con eliminación
gaussiana sobre `[A|b]` tenés que rehacer todo cada vez:

| Enfoque | m lados derechos |
|---|---|
| Eliminación gaussiana sobre `[A\|b]` | `m · (2/3)n³` |
| LU una vez + m pares de sustituciones | `(2/3)n³ + m · 2n²` |

Para `n = 1000` y `m = 100`: la primera opción son `6.7 × 10¹⁰` operaciones,
la segunda `6.9 × 10⁸`. **Cien veces menos.**

Además L y U sirven para otras cosas casi gratis:

- `det(A) = (−1)^(número de intercambios) · u₀₀ · u₁₁ · ⋯ · u₍ₙ₋₁₎₍ₙ₋₁₎`
- Estimación del número de condición `κ(A)`
- Refinamiento iterativo de la solución

---

## 8. Cómo lo hace nuestro código

En `src/lu_serial.c`, la función `lu_factor()` implementa exactamente esto.
Mirá el lazo:

```c
for (int k = 0; k < n; k++) {
    ...
    double inv_akk = 1.0 / rowk[k];
    for (int i = k + 1; i < n; i++) {
        double l = rowi[k] * inv_akk;   /* multiplicador     */
        rowi[k] = l;                    /* se guarda in-place */
        for (int j = k + 1; j < n; j++)
            rowi[j] -= l * rowk[j];     /* actualización      */
    }
}
```

Tres lazos anidados: `k`, `i`, `j`. De ahí sale el costo `n³` que analizamos en
el [documento 03](03-costo-y-por-que-paralelizar.md).

Probalo con el ejemplo de esta guía:

```bash
make
./bin/lu_serial -n 3 -k random
```

---

## 9. Errores frecuentes

| Error | Por qué está mal |
|---|---|
| Calcular `A⁻¹` y multiplicar | Tres veces más caro y numéricamente peor |
| Olvidar que la diagonal de L es de unos | Si la guardás explícitamente, desperdiciás memoria y podés confundirte al leer el arreglo in-place |
| Aplicar la sustitución hacia atrás sobre `b` en vez de sobre `y` | Hay que hacer las DOS sustituciones, en orden |
| Pensar que `L·U = A` cuando hubo intercambios | Con pivoteo lo que vale es `P·A = L·U`. Hay que aplicar P a `b` también |
| Dividir dentro del lazo interno | `1.0/a_kk` se calcula una vez afuera: una división cuesta ~5× lo que una multiplicación |

---

## 10. Preguntas típicas de parcial

1. Factorizá a mano la matriz `[[4,3],[6,3]]` y resolvé `Ax = [10, 12]ᵀ`.
2. ¿Por qué la diagonal de L es de unos? ¿Es una elección o una consecuencia?
   *(Es una convención — la factorización de Doolittle. La alternativa, poner
   los unos en U, se llama factorización de Crout. Ambas existen; hay que
   fijar una para que la factorización sea única.)*
3. ¿Cuántas operaciones cuesta resolver un segundo lado derecho una vez que ya
   tenés L y U?
4. Dada una matriz invertible, ¿siempre existe `A = LU`? Dar un contraejemplo.
5. ¿Cómo calcularías el determinante a partir de la factorización?

---

**Siguiente:** [02 · Pivoteo y estabilidad numérica](02-pivoteo-y-estabilidad.md)
