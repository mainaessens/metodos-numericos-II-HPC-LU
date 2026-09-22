# 02 · Pivoteo y estabilidad numérica

> Prerrequisitos: [documento 01](01-que-es-la-factorizacion-lu.md).
> Este es el documento que más preguntas de defensa genera.

---

## 1. Qué problema intenta resolver

El algoritmo del documento anterior tiene dos agujeros:

**Agujero 1 — el pivote puede ser cero.**

```
A = ⎡ 0  1 ⎤
    ⎣ 1  0 ⎦
```

Es invertible (determinante −1), pero `l₁₀ = 1/0`. El algoritmo explota en el
primer paso.

**Agujero 2 — el pivote puede ser muy chico, que es peor.**

Peor porque no explota: devuelve un número, y ese número está mal. Un error
silencioso es más peligroso que un crash.

---

## 2. Idea intuitiva

El problema es que estamos **dividiendo** por el pivote. Dividir por un número
chico produce un número enorme, y sumar un número enorme a uno normal en
aritmética de punto flotante **borra** el número normal: los bits del chico no
entran en la mantisa.

La solución es de una simpleza notable:

> **Antes de dividir, buscá el elemento más grande de la columna y traelo
> arriba intercambiando filas.**

Intercambiar dos ecuaciones de un sistema no cambia la solución. Es
matemáticamente gratis. Y garantiza que todos los multiplicadores cumplan
`|l_ik| ≤ 1`, lo que impide que los números crezcan descontroladamente.

Eso se llama **pivoteo parcial**.

---

## 3. El ejemplo que hay que saber contar

Este ejemplo tiene que estar en el informe y en las slides. Es el argumento
completo en cuatro líneas de cuentas.

Sistema:

```
0.0001·x + 1·y = 1
     1·x + 1·y = 2
```

La solución exacta es `x = 1.00010001…`, `y = 0.99989999…`, o sea
aproximadamente `x = 1`, `y = 1`.

Vamos a hacer las cuentas con **3 dígitos significativos**, para simular en
chiquito lo que le pasa a una computadora con 16 dígitos.

### SIN pivoteo

El pivote es `0.0001`. El multiplicador:

```
l = 1 / 0.0001 = 10000
```

Actualizamos la segunda fila:

```
y:  1 − 10000·1 = −9999   →  a 3 dígitos:  −10000
b:  2 − 10000·1 = −9998   →  a 3 dígitos:  −10000
```

Fijate lo que pasó: el `1` original y el `2` original **desaparecieron**. Se
los comió el redondeo frente a 10000. Toda la información de la segunda
ecuación se perdió.

Sustitución hacia atrás:

```
y = −10000 / −10000 = 1.00        ← este está bien, por casualidad
x = (1 − 1·1.00) / 0.0001 = 0 / 0.0001 = 0
```

**`x = 0`.** La respuesta correcta es `x ≈ 1`. Cien por ciento de error.

### CON pivoteo

Intercambiamos las filas porque `|1| > |0.0001|`:

```
     1·x + 1·y = 2
0.0001·x + 1·y = 1
```

Multiplicador:

```
l = 0.0001 / 1 = 0.0001
```

Actualización:

```
y:  1 − 0.0001·1 = 0.9999   →  a 3 dígitos:  1.00
b:  1 − 0.0001·2 = 0.9998   →  a 3 dígitos:  1.00
```

Sustitución hacia atrás:

```
y = 1.00 / 1.00 = 1.00
x = (2 − 1·1.00) / 1 = 1.00
```

**`x = 1`, `y = 1`.** Correcto.

### La moraleja

Misma matriz, mismo algoritmo, misma aritmética. Lo único que cambió fue el
**orden de las filas**. Sin pivoteo el resultado es basura; con pivoteo es
exacto hasta la precisión disponible.

**Y la matriz de este ejemplo está perfectamente bien condicionada.** No es un
problema difícil: es un algoritmo mal implementado. Esa distinción —problema
mal condicionado vs. algoritmo inestable— es exactamente lo que van a querer
escucharles decir en la defensa.

---

## 4. Fundamento teórico

### El teorema

> Para toda matriz invertible `A` existe una matriz de permutación `P` tal que
> ```
> P · A = L · U
> ```
> donde L es triangular inferior con diagonal unitaria y todos sus elementos
> cumplen `|l_ij| ≤ 1`, y U es triangular superior.

La cota `|l_ij| ≤ 1` es consecuencia directa de elegir siempre el máximo: si
`|a_kk|` es el mayor de la columna, entonces `|a_ik / a_kk| ≤ 1` para todo `i`.

### El factor de crecimiento

La cantidad que gobierna la estabilidad es

```
        max_ij |u_ij|
   ρ = ───────────────
        max_ij |a_ij|
```

Mide cuánto crecieron los elementos durante la eliminación. La cota clásica
del error backward de la eliminación gaussiana es proporcional a `ρ`:

```
   ‖ΔA‖ ≲ c · n · ρ · ε · ‖A‖
```

donde `ε ≈ 2.2·10⁻¹⁶` es el epsilon de la máquina en doble precisión.

- **Sin pivoteo**, `ρ` no está acotado: puede ser arbitrariamente grande. El
  algoritmo es **inestable**.
- **Con pivoteo parcial**, la cota teórica es `ρ ≤ 2^(n−1)`. Para n = 100 eso
  es `10³⁰`, que es una cota inútil. **Pero en la práctica `ρ` casi nunca pasa
  de 10 o 20.**

Ese hueco entre la teoría y la práctica es un problema abierto clásico del
análisis numérico (Higham, 2002, cap. 9). Se conocen matrices patológicas que
alcanzan la cota `2^(n−1)` —la matriz de Wilkinson es el ejemplo estándar—
pero son de medida cero: no aparecen en problemas reales.

**Cómo decirlo en la defensa:** "El pivoteo parcial no tiene una garantía
teórica buena, pero tiene 70 años de evidencia empírica y es lo que usan
LAPACK y todas las bibliotecas serias. Nuestro código mide el factor de
crecimiento en cada corrida para verificarlo."

---

## 5. Pivoteo parcial vs. pivoteo completo

| | Parcial | Completo |
|---|---|---|
| Qué busca | el máximo de la **columna k**, de la fila k para abajo | el máximo de **toda la submatriz** activa |
| Costo de la búsqueda | `O(n²)` comparaciones en total | `O(n³)` comparaciones en total |
| Qué permuta | solo filas: `PA = LU` | filas y columnas: `PAQ = LU` |
| Cota de crecimiento | `2^(n−1)` (pesimista) | `O(n^(½ log n))` (mucho mejor) |
| ¿Se usa? | **sí, en todos lados** | casi nunca |

El pivoteo completo es más estable pero su búsqueda cuesta lo mismo que el
algoritmo entero: duplicarías el tiempo de ejecución para protegerte de un
caso que en la práctica no ocurre. Además permutar columnas reordena las
incógnitas, lo que agrega contabilidad.

**Conclusión para el informe:** el pivoteo parcial es el compromiso adoptado
por LAPACK y por nosotros.

---

## 6. Cómo se representa P

Nunca como una matriz. Una matriz de permutación de `n × n` tendría `n²`
elementos de los cuales solo `n` son distintos de cero: sería un desperdicio
absurdo de memoria y de tiempo.

Se guarda como un **vector de enteros** `piv[n]`, donde `piv[k]` es el índice
de la fila con la que se intercambió la fila `k` en el paso `k`.

Para aplicar `P` a un vector `b`, se repiten los mismos intercambios en el
mismo orden:

```c
for (int k = 0; k < n; k++) {
    int p = piv[k];
    if (p != k) { double t = b[k]; b[k] = b[p]; b[p] = t; }
}
```

**Detalle que se pasa por alto:** al intercambiar dos filas de la matriz hay
que intercambiar la **fila completa**, no solo la parte que todavía es U. La
parte izquierda contiene multiplicadores de L que ya calculamos, y esos
también tienen que viajar: `L` corresponde a la matriz *permutada*, no a la
original. En `lu_serial.c`:

```c
if (p != k) {
    double *rk = &A[(size_t)k * n];
    double *rp = &A[(size_t)p * n];
    for (int j = 0; j < n; j++) {      /* j desde 0, no desde k */
        double t = rk[j]; rk[j] = rp[j]; rp[j] = t;
    }
}
```

---

## 7. Residuo, error y condicionamiento

Esta es la sección que hay que tener clarísima, porque es la que confunde a
todo el mundo.

### El residuo

```
   r = A·x̂ − b
```

donde `x̂` es la solución que calculó la computadora. Lo normalizamos:

```
        ‖A·x̂ − b‖∞
   r = ──────────────
        ‖A‖∞ · ‖x̂‖∞
```

En nuestro código lo calcula `residual_rel()` en `src/verify.c`.

### Qué significa un residuo chico

Un residuo del orden de `ε ≈ 10⁻¹⁶` significa que **`x̂` es la solución exacta
de un sistema levemente perturbado**:

```
   (A + ΔA) · x̂ = b     con    ‖ΔA‖ ≈ r · ‖A‖
```

Eso se llama **estabilidad backward**: el algoritmo no agregó más error del
que ya tendrías si hubieras medido los datos de entrada con precisión de
máquina. Es lo máximo que se le puede pedir a un algoritmo.

### Lo que un residuo chico NO significa

**No significa que `x̂` esté cerca de la solución verdadera.** La relación es:

```
   ‖x̂ − x‖                      ‖A·x̂ − b‖
   ─────────   ≲   κ(A) ·  ──────────────────
    ‖x‖                       ‖A‖ · ‖x̂‖
```

donde `κ(A) = ‖A‖·‖A⁻¹‖` es el **número de condición**.

- Error **backward** = el residuo. Es culpa del **algoritmo**.
- Error **forward** = qué tan lejos está `x̂` de `x`. Es culpa del
  **algoritmo × el problema**.
- `κ(A)` = qué tan sensible es el **problema**. No es culpa de nadie: es una
  propiedad de la matriz.

**Regla práctica:** con `κ(A) ≈ 10^s`, perdés aproximadamente `s` dígitos
significativos. En doble precisión arrancás con 16.

### La matriz de Hilbert

```
   a_ij = 1 / (i + j + 1)
```

Su número de condición crece como `e^(3.5n)`:

| n | κ(A) aproximado | dígitos correctos esperables |
|---|---|---|
| 5 | 4.8 × 10⁵ | ~10 |
| 10 | 1.6 × 10¹³ | ~3 |
| 13 | 3 × 10¹⁸ | 0 |

Para `n = 13` en doble precisión ya no queda **ningún** dígito confiable en la
solución… **aunque el residuo siga siendo del orden de 10⁻¹⁶.**

Ese es el experimento E6 del proyecto, y es lo que hay que mostrar en el
informe: la gráfica del residuo se mantiene plana en `ε` para todos los tipos
de matriz, pero eso no quiere decir que todas las soluciones sean igual de
buenas.

Probalo:

```bash
./bin/lu_serial -n 12 -k hilbert
./bin/lu_serial -n 12 -k diagdom
```

---

## 8. Las matrices de prueba de nuestro código

`src/matgen.c` genera cuatro tipos, cada uno para responder una pregunta
distinta:

| `-k` | Qué es | Para qué la usamos |
|---|---|---|
| `diagdom` | diagonal estrictamente dominante: `\|a_ii\| > Σ_{j≠i} \|a_ij\|` | El caso de referencia. Se demuestra que **no necesita pivoteo**. Todos los experimentos de escalabilidad usan esta. |
| `random` | uniforme en [−1,1) | Caso genérico, condicionamiento medio |
| `hilbert` | `1/(i+j+1)` | Mal condicionada: residuo vs. error |
| `zerodiag` | aleatoria con toda la diagonal en cero | **Sin pivoteo falla en el paso 0.** Es la demostración |

La demostración de una línea:

```bash
./bin/lu_serial -n 200 -k zerodiag             # anda: residuo ~1e-16
./bin/lu_serial -n 200 -k zerodiag --nopivot   # FALLO: pivote nulo en k=0
```

---

## 9. Errores frecuentes

| Error | Consecuencia |
|---|---|
| Buscar el máximo sin valor absoluto | Elegís un número grande negativo… o uno chico positivo. Hay que comparar `fabs()` |
| Intercambiar solo desde la columna k | La parte de L se desincroniza y `PA ≠ LU`. Los resultados dan mal sin mensaje de error |
| Olvidar aplicar P a `b` | Resolvés un sistema con las ecuaciones en otro orden |
| Confundir residuo con error | "El residuo dio 1e-16 así que la solución es exacta" es falso |
| Decir que Hilbert falla por culpa del algoritmo | Falla por culpa del **problema**: κ(A) enorme. El algoritmo se portó bien (residuo chico) |

---

## 10. Preguntas típicas de parcial y defensa

1. ¿Por qué el pivoteo parcial y no el completo?
2. Mostrá con un ejemplo de 2×2 qué pasa sin pivoteo.
3. Si el residuo es `10⁻¹⁵`, ¿la solución es correcta? Justificá.
4. ¿Qué es el factor de crecimiento y por qué importa?
5. ¿Cómo se representa la matriz de permutación en memoria y por qué?
6. Una matriz diagonal dominante, ¿necesita pivoteo? ¿Por qué?
7. Diferencia entre un problema mal condicionado y un algoritmo inestable.
8. ¿Cuántos dígitos correctos esperás si `κ(A) = 10⁸`?

---

**Siguiente:** [03 · Costo y por qué paralelizar](03-costo-y-por-que-paralelizar.md)
