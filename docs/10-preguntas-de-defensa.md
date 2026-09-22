# 10 · Banco de preguntas de defensa

> Todas las preguntas, con respuestas trabajadas. Úsenlo en modo examen: una
> persona pregunta, otra contesta **en voz alta**, la tercera controla contra
> el texto. Si la respuesta no sale fluida, ese es el tema a repasar.
>
> Regla del ensayo: **cada persona contesta las preguntas de la parte que no
> programó.**

---

## Bloque A — Fundamentos de LU

**A1. ¿Qué es la factorización LU y para qué sirve?**

Descompone una matriz cuadrada invertible `A` en el producto de una triangular
inferior `L` (con diagonal de unos) por una triangular superior `U`. Sirve para
resolver `Ax = b`: el sistema se parte en dos triangulares, `Ly = b` hacia
adelante y `Ux = y` hacia atrás, que cuestan `n²` cada uno. No es un método
nuevo: es la eliminación gaussiana guardando los multiplicadores en vez de
descartarlos.

**A2. ¿Por qué no calcular `A⁻¹` y multiplicar?**

Cuesta unas tres veces más y es numéricamente peor: introduce más error de
redondeo. Invertir una matriz para resolver un sistema es el error clásico.

**A3. ¿Por qué factorizar en vez de hacer eliminación sobre `[A|b]`?**

Porque la factorización se reutiliza. Con `m` lados derechos, la eliminación
sobre la aumentada cuesta `m·(2/3)n³` y LU cuesta `(2/3)n³ + m·2n²`. Con
`n = 1000` y `m = 100`, es cien veces menos.

**A4. ¿Por qué la diagonal de `L` es de unos?**

Es una convención (factorización de Doolittle) necesaria para que la
factorización sea **única**: sin fijar una normalización, hay infinitas
factorizaciones. La alternativa —unos en `U`— se llama factorización de Crout.
Además permite guardar `L` y `U` in-place en un solo arreglo.

**A5. ¿Toda matriz invertible admite `A = LU`?**

No. `[[0,1],[1,0]]` es invertible (determinante −1) pero su primer pivote es
cero. Con permutaciones sí: **toda matriz invertible admite `PA = LU`**.

**A6. ¿Cuántas operaciones cuesta?**

`(2/3)n³` para la factorización, `2n²` para las dos sustituciones. Se deduce
sumando `2(n−k−1)²` sobre `k`, que es `2·Σm²  ≈ (2/3)n³`.

---

## Bloque B — Pivoteo y estabilidad

**B1. ¿Por qué el pivoteo parcial y no el completo?**

El completo busca el máximo en toda la submatriz: `O(n³)` comparaciones, del
mismo orden que el algoritmo entero, y además permuta columnas, lo que
reordena las incógnitas. El parcial cuesta `O(n²)` y en la práctica mantiene
el factor de crecimiento acotado. Es el compromiso que adopta LAPACK.

**B2. Mostrá con un ejemplo qué pasa sin pivoteo.**

*(Escribir en el pizarrón el sistema `0.0001x + y = 1`, `x + y = 2` con 3
dígitos significativos. Ver [doc 02 §3](02-pivoteo-y-estabilidad.md). Sin
pivoteo da `x = 0`; con pivoteo da `x = 1`, que es lo correcto. Misma matriz,
mismo algoritmo: lo único que cambió es el orden de las filas.)*

**B3. ¿Qué es el factor de crecimiento?**

`ρ = max|u_ij| / max|a_ij|`. Mide cuánto crecieron los elementos durante la
eliminación, y es la cantidad que acota el error backward. Sin pivoteo no está
acotado. Con pivoteo parcial la cota teórica es `2^(n−1)` —inútil— pero en la
práctica ρ casi nunca pasa de 10 o 20. Nuestro código lo mide en cada corrida.

**B4. Si el residuo es `10⁻¹⁵`, ¿la solución es correcta?**

No necesariamente. Un residuo chico significa **estabilidad backward**:
resolvimos exactamente un sistema cercano al nuestro. El error en la solución
está acotado por `κ(A)` por el residuo relativo. Con una matriz de Hilbert de
`n = 12`, `κ(A) ≈ 10¹³` y la solución puede tener 3 dígitos correctos aunque
el residuo sea mínimo. Lo mostramos en el experimento E6.

**B5. Diferencia entre problema mal condicionado y algoritmo inestable.**

El condicionamiento es una propiedad del **problema**: cuánto amplifica los
errores de entrada. La estabilidad es una propiedad del **algoritmo**: cuánto
error agrega él. LU sin pivoteo es inestable incluso con problemas bien
condicionados (ejemplo B2). Hilbert está mal condicionada aunque el algoritmo
se porte perfecto.

**B6. ¿Cómo representan `P` en memoria?**

Como un vector de enteros `piv[n]`, donde `piv[k]` es la fila con la que se
intercambió la fila `k`. Guardar una matriz de `n²` para representar `n`
intercambios sería un desperdicio.

**B7. Una matriz diagonal dominante, ¿necesita pivoteo?**

No. Se demuestra que en ese caso el pivote natural nunca es cero y el factor
de crecimiento queda acotado. Por eso todos nuestros experimentos de
escalabilidad usan matrices diagonal dominantes: aislan el efecto de la
paralelización del efecto del pivoteo.

---

## Bloque C — Paralelización

**C1. ¿Dónde está el paralelismo en LU?**

En la actualización de la submatriz dentro de cada paso `k`: son `(n−k−1)²`
operaciones completamente independientes entre sí. El lazo sobre `k` **no** se
puede paralelizar: el paso `k+1` necesita la matriz actualizada por el `k`.
Son `n` barreras de sincronización obligatorias.

**C2. ¿Por qué distribución cíclica y no en bloques?**

Porque la eliminación consume la matriz de arriba hacia abajo. Con bloques
contiguos, los procesos dueños de las primeras filas se quedan sin trabajo
activo a mitad de camino: la fracción media de procesos ocupados es 1/2, o sea
que la eficiencia máxima alcanzable es 50 % aunque la comunicación fuera
gratis. Con reparto cíclico, la diferencia de carga entre procesos nunca supera
una fila. El experimento E5 lo confirma.

**C3. ¿Qué se comunica en cada paso?**

Dos colectivas y a veces un punto a punto: (1) `MPI_Allreduce` con `MPI_MAXLOC`
para encontrar el pivote, porque nadie tiene la columna entera; (2) un
`MPI_Sendrecv_replace` para intercambiar las filas, solo si viven en procesos
distintos; (3) `MPI_Bcast` de la fila pivote desde su dueño. La actualización
es 100 % local.

**C4. ¿Por qué `MPI_MAXLOC` y no `MPI_MAX`?**

Porque no alcanza con el valor máximo: necesitamos saber **qué fila** lo
contiene para intercambiarla. `MAXLOC` sobre `MPI_DOUBLE_INT` devuelve el par
(valor, índice) en una sola colectiva; con `MAX` haría falta una segunda
comunicación.

**C5. ¿Por qué `Allreduce` y no `Reduce`?**

Porque **todos** los procesos necesitan el resultado: el dueño de la fila `k`
para saber con quién intercambiar, y todos para saber si la operación les toca.
`Allreduce` es más eficiente que `Reduce` + `Bcast`.

**C6. ¿Por qué difunden solo desde la columna `k`?**

Porque los elementos a la izquierda de `k` ya son parte de `L` y nadie los va a
usar en la actualización. Difundir `n−k` en vez de `n` reduce el volumen total
de comunicación de `n²` a `n²/2`.

**C7. ¿Por qué `MPI_Sendrecv` y no `Send` + `Recv`?**

Porque dos `MPI_Send` simétricos con mensajes grandes son un deadlock. `Send`
solo retorna inmediatamente si el mensaje entra en el buffer interno del
sistema; con `n` grande espera al receptor, y si los dos esperan, se cuelgan.
Es el mismo problema de `anillo.c` del TP1 con `N = 1000000`.

**C8. ¿Qué pasa si `p` no divide a `n`?**

Con reparto cíclico no pasa nada especial: el proceso `r` recibe las filas
`r, r+p, r+2p…` y algunos quedan con una fila más. Lo que hay que cuidar es
alocar `⌈n/p⌉` filas por proceso y separar `nloc` de `maxloc`. Lo verificamos
corriendo con `p = 3, 5, 7`, que son los que rompen los códigos mal escritos.

**C9. Escribí el modelo de costos.**

`T_comp ≈ (2/3)n³/(p·γ)`, `T_com ≈ n·log(p)·α + (n²/2)·log(p)·β`. El cociente
es `T_com/T_comp ~ p/n`. Predice que la eficiencia **mejora** al crecer `N` a
`p` fijo, y **empeora** al crecer `p` a `N` fijo. El parámetro que gobierna
todo es `n/p`: cuántas filas tiene cada proceso.

---

## Bloque D — Resultados y métricas

**D1. ¿Por qué `T₁` es el serial puro?**

Porque si usáramos la versión MPI con `p = 1`, el mismo overhead estaría
arriba y abajo de la fracción y se cancelaría, inflando el speedup. `T₁` es el
tiempo del **algoritmo secuencial**, compilado con las mismas banderas.

**D2. ¿Por qué el speedup no es lineal?**

Por tres causas simultáneas: la parte serial que no se paraleliza (Amdahl);
el costo de comunicación, que son `O(n)` colectivas de `O(log p)` etapas y
crece con `p`; y el tiempo ocioso por desbalance, porque en el paso `k` solo
trabajan los procesos con filas `> k`. Nuestra descomposición del tiempo (E4)
cuantifica cuánto pesa cada una.

**D3. ¿Por qué la eficiencia mejora al agrandar `N`?**

Porque el cómputo crece como `n³` y la comunicación como `n²`. A `p` fijo,
duplicar `n` reduce a la mitad el peso relativo de la comunicación. Es la razón
por la que HPC tiene sentido para problemas grandes y no para chicos.

**D4. ¿Por qué la mediana y no el promedio?**

Porque el cluster no es dedicado: la distribución de tiempos tiene una cola
larga hacia arriba (corridas contaminadas por otros jobs) pero no hacia abajo.
El promedio se deja arrastrar por esa cola; la mediana es robusta. Reportamos
también mínimo y máximo para mostrar cuánto ruido había.

**D5. ¿Puede haber speedup superlineal?**

Sí, y hay que explicarlo, no esconderlo. Suele ser efecto de cache: al
repartir la matriz, la porción de cada proceso entra en niveles de cache más
rápidos que la matriz completa en la corrida serial. No viola la teoría: es
que `T₁` está medido con una jerarquía de memoria peor.

**D6. Diferencia entre escalabilidad fuerte y débil.**

Fuerte: `N` fijo, `p` creciente. Responde "¿resuelvo este problema más
rápido?". Su límite es Amdahl. Débil: el trabajo **por proceso** se mantiene
constante. Responde "¿puedo resolver problemas más grandes?". En LU el trabajo
crece como `n³`, así que `n` tiene que crecer como `p^(1/3)`, no como `p`.

**D7. ¿Por qué su residuo no cambia con `p`?**

Porque el orden de las operaciones aritméticas es idéntico: cada `a_ij` se
actualiza con la misma secuencia de restas sin importar quién sea su dueño. Lo
único que cambia es qué proceso las ejecuta. Si cambiara con `p`, sería un bug.

**D8. Si la eficiencia con 16 procesos es 40 %, ¿qué significa?**

Que de los 16 procesadores, el equivalente a 6.4 estuvo haciendo trabajo útil
y 9.6 se consumieron en overhead. Es información para decidir: si el cluster
se cobra por core-hora, conviene correr con menos procesos.

---

## Bloque E — Contexto y estado del arte

**E1. ¿Por qué su implementación es más lenta que LAPACK?**

Porque la nuestra es escalar, basada en operaciones vector-matriz (BLAS nivel
2), con razón flops/accesos-a-memoria baja. LAPACK usa un algoritmo **por
bloques**: factoriza un panel angosto y aplica la actualización como producto
matriz-matriz (`dgemm`, BLAS nivel 3), que reutiliza cada dato cargado en cache
`O(b)` veces. La diferencia es de uso de la jerarquía de memoria, no de
cantidad de operaciones.

**E2. ¿Qué cambiarían para correr en cientos de procesos?**

Tres cosas: (1) distribución **2D block-cyclic** como ScaLAPACK — con reparto
1D la difusión del panel involucra a los `p` procesos, mientras que en una
grilla `√p × √p` involucra solo a `√p`; (2) algoritmo **por bloques** para
usar BLAS3; (3) **pivoteo por torneo** (CALU) para reducir el número de
reducciones globales de `n` a `n/b`.

**E3. ¿Qué es CALU?**

*Communication-Avoiding LU*. Usa pivoteo por torneo: en vez de una reducción
global por columna, elige buenos pivotes para un panel entero de `b` columnas
con una sola reducción en forma de torneo. Se demuestra que alcanza las cotas
inferiores teóricas de comunicación para LU (Grigori, Demmel & Xiang, 2011).

**E4. ¿Qué se usa hoy en producción?**

ScaLAPACK sigue siendo el estándar de facto en memoria distribuida. Para
sistemas con GPU, SLATE —desarrollada en el marco del Exascale Computing
Project del Departamento de Energía de EE.UU.— reemplaza el pivoteo parcial
clásico por variantes que mueven menos datos: CALU con pivoteo por torneo,
pivoteo por umbral, transformada de mariposa aleatoria (RBT). El punto en
común es mover la factorización del panel a la GPU para que la CPU deje de ser
el cuello de botella (Gates et al., 2025).

**E5. ¿Por qué LU aparece en el TOP500?**

Porque el benchmark HPL (*High Performance Linpack*), que es el que ordena la
lista, resuelve un sistema denso `Ax = b` por factorización LU con pivoteo
parcial. Es un benchmark de LU por bloques distribuido, y el `Rmax` que se
publica es su rendimiento sostenido en GFLOP/s.

**E6. ¿Dónde aparece LU en problemas reales?**

En el paso lineal interno de casi todo: métodos de elementos finitos, análisis
de circuitos, optimización (el paso de Newton resuelve un sistema lineal),
mínimos cuadrados, métodos implícitos para ecuaciones diferenciales. Es una de
las operaciones más ejecutadas de la computación científica.

---

## Bloque F — Las que descolocan

**F1. ¿Por qué eligieron C y no Python?**

Porque HPC requiere control sobre el layout de memoria y sobre cuándo se
comunica. Python tiene gestión automática de memoria y un intérprete que
agrega overhead por operación; con `mpi4py` los tiempos serían órdenes de
magnitud peores y el análisis de escalabilidad perdería sentido. Además `mpicc`
es lo que hay en cualquier cluster.

**F2. ¿Podrían usar OpenMP en vez de MPI?**

Para un solo nodo, sí, y sería más simple porque hay memoria compartida. Pero
OpenMP no cruza nodos: no serviría para un cluster. Lo ideal a gran escala es
híbrido —MPI entre nodos, OpenMP dentro— pero eso complica mucho la medición y
no era necesario para responder la consigna.

**F3. ¿Cuál fue la parte más difícil?**

*(Contesten con la verdad, no con una frase hecha. Candidatos honestos: el
mapeo global↔local con `p` que no divide a `n`; separar el tiempo de
comunicación del tiempo de espera por desbalance; conseguir que el residuo diera
idéntico para todo `p`.)*

**F4. ¿Qué harían distinto si empezaran de nuevo?**

*(También con la verdad. Buenas respuestas: instrumentar los tiempos desde el
primer día; empezar por el diseño de la distribución de datos antes que por el
código; automatizar la batería de corridas antes de tener acceso al cluster.)*

**F5. ¿Cómo saben que su resultado está bien?**

Tres verificaciones independientes: el residuo `‖Ax−b‖` está en el orden de
`ε` de máquina; el resultado del paralelo coincide **bit a bit** con el del
serial para `p = 1, 2, 3, 4, 5, 7`; y el factor de crecimiento se mantiene
cerca de 1, lo que confirma que el pivoteo está funcionando.

---

## Referencias citadas

- Gates, M., Abdelfattah, A., Akbudak, K., Al Farhan, M., Alomairy, R.,
  Bielich, D., Burgess, T., Cayrols, S., Lindquist, N., Sukkari, D., &
  YarKhan, A. (2025). Evolution of the SLATE linear algebra library.
  *The International Journal of High Performance Computing Applications,
  39*(1), 3–17. https://doi.org/10.1177/10943420241286531
- Grigori, L., Demmel, J. W., & Xiang, H. (2011). CALU: A communication optimal
  LU factorization algorithm. *SIAM Journal on Matrix Analysis and
  Applications, 32*(4). https://doi.org/10.1137/100788926
- Higham, N. J. (2002). *Accuracy and stability of numerical algorithms*
  (2ª ed.). SIAM.
- Golub, G. H., & Van Loan, C. F. (2013). *Matrix computations* (4ª ed.).
  Johns Hopkins University Press.

La lista completa está en [`informe/referencias.md`](../informe/referencias.md).
