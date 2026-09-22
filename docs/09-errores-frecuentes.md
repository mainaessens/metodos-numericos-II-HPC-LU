# 09 · Errores frecuentes

> Lista consolidada de las trampas de este proyecto, ordenadas por lo que
> cuesta descubrirlas tarde. Conviene leerla entera una vez antes de empezar a
> programar, y volver a ella cuando algo no anda.

---

## Nivel 1 — Los que arruinan el proyecto

### 1. Deadlock con `MPI_Send` bloqueante

**Síntoma:** el programa anda con `n = 500` y se cuelga con `n = 5000`. Sin
mensaje de error, sin crash: simplemente se queda ahí.

**Causa:** `MPI_Send` copia a un buffer interno solo si el mensaje es chico.
Con mensajes grandes espera al receptor, y si los dos procesos hacen `Send`
primero, ninguno llega nunca al `Recv`.

Es exactamente lo que vieron en `anillo.c` del TP1 al pasar `N` de 10 a
1 000 000.

**Solución:** `MPI_Sendrecv` o `MPI_Sendrecv_replace`. Nunca `Send` + `Recv`
simétricos.

> **No confíen en que "anduvo con n chico".** Este bug está latente hasta el
> día de la corrida grande.

### 2. Replicar la matriz completa en cada proceso

**Síntoma:** anda con `n = 1000`, se queda sin memoria con `n = 8000`, o el
sistema empieza a hacer swap y los tiempos se vuelven absurdos.

**Causa:** alocar `n²` doubles por proceso, o hacer `MPI_Scatter` desde el
rank 0 (que necesita `n²` para tener qué repartir).

**Solución:** cada proceso aloca `⌈n/p⌉ · n` y **genera** sus propias filas.
Nuestro `mat_generate_row()` existe exactamente para esto.

Esta decisión hay que tomarla **el primer día**. Cambiarla después implica
reescribir el programa entero.

### 3. Medir `T₁` con la versión MPI corrida con un proceso

**Síntoma:** el speedup da sospechosamente bueno.

**Causa:** el numerador y el denominador arrastran el mismo overhead, así que
se cancela. Estás midiendo el programa contra sí mismo.

**Solución:** `T₁` es `bin/lu_serial`, compilado con las mismas banderas.

> Es lo primero que revisa un docente al leer un informe de HPC.

### 4. No separar `T_cómputo` de `T_comunicación`

**Síntoma:** cuando llega el momento de escribir el informe, no se puede
responder la parte de la consigna que pregunta cómo escala la comunicación.

**Causa:** instrumentar solo el tiempo total.

**Solución:** instrumentar desde la Fase 2. Reinstrumentar después de haber
hecho las corridas significa **rehacerlas todas**, y eso son días de cluster.

---

## Nivel 2 — Los que dan resultados mal sin avisar

### 5. Intercambiar solo la parte derecha de las filas

**Síntoma:** el residuo da `10⁻²` en vez de `10⁻¹⁵`. Ningún error, ningún
crash.

**Causa:** al intercambiar las filas `k` y `p`, se intercambia solo desde la
columna `k` en adelante. Pero la parte izquierda contiene multiplicadores de
`L` que también tienen que viajar: `L` corresponde a la matriz **permutada**.

**Solución:** intercambiar la fila completa, `j` desde 0.

### 6. Buscar el pivote sin valor absoluto

**Síntoma:** residuo degradado, sobre todo con matrices con entradas negativas
grandes.

**Causa:** comparar `a_ik > max` en vez de `fabs(a_ik) > max`. Un `−100` es
un pivote excelente y un `+0.001` es pésimo.

**Solución:** `fabs()` siempre.

### 7. Olvidar aplicar la permutación a `b`

**Síntoma:** la solución está mal, pero L y U están bien.

**Causa:** `PA = LU`, entonces `LUx = Pb`. Si resolvés `LUx = b` estás
resolviendo un sistema con las ecuaciones en otro orden.

**Solución:** aplicar los mismos intercambios a `b`, en el mismo orden.

### 8. Suponer que `n` es múltiplo de `p`

**Síntoma:** anda perfecto con `p = 1, 2, 4, 8` y se rompe con `p = 3, 5, 7`.
O peor: no se rompe, pero da mal.

**Causa:** `nloc = n/p` en vez de la fórmula correcta; usar `MPI_Gather` donde
hace falta `MPI_Gatherv`.

**Solución:** `nloc = n/p + (rank < n % p ? 1 : 0)`, separar `nloc` de
`maxloc`, y **probar con `p` primo siempre**. `make test` corre 1, 2, 3, 4, 5, 7
por esta razón.

### 9. Usar `MPI_MAX` en vez de `MPI_MAXLOC`

**Síntoma:** funciona pero necesitás una segunda comunicación para averiguar
qué proceso tenía el máximo, o peor, elegís mal la fila.

**Solución:** `MPI_MAXLOC` con `MPI_DOUBLE_INT` devuelve valor e índice juntos.

### 10. Difundir la fila pivote antes de intercambiarla

**Síntoma:** el residuo cambia con `p`.

**Causa:** orden de operaciones. Primero se elige el pivote, después se
intercambia, y **recién después** se difunde la fila que quedó en la posición
`k`.

---

## Nivel 3 — Los que arruinan las mediciones

### 11. Una sola corrida por punto

Las curvas quedan con dientes de sierra que no significan nada. Mínimo 5
repeticiones, y se reporta la **mediana**.

### 12. Incluir la generación de la matriz en el tiempo medido

Distorsiona el speedup, y cada vez más a medida que crece `n`. El cronómetro
arranca después de un `MPI_Barrier` y para antes de escribir resultados.

### 13. Banderas de compilación distintas entre serial y MPI

Si el serial va con `-O3` y el MPI con `-O0`, el speedup es ficción. El
`Makefile` usa la misma variable para los dos: no la toquen a mano.

### 14. Reportar el tiempo del rank 0

Un programa paralelo termina cuando termina el **último** proceso.
`MPI_Allreduce` con `MPI_MAX` sobre los tiempos.

### 15. Medir con matrices chicas

Con `n = 500` y `p = 16` hay 31 filas por proceso: se está midiendo latencia
de red, no el algoritmo. No es un bug, es el régimen equivocado. Para conclusiones
sobre escalabilidad hacen falta al menos `n = 2000`.

### 16. Usar las corridas con `--split-idle` para el speedup

Las barreras agregan overhead. Esas corridas sirven solo para el gráfico de
composición del tiempo.

### 17. Escalabilidad débil con `n = n₀·p`

El trabajo de LU crece como `n³`. Con `n = n₀·p`, el trabajo por proceso se
multiplica por `p²` y la curva da un desastre que no significa nada. Lo
correcto es `n = n₀·p^(1/3)`.

---

## Nivel 4 — Los de C que cuestan una tarde

### 18. Desbordamiento de enteros en los índices

Con `n = 50000`, `i*n` pasa de 2³¹ y desborda el `int`. Segfault, o peor,
corrupción silenciosa. Siempre `(size_t)i * n`.

### 19. `double **A` en vez de `double *A`

MPI necesita buffers **contiguos**. Con `double**` cada fila puede estar en
cualquier lado del heap, y además se pierde la localidad de cache.

### 20. Dividir dentro del lazo interno

`1.0/a_kk` se calcula una vez por paso, afuera. Una división cuesta 4–5 veces
lo que una multiplicación.

### 21. `clock_gettime` no compila con `-std=c11`

Es POSIX, no C estándar. `#define _POSIX_C_SOURCE 199309L` **antes** de los
`#include`. (Ya está en `verify.c`.)

### 22. Una colectiva llamada solo por algunos procesos

`MPI_Bcast`, `MPI_Allreduce`, `MPI_Barrier`: **todos** los procesos del
comunicador tienen que llamarlas. Un `if (rank == 0) MPI_Bcast(...)` es un
cuelgue garantizado.

---

## Nivel 5 — Los de gestión del proyecto

### 23. Dividir el código en tres pedazos, uno por persona

El lazo principal de LU no se puede partir en tres. Repartan por **rol**
(núcleo numérico / paralelización / infraestructura y análisis), no por
líneas de código.

### 24. Dejar el informe para el final

La introducción y la sección de métodos se pueden escribir en las primeras dos
semanas. Solo resultados y conclusiones dependen de las mediciones.

### 25. Que cada uno entienda solo su parte

En la exposición oral pueden preguntarle a cualquiera sobre cualquier cosa.
Ensayo cruzado: cada persona contesta preguntas de lo que **no** programó.

### 26. Conseguir el acceso al cluster tarde

Sin acceso no hay mediciones, y sin mediciones no hay informe. Es lo primero
que hay que destrabar.

---

**Siguiente:** [10 · Preguntas de defensa](10-preguntas-de-defensa.md)
