# Prompt para redactar el informe final — Proyecto 6: Factorización LU Paralela

> **Cómo usarlo:** copiá todo lo que está debajo de la línea en un chat nuevo y
> adjuntá:
> 1. Las figuras: `figuras/cluster-maia/fig1_speedup.png`,
>    `fig2_eficiencia.png`, `fig6_residuo.png`, y
>    `figuras/comparacion/figC1_cluster.png`, `figC2_pc.png`,
>    `figC3_descomposicion.png`.
> 2. `results/comparacion.md` (todas las tablas).
> 3. Los dos informes de referencia de años anteriores (Monte Carlo y C-14) y
>    `Formato de informe.docx`.
>
> Antes de pegarlo, completá lo que está entre `[CORCHETES]`.

---

Sos un docente de Métodos Numéricos y Computación de Alto Rendimiento que ayuda a dos estudiantes de Ingeniería a redactar el **informe final de su proyecto**. El informe se entrega a la cátedra de **Métodos Numéricos II y Computación Científica** de la **Facultad de Ciencias Exactas y Tecnología (FACET), Universidad Nacional de Tucumán**, en 2026. Después hay una **defensa oral**: todo lo que escribas los estudiantes lo tienen que poder explicar y justificar.

Tu trabajo es redactar el informe **completo, listo para entregar**, en **español rioplatense formal** (registro académico, sin voseo en el texto del informe), siguiendo el formato de la cátedra y el estilo de los informes de referencia que te adjunto. **Usá exclusivamente los datos que figuran en este prompt y en los archivos adjuntos. No inventes números, experimentos, hardware ni referencias.** Si algo te falta, dejalo marcado como `[COMPLETAR: …]` en vez de suponerlo.

---

## 1. Datos de la portada

- **Materia:** Métodos Numéricos II y Computación Científica
- **Carrera:** Ingeniería en Informática `[CONFIRMAR la carrera de Ezequiel]`
- **Integrantes:** Naessens, Maia · Martínez, Ezequiel `[CONFIRMAR apellido]`
- **Título:** *Factorización LU paralela con MPI: escalabilidad del cómputo y la comunicación, y precisión del residuo*
- **Año:** 2026
- **Repositorio (obligatorio en Resultados):** https://github.com/mainaessens/metodos-numericos-II-HPC-LU

---

## 2. Formato obligatorio de la cátedra

Respetá **estas secciones, en este orden, con estos límites de extensión**:

| # | Sección | Extensión | Qué tiene que tener |
|---|---|---|---|
| 1 | **Portada** | — | Materia, carrera, integrantes, título |
| — | **Índice** | — | Numerado con subsecciones, como el informe de C-14 |
| 2 | **Resumen** | **hasta 200 palabras** | De qué trata el trabajo **y los resultados obtenidos, con números** |
| 3 | **Introducción** | **hasta 1 carilla** | Explicación teórica, **estado del arte** y **objetivos** del trabajo |
| 4 | **Datos y Métodos** | **máximo 1½ carillas** | Datos utilizados y métodos empleados |
| 5 | **Resultados** | **máximo 2 carillas** | Gráficas, tablas y sus explicaciones + **link al repositorio** de código y datos |
| 6 | **Conclusiones** | **1 carilla** | Cómo se respondió **cada objetivo** planteado |
| 7 | **Referencias** | — | **Normas APA 7** |

Los límites son estrictos: priorizá lo que se mide y se interpreta por sobre la teoría general. Si una sección se pasa, recortá teoría, no resultados.

---

## 3. Estilo, tomado de los informes de referencia

Los dos informes adjuntos (*Aproximación de π por Monte Carlo en paralelo* y *Estimación de la edad de muestras por C-14*) aprobaron en años anteriores. Imitá estos rasgos:

1. **Secciones y subsecciones numeradas** (3.1, 3.2, 3.2.1…) con índice al principio, como el de C-14.
2. **Tablas con título arriba** del tipo *"Tabla 1: Speedup y eficiencia en el cluster para N = 4000"*. Encabezados claros, unidades entre corchetes y valores con 2 decimales.
3. **Gráficos con título numerado** (*"Gráfico 1: …"*) y, **inmediatamente debajo, un párrafo que lo interprete**. En las referencias se usa el formato "Se puede observar: …" seguido de las observaciones. Ningún gráfico ni tabla queda sin explicar.
4. **Las funciones MPI se explican una por una** (nombre, qué hace y para qué se usa en *este* código), como en el informe de Monte Carlo, pero **solo las que el código realmente usa** (lista en §6.3).
5. **Comparación PC vs cluster**, con análisis de "hasta qué punto conviene paralelizar" y cuál es el número de procesos más conveniente, como en Monte Carlo.
6. **Conclusiones que retoman cada objetivo** y cierran con un bloque *"En términos prácticos:"* en viñetas, como en C-14.
7. **Referencias en una lista numerada en APA.** Si se usó un asistente de IA, se declara como referencia, como hizo el informe de C-14: *Anthropic. (2026). Claude [Asistente conversacional]. Consultado para apoyo en la implementación, la organización de las mediciones y la redacción.*
8. **Tono:** explicativo, de estudiante que entiende lo que hizo. Frases claras y sin relleno ni adjetivos de marketing.

**No copies** las debilidades de las referencias: el informe de Monte Carlo usa el promedio de una sola corrida y llama "Gráfico 4: Tiempo vs Speed Up" a un gráfico ambiguo. Nosotros tenemos 5 repeticiones con mediana y T₁ serial puro, y eso se tiene que notar.

---

## 4. Consigna del proyecto (Proyecto 6)

> Investigar la implementación paralela de la eliminación gaussiana (factorización A = LU). Analizar cómo escalan el cómputo y la comunicación al incrementar N, evaluar el speedup alcanzado y la precisión del error residual ‖Ax − b‖.

**Objetivos del trabajo.** Ponelos numerados en la Introducción y respondelos uno por uno en las Conclusiones:

1. Implementar la factorización LU con pivoteo parcial en versión serial y en versión paralela con MPI (memoria distribuida).
2. Verificar la corrección numérica mediante el residuo relativo ‖Ax − b‖∞ / (‖A‖∞‖x‖∞).
3. Medir la escalabilidad fuerte (speedup y eficiencia) en dos tipos de sistema: computadoras personales y el cluster.
4. Cuantificar por separado cómo escalan el cómputo, la comunicación y el tiempo ocioso por desbalance al aumentar el número de procesos y el tamaño N.
5. Comparar el comportamiento entre sistemas: la reproducibilidad en el cluster y dos notebooks con arquitecturas distintas.

---

## 5. Introducción — contenido técnico

**El problema (breve).** La factorización PA = LU descompone A en una triangular inferior L (con unos en la diagonal) y una triangular superior U, con P una matriz de permutación del pivoteo parcial. Resolver Ax = b se reduce a Ly = Pb (sustitución hacia adelante) y Ux = y (hacia atrás). Se factoriza una vez y se reutiliza para muchos b. Es la base de la resolución de sistemas densos.

**Por qué merece HPC.** El costo es (2/3)N³ operaciones de punto flotante, frente a O(N²) de las sustituciones. Duplicar N multiplica el trabajo por 8. La memoria crece como N²·8 bytes (N = 8000 ocupa 512 MB). Con N = 8000, la versión serial tardó 82 s en el cluster.

**Pivoteo parcial.** En cada paso k se elige como pivote el elemento de mayor módulo de la columna k (desde la fila k hacia abajo). Sin pivoteo, una matriz con un cero en la diagonal falla en el paso 0. Con pivoteo, el método es estable en la práctica (Higham, 2002).

**Estado del arte** (un párrafo, citando):
- **ScaLAPACK** (Blackford et al., 1997): distribución **2D en bloques cíclicos** y algoritmos por bloques que usan BLAS-3.
- **HPL / LINPACK** (Petitet et al.; Dongarra et al., 2003): la factorización LU es el **benchmark del TOP500**, el ranking de supercomputadoras se mide resolviendo exactamente este problema.
- **CALU** (Grigori et al., 2011): pivoteo "por torneo" que reduce la comunicación al mínimo teórico.
- **SLATE** (Gates et al., 2025): sucesor de ScaLAPACK para arquitecturas exascale con GPU.
- **MPI** (Message Passing Interface Forum, 2025): estándar de pasaje de mensajes, versión actual 5.0. Implementación usada: Open MPI 4.1.6.

Aclarar que nuestra implementación es la versión **didáctica** del algoritmo: distribución **1D cíclica por filas**, sin bloques. Sirve para estudiar con claridad la relación entre cómputo y comunicación.

---

## 6. Datos y Métodos — contenido técnico

### 6.1 Datos (matrices de prueba)
Las matrices se generan sintéticamente con un **generador sin estado** (función hash *splitmix64* sobre (semilla, i, j)): cada elemento a_ij depende solo de la semilla y su posición. Esto permite que **cada proceso genere solo sus filas**, sin que ninguno aloque la matriz completa, y garantiza que la versión serial y la paralela resuelvan **exactamente el mismo sistema**. Semilla 42 en todas las corridas. Cuatro tipos:

| Tipo | Para qué sirve |
|---|---|
| Diagonal dominante (`diagdom`) | Bien condicionada. Es la que se usa en todos los experimentos de tiempo |
| Aleatoria (`random`) | Caso general |
| Diagonal nula (`zerodiag`) | Demuestra que el pivoteo es necesario: sin pivoteo falla en el paso k = 0 |
| Hilbert (`hilbert`) | Muy mal condicionada (κ(A) ≈ 10²⁰ para N grande). Muestra que residuo chico ≠ solución exacta |

### 6.2 Algoritmo
Eliminación gaussiana *right-looking* con pivoteo parcial, almacenamiento *in-place* (L y U sobreescriben A) y permutación guardada como vector de enteros.

### 6.3 Paralelización con MPI
- **Distribución cíclica por filas:** la fila i pertenece al proceso i mod p. Como la eliminación "consume" la matriz de arriba hacia abajo, el reparto cíclico mantiene a todos los procesos con trabajo hasta el final: la diferencia de carga entre procesos nunca supera una fila. Con bloques contiguos, los dueños de las primeras filas se quedarían sin trabajo. *(Incluir una figura o tabla chica que muestre el reparto de 8 filas entre 4 procesos.)*
- **El paso k tiene cuatro operaciones.** Explicá cada función MPI en este formato, que es el del informe de Monte Carlo:
  1. **Búsqueda del pivote:** `MPI_Allreduce` con `MPI_MAXLOC` sobre `MPI_DOUBLE_INT`. Cada proceso busca su máximo local en la columna k y la reducción devuelve a todos el valor máximo **y** su fila en una sola operación (comunicación **colectiva**).
  2. **Intercambio de filas:** si las dos filas están en el mismo proceso, es un intercambio local. Si no, `MPI_Sendrecv_replace` entre los dos procesos (comunicación **punto a punto**). Se usa Sendrecv en lugar de Send + Recv para evitar el **deadlock** de dos envíos simétricos.
  3. **Difusión de la fila pivote:** `MPI_Bcast` desde el dueño de la fila k, solo de los elementos k…N−1 (comunicación **colectiva**).
  4. **Actualización:** cada proceso actualiza sus filas por debajo de k. Es 100 % local y sin comunicación, el O(N³/p) del trabajo.
- Otras funciones: `MPI_Init`, `MPI_Comm_rank`, `MPI_Comm_size`, `MPI_Barrier` (largada pareja antes de cronometrar), `MPI_Wtime`, `MPI_Reduce` (junta los tiempos de cada proceso en el rango 0) y `MPI_Finalize`.
- **Modelo de costos:** T_comp ≈ (2/3)N³/p, y T_com ≈ 2N colectivas de O(log p) etapas. Por lo tanto **T_com / T_comp ~ p / N**: la eficiencia mejora al crecer N con p fijo y empeora al crecer p con N fijo. Es la hipótesis que los resultados tienen que verificar.

### 6.4 Métricas
- **Speedup** S(p) = T₁ / T(p), con **T₁ = la versión serial pura** (no la MPI con p = 1, que inflaría el speedup).
- **Eficiencia** E(p) = S(p) / p.
- **GFLOP/s** = (2/3)N³ / T / 10⁹ (convención de HPL).
- **Residuo relativo** r = ‖Ax − b‖∞ / (‖A‖∞ ‖x‖∞). Lo esperable es O(ε) ≈ 10⁻¹⁶.
- **Descomposición del tiempo:** cómputo, comunicación y **ocioso**. Con la opción `--split-idle` se pone una barrera antes de cada colectiva, de modo que la espera por desbalance se mide aparte y no se confunde con comunicación. Esas corridas agregan overhead por las barreras y **no se usan para el speedup**.

### 6.5 Protocolo de medición
5 repeticiones por punto y se reporta la **mediana** (es robusta a corridas lentas por ruido del sistema; el promedio no). Solo se cronometra el lazo de factorización: quedan afuera la generación de la matriz, las sustituciones y el residuo. Serial y MPI se compilan con las mismas banderas (`gcc -O3 -march=native`).

### 6.6 Hardware y software

| | PC 1 (Ezequiel) | PC 2 (Maia) | Cluster |
|---|---|---|---|
| CPU | AMD Ryzen 7 7735HS | Intel Core i7-1355U (13.ª gen.) | AMD EPYC 7B12 (2 sockets × 64 núcleos) |
| Núcleos | 8 núcleos iguales / 16 hilos | **Híbrida:** 2 núcleos de rendimiento (P) + 8 de eficiencia (E), 12 hilos. WSL los expone como 6 núcleos × 2 hilos | 48 núcleos asignados, sin hyperthreading |
| Caché L3 | 16 MiB compartida | 12 MiB compartida | 16 MiB por grupo de 4 núcleos (CCX); 512 MiB en total |
| RAM visible | 6,6 GiB (WSL2) | 7,6 GiB (WSL2) | 192 GB asignados |
| Entorno | Ubuntu 24.04 en WSL2 | Ubuntu en WSL2 | Ubuntu 24.04; entorno web (VS Code) del CCAD-UNC, lab.ccad.unc.edu.ar |
| p usados | 2, 4, 8 | 2, 4, 8 (con `--use-hwthread-cpus`, ver nota) | 2, 4, 8, 16, 32 |
| N usados | 1000, 2000, 4000 | 1000, 2000, 4000 | 1000, 2000, 4000, 8000 |

Compilador gcc 13.3.0 y Open MPI 4.1.6 en todos los casos. **Nota PC 2:** como WSL expone 6 núcleos, Open MPI no permitía lanzar 8 procesos. Se usó `--use-hwthread-cpus` para habilitar los 12 hilos, en todas las corridas de esa máquina para mantener las condiciones iguales.

**Experimentos realizados:**

| Experimento | Dónde | Qué responde |
|---|---|---|
| Validación (residuo, serial vs paralelo) | Todos | ¿El paralelo da lo mismo que el serial? |
| Escalabilidad fuerte | Todos | Speedup y eficiencia |
| Descomposición del tiempo con `--split-idle` (N = 4000) | Cluster (Maia), PC 2 | ¿En qué se va el tiempo? |
| Precisión (4 tipos de matriz, N = 100 a 1000) | Todos | ¿El residuo está en el orden de ε? ¿Depende de p o del hardware? |

El cluster se midió **dos veces de forma independiente** (una corrida de cada integrante, con la misma configuración). **No** se hicieron escalabilidad débil, comparación cíclico vs bloques ni comparación contra LAPACK: van en trabajo futuro, no como resultados.

---

## 7. Resultados — DATOS REALES (usar exactamente estos)

### 7.1 Escalabilidad fuerte en el cluster (corrida de Maia, mediana de 5)

| N | T₁ serial [s] | p = 2 | p = 4 | p = 8 | p = 16 | p = 32 |
|---|---|---|---|---|---|---|
| 1000 | 0,074 | 1,77 (88 %) | 2,65 (66 %) | **3,26 (41 %)** | 2,94 (18 %) | 2,67 (8 %) |
| 2000 | 1,131 | 1,40 (70 %) | 5,63 (141 %) | 8,40 (105 %) | 11,03 (69 %) | 12,94 (40 %) |
| 4000 | 10,441 | 1,23 (62 %) | 4,35 (109 %) | 11,28 (141 %) | 15,33 (96 %) | 25,57 (80 %) |
| 8000 | 82,421 | 1,19 (59 %) | 3,73 (93 %) | 5,65 (71 %) | 10,67 (67 %) | **34,77 (109 %)** |

*Formato: speedup (eficiencia).* Máximo rendimiento: **143,98 GFLOP/s** con N = 8000 y p = 32, contra ~4,1 GFLOP/s del serial con el mismo N.

### 7.2 Reproducibilidad: dos corridas independientes en el cluster

| N | p | S (Ezequiel) | S (Maia) |
|---|---|---|---|
| 2000 | 8 | 8,75 | 8,40 |
| 2000 | 32 | 11,90 | 12,94 |
| 4000 | 8 | 10,87 | 11,28 |
| 4000 | 32 | 24,53 | 25,57 |
| 8000 | 32 | 34,15 | 34,77 |

T₁ serial prácticamente idéntico en las dos corridas: 10,446 s y 10,441 s con N = 4000; 82,83 s y 82,42 s con N = 8000. Para p ≥ 4 las diferencias de speedup son menores al ~10 %, y en la mayoría de los puntos menores al 5 %. **Única discrepancia sistemática: p = 2.** La corrida de Maia dio entre un 10 y un 15 % menos (por ejemplo 1,23 contra 1,40 con N = 4000). En ambas corridas p = 2 queda muy por debajo del resto. Mencionalo con honestidad como variación entre sesiones del cluster (no es un nodo dedicado). No le inventes una causa.

### 7.3 Descomposición del tiempo en el cluster (N = 4000, `--split-idle`, corrida de Maia)

| p | Cómputo | Comunicación | Ocioso |
|---|---|---|---|
| 2 | 98,9 % | 0,8 % | 0,4 % |
| 4 | 92,6 % | 3,9 % | 3,5 % |
| 8 | 76,2 % | 10,9 % | 13,0 % |
| 16 | 53,8 % | 19,3 % | 26,9 % |
| 32 | 36,1 % | 36,1 % | 27,9 % |

Con barreras, T(p = 32) = 0,440 s contra 0,408 s sin barreras: el overhead de `--split-idle` es de ~8 %, lo que justifica no usar esas corridas para el speedup.

### 7.4 Escalabilidad fuerte en las dos notebooks

| N | p | S — PC 1 Ryzen (8 núcleos iguales) | S — PC 2 i7 híbrida |
|---|---|---|---|
| 1000 | 2 | 1,44 | 2,15 |
| 1000 | 4 | 1,79 | 1,84 |
| 1000 | 8 | 2,41 | 1,86 |
| 2000 | 2 | 1,49 | 1,63 |
| 2000 | 4 | 2,10 | 1,61 |
| 2000 | 8 | 2,04 | 2,01 |
| 4000 | 2 | 1,41 | 1,59 |
| 4000 | 4 | 1,77 | **1,28** |
| 4000 | 8 | 1,88 | 1,67 |

T₁ serial: PC 1 = 0,110 / 1,834 / 17,995 s; PC 2 = 0,074 / 1,124 / 11,270 s (N = 1000 / 2000 / 4000). El serial de la PC 2 es más rápido (núcleo P de alta frecuencia), pero escala peor.

### 7.5 Descomposición del tiempo en la PC 2 (N = 4000, `--split-idle`)

| p | Cómputo | Comunicación | Ocioso |
|---|---|---|---|
| 2 | 97,5 % | 0,9 % | 1,6 % |
| 4 | 73,3 % | 1,2 % | **25,5 %** |
| 8 | 82,1 % | 2,5 % | 15,4 % |

### 7.6 Precisión

Residuo relativo (serial, semilla 42), **idéntico en los cuatro sistemas** hasta el último dígito reportado:

| Matriz | N = 100 | N = 200 | N = 500 | N = 1000 |
|---|---|---|---|---|
| Diagonal dominante | 5,03e-16 | 8,32e-16 | 1,66e-15 | 2,42e-15 |
| Aleatoria | 1,14e-16 | 2,37e-16 | 3,80e-16 | 6,81e-16 |
| Diagonal nula | 1,24e-16 | 2,96e-16 | 4,34e-16 | 1,07e-15 |
| Hilbert | 1,60e-18 | 1,29e-18 | 1,54e-18 | 7,03e-19 |

Diagonal dominante con N mayores, en el cluster: 3,38e-15 (N = 2000), 4,68e-15 (N = 4000), 6,82e-15 (N = 8000). El residuo **no cambia con p**: la versión MPI da el mismo valor que la serial. Sin pivoteo, la matriz de diagonal nula **falla en el paso k = 0** (pivote nulo).

---

## 8. Interpretaciones que el informe TIENE que desarrollar

Usá estas explicaciones, que están respaldadas por los datos. Cada una va en el párrafo del gráfico o tabla que corresponde.

1. **Las curvas de speedup se separan por N** (Gráfico de speedup del cluster). Con N = 1000 el máximo es 3,26 en p = 8 y **después baja**: con 125 filas por proceso, las ~2N colectivas por factorización no se achican con p, mientras que el cómputo por proceso sí. Con N = 8000 se llega a 34,77 con p = 32. Es la **verificación experimental de T_com/T_comp ~ p/N**, la conclusión central del trabajo.

2. **Speedup superlineal (E > 100 %) por efecto de caché.** Ocurre cuando el bloque de filas de cada proceso (filas × N × 8 bytes) entra en la L3 de 16 MiB de su grupo de núcleos, mientras que el serial trabaja desde la RAM:

   | Caso | Datos por proceso | Eficiencia |
   |---|---|---|
   | N = 2000, p = 4 | 500 × 2000 × 8 B = 8 MB | 141 % |
   | N = 4000, p = 8 | 500 × 4000 × 8 B = 16 MB | 141 % |
   | N = 8000, p = 32 | 250 × 8000 × 8 B = 16 MB | 109 % |

   Con p = 2 los bloques no entran en ningún caso y el speedup queda en ~1,2–1,5. **No viola nada:** T₁ se midió con una jerarquía de memoria peor. Hay que explicarlo, no esconderlo, y se repitió en las dos corridas.

3. **¿Por qué se aplana el speedup?** Según la descomposición del cluster, con p = 16 más de un cuarto del tiempo (26,9 %) es **espera por desbalance**, no comunicación. Con p = 32, comunicación y cómputo empatan (36 % cada uno). Las tres causas son Amdahl (la búsqueda del pivote y la difusión son seriales en cada paso), el costo de comunicación que crece con p y el desbalance (en el paso k solo trabajan las filas por debajo de k). La figura de descomposición **cuantifica cuánto pesa cada una**. Sin `--split-idle`, la espera quedaba contada como comunicación: es una mejora metodológica respecto de una primera medición.

4. **En las notebooks el límite es la memoria, no la comunicación.** La comunicación es menor al 3 % en la PC 2 (memoria compartida: los mensajes son copias dentro de la RAM), y aun así el speedup no pasa de ~2. Todos los procesos comparten **una sola L3 y el mismo bus de memoria**: agregar procesos no agrega ancho de banda. El serial rinde ~9 GFLOP/s con N = 1000, cuando la matriz de 8 MB entra en la L3, y cae a ~2–4 GFLOP/s con N = 4000 (128 MB).

5. **La anomalía de p = 4 en la PC 2 se explica por hardware heterogéneo.** Con N = 4000, p = 4 (S = 1,28) es **más lento que p = 2 (1,59) y que p = 8 (1,67)**. Se descartó que fuera un error de medición: las 5 repeticiones fueron muy parejas (8,72–8,84 s) y se probaron tres configuraciones de afinidad de procesos (`--bind-to` hwthread, none y la configuración por defecto), que dieron todas ~9 s. La descomposición muestra **25,5 % de tiempo ocioso con p = 4**, contra 1,6 % con p = 2. La explicación es que el i7-1355U es **híbrido** (2 núcleos P rápidos + 8 núcleos E lentos) y tiene un límite de potencia bajo: con p = 2 ambos procesos corren en núcleos P; con p = 4 algunos caen en núcleos E. Como **cada paso k termina en una colectiva**, todos esperan al proceso más lento, y el desbalance ya no es de trabajo sino de **velocidad de los núcleos**. Con p = 8 cada proceso tiene la mitad de filas y la espera relativa baja (15,4 %). La PC 1, con 8 núcleos iguales, no muestra la anomalía. **Es uno de los resultados más interesantes del trabajo:** el reparto cíclico equilibra filas, pero no puede equilibrar núcleos de distinta velocidad. Presentalo como explicación coherente con los datos y con la arquitectura documentada del procesador. No afirmes que se verificó en qué núcleo corrió cada proceso, porque no se midió.

6. **Precisión.** El residuo está en el orden de ε (10⁻¹⁶–10⁻¹⁵) para las matrices bien condicionadas y crece levemente con N, como predice la teoría (∝ N·ε). Es **idéntico entre serial y paralelo, para todo p, y en las cuatro máquinas**, lo que confirma que la paralelización no introduce error y que el resultado depende solo del algoritmo y de la aritmética IEEE 754, no del hardware. **Hilbert:** el residuo es diminuto (~10⁻¹⁸), pero con κ(A) ≈ 10²⁰ la solución calculada no sirve: un residuo chico significa que x es la solución exacta de un problema *cercano*, no que x sea exacta. **Diagonal nula:** se resuelve gracias al pivoteo y sin pivoteo falla en el paso 0.

7. **¿Cuántos procesos conviene usar?** (estilo Monte Carlo). En el cluster depende de N:
   - **N = 1000:** p = 8 da el máximo speedup, y más procesos lo empeoran.
   - **N = 2000:** p = 8, con eficiencia ~105 %.
   - **N = 4000:** p = 16 mantiene ~96 % de eficiencia con S ≈ 15.
   - **N = 8000:** p = 32 sigue siendo eficiente.

   En las notebooks no conviene pasar de p = 2 para N grande: es el mejor equilibrio entre speedup y eficiencia, ~80 % en la PC 2.

---

## 9. Figuras a incluir (archivos adjuntos)

| Gráfico | Archivo | Va en |
|---|---|---|
| Distribución cíclica de filas (dibujarla como tabla o esquema simple) | — | Datos y Métodos |
| Speedup en el cluster | `figuras/cluster-maia/fig1_speedup.png` | Resultados |
| Eficiencia en el cluster | `figuras/cluster-maia/fig2_eficiencia.png` | Resultados |
| Reproducibilidad (dos corridas) | `figuras/comparacion/figC1_cluster.png` | Resultados |
| Descomposición del tiempo, PC 2 vs cluster | `figuras/comparacion/figC3_descomposicion.png` | Resultados |
| Notebooks: homogénea vs híbrida | `figuras/comparacion/figC2_pc.png` | Resultados |
| Residuo vs N por tipo de matriz | `figuras/cluster-maia/fig6_residuo.png` | Resultados |

Las 2 carillas de Resultados son poco espacio. Si no entra todo, priorizá **speedup, descomposición, notebooks y residuo**, y remití el resto al repositorio con una frase ("las figuras restantes y los datos crudos están en el repositorio"), como hizo el informe de C-14.

---

## 10. Conclusiones — qué tienen que decir

Respondé **objetivo por objetivo** (los 5 de §4). Después:
- **El límite encontrado:** con N = 1000 la paralelización deja de convenir a partir de p = 8. La causa, según la descomposición, es que la comunicación y la espera por desbalance superan al cómputo.
- **Lo que se confirmó del modelo:** la separación de las curvas por N verifica T_com/T_comp ~ p/N.
- **El hallazgo de hardware:** en máquinas con núcleos heterogéneos, un algoritmo sincronizado paso a paso queda limitado por el núcleo más lento.
- **Trabajo futuro**, una frase cada uno: distribución 2D en bloques cíclicos (como ScaLAPACK); algoritmo por bloques con BLAS-3 para aprovechar la caché; pivoteo por torneo (CALU) para reducir colectivas; escalabilidad débil; comparación cíclico vs bloques; comparación contra LAPACK/HPL.
- Bloque final **"En términos prácticos:"** con 4 o 5 viñetas.

---

## 11. Referencias (APA 7, numeradas; citarlas en el texto)

1. Blackford, L. S., Choi, J., Cleary, A., D'Azevedo, E., Demmel, J., Dhillon, I., Dongarra, J., Hammarling, S., Henry, G., Petitet, A., Stanley, K., Walker, D., & Whaley, R. C. (1997). *ScaLAPACK users' guide*. Society for Industrial and Applied Mathematics. https://www.netlib.org/scalapack/slug/
2. Burden, R. L., Faires, J. D., & Burden, A. M. (2015). *Numerical analysis* (10.ª ed.). Cengage Learning.
3. Cátedra de Métodos Numéricos II y Computación Científica. (2026). *Computación de alto rendimiento (HPC)* [Material de clase]. Facultad de Ciencias Exactas y Tecnología, Universidad Nacional de Tucumán.
4. Cátedra de Métodos Numéricos II y Computación Científica. (2026). *Trabajo Práctico 1: HPC* [Material de clase]. Facultad de Ciencias Exactas y Tecnología, Universidad Nacional de Tucumán.
5. Dongarra, J. J., Luszczek, P., & Petitet, A. (2003). The LINPACK benchmark: Past, present and future. *Concurrency and Computation: Practice and Experience, 15*(9), 803–820. https://doi.org/10.1002/cpe.728
6. Gates, M., Abdelfattah, A., Akbudak, K., Al Farhan, M., Alomairy, R., Bielich, D., Burgess, T., Cayrols, S., Lindquist, N., Sukkari, D., & YarKhan, A. (2025). Evolution of the SLATE linear algebra library. *The International Journal of High Performance Computing Applications, 39*(1), 3–17. https://doi.org/10.1177/10943420241286531
7. Golub, G. H., & Van Loan, C. F. (2013). *Matrix computations* (4.ª ed.). Johns Hopkins University Press.
8. Grigori, L., Demmel, J. W., & Xiang, H. (2011). CALU: A communication optimal LU factorization algorithm. *SIAM Journal on Matrix Analysis and Applications, 32*(4), 1317–1350. https://doi.org/10.1137/100788926
9. Higham, N. J. (2002). *Accuracy and stability of numerical algorithms* (2.ª ed.). Society for Industrial and Applied Mathematics.
10. Message Passing Interface Forum. (2025). *MPI: A message-passing interface standard* (Versión 5.0). https://www.mpi-forum.org/docs/mpi-5.0/mpi50-report.pdf
11. Naessens, M., & Martínez, E. (2026). *Factorización LU paralela con MPI* [Repositorio de código y datos]. GitHub. https://github.com/mainaessens/metodos-numericos-II-HPC-LU
12. Petitet, A., Whaley, R. C., Dongarra, J., & Cleary, A. (s. f.). *HPL: A portable implementation of the high-performance Linpack benchmark for distributed-memory computers*. Netlib. https://www.netlib.org/benchmark/hpl/
13. The Open MPI Project. (2026). *Open MPI: Open source high performance computing* [Software]. https://www.open-mpi.org/
14. Anthropic. (2026). *Claude* [Asistente conversacional]. Consultado para apoyo en la implementación, la organización de las mediciones y la redacción.

Toda referencia de la lista tiene que estar citada en el texto, y toda cita del texto tiene que estar en la lista.

---

## 12. Entrega

1. Primero mostrame un **esquema** con el índice y, en una línea, qué va en cada sección, para que lo apruebe.
2. Después redactá el informe completo, sección por sección, respetando los límites de extensión. Marcá dónde va cada gráfico como `[Gráfico N: archivo.png]`.
3. Al final agregá:
   - un **conteo aproximado de palabras del Resumen** (tiene que ser ≤ 200);
   - un **checklist**: ¿cada objetivo tiene respuesta en Conclusiones?, ¿cada gráfico o tabla tiene su párrafo?, ¿están todos los `[COMPLETAR]` señalados?, ¿está el link al repositorio?, ¿todas las referencias están citadas?
   - **10 preguntas probables de la defensa oral** con una respuesta corta de 2 o 3 líneas cada una, basadas en este informe. Por ejemplo: por qué T₁ es el serial puro, por qué la mediana, qué es el speedup superlineal, por qué p = 4 es lento en la PC 2, por qué Sendrecv y no Send + Recv, qué significa que el residuo de Hilbert sea chico.
