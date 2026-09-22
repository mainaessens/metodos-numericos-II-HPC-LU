# Guía de estudio — Factorización LU Paralela

**Proyecto 6 · Métodos Numéricos II y Computación Científica · FACET-UNT · 2026**

Esta carpeta es el material de estudio del proyecto. Está escrita asumiendo
que **no sabés qué es la factorización LU** y que **nunca escribiste un
programa con MPI**. Cada documento arranca de cero y va subiendo.

---

## Cómo usar esta guía

Leelos en orden. Cada uno depende del anterior.

| # | Documento | Qué vas a poder hacer al terminarlo |
|---|---|---|
| 01 | [¿Qué es la factorización LU?](01-que-es-la-factorizacion-lu.md) | Factorizar una matriz 3×3 a mano y explicar para qué sirve |
| 02 | [Pivoteo y estabilidad numérica](02-pivoteo-y-estabilidad.md) | Explicar por qué sin pivoteo el método se rompe, con un ejemplo concreto |
| 03 | [Costo y por qué paralelizar](03-costo-y-por-que-paralelizar.md) | Contar operaciones y justificar por qué este problema merece HPC |
| 04 | [MPI desde cero](04-mpi-desde-cero.md) | Entender rank, size, Send/Recv y las colectivas que usamos |
| 05 | [LU en paralelo](05-lu-en-paralelo.md) | Explicar el algoritmo distribuido y el modelo de costos |
| 06 | [El código, explicado](06-el-codigo-explicado.md) | Leer `lu_serial.c` y `lu_mpi.c` línea por línea |
| 07 | [Experimentos y métricas](07-experimentos-y-metricas.md) | Diseñar y ejecutar las mediciones, y calcular speedup y eficiencia |
| 08 | [Guía del cluster](08-guia-del-cluster.md) | Compilar y lanzar jobs en el servidor de la cátedra |
| 09 | [Errores frecuentes](09-errores-frecuentes.md) | Evitar las trampas que hacen perder días |
| 10 | [Preguntas de defensa](10-preguntas-de-defensa.md) | Contestar cualquier cosa que les pregunten en la exposición |

---

## El proyecto en tres frases

La **factorización LU** descompone una matriz cuadrada invertible `A` en el
producto de una triangular inferior `L` por una triangular superior `U`. Es la
base del método directo para resolver sistemas de ecuaciones lineales `Ax = b`,
y el 99 % de su costo está en una sola operación que se repite `n` veces y que
se puede repartir entre varios procesos.

Nuestro trabajo consiste en implementar esa factorización en paralelo con
**MPI**, medir cómo escalan el cómputo y la comunicación al crecer `N`, y
evaluar el **speedup** y la **precisión del residuo `‖Ax − b‖`**.

---

## Mapa del repositorio

```
mnii-hpc-lu/
├── include/lu.h            declaraciones comunes
├── src/
│   ├── matgen.c            generación reproducible de matrices de prueba
│   ├── verify.c            normas, residuo, factor de crecimiento, reloj
│   ├── lu_serial.c         LU serial con pivoteo parcial  (la referencia)
│   └── lu_mpi.c            LU paralela con MPI            (el proyecto)
├── scripts/
│   ├── run_experiments.sh  batería completa de corridas
│   ├── job.slurm           plantilla de job para el cluster
│   └── plots.py            estadística y figuras del informe
├── docs/                   ESTA carpeta
├── results/                CSV crudos de las corridas (se versionan)
├── figuras/                las figuras que van al informe
└── informe/                borrador del informe
```

---

## Arranque rápido

```bash
make                  # compila bin/lu_serial y bin/lu_mpi
make test             # validación: serial vs MPI con 1,2,3,4,5,7 procesos

./bin/lu_serial -n 1000
mpirun -n 4 ./bin/lu_mpi -n 1000
```

Si no tenés 4 cores libres, agregá `--oversubscribe` a `mpirun`. Los tiempos
no van a servir para nada, pero el resultado numérico sí.

---

## Una advertencia sobre cómo estudiar esto

En la defensa van a preguntar **por qué**, no **qué**. "Usamos distribución
cíclica" no es una respuesta; "usamos distribución cíclica porque la
eliminación gaussiana consume la matriz de arriba hacia abajo y con bloques
contiguos los primeros procesos quedan ociosos a mitad de camino" sí lo es.

Cada documento de esta carpeta tiene al final una sección de **preguntas
típicas**. Usenlas.
