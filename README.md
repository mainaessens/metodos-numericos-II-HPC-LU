# Factorización LU Paralela con MPI

**Proyecto 6 · Métodos Numéricos II y Computación Científica**
Facultad de Ciencias Exactas y Tecnología — Universidad Nacional de Tucumán — 2026

---

## El problema

La factorización LU descompone una matriz cuadrada invertible `A` en el
producto de una triangular inferior `L` por una triangular superior `U`,
constituyendo la base del cálculo de sistemas de ecuaciones lineales `Ax = b`.

Este proyecto implementa la **eliminación gaussiana con pivoteo parcial en
paralelo** sobre memoria distribuida (MPI), y analiza:

- cómo escalan el **tiempo de cómputo** y el de **comunicación** al incrementar N,
- el **speedup** y la **eficiencia** alcanzados,
- la **precisión numérica** medida por el residuo `‖Ax − b‖`.

---

## Arranque rápido

```bash
make                              # compila bin/lu_serial y bin/lu_mpi
make test                         # validación: serial vs MPI con p = 1,2,3,4,5,7

./bin/lu_serial -n 2000
mpirun -n 8 ./bin/lu_mpi -n 2000
```

Requiere `gcc`, una implementación de MPI (`mpicc` / `mpirun`), y —para las
figuras— Python con `pandas` y `matplotlib`.

### Batería completa de experimentos

```bash
./scripts/run_experiments.sh all      # produce results/raw.csv
python3 scripts/plots.py              # produce figuras/ y results/resumen.md
```

En el cluster, lanzar desde un job: `sbatch scripts/job.slurm`
(editar antes con los valores reales del cluster).

---

## Opciones

### `bin/lu_serial`

| Opción | Qué hace |
|---|---|
| `-n <int>` | dimensión de la matriz (default 1000) |
| `-k <tipo>` | `random` · `diagdom` · `hilbert` · `zerodiag` (default `diagdom`) |
| `-s <int>` | semilla (default 42) |
| `-r <int>` | repeticiones cronometradas |
| `--nopivot` | factorizar **sin** pivoteo (experimento de estabilidad) |
| `--csv` | una línea CSV en vez del reporte legible |

### `bin/lu_mpi`

Las mismas, más:

| Opción | Qué hace |
|---|---|
| `-d <reparto>` | `cyclic` (default) · `block` |
| `--split-idle` | barreras para medir el tiempo ocioso por separado |

---

## Estructura

```
├── include/lu.h            declaraciones comunes
├── src/
│   ├── matgen.c            generación reproducible de matrices
│   ├── verify.c            normas, residuo, factor de crecimiento, reloj
│   ├── lu_serial.c         LU serial con pivoteo parcial (la referencia)
│   └── lu_mpi.c            LU paralela con MPI (el proyecto)
├── scripts/
│   ├── run_experiments.sh  batería completa de corridas
│   ├── job.slurm           plantilla de job para el cluster
│   └── plots.py            estadística y figuras
├── docs/                   ← GUÍA DE ESTUDIO, empezar por acá
├── results/                CSV crudos de las corridas
├── figuras/                figuras del informe
└── informe/                estructura y referencias
```

---

## Documentación

La carpeta [`docs/`](docs/) es una guía de estudio completa, escrita asumiendo
que no sabés qué es la factorización LU ni cómo se programa con MPI.

| # | Documento |
|---|---|
| 01 | [¿Qué es la factorización LU?](docs/01-que-es-la-factorizacion-lu.md) |
| 02 | [Pivoteo y estabilidad numérica](docs/02-pivoteo-y-estabilidad.md) |
| 03 | [Costo y por qué paralelizar](docs/03-costo-y-por-que-paralelizar.md) |
| 04 | [MPI desde cero](docs/04-mpi-desde-cero.md) |
| 05 | [LU en paralelo](docs/05-lu-en-paralelo.md) |
| 06 | [El código, explicado](docs/06-el-codigo-explicado.md) |
| 07 | [Experimentos y métricas](docs/07-experimentos-y-metricas.md) |
| 08 | [Guía del cluster](docs/08-guia-del-cluster.md) |
| 09 | [Errores frecuentes](docs/09-errores-frecuentes.md) |
| 10 | [Preguntas de defensa](docs/10-preguntas-de-defensa.md) |

---

## Decisiones de diseño

| Decisión | Elección | Por qué |
|---|---|---|
| Distribución de datos | **Cíclica por filas** | La eliminación consume la matriz de arriba hacia abajo; con bloques contiguos la eficiencia máxima alcanzable es 50 % |
| Variante del algoritmo | **Right-looking** (kij) | Expone la actualización de la submatriz como un bloque grande de trabajo independiente |
| Pivoteo | **Parcial** | Estabilidad a costo `O(n²)`; el completo cuesta `O(n³)` |
| Búsqueda del pivote | **`MPI_Allreduce` + `MPI_MAXLOC`** | Nadie tiene la columna entera; MAXLOC devuelve valor e índice en una colectiva |
| Difusión del panel | **`MPI_Bcast`**, solo desde la columna `k` | Colectiva optimizada; difundir desde `k` reduce el volumen a la mitad |
| Intercambio de filas | **`MPI_Sendrecv_replace`** | `Send` + `Recv` simétricos son deadlock con mensajes grandes |
| Almacenamiento | **In-place** sobre A | Mitad de memoria; la diagonal de L es de unos y no se guarda |
| Layout de memoria | **1D contiguo** `A[i*n+j]` | MPI necesita buffers contiguos; además respeta la localidad de cache |
| Generación de matrices | **Hash posicional sin estado** | Cada proceso genera sus filas sin que nadie aloque `n²`, y serial y paralelo resuelven el mismo sistema |

El razonamiento completo está en [`docs/05`](docs/05-lu-en-paralelo.md).

---

## Verificación

```bash
make test
```

El residuo relativo `‖Ax−b‖∞ / (‖A‖∞‖x‖∞)` del paralelo debe ser **idéntico**
al del serial para cualquier número de procesos, incluidos los que no dividen
a `n`:

```
serial : 1.232598e-15
p = 1  : 1.232598e-15
p = 2  : 1.232598e-15
p = 3  : 1.232598e-15
p = 4  : 1.232598e-15
p = 5  : 1.232598e-15
p = 7  : 1.232598e-15
```

Y la demostración de por qué el pivoteo no es opcional:

```bash
./bin/lu_serial -n 200 -k zerodiag              # residuo ~1e-16
./bin/lu_serial -n 200 -k zerodiag --nopivot    # FALLO: pivote nulo en k=0
```

---

## Equipo

- *(completar)*
- *(completar)*
- *(completar)*

---

## Referencias

En [`informe/referencias.md`](informe/referencias.md).
