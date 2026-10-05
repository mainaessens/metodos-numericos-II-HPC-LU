# Factorización LU paralela (A = LU) — versión simple

Proyecto 6 de Métodos Numéricos II / Computación Científica (FACET-UNT).
Maia Naessens y Ezequiel Martínez.

Hay solo dos programas, cada uno en un único archivo `.c`, sin bibliotecas propias ni `.h`:

| Archivo | Qué hace |
|---|---|
| `LU_Serial.c` | Factoriza A = LU en un solo proceso. Es la referencia para comparar. |
| `LU_MPI.c` | Hace la misma factorización repartiendo las filas entre P procesos con MPI. |

Los dos generan **la misma matriz** (misma semilla), así que tienen que dar el mismo
residuo y el mismo error. Esa es la primera forma de comprobar que la versión paralela está bien.

## Qué hacen los programas

1. Generan una matriz A de N × N **diagonalmente dominante** (a cada elemento de la
   diagonal se le suma N) y un vector b = A · (1, 1, …, 1). Así la solución exacta es x = 1.
2. Factorizan A = LU por eliminación gaussiana **sin pivoteo**. Como A es
   diagonalmente dominante, no hace falta pivotear: nunca aparece un pivote cero o muy chico.
   L y U se guardan en la misma matriz (L debajo de la diagonal, con 1 en la diagonal; U en
   la diagonal y arriba).
3. Resuelven L y = b (hacia adelante) y U x = y (hacia atrás).
4. Muestran el tiempo de la factorización, el residuo ‖Ax − b‖ y el error máximo |x − 1|.

### Cómo se paraleliza (`LU_MPI.c`)

- **Reparto cíclico por filas:** la fila i la tiene el proceso `i % P`.
  Con P = 4 y N = 8:

  | Proceso | Filas globales |
  |---|---|
  | 0 | 0, 4 |
  | 1 | 1, 5 |
  | 2 | 2, 6 |
  | 3 | 3, 7 |

- El proceso 0 genera A, la ordena y la reparte con `MPI_Scatter`.
- En cada paso k, el dueño de la fila k (proceso `k % P`) la envía a todos con `MPI_Bcast`
  (solo las columnas k…N−1), y cada proceso actualiza sus filas que están debajo de k.
- Al final, el proceso 0 junta todo con `MPI_Gather`, resuelve el sistema y verifica.
- Los tiempos se miden con `MPI_Wtime` y se toma el del proceso más lento con `MPI_Reduce` (`MPI_MAX`).

Se usa el reparto cíclico y no por bloques para **balancear la carga**: a medida que k avanza,
las filas de arriba ya terminaron. Con bloques, los primeros procesos se quedarían sin
trabajo; con el reparto cíclico todos siguen teniendo filas por debajo de k.

Funciones MPI usadas (todas vistas en la teoría o en el TP1): `MPI_Init`, `MPI_Comm_rank`,
`MPI_Comm_size`, `MPI_Finalize`, `MPI_Wtime`, `MPI_Barrier`, `MPI_Scatter`, `MPI_Bcast`,
`MPI_Gather`, `MPI_Reduce`.

## Cómo compilar y correr en la PC

Hace falta Linux o WSL (Ubuntu) con `gcc` y OpenMPI. Si no están instalados:

```bash
sudo apt update
sudo apt install build-essential openmpi-bin libopenmpi-dev
```

Compilar (desde la carpeta del repo):

```bash
gcc -O2 -o LU_Serial LU_Serial.c -lm
mpicc -O2 -o LU_MPI LU_MPI.c -lm
```

Correr (el número es N, el tamaño de la matriz; si no se pone, usa N = 1000):

```bash
./LU_Serial 1000
mpirun -np 4 ./LU_MPI 1000
```

- **N tiene que ser múltiplo de la cantidad de procesos** (por ejemplo N = 1200 sirve
  para 1, 2, 3, 4, 6 y 8 procesos). Si no, el programa avisa y termina.
- Si `mpirun` dice que no hay suficientes "slots" (pediste más procesos que núcleos),
  agregá `--oversubscribe`: `mpirun --oversubscribe -np 8 ./LU_MPI 1200`.

## Cómo correr en el cluster (CCAD-UNC)

1. Entrar a `lab.ccad.unc.edu.ar` y abrir una terminal.
2. Subir `LU_Serial.c` y `LU_MPI.c` (o clonar el repo con `git clone`).
3. Si `mpicc` no se encuentra, cargar los módulos: `module load gcc openmpi`
   (con `module avail` se ve qué hay disponible).
4. Compilar y correr igual que en la PC:

```bash
gcc -O2 -o LU_Serial LU_Serial.c -lm
mpicc -O2 -o LU_MPI LU_MPI.c -lm
./LU_Serial 2000
mpirun -np 8 ./LU_MPI 2000
```

## Cómo leer la salida

```
LU MPI  N = 1200  procesos = 4
Tiempo de factorizacion: 0.267761 s
Tiempo total (con Scatter y Gather): 0.274157 s
Residuo ||Ax - b||:      1.150686e-10
Error maximo |x - 1|:    9.325873e-15
```

- **Tiempo de factorización:** solo el bucle de eliminación (cálculo + `MPI_Bcast` de cada fila pivote).
- **Tiempo total:** agrega el reparto inicial (`MPI_Scatter`) y la recolección final (`MPI_Gather`).
  La diferencia entre los dos es el costo de mover la matriz.
- **Residuo y error:** tienen que ser chicos (del orden de 1e-10 o menos) e **iguales** a los
  de `LU_Serial` con el mismo N. Si dan distinto, la versión paralela tiene un error.

Speedup: S(P) = tiempo de `LU_Serial` / tiempo de `LU_MPI` con P procesos.
