# Factorización LU paralela (A = LU)

Proyecto 6 de Métodos Numéricos II / Computación Científica (FACET-UNT).
Maia Naessens y Ezequiel Martínez.

> **Consigna.** La factorización LU descompone una matriz cuadrada invertible A en el
> producto de una matriz triangular inferior L y una triangular superior U, constituyendo
> la base del cálculo de sistemas de ecuaciones lineales Ax = b. Investigue la
> implementación paralela de la eliminación gaussiana. Analice cómo escala el tiempo de
> cómputo y la comunicación al incrementar N, evaluando el speedup y la precisión numérica
> del error residual ‖Ax − b‖.

Hay solo dos programas, cada uno en un único archivo `.c`, sin bibliotecas propias ni `.h`:

| Archivo | Qué hace |
|---|---|
| `LU_Serial.c` | Factoriza A = LU en un solo proceso. Es la referencia para el speedup. |
| `LU_MPI.c` | Hace la misma factorización repartiendo las filas entre P procesos con MPI. |

Los dos generan **la misma matriz** (misma semilla), así que tienen que dar el mismo
residuo y el mismo error. Esa es la primera forma de comprobar que la versión paralela está bien.

## Dónde se cumple cada parte de la consigna

| Lo que pide la consigna | Dónde está |
|---|---|
| A = L·U con A cuadrada invertible | A es diagonalmente dominante ⇒ invertible. `factorizar_lu` en `LU_Serial.c`; el bucle de `k` en `LU_MPI.c`. |
| Base del cálculo de Ax = b | `resolver`: sustitución hacia adelante (Ly = b) y hacia atrás (Ux = y), en los dos archivos. |
| Implementación paralela de la eliminación gaussiana | `LU_MPI.c`: reparto cíclico de filas + `MPI_Bcast` de la fila pivote en cada paso. |
| Cómo escala el **tiempo de cómputo** al incrementar N | Los dos programas imprimen *Tiempo de cómputo* por separado. |
| Cómo escala la **comunicación** al incrementar N | `LU_MPI.c` acumula aparte el tiempo de `MPI_Scatter`, `MPI_Bcast` y `MPI_Gather` e imprime *Tiempo de comunicación*. |
| **Speedup** | Tiempo total de `LU_Serial` ÷ tiempo total de `LU_MPI` con P procesos. Los dos usan reloj de pared, así que son comparables. |
| **Precisión del error residual ‖Ax − b‖** | Los dos imprimen el residuo, el residuo relativo ‖Ax − b‖/‖b‖ y el error máximo \|x − 1\|. |

## Qué hacen los programas

1. Generan una matriz A de N × N **diagonalmente dominante** (a cada elemento de la
   diagonal se le suma N) y un vector b = A · (1, 1, …, 1). Así la solución exacta es x = 1
   y se puede medir el error directamente.
2. Factorizan A = LU por eliminación gaussiana **sin pivoteo**. Como A es
   diagonalmente dominante, no hace falta pivotear: nunca aparece un pivote cero o muy chico.
   L y U se guardan en la misma matriz (L debajo de la diagonal, con 1 en la diagonal; U en
   la diagonal y arriba), así no se usa memoria de más.
3. Resuelven L y = b (hacia adelante) y U x = y (hacia atrás).
4. Muestran los tiempos, el residuo ‖Ax − b‖, el residuo relativo y el error máximo |x − 1|.

El costo de la factorización es aproximadamente **(2/3)·N³** operaciones; por eso al duplicar
N el tiempo de cómputo se multiplica por ~8.

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
  (solo las columnas k…N−1, lo anterior ya no se usa), y cada proceso actualiza sus filas
  que están debajo de k.
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
./LU_Serial 1200
mpirun -np 4 ./LU_MPI 1200
```

- **N tiene que ser múltiplo de la cantidad de procesos** (por ejemplo N = 1200 sirve
  para 1, 2, 3, 4, 6 y 8 procesos). Si no, el programa avisa y termina.
- Si `mpirun` dice que no hay suficientes "slots" (pediste más procesos que núcleos),
  agregá `--oversubscribe`. Ojo: con más procesos que núcleos los tiempos no sirven
  para medir speedup, porque los procesos se pelean por el mismo procesador.

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
LU MPI  N = 1200  procesos = 2
Tiempo de computo:        0.177995 s
Tiempo de comunicacion:   0.039311 s
Tiempo total:             0.199276 s
Residuo ||Ax - b||:       1.150686e-10
Residuo relativo:         1.845492e-15
Error maximo |x - 1|:     9.325873e-15
CSV,1200,2,0.177995,0.039311,0.199276,1.150686e-10,1.845492e-15
```

- **Tiempo de cómputo:** solo las cuentas (los multiplicadores y la actualización de las filas).
- **Tiempo de comunicación:** `MPI_Scatter` + todos los `MPI_Bcast` + `MPI_Gather`.
- **Tiempo total:** desde que arrancan todos los procesos hasta que termina el `Gather`.
- **Residuo relativo:** es el que conviene comparar entre distintos N, porque el residuo
  absoluto crece con N aunque la precisión sea la misma.
- La última línea (`CSV,…`) tiene los mismos números separados por comas, en el orden
  **N, P, cómputo, comunicación, total, residuo, residuo relativo**, para copiar a una planilla.

Dos aclaraciones para la defensa:

1. **Cómputo + comunicación no da exactamente el total**, porque cada máximo puede venir
   de un proceso distinto.
2. **El tiempo de `MPI_Bcast` incluye la espera**: un proceso que llega antes se queda
   esperando a los demás. Por eso ese tiempo mezcla comunicación real con desbalance de carga,
   y crece cuando hay más procesos que núcleos.

### Qué se espera al analizar los resultados

- **Cómputo:** crece como N³ y se reparte entre los procesos, así que baja aproximadamente
  a la mitad al duplicar P.
- **Comunicación:** crece como N² (se mandan filas, no la matriz entera) y **aumenta** con P.
- Por eso el speedup mejora hasta cierto punto y después se achata: llega un momento en que
  la comunicación pesa más que lo que se gana repartiendo las cuentas.
- **Speedup:** S(P) = tiempo total de `LU_Serial` / tiempo total de `LU_MPI` con P procesos.
  **Eficiencia:** E(P) = S(P) / P.
- **Precisión:** el residuo relativo tiene que quedar cerca de 1e-15 y **no depender de P**.
  Si cambia al cambiar la cantidad de procesos, hay un error en la versión paralela.
