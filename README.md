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
| `LU_Serial.c` | Factoriza P·A = L·U en un solo proceso. Es la referencia para el speedup. |
| `LU_MPI.c` | Hace la misma factorización repartiendo las filas entre P procesos con MPI. |

Los dos generan **la misma matriz y el mismo vector** (misma semilla), así que tienen que
dar **exactamente el mismo residuo**. Esa es la primera forma de comprobar que la versión
paralela está bien.

## Dónde se cumple cada parte de la consigna

| Lo que pide la consigna | Dónde está |
|---|---|
| A = L·U con A cuadrada invertible | `factorizar_lu` en `LU_Serial.c`; el bucle de `k` en `LU_MPI.c`. A se genera al azar: es invertible con probabilidad 1. |
| Base del cálculo de Ax = b | `resolver`: sustitución hacia adelante (Ly = Pb) y hacia atrás (Ux = y), en los dos archivos. |
| Implementación paralela de la eliminación gaussiana | `LU_MPI.c`: reparto cíclico de filas + pivoteo parcial distribuido + `MPI_Bcast` de la fila pivote en cada paso. |
| Cómo escala el **tiempo de cómputo** al incrementar N | Los dos programas imprimen *Tiempo de cómputo* por separado. |
| Cómo escala la **comunicación** al incrementar N | `LU_MPI.c` acumula aparte el tiempo de todas las llamadas MPI e imprime *Tiempo de comunicación*. |
| **Speedup** | Tiempo total de `LU_Serial` ÷ tiempo total de `LU_MPI` con P procesos. Los dos usan reloj de pared, así que son comparables. |
| **Precisión del error residual ‖Ax − b‖** | Los dos imprimen el residuo ‖Ax − b‖ y el residuo relativo ‖Ax − b‖/‖b‖. |

## Qué hacen los programas

1. Generan **al azar** una matriz A de N × N y un vector b, con valores uniformes
   en [−1, 1] y semilla fija.
2. Factorizan **P·A = L·U** por eliminación gaussiana **con pivoteo parcial**.
   L y U se guardan en la misma matriz (L debajo de la diagonal, con 1 en la diagonal;
   U en la diagonal y arriba), así no se usa memoria de más. La permutación P se guarda
   como un vector de índices `perm`, no como matriz.
3. Resuelven L y = P·b (hacia adelante) y U x = y (hacia atrás).
4. Muestran los tiempos, el residuo ‖Ax − b‖ y el residuo relativo.

El costo de la factorización es aproximadamente **(2/3)·N³** operaciones; por eso al duplicar
N el tiempo de cómputo se multiplica por ~8.

### Por qué hace falta pivotear

Como A se genera al azar, **no** es diagonalmente dominante: en la columna k puede aparecer
un pivote muy chico (o cero), y dividir por él amplifica los errores de redondeo.

El **pivoteo parcial** busca, en la columna k de la fila k para abajo, el elemento de mayor
valor absoluto, e intercambia esa fila con la fila k. Así todos los multiplicadores cumplen
|l_ik| ≤ 1 y el error no se amplifica. Medido sobre la misma matriz, con N = 1200:

| | Residuo relativo ‖Ax − b‖/‖b‖ |
|---|---|
| Sin pivoteo | 5.85 × 10⁻¹⁰ |
| **Con pivoteo parcial** | **2.55 × 10⁻¹²** |

Es decir, unas **230 veces más preciso** por el solo hecho de elegir bien el pivote.

### Cómo se paraleliza (`LU_MPI.c`)

- **Reparto cíclico por filas:** la fila i la tiene el proceso `i % P`.
  Con P = 4 y N = 8:

  | Proceso | Filas globales |
  |---|---|
  | 0 | 0, 4 |
  | 1 | 1, 5 |
  | 2 | 2, 6 |
  | 3 | 3, 7 |

  Se usa el reparto cíclico y no por bloques para **balancear la carga**: a medida que k
  avanza, las filas de arriba ya terminaron. Con bloques, los primeros procesos se quedarían
  sin trabajo; con el reparto cíclico todos siguen teniendo filas por debajo de k.

- El proceso 0 genera A, la ordena y la reparte con `MPI_Scatter`.
- En **cada paso k** se hacen cuatro cosas:

  1. **Búsqueda local del pivote** (cómputo): cada proceso busca, entre *sus* filas i ≥ k,
     la de mayor |a_ik|.
  2. **Elección global del pivote** (comunicación): el proceso 0 junta los P candidatos con
     `MPI_Gather`, elige el mayor de todos y avisa a todos cuál es la fila pivote con
     `MPI_Bcast`.
  3. **Intercambio de filas** (comunicación): la fila k y la fila pivote se intercambian
     enteras. Si están en el mismo proceso es un intercambio local; si están en procesos
     distintos, los dos dueños se las mandan con `MPI_Send` / `MPI_Recv`. Para que no se
     queden los dos esperando a la vez (*deadlock*), el de rank menor manda primero y
     después recibe, y el otro al revés.
  4. **Eliminación**: el dueño de la fila k la difunde con `MPI_Bcast` (solo las columnas
     k…N−1, lo anterior ya no se usa) y cada proceso actualiza sus filas que están debajo de k.

- Al final, el proceso 0 junta todo con `MPI_Gather`, aplica la permutación a b, resuelve el
  sistema y calcula el residuo contra la A **original**.
- Los tiempos se miden con `MPI_Wtime`; el proceso 0 los junta con `MPI_Gather` y se queda
  con el del proceso más lento.

Funciones MPI usadas (todas vistas en la teoría / TP1): `MPI_Init`, `MPI_Comm_rank`,
`MPI_Comm_size`, `MPI_Finalize`, `MPI_Wtime`, `MPI_Barrier`, `MPI_Scatter`, `MPI_Bcast`,
`MPI_Gather`, `MPI_Send`, `MPI_Recv`.

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
LU MPI (con pivoteo parcial)  N = 1200  procesos = 2
Tiempo de computo:        0.112541 s
Tiempo de comunicacion:   0.022785 s
Tiempo total:             0.134512 s
Residuo ||Ax - b||:       5.125055e-11
Residuo relativo:         2.553602e-12
CSV,1200,2,0.112541,0.022785,0.134512,5.125055e-11,2.553602e-12
```

- **Tiempo de cómputo:** las cuentas (búsqueda local del pivote, multiplicadores y
  actualización de las filas).
- **Tiempo de comunicación:** `MPI_Scatter` + las colectivas del pivoteo + los
  `MPI_Send`/`MPI_Recv` de los intercambios + los `MPI_Bcast` de la fila pivote + `MPI_Gather`.
- **Tiempo total:** desde que arrancan todos los procesos hasta que termina el `Gather`.
- **Residuo relativo:** es el que conviene comparar entre distintos N, porque el residuo
  absoluto crece con N aunque la precisión sea la misma.
- La última línea (`CSV,…`) tiene los mismos números separados por comas, en el orden
  **N, P, cómputo, comunicación, total, residuo, residuo relativo**, para copiar a una planilla.

Dos aclaraciones para la defensa:

1. **Cómputo + comunicación no da exactamente el total**, porque cada máximo puede venir
   de un proceso distinto.
2. **El tiempo de las colectivas incluye la espera**: un proceso que llega antes se queda
   esperando a los demás. Por eso ese tiempo mezcla comunicación real con desbalance de carga,
   y crece cuando hay más procesos que núcleos.

### Qué se espera al analizar los resultados

- **Cómputo:** crece como N³ y se reparte entre los procesos, así que baja aproximadamente
  a la mitad al duplicar P.
- **Comunicación:** crece como N² en volumen (se mandan filas, no la matriz entera), pero
  además **el pivoteo agrega dos colectivas chiquitas por paso** (`Gather` + `Bcast`), o sea
  ~2N colectivas en total. Esas son de *latencia*, no de volumen: casi no mandan datos pero
  cuestan tiempo, y ese costo **crece con P**. Es el precio de la estabilidad numérica.
- Por eso el speedup mejora hasta cierto punto y después se achata: llega un momento en que
  la comunicación pesa más que lo que se gana repartiendo las cuentas.
- **Speedup:** S(P) = tiempo total de `LU_Serial` / tiempo total de `LU_MPI` con P procesos.
  **Eficiencia:** E(P) = S(P) / P.
- **Precisión:** el residuo relativo **no tiene que depender de P**. Como los dos programas
  generan la misma matriz y eligen los mismos pivotes, el residuo de `LU_MPI` tiene que dar
  **idéntico** al de `LU_Serial` para cualquier cantidad de procesos. Si cambia al cambiar P,
  hay un error en la versión paralela.
