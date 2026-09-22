# 04 · MPI desde cero

> Prerrequisitos: saber C. Nada de MPI.
> Este documento cubre exactamente lo que necesitamos para el proyecto, ni más
> ni menos. Se apoya en los programas del TP1 de la cátedra (`hola.c`,
> `recibido.c`, `anillo.c`, `broadcast.c`, `scatter_gather.c`).

---

## 1. Qué problema intenta resolver

Una máquina con memoria compartida (tu notebook) tiene un límite: la cantidad
de cores y la RAM de esa máquina. Un **cluster** son muchas máquinas
conectadas por red. Cada una tiene **su propia memoria**, y ninguna puede leer
la memoria de otra.

Ese es el modelo de **memoria distribuida**. Si el proceso A tiene un dato que
el proceso B necesita, alguien tiene que **mandarlo por la red explícitamente**.
No hay variables compartidas. No hay magia.

**MPI** (*Message Passing Interface*) es el estándar para hacer eso. No es un
lenguaje ni un compilador: es una **especificación de biblioteca**. Las
implementaciones más comunes son OpenMPI (la que usamos) y MPICH.

---

## 2. Idea intuitiva: SPMD

MPI usa el modelo **SPMD** (*Single Program, Multiple Data*):

> **Se ejecuta el MISMO programa en todos los procesos, al mismo tiempo.**
> Cada proceso sabe su número, y usa ese número para decidir qué parte del
> trabajo hacer.

Es como repartir el mismo examen a 8 personas diciéndoles "el que tenga el
número 3 resuelve el ejercicio 3". Mismo enunciado, distinta tarea.

Cuando corrés:

```bash
mpirun -n 4 ./mi_programa
```

`mpirun` arranca **cuatro copias** de `mi_programa`. Son cuatro procesos del
sistema operativo, independientes, con memoria separada. Podrían estar en
cuatro máquinas distintas.

---

## 3. El programa mínimo

```c
#include <mpi.h>
#include <stdio.h>

int main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);           /* arranca el entorno MPI */

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);   /* ¿quién soy yo?      */
    MPI_Comm_size(MPI_COMM_WORLD, &size);   /* ¿cuántos somos?     */

    printf("Hola, soy el proceso %d de %d\n", rank, size);

    MPI_Finalize();                   /* cierra el entorno MPI */
    return 0;
}
```

### Los cuatro conceptos que hay que tener claros

| Concepto | Qué es |
|---|---|
| **`MPI_COMM_WORLD`** | El **comunicador**: el grupo de todos los procesos del programa. Es el "universo". Se pueden crear subgrupos, pero no lo necesitamos |
| **`rank`** | El identificador de cada proceso dentro del comunicador. Va de `0` a `size−1`. **Es la variable más importante de todo programa MPI** |
| **`size`** | Cuántos procesos hay. Es el `p` de nuestras fórmulas |
| **`MPI_Init` / `MPI_Finalize`** | Delimitan la zona donde se puede usar MPI. Antes de `Init` y después de `Finalize` no hay MPI |

### La pregunta del TP1: ¿en qué orden sale el `printf`?

Con `p = 5`, la salida puede ser `0 1 2 3 4`, o `4 3 2 1 0`, o
`1 3 2 0 4`, o cualquier otro orden. **No hay ningún orden garantizado.**

Los cinco procesos corren de verdad al mismo tiempo y escriben a la misma
terminal; quién llega primero depende del planificador del sistema operativo,
de la carga de la máquina, de la suerte. Esto se llama **ejecución no
determinística** y es una de las razones por las que depurar programas
paralelos es más difícil.

Si querés salida ordenada, tenés que ordenarla vos: que todos manden su
mensaje al rank 0 y que el rank 0 imprima. O simplemente que imprima solo el
rank 0, que es lo que hace nuestro `lu_mpi.c`.

**Y si el `printf` estuviera debajo de `MPI_Finalize()`?** Se ejecutaría
igual, las mismas 5 veces: `MPI_Finalize` cierra MPI, no termina el proceso.
Pero `rank` y `size` ya fueron leídos antes, así que imprimiría lo mismo. Lo
que **no** podés hacer después de `Finalize` es llamar a cualquier otra
función `MPI_*`.

---

## 4. Comunicación punto a punto

Un proceso manda, otro recibe. Los dos tienen que participar.

```c
MPI_Send(buffer, count, tipo, destino, tag, comunicador);
MPI_Recv(buffer, count, tipo, origen,  tag, comunicador, &status);
```

| Argumento | Qué es |
|---|---|
| `buffer` | puntero a los datos (contiguos en memoria) |
| `count` | cuántos elementos |
| `tipo` | `MPI_DOUBLE`, `MPI_INT`, `MPI_CHAR`, … |
| `destino` / `origen` | el rank del otro |
| `tag` | una etiqueta entera, para distinguir mensajes distintos entre el mismo par |
| `comunicador` | normalmente `MPI_COMM_WORLD` |

**El tipo tiene que coincidir en los dos lados.** Si uno manda `MPI_DOUBLE` y
el otro recibe `MPI_INT`, los bytes se reinterpretan y salen números sin
sentido — sin ningún mensaje de error.

Ejemplo:

```c
if (rank == 0) {
    double dato = 3.14;
    MPI_Send(&dato, 1, MPI_DOUBLE, 1, 0, MPI_COMM_WORLD);
} else if (rank == 1) {
    double recibido;
    MPI_Recv(&recibido, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    printf("Recibi %f\n", recibido);
}
```

---

## 5. Deadlock: el error que vieron en `anillo.c`

Esta sección importa porque en el TP1 ya se toparon con esto, y porque nuestro
`lu_mpi.c` está escrito específicamente para evitarlo.

### El problema

`MPI_Send` es **bloqueante**: no retorna hasta que es seguro reusar el buffer.
Pero hay una sutileza enorme:

- Si el mensaje es **chico**, MPI lo copia a un buffer interno del sistema y
  `MPI_Send` retorna enseguida, sin esperar a que el otro reciba. Esto se
  llama modo **eager**.
- Si el mensaje es **grande** (típicamente más de 64 KB, depende de la
  implementación), no hay buffer interno suficiente: `MPI_Send` **espera** a
  que el receptor haga su `MPI_Recv`. Modo **rendezvous**.

Entonces este código:

```c
MPI_Send(datos, N, MPI_DOUBLE, vecino, 0, MPI_COMM_WORLD);
MPI_Recv(datos, N, MPI_DOUBLE, vecino, 0, MPI_COMM_WORLD, &st);
```

**anda perfecto con `N = 10` y se cuelga para siempre con `N = 1000000`.**

Con `N = 10`: los dos procesos hacen `Send`, los dos retornan inmediatamente
gracias al buffer del sistema, los dos hacen `Recv`, todo bien.

Con `N = 1000000`: los dos procesos hacen `Send` y los dos se quedan esperando
a que el otro reciba. Pero ninguno llega nunca al `Recv`, porque está
bloqueado en el `Send`. **Deadlock.** El programa no da error: simplemente se
queda ahí colgado.

> Esto es exactamente lo que pasa en el punto 1c del TP1 al cambiar
> `#define N 10` por `#define N 1000000`.

### Las tres soluciones

**(a) `MPI_Sendrecv` — la que usamos.**

Manda y recibe en una sola llamada. MPI se encarga de ordenar las operaciones
internamente para que no haya deadlock.

```c
MPI_Sendrecv(envio,   N, MPI_DOUBLE, destino, tag,
             recepcion, N, MPI_DOUBLE, origen,  tag,
             MPI_COMM_WORLD, &status);
```

Y si el buffer de envío y el de recepción son el mismo (un intercambio puro),
existe la variante que usamos en `lu_mpi.c`:

```c
MPI_Sendrecv_replace(fila, n, MPI_DOUBLE, otro, TAG, otro, TAG,
                     MPI_COMM_WORLD, &status);
```

**(b) Alternar el orden por paridad.**

```c
if (rank % 2 == 0) { MPI_Send(...); MPI_Recv(...); }
else               { MPI_Recv(...); MPI_Send(...); }
```

Funciona, pero es frágil: hay que pensarlo de nuevo cada vez.

**(c) Comunicación no bloqueante.**

```c
MPI_Request req[2];
MPI_Irecv(rx, N, MPI_DOUBLE, origen,  0, MPI_COMM_WORLD, &req[0]);
MPI_Isend(tx, N, MPI_DOUBLE, destino, 0, MPI_COMM_WORLD, &req[1]);
/* acá se puede hacer cómputo mientras la comunicación viaja */
MPI_Waitall(2, req, MPI_STATUSES_IGNORE);
```

La `I` es de *immediate*: la función retorna enseguida y la comunicación sigue
en segundo plano. Es la más potente (permite **solapar comunicación y
cómputo**) y la más fácil de romper: no podés tocar el buffer hasta después
del `Wait`.

---

## 6. Comunicación colectiva

Una **colectiva** es una operación en la que participan **todos** los procesos
del comunicador. Todos tienen que llamarla, y todos con los mismos argumentos
de control.

Son más eficientes que hacerlo a mano con `Send`/`Recv` porque la
implementación usa árboles: un broadcast a `p` procesos cuesta `O(log p)`
etapas, no `p`.

### `MPI_Bcast` — uno le manda a todos

```c
MPI_Bcast(buffer, count, tipo, root, comunicador);
```

El proceso `root` tiene el dato; al volver, todos lo tienen.

```
  antes:   P0:[X]  P1:[ ]  P2:[ ]  P3:[ ]
  después: P0:[X]  P1:[X]  P2:[X]  P3:[X]
```

**En nuestro proyecto:** difundir la fila pivote en cada paso `k`.

### `MPI_Reduce` y `MPI_Allreduce` — todos combinan en uno

```c
MPI_Reduce(&local, &global, count, tipo, operacion, root, comunicador);
MPI_Allreduce(&local, &global, count, tipo, operacion, comunicador);
```

Combina un valor de cada proceso con una operación (`MPI_SUM`, `MPI_MAX`,
`MPI_MIN`, `MPI_MAXLOC`, …).

```
  MPI_Reduce con SUM:      MPI_Allreduce con SUM:
  P0:[2] ┐                 P0:[2] ┐        ┌ P0:[10]
  P1:[3] ├→ P0:[10]        P1:[3] ├→ 10 ──┼ P1:[10]
  P2:[1] │                 P2:[1] │        ├ P2:[10]
  P3:[4] ┘                 P3:[4] ┘        └ P3:[10]
```

La diferencia: `Reduce` deja el resultado solo en `root`; `Allreduce` se lo da
a todos. `Allreduce` = `Reduce` + `Bcast`, pero implementado de forma más
eficiente que hacer las dos por separado.

### `MPI_MAXLOC` — la operación estrella del proyecto

Esta merece su propia sección porque es la que resuelve el pivoteo.

Queremos el máximo de `|a_ik|` sobre todas las filas, **pero también necesitamos
saber qué fila lo tiene** (para intercambiarla). Con `MPI_MAX` obtendríamos
solo el valor, y haría falta una segunda comunicación para averiguar el dueño.

`MPI_MAXLOC` devuelve el par (valor, índice) de una sola vez:

```c
struct { double val; int idx; } local, global;

local.val = mi_maximo;
local.idx = indice_global_de_ese_maximo;

MPI_Allreduce(&local, &global, 1, MPI_DOUBLE_INT, MPI_MAXLOC, MPI_COMM_WORLD);

/* ahora global.val es el máximo de todos
   y global.idx dice de qué fila salió     */
```

El tipo `MPI_DOUBLE_INT` corresponde exactamente a esa struct de
`{double; int;}`. En caso de empate, `MAXLOC` devuelve el índice **menor**, lo
que hace que el algoritmo sea determinista.

### `MPI_Scatter` y `MPI_Gather` — repartir y juntar

```
  Scatter:  P0:[a b c d] → P0:[a] P1:[b] P2:[c] P3:[d]
  Gather:   P0:[a] P1:[b] P2:[c] P3:[d] → P0:[a b c d]
```

Las variantes `Scatterv` / `Gatherv` (con `v` de *variable*) permiten que cada
proceso reciba o mande una cantidad **distinta** de elementos. Son las que hay
que usar cuando `n` no es múltiplo de `p`.

> **En nuestro `lu_mpi.c` NO usamos Scatter.** Cada proceso **genera** sus
> propias filas con `mat_generate_row()`. Así nunca existe una copia completa
> de la matriz en ningún proceso, ni siquiera un instante. Es la decisión que
> permite correr con `n = 8000`.

### `MPI_Barrier` — sincronizar

```c
MPI_Barrier(MPI_COMM_WORLD);
```

Ningún proceso pasa hasta que llegaron todos. **No transmite datos.** Lo
usamos para dos cosas:

1. Arrancar el cronómetro parejo en todos los procesos.
2. Con `--split-idle`, separar el tiempo de espera por desbalance del tiempo
   de comunicación real.

**No lo pongas "por las dudas":** cada barrera cuesta, y de más obliga a todos
a esperar al más lento sin necesidad.

---

## 7. Medición de tiempos

```c
double t0 = MPI_Wtime();
/* ... lo que se quiere medir ... */
double t1 = MPI_Wtime();
double dt = t1 - t0;
```

`MPI_Wtime()` devuelve segundos de reloj de pared, en doble precisión. Es
portable y tiene buena resolución (`MPI_Wtick()` te dice cuál).

### Dos reglas

**(1) Barrera antes de arrancar.** Si no, los procesos que arrancaron primero
miden tiempo de espera de los que arrancaron después:

```c
MPI_Barrier(MPI_COMM_WORLD);
double t0 = MPI_Wtime();
```

**(2) El tiempo del programa es el del proceso MÁS LENTO.** Un programa
paralelo termina cuando termina el último:

```c
double dt_max;
MPI_Allreduce(&dt, &dt_max, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
```

Reportar el tiempo del rank 0, o el promedio, es un error clásico que infla
los resultados.

---

## 8. Compilar y ejecutar

```bash
mpicc -O3 -Wall programa.c -o programa -lm     # mpicc es un wrapper de gcc
mpirun -n 4 ./programa                          # 4 procesos
```

Si pedís más procesos que cores físicos, OpenMPI se queja. Para desarrollo:

```bash
mpirun --oversubscribe -n 8 ./programa
```

Los resultados numéricos van a ser correctos; los tiempos, basura.

Para ver cuántos cores tenés: `nproc` o `lscpu`.

---

## 9. Errores frecuentes

| Error | Síntoma | Solución |
|---|---|---|
| `MPI_Send` simétrico con mensajes grandes | El programa se cuelga sin error | `MPI_Sendrecv` |
| Una colectiva llamada solo por algunos procesos | Cuelgue, o error críptico | **Todos** llaman a la colectiva, siempre |
| Tipos distintos en `Send` y `Recv` | Números sin sentido, sin error | Verificar que coincidan |
| Usar el buffer de un `Isend` antes del `Wait` | Datos corruptos aleatorios | No tocarlo hasta después del `Wait` |
| Suponer un orden en la salida de `printf` | Salida desordenada | Que imprima solo el rank 0 |
| Medir con el reloj del rank 0 | Tiempos optimistas | `Allreduce` con `MPI_MAX` |
| Olvidar `MPI_Finalize()` | Warnings o cuelgues al terminar | Ponerlo siempre |
| Suponer que `n` es múltiplo de `p` | Se rompe con p = 3, 5, 7 | Usar `Gatherv`, y probar con primos |

---

## 10. Preguntas típicas de parcial y defensa

1. ¿Qué es el rank y para qué sirve?
2. Diferencia entre comunicación punto a punto y colectiva. Dar un ejemplo de
   cada una en nuestro código.
3. Explicá por qué `anillo.c` anda con N=10 y se cuelga con N=1000000.
4. ¿Cuál es la diferencia entre `MPI_Reduce` y `MPI_Allreduce`?
5. ¿Por qué usamos `MPI_MAXLOC` y no `MPI_MAX` para el pivoteo?
6. ¿Por qué un `MPI_Bcast` a `p` procesos cuesta `O(log p)` y no `O(p)`?
7. ¿Por qué hay que hacer `Allreduce` con `MPI_MAX` sobre los tiempos?
8. ¿Qué hace `MPI_Barrier` y cuándo conviene usarlo?

---

**Siguiente:** [05 · LU en paralelo](05-lu-en-paralelo.md)
