# 08 · Guía del cluster

> Este documento es una guía práctica. Como todavía no conocemos el cluster de
> la cátedra, está escrito como una **lista de cosas que hay que averiguar** y
> los comandos para averiguarlas. Actualícenlo apenas tengan el acceso: va a
> ser la sección de "hardware utilizado" del informe.

---

## 1. El primer día: qué averiguar

Correr estos comandos y **guardar la salida en un archivo**. Después la
necesitan para el informe.

### ¿Qué hardware hay?

```bash
lscpu                       # modelo de CPU, cores, sockets, caches
free -h                     # memoria RAM
nproc                       # cores visibles en este nodo
cat /proc/cpuinfo | grep "model name" | head -1
```

Lo que hay que anotar para el informe:

- Modelo exacto de CPU y frecuencia
- Cores físicos por nodo (¡ojo con el hyperthreading: `Thread(s) per core`)
- Sockets por nodo
- Tamaño de L3
- RAM por nodo
- Cantidad de nodos disponibles
- Tipo de red entre nodos (Ethernet 1G / 10G / InfiniBand). Esto cambia
  **muchísimo** los resultados de comunicación

### ¿Qué software hay?

```bash
which gcc mpicc mpirun
gcc --version
mpirun --version            # OpenMPI o MPICH, y qué versión
module avail                # si el cluster usa "environment modules"
```

Si aparece una lista con `module avail`, casi seguro hay que cargar el
compilador y MPI antes de usarlos:

```bash
module load gcc
module load openmpi
module list                 # ver qué está cargado
```

Anotar los nombres exactos: van en `scripts/job.slurm`.

### ¿Cómo se lanzan trabajos?

```bash
which sbatch squeue sinfo   # SLURM
which qsub qstat            # PBS / Torque
```

Si es SLURM (lo más común):

```bash
sinfo                       # particiones, nodos, estado
sinfo -o "%P %D %c %m %l"   # partición, nodos, cores, memoria, tiempo límite
scontrol show partition     # detalle de límites
squeue -u $USER             # mis trabajos en cola
```

Anotar: nombre de la partición, límite de tiempo por job, límite de nodos por
usuario, si se puede pedir `--exclusive`.

---

## 2. Preguntas para hacerle a la cátedra

Vale la pena mandar un mail con esto el primer día:

1. ¿Cuántos nodos y cuántos cores por nodo podemos usar?
2. ¿Qué gestor de colas hay y cuál es el límite de tiempo por job?
3. ¿Se pueden pedir nodos **exclusivos**? *(Importa mucho: sin exclusividad,
   otro trabajo corriendo en el mismo nodo contamina las mediciones)*
4. ¿Qué red conecta los nodos? *(Ethernet vs InfiniBand cambia el análisis de
   comunicación por completo)*
5. ¿Hay cuota de disco? *(Nuestros CSV son chicos, pero conviene saber)*
6. ¿Hay LAPACK / BLAS optimizado instalado, para el experimento E7?

---

## 3. Flujo de trabajo

### Subir el código

Lo más limpio es clonar el repositorio directamente en el cluster:

```bash
ssh usuario@cluster
git clone https://github.com/mainaessens/mnii-hpc-lu.git
cd mnii-hpc-lu
```

Así lo que corre en el cluster es exactamente lo que está versionado, sin
copias sueltas que se desincronizan. Si el cluster no tiene salida a internet:

```bash
# desde tu máquina
rsync -avz --exclude bin --exclude build . usuario@cluster:~/mnii-hpc-lu/
```

### Compilar

```bash
module load gcc openmpi      # si hace falta
make clean && make
```

Si `-march=native` da problemas:

```bash
make ARCH=                   # compila sin optimizaciones específicas de CPU
```

### Probar antes de lanzar la batería completa

**Nunca** correr `run_experiments.sh all` de entrada. Primero:

```bash
make test                                    # validación numérica
mpirun -n 4 ./bin/lu_mpi -n 1000             # una corrida chica
```

Verificar que el residuo esté en `~1e-15`. Después sí, el job completo.

### Lanzar el job

Editar `scripts/job.slurm` con los valores reales del cluster, y:

```bash
sbatch scripts/job.slurm
squeue -u $USER              # ver el estado
tail -f logs/lu-<jobid>.out  # seguir la salida
scancel <jobid>              # cancelar si hace falta
```

### Traer los resultados

```bash
# desde tu máquina
scp usuario@cluster:~/mnii-hpc-lu/results/raw.csv results/
python3 scripts/plots.py
```

O, mejor, hacer commit desde el cluster:

```bash
git add results/raw.csv
git commit -m "resultados: corrida del <fecha>, <p> procesos"
git push
```

Así los datos quedan versionados junto con el código que los produjo, que es
justo lo que la consigna pide entregar.

---

## 4. Reglas de convivencia

Un cluster es un recurso compartido. Estas no son formalidades:

- **No correr nada pesado en el nodo de login.** El nodo de login es para
  editar, compilar y lanzar jobs. Correr ahí un `mpirun -n 16` molesta a todo
  el mundo y en algunos clusters es causa de suspensión de la cuenta.
- **Pedir solo lo que se va a usar.** Un job que reserva 32 cores y usa 4
  bloquea a los demás.
- **Estimar bien el tiempo.** Si pedís 24 h para algo que tarda 20 minutos, el
  planificador te va a poner al final de la cola.
- **Probar chico antes de correr grande.** Un bug que se manifiesta después de
  2 horas de cómputo cuesta 2 horas de cluster a todos.

---

## 5. Problemas típicos y qué hacer

| Síntoma | Causa probable | Solución |
|---|---|---|
| `mpicc: command not found` | MPI no cargado | `module load openmpi` |
| `Illegal instruction` en el nodo de cómputo | `-march=native` compilado en un nodo distinto | Compilar dentro del job, o `make ARCH=` |
| `There are not enough slots available` | Pedís más procesos que cores asignados | Ajustar `-n`, o `--oversubscribe` (solo para probar) |
| El job queda en `PD` (pending) mucho tiempo | Cola llena, o pedís más recursos de los que hay | `squeue --start` para ver la estimación; pedir menos |
| El job muere con `OOM` | Memoria insuficiente | Bajar `n`, o subir `--mem` en el sbatch |
| El programa se cuelga | Deadlock MPI | Ver el [doc 04](04-mpi-desde-cero.md) §5 |
| Tiempos con enorme variabilidad entre repeticiones | Nodo compartido con otro job | Pedir `--exclusive`; si no se puede, subir el número de repeticiones y **reportarlo** |
| Speedup ridículo con muchos procesos | Estás corriendo con `--oversubscribe` | Usar solo tantos procesos como cores reales |

---

## 6. Qué va en el informe

La sección de "Datos y Métodos" tiene que incluir la descripción del hardware.
Algo como:

> Las mediciones se realizaron en el cluster \<nombre\> de la FACET-UNT,
> compuesto por \<N\> nodos con \<modelo de CPU\> (\<cores\> cores por nodo,
> \<frecuencia\> GHz, \<L3\> MB de caché L3) y \<RAM\> GB de memoria por nodo,
> interconectados por \<red\>. El código se compiló con \<compilador y
> versión\> usando \<banderas\> y \<implementación MPI y versión\>. Los
> trabajos se lanzaron mediante \<gestor de colas\> \<con/sin\> asignación
> exclusiva de nodos. Cada punto experimental corresponde a la mediana de
> \<R\> repeticiones.

Completen los corchetes con los datos reales del primer día. **Un informe de
HPC sin descripción del hardware no se puede evaluar**: los tiempos no
significan nada sin saber en qué corrieron.

---

**Siguiente:** [09 · Errores frecuentes](09-errores-frecuentes.md)
