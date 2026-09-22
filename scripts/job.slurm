#!/bin/bash
# ============================================================================
#  job.slurm — Script de job para SLURM
#
#  Proyecto 6 · Factorización LU Paralela · FACET-UNT 2026
#
#  ATENCIÓN: este archivo es una PLANTILLA. Antes de usarlo hay que ajustar
#  el nombre de la partición, la cantidad de nodos/tareas y el tiempo límite
#  a lo que realmente tenga el cluster de la cátedra. Preguntar:
#
#      sinfo                      -> qué particiones hay y cuántos nodos
#      scontrol show partition    -> límites de tiempo por partición
#      lscpu                      -> cores por nodo
#      module avail               -> qué compiladores y MPI hay disponibles
#
#  Enviar el job:     sbatch scripts/job.slurm
#  Ver la cola:       squeue -u $USER
#  Cancelar:          scancel <jobid>
#  Ver la salida:     tail -f logs/lu-<jobid>.out
# ============================================================================

#SBATCH --job-name=lu-hpc
#SBATCH --output=logs/lu-%j.out
#SBATCH --error=logs/lu-%j.err
#SBATCH --nodes=1                 # <-- ajustar
#SBATCH --ntasks=16               # <-- ajustar: cantidad máxima de procesos MPI
#SBATCH --cpus-per-task=1
#SBATCH --time=02:00:00           # <-- ajustar al límite de la partición
# #SBATCH --partition=NOMBRE      # <-- descomentar y poner el nombre real
# #SBATCH --exclusive             # <-- si el cluster lo permite: nodo dedicado.
#                                 #     Mejora MUCHÍSIMO la reproducibilidad de
#                                 #     los tiempos. Vale la pena pedirlo.

set -e
mkdir -p logs results

echo "=============================================="
echo " Job        : $SLURM_JOB_ID"
echo " Nodos      : $SLURM_JOB_NUM_NODES"
echo " Tareas     : $SLURM_NTASKS"
echo " Lista nodos: $SLURM_JOB_NODELIST"
echo " Fecha      : $(date)"
echo "=============================================="

# --- Entorno ---------------------------------------------------------------
# Muchos clusters usan "modules". Descomentar y ajustar a lo que muestre
# `module avail`. Si el cluster no usa modules, borrar estas líneas.
# module load gcc
# module load openmpi

echo "--- entorno ---"
which mpicc mpirun gcc || true
mpirun --version | head -1 || true
lscpu | grep -E "Model name|^CPU\(s\)|Thread|Core|Socket|L3" || true
echo "---------------"

# --- Compilación -----------------------------------------------------------
# Compilamos DENTRO del job para que el binario se genere en el mismo tipo de
# nodo donde va a correr. Si el nodo de login tiene otra CPU que los de
# cómputo, -march=native compilado afuera puede generar instrucciones
# ilegales ("Illegal instruction") en el nodo de cómputo.
make clean
make

# --- Experimentos ----------------------------------------------------------
# PROCS tiene que ser una lista de valores <= SLURM_NTASKS.
export PROCS="1 2 4 8 16"
export REPS=5
export N_STRONG="1000 2000 4000"

# En SLURM, `srun` conoce la asignación del job y no hace falta pasarle -n.
# Pero nuestro script usa mpirun -n <p> para variar p, así que dejamos mpirun.
export MPIRUN="mpirun"
export MPIFLAGS=""

./scripts/run_experiments.sh all

echo "=============================================="
echo " Terminado: $(date)"
echo " Resultados en results/raw.csv"
echo "=============================================="
