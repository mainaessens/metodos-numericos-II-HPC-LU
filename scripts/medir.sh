#!/usr/bin/env bash
# ============================================================================
#  medir.sh — Corre TODA la batería acordada en un sistema y deja todo listo
#             para comparar: datos crudos, info de la máquina, resumen y figuras.
#
#  Proyecto 6 · Factorización LU Paralela · FACET-UNT 2026
#
#  Uso:
#      ./scripts/medir.sh <sistema> <pc|cluster>
#
#  Ejemplos:
#      ./scripts/medir.sh pc-maia      pc        # en la notebook (WSL)
#      ./scripts/medir.sh cluster-maia cluster   # en el cluster (CCAD)
#
#  Qué corre (mismos parámetros que las corridas de Ezequiel + E4):
#      E2  escalabilidad fuerte  (serial + MPI, diagdom, cíclico, seed 42, 5 reps)
#      E4  descomposición del tiempo con --split-idle (N = 4000)
#      E6  precisión (4 tipos de matriz, N = 100..1000, serial y p = 4)
#
#  Deja:
#      results/<sistema>/raw.csv        todas las corridas
#      results/<sistema>/info_<sistema>.txt   hardware y software
#      results/<sistema>/corrida.log    log completo de la corrida
#      results/<sistema>/resumen.md     (lo genera plots.py)
#      figuras/<sistema>/*.png          (lo genera plots.py)
#
#  Los parámetros se pueden pisar desde el entorno, p.ej.:
#      PROCS="2 4 8" REPS=7 ./scripts/medir.sh pc-maia pc
# ============================================================================

set -u

SIS=${1:-}
TIPO=${2:-}
if [ -z "$SIS" ] || { [ "$TIPO" != "pc" ] && [ "$TIPO" != "cluster" ]; }; then
    echo "uso: $0 <sistema> <pc|cluster>      ej: $0 pc-maia pc"
    exit 1
fi

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
DIR="$ROOT/results/$SIS"
mkdir -p "$DIR"

# --- parámetros: iguales a los de Ezequiel ----------------------------------
if [ "$TIPO" = "pc" ]; then
    export PROCS=${PROCS:-"2 4 8"}
    export N_STRONG=${N_STRONG:-"1000 2000 4000"}
else
    export PROCS=${PROCS:-"2 4 8 16 32"}
    export N_STRONG=${N_STRONG:-"1000 2000 4000 8000"}
fi
export REPS=${REPS:-5}
export SEED=${SEED:-42}
export N_BREAK=${N_BREAK:-4000}
export OUT="$DIR/raw.csv"

if [ -s "$OUT" ]; then
    echo "ATENCION: $OUT ya existe y las corridas nuevas se AGREGARIAN."
    echo "Si es una corrida nueva, renombralo o borralo antes. Abortando."
    exit 1
fi

LOG="$DIR/corrida.log"
exec > >(tee -a "$LOG") 2>&1

echo "=============================================="
echo " Sistema   : $SIS ($TIPO)"
echo " PROCS     : $PROCS"
echo " N_STRONG  : $N_STRONG"
echo " N_BREAK   : $N_BREAK   REPS: $REPS   SEED: $SEED"
echo " Inicio    : $(date)"
echo "=============================================="

# --- 1. info de la máquina (va al informe) -----------------------------------
INFO="$DIR/info_$SIS.txt"
{
    echo "# Sistema: $SIS ($TIPO) — $(date)"
    echo; echo "## lscpu";    lscpu
    echo; echo "## memoria";  free -h
    echo; echo "## nproc";    nproc
    echo; echo "## compilador y MPI"
    gcc --version | head -1
    mpirun --version 2>/dev/null | head -1
    echo; echo "## sistema operativo"
    (grep PRETTY_NAME /etc/os-release || uname -a) 2>/dev/null
    uname -r
} > "$INFO" 2>&1
echo "Info de la máquina en $INFO"

# --- 2. compilar limpio (mismas banderas para serial y MPI) ------------------
(cd "$ROOT" && make clean >/dev/null && make) || { echo "Fallo make"; exit 1; }

# --- 3. validación rápida antes de gastar tiempo -----------------------------
echo "--- validación (N=1000, p=4) ---"
mpirun ${MPIFLAGS:-} -n 4 "$ROOT/bin/lu_mpi" -n 1000 -k diagdom -s "$SEED" | grep -E "residuo|OK|ATENCION"

# --- 4. experimentos ---------------------------------------------------------
"$ROOT/scripts/run_experiments.sh" strong
"$ROOT/scripts/run_experiments.sh" breakdown
"$ROOT/scripts/run_experiments.sh" precision

# --- 5. figuras y resumen ----------------------------------------------------
echo "--- figuras ---"
python3 "$ROOT/scripts/plots.py" "$SIS" || \
    echo "(plots.py falló: faltan pandas/matplotlib. Los datos están en $OUT; se puede graficar en otra máquina)"

echo "=============================================="
echo " Fin: $(date)"
echo " Datos   : results/$SIS/"
echo " Figuras : figuras/$SIS/"
echo "=============================================="
