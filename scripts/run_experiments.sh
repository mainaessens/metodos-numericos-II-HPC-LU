#!/usr/bin/env bash
# ============================================================================
#  run_experiments.sh — Batería completa de experimentos del proyecto
#
#  Proyecto 6 · Factorización LU Paralela · FACET-UNT 2026
#
#  Produce results/raw.csv con TODAS las corridas. Cada punto se repite
#  REPS veces y se guardan TODAS las repeticiones: la estadística (mediana)
#  la hace después scripts/plots.py, que es donde se puede revisar.
#
#  Uso:
#      ./scripts/run_experiments.sh                  # todo
#      ./scripts/run_experiments.sh strong           # solo escalabilidad fuerte
#      PROCS="1 2 4 8" ./scripts/run_experiments.sh  # sobreescribir procesos
#
#  En el cluster conviene lanzarlo desde un job (ver scripts/job.slurm),
#  no en el nodo de login.
# ============================================================================

set -u

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BIN="$ROOT/bin"
OUT="$ROOT/results/raw.csv"

# --- parámetros (se pueden pisar desde el entorno) --------------------------
REPS=${REPS:-5}                       # repeticiones por punto: nunca menos de 5
SEED=${SEED:-42}
PROCS=${PROCS:-"1 2 4 8 16"}          # ajustar a los cores del cluster
N_STRONG=${N_STRONG:-"1000 2000 4000"}
MPIRUN=${MPIRUN:-mpirun}
MPIFLAGS=${MPIFLAGS:-""}              # p.ej. "--oversubscribe" en una notebook

WHICH=${1:-all}

# --- preparación ------------------------------------------------------------
mkdir -p "$ROOT/results"
if [ ! -x "$BIN/lu_mpi" ]; then
    echo "No estan los binarios. Corriendo make..."
    (cd "$ROOT" && make) || exit 1
fi

# Cabecera del CSV. Si el archivo ya existe se AGREGA, no se pisa:
# así se pueden acumular corridas de días distintos.
if [ ! -f "$OUT" ]; then
    echo "impl,n,p,dist,kind,variant,t_total,t_comp,t_comm,t_idle,gflops,residual,growth" > "$OUT"
fi

log() { echo "[$(date +%H:%M:%S)] $*"; }

# ============================================================================
#  E1 + E2 — Escalabilidad fuerte
#  N fijo, p creciente. Es LA figura del informe.
#  Incluye el serial puro, que es el T1 del speedup.
# ============================================================================
exp_strong() {
    log "E2: escalabilidad fuerte"
    for n in $N_STRONG; do
        for r in $(seq 1 "$REPS"); do
            log "  serial n=$n rep=$r"
            "$BIN/lu_serial" -n "$n" -k diagdom -s "$SEED" --csv >> "$OUT"
        done
        for p in $PROCS; do
            for r in $(seq 1 "$REPS"); do
                log "  mpi n=$n p=$p rep=$r"
                $MPIRUN $MPIFLAGS -n "$p" "$BIN/lu_mpi" \
                    -n "$n" -k diagdom -s "$SEED" --csv >> "$OUT"
            done
        done
    done
}

# ============================================================================
#  E3 — Escalabilidad débil
#  El trabajo por proceso se mantiene constante. Como el trabajo de LU crece
#  como n^3, para que n^3/p sea constante hay que tomar n = n0 * p^(1/3).
#  (Este es el error clásico: usar n = n0*p, que multiplica el trabajo por
#  proceso por p^2 y da una curva que parece un desastre sin serlo.)
# ============================================================================
exp_weak() {
    log "E3: escalabilidad debil (n = 1000 * p^(1/3))"
    for p in $PROCS; do
        n=$(python3 -c "print(int(round(1000*($p)**(1/3))))")
        for r in $(seq 1 "$REPS"); do
            log "  mpi n=$n p=$p rep=$r"
            $MPIRUN $MPIFLAGS -n "$p" "$BIN/lu_mpi" \
                -n "$n" -k diagdom -s "$SEED" --csv >> "$OUT"
        done
    done
}

# ============================================================================
#  E4 — Descomposición del tiempo (cómputo / comunicación / ocioso)
#  Con --split-idle metemos una barrera antes de cada colectiva, de modo que
#  la espera por desbalance se contabiliza como T_ocioso y no se disfraza de
#  comunicación. Las barreras agregan overhead, así que estas corridas NO se
#  usan para el speedup: solo para el gráfico de composición.
# ============================================================================
exp_breakdown() {
    log "E4: descomposicion del tiempo"
    for p in $PROCS; do
        for r in $(seq 1 "$REPS"); do
            log "  mpi n=4000 p=$p rep=$r (split-idle)"
            $MPIRUN $MPIFLAGS -n "$p" "$BIN/lu_mpi" \
                -n 4000 -k diagdom -s "$SEED" --split-idle --csv >> "$OUT"
        done
    done
}

# ============================================================================
#  E5 — Reparto cíclico vs bloques contiguos
#  La evidencia experimental de la figura de distribución de datos.
# ============================================================================
exp_dist() {
    log "E5: ciclico vs bloques"
    for p in $PROCS; do
        [ "$p" -eq 1 ] && continue          # con p=1 los dos repartos coinciden
        for d in cyclic block; do
            for r in $(seq 1 "$REPS"); do
                log "  mpi n=2000 p=$p dist=$d rep=$r"
                $MPIRUN $MPIFLAGS -n "$p" "$BIN/lu_mpi" \
                    -n 2000 -d "$d" -k diagdom -s "$SEED" --csv >> "$OUT"
            done
        done
    done
}

# ============================================================================
#  E6 — Precisión y condicionamiento
#  Barremos tipos de matriz. Hilbert se corta en n=1000 porque por encima el
#  número de condición ya desbordó cualquier interpretación razonable.
# ============================================================================
exp_precision() {
    log "E6: precision"
    for k in diagdom random hilbert zerodiag; do
        for n in 100 200 500 1000; do
            "$BIN/lu_serial" -n "$n" -k "$k" -s "$SEED" --csv >> "$OUT" 2>/dev/null
            $MPIRUN $MPIFLAGS -n 4 "$BIN/lu_mpi" -n "$n" -k "$k" -s "$SEED" --csv >> "$OUT" 2>/dev/null
        done
    done
    log "  y la demostracion de por que hace falta pivotear:"
    "$BIN/lu_serial" -n 200 -k zerodiag --nopivot --csv >> "$OUT" 2>&1 || \
        log "  (fallo como se esperaba: pivote nulo sin pivoteo)"
}

# --- despacho ---------------------------------------------------------------
case "$WHICH" in
    strong)     exp_strong ;;
    weak)       exp_weak ;;
    breakdown)  exp_breakdown ;;
    dist)       exp_dist ;;
    precision)  exp_precision ;;
    all)        exp_strong; exp_weak; exp_breakdown; exp_dist; exp_precision ;;
    *)          echo "uso: $0 [all|strong|weak|breakdown|dist|precision]"; exit 1 ;;
esac

log "Listo. Resultados en $OUT"
log "Siguiente paso:  python3 scripts/plots.py"
