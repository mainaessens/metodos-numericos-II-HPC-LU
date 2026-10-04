#!/usr/bin/env python3
# ============================================================================
#  plots.py — Análisis estadístico y figuras del informe
#
#  Proyecto 6 · Factorización LU Paralela · FACET-UNT 2026
#
#  Lee results/raw.csv (todas las repeticiones crudas), calcula la MEDIANA
#  por punto experimental y genera las figuras en figuras/.
#
#  Requisitos:  pip install pandas matplotlib
#  Uso:         python3 scripts/plots.py              (results/raw.csv -> figuras/)
#               python3 scripts/plots.py pc-maia      (results/pc-maia/raw.csv
#                                                      -> figuras/pc-maia/)
#
#  ---------------------------------------------------------------------------
#  POR QUÉ LA MEDIANA Y NO EL PROMEDIO
#
#  El cluster no es dedicado: hay otros jobs, ruido del sistema operativo,
#  contención de red. Eso produce corridas ocasionalmente MUY lentas, pero
#  nunca corridas anormalmente rápidas: la distribución de tiempos tiene una
#  cola larga hacia arriba. El promedio se deja arrastrar por esa cola; la
#  mediana no. Reportamos además la dispersión (mín, máx) para que se vea
#  cuánto ruido había.
# ============================================================================

import os
import sys

try:
    import pandas as pd
    import matplotlib
    matplotlib.use("Agg")          # backend sin ventana: funciona por SSH
    import matplotlib.pyplot as plt
except ImportError:
    sys.exit("Faltan dependencias.  pip install pandas matplotlib")

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
# Opcional: nombre del sistema (pc, cluster, pc-maia, cluster-maia, ...).
# Con nombre, lee results/<sistema>/raw.csv y escribe en figuras/<sistema>/.
SISTEMA = sys.argv[1].strip("/") if len(sys.argv) > 1 else ""
RES  = os.path.join(ROOT, "results", SISTEMA)
RAW  = os.path.join(RES, "raw.csv")
FIGS = os.path.join(ROOT, "figuras", SISTEMA)
os.makedirs(FIGS, exist_ok=True)

# ---------------------------------------------------------------------------
#  PALETA
#
#  Orden fijo: el color identifica la SERIE, nunca su posición en el ranking.
#  Si se filtra una serie, las demás conservan su color.
#  Validada para daltonismo (protanopía/deuteranopía/tritanopía) y para
#  contraste sobre fondo claro.
# ---------------------------------------------------------------------------
C = ["#1F6FB2",   # azul     - slot 1
     "#D1651A",   # naranja  - slot 2
     "#8B4FA8",   # violeta  - slot 3
     "#1B7F5A"]   # verde    - slot 4

INK   = "#1A1D24"      # texto principal
INK2  = "#5B6472"      # texto secundario
GRID  = "#DDE1E7"      # grilla: recesiva, nunca compite con los datos
IDEAL = "#9AA3B2"      # líneas de referencia (ideal, teórico)

plt.rcParams.update({
    "figure.dpi":        130,
    "savefig.dpi":       200,
    "savefig.bbox":      "tight",
    "font.size":         10,
    "axes.titlesize":    12,
    "axes.titleweight":  "bold",
    "axes.labelsize":    10,
    "axes.edgecolor":    GRID,
    "axes.labelcolor":   INK2,
    "axes.titlecolor":   INK,
    "text.color":        INK,
    "xtick.color":       INK2,
    "ytick.color":       INK2,
    "axes.spines.top":   False,
    "axes.spines.right": False,
    "grid.color":        GRID,
    "grid.linewidth":    0.8,
    "legend.frameon":    False,
    "legend.fontsize":   9,
})


def load():
    if not os.path.exists(RAW):
        sys.exit(f"No existe {RAW}.\nCorré primero:  ./scripts/run_experiments.sh")
    # on_bad_lines="skip": descarta lineas basura (p.ej. mensajes de error que
    # se hayan colado en el CSV) en vez de romper la lectura.
    df = pd.read_csv(RAW, on_bad_lines="skip")
    for c in ["n", "p", "t_total", "t_comp", "t_comm", "t_idle", "gflops", "residual"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(subset=["n", "p", "t_total"])
    df["n"] = df["n"].astype(int)
    df["p"] = df["p"].astype(int)
    df["t_idle"] = df["t_idle"].fillna(0)
    return df


def med(df, by, val="t_total"):
    """Mediana por grupo + dispersión, para reportar honestamente el ruido."""
    g = df.groupby(by)[val]
    out = g.median().reset_index().rename(columns={val: "med"})
    out["lo"]  = g.min().values
    out["hi"]  = g.max().values
    out["reps"] = g.count().values
    return out


def style(ax, xlabel, ylabel, title, subtitle=None):
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title, loc="left", pad=16 if subtitle else 10)
    if subtitle:
        ax.text(0, 1.02, subtitle, transform=ax.transAxes,
                fontsize=9, color=INK2, va="bottom")
    ax.grid(True, axis="y", alpha=0.9)
    ax.set_axisbelow(True)


def save(fig, name):
    path = os.path.join(FIGS, name)
    fig.savefig(path)
    plt.close(fig)
    print(f"  -> {os.path.relpath(path, ROOT)}")


# ===========================================================================
#  FIG 1 y 2 — Escalabilidad fuerte: speedup y eficiencia
#
#  Forma elegida: líneas, porque el eje x (cantidad de procesos) es una
#  progresión ordenada y lo que importa es la TENDENCIA, no comparar
#  categorías sueltas.
#
#  T1 sale de las corridas SERIALES, no de mpi con p=1. Es la diferencia
#  entre medir el speedup y medirse a uno mismo.
# ===========================================================================
def fig_strong(df):
    ser = df[(df.impl == "serial") & (df.variant == "pivot") & (df.kind == "diagdom")]
    # t_idle == 0 excluye las corridas con --split-idle (E4): sus barreras
    # agregan overhead y NO deben usarse para el speedup.
    par = df[(df.impl == "mpi") & (df.dist == "cyclic") & (df.kind == "diagdom")
             & (df.t_idle == 0)]
    if ser.empty or par.empty:
        print("  (sin datos de escalabilidad fuerte)"); return

    t1 = med(ser, ["n"]).set_index("n")["med"].to_dict()
    mp = med(par, ["n", "p"])

    # Quedarnos SOLO con los N que son puntos de escalabilidad fuerte: los que
    # se corrieron con al menos 3 valores distintos de p. Así se descartan las
    # corridas de los experimentos de precisión (un solo p) y las de
    # escalabilidad débil (un N distinto por cada p), que si no se colarían
    # como series falsas en este gráfico.
    cuenta_p = mp.groupby("n")["p"].nunique()
    ns = [int(n) for n in sorted(set(mp.n) & set(t1.keys())) if cuenta_p.get(n, 0) >= 3]

    # La paleta tiene 4 slots y NO se recicla: si hubiera más de 4 N, nos
    # quedamos con los 4 más grandes, que son los que importan para el análisis.
    if len(ns) > len(C):
        ns = ns[-len(C):]

    if not ns:
        print("  (no hay N con serial y paralelo: no se puede calcular speedup)"); return

    for fname, ylab, title, sub, transform in [
        ("fig1_speedup.png", "Speedup  S(p) = T₁ / T(p)",
         "Escalabilidad fuerte: speedup",
         "Matriz diagonal dominante · reparto cíclico · mediana de las repeticiones",
         lambda s, p: s),
        ("fig2_eficiencia.png", "Eficiencia  E(p) = S(p) / p",
         "Escalabilidad fuerte: eficiencia",
         "El mismo dato que el speedup, pero normalizado: es el número honesto",
         lambda s, p: s / p),
    ]:
        fig, ax = plt.subplots(figsize=(7.2, 4.6))
        pmax = int(mp.p.max())

        if fname.startswith("fig1"):
            ps = sorted(mp.p.unique())
            ax.plot(ps, ps, "--", color=IDEAL, lw=1.4, zorder=1)
            mid = ps[len(ps) // 2]
            ax.annotate("speedup lineal (ideal)", xy=(mid, mid),
                        xytext=(8, -14), textcoords="offset points",
                        ha="left", fontsize=8.5, color=IDEAL)
        else:
            ax.axhline(1.0, ls="--", color=IDEAL, lw=1.4, zorder=1)
            ax.annotate("100 % de eficiencia", xy=(pmax, 1.0),
                        xytext=(-6, 6), textcoords="offset points",
                        ha="right", fontsize=8.5, color=IDEAL)

        for i, n in enumerate(ns):
            sub_df = mp[mp.n == n].sort_values("p")
            n = int(n)
            sp = t1[n] / sub_df["med"].values
            y  = transform(sp, sub_df["p"].values)
            ax.plot(sub_df.p, y, "-o", color=C[i % len(C)], lw=1.8, ms=6,
                    label=f"N = {n}", zorder=3)
            # Etiqueta directa en el extremo: identidad sin depender del color.
            ax.annotate(f"N={n}", xy=(sub_df.p.values[-1], y[-1]),
                        xytext=(7, 0), textcoords="offset points",
                        fontsize=9, color=C[i % len(C)], va="center")

        ax.set_xscale("log", base=2)
        ax.set_xticks(sorted(mp.p.unique()))
        ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
        ax.margins(x=0.12)
        style(ax, "Cantidad de procesos  p", ylab, title, sub)
        ax.legend(loc="upper left")
        save(fig, fname)


# ===========================================================================
#  FIG 3 — Descomposición del tiempo
#
#  Forma elegida: barras apiladas en FRACCIÓN del total, porque la pregunta
#  no es "cuánto tardó" sino "en qué se fue el tiempo". Las tres partes son
#  componentes de un mismo todo: apilarlas es la forma correcta.
# ===========================================================================
def fig_breakdown(df):
    # Se usan SOLO las corridas con --split-idle (t_idle > 0). Si no hay,
    # se cae a las corridas comunes, donde la espera queda dentro de "comunicacion".
    d = df[(df.impl == "mpi") & (df.dist == "cyclic") & (df.t_idle > 0)]
    if d.empty:
        print("  (aviso: no hay corridas con --split-idle; el ocioso queda en 0)")
        d = df[(df.impl == "mpi") & (df.dist == "cyclic")]
    if d.empty:
        print("  (sin datos de descomposicion)"); return

    n = int(d.n.mode().iloc[0])
    d = d[d.n == n]
    g = d.groupby("p")[["t_comp", "t_comm", "t_idle"]].median().reset_index().sort_values("p")
    tot = (g.t_comp + g.t_comm + g.t_idle).replace(0, 1)

    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    x = range(len(g))
    comp = 100 * g.t_comp / tot
    comm = 100 * g.t_comm / tot
    idle = 100 * g.t_idle / tot

    # El borde blanco de 2px es el "spacer" que separa segmentos contiguos:
    # sin él, dos colores vecinos se leen como uno solo.
    ax.bar(x, comp, color=C[0], label="Cómputo local",  edgecolor="white", linewidth=2)
    ax.bar(x, comm, bottom=comp, color=C[1], label="Comunicación", edgecolor="white", linewidth=2)
    ax.bar(x, idle, bottom=comp + comm, color=C[2], label="Ocioso (desbalance)",
           edgecolor="white", linewidth=2)

    for xi, c in zip(x, comp):
        if c > 12:
            ax.text(xi, c / 2, f"{c:.0f}%", ha="center", va="center",
                    color="white", fontsize=9, fontweight="bold")

    ax.set_xticks(list(x))
    ax.set_xticklabels([str(int(v)) for v in g.p])
    ax.set_ylim(0, 100)
    style(ax, "Cantidad de procesos  p", "Porcentaje del tiempo total",
          f"¿En qué se va el tiempo?  (N = {n})",
          "La fracción de comunicación crece con p: es la razón por la que el speedup se aplana")
    ax.legend(loc="lower right", ncol=1)
    save(fig, "fig3_descomposicion.png")


# ===========================================================================
#  FIG 4 — Reparto cíclico vs bloques contiguos
#
#  Forma elegida: barras agrupadas. Son dos categorías discretas comparadas
#  para cada p: la comparación directa de alturas es lo que se quiere leer.
# ===========================================================================
def fig_dist(df):
    d = df[(df.impl == "mpi") & (df.p > 1)]
    if d.empty or d.dist.nunique() < 2:
        print("  (sin comparacion de repartos)"); return

    n = int(d[d.dist == "block"].n.mode().iloc[0])
    d = d[d.n == n]
    g = med(d, ["dist", "p"]).pivot(index="p", columns="dist", values="med").sort_index()
    if "block" not in g or "cyclic" not in g:
        print("  (faltan corridas de uno de los repartos)"); return

    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    x = range(len(g)); w = 0.38
    ax.bar([i - w/2 for i in x], g["cyclic"], w, color=C[0],
           label="Cíclico por filas", edgecolor="white", linewidth=2)
    ax.bar([i + w/2 for i in x], g["block"], w, color=C[1],
           label="Bloques contiguos", edgecolor="white", linewidth=2)

    for i, (cy, bl) in enumerate(zip(g["cyclic"], g["block"])):
        if cy > 0:
            ax.text(i + w/2, bl, f"+{100*(bl-cy)/cy:.0f}%", ha="center",
                    va="bottom", fontsize=9, color=C[1], fontweight="bold")

    ax.set_xticks(list(x))
    ax.set_xticklabels([str(int(v)) for v in g.index])
    style(ax, "Cantidad de procesos  p", "Tiempo de factorización [s]",
          f"El reparto de datos decide la escalabilidad  (N = {n})",
          "Mismo algoritmo, misma máquina: lo único que cambia es quién tiene qué fila")
    ax.legend(loc="upper right")
    save(fig, "fig4_reparto.png")


# ===========================================================================
#  FIG 5 — Escalabilidad débil
#  Una sola serie: no lleva caja de leyenda, el título ya la nombra.
# ===========================================================================
def fig_weak(df):
    d = df[(df.impl == "mpi") & (df.dist == "cyclic") & (df.kind == "diagdom")].copy()
    if d.empty: return
    # Los puntos de escalabilidad débil son los que cumplen n ≈ 1000·p^(1/3).
    d["esperado"] = (1000 * d.p ** (1/3)).round()
    d = d[(d.n - d.esperado).abs() <= 2]
    if d.p.nunique() < 2:          # hace falta al menos 2 valores de p
        print("  (sin datos de escalabilidad debil)"); return

    g = med(d, ["p"]).sort_values("p")
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    ax.plot(g.p, g["med"], "-o", color=C[0], lw=1.8, ms=6, zorder=3)
    ax.fill_between(g.p, g.lo, g.hi, color=C[0], alpha=0.14, lw=0, zorder=2)
    ax.axhline(g["med"].iloc[0], ls="--", color=IDEAL, lw=1.4, zorder=1)
    ax.annotate("escalabilidad débil perfecta", xy=(g.p.max(), g["med"].iloc[0]),
                xytext=(-6, 6), textcoords="offset points", ha="right",
                fontsize=8.5, color=IDEAL)

    ax.set_xscale("log", base=2)
    ax.set_xticks(list(g.p))
    ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
    ax.margins(x=0.12)
    style(ax, "Cantidad de procesos  p", "Tiempo de factorización [s]",
          "Escalabilidad débil: trabajo constante por proceso",
          "N crece como p^(1/3) porque el trabajo de LU crece como N³ · banda = mín–máx observado")
    save(fig, "fig5_debil.png")


# ===========================================================================
#  FIG 6 — Precisión: residuo vs tamaño, por tipo de matriz
#  Escala logarítmica en y porque los residuos abarcan varios órdenes.
# ===========================================================================
def fig_precision(df):
    d = df[(df.impl == "serial") & (df.variant == "pivot")]
    if d.empty: return
    g = d.groupby(["kind", "n"])["residual"].median().reset_index()
    kinds = [k for k in ["diagdom", "random", "zerodiag", "hilbert"] if k in set(g.kind)]
    if not kinds: return

    nombres = {"diagdom": "Diagonal dominante", "random": "Aleatoria",
               "hilbert": "Hilbert (mal condicionada)", "zerodiag": "Diagonal nula"}

    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    for i, k in enumerate(kinds):
        s = g[g.kind == k].sort_values("n")
        ax.plot(s.n, s.residual, "-o", color=C[i % len(C)], lw=1.8, ms=6,
                label=nombres.get(k, k), zorder=3)

    eps = 2.22e-16
    ax.axhline(eps, ls="--", color=IDEAL, lw=1.4, zorder=1)
    ax.annotate("ε de máquina ≈ 2.2e-16", xy=(g.n.max(), eps),
                xytext=(-6, 6), textcoords="offset points", ha="right",
                fontsize=8.5, color=IDEAL)

    ax.set_yscale("log")
    style(ax, "Dimensión de la matriz  N", "Residuo relativo  ‖Ax−b‖∞ / (‖A‖∞‖x‖∞)",
          "Precisión: el residuo se mantiene en el orden de ε",
          "Con pivoteo parcial. Un residuo chico NO garantiza solución exacta: ver κ(A)")
    ax.legend(loc="best")
    save(fig, "fig6_residuo.png")


# ===========================================================================
#  Tabla resumen en texto, para pegar en el informe.
# ===========================================================================
def tabla(df):
    ser = df[(df.impl == "serial") & (df.variant == "pivot") & (df.kind == "diagdom")]
    # t_idle == 0 excluye las corridas con --split-idle (E4): sus barreras
    # agregan overhead y NO deben usarse para el speedup.
    par = df[(df.impl == "mpi") & (df.dist == "cyclic") & (df.kind == "diagdom")
             & (df.t_idle == 0)]
    if ser.empty or par.empty: return

    t1 = med(ser, ["n"]).set_index("n")["med"].to_dict()
    mp = med(par, ["n", "p"])

    path = os.path.join(RES, "resumen.md")
    with open(path, "w", encoding="utf-8") as f:
        f.write("# Resumen de resultados\n\n")
        f.write("Mediana de las repeticiones. T1 = version serial pura.\n\n")
        f.write("| N | p | T(p) [s] | min-max [s] | reps | Speedup | Eficiencia | GFLOP/s |\n")
        f.write("|---|---|---|---|---|---|---|---|\n")
        # Solo los N de escalabilidad fuerte (2 o mas valores de p): asi no se
        # cuelan las corridas sueltas del experimento de precision.
        cuenta_p = mp.groupby("n")["p"].nunique()
        for n in sorted(set(mp.n) & set(t1)):
            if cuenta_p.get(n, 0) < 2:
                continue
            f.write(f"| {n} | 1 (serial) | {t1[n]:.4f} | - | - | 1.00 | 100% | - |\n")
            for _, r in mp[mp.n == n].sort_values("p").iterrows():
                sp = t1[n] / r["med"]
                f.write(f"| {n} | {int(r.p)} | {r['med']:.4f} | "
                        f"{r.lo:.4f}–{r.hi:.4f} | {int(r.reps)} | "
                        f"{sp:.2f} | {100*sp/r.p:.0f}% | "
                        f"{(2/3)*n**3/r['med']/1e9:.2f} |\n")
        f.write("\n")
    print(f"  -> {os.path.relpath(path, ROOT)}")


def main():
    print(f"Leyendo {os.path.relpath(RAW, ROOT)} ...")
    df = load()
    print(f"  {len(df)} corridas\n")
    print("Generando figuras:")
    fig_strong(df)
    fig_breakdown(df)
    fig_dist(df)
    fig_weak(df)
    fig_precision(df)
    tabla(df)
    print("\nListo.")


if __name__ == "__main__":
    main()
