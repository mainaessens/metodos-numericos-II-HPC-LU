#!/usr/bin/env python3
# ============================================================================
#  comparar.py — Compara las corridas de todos los sistemas medidos
#
#  Proyecto 6 · Factorización LU Paralela · FACET-UNT 2026
#
#  Lee results/<sistema>/raw.csv para cada sistema y genera:
#     results/comparacion.md                 tablas lado a lado
#     figuras/comparacion/figC1_cluster.png  reproducibilidad en el cluster
#     figuras/comparacion/figC2_pc.png       PC Ryzen (8 núcleos iguales) vs
#                                            PC i7 híbrida (2P + 8E)
#     figuras/comparacion/figC3_descomposicion.png  E4 con --split-idle,
#                                            cluster y PC lado a lado
#
#  Uso:  python3 scripts/comparar.py
# ============================================================================

import os
import sys

sys.argv = sys.argv[:1]              # plots.py lee argv[1]; acá no aplica
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import plots as P                    # reutiliza estilo, paleta y helpers
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib

ROOT = P.ROOT
OUT_FIG = os.path.join(ROOT, "figuras", "comparacion")
os.makedirs(OUT_FIG, exist_ok=True)

SISTEMAS = {
    "pc":           "PC Ezequiel (Ryzen 7 7735HS)",
    "pc-maia":      "PC Maia (i7-1355U)",
    "cluster":      "Cluster — corrida Ezequiel",
    "cluster-maia": "Cluster — corrida Maia",
}


def cargar(s):
    path = os.path.join(ROOT, "results", s, "raw.csv")
    if not os.path.exists(path):
        return None
    df = pd.read_csv(path, on_bad_lines="skip")
    for c in ["n", "p", "t_total", "t_comp", "t_comm", "t_idle", "residual"]:
        df[c] = pd.to_numeric(df[c], errors="coerce")
    df = df.dropna(subset=["n", "p", "t_total"])
    df["n"] = df.n.astype(int); df["p"] = df.p.astype(int)
    df["t_idle"] = df.t_idle.fillna(0)
    return df


def speedups(df):
    """Tabla N, p, T, S, E con T1 = serial puro; excluye corridas split-idle."""
    ser = df[(df.impl == "serial") & (df.variant == "pivot") & (df.kind == "diagdom")]
    par = df[(df.impl == "mpi") & (df.dist == "cyclic") & (df.kind == "diagdom") & (df.t_idle == 0)]
    t1 = ser.groupby("n").t_total.median()
    g = par.groupby(["n", "p"]).t_total.median().reset_index()
    g = g[g.n.map(g.groupby("n").p.nunique()) >= 2]       # solo escalabilidad fuerte
    g["T1"] = g.n.map(t1)
    g["S"] = g.T1 / g.t_total
    g["E"] = g.S / g.p
    return g.dropna()


def descomp(df):
    d = df[(df.impl == "mpi") & (df.t_idle > 0)]
    if d.empty:
        return None
    g = d.groupby("p")[["t_comp", "t_comm", "t_idle"]].median()
    tot = g.sum(axis=1)
    return (100 * g.div(tot, axis=0)).round(1)


data = {s: cargar(s) for s in SISTEMAS}
data = {s: d for s, d in data.items() if d is not None}
sp = {s: speedups(d) for s, d in data.items()}

# ---------------------------------------------------------------- tablas ---
L = ["# Comparación entre sistemas\n",
     "Mediana de 5 repeticiones. T₁ = serial puro. Las corridas con `--split-idle` "
     "no entran en el speedup.\n"]

for grupo, sis in [("Cluster (AMD EPYC 7B12, 48 núcleos asignados)", ["cluster", "cluster-maia"]),
                   ("PC", ["pc", "pc-maia"])]:
    sis = [s for s in sis if s in sp]
    if not sis:
        continue
    L.append(f"\n## {grupo}\n")
    L.append("| N | p | " + " | ".join(f"S — {SISTEMAS[s]}" for s in sis)
             + " | " + " | ".join(f"E — {SISTEMAS[s]}" for s in sis) + " |")
    L.append("|" + "---|" * (2 + 2 * len(sis)))
    keys = sorted(set().union(*[set(zip(sp[s].n, sp[s].p)) for s in sis]))
    for n, p in keys:
        row = [str(n), str(p)]
        for col in ["S", "E"]:
            for s in sis:
                r = sp[s][(sp[s].n == n) & (sp[s].p == p)]
                if r.empty: row.append("–")
                else: row.append(f"{r[col].iloc[0]:.2f}" if col == "S" else f"{100*r[col].iloc[0]:.0f}%")
        L.append("| " + " | ".join(row) + " |")
    # T1 seriales
    L.append("\n| N | " + " | ".join(f"T₁ [s] — {SISTEMAS[s]}" for s in sis) + " |")
    L.append("|" + "---|" * (1 + len(sis)))
    for n in sorted(set().union(*[set(sp[s].n) for s in sis])):
        row = [str(n)]
        for s in sis:
            r = sp[s][sp[s].n == n]
            row.append(f"{r.T1.iloc[0]:.3f}" if not r.empty else "–")
        L.append("| " + " | ".join(row) + " |")

L.append("\n## Descomposición del tiempo (E4, `--split-idle`, N = 4000)\n")
for s, d in data.items():
    t = descomp(d)
    if t is None:
        continue
    L.append(f"\n**{SISTEMAS[s]}**\n")
    L.append("| p | Cómputo | Comunicación | Ocioso |\n|---|---|---|---|")
    for p, r in t.iterrows():
        L.append(f"| {p} | {r.t_comp}% | {r.t_comm}% | {r.t_idle}% |")

L.append("\n## Precisión (residuo relativo, serial, mediana)\n")
L.append("| Matriz | N | " + " | ".join(SISTEMAS[s] for s in data) + " |")
L.append("|" + "---|" * (2 + len(data)))
res = {s: d[(d.impl == "serial") & (d.variant == "pivot")].groupby(["kind", "n"]).residual.median()
       for s, d in data.items()}
idx = sorted(set().union(*[set(r.index) for r in res.values()]))
for k, n in idx:
    if n > 1000 and k == "diagdom":
        continue
    L.append(f"| {k} | {n} | " + " | ".join(
        f"{res[s].get((k, n)):.3e}" if (k, n) in res[s] else "–" for s in data) + " |")

with open(os.path.join(ROOT, "results", "comparacion.md"), "w", encoding="utf-8") as f:
    f.write("\n".join(L) + "\n")
print("  -> results/comparacion.md")


# --------------------------------------------------------------- figuras ---
def fig_pares(sis, ns, fname, titulo, sub):
    sis = [s for s in sis if s in sp]
    if len(sis) < 2:
        return
    fig, ax = plt.subplots(figsize=(7.2, 4.6))
    ps = sorted(set().union(*[set(sp[s].p) for s in sis]))
    ax.plot(ps, ps, "--", color=P.IDEAL, lw=1.4, zorder=1)
    ax.annotate("speedup lineal (ideal)", xy=(ps[len(ps)//2], ps[len(ps)//2]),
                xytext=(8, -14), textcoords="offset points", fontsize=8.5, color=P.IDEAL)
    estilos = ["-o", "--s"]
    for i, n in enumerate(ns):
        for j, s in enumerate(sis):
            d = sp[s][sp[s].n == n].sort_values("p")
            if d.empty: continue
            ax.plot(d.p, d.S, estilos[j], color=P.C[i % len(P.C)], lw=1.8, ms=6,
                    mfc="white" if j else P.C[i % len(P.C)], zorder=3,
                    label=f"N = {n} · {SISTEMAS[s]}")
    # Ejes log-log en base 2: así la recta ideal S = p se ve RECTA.
    ax.set_xscale("log", base=2); ax.set_xticks(ps)
    ax.set_yscale("log", base=2)
    ymax = ps[-1] * 1.3
    ax.set_yticks([t for t in [1, 2, 4, 8, 16, 32] if t <= ymax]); ax.set_ylim(0.8, ymax)
    for a in (ax.xaxis, ax.yaxis):
        a.set_major_formatter(matplotlib.ticker.ScalarFormatter())
        a.set_minor_formatter(matplotlib.ticker.NullFormatter())
    ax.margins(x=0.08)
    P.style(ax, "Cantidad de procesos  p", "Speedup  S(p) = T₁ / T(p)", titulo, sub)
    ax.legend(loc="upper left", fontsize=8)
    fig.savefig(os.path.join(OUT_FIG, fname)); plt.close(fig)
    print(f"  -> figuras/comparacion/{fname}")


fig_pares(["cluster", "cluster-maia"], [2000, 4000, 8000], "figC1_cluster.png",
          "Reproducibilidad en el cluster: dos corridas independientes",
          "Línea llena: corrida de Ezequiel · línea punteada: corrida de Maia")
fig_pares(["pc", "pc-maia"], [1000, 4000], "figC2_pc.png",
          "Dos notebooks: núcleos homogéneos vs híbridos",
          "Ryzen 7 7735HS (8 núcleos iguales) vs i7-1355U (2 rápidos + 8 eficientes)")

# E4 lado a lado
desc = [(s, descomp(d)) for s, d in data.items()]
desc = [(s, t) for s, t in desc if t is not None]
if desc:
    fig, axs = plt.subplots(1, len(desc), figsize=(4.2 * len(desc), 4.4), sharey=True,
                            squeeze=False)
    for ax, (s, t) in zip(axs[0], desc):
        x = range(len(t))
        ax.bar(x, t.t_comp, color=P.C[0], edgecolor="white", linewidth=2, label="Cómputo local")
        ax.bar(x, t.t_comm, bottom=t.t_comp, color=P.C[1], edgecolor="white", linewidth=2,
               label="Comunicación")
        ax.bar(x, t.t_idle, bottom=t.t_comp + t.t_comm, color=P.C[2], edgecolor="white",
               linewidth=2, label="Ocioso (desbalance)")
        ax.set_xticks(list(x)); ax.set_xticklabels([str(p) for p in t.index])
        ax.set_ylim(0, 100)
        P.style(ax, "Procesos  p", "% del tiempo total" if ax is axs[0][0] else "",
                SISTEMAS[s])
    h, l = axs[0][0].get_legend_handles_labels()
    fig.legend(h, l, loc="lower center", ncol=3, fontsize=9, bbox_to_anchor=(0.5, -0.04))
    fig.suptitle("¿En qué se va el tiempo?  N = 4000, con --split-idle", x=0.02, ha="left",
                 fontweight="bold", color=P.INK)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.savefig(os.path.join(OUT_FIG, "figC3_descomposicion.png")); plt.close(fig)
    print("  -> figuras/comparacion/figC3_descomposicion.png")
