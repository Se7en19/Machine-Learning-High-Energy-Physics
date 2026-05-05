"""
=============================================================
ANÁLISIS CORREGIDO — Bar Strip Detector
=============================================================
Correcciones respecto a scripts originales:
  1) Bug crítico: ang_list vs mu_angles en analisis_separado.py
  2) Ahora se usa ConeAngle (col 14) en lugar de reconstrucción
  3) Barras de error binomial en TODAS las gráficas
  4) Cálculo de eficiencia robusto

3 gráficas:
  1) landau_corregida.png    — Distribución de Landau (run 45, ~1 GeV/c)
  2) eff_momento_corregida.png — Eficiencia vs momento (GeV/c)
  3) eff_angulo_corregida.png  — Eficiencia vs ángulo del cono
=============================================================
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import uproot, glob, os, re

MIXED_DIR = "../Classifier/data/mixed"
plt.rcParams.update({
    "font.size": 12, "axes.titlesize": 14, "axes.labelsize": 12,
    "xtick.labelsize": 10, "ytick.labelsize": 10, "legend.fontsize": 11,
    "figure.dpi": 150, "figure.facecolor": "white",
})

def best_cycle(f, tree_name="Hits"):
    keys = [k for k in f.keys(cycle=True) if k.split(";")[0] == tree_name]
    if not keys: return None
    best = max(keys, key=lambda k: f[k].num_entries)
    return f[best] if f[best].num_entries > 0 else None

def run_num(f):
    m = re.search(r"output_run(\d+)", f)
    return int(m.group(1)) if m else 0

files = sorted(glob.glob(os.path.join(MIXED_DIR, "output_run*.root")), key=run_num)
N_RUNS = len(files)
info_str = r"Geant4  |  Bar Strip Detector  |  70 cm Fe + BC404  |  FTFP\_BERT"

# ===========================================================================
# 1) LANDAU — run 45 (p0 ≈ 1.0225 GeV/c)
# ===========================================================================
fpath = os.path.join(MIXED_DIR, "output_run45.root")
with uproot.open(fpath) as f:
    tree = best_cycle(f)
    d45 = tree.arrays(["fdEdx", "particleID", "fEvent", "layerID", "Momentum"], library="np")

mu_mask = d45["particleID"] == 0
pi_mask = d45["particleID"] == 1

mu_dedx   = d45["fdEdx"][mu_mask]
pi_dedx   = d45["fdEdx"][pi_mask]
mu_events = len(np.unique(d45["fEvent"][mu_mask]))
pi_events = len(np.unique(d45["fEvent"][pi_mask]))
tot_events = len(np.unique(d45["fEvent"]))

mu_hits_tot = len(mu_dedx)
pi_hits_tot = len(pi_dedx)

mu_mean, mu_med = np.mean(mu_dedx), np.median(mu_dedx)
pi_mean, pi_med = np.mean(pi_dedx), np.median(pi_dedx)

# Fracciones sobre umbrales
mu_above = np.sum(mu_dedx > 0.5)
pi_above = np.sum(pi_dedx > 0.5)

fig1, ax1 = plt.subplots(figsize=(11, 7))
bins = np.linspace(0.01, 5.0, 100)

ax1.hist(mu_dedx, bins=bins, histtype="stepfilled",
         color="steelblue", alpha=0.5, density=True,
         label=rf"$\mu^+$  ({mu_events} eventos, {mu_hits_tot} hits)")
ax1.hist(mu_dedx, bins=bins, histtype="step",
         color="steelblue", lw=2, density=True)
ax1.hist(pi_dedx, bins=bins, histtype="stepfilled",
         color="tomato", alpha=0.4, density=True,
         label=rf"$\pi^+$  ({pi_events} eventos, {pi_hits_tot} hits)")
ax1.hist(pi_dedx, bins=bins, histtype="step",
         color="tomato", lw=2, density=True)

ax1.axvline(0.5, color="green", lw=2.5, ls="--",
            label=r"Umbral dE/dx = 0.5 MeV/mm")
ax1.set_xlabel(r"$dE/dx$ por paso  (MeV/mm)")
ax1.set_ylabel("Densidad (normalizada)")
ax1.set_yscale("log")
ax1.set_xlim(0.01, 5.0)
ax1.set_ylim(1e-3, 5)
ax1.legend(fontsize=10, loc="upper right")
ax1.set_title(r"Distribución de Landau — haz mixto $\mu^+$/$\pi^+$  ($p_0 \approx 1$ GeV/c)", fontsize=13)
ax1.grid(True, alpha=0.3, which="both")

stats_text = (
    f"Estadísticas de dE/dx (MeV/mm):\n"
    f"{'':>8} {'μ⁺':>10} {'π⁺':>10}\n"
    f"{'Media':>8} {mu_mean:>8.4f}  {pi_mean:>8.4f}\n"
    f"{'Mediana':>8} {mu_med:>8.4f}  {pi_med:>8.4f}\n"
    f"{'Hits totales':>8} {mu_hits_tot:>8d}  {pi_hits_tot:>8d}"
)
ax1.text(0.55, 0.65, stats_text, transform=ax1.transAxes, fontsize=8.5,
         family="monospace", va="top", ha="left",
         bbox=dict(boxstyle="round,pad=0.5", facecolor="white", alpha=0.85, edgecolor="gray"))

info1 = (f"Geant4  |  ~1000 μ⁺ + ~1000 π⁺ generados  |  "
         f"{mu_events} μ⁺ detectados + {pi_events} π⁺ detectados = {tot_events} total")
ax1.text(0.01, 1.012, info1, transform=ax1.transAxes,
         fontsize=8.5, va="bottom", ha="left", color="#333333", style="italic")

fig1.tight_layout()
fig1.savefig("img/landau_corregida.png", dpi=150, bbox_inches="tight")
print("Saved: img/landau_corregida.png")
plt.close(fig1)

# ===========================================================================
# 2) EFICIENCIA VS MOMENTO (GeV/c)
# ===========================================================================
p0_vals_GeV = np.logspace(np.log10(0.05), np.log10(10.0), 80)
records_mu, records_pi = [], []

for fpath in files:
    rn = run_num(fpath)
    if rn >= 80: continue
    p0 = p0_vals_GeV[rn]
    with uproot.open(fpath) as f:
        tree = best_cycle(f)
        if tree is None: continue
        d = tree.arrays(["particleID", "fEvent"], library="np")

    for pid, rec_list in [(0, records_mu), (1, records_pi)]:
        mask = d["particleID"] == pid
        n = mask.sum()
        if n == 0:
            rec_list.append({"p0": p0, "eff": 0, "n_det": 0, "run": rn})
        else:
            n_det = len(np.unique(d["fEvent"][mask]))
            eff = n_det / 1000.0
            rec_list.append({"p0": p0, "eff": eff, "n_det": n_det, "run": rn})

mu_recs = sorted(records_mu, key=lambda r: r["p0"])
pi_recs = sorted(records_pi, key=lambda r: r["p0"])

p_mu = np.array([r["p0"] for r in mu_recs])
e_mu = np.array([r["eff"] for r in mu_recs])
n_mu = np.array([r["n_det"] for r in mu_recs])
p_pi = np.array([r["p0"] for r in pi_recs])
e_pi = np.array([r["eff"] for r in pi_recs])
n_pi = np.array([r["n_det"] for r in pi_recs])

# Barras de error binomial
err_mu = np.sqrt(e_mu * (1 - e_mu) / 1000)
err_pi = np.sqrt(e_pi * (1 - e_pi) / 1000)

fig2, ax2 = plt.subplots(figsize=(11, 7))

ax2.errorbar(p_mu, e_mu*100, yerr=err_mu*100, fmt="o-", color="steelblue",
             lw=2, ms=5, capsize=3, capthick=1, label=r"$\mu^+$", zorder=3)
ax2.errorbar(p_pi, e_pi*100, yerr=err_pi*100, fmt="s-", color="tomato",
             lw=2, ms=5, capsize=3, capthick=1, label=r"$\pi^+$", zorder=3)

ax2.axhline(50, color="gray", ls=":", lw=1, alpha=0.5)
ax2.set_xlabel(r"Momento inicial $p_0$  (GeV/c)")
ax2.set_ylabel(r"Eficiencia $\varepsilon$ (%)")
ax2.set_ylim(-3, 105)
ax2.set_xlim(0.04, 11)
ax2.legend(fontsize=12, framealpha=0.85, loc="upper left")
ax2.set_title(r"Eficiencia de detección $\varepsilon = N_{\mathrm{det}}/N_{\mathrm{gen}}$ vs $p_0$", fontsize=13)
ax2.grid(True, alpha=0.3, which="both")

# Regímenes
ax2.axvspan(0.04, 0.5, alpha=0.06, color="gray")
ax2.text(0.15, 52, "Régimen I:\nμ⁺ no penetran\n70 cm Fe\n(R μ⁺ < 70 cm)", fontsize=9,
         color="gray", ha="center", va="bottom", alpha=0.7)
ax2.axvspan(0.5, 1.1, alpha=0.06, color="steelblue")
ax2.text(0.74, 52, "Régimen II:\ntransición μ⁺\n(R μ⁺ ≈ 70 cm)\nsubida sigmoidea",
         fontsize=9, color="steelblue", ha="center", va="bottom", alpha=0.7)
ax2.axvspan(1.1, 10, alpha=0.06, color="green")
ax2.text(4.0, 52, "Régimen III:\nmeseta μ⁺\n(~89±1%)\ntodos penetran",
         fontsize=9, color="green", ha="center", va="bottom", alpha=0.7)

ax2.errorbar([], [], fmt="none",
             label=r"Barras: $\sigma_\varepsilon = \sqrt{\varepsilon(1-\varepsilon)/1000}$")
ax2.legend(fontsize=10, framealpha=0.85, loc="upper left")

ax2.text(0.01, 1.012, info_str, transform=ax2.transAxes,
         fontsize=8.5, va="bottom", ha="left", color="#333333", style="italic")

fig2.savefig("img/eff_momento_corregida.png", dpi=150, bbox_inches="tight")
print("Saved: img/eff_momento_corregida.png")
plt.close(fig2)

# ===========================================================================
# 3) EFICIENCIA VS ÁNGULO DEL CONO
#    CORREGIDO: Usa ConeAngle (col 14) y ang_list en lugar de mu_angles
# ===========================================================================
L_src = 2000.0
mu_angles_correct, pi_angles_correct = [], []

for fpath in files:
    rn = run_num(fpath)
    if rn >= 80: continue
    with uproot.open(fpath) as f:
        tree = best_cycle(f)
        if tree is None: continue
        # Leer ConeAngle (col 14) + particleID + fEvent
        d = tree.arrays(["fEvent", "particleID", "ConeAngle", "layerID"], library="np")

    for pid_val, ang_list in [(0, mu_angles_correct), (1, pi_angles_correct)]:
        mask_pid = d["particleID"] == pid_val
        if mask_pid.sum() == 0: continue

        ev_p = d["fEvent"][mask_pid]
        cone_p = d["ConeAngle"][mask_pid]
        lay_p = d["layerID"][mask_pid]

        for ev in np.unique(ev_p):
            m_ev = ev_p == ev
            lay_ev = lay_p[m_ev]
            cone_ev = cone_p[m_ev]

            m0 = lay_ev == 0
            m1 = lay_ev == 1
            if m0.sum() == 0 or m1.sum() == 0:
                continue

            cone_med = np.median(cone_ev)
            ang_list.append(np.degrees(cone_med))

theta_max = np.degrees(np.arctan2(np.sqrt(350**2 + 350**2), L_src))
n_bins = 18
bin_edges = np.linspace(0, theta_max, n_bins + 1)
bc = 0.5 * (bin_edges[:-1] + bin_edges[1:])

# N_gen teórico
rng = np.random.default_rng(42)
N_mc = 2_000_000
tx_mc = rng.uniform(-350, 350, N_mc)
ty_mc = rng.uniform(-350, 350, N_mc)
theta_mc = np.degrees(np.arctan2(np.sqrt(tx_mc**2 + ty_mc**2), L_src))
h_mc, _ = np.histogram(theta_mc, bins=bin_edges)
p_theory = h_mc / h_mc.sum()
n_total_per_species = N_RUNS * 1000.0
n_gen = p_theory * n_total_per_species

fig3, ax3 = plt.subplots(figsize=(11, 7))

for angles, color, label, marker in [
    (mu_angles_correct, "steelblue", r"$\mu^+$", "o"),
    (pi_angles_correct, "tomato", r"$\pi^+$", "s"),
]:
    if not angles: continue
    h_det, _ = np.histogram(angles, bins=bin_edges)
    eff = np.where(n_gen > 10, h_det / n_gen, np.nan)
    err = np.where(n_gen > 10, np.sqrt(eff * (1 - eff) / n_gen), np.nan)
    eff = np.clip(eff, 0, 1.05)
    ax3.errorbar(bc, eff*100, yerr=err*100, fmt=marker+"-", color=color,
                 lw=2, ms=6, capsize=3, capthick=1, label=label)

ax3.set_xlabel(r"Ángulo del cono $\theta$ (°)")
ax3.set_ylabel(r"Eficiencia $\varepsilon$ (%)")
ax3.set_xlim(0, theta_max + 0.3)
ax3.set_ylim(-3, 110)
ax3.axhline(50, color="gray", ls=":", lw=1, alpha=0.5)
ax3.legend(fontsize=12, framealpha=0.85)
ax3.set_title(r"Eficiencia vs ángulo del cono $\theta$ — Bar Strip Detector", fontsize=13)
ax3.grid(True, alpha=0.3)

# Límites físicos
theta_lateral = np.degrees(np.arctan(35.0 / 270.0))  # ~7.4°
tx_acc_cm = 50.0 * 200.0 / 300.5
theta_geom = np.degrees(np.arctan(tx_acc_cm / 200.0))  # ~9.4°

ax3.axvline(theta_lateral, color="goldenrod", ls=":", lw=2, alpha=0.8, zorder=1)
ax3.axvline(theta_geom, color="#8B008B", ls=":", lw=2, alpha=0.8, zorder=1)

ax3.annotate("Partícula sale\npor cara lateral\ndel Fe",
             xy=(theta_lateral, 50), xytext=(theta_lateral + 1.5, 70),
             fontsize=9, color="goldenrod",
             arrowprops=dict(arrowstyle="->", color="goldenrod", lw=1.5),
             bbox=dict(boxstyle="round", facecolor="#fff8dc", alpha=0.8))

ax3.annotate("Límite geométrico:\nbarras ±50 cm",
             xy=(theta_geom, 50), xytext=(theta_geom - 3, 70),
             fontsize=9, color="#8B008B",
             arrowprops=dict(arrowstyle="->", color="#8B008B", lw=1.5),
             bbox=dict(boxstyle="round", facecolor="#f0e6f6", alpha=0.8))

ax3.text(0.03, 0.05,
         "θ medido con columna ConeAngle\n"
         "(ángulo inicial traza vs eje z):\n"
         "  • Más preciso que reconstrucción\n"
         "    desde posiciones de barra\n"
         "  • Resolución limitada por paso\n"
         "    de barra de 5 cm (~0.7°)",
         transform=ax3.transAxes, fontsize=8.5, color="gray",
         va="bottom", ha="left", style="italic",
         bbox=dict(boxstyle="round", facecolor="white", alpha=0.8, edgecolor="gray"))

ax3.text(0.01, 1.012, info_str, transform=ax3.transAxes,
         fontsize=8.5, va="bottom", ha="left", color="#333333", style="italic")

fig3.savefig("img/eff_angulo_corregida.png", dpi=150, bbox_inches="tight")
print("Saved: img/eff_angulo_corregida.png")
plt.close(fig3)

# ===========================================================================
# RESUMEN NUMÉRICO COMPLETO
# ===========================================================================
print("\n" + "="*70)
print("RESUMEN COMPLETO — ANÁLISIS CORREGIDO")
print("="*70)
print(f"\nArchivos procesados: {N_RUNS} runs (run 0 a {N_RUNS-1})")
print(f"Eventos por run: 2000 (~1000 μ⁺ + ~1000 π⁺)")
print(f"Rango de momento: {p0_vals_GeV[0]:.4f} – {p0_vals_GeV[-1]:.4f} GeV/c")
print(f"Total eventos simulados: {N_RUNS * 2000}")

print(f"\n--- LANDAU (run 45, p0=1.0225 GeV/c) ---")
print(f"  Detectados: {tot_events} ({mu_events} μ⁺ + {pi_events} π⁺)")
print(f"  dE/dx μ⁺: media={mu_mean:.4f}, mediana={mu_med:.4f} MeV/mm")
print(f"  dE/dx π⁺: media={pi_mean:.4f}, mediana={pi_med:.4f} MeV/mm")

print(f"\n--- EFICIENCIA VS MOMENTO (selección) ---")
print(f"{'Run':>4} {'p0(GeV/c)':>12} {'ε(μ⁺)%':>10} {'Nμ':>5} {'ε(π⁺)%':>10} {'Nπ':>5}")
step = max(1, len(mu_recs) // 15)
for rn in sorted(set(r["run"] for r in mu_recs))[::step]:
    p = p0_vals_GeV[rn]
    mu_r = next(x for x in mu_recs if x["run"] == rn)
    pi_r = next(x for x in pi_recs if x["run"] == rn)
    print(f"{rn:4d} {p:12.4f} {mu_r['eff']*100:9.1f}% {mu_r['n_det']:5d} "
          f"{pi_r['eff']*100:9.1f}% {pi_r['n_det']:5d}")

mu_plateau = [r for r in mu_recs if r["p0"] > 1.1]
pi_plateau = [r for r in pi_recs if r["p0"] > 1.1]
if mu_plateau:
    mu_effs = np.array([r["eff"] for r in mu_plateau])
    print(f"\nμ⁺ meseta (>1.1 GeV/c): ε = {np.mean(mu_effs)*100:.1f}% ± {np.std(mu_effs)*100:.1f}%")
if pi_plateau:
    pi_effs = np.array([r["eff"] for r in pi_plateau])
    print(f"π⁺ meseta (>1.1 GeV/c): ε = {np.mean(pi_effs)*100:.1f}% ± {np.std(pi_effs)*100:.1f}%")

print("\n" + "="*70)
