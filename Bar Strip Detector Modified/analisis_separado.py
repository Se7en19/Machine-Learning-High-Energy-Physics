"""
Tres gráficas separadas — Bar Strip Detector
=============================================
1) landau_separada.png   — Distribución de Landau (run 45, ~1 GeV/c) con umbral 0.5
2) eff_momento_separada.png — Eficiencia vs momento (GeV/c)
3) eff_angulo_separada.png   — Eficiencia vs ángulo del cono
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
info_str = r"Geant4  |  Bar Strip Detector  |  70 cm Fe + BC404  |  FTFP\_BERT"

# ============================================================
# 1) LANDAU — run 45 (p0 ≈ 1023 MeV/c, 2000 events mixed)
# ============================================================
fpath = os.path.join(MIXED_DIR, "output_run45.root")
with uproot.open(fpath) as f:
    tree = best_cycle(f)
    d = tree.arrays(["fdEdx", "particleID", "fEvent"], library="np")

mu_mask = d["particleID"] == 0
pi_mask = d["particleID"] == 1
mu_dedx = d["fdEdx"][mu_mask]
pi_dedx = d["fdEdx"][pi_mask]

mu_events = len(np.unique(d["fEvent"][mu_mask]))
pi_events = len(np.unique(d["fEvent"][pi_mask]))

fig1, ax1 = plt.subplots(figsize=(10, 6))
bins = np.linspace(0.01, 5.0, 100)

ax1.hist(mu_dedx, bins=bins, histtype="stepfilled",
         color="steelblue", alpha=0.5, density=True,
         label=rf"$\mu^+$  ({mu_events} eventos, 2000 generados)")
ax1.hist(mu_dedx, bins=bins, histtype="step",
         color="steelblue", lw=2, density=True)
ax1.hist(pi_dedx, bins=bins, histtype="stepfilled",
         color="tomato", alpha=0.4, density=True,
         label=rf"$\pi^+$  ({pi_events} eventos, 2000 generados)")
ax1.hist(pi_dedx, bins=bins, histtype="step",
         color="tomato", lw=2, density=True)

ax1.axvline(0.5, color="green", lw=2.5, ls="--",
            label=r"Umbral dE/dx = 0.5 MeV/mm")
ax1.set_xlabel(r"$dE/dx$ por paso  (MeV/mm)")
ax1.set_ylabel("Densidad (normalizada)")
ax1.set_yscale("log")
ax1.set_xlim(0.01, 5.0)
ax1.legend(fontsize=10, loc="upper right")
ax1.set_title(r"Distribución de Landau — haz mixto $\mu^+$/$\pi^+$  ($p_0 \approx 1$ GeV/c)")

mu_pct = np.sum(mu_dedx > 0.5) / len(mu_dedx) * 100
pi_pct = np.sum(pi_dedx > 0.5) / len(pi_dedx) * 100

ax1.text(0.55, 0.82,
         "Línea verde: dE/dx = 0.5 MeV/mm\n"
         f"{mu_pct:.1f}% hits μ⁺ y {pi_pct:.1f}% hits π⁺\n"
         "superan este umbral (δ-rays)",
         transform=ax1.transAxes, fontsize=9, color="green",
         va="top", ha="left",
         bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))

info1 = (f"Geant4  |  ~1000 μ⁺ + ~1000 π⁺ generados  |  "
         f"{mu_events} μ⁺ detectados + {pi_events} π⁺ detectados = {mu_events+pi_events} total  |  "
         f"$p_0 \\approx 1.0$ GeV/c")
ax1.text(0.01, 1.008, info1, transform=ax1.transAxes,
         fontsize=8.5, va="bottom", ha="left", color="#333333", style="italic")

fig1.tight_layout()
fig1.savefig("img/landau_separada.png", dpi=150, bbox_inches="tight")
print("Saved: img/landau_separada.png")
plt.close(fig1)

# ============================================================
# 2) EFICIENCIA VS MOMENTO (GeV/c)
# ============================================================
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
            rec_list.append({"p0": p0, "eff": 0, "n_det": 0})
        else:
            n_det = len(np.unique(d["fEvent"][mask]))
            eff = n_det / 1000.0
            rec_list.append({"p0": p0, "eff": eff, "n_det": n_det})

mu_recs = sorted(records_mu, key=lambda r: r["p0"])
pi_recs = sorted(records_pi, key=lambda r: r["p0"])

N_gen = 1000.0
p_mu = np.array([r["p0"] for r in mu_recs])
e_mu = np.array([r["eff"] for r in mu_recs])
n_mu = np.array([r["n_det"] for r in mu_recs])
p_pi = np.array([r["p0"] for r in pi_recs])
e_pi = np.array([r["eff"] for r in pi_recs])
n_pi = np.array([r["n_det"] for r in pi_recs])

# Binomial error bars
err_mu = np.sqrt(e_mu * (1 - e_mu) / N_gen)
err_pi = np.sqrt(e_pi * (1 - e_pi) / N_gen)

fig2, ax2 = plt.subplots(figsize=(10, 6))

ax2.errorbar(p_mu, e_mu, yerr=err_mu, fmt="o-", color="steelblue",
             lw=2, ms=4, capsize=3, capthick=1, label=r"$\mu^+$")
ax2.errorbar(p_pi, e_pi, yerr=err_pi, fmt="s-", color="tomato",
             lw=2, ms=5, capsize=3, capthick=1, label=r"$\pi^+$")

ax2.axhline(0.5, color="gray", ls=":", lw=1, alpha=0.5)
ax2.set_xlabel(r"Momento inicial $p_0$  (GeV/c)")
ax2.set_ylabel(r"$\varepsilon = N_{\mathrm{det}}\,/\,N_{\mathrm{gen}}$")
ax2.set_ylim(-0.05, 1.10)
ax2.set_xlim(0.04, 11)
ax2.legend(fontsize=11, framealpha=0.85, loc="upper left")
ax2.set_title(r"Eficiencia de detección vs momento — Bar Strip Detector")
ax2.grid(True, alpha=0.3)

# Regímenes
ax2.axvspan(0.04, 0.5, alpha=0.08, color="gray")
ax2.text(0.15, 0.60, "Régimen I:\nμ⁺ no penetran\n70 cm Fe", fontsize=8.5, color="gray", ha="center", alpha=0.7)
ax2.axvspan(0.5, 1.1, alpha=0.08, color="steelblue")
ax2.text(0.74, 0.60, "Régimen II:\ntransición μ⁺\n(subida sigmoidea)", fontsize=8.5, color="steelblue", ha="center", alpha=0.7)
ax2.axvspan(1.1, 10, alpha=0.08, color="green")
ax2.text(4.0, 0.60, "Régimen III:\nmeseta μ⁺\n(~89%)", fontsize=8.5, color="green", ha="center", alpha=0.7)

# Inset zoom for pions
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
ax_inset = inset_axes(ax2, width="30%", height="30%", loc="lower right",
                      bbox_to_anchor=(0.02, 0.02, 0.5, 0.5),
                      bbox_transform=ax2.transAxes)
ax_inset.errorbar(p_pi, e_pi*100, yerr=err_pi*100, fmt="s-", color="tomato",
                  lw=1.5, ms=3, capsize=2, capthick=1)
ax_inset.set_xscale("log")
ax_inset.set_xlim(0.04, 11)
ax_inset.set_ylim(0, 20)
ax_inset.set_ylabel(r"$\varepsilon_{\pi^+}$ (%)", fontsize=8)
ax_inset.set_xlabel(r"$p_0$ (GeV/c)", fontsize=8)
ax_inset.tick_params(labelsize=7)
ax_inset.grid(True, alpha=0.3)
ax_inset.axhline(7.5, color="gray", ls="--", lw=0.8, alpha=0.5)
ax_inset.text(0.1, 8.5, "media ~7.5%", fontsize=6, color="gray", alpha=0.7)

ax2.text(0.01, 1.008, info_str, transform=ax2.transAxes,
         fontsize=8.5, va="bottom", ha="left", color="#333333", style="italic")

fig2.tight_layout()
fig2.savefig("img/eff_momento_separada.png", dpi=150, bbox_inches="tight")
print("Saved: img/eff_momento_separada.png")
plt.close(fig2)

# ============================================================
# 3) EFICIENCIA VS ÁNGULO DEL CONO
# ============================================================
L_src = 2000.0
mu_angles, pi_angles = [], []
z_L1, z_L2 = 1005.0, 1035.0
scale_L1 = (z_L1 + L_src) / L_src
scale_L2 = (z_L2 + L_src) / L_src

for fpath in files:
    rn = run_num(fpath)
    if rn >= 80: continue
    with uproot.open(fpath) as f:
        tree = best_cycle(f)
        if tree is None: continue
        d = tree.arrays(["fEvent", "fX", "fY", "layerID", "particleID"], library="np")

    for pid_val, ang_list in [(0, mu_angles), (1, pi_angles)]:
        mask_pid = d["particleID"] == pid_val
        if mask_pid.sum() == 0: continue
        fX_p, fY_p = d["fX"][mask_pid], d["fY"][mask_pid]
        lay_p, ev_p = d["layerID"][mask_pid], d["fEvent"][mask_pid]
        for ev in np.unique(ev_p):
            m_ev = ev_p == ev
            lay_ev, fX_ev, fY_ev = lay_p[m_ev], fX_p[m_ev], fY_p[m_ev]
            m0, m1 = lay_ev == 0, lay_ev == 1
            if m0.sum() == 0 or m1.sum() == 0:
                continue
            ty_mm = np.median(fY_ev[m0]) / scale_L1
            tx_mm = np.median(fX_ev[m1]) / scale_L2
            r_face = np.sqrt(tx_mm**2 + ty_mm**2)
            mu_angles.append(np.degrees(np.arctan2(r_face, L_src)))

theta_max = np.degrees(np.arctan2(np.sqrt(350**2 + 350**2), L_src))
n_bins = 18
bin_edges = np.linspace(0, theta_max, n_bins + 1)
bc = 0.5 * (bin_edges[:-1] + bin_edges[1:])

rng = np.random.default_rng(42)
N_mc = 2_000_000
tx_mc = rng.uniform(-350, 350, N_mc)
ty_mc = rng.uniform(-350, 350, N_mc)
theta_mc = np.degrees(np.arctan2(np.sqrt(tx_mc**2 + ty_mc**2), L_src))
h_mc, _ = np.histogram(theta_mc, bins=bin_edges)
p_theory = h_mc / h_mc.sum()
n_runs = len(files)
n_total_per_species = n_runs * 1000.0
n_gen = p_theory * n_total_per_species

fig3, ax3 = plt.subplots(figsize=(10, 6))

for angles, color, label in [
    (mu_angles, "steelblue", r"$\mu^+$"),
    (pi_angles, "tomato", r"$\pi^+$"),
]:
    if not angles: continue
    h_det, _ = np.histogram(angles, bins=bin_edges)
    eff = np.where(n_gen > 10, h_det / n_gen, np.nan)
    eff = np.clip(eff, 0, 1.05)
    ax3.plot(bc, eff, "o-", color=color, lw=2, ms=5, label=label)

ax3.set_xlabel(r"Ángulo del cono $\theta$ (°)")
ax3.set_ylabel(r"$\varepsilon = N_{\mathrm{det}}\,/\,N_{\mathrm{gen}}$")
ax3.set_xlim(0, theta_max + 0.3)
ax3.set_ylim(-0.05, 1.10)
ax3.axhline(0.5, color="gray", ls=":", lw=1, alpha=0.5)
ax3.legend(fontsize=11, framealpha=0.85)
ax3.set_title(r"Eficiencia vs ángulo del cono — Bar Strip Detector")

theta_lateral = np.degrees(np.arctan(35.0 / 270.0))
tx_acc_cm = 50.0 * 200.0 / 300.5
theta_geom = np.degrees(np.arctan(tx_acc_cm / 200.0))

ax3.axvline(theta_lateral, color="goldenrod", ls=":", lw=1.8, alpha=0.9, zorder=1)
ax3.axvline(theta_geom, color="#8B008B", ls=":", lw=1.8, alpha=0.9, zorder=1)
ax3.text(theta_lateral + 0.15, 1.09, f"cara lateral\n({theta_lateral:.1f}°)",
         fontsize=7.5, color="goldenrod", va="top", ha="left",
         transform=ax3.get_xaxis_transform())
ax3.text(theta_geom + 0.15, 1.09, f"límite geom.\n({theta_geom:.1f}°)",
         fontsize=7.5, color="#8B008B", va="top", ha="left",
         transform=ax3.get_xaxis_transform())

ax3.text(0.03, 0.06,
         r"$\theta$ reconstruido de posiciones de barras" + "\n"
         r"Capa 1 (Y) + Capa 2 (X) $ \rightarrow$ (t$_x$, t$_y$)" + "\n"
         r"resolución $\approx 0.7^\circ$ por coordenada",
         transform=ax3.transAxes, fontsize=8, color="gray",
         va="bottom", style="italic")

ax3.text(0.01, 1.008, info_str, transform=ax3.transAxes,
         fontsize=8.5, va="bottom", ha="left", color="#333333", style="italic")

fig3.tight_layout()
fig3.savefig("img/eff_angulo_separada.png", dpi=150, bbox_inches="tight")
print("Saved: img/eff_angulo_separada.png")
plt.close(fig3)

# ============================================================
print("\n" + "="*60)
print("RESUMEN DE EFICIENCIAS")
print("="*60)
for r in mu_recs:
    if 0.95 <= r["p0"] <= 1.05:
        print(f"μ⁺ @ p0={r['p0']:.3f} GeV/c: ε = {r['eff']:.1%}  ({r['n_det']} eventos)")
        break
for r in pi_recs:
    if 0.95 <= r["p0"] <= 1.05:
        print(f"π⁺ @ p0={r['p0']:.3f} GeV/c: ε = {r['eff']:.1%}  ({r['n_det']} eventos)")
        break
mu_plateau = [r for r in mu_recs if r["p0"] > 1.1]
pi_plateau = [r for r in pi_recs if r["p0"] > 1.1]
if mu_plateau:
    print(f"μ⁺ meseta (>1.1 GeV/c): ε media = {np.mean([r['eff'] for r in mu_plateau]):.1%}")
if pi_plateau:
    print(f"π⁺ meseta (>1.1 GeV/c): ε media = {np.mean([r['eff'] for r in pi_plateau]):.1%}")
