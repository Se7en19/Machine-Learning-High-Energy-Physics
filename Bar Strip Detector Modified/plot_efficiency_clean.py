"""
Efficiency plots — Bar Strip Detector (mixed beam, 70 cm Fe + BC404 bars)
Efficiency vs initial momentum + Efficiency vs cone angle
Consistent with Landau plot (run 45, p0 = 1023 MeV/c)
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import uproot, glob, os, re

MIXED_DIR = "../Classifier/data/mixed"

plt.rcParams.update({
    "font.size": 12, "axes.titlesize": 13, "axes.labelsize": 12,
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

# Initial momenta: 80 runs log from 0.05 to 10.0 GeV/c
p0_vals_GeV = np.logspace(np.log10(0.05), np.log10(10.0), 80)

# ============================================================================
# EFFICIENCY VS MOMENTUM
# ============================================================================
files = sorted(glob.glob(os.path.join(MIXED_DIR, "output_run*.root")), key=run_num)

records_mu, records_pi = [], []
L_src = 2000.0  # mm

for fpath in files:
    rn = run_num(fpath)
    if rn >= 80: continue
    p0 = p0_vals_GeV[rn] * 1000  # MeV/c
    with uproot.open(fpath) as f:
        tree = best_cycle(f)
        if tree is None: continue
        d = tree.arrays(["particleID", "fEvent", "Momentum", "fX", "fY", "layerID"], library="np")

    for pid, name, rec_list in [(0, "mu+", records_mu), (1, "pi+", records_pi)]:
        mask = d["particleID"] == pid
        n = mask.sum()
        if n == 0:
            rec_list.append({"p0": p0, "eff": 0, "n_det": 0})
        else:
            n_det = len(np.unique(d["fEvent"][mask]))
            eff = n_det / 1000.0  # 1000 per species per run
            rec_list.append({"p0": p0, "eff": eff, "n_det": n_det})

# Sort by p0
mu_recs = sorted(records_mu, key=lambda r: r["p0"])
pi_recs = sorted(records_pi, key=lambda r: r["p0"])

fig1, ax1 = plt.subplots(figsize=(10, 6))

p_mu = np.array([r["p0"] for r in mu_recs])
e_mu = np.array([r["eff"] for r in mu_recs])
p_pi = np.array([r["p0"] for r in pi_recs])
e_pi = np.array([r["eff"] for r in pi_recs])

ax1.semilogx(p_mu, e_mu, "o-", color="steelblue", lw=2, ms=4, label=r"$\mu^+$")
ax1.semilogx(p_pi, e_pi, "s-", color="tomato", lw=2, ms=4, label=r"$\pi^+$")

ax1.axhline(0.5, color="gray", ls=":", lw=1, alpha=0.5)
ax1.set_xlabel(r"Momento inicial $p_0$  (MeV/c)", fontsize=12)
ax1.set_ylabel(r"$\varepsilon = N_{\mathrm{det}}\,/\,N_{\mathrm{gen}}$", fontsize=12)
ax1.set_ylim(-0.05, 1.10)
ax1.set_xlim(40, 11000)
ax1.legend(fontsize=11, framealpha=0.85)
ax1.set_title(r"Eficiencia de detección vs momento inicial — Bar Strip Detector", fontsize=13)

# Annotations for regimes
ax1.axvspan(40, 500, alpha=0.08, color="gray", label=None)
ax1.text(150, 0.60, "Muones: no penetran\n70 cm de Fe", fontsize=8, color="gray", ha="center", alpha=0.7)
ax1.axvspan(500, 1100, alpha=0.08, color="steelblue", label=None)
ax1.text(700, 0.60, "Régimen de\ntransición μ⁺", fontsize=8, color="steelblue", ha="center", alpha=0.7)
ax1.axvspan(1100, 10000, alpha=0.08, color="green", label=None)
ax1.text(5000, 0.60, "Meseta μ⁺\n(~85-92%)", fontsize=8, color="green", ha="center", alpha=0.7)

info = r"Geant4  |  Bar Strip Detector  |  70 cm Fe + BC404  |  FTFP\_BERT"
ax1.text(0.01, 1.008, info, transform=ax1.transAxes,
         fontsize=8.5, va="bottom", ha="left", color="#333333", style="italic")

fig1.tight_layout()
fig1.savefig("img/efficiency_vs_momentum_clean.png", dpi=150, bbox_inches="tight")
print("Saved: img/efficiency_vs_momentum_clean.png")
plt.close(fig1)

# ============================================================================
# EFFICIENCY VS CONE ANGLE
# ============================================================================
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
            ang_list.append(np.degrees(np.arctan2(r_face, L_src)))

theta_max = np.degrees(np.arctan2(np.sqrt(350**2 + 350**2), L_src))
n_bins = 18
bin_edges = np.linspace(0, theta_max, n_bins + 1)
bc = 0.5 * (bin_edges[:-1] + bin_edges[1:])

# Theoretical N_gen
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

fig2, ax2 = plt.subplots(figsize=(10, 6))

for angles, color, label in [
    (mu_angles, "steelblue", r"$\mu^+$"),
    (pi_angles, "tomato", r"$\pi^+$"),
]:
    if not angles: continue
    h_det, _ = np.histogram(angles, bins=bin_edges)
    eff = np.where(n_gen > 10, h_det / n_gen, np.nan)
    eff = np.clip(eff, 0, 1.05)
    ax2.plot(bc, eff, "o-", color=color, lw=2, ms=5, label=label)

ax2.set_xlabel(r"Ángulo del cono $\theta$ (°)", fontsize=12)
ax2.set_ylabel(r"$\varepsilon = N_{\mathrm{det}}\,/\,N_{\mathrm{gen}}$", fontsize=12)
ax2.set_xlim(0, theta_max + 0.3)
ax2.set_ylim(-0.05, 1.10)
ax2.axhline(0.5, color="gray", ls=":", lw=1, alpha=0.5)
ax2.legend(fontsize=11, framealpha=0.85)
ax2.set_title(r"Eficiencia vs ángulo del cono — Bar Strip Detector", fontsize=13)

# Regime lines
theta_lateral = np.degrees(np.arctan(35.0 / 270.0))
tx_acc_cm = 50.0 * 200.0 / 300.5
theta_geom = np.degrees(np.arctan(tx_acc_cm / 200.0))

ax2.axvline(theta_lateral, color="goldenrod", ls=":", lw=1.8, alpha=0.9, zorder=1)
ax2.axvline(theta_geom, color="#8B008B", ls=":", lw=1.8, alpha=0.9, zorder=1)
ax2.text(theta_lateral + 0.15, 1.09, f"cara lateral\n({theta_lateral:.1f}°)",
         fontsize=7.5, color="goldenrod", va="top", ha="left",
         transform=ax2.get_xaxis_transform())
ax2.text(theta_geom + 0.15, 1.09, f"límite geom.\n({theta_geom:.1f}°)",
         fontsize=7.5, color="#8B008B", va="top", ha="left",
         transform=ax2.get_xaxis_transform())

ax2.text(0.03, 0.06,
         r"$\theta$ reconstruido de posiciones de barras" + "\n"
         r"Capa 1 (Y) + Capa 2 (X) $ \rightarrow$ (t$_x$, t$_y$)" + "\n"
         r"resolución $\approx 0.7^\circ$ por coordenada",
         transform=ax2.transAxes, fontsize=8, color="gray",
         va="bottom", style="italic")
ax2.text(0.01, 1.008, info, transform=ax2.transAxes,
         fontsize=8.5, va="bottom", ha="left", color="#333333", style="italic")

fig2.tight_layout()
fig2.savefig("img/efficiency_vs_angle_clean.png", dpi=150, bbox_inches="tight")
print("Saved: img/efficiency_vs_angle_clean.png")
plt.close(fig2)

print("\n--- Resumen de eficiencias ---")
for r in mu_recs:
    if r["p0"] >= 1000 and r["p0"] <= 1050:
        print(f"μ⁺ @ p0={r['p0']:.0f} MeV/c: ε = {r['eff']:.1%}  ({r['n_det']} eventos detectados)")
        break
for r in pi_recs:
    if r["p0"] >= 1000 and r["p0"] <= 1050:
        print(f"π⁺ @ p0={r['p0']:.0f} MeV/c: ε = {r['eff']:.1%}  ({r['n_det']} eventos detectados)")
        break
mu_plateau = [r for r in mu_recs if r["p0"] > 1100]
pi_plateau = [r for r in pi_recs if r["p0"] > 1100]
if mu_plateau:
    print(f"μ⁺ meseta (>1100 MeV/c): ε media = {np.mean([r['eff'] for r in mu_plateau]):.1%}")
if pi_plateau:
    print(f"π⁺ meseta (>1100 MeV/c): ε media = {np.mean([r['eff'] for r in pi_plateau]):.1%}")
