"""
Análisis completo — Bar Strip Detector
========================================
1) Distribución de Landau (run 45, ~1 GeV/c) con umbral dE/dx = 0.5 MeV/mm
2) Eficiencia vs momento (en GeV/c) para μ⁺ y π⁺
3) Eficiencia vs ángulo del cono
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

fig, axes = plt.subplots(1, 3, figsize=(18, 6))

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

ax1 = axes[0]
bins = np.linspace(0.01, 5.0, 100)
ax1.hist(mu_dedx, bins=bins, histtype="stepfilled",
         color="steelblue", alpha=0.5, density=True,
         label=rf"$\mu^+$  ({mu_events} eventos)")
ax1.hist(mu_dedx, bins=bins, histtype="step",
         color="steelblue", lw=2, density=True)
ax1.hist(pi_dedx, bins=bins, histtype="stepfilled",
         color="tomato", alpha=0.4, density=True,
         label=rf"$\pi^+$  ({pi_events} eventos)")
ax1.hist(pi_dedx, bins=bins, histtype="step",
         color="tomato", lw=2, density=True)

ax1.axvline(0.5, color="green", lw=2.5, ls="--",
            label=r"Umbral dE/dx = 0.5 MeV/mm")
ax1.set_xlabel(r"$dE/dx$ por paso  (MeV/mm)")
ax1.set_ylabel("Densidad (normalizada)")
ax1.set_yscale("log")
ax1.set_xlim(0.01, 5.0)
ax1.legend(fontsize=9, loc="upper right")
ax1.set_title(r"Distribución de Landau ($p_0 \approx 1$ GeV/c)")

info = (r"Geant4  |  Bar Strip Detector  |  "
        r"mixed $\mu^+$/$\pi^+$  |  $p_0 \approx 1.0$ GeV/c  |  2000 eventos")
ax1.text(0.01, 1.008, info, transform=ax1.transAxes,
         fontsize=8, va="bottom", ha="left", color="#333333", style="italic")

# ============================================================
# 2) EFICIENCIA VS MOMENTO (en GeV/c)
# ============================================================
files = sorted(glob.glob(os.path.join(MIXED_DIR, "output_run*.root")), key=run_num)
p0_vals_GeV = np.logspace(np.log10(0.05), np.log10(10.0), 80)

records_mu, records_pi = [], []

for fpath in files:
    rn = run_num(fpath)
    if rn >= 80: continue
    p0 = p0_vals_GeV[rn]  # GeV/c
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

ax2 = axes[1]
p_mu = np.array([r["p0"] for r in mu_recs])
e_mu = np.array([r["eff"] for r in mu_recs])
p_pi = np.array([r["p0"] for r in pi_recs])
e_pi = np.array([r["eff"] for r in pi_recs])

ax2.semilogx(p_mu, e_mu, "o-", color="steelblue", lw=2, ms=4, label=r"$\mu^+$")
ax2.semilogx(p_pi, e_pi, "s-", color="tomato", lw=2, ms=4, label=r"$\pi^+$")

ax2.axhline(0.5, color="gray", ls=":", lw=1, alpha=0.5)
ax2.set_xlabel(r"Momento inicial $p_0$  (GeV/c)")
ax2.set_ylabel(r"$\varepsilon = N_{\mathrm{det}}\,/\,N_{\mathrm{gen}}$")
ax2.set_ylim(-0.05, 1.10)
ax2.set_xlim(0.04, 11)
ax2.legend(fontsize=11, framealpha=0.85)
ax2.set_title(r"Eficiencia vs momento — Bar Strip Detector")

# Regímenes
ax2.axvspan(0.04, 0.5, alpha=0.08, color="gray")
ax2.text(0.15, 0.60, "Muones: no penetran\n70 cm de Fe", fontsize=8, color="gray", ha="center", alpha=0.7)
ax2.axvspan(0.5, 1.1, alpha=0.08, color="steelblue")
ax2.text(0.7, 0.60, "Régimen de\ntransición μ⁺", fontsize=8, color="steelblue", ha="center", alpha=0.7)
ax2.axvspan(1.1, 10, alpha=0.08, color="green")
ax2.text(5.0, 0.60, "Meseta μ⁺\n(~85-92%)", fontsize=8, color="green", ha="center", alpha=0.7)

info = r"Geant4  |  Bar Strip Detector  |  70 cm Fe + BC404  |  FTFP\_BERT"
ax2.text(0.01, 1.008, info, transform=ax2.transAxes,
         fontsize=8.5, va="bottom", ha="left", color="#333333", style="italic")

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
            ang_list.append(np.degrees(np.arctan2(r_face, L_src)))

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

ax3 = axes[2]
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

ax3.text(0.01, 1.008, info, transform=ax3.transAxes,
         fontsize=8.5, va="bottom", ha="left", color="#333333", style="italic")

# ============================================================
fig.tight_layout()
fig.savefig("img/analisis_completo.png", dpi=150, bbox_inches="tight")
print("Saved: img/analisis_completo.png")
plt.close(fig)

# ============================================================
# RESUMEN
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
