"""
Landau distribution — Bar Strip Detector
Run 45 (p0 = 1023 MeV/c, ~1 GeV/c): 2000 events (~1000 mu+ + ~1000 pi+ generated)
Green line at dE/dx = 0.5 MeV/mm
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import uproot, os, re

MIXED_DIR = "../Classifier/data/mixed"

plt.rcParams.update({
    "font.size": 12, "axes.titlesize": 14, "axes.labelsize": 12,
    "xtick.labelsize": 10, "ytick.labelsize": 10, "legend.fontsize": 10,
    "figure.dpi": 150, "figure.facecolor": "white",
})

def best_cycle(f, tree_name="Hits"):
    keys = [k for k in f.keys(cycle=True) if k.split(";")[0] == tree_name]
    if not keys: return None
    best = max(keys, key=lambda k: f[k].num_entries)
    return f[best] if f[best].num_entries > 0 else None

# Run 45 → p0 ≈ 1023 MeV/c (closest to 1 GeV/c in the sweep)
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

print(f"=== Run 45: p0 = 1023 MeV/c (2000 events = ~1000 mu+ + ~1000 pi+ generated) ===")
print(f"Mu+: {mu_events} eventos detectados, {len(mu_dedx)} hits")
print(f"Pi+: {pi_events} eventos detectados, {len(pi_dedx)} hits")
print(f"dE/dx mu+ medio={np.mean(mu_dedx):.4f}, mediana={np.median(mu_dedx):.4f}")
print(f"dE/dx pi+ medio={np.mean(pi_dedx):.4f}, mediana={np.median(pi_dedx):.4f}")

fig, ax = plt.subplots(figsize=(10, 6))

bins = np.linspace(0.01, 5.0, 100)

ax.hist(mu_dedx, bins=bins, histtype="stepfilled",
        color="steelblue", alpha=0.5, density=True,
        label=rf"$\mu^+$  ({mu_events} eventos, {len(mu_dedx)} hits)")
ax.hist(mu_dedx, bins=bins, histtype="step",
        color="steelblue", lw=2, density=True)
ax.hist(pi_dedx, bins=bins, histtype="stepfilled",
        color="tomato", alpha=0.4, density=True,
        label=rf"$\pi^+$  ({pi_events} eventos, {len(pi_dedx)} hits)")
ax.hist(pi_dedx, bins=bins, histtype="step",
        color="tomato", lw=2, density=True)

ax.axvline(0.5, color="green", lw=2.5, ls="--",
           label=r"Umbral dE/dx = 0.5 MeV/mm")

# MPV markers
mu_mpv = bins[np.argmax(np.histogram(mu_dedx, bins=bins, density=True)[0])]
pi_mpv = bins[np.argmax(np.histogram(pi_dedx, bins=bins, density=True)[0])]
ax.axvline(mu_mpv, color="steelblue", lw=1, ls=":", alpha=0.5)
ax.axvline(pi_mpv, color="tomato", lw=1, ls=":", alpha=0.5)
ax.text(mu_mpv + 0.02, 0.5, f"MPV$_{{\\mu}}$={mu_mpv:.3f}", fontsize=8, color="steelblue", rotation=90)
ax.text(pi_mpv + 0.02, 0.3, f"MPV$_{{\\pi}}$={pi_mpv:.3f}", fontsize=8, color="tomato", rotation=90)

ax.set_xlabel(r"$dE/dx$ por paso  (MeV/mm)", fontsize=12)
ax.set_ylabel("Densidad (normalizada)", fontsize=12)
ax.set_yscale("log")
ax.set_xlim(0.01, 5.0)
ax.legend(fontsize=9, loc="upper right")

ax.text(0.55, 0.82,
        "Línea verde: dE/dx = 0.5 MeV/mm\n"
        "Eventos a la derecha = pasos con\n"
        r"rayos $\delta$ energéticos (barras 1 cm)",
        transform=ax.transAxes, fontsize=9, color="green",
        va="top", ha="left",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))

info = (r"Geant4  |  Bar Strip Detector  |  "
        r"mixed $\mu^+$/$\pi^+$  |  "
        r"FTFP\_BERT  |  70 cm Fe  |  "
        r"$p_0$ = 1.0 GeV/c")
ax.text(0.01, 1.008, info, transform=ax.transAxes,
        fontsize=8.5, va="bottom", ha="left", color="#333333", style="italic")

fig.tight_layout()
fig.savefig("img/landau_fresh_threshold.png", dpi=150, bbox_inches="tight")
print("\nSaved: img/landau_fresh_threshold.png")
plt.close(fig)
