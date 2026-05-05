"""
Landau distribution — Bar Strip Detector (mixed μ⁺/π⁺ beam, 70cm Fe + BC404 bars)
1k muon events + 1k pion events. Green vertical line at dE/dx = 0.5 MeV/mm (threshold).
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import uproot

MIXED_FILE = "../Classifier/data/mixed/output_run0.root"

plt.rcParams.update({
    "font.size": 12,
    "axes.titlesize": 14,
    "axes.labelsize": 12,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 10,
    "figure.dpi": 150,
    "figure.facecolor": "white",
})

with uproot.open(MIXED_FILE) as f:
    tree = f["Hits"]
    d = tree.arrays(["fdEdx", "particleID", "Momentum", "layerID", "fEvent"], library="np")

mu_fdEdx = d["fdEdx"][d["particleID"] == 0]
pi_fdEdx = d["fdEdx"][d["particleID"] == 1]

DEDX_MIN, DEDX_MAX = 0.01, 5.0
mask_mu = (mu_fdEdx >= DEDX_MIN) & (mu_fdEdx <= DEDX_MAX)
mask_pi = (pi_fdEdx >= DEDX_MIN) & (pi_fdEdx <= DEDX_MAX)
mu_fdEdx = mu_fdEdx[mask_mu]
pi_fdEdx = pi_fdEdx[mask_pi]

n_mu_events = len(np.unique(d["fEvent"][d["particleID"] == 0]))
n_pi_events = len(np.unique(d["fEvent"][d["particleID"] == 1]))

print(f"Eventos generados: ~1000 mu+ + ~1000 pi+")
print(f"Eventos con hits en detector: {n_mu_events} mu+, {n_pi_events} pi+")
print(f"Hits totales tras filtro: {len(mu_fdEdx)} mu+, {len(pi_fdEdx)} pi+")
print(f"dE/dx medio: mu+={np.mean(mu_fdEdx):.4f}, pi+={np.mean(pi_fdEdx):.4f} MeV/mm")
print(f"dE/dx mediana: mu+={np.median(mu_fdEdx):.4f}, pi+={np.median(pi_fdEdx):.4f} MeV/mm")

fig, ax = plt.subplots(figsize=(10, 6))

bins = np.linspace(DEDX_MIN, DEDX_MAX, 100)

ax.hist(mu_fdEdx, bins=bins, histtype="stepfilled",
        color="steelblue", alpha=0.5, density=True,
        label=rf"$\mu^+$  ({n_mu_events} events, {len(mu_fdEdx)} hits)")
ax.hist(mu_fdEdx, bins=bins, histtype="step",
        color="steelblue", lw=2, density=True)
ax.hist(pi_fdEdx, bins=bins, histtype="stepfilled",
        color="tomato", alpha=0.4, density=True,
        label=rf"$\pi^+$  ({n_pi_events} events, {len(pi_fdEdx)} hits)")
ax.hist(pi_fdEdx, bins=bins, histtype="step",
        color="tomato", lw=2, density=True)

THRESHOLD = 0.5
ax.axvline(THRESHOLD, color="green", lw=2.5, ls="--",
           label=rf"Umbral dE/dx = {THRESHOLD} MeV/mm")

# MPV markers
mu_mpv = bins[np.argmax(np.histogram(mu_fdEdx, bins=bins, density=True)[0])]
pi_mpv = bins[np.argmax(np.histogram(pi_fdEdx, bins=bins, density=True)[0])]
ax.axvline(mu_mpv, color="steelblue", lw=1, ls=":", alpha=0.5)
ax.axvline(pi_mpv, color="tomato", lw=1, ls=":", alpha=0.5)
ax.text(mu_mpv + 0.02, 0.5, f"MPV$_{{\\mu}}$={mu_mpv:.3f}", fontsize=8, color="steelblue", rotation=90)
ax.text(pi_mpv + 0.02, 0.3, f"MPV$_{{\\pi}}$={pi_mpv:.3f}", fontsize=8, color="tomato", rotation=90)

ax.set_xlabel(r"$dE/dx$ por paso  (MeV/mm)", fontsize=12)
ax.set_ylabel("Densidad (normalizada)", fontsize=12)
ax.set_yscale("log")
ax.set_xlim(DEDX_MIN, DEDX_MAX)
ax.legend(fontsize=9, loc="upper right")

# Annotation about the threshold
ax.text(0.55, 0.82,
        "Línea verde: dE/dx = 0.5 MeV/mm\n"
        "Eventos a la derecha = pasos con\n"
        r"rayos $\delta$ energéticos en barras de 1 cm",
        transform=ax.transAxes, fontsize=9, color="green",
        va="top", ha="left",
        bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))

# Info
info = (r"Geant4  |  Bar Strip Detector  |  "
        r"$\mu^+$/$\pi^+$ mixed beam  |  "
        r"FTFP\_BERT  |  70 cm Fe  |  "
        r"$p_0$ = 1.0 GeV/c")
ax.text(0.01, 1.008, info, transform=ax.transAxes,
        fontsize=8.5, va="bottom", ha="left", color="#333333", style="italic")

fig.tight_layout()
fig.savefig("img/landau_bars_threshold.png", dpi=150, bbox_inches="tight")
print("\nSaved: img/landau_bars_threshold.png")
plt.close(fig)
