"""
Bethe-Bloch curve plotter — mu+ vs pi+  (Geant4, BC404 bar strip detector)

Genera 7 gráficas:
  1. dE/dx vs βγ          (histograma 2D log-log, mu+ vs pi+ lado a lado)
  2. dE/dx vs β           (dependencia en velocidad)
  3. dE/dx vs p           (GeV/c, estilo PID, Y lineal)
  4. Distribución de Landau  (dE/dx por paso)
  5. Overlay βγ           (mediana dE/dx, mu+ vs pi+ en un panel)
  6. PID combinado        (histograma 2D mu+ vs pi+ en momento, panel único)
  7. Hits por capa        (distribución de hits en Capa 1 vs Capa 2)

Uso:
    python plot_bethe_bloch.py \
        --muon "../Classifier/data/muon_bars/output_run*.root" \
        --pion "../Classifier/data/pion_bars/output_run*.root" \
        --out  img/
"""

import argparse
import glob
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import uproot

# ── Particle masses  (MeV/c²) ───────────────────────────────────────────────
M_MU  = 105.6583755
M_PI  = 139.5702
M_E   = 0.51099895

# ── BC404 plastic scintillator (G4_PLASTIC_SC_VINYLTOLUENE) ─────────────────
# C9H10: rho=1.032 g/cm³, I=64.7 eV, Z/A=0.5424
AIR = dict(Z_over_A=0.5424, I=64.7e-6, rho=1.032)
K   = 0.307075   # MeV·cm²/mol

plt.rcParams.update({
    "font.size": 12,
    "axes.titlesize": 13,
    "axes.labelsize": 12,
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
    "legend.fontsize": 10,
    "figure.dpi": 150,
    "axes.facecolor": "white",
    "figure.facecolor": "white",
    "axes.grid": False,
})

CMAP = "jet"

DEDX_MIN = 0.01   # MeV/mm
DEDX_MAX = 5.0    # MeV/mm


# ============================================================================
# Sternheimer density-effect correction δ(βγ) for BC404
# ============================================================================
_STERN_BC404 = dict(C=-3.7936, x0=0.1496, x1=2.4815, a=0.15018, m=3.4083, d0=0.00)

def density_effect(bg, stern=_STERN_BC404):
    bg  = np.asarray(bg, dtype=float)
    x   = np.log10(bg)
    C, x0, x1, a, m, d0 = (stern[k] for k in ("C","x0","x1","a","m","d0"))
    delta = np.where(
        x >= x1,
        2.0 * np.log(10) * x + C,
        np.where(
            x >= x0,
            2.0 * np.log(10) * x + C + a * (x1 - x)**m,
            d0 * 10.0**(2.0 * (x - x0))
        )
    )
    return delta


# ============================================================================
# Analytic curves  [MeV/mm] vs βγ
# ============================================================================
def bethe_bloch(bg, mass=M_MU, mat=AIR):
    bg   = np.asarray(bg, dtype=float)
    b2   = bg**2 / (1.0 + bg**2)
    g    = np.sqrt(1.0 + bg**2)
    Tmax = (2.0 * M_E * bg**2) / (1.0 + 2.0 * g * M_E / mass + (M_E / mass)**2)
    logA = np.where(Tmax > 0, 2.0 * M_E * b2 * g**2 * Tmax / mat["I"]**2, np.nan)
    delta = density_effect(bg)
    dedx = K * mat["Z_over_A"] / b2 * (0.5 * np.log(logA) - b2 - delta / 2.0)
    dedx = np.where(dedx > 0, dedx, np.nan)
    return dedx * mat["rho"] / 10.0   # MeV/mm


def landau_mpv(bg, x_mm=10.0, mat=AIR):
    """
    Landau MPV of dE/dx for a step of x_mm (mm).
    Para las barras de 1 cm de grosor, x_mm = 10.0 (incidencia normal).
    """
    bg   = np.asarray(bg, dtype=float)
    b2   = bg**2 / (1.0 + bg**2)
    bg2  = bg**2
    x_cm = x_mm / 10.0
    xi   = 0.5 * K * mat["Z_over_A"] * mat["rho"] * x_cm / b2
    delta = density_effect(bg)
    mpv = (xi / x_cm) * (np.log(2.0 * M_E * bg2 / mat["I"])
                         + np.log(xi / mat["I"])
                         + 0.2 - b2 - delta)
    mpv /= 10.0   # MeV/mm
    return np.where(mpv > 0, mpv, np.nan)


# ============================================================================
# Data loading
# ============================================================================
def _best_cycle(f, tree_name="Hits"):

    keys = [k for k in f.keys(cycle=True) if k.split(";")[0] == tree_name]
    if not keys:
        return None
    best_key = max(keys, key=lambda k: f[k].num_entries)
    return f[best_key] if f[best_key].num_entries > 0 else None


def load_all_hits(glob_pattern: str, mass: float, label: str):
    """
    Carga hits de todos los ROOT files que coincidan con glob_pattern.
    Incluye las columnas layerID y barID propias del detector de barras.
    """
    files = sorted(glob.glob(glob_pattern))
    if not files:
        raise FileNotFoundError(f"[ERROR] No files match: {glob_pattern}")

    print(f"\n--- {label} : {len(files)} file(s) matching '{glob_pattern}' ---")
    branches = ["fdEdx", "Momentum", "Ekin", "layerID", "barID"]

    parts = {b: [] for b in branches}
    total = 0
    for fpath in files:
        with uproot.open(fpath) as f:
            tree = _best_cycle(f)
            if tree is None:
                print(f"  SKIP {os.path.basename(fpath)} (no entries)")
                continue
            n = tree.num_entries
            d = tree.arrays(branches, library="np")
            for b in branches:
                parts[b].append(d[b])
            total += n
            print(f"  {os.path.basename(fpath):30s}  {n:7,} hits")

    print(f"  Total hits loaded: {total:,}")
    if total == 0:
        raise RuntimeError(f"[ERROR] No entries found for '{glob_pattern}'.")

    data = {b: np.concatenate(parts[b]) for b in branches}

    mask = (data["fdEdx"] >= DEDX_MIN) & (data["fdEdx"] <= DEDX_MAX) & (data["Momentum"] > 0)
    data = {b: data[b][mask] for b in branches}

    data["bg"]   = data["Momentum"] / mass
    data["beta"] = data["Momentum"] / np.sqrt(data["Momentum"]**2 + mass**2)

    print(f"  Hits en rango dEdx [{DEDX_MIN:.0e}, {DEDX_MAX:.0e}]: {mask.sum():,}")
    return data


def load_mixed_hits(glob_pattern: str):
    """
    Carga el dataset mixto (particleID: 0=mu+, 1=pi+) y lo separa por especie.
    """
    files = sorted(glob.glob(glob_pattern))
    if not files:
        raise FileNotFoundError(f"[ERROR] No files match: {glob_pattern}")

    print(f"\n--- mixed : {len(files)} file(s) matching '{glob_pattern}' ---")
    branches = ["fdEdx", "Momentum", "Ekin", "layerID", "barID", "particleID"]

    parts = {b: [] for b in branches}
    total = 0
    for fpath in files:
        with uproot.open(fpath) as f:
            tree = _best_cycle(f)
            if tree is None:
                print(f"  SKIP {os.path.basename(fpath)} (no entries)")
                continue
            n = tree.num_entries
            d = tree.arrays(branches, library="np")
            for b in branches:
                parts[b].append(d[b])
            total += n
            print(f"  {os.path.basename(fpath):30s}  {n:7,} hits")

    print(f"  Total hits cargados: {total:,}")
    if total == 0:
        raise RuntimeError(f"[ERROR] No entries found for '{glob_pattern}'.")

    raw = {b: np.concatenate(parts[b]) for b in branches}

    def _make_species(pid_val, mass, label):
        mask_pid = raw["particleID"] == pid_val
        d = {b: raw[b][mask_pid] for b in branches if b != "particleID"}
        mask_range = (d["fdEdx"] >= DEDX_MIN) & (d["fdEdx"] <= DEDX_MAX) & (d["Momentum"] > 0)
        d = {b: d[b][mask_range] for b in d}
        d["bg"]   = d["Momentum"] / mass
        d["beta"] = d["Momentum"] / np.sqrt(d["Momentum"]**2 + mass**2)
        print(f"  {label}: {mask_pid.sum():,} hits totales → {mask_range.sum():,} en rango dEdx")
        return d

    mu = _make_species(0, M_MU, "mu+")
    pi = _make_species(1, M_PI, "pi+")
    return mu, pi


# ============================================================================
# Plot helpers
# ============================================================================
def _add_info(ax, particle=""):
    label = (r"Geant4  |  Bar Strip Detector  |  "
             + (particle + r"  in BC404  |  " if particle else "")
             + r"FTFP\_BERT  |  70cm Fe absorber")
    ax.text(0.01, 1.008, label, transform=ax.transAxes, fontsize=8.5,
            va="bottom", ha="left", color="#333333", style="italic")


def _save(fig, out_dir, filename):
    path = os.path.join(out_dir, filename)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {path}")


def _bethe_line(ax, bg_range, mass, color="black", lw=2, label="Landau MPV (BB)"):
    bg_th = np.logspace(np.log10(bg_range[0]), np.log10(bg_range[1]), 1000)
    mpv   = landau_mpv(bg_th, x_mm=10.0)
    v     = ~np.isnan(mpv)
    ax.plot(bg_th[v], mpv[v], color=color, lw=lw, ls="-", label=label, zorder=5)


# ============================================================================
# PLOT 1: dE/dx vs βγ
# ============================================================================
def plot_dedx_vs_bg(mu, pi, out_dir):
    fig, axes = plt.subplots(1, 2, figsize=(16, 6), sharey=True)
    fig.subplots_adjust(wspace=0.05)

    bg_range = (0.3, 100)
    norm     = mcolors.LogNorm(vmin=1, vmax=None)

    for ax, data, mass, particle_label in [
        (axes[0], mu, M_MU, r"$\mu^+$"),
        (axes[1], pi, M_PI, r"$\pi^+$"),
    ]:
        h = ax.hist2d(data["bg"], data["fdEdx"],
                      bins=[200, 200],
                      range=[list(bg_range), [DEDX_MIN, DEDX_MAX]],
                      norm=norm, cmap=CMAP)
        cb = plt.colorbar(h[3], ax=ax)
        cb.set_label("Counts", fontsize=10)

        _bethe_line(ax, bg_range, mass)

        bg_mip = 3.5 * (M_MU / mass)
        ax.axvline(bg_mip, color="gold", ls="--", lw=1.5, alpha=0.9)
        ax.text(bg_mip * 1.06, DEDX_MIN * 1.5, "MIP",
                color="goldenrod", fontsize=9, va="bottom")

        ax.set_xscale("log"); ax.set_yscale("log")
        ax.set_xlabel(r"$\beta\gamma = p\,/\,mc$", fontsize=12)
        if ax is axes[0]:
            ax.set_ylabel(r"$dE/dx$  (MeV/mm)", fontsize=12)
        ax.set_xlim(bg_range); ax.set_ylim(DEDX_MIN, DEDX_MAX)
        ax.legend(loc="upper right", fontsize=9, framealpha=0.7)
        _add_info(ax, particle=particle_label)

    fig.tight_layout()
    _save(fig, out_dir, "bethe_bloch_bg.png")


# ============================================================================
# PLOT 2: dE/dx vs β
# ============================================================================
def plot_dedx_vs_beta(mu, pi, out_dir):
    fig, axes = plt.subplots(1, 2, figsize=(16, 6), sharey=True)
    fig.subplots_adjust(wspace=0.05)

    beta_range = (0.1, 1.0)
    norm = mcolors.LogNorm(vmin=1, vmax=None)

    for ax, data, mass, particle_label in [
        (axes[0], mu, M_MU, r"$\mu^+$"),
        (axes[1], pi, M_PI, r"$\pi^+$"),
    ]:
        h = ax.hist2d(data["beta"], data["fdEdx"],
                      bins=[200, 200],
                      range=[list(beta_range), [DEDX_MIN, DEDX_MAX]],
                      norm=norm, cmap=CMAP)
        cb = plt.colorbar(h[3], ax=ax)
        cb.set_label("Counts", fontsize=10)

        beta_th = np.linspace(beta_range[0], 0.9999, 2000)
        bg_th   = beta_th / np.sqrt(1.0 - beta_th**2)
        mpv = landau_mpv(bg_th, x_mm=10.0)
        v   = ~np.isnan(mpv)
        ax.plot(beta_th[v], mpv[v], color="black", lw=2, ls="-", label="Landau MPV (BB)")

        ax.set_yscale("log")
        ax.set_xlabel(r"$\beta = v/c$", fontsize=12)
        if ax is axes[0]:
            ax.set_ylabel(r"$dE/dx$  (MeV/mm)", fontsize=12)
        ax.set_xlim(beta_range); ax.set_ylim(DEDX_MIN, DEDX_MAX)
        ax.legend(loc="upper right", fontsize=9, framealpha=0.7)
        _add_info(ax, particle=particle_label)

    fig.tight_layout()
    _save(fig, out_dir, "dedx_vs_beta.png")


# ============================================================================
# PLOT 3: dE/dx vs p (GeV/c)
# ============================================================================
def plot_dedx_vs_p(mu, pi, out_dir):
    fig, axes = plt.subplots(1, 2, figsize=(16, 6), sharey=True)
    fig.subplots_adjust(wspace=0.05)

    p_range_GeV  = (0.03, 15.0)
    DEDX_LIN_MAX = 3.0
    norm    = mcolors.LogNorm(vmin=1, vmax=None)
    p_bins  = np.logspace(np.log10(p_range_GeV[0]), np.log10(p_range_GeV[1]), 200)
    dedx_bins = np.linspace(0, DEDX_LIN_MAX, 200)

    for ax, data, mass, particle_label in [
        (axes[0], mu, M_MU, r"$\mu^+$"),
        (axes[1], pi, M_PI, r"$\pi^+$"),
    ]:
        p_GeV = data["Momentum"] / 1000.0

        h = ax.hist2d(p_GeV, data["fdEdx"],
                      bins=[p_bins, dedx_bins], norm=norm, cmap=CMAP)
        cb = plt.colorbar(h[3], ax=ax)
        cb.set_label("Counts", fontsize=10)

        p_th_GeV = np.logspace(np.log10(p_range_GeV[0]), np.log10(p_range_GeV[1]), 2000)
        bg_th    = (p_th_GeV * 1000) / mass
        mpv      = landau_mpv(bg_th, x_mm=10.0)
        v        = ~np.isnan(mpv) & (mpv < DEDX_LIN_MAX)
        ax.plot(p_th_GeV[v], mpv[v], color="black", lw=2, ls="-", label="Landau MPV (BB)")

        p_mip_GeV = 3.5 * mass / 1000.0
        ax.axvline(p_mip_GeV, color="gold", ls="--", lw=1.5, alpha=0.9)
        ax.text(p_mip_GeV * 1.1, 0.05, "MIP", color="goldenrod", fontsize=9, va="bottom")

        ax.set_xscale("log")
        ax.set_xlabel(r"$p$  (GeV/c)", fontsize=12)
        if ax is axes[0]:
            ax.set_ylabel(r"$dE/dx$  (MeV/mm)", fontsize=12)
        ax.set_xlim(p_range_GeV); ax.set_ylim(0, DEDX_LIN_MAX)
        ax.legend(loc="upper right", fontsize=9, framealpha=0.7)
        _add_info(ax, particle=particle_label)

    fig.tight_layout()
    _save(fig, out_dir, "dedx_vs_momentum.png")


# ============================================================================
# PLOT 4: Distribución de Landau
# ============================================================================
def plot_landau(mu, pi, out_dir):
    fig, ax = plt.subplots(figsize=(10, 6))

    bins = np.linspace(DEDX_MIN, DEDX_MAX, 150)

    ax.hist(mu["fdEdx"], bins=bins, histtype="stepfilled",
            color="steelblue", alpha=0.5, density=True, label=r"$\mu^+$")
    ax.hist(mu["fdEdx"], bins=bins, histtype="step",
            color="steelblue", lw=2, density=True)
    ax.hist(pi["fdEdx"], bins=bins, histtype="stepfilled",
            color="tomato", alpha=0.4, density=True, label=r"$\pi^+$")
    ax.hist(pi["fdEdx"], bins=bins, histtype="step",
            color="tomato", lw=2, density=True)

    ax.set_xlabel(r"$dE/dx$ por paso  (MeV/mm)", fontsize=12)
    ax.set_ylabel("Densidad (normalizada)", fontsize=12)
    ax.set_yscale("log")
    ax.legend(fontsize=10)
    ax.text(0.55, 0.88,
            r"Cola asimétrica: $\delta$-rays energéticos" + "\n(fluctuaciones de Landau, barras de 1 cm)",
            transform=ax.transAxes, fontsize=9, color="gray", va="top", ha="left",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
    _add_info(ax)
    fig.tight_layout()
    _save(fig, out_dir, "landau_distribution.png")


# ============================================================================
# PLOT 5: Overlay μ⁺ vs π⁺ — mediana dE/dx vs βγ
# ============================================================================
def plot_overlay(mu, pi, out_dir):
    fig, ax = plt.subplots(figsize=(10, 7))

    bg_range = (0.3, 100)
    bg_bins  = np.logspace(np.log10(bg_range[0]), np.log10(bg_range[1]), 80)
    bg_centers = 0.5 * (bg_bins[:-1] + bg_bins[1:])

    def profile_median(bg_arr, dedx_arr, bins):
        med, lo, hi = [], [], []
        for lo_e, hi_e in zip(bins[:-1], bins[1:]):
            mask = (bg_arr >= lo_e) & (bg_arr < hi_e)
            if mask.sum() < 20:
                med.append(np.nan); lo.append(np.nan); hi.append(np.nan)
                continue
            d = dedx_arr[mask]
            med.append(np.median(d))
            lo.append(np.percentile(d, 25))
            hi.append(np.percentile(d, 75))
        return np.array(med), np.array(lo), np.array(hi)

    mu_med, mu_lo, mu_hi = profile_median(mu["bg"], mu["fdEdx"], bg_bins)
    pi_med, pi_lo, pi_hi = profile_median(pi["bg"], pi["fdEdx"], bg_bins)

    v_mu = ~np.isnan(mu_med)
    v_pi = ~np.isnan(pi_med)
    ax.fill_between(bg_centers[v_mu], mu_lo[v_mu], mu_hi[v_mu], color="steelblue", alpha=0.25)
    ax.fill_between(bg_centers[v_pi], pi_lo[v_pi], pi_hi[v_pi], color="tomato",   alpha=0.25)
    ax.plot(bg_centers[v_mu], mu_med[v_mu], color="steelblue", lw=2, label=r"$\mu^+$  mediana dE/dx")
    ax.plot(bg_centers[v_pi], pi_med[v_pi], color="tomato",   lw=2, label=r"$\pi^+$  mediana dE/dx")

    bg_th = np.logspace(np.log10(bg_range[0]), np.log10(bg_range[1]), 2000)
    mpv   = landau_mpv(bg_th, x_mm=10.0, mat=AIR)
    v = ~np.isnan(mpv)
    ax.plot(bg_th[v], mpv[v], color="black", lw=2, ls="--",
            label=r"Landau MPV (BB) — universal $\mu^+\!/\pi^+$")

    ax.axvline(3.5, color="gray", ls=":", lw=1.5, alpha=0.6)
    ax.text(3.5 * 1.05, DEDX_MIN * 1.3, "MIP\n(βγ≈3.5)", color="gray", fontsize=9, va="bottom")

    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel(r"$\beta\gamma = p\,/\,mc$", fontsize=13)
    ax.set_ylabel(r"$dE/dx$  (MeV/mm)", fontsize=13)
    ax.set_xlim(bg_range); ax.set_ylim(DEDX_MIN, DEDX_MAX)
    ax.legend(fontsize=10, framealpha=0.8)
    _add_info(ax)
    fig.suptitle(r"Bethe-Bloch overlay: $\mu^+$ vs $\pi^+$ en BC404 (barras 1 cm)", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    _save(fig, out_dir, "bethe_bloch_overlay.png")


# ============================================================================
# PLOT 6: PID combinado
# ============================================================================
def plot_pid_combined(mu, pi, out_dir):
    fig, ax = plt.subplots(figsize=(10, 7))

    p_range_GeV  = (0.03, 15.0)
    DEDX_LIN_MAX = 3.0
    p_bins    = np.logspace(np.log10(p_range_GeV[0]), np.log10(p_range_GeV[1]), 200)
    dedx_bins = np.linspace(0, DEDX_LIN_MAX, 200)

    def hist2d_arrays(p_MeV, dedx):
        H, xedges, yedges = np.histogram2d(p_MeV / 1000.0, dedx, bins=[p_bins, dedx_bins])
        return H, xedges, yedges

    H_mu, xe, ye = hist2d_arrays(mu["Momentum"], mu["fdEdx"])
    H_pi, _,  _  = hist2d_arrays(pi["Momentum"], pi["fdEdx"])

    def col_norm(H):
        s = H.sum(axis=1, keepdims=True)
        s[s == 0] = 1
        return H / s

    H_mu_n = col_norm(H_mu)
    H_pi_n = col_norm(H_pi)

    Xc = 0.5 * (xe[:-1] + xe[1:])
    Yc = 0.5 * (ye[:-1] + ye[1:])
    X, Y = np.meshgrid(Xc, Yc, indexing="ij")
    vmax = max(H_mu_n.max(), H_pi_n.max()) * 0.6

    ax.pcolormesh(X, Y, H_mu_n, cmap="Blues", vmin=0, vmax=vmax, alpha=0.85, shading="auto")
    ax.pcolormesh(X, Y, H_pi_n, cmap="Reds",  vmin=0, vmax=vmax, alpha=0.65, shading="auto")

    p_th_GeV = np.logspace(np.log10(p_range_GeV[0]), np.log10(p_range_GeV[1]), 2000)
    for mass, color, label in [
        (M_MU, "royalblue", r"Landau MPV — $\mu^+$"),
        (M_PI, "firebrick", r"Landau MPV — $\pi^+$"),
    ]:
        bg_th = (p_th_GeV * 1000) / mass
        mpv   = landau_mpv(bg_th, x_mm=10.0)
        v     = ~np.isnan(mpv) & (mpv < DEDX_LIN_MAX)
        ax.plot(p_th_GeV[v], mpv[v], color=color, lw=2, ls="-", label=label, zorder=5)

    from matplotlib.patches import Patch
    legend_handles = [
        Patch(facecolor="steelblue", alpha=0.8, label=r"$\mu^+$ datos"),
        Patch(facecolor="tomato",    alpha=0.8, label=r"$\pi^+$ datos"),
        plt.Line2D([0], [0], color="royalblue", lw=2, label=r"Landau MPV — $\mu^+$"),
        plt.Line2D([0], [0], color="firebrick", lw=2, label=r"Landau MPV — $\pi^+$"),
    ]
    ax.legend(handles=legend_handles, fontsize=10, framealpha=0.85)

    p_mip_pi = 3.5 * M_PI / 1000.0
    ax.axvline(p_mip_pi, color="gold", ls="--", lw=1.2, alpha=0.8)
    ax.text(p_mip_pi * 1.08, 0.05, r"MIP ($\pi^+$)", color="goldenrod", fontsize=8, va="bottom")

    ax.set_xscale("log")
    ax.set_xlabel(r"$p$  (GeV/c)", fontsize=13)
    ax.set_ylabel(r"$dE/dx$  (MeV/mm)", fontsize=13)
    ax.set_xlim(p_range_GeV); ax.set_ylim(0, DEDX_LIN_MAX)
    _add_info(ax)
    fig.suptitle(r"PID: $\mu^+$ vs $\pi^+$ en BC404 bar strip detector", fontsize=13)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    _save(fig, out_dir, "pid_combined.png")


# ============================================================================
# PLOT 7: Distribución de hits por capa  (específico del bar strip detector)
# ============================================================================
def plot_layer_hits(mu, pi, out_dir):
    """
    Muestra cuántos hits registra cada capa del detector.
    La capa 1 (barras en X) y capa 2 (barras en Y) deben tener recuentos similares
    si la partícula pasa recto. Diferencias indican dispersión o absorción parcial.
    """
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    for ax, data, particle_label, color in [
        (axes[0], mu, r"$\mu^+$", "steelblue"),
        (axes[1], pi, r"$\pi^+$", "tomato"),
    ]:
        layer0 = (data["layerID"] == 0).sum()
        layer1 = (data["layerID"] == 1).sum()
        bars   = ax.bar(["Capa 1\n(barras X)", "Capa 2\n(barras Y)"],
                        [layer0, layer1], color=[color, color], alpha=0.7,
                        edgecolor="black")
        for rect, val in zip(bars, [layer0, layer1]):
            ax.text(rect.get_x() + rect.get_width() / 2, rect.get_height() * 1.01,
                    f"{val:,}", ha="center", va="bottom", fontsize=10)
        ax.set_ylabel("Número de hits", fontsize=11)
        ax.set_ylim(0, max(layer0, layer1) * 1.15)
        ax.text(0.5, 0.97, f"Hits por capa — {particle_label}",
                transform=ax.transAxes, fontsize=12, ha="center", va="top", fontweight="bold")
        _add_info(ax, particle=particle_label)

    fig.tight_layout()
    _save(fig, out_dir, "layer_hits.png")


# ============================================================================
# PLOT 8: Eficiencia de detección vs momento
# ============================================================================
def load_efficiency_mixed(glob_pattern: str, n_total_per_species: int = 1000):
    """
    Por cada run (archivo ROOT), cuenta cuántos eventos únicos de cada especie
    dejaron al menos un hit en el centellador.

    e = (n_total - n_detectado) / n_total   (fracción no detectada / absorbida)
    eff = n_detectado / n_total              (fracción detectada, eficiencia estándar)

    n_total_per_species: eventos disparados por especie por run (= 2000/2 = 1000).
    """
    files = sorted(glob.glob(glob_pattern))
    if not files:
        raise FileNotFoundError(f"[ERROR] No files match: {glob_pattern}")

    records = []
    for fpath in files:
        with uproot.open(fpath) as f:
            tree = _best_cycle(f)
            if tree is None:
                continue
            d = tree.arrays(["fEvent", "particleID", "Momentum"], library="np")

        for pid, mass in [(0, M_MU), (1, M_PI)]:
            mask = d["particleID"] == pid
            if mask.sum() == 0:
                n_det = 0
                p_med = np.nan
            else:
                n_det = len(np.unique(d["fEvent"][mask]))
                p_med = np.median(d["Momentum"][mask]) / 1000.0   # MeV/c → GeV/c

            eff    = min(n_det / n_total_per_species, 1.0)
            e_loss = 1.0 - eff
            records.append({"pid": pid, "p_GeV": p_med, "n_det": n_det,
                            "eff": eff, "e_loss": e_loss})

    return records


def plot_efficiency(glob_pattern: str, out_dir: str, n_total: int = 1000):
    print("\n--- Calculando eficiencia por run ---")
    records = load_efficiency_mixed(glob_pattern, n_total_per_species=n_total)

    mu_recs = sorted([r for r in records if r["pid"] == 0 and not np.isnan(r["p_GeV"])],
                     key=lambda r: r["p_GeV"])
    pi_recs = sorted([r for r in records if r["pid"] == 1 and not np.isnan(r["p_GeV"])],
                     key=lambda r: r["p_GeV"])

    mu_p   = np.array([r["p_GeV"]  for r in mu_recs])
    mu_eff = np.array([r["eff"]    for r in mu_recs])
    mu_e   = np.array([r["e_loss"] for r in mu_recs])

    pi_p   = np.array([r["p_GeV"]  for r in pi_recs])
    pi_eff = np.array([r["eff"]    for r in pi_recs])
    pi_e   = np.array([r["e_loss"] for r in pi_recs])

    fig, ax = plt.subplots(figsize=(10, 6))
    fig.suptitle(r"Eficiencia del detector vs momento  —  Bar Strip Detector (70 cm Fe + BC404)",
                 fontsize=13)

    ax.semilogx(mu_p, mu_eff, 'o-', color="steelblue", lw=2, ms=4,
                label=r"$\mu^+$")
    ax.semilogx(pi_p, pi_eff, 's-', color="tomato",    lw=2, ms=4,
                label=r"$\pi^+$")
    ax.axhline(0.5, color="gray", ls=":", lw=1, alpha=0.5)
    ax.set_xlabel(r"$p$  (GeV/c)", fontsize=12)
    ax.set_ylabel(r"$\varepsilon = N_\mathrm{det}\,/\,N_\mathrm{total}$", fontsize=12)
    ax.set_ylim(-0.05, 1.10)
    ax.legend(fontsize=11, framealpha=0.85)
    _add_info(ax)

    fig.tight_layout(rect=[0, 0, 1, 0.94])
    _save(fig, out_dir, "efficiency_vs_momentum.png")


# ============================================================================
# Main
# ============================================================================
def main(mu_path, pi_path, mixed_path, out_dir):
    os.makedirs(out_dir, exist_ok=True)

    if mixed_path:
        mu, pi = load_mixed_hits(mixed_path)
    else:
        mu = load_all_hits(mu_path, M_MU, "mu+")
        pi = load_all_hits(pi_path, M_PI, "pi+")

    plot_dedx_vs_bg(mu, pi, out_dir)
    plot_dedx_vs_beta(mu, pi, out_dir)
    plot_dedx_vs_p(mu, pi, out_dir)
    plot_landau(mu, pi, out_dir)
    plot_overlay(mu, pi, out_dir)
    plot_pid_combined(mu, pi, out_dir)
    plot_layer_hits(mu, pi, out_dir)

    if mixed_path:
        plot_efficiency(mixed_path, out_dir, n_total=1000)

    print(f"\nAll plots saved to: {out_dir}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--mixed", default=None,
                        help='Glob para ROOT files del dataset mixto (particleID col), '
                             'ej. "../Classifier/data/mixed/output_run*.root"')
    parser.add_argument("--muon", default=None,
                        help='Glob para ROOT files de muones (dataset separado)')
    parser.add_argument("--pion", default=None,
                        help='Glob para ROOT files de piones (dataset separado)')
    parser.add_argument("--out", default="img")
    args = parser.parse_args()

    if not args.mixed and not (args.muon and args.pion):
        parser.error("Usa --mixed <glob>  o bien  --muon <glob> --pion <glob>")

    main(args.muon, args.pion, args.mixed, args.out)
