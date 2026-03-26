"""
Bethe-Bloch curve plotter — mu+ vs pi+  (Geant4, iron detector — G4_Fe 7m)

Inspired by Muon_dEdx_Simulation/analysis/plot_dedx.C

Generates 4 publication-quality plots:
  1. dE/dx vs βγ          (canonical Bethe-Bloch, log-log, mu+ vs pi+ overlay)
  2. dE/dx vs β           (velocity dependence)
  3. dE/dx vs p           (momentum, log x)
  4. Landau distribution  (dE/dx per step, density comparison)

Usage:
    python plot_bethe_bloch.py \
        --muon "simulation_mu/build/output_run*.root" \
        --pion "simulation_pi/build/output_run*.root" \
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
AIR = dict(Z_over_A=0.5424, I=64.7e-6, rho=1.032)   # BC404 plastic scintillator
K   = 0.307075   # MeV·cm²/mol

# ── Plot style  (matches Muon_dEdx_Simulation/img/ reference) ────────────────
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

CMAP = "jet"   # blue→cyan→green→yellow→red, matching ROOT rainbow palette

# ── dE/dx range para BC404 (rho=1.032 g/cm³) ────────────────────────────────
# MIP en plástico ≈ 0.20 MeV/mm; rango: captura la distribución de Landau completa
DEDX_MIN = 0.01   # MeV/mm  (por debajo del MIP, excluye pasos con edep≈0)
DEDX_MAX = 5.0    # MeV/mm  (captura la cola de Landau en plástico)


# ============================================================================
# Sternheimer density-effect correction δ(βγ) for Fe (PDG parameters)
# ============================================================================
# Parámetros para G4_Fe de la tabla de Sternheimer (PDG / NIST ESTAR):
_STERN_Fe = dict(C=-3.7936, x0=0.1496, x1=2.4815, a=0.15018, m=3.4083, d0=0.00)

def density_effect(bg, stern=_STERN_Fe):
    """Corrección por efecto de densidad δ(βγ) de Sternheimer para Fe."""
    bg  = np.asarray(bg, dtype=float)
    x   = np.log10(bg)                         # x = log10(βγ)
    C, x0, x1, a, m, d0 = (stern[k] for k in ("C","x0","x1","a","m","d0"))
    delta = np.where(
        x >= x1,
        2.0 * np.log(10) * x + C,
        np.where(
            x >= x0,
            2.0 * np.log(10) * x + C + a * (x1 - x)**m,
            d0 * 10.0**(2.0 * (x - x0))        # conductores: δ₀≠0 por debajo de x0
        )
    )
    return delta


# ============================================================================
# Analytic curves  [MeV/mm] vs βγ
# ============================================================================
def bethe_bloch(bg, mass=M_MU, mat=AIR):
    """Bethe-Bloch MEAN energy loss (MeV/mm) con corrección de densidad."""
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
    Landau Most Probable Value of dE/dx for a step of x_mm (mm).

    PDG formula:
        Δp/x = ξ/x · [ln(2mₑβ²γ²/I) + ln(ξ/I) + 0.2 − β² − δ(βγ)]

    El término δ (density effect) suprime la subida relativista en materiales
    densos (hierro), produciendo el plateau de Fermi observado en los datos
    de Geant4.
    """
    bg   = np.asarray(bg, dtype=float)
    b2   = bg**2 / (1.0 + bg**2)
    bg2  = bg**2                            # = β²γ²

    x_cm = x_mm / 10.0
    xi   = 0.5 * K * mat["Z_over_A"] * mat["rho"] * x_cm / b2  # MeV
    delta = density_effect(bg)

    mpv = (xi / x_cm) * (np.log(2.0 * M_E * bg2 / mat["I"])
                         + np.log(xi / mat["I"])
                         + 0.2 - b2 - delta)   # MeV/cm
    mpv /= 10.0                                 # → MeV/mm
    return np.where(mpv > 0, mpv, np.nan)


# ============================================================================
# Data loading — per-run ROOT files, robustly reads the best cycle per file
# ============================================================================
def _best_cycle(f, tree_name="Hits"):
    """
    From an open uproot file, return the TTree cycle with the most entries.
    This handles the case where G4AnalysisManager appended empty cycles to an
    existing file (UPDATE mode on re-run), so Hits;1 has real data and later
    cycles are empty.
    """
    keys = [k for k in f.keys(cycle=True) if k.split(";")[0] == tree_name]
    if not keys:
        return None
    best_key = max(keys, key=lambda k: f[k].num_entries)
    return f[best_key] if f[best_key].num_entries > 0 else None


def load_all_hits(glob_pattern: str, mass: float, label: str):
    """
    Load hits from all per-run ROOT files matching *glob_pattern*.
    From each file the cycle with the most entries is used (robustly handles
    stale empty cycles left by G4AnalysisManager UPDATE-mode re-runs).

    Returns dict with arrays: fdEdx, Momentum, Ekin, bg, beta.
    """
    files = sorted(glob.glob(glob_pattern))
    if not files:
        raise FileNotFoundError(f"[ERROR] No files match: {glob_pattern}")

    print(f"\n--- {label} : {len(files)} file(s) matching '{glob_pattern}' ---")
    branches = ["fdEdx", "Momentum", "Ekin"]

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
        raise RuntimeError(f"[ERROR] No valid entries found for '{glob_pattern}'. Rerun the simulation.")

    data = {b: np.concatenate(parts[b]) for b in branches}

    # Keep only hits with valid dEdx in display range
    mask = (data["fdEdx"] >= DEDX_MIN) & (data["fdEdx"] <= DEDX_MAX) & (data["Momentum"] > 0)
    data = {b: data[b][mask] for b in branches}

    data["bg"]   = data["Momentum"] / mass
    data["beta"] = data["Momentum"] / np.sqrt(data["Momentum"]**2 + mass**2)

    print(f"  Hits in dEdx range [{DEDX_MIN:.0e}, {DEDX_MAX:.0e}]: {mask.sum():,}")
    return data


# ============================================================================
# Plot helpers
# ============================================================================
def _add_info(ax, particle=""):
    """Top-left italic label — mirrors the ROOT title in the reference images."""
    label = r"Geant4  |  " + (particle + r"  in BC404  |  " if particle else "") + r"FTFP\_BERT  |  $10^3$ events/beam"
    ax.text(0.01, 1.008, label,
            transform=ax.transAxes, fontsize=8.5,
            va="bottom", ha="left", color="#333333", style="italic")


def _save(fig, out_dir, filename):
    path = os.path.join(out_dir, filename)
    fig.savefig(path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {path}")


def _bethe_line(ax, bg_range, mass, color="black", lw=2,
                label="Landau MPV (BB)"):
    """Overlay the Landau MPV analytic curve (matches the data mode)."""
    bg_th = np.logspace(np.log10(bg_range[0]), np.log10(bg_range[1]), 1000)
    mpv   = landau_mpv(bg_th, x_mm=10.0)   # 10-mm cell, normal incidence
    v     = ~np.isnan(mpv)
    ax.plot(bg_th[v], mpv[v], color=color, lw=lw, ls="-", label=label, zorder=5)


# ============================================================================
# PLOT 1: dE/dx vs βγ   (side-by-side, jet colormap, reference style)
# ============================================================================
def plot_dedx_vs_bg(mu, pi, out_dir):
    fig, axes = plt.subplots(1, 2, figsize=(16, 6), sharey=True)
    fig.subplots_adjust(wspace=0.05)

    bg_range = (0.3, 100)  # cubre rising edge (βγ<1), MIP (~3.5) y subida relativista (>10)
    norm     = mcolors.LogNorm(vmin=1, vmax=None)

    for ax, data, mass, label, particle_label in [
        (axes[0], mu, M_MU, r"$\mu^+$",  r"$\mu^+$"),
        (axes[1], pi, M_PI, r"$\pi^+$",  r"$\pi^+$"),
    ]:
        h = ax.hist2d(
            data["bg"], data["fdEdx"],
            bins=[200, 200],
            range=[list(bg_range), [DEDX_MIN, DEDX_MAX]],
            norm=norm,
            cmap=CMAP,
        )
        cb = plt.colorbar(h[3], ax=ax)
        cb.set_label("Counts", fontsize=10)

        _bethe_line(ax, bg_range, mass, color="black", lw=2, label="Landau MPV (BB)")

        # MIP marker (dashed yellow, like reference)
        bg_mip = 3.5 * (M_MU / mass)
        ax.axvline(bg_mip, color="gold", ls="--", lw=1.5, alpha=0.9)
        ax.text(bg_mip * 1.06, DEDX_MIN * 1.5, "MIP",
                color="goldenrod", fontsize=9, va="bottom")

        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel(r"$\beta\gamma = p\,/\,mc$", fontsize=12)
        if ax is axes[0]:
            ax.set_ylabel(r"$dE/dx$  (MeV/mm)", fontsize=12)
        ax.set_xlim(bg_range)
        ax.set_ylim(DEDX_MIN, DEDX_MAX)
        ax.legend(loc="upper right", fontsize=9, framealpha=0.7)
        _add_info(ax, particle=particle_label)

    fig.tight_layout()
    _save(fig, out_dir, "bethe_bloch_bg.png")


# ============================================================================
# PLOT 2: dE/dx vs β  (side-by-side, reference style)
# ============================================================================
def plot_dedx_vs_beta(mu, pi, out_dir):
    fig, axes = plt.subplots(1, 2, figsize=(16, 6), sharey=True)
    fig.subplots_adjust(wspace=0.05)

    beta_range = (0.1, 1.0)  # βγ=0.3 → β≈0.287; ampliar para ver el rising edge
    norm = mcolors.LogNorm(vmin=1, vmax=None)

    for ax, data, mass, particle_label in [
        (axes[0], mu, M_MU, r"$\mu^+$"),
        (axes[1], pi, M_PI, r"$\pi^+$"),
    ]:
        h = ax.hist2d(
            data["beta"], data["fdEdx"],
            bins=[200, 200],
            range=[list(beta_range), [DEDX_MIN, DEDX_MAX]],
            norm=norm, cmap=CMAP,
        )
        cb = plt.colorbar(h[3], ax=ax)
        cb.set_label("Counts", fontsize=10)

        # Landau MPV overlay
        beta_th = np.linspace(beta_range[0], 0.9999, 2000)
        bg_th   = beta_th / np.sqrt(1.0 - beta_th**2)
        mpv = landau_mpv(bg_th, x_mm=10.0)
        v   = ~np.isnan(mpv)
        ax.plot(beta_th[v], mpv[v], color="black", lw=2, ls="-", label="Landau MPV (BB)")

        # MIP line
        beta_mip = 3.5 / np.sqrt(1 + 3.5**2)
        ax.axvline(beta_mip * (M_MU / mass) / np.sqrt(1 + (3.5 * M_MU / mass)**2) *
                   np.sqrt(1 + (3.5 * M_MU / mass)**2),
                   color="gold", ls="--", lw=1.5, alpha=0.9)

        ax.set_yscale("log")
        ax.set_xlabel(r"$\beta = v/c$", fontsize=12)
        if ax is axes[0]:
            ax.set_ylabel(r"$dE/dx$  (MeV/mm)", fontsize=12)
        ax.set_xlim(beta_range)
        ax.set_ylim(DEDX_MIN, DEDX_MAX)
        ax.legend(loc="upper right", fontsize=9, framealpha=0.7)
        _add_info(ax, particle=particle_label)

    fig.tight_layout()
    _save(fig, out_dir, "dedx_vs_beta.png")


# ============================================================================
# PLOT 3: dE/dx vs momentum  (side-by-side, reference style)
# ============================================================================
def plot_dedx_vs_p(mu, pi, out_dir):
    """PID-style plot: dE/dx vs momentum (GeV/c), linear y-axis — estilo ALICE/LHCb."""
    fig, axes = plt.subplots(1, 2, figsize=(16, 6), sharey=True)
    fig.subplots_adjust(wspace=0.05)

    # Momentum en GeV/c (columna Momentum está en MeV/c → dividir entre 1000)
    p_range_GeV = (0.03, 15.0)   # GeV/c
    DEDX_LIN_MAX = 3.0            # MeV/mm — rango lineal para plástico
    norm    = mcolors.LogNorm(vmin=1, vmax=None)
    p_bins  = np.logspace(np.log10(p_range_GeV[0]), np.log10(p_range_GeV[1]), 200)
    dedx_bins = np.linspace(0, DEDX_LIN_MAX, 200)

    for ax, data, mass, particle_label in [
        (axes[0], mu, M_MU, r"$\mu^+$"),
        (axes[1], pi, M_PI, r"$\pi^+$"),
    ]:
        p_GeV = data["Momentum"] / 1000.0   # MeV/c → GeV/c

        h = ax.hist2d(
            p_GeV, data["fdEdx"],
            bins=[p_bins, dedx_bins],
            norm=norm, cmap=CMAP,
        )
        cb = plt.colorbar(h[3], ax=ax)
        cb.set_label("Counts", fontsize=10)

        # Landau MPV overlay
        p_th_GeV = np.logspace(np.log10(p_range_GeV[0]), np.log10(p_range_GeV[1]), 2000)
        bg_th    = (p_th_GeV * 1000) / mass   # convertir a MeV/c para βγ
        mpv      = landau_mpv(bg_th, x_mm=10.0)
        v        = ~np.isnan(mpv) & (mpv < DEDX_LIN_MAX)
        ax.plot(p_th_GeV[v], mpv[v], color="black", lw=2, ls="-", label="Landau MPV (BB)")

        # MIP line at βγ≈3.5
        p_mip_GeV = 3.5 * mass / 1000.0
        ax.axvline(p_mip_GeV, color="gold", ls="--", lw=1.5, alpha=0.9)
        ax.text(p_mip_GeV * 1.1, 0.05, "MIP", color="goldenrod", fontsize=9, va="bottom")

        ax.set_xscale("log")
        # eje Y lineal — igual que plots de PID experimentales (ALICE, LHCb, etc.)
        ax.set_xlabel(r"$p$  (GeV/c)", fontsize=12)
        if ax is axes[0]:
            ax.set_ylabel(r"$dE/dx$  (MeV/mm)", fontsize=12)
        ax.set_xlim(p_range_GeV)
        ax.set_ylim(0, DEDX_LIN_MAX)
        ax.legend(loc="upper right", fontsize=9, framealpha=0.7)
        _add_info(ax, particle=particle_label)

    fig.tight_layout()
    _save(fig, out_dir, "dedx_vs_momentum.png")


# ============================================================================
# PLOT 4: Landau distribution
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

    ax.set_xlabel(r"$dE/dx$ per step  (MeV/mm)", fontsize=12)
    ax.set_ylabel("Density (normalized)", fontsize=12)
    ax.set_yscale("log")
    ax.legend(fontsize=10)
    ax.text(0.55, 0.88,
            r"Asymmetric tail: high-energy $\delta$-rays" + "\n(Landau fluctuations in thin layers)",
            transform=ax.transAxes, fontsize=9, color="gray",
            va="top", ha="left",
            bbox=dict(boxstyle="round,pad=0.3", facecolor="white", alpha=0.8))
    _add_info(ax)
    fig.tight_layout()
    _save(fig, out_dir, "landau_distribution.png")


# ============================================================================
# PLOT 5: Overlay μ⁺ vs π⁺  —  mediana dE/dx por bin de βγ + curvas analíticas
# ============================================================================
def plot_overlay(mu, pi, out_dir):
    fig, ax = plt.subplots(figsize=(10, 7))

    bg_range = (0.3, 100)
    bg_bins  = np.logspace(np.log10(bg_range[0]), np.log10(bg_range[1]), 80)
    bg_centers = 0.5 * (bg_bins[:-1] + bg_bins[1:])

    def profile_median(bg_arr, dedx_arr, bins):
        """Mediana de dE/dx en cada bin de βγ (mínimo 20 hits para fiabilidad)."""
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

    # Banda IQR (25–75 percentil)
    v_mu = ~np.isnan(mu_med)
    v_pi = ~np.isnan(pi_med)
    ax.fill_between(bg_centers[v_mu], mu_lo[v_mu], mu_hi[v_mu],
                    color="steelblue", alpha=0.25, label=None)
    ax.fill_between(bg_centers[v_pi], pi_lo[v_pi], pi_hi[v_pi],
                    color="tomato",   alpha=0.25, label=None)

    # Mediana de datos
    ax.plot(bg_centers[v_mu], mu_med[v_mu], color="steelblue", lw=2,
            label=r"$\mu^+$  mediana dE/dx")
    ax.plot(bg_centers[v_pi], pi_med[v_pi], color="tomato",   lw=2,
            label=r"$\pi^+$  mediana dE/dx")

    # Curva analítica Landau MPV — universal en βγ para partículas cargadas pesadas
    # (μ y π ambas tienen M >> mₑ, por lo que Tmax ≈ 2mₑβ²γ² y la curva es idéntica)
    bg_th = np.logspace(np.log10(bg_range[0]), np.log10(bg_range[1]), 2000)
    mpv   = landau_mpv(bg_th, x_mm=10.0, mat=AIR)
    v = ~np.isnan(mpv)
    ax.plot(bg_th[v], mpv[v], color="black", lw=2, ls="--",
            label=r"Landau MPV (BB) — universal $\mu^+\!/\pi^+$")

    ax.axvline(3.5, color="gray", ls=":", lw=1.5, alpha=0.6)
    ax.text(3.5 * 1.05, DEDX_MIN * 1.3, "MIP\n(βγ≈3.5)",
            color="gray", fontsize=9, va="bottom")

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(r"$\beta\gamma = p\,/\,mc$", fontsize=13)
    ax.set_ylabel(r"$dE/dx$  (MeV/mm)", fontsize=13)
    ax.set_title(r"Bethe-Bloch overlay: $\mu^+$ vs $\pi^+$ en BC404 scintillator", fontsize=13)
    ax.set_xlim(bg_range)
    ax.set_ylim(DEDX_MIN, DEDX_MAX)
    ax.legend(fontsize=10, framealpha=0.8)
    _add_info(ax)

    fig.tight_layout()
    _save(fig, out_dir, "bethe_bloch_overlay.png")


# ============================================================================
# PLOT 6: Combined PID plot — μ⁺ y π⁺ en un solo panel (p GeV/c, Y lineal)
# ============================================================================
def plot_pid_combined(mu, pi, out_dir):
    """
    Histograma 2D de μ⁺ (azul) y π⁺ (rojo) superpuestos en un solo panel.
    Eje X: momentum en GeV/c (log).  Eje Y: dE/dx en MeV/mm (lineal).
    Estilo PID experimental — muestra la separación entre species.
    """
    fig, ax = plt.subplots(figsize=(10, 7))

    p_range_GeV  = (0.03, 15.0)
    DEDX_LIN_MAX = 3.0
    p_bins    = np.logspace(np.log10(p_range_GeV[0]), np.log10(p_range_GeV[1]), 200)
    dedx_bins = np.linspace(0, DEDX_LIN_MAX, 200)

    # Calcular histogramas 2D para cada partícula
    def hist2d_arrays(p_MeV, dedx):
        H, xedges, yedges = np.histogram2d(
            p_MeV / 1000.0, dedx,
            bins=[p_bins, dedx_bins]
        )
        return H, xedges, yedges

    H_mu, xe, ye = hist2d_arrays(mu["Momentum"], mu["fdEdx"])
    H_pi, _,  _  = hist2d_arrays(pi["Momentum"], pi["fdEdx"])

    # Normalizar por columna (perfil de densidad) para que ambas species sean visibles
    # aunque tengan distinto número de hits
    def col_norm(H):
        col_sum = H.sum(axis=1, keepdims=True)
        col_sum[col_sum == 0] = 1
        return H / col_sum

    H_mu_n = col_norm(H_mu)
    H_pi_n = col_norm(H_pi)

    # Mostrar como pcolormesh con colormaps distintos y alpha
    Xc = 0.5 * (xe[:-1] + xe[1:])
    Yc = 0.5 * (ye[:-1] + ye[1:])
    X, Y = np.meshgrid(Xc, Yc, indexing="ij")

    vmax = max(H_mu_n.max(), H_pi_n.max()) * 0.6   # saturar un poco para contraste

    ax.pcolormesh(X, Y, H_mu_n,
                  cmap="Blues", vmin=0, vmax=vmax, alpha=0.85, shading="auto")
    ax.pcolormesh(X, Y, H_pi_n,
                  cmap="Reds",  vmin=0, vmax=vmax, alpha=0.65, shading="auto")

    # Curvas de Bethe-Bloch para cada partícula
    p_th_GeV = np.logspace(np.log10(p_range_GeV[0]), np.log10(p_range_GeV[1]), 2000)
    for mass, color, label in [
        (M_MU, "royalblue", r"Landau MPV — $\mu^+$"),
        (M_PI, "firebrick", r"Landau MPV — $\pi^+$"),
    ]:
        bg_th = (p_th_GeV * 1000) / mass
        mpv   = landau_mpv(bg_th, x_mm=10.0)
        v     = ~np.isnan(mpv) & (mpv < DEDX_LIN_MAX)
        ax.plot(p_th_GeV[v], mpv[v], color=color, lw=2, ls="-", label=label, zorder=5)

    # Parches para la leyenda del histograma
    from matplotlib.patches import Patch
    legend_handles = [
        Patch(facecolor="steelblue", alpha=0.8, label=r"$\mu^+$ datos"),
        Patch(facecolor="tomato",    alpha=0.8, label=r"$\pi^+$ datos"),
        plt.Line2D([0], [0], color="royalblue", lw=2, label=r"Landau MPV — $\mu^+$"),
        plt.Line2D([0], [0], color="firebrick", lw=2, label=r"Landau MPV — $\pi^+$"),
    ]
    ax.legend(handles=legend_handles, fontsize=10, framealpha=0.85)

    # MIP marker (pion, el más pesado → MIP a mayor p)
    p_mip_pi = 3.5 * M_PI / 1000.0
    ax.axvline(p_mip_pi, color="gold", ls="--", lw=1.2, alpha=0.8)
    ax.text(p_mip_pi * 1.08, 0.05, r"MIP ($\pi^+$)",
            color="goldenrod", fontsize=8, va="bottom")

    ax.set_xscale("log")
    ax.set_xlabel(r"$p$  (GeV/c)", fontsize=13)
    ax.set_ylabel(r"$dE/dx$  (MeV/mm)", fontsize=13)
    ax.set_title(r"PID: $\mu^+$ vs $\pi^+$ en BC404 scintillator", fontsize=13)
    ax.set_xlim(p_range_GeV)
    ax.set_ylim(0, DEDX_LIN_MAX)
    _add_info(ax)

    fig.tight_layout()
    _save(fig, out_dir, "pid_combined.png")


# ============================================================================
# Main
# ============================================================================
def main(mu_path, pi_path, out_dir):
    os.makedirs(out_dir, exist_ok=True)

    mu = load_all_hits(mu_path, M_MU, "mu+")
    pi = load_all_hits(pi_path, M_PI, "pi+")

    plot_dedx_vs_bg(mu, pi, out_dir)
    plot_dedx_vs_beta(mu, pi, out_dir)
    plot_dedx_vs_p(mu, pi, out_dir)
    plot_landau(mu, pi, out_dir)
    plot_overlay(mu, pi, out_dir)
    plot_pid_combined(mu, pi, out_dir)

    print(f"\nAll plots saved to: {out_dir}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--muon", required=True,
                        help='Glob pattern for muon ROOT files, e.g. "simulation_mu/build/output_run*.root"')
    parser.add_argument("--pion", required=True,
                        help='Glob pattern for pion ROOT files, e.g. "simulation_pi/build/output_run*.root"')
    parser.add_argument("--out",  default="img")
    args = parser.parse_args()
    main(args.muon, args.pion, args.out)
