"""
================================================================================
plot_all.py — Bar Strip Detector (unified plotter)
================================================================================
Genera 10 imágenes de un haz mixto mu+/pi+ atravesando 70 cm de Fe +
dos capas de barras centelladoras BC404.

Uso:
    python plot_all.py --mixed "../Classifier/data/mixed/output_run*.root"
    python plot_all.py --mixed "ruta/output_run*.root" --out img/

Dependencias:
    pip install numpy matplotlib uproot

Imágenes generadas:
  Bethe-Bloch / PID:
    1. bethe_bloch_bg.png         — dE/dx vs beta-gamma (hist 2D, lado a lado)
    2. dedx_vs_beta.png           — dE/dx vs beta (velocidad)
    3. dedx_vs_momentum.png       — dE/dx vs p (GeV/c, estilo PID)
    4. bethe_bloch_overlay.png    — Overlay mediana dE/dx mu+ vs pi+
    5. pid_combined.png           — Histograma 2D combinado mu+/pi+
    6. detector_layout.png        — Diagrama esquematico del detector
    7. layer_hits.png             — Hits por capa (Capa 1 vs Capa 2)

  Landau:
    8. landau_corregida.png       — Landau run 45 con tabla de estadisticas

  Eficiencia:
    9. eff_momento_corregida.png  — Eficiencia vs momento (barras error, regimenes)
   10. eff_angulo_corregida.png   — Eficiencia vs angulo ConeAngle (barras error)
================================================================================
"""

import argparse
import glob
import os
import re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import uproot

# ============================================================================
# Constantes fisicas
# ============================================================================
M_MU  = 105.6583755   # MeV/c^2
M_PI  = 139.5702      # MeV/c^2
M_E   = 0.51099895    # MeV/c^2

# BC404 plastic scintillator (G4_PLASTIC_SC_VINYLTOLUENE)
# C9H10: rho=1.032 g/cm^3, I=64.7 eV, Z/A=0.5424
AIR = dict(Z_over_A=0.5424, I=64.7e-6, rho=1.032)
K   = 0.307075        # MeV*cm^2/mol

DEDX_MIN = 0.01       # MeV/mm
DEDX_MAX = 5.0        # MeV/mm
CMAP = "viridis"

# Color palette — tableau-inspired, colorblind-safe
COLOR_MU = "#4C72B0"   # blue for mu+
COLOR_PI = "#DD8452"   # warm orange for pi+
COLOR_MU_LIGHT = "#4C72B0"
COLOR_PI_LIGHT = "#DD8452"

# Sternheimer density-effect correction for BC404
_STERN_BC404 = dict(C=-3.7936, x0=0.1496, x1=2.4815, a=0.15018, m=3.4083, d0=0.00)

# HEP paper-style rcParams (serif, CM math, inward ticks, subtle grid)
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["DejaVu Serif", "Times New Roman", "Liberation Serif"],
    "mathtext.fontset": "cm",
    "font.size": 11,
    "axes.titlesize": 14,
    "axes.titleweight": "bold",
    "axes.labelsize": 12,
    "xtick.labelsize": 11,
    "ytick.labelsize": 11,
    "legend.fontsize": 10,
    "figure.dpi": 200,
    "savefig.dpi": 200,
    "axes.facecolor": "white",
    "figure.facecolor": "white",
    "axes.grid": True,
    "grid.alpha": 0.25,
    "grid.linestyle": "--",
    "grid.linewidth": 0.7,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "xtick.minor.visible": True,
    "ytick.minor.visible": True,
})

INFO_STR = r"Geant4  |  Bar Strip Detector  |  70 cm Fe + BC404  |  FTFP\_BERT  |  $B_x = 1.0\,\mathrm{T}$"

# ============================================================================
# Fisica: Bethe-Bloch, Landau MPV, Sternheimer
# ============================================================================
def density_effect(bg, stern=_STERN_BC404):
    bg = np.asarray(bg, dtype=float)
    x = np.log10(bg)
    C, x0, x1, a, m, d0 = (stern[k] for k in ("C", "x0", "x1", "a", "m", "d0"))
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


def bethe_bloch(bg, mass=M_MU, mat=AIR):
    bg = np.asarray(bg, dtype=float)
    b2 = bg**2 / (1.0 + bg**2)
    g = np.sqrt(1.0 + bg**2)
    Tmax = (2.0 * M_E * bg**2) / (1.0 + 2.0 * g * M_E / mass + (M_E / mass)**2)
    logA = np.where(Tmax > 0, 2.0 * M_E * b2 * g**2 * Tmax / mat["I"]**2, np.nan)
    delta = density_effect(bg)
    dedx = K * mat["Z_over_A"] / b2 * (0.5 * np.log(logA) - b2 - delta / 2.0)
    dedx = np.where(dedx > 0, dedx, np.nan)
    return dedx * mat["rho"] / 10.0


def landau_mpv(bg, x_mm=10.0, mat=AIR):
    bg = np.asarray(bg, dtype=float)
    b2 = bg**2 / (1.0 + bg**2)
    bg2 = bg**2
    x_cm = x_mm / 10.0
    xi = 0.5 * K * mat["Z_over_A"] * mat["rho"] * x_cm / b2
    delta = density_effect(bg)
    mpv = (xi / x_cm) * (np.log(2.0 * M_E * bg2 / mat["I"])
                         + np.log(xi / mat["I"])
                         + 0.2 - b2 - delta)
    mpv /= 10.0
    return np.where(mpv > 0, mpv, np.nan)


# ============================================================================
# Helpers
# ============================================================================
def _best_cycle(f, tree_name="Hits"):
    keys = [k for k in f.keys(cycle=True) if k.split(";")[0] == tree_name]
    if not keys:
        return None
    best_key = max(keys, key=lambda k: f[k].num_entries)
    return f[best_key] if f[best_key].num_entries > 0 else None


def run_num(f):
    m = re.search(r"output_run(\d+)", f)
    return int(m.group(1)) if m else 0


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
    mpv = landau_mpv(bg_th, x_mm=10.0)
    v = ~np.isnan(mpv)
    ax.plot(bg_th[v], mpv[v], color=color, lw=lw, ls="-", label=label, zorder=5)


def _setup_ax(ax, log_x=False, log_y=False):
    """Apply consistent HEP-style formatting to an axis."""
    ax.minorticks_on()
    ax.tick_params(which="both", direction="in", top=True, right=True)
    if log_x:
        ax.set_xscale("log")
    if log_y:
        ax.set_yscale("log")


def _legend_kwargs(**overrides):
    """Default legend kwargs with HEP styling."""
    kw = dict(framealpha=0.9, edgecolor="lightgray", fancybox=True)
    kw.update(overrides)
    return kw


def _info_box_kwargs(**overrides):
    """Default info/annotation box kwargs."""
    kw = dict(boxstyle="round,pad=0.4", facecolor="white", edgecolor="#CCCCCC", alpha=0.8)
    kw.update(overrides)
    return kw


# ============================================================================
# Data loading
# ============================================================================
def load_mixed_hits(glob_pattern: str):
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
        d["bg"] = d["Momentum"] / mass
        d["beta"] = d["Momentum"] / np.sqrt(d["Momentum"]**2 + mass**2)
        print(f"  {label}: {mask_pid.sum():,} hits totales -> {mask_range.sum():,} en rango dEdx")
        return d

    mu = _make_species(0, M_MU, "mu+")
    pi = _make_species(1, M_PI, "pi+")
    return mu, pi


# ============================================================================
# PLOT 1: dE/dx vs beta-gamma
# ============================================================================
def plot_dedx_vs_bg(mu, pi, out_dir):
    fig, axes = plt.subplots(1, 2, figsize=(16, 6), sharey=True)
    fig.subplots_adjust(wspace=0.05)

    bg_range = (0.3, 100)
    norm = mcolors.LogNorm(vmin=1, vmax=None)

    for ax, data, mass, particle_label in [
        (axes[0], mu, M_MU, r"$\mu^+$"),
        (axes[1], pi, M_PI, r"$\pi^+$"),
    ]:
        h = ax.hist2d(data["bg"], data["fdEdx"],
                      bins=[200, 200],
                      range=[list(bg_range), [DEDX_MIN, DEDX_MAX]],
                      norm=norm, cmap=CMAP)
        cb = plt.colorbar(h[3], ax=ax, shrink=0.82, aspect=14)
        cb.set_label("Counts", fontsize=10)
        cb.ax.minorticks_on()
        cb.ax.tick_params(direction="in")

        _bethe_line(ax, bg_range, mass)

        bg_mip = 3.5 * (M_MU / mass)
        ax.axvline(bg_mip, color="gold", ls="--", lw=1.5, alpha=0.9)
        ax.text(bg_mip * 1.06, DEDX_MIN * 1.5, "MIP",
                color="goldenrod", fontsize=9, va="bottom")

        _setup_ax(ax, log_x=True, log_y=True)
        ax.set_xlabel(r"$\beta\gamma = p\,/\,mc$", fontsize=12)
        if ax is axes[0]:
            ax.set_ylabel(r"$dE/dx$  (MeV/mm)", fontsize=12)
        ax.set_xlim(bg_range); ax.set_ylim(DEDX_MIN, DEDX_MAX)
        ax.legend(loc="upper right", fontsize=9, **_legend_kwargs())
        _add_info(ax, particle=particle_label)

    fig.tight_layout()
    _save(fig, out_dir, "bethe_bloch_bg.png")


# ============================================================================
# PLOT 2: dE/dx vs beta
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
        cb = plt.colorbar(h[3], ax=ax, shrink=0.82, aspect=14)
        cb.set_label("Counts", fontsize=10)
        cb.ax.minorticks_on()
        cb.ax.tick_params(direction="in")

        beta_th = np.linspace(beta_range[0], 0.9999, 2000)
        bg_th = beta_th / np.sqrt(1.0 - beta_th**2)
        mpv = landau_mpv(bg_th, x_mm=10.0)
        v = ~np.isnan(mpv)
        ax.plot(beta_th[v], mpv[v], color="black", lw=2, ls="-", label="Landau MPV (BB)")

        _setup_ax(ax, log_y=True)
        ax.set_xlabel(r"$\beta = v/c$", fontsize=12)
        if ax is axes[0]:
            ax.set_ylabel(r"$dE/dx$  (MeV/mm)", fontsize=12)
        ax.set_xlim(beta_range); ax.set_ylim(DEDX_MIN, DEDX_MAX)
        ax.legend(loc="upper right", fontsize=9, **_legend_kwargs())
        _add_info(ax, particle=particle_label)

    fig.tight_layout()
    _save(fig, out_dir, "dedx_vs_beta.png")


# ============================================================================
# PLOT 3: dE/dx vs p (GeV/c)
# ============================================================================
def plot_dedx_vs_p(mu, pi, out_dir):
    fig, axes = plt.subplots(1, 2, figsize=(16, 6), sharey=True)
    fig.subplots_adjust(wspace=0.05)

    p_range_GeV = (0.03, 15.0)
    DEDX_LIN_MAX = 3.0
    norm = mcolors.LogNorm(vmin=1, vmax=None)
    p_bins = np.logspace(np.log10(p_range_GeV[0]), np.log10(p_range_GeV[1]), 200)
    dedx_bins = np.linspace(0, DEDX_LIN_MAX, 200)

    for ax, data, mass, particle_label in [
        (axes[0], mu, M_MU, r"$\mu^+$"),
        (axes[1], pi, M_PI, r"$\pi^+$"),
    ]:
        p_GeV = data["Momentum"] / 1000.0

        h = ax.hist2d(p_GeV, data["fdEdx"],
                      bins=[p_bins, dedx_bins], norm=norm, cmap=CMAP)
        cb = plt.colorbar(h[3], ax=ax, shrink=0.82, aspect=14)
        cb.set_label("Counts", fontsize=10)
        cb.ax.minorticks_on()
        cb.ax.tick_params(direction="in")

        p_th_GeV = np.logspace(np.log10(p_range_GeV[0]), np.log10(p_range_GeV[1]), 2000)
        bg_th = (p_th_GeV * 1000) / mass
        mpv = landau_mpv(bg_th, x_mm=10.0)
        v = ~np.isnan(mpv) & (mpv < DEDX_LIN_MAX)
        ax.plot(p_th_GeV[v], mpv[v], color="black", lw=2, ls="-", label="Landau MPV (BB)")

        p_mip_GeV = 3.5 * mass / 1000.0
        ax.axvline(p_mip_GeV, color="gold", ls="--", lw=1.5, alpha=0.9)
        ax.text(p_mip_GeV * 1.1, 0.05, "MIP", color="goldenrod", fontsize=9, va="bottom")

        _setup_ax(ax, log_x=True)
        ax.set_xlabel(r"$p$  (GeV/c)", fontsize=12)
        if ax is axes[0]:
            ax.set_ylabel(r"$dE/dx$  (MeV/mm)", fontsize=12)
        ax.set_xlim(p_range_GeV); ax.set_ylim(0, DEDX_LIN_MAX)
        ax.legend(loc="upper right", fontsize=9, **_legend_kwargs())
        _add_info(ax, particle=particle_label)

    fig.tight_layout()
    _save(fig, out_dir, "dedx_vs_momentum.png")


# ============================================================================
# PLOT 4: Detector layout (schematic)
# ============================================================================
def plot_detector_layout(out_dir, seed=42):
    from matplotlib.patches import Rectangle, Patch
    from matplotlib.lines import Line2D

    rng = np.random.default_rng(seed)

    z_src = -200.0
    z_fe0 = 0.0
    z_fe1 = 70.0
    x_fe = 35.0
    z_L1 = 100.5
    z_L2 = 103.5
    x_bar = 50.0
    x_bar_c = 47.5
    dist_L1 = z_L1 - z_src
    dist_L2 = z_L2 - z_src
    footprint_L1 = round(x_fe * dist_L1 / (-z_src), 1)
    footprint_L2 = round(x_fe * dist_L2 / (-z_src), 1)
    tx_acc = x_bar * (-z_src) / dist_L1
    theta_acc = np.degrees(np.arctan(tx_acc / (-z_src)))
    tx_lat = x_fe / (1.0 + z_fe1 / (-z_src))
    theta_lat = np.degrees(np.arctan(tx_lat / (-z_src)))

    fig = plt.figure(figsize=(21, 7))
    fig.suptitle(
        r"Bar Strip Detector + $B_x = 1.0\ \mathrm{T}$ — Configuración del sistema" + "\n"
        r"Absorbedor: G4\_Fe  70$\times$70$\times$70 cm  |  "
        r"Centellador: BC404  |  "
        r"Campo magnético: $B_x = 1.0\ \mathrm{T}$ en región vacío (fuente → Fe)",
        fontsize=12, fontweight="bold", y=0.99, va="top")

    gs = fig.add_gridspec(1, 3, width_ratios=[3.2, 2, 2],
                          wspace=0.22, left=0.04, right=0.98,
                          top=0.87, bottom=0.12)
    ax_s = fig.add_subplot(gs[0])
    ax_1 = fig.add_subplot(gs[1])
    ax_2 = fig.add_subplot(gs[2])

    ax = ax_s
    ax.set_title("Vista lateral  (plano Z – X,  corte en Y = 0)", fontsize=10, pad=4)
    ax.set_xlabel("z (cm)", fontsize=10)
    ax.set_ylabel("x (cm)", fontsize=10)

    x_at_L1 = x_fe * dist_L1 / (-z_src)
    ax.fill(
        [z_src, z_fe0, z_L1 + 8, z_L1 + 8, z_fe0, z_src],
        [0, x_fe, x_at_L1 + 2, -(x_at_L1 + 2), -x_fe, 0],
        color="orange", alpha=0.13, zorder=0)

    ax.add_patch(Rectangle((z_fe0, -x_fe), z_fe1 - z_fe0, 2 * x_fe,
                            facecolor="#C8860A", edgecolor="black",
                            lw=1.5, alpha=0.72, zorder=1))
    ax.text((z_fe0 + z_fe1) / 2, 0,
            r"G4\_Fe" + "\n70×70×70 cm\n= 4.17 λ_I",
            ha="center", va="center", fontsize=9,
            fontweight="bold", color="white", zorder=2)

    ax.add_patch(Rectangle((z_L1 - 0.5, -(x_bar_c + 3)), 1.0, 2 * (x_bar_c + 3),
                            facecolor="royalblue", edgecolor="royalblue",
                            alpha=0.80, zorder=2))
    ax.text(z_L1, x_bar_c + 5,
            "Capa 1\n(barras X)", ha="center", va="bottom",
            fontsize=8.5, color="royalblue", fontweight="bold")

    ax.add_patch(Rectangle((z_L2 - 0.5, -(x_bar_c + 3)), 1.0, 2 * (x_bar_c + 3),
                            facecolor="firebrick", edgecolor="firebrick",
                            alpha=0.80, zorder=2))
    ax.text(z_L2, -(x_bar_c + 5),
            "Capa 2\n(barras Y)", ha="center", va="top",
            fontsize=8.5, color="firebrick", fontweight="bold")

    # ── Región del campo magnético Bz = 0.5 T (z = -200 cm a z = 0) ──
    ax.fill([z_src, z_fe0, z_fe0, z_src],
            [-x_fe - 4, -x_fe - 4, x_fe + 4, x_fe + 4],
            color="deepskyblue", alpha=0.08, zorder=0)
    ax.text((z_src + z_fe0) / 2, x_fe + 2,
            r"$\vec{B}_x = 1.0\ \mathrm{T}$",
            ha="center", va="bottom", fontsize=10,
            color="deepskyblue", fontweight="bold", zorder=4)

    # Flechas de campo magnético (→ hacia +x, representando Bx) — 4 flechas uniformes
    for zb in np.linspace(z_src + 30, z_fe0 - 15, 4):
        ax.annotate("", xy=(zb, 18), xytext=(zb, 3),
                    arrowprops=dict(arrowstyle="-|>", color="deepskyblue",
                                    lw=1.5, mutation_scale=12),
                    zorder=4)

    for tx in np.linspace(-28, 28, 6):
        x_end = tx * (z_L2 + 6 - z_src) / (-z_src)
        ax.plot([z_src, z_L2 + 6], [0, x_end],
                color="darkorange", lw=0.9, alpha=0.55, zorder=1)
        x_hit = tx * dist_L1 / (-z_src)
        if abs(x_hit) < x_bar:
            ax.plot([z_L1], [x_hit], "o", color="darkorange", ms=3.5, zorder=3)

    for tx, z_stop in [(-15, 18), (8, 42), (22, 25), (-28, 55)]:
        x_stop = tx * (z_stop - z_src) / (-z_src)
        ax.plot([z_src, z_stop], [0, x_stop],
                color="steelblue", lw=1.0, ls="--", alpha=0.75, zorder=1)
        ax.plot([z_stop], [x_stop], "x",
                color="steelblue", ms=7, mew=1.8, zorder=3)

    for tx_pi in [-6, 20, -18]:
        x_hit_pi = tx_pi * dist_L1 / (-z_src)
        if abs(x_hit_pi) < x_bar:
            ax.plot([z_L1], [x_hit_pi], "o", color="tomato", ms=5, zorder=3)

    z_ext = z_L2 + 10
    x_ext = x_bar * (z_ext - z_src) / dist_L1
    ax.plot([z_src, z_ext], [0, x_ext], color="#8B008B",
            ls="--", lw=1.3, alpha=0.82, zorder=2)
    ax.plot([z_src, z_ext], [0, -x_ext], color="#8B008B",
            ls="--", lw=1.3, alpha=0.82, zorder=2)
    ax.text(20, 43, rf"$\theta_{{acc}}\approx{theta_acc:.1f}°$",
            fontsize=9, color="#8B008B", va="bottom", ha="left",
            rotation=18, rotation_mode='anchor', zorder=5)

    x_ext_lat = tx_lat * (z_ext - z_src) / (-z_src)
    ax.plot([z_src, z_ext], [0, x_ext_lat], color="goldenrod",
            ls="--", lw=1.3, alpha=0.82, zorder=2)
    ax.plot([z_src, z_ext], [0, -x_ext_lat], color="goldenrod",
            ls="--", lw=1.3, alpha=0.82, zorder=2)
    ax.text(10, 28, rf"$\theta_{{lat}}\approx{theta_lat:.1f}°$",
            fontsize=9, color="goldenrod", va="bottom", ha="left",
            rotation=14, rotation_mode='anchor', zorder=5)

    ax.plot([z_src], [0], "*", color="darkred", ms=13, zorder=5)
    ax.annotate("Fuente\n(0, 0, −2 m)\n→ dirección +z",
                xy=(z_src, 0), xytext=(z_src - 8, -40),
                fontsize=8.5, fontweight="bold", color="darkred", zorder=5,
                arrowprops=dict(arrowstyle="->", color="darkred", lw=1.2),
                bbox=_info_box_kwargs(edgecolor="darkred", alpha=0.92))

    y_arr = -x_fe - 7
    ax.annotate("", xy=(z_fe1, y_arr), xytext=(z_fe0, y_arr),
                arrowprops=dict(arrowstyle="<->", color="black", lw=1.2))
    ax.text((z_fe0 + z_fe1) / 2, y_arr - 1.2,
            "70 cm", ha="center", va="top", fontsize=9)
    ax.annotate("", xy=(z_L1, y_arr), xytext=(z_fe1, y_arr),
                arrowprops=dict(arrowstyle="<->", color="black", lw=1.2))
    ax.text((z_fe1 + z_L1) / 2, y_arr - 1.2,
            "30 cm", ha="center", va="top", fontsize=9)

    legend_elems = [
        Patch(facecolor="deepskyblue", alpha=0.25,
              label=r"Región $B_x = 1.0\ \mathrm{T}$"),
        Line2D([0], [0], marker="*", color="darkred", linestyle="None", ms=11,
               label="fuente puntual"),
        Line2D([0], [0], color="darkorange", lw=2,
               label=r"$\mu^+$ (atraviesa el Fe)"),
        Line2D([0], [0], marker="o", color="darkorange", linestyle="None",
               ms=5, label=r"impacto $\mu^+$ en centellador"),
        Line2D([0], [0], color="steelblue", lw=1.5, ls="--",
               label=r"$\pi^+$ (absorbido en Fe)"),
        Line2D([0], [0], marker="x", color="steelblue", linestyle="None",
               ms=7, mew=1.8, label="absorción hadrónica"),
        Line2D([0], [0], marker="o", color="tomato", linestyle="None",
               ms=5, label=r"impacto $\pi^+$ (≈10%)"),
        Line2D([0], [0], color="#8B008B", lw=1.5, ls="--",
               label=rf"límite geom. $\theta_{{acc}}\approx{theta_acc:.1f}°$"),
        Line2D([0], [0], color="goldenrod", lw=1.5, ls="--",
               label=rf"cara lateral $\theta_{{lat}}\approx{theta_lat:.1f}°$"),
    ]
    ax.legend(handles=legend_elems, fontsize=8.5, loc="upper left", **_legend_kwargs(framealpha=0.88))
    ax.set_xlim(-218, 118)
    ax.set_ylim(-63, 63)
    ax.tick_params(labelsize=9)

    def _layer_panel(ax, horizontal_bars, title, xlabel, ylabel,
                     color, footprint, x_hits, y_hits):
        ax.set_title(title, fontsize=9, pad=4)
        ax.set_xlabel(xlabel, fontsize=9)
        ax.set_ylabel(ylabel, fontsize=9)

        bar_c = np.linspace(-47.5, 47.5, 20)
        if horizontal_bars:
            for xc in bar_c:
                ax.add_patch(Rectangle((xc - 2.5, -50), 5, 100,
                                       facecolor=color, edgecolor=color,
                                       alpha=0.18, lw=0.4, zorder=0))
                ax.axvline(xc + 2.5, color=color, lw=0.3, alpha=0.4)
        else:
            for yc in bar_c:
                ax.add_patch(Rectangle((-50, yc - 2.5), 100, 5,
                                       facecolor=color, edgecolor=color,
                                       alpha=0.18, lw=0.4, zorder=0))
                ax.axhline(yc + 2.5, color=color, lw=0.3, alpha=0.4)

        ax.add_patch(Rectangle((-35, -35), 70, 70,
                               facecolor="none", edgecolor="gray",
                               lw=1.5, ls="-", zorder=2))
        ax.text(0, 0, "huella del Fe  70×70 cm",
                ha="center", va="center", fontsize=7.5, color="gray",
                alpha=0.80, zorder=2)

        ax.add_patch(Rectangle((-footprint, -footprint),
                               2 * footprint, 2 * footprint,
                               facecolor="none", edgecolor="darkorange",
                               lw=1.3, ls="--", zorder=3))
        ax.text(-footprint + 2, -footprint + 1.8,
                rf"huella del haz  $\pm${footprint:.1f} cm",
                ha="left", va="bottom", fontsize=7.5,
                color="darkorange", alpha=0.92, zorder=3,
                bbox=_info_box_kwargs(alpha=0.7, edgecolor="darkorange"))

        mask = (np.abs(x_hits) < x_bar) & (np.abs(y_hits) < x_bar)
        ax.scatter(x_hits[mask], y_hits[mask], color="darkorange",
                   s=10, alpha=0.60, zorder=4,
                   label="impactos (dist. uniforme\nen cara del Fe)")
        ax.legend(fontsize=7.5, loc="lower right", **_legend_kwargs(framealpha=0.85))
        ax.set_xlim(-55, 55)
        ax.set_ylim(-55, 55)
        ax.set_aspect("equal")
        ax.tick_params(labelsize=9)

    n_dots = 220
    tx1 = rng.uniform(-x_fe, x_fe, n_dots)
    ty1 = rng.uniform(-x_fe, x_fe, n_dots)
    x1 = tx1 * dist_L1 / (-z_src)
    y1 = ty1 * dist_L1 / (-z_src)
    _layer_panel(ax_1, horizontal_bars=False,
                 title="Capa 1 — 20 barras a lo largo de X\n"
                       r"z = 100.5 cm  →  detecta posición Y",
                 xlabel="x (cm)  [largo de barra]",
                 ylabel="y (cm)  [posición medida]",
                 color="royalblue", footprint=footprint_L1,
                 x_hits=x1, y_hits=y1)

    tx2 = rng.uniform(-x_fe, x_fe, n_dots)
    ty2 = rng.uniform(-x_fe, x_fe, n_dots)
    x2 = tx2 * dist_L2 / (-z_src)
    y2 = ty2 * dist_L2 / (-z_src)
    _layer_panel(ax_2, horizontal_bars=True,
                 title="Capa 2 — 20 barras a lo largo de Y\n"
                       r"z = 103.5 cm  →  detecta posición X",
                 xlabel="x (cm)  [posición medida]",
                 ylabel="y (cm)  [largo de barra]",
                 color="firebrick", footprint=footprint_L2,
                 x_hits=x2, y_hits=y2)

    fig.text(0.5, 0.005,
             r"Fuente puntual en (0, 0, $-$2 m).  Dirección por evento: apunta a punto uniforme "
             r"en cara del Fe (70×70 cm)  |  Capa 1 + Capa 2 $\rightarrow$ plano sensible X-Y de 1 m × 1 m",
             ha="center", va="bottom", fontsize=8.5,
             bbox=_info_box_kwargs(facecolor="lightyellow", edgecolor="gray", alpha=0.85, pad=3))

    _save(fig, out_dir, "detector_layout.png")


# ============================================================================
# PLOT 5: Overlay mu+ vs pi+ — mediana dE/dx vs beta-gamma
# ============================================================================
def plot_overlay(mu, pi, out_dir):
    fig, ax = plt.subplots(figsize=(10, 7))

    bg_range = (0.3, 100)
    bg_bins = np.logspace(np.log10(bg_range[0]), np.log10(bg_range[1]), 80)
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
    ax.fill_between(bg_centers[v_mu], mu_lo[v_mu], mu_hi[v_mu], color=COLOR_MU, alpha=0.15)
    ax.fill_between(bg_centers[v_pi], pi_lo[v_pi], pi_hi[v_pi], color=COLOR_PI, alpha=0.15)
    ax.plot(bg_centers[v_mu], mu_med[v_mu], color=COLOR_MU, lw=2.5, label=r"$\mu^+$  mediana dE/dx")
    ax.plot(bg_centers[v_pi], pi_med[v_pi], color=COLOR_PI, lw=2.5, label=r"$\pi^+$  mediana dE/dx")

    bg_th = np.logspace(np.log10(bg_range[0]), np.log10(bg_range[1]), 2000)
    mpv = landau_mpv(bg_th, x_mm=10.0, mat=AIR)
    v = ~np.isnan(mpv)
    ax.plot(bg_th[v], mpv[v], color="black", lw=2, ls="--",
            label=r"Landau MPV (BB) — universal $\mu^+\!/\pi^+$")

    ax.axvline(3.5, color="gray", ls=":", lw=1.5, alpha=0.6)
    ax.text(3.5 * 1.05, DEDX_MIN * 1.3, "MIP\n(βγ≈3.5)", color="gray", fontsize=9, va="bottom")

    _setup_ax(ax, log_x=True, log_y=True)
    ax.set_xlabel(r"$\beta\gamma = p\,/\,mc$", fontsize=13)
    ax.set_ylabel(r"$dE/dx$  (MeV/mm)", fontsize=13)
    ax.set_xlim(bg_range); ax.set_ylim(DEDX_MIN, DEDX_MAX)
    ax.legend(fontsize=10, **_legend_kwargs())
    _add_info(ax)
    fig.suptitle(r"Bethe-Bloch overlay: $\mu^+$ vs $\pi^+$ en BC404 (barras 1 cm)", fontsize=14, fontweight="bold", y=0.88)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    _save(fig, out_dir, "bethe_bloch_overlay.png")


# ============================================================================
# PLOT 6: PID combinado
# ============================================================================
def plot_pid_combined(mu, pi, out_dir):
    fig, ax = plt.subplots(figsize=(10, 7))

    p_range_GeV = (0.03, 15.0)
    DEDX_LIN_MAX = 3.0
    p_bins = np.logspace(np.log10(p_range_GeV[0]), np.log10(p_range_GeV[1]), 200)
    dedx_bins = np.linspace(0, DEDX_LIN_MAX, 200)

    def hist2d_arrays(p_MeV, dedx):
        H, xedges, yedges = np.histogram2d(p_MeV / 1000.0, dedx, bins=[p_bins, dedx_bins])
        return H, xedges, yedges

    H_mu, xe, ye = hist2d_arrays(mu["Momentum"], mu["fdEdx"])
    H_pi, _, _ = hist2d_arrays(pi["Momentum"], pi["fdEdx"])

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
    ax.pcolormesh(X, Y, H_pi_n, cmap="Reds", vmin=0, vmax=vmax, alpha=0.65, shading="auto")

    p_th_GeV = np.logspace(np.log10(p_range_GeV[0]), np.log10(p_range_GeV[1]), 2000)
    for mass, color, label in [
        (M_MU, COLOR_MU, r"Landau MPV — $\mu^+$"),
        (M_PI, COLOR_PI, r"Landau MPV — $\pi^+$"),
    ]:
        bg_th = (p_th_GeV * 1000) / mass
        mpv = landau_mpv(bg_th, x_mm=10.0)
        v = ~np.isnan(mpv) & (mpv < DEDX_LIN_MAX)
        ax.plot(p_th_GeV[v], mpv[v], color=color, lw=2.5, ls="-", label=label, zorder=5)

    from matplotlib.patches import Patch
    legend_handles = [
        Patch(facecolor="steelblue", alpha=0.8, label=r"$\mu^+$ datos"),
        Patch(facecolor="tomato", alpha=0.8, label=r"$\pi^+$ datos"),
        plt.Line2D([0], [0], color=COLOR_MU, lw=2.5, label=r"Landau MPV — $\mu^+$"),
        plt.Line2D([0], [0], color=COLOR_PI, lw=2.5, label=r"Landau MPV — $\pi^+$"),
    ]
    ax.legend(handles=legend_handles, fontsize=10, **_legend_kwargs())

    p_mip_pi = 3.5 * M_PI / 1000.0
    ax.axvline(p_mip_pi, color="gold", ls="--", lw=1.2, alpha=0.8)
    ax.text(p_mip_pi * 1.08, 0.05, r"MIP ($\pi^+$)", color="goldenrod", fontsize=8, va="bottom")

    _setup_ax(ax, log_x=True)
    ax.set_xlabel(r"$p$  (GeV/c)", fontsize=13)
    ax.set_ylabel(r"$dE/dx$  (MeV/mm)", fontsize=13)
    ax.set_xlim(p_range_GeV); ax.set_ylim(0, DEDX_LIN_MAX)
    _add_info(ax)
    fig.suptitle(r"PID: $\mu^+$ vs $\pi^+$ en BC404 bar strip detector", fontsize=14, fontweight="bold", y=0.88)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    _save(fig, out_dir, "pid_combined.png")


# ============================================================================
# PLOT 7: Hits por capa
# ============================================================================
def plot_layer_hits(mu, pi, out_dir):
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))

    for ax, data, particle_label, color in [
        (axes[0], mu, r"$\mu^+$", COLOR_MU),
        (axes[1], pi, r"$\pi^+$", COLOR_PI),
    ]:
        layer0 = (data["layerID"] == 0).sum()
        layer1 = (data["layerID"] == 1).sum()
        bars = ax.bar(["Capa 1\n(barras X)", "Capa 2\n(barras Y)"],
                      [layer0, layer1], color=[color, color], alpha=0.7,
                      edgecolor="black", width=0.6)
        for rect, val in zip(bars, [layer0, layer1]):
            ax.text(rect.get_x() + rect.get_width() / 2, rect.get_height() * 1.01,
                    f"{val:,}", ha="center", va="bottom", fontsize=10, fontweight="bold")
        ax.set_ylabel("Número de hits", fontsize=11)
        ax.set_ylim(0, max(layer0, layer1) * 1.15)
        ax.text(0.5, 0.97, f"Hits por capa — {particle_label}",
                transform=ax.transAxes, fontsize=12, ha="center", va="top", fontweight="bold")
        _setup_ax(ax)
        _add_info(ax, particle=particle_label)

    fig.tight_layout()
    _save(fig, out_dir, "layer_hits.png")


# ============================================================================
# PLOT 8: Landau corregida — run 45 con tabla de estadisticas
# ============================================================================
def plot_landau_corregida(mixed_path, out_dir):
    files = sorted(glob.glob(mixed_path), key=run_num)
    fpath = None
    for f in files:
        if "output_run45" in f:
            fpath = f
            break
    if fpath is None:
        print("  SKIP landau_corregida: no output_run45.root found")
        return

    with uproot.open(fpath) as f:
        tree = _best_cycle(f)
        if tree is None:
            print("  SKIP landau_corregida: no tree in output_run45.root")
            return
        d45 = tree.arrays(["fdEdx", "particleID", "fEvent", "layerID", "Momentum"], library="np")

    mu_mask = d45["particleID"] == 0
    pi_mask = d45["particleID"] == 1

    mu_dedx = d45["fdEdx"][mu_mask]
    pi_dedx = d45["fdEdx"][pi_mask]
    mu_events = len(np.unique(d45["fEvent"][mu_mask]))
    pi_events = len(np.unique(d45["fEvent"][pi_mask]))
    tot_events = len(np.unique(d45["fEvent"]))
    mu_hits_tot = len(mu_dedx)
    pi_hits_tot = len(pi_dedx)

    mu_mean, mu_med = np.mean(mu_dedx), np.median(mu_dedx)
    pi_mean, pi_med = np.mean(pi_dedx), np.median(pi_dedx)

    fig, ax = plt.subplots(figsize=(11, 7))
    bins = np.linspace(0.01, 5.0, 100)

    ax.hist(mu_dedx, bins=bins, histtype="stepfilled",
            color=COLOR_MU, alpha=0.5, density=True,
            label=rf"$\mu^+$  ({mu_events} eventos, {mu_hits_tot} hits)")
    ax.hist(mu_dedx, bins=bins, histtype="step",
            color=COLOR_MU, lw=2, density=True)
    ax.hist(pi_dedx, bins=bins, histtype="stepfilled",
            color=COLOR_PI, alpha=0.4, density=True,
            label=rf"$\pi^+$  ({pi_events} eventos, {pi_hits_tot} hits)")
    ax.hist(pi_dedx, bins=bins, histtype="step",
            color=COLOR_PI, lw=2, density=True)

    ax.axvline(0.5, color="green", lw=2.5, ls="--",
               label=r"Umbral dE/dx = 0.5 MeV/mm")
    ax.set_xlabel(r"$dE/dx$ por paso  (MeV/mm)")
    ax.set_ylabel("Densidad (normalizada)")
    ax.set_yscale("log")
    ax.set_xlim(0.01, 5.0)
    ax.set_ylim(1e-3, 5)
    ax.legend(fontsize=10, loc="upper right", **_legend_kwargs())
    ax.set_title(r"Distribución de Landau — haz mixto $\mu^+$/$\pi^+$  ($p_0 \approx 1$ GeV/c)", fontsize=14, fontweight="bold", pad=16)
    _setup_ax(ax, log_y=True)

    stats_text = (
        f"Estadísticas de dE/dx (MeV/mm):\n"
        f"{'':>8} {'μ⁺':>10} {'π⁺':>10}\n"
        f"{'Media':>8} {mu_mean:>8.4f}  {pi_mean:>8.4f}\n"
        f"{'Mediana':>8} {mu_med:>8.4f}  {pi_med:>8.4f}\n"
        f"{'Hits totales':>8} {mu_hits_tot:>8d}  {pi_hits_tot:>8d}"
    )
    ax.text(0.02, 0.98, stats_text, transform=ax.transAxes, fontsize=8.5,
            family="monospace", va="top", ha="left",
            bbox=_info_box_kwargs())

    info1 = (f"Geant4  |  ~1000 μ⁺ + ~1000 π⁺ generados  |  "
             f"{mu_events} μ⁺ detectados + {pi_events} π⁺ detectados = {tot_events} total")
    ax.text(0.01, 1.005, info1, transform=ax.transAxes,
            fontsize=8.5, va="bottom", ha="left", color="#333333", style="italic")

    fig.tight_layout(rect=[0, 0, 1, 0.90])
    _save(fig, out_dir, "landau_corregida.png")



# ============================================================================
# PLOT 10: Eficiencia vs momento corregida (barras error + regimenes)
# ============================================================================
def plot_eff_momento_corregida(mixed_path, out_dir):
    files = sorted(glob.glob(mixed_path), key=run_num)
    p0_vals_GeV = np.logspace(np.log10(0.05), np.log10(10.0), 80)
    records_mu, records_pi = [], []

    for fpath in files:
        rn = run_num(fpath)
        if rn >= 80:
            continue
        p0 = p0_vals_GeV[rn]
        with uproot.open(fpath) as f:
            tree = _best_cycle(f)
            if tree is None:
                continue
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
    p_pi = np.array([r["p0"] for r in pi_recs])
    e_pi = np.array([r["eff"] for r in pi_recs])

    err_mu = np.sqrt(e_mu * (1 - e_mu) / 1000)
    err_pi = np.sqrt(e_pi * (1 - e_pi) / 1000)

    fig, ax = plt.subplots(figsize=(11, 7))

    ax.errorbar(p_mu, e_mu * 100, yerr=err_mu * 100, fmt="o-", color=COLOR_MU,
                lw=2, ms=5, capsize=3, capthick=1, label=r"$\mu^+$", zorder=3)
    ax.errorbar(p_pi, e_pi * 100, yerr=err_pi * 100, fmt="s-", color=COLOR_PI,
                lw=2, ms=5, capsize=3, capthick=1, label=r"$\pi^+$", zorder=3)

    ax.axhline(50, color="gray", ls=":", lw=1, alpha=0.5)
    ax.set_xlabel(r"Momento inicial $p_0$  (GeV/c)")
    ax.set_ylabel(r"Eficiencia $\varepsilon$ (%)")
    ax.set_ylim(-3, 105)
    ax.set_xlim(0.04, 11)
    ax.set_title(r"Eficiencia de detección vs $p_0$", fontsize=14, fontweight="bold", pad=16)
    _setup_ax(ax)

    ax.axvline(0.5, color="gray", ls=":", lw=1.2, alpha=0.6, zorder=0)
    ax.axvline(1.1, color=COLOR_MU, ls=":", lw=1.2, alpha=0.6, zorder=0)

    ax.text(0.085, 55, "Régimen I:\nμ⁺ no penetran\n70 cm Fe",
            fontsize=8.5, color="gray", ha="center", va="bottom",
            alpha=0.8, style="italic")
    ax.text(0.74, 55, "Régimen II:\ntransición μ⁺\nsigmoidal",
            fontsize=8.5, color=COLOR_MU, ha="center", va="bottom",
            alpha=0.8, style="italic")
    ax.text(3.5, 55, "Régimen III:\nmeseta ~89%",
            fontsize=8.5, color="green", ha="center", va="bottom",
            alpha=0.8, style="italic")

    ax.errorbar([], [], fmt="none",
                label=r"Barras: $\sigma_\varepsilon = \sqrt{\varepsilon(1-\varepsilon)/1000}$")
    ax.legend(fontsize=10, **_legend_kwargs(), loc="upper left")

    ax.text(0.01, 1.005, INFO_STR, transform=ax.transAxes,
            fontsize=8.5, va="bottom", ha="left", color="#333333", style="italic")

    fig.tight_layout(rect=[0, 0, 1, 0.90])
    _save(fig, out_dir, "eff_momento_corregida.png")


# ============================================================================
# PLOT 11: Eficiencia vs angulo corregida (ConeAngle + barras error)
# ============================================================================
def plot_eff_angulo_corregida(mixed_path, out_dir):
    files = sorted(glob.glob(mixed_path), key=run_num)
    n_runs = len([f for f in files if run_num(f) < 80])

    L_src = 2000.0
    mu_angles, pi_angles = [], []

    for fpath in files:
        rn = run_num(fpath)
        if rn >= 80:
            continue
        with uproot.open(fpath) as f:
            tree = _best_cycle(f)
            if tree is None:
                continue
            d = tree.arrays(["fEvent", "particleID", "ConeAngle", "layerID"], library="np")

        for pid_val, ang_list in [(0, mu_angles), (1, pi_angles)]:
            mask_pid = d["particleID"] == pid_val
            if mask_pid.sum() == 0:
                continue

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

    rng = np.random.default_rng(42)
    N_mc = 2_000_000
    tx_mc = rng.uniform(-350, 350, N_mc)
    ty_mc = rng.uniform(-350, 350, N_mc)
    theta_mc = np.degrees(np.arctan2(np.sqrt(tx_mc**2 + ty_mc**2), L_src))
    h_mc, _ = np.histogram(theta_mc, bins=bin_edges)
    p_theory = h_mc / h_mc.sum()
    n_total_per_species = n_runs * 1000.0
    n_gen = p_theory * n_total_per_species

    fig, ax = plt.subplots(figsize=(11, 7))

    for angles, color, label, marker in [
        (mu_angles, COLOR_MU, r"$\mu^+$", "o"),
        (pi_angles, COLOR_PI, r"$\pi^+$", "s"),
    ]:
        if not angles:
            continue
        h_det, _ = np.histogram(angles, bins=bin_edges)
        eff = np.where(n_gen > 10, h_det / n_gen, np.nan)
        err = np.where(n_gen > 10, np.sqrt(eff * (1 - eff) / n_gen), np.nan)
        eff = np.clip(eff, 0, 1.05)
        ax.errorbar(bc, eff * 100, yerr=err * 100, fmt=marker + "-", color=color,
                    lw=2, ms=6, capsize=3, capthick=1, label=label)

    ax.set_xlabel(r"Ángulo del cono $\theta$ (°)")
    ax.set_ylabel(r"Eficiencia $\varepsilon$ (%)")
    ax.set_xlim(0, theta_max + 0.3)
    ax.set_ylim(-3, 110)
    ax.axhline(50, color="gray", ls=":", lw=1, alpha=0.5)
    ax.legend(fontsize=12, **_legend_kwargs())
    ax.set_title(r"Eficiencia vs ángulo del cono $\theta$", fontsize=14, fontweight="bold", pad=16)
    _setup_ax(ax)

    theta_lateral = np.degrees(np.arctan(35.0 / 270.0))
    tx_acc_cm = 50.0 * 200.0 / 300.5
    theta_geom = np.degrees(np.arctan(tx_acc_cm / 200.0))

    ax.axvline(theta_lateral, color="goldenrod", ls=":", lw=1.5, alpha=0.8, zorder=1)
    ax.axvline(theta_geom, color="#8B008B", ls=":", lw=1.5, alpha=0.8, zorder=1)

    ax.annotate("Sale por cara\nlateral del Fe",
                xy=(theta_lateral, 30), xytext=(theta_lateral - 2.5, 75),
                fontsize=8.5, color="goldenrod",
                arrowprops=dict(arrowstyle="->", color="goldenrod", lw=1.5),
                bbox=dict(boxstyle="round", facecolor="#fff8dc", alpha=0.85, edgecolor="goldenrod"))

    ax.annotate("Límite geométrico\nbarras ±50 cm",
                xy=(theta_geom, 20), xytext=(theta_geom + 0.5, 75),
                fontsize=8.5, color="#8B008B",
                arrowprops=dict(arrowstyle="->", color="#8B008B", lw=1.5),
                bbox=dict(boxstyle="round", facecolor="#f0e6f6", alpha=0.85, edgecolor="#8B008B"))

    ax.text(0.02, 0.02,
            "θ medido con columna ConeAngle\n"
            "(ángulo inicial traza vs eje z):\n"
            "  · Más preciso que reconstrucción\n"
            "    desde posiciones de barra\n"
            "  · Resolución limitada por paso\n"
            "    de barra de 5 cm (~0.7°)",
            transform=ax.transAxes, fontsize=8, color="gray",
            va="bottom", ha="left", style="italic",
            bbox=dict(boxstyle="round", facecolor="white", alpha=0.8, edgecolor="lightgray"))

    ax.text(0.01, 1.005, INFO_STR, transform=ax.transAxes,
            fontsize=8.5, va="bottom", ha="left", color="#333333", style="italic")

    fig.tight_layout(rect=[0, 0, 1, 0.90])
    _save(fig, out_dir, "eff_angulo_corregida.png")


# ============================================================================
# Main
# ============================================================================
def main(mixed_path, out_dir):
    os.makedirs(out_dir, exist_ok=True)

    print("=" * 60)
    print("Bar Strip Detector — Generating 10 plots")
    print("=" * 60)

    # Load mixed hits for Bethe-Bloch / PID plots
    print("\nLoading mixed hits data...")
    mu, pi = load_mixed_hits(mixed_path)

    # Bethe-Bloch / PID plots (use mu, pi data)
    print("\n[1/10] bethe_bloch_bg.png")
    plot_dedx_vs_bg(mu, pi, out_dir)

    print("[2/10] dedx_vs_beta.png")
    plot_dedx_vs_beta(mu, pi, out_dir)

    print("[3/10] dedx_vs_momentum.png")
    plot_dedx_vs_p(mu, pi, out_dir)

    print("[4/10] detector_layout.png")
    plot_detector_layout(out_dir)

    print("[5/10] bethe_bloch_overlay.png")
    plot_overlay(mu, pi, out_dir)

    print("[6/10] pid_combined.png")
    plot_pid_combined(mu, pi, out_dir)

    print("[7/10] layer_hits.png")
    plot_layer_hits(mu, pi, out_dir)

    # Landau plot (read run 45 directly)
    print("\n[8/10] landau_corregida.png")
    plot_landau_corregida(mixed_path, out_dir)

    # Efficiency plots (read all runs directly)
    print("\n[9/10] eff_momento_corregida.png")
    plot_eff_momento_corregida(mixed_path, out_dir)

    print("[10/10] eff_angulo_corregida.png")
    plot_eff_angulo_corregida(mixed_path, out_dir)

    print("\n" + "=" * 60)
    print(f"All 10 plots saved to: {out_dir}")
    print("=" * 60)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Bar Strip Detector — Generate 10 analysis plots from mixed beam simulation")
    parser.add_argument("--mixed", required=True,
                        help='Glob para ROOT files del dataset mixto, '
                             'ej. "../Classifier/data/mixed_Bfield05T/output_run*.root"')
    parser.add_argument("--out", default="img",
                        help='Directorio de salida para las imágenes (default: img/)')
    args = parser.parse_args()

    main(args.mixed, args.out)
