"""Publication-grade EOS figure for NPT scan results.

Reads a DataFrame produced by :meth:`ElementEOSCalculator.analyze_npt`
and renders a 2×3 multi-panel figure:

    +---------+---------+---------+
    |  ① ρ-T  |  ② V-T  |  ③ P-T  |   EOS row
    +---------+---------+---------+
    |  ④ H-T  | ⑤ MSD-T |  ⑥ D-T  |   thermo + dynamics row
    +---------+---------+---------+

Panels ①-⑤ use the mean/std already stored in the DataFrame (so a
single CSV is enough).  Panel ⑥ fits the diffusion coefficient D from
the raw ``thermo.dat`` time series of each *liquid* directory (Einstein
relation, D = slope(MSD vs t) / 6); this is the only place where raw
files are read.

All I/O is concentrated in :func:`_collect_diffusions`; the six
:func:`_panel_*` helpers and :func:`plot_eos_figure` are pure-data and
unit-testable.
"""

import os

import numpy as np


# ---- style (Nature / matplotlib quick-start from nature-figure) -----------

PALETTE = {
    'solid':        '#0F4D92',   # blue_main
    'liquid':       '#B64342',   # red_strong
    'tm_line':      '#767676',   # neutral_mid
    'reference':    '#CFCECE',   # neutral_light
    'fit':          '#4D4D4D',   # neutral_dark
}

_RC = {
    'font.family':       'sans-serif',
    'font.sans-serif':   ['Arial', 'Helvetica', 'DejaVu Sans', 'sans-serif'],
    'svg.fonttype':      'none',     # editable text in SVG
    'pdf.fonttype':      42,         # editable TrueType text in PDF
    'font.size':         7,
    'axes.linewidth':    0.8,
    'axes.spines.right': False,
    'axes.spines.top':   False,
    'legend.frameon':   False,
    'xtick.major.size': 3,
    'ytick.major.size': 3,
    'xtick.major.width': 0.8,
    'ytick.major.width': 0.8,
    'lines.linewidth':  0.8,
}


def _apply_rc(plt):
    plt.rcParams.update(_RC)


# ---- column-name conventions ----------------------------------------------
# The LAMMPS thermo header written by npt-*.j2 templates.  Centralised
# here so a future header change only edits one mapping.

COL = {
    'temp':     'temperature',
    'phase':    'phase',
    'natoms':   'natoms',
    'T':        'temp[K]_mean',
    'T_std':    'temp[K]_std',
    'P':        'press[bars]_mean',
    'P_std':    'press[bars]_std',
    'V':        'vol[A^3]_mean',
    'V_std':    'vol[A^3]_std',
    'rho':      'density[gcc]_mean',
    'rho_std':  'density[gcc]_std',
    'E':        'etotal[eV]_mean',
    'H':        'enthalpy[eV]_mean',
    'H_std':    'enthalpy[eV]_std',
    'MSD':      'MSD[ang^2/ps]_mean',
    'MSD_std':  'MSD[ang^2/ps]_std',
    # raw thermo.dat column names (for D fit)
    'raw_step': 'step',
    'raw_msd':  'MSD[ang^2/ps]',
}


# ---- I/O: only place that reads raw thermo.dat ---------------------------

def _collect_diffusions(df, work_root, element, structure=None,
                        dt=0.001, skip_ratio=0.3):
    """Fit D (cm²/s) for every liquid row from raw ``thermo.dat``.

    Walks the ``work_root/{element}/npt/`` directory, matches each
    liquid row of ``df`` to its ``{T}k*_{phase}/`` folder, reads the
    raw (step, MSD) time series, and applies the Einstein relation
    via :func:`extrempy.md.thermo.calculate_liquid_diffusion_coefficient`.

    Returns
    -------
    dict
        ``{temperature: {'D': float, 'r2': float}}`` keyed by the
        liquid temperatures present in ``df``.  Rows without a
        matching directory or with too few points are silently
        skipped.
    """
    from extrempy.md.thermo import calculate_liquid_diffusion_coefficient

    liquid = df[df[COL['phase']] == 'liquid']
    if liquid.empty:
        return {}

    npt_dir = os.path.join(work_root, element, 'npt')
    out = {}
    for _, row in liquid.iterrows():
        T = int(row[COL['temp']])
        # match {T}k*_{phase}/  (phase is always 'liquid' here)
        match = None
        if os.path.isdir(npt_dir):
            for name in os.listdir(npt_dir):
                if name.startswith(f'{T}k') and name.endswith('_liquid'):
                    match = os.path.join(npt_dir, name)
                    break
        if match is None:
            continue
        thermo_path = os.path.join(match, 'thermo.dat')
        if not os.path.exists(thermo_path):
            continue
        steps, msd = _read_msd_series(thermo_path)
        if steps is None or len(steps) < 10:
            continue
        D, _slope, r2 = calculate_liquid_diffusion_coefficient(
            msd, steps, time_step=dt)
        out[T] = {'D': D, 'r2': r2}
    return out


def _read_msd_series(thermo_path):
    """Read (step, MSD) time series from a raw thermo.dat."""
    import pandas as pd
    try:
        data = pd.read_csv(thermo_path, sep=r'\s+', comment='#',
                           header=None, skipinitialspace=True)
    except Exception:
        return None, None
    # column names from the leading '# ...' line
    with open(thermo_path) as f:
        first = f.readline().strip()
    names = first.lstrip('#').split() if first.startswith('#') else None
    if names is None or len(names) != data.shape[1]:
        return None, None
    data.columns = names
    if 'step' not in data.columns or COL['raw_msd'] not in data.columns:
        return None, None
    return data['step'].values, data[COL['raw_msd']].values


# ---- per-panel plotters (pure data → ax) ----------------------------------

def _panel_density(ax, df, Tm):
    """①  ρ vs T."""
    for phase, color in [('solid', PALETTE['solid']),
                         ('liquid', PALETTE['liquid'])]:
        sub = df[df[COL['phase']] == phase]
        if sub.empty:
            continue
        ax.errorbar(sub[COL['temp']], sub[COL['rho']],
                    yerr=sub[COL['rho_std']], fmt='o', ms=4, mew=0.6,
                    color=color, mfc='white', lw=0.8, capsize=1.5,
                    capthick=0.6, label=phase)
        _label_end(ax, sub[COL['temp']], sub[COL['rho']], phase, color)
    ax.axvline(Tm, color=PALETTE['tm_line'], lw=0.6, ls='--', zorder=0)
    ax.set_xlabel('Temperature (K)')
    ax.set_ylabel(r'Density $\rho$ (g cm$^{-3}$)')


def _panel_volume(ax, df, Tm, fit_alpha=True):
    """②  V/N vs T."""
    import numpy as np
    for phase, color in [('solid', PALETTE['solid']),
                         ('liquid', PALETTE['liquid'])]:
        sub = df[df[COL['phase']] == phase].dropna(subset=[COL['natoms']])
        if sub.empty:
            continue
        natoms = sub[COL['natoms']].iloc[0]
        V_atom = sub[COL['V']] / natoms
        V_err = sub[COL['V_std']] / natoms
        ax.errorbar(sub[COL['temp']], V_atom, yerr=V_err, fmt='o', ms=4,
                    mew=0.6, color=color, mfc='white', lw=0.8,
                    capsize=1.5, capthick=0.6, label=phase)
        if fit_alpha and len(sub) >= 2:
            _plot_linear_fit(ax, sub[COL['temp']], V_atom, color, ls=':')
        _label_end(ax, sub[COL['temp']], V_atom, phase, color)
    ax.axvline(Tm, color=PALETTE['tm_line'], lw=0.6, ls='--', zorder=0)
    ax.set_xlabel('Temperature (K)')
    ax.set_ylabel(r'Volume per atom $V/N$ ($\rm \AA^3$)')


def _panel_pressure(ax, df, Tm):
    """③  P vs T (sanity: all points near target pressure)."""
    target = 1.0  # bar — set by melt.py default pressure=0.0001 kbar
    for phase, color in [('solid', PALETTE['solid']),
                         ('liquid', PALETTE['liquid'])]:
        sub = df[df[COL['phase']] == phase]
        if sub.empty:
            continue
        ax.errorbar(sub[COL['temp']], sub[COL['P']], yerr=sub[COL['P_std']],
                    fmt='o', ms=4, mew=0.6, color=color, mfc='white',
                    lw=0.8, capsize=1.5, capthick=0.6, label=phase)
    ax.axhline(target, color=PALETTE['reference'], lw=0.6, ls='--', zorder=0)
    ax.axvline(Tm, color=PALETTE['tm_line'], lw=0.6, ls='--', zorder=0)
    ax.set_xlabel('Temperature (K)')
    ax.set_ylabel('Pressure (bar)')


def _panel_enthalpy(ax, df, Tm, dh_f=None):
    """④  H/N vs T, with latent heat ΔH_f annotation."""
    for phase, color in [('solid', PALETTE['solid']),
                         ('liquid', PALETTE['liquid'])]:
        sub = df[df[COL['phase']] == phase].dropna(subset=[COL['natoms']])
        if sub.empty:
            continue
        natoms = sub[COL['natoms']].iloc[0]
        H_atom = sub[COL['H']] / natoms
        H_err = sub[COL['H_std']] / natoms
        ax.errorbar(sub[COL['temp']], H_atom, yerr=H_err, fmt='o', ms=4,
                    mew=0.6, color=color, mfc='white', lw=0.8,
                    capsize=1.5, capthick=0.6, label=phase)
        _plot_linear_fit(ax, sub[COL['temp']], H_atom, color, ls=':')
    ax.axvline(Tm, color=PALETTE['tm_line'], lw=0.6, ls='--', zorder=0)
    if dh_f is not None:
        ax.annotate(rf'$\Delta H_{{\rm f}}$ = {dh_f * 1000:.1f} meV/atom',
                    xy=(0.55, 0.05), xycoords='axes fraction', fontsize=6,
                    color=PALETTE['fit'])
    ax.set_xlabel('Temperature (K)')
    ax.set_ylabel(r'Enthalpy per atom $H/N$ (eV)')


def _panel_msd(ax, df, Tm):
    """⑤  log10(MSD) vs T."""
    for phase, color in [('solid', PALETTE['solid']),
                         ('liquid', PALETTE['liquid'])]:
        sub = df[df[COL['phase']] == phase]
        if sub.empty:
            continue
        msd = sub[COL['MSD']].clip(lower=1e-4)   # avoid log10(0)
        ax.plot(sub[COL['temp']], np.log10(msd), 'o', ms=4, mew=0.6,
                color=color, mfc='white', lw=0.8, label=phase)
        _label_end(ax, sub[COL['temp']], np.log10(msd), phase, color)
    ax.axvline(Tm, color=PALETTE['tm_line'], lw=0.6, ls='--', zorder=0)
    ax.set_xlabel('Temperature (K)')
    ax.set_ylabel(r'$\log_{10}$ MSD ($\rm \AA^2$)')


def _panel_arrhenius(ax, diffusions, ea=None):
    """⑥  ln(D) vs 1000/T for liquid, with Arrhenius fit."""
    from extrempy.md.thermo import fit_arrhenius
    if not diffusions:
        ax.text(0.5, 0.5, 'no liquid data', ha='center', va='center',
                transform=ax.transAxes, fontsize=7, color='gray')
        ax.set_xlabel(r'$1000/T$ (K$^{-1}$)')
        ax.set_ylabel(r'$\ln D$ (cm$^2$ s$^{-1}$)')
        return
    Ts = sorted(diffusions.keys())
    Ds = [diffusions[T]['D'] for T in Ts]
    inv = [1000.0 / T for T in Ts]
    lnD = [np.log(D) for D in Ds if D > 0]
    inv_plot = [1000.0 / T for T, D in zip(Ts, Ds) if D > 0]
    ax.plot(inv_plot, lnD, 'o', ms=4, mew=0.6, color=PALETTE['liquid'],
            mfc='white', lw=0.8)
    fit = fit_arrhenius(Ts, Ds)
    if fit['Ea_eV'] is not None and fit['r_squared'] > 0.5:
        x = np.linspace(min(inv_plot), max(inv_plot), 50)
        # ln D = ln D0 - Ea/(k T) = ln D0 - (Ea/k) * (1/T)
        #       = ln D0 - (Ea/k) * (x/1000)
        y = np.log(fit['D0']) - (fit['Ea_eV'] / 8.617333262e-5) * (x / 1000.0)
        ax.plot(x, y, '-', color=PALETTE['fit'], lw=0.8)
        ax.annotate(rf'$E_{{\rm a}}$ = {fit["Ea_eV"]:.2f} eV'
                    '\n' + rf'$R^2$ = {fit["r_squared"]:.3f}',
                    xy=(0.55, 0.05), xycoords='axes fraction',
                    fontsize=6, color=PALETTE['fit'])
    ax.set_xlabel(r'$1000/T$ (K$^{-1}$)')
    ax.set_ylabel(r'$\ln D$ (cm$^2$ s$^{-1}$)')


# ---- small helpers --------------------------------------------------------

def _label_end(ax, x_series, y_series, label, color):
    """Direct-label the last point of a series (avoid legend eye travel)."""
    ax.annotate(label, xy=(x_series.iloc[-1], y_series.iloc[-1]),
                xytext=(4, 0), textcoords='offset points',
                fontsize=6, color=color, va='center')


def _plot_linear_fit(ax, x, y, color, ls='--', n=50):
    """Overlay a thin dashed linear fit line."""
    import numpy as np
    if len(x) < 2:
        return
    slope, intercept = np.polyfit(x, y, 1)
    xs = np.linspace(x.min(), x.max(), n)
    ax.plot(xs, slope * xs + intercept, ls=ls, lw=0.5,
            color=color, alpha=0.6, zorder=0)


# ---- main entry point -----------------------------------------------------

def plot_eos_figure(df, element, Tm, out_dir, work_root=None,
                    structure=None, dt=0.001, dpi=300):
    """Render a 2×3 EOS figure and save SVG / PDF / PNG.

    Parameters
    ----------
    df : pandas.DataFrame
        Output of :meth:`ElementEOSCalculator.analyze_npt`.  Must
        contain ``_mean`` and ``_std`` columns plus ``natoms``.
    element : str
    Tm : float
        Melting point (K) for the vertical reference line.
    out_dir : str
        Output directory (created if missing).
    work_root : str or None
        Required only for panel ⑥ (diffusion fit): used to locate the
        raw ``thermo.dat`` files under ``{work_root}/{element}/npt/``.
        If ``None``, panel ⑥ is left empty.
    structure : str or None
        Crystal-structure tag (e.g. ``'fcc'``) included in the output
        filename.  Defaults to a per-row column lookup.
    dt : float
        MD timestep (ps) — used by the diffusion fit.
    dpi : int
        PNG export resolution.
    """
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    _apply_rc(plt)

    # ΔH_f via linear extrapolation to Tm
    from extrempy.md.thermo import extrapolate_enthalpy_to_tm
    dh = extrapolate_enthalpy_to_tm(df, Tm)
    dh_f = dh['dH_f']

    # Diffusion coefficients from raw thermo.dat (liquid only)
    diffusions = (_collect_diffusions(df, work_root, element, structure, dt)
                  if work_root else {})

    fig, axes = plt.subplots(2, 3, figsize=(7.09, 4.5),
                             constrained_layout=True)

    _panel_density(axes[0, 0], df, Tm)
    _panel_volume(axes[0, 1], df, Tm)
    _panel_pressure(axes[0, 2], df, Tm)
    _panel_enthalpy(axes[1, 0], df, Tm, dh_f=dh_f)
    _panel_msd(axes[1, 1], df, Tm)
    _panel_arrhenius(axes[1, 2], diffusions)

    suffix = f'_{structure}' if structure else ''
    base = os.path.join(out_dir, f'{element}_eos{suffix}')
    os.makedirs(out_dir, exist_ok=True)
    for ext in ('svg', 'pdf', 'png'):
        fig.savefig(f'{base}.{ext}', dpi=dpi if ext == 'png' else None,
                    bbox_inches='tight')
    plt.close(fig)

    # also dump the fitted physical quantities
    import json
    summary = {
        'element': element,
        'structure': structure,
        'Tm': Tm,
        'dH_f_eV_per_atom': dh_f,
        'H_solid_at_Tm': dh['H_solid_Tm'],
        'H_liquid_at_Tm': dh['H_liquid_Tm'],
        'solid_Cp_slope_eV_per_atom_K': dh['solid_slope'],
        'liquid_Cp_slope_eV_per_atom_K': dh['liquid_slope'],
        'diffusions': diffusions,
    }
    with open(f'{base}_fit.json', 'w') as f:
        json.dump(summary, f, indent=2)
    return summary
