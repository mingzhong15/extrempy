import pandas as pd
import numpy as np
import glob
import os
import re


def read_thermo_dat(filepath, skip_unstable=True, stable_ratio=0.3, stability_window=300):
    """
    Read thermo.dat file and calculate averaged properties.

    Parameters
    ----------
    filepath : str
        Path to thermo.dat file
    skip_unstable : bool
        Whether to skip unstable initial part
    stable_ratio : float
        Ratio of data to skip from the beginning
    stability_window : int
        Window size for stability check based on temperature/pressure

    Returns
    -------
    data : pandas.DataFrame or None
    averages : dict or None
    stable_start : int or None
    """
    try:
        with open(filepath, 'r') as f:
            first_line = f.readline().strip()
            if first_line.startswith('#'):
                header_line = first_line[1:].strip()
                column_names = header_line.split()
            else:
                f.seek(0)
                column_names = None

        data = pd.read_csv(filepath, sep=r'\s+', comment='#', header=None,
                           skipinitialspace=True)

        if column_names is not None:
            if len(column_names) == len(data.columns):
                data.columns = column_names
            else:
                data = pd.read_csv(filepath, sep=r'\s+', skiprows=1,
                                   header=0, skipinitialspace=True)
    except Exception as e:
        print(f"Error reading {filepath}: {e}")
        return None, None, None

    if skip_unstable:
        n_total = len(data)
        skip_steps = int(n_total * stable_ratio)

        if n_total > stability_window * 2:
            if 'temp[K]' in data.columns:
                temp_std = data['temp[K]'].rolling(window=stability_window, min_periods=100).std()
                temp_threshold = temp_std.quantile(0.1) * 2
                stable_idx = np.where(temp_std < temp_threshold)[0]
                if len(stable_idx) > 0:
                    stable_start = max(skip_steps, stable_idx[0])
                else:
                    stable_start = skip_steps
            else:
                stable_start = skip_steps
        else:
            stable_start = skip_steps
    else:
        stable_start = 0

    data_stable = data.iloc[stable_start:].copy()

    averages = {}
    for col in data_stable.columns:
        if col != 'step':
            col_data = data_stable[col].values
            averages[col] = {
                'mean': np.mean(col_data),
                'std': np.std(col_data),
                'min': np.min(col_data),
                'max': np.max(col_data)
            }

    return data, averages, stable_start


def detect_phase_from_msd(msd_data, steps, time_step=0.001, slope_threshold=0.001):
    """
    Automatically detect phase (solid or liquid) from MSD data.

    Parameters
    ----------
    msd_data : array
        MSD values over time
    steps : array
        Step numbers
    time_step : float
        Time step in ps per step
    slope_threshold : float
        Threshold for slope (ang^2/ps^2)

    Returns
    -------
    phase : str or None ('solid', 'liquid', or None)
    slope : float
    r_squared : float
    """
    skip_ratio = 0.3
    skip_idx = int(len(msd_data) * skip_ratio)

    if len(msd_data) - skip_idx < 100:
        return None, 0.0, 0.0

    time = steps[skip_idx:] * time_step
    msd = msd_data[skip_idx:]

    try:
        from scipy import stats
        slope, intercept, r_value, p_value, std_err = stats.linregress(time, msd)
        r_squared = r_value ** 2
    except Exception:
        coeffs = np.polyfit(time, msd, 1)
        slope = coeffs[0]
        msd_pred = slope * time + coeffs[1]
        ss_res = np.sum((msd - msd_pred) ** 2)
        ss_tot = np.sum((msd - np.mean(msd)) ** 2)
        r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0.0

    if slope < slope_threshold:
        phase = 'solid'
    else:
        phase = 'liquid'

    return phase, slope, r_squared


def calculate_solid_msd_plateau(msd_data, steps, window_size=1000, min_plateau_length=500):
    """
    Calculate plateau average MSD for solid phase.

    Returns
    -------
    plateau_mean : float
    plateau_std : float
    plateau_start_idx : int
    """
    df = pd.DataFrame({'step': steps, 'msd': msd_data})
    df['rolling_mean'] = df['msd'].rolling(window=window_size, min_periods=100).mean()
    df['rolling_std'] = df['msd'].rolling(window=window_size, min_periods=100).std()

    skip_ratio = 0.3
    skip_idx = int(len(df) * skip_ratio)

    df_valid = df.iloc[skip_idx:].copy()
    df_valid['rel_std'] = df_valid['rolling_std'] / (df_valid['rolling_mean'] + 1e-10)

    min_rel_std_idx = df_valid['rel_std'].idxmin()
    plateau_start_idx = max(skip_idx, min_rel_std_idx - min_plateau_length // 2)
    plateau_end_idx = min(len(df), plateau_start_idx + min_plateau_length)

    plateau_data = msd_data[plateau_start_idx:plateau_end_idx]
    plateau_mean = np.mean(plateau_data)
    plateau_std = np.std(plateau_data)

    return plateau_mean, plateau_std, plateau_start_idx


def calculate_liquid_diffusion_coefficient(msd_data, steps, time_step=0.001):
    """
    Calculate diffusion coefficient from MSD linear fit.

    Returns
    -------
    D : float
        Diffusion coefficient in cm^2/s
    slope : float
    r_squared : float
    """
    skip_ratio = 0.3
    skip_idx = int(len(msd_data) * skip_ratio)

    time = steps[skip_idx:] * time_step
    msd = msd_data[skip_idx:]

    if len(time) < 2:
        return 0.0, 0.0, 0.0

    try:
        from scipy import stats
        slope, intercept, r_value, p_value, std_err = stats.linregress(time, msd)
        r_squared = r_value ** 2
    except Exception:
        coeffs = np.polyfit(time, msd, 1)
        slope = coeffs[0]
        msd_pred = slope * time + coeffs[1]
        ss_res = np.sum((msd - msd_pred) ** 2)
        ss_tot = np.sum((msd - np.mean(msd)) ** 2)
        r_squared = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0.0

    D_ang2_ps = slope / 6.0
    D_cm2_s = D_ang2_ps * 1e-4

    return D_cm2_s, slope, r_squared


def find_first_peak_rdf(r, g_r, r_min=0.5, r_max=5.0):
    """Find the first peak position in RDF."""
    if r is None or g_r is None or len(r) == 0:
        return None, None

    mask = (r >= r_min) & (r <= r_max)
    r_filtered = r[mask]
    g_filtered = g_r[mask]

    if len(r_filtered) == 0:
        return None, None

    peak_idx = np.argmax(g_filtered)
    peak_r = r_filtered[peak_idx]
    peak_g = g_filtered[peak_idx]

    return peak_r, peak_g


def process_npt_directories(base_dir, element, phases=None, skip_unstable=True,
                            stable_ratio=0.3, structure='fcc'):
    """
    Process all NPT simulation directories and extract averaged properties.

    Returns
    -------
    summary : pandas.DataFrame
    all_data : dict
    """
    if phases is None:
        phases = ['solid', 'liquid']

    summary_list = []
    all_data = {}

    for phase in phases:
        dirs = []
        dir_pattern = os.path.join(base_dir, f'{element}_*k_npt_{phase}')
        found_dirs = glob.glob(dir_pattern)
        dirs.extend(found_dirs)
        dirs = sorted(list(set(dirs)))

        if len(dirs) == 0:
            print(f"Warning: No directories found for pattern {element}_*k_npt_{phase}")
            continue

        print(f"Found {len(dirs)} directories for {phase} phase")

        for dir_path in dirs:
            thermo_file = os.path.join(dir_path, 'thermo.dat')
            if not os.path.exists(thermo_file):
                print(f"Warning: {thermo_file} not found, skipping...")
                continue

            match = re.search(r'(\d+\.?\d*)k', os.path.basename(dir_path), re.IGNORECASE)
            if not match:
                match = re.search(r'(\d+\.?\d*)\.k', os.path.basename(dir_path), re.IGNORECASE)
            if match:
                temperature = float(match.group(1))
            else:
                print(f"Warning: Could not extract temperature from {dir_path}, skipping...")
                continue

            data, averages, stable_start = read_thermo_dat(
                thermo_file, skip_unstable=skip_unstable, stable_ratio=stable_ratio
            )

            if data is None or averages is None:
                print(f"Warning: Failed to read {thermo_file}")
                continue

            detected_phase = None
            msd_col = None
            for col in data.columns:
                if 'MSD' in col.upper():
                    msd_col = col
                    break

            if msd_col:
                steps = data['step'].values
                msd_data = data[msd_col].values
                detected_phase, msd_slope, msd_r2 = detect_phase_from_msd(msd_data, steps)

                if detected_phase is not None and detected_phase != phase:
                    print(f"Warning: Directory '{phase}' but MSD suggests '{detected_phase}', skipping.")
                    continue

            if detected_phase is None:
                detected_phase = phase

            if msd_col:
                steps = data['step'].values
                msd_data = data[msd_col].values

                if detected_phase == 'solid':
                    plateau_mean, plateau_std, _ = calculate_solid_msd_plateau(msd_data, steps)
                    averages['MSD_plateau'] = {
                        'mean': plateau_mean, 'std': plateau_std,
                        'min': plateau_mean - plateau_std,
                        'max': plateau_mean + plateau_std
                    }
                elif detected_phase == 'liquid':
                    D, slope, r2 = calculate_liquid_diffusion_coefficient(msd_data, steps)
                    averages['diffusion_coefficient'] = {'mean': D, 'std': 0, 'min': D, 'max': D}
                    averages['MSD_slope'] = {'mean': slope, 'std': 0, 'min': slope, 'max': slope}
                    averages['MSD_fit_r2'] = {'mean': r2, 'std': 0, 'min': r2, 'max': r2}

            summary_row = {'temperature': temperature, 'phase': detected_phase}
            for col_name, stats_dict in averages.items():
                for stat_name, stat_value in stats_dict.items():
                    summary_row[f'{col_name}_{stat_name}'] = stat_value

            summary_list.append(summary_row)
            all_data[(temperature, detected_phase)] = {'data': data, 'averages': averages}

    summary = pd.DataFrame(summary_list)

    if len(summary) > 0:
        summary = summary.sort_values('temperature').reset_index(drop=True)

    return summary, all_data


def read_dicts_from_file(filename):
    """Read Python dict lines from a data.log-style file."""
    import ast
    dict_list = []
    with open(filename, 'r', encoding='utf-8') as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()
            if line:
                try:
                    data_dict = ast.literal_eval(line)
                    dict_list.append(data_dict)
                except (SyntaxError, ValueError) as e:
                    print(f"Parse error at line {line_num}: {line[:50]}...")
    return dict_list


def write_dict_to_file(filename, data_dict):
    """Append a Python dict as a single line to a data.log-style file."""
    with open(filename, 'a', encoding='utf-8') as f:
        f.write(str(data_dict) + '\n')


def plot_thermo_summary(summary, element, ax=None):
    """Plot temperature vs volume, density, and energy from NPT summary.

    Parameters
    ----------
    summary : pandas.DataFrame
        Output of ``process_npt_directories()``.
    element : str
        Element label for titles.
    ax : array of matplotlib.axes.Axes or None
        If None, create new figure with 3 subplots.
    """
    import matplotlib.pyplot as plt

    if ax is None:
        fig, ax = plt.subplots(1, 3, figsize=(12, 3.5))
    if not isinstance(ax, (list, np.ndarray)):
        ax = [ax]

    cols = summary.columns
    vol_col = next((c for c in cols if 'VOL' in c or 'vol' in c), None)
    rho_col = next((c for c in cols if 'RHO' in c or 'rho' in c or 'density' in c), None)
    ener_col = next((c for c in cols if 'ETOTAL' in c or 'etotal' in c or 'energy' in c), None)

    for i, (ycol, ylabel) in enumerate(zip(
            [vol_col, rho_col, ener_col],
            ['Volume (A^3)', 'Density (g/cc)', 'Energy (eV)'])):
        if ycol is None or ycol not in cols:
            continue
        mean_col = ycol if ycol.endswith('_mean') else f'{ycol}_mean'
        std_col = ycol.replace('_mean', '_std') if '_mean' in ycol else f'{ycol}_std'

        if mean_col in cols:
            for phase in summary['phase'].unique():
                mask = summary['phase'] == phase
                x = summary.loc[mask, 'temperature']
                y = summary.loc[mask, mean_col]
                yerr = summary.loc[mask, std_col] if std_col in cols else None
                ax[i].errorbar(x, y, yerr=yerr, fmt='o-', label=phase, capsize=3)

        ax[i].set_xlabel('Temperature (K)')
        ax[i].set_ylabel(ylabel)
        ax[i].set_title(f'{element} — {ylabel}')
        ax[i].legend()
        ax[i].grid(True, alpha=0.3)

    plt.tight_layout()


# Predefined experimental and DP melting points (from notebook melt_list)
MELT_DATA = {
    'fcc': {
        'Al':  {'Tm': 933.0,  'Tm_PBEsol_DP': 855.0,  'Tm_PBEsol_DP_error': 20.0},
        'Au':  {'Tm': 1337.0, 'Tm_PBEsol_DP': 1480.0, 'Tm_PBEsol_DP_error': 30.0},
        'Cu':  {'Tm': 1358.0, 'Tm_PBEsol_DP': 1220.0, 'Tm_PBEsol_DP_error': 25.0},
        'Ag':  {'Tm': 1235.0, 'Tm_PBEsol_DP': 1120.0, 'Tm_PBEsol_DP_error': 25.0},
        'Pb':  {'Tm': 601.0,  'Tm_PBEsol_DP': 460.0,  'Tm_PBEsol_DP_error': 20.0},
        'Ni':  {'Tm': 1728.0, 'Tm_PBEsol_DP': 1650.0, 'Tm_PBEsol_DP_error': 30.0},
    },
    'bcc': {
        'Li':  {'Tm': 454.0,  'Tm_PBEsol_DP': 500.0,  'Tm_PBEsol_DP_error': 30.0},
        'Na':  {'Tm': 371.0,  'Tm_PBEsol_DP': 440.0,  'Tm_PBEsol_DP_error': 30.0},
        'K':   {'Tm': 337.0,  'Tm_PBEsol_DP': 400.0,  'Tm_PBEsol_DP_error': 30.0},
        'V':   {'Tm': 2183.0, 'Tm_PBEsol_DP': 2120.0, 'Tm_PBEsol_DP_error': 40.0},
        'Nb':  {'Tm': 2750.0, 'Tm_PBEsol_DP': 2680.0, 'Tm_PBEsol_DP_error': 40.0},
        'Ta':  {'Tm': 3290.0, 'Tm_PBEsol_DP': 3200.0, 'Tm_PBEsol_DP_error': 50.0},
        'Cr':  {'Tm': 2180.0, 'Tm_PBEsol_DP': 2150.0, 'Tm_PBEsol_DP_error': 35.0},
        'Mo':  {'Tm': 2896.0, 'Tm_PBEsol_DP': 2820.0, 'Tm_PBEsol_DP_error': 40.0},
        'W':   {'Tm': 3687.0, 'Tm_PBEsol_DP': 3650.0, 'Tm_PBEsol_DP_error': 50.0},
        'Fe':  {'Tm': 1811.0, 'Tm_PBEsol_DP': 1700.0, 'Tm_PBEsol_DP_error': 40.0},
    },
    'hcp': {
        'Be':  {'Tm': 1560.0, 'Tm_PBEsol_DP': 1650.0, 'Tm_PBEsol_DP_error': 40.0},
        'Mg':  {'Tm': 923.0,  'Tm_PBEsol_DP': 940.0,  'Tm_PBEsol_DP_error': 25.0},
        'Ti':  {'Tm': 1941.0, 'Tm_PBEsol_DP': 1860.0, 'Tm_PBEsol_DP_error': 40.0},
        'Zr':  {'Tm': 2128.0, 'Tm_PBEsol_DP': 2120.0, 'Tm_PBEsol_DP_error': 50.0},
        'Zn':  {'Tm': 693.0,  'Tm_PBEsol_DP': 610.0,  'Tm_PBEsol_DP_error': 20.0},
    },
    'diamond': {
        'Si':  {'Tm': 1687.0, 'Tm_PBEsol_DP': 1890.0, 'Tm_PBEsol_DP_error': 40.0},
        'Ge':  {'Tm': 1211.0, 'Tm_PBEsol_DP': 1300.0, 'Tm_PBEsol_DP_error': 30.0},
    },
}

MELT_LIST = {}
for phase, phase_data in MELT_DATA.items():
    for elem, data in phase_data.items():
        MELT_LIST[elem] = {**data, 'phase': phase}
