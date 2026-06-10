import numpy as np

try:
    from extrempy.constant import *
except ImportError:
    from ..constant import *


def read_dump_file(dump_path):
    """Read LAMMPS dump file using ASE."""
    import ase.io
    try:
        atoms = ase.io.read(dump_path, format='lammps-dump-text')
        return atoms
    except Exception as e:
        print(f"Error reading {dump_path}: {e}")
        return None


def calculate_rdf(atoms, r_max=10.0, bins=100, cutoff=None):
    """
    Calculate radial distribution function (RDF).

    Returns
    -------
    r : ndarray
    g_r : ndarray
    """
    from ase.neighborlist import NeighborList

    if cutoff is None:
        cutoff = r_max

    positions = atoms.get_positions()
    cell = atoms.get_cell()
    volume = atoms.get_volume()
    n_atoms = len(atoms)

    r = np.linspace(0, r_max, bins)
    dr = r[1] - r[0]

    nl = NeighborList([cutoff / 2] * n_atoms, self_interaction=False,
                      bothways=True, skin=0.0)
    nl.update(atoms)

    distances = []
    for i in range(n_atoms):
        indices, offsets = nl.get_neighbors(i)
        if len(indices) > 0:
            pos_i = positions[i]
            pos_js = positions[indices] + np.dot(offsets, cell)
            vecs = pos_js - pos_i
            dists = np.linalg.norm(vecs, axis=1)
            valid = (dists < r_max) & (indices > i)
            distances.extend(dists[valid].tolist())

    if len(distances) > 0:
        distances = np.array(distances)
        hist, bin_edges = np.histogram(distances, bins=bins, range=(0, r_max))
        rho = n_atoms / volume
        r_centers = (bin_edges[:-1] + bin_edges[1:]) / 2.0
        r_shells = 4 * np.pi * r_centers ** 2 * dr
        g_r = hist / (r_shells * rho * n_atoms)
        g_r[0] = 0.0
    else:
        g_r = np.zeros_like(r)

    return r, g_r


def find_first_minimum_rdf(r, g_r, r_min=2.0, r_max=6.0):
    """
    Find the first minimum in RDF for Q4/Q6 cutoff determination.

    Returns
    -------
    cutoff : float
    r_min_idx : int or None
    """
    from scipy.signal import argrelextrema

    mask = (r >= r_min) & (r <= r_max)
    if not np.any(mask):
        return 3.5, None

    r_search = r[mask]
    g_r_search = g_r[mask]

    min_indices = argrelextrema(g_r_search, np.less, order=3)[0]

    if len(min_indices) > 0:
        first_min_idx = min_indices[0]
        cutoff = r_search[first_min_idx]
        original_indices = np.where(mask)[0]
        r_min_idx = original_indices[first_min_idx]
        return cutoff, r_min_idx
    else:
        return 3.5, None


def calculate_coordination_number(atoms, cutoff=None):
    """
    Calculate coordination number for each atom.

    Returns
    -------
    coord_numbers : ndarray
    avg_coord : float
    """
    from ase.neighborlist import NeighborList

    n_atoms = len(atoms)

    if cutoff is None:
        cutoff = 3.5

    nl = NeighborList([cutoff / 2] * n_atoms, self_interaction=False,
                      bothways=True, skin=0.0)
    nl.update(atoms)

    coord_numbers = np.array([len(nl.get_neighbors(i)[0]) for i in range(n_atoms)])
    avg_coord = np.mean(coord_numbers)

    return coord_numbers, avg_coord


def calculate_q4_q6(atoms, cutoff=None):
    """
    Calculate Q4 and Q6 bond-orientational order parameters.

    Returns
    -------
    q4 : float
    q6 : float
    q4_per_atom : ndarray
    q6_per_atom : ndarray
    """
    from ase.neighborlist import NeighborList
    try:
        from scipy.special import sph_harm_y as sph_harm
    except ImportError:
        from scipy.special import sph_harm

    n_atoms = len(atoms)

    if cutoff is None:
        cutoff = 3.5

    positions = atoms.get_positions()
    cell = atoms.get_cell()

    nl = NeighborList([cutoff / 2] * n_atoms, self_interaction=False,
                      bothways=True, skin=0.0)
    nl.update(atoms)

    q4_per_atom = np.zeros(n_atoms)
    q6_per_atom = np.zeros(n_atoms)

    m4_array = np.arange(-4, 5)
    m6_array = np.arange(-6, 7)

    for i in range(n_atoms):
        indices, offsets = nl.get_neighbors(i)
        if len(indices) == 0:
            continue

        pos_i = positions[i]
        pos_js = positions[indices] + np.dot(offsets, cell)
        vecs = pos_js - pos_i
        dists = np.linalg.norm(vecs, axis=1, keepdims=True)
        valid = dists.flatten() > 1e-10
        if not np.any(valid):
            continue

        vecs = vecs[valid]
        dists = dists[valid]
        bond_vectors = vecs / dists

        theta = np.arccos(bond_vectors[:, 2])
        phi = np.arctan2(bond_vectors[:, 1], bond_vectors[:, 0])

        nb = len(theta)
        q4_m = []
        for m in m4_array:
            ylm = sph_harm(m, 4, phi, theta)
            q4_m.append(np.sum(ylm) / nb)
        q4_atom = np.sqrt(4 * np.pi / 9 * np.sum(np.abs(np.array(q4_m)) ** 2))
        q4_per_atom[i] = q4_atom

        q6_m = []
        for m in m6_array:
            ylm = sph_harm(m, 6, phi, theta)
            q6_m.append(np.sum(ylm) / nb)
        q6_atom = np.sqrt(4 * np.pi / 13 * np.sum(np.abs(np.array(q6_m)) ** 2))
        q6_per_atom[i] = q6_atom

    q4 = np.mean(q4_per_atom)
    q6 = np.mean(q6_per_atom)

    return q4, q6, q4_per_atom, q6_per_atom


def diagnose_structure_split_z(atoms, cutoff=None, z_mid=None):
    """
    Diagnose solid/liquid phase by splitting along z-axis.

    Parameters
    ----------
    atoms : ase.Atoms
    cutoff : float or None
        Cutoff radius for Q4/Q6, auto-detected from RDF if None
    z_mid : float or None
        Midpoint z for splitting; auto-detected as mean if None

    Returns
    -------
    result : dict
        Keys: 'upper_phase', 'lower_phase', 'upper_q4', 'lower_q4', ...
    """
    positions = atoms.get_positions()
    cell = atoms.get_cell()

    if z_mid is None:
        z_mid = np.mean(positions[:, 2])

    upper_mask = positions[:, 2] > z_mid
    lower_mask = positions[:, 2] <= z_mid

    if cutoff is None:
        r, g_r = calculate_rdf(atoms)
        cutoff, _ = find_first_minimum_rdf(r, g_r)

    result = {'cutoff': cutoff, 'z_mid': z_mid}

    for label, mask in [('upper', upper_mask), ('lower', lower_mask)]:
        if np.sum(mask) < 4:
            result[f'{label}_phase'] = 'unknown'
            continue

        indices = np.where(mask)[0]
        sub_atoms = atoms[indices]
        sub_atoms.set_cell(cell)
        sub_atoms.set_pbc(atoms.get_pbc())

        q4, q6, q4pa, q6pa = calculate_q4_q6(sub_atoms, cutoff)
        result[f'{label}_q4'] = q4
        result[f'{label}_q6'] = q6
        result[f'{label}_q4_per_atom'] = q4pa
        result[f'{label}_q6_per_atom'] = q6pa

        if q4 > 0.1 and q6 > 0.3:
            result[f'{label}_phase'] = 'solid (FCC-like)'
        elif q4 > 0.05:
            result[f'{label}_phase'] = 'partial order'
        else:
            result[f'{label}_phase'] = 'liquid'

    return result


def _gaussian(x, A, x0, sigma):
    return A * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))


def find_first_peak_rdf_gaussian(r, g_r, r_min=0.5, r_max=5.0):
    """
    Fit Gaussian to the first RDF peak.

    Returns
    -------
    result : dict
        Keys: 'peak_r', 'peak_g', 'FWHM', 'sigma', 'fit_r2'
    """
    from scipy.optimize import curve_fit

    if r is None or g_r is None or len(r) == 0:
        return None

    mask = (r >= r_min) & (r <= r_max)
    r_fit = r[mask]
    g_fit = g_r[mask]

    if len(r_fit) < 5:
        return None

    peak_idx = np.argmax(g_fit)
    peak_r_guess = r_fit[peak_idx]
    peak_g_guess = g_fit[peak_idx]

    try:
        popt, _ = curve_fit(_gaussian, r_fit, g_fit,
                             p0=[peak_g_guess, peak_r_guess, 0.3],
                             maxfev=10000)
        A, x0, sigma = popt
        g_pred = _gaussian(r_fit, *popt)
        ss_res = np.sum((g_fit - g_pred) ** 2)
        ss_tot = np.sum((g_fit - np.mean(g_fit)) ** 2)
        r2 = 1 - (ss_res / ss_tot) if ss_tot > 0 else 0.0

        FWHM = 2 * np.sqrt(2 * np.log(2)) * sigma

        return {
            'peak_r': x0,
            'peak_g': A,
            'FWHM': FWHM,
            'sigma': sigma,
            'fit_r2': r2
        }
    except Exception as e:
        return {
            'peak_r': peak_r_guess,
            'peak_g': peak_g_guess,
            'FWHM': None,
            'sigma': None,
            'fit_r2': None,
            'error': str(e)
        }
