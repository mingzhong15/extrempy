import numpy as np
import os

from ase.build import bulk, make_supercell
from ase.io import write, read
from ase import Atoms

from .lazy.lib import (_get_lattice_from_data, ELEMENT_PHASE_DATA,
                       ELEMENT_PRIMITIVE_ATOMS, LATTICE_CONSTANTS,
                       SUPPORTED_STRUCTURES, get_phase_segments)


# ──────────────────────────────────────────────
#  Internal helpers
# ──────────────────────────────────────────────

def _generate_atoms(element, structure_type, supercell=None, target_atoms=100):
    """Generate ASE Atoms for standard structure types without writing files.

    Returns (atoms, n_atoms).
    """
    if supercell is None:
        use_supercell = None
    elif isinstance(supercell, int):
        use_supercell = (supercell, supercell, supercell)
    elif isinstance(supercell, (tuple, list)) and len(supercell) == 3:
        use_supercell = tuple(supercell)
    else:
        raise ValueError(
            f"supercell must be int or tuple of length 3, got {supercell}"
        )

    lattice_key = (element, structure_type)
    has_custom_lattice = lattice_key in LATTICE_CONSTANTS

    if structure_type == 'dhcp':
        if has_custom_lattice:
            a, c4 = LATTICE_CONSTANTS[lattice_key]
            c = c4 / 2
        else:
            rt_phase = ELEMENT_PHASE_DATA.get(element, {}).get('phases', [{}])[0]
            a = rt_phase.get('a', 3.5)
            c4 = rt_phase.get('c', 11.8)
            c = c4 / 2

        atoms = bulk(element, 'hcp', a=a, c=c)

        supercell_matrix = np.eye(3, dtype=int)
        supercell_matrix[2, 2] = 2
        atoms = make_supercell(atoms, supercell_matrix)

        pos = atoms.get_positions()
        cell = atoms.get_cell()
        frac = pos @ np.linalg.inv(cell)
        atoms.set_positions(frac @ cell)

        if use_supercell is not None:
            nx, ny, nz = use_supercell
        else:
            _prim = 4
            n = int(np.floor((target_atoms / _prim) ** (1 / 3)))
            n = max(1, n)
            nx = ny = n
            nz = max(1, n // 2)
        s_mat = np.diag([nx, ny, nz])
        atoms = make_supercell(atoms, s_mat)

    elif structure_type in ('fcc', 'bcc', 'diamond', 'sc'):
        if has_custom_lattice:
            a = LATTICE_CONSTANTS[lattice_key][0]
            atoms = bulk(element, structure_type, a=a, cubic=True)
        else:
            a = _get_lattice_from_data(element, structure_type)[0]
            atoms = bulk(element, structure_type, a=a, cubic=True)
    elif structure_type == 'hcp':
        if has_custom_lattice:
            a, c = LATTICE_CONSTANTS[lattice_key]
            atoms = bulk(element, 'hcp', a=a, c=c)
        else:
            a, c = _get_lattice_from_data(element, "hcp")
            atoms = bulk(element, "hcp", a=a, c=c)
    elif structure_type == 'bct':
        if has_custom_lattice:
            a, c = LATTICE_CONSTANTS[lattice_key]
        else:
            a, c = 3.252, 4.946
        positions = np.array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]])
        cell = np.array([[a, 0.0, 0.0], [0.0, a, 0.0], [0.0, 0.0, c]])
        cart_positions = positions @ cell
        atoms = Atoms(element * 2, positions=cart_positions, cell=cell, pbc=True)
    else:
        raise ValueError(f"Unsupported structure type: {structure_type}")

    if use_supercell is not None and structure_type != 'dhcp':
        _supercell_matrix = np.diag(use_supercell)
        result = make_supercell(atoms, _supercell_matrix)
    else:
        _prim = ELEMENT_PRIMITIVE_ATOMS.get(structure_type, 2)
        n = int(np.floor((target_atoms / _prim) ** (1 / 3)))
        n = max(1, n)
        _supercell_matrix = np.eye(3, dtype=int) * n
        result = make_supercell(atoms, _supercell_matrix)
        n_atoms = len(result)
        if n_atoms > target_atoms and structure_type in ('fcc', 'bcc', 'diamond', 'sc', 'bct'):
            n = int(np.floor((target_atoms / _prim) ** (1 / 3)))
            n = max(1, n)
            _supercell_matrix = np.eye(3, dtype=int) * n
            result = make_supercell(atoms, _supercell_matrix)

    return result, len(result)


def _to_supercell_matrix(supercell):
    if isinstance(supercell, int):
        return np.diag([supercell, supercell, supercell])
    if isinstance(supercell, (tuple, list)) and len(supercell) == 3:
        return np.diag(supercell)
    raise ValueError(
        f"supercell must be int or tuple of length 3, got {supercell}"
    )


def calculate_supercell(n_atoms_cell, target_atoms=100):
    """
    Calculate optimal isotropic supercell (nx, ny, nz).

    Parameters
    ----------
    n_atoms_cell : int
        Number of atoms in the primitive cell.
    target_atoms : int
        Target total number of atoms.

    Returns
    -------
    (int, int, int)
    """
    ratio = target_atoms / n_atoms_cell
    n = max(1, int(np.round(ratio ** (1 / 3))))
    best, best_diff = None, float("inf")
    for nx in range(max(1, n - 1), n + 2):
        for ny in range(max(1, n - 1), n + 2):
            for nz in range(max(1, n - 1), n + 2):
                n_total = nx * ny * nz * n_atoms_cell
                diff = abs(n_total - target_atoms)
                if diff < best_diff:
                    best_diff, best = diff, (nx, ny, nz)
    return best


# Default supercell sizes per standard structure type (used by ase_source
# when neither `supercell` nor `target_atoms` forces a different choice).
DEFAULT_SUPERCELL = {
    'fcc': (2, 2, 2),
    'bcc': (3, 3, 3),
    'hcp': (3, 3, 4),
    'diamond': (2, 2, 2),
    'dhcp': (3, 3, 2),
    'sc': (3, 3, 3),
    'bct': (3, 3, 3),
}


# ──────────────────────────────────────────────
#  Structure sources — factories returning () -> ase.Atoms
# ──────────────────────────────────────────────

def ase_source(element, structure_type=None, supercell=None, target_atoms=100):
    """Return a callable that builds an ase.Atoms via ASE ``bulk()``.

    When *structure_type* is None, auto-detects from
    ``ELEMENT_PHASE_DATA[element]['rt_structure']``.
    """
    def _make():
        st = structure_type
        if st is None:
            rt = ELEMENT_PHASE_DATA.get(element, {}).get('rt_structure')
            st = rt or 'fcc'
        sc = supercell or DEFAULT_SUPERCELL.get(st)
        atoms, _ = _generate_atoms(element, st, supercell=sc,
                                    target_atoms=target_atoms)
        return atoms
    return _make


def mc3d_source(structure_uuid, supercell=None, target_atoms=100,
                method='pbesol-v2'):
    """Return a callable that downloads a structure from MC3D."""
    def _make():
        from .lazy.mc3d import download_atoms
        atoms = download_atoms(structure_uuid, method=method)
        if supercell is not None:
            atoms = make_supercell(atoms, _to_supercell_matrix(supercell))
        else:
            sc = calculate_supercell(len(atoms), target_atoms)
            atoms = make_supercell(atoms, np.diag(sc))
        return atoms
    return _make


# ──────────────────────────────────────────────
#  Unified resolver
# ──────────────────────────────────────────────

def resolve_poscar(element, label, *, confs_dir, source,
                   force=False, verbose=True):
    """Resolve a structure to ``{confs_dir}/{label}.POSCAR``.

    This is the single entry point used by both :class:`DPBuilder` and
    :class:`EOSCalculator`.  It handles only the "already exists /
    LIQ placeholder / write file" concerns; the actual structure
    production is delegated to *source*.

    Parameters
    ----------
    element : str
    label : str
        Output filename stem.  Labels ending in ``'-LIQ'`` are treated
        as liquid placeholders and skipped (returns ``None``).
    confs_dir : str
    source : callable or ase.Atoms or str
        - callable: ``() -> ase.Atoms`` (e.g. from
          :func:`ase_source` / :func:`mc3d_source`)
        - ``ase.Atoms``: used directly
        - ``str``: path to an existing POSCAR/CIF/vasp file, read in
    force : bool
        Overwrite an existing POSCAR at the target path.
    verbose : bool

    Returns
    -------
    str or None
        Path to the POSCAR, or ``None`` for LIQ placeholders.
    """
    out_path = os.path.join(confs_dir, f'{label}.POSCAR')

    if label.endswith('-LIQ'):
        if verbose:
            print(f'  - {label}: LIQ (placeholder, from AIMD CONTCAR)')
        return None

    if os.path.exists(out_path) and not force:
        if verbose:
            try:
                n = len(read(out_path, format='vasp'))
            except Exception:
                n = 0
            print(f'  - {label}: exists ({n} atoms, skip)')
        return out_path

    # Obtain atoms from source.
    if callable(source):
        atoms = source()
    elif hasattr(source, 'get_positions'):
        atoms = source
    elif isinstance(source, str):
        if not os.path.exists(source):
            raise FileNotFoundError(f'source path not found: {source}')
        atoms = read(source, format='vasp')
    else:
        raise TypeError(
            f'source must be callable, ase.Atoms, or path str; '
            f'got {type(source)}')

    os.makedirs(confs_dir, exist_ok=True)
    write(out_path, atoms, format='vasp', direct=True)
    if verbose:
        print(f'  \u2713 {label}: {len(atoms)} atoms -> {out_path}')
    return out_path


def _file_source_resolve(lib_path, label):
    """Try to load a structure from a directory.

    Looks for {label}.POSCAR, {label}.cif, or {label}.vasp (in that order).
    Returns ase.Atoms or None.
    """
    for ext, fmt in [('.POSCAR', 'vasp'), ('.cif', 'cif'), ('.vasp', 'vasp')]:
        path = os.path.join(lib_path, f'{label}{ext}')
        if os.path.exists(path):
            return read(path, format=fmt)
    return None


def prepare_confs(element, work_root,
                  structure_sources=None,
                  segs=None,
                  supercell=None,
                  force=False):
    """Prepare POSCAR files for all phases of *element* into
    ``{work_root}/{element}/confs/``.

    Thin batch wrapper around :func:`resolve_poscar`.  For every
    phase segment (from :func:`get_phase_segments` or passed
    explicitly via *segs*) a structure source is constructed and
    passed to ``resolve_poscar``.

    Source selection per seg:

    1. If ``seg['structure']`` is in :data:`SUPPORTED_STRUCTURES` →
       :func:`ase_source`.
    2. Else if ``structure_sources`` is given, try each in order:

       - ``str`` → directory; looked up via
         :func:`_file_source_resolve` (supports ``.POSCAR`` /
         ``.cif`` / ``.vasp``).  The resolved file path is passed
         directly to ``resolve_poscar``.
       - ``callable(element, label, structure_type)`` → must return
         an ``ase.Atoms``, a file path (``str``), or ``None``.

    Returns
    -------
    dict
        ``{label: {'path': ..., 'status': 'ok'|'skipped'|'skipped_liq'|'error', ...}}``

    Parameters
    ----------
    element : str
    work_root : str
    structure_sources : list[str | callable], optional
    segs : list[dict], optional
        If ``None``, auto-detected via
        ``get_phase_segments(..., skip_unsupported=False)``.
    supercell : int or (int, int, int), optional
        Applied to ASE-generated structures (passed through to
        :func:`ase_source`).
    force : bool
        Overwrite existing POSCAR files.
    """
    confs_dir = os.path.join(work_root, element, 'confs')
    if segs is None:
        segs = get_phase_segments(element, skip_unsupported=False)

    results = {}
    for seg in segs:
        label = seg['label']
        st = seg.get('structure')
        target = os.path.join(confs_dir, f'{label}.POSCAR')

        # Build the source for this seg.
        if st in SUPPORTED_STRUCTURES:
            src = ase_source(element, structure_type=st,
                             supercell=supercell)
        elif structure_sources:
            src = None
            for s in structure_sources:
                try:
                    if isinstance(s, str):
                        found = _file_source_resolve(s, label)
                        if found is not None:
                            src = found   # an Atoms object, passed as-is
                    elif callable(s):
                        r = s(element, label, st)
                        if r is not None:
                            if isinstance(r, str):
                                src = r    # path str
                            elif hasattr(r, 'get_positions'):
                                src = r    # Atoms
                    if src is not None:
                        break
                except Exception as e:
                    print(f'  - {label}: source error ({e})')
            if src is None:
                print(f'  - {label}: no source resolved')
                results[label] = {'status': 'error',
                                  'error': 'no source resolved'}
                continue
            # Apply supercell to file/callable sources (ASE source already
            # applied it internally via _generate_atoms).
            if supercell is not None and hasattr(src, 'get_positions'):
                src = make_supercell(src, _to_supercell_matrix(supercell))
        else:
            print(f'  - {label}: structure type {st!r} not supported '
                  f'and no structure_sources given')
            results[label] = {'status': 'error',
                              'error': f'unsupported structure {st!r}'}
            continue

        try:
            path = resolve_poscar(element, label, confs_dir=confs_dir,
                                  source=src, force=force, verbose=True)
            if path is None:
                results[label] = {'status': 'skipped_liq'}
            elif path == target and not force and os.path.exists(target):
                # resolve_poscar printed "exists"; we still record n_atoms
                try:
                    n = len(read(path, format='vasp'))
                except Exception:
                    n = 0
                results[label] = {'path': path, 'n_atoms': n,
                                  'status': 'skipped'}
            else:
                try:
                    n = len(read(path, format='vasp'))
                except Exception:
                    n = 0
                results[label] = {'path': path, 'n_atoms': n,
                                  'status': 'ok'}
        except Exception as e:
            print(f'  - {label}: resolve_poscar failed ({e})')
            results[label] = {'status': 'error', 'error': str(e)}

    return results


# ──────────────────────────────────────────────
#  Legacy single-element helper (kept as a simpler alternative to
#  resolve_poscar when batch/segment logic is not needed)
# ──────────────────────────────────────────────

def generate_element_structure(element,
                                output_dir="structures",
                                target_atoms=100,
                                supercell=None,
                                structure_type=None,
                                verbose=True):
    """
    Generate POSCAR file for element with automatic supercell expansion.

    Parameters
    ----------
    element : str
    output_dir : str
    target_atoms : int
        Target number of atoms (used when supercell is None).
    supercell : int or tuple or None
        Supercell dimensions. If int, use as n x n x n.
        If tuple (nx, ny, nz), use as anisotropic supercell.
        If None, compute from target_atoms.
    structure_type : str or None
        One of: 'fcc', 'bcc', 'hcp', 'diamond', 'sc', 'bct', 'dhcp'
        If None, auto-detect from ELEMENT_PHASE_DATA rt_structure
    verbose : bool

    Returns
    -------
    atoms : ase.Atoms
    poscar_path : str
    """
    if structure_type is None:
        rt = ELEMENT_PHASE_DATA.get(element, {}).get('rt_structure')
        if rt:
            structure_type = rt
        else:
            structure_type = 'fcc'

    atoms, n_atoms = _generate_atoms(
        element, structure_type,
        supercell=supercell, target_atoms=target_atoms
    )

    os.makedirs(output_dir, exist_ok=True)
    poscar_path = os.path.join(
        output_dir, f"{element}-{structure_type.upper()}.POSCAR"
    )
    write(poscar_path, atoms, format='vasp', direct=True)

    if verbose:
        print(
            f"Generated {element} ({structure_type.upper()}): "
            f"{n_atoms} atoms -> {poscar_path}"
        )

    return atoms, poscar_path
