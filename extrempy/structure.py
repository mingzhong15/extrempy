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


# ──────────────────────────────────────────────
#  Public API – existing
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


# ──────────────────────────────────────────────
#  Public API – new
# ──────────────────────────────────────────────

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

    For every phase segment (from ``get_phase_segments`` or passed
    explicitly via *segs*):

    1. Skip LIQ phases (they are populated later from AIMD CONTCAR).
    2. If the target POSCAR already exists and *force* is False → skip.
    3. Try ASE ``_generate_atoms`` for standard structure types
       (fcc, bcc, hcp, dhcp, diamond, sc, bct).
    4. Try each entry in *structure_sources* in order:

       - ``str`` → treated as a directory path; looked up via
         ``_file_source_resolve`` (supports ``.POSCAR``, ``.cif``,
         ``.vasp``).
       - ``callable(element, label, structure_type)`` → must return
         an ``ase.Atoms`` object, a file path (``str``), or ``None``.

    Returns a dict::

        {label: {'path': ..., 'n_atoms': ..., 'source': ...,
                 'status': 'ok'|'skipped'|'skipped_liq'|'error'}}

    Parameters
    ----------
    element : str
        Element symbol (e.g. ``'Ga'``, ``'Bi'``).
    work_root : str
        Root working directory. The confs directory is
        ``{work_root}/{element}/confs/``.
    structure_sources : list[str | callable], optional
        Additional structure sources tried after ASE generation fails
        (or for structure types not supported by ASE). Each entry:

        - ``str`` – path to a directory containing
          ``{label}.POSCAR`` / ``.cif`` / ``.vasp`` files.
        - ``callable(element, label, structure_type)`` – returns
          ``ase.Atoms``, a file path, or ``None``.
    segs : list[dict], optional
        Phase segment list. If ``None``, auto-detected via
        ``get_phase_segments(..., skip_unsupported=False)``.
    supercell : int or (int, int, int), optional
        Supercell expansion applied to structures obtained from
        *structure_sources* (ASE-generated structures already have
        supercell applied internally).  ``None`` → no expansion.
    force : bool
        Overwrite existing POSCAR files (default ``False``).
    """
    confs_dir = os.path.join(work_root, element, 'confs')
    os.makedirs(confs_dir, exist_ok=True)

    if segs is None:
        segs = get_phase_segments(element, skip_unsupported=False)

    results = {}
    for seg in segs:
        label = seg['label']
        st = seg.get('structure')
        target = os.path.join(confs_dir, f'{label}.POSCAR')

        if label.endswith('-LIQ'):
            print(f"  - {label}: LIQ (placeholder, from AIMD CONTCAR)")
            results[label] = {'status': 'skipped_liq'}
            continue

        if os.path.exists(target) and not force:
            try:
                n_atoms = len(read(target, format='vasp'))
            except Exception:
                n_atoms = 0
            print(f"  - {label}: exists ({n_atoms} atoms, skip)")
            results[label] = {
                'path': target, 'n_atoms': n_atoms, 'status': 'skipped'
            }
            continue

        atoms = None
        source_desc = 'none'

        if st in SUPPORTED_STRUCTURES:
            try:
                atoms, n = _generate_atoms(element, st, supercell=None)
                source_desc = 'ASE'
                print(f"  - {label}: ASE ({n} atoms)")
            except Exception as e:
                print(f"  - {label}: ASE failed ({e})")

        if atoms is None and structure_sources:
            for src_idx, src in enumerate(structure_sources):
                try:
                    if isinstance(src, str):
                        result = _file_source_resolve(src, label)
                        if result is not None:
                            atoms = result
                            source_desc = f'file[{src}]'
                    elif callable(src):
                        result = src(element, label, st)
                        if result is not None:
                            if isinstance(result, str):
                                atoms = read(result, format='vasp')
                            elif hasattr(result, 'get_positions'):
                                atoms = result
                            else:
                                continue
                            source_desc = f'callable[{src_idx}]'
                    if atoms is not None:
                        break
                except Exception as e:
                    print(f"  - {label}: source[{src_idx}] error ({e})")
                    continue

        if atoms is None:
            msg = (
                f"No structure available for '{label}' (type='{st}').\n"
                f"  Place a POSCAR manually at:\n"
                f"    {target}\n"
                f"  or use prepare_confs({element!r}, work_root=...,\n"
                f"    structure_sources=[...]) with a suitable source."
            )
            print(f"  - {label}: {msg}")
            results[label] = {'status': 'error', 'error': msg}
            continue

        if supercell is not None and source_desc != 'ASE':
            atoms = make_supercell(atoms, _to_supercell_matrix(supercell))

        write(target, atoms, format='vasp', direct=True)
        n_atoms = len(atoms)
        print(f"  \u2713 {label}: {n_atoms} atoms -> {target} ({source_desc})")
        results[label] = {
            'path': target, 'n_atoms': n_atoms,
            'source': source_desc, 'status': 'ok',
        }

    return results


# ──────────────────────────────────────────────
#  Legacy batch helpers (unchanged)
# ──────────────────────────────────────────────

def batch_generate_structures(elements, output_dir="structures",
                              target_atoms=100, structure_type=None):
    """Generate POSCAR files for multiple elements."""
    results = {}
    for element in elements:
        try:
            atoms, poscar_path = generate_element_structure(
                element, output_dir, target_atoms,
                structure_type=structure_type, verbose=True
            )
            results[element] = {
                'atoms': atoms, 'path': poscar_path, 'n_atoms': len(atoms)
            }
        except Exception as e:
            print(f"Error generating structure for {element}: {e}")
            results[element] = None
    return results


def generate_all_typical_elements(output_dir="structures", target_atoms=100):
    """Generate POSCAR files for all typical metallic elements."""
    results = {}
    total_count = 0
    success_count = 0

    elements_to_generate = {}
    for sym, data in ELEMENT_PHASE_DATA.items():
        rt = data.get('rt_structure')
        if rt and rt in ('fcc', 'bcc', 'hcp', 'diamond', 'sc', 'bct', 'dhcp'):
            if rt not in elements_to_generate:
                elements_to_generate[rt] = []
            elements_to_generate[rt].append(sym)

    print(
        f"Generating POSCAR files for typical metallic elements "
        f"(target: {target_atoms} atoms)..."
    )
    print("=" * 70)

    for struct_type, elements_list in elements_to_generate.items():
        print(f"\n{struct_type.upper()} structures:")
        for element in elements_list:
            total_count += 1
            try:
                atoms, poscar_path = generate_element_structure(
                    element, output_dir=output_dir,
                    target_atoms=target_atoms,
                    structure_type=struct_type, verbose=True
                )
                results[element] = {
                    'structure': struct_type,
                    'atoms': atoms, 'path': poscar_path,
                    'n_atoms': len(atoms),
                }
                success_count += 1
            except Exception as e:
                print(f"  ERROR: Failed to generate {element} ({struct_type}): {e}")
                results[element] = None

    print("\n" + "=" * 70)
    print(f"Generation complete: {success_count}/{total_count} successful")
    print(f"Output directory: {output_dir}")
