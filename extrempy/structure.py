import numpy as np
import os

from ase.build import bulk, make_supercell
from ase.io import write
from ase import Atoms

from .lazy.lib import (_get_lattice_from_data, ELEMENT_PHASE_DATA,
                       ELEMENT_PRIMITIVE_ATOMS, LATTICE_CONSTANTS)

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

    if supercell is None:
        use_supercell = None
    elif isinstance(supercell, int):
        use_supercell = (supercell, supercell, supercell)
    elif isinstance(supercell, (tuple, list)) and len(supercell) == 3:
        use_supercell = tuple(supercell)
    else:
        raise ValueError(f"supercell must be int or tuple of length 3, got {supercell}")

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
        for i in range(len(atoms)):
            if 0.4 < frac[i, 2] < 0.6:
                pass
        atoms.set_positions(frac @ cell)

        if use_supercell is not None:
            nx, ny, nz = use_supercell
        else:
            primitive_atoms = 4
            n = int(np.floor((target_atoms / primitive_atoms) ** (1 / 3)))
            n = max(1, n)
            nx = ny = n
            nz = max(1, n // 2)
        supercell_matrix = np.diag([nx, ny, nz])
        atoms = make_supercell(atoms, supercell_matrix)

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
        supercell = make_supercell(atoms, _supercell_matrix)
    else:
        primitive_atoms = ELEMENT_PRIMITIVE_ATOMS.get(structure_type, 2)

        n = int(np.floor((target_atoms / primitive_atoms) ** (1 / 3)))
        n = max(1, n)
        _supercell_matrix = np.eye(3, dtype=int) * n

        supercell = make_supercell(atoms, _supercell_matrix)
        n_atoms = len(supercell)

        if n_atoms > target_atoms and structure_type in ('fcc', 'bcc', 'diamond', 'sc', 'bct'):
            n = int(np.floor((target_atoms / primitive_atoms) ** (1 / 3)))
            n = max(1, n)
            _supercell_matrix = np.eye(3, dtype=int) * n
            supercell = make_supercell(atoms, _supercell_matrix)

    n_atoms = len(supercell)

    os.makedirs(output_dir, exist_ok=True)
    poscar_path = os.path.join(output_dir, f"{element}-{structure_type.upper()}.POSCAR")
    write(poscar_path, supercell, format='vasp', direct=True)

    if verbose:
        print(f"Generated {element} ({structure_type.upper()}): {n_atoms} atoms -> {poscar_path}")

    return supercell, poscar_path


def batch_generate_structures(elements, output_dir="structures", target_atoms=100, structure_type=None):
    """Generate POSCAR files for multiple elements."""
    results = {}
    for element in elements:
        try:
            atoms, poscar_path = generate_element_structure(
                element, output_dir, target_atoms,
                structure_type=structure_type, verbose=True
            )
            results[element] = {'atoms': atoms, 'path': poscar_path, 'n_atoms': len(atoms)}
        except Exception as e:
            print(f"Error generating structure for {element}: {e}")
            results[element] = None
    return results


def generate_all_typical_elements(output_dir="structures", target_atoms=100):
    """Generate POSCAR files for all typical metallic elements."""
    results = {}
    total_count = 0
    success_count = 0

    # Only generate for elements with 'rt_structure' set
    elements_to_generate = {}
    for sym, data in ELEMENT_PHASE_DATA.items():
        rt = data.get('rt_structure')
        if rt and rt in ('fcc', 'bcc', 'hcp', 'diamond', 'sc', 'bct', 'dhcp'):
            if rt not in elements_to_generate:
                elements_to_generate[rt] = []
            elements_to_generate[rt].append(sym)

    print(f"Generating POSCAR files for typical metallic elements (target: {target_atoms} atoms)...")
    print("=" * 70)

    for structure_type, elements_list in elements_to_generate.items():
        print(f"\n{structure_type.upper()} structures:")
        for element in elements_list:
            total_count += 1
            try:
                atoms, poscar_path = generate_element_structure(
                    element, output_dir=output_dir, target_atoms=target_atoms,
                    structure_type=structure_type, verbose=True
                )
                results[element] = {
                    'structure': structure_type,
                    'atoms': atoms, 'path': poscar_path, 'n_atoms': len(atoms)
                }
                success_count += 1
            except Exception as e:
                print(f"  ERROR: Failed to generate {element} ({structure_type}): {e}")
                results[element] = None

    print("\n" + "=" * 70)
    print(f"Generation complete: {success_count}/{total_count} successful")
    print(f"Output directory: {output_dir}")
