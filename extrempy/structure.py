import numpy as np
import os

from ase.build import bulk, make_supercell
from ase.io import write
from ase.data import chemical_symbols, atomic_numbers
from ase import Atoms


ELEMENTS_BY_STRUCTURE = {
    'fcc': ['Au', 'Cu', 'Ag', 'Al', 'Ni', 'Pt', 'Pd', 'Rh', 'Ir', 'Pb', 'Th', 'Ca', 'Sr', 'Yb', 'Ce'],
    'bcc': ['Fe', 'W', 'Mo', 'V', 'Cr', 'Nb', 'Ta', 'Ba', 'K', 'Na', 'Li', 'Cs', 'Eu'],
    'hcp': ['Mg', 'Zn', 'Ti', 'Zr', 'Co', 'Be', 'Cd', 'Re', 'Os', 'Ru', 'Y', 'Sc',
            'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Lu', 'Hf', 'Tl', 'Tc'],
    'diamond': ['Si', 'Ge', 'C', 'Sn'],
    'sc': ['Po'],
    'bct': ['In'],
}

LATTICE_CONSTANTS = {
    ('Sn', 'diamond'): (6.489,),
    ('In', 'bct'): (3.251, 4.947),
}

ELEMENT_PRIMITIVE_ATOMS = {
    'fcc': 4,
    'bcc': 2,
    'hcp': 2,
    'diamond': 8,
    'sc': 1,
    'bct': 2,
}


def generate_element_structure(element,
                               output_dir="structures",
                               target_atoms=100,
                               structure_type=None,
                               verbose=True):
    """
    Generate POSCAR file for element with automatic supercell expansion.

    Parameters
    ----------
    element : str
        Element symbol (e.g., 'Au', 'Si', 'Fe')
    output_dir : str
        Output directory for POSCAR file
    target_atoms : int
        Target number of atoms (default: 100)
    structure_type : str or None
        Structure type: 'fcc', 'bcc', 'hcp', 'diamond', 'sc', 'bct'
        If None, auto-detect from ELEMENTS_BY_STRUCTURE
    verbose : bool
        Print progress information

    Returns
    -------
    atoms : ase.Atoms
    poscar_path : str
    """
    if structure_type is None:
        structure_type = 'fcc'
        for stype, elements_list in ELEMENTS_BY_STRUCTURE.items():
            if element in elements_list:
                structure_type = stype
                break

    lattice_key = (element, structure_type)
    has_custom_lattice = lattice_key in LATTICE_CONSTANTS

    if structure_type in ('fcc', 'bcc', 'diamond', 'sc'):
        if has_custom_lattice:
            a = LATTICE_CONSTANTS[lattice_key][0]
            atoms = bulk(element, structure_type, a=a, cubic=True)
        else:
            atoms = bulk(element, structure_type, cubic=True)
    elif structure_type == 'hcp':
        if has_custom_lattice:
            a, c = LATTICE_CONSTANTS[lattice_key]
            atoms = bulk(element, 'hcp', a=a, c=c)
        else:
            atoms = bulk(element, 'hcp')
    elif structure_type == 'bct':
        if has_custom_lattice:
            a, c = LATTICE_CONSTANTS[lattice_key]
        else:
            a, c = 3.251, 4.947
        positions = np.array([[0.0, 0.0, 0.0], [0.5, 0.5, 0.5]])
        cell = np.array([[a, 0.0, 0.0], [0.0, a, 0.0], [0.0, 0.0, c]])
        cart_positions = positions @ cell
        atoms = Atoms(element * 2, positions=cart_positions, cell=cell, pbc=True)
    else:
        raise ValueError(f"Unsupported structure type: {structure_type}")

    primitive_atoms = ELEMENT_PRIMITIVE_ATOMS.get(structure_type, 2)

    n = int(np.floor((target_atoms / primitive_atoms) ** (1 / 3)))
    n = max(1, n)
    supercell_matrix = np.eye(3, dtype=int) * n

    supercell = make_supercell(atoms, supercell_matrix)
    n_atoms = len(supercell)

    if n_atoms > target_atoms and structure_type in ('fcc', 'bcc', 'diamond', 'sc', 'bct'):
        n = int(np.floor((target_atoms / primitive_atoms) ** (1 / 3)))
        n = max(1, n)
        supercell_matrix = np.eye(3, dtype=int) * n
        supercell = make_supercell(atoms, supercell_matrix)
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
    """Generate POSCAR files for all typical elements with their standard structures."""
    results = {}
    total_count = 0
    success_count = 0

    print(f"Generating POSCAR files for all typical elements (target: {target_atoms} atoms)...")
    print("=" * 70)

    for structure_type, elements_list in ELEMENTS_BY_STRUCTURE.items():
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
    return results
