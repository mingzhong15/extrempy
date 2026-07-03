"""Surface slab generation from bulk crystal structures.

A thin, pure-ASE re-implementation of the slab-cutting logic from
``dpgen/data/surf.py``.  No pymatgen dependency.

Public API
----------
make_slab   — cut a single surface slab from bulk Atoms
make_slabs  — batch-generate slabs for multiple Miller indices
"""

import os

import numpy as np
from ase.build import general_surface, make_supercell
from ase.io import read, write


def _to_surface_supercell_matrix(supercell):
    """Convert a supercell spec to a 3×3 integer matrix.

    For surfaces an *int* means in-plane *n*×*n* expansion (z unchanged),
    so that the slab thickness is controlled solely by ``layers``/``z_min``.
    """
    if isinstance(supercell, int):
        return np.diag([supercell, supercell, 1])
    if isinstance(supercell, (tuple, list)):
        if len(supercell) == 2:
            return np.diag([supercell[0], supercell[1], 1])
        if len(supercell) == 3:
            return np.diag(supercell)
    raise ValueError(
        f"supercell must be int or tuple of length 2/3, got {supercell}"
    )


def make_slab(bulk_atoms, miller,
              layers=None, z_min=None,
              vacuum_min=None, supercell=(1, 1, 1)):
    """Cut a surface slab from bulk atoms.

    Parameters
    ----------
    bulk_atoms : ase.Atoms
        Bulk crystal structure.
    miller : list[int] or tuple[int]
        Miller index [h, k, l].
    layers : int, optional
        Number of atomic layers in the slab.  Mutually exclusive with
        *z_min*.
    z_min : float, optional
        Minimum slab thickness in Angstrom (excluding vacuum).
        Mutually exclusive with *layers*.  Automatically increments the
        layer count from 1 up to 50 until the slab reaches this thickness.
    vacuum_min : float, optional
        Vacuum thickness in Angstrom.  Defaults to 2× the largest atomic
        radius in the structure.
    supercell : int or tuple[int], optional
        Supercell expansion.  *int* means in-plane *n*×*n* (z unchanged);
        *(nx, ny)* means in-plane only; *(nx, ny, nz)* gives full control
        (nz stacks slabs along the surface normal).

    Returns
    -------
    ase.Atoms
        Surface slab structure.
    """
    if layers is None and z_min is None:
        raise ValueError("must specify either layers or z_min")
    if layers is not None and z_min is not None:
        raise ValueError("layers and z_min are mutually exclusive")

    if vacuum_min is None:
        from ase.data import covalent_radii
        vacuum_min = 2 * max(covalent_radii[bulk_atoms.numbers])

    if layers is not None:
        slab = general_surface.surface(
            bulk_atoms, indices=miller, vacuum=vacuum_min, layers=layers
        )
    else:
        slab = None
        for n in range(1, 51):
            s = general_surface.surface(
                bulk_atoms, indices=miller, vacuum=vacuum_min, layers=n
            )
            if s.cell.lengths()[-1] >= z_min:
                slab = s
                break
        if slab is None:
            raise ValueError(
                f"could not reach z_min={z_min} within 50 layers "
                f"(miller={miller})"
            )

    smat = _to_surface_supercell_matrix(supercell)
    if not np.array_equal(smat, np.eye(3, dtype=int)):
        slab = make_supercell(slab, smat)

    return slab


def _resolve_source(source):
    """Resolve *source* to an ``ase.Atoms`` instance.

    Accepts an ``Atoms`` object, a file path (POSCAR/CIF), or a
    callable that returns ``Atoms``.
    """
    if callable(source):
        return source()
    if hasattr(source, "get_positions"):
        return source
    if isinstance(source, str):
        if source.lower().endswith(".cif"):
            return read(source, format="cif")
        return read(source, format="vasp")
    raise TypeError(
        f"unsupported source type: {type(source)} "
        "(expected ase.Atoms, file path str, or callable)"
    )


def make_slabs(bulk_source, millers, out_dir=None,
               layers=None, z_min=None, vacuum_min=None,
               supercell=(1, 1, 1)):
    """Batch-generate surface slabs for multiple Miller indices.

    Parameters
    ----------
    bulk_source : ase.Atoms / str / callable
        Bulk crystal structure.  Can be an ``Atoms`` object, a
        POSCAR/CIF file path, or a callable returning ``Atoms``.
    millers : list[list[int]]
        List of Miller indices, e.g. ``[[1,0,0], [1,1,1]]``.
    out_dir : str, optional
        Output directory.  If given, writes
        ``{out_dir}/surf-{hkl}.POSCAR`` for each surface.  If ``None``,
        returns ``Atoms`` objects instead.
    layers : int, optional
        Number of slab layers.  Mutually exclusive with *z_min*.
    z_min : float, optional
        Minimum slab thickness in Angstrom.  Mutually exclusive with
        *layers*.
    vacuum_min : float, optional
        Vacuum thickness in Angstrom.
    supercell : int or tuple[int], optional
        Supercell expansion (see :func:`make_slab`).

    Returns
    -------
    dict
        If *out_dir* is given: ``{label: path_str}``.
        Otherwise: ``{label: ase.Atoms}``.
    """
    bulk = _resolve_source(bulk_source)

    results = {}
    for miller in millers:
        slab = make_slab(
            bulk, miller,
            layers=layers, z_min=z_min,
            vacuum_min=vacuum_min, supercell=supercell,
        )
        hkl = "".join(str(i) for i in miller)
        label = f"surf-{hkl}"

        if out_dir is not None:
            os.makedirs(out_dir, exist_ok=True)
            path = os.path.join(out_dir, f"{label}.POSCAR")
            write(path, slab, format="vasp", direct=True)
            results[label] = path
        else:
            results[label] = slab

    return results
