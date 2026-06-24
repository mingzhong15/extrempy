import re
import requests
import numpy as np
from io import StringIO
from ase.io import read
from ase import Atoms

MC3D_BASE = "https://mcxd-api.materialscloud.org/mc3d"
AIIDA_BASE = "https://aiida.materialscloud.org"

_overview_cache = {}
_detail_cache = {}


def _parse_atoms(formula):
    elems = re.findall(r"[A-Z][a-z]?", formula)
    counts = re.findall(r"\d+", formula) if re.search(r"\d", formula) else []
    if counts:
        total = sum(int(c) for _, c in zip(elems, counts))
    else:
        total = len(elems)
    return total, set(elems)


def _get_overview(method):
    if method not in _overview_cache:
        url = f"{MC3D_BASE}/{method}/overview"
        _overview_cache[method] = requests.get(url, timeout=120).json()
    return _overview_cache[method]


def _get_detail(entry_id, method):
    key = f"{method}/{entry_id}"
    if key not in _detail_cache:
        url = f"{MC3D_BASE}/{method}/core_base/{entry_id}"
        _detail_cache[key] = requests.get(url, timeout=30).json()
    return _detail_cache[key]


def get_phases(element, method="pbesol-v2", mode="ambient"):
    """
    Query MC3D for all polymorphs of a pure element.

    Parameters
    ----------
    element : str
    method : str
    mode : str
        'ambient' — skip theoretical, high_pressure, high_temperature.
        'experimental' — skip theoretical only.
        'all' — everything.
    """
    overview = _get_overview(method)
    entries = [
        e for e in overview
        if _parse_atoms(e["formula"])[1] == {element}
    ]
    if not entries:
        raise ValueError(f"No entries for {element} in MC3D/{method}")

    if mode == "ambient":
        entries = [e for e in entries
                   if e.get("th") is not True
                   and e.get("hp") is not True
                   and e.get("ht") is not True]
    elif mode == "experimental":
        entries = [e for e in entries if e.get("th") is not True]

    results = []
    for e in entries:
        d = _get_detail(e["id"], method)
        g = d["general"]
        p = d["properties"]
        n_atoms, _ = _parse_atoms(g["formula_hill"])
        te = p.get("total_energy", {}).get("value")
        vol = p.get("cell_volume", {}).get("value")

        if e.get("th") is True:
            ptype = "theoretical"
        elif e.get("hp") is True:
            ptype = "high_pressure"
        elif e.get("ht") is True:
            ptype = "high_temperature"
        else:
            ptype = "ambient"

        results.append({
            "id": e["id"],
            "sg": g["spacegroup_number"],
            "bravais": g["bravais_lattice"],
            "spg_intl": g["spacegroup_international"],
            "energy_per_atom": te / n_atoms if te is not None else None,
            "cell_volume": vol,
            "phase_type": ptype,
            "sdb": e.get("sdb"),
            "source_id": e.get("sid"),
            "n_atoms_cell": n_atoms,
            "structure_uuid": g["structure_uuid"],
        })

    results.sort(key=lambda x: (x["energy_per_atom"] is None,
                                  x["energy_per_atom"] or 0))
    return results


def list_phases(element, method="pbesol-v2", mode="ambient"):
    """
    Print a formatted table of MC3D phases.

    Returns list[dict], same as get_phases().
    """
    phases = get_phases(element, method=method, mode=mode)

    print()
    print(f"─── MC3D Phases: {element}  "
          f"(method={method}, mode={mode}) ───")
    header = (f"  {'#':<4s} {'MC3D ID':<16s} {'SG':>4s}  "
              f"{'Internat.':<10s}  {'eV/atom':>12s}  "
              f"{'Vol (A^3)':>10s}  {'Type':<16s}  {'Src':>5s}")
    print(header)
    print(f"  {'─' * 89}")

    for i, p in enumerate(phases):
        e_pa = (f"{p['energy_per_atom']:12.4f}"
                if p["energy_per_atom"] is not None else "           --")
        print(f"  {i+1:<4d} {p['id']:<16s} {p['sg']:4d}  "
              f"{p['spg_intl']:<10s}  {e_pa}  "
              f"{p['cell_volume']:10.2f}  "
              f"{p['phase_type']:<16s}  {p['sdb']:>5s}")

    print(f"  {'─' * 89}")
    n = len(phases)
    if n:
        best = phases[0]
        print(f"  {n} phase(s).  Most stable: "
              f"#{1} {best['id']} ({best['spg_intl']}, "
              f"{best['energy_per_atom']:.4f} eV/atom)")
    print()
    return phases


def download_atoms(structure_uuid, method="pbesol-v2"):
    url = (f"{AIIDA_BASE}/mc3d-{method}/api/v4/nodes/"
           f"{structure_uuid}/download?download_format=cif")
    resp = requests.get(url, timeout=60)
    resp.raise_for_status()
    return read(StringIO(resp.text), format="cif")


def make_phase_segments(element, method="pbesol-v2", mode="ambient",
                        Tm=None, T_min=300, T_max_factor=2):
    """
    Generate DPBuilder-compatible phase segments from MC3D.

    Solid phases evenly split Tm range; a LIQ placeholder is appended.
    """
    from .lib import ELEMENT_PHASE_DATA

    if Tm is None:
        Tm = ELEMENT_PHASE_DATA.get(element, {}).get("Tm", 1000)

    phases = get_phases(element, method=method, mode=mode)
    if not phases:
        raise ValueError(f"No MC3D phases for {element} (mode={mode})")

    n = len(phases)
    segs = []
    for i, p in enumerate(phases):
        segs.append({
            "label": f"{element}-{p['sg']}",
            "structure": "mc3d",
            "T_core": (Tm * i / n, Tm * (i + 1) / n),
            "T_explore": (T_min, T_max_factor * Tm),
            "mc3d_id": p["id"],
            "structure_uuid": p["structure_uuid"],
            "sg": p["sg"],
            "energy_per_atom": p["energy_per_atom"],
            "phase_type": p["phase_type"],
        })

    # deduplicate labels matching lib.py's _deduplicate_labels
    seen = {}
    for s in segs:
        lbl = s["label"]
        if lbl not in seen:
            seen[lbl] = 0
        else:
            seen[lbl] += 1
            s["label"] = f"{lbl}-{seen[lbl] + 1}"

    segs.append({
        "label": f"{element}-LIQ",
        "structure": "mc3d",
        "T_core": (Tm, T_max_factor * Tm),
        "T_explore": (Tm, T_max_factor * Tm),
    })
    return segs
