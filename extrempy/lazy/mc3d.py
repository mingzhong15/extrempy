import re
import requests
import numpy as np
from io import StringIO
from ase.io import read
from ase import Atoms

from extrempy.lazy.minerals import parse_formula as _parse_formula

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

    # Explicitly sort by energy_per_atom ascending (lowest energy first);
    # None values sort last.  Makes the output contract self-evident and
    # independent of get_phases' internal ordering.
    phases.sort(key=lambda p: (p["energy_per_atom"] is None,
                                p["energy_per_atom"] or 0))

    n = len(phases)
    segs = []
    for i, p in enumerate(phases):
        # Sanitize spg_intl for filename safety: remove spaces, replace
        # '/' with '-' (e.g. 'C 222' -> 'C222', 'P21/c' -> 'P21-c').
        spg_safe = p['spg_intl'].replace(' ', '').replace('/', '-')
        segs.append({
            "label": f"{element}-SG{p['sg']}-{spg_safe}",
            "short_name": spg_safe,
            "structure": "mc3d",
            "T_core": (Tm * i / n, Tm * (i + 1) / n),
            "T_explore": (T_min, T_max_factor * Tm),
            "mc3d_id": p["id"],
            "structure_uuid": p["structure_uuid"],
            "sg": p["sg"],
            "spg_intl": p["spg_intl"],
            "n_atoms_cell": p["n_atoms_cell"],
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
        "short_name": "LIQ",
        "structure": "mc3d",
        "T_core": (Tm, T_max_factor * Tm),
        "T_explore": (Tm, T_max_factor * Tm),
    })
    return segs


# ================================================================
#  Compound query (multi-element, includes HP/HT phases)
# ================================================================

def get_phases_compound(formula, method="pbesol-v2", mode="all"):
    """Query MC3D for all polymorphs of a compound formula.

    Unlike :func:`get_phases` (single element), this matches by exact
    formula composition and includes high-pressure / high-temperature /
    theoretical phases by default (mode='all').

    Parameters
    ----------
    formula : str
        Chemical formula (e.g. 'MgO', 'MgSiO3', 'Mg2SiO4').
    method : str
        MC3D method tag.
    mode : str
        'ambient' — skip theoretical, high_pressure, high_temperature.
        'experimental' — skip theoretical only.
        'all' — everything (default).

    Returns
    -------
    list[dict]
        Same fields as :func:`get_phases`, plus ``formula`` and ``phase_type``.
    """
    target_composition = _parse_formula(formula)
    overview = _get_overview(method)

    # Match by composition: same elements, same counts (order-independent)
    entries = []
    for e in overview:
        n_atoms, elems = _parse_atoms(e["formula"])
        if elems != set(target_composition.keys()):
            continue
        # Verify counts match (Hill notation may differ)
        try:
            e_counts = _parse_formula(e["formula"])
        except ValueError:
            continue
        if e_counts != target_composition:
            continue
        entries.append(e)

    if not entries:
        raise ValueError(
            f"No MC3D entries for formula '{formula}' in method={method}")

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
            "formula": formula,
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


def list_phases_compound(formula, method="pbesol-v2", mode="all"):
    """Print a formatted table of MC3D phases for a compound formula.

    Returns list[dict], same as :func:`get_phases_compound`.
    """
    phases = get_phases_compound(formula, method=method, mode=mode)

    print()
    print(f"─── MC3D Phases: {formula}  "
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


def backfill_mc3d_uuid(db, method="pbesol-v2"):
    """Fill missing mc3d_uuid by matching (formula, spacegroup) against MC3D.

    For each phase in *db* with ``structure_source='mc3d'``,
    ``mc3d_uuid=None``, and a non-null ``spacegroup``:

      1. Query MC3D by the compound's formula → all polymorph candidates.
      2. Filter candidates by ``spacegroup_number == phase.spacegroup``.
      3. Disposition:
         - **1 match** → fill ``mc3d_uuid`` (recorded in ``matched``).
         - **N matches** → leave ``mc3d_uuid=None``; record all N in
           ``multiple_candidates`` for manual review.  We do NOT
           auto-select, to avoid silent mismatches between polymorphs
           that share a spacegroup.
         - **0 matches** → mark ``structure_source='manual'`` (so the
           phase is no longer claimed as MC3D-sourced); record in
           ``no_match``.

    Phases without ``spacegroup`` (LIQ phases, or any mc3d phase missing
    the field) are skipped and listed in ``skipped``.

    Parameters
    ----------
    db : dict
        Loaded mineral database (modified in place).
    method : str
        MC3D method tag (default 'pbesol-v2').

    Returns
    -------
    dict
        ``{
          'matched': [label, ...],
          'multiple_candidates': {label: [{uuid, sg, spg_intl, energy_per_atom, mc3d_id}, ...]},
          'no_match': [label, ...],
          'skipped': [label, ...],
          'api_errors': [label, ...],
        }``
    """
    report = {
        'matched': [],
        'multiple_candidates': {},
        'no_match': [],
        'skipped': [],
        'api_errors': [],
    }

    # Cache MC3D results per formula (avoid re-querying for multiple
    # phases of the same compound).
    mc3d_cache = {}

    for formula, compound in db['compounds'].items():
        for ph in compound['phases']:
            label = ph.get('label', '<no-label>')
            is_liq = label.endswith('-LIQ')

            # Only act on mc3d-sourced solid phases missing uuid.
            if ph.get('structure_source') != 'mc3d':
                continue
            if ph.get('mc3d_uuid'):
                continue  # already filled
            if is_liq:
                continue  # LIQ has no crystal structure

            sg = ph.get('spacegroup')
            if sg is None:
                report['skipped'].append(
                    f"{formula}/{label}: no spacegroup field")
                continue

            # Query MC3D for this formula (cached).
            if formula not in mc3d_cache:
                try:
                    mc3d_cache[formula] = get_phases_compound(
                        formula, method=method, mode='all')
                except Exception as e:
                    mc3d_cache[formula] = None
                    report['api_errors'].append(
                        f"{formula}/{label}: MC3D query failed ({e})")
                    continue

            mc3d_phases = mc3d_cache[formula]
            if mc3d_phases is None:
                report['api_errors'].append(
                    f"{formula}/{label}: MC3D query failed (cached None)")
                continue

            # Filter candidates by spacegroup.
            candidates = [p for p in mc3d_phases if p['sg'] == sg]

            if len(candidates) == 1:
                ph['mc3d_uuid'] = candidates[0]['structure_uuid']
                report['matched'].append(label)
            elif len(candidates) > 1:
                # Multiple MC3D entries with same sg → record for review.
                # Do NOT auto-fill; manual selection needed.
                report['multiple_candidates'][label] = [
                    {
                        'uuid': c['structure_uuid'],
                        'mc3d_id': c['id'],
                        'sg': c['sg'],
                        'spg_intl': c['spg_intl'],
                        'energy_per_atom': c['energy_per_atom'],
                        'phase_type': c['phase_type'],
                    }
                    for c in candidates
                ]
            else:
                # No MC3D entry with this sg → mark as manual.
                ph['structure_source'] = 'manual'
                report['no_match'].append(
                    f"{formula}/{label}: sg={sg} not in MC3D "
                    f"(marked structure_source='manual')")

    return report
