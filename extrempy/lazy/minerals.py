"""Deep-Earth mineral thermodynamic database.

Loads `extrempy/data/mineral_phases.json` and provides:
  - constants: EARTH_LAYERS, DEEP_EARTH_ELEMENTS
  - query: list_compounds, get_compound, get_compound_phases
  - segment construction: make_extreme_segments (2D P-T aware)
  - validation: validate_db

Database schema (per phase, 8 fields):
  label, spacegroup, spg_intl, mc3d_uuid, structure_source,
  P_range_GPa, T_range_K, refs

LIQ phases have spacegroup=null and spg_intl=null.

Derived quantities (NOT stored):
  P_explore = P_range +/- p_margin (default 0.15)
  T_explore = T_range +/- t_margin (default 300 K)
  earth_layers = lookup from P_range against EARTH_LAYERS
"""

import json
import os
import re

try:
    from importlib.resources import files as _resource_files
except ImportError:
    _resource_files = None


# ================================================================
#  Constants
# ================================================================

EARTH_LAYERS = {
    'upper_mantle':    {'P_GPa': (0, 13),   'T_K': (300, 2000)},
    'transition_zone': {'P_GPa': (13, 23),  'T_K': (500, 2200)},
    'lower_mantle':    {'P_GPa': (23, 125), 'T_K': (1500, 4000)},
    'd_prime':         {'P_GPa': (125, 135),'T_K': (2500, 4000)},
    'outer_core':      {'P_GPa': (135, 330),'T_K': (4000, 6000)},
    'inner_core':      {'P_GPa': (330, 400),'T_K': (5000, 10000)},
}

DEEP_EARTH_ELEMENTS = [
    'Mg', 'Si', 'O', 'Fe', 'Ca', 'Ni', 'H', 'C',
    'N', 'Al', 'He', 'Na', 'P', 'S', 'K', 'Ti',
]

GLOBAL_P_RANGE_GPa = (0, 400)
GLOBAL_T_RANGE_K = (300, 10000)

DEFAULT_P_MARGIN = 0.15
DEFAULT_T_MARGIN = 300


# ================================================================
#  Database loading
# ================================================================

def _default_db_path():
    """Locate the packaged mineral_phases.json via importlib.resources."""
    if _resource_files is not None:
        try:
            return str(_resource_files('extrempy.data') / 'mineral_phases.json')
        except (ModuleNotFoundError, FileNotFoundError, TypeError):
            pass
    # Fallback: relative to this file
    here = os.path.dirname(os.path.abspath(__file__))
    return os.path.normpath(os.path.join(here, '..', 'data', 'mineral_phases.json'))


def load_mineral_db(json_path=None, *, reload=False):
    """Load the mineral thermodynamic database from JSON.

    Always reads the file fresh — no caching. For a 10 KB JSON file the
    IO cost is negligible (<1 ms), and always-fresh reads make manual
    edits during a session immediately visible.

    The ``reload`` parameter is kept only for backward compatibility
    with older callers/tests and has no effect.

    Returns
    -------
    dict
        Mirrors `mineral_phases.json` top-level structure.
    """
    path = json_path or _default_db_path()
    if not os.path.exists(path):
        raise FileNotFoundError(
            f"mineral_phases.json not found at {path}. "
            f"If running from source, ensure extrempy/data/ exists.")
    with open(path, 'r', encoding='utf-8') as f:
        db = json.load(f)
    if 'compounds' not in db:
        raise KeyError(f"JSON at {path} missing top-level 'compounds' key")
    return db


# ================================================================
#  Formula parsing
# ================================================================

def parse_formula(formula):
    """Parse a chemical formula into {element: count}.

    Examples:
        'Mg' -> {'Mg': 1}
        'MgO' -> {'Mg': 1, 'O': 1}
        'MgSiO3' -> {'Mg': 1, 'Si': 1, 'O': 3}
        'Mg2SiO4' -> {'Mg': 2, 'Si': 1, 'O': 4}
        'NaAlSi3O8' -> {'Na': 1, 'Al': 1, 'Si': 3, 'O': 8}

    Returns
    -------
    dict[str, int]
    """
    # Pattern: element symbol (capital + optional lowercase) + optional count
    tokens = re.findall(r'([A-Z][a-z]?)(\d*)', formula)
    result = {}
    consumed = 0
    for elem, count_str in tokens:
        consumed += len(elem) + len(count_str)
        count = int(count_str) if count_str else 1
        result[elem] = result.get(elem, 0) + count
    if consumed != len(formula):
        raise ValueError(
            f"Failed to parse formula '{formula}': "
            f"consumed {consumed} of {len(formula)} chars")
    if not result:
        raise ValueError(f"Empty formula: '{formula}'")
    return result


def formula_elements(formula):
    """Return sorted unique element list of a formula."""
    return sorted(parse_formula(formula).keys())


# ================================================================
#  Query
# ================================================================

def list_compounds(db=None):
    """List all compound formulas in the database."""
    if db is None:
        db = load_mineral_db()
    return sorted(db['compounds'].keys())


def get_compound(formula, db=None):
    """Get a full compound entry: {Tm_at_1bar_K, phases: [...]}.

    Raises KeyError if formula not in database.
    """
    if db is None:
        db = load_mineral_db()
    if formula not in db['compounds']:
        raise KeyError(
            f"Compound '{formula}' not in mineral database. "
            f"Available: {list_compounds(db)}")
    return db['compounds'][formula]


def get_compound_phases(formula, db=None):
    """Get the phases list for a compound."""
    return get_compound(formula, db)['phases']


# ================================================================
#  Earth-layer lookup (derived)
# ================================================================

def lookup_earth_layers(P_range_GPa):
    """Return earth-layer names whose P range overlaps with the given range."""
    p_lo, p_hi = P_range_GPa
    matched = []
    for name, info in EARTH_LAYERS.items():
        layer_lo, layer_hi = info['P_GPa']
        if p_lo < layer_hi and p_hi > layer_lo:
            matched.append(name)
    return matched


# ================================================================
#  Segment construction (2D P-T aware)
# ================================================================

def _derive_explore(P_range, T_range, p_margin, t_margin):
    """Compute P_explore and T_explore from ranges + margins."""
    p_lo, p_hi = P_range
    p_span = p_hi - p_lo
    p_explore = (max(0, p_lo - p_margin * p_span),
                 p_hi + p_margin * p_span)

    t_lo, t_hi = T_range
    t_span = t_hi - t_lo
    t_explore = (max(300, t_lo - t_margin),
                 t_hi + t_margin)
    return p_explore, t_explore


def make_extreme_segments(formula, db=None, *,
                          p_margin=DEFAULT_P_MARGIN,
                          t_margin=DEFAULT_T_MARGIN):
    """Build 2D (P,T) segments for a compound from the mineral database.

    Each returned seg dict contains:
        label, spacegroup, spg_intl, mc3d_uuid, structure_source,
        T_core, T_explore, P_core, P_explore (all in GPa / K),
        earth_layers, refs

    Parameters
    ----------
    formula : str
        Compound formula (e.g. 'MgO', 'MgSiO3', or pure element 'Fe').
    db : dict or None
        Loaded database; None loads the default JSON.
    p_margin : float
        Fractional expansion of P_range to get P_explore (default 0.15).
    t_margin : float
        Absolute expansion (in K) of T_range to get T_explore (default 300).

    Returns
    -------
    list[dict]
    """
    phases = get_compound_phases(formula, db)
    segs = []
    for ph in phases:
        P_range = ph['P_range_GPa']
        T_range = ph['T_range_K']
        P_explore, T_explore = _derive_explore(P_range, T_range,
                                               p_margin, t_margin)
        seg = {
            'label': ph['label'],
            'spacegroup': ph.get('spacegroup'),
            'spg_intl': ph.get('spg_intl'),
            'mc3d_uuid': ph.get('mc3d_uuid'),
            'structure_source': ph.get('structure_source', 'mc3d'),
            'P_core': tuple(P_range),
            'T_core': tuple(T_range),
            'P_explore': P_explore,
            'T_explore': T_explore,
            'earth_layers': lookup_earth_layers(P_range),
            'refs': ph.get('refs', []),
        }
        segs.append(seg)
    return segs


# ================================================================
#  Validation
# ================================================================

def validate_db(db=None):
    """Return a list of warnings about incomplete or malformed entries.

    Empty list = no warnings.

    Checks performed:
      - Compound has Tm_at_1bar_K and at least one phase
      - No duplicate phase labels within a compound
      - Every phase has P_range_GPa and T_range_K (well-formed)
      - Every non-LIQ phase with structure_source='mc3d' has spacegroup
        (needed for backfill_mc3d_uuid matching)
      - Every non-LIQ mc3d phase has mc3d_uuid (or warns to run backfill)
      - Every phase has at least one ref (confidence low otherwise)
    """
    if db is None:
        db = load_mineral_db()
    warnings = []
    for formula, compound in db['compounds'].items():
        if 'Tm_at_1bar_K' not in compound:
            warnings.append(f"{formula}: missing Tm_at_1bar_K")
        if 'phases' not in compound or not compound['phases']:
            warnings.append(f"{formula}: no phases")
            continue
        labels = [ph.get('label') for ph in compound['phases']]
        if len(labels) != len(set(labels)):
            warnings.append(f"{formula}: duplicate phase labels {labels}")
        for ph in compound['phases']:
            lbl = ph.get('label', '<no-label>')
            is_liq = lbl.endswith('-LIQ')
            for field in ('P_range_GPa', 'T_range_K'):
                if field not in ph:
                    warnings.append(f"{formula}/{lbl}: missing {field}")
                elif not (isinstance(ph[field], list)
                          and len(ph[field]) == 2):
                    warnings.append(
                        f"{formula}/{lbl}: malformed {field}={ph.get(field)}")
            if not ph.get('refs'):
                warnings.append(f"{formula}/{lbl}: no refs (confidence low)")
            src = ph.get('structure_source', 'mc3d')
            # Non-LIQ mc3d phases need spacegroup for backfill matching
            if src == 'mc3d' and not is_liq:
                if ph.get('spacegroup') is None:
                    warnings.append(
                        f"{formula}/{lbl}: structure_source='mc3d' but "
                        f"spacegroup is null (needed for backfill_mc3d_uuid)")
                if not ph.get('mc3d_uuid'):
                    warnings.append(
                        f"{formula}/{lbl}: structure_source='mc3d' but "
                        f"mc3d_uuid is null (run backfill_mc3d_uuid)")
    return warnings


# ================================================================
#  Printing helpers
# ================================================================

def print_compound_summary(formula, db=None):
    """Pretty-print a compound's phases with P-T ranges and refs."""
    segs = make_extreme_segments(formula, db)
    print(f"\n─── {formula}: {len(segs)} phases ───")
    hdr = (f"  {'#':<3} {'Label':<28} {'P (GPa)':<16} "
           f"{'T (K)':<16} {'Layers':<20} {'Refs'}")
    print(hdr)
    print(f"  {'-' * 90}")
    for i, s in enumerate(segs):
        p_str = f"{s['P_core'][0]:.0f}-{s['P_core'][1]:.0f}"
        t_str = f"{s['T_core'][0]:.0f}-{s['T_core'][1]:.0f}"
        layers = ','.join(s['earth_layers']) if s['earth_layers'] else '-'
        refs = '; '.join(s['refs']) if s['refs'] else '(none)'
        print(f"  {i:<3} {s['label']:<28} {p_str:<16} "
              f"{t_str:<16} {layers:<20} {refs}")
    print()
