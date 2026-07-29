"""Two-phase melt diagnosis from LAMMPS ``fix ave/chunk`` output.

This module parses ``chunk.profile`` (layered density / temperature / Q4 / Q6
along z) produced by the ``two-phase.j2`` template, and diagnoses each
candidate temperature as ``solid`` / ``liquid`` / ``coexist`` / ``unknown``
by comparing the upper half (liquid candidate) vs lower half (solid
candidate) of the simulation box.

The parser is data-driven: it reads the LAMMPS header line to build a
column-name → index map, so changes to column order or added columns
(e.g. extra computes) do not break it.

Cross-validation with ``rdf_top.txt`` is available via
:func:`cross_validate_with_rdf`, which does **not** mutate the chunk
verdict — it only annotates a ``confidence`` flag in the details.
"""

import os

import numpy as np


# ---- tunable thresholds (module-level constants) ------------------------

Q4_SOLID = 0.10          # Q4 above this → solid-like
Q6_SOLID = 0.30          # Q6 above this → solid-like (FCC theory ~0.57)
Q6_LIQUID = 0.15         # Q6 above this → not liquid-like (half of Q6_SOLID);
                          # used jointly with Q4_SOLID to confirm solid order
DELTA_Q4 = 0.05          # upper/lower Q4 difference threshold for stratification
DELTA_RHO = 0.02         # relative density difference threshold
RDF_PEAK_THRESHOLD = 1.5
RDF_SOLID_MIN_PEAKS = 3
RDF_R_MIN, RDF_R_MAX = 2.5, 6.0


# ---- parsing ------------------------------------------------------------

def read_chunk_profile(path):
    """Parse a LAMMPS ``fix ave/chunk`` output file.

    The file contains repeated blocks, each preceded by a ``#`` header
    line listing the column names (e.g.
    ``# Chunk Coord1 Ncount density/mass temp v_virial_atom c_Q[2] c_Q[3]``)
    and a timestep line.  This function is data-driven: it reads the
    header to build a column-name → index map, so column reordering or
    additions do not break it.

    Parameters
    ----------
    path : str

    Returns
    -------
    blocks : list[dict]
        Each dict has ``timestep`` (int or None) and ``chunks``
        (list[dict] keyed by the normalised column names —
        ``c_Q[2]`` is mapped to ``q4``, ``c_Q[3]`` to ``q6``,
        ``density/mass`` to ``density``).
        Returns ``[]`` if the file is missing or unparseable.
    """
    if not os.path.exists(path):
        return []

    with open(path) as f:
        lines = f.readlines()

    blocks = []
    cur_header = None
    cur_timestep = None
    cur_rows = []

    for raw in lines:
        line = raw.strip()
        if not line:
            continue

        if line.startswith('#'):
            # Header line: "# Chunk Coord1 Ncount density/mass temp ..."
            tokens = line.lstrip('#').split()
            # Recognise the column-definition line by the presence of
            # canonical chunk/coord keywords (LAMMPS version-agnostic).
            # Older LAMMPS may lead with "Chunk", newer ones with
            # "Coordination" or "Coord1" — match any of these.
            if tokens and any(t.lower() in ('chunk', 'coord1',
                                            'coord', 'ncount')
                              for t in tokens[:3]):
                cur_header = _normalise_header(tokens)
            # else: ignore other comment lines (e.g. file metadata)
            continue

        # Non-comment line: could be a timestep marker or data row.
        parts = line.split()
        # Timestep marker: a single integer on its own line.
        if len(parts) == 1:
            try:
                ts = int(parts[0])
            except ValueError:
                continue
            # Flush previous block
            if cur_rows and cur_header is not None:
                blocks.append(_build_block(cur_timestep, cur_header, cur_rows))
            cur_timestep = ts
            cur_rows = []
            continue

        # Data row
        if cur_header is not None:
            cur_rows.append(parts)

    # Flush final block
    if cur_rows and cur_header is not None:
        blocks.append(_build_block(cur_timestep, cur_header, cur_rows))

    return blocks


def _normalise_header(tokens):
    """Map raw LAMMPS column names to canonical keys.

    Returns a list of canonical names aligned with ``tokens``.
    Unknown names are kept as-is.
    """
    mapping = {
        'chunk': 'id',
        'coordination': 'id',
        'coord1': 'coord',
        'ncount': 'ncount',
        'density/mass': 'density',
        'temp': 'temp',
        'v_virial_atom': 'virial',
        'c_q[2]': 'q4',
        'c_q[3]': 'q6',
    }
    return [mapping.get(t.lower(), t) for t in tokens]


def _build_block(timestep, header, rows):
    """Construct a block dict from parsed rows."""
    chunks = []
    for r in rows:
        if len(r) < len(header):
            continue
        d = {}
        for name, val in zip(header, r):
            try:
                d[name] = float(val)
            except ValueError:
                d[name] = val
        # Preserve integer id if present
        if 'id' in d:
            try:
                d['id'] = int(d['id'])
            except (ValueError, TypeError):
                pass
        chunks.append(d)
    return {'timestep': timestep, 'chunks': chunks}


# ---- per-block diagnosis ------------------------------------------------

def diagnose_block(block, n_half=None):
    """Diagnose a single chunk block by splitting it into upper/lower halves.

    Parameters
    ----------
    block : dict
        Output of :func:`read_chunk_profile` element.
    n_half : int or None
        Number of chunks per half.  Defaults to ``len(chunks) // 2`` so
        any layer count works.

    Returns
    -------
    dict
        Keys: ``lower_q4``, ``upper_q4``, ``lower_q6``, ``upper_q6``,
        ``lower_rho``, ``upper_rho``, ``delta_q4``, ``delta_rho``,
        ``verdict``.
        ``verdict`` is ``'solid'`` / ``'liquid'`` / ``'coexist'`` /
        ``'unknown'``.
    """
    chunks = block.get('chunks', [])
    if not chunks:
        return _empty_verdict()

    if n_half is None:
        n_half = len(chunks) // 2
    if n_half < 1:
        return _empty_verdict()

    lower = chunks[:n_half]
    upper = chunks[n_half:n_half * 2] or chunks[n_half:]

    lower_q4 = np.mean([c.get('q4', 0.0) for c in lower])
    upper_q4 = np.mean([c.get('q4', 0.0) for c in upper])
    lower_q6 = np.mean([c.get('q6', 0.0) for c in lower])
    upper_q6 = np.mean([c.get('q6', 0.0) for c in upper])
    lower_rho = np.mean([c.get('density', 0.0) for c in lower])
    upper_rho = np.mean([c.get('density', 0.0) for c in upper])

    delta_q4 = upper_q4 - lower_q4
    # Relative density difference (guard against zero).
    rho_ref = lower_rho if lower_rho else 1.0
    delta_rho = abs(upper_rho - lower_rho) / rho_ref

    lower_solid = (lower_q4 > Q4_SOLID) and (lower_q6 > Q6_LIQUID)
    upper_solid = (upper_q4 > Q4_SOLID) and (upper_q6 > Q6_LIQUID)
    stratified = (abs(delta_q4) > DELTA_Q4) and (delta_rho > DELTA_RHO)

    if lower_solid and not upper_solid and stratified:
        verdict = 'coexist'
    elif lower_solid and upper_solid:
        verdict = 'solid'
    elif not lower_solid and not upper_solid:
        verdict = 'liquid'
    else:
        verdict = 'unknown'

    return dict(
        lower_q4=float(lower_q4), upper_q4=float(upper_q4),
        lower_q6=float(lower_q6), upper_q6=float(upper_q6),
        lower_rho=float(lower_rho), upper_rho=float(upper_rho),
        delta_q4=float(delta_q4), delta_rho=float(delta_rho),
        verdict=verdict,
    )


def _empty_verdict():
    return dict(
        lower_q4=0.0, upper_q4=0.0, lower_q6=0.0, upper_q6=0.0,
        lower_rho=0.0, upper_rho=0.0, delta_q4=0.0, delta_rho=0.0,
        verdict='unknown',
    )


# ---- RDF cross-validation ----------------------------------------------

def diagnose_rdf(r, g_r, r_min=RDF_R_MIN, r_max=RDF_R_MAX,
                 peak_threshold=RDF_PEAK_THRESHOLD,
                 min_peaks=RDF_SOLID_MIN_PEAKS):
    """Classify an RDF curve as ``'solid'`` / ``'liquid'`` / ``'unknown'``.

    A solid shows multiple pronounced peaks (long-range order); a
    liquid shows only a first peak followed by dampened oscillations.

    Parameters
    ----------
    r, g_r : array-like
    r_min, r_max : float
        Range to search for secondary peaks (Å).
    peak_threshold : float
        Minimum ``g(r)`` value to count as a peak.
    min_peaks : int
        Number of peaks above ``peak_threshold`` required for ``'solid'``.
    """
    if r is None or g_r is None:
        return 'unknown'
    r = np.asarray(r)
    g_r = np.asarray(g_r)
    if len(r) == 0:
        return 'unknown'

    mask = (r >= r_min) & (r <= r_max)
    if not np.any(mask):
        return 'unknown'
    g_search = g_r[mask]
    n_peaks = int(np.sum(g_search > peak_threshold))
    return 'solid' if n_peaks >= min_peaks else 'liquid'


def cross_validate_with_rdf(chunk_verdict, rdf_verdict):
    """Annotate a chunk-based verdict with RDF confidence.

    Does **not** mutate ``chunk_verdict``; only returns a confidence
    flag so the caller can decide how much to trust the result.

    Returns
    -------
    dict
        ``{'rdf_verdict': str, 'confidence': str}`` where confidence is
        ``'high'`` (chunk coexist + RDF liquid), ``'low'`` (chunk
        coexist + RDF solid, i.e. upper half may not have truly melted),
        or ``'medium'`` otherwise.
    """
    if chunk_verdict == 'coexist' and rdf_verdict == 'solid':
        confidence = 'low'
    elif chunk_verdict == 'coexist' and rdf_verdict == 'liquid':
        confidence = 'high'
    else:
        confidence = 'medium'
    return {'rdf_verdict': rdf_verdict, 'confidence': confidence}


# ---- case-level aggregation --------------------------------------------

def diagnose_case(chunk_path, rdf_top_path=None, vote_window=5):
    """Diagnose one temperature directory.

    Reads ``chunk.profile``, votes over the last ``adaptive_window``
    blocks (``adaptive_window = min(vote_window, max(1, n // 4))`` so
    short trajectories do not over-sample), and optionally cross-checks
    with ``rdf_top.txt``.

    Parameters
    ----------
    chunk_path : str
    rdf_top_path : str or None
        Path to ``rdf_top.txt`` (upper-half RDF).  If given, adds
        ``rdf_top_verdict`` and ``confidence`` to the result.
    vote_window : int

    Returns
    -------
    dict
        Keys: ``verdict``, ``lower_q4``, ``upper_q4``, ``lower_q6``,
        ``upper_q6``, ``lower_rho``, ``upper_rho``, ``n_blocks_voted``,
        ``rdf_top_verdict`` (or None), ``confidence`` (or None).
    """
    blocks = read_chunk_profile(chunk_path)
    if not blocks:
        return _empty_case()

    adaptive_window = min(vote_window, max(1, len(blocks) // 4))
    voted = blocks[-adaptive_window:]

    # Diagnose each block ONCE; reuse for both voting and averaging.
    block_diags = [diagnose_block(b) for b in voted]
    verdicts = [d['verdict'] for d in block_diags]
    # Majority vote; tie → unknown.
    counts = {}
    for v in verdicts:
        counts[v] = counts.get(v, 0) + 1
    winner = max(counts, key=counts.get)
    top_count = counts[winner]
    # Tie check
    tied = [v for v, c in counts.items() if c == top_count]
    verdict = winner if len(tied) == 1 else 'unknown'

    # Average diagnostics over voted blocks.
    diag_avg = _average_diag(block_diags)

    result = dict(verdict=verdict, n_blocks_voted=len(voted), **diag_avg)

    # RDF cross-validation (annotation only, does not change verdict).
    if rdf_top_path and os.path.exists(rdf_top_path):
        r, g_r = _read_rdf_file(rdf_top_path)
        rdf_v = diagnose_rdf(r, g_r)
        cv = cross_validate_with_rdf(verdict, rdf_v)
        result['rdf_top_verdict'] = rdf_v
        result['confidence'] = cv['confidence']
    else:
        result['rdf_top_verdict'] = None
        result['confidence'] = None

    return result


def _read_rdf_file(filepath):
    """Minimal LAMMPS fix ave/time RDF parser (self-contained, numpy only).

    Returns ``(r, g_r)`` as 1-D arrays averaged over all time frames, or
    ``(None, None)`` if unparseable.  Kept inline so this module stays
    independent of the heavier ``extrempy.md`` package (which pulls in
    polars / ase / scipy at import time).
    """
    try:
        with open(filepath) as f:
            lines = f.readlines()
    except OSError:
        return None, None

    nbins = 0
    data = []
    for line in lines:
        parts = line.split()
        if len(parts) == 2 and nbins == 0:
            try:
                nbins = int(parts[-1])
            except ValueError:
                continue
        if len(parts) == 4:
            try:
                data.append([float(parts[1]), float(parts[2])])
            except (ValueError, IndexError):
                continue
    if nbins == 0 or not data:
        return None, None
    arr = np.array(data).reshape(-1, nbins, 2)
    averaged = arr.mean(axis=0)
    return averaged[:, 0], averaged[:, 1]


def _average_diag(diags):
    """Average numeric diagnostic fields across multiple block verdicts."""
    keys = ['lower_q4', 'upper_q4', 'lower_q6', 'upper_q6',
            'lower_rho', 'upper_rho', 'delta_q4', 'delta_rho']
    out = {}
    for k in keys:
        vals = [d[k] for d in diags if k in d]
        out[k] = float(np.mean(vals)) if vals else 0.0
    return out


def _empty_case():
    return dict(
        verdict='unknown', n_blocks_voted=0,
        lower_q4=0.0, upper_q4=0.0, lower_q6=0.0, upper_q6=0.0,
        lower_rho=0.0, upper_rho=0.0, delta_q4=0.0, delta_rho=0.0,
        rdf_top_verdict=None, confidence=None,
    )
