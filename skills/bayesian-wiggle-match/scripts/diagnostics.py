"""Annual posterior concentration and candidate peaks; never a uniqueness pass."""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path

import numpy as np

VERSION = '0.3.0'
WINDOWS = (1, 2, 5, 10, 20, 50)


def ranges(values):
    out = []
    for value in values:
        value = int(value)
        if out and value == out[-1][1] + 1:
            out[-1][1] = value
        else:
            out.append([value, value])
    return out


def integer(value, name):
    try:
        v = float(value)
    except (TypeError, ValueError):
        raise ValueError(name + ' must be a finite integer') from None
    if isinstance(value, bool) or not math.isfinite(v) or v != int(v):
        raise ValueError(name + ' must be a finite integer')
    return int(v)


def analyze(years, probability, start_min, start_max):
    """Describe a normalized posterior over assessed integer calendar positions.

    Missing positions never become rejected dates. Peak basin masses use a
    documented bookkeeping partition, not a validated scientific mode classifier.
    """
    lo = integer(start_min, 'start_min'); hi = integer(start_max, 'start_max')
    if hi < lo or hi - lo > 100000:
        raise ValueError('Ordered search bounds may span at most 100001 years')
    y = np.array([integer(v, 'year') for v in years], dtype=np.int64)
    p = np.array(probability, dtype=float)
    if y.ndim != 1 or p.ndim != 1 or len(y) == 0 or len(y) != len(p):
        raise ValueError('Nonempty one-dimensional years and probability must have equal length')
    if np.any(np.diff(y) <= 0) or y[0] < lo or y[-1] > hi:
        raise ValueError('Years must be strictly increasing, unique and within search bounds')
    if np.any(~np.isfinite(p)) or np.any(p < 0) or abs(float(p.sum()) - 1) > 1e-10:
        raise ValueError('Probabilities must be finite, nonnegative and normalized')
    dense = np.zeros(hi-lo+1)
    assessed = np.zeros(len(dense), dtype=np.int64)
    dense[y-lo] = p; assessed[y-lo] = 1
    conditional = len(y) != len(dense)
    windows = {}
    for width in WINDOWS:
        item = dict(width_years=width, conditional_on_reference_coverage=conditional,
                    numerical_tie_absolute_tolerance=1e-12)
        if width > len(dense):
            item.update(status='wider_than_search', conditional_best_mass=None,
                        maximizing_start_ranges=[], fully_assessed_best_mass=None,
                        fully_assessed_maximizing_start_ranges=[])
        else:
            sums = np.convolve(dense, np.ones(width), mode='valid')
            count = np.convolve(assessed, np.ones(width, dtype=np.int64), mode='valid')
            best = float(np.max(sums))
            ties = np.flatnonzero(np.abs(sums-best) <= 1e-12) + lo
            full = np.flatnonzero(count == width)
            full_best = float(np.max(sums[full])) if len(full) else None
            full_ties = full[np.abs(sums[full]-full_best) <= 1e-12] + lo if len(full) else []
            item.update(status='conditional_description' if conditional else 'assessed',
                        conditional_best_mass=min(1., best), maximizing_start_ranges=ranges(ties),
                        fully_assessed_best_mass=min(1., full_best) if full_best is not None else None,
                        fully_assessed_maximizing_start_ranges=ranges(full_ties))
        windows[str(width)] = item
    # Split on unassessed years, not on zero-probability assessed positions.
    cuts = [0] + (np.flatnonzero(np.diff(y) > 1) + 1).tolist() + [len(y)]
    peaks = []; valleys = []
    for block, (begin, end) in enumerate(zip(cuts[:-1], cuts[1:])):
        candidates = []
        i = begin
        while i < end:
            j = i
            while j + 1 < end and p[j+1] == p[i]:
                j += 1
            if p[i] > 0 and (i == begin or p[i] > p[i-1]) and (j+1 == end or p[j] > p[j+1]):
                candidates.append((i, j))
            i = j+1
        boundaries = []
        for (left_start, left_end), (right_start, right_end) in zip(candidates[:-1], candidates[1:]):
            segment = p[left_end+1:right_start]
            minimum = float(segment.min())
            ties = np.flatnonzero(segment == minimum) + left_end+1
            split = (int(ties[0]) + int(ties[-1])) // 2
            boundaries.append(split)
            valleys.append(dict(block=block, left_peak_start=int(y[left_start]),
                                right_peak_start=int(y[right_start]), valley_start=int(y[ties[0]]),
                                valley_end=int(y[ties[-1]]), probability=minimum,
                                split_after_year=int(y[split]),
                                plateau_center_separation_years=float((y[right_start]+y[right_end]-y[left_start]-y[left_end])/2)))
        for k, (a, b) in enumerate(candidates):
            basin_left = begin if k == 0 else boundaries[k-1]+1
            basin_right = end-1 if k == len(candidates)-1 else boundaries[k]
            edge = a == begin or b == end-1
            prominence = None
            if not edge:
                base = []
                for indices in [range(a-1, begin-1, -1), range(b+1, end)]:
                    minimum = float(p[a])
                    for idx in indices:
                        if p[idx] > p[a]:
                            break
                        minimum = min(minimum, float(p[idx]))
                    base.append(minimum)
                prominence = float(p[a] - max(base))
            peaks.append(dict(block=block, plateau_start=int(y[a]), plateau_end=int(y[b]),
                              height=float(p[a]), plateau_mass=float(p[a:b+1].sum()),
                              prominence=prominence,
                              relative_prominence=prominence/float(p[a]) if prominence is not None else None,
                              basin_start=int(y[basin_left]), basin_end=int(y[basin_right]),
                              basin_mass=float(p[basin_left:basin_right+1].sum()),
                              touches_search_boundary=bool(y[a] == lo or y[b] == hi),
                              touches_unassessed_gap=bool((a == begin and y[a] > lo) or (b == end-1 and y[b] < hi))))
    mass = sum(peak['basin_mass'] for peak in peaks)
    if abs(mass - float(p.sum())) > 1e-10:
        raise ValueError('Candidate basin partition did not conserve posterior probability')
    flat = bool(np.all(p == p[0]))
    return dict(diagnostics_version=VERSION, requested_start_min=lo, requested_start_max=hi,
                assessed_positions=len(y), assessed_fraction=len(y)/len(dense),
                conditional_on_reference_coverage=conditional,
                unassessed_start_ranges=ranges(np.flatnonzero(assessed == 0)+lo),
                boundary_five_year_mass=float(p[(y-lo < 5) | (hi-y < 5)].sum()),
                structure_status='reference_incomplete' if conditional else ('flat_on_assessed_support' if flat else 'candidate_modes_reported_requires_scale_judgment'),
                uniqueness_status='not_automatically_evaluated', windows=windows,
                candidate_peak_count=len(peaks), candidate_peaks=peaks, valleys=valleys,
                candidate_basin_mass_sum=mass,
                algorithm=dict(smoothing='none', plateau_equality='exact stored probability equality',
                               gap_policy='separate assessed blocks; unknown dates are not excluded',
                               basin_boundary='midpoint of minimum inter-peak valley plateau, floor; split year assigned left',
                               prominence='peak minus larger of two side minima up to a strictly higher point or block end; null at block edges'),
                interpretation='Candidate local extrema and bookkeeping basin masses, not proven calendar modes, significant peaks, accuracy, or a uniqueness pass. Interpret full posterior, reference coverage, priors and research timescale together.')


def write_outputs(out, result, create=True):
    out = Path(out).resolve()
    if os.name == 'nt' and out.drive.lower() == 'c:':
        raise ValueError('Output on C: is prohibited')
    if create:
        out.mkdir(parents=True, exist_ok=False)
    names = ['diagnostics.json', 'candidate_peaks.csv', 'peak_valleys.csv']
    if any((out/name).exists() for name in names):
        raise ValueError('Diagnostics output already exists; no overwrite permitted')
    (out/names[0]).write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False), encoding='utf-8')
    for name, rows, fallback in [(names[1], result['candidate_peaks'], ['plateau_start','plateau_end']),
                                 (names[2], result['valleys'], ['left_peak_start','right_peak_start'])]:
        with (out/name).open('x', encoding='utf-8', newline='') as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0]) if rows else fallback)
            w.writeheader(); w.writerows(rows)


def analyze_saved(posterior_path, summary_path, out):
    with Path(posterior_path).open(encoding='utf-8-sig', newline='') as f:
        rows = list(csv.DictReader(f))
    summary = json.loads(Path(summary_path).read_text(encoding='utf-8-sig'))
    result = analyze([r['astro_year'] for r in rows], [r['probability'] for r in rows],
                     summary['requested_start_min'], summary['requested_start_max'])
    if result['assessed_positions'] != summary['assessed_positions'] or result['conditional_on_reference_coverage'] != summary['conditional_on_reference_coverage']:
        raise ValueError('Saved summary and posterior coverage disagree')
    supplied = sorted(integer(v, 'unassessed year') for v in summary['unassessed_start_years'])
    expected = [y for a,b in result['unassessed_start_ranges'] for y in range(a,b+1)]
    if supplied != expected:
        raise ValueError('Saved summary unassessed years disagree with posterior')
    result['sources'] = [dict(path=str(Path(p).resolve()), sha256=hashlib.sha256(Path(p).read_bytes()).hexdigest())
                         for p in [posterior_path, summary_path]]
    result['diagnostics_sha256'] = hashlib.sha256(Path(__file__).read_bytes()).hexdigest()
    write_outputs(out, result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--posterior', required=True)
    parser.add_argument('--summary', required=True)
    parser.add_argument('--out', required=True)
    args = parser.parse_args()
    try:
        result = analyze_saved(args.posterior, args.summary, args.out)
    except (ValueError, KeyError, OSError) as e:
        parser.exit(2, 'Diagnostics error: ' + str(e) + '\n')
    print(json.dumps({k:result[k] for k in ['diagnostics_version','candidate_peak_count','structure_status','uniqueness_status']}))


if __name__ == '__main__':
    main()
