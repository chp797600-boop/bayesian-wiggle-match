"""Compare explicit Gaussian offset priors without selecting a winning model."""
import argparse
import hashlib
import itertools
import json
import math
import os
from pathlib import Path

try:
    from .core import match, number
    from .prepare import read_csv, write_csv
    from .run import execute, tokens
except ImportError:
    from core import match, number
    from prepare import read_csv, write_csv
    from run import execute, tokens


def evaluate(samples_path,reference_path,out,start_min,start_max,offset_sds,offset_mean=0,
             allow_incomplete=False,truth_start=None):
    out=Path(out).resolve()
    if out.exists():raise ValueError('Output directory already exists; choose a new directory')
    if os.name=='nt' and out.drive.lower()=='c:':raise ValueError('Output on C: is prohibited')
    priors=[number(x,'offset_sd') for x in offset_sds]
    if len(priors)<2 or len(set(priors))!=len(priors):
        raise ValueError('Supply at least two distinct explicit offset SD values')
    if any(x<0 or not math.isfinite(x*x) for x in priors):
        raise ValueError('Offset SD must be nonnegative with finite variance')
    if truth_start is not None:truth_start=number(truth_start,'truth_start',integer=True)
    samples=read_csv(samples_path);curve=read_csv(reference_path)
    if tokens(samples,'measurement_ids') & tokens(curve,'member_ids') or tokens(samples,'tree_id') & tokens(curve,'tree_ids'):
        raise ValueError('Reference/target overlap in declared measurement IDs or tree IDs')
    # Validate all models before creating a directory. Truth never enters match().
    for sd in priors:
        check=match(samples,curve,start_min,start_max,offset_mean,sd)
        if check['conditional_on_reference_coverage'] and not allow_incomplete:
            raise ValueError('Reference coverage is incomplete; explicit allow_incomplete is required')
    inputs=[dict(path=str(Path(p).resolve()),sha256=hashlib.sha256(Path(p).read_bytes()).hexdigest())
            for p in [samples_path,reference_path]]
    out.mkdir(parents=True,exist_ok=False)
    rows=[]; distributions=[]
    for sd in priors:
        label='offset_sd_'+format(sd,'.15g')
        # Distinct finite floats may share a short decimal representation.
        if (out/label).exists():label='offset_sd_'+repr(sd)
        r=execute(samples_path,reference_path,out/label,start_min,start_max,offset_mean,sd,allow_incomplete)
        intervals=r['hpd']['0.95']['intervals']
        probability={x['astro_year']:x['probability'] for x in r['posterior']}
        status='not_supplied';covered=None
        if truth_start is not None:
            if not start_min<=truth_start<=start_max:status='outside_search_window'
            elif truth_start not in probability:status='unassessed'
            else:
                covered=any(a<=truth_start<=b for a,b in intervals)
                status='covered' if covered else 'missed'
        rows.append(dict(offset_sd_bp=sd,map_start_astro_year=r['map_start_astro_year'],
                         hpd95_intervals=json.dumps(intervals),hpd95_year_count=sum(b-a+1 for a,b in intervals),
                         hpd95_interval_count=len(intervals),hpd95_actual_mass=r['hpd']['0.95']['actual_mass'],
                         assessed_positions=r['assessed_positions'],
                         assessed_fraction=r['assessed_positions']/(start_max-start_min+1),
                         conditional_on_reference_coverage=r['conditional_on_reference_coverage'],
                         boundary_five_year_mass=r['boundary_five_year_mass'],
                         truth_status=status,truth_in_95_hpd=covered,run_directory=label))
        distributions.append(probability)
    for record in inputs:
        if hashlib.sha256(Path(record['path']).read_bytes()).hexdigest()!=record['sha256']:
            raise ValueError('An input changed during sensitivity evaluation; do not interpret this run')
    distances=[]
    for a,b in itertools.combinations(range(len(rows)),2):
        pa=distributions[a];pb=distributions[b]
        if pa.keys()!=pb.keys():raise ValueError('Assessable years changed across priors')
        tv=.5*sum(abs(pa[t]-pb[t]) for t in pa)
        distances.append(dict(offset_sd_a_bp=priors[a],offset_sd_b_bp=priors[b],total_variation=tv))
    report=dict(tool_version='0.4.0',inputs=inputs,
                entrypoint_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                parameters=dict(start_min=start_min,start_max=start_max,offset_mean_bp=offset_mean,
                                offset_sds_bp=priors,allow_incomplete=allow_incomplete),
                truth_start_astro_year=truth_start,truth_usage='Post-hoc scoring only; never used in model fitting',
                runs=rows,pairwise_total_variation=distances,
                max_pairwise_total_variation=max(x['total_variation'] for x in distances),
                map_span_years=max(x['map_start_astro_year'] for x in rows)-min(x['map_start_astro_year'] for x in rows),
                interpretation='Sensitivity on the same assessable positions; no best prior, accuracy estimate, or model averaging selected. A broad or fragmented HPD containing truth is not single-year recovery.')
    write_csv(out/'sensitivity.csv',rows)
    (out/'sensitivity.json').write_text(json.dumps(report,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--samples',required=True);p.add_argument('--reference',required=True)
    p.add_argument('--out',required=True);p.add_argument('--start-min',required=True,type=int)
    p.add_argument('--start-max',required=True,type=int)
    p.add_argument('--offset-sds',required=True,nargs='+',type=float)
    p.add_argument('--offset-mean',type=float,default=0)
    p.add_argument('--allow-incomplete',action='store_true')
    p.add_argument('--truth-start',type=int,help='Known-date validation only; never used to fit or select a prior')
    args=p.parse_args()
    try:
        r=evaluate(args.samples,args.reference,args.out,args.start_min,args.start_max,args.offset_sds,
                   args.offset_mean,args.allow_incomplete,args.truth_start)
    except (ValueError,KeyError,FileNotFoundError) as e:p.exit(2,'Input/run error: '+str(e)+'\n')
    print(json.dumps({k:r[k] for k in ['tool_version','map_span_years','max_pairwise_total_variation']},ensure_ascii=False))


if __name__=='__main__':main()
