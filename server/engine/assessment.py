"""Traceable model assessment of a saved annual Bayesian wiggle match."""
import argparse
import csv
import hashlib
import json
import math
import os
from pathlib import Path

import numpy as np

try:
    from .core import match,number
    from .prepare import read_csv
    from .run import tokens
    from .model_checks import fixed_agreement,offset_summary
    from .predictive import posterior_predictive,leave_one_out,compare_positions,design
except ImportError:
    from core import match,number
    from prepare import read_csv
    from run import tokens
    from model_checks import fixed_agreement,offset_summary
    from predictive import posterior_predictive,leave_one_out,compare_positions,design

VERSION='0.4.0'


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def load(path):return json.loads(Path(path).read_text(encoding='utf-8-sig'))
def save(path,value):Path(path).write_text(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')


def table(path,rows):
    if not rows:return
    keys=list(dict.fromkeys(k for r in rows for k in r))
    with Path(path).open('x',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=keys);w.writeheader()
        for r in rows:w.writerow({k:json.dumps(v,ensure_ascii=False) if isinstance(v,(dict,list)) else v for k,v in r.items()})


def no_overlap(samples,curve):
    if tokens(samples,'measurement_ids')&tokens(curve,'member_ids') or tokens(samples,'tree_id')&tokens(curve,'tree_ids'):
        raise ValueError('Target and reference have declared measurement/tree overlap')


def known_truth_score(result,truth):
    if truth is None:return dict(status='not_applicable_no_independent_truth')
    truth=number(truth,'truth_start',integer=True)
    if truth<result['requested_start_min'] or truth>result['requested_start_max']:
        return dict(status='outside_search_window',truth_start=truth)
    if truth in result['unassessed_start_years']:return dict(status='unassessed',truth_start=truth)
    intervals=result['hpd']['0.95']['intervals']
    return dict(status='covered' if any(a<=truth<=b for a,b in intervals) else 'missed',
                truth_start=truth,map_error_years=result['map_start_astro_year']-truth,
                hpd95_year_count=sum(b-a+1 for a,b in intervals),
                truth_independence='Supplied by caller; must be verified externally',
                truth_usage='Post-hoc score only, never used for model selection or fitting')


def wilson(success,n):
    z=1.959963984540054; p=success/n; den=1+z*z/n
    middle=(p+z*z/(2*n))/den; half=z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return [max(0.,middle-half),min(1.,middle+half)]


def simulation_coverage(samples,curve,baseline,replicates,seed):
    replicates=number(replicates,'simulation_replicates',positive=True,integer=True)
    seed=number(seed,'seed',integer=True)
    if seed<0:raise ValueError('Seed must be nonnegative')
    eligible=np.array([int(r['astro_year']) for r in baseline['posterior']])
    means,variances=design(samples,curve,eligible)
    rng=np.random.default_rng(seed); rows=[]
    for i in range(replicates):
        index=int(rng.integers(len(eligible))); truth=int(eligible[index])
        delta=float(rng.normal(baseline['offset_prior_mean_bp'],baseline['offset_prior_sd_bp']))
        y=means[index]+delta+rng.normal(size=len(samples))*np.sqrt(variances[index])
        generated=[dict(s,c14_age_bp=float(v)) for s,v in zip(samples,y)]
        # Neither truth nor delta is passed to the inference function.
        fit=match(generated,curve,baseline['requested_start_min'],baseline['requested_start_max'],
                  baseline['offset_prior_mean_bp'],baseline['offset_prior_sd_bp'])
        row=dict(replicate=i,data_kind='synthetic_model_check',true_start=truth,true_offset_bp=delta,
                 synthetic_c14_age_bp=y.tolist(),map_start=fit['map_start_astro_year'],
                 map_error_years=fit['map_start_astro_year']-truth)
        for label in ['0.95','0.9545']:
            h=fit['hpd'][label]
            row['covered_'+label]=any(a<=truth<=b for a,b in h['intervals'])
            row['hpd_year_count_'+label]=sum(b-a+1 for a,b in h['intervals'])
            row['hpd_actual_mass_'+label]=h['actual_mass']
        rows.append(row)
    coverage={}
    for label in ['0.95','0.9545']:
        count=sum(r['covered_'+label] for r in rows)
        coverage[label]=dict(covered=count,replicates=replicates,coverage_fraction=count/replicates,
                             wilson95_interval=wilson(count,replicates),
                             mean_posterior_actual_hpd_mass=float(np.mean([r['hpd_actual_mass_'+label] for r in rows])),
                             median_hpd_year_count=float(np.median([r['hpd_year_count_'+label] for r in rows])))
    return dict(replicates=replicates,seed=seed,coverage=coverage,runs=rows,
                position_generation='Uniform on assessable candidate positions, NOT baseline posterior probabilities',
                offset_generation='Declared Gaussian prior; independent effective measurement/reference errors',
                conditional_on_reference_coverage=baseline['conditional_on_reference_coverage'],
                truth_usage='Scoring only; not passed to match or used to select settings',
                interpretation='Synthetic within-model coverage check, not new measured data, independent archaeological accuracy, or a full SBC rank-uniformity test')


def execute(run_dir,out,seed,draws=2048,expanded_bounds=None,alternative_references=None,
            simulation_replicates=0,truth_start=None):
    out=Path(out).resolve(); run_dir=Path(run_dir).resolve()
    if out.exists():raise ValueError('Assessment output already exists; choose a new directory')
    if os.name=='nt' and out.drive.lower()=='c:':raise ValueError('Output on C: is prohibited')
    seed=number(seed,'seed',integer=True); draws=number(draws,'draws',positive=True,integer=True)
    count=number(simulation_replicates,'simulation_replicates',integer=True)
    if seed<0 or draws<2 or count<0:raise ValueError('Use nonnegative seed/replicates and at least two predictive draws')
    manifest=load(run_dir/'manifest.json'); summary=load(run_dir/'summary.json')
    if len(manifest['inputs'])!=2:raise ValueError('Expected a saved run with samples and reference inputs')
    inputs=manifest['inputs']
    for row in inputs:
        if sha(row['path'])!=row['sha256']:raise ValueError('Saved input hash changed: '+row['path'])
    core_hash=sha(Path(__file__).with_name('core.py'))
    if manifest['engine_sha256']!=core_hash:raise ValueError('Saved kernel identity differs from installed mathematical kernel')
    samples=read_csv(inputs[0]['path']); curve=read_csv(inputs[1]['path']); no_overlap(samples,curve)
    params=manifest['parameters']
    result=match(samples,curve,params['start_min'],params['start_max'],params['offset_mean_bp'],params['offset_sd_bp'])
    saved=read_csv(run_dir/'posterior.csv')
    if [int(r['astro_year']) for r in saved]!=[r['astro_year'] for r in result['posterior']]:
        raise ValueError('Saved posterior support differs from declared final model')
    error=max(abs(float(a['probability'])-b['probability']) for a,b in zip(saved,result['posterior']))
    conditional_error=max(max(abs(float(a[k])-b[k]) for k in ['offset_conditional_mean_bp','offset_conditional_sd_bp']) for a,b in zip(saved,result['posterior']))
    if error>1e-12 or conditional_error>1e-9:raise ValueError('Saved final posterior does not match input/model identity')
    for key in ['requested_start_min','requested_start_max','offset_prior_mean_bp','offset_prior_sd_bp','unassessed_start_years']:
        if summary[key]!=result[key]:raise ValueError('Saved summary and declared model disagree: '+key)
    sources=[dict(path=str(run_dir/name),sha256=sha(run_dir/name)) for name in ['summary.json','posterior.csv','manifest.json']]
    comparisons=[]; extra_sources=[]
    if expanded_bounds is not None:
        lo,hi=[number(v,'expanded_bound',integer=True) for v in expanded_bounds]
        if lo>params['start_min'] or hi<params['start_max'] or (lo,hi)==(params['start_min'],params['start_max']):
            raise ValueError('Expanded bounds must contain and enlarge the original search window')
        other=match(samples,curve,lo,hi,params['offset_mean_bp'],params['offset_sd_bp'])
        outside=sum(r['probability'] for r in other['posterior'] if r['astro_year']<params['start_min'] or r['astro_year']>params['start_max'])
        comparisons.append(dict(kind='expanded_search',start_min=lo,start_max=hi,
                                posterior_mass_outside_original_window=outside,
                                coverage_conditional=other['conditional_on_reference_coverage'],
                                comparison=compare_positions(result,other),summary={k:v for k,v in other.items() if k not in ['posterior','fit_at_map']},posterior=other['posterior']))
    for path in alternative_references or []:
        path=Path(path).resolve(); ref=read_csv(path); no_overlap(samples,ref)
        identity=dict(path=str(path),sha256=sha(path)); extra_sources.append(identity)
        try:
            other=match(samples,ref,params['start_min'],params['start_max'],params['offset_mean_bp'],params['offset_sd_bp'])
        except ValueError as e:
            if 'No assessable candidate positions' not in str(e):raise
            comparisons.append(dict(kind='alternative_reference',reference=identity,status='not_evaluated_no_assessable_positions'));continue
        comparisons.append(dict(kind='alternative_reference',reference=identity,status='computed',
                                coverage_conditional=other['conditional_on_reference_coverage'],
                                comparison=compare_positions(result,other),
                                summary={k:v for k,v in other.items() if k not in ['posterior','fit_at_map']},posterior=other['posterior']))
    report=dict(tool_version=VERSION,source_run=str(run_dir),
                model_identity=dict(kernel_sha256=core_hash,posterior_max_probability_difference=error,
                                    offset_conditional_max_difference=conditional_error,parameters=params,
                                    sources=sources,inputs=inputs,alternative_inputs=extra_sources),
                statistical_consistency=fixed_agreement(samples,curve,result),offset_joint=offset_summary(result),
                leave_one_out=leave_one_out(samples,curve,result),
                posterior_predictive=posterior_predictive(samples,curve,result,draws,seed),
                position_reference_comparisons=comparisons,
                comparison_status='computed_explicit_comparisons' if comparisons else 'not_requested_no_comparison_settings',
                simulation_coverage=simulation_coverage(samples,curve,result,count,seed+1) if count else dict(status='not_requested'),
                independent_truth=known_truth_score(result,truth_start),
                limitations=['No automatic five-gate acceptance or scientific uniqueness classification',
                             'Shared-offset final Acomb is not supplied by a fixed-offset formula',
                             'Cross-year/batch covariance, uncertain gaps and multi-ring samples remain unsupported',
                             'Physical offset explanation requires independent external evidence'])
    for row in inputs+extra_sources+sources:
        if sha(row['path'])!=row['sha256']:raise ValueError('Input or saved result changed during assessment')
    report['assessment_code_sha256']={name:sha(Path(__file__).with_name(name)) for name in ['assessment.py','model_checks.py','predictive.py']}
    out.mkdir(parents=True,exist_ok=False)
    save(out/'assessment.json',report)
    table(out/'leave_one_out.csv',[{k:v for k,v in r.items() if k!='posterior'} for r in report['leave_one_out']])
    table(out/'predictive_marginals.csv',report['posterior_predictive']['marginals'])
    table(out/'predictive_statistics.csv',report['posterior_predictive']['replicate_statistics'])
    table(out/'simulation_runs.csv',report['simulation_coverage'].get('runs',[]))
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run-dir',required=True);p.add_argument('--out',required=True)
    p.add_argument('--seed',required=True,type=int);p.add_argument('--draws',type=int,default=2048)
    p.add_argument('--expanded-start-min',type=int);p.add_argument('--expanded-start-max',type=int)
    p.add_argument('--reference-alternative',action='append',default=[])
    p.add_argument('--simulation-replicates',type=int,default=0)
    p.add_argument('--truth-start',type=int,help='Independent known-year validation only; never selects settings')
    a=p.parse_args()
    if (a.expanded_start_min is None)!=(a.expanded_start_max is None):p.error('Supply both expanded bounds')
    bounds=None if a.expanded_start_min is None else (a.expanded_start_min,a.expanded_start_max)
    try:r=execute(a.run_dir,a.out,a.seed,a.draws,bounds,a.reference_alternative,a.simulation_replicates,a.truth_start)
    except (ValueError,KeyError,OSError) as e:p.exit(2,'Assessment error: '+str(e)+'\n')
    print(json.dumps(dict(tool_version=VERSION,consistency_status=r['statistical_consistency']['status'],
                         loo_records=len(r['leave_one_out']),predictive_draws=r['posterior_predictive']['draws'],comparison_status=r['comparison_status'])))


if __name__=='__main__':main()
