"""Run annual Bayesian matching from explicit, traceable CSV inputs."""
import argparse
import hashlib
import json
import os
import platform
import sys
from pathlib import Path

import numpy as np

try:
    from .core import match
    from .prepare import read_csv,write_csv
    from .diagnostics import analyze,write_outputs
    from .model_checks import fixed_agreement,offset_summary
except ImportError:
    from core import match
    from prepare import read_csv,write_csv
    from diagnostics import analyze,write_outputs
    from model_checks import fixed_agreement,offset_summary


def tokens(rows,field):
    result=set()
    for row in rows:
        value=row.get(field,'')
        if not value.strip():raise ValueError('Missing traceability field '+field)
        result.update(x.strip() for x in value.split(';') if x.strip())
    return result


def execute(samples_path,reference_path,out,start_min,start_max,offset_mean=0,offset_sd=None,allow_incomplete=False):
    if offset_sd is None:
        raise ValueError('An explicit offset_sd in 14C BP is required; no default scientific prior is selected')
    out=Path(out).resolve()
    if out.exists():raise ValueError('Output directory already exists; choose a new directory')
    if os.name=='nt' and out.drive.lower()=='c:':raise ValueError('Output on C: is prohibited')
    samples=read_csv(samples_path); curve=read_csv(reference_path)
    measurement_overlap=tokens(samples,'measurement_ids') & tokens(curve,'member_ids')
    tree_overlap=tokens(samples,'tree_id') & tokens(curve,'tree_ids')
    if measurement_overlap or tree_overlap:
        raise ValueError('Reference/target overlap in declared measurement IDs or tree IDs; cannot establish independent matching')
    result=match(samples,curve,start_min,start_max,offset_mean,offset_sd)
    if result['conditional_on_reference_coverage'] and not allow_incomplete:
        raise ValueError('Reference leaves '+str(len(result['unassessed_start_years']))+
                         ' candidate positions unassessed. Supply adequate coverage, or explicitly use --allow-incomplete for conditional exploratory output.')
    result['result_status']='conditional_exploratory' if result['conditional_on_reference_coverage'] else 'window_conditional_research_result'
    result['posterior_diagnostics']=analyze([r['astro_year'] for r in result['posterior']],
                                          [r['probability'] for r in result['posterior']],start_min,start_max)
    result['statistical_consistency']=fixed_agreement(samples,curve,result)
    result['offset_joint']=offset_summary(result)
    # Finish calculation and validation before creating any output directory.
    out.mkdir(parents=True,exist_ok=False)
    write_csv(out/'posterior.csv',result['posterior'])
    write_csv(out/'fit_at_map.csv',result['fit_at_map'])
    write_outputs(out,result['posterior_diagnostics'],create=False)
    (out/'model_checks.json').write_text(json.dumps({k:result[k] for k in ['statistical_consistency','offset_joint']},ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
    summary={k:v for k,v in result.items() if k not in ['posterior','fit_at_map']}
    (out/'summary.json').write_text(json.dumps(summary,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
    inputs=[]
    for p in [samples_path,reference_path]:
        p=Path(p).resolve(); inputs.append(dict(path=str(p),sha256=hashlib.sha256(p.read_bytes()).hexdigest()))
    manifest=dict(inputs=inputs,engine_version='0.4.0',mathematical_kernel_version='0.1.0',python=platform.python_version(),numpy=np.__version__,
                  engine_sha256=hashlib.sha256(Path(__file__).with_name('core.py').read_bytes()).hexdigest(),
                  entrypoint_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                  diagnostics_sha256=hashlib.sha256(Path(__file__).with_name('diagnostics.py').read_bytes()).hexdigest(),
                  model_checks_sha256=hashlib.sha256(Path(__file__).with_name('model_checks.py').read_bytes()).hexdigest(),
                  parameters=dict(start_min=start_min,start_max=start_max,offset_mean_bp=offset_mean,offset_sd_bp=offset_sd,allow_incomplete=allow_incomplete),
                  implementation='Independent analytic Gaussian model; not OxCal or reproduction of Pearson2026 code',
                  assumptions=['Exact integer ring gaps; one year per sample; one target per ring offset',
                               'Independent measurement errors and independent reference errors between years',
                               'One shared constant radiocarbon offset; no automatic outlier rejection',
                               'Uniform position prior within stated window; probabilities conditional on assessable positions'],
                  overlap_check='No overlap in supplied measurement and tree IDs; accuracy of external provenance remains user responsibility')
    (out/'manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2),encoding='utf-8')
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--samples',required=True); parser.add_argument('--reference',required=True)
    parser.add_argument('--out',required=True); parser.add_argument('--start-min',required=True,type=int)
    parser.add_argument('--start-max',required=True,type=int)
    parser.add_argument('--offset-mean',type=float,default=0,help='sample minus reference, 14C BP')
    parser.add_argument('--offset-sd',type=float,required=True,help='Explicit scientific prior SD in 14C BP; 0 means fixed offset')
    parser.add_argument('--allow-incomplete',action='store_true',help='Explicitly permit coverage-conditioned exploratory results')
    args=parser.parse_args()
    try:
        r=execute(args.samples,args.reference,args.out,args.start_min,args.start_max,args.offset_mean,args.offset_sd,args.allow_incomplete)
    except (ValueError,KeyError,FileNotFoundError) as e:
        parser.exit(2,'Input/run error: '+str(e)+'\n')
    print(json.dumps({k:r[k] for k in ['map_start_astro_year','map_probability','offset_posterior_mean_bp','offset_posterior_sd_bp','warnings']},ensure_ascii=False))


if __name__=='__main__':main()
