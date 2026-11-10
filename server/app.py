"""Project-local tools around the frozen annual Gaussian matching engine."""
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import sys

PACKAGE=Path(__file__).resolve().parent
if str(PACKAGE) not in sys.path:sys.path.insert(0,str(PACKAGE))
from engine import run,assessment
from engine.core import number


class Backend:
    def __init__(self,workspace):
        self.workspace=Path(workspace).resolve()
        if self.workspace.drive.lower()=='c:':raise ValueError('Workspace must not be on C:')
        self.output_root=self.workspace/'output/wiggle_mcp_runs'
        self.catalog=json.loads((PACKAGE/'data/catalog.json').read_text(encoding='utf-8'))
        for name,expected in json.loads((PACKAGE/'engine_manifest.json').read_text(encoding='utf-8')).items():
            if hashlib.sha256((PACKAGE/'engine'/name).read_bytes()).hexdigest()!=expected:
                raise ValueError('Frozen calculation engine hash mismatch')

    def _input(self,path):
        p=Path(path)
        if not p.is_absolute():p=self.workspace/p
        p=p.resolve()
        if not p.is_relative_to(self.workspace) or p.drive.lower()=='c:':
            raise ValueError('Input must be a file inside the configured E-drive project workspace')
        if not p.is_file():raise ValueError('Input file does not exist')
        return p

    def _output_path(self,run_id):
        if not isinstance(run_id,str) or not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,63}',run_id):
            raise ValueError('Use a new run_id containing only letters, numbers, underscore or hyphen')
        p=(self.output_root/run_id).resolve()
        if not p.is_relative_to(self.workspace) or not p.is_relative_to(self.output_root.resolve()):
            raise ValueError('Output must remain inside project output root')
        if p.exists():raise ValueError('Output exists; use a new run_id')
        return p

    def _reference_path(self,row):
        p=(PACKAGE/'data'/row['file']).resolve()
        if not p.is_relative_to(PACKAGE):raise ValueError('Invalid bundled reference path')
        if hashlib.sha256(p.read_bytes()).hexdigest()!=row['sha256']:
            raise ValueError('Reference hash mismatch; do not interpret altered reference data')
        return p

    def get_reference_info(self,reference_id):
        row=next((x for x in self.catalog['references'] if x['reference_id']==reference_id),None)
        if row is None:raise ValueError('Unknown reference_id; call list_reference_data')
        self._reference_path(row)
        return dict(row)

    def list_reference_data(self,start_astro_year=None,end_astro_year=None,include_demo=False):
        if (start_astro_year is None)!=(end_astro_year is None):raise ValueError('Supply both date bounds or neither')
        if start_astro_year is not None:
            start_astro_year=number(start_astro_year,'start',integer=True)
            end_astro_year=number(end_astro_year,'end',integer=True)
            if end_astro_year<start_astro_year:raise ValueError('Date bounds reversed')
        refs=[]
        for row in self.catalog['references']:
            if row['is_demo'] and not include_demo:continue
            if start_astro_year is not None and (row['end_astro_year']<start_astro_year or row['start_astro_year']>end_astro_year):continue
            self._reference_path(row);refs.append(dict(row))
        return dict(references=refs,calendar=self.catalog['calendar'],
                    interpretation='Temporal overlap is not a scientific suitability decision; library coverage is incomplete')

    def validate_floating_series(self,samples_path):
        path=self._input(samples_path)
        with path.open(encoding='utf-8-sig',newline='') as f:rows=list(csv.DictReader(f))
        if not 2<=len(rows)<=2000:raise ValueError('Supply 2 to 2000 single-ring sample rows')
        offsets=[];seen_measurements=set()
        for row in rows:
            if number(row.get('sample_span_years'),'sample_span_years')!=1:raise ValueError('Only single annual rings supported')
            offsets.append(number(row.get('ring_offset'),'ring_offset',integer=True))
            number(row.get('c14_age_bp'),'c14_age_bp');number(row.get('sigma_bp'),'sigma_bp',positive=True)
            for key in ['sample_id','measurement_ids','tree_id']:
                if not str(row.get(key,'')).strip():raise ValueError('Missing traceability field '+key)
            identities={x.strip() for x in row['measurement_ids'].split(';') if x.strip()}
            if not identities or seen_measurements & identities:
                raise ValueError('measurement_ids must identify distinct measurements across ring offsets')
            seen_measurements.update(identities)
        if min(offsets)!=0 or len(offsets)!=len(set(offsets)):raise ValueError('Offsets must be distinct, nonnegative, earliest zero')
        if len({r['sample_id'] for r in rows})!=len(rows):raise ValueError('sample_id must be unique')
        return dict(status='input_structure_valid',annual_samples=len(rows),span_years=max(offsets)+1,
                    continuous_measurements=len(rows)==max(offsets)+1,input_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),
                    scientific_acceptance='not_automatically_evaluated',
                    limitations=['Relative ring gaps and source identities remain externally supplied facts',
                                 'A valid input is not a reliable-date acceptance'])

    def run_wiggle_match(self,samples_path,reference_id,start_min,start_max,offset_mean_bp,offset_sd_bp,run_id,allow_incomplete=False):
        self.validate_floating_series(samples_path)
        info=self.get_reference_info(reference_id);ref=self._reference_path(info)
        samples=self._input(samples_path)
        target=run.read_csv(samples)
        aliases={x.casefold() for x in info['physical_tree_aliases']}
        supplied={x.strip().casefold() for r in target for x in r['tree_id'].split(';')}
        if aliases & supplied:raise ValueError('Target and reference identify the same physical tree')
        out=self._output_path(run_id)
        out.parent.mkdir(parents=True,exist_ok=True)
        result=run.execute(samples,ref,out,start_min,start_max,offset_mean_bp,offset_sd_bp,allow_incomplete)
        summary={k:v for k,v in result.items() if k not in ['posterior','fit_at_map']}
        wrapper=dict(mcp_version='0.1.0',reference=info,run_id=run_id,
          scientific_acceptance='not_automatically_evaluated',
          output_directory=str(out),summary=summary,
          limitations=['All probabilities condition on declared search and supported reference positions',
           'Years outside this measured segment are not excluded by the data library',
           'Shared-unknown-offset final Acomb and automatic five-gate acceptance remain unavailable'])
        (out/'mcp_context.json').write_text(json.dumps({k:v for k,v in wrapper.items() if k!='summary'},ensure_ascii=False,indent=2),encoding='utf-8')
        return wrapper

    def assess_match(self,run_id,assessment_id,seed,draws=256):
        if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_-]{0,63}',run_id):raise ValueError('Invalid run_id')
        source=(self.output_root/run_id).resolve()
        if not source.is_relative_to(self.workspace) or not source.is_relative_to(self.output_root.resolve()):raise ValueError('Invalid run path')
        if not (source/'mcp_context.json').is_file():raise ValueError('No MCP matching run with this ID')
        draws=number(draws,'draws',integer=True);seed=number(seed,'seed',integer=True)
        if not 32<=draws<=4096:raise ValueError('draws must be 32 to 4096; report Monte Carlo uncertainty')
        out=self._output_path(assessment_id)
        report=assessment.execute(source,out,seed,draws)
        return dict(status='assessment_completed',output_directory=str(out),
           consistency=report['statistical_consistency'],leave_one_out_records=len(report['leave_one_out']),
           predictive_draws=report['posterior_predictive']['draws'],
           scientific_acceptance='not_automatically_evaluated',
           independent_truth='not_supplied_for_unknown_floating_series',
           limitations=report['limitations'])
