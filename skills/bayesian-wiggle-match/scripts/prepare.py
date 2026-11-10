"""Source-preserving import and explicit annual reference construction."""
import argparse
import csv
import hashlib
import io
import json
import math
from collections import defaultdict
from pathlib import Path

try:
    from .core import number
except ImportError:
    from core import number


def write_csv(path,rows):
    rows=list(rows)
    if not rows:
        raise ValueError('Refusing to write an empty data table: '+str(path))
    with Path(path).open('x',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)


def read_csv(path):
    with Path(path).open(encoding='utf-8-sig',newline='') as f:
        return list(csv.DictReader(f))


def import_collection(raw_root):
    root=Path(raw_root).resolve(); kept={}; exclusions=[]; sources=[]; duplicate_count=0
    def add(row,origin):
        nonlocal duplicate_count
        key=row['measurement_id']
        if key in kept:
            old=kept[key]
            for k in ['astro_year','c14_age_bp','sigma_bp','tree_id']:
                if old[k]!=row[k]:
                    raise ValueError('Conflicting duplicate '+key+' field '+k)
            duplicate_count+=1
            old['origins']=json.dumps(json.loads(old['origins'])+[origin],ensure_ascii=False)
            return
        row['origins']=json.dumps([origin],ensure_ascii=False)
        kept[key]=row
    for path in sorted(root.glob('*_IntChron/*.json')):
        source_hash=hashlib.sha256(path.read_bytes()).hexdigest()
        sources.append(dict(path=str(path),sha256=source_hash))
        obj=json.loads(path.read_text(encoding='utf-8'))
        for rec in obj.get('records',[]):
            fd=rec.get('file_data',{})
            for series in fd.get('series_list',[]):
                data=series.get('data',{})
                if not isinstance(data,dict) or 'r_date' not in data:
                    continue
                for i in range(len(data['r_date'])):
                    origin=dict(file=str(path.relative_to(root)),sha256=source_hash,record=rec['record'],series=series['series'],row_index_0=i)
                    def get(k,default=None):
                        a=data.get(k)
                        return a[i] if isinstance(a,list) and i<len(a) else default
                    try:
                        header=fd.get('header',{})
                        if header.get('site_type')!='Dendrochronological' or number(header.get('latitude'),'latitude')<=0:
                            raise ValueError('not a confirmed Northern Hemisphere tree-ring record')
                        if series.get('intcal_set_type')!='Int':
                            raise ValueError('not Northern Hemisphere Int dataset')
                        if get('r_date_sigma_batch') not in [None,0]:
                            raise ValueError('declared correlated batch uncertainty not supported by v0.1')
                        if get('r_date_sigma_extra',0) not in [None,0] or get('r_date_sigma_mult',1) not in [None,1]:
                            raise ValueError('additional or multiplicative uncertainty requires separate interpretation')
                        if get('intcal_reservoir',0) not in [None,0] or get('intcal_reservoir_sigma',0) not in [None,0]:
                            raise ValueError('reservoir correction not supported')
                        if number(get('calage_range'),'sample width')!=1:
                            raise ValueError('not a one-year sample')
                        for name in ['calage_sigmaI','calage_sigmaD','calage_sigmaC']:
                            if number(get(name),'calendar uncertainty')!=0:
                                raise ValueError('nonzero or undefined calendar uncertainty')
                        year=1950-number(get('calage'),'calendar BP',integer=True)
                        identity=get('intcal_data_id')
                        if identity is None:
                            raise ValueError('missing stable IntCal measurement ID')
                        row=dict(measurement_id='intchron:'+str(identity),
                                 tree_id='intchron:'+rec['record']+':'+str(series.get('series_sample')),
                                 series_id=series['series'],labcode=str(get('labcode','')),
                                 astro_year=year,c14_age_bp=number(get('r_date'),'14C age'),
                                 sigma_bp=number(get('r_date_sigma'),'14C sigma',positive=True),
                                 sample_span_years=1,calendar_basis='IntChron calendar BP; all three calendar sigma fields zero',
                                 source_dois=';'.join(series.get('refs',[])))
                        add(row,origin)
                    except ValueError as e:
                        if str(e).startswith('Conflicting'):
                            raise
                        exclusions.append(dict(**origin,reason=str(e)))
    pangaea=root/'Brehm2021_PANGAEA921808/extracted/datasets'
    for path in sorted(pangaea.glob('*.tab')):
        source_hash=hashlib.sha256(path.read_bytes()).hexdigest()
        sources.append(dict(path=str(path),sha256=source_hash))
        raw=path.read_text(encoding='utf-8-sig'); body=raw.split('*/',1)[-1].strip()
        for i,r in enumerate(csv.DictReader(io.StringIO(body),delimiter='\t')):
            origin=dict(file=str(path.relative_to(root)),sha256=source_hash,row_index_0=i)
            try:
                if r.get('Interval comm (block interval)')!='EL':
                    raise ValueError('Unrecognized annual ring segment')
                bp=number(r['Age [ka BP]'],'calendar ka BP')*1000
                if not math.isclose(bp,round(bp),abs_tol=1e-7):
                    raise ValueError('Noninteger calendar year')
                target=r.get('Sample ID') or r.get('ID (measurement)')
                if not target:
                    raise ValueError('Missing stable target ID')
                row=dict(measurement_id='pangaea921808:'+target,
                         tree_id='pangaea921808:'+r['Name'],series_id=r['Name'],labcode=r.get('Lab label',''),
                         astro_year=1950-round(bp),c14_age_bp=number(r['Age dated [ka]'],'14C ka')*1000,
                         sigma_bp=number(r['Age dated std dev [±]'],'14C sigma ka',positive=True)*1000,
                         sample_span_years=1,calendar_basis='published absolute annual tree chronology; EL=annual early+late wood; Age ka BP to year',
                         source_dois='10.1594/PANGAEA.921808;10.1038/s41561-020-00674-0')
                add(row,origin)
            except ValueError as e:
                if str(e).startswith('Conflicting'):
                    raise
                exclusions.append(dict(**origin,reason=str(e)))
    rows=sorted(kept.values(),key=lambda r:(r['astro_year'],r['tree_id'],r['measurement_id']))
    return rows,dict(sources=sources,excluded=exclusions,duplicate_import_rows=duplicate_count,
                     imported_measurements=len(rows),scope='Six IntChron exports plus Brehm2021 PANGAEA annual records; other downloaded formats not yet normalized',
                     calendar_conversion='astronomical year = 1950 - integer calendar BP; 0=1 BCE')


def build_reference(rows,exclude_trees=()):
    groups=defaultdict(list); seen={}
    for r in rows:
        if r['tree_id'] in set(exclude_trees):continue
        identity=r['measurement_id']
        values=(number(r['astro_year'],'astro_year',integer=True),
                number(r['c14_age_bp'],'14C age'),number(r['sigma_bp'],'sigma',positive=True),r['tree_id'])
        if identity in seen:
            if seen[identity]!=values:raise ValueError('Conflicting duplicate measurement '+identity)
            continue
        seen[identity]=values
        if r.get('sample_span_years') is not None and number(r['sample_span_years'],'span')!=1:
            raise ValueError('Reference contains nonannual sample')
        groups[values[0]].append((r,values))
    output=[]
    for year,group in sorted(groups.items()):
        weights=[1/v[2]**2 for r,v in group]; W=sum(weights)
        mean=sum(w*v[1] for w,(r,v) in zip(weights,group))/W
        q=sum(w*(v[1]-mean)**2 for w,(r,v) in zip(weights,group))
        # Explicit conservative precision convention; no observations removed.
        birge=max(1.,q/(len(group)-1)) if len(group)>1 else 1.
        output.append(dict(astro_year=year,c14_age_bp=mean,sigma_bp=math.sqrt(birge/W),
                           n_measurements=len(group),birge_variance_factor=birge,
                           member_ids=';'.join(r['measurement_id'] for r,v in group),
                           tree_ids=';'.join(sorted(set(r['tree_id'] for r,v in group)))))
    return output


def main():
    p=argparse.ArgumentParser(description=__doc__)
    sub=p.add_subparsers(dest='command',required=True)
    imp=sub.add_parser('import'); imp.add_argument('--raw-root',required=True); imp.add_argument('--out',required=True)
    ref=sub.add_parser('reference'); ref.add_argument('--measurements',required=True); ref.add_argument('--out',required=True)
    ref.add_argument('--tree',action='append',required=True,help='Explicit source tree; repeat to select several')
    args=p.parse_args()
    out=Path(args.out).resolve()
    if out.exists():raise ValueError('Output directory already exists')
    if args.command=='import':
        rows,report=import_collection(args.raw_root)
        out.mkdir(parents=True); write_csv(out/'measurements.csv',rows)
        (out/'provenance.json').write_text(json.dumps(report,ensure_ascii=False,indent=2),encoding='utf-8')
        print(json.dumps({'measurements':len(rows),'duplicates':report['duplicate_import_rows'],'excluded':len(report['excluded'])}))
    else:
        all_rows=read_csv(args.measurements)
        available={r['tree_id'] for r in all_rows}
        if not set(args.tree)<=available:raise ValueError('Unknown tree IDs')
        selected=[r for r in all_rows if r['tree_id'] in args.tree]
        rows=build_reference(selected)
        out.mkdir(parents=True); write_csv(out/'reference.csv',rows)
        (out/'provenance.json').write_text(json.dumps(dict(trees=args.tree,measurement_file=str(Path(args.measurements).resolve()),measurement_sha256=hashlib.sha256(Path(args.measurements).read_bytes()).hexdigest(),aggregation='inverse variance annual mean; error multiplied by sqrt(max(1,Q/(n-1))); independent measurement errors assumed'),ensure_ascii=False,indent=2),encoding='utf-8')
        print(json.dumps({'reference_years':len(rows)}))


if __name__=='__main__':main()
