"""Exact Gaussian marginalization of one shared radiocarbon offset.

Independent reimplementation; not OxCal and not the Pearson 2026 author code.
Calendar integers are astronomical years (0 = 1 BCE). Ages/errors are 14C BP.
"""
import math

import numpy as np


def number(value, name, positive=False, integer=False):
    try:
        value=float(value)
    except (TypeError,ValueError) as e:
        raise ValueError(name+' must be numeric') from e
    if not math.isfinite(value) or (positive and value<=0):
        raise ValueError(name+' must be finite'+(' and positive' if positive else ''))
    if integer:
        if not value.is_integer():
            raise ValueError(name+' must be an integer')
        return int(value)
    return value


def hpd_sets(years, probabilities, levels=(.6827,.9545,.95,.99)):
    years=np.asarray(years,dtype=int)
    p=np.asarray(probabilities,dtype=float)
    if len(years)!=len(p) or not len(p) or not np.all(np.isfinite(p)) or np.any(p<0) or not np.isclose(p.sum(),1):
        raise ValueError('HPD requires normalized finite probabilities')
    out={}
    for level in levels:
        if not 0<float(level)<=1:
            raise ValueError('HPD mass must be in (0,1]')
        ranked=np.sort(p)[::-1]
        k=min(int(np.searchsorted(np.cumsum(ranked),level)),len(p)-1)
        cutoff=ranked[k]
        # Include exact probability ties rather than arbitrarily selecting a year.
        mask=p>=cutoff
        selected=sorted(int(x) for x in years[mask])
        intervals=[]
        for year in selected:
            if intervals and year==intervals[-1][1]+1:
                intervals[-1][1]=year
            else:
                intervals.append([year,year])
        out[str(level)]={'intervals':intervals,'actual_mass':float(p[mask].sum()),'threshold':float(cutoff)}
    return out


def match(samples,curve,start_min,start_max,offset_mean=0.,offset_sd=50.):
    """Match annual samples at exact offsets to a discrete annual reference.

    Position prior: uniform over the explicitly requested integer search years.
    Missing reference windows are UNASSESSED, not scientifically excluded.
    Returned probabilities condition on the subset of assessable positions.
    Different reference years and measurement errors are assumed independent.
    """
    start_min=number(start_min,'start_min',integer=True)
    start_max=number(start_max,'start_max',integer=True)
    if start_max<start_min or start_max-start_min>100000:
        raise ValueError('Search window must be ordered and at most 100001 years')
    mu=number(offset_mean,'offset_mean')
    tau=number(offset_sd,'offset_sd')
    if tau<0 or not math.isfinite(tau*tau):
        raise ValueError('offset_sd must be nonnegative with finite variance')
    if len(samples)<2:
        raise ValueError('At least two annual samples are required')
    clean=[]
    for row in samples:
        if number(row.get('sample_span_years'),'sample_span_years')!=1:
            raise ValueError('Only single-year samples supported; do not assign multi-ring averages to midpoints')
        clean.append(dict(sample_id=str(row.get('sample_id','')),
                          ring_offset=number(row.get('ring_offset'),'ring_offset',integer=True),
                          c14_age_bp=number(row.get('c14_age_bp'),'c14_age_bp'),
                          sigma_bp=number(row.get('sigma_bp'),'sigma_bp',positive=True)))
    offsets=[s['ring_offset'] for s in clean]
    if len(offsets)!=len(set(offsets)):
        raise ValueError('Repeated ring offsets share reference error; combine justified replicates upstream')
    if min(offsets)!=0:
        raise ValueError('ring_offset must be nonnegative, with earliest sample at offset 0')
    byyear={}
    for row in curve:
        year=number(row.get('astro_year'),'astro_year',integer=True)
        if year in byyear:
            raise ValueError('Reference must have exactly one mean and error per year')
        byyear[year]=dict(c14_age_bp=number(row.get('c14_age_bp'),'reference age'),
                          sigma_bp=number(row.get('sigma_bp'),'reference error',positive=True))
    y=np.array([s['c14_age_bp'] for s in clean])
    measurement_var=np.array([s['sigma_bp']**2 for s in clean])
    positions=[]; unassessed=[]; logs=[]; conditional=[]
    for start in range(start_min,start_max+1):
        years=[start+d for d in offsets]
        if any(t not in byyear for t in years):
            unassessed.append(start)
            continue
        ref=np.array([byyear[t]['c14_age_bp'] for t in years])
        v=measurement_var+np.array([byyear[t]['sigma_bp']**2 for t in years])
        if not np.all(np.isfinite(v)):
            raise ValueError('Age-error variance overflow')
        w=1/v; residual=y-ref; W=float(w.sum()); tau2=tau*tau
        if tau2==0:
            cm=mu; cv=0.; q=float(np.sum(w*(residual-mu)**2)); logdet=float(np.log(v).sum())
        else:
            cv=1/(W+1/tau2)
            cm=mu+cv*float(np.sum(w*(residual-mu)))
            q=float(np.sum(w*(residual-cm)**2)+(cm-mu)**2/tau2)
            logdet=float(np.log(v).sum()+np.log1p(tau2*W))
        logs.append(-.5*(len(clean)*math.log(2*math.pi)+logdet+q))
        positions.append(start); conditional.append((cm,cv))
    if not positions:
        raise ValueError('No assessable candidate positions: reference coverage is insufficient')
    log=np.array(logs)
    if not np.all(np.isfinite(log)):
        raise ValueError('Likelihood is nonfinite; inspect input magnitude and errors')
    p=np.exp(log-log.max()); p/=p.sum()
    m=np.array([x[0] for x in conditional]); variance=np.array([x[1] for x in conditional])
    mean=float(p@m)
    # Centered mixture variance avoids cancellation for large common offsets.
    sd=float(np.sqrt(p@(variance+(m-mean)**2)))
    imap=int(np.argmax(p)); mapyear=positions[imap]
    edge_mass=float(sum(v for t,v in zip(positions,p) if t-start_min<5 or start_max-t<5))
    warnings=[]
    if unassessed:
        warnings.append('Some positions are unassessed. Posterior is conditional on reference coverage; missing windows are not excluded dates.')
    if edge_mass>.01:
        warnings.append('More than 1% posterior mass lies within five years of the requested search boundaries; review/expand the search window.')
    if len(positions)<10:
        warnings.append('Fewer than ten assessable years: high posterior probability can reflect the restricted candidate set rather than uniquely diagnostic data.')
    fit=[]
    for sample in clean:
        year=mapyear+sample['ring_offset']; ref=byyear[year]
        fit.append({**sample,'astro_year_at_map':year,'reference_c14_age_bp':ref['c14_age_bp'],
                    'reference_sigma_bp':ref['sigma_bp'],
                    'residual_after_conditional_offset_bp':sample['c14_age_bp']-ref['c14_age_bp']-conditional[imap][0]})
    return dict(model='annual_known_gap_shared_gaussian_offset_v0.1',
                position_prior='uniform over integer start years in requested window',
                calendar_scale='astronomical integer year: 0=1 BCE, -1=2 BCE',
                offset_sign='positive means sample radiocarbon age is older than reference',
                offset_prior_mean_bp=mu,offset_prior_sd_bp=tau,
                requested_start_min=start_min,requested_start_max=start_max,
                assessed_positions=len(positions),unassessed_start_years=unassessed,
                conditional_on_reference_coverage=bool(unassessed),
                map_start_astro_year=mapyear,map_probability=float(p[imap]),
                offset_posterior_mean_bp=mean,offset_posterior_sd_bp=sd,
                hpd=hpd_sets(positions,p),boundary_five_year_mass=edge_mass,warnings=warnings,
                posterior=[dict(astro_year=t,probability=float(prob),offset_conditional_mean_bp=cm,
                                offset_conditional_sd_bp=math.sqrt(cv),log_marginal_likelihood=ll)
                           for t,prob,(cm,cv),ll in zip(positions,p,conditional,logs)],fit_at_map=fit)
