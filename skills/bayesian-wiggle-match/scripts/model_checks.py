"""Model-matched overlap indices and shared-offset mixture diagnostics."""
import math

import numpy as np


def overlap_agreement(profiles):
    """Apply the published discrete overlap definition to normalized profiles.

    Each tuple contains a single-date likelihood and its modeled posterior on
    the same uniform grid. This is not a multi-parameter joint agreement index.
    """
    factors = []
    for likelihood, posterior in profiles:
        f = np.asarray(likelihood, dtype=float); q = np.asarray(posterior, dtype=float)
        if f.ndim != 1 or f.shape != q.shape or not len(f):
            raise ValueError('Agreement profiles require equal nonempty one-dimensional grids')
        if np.any(~np.isfinite(f)) or np.any(~np.isfinite(q)) or np.any(f < 0) or np.any(q < 0):
            raise ValueError('Agreement profiles must be finite and nonnegative')
        if abs(float(f.sum())-1)>1e-10 or abs(float(q.sum())-1)>1e-10:
            raise ValueError('Agreement profiles must be normalized probability masses')
        factors.append(float(f@q/(f@f)))
    if not factors:
        raise ValueError('At least one agreement profile is required')
    n = len(factors)
    acomb = 100*math.exp(sum(math.log(f) for f in factors)/math.sqrt(n)) if min(factors)>0 else 0.
    an = 100/math.sqrt(2*n)
    return dict(n=n, individual_percent=[100*f for f in factors],acomb_percent=acomb,
                an_percent=an,passes_combination_threshold=acomb>=an)


def fixed_agreement(samples, curve, result):
    base=dict(implementation='Published overlap-integral definition on the supplied annual reference domain',
              oxcal_equivalence_claim=None, threshold_is_not_accuracy_probability=True)
    if float(result['offset_prior_sd_bp']) != 0:
        return dict(base,status='not_evaluated_shared_offset',reason='Acomb for the fixed one-position combination is not transferred to a shared-unknown-offset model')
    years = sorted(int(r['astro_year']) for r in curve)
    if len(years)!=len(set(years)) or years[-1]-years[0]+1!=len(years):
        return dict(base,status='not_evaluated_reference_gaps',reason='Missing reference years are not zero likelihood')
    if result.get('conditional_on_reference_coverage',False):
        return dict(base,status='not_evaluated_incomplete_search',reason='The final posterior is conditional on an incompletely assessed search')
    byyear={int(r['astro_year']):r for r in curve}
    means=np.array([float(byyear[y]['c14_age_bp']) for y in years])
    curve_var=np.array([float(byyear[y]['sigma_bp'])**2 for y in years])
    mu=float(result['offset_prior_mean_bp'])
    profiles=[]; edge=[]
    for sample in samples:
        v=curve_var+float(sample['sigma_bp'])**2
        log=-.5*((float(sample['c14_age_bp'])-mu-means)**2/v+np.log(v))
        f=np.exp(log-log.max()); f/=f.sum()
        q=np.zeros(len(years))
        for row in result['posterior']:
            idx=int(row['astro_year'])+int(sample['ring_offset'])-years[0]
            if not 0<=idx<len(q):
                raise ValueError('Final posterior uses dates outside the reference domain')
            q[idx]=float(row['probability'])
        profiles.append((f,q))
        edge.append(float(f[:min(5,len(f))].sum()+f[max(5,len(f)-5):].sum()))
    computed=overlap_agreement(profiles)
    return dict(base,**computed,status='computed_fixed_offset',
                reference_start_year=years[0],reference_end_year=years[-1],
                likelihood_edge_five_year_masses=edge,
                search_start_year=result['requested_start_min'],search_end_year=result['requested_start_max'],
                fixed_offset_bp=mu,
                interpretation='Conditioned on supplied reference domain and declared position prior; finite reference tails or search truncation can change the index. Not a portable OxCal-equivalence or five-gate pass.')


def posterior_arrays(result):
    rows=result['posterior']
    years=np.array([int(r['astro_year']) for r in rows])
    p=np.array([float(r['probability']) for r in rows])
    m=np.array([float(r['offset_conditional_mean_bp']) for r in rows])
    sd=np.array([float(r['offset_conditional_sd_bp']) for r in rows])
    return years,p,m,sd


def mixture_cdf(value,p,means,sds):
    return float(sum(float(w)*(.5*(1+math.erf((value-float(m))/(math.sqrt(2)*float(s)))) if s>0 else float(value>=m))
                     for w,m,s in zip(p,means,sds)))


def mixture_quantile(q,p,means,sds):
    if not 0<q<1:raise ValueError('Quantile must be in (0,1)')
    lo=float(np.min(means-12*sds))-1; hi=float(np.max(means+12*sds))+1
    for _ in range(65):
        mid=(lo+hi)/2
        if mixture_cdf(mid,p,means,sds)<q:lo=mid
        else:hi=mid
    return float((lo+hi)/2)


def offset_summary(result):
    years,p,m,sd=posterior_arrays(result)
    mean=float(p@m); date_mean=float(p@years)
    between=float(p@((m-mean)**2)); variance=float(p@(sd**2))+between
    date_var=float(p@((years-date_mean)**2))
    cov=float(p@((years-date_mean)*(m-mean)))
    tau=float(result['offset_prior_sd_bp'])
    if tau==0:
        # A fixed parameter is a point mass even when stored weights sum to
        # 1 +/- floating-point rounding; do not invent a variance/correlation.
        mean=float(result['offset_prior_mean_bp']); between=0.; variance=0.; cov=0.
    quantiles=[mean,mean] if variance==0 else [mixture_quantile(q,p,m,sd) for q in [.025,.975]]
    return dict(mean_bp=mean,sd_bp=math.sqrt(variance),equal_tail95_bp=quantiles,
                sd_contraction=tau/math.sqrt(variance) if tau>0 and variance>0 else None,
                date_offset_covariance=cov,
                date_offset_correlation=max(-1.,min(1.,cov/math.sqrt(date_var*variance))) if date_var>0 and variance>0 else None,
                between_date_variance_fraction=between/variance if variance>0 else None,
                fixed_offset=tau==0, physical_explanation_status='not_independently_verified',
                interpretation='Exact mixture moments and equal-tail interval, not offset HPD or proof of a physical cause')
