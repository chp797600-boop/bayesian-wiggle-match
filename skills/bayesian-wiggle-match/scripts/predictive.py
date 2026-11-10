"""Exact leave-one-out refits and replicated-data posterior predictive checks."""
import math

import numpy as np

try:
    from .core import match,number
    from .model_checks import posterior_arrays,mixture_cdf,mixture_quantile
except ImportError:
    from core import match,number
    from model_checks import posterior_arrays,mixture_cdf,mixture_quantile


def compare_positions(baseline,comparison):
    a={int(r['astro_year']):float(r['probability']) for r in baseline['posterior']}
    b={int(r['astro_year']):float(r['probability']) for r in comparison['posterior']}
    common=sorted(a.keys()&b.keys())
    ma=sum(a[y] for y in common); mb=sum(b[y] for y in common)
    tv=.5*sum(abs(a[y]/ma-b[y]/mb) for y in common) if ma>0 and mb>0 else None
    same=a.keys()==b.keys()
    return dict(same_assessed_support=same,common_assessed_positions=len(common),
                baseline_mass_on_common_support=ma,comparison_mass_on_common_support=mb,
                total_variation_same_support=.5*sum(abs(a[y]-b[y]) for y in common) if same else None,
                conditional_common_support_tv=tv,
                map_shift_years=comparison['map_start_astro_year']-baseline['map_start_astro_year'],
                interpretation='Missing support is unassessed, never an excluded zero; restricted TV must be read with common-support masses')


def design(samples,curve,years):
    byyear={int(r['astro_year']):r for r in curve}
    means=np.array([[float(byyear[int(t)+int(s['ring_offset'])]['c14_age_bp']) for s in samples] for t in years])
    variances=np.array([[float(byyear[int(t)+int(s['ring_offset'])]['sigma_bp'])**2+float(s['sigma_bp'])**2 for s in samples] for t in years])
    return means,variances


def statistics(z,offsets):
    centered=np.asarray(offsets,dtype=float); centered-=centered.mean()
    norm=float(np.sqrt(centered@centered))
    trend=np.abs(z@centered/norm) if norm>0 else np.zeros(len(z))
    order=np.argsort(offsets)
    run=np.zeros(len(z),dtype=int); longest=np.zeros(len(z),dtype=int); previous=np.zeros(len(z))
    for col in order:
        sign=np.sign(z[:,col])
        run=np.where(sign==0,0,np.where(sign==previous,run+1,1))
        longest=np.maximum(longest,run); previous=sign
    return dict(squared_residual_sum=np.sum(z*z,axis=1),max_abs_residual=np.max(np.abs(z),axis=1),
                absolute_trend=trend,longest_same_sign_run=longest)


def posterior_predictive(samples,curve,result,draws,seed):
    draws=number(draws,'draws',positive=True,integer=True); seed=number(seed,'seed',integer=True)
    if draws<2 or seed<0:raise ValueError('At least two predictive draws and a nonnegative explicit seed are required')
    years,p,m,sd=posterior_arrays(result)
    means,variances=design(samples,curve,years)
    observed=np.array([float(s['c14_age_bp']) for s in samples])
    rng=np.random.default_rng(seed)
    index=rng.choice(len(years),size=draws,p=p)
    delta=rng.normal(m[index],sd[index])
    centers=means[index]+delta[:,None]
    scale=np.sqrt(variances[index])
    replicate=centers+rng.normal(size=centers.shape)*scale
    obs_stats=statistics((observed[None,:]-centers)/scale,[int(s['ring_offset']) for s in samples])
    rep_stats=statistics((replicate-centers)/scale,[int(s['ring_offset']) for s in samples])
    checks={}
    for key in obs_stats:
        frequency=float(np.mean(rep_stats[key]>=obs_stats[key]))
        checks[key]=dict(posterior_predictive_upper_tail_frequency=frequency,
                         monte_carlo_se=math.sqrt(frequency*(1-frequency)/draws),
                         observed_statistic_posterior_mean=float(np.mean(obs_stats[key])),
                         replicated_statistic_mean=float(np.mean(rep_stats[key])))
    marginals=[]
    for i,sample in enumerate(samples):
        center=means[:,i]+m; conditional_sd=np.sqrt(variances[:,i]+sd**2)
        mean=float(p@center)
        sigma=math.sqrt(float(p@(conditional_sd**2+(center-mean)**2)))
        marginals.append(dict(sample_id=sample['sample_id'],predictive_mean_bp=mean,predictive_sd_bp=sigma,
                              equal_tail95_bp=[mixture_quantile(q,p,center,conditional_sd) for q in [.025,.975]],
                              observed_c14_age_bp=float(sample['c14_age_bp'])))
    rows=[dict(draw=i,start_astro_year=int(years[index[i]]),shared_offset_bp=float(delta[i]),
               **{'observed_'+k:float(obs_stats[k][i]) for k in obs_stats},
               **{'replicated_'+k:float(rep_stats[k][i]) for k in rep_stats}) for i in range(draws)]
    return dict(draws=draws,seed=seed,checks=checks,marginals=marginals,replicate_statistics=rows,
                replicate_measurements=replicate.tolist(),
                conditional_on_reference_coverage=result.get('conditional_on_reference_coverage',False),
                interpretation='Replicate the fitted effective Gaussian model using one shared offset per sequence; not repeated assays of fixed physical rings. Predictive frequencies are not frequentist p-values, Acomb, or automatic acceptance.')


def heldout_prediction(sample,curve,result):
    byyear={int(r['astro_year']):r for r in curve}
    rows=[r for r in result['posterior'] if int(r['astro_year'])+int(sample['ring_offset']) in byyear]
    mass=sum(float(r['probability']) for r in rows)
    if mass<=0:return dict(heldout_predictive_assessed_mass=mass,heldout_predictive_status='unassessed')
    p=np.array([float(r['probability'])/mass for r in rows])
    center=np.array([float(byyear[int(r['astro_year'])+int(sample['ring_offset'])]['c14_age_bp'])+float(r['offset_conditional_mean_bp']) for r in rows])
    sd=np.sqrt(np.array([float(sample['sigma_bp'])**2+float(byyear[int(r['astro_year'])+int(sample['ring_offset'])]['sigma_bp'])**2+float(r['offset_conditional_sd_bp'])**2 for r in rows]))
    y=float(sample['c14_age_bp'])
    # Compute both tails directly: 1 - CDF loses small upper tails to rounding.
    cdf=float(sum(float(w)*.5*math.erfc((float(c)-y)/(math.sqrt(2)*float(s))) for w,c,s in zip(p,center,sd)))
    survival=float(sum(float(w)*.5*math.erfc((y-float(c))/(math.sqrt(2)*float(s))) for w,c,s in zip(p,center,sd)))
    log=-.5*((y-center)/sd)**2-np.log(sd)-.5*math.log(2*math.pi)
    positive=p>0
    weighted_log=log[positive]+np.log(p[positive])
    best=float(weighted_log.max()); density_log=best+math.log(float(np.exp(weighted_log-best).sum()))
    return dict(heldout_predictive_assessed_mass=mass,
                heldout_predictive_status='assessed' if abs(mass-1)<1e-10 else 'conditional_on_predictable_reference',
                heldout_predictive_cdf=cdf,heldout_predictive_survival=survival,
                heldout_predictive_two_tail_area=min(1.,2*min(cdf,survival)),
                heldout_predictive_log_density=density_log,
                heldout_predictive_mean_bp=float(p@center))


def leave_one_out(samples,curve,result):
    if len(samples)<3:return [dict(status='not_evaluated_insufficient_samples',reason='Core requires at least two retained annual samples')]
    out=[]
    for omitted,sample in enumerate(samples):
        retained=[dict(s) for i,s in enumerate(samples) if i!=omitted]
        shift=min(int(s['ring_offset']) for s in retained)
        for s in retained:s['ring_offset']=int(s['ring_offset'])-shift
        fit=match(retained,curve,result['requested_start_min']+shift,result['requested_start_max']+shift,
                  result['offset_prior_mean_bp'],result['offset_prior_sd_bp'])
        fit['map_start_astro_year']-=shift
        for row in fit['posterior']:row['astro_year']-=shift
        intervals=[[a-shift,b-shift] for a,b in fit['hpd']['0.95']['intervals']]
        out.append(dict(status='computed',heldout_sample_id=sample['sample_id'],omitted_offset=int(sample['ring_offset']),
                        rebased_by_years=shift,map_start_astro_year=fit['map_start_astro_year'],
                        hpd95_intervals=intervals,hpd95_year_count=sum(b-a+1 for a,b in intervals),
                        offset_mean_bp=fit['offset_posterior_mean_bp'],offset_sd_bp=fit['offset_posterior_sd_bp'],
                        comparison=compare_positions(result,fit),**heldout_prediction(sample,curve,fit),
                        posterior=fit['posterior'],
                        interpretation='Influence is not an outlier verdict or authorization to delete a sample'))
    return out
