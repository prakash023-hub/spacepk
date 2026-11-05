"""
================================================================================
LAYER 3: BAYESIAN PK ESTIMATION (pooled, study-level summary data)
SpacePK Framework
PyMC — posterior estimation of paracetamol PK parameters, Earth vs spaceflight
================================================================================

What this model is (and is not)
-------------------------------
* It is a POOLED Bayesian PK model fitted to published study-level summary
  statistics (mean Cmax, Tmax, AUC0-inf per study arm). One parameter set is
  estimated per condition (Earth ground control vs real spaceflight).
* It is NOT a hierarchical population-PK (PopPK) model: there are no
  subject-level random effects, because individual-subject data are not
  available. The next step is a hierarchical model with study- and
  subject-level random effects once individual data can be obtained.

Structural model: one-compartment, first-order oral absorption
---------------------------------------------------------------
    C(t)  = F·D·ka / (V·(ka−ke)) · (e^(−ke·t) − e^(−ka·t))
    tmax  = ln(ka/ke) / (ka − ke)            (→ 1/ke when ka → ke)
    Cmax  = C(tmax) = (F·D/V) · e^(−ke·tmax) (exact, algebraically equivalent)
    AUC∞  = F·D / (V·ke)

Units: D in mg, V in L  →  mg/L = µg/mL; AUC in µg·h/mL.

Identifiability note: Cmax, Tmax and AUC depend on F and V only through the
ratio V/F (apparent volume). F and V are therefore NOT separately identified
by these data — their individual posteriors are prior-driven. The quantities
the data actually inform are ka, ke, V/F and CL/F, and these are reported.

Assumptions (stated explicitly)
-------------------------------
1. Linear PK between 500 and 1000 mg (each observation is fitted at its own
   administered dose; Cmax and AUC scale with dose). Reasonable for
   paracetamol at therapeutic doses, but an assumption.
2. Saliva and plasma paracetamol concentrations are treated as interchangeable
   (saliva:plasma ratio ≈ 1 for paracetamol). Matrix is recorded per row.
3. Multiplicative (log-normal) residual error, which absorbs between-study
   heterogeneity (formulation, assay, sampling schedule).
4. Only AUC0-inf is used (AUC0-t is not the same quantity as F·D/(V·ke)).
5. Head-down-tilt bed rest (Gandia 2003 HDT days 1–80, Rumble, Renwick) is a
   ground ANALOG, not Earth baseline and not spaceflight — it is excluded from
   both groups here. Renwick 1992 is also excluded (nifedipine co-dosing).
6. Small data: 6 Earth / 4 space study arms. Results are preliminary and ka
   is weakly identifiable.
================================================================================
"""

import argparse
import numpy as np
import pymc as pm
import pytensor.tensor as pt
import arviz as az
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

# ── Priors ────────────────────────────────────────────────────────────────────
# ONE shared set of weakly-informative, textbook paracetamol priors is used for
# BOTH conditions. Any Earth–space difference in the posterior is therefore
# driven by the data, not by the prior.
#
# Previous version problems (fixed):
#  * The space F prior had mean 1.08 (=126.72% relative bioavailability × 0.85)
#    but F was truncated at 0.99, so most prior mass sat on the boundary.
#    126.72% (Kovachevich 2009, tablets) is a RELATIVE bioavailability ratio,
#    not an absolute F, and it conflicts with Polyakov 2021 (48.02%, reduced)
#    and Kovachevich capsules (no significant change). Any real change in
#    exposure is captured through V/F, which the data identify.
#  * Space priors were built from the same in-flight values used in the
#    likelihood (double use of data). Shared priors remove this.
F_MEAN, F_SD = 0.85, 0.08            # absolute oral bioavailability

def _beta_params(mean, sd):
    """Beta(α, β) with the given mean and SD — keeps F strictly inside (0, 1)."""
    k = mean * (1 - mean) / sd**2 - 1
    return mean * k, (1 - mean) * k

PRIORS = {
    'F_alpha': _beta_params(F_MEAN, F_SD)[0],
    'F_beta':  _beta_params(F_MEAN, F_SD)[1],
    'ka_mean': 1.50,  'ka_sd': 0.75,   # 1/h  (wide: ka weakly identifiable)
    'ke_mean': 0.277, 'ke_sd': 0.07,   # 1/h  (t½ ≈ 2.5 h)
    'Vd_mean': 67.5,  'Vd_sd': 15.0,   # L    (~0.9 L/kg × 75 kg)
}

# ── Observed study-arm data (traceable to data/paracetamol_pk_clean.csv) ──────
# Each row = one study arm. AUC = AUC0-inf (None if not reported).
# Units harmonised: Cmax µg/mL (Cintron 10780 ng/mL → 10.78), Tmax h
# (Rumble 35.6 min → 0.593 h).
#
# Previous version problems (fixed):
#  * Earth group contained HDT bed-rest days 1/18/80 (simulated microgravity),
#    plus session subsets AND their global means (double counting).
#  * Space group contained three GROUND-CONTROL values (Kovachevich 5.13,
#    Polyakov 11.5, Cintron 10.78 Cmax; Kovachevich 16.21 and Polyakov 45.5
#    AUC) mislabelled as spaceflight, and omitted the Kovachevich capsule arm.
#  * All observations were fitted at 500 mg although Gandia used 1000 mg and
#    Cintron/Polyakov used 650/625 mg. Each arm now uses its own dose.
STUDY_ARMS = [
    # source,                       condition, dose, Cmax,  Tmax,  AUCinf, matrix
    ('Gandia2003 HDT0 (ambulatory)',  'earth', 1000, 8.04,  2.50,  None,   'plasma'),
    ('Kovachevich2009 tablets',       'earth',  500, 5.13,  1.12, 16.21,   'saliva'),
    ('Kovachevich2009 capsules',      'earth',  500, 5.00,  0.90, 14.81,   'saliva'),
    ('Cintron1987a (DelloRusso2022)', 'earth',  650, 10.78, 0.50,  None,   'saliva'),
    ('Polyakov2021 (DelloRusso2022)', 'earth',  625, 11.50, 0.78, 45.50,   'saliva'),
    ('Rumble1991 ambulatory morning', 'earth',  500, 9.47,  0.593, None,   'plasma'),

    ('Kovachevich2009 tablets ISS',   'space',  500, 4.80,  1.80, 19.79,   'saliva'),
    ('Kovachevich2009 capsules ISS',  'space',  500, 4.17,  0.60, 17.23,   'saliva'),
    ('Cintron1987a Shuttle',          'space',  650, 11.40, 0.70,  None,   'saliva'),
    ('Polyakov2021 MIR/ISS',          'space',  625, 5.40,  1.80, 19.80,   'saliva'),
]


def get_condition_data(condition):
    rows = [r for r in STUDY_ARMS if r[1] == condition]
    dose = np.array([r[2] for r in rows], dtype=float)
    cmax = np.array([r[3] for r in rows], dtype=float)
    tmax = np.array([r[4] for r in rows], dtype=float)
    auc_rows = [r for r in rows if r[5] is not None]
    auc_dose = np.array([r[2] for r in auc_rows], dtype=float)
    auc = np.array([r[5] for r in auc_rows], dtype=float)
    return {'dose': dose, 'Cmax': cmax, 'Tmax': tmax,
            'AUC_dose': auc_dose, 'AUC': auc, 'n_arms': len(rows)}


# ── One-compartment oral PK (exact) ───────────────────────────────────────────
_EQ_TOL = 1e-6

def tmax_1cmt(ka, ke, lib=np):
    """tmax = ln(ka/ke)/(ka−ke); limit 1/ke when ka == ke (works for flip-flop)."""
    diff = ka - ke
    if lib is pt:
        near = pt.lt(pt.abs(diff), _EQ_TOL)
        safe = pt.switch(near, 1.0, diff)
        return pt.switch(near, 1.0 / ke, pt.log(ka / ke) / safe)
    diff = np.asarray(diff, dtype=float)
    near = np.abs(diff) < _EQ_TOL
    safe = np.where(near, 1.0, diff)
    return np.where(near, 1.0 / ke, np.log(ka / ke) / safe)


def cmax_1cmt(F, dose_mg, ka, ke, V, lib=np):
    """Exact one-compartment oral Cmax (µg/mL) = (F·D/V)·e^(−ke·tmax)."""
    tmax = tmax_1cmt(ka, ke, lib=lib)
    return (F * dose_mg / V) * lib.exp(-ke * tmax)


def auc_1cmt(F, dose_mg, ke, V):
    """AUC0-inf (µg·h/mL) = F·D / (V·ke)."""
    return F * dose_mg / (V * ke)


# ── Model ─────────────────────────────────────────────────────────────────────
def _fit_condition(data, n_draws, n_tune, seed):
    p = PRIORS
    with pm.Model():
        F  = pm.Beta('F', alpha=p['F_alpha'], beta=p['F_beta'])
        ka = pm.TruncatedNormal('ka', mu=p['ka_mean'], sigma=p['ka_sd'], lower=0.1)
        ke = pm.TruncatedNormal('ke', mu=p['ke_mean'], sigma=p['ke_sd'], lower=0.05)
        Vd = pm.TruncatedNormal('Vd', mu=p['Vd_mean'], sigma=p['Vd_sd'], lower=10.0)

        # Log-scale residual SDs (multiplicative error, ~CV)
        s_cmax = pm.HalfNormal('sigma_Cmax', sigma=0.5)
        s_tmax = pm.HalfNormal('sigma_Tmax', sigma=0.7)
        s_auc  = pm.HalfNormal('sigma_AUC',  sigma=0.5)

        cmax_pred = cmax_1cmt(F, data['dose'], ka, ke, Vd, lib=pt)
        tmax_pred = tmax_1cmt(ka, ke, lib=pt)
        auc_pred  = auc_1cmt(F, data['AUC_dose'], ke, Vd)

        pm.LogNormal('Cmax_obs', mu=pt.log(cmax_pred), sigma=s_cmax,
                     observed=data['Cmax'])
        pm.LogNormal('Tmax_obs', mu=pt.log(tmax_pred), sigma=s_tmax,
                     observed=data['Tmax'])
        pm.LogNormal('AUC_obs',  mu=pt.log(auc_pred),  sigma=s_auc,
                     observed=data['AUC'])

        pm.Deterministic('V_F',  Vd / F)          # apparent volume (identified)
        pm.Deterministic('CL_F', Vd * ke / F)     # apparent clearance (identified)

        trace = pm.sample(n_draws, tune=n_tune, chains=4, cores=1,
                          progressbar=False, target_accept=0.99,
                          return_inferencedata=True, random_seed=seed)
    return trace


REPORT_PARAMS = ['F', 'ka', 'ke', 'Vd', 'V_F', 'CL_F']


def _summarise(trace):
    out, post = {}, {}
    for prm in REPORT_PARAMS:
        x = trace.posterior[prm].values.flatten()
        post[prm] = x
        out[prm] = {'mean': x.mean(), 'sd': x.std(),
                    '2.5%': np.percentile(x, 2.5),
                    '97.5%': np.percentile(x, 97.5)}
    out['trace'] = trace
    out['posteriors'] = post
    out['diagnostics'] = az.summary(trace, var_names=['F', 'ka', 'ke', 'Vd'],
                                    kind='diagnostics')
    return out


def run_bayesian_popPK(drug='Paracetamol', n_draws=500, n_tune=300):
    """Pooled Bayesian PK estimation for Earth and spaceflight.
    (Function name kept for backward compatibility.)"""
    if drug != 'Paracetamol':
        raise ValueError('Only paracetamol has real in-flight PK data.')
    print(f"\n{'='*64}")
    print(f"  BAYESIAN PK ESTIMATION (pooled): {drug}")
    print(f"  n_draws={n_draws}, n_tune={n_tune}, chains=4")
    print(f"{'='*64}")

    results = {}
    for cond, seed in (('earth', 42), ('space', 43)):
        data = get_condition_data(cond)
        print(f"\n  Fitting {cond}: {data['n_arms']} study arms "
              f"({len(data['Cmax'])} Cmax, {len(data['Tmax'])} Tmax, "
              f"{len(data['AUC'])} AUC0-inf)")
        results[cond] = _summarise(_fit_condition(data, n_draws, n_tune, seed))

    # Convergence check
    for cond in ('earth', 'space'):
        d = results[cond]['diagnostics']
        rhat, ess = d['r_hat'].max(), d['ess_bulk'].min()
        flag = 'OK' if rhat < 1.01 and ess > 400 else 'CHECK'
        print(f"  {cond:<6} max R-hat={rhat:.3f}  min ESS_bulk={ess:.0f}  [{flag}]")

    print(f"\n  {'Param':<6} {'Earth mean':>11} {'Earth 95% CrI':>20} "
          f"{'Space mean':>11} {'Space 95% CrI':>20} {'Change':>8}")
    print("  " + "-"*80)
    for prm in REPORT_PARAMS:
        re_, rs = results['earth'][prm], results['space'][prm]
        chg = (rs['mean'] - re_['mean']) / re_['mean'] * 100
        ci_e = f"[{re_['2.5%']:.3f}, {re_['97.5%']:.3f}]"
        ci_s = f"[{rs['2.5%']:.3f}, {rs['97.5%']:.3f}]"
        print(f"  {prm:<6} {re_['mean']:>11.3f} {ci_e:>20} "
              f"{rs['mean']:>11.3f} {ci_s:>20} {chg:>+7.1f}%")
    print("  Note: F and Vd are prior-driven (only V/F is identified); "
          "interpret ka, ke, V/F, CL/F.")
    return results


def compute_bayesian_dose_ci(results, earth_dose_mg=500.0):
    """Space dose that matches Earth Cmax (primary) and AUC (secondary),
    using the exact one-compartment Cmax. Linear PK → factor = ratio of
    predicted exposures at the same dose."""
    e, s = results['earth']['posteriors'], results['space']['posteriors']
    n = min(len(e['F']), len(s['F']))
    rng = np.random.default_rng(0)
    ie, is_ = rng.permutation(len(e['F']))[:n], rng.permutation(len(s['F']))[:n]

    cmax_e = cmax_1cmt(e['F'][ie], earth_dose_mg, e['ka'][ie], e['ke'][ie], e['Vd'][ie])
    cmax_s = cmax_1cmt(s['F'][is_], earth_dose_mg, s['ka'][is_], s['ke'][is_], s['Vd'][is_])
    auc_e = auc_1cmt(e['F'][ie], earth_dose_mg, e['ke'][ie], e['Vd'][ie])
    auc_s = auc_1cmt(s['F'][is_], earth_dose_mg, s['ke'][is_], s['Vd'][is_])

    r_c, r_a = cmax_e / cmax_s, auc_e / auc_s
    q = lambda x, p: float(np.percentile(x, p))
    return {
        'earth_dose_mg': earth_dose_mg,
        'space_dose_mean_mg': round(float(np.median(r_c) * earth_dose_mg), 1),
        'space_dose_lo_mg': round(q(r_c, 2.5) * earth_dose_mg, 1),
        'space_dose_hi_mg': round(q(r_c, 97.5) * earth_dose_mg, 1),
        'factor_mean': round(float(np.median(r_c)), 3),
        'factor_lo': round(q(r_c, 2.5), 3),
        'factor_hi': round(q(r_c, 97.5), 3),
        'factor_auc_median': round(float(np.median(r_a)), 3),
        'factor_auc_lo': round(q(r_a, 2.5), 3),
        'factor_auc_hi': round(q(r_a, 97.5), 3),
    }


# ── Plot posteriors ───────────────────────────────────────────────────────────
def plot_posteriors(results, drug='Paracetamol', out_dir=None):
    from pathlib import Path
    if out_dir is None:
        out_dir = Path(__file__).resolve().parent.parent / "figures"
    out_dir = Path(out_dir)
    out_dir.mkdir(exist_ok=True)

    params = ['ka', 'ke', 'V_F', 'CL_F']
    labels = ['ka (1/h)', 'ke (1/h)', 'V/F (L)', 'CL/F (L/h)']

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    fig.patch.set_facecolor('#0a0a1a')
    for ax, prm, lab in zip(axes.flat, params, labels):
        ax.set_facecolor('#0d1117')
        e_post = results['earth']['posteriors'][prm]
        s_post = results['space']['posteriors'][prm]
        ax.hist(e_post, bins=40, alpha=0.6, color='#2196F3',
                label=f'Earth (median={np.median(e_post):.3f})', density=True)
        ax.hist(s_post, bins=40, alpha=0.6, color='#F44336',
                label=f'Space (median={np.median(s_post):.3f})', density=True)
        ax.axvline(np.median(e_post), color='#2196F3', lw=2, ls='--')
        ax.axvline(np.median(s_post), color='#F44336', lw=2, ls='--')
        ax.set(xlabel=lab, ylabel='Density', title=f'Posterior: {lab}')
        ax.legend(facecolor='#1a1a2e', labelcolor='white', fontsize=8)
        ax.tick_params(colors='white')
        ax.xaxis.label.set_color('white')
        ax.yaxis.label.set_color('white')
        ax.title.set_color('white')
        ax.grid(True, alpha=0.2)

    plt.suptitle(f'Bayesian PK (pooled) posteriors — {drug}: Earth vs Space\n'
                 'Preliminary: 6 Earth / 4 space study arms',
                 color='white', fontsize=12, fontweight='bold')
    plt.tight_layout()
    plt.savefig(out_dir / "bayesian_posteriors.png",
                dpi=150, bbox_inches='tight', facecolor='#0a0a1a')
    print("  Posterior plot saved")


# ── Self-check of the PK formulas ─────────────────────────────────────────────
def _self_check():
    from scipy.integrate import quad
    F, D, ka, ke, V = 0.85, 500.0, 1.5, 0.277, 67.5
    t = np.linspace(0, 48, 200001)
    C = F * D * ka / (V * (ka - ke)) * (np.exp(-ke * t) - np.exp(-ka * t))
    assert abs(cmax_1cmt(F, D, ka, ke, V) - C.max()) < 1e-3
    assert abs(tmax_1cmt(ka, ke) - t[C.argmax()]) < 1e-3
    auc_num = quad(lambda x: F*D*ka/(V*(ka-ke))*(np.exp(-ke*x)-np.exp(-ka*x)), 0, np.inf)[0]
    assert abs(auc_1cmt(F, D, ke, V) - auc_num) < 1e-6
    # ka == ke limit: Cmax = F·D/(V·e)
    assert abs(cmax_1cmt(F, D, 0.3, 0.3, V) - F * D / (V * np.e)) < 1e-9
    print(f"  Formula self-check passed (Earth-prior Cmax = "
          f"{cmax_1cmt(F, D, ka, ke, V):.2f} µg/mL)")


# ── Run ───────────────────────────────────────────────────────────────────────
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='SpacePK Bayesian PK (Layer 3)')
    parser.add_argument('--fast', action='store_true', help='Quick demo (1000 draws)')
    parser.add_argument('--full', action='store_true', help='Publication quality (3000 draws)')
    args = parser.parse_args()
    n_draws, n_tune = (3000, 2000) if args.full else (1000, 1000)

    print("SPACE PK PIPELINE — LAYER 3: BAYESIAN PK ESTIMATION (pooled)\n")
    _self_check()
    results = run_bayesian_popPK('Paracetamol', n_draws=n_draws, n_tune=n_tune)
    d = compute_bayesian_dose_ci(results)
    print("\n  MODEL-DERIVED DOSE EVALUATION (posterior median, 95% CrI):")
    print(f"    Earth dose            : {d['earth_dose_mg']:.0f} mg")
    print(f"    Space dose (Cmax-match): {d['space_dose_mean_mg']:.0f} mg "
          f"[{d['space_dose_lo_mg']:.0f}–{d['space_dose_hi_mg']:.0f}]")
    print(f"    Factor (Cmax, primary) : ×{d['factor_mean']:.2f} "
          f"[{d['factor_lo']:.2f}–{d['factor_hi']:.2f}]")
    print(f"    Factor (AUC, secondary): ×{d['factor_auc_median']:.2f} "
          f"[{d['factor_auc_lo']:.2f}–{d['factor_auc_hi']:.2f}]")
    plot_posteriors(results)
    print("\n✅ Layer 3 complete (preliminary — not for clinical use)")
