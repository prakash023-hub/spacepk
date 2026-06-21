"""
================================================================================
LAYER 3: BAYESIAN POPULATION PHARMACOKINETICS
SpacePK Framework — CPT:PSP Paper
PyMC Bayesian PopPK — Individual astronaut variability
Posterior space parameter estimation
================================================================================
"""

import argparse
import numpy as np
import pymc as pm
import pytensor.tensor as pt
from scipy.integrate import trapezoid
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')

# ── Literature data for Bayesian priors ───────────────────────────────────────
# From Gandia 2003, Kovachevich 2009, Polyakov 2021

LITERATURE_PRIORS = {
    'Paracetamol': {
        # Earth priors (mean, sigma) — Gandia 2003 Table II
        'earth': {
            'Cmax_mean': 8.04, 'Cmax_sd': 2.72,  # µg/mL
            'Tmax_mean': 2.50, 'Tmax_sd': 1.20,  # h
            'AUC_mean': 29.17, 'AUC_sd': 10.27,  # µg·h/mL
            't12_mean': 3.17,  't12_sd': 1.58,   # h
            'F_mean': 0.85,    'F_sd': 0.08,
            'ka_mean': 1.50,   'ka_sd': 0.50,
            'ke_mean': 0.277,  'ke_sd': 0.05,
            'Vd_mean': 67.5,   'Vd_sd': 15.0,
        },
        # Space priors (Kovachevich 2009, Polyakov 2021)
        'space': {
            'Cmax_mean': 5.40, 'Cmax_sd': 1.20,
            'Tmax_mean': 1.80, 'Tmax_sd': 0.64,
            'AUC_mean': 19.79, 'AUC_sd': 3.15,
            't12_mean': 3.24,  't12_sd': 1.04,
            'F_mean': 1.08,    'F_sd': 0.24,    # 126.72% × 0.85
            'ka_mean': 0.90,   'ka_sd': 0.40,
            'ke_mean': 0.266,  'ke_sd': 0.06,
            'Vd_mean': 65.5,   'Vd_sd': 14.0,
        }
    }
}

# ── Observed data (from our dataset) ─────────────────────────────────────────
OBSERVED_DATA = {
    'Paracetamol_Earth': {
        'Cmax': [6.24, 9.84, 8.04, 9.01, 12.21, 10.43, 11.81, 11.95, 11.89,
                 15.48, 13.69, 14.74],
        'Tmax': [3.11, 1.89, 2.50, 1.47, 1.33, 1.41, 0.82, 1.25, 1.05,
                 1.22, 0.48, 0.91],
        'AUC':  [21.18, 37.16, 29.17, 26.68, 39.77, 32.50, 32.22, 38.04,
                 35.33, 29.42, 41.12, 34.24],
        'n': 18,
        'condition': 'Earth/Simulated'
    },
    'Paracetamol_Space': {
        'Cmax': [5.13, 4.80, 11.5, 5.40, 10.78, 11.40],
        'Tmax': [1.12, 1.80, 0.78, 1.80, 0.50, 0.70],
        'AUC':  [16.21, 19.79, 45.5, 19.8],
        'n': 5,
        'condition': 'Real spaceflight'
    }
}

# ── Bayesian PopPK Model ───────────────────────────────────────────────────────
def run_bayesian_popPK(drug='Paracetamol', n_draws=500, n_tune=300):
    """
    Bayesian Population PK estimation
    Estimates posterior distributions of PK parameters
    for Earth and Space conditions
    """
    print(f"\n{'='*60}")
    print(f"  BAYESIAN POPUPK: {drug}")
    print(f"  n_draws={n_draws}, n_tune={n_tune}")
    print(f"{'='*60}")

    priors = LITERATURE_PRIORS[drug]
    obs    = OBSERVED_DATA

    earth_Cmax = np.array(obs[f'{drug}_Earth']['Cmax'])
    earth_Tmax = np.array(obs[f'{drug}_Earth']['Tmax'])
    earth_AUC  = np.array(obs[f'{drug}_Earth']['AUC'])
    space_Cmax = np.array(obs[f'{drug}_Space']['Cmax'])
    space_AUC  = np.array(obs[f'{drug}_Space']['AUC'])

    results = {}

    # ── Earth Model ───────────────────────────────────────────────────────────
    print("\n  Running Earth Bayesian model...")
    with pm.Model() as earth_model:
        # Priors — based on literature
        F_e  = pm.TruncatedNormal('F',
               mu=priors['earth']['F_mean'],
               sigma=priors['earth']['F_sd'],
               lower=0.1, upper=0.99)
        ka_e = pm.TruncatedNormal('ka',
               mu=priors['earth']['ka_mean'],
               sigma=priors['earth']['ka_sd'],
               lower=0.1)
        ke_e = pm.TruncatedNormal('ke',
               mu=priors['earth']['ke_mean'],
               sigma=priors['earth']['ke_sd'],
               lower=0.01)
        Vd_e = pm.TruncatedNormal('Vd',
               mu=priors['earth']['Vd_mean'],
               sigma=priors['earth']['Vd_sd'],
               lower=10.0)
        sigma_Cmax = pm.HalfNormal('sigma_Cmax', sigma=3.0)
        sigma_AUC  = pm.HalfNormal('sigma_AUC', sigma=10.0)

        # Predicted Cmax = F*D*ka / (Vd*(ka-ke)) * ...
        # Simplified: Cmax_pred = F * dose / Vd (at peak)
        dose = 500.0
        Cmax_pred = (F_e * dose * 1000 * ka_e) / \
                    (Vd_e * 1000 * (ka_e - ke_e + 1e-6))
        AUC_pred  = (F_e * dose * 1000) / (Vd_e * 1000 * ke_e)

        # Likelihood
        pm.Normal('Cmax_obs', mu=Cmax_pred, sigma=sigma_Cmax,
                  observed=earth_Cmax)
        pm.Normal('AUC_obs', mu=AUC_pred, sigma=sigma_AUC,
                  observed=earth_AUC)

        # Sample
        trace_e = pm.sample(n_draws, tune=n_tune, progressbar=False,
                            target_accept=0.9, return_inferencedata=True,
                            random_seed=42)

    # Extract posteriors
    F_post_e  = trace_e.posterior['F'].values.flatten()
    ka_post_e = trace_e.posterior['ka'].values.flatten()
    ke_post_e = trace_e.posterior['ke'].values.flatten()
    Vd_post_e = trace_e.posterior['Vd'].values.flatten()

    results['earth'] = {
        'F':  {'mean': F_post_e.mean(),  'sd': F_post_e.std(),
               '2.5%': np.percentile(F_post_e, 2.5),
               '97.5%': np.percentile(F_post_e, 97.5)},
        'ka': {'mean': ka_post_e.mean(), 'sd': ka_post_e.std(),
               '2.5%': np.percentile(ka_post_e, 2.5),
               '97.5%': np.percentile(ka_post_e, 97.5)},
        'ke': {'mean': ke_post_e.mean(), 'sd': ke_post_e.std(),
               '2.5%': np.percentile(ke_post_e, 2.5),
               '97.5%': np.percentile(ke_post_e, 97.5)},
        'Vd': {'mean': Vd_post_e.mean(), 'sd': Vd_post_e.std(),
               '2.5%': np.percentile(Vd_post_e, 2.5),
               '97.5%': np.percentile(Vd_post_e, 97.5)},
        'trace': trace_e,
        'posteriors': {'F': F_post_e, 'ka': ka_post_e,
                       'ke': ke_post_e, 'Vd': Vd_post_e}
    }

    # ── Space Model ───────────────────────────────────────────────────────────
    print("  Running Space Bayesian model...")
    with pm.Model() as space_model:
        F_s  = pm.TruncatedNormal('F',
               mu=priors['space']['F_mean'],
               sigma=priors['space']['F_sd'],
               lower=0.1, upper=0.99)
        ka_s = pm.TruncatedNormal('ka',
               mu=priors['space']['ka_mean'],
               sigma=priors['space']['ka_sd'],
               lower=0.1)
        ke_s = pm.TruncatedNormal('ke',
               mu=priors['space']['ke_mean'],
               sigma=priors['space']['ke_sd'],
               lower=0.01)
        Vd_s = pm.TruncatedNormal('Vd',
               mu=priors['space']['Vd_mean'],
               sigma=priors['space']['Vd_sd'],
               lower=10.0)
        sigma_Cmax = pm.HalfNormal('sigma_Cmax', sigma=2.0)
        sigma_AUC  = pm.HalfNormal('sigma_AUC', sigma=8.0)

        Cmax_pred_s = (F_s * 500.0 * 1000 * ka_s) / \
                      (Vd_s * 1000 * (ka_s - ke_s + 1e-6))
        AUC_pred_s  = (F_s * 500.0 * 1000) / (Vd_s * 1000 * ke_s)

        pm.Normal('Cmax_obs', mu=Cmax_pred_s, sigma=sigma_Cmax,
                  observed=space_Cmax)
        pm.Normal('AUC_obs', mu=AUC_pred_s, sigma=sigma_AUC,
                  observed=space_AUC)

        trace_s = pm.sample(n_draws, tune=n_tune, progressbar=False,
                            target_accept=0.9, return_inferencedata=True,
                            random_seed=42)

    F_post_s  = trace_s.posterior['F'].values.flatten()
    ka_post_s = trace_s.posterior['ka'].values.flatten()
    ke_post_s = trace_s.posterior['ke'].values.flatten()
    Vd_post_s = trace_s.posterior['Vd'].values.flatten()

    results['space'] = {
        'F':  {'mean': F_post_s.mean(),  'sd': F_post_s.std(),
               '2.5%': np.percentile(F_post_s, 2.5),
               '97.5%': np.percentile(F_post_s, 97.5)},
        'ka': {'mean': ka_post_s.mean(), 'sd': ka_post_s.std(),
               '2.5%': np.percentile(ka_post_s, 2.5),
               '97.5%': np.percentile(ka_post_s, 97.5)},
        'ke': {'mean': ke_post_s.mean(), 'sd': ke_post_s.std(),
               '2.5%': np.percentile(ke_post_s, 2.5),
               '97.5%': np.percentile(ke_post_s, 97.5)},
        'Vd': {'mean': Vd_post_s.mean(), 'sd': Vd_post_s.std(),
               '2.5%': np.percentile(Vd_post_s, 2.5),
               '97.5%': np.percentile(Vd_post_s, 97.5)},
        'trace': trace_s,
        'posteriors': {'F': F_post_s, 'ka': ka_post_s,
                       'ke': ke_post_s, 'Vd': Vd_post_s}
    }

    # ── Print results ─────────────────────────────────────────────────────────
    print(f"\n  {'Parameter':<10} {'Earth Mean':>12} {'Earth 95%CI':>18} "
          f"{'Space Mean':>12} {'Space 95%CI':>18} {'Change':>8}")
    print("  " + "-"*82)
    for param in ['F', 'ka', 'ke', 'Vd']:
        re = results['earth'][param]
        rs = results['space'][param]
        chg = ((rs['mean'] - re['mean']) / re['mean']) * 100
        ci_e = f"[{re['2.5%']:.3f},{re['97.5%']:.3f}]"
        ci_s = f"[{rs['2.5%']:.3f},{rs['97.5%']:.3f}]"
        print(f"  {param:<10} {re['mean']:>12.4f} {ci_e:>18} "
              f"{rs['mean']:>12.4f} {ci_s:>18} {chg:>+7.1f}%")

    return results


def compute_bayesian_dose_ci(results, earth_dose_mg=500.0):
    """Posterior dose factor to match Earth Cmax in space (ratio of predicted Cmax)."""
    e = results['earth']['posteriors']
    s = results['space']['posteriors']
    dose = earth_dose_mg
    cmax_e = (e['F'] * dose * e['ka']) / (e['Vd'] * (e['ka'] - e['ke'] + 1e-6))
    cmax_s = (s['F'] * dose * s['ka']) / (s['Vd'] * (s['ka'] - s['ke'] + 1e-6))
    ratio = cmax_e / (cmax_s + 1e-6)
    return {
        'earth_dose_mg': earth_dose_mg,
        'space_dose_mean_mg': round(float(ratio.mean() * earth_dose_mg), 1),
        'space_dose_lo_mg': round(float(np.percentile(ratio, 2.5) * earth_dose_mg), 1),
        'space_dose_hi_mg': round(float(np.percentile(ratio, 97.5) * earth_dose_mg), 1),
        'factor_mean': round(float(ratio.mean()), 3),
        'factor_lo': round(float(np.percentile(ratio, 2.5)), 3),
        'factor_hi': round(float(np.percentile(ratio, 97.5)), 3),
    }

# ── Plot posteriors ────────────────────────────────────────────────────────────
def plot_posteriors(results, drug='Paracetamol', out_dir=None):
    from pathlib import Path
    if out_dir is None:
        out_dir = Path(__file__).resolve().parent.parent / "figures"
    out_dir = Path(out_dir)
    out_dir.mkdir(exist_ok=True)

    params = ['F', 'ka', 'ke', 'Vd']
    units  = ['(fraction)', '(/h)', '(/h)', '(L)']

    fig, axes = plt.subplots(2, 2, figsize=(12, 8))
    fig.patch.set_facecolor('#0a0a1a')

    for ax, param, unit in zip(axes.flat, params, units):
        ax.set_facecolor('#0d1117')
        e_post = results['earth']['posteriors'][param]
        s_post = results['space']['posteriors'][param]

        ax.hist(e_post, bins=40, alpha=0.6, color='#2196F3',
                label=f'Earth (μ={e_post.mean():.3f})', density=True)
        ax.hist(s_post, bins=40, alpha=0.6, color='#F44336',
                label=f'Space (μ={s_post.mean():.3f})', density=True)
        ax.axvline(e_post.mean(), color='#2196F3', lw=2, ls='--')
        ax.axvline(s_post.mean(), color='#F44336', lw=2, ls='--')
        ax.set(xlabel=f'{param} {unit}', ylabel='Density',
               title=f'Posterior: {param}')
        ax.legend(facecolor='#1a1a2e', labelcolor='white', fontsize=8)
        ax.tick_params(colors='white')
        ax.xaxis.label.set_color('white')
        ax.yaxis.label.set_color('white')
        ax.title.set_color('white')
        ax.grid(True, alpha=0.2)

    plt.suptitle(f'Bayesian PopPK Posterior Distributions — {drug}\nEarth vs Space',
                 color='white', fontsize=12, fontweight='bold')
    plt.tight_layout()
    plt.savefig(out_dir / "bayesian_posteriors.png",
                dpi=150, bbox_inches='tight', facecolor='#0a0a1a')
    print("  Posterior plot saved")

# ── Run ────────────────────────────────────────────────────────────────────────
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='SpacePK Bayesian PopPK')
    parser.add_argument('--fast', action='store_true', help='Quick demo (500 draws)')
    parser.add_argument('--full', action='store_true', help='Publication quality (2000 draws)')
    args = parser.parse_args()
    n_draws, n_tune = (2000, 1000) if args.full else (500, 300)

    print("SPACE PK PIPELINE — LAYER 3: BAYESIAN POPUPK")
    print("CPT: Pharmacometrics & Systems Pharmacology\n")
    results = run_bayesian_popPK('Paracetamol', n_draws=n_draws, n_tune=n_tune)
    dose_ci = compute_bayesian_dose_ci(results)
    print(f"\n  BAYESIAN DOSE RECOMMENDATION (95% CI):")
    print(f"    Earth dose     : {dose_ci['earth_dose_mg']:.0f} mg")
    print(f"    Space dose     : {dose_ci['space_dose_mean_mg']:.0f} mg "
          f"[{dose_ci['space_dose_lo_mg']:.0f}–{dose_ci['space_dose_hi_mg']:.0f}]")
    print(f"    Factor         : ×{dose_ci['factor_mean']:.2f} "
          f"[{dose_ci['factor_lo']:.2f}–{dose_ci['factor_hi']:.2f}]")
    plot_posteriors(results)
    print("\n✅ Layer 3 complete")
