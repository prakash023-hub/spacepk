"""
================================================================================
SpacePK — SINGLE SOURCE OF TRUTH
================================================================================
Every figure, table, and the Streamlit app import the PK engine from HERE.
Nothing else in the repo may define pbpk_odes / extract_pk / recommend_dose.
This is what guarantees the app, the figures, and the manuscript can never
disagree again.

Contents
  1. get_physiology()  — 7-compartment volumes/flows + mission-phase modifiers
  2. calibrate_kp()    — tissue:plasma Kp calibrated so Vss ≈ literature Vd
                         (Rodgers & Rowland / Poulin & Theil practice)
  3. pbpk_odes()       — 7-compartment ODE system (GI→Portal→Liver→Venous→
                         Arterial→Tissue→Kidney). Literature F already includes
                         first-pass, so absorbed drug enters VENOUS blood
                         (no double first-pass extraction).
  4. extract_pk()      — Cmax/Tmax/AUC/t½ (NO rounding here; round at display)
  5. recommend_dose()  — ONE primary number, criterion='cmax' by default;
                         AUC-matched reported as secondary sensitivity only.
  6. run_pbpk_analysis()— Earth vs Space simulation + dose recommendation.

Validation (paracetamol, Earth, dose-normalised anchors — see dataset_census):
  Cmax ✓ (<10%), t½ ✓ (<6%), AUC within AUC0t→AUC0inf range.
  Tmax runs fast (single first-order absorption) — disclosed limitation.

NEVER rescale a simulated curve to a literature value. If validation fails,
fix the model — do not normalise the output.
================================================================================
"""

from __future__ import annotations

import numpy as np
from scipy.integrate import solve_ivp, trapezoid


# ── 1. Physiology + mission-phase modifiers ──────────────────────────────────
def get_physiology(body_weight=75, mission_days=0):
    physio = {
        'Q_portal': 72.0, 'Q_liver': 90.0, 'Q_kidney': 72.0, 'Q_tissue': 228.0,
        'V_gi': 1.65, 'V_portal': 0.99, 'V_liver': 1.65, 'V_arterial': 1.65,
        'V_venous': 3.30, 'V_tissue': 22.5, 'V_kidney': 0.28,
        'body_weight': body_weight, 'plasma_volume': 3.0,
        'total_body_water': body_weight * 0.60,
        'mission_days': mission_days,
        'gravity': 9.8 if mission_days == 0 else 0.0,
    }
    if mission_days > 0:
        if mission_days <= 3:
            physio['Q_portal'] *= 1.15
            physio['Q_liver'] *= 1.10
            physio['V_tissue'] *= 0.99
            physio['total_body_water'] *= 0.99
            physio['phase'] = 'acute'
        elif mission_days <= 14:
            physio['Q_liver'] *= 1.05
            physio['V_tissue'] *= 0.98
            physio['plasma_volume'] *= 0.97
            physio['phase'] = 'adaptation'
        else:
            physio['Q_liver'] *= 0.95
            physio['Q_kidney'] *= 0.98
            physio['V_tissue'] *= 0.97
            physio['plasma_volume'] *= 0.95
            physio['total_body_water'] *= 0.97
            physio['phase'] = 'chronic'
    else:
        physio['phase'] = 'earth'
    return physio


# ── 2. Vd-calibrated tissue:plasma partition ─────────────────────────────────
def calibrate_kp(Vd, physio, logP=0.0):
    """Choose Kp so Vss ≈ literature Vd: Vss ≈ V_blood + Kp·V_perfused."""
    V_blood = physio['V_portal'] + physio['V_arterial'] + physio['V_venous']
    V_perf = physio['V_liver'] + physio['V_tissue'] + physio['V_kidney']
    kp_guess = 1.0 + max(logP, -1.0) * 0.15
    kp_vd = (max(Vd, V_blood + 0.5) - V_blood) / max(V_perf, 0.1)
    Kp = 0.85 * kp_vd + 0.15 * kp_guess
    return float(np.clip(Kp, 0.05, 80.0))


# ── 3. 7-compartment ODE system ──────────────────────────────────────────────
def pbpk_odes(t, y, drug_params, physio):
    """
    State (conc mg/L except A_gi, A_excr in mg):
      y0 A_gi, y1 C_portal, y2 C_liver, y3 C_arterial,
      y4 C_venous, y5 C_tissue, y6 C_kidney, y7 A_excreted
    """
    A_gi, C_portal, C_liver, C_arterial, C_venous, C_tissue, C_kidney, A_excr = y

    ka = max(drug_params['ka'], 1e-8)
    F = float(np.clip(drug_params['F'], 0.0, 1.0))
    CLh = max(drug_params['CLh'], 0.0)
    CLr = max(drug_params['CLr'], 0.0)
    Kp = max(drug_params['Kp'], 0.05)
    Kp_l = max(drug_params.get('Kp_liver', Kp * 0.9), 0.05)
    Kp_k = max(drug_params.get('Kp_kidney', Kp * 0.5), 0.05)

    Q_p = physio['Q_portal']
    Q_l = physio['Q_liver']
    Q_k = physio['Q_kidney']
    Q_t = physio['Q_tissue']
    V_po = physio['V_portal']
    V_l = physio['V_liver']
    V_a = physio['V_arterial']
    V_v = physio['V_venous']
    V_ti = physio['V_tissue']
    V_k = physio['V_kidney']

    # Literature F already includes first-pass loss (F_abs·F_g·F_h).
    # Absorbed drug enters systemic VENOUS blood so hepatic extraction is NOT
    # applied a second time. Portal/liver handle recirculating drug only.
    dA_gi = -ka * A_gi
    absorption_rate = ka * A_gi * F

    dC_portal = (Q_p * C_arterial - Q_p * C_portal) / V_po

    Q_ha = max(Q_l - Q_p, 1.0)
    liver_in = Q_p * C_portal + Q_ha * C_arterial
    liver_out = Q_l * (C_liver / Kp_l)
    hepatic_clearance = CLh * (C_liver / Kp_l)
    dC_liver = (liver_in - liver_out - hepatic_clearance) / V_l

    dC_tissue = (Q_t * C_arterial - Q_t * (C_tissue / Kp)) / V_ti

    renal_clearance = CLr * (C_kidney / Kp_k)
    dC_kidney = (Q_k * C_arterial - Q_k * (C_kidney / Kp_k) - renal_clearance) / V_k

    venous_in = (
        Q_t * (C_tissue / Kp)
        + Q_k * (C_kidney / Kp_k)
        + Q_l * (C_liver / Kp_l)
        + absorption_rate
    )
    Q_return = Q_t + Q_k + Q_l
    dC_venous = (venous_in - Q_return * C_venous) / V_v

    dC_arterial = (Q_return * C_venous - Q_return * C_arterial) / V_a
    dA_excr = renal_clearance

    return [dA_gi, dC_portal, dC_liver, dC_arterial, dC_venous,
            dC_tissue, dC_kidney, dA_excr]


# ── 4. PK metric extraction (NO rounding here) ───────────────────────────────
def extract_pk(t, C):
    """Return raw floats. Rounding happens ONLY at the display layer.
    (Rounding here previously collapsed low-conc drugs like scopolamine.)"""
    C = np.clip(np.asarray(C, dtype=float), 0, None)
    if C.max() <= 0:
        return {'Cmax': 0.0, 'Tmax': 0.0, 'AUC': 0.0, 't12': float('nan')}
    idx = int(np.argmax(C))
    Cmax = float(C[idx])
    Tmax = float(t[idx])
    AUC = float(trapezoid(C, t))
    elim, t_el = C[idx:], t[idx:]
    mask = elim > Cmax * 0.05
    t12 = float('nan')
    if mask.sum() > 5:
        try:
            sl, _ = np.polyfit(t_el[mask], np.log(elim[mask] + 1e-12), 1)
            if sl < 0:
                t12 = 0.693 / (-sl)
        except Exception:
            pass
    return {'Cmax': Cmax, 'Tmax': Tmax, 'AUC': AUC, 't12': t12}


# ── 5. Dose recommendation (ONE primary number) ──────────────────────────────
def recommend_dose(drug_profile, pk_earth, pk_space, criterion='cmax'):
    """
    Primary criterion defaults to 'cmax' (preserve peak exposure).
    AUC-matched dose is returned as a SECONDARY sensitivity value only.
    """
    earth_dose = drug_profile['dose_mg']
    cmax_e = pk_earth.get('Cmax', 0) or 0
    cmax_s = pk_space.get('Cmax', 0) or 0
    auc_e = pk_earth.get('AUC', 0) or 0
    auc_s = pk_space.get('AUC', 0) or 0

    if drug_profile.get('is_iv'):
        ke_e = drug_profile['ke_earth']
        ke_s = max(drug_profile['ke_space'], 1e-6)
        factor = ke_e / ke_s
        return {
            'earth_dose_mg': earth_dose,
            'space_dose_mg': round(earth_dose * factor, 1),
            'space_dose_auc_mg': round(earth_dose * factor, 1),
            'adjustment_factor': round(factor, 4),
            'adjustment_pct': round((factor - 1) * 100, 1),
            'criterion': 'IV clearance-matched (ke)',
            'rationale': 'IV: clearance-matched dose (ke Earth / ke Space)',
        }

    if cmax_e <= 0 or cmax_s <= 0:
        return None

    factor_cmax = cmax_e / cmax_s
    factor_auc = (auc_e / auc_s) if auc_s > 0 else factor_cmax
    primary = factor_cmax if criterion == 'cmax' else factor_auc
    return {
        'earth_dose_mg': earth_dose,
        'space_dose_mg': round(earth_dose * primary, 1),
        'space_dose_auc_mg': round(earth_dose * factor_auc, 1),
        'adjustment_factor': round(primary, 4),
        'adjustment_pct': round((primary - 1) * 100, 1),
        'factor_auc': round(factor_auc, 4),
        'factor_cmax': round(factor_cmax, 4),
        'criterion': f'{criterion.upper()}-matched (primary)',
        'rationale': (
            f'{criterion.upper()}-matched dose to preserve target exposure '
            f'(AUC-matched sensitivity: {earth_dose * factor_auc:.1f} mg)'
        ),
    }


# Backward-compatible alias
compute_dose_recommendation = recommend_dose


# ── 6. Simulation driver ─────────────────────────────────────────────────────
def _build_params(drug_profile, physio, earth=True):
    prefix = 'earth' if earth else 'space'
    ke = max(drug_profile[f'ke_{prefix}'], 1e-6)
    ka = drug_profile.get(f'ka_{prefix}')
    if ka is None:
        ka = drug_profile.get('ka_earth') or 1.5
        if not earth and drug_profile.get('space_modifiers'):
            ka = ka * drug_profile['space_modifiers'].get('ka_factor', 1.0)
    Vd = max(drug_profile[f'Vd_{prefix}'], 1.0)
    F = drug_profile[f'F_{prefix}']
    CL = ke * Vd
    Kp = calibrate_kp(Vd, physio, logP=drug_profile.get('logP', 0.0))
    return {
        'drug_name': drug_profile['drug_name'],
        'dose_mg': drug_profile['dose_mg'],
        'F': F, 'Vd': Vd, 'ka': max(ka, 1e-6), 'ke': ke,
        'CLh': CL * 0.7, 'CLr': CL * 0.3,
        'Kp': Kp, 'Kp_liver': Kp * 0.9, 'Kp_kidney': max(Kp * 0.4, 0.05),
        'route': 'iv' if drug_profile.get('is_iv') else 'oral',
    }


def simulate_pbpk(params_e, params_s, physio_e, physio_s, t_end=24):
    t_eval = np.linspace(0, t_end, 600)
    dose = params_e['dose_mg']

    def y0(params, physio):
        if params.get('route') == 'iv':
            return [0, 0, 0, 0, dose / max(physio['V_venous'], 0.1), 0, 0, 0]
        return [dose, 0, 0, 0, 0, 0, 0, 0]

    sol_e = solve_ivp(
        pbpk_odes, (0, t_end), y0(params_e, physio_e), t_eval=t_eval,
        args=(params_e, physio_e), method='RK45', rtol=1e-6, atol=1e-9,
    )
    sol_s = solve_ivp(
        pbpk_odes, (0, t_end), y0(params_s, physio_s), t_eval=t_eval,
        args=(params_s, physio_s), method='RK45', rtol=1e-6, atol=1e-9,
    )
    return t_eval, sol_e, sol_s


def run_pbpk_analysis(drug_profile, mission_days=30, body_weight=75,
                      criterion='cmax', verbose=True):
    t_end = max(24.0, float(drug_profile.get('t12_earth', 2) or 2) * 6)
    physio_e = get_physiology(body_weight, 0)
    physio_s = get_physiology(body_weight, mission_days)
    params_e = _build_params(drug_profile, physio_e, earth=True)
    params_s = _build_params(drug_profile, physio_s, earth=False)

    t, sol_e, sol_s = simulate_pbpk(params_e, params_s, physio_e, physio_s, t_end=t_end)

    C_earth = sol_e.y[3]
    C_space = sol_s.y[3]
    pk_earth = extract_pk(t, C_earth)
    pk_space = extract_pk(t, C_space)
    dose_rec = recommend_dose(drug_profile, pk_earth, pk_space, criterion=criterion)

    if verbose:
        print(f"\n{'='*60}")
        print(f"  7-COMPARTMENT PBPK: {drug_profile['drug_name']}")
        print(f"  Mission day {mission_days} | Kp_earth={params_e['Kp']:.2f} Kp_space={params_s['Kp']:.2f}")
        print(f"{'='*60}")
        print(f"  Earth: F={params_e['F']:.3f} ka={params_e['ka']:.3f} ke={params_e['ke']:.4f} Vd={params_e['Vd']:.1f}L")
        print(f"  Space: F={params_s['F']:.3f} ka={params_s['ka']:.3f} ke={params_s['ke']:.4f} Vd={params_s['Vd']:.1f}L")
        for param in ['Cmax', 'Tmax', 'AUC', 't12']:
            e, s = pk_earth[param], pk_space[param]
            chg = f"{((s-e)/e)*100:+.1f}%" if isinstance(e, float) and e > 0 else 'NA'
            print(f"    {param}: Earth={e:.4g}  Space={s:.4g}  ({chg})")
        if dose_rec:
            print(f"  Dose {dose_rec['criterion']}: {dose_rec['earth_dose_mg']:.0f} → "
                  f"{dose_rec['space_dose_mg']:.0f} mg ({dose_rec['adjustment_pct']:+.1f}%)")
            print(f"  Dose AUC-matched (sensitivity): {dose_rec['earth_dose_mg']:.0f} → "
                  f"{dose_rec['space_dose_auc_mg']:.0f} mg")

    return {
        't': t,
        'C_earth': C_earth, 'C_space': C_space,
        'C_liver_e': sol_e.y[2], 'C_liver_s': sol_s.y[2],
        'C_tissue_e': sol_e.y[5], 'C_tissue_s': sol_s.y[5],
        'pk_earth': pk_earth, 'pk_space': pk_space,
        'dose_recommendation': dose_rec,
        'drug_params_earth': params_e,
        'drug_params_space': params_s,
        'physio_earth': physio_e,
        'physio_space': physio_s,
        'model': '7-compartment PBPK (Vd-calibrated Kp)',
    }


# Backward-compatible alias
extract_pk_pbpk = extract_pk
