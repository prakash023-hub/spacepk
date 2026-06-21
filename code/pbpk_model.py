"""
================================================================================
LAYER 2: PBPK 7-COMPARTMENT WHOLE BODY MODEL
SpacePK Framework — CPT:PSP Paper
Physiologically-Based Pharmacokinetic Model with Space Modifications
================================================================================
Compartments:
    1. GI tract (absorption site)
    2. Portal blood
    3. Liver (first-pass + metabolism)
    4. Arterial blood
    5. Venous blood
    6. Peripheral tissue
    7. Kidney (excretion)
================================================================================
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.integrate import trapezoid
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# ── PHYSIOLOGICAL PARAMETERS ──────────────────────────────────────────────────
# Source: Brown et al. 1997, Williams & Leggett 1981
# Space modifications: Leach 1981, Dello Russo 2022

def get_physiology(body_weight=75, mission_days=0):
    """
    Returns physiological parameters
    mission_days=0 → Earth; >0 → Space modified
    """
    # Earth baseline (70-75kg adult)
    physio = {
        # Blood flows (L/h) — cardiac output distributed
        'Q_total': 390,       # Total cardiac output L/h
        'Q_portal': 72,       # Portal blood flow L/h
        'Q_liver': 90,        # Hepatic blood flow L/h (portal + arterial)
        'Q_kidney': 72,       # Renal blood flow L/h
        'Q_tissue': 228,      # Peripheral tissue flow L/h

        # Compartment volumes (L)
        'V_gi': 1.65,         # GI lumen volume
        'V_portal': 0.99,     # Portal blood volume
        'V_liver': 1.65,      # Liver volume
        'V_arterial': 1.65,   # Arterial blood
        'V_venous': 3.30,     # Venous blood
        'V_tissue': 22.5,     # Peripheral tissue (muscle+fat)
        'V_kidney': 0.28,     # Kidney volume

        # Body composition
        'body_weight': body_weight,
        'plasma_volume': 3.0, # L
        'total_body_water': body_weight * 0.60,  # L
        
        # Mission parameters
        'mission_days': mission_days,
        'gravity': 9.8 if mission_days == 0 else 0.0,
    }

    # ── Space modifications ───────────────────────────────────────────────────
    if mission_days > 0:
        if mission_days <= 3:
            # Acute phase: headward fluid shift
            physio['Q_portal']  *= 1.15  # increased GI blood flow (Gandia 2003)
            physio['Q_liver']   *= 1.10
            physio['V_tissue']  *= 0.99
            physio['total_body_water'] *= 0.99
        elif mission_days <= 14:
            # Adaptation phase
            physio['Q_liver']   *= 1.05
            physio['V_tissue']  *= 0.98
            physio['plasma_volume'] *= 0.97
        else:
            # Chronic phase (Leach 1981, Dello Russo 2022)
            physio['Q_liver']   *= 0.95  # reduced hepatic flow
            physio['Q_kidney']  *= 0.98
            physio['V_tissue']  *= 0.97
            physio['plasma_volume'] *= 0.95  # -5% plasma volume
            physio['total_body_water'] *= 0.97  # Leach 1981: -3%

    return physio

# ── PBPK ODE SYSTEM ───────────────────────────────────────────────────────────
def pbpk_odes(t, y, drug_params, physio, route='oral'):
    """
    7-compartment PBPK differential equations
    
    State vector y:
        y[0]: A_gi      — amount in GI tract (mg)
        y[1]: C_portal  — concentration in portal blood (mg/L)
        y[2]: C_liver   — concentration in liver (mg/L)
        y[3]: C_arterial— concentration in arterial blood (mg/L)
        y[4]: C_venous  — concentration in venous blood (mg/L)
        y[5]: C_tissue  — concentration in peripheral tissue (mg/L)
        y[6]: C_kidney  — concentration in kidney (mg/L)
        y[7]: A_excreted— cumulative amount excreted (mg)
    """
    A_gi, C_portal, C_liver, C_arterial, C_venous, C_tissue, C_kidney, A_excr = y

    # Drug parameters
    ka   = drug_params['ka']       # absorption rate constant (1/h)
    F    = drug_params['F']        # bioavailability
    Vd   = drug_params['Vd']       # Vd in L
    ke   = drug_params['ke']       # elimination rate constant (1/h)
    CLh  = drug_params['CLh']      # hepatic clearance (L/h)
    CLr  = drug_params['CLr']      # renal clearance (L/h)
    Kp   = drug_params.get('Kp', 1.0)  # tissue:plasma partition coefficient
    dose = drug_params['dose_mg']

    # Physiology
    Q_p  = physio['Q_portal']
    Q_l  = physio['Q_liver']
    Q_k  = physio['Q_kidney']
    Q_t  = physio['Q_tissue']
    V_gi = physio['V_gi']
    V_po = physio['V_portal']
    V_l  = physio['V_liver']
    V_a  = physio['V_arterial']
    V_v  = physio['V_venous']
    V_ti = physio['V_tissue']
    V_k  = physio['V_kidney']

    # ── GI absorption ─────────────────────────────────────────────────────────
    # Lag time handled by initial conditions
    dA_gi = -ka * A_gi

    # Amount absorbed into portal blood (mg/h)
    absorption_rate = ka * A_gi * F

    # ── Portal blood ──────────────────────────────────────────────────────────
    dC_portal = (absorption_rate - Q_p * C_portal) / V_po

    # ── Liver (first-pass + hepatic metabolism) ───────────────────────────────
    # Hepatic extraction
    liver_in  = Q_p * C_portal + (Q_l - Q_p) * C_arterial
    liver_out = Q_l * C_liver
    hepatic_clearance = CLh * C_liver
    dC_liver  = (liver_in - liver_out - hepatic_clearance) / V_l

    # ── Arterial blood ────────────────────────────────────────────────────────
    # Receives from: lung (venous)
    # Sends to: liver, kidney, tissue
    dC_arterial = (Q_t * C_venous/Kp - Q_t * C_arterial) / V_a

    # ── Venous blood ──────────────────────────────────────────────────────────
    # Collects from: tissue, kidney
    venous_in = Q_t * C_tissue/Kp + Q_k * C_kidney + Q_l * C_liver
    dC_venous = (venous_in - (Q_t + Q_k + Q_l) * C_venous) / V_v

    # ── Peripheral tissue ─────────────────────────────────────────────────────
    dC_tissue = (Q_t * C_arterial - Q_t * C_tissue/Kp) / V_ti

    # ── Kidney ────────────────────────────────────────────────────────────────
    renal_clearance = CLr * C_kidney
    dC_kidney = (Q_k * C_arterial - Q_k * C_kidney - renal_clearance) / V_k

    # ── Excretion ─────────────────────────────────────────────────────────────
    dA_excr = renal_clearance * V_k

    return [dA_gi, dC_portal, dC_liver, dC_arterial, dC_venous,
            dC_tissue, dC_kidney, dA_excr]

# ── SIMULATE DRUG ──────────────────────────────────────────────────────────────
def simulate_pbpk(drug_name, drug_params_earth, drug_params_space,
                  physio_earth, physio_space, t_end=12):
    """
    Run PBPK simulation for Earth and Space conditions
    Returns time array and concentration profiles
    """
    t_span = (0, t_end)
    t_eval = np.linspace(0, t_end, 500)

    # Initial conditions — dose in GI tract
    dose = drug_params_earth['dose_mg']
    y0 = [dose, 0, 0, 0, 0, 0, 0, 0]

    # ── Earth simulation ──────────────────────────────────────────────────────
    sol_e = solve_ivp(
        pbpk_odes, t_span, y0, t_eval=t_eval,
        args=(drug_params_earth, physio_earth),
        method='RK45', rtol=1e-6, atol=1e-8
    )

    # ── Space simulation ──────────────────────────────────────────────────────
    sol_s = solve_ivp(
        pbpk_odes, t_span, y0, t_eval=t_eval,
        args=(drug_params_space, physio_space),
        method='RK45', rtol=1e-6, atol=1e-8
    )

    return t_eval, sol_e, sol_s

# ── EXTRACT PK PARAMETERS ─────────────────────────────────────────────────────
def extract_pk_pbpk(t, C):
    """Extract Cmax, Tmax, AUC, t½ from PBPK output"""
    C = np.array(C)
    C = np.clip(C, 0, None)
    
    if C.max() == 0:
        return {'Cmax': 0, 'Tmax': 0, 'AUC': 0, 't12': 0}
    
    idx  = np.argmax(C)
    Cmax = C[idx]
    Tmax = t[idx]
    AUC  = trapezoid(C, t)

    # t½ from elimination phase
    elim = C[idx:]
    t_el = t[idx:]
    mask = elim > Cmax * 0.05
    if mask.sum() > 5:
        try:
            sl, _ = np.polyfit(t_el[mask], np.log(elim[mask] + 1e-10), 1)
            t12 = 0.693 / (-sl) if sl < 0 else np.nan
        except:
            t12 = np.nan
    else:
        t12 = np.nan

    return {'Cmax': round(Cmax,3), 'Tmax': round(Tmax,3),
            'AUC': round(AUC,3), 't12': round(t12,3) if not np.isnan(t12) else 'NA'}

def compute_dose_recommendation(drug_profile, pk_earth, pk_space):
    """Match Earth Cmax exposure with adjusted space dose."""
    earth_dose = drug_profile['dose_mg']
    cmax_e = pk_earth.get('Cmax', 0)
    cmax_s = pk_space.get('Cmax', 0)
    auc_e = pk_earth.get('AUC', 0)
    auc_s = pk_space.get('AUC', 0)

    if not cmax_e or not cmax_s or cmax_s <= 0:
        return None

    factor_cmax = cmax_e / cmax_s
    factor_auc = (auc_e / auc_s) if auc_s and auc_s > 0 else factor_cmax
    factor = (factor_cmax + factor_auc) / 2

    return {
        'earth_dose_mg': earth_dose,
        'space_dose_mg': round(earth_dose * factor_cmax, 1),
        'space_dose_auc_mg': round(earth_dose * factor_auc, 1),
        'adjustment_factor': round(factor_cmax, 3),
        'adjustment_pct': round((factor_cmax - 1) * 100, 1),
        'rationale': 'Cmax-matched dose to preserve peak exposure',
    }

# ── MAIN ──────────────────────────────────────────────────────────────────────
def run_pbpk_analysis(drug_profile, mission_days=30, body_weight=75, verbose=True):
    """
    Complete PBPK analysis for a drug
    drug_profile: output from drug_properties.analyze_drug()
    """
    drug_name = drug_profile['drug_name']
    
    if verbose:
        print(f"\n{'='*60}")
        print(f"  PBPK SIMULATION: {drug_name}")
        print(f"  Mission: {mission_days} days | Weight: {body_weight}kg")
        print(f"{'='*60}")

    # ── Earth drug parameters ─────────────────────────────────────────────────
    ke_e = drug_profile['ke_earth']
    ka_e = drug_profile['ka_earth'] or 1.5
    Vd_e = drug_profile['Vd_earth']
    F_e  = drug_profile['F_earth']
    
    # Clearance = ke × Vd
    CL_e = ke_e * Vd_e
    CLh_e = CL_e * 0.7   # 70% hepatic clearance (typical)
    CLr_e = CL_e * 0.3   # 30% renal clearance

    drug_params_earth = {
        'drug_name': drug_name,
        'dose_mg': drug_profile['dose_mg'],
        'F': F_e, 'Vd': Vd_e,
        'ka': ka_e, 'ke': ke_e,
        'CLh': CLh_e, 'CLr': CLr_e,
        'Kp': 1.0 + drug_profile['logP'] * 0.2,
    }

    # ── Space drug parameters ─────────────────────────────────────────────────
    ke_s = drug_profile['ke_space']
    ka_s = drug_profile['ka_space'] or ka_e * drug_profile['space_modifiers']['ka_factor']
    Vd_s = drug_profile['Vd_space']
    F_s  = drug_profile['F_space']
    
    CL_s  = ke_s * Vd_s
    CLh_s = CL_s * 0.7
    CLr_s = CL_s * 0.3

    drug_params_space = {
        'drug_name': drug_name,
        'dose_mg': drug_profile['dose_mg'],
        'F': F_s, 'Vd': Vd_s,
        'ka': ka_s, 'ke': ke_s,
        'CLh': CLh_s, 'CLr': CLr_s,
        'Kp': 1.0 + drug_profile['logP'] * 0.2,
    }

    # ── Physiology ────────────────────────────────────────────────────────────
    physio_earth = get_physiology(body_weight, mission_days=0)
    physio_space = get_physiology(body_weight, mission_days=mission_days)

    # ── Run simulation ────────────────────────────────────────────────────────
    t, sol_e, sol_s = simulate_pbpk(
        drug_name, drug_params_earth, drug_params_space,
        physio_earth, physio_space
    )

    # Venous blood concentration (systemic) = compartment [4]
    C_earth  = sol_e.y[4]   # venous = systemic
    C_space  = sol_s.y[4]
    C_liver_e = sol_e.y[2]
    C_liver_s = sol_s.y[2]
    C_tissue_e = sol_e.y[5]
    C_tissue_s = sol_s.y[5]

    # ── PK parameters ─────────────────────────────────────────────────────────
    pk_earth = extract_pk_pbpk(t, C_earth)
    pk_space = extract_pk_pbpk(t, C_space)

    if verbose:
        print(f"\n{'Parameter':<15} {'Earth':>12} {'Space':>12} {'Change':>10}")
        print("-"*52)
        for param in ['Cmax', 'Tmax', 'AUC', 't12']:
            val_e = pk_earth[param]
            val_s = pk_space[param]
            if isinstance(val_e, (int,float)) and isinstance(val_s, (int,float)) and val_e > 0:
                change = f"{((val_s-val_e)/val_e)*100:+.1f}%"
            else:
                change = "NA"
            unit = 'µg/mL' if param=='Cmax' else 'h' if param in ['Tmax','t12'] else 'µg·h/mL'
            print(f"{param+' ('+unit+')':<15} {str(val_e):>12} {str(val_s):>12} {change:>10}")

        # Dose recommendation
        dose_rec = compute_dose_recommendation(drug_profile, pk_earth, pk_space)
        if dose_rec:
            print(f"\n🎯 DOSE RECOMMENDATION:")
            print(f"   Earth dose  : {dose_rec['earth_dose_mg']:.0f} mg")
            print(f"   Space dose  : {dose_rec['space_dose_mg']:.0f} mg (×{dose_rec['adjustment_factor']:.2f})")
            print(f"   Adjustment  : {dose_rec['adjustment_pct']:+.1f}%")

    dose_rec = compute_dose_recommendation(drug_profile, pk_earth, pk_space)

    return {
        't': t,
        'C_earth': C_earth, 'C_space': C_space,
        'C_liver_e': C_liver_e, 'C_liver_s': C_liver_s,
        'C_tissue_e': C_tissue_e, 'C_tissue_s': C_tissue_s,
        'pk_earth': pk_earth, 'pk_space': pk_space,
        'dose_recommendation': dose_rec,
        'drug_params_earth': drug_params_earth,
        'drug_params_space': drug_params_space,
        'physio_earth': physio_earth,
        'physio_space': physio_space,
    }

if __name__ == '__main__':
    from pathlib import Path
    from drug_properties import analyze_drug

    ROOT = Path(__file__).resolve().parent.parent
    OUT = ROOT / "figures"
    OUT.mkdir(exist_ok=True)

    print("SPACE PK PIPELINE — LAYER 2: PBPK MODEL")
    print("CPT: Pharmacometrics & Systems Pharmacology\n")

    # Run for paracetamol
    drug = analyze_drug('Paracetamol', mission_days=30)
    result = run_pbpk_analysis(drug, mission_days=30)

    # Quick plot
    fig, axes = plt.subplots(1, 3, figsize=(15, 5))
    fig.patch.set_facecolor('#0a0a1a')
    titles = ['Venous Blood (Systemic)', 'Liver', 'Peripheral Tissue']
    data_pairs = [
        (result['C_earth'], result['C_space']),
        (result['C_liver_e'], result['C_liver_s']),
        (result['C_tissue_e'], result['C_tissue_s']),
    ]
    for ax, title, (Ce, Cs) in zip(axes, titles, data_pairs):
        ax.set_facecolor('#0d1117')
        ax.plot(result['t'], Ce, '#2196F3', lw=2.5, label='Earth')
        ax.plot(result['t'], Cs, '#F44336', lw=2.5, label='Space')
        ax.fill_between(result['t'], Ce, alpha=0.15, color='#2196F3')
        ax.fill_between(result['t'], Cs, alpha=0.15, color='#F44336')
        ax.set(xlabel='Time (h)', ylabel='Conc (mg/L)', title=title)
        ax.legend(facecolor='#1a1a2e', labelcolor='white', fontsize=9)
        ax.tick_params(colors='white')
        ax.xaxis.label.set_color('white')
        ax.yaxis.label.set_color('white')
        ax.title.set_color('white')
        ax.grid(True, alpha=0.2)

    plt.suptitle('PBPK: Paracetamol — Earth vs Space (30-day mission)',
                 color='white', fontsize=12, fontweight='bold')
    plt.tight_layout()
    plt.savefig(OUT / "pbpk_layer2.png",
                dpi=150, bbox_inches='tight', facecolor='#0a0a1a')
    print("\n✅ Layer 2 complete — PBPK plot saved")
