"""
================================================================================
LAYER 1: DRUG PROPERTIES PIPELINE
SpacePK Framework — CPT:PSP Paper
RDKit + BioPython molecular descriptor extraction
Auto-parameterization for PBPK model
================================================================================
"""

import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors, Crippen
from rdkit.Chem import QED

# ── Drug Database (SMILES + literature PK) ───────────────────────────────────
DRUG_DATABASE = {
    'Paracetamol': {
        'smiles': 'CC(=O)Nc1ccc(O)cc1',
        'dose_mg': 500,
        'lit_F': 0.85,
        'lit_Vd_Lkg': 0.90,
        'lit_t12_h': 2.5,
        'lit_ka_h': 1.5,
        'lit_ke_h': 0.277,
        'protein_binding': 0.15,
        'bcs_class': 1,
        'space_data': True,
        'cyp': 'CYP2E1/CYP3A4',
        'transporter': 'None'
    },
    'Ibuprofen': {
        'smiles': 'CC(C)Cc1ccc(C(C)C(=O)O)cc1',
        'dose_mg': 600,
        'lit_F': 0.80,
        'lit_Vd_Lkg': 0.15,
        'lit_t12_h': 2.0,
        'lit_ka_h': 2.0,
        'lit_ke_h': 0.347,
        'protein_binding': 0.99,
        'bcs_class': 2,
        'space_data': True,
        'cyp': 'CYP2C9',
        'transporter': 'OATP1B1'
    },
    'Ciprofloxacin': {
        'smiles': 'O=C(O)c1cn(C2CC2)c2cc(N3CCNCC3)c(F)cc2c1=O',
        'dose_mg': 250,
        'lit_F': 0.70,
        'lit_Vd_Lkg': 2.50,
        'lit_t12_h': 4.0,
        'lit_ka_h': 0.9,
        'lit_ke_h': 0.173,
        'protein_binding': 0.30,
        'bcs_class': 4,
        'space_data': True,
        'cyp': 'CYP1A2',
        'transporter': 'P-gp'
    },
    'Promethazine': {
        'smiles': 'CN(C)CCCN1c2ccccc2Sc2ccccc21',
        'dose_mg': 25,
        'lit_F': 0.25,
        'lit_Vd_Lkg': 13.0,
        'lit_t12_h': 12.0,
        'lit_ka_h': 0.5,
        'lit_ke_h': 0.058,
        'protein_binding': 0.76,
        'bcs_class': 1,
        'space_data': True,
        'cyp': 'CYP2D6',
        'transporter': 'None'
    },
    'Scopolamine': {
        'smiles': 'CN1C2CC(OC(=O)C(CO)c3ccccc3)CC1C2',
        'dose_mg': 0.4,
        'lit_F': 0.60,
        'lit_Vd_Lkg': 1.35,
        'lit_t12_h': 8.0,
        'lit_ka_h': 0.8,
        'lit_ke_h': 0.087,
        'protein_binding': 0.60,
        'bcs_class': 1,
        'space_data': True,
        'cyp': 'CYP3A4',
        'transporter': 'None'
    },
    'Lidocaine': {
        'smiles': 'CCN(CC)CC(=O)Nc1c(C)cccc1C',
        'dose_mg': 75,
        'lit_F': 0.35,
        'lit_Vd_Lkg': 1.30,
        'lit_t12_h': 1.8,
        'lit_ka_h': 'IV',
        'lit_ke_h': 0.385,
        'protein_binding': 0.70,
        'bcs_class': 1,
        'space_data': True,
        'cyp': 'CYP3A4/CYP1A2',
        'transporter': 'None'
    },
}

# ── BCS Classification thresholds ────────────────────────────────────────────
def classify_bcs(logP, mw, solubility_est):
    """
    BCS Classification:
    Class I: High solubility, High permeability
    Class II: Low solubility, High permeability
    Class III: High solubility, Low permeability
    Class IV: Low solubility, Low permeability
    """
    high_perm = logP > 1.72  # Caco-2 threshold
    high_sol  = solubility_est > 1.0  # mg/mL threshold (simplified)
    if high_perm and high_sol: return 'I'
    elif high_perm and not high_sol: return 'II'
    elif not high_perm and high_sol: return 'III'
    else: return 'IV'

# ── P-gp substrate prediction (simplified Lipinski-based) ───────────────────
def predict_pgp(mol, mw, logP):
    hbd = rdMolDescriptors.CalcNumHBD(mol)
    hba = rdMolDescriptors.CalcNumHBA(mol)
    # Seelig rules for P-gp substrate
    if mw > 400 and logP > 2 and hba > 4:
        return True
    return False

# ── Estimate bioavailability from molecular properties ───────────────────────
def estimate_bioavailability(logP, mw, hbd, hba, psa):
    """
    Lipinski rule-based F estimation
    Not precise but sufficient for parameterization
    """
    score = 1.0
    if mw > 500: score *= 0.7
    if logP > 5 or logP < 0: score *= 0.8
    if hbd > 5: score *= 0.75
    if hba > 10: score *= 0.75
    if psa > 140: score *= 0.6
    return round(min(score, 0.95), 2)

# ── Estimate Vd from logP ─────────────────────────────────────────────────────
def estimate_vd(logP, protein_binding=0.5):
    """
    Oie-Tozer equation simplified
    Vd ~ function of logP and protein binding
    """
    # Empirical: higher logP → higher Vd (tissue distribution)
    vd_base = 0.5 + (logP * 0.8)
    vd_adj  = vd_base * (1 - protein_binding * 0.3)
    return round(max(0.1, min(vd_adj, 20.0)), 2)

# ── Space modification factors ────────────────────────────────────────────────
def get_space_modifiers(mission_days=30, drug_bcs='I'):
    """
    Space physiology modifiers based on literature:
    Gandia 2003, Kovachevich 2009, Dello Russo 2022,
    Polyakov 2021, Leach 1981
    """
    if mission_days == 0:
        return {
            'phase': 'earth',
            'mission_days': 0,
            'ka_factor': 1.0,
            'F_factor': 1.0,
            'Vd_factor': 1.0,
            'ke_factor': 1.0,
        }

    # Phase-specific modifications
    if mission_days <= 3:
        phase = 'acute'
        # Fluid shift phase — faster gastric emptying (Gandia 2003)
        ka_factor  = 1.50   # absorption faster
        F_factor   = 1.10
        Vd_factor  = 0.99   # minimal change
        ke_factor  = 1.00
    elif mission_days <= 14:
        phase = 'adaptation'
        # Mixed phase — variable (Cintron 1987)
        ka_factor  = 1.00
        F_factor   = 1.15
        Vd_factor  = 0.98
        ke_factor  = 0.98
    else:
        phase = 'chronic'
        # Long-term — delayed absorption (Kovachevich 2009, Polyakov 2021)
        ka_factor  = 0.64   # Kovachevich: absorption rate 124.45%±24.27 → 36% reduction
        F_factor   = 1.267  # Kovachevich: bioavail 126.72%±24.04
        Vd_factor  = 0.97   # Leach 1981: body water -3%
        ke_factor  = 0.96   # Conservative estimate
    
    # BCS class adjustment
    if drug_bcs == 'II':
        # Low solubility drugs — dissolution limited
        # Slower gastric transit actually helps dissolution
        ka_factor *= 1.2
    elif drug_bcs == 'IV':
        # Worst case — both solubility and permeability affected
        ka_factor *= 0.8
        F_factor  *= 0.9

    return {
        'phase': phase,
        'mission_days': mission_days,
        'ka_factor': round(ka_factor, 3),
        'F_factor': round(F_factor, 3),
        'Vd_factor': round(Vd_factor, 3),
        'ke_factor': round(ke_factor, 3)
    }

def _build_profile(drug_name, smiles, dose_mg, body_weight_kg, mission_days,
                   mw, logP, hbd, hba, psa, rotb, arom, qed, lipinski, bcs, pgp,
                   pb, F, Vd_per_kg, Vd, t12, ke, ka, space,
                   F_space, Vd_space, ke_space, ka_space):
    return {
        'drug_name': drug_name,
        'smiles': smiles,
        'dose_mg': dose_mg,
        'body_weight_kg': body_weight_kg,
        'mission_days': mission_days,
        'MW': mw, 'logP': logP, 'HBD': hbd, 'HBA': hba,
        'PSA': psa, 'QED': qed, 'BCS': bcs, 'PgP': pgp,
        'lipinski_ok': lipinski,
        'F_earth': F, 'Vd_earth': Vd, 'ke_earth': ke,
        'ka_earth': ka, 't12_earth': t12, 'protein_binding': pb,
        'F_space': F_space, 'Vd_space': Vd_space,
        'ke_space': ke_space, 'ka_space': ka_space,
        'space_phase': space['phase'],
        'space_modifiers': space,
    }

# Drugs with literature space-PK data or ISS-relevant profiles
DRUG_LIST = list(DRUG_DATABASE.keys())

# ── MAIN FUNCTION: Analyze any drug ──────────────────────────────────────────
def analyze_drug(drug_name, smiles=None, dose_mg=None,
                 mission_days=30, body_weight_kg=75, verbose=True):
    """
    Complete drug property analysis for Space PK simulation
    
    Parameters:
        drug_name: str — drug name (checks database first)
        smiles: str — SMILES string (if not in database)
        dose_mg: float — oral dose in mg
        mission_days: int — spaceflight duration
        body_weight_kg: float — astronaut body weight
    
    Returns:
        dict — complete drug profile + space modifiers
    """
    if verbose:
        print(f"\n{'='*60}")
        print(f"  DRUG ANALYSIS: {drug_name}")
        print(f"{'='*60}")

    # Get from database or use provided SMILES
    if drug_name in DRUG_DATABASE:
        db = DRUG_DATABASE[drug_name]
        smiles  = db['smiles']
        dose_mg = dose_mg or db['dose_mg']
        lit_data = db
    else:
        if not smiles:
            raise ValueError(f"Drug '{drug_name}' not in database. Provide SMILES.")
        lit_data = {}
        dose_mg = dose_mg or 500

    # Parse molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")

    # ── RDKit descriptors ────────────────────────────────────────────────────
    mw     = Descriptors.MolWt(mol)
    logP   = Crippen.MolLogP(mol)
    hbd    = rdMolDescriptors.CalcNumHBD(mol)
    hba    = rdMolDescriptors.CalcNumHBA(mol)
    psa    = rdMolDescriptors.CalcTPSA(mol)
    rotb   = rdMolDescriptors.CalcNumRotatableBonds(mol)
    rings  = rdMolDescriptors.CalcNumRings(mol)
    arom   = rdMolDescriptors.CalcNumAromaticRings(mol)
    qed    = QED.qed(mol)

    # Lipinski compliance
    lipinski = all([mw<=500, logP<=5, hbd<=5, hba<=10])

    # Solubility estimate (ESOL model simplified)
    sol_est = 10 ** (0.5 - 0.01*(mw-350)/50 - 0.5*logP)

    # BCS Classification
    bcs = classify_bcs(logP, mw, sol_est)

    # P-gp prediction
    pgp = predict_pgp(mol, mw, logP)

    # PK parameter estimates
    pb  = lit_data.get('protein_binding', 0.5)
    F   = lit_data.get('lit_F', estimate_bioavailability(logP,mw,hbd,hba,psa))
    Vd_per_kg = lit_data.get('lit_Vd_Lkg', estimate_vd(logP, pb))
    Vd  = Vd_per_kg * body_weight_kg
    t12 = lit_data.get('lit_t12_h', 2.0)
    ke  = 0.693 / t12
    ka  = lit_data.get('lit_ka_h', 1.5) if lit_data.get('lit_ka_h','IV') != 'IV' else None

    # Space modifiers
    space = get_space_modifiers(mission_days, bcs)

    # Space-adjusted parameters
    F_space  = min(F  * space['F_factor'], 0.99)
    Vd_space = Vd * space['Vd_factor']
    ke_space = ke * space['ke_factor']
    ka_space = (ka * space['ka_factor']) if ka else None

    # ── Print results ────────────────────────────────────────────────────────
    if not verbose:
        return _build_profile(
            drug_name, smiles, dose_mg, body_weight_kg, mission_days,
            mw, logP, hbd, hba, psa, rotb, arom, qed, lipinski, bcs, pgp,
            pb, F, Vd_per_kg, Vd, t12, ke, ka, space,
            F_space, Vd_space, ke_space, ka_space,
        )

    print(f"\n📊 MOLECULAR PROPERTIES (RDKit)")
    print(f"  MW              : {mw:.2f} g/mol")
    print(f"  LogP            : {logP:.2f}")
    print(f"  H-bond donors   : {hbd}")
    print(f"  H-bond acceptors: {hba}")
    print(f"  PSA             : {psa:.1f} Å²")
    print(f"  Rotatable bonds : {rotb}")
    print(f"  Aromatic rings  : {arom}")
    print(f"  QED score       : {qed:.3f}")
    print(f"  Lipinski OK     : {'✅' if lipinski else '❌'}")
    print(f"  BCS Class       : {bcs}")
    print(f"  P-gp substrate  : {'Yes ⚠️' if pgp else 'No ✅'}")

    print(f"\n💊 PK PARAMETERS (Literature + Estimated)")
    print(f"  Dose            : {dose_mg} mg")
    print(f"  Bioavailability : {F:.3f} ({F*100:.0f}%)")
    print(f"  Vd (Earth)      : {Vd:.1f} L ({Vd_per_kg:.2f} L/kg)")
    print(f"  t½ (Earth)      : {t12:.2f} h")
    print(f"  ke (Earth)      : {ke:.4f} /h")
    if ka: print(f"  ka (Earth)      : {ka:.4f} /h")
    print(f"  Protein binding : {pb*100:.0f}%")

    print(f"\n🚀 SPACE MODIFIERS (Mission day {mission_days} — {space['phase']} phase)")
    print(f"  ka factor       : ×{space['ka_factor']} (absorption rate)")
    print(f"  F factor        : ×{space['F_factor']} (bioavailability)")
    print(f"  Vd factor       : ×{space['Vd_factor']} (distribution)")
    print(f"  ke factor       : ×{space['ke_factor']} (elimination)")

    print(f"\n🛸 SPACE-ADJUSTED PK PARAMETERS")
    print(f"  F (Space)       : {F_space:.3f} ({F_space*100:.0f}%)")
    print(f"  Vd (Space)      : {Vd_space:.1f} L")
    print(f"  ke (Space)      : {ke_space:.4f} /h")
    if ka_space: print(f"  ka (Space)      : {ka_space:.4f} /h")

    return _build_profile(
        drug_name, smiles, dose_mg, body_weight_kg, mission_days,
        mw, logP, hbd, hba, psa, rotb, arom, qed, lipinski, bcs, pgp,
        pb, F, Vd_per_kg, Vd, t12, ke, ka, space,
        F_space, Vd_space, ke_space, ka_space,
    )

# ── Run all drugs ─────────────────────────────────────────────────────────────
if __name__ == '__main__':
    print("SPACE PK PIPELINE — LAYER 1: DRUG PROPERTIES")
    print("CPT: Pharmacometrics & Systems Pharmacology")
    
    results = {}
    for drug in ['Paracetamol', 'Ibuprofen', 'Ciprofloxacin', 'Promethazine']:
        results[drug] = analyze_drug(drug, mission_days=30)
    
    print(f"\n{'='*60}")
    print("  MULTI-DRUG SPACE PARAMETER SUMMARY")
    print(f"{'='*60}")
    print(f"{'Drug':<15} {'F Earth':>8} {'F Space':>8} {'Vd Earth':>10} {'Vd Space':>10} {'t½ Earth':>9} {'BCS':>4}")
    print("-"*70)
    for drug, r in results.items():
        t12_e = round(0.693/r['ke_earth'], 2)
        print(f"{drug:<15} {r['F_earth']:>8.3f} {r['F_space']:>8.3f} "
              f"{r['Vd_earth']:>10.1f} {r['Vd_space']:>10.1f} "
              f"{t12_e:>9.2f} {r['BCS']:>4}")
    
    print("\n✅ Layer 1 complete")
