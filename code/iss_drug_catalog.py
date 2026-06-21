"""
ISS / NASA spaceflight medical kit + literature-validated PK catalog.
Merged into DRUG_DATABASE at import time.
"""

# category, space_data = studied in HDT/spaceflight literature
ISS_DRUG_CATALOG = {
    # ── Analgesics / Antipyretics ─────────────────────────────────────────────
    'Paracetamol': {
        'aliases': ['Acetaminophen'],
        'category': 'Analgesic',
        'iss_kit': True,
        'smiles': 'CC(=O)Nc1ccc(O)cc1',
        'dose_mg': 500,
        'lit_F': 0.85, 'lit_Vd_Lkg': 0.90, 'lit_t12_h': 2.5,
        'lit_ka_h': 1.5, 'lit_ke_h': 0.277,
        'protein_binding': 0.15, 'bcs_class': 1,
        'space_data': True, 'cyp': 'CYP2E1/CYP3A4', 'transporter': 'None',
    },
    'Aspirin': {
        'category': 'Analgesic / Antiplatelet',
        'iss_kit': True,
        'smiles': 'CC(=O)Oc1ccccc1C(=O)O',
        'dose_mg': 325,
        'lit_F': 0.68, 'lit_Vd_Lkg': 0.15, 'lit_t12_h': 0.25,
        'lit_ka_h': 3.5, 'lit_ke_h': 2.77,
        'protein_binding': 0.90, 'bcs_class': 2,
        'space_data': False, 'cyp': 'CYP2C9', 'transporter': 'None',
    },
    'Ibuprofen': {
        'category': 'NSAID',
        'iss_kit': True,
        'smiles': 'CC(C)Cc1ccc(C(C)C(=O)O)cc1',
        'dose_mg': 600,
        'lit_F': 0.80, 'lit_Vd_Lkg': 0.15, 'lit_t12_h': 2.0,
        'lit_ka_h': 2.0, 'lit_ke_h': 0.347,
        'protein_binding': 0.99, 'bcs_class': 2,
        'space_data': True, 'cyp': 'CYP2C9', 'transporter': 'OATP1B1',
    },
    'Tramadol': {
        'category': 'Analgesic (Opioid)',
        'iss_kit': True,
        'smiles': 'COC1=CC=CC(=C1)C2(CCC(C2)(C(=O)O)O)N(C)C',
        'dose_mg': 50,
        'lit_F': 0.75, 'lit_Vd_Lkg': 2.6, 'lit_t12_h': 6.0,
        'lit_ka_h': 1.2, 'lit_ke_h': 0.116,
        'protein_binding': 0.20, 'bcs_class': 1,
        'space_data': False, 'cyp': 'CYP2D6/CYP3A4', 'transporter': 'None',
    },
    # ── Antibiotics ────────────────────────────────────────────────────────────
    'Ciprofloxacin': {
        'category': 'Antibiotic (Fluoroquinolone)',
        'iss_kit': True,
        'smiles': 'O=C(O)c1cn(C2CC2)c2cc(N3CCNCC3)c(F)cc2c1=O',
        'dose_mg': 250,
        'lit_F': 0.70, 'lit_Vd_Lkg': 2.50, 'lit_t12_h': 4.0,
        'lit_ka_h': 0.9, 'lit_ke_h': 0.173,
        'protein_binding': 0.30, 'bcs_class': 4,
        'space_data': True, 'cyp': 'CYP1A2', 'transporter': 'P-gp',
    },
    'Amoxicillin': {
        'category': 'Antibiotic (Penicillin)',
        'iss_kit': True,
        'smiles': 'CC1(C(N2C(C(=O)NC2=O)SC1)C(=O)O)NC(=O)C3=CC=CC=C3',
        'dose_mg': 500,
        'lit_F': 0.90, 'lit_Vd_Lkg': 0.25, 'lit_t12_h': 1.2,
        'lit_ka_h': 2.5, 'lit_ke_h': 0.578,
        'protein_binding': 0.18, 'bcs_class': 3,
        'space_data': False, 'cyp': 'None', 'transporter': 'None',
    },
    'Azithromycin': {
        'category': 'Antibiotic (Macrolide)',
        'iss_kit': True,
        'smiles': 'CC[C@H]1OC(=O)[C@H](C)[C@@H](O[C@H]2C[C@@](C)(OC)[C@@H](O)[C@H](C)O2)[C@H](C)[C@@H](O[C@@H]2O[C@H](C)C[C@H](N(C)C)[C@H]2O)C(=O)[C@H](C)C(=O)[C@H](C)[C@@H](O)[C@]1(C)O',
        'dose_mg': 250,
        'lit_F': 0.38, 'lit_Vd_Lkg': 31.0, 'lit_t12_h': 68.0,
        'lit_ka_h': 0.5, 'lit_ke_h': 0.010,
        'protein_binding': 0.51, 'bcs_class': 4,
        'space_data': False, 'cyp': 'CYP3A4', 'transporter': 'P-gp',
    },
    'Doxycycline': {
        'category': 'Antibiotic (Tetracycline)',
        'iss_kit': True,
        'smiles': 'C[C@@H]1[C@H]2[C@H]3C(=C(O)[C@]1(O)C(=O)C2=O)C(=O)C4=C(C3=O)C(=O)C=C(C4=O)N(C)C',
        'dose_mg': 100,
        'lit_F': 0.95, 'lit_Vd_Lkg': 0.95, 'lit_t12_h': 18.0,
        'lit_ka_h': 0.8, 'lit_ke_h': 0.039,
        'protein_binding': 0.90, 'bcs_class': 4,
        'space_data': False, 'cyp': 'CYP3A4', 'transporter': 'None',
    },
    # ── Anti-nausea / Motion sickness ──────────────────────────────────────────
    'Promethazine': {
        'category': 'Antiemetic / Antihistamine',
        'iss_kit': True,
        'smiles': 'CN(C)CCCN1c2ccccc2Sc2ccccc21',
        'dose_mg': 25,
        'lit_F': 0.25, 'lit_Vd_Lkg': 13.0, 'lit_t12_h': 12.0,
        'lit_ka_h': 0.5, 'lit_ke_h': 0.058,
        'protein_binding': 0.76, 'bcs_class': 1,
        'space_data': True, 'cyp': 'CYP2D6', 'transporter': 'None',
    },
    'Scopolamine': {
        'category': 'Antiemetic / Motion sickness',
        'iss_kit': True,
        'smiles': 'CN1C2CC(OC(=O)C(CO)c3ccccc3)CC1C2',
        'dose_mg': 0.4,
        'lit_F': 0.60, 'lit_Vd_Lkg': 1.35, 'lit_t12_h': 8.0,
        'lit_ka_h': 0.8, 'lit_ke_h': 0.087,
        'protein_binding': 0.60, 'bcs_class': 1,
        'space_data': True, 'cyp': 'CYP3A4', 'transporter': 'None',
    },
    'Meclizine': {
        'category': 'Antiemetic / Motion sickness',
        'iss_kit': True,
        'smiles': 'CN1CCN(CC1)C(C2=CC=CC=C2)C3=CC=CC=C3Cl',
        'dose_mg': 25,
        'lit_F': 0.50, 'lit_Vd_Lkg': 6.0, 'lit_t12_h': 6.0,
        'lit_ka_h': 0.6, 'lit_ke_h': 0.116,
        'protein_binding': 0.90, 'bcs_class': 2,
        'space_data': False, 'cyp': 'CYP2D6', 'transporter': 'None',
    },
    'Ondansetron': {
        'category': 'Antiemetic (5-HT3)',
        'iss_kit': True,
        'smiles': 'CN1C2=C(C(=O)N(C1=O)C3=CC=CC=C3)N(C=N2)C4=CC=C(C=C4)C(F)(F)F',
        'dose_mg': 8,
        'lit_F': 0.60, 'lit_Vd_Lkg': 2.5, 'lit_t12_h': 3.5,
        'lit_ka_h': 1.8, 'lit_ke_h': 0.198,
        'protein_binding': 0.70, 'bcs_class': 3,
        'space_data': False, 'cyp': 'CYP3A4/CYP1A2', 'transporter': 'P-gp',
    },
    'Metoclopramide': {
        'category': 'Antiemetic / Prokinetic',
        'iss_kit': True,
        'smiles': 'CCN(CC)CCNC(=O)C1=CC(=C(C=C1OC)Cl)N',
        'dose_mg': 10,
        'lit_F': 0.80, 'lit_Vd_Lkg': 3.5, 'lit_t12_h': 5.0,
        'lit_ka_h': 1.5, 'lit_ke_h': 0.139,
        'protein_binding': 0.30, 'bcs_class': 3,
        'space_data': False, 'cyp': 'CYP2D6', 'transporter': 'None',
    },
    'Diphenhydramine': {
        'category': 'Antihistamine / Sedative',
        'iss_kit': True,
        'smiles': 'CN(C)CCOC(c1ccccc1)c2ccccc2',
        'dose_mg': 25,
        'lit_F': 0.70, 'lit_Vd_Lkg': 3.3, 'lit_t12_h': 4.0,
        'lit_ka_h': 1.2, 'lit_ke_h': 0.173,
        'protein_binding': 0.98, 'bcs_class': 1,
        'space_data': False, 'cyp': 'CYP2D6', 'transporter': 'None',
    },
    # ── Cardiovascular ─────────────────────────────────────────────────────────
    'Nifedipine': {
        'category': 'Cardiovascular (CCB)',
        'iss_kit': False,
        'smiles': 'COC(=O)C1=C(C)NC(C)=C(C1c2ccccc2[N+](=O)[O-])C(=O)OC',
        'dose_mg': 20,
        'lit_F': 0.65, 'lit_Vd_Lkg': 1.5, 'lit_t12_h': 2.0,
        'lit_ka_h': 2.5, 'lit_ke_h': 0.347,
        'protein_binding': 0.95, 'bcs_class': 2,
        'space_data': True, 'cyp': 'CYP3A4', 'transporter': 'None',
    },
    'Verapamil': {
        'category': 'Cardiovascular (CCB)',
        'iss_kit': False,
        'smiles': 'COc1ccc(CCN(C)CCc2ccc(OC)c(OC)c2)cc1OC',
        'dose_mg': 80,
        'lit_F': 0.22, 'lit_Vd_Lkg': 3.5, 'lit_t12_h': 6.0,
        'lit_ka_h': 1.0, 'lit_ke_h': 0.116,
        'protein_binding': 0.90, 'bcs_class': 2,
        'space_data': True, 'cyp': 'CYP3A4', 'transporter': 'P-gp',
    },
    'Propranolol': {
        'category': 'Cardiovascular (Beta-blocker)',
        'iss_kit': False,
        'smiles': 'CC(C)NCC(O)c1ccc(O)cc1',
        'dose_mg': 80,
        'lit_F': 0.25, 'lit_Vd_Lkg': 4.0, 'lit_t12_h': 4.0,
        'lit_ka_h': 1.5, 'lit_ke_h': 0.173,
        'protein_binding': 0.90, 'bcs_class': 1,
        'space_data': True, 'cyp': 'CYP2D6/CYP1A2', 'transporter': 'None',
    },
    'Metoprolol': {
        'category': 'Cardiovascular (Beta-blocker)',
        'iss_kit': True,
        'smiles': 'COCCc1ccc(OCC(O)CNC(C)C)cc1',
        'dose_mg': 50,
        'lit_F': 0.50, 'lit_Vd_Lkg': 4.0, 'lit_t12_h': 3.5,
        'lit_ka_h': 1.8, 'lit_ke_h': 0.198,
        'protein_binding': 0.12, 'bcs_class': 1,
        'space_data': False, 'cyp': 'CYP2D6', 'transporter': 'None',
    },
    'Furosemide': {
        'category': 'Diuretic',
        'iss_kit': True,
        'smiles': 'NS(=O)(=O)c1cc(C(=O)Nc2ccccc2Cl)c2c(c1)NC(=O)CC2',
        'dose_mg': 40,
        'lit_F': 0.65, 'lit_Vd_Lkg': 0.12, 'lit_t12_h': 2.0,
        'lit_ka_h': 2.0, 'lit_ke_h': 0.347,
        'protein_binding': 0.98, 'bcs_class': 4,
        'space_data': True, 'cyp': 'None', 'transporter': 'OAT1/OAT3',
    },
    'Enalapril': {
        'category': 'Cardiovascular (ACE inhibitor)',
        'iss_kit': True,
        'smiles': 'CCOC(=O)C(CCC1=CC=CC=C1)N[C@@H](C)C(=O)N2CCC[C@H]2C(=O)O',
        'dose_mg': 10,
        'lit_F': 0.60, 'lit_Vd_Lkg': 1.8, 'lit_t12_h': 11.0,
        'lit_ka_h': 0.8, 'lit_ke_h': 0.063,
        'protein_binding': 0.60, 'bcs_class': 3,
        'space_data': False, 'cyp': 'None', 'transporter': 'None',
    },
    # ── Sedatives / Anxiolytics / Stimulants ───────────────────────────────────
    'Lorazepam': {
        'category': 'Anxiolytic / Sedative',
        'iss_kit': True,
        'smiles': 'OC1N=C(c2ccccc2Cl)c2cc(Cl)ccc2N1C1OC(CO)C(O)C1O',
        'dose_mg': 1,
        'lit_F': 0.90, 'lit_Vd_Lkg': 1.3, 'lit_t12_h': 12.0,
        'lit_ka_h': 1.0, 'lit_ke_h': 0.058,
        'protein_binding': 0.85, 'bcs_class': 2,
        'space_data': False, 'cyp': 'CYP3A4/UGT', 'transporter': 'None',
    },
    'Diazepam': {
        'category': 'Anxiolytic / Sedative',
        'iss_kit': True,
        'smiles': 'CN1C(=O)CN=C(C2=C1C=CC(=C2)Cl)C3=CC=CC=C3',
        'dose_mg': 5,
        'lit_F': 1.00, 'lit_Vd_Lkg': 1.1, 'lit_t12_h': 48.0,
        'lit_ka_h': 0.8, 'lit_ke_h': 0.014,
        'protein_binding': 0.98, 'bcs_class': 2,
        'space_data': False, 'cyp': 'CYP3A4/CYP2C19', 'transporter': 'None',
    },
    'Midazolam': {
        'category': 'Sedative (Benzodiazepine)',
        'iss_kit': True,
        'smiles': 'CC1=NC=C2C(=N1)C(=O)N(C(=O)N2C)C3=CC=CC=C3Cl',
        'dose_mg': 5,
        'lit_F': 0.40, 'lit_Vd_Lkg': 1.5, 'lit_t12_h': 2.5,
        'lit_ka_h': 1.5, 'lit_ke_h': 0.277,
        'protein_binding': 0.97, 'bcs_class': 2,
        'space_data': False, 'cyp': 'CYP3A4', 'transporter': 'P-gp',
    },
    'Dextroamphetamine': {
        'category': 'Stimulant (CNS)',
        'iss_kit': True,
        'smiles': 'CC(N)Cc1ccccc1',
        'dose_mg': 5,
        'lit_F': 0.75, 'lit_Vd_Lkg': 4.0, 'lit_t12_h': 10.0,
        'lit_ka_h': 1.2, 'lit_ke_h': 0.069,
        'protein_binding': 0.20, 'bcs_class': 1,
        'space_data': True, 'cyp': 'CYP2D6', 'transporter': 'None',
    },
    'Modafinil': {
        'category': 'Stimulant / Wakefulness',
        'iss_kit': False,
        'smiles': 'CC(C)S(=O)(=O)C1=CC=C(C=C1)C(=O)C2=CC=CC=N2',
        'dose_mg': 200,
        'lit_F': 0.62, 'lit_Vd_Lkg': 0.9, 'lit_t12_h': 15.0,
        'lit_ka_h': 0.5, 'lit_ke_h': 0.046,
        'protein_binding': 0.60, 'bcs_class': 2,
        'space_data': False, 'cyp': 'CYP3A4', 'transporter': 'None',
    },
    # ── GI ─────────────────────────────────────────────────────────────────────
    'Famotidine': {
        'category': 'GI (H2 blocker)',
        'iss_kit': True,
        'smiles': 'NC(=N)NC1=NC(CSC2=NC=CN2)=CS1',
        'dose_mg': 20,
        'lit_F': 0.50, 'lit_Vd_Lkg': 1.3, 'lit_t12_h': 3.0,
        'lit_ka_h': 1.5, 'lit_ke_h': 0.231,
        'protein_binding': 0.15, 'bcs_class': 3,
        'space_data': False, 'cyp': 'CYP1A2', 'transporter': 'None',
    },
    'Omeprazole': {
        'category': 'GI (PPI)',
        'iss_kit': True,
        'smiles': 'COc1ccc2[nH]c(S(=O)Cc3ncc(C)c(OC)c3C)nc2c1',
        'dose_mg': 20,
        'lit_F': 0.40, 'lit_Vd_Lkg': 0.3, 'lit_t12_h': 1.0,
        'lit_ka_h': 1.0, 'lit_ke_h': 0.693,
        'protein_binding': 0.95, 'bcs_class': 2,
        'space_data': False, 'cyp': 'CYP2C19/CYP3A4', 'transporter': 'None',
    },
    'Loperamide': {
        'category': 'GI (Antidiarrheal)',
        'iss_kit': True,
        'smiles': 'CC(C)(C)C(=O)N(CCC(O)(C1CCCCC1)C2CCCCC2)CC3=CC=C(C=C3)Cl',
        'dose_mg': 2,
        'lit_F': 0.40, 'lit_Vd_Lkg': 18.0, 'lit_t12_h': 11.0,
        'lit_ka_h': 0.8, 'lit_ke_h': 0.063,
        'protein_binding': 0.97, 'bcs_class': 2,
        'space_data': False, 'cyp': 'CYP3A4/CYP2C8', 'transporter': 'P-gp',
    },
    'Bisacodyl': {
        'category': 'GI (Laxative)',
        'iss_kit': True,
        'smiles': 'CC(=O)Oc1cc(C(=O)Oc2ccccc2)ccc1C(=O)O',
        'dose_mg': 5,
        'lit_F': 0.05, 'lit_Vd_Lkg': 0.5, 'lit_t12_h': 8.0,
        'lit_ka_h': 0.3, 'lit_ke_h': 0.087,
        'protein_binding': 0.95, 'bcs_class': 4,
        'space_data': False, 'cyp': 'None', 'transporter': 'None',
    },
    # ── Respiratory ────────────────────────────────────────────────────────────
    'Albuterol': {
        'category': 'Respiratory (Bronchodilator)',
        'iss_kit': True,
        'smiles': 'CC(C)(C)NCC(O)c1ccc(O)c(CO)c1',
        'dose_mg': 2.5,
        'lit_F': 0.10, 'lit_Vd_Lkg': 2.5, 'lit_t12_h': 4.0,
        'lit_ka_h': 2.0, 'lit_ke_h': 0.173,
        'protein_binding': 0.10, 'bcs_class': 3,
        'space_data': False, 'cyp': 'None', 'transporter': 'None',
    },
    'Pseudoephedrine': {
        'category': 'Respiratory (Decongestant)',
        'iss_kit': True,
        'smiles': 'CNCC(O)c1ccccc1',
        'dose_mg': 60,
        'lit_F': 1.00, 'lit_Vd_Lkg': 2.6, 'lit_t12_h': 6.0,
        'lit_ka_h': 2.0, 'lit_ke_h': 0.116,
        'protein_binding': 0.70, 'bcs_class': 1,
        'space_data': False, 'cyp': 'None', 'transporter': 'None',
    },
    'Fluticasone': {
        'category': 'Respiratory (Corticosteroid)',
        'iss_kit': True,
        'smiles': 'CC12C(C(C(O1)C(C)(F)F)C(C2=O)Cl)C(=O)SCF',
        'dose_mg': 0.1,
        'lit_F': 0.01, 'lit_Vd_Lkg': 4.0, 'lit_t12_h': 7.0,
        'lit_ka_h': 0.5, 'lit_ke_h': 0.099,
        'protein_binding': 0.99, 'bcs_class': 4,
        'space_data': False, 'cyp': 'CYP3A4', 'transporter': 'P-gp',
    },
    # ── Emergency / IV / Local ─────────────────────────────────────────────────
    'Lidocaine': {
        'category': 'Local anesthetic / Antiarrhythmic',
        'iss_kit': True,
        'smiles': 'CCN(CC)CC(=O)Nc1c(C)cccc1C',
        'dose_mg': 75,
        'lit_F': 0.35, 'lit_Vd_Lkg': 1.30, 'lit_t12_h': 1.8,
        'lit_ka_h': 'IV', 'lit_ke_h': 0.385,
        'protein_binding': 0.70, 'bcs_class': 1,
        'space_data': True, 'cyp': 'CYP3A4/CYP1A2', 'transporter': 'None',
    },
    'Epinephrine': {
        'category': 'Emergency (Vasopressor)',
        'iss_kit': True,
        'smiles': 'CNCC(O)c1ccc(O)c(O)c1',
        'dose_mg': 0.3,
        'lit_F': 0.05, 'lit_Vd_Lkg': 0.5, 'lit_t12_h': 0.05,
        'lit_ka_h': 'IV', 'lit_ke_h': 13.86,
        'protein_binding': 0.00, 'bcs_class': 3,
        'space_data': False, 'cyp': 'COMT/MAO', 'transporter': 'None',
    },
    'Atropine': {
        'category': 'Emergency (Anticholinergic)',
        'iss_kit': True,
        'smiles': 'CN1[C@H]2CC[C@@H]1CC(C2)OC(=O)C(CO)c3ccccc3',
        'dose_mg': 0.5,
        'lit_F': 0.50, 'lit_Vd_Lkg': 2.5, 'lit_t12_h': 3.0,
        'lit_ka_h': 'IV', 'lit_ke_h': 0.231,
        'protein_binding': 0.18, 'bcs_class': 1,
        'space_data': False, 'cyp': 'None', 'transporter': 'None',
    },
    'Dexamethasone': {
        'category': 'Corticosteroid',
        'iss_kit': True,
        'smiles': 'CC1CC2C3CCC4=CC(=O)C=CC4(C)C3(F)C(O)CC2(C)C1(O)C(=O)CO',
        'dose_mg': 4,
        'lit_F': 0.80, 'lit_Vd_Lkg': 0.8, 'lit_t12_h': 4.0,
        'lit_ka_h': 1.5, 'lit_ke_h': 0.173,
        'protein_binding': 0.77, 'bcs_class': 2,
        'space_data': False, 'cyp': 'CYP3A4', 'transporter': 'P-gp',
    },
    'Prednisone': {
        'category': 'Corticosteroid',
        'iss_kit': True,
        'smiles': 'CC1=CC(=O)C2C(C1=O)C=CC3C2CCC4(C3CCC4(O)C(=O)CO)C',
        'dose_mg': 20,
        'lit_F': 0.80, 'lit_Vd_Lkg': 0.9, 'lit_t12_h': 3.0,
        'lit_ka_h': 1.2, 'lit_ke_h': 0.231,
        'protein_binding': 0.90, 'bcs_class': 2,
        'space_data': False, 'cyp': 'CYP3A4', 'transporter': 'None',
    },
    'Acetazolamide': {
        'category': 'Diuretic / Altitude sickness',
        'iss_kit': True,
        'smiles': 'CC(=O)Nc1nnc(S(N)(=O)=O)s1',
        'dose_mg': 250,
        'lit_F': 0.95, 'lit_Vd_Lkg': 0.15, 'lit_t12_h': 6.0,
        'lit_ka_h': 1.5, 'lit_ke_h': 0.116,
        'protein_binding': 0.90, 'bcs_class': 3,
        'space_data': False, 'cyp': 'None', 'transporter': 'None',
    },
    'Melatonin': {
        'category': 'Sleep / Circadian',
        'iss_kit': True,
        'smiles': 'CC(=O)NCCc1c[nH]c2ccc(O)cc12',
        'dose_mg': 3,
        'lit_F': 0.15, 'lit_Vd_Lkg': 1.1, 'lit_t12_h': 0.75,
        'lit_ka_h': 2.0, 'lit_ke_h': 0.924,
        'protein_binding': 0.60, 'bcs_class': 1,
        'space_data': False, 'cyp': 'CYP1A2', 'transporter': 'None',
    },
    'Levetiracetam': {
        'category': 'Anticonvulsant',
        'iss_kit': True,
        'smiles': 'CC[C@H](N)C(=O)N1CCCC1=O',
        'dose_mg': 500,
        'lit_F': 1.00, 'lit_Vd_Lkg': 0.7, 'lit_t12_h': 7.0,
        'lit_ka_h': 2.0, 'lit_ke_h': 0.099,
        'protein_binding': 0.10, 'bcs_class': 1,
        'space_data': False, 'cyp': 'None', 'transporter': 'None',
    },
    'Ranitidine': {
        'category': 'GI (H2 blocker)',
        'iss_kit': True,
        'smiles': 'CN(C)C/C=C/CSC1=NC2=C(N1)SC(=N2)N',
        'dose_mg': 150,
        'lit_F': 0.50, 'lit_Vd_Lkg': 1.4, 'lit_t12_h': 2.5,
        'lit_ka_h': 1.8, 'lit_ke_h': 0.277,
        'protein_binding': 0.15, 'bcs_class': 3,
        'space_data': False, 'cyp': 'CYP1A2', 'transporter': 'None',
    },
}


def catalog_as_database():
    """Strip UI metadata for PBPK engine."""
    db = {}
    for name, entry in ISS_DRUG_CATALOG.items():
        row = {k: v for k, v in entry.items()
               if k not in ('category', 'iss_kit', 'aliases')}
        db[name] = row
    return db


def get_catalog_table():
    """Rows for Streamlit drug browser."""
    rows = []
    for name, e in ISS_DRUG_CATALOG.items():
        rows.append({
            'Drug': name,
            'Category': e.get('category', 'Other'),
            'ISS Kit': 'Yes' if e.get('iss_kit') else 'No',
            'Space PK Data': 'Yes' if e.get('space_data') else 'Estimated',
            'Dose (mg)': e.get('dose_mg'),
            'Route': 'IV' if e.get('lit_ka_h') == 'IV' else 'Oral',
            't½ (h)': e.get('lit_t12_h'),
            'BCS': e.get('bcs_class'),
        })
    return rows


def get_categories():
    cats = sorted({e.get('category', 'Other') for e in ISS_DRUG_CATALOG.values()})
    return cats


def search_drugs(query='', category=None, iss_only=False, space_data_only=False):
    q = query.lower().strip()
    out = []
    for name, e in ISS_DRUG_CATALOG.items():
        if category and e.get('category') != category:
            continue
        if iss_only and not e.get('iss_kit'):
            continue
        if space_data_only and not e.get('space_data'):
            continue
        aliases = e.get('aliases', [])
        hay = ' '.join([name.lower()] + [a.lower() for a in aliases])
        if q and q not in hay:
            continue
        out.append(name)
    return sorted(out)
