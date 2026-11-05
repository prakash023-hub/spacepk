# SpacePK

[![Streamlit App](https://static.streamlit.io/badges/streamlit_badge_black_white.svg)](https://prakash023-hub-spacepk-app-kua6sq.streamlit.app)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)

**An open, reproducible cheminformatics → PBPK → Bayesian framework for comparing Earth and spaceflight pharmacokinetics.**

> Research prototype and hypothesis-generating tool. **Not a clinical dosing system.** Not validated in astronauts; not approved by any space agency or medical authority.

## Live demo
https://prakash023-hub-spacepk-app-kua6sq.streamlit.app

## What it does
Select a drug + mission day + body weight → get Earth vs space concentration–time curves, mission-phase classification, and a model-derived (Cmax-matched) dose evaluation with uncertainty.

## Pipeline
- **Layer 1 — Cheminformatics:** RDKit molecular descriptors → BCS class, P-gp flag, Lipinski.
- **Layer 2 — 7-compartment PBPK:** GI → portal → liver → arterial → venous → tissue → kidney, with mission-phase (acute/adaptation/chronic) physiological modifiers. Vd-calibrated Kp; no double first-pass.
- **Layer 3 — Bayesian PK estimation (PyMC, pooled):** exact one-compartment oral model fitted to study-arm Cmax, Tmax and AUC0-inf, each at its own dose (assumes linear PK 500–1000 mg). Reports ka, ke, V/F, CL/F with 95% credible intervals (F and V are not separately identifiable from these data). Shared priors for Earth and space so differences come from data. Pooled, not hierarchical PopPK — subject-level random effects are the next step. **Preliminary** (6 Earth / 4 in-flight study arms; ka weakly identifiable).
- **Dose evaluation:** Cmax-matched primary, AUC-matched secondary.

## Dataset (honest census)
`data/space_pk_master_clean.csv` — cleaned analysis dataset. See `data/dataset_census.json` for reproducible counts.

- **11 peer-reviewed studies**, **10 drugs with observed PK**, **190 observed parameter values**.
- Real in-flight PK data exist for **paracetamol** (Shuttle/ISS/MIR) plus **one scopolamine/dextroamphetamine** Cmax (Shuttle).
- All other observed drugs are **HDT / bed-rest analogue**, identifiable from the `Gravity_Model` / `Condition` columns (e.g. `Real ISS`, `Real spaceflight` vs `HDT -6 degrees`).
- The wider ISS-formulary catalogue in the app is **model-predicted** and labelled as such — it is not observed data.

Evidence tiers used above: `earth_observed`, `spaceflight_real`, `analog_hdt_bedrest`, `model_predicted` — a narrative classification of the studies, not a column in the CSV.

## Status of results
- **Earth calibration (paracetamol 500 mg):** model Cmax 4.65 vs 5.13 µg/mL literature (−9.4%); t½ +1%; AUC in range; Tmax predicted early (disclosed limitation).
- **Earth vs space (dose-normalised):** paracetamol space exposure **lower** than Earth (directional; underlying data sparse and inconsistent).
- **Bayesian layer (500 mg reference):** Cmax-matched space dose factor ×1.27 (95% CrI 0.86–1.84); AUC-matched ×1.31 (0.81–2.06). Credible intervals include 1, so no dose change is supported by current data.
- Bayesian posteriors and multi-drug tables are regenerated from source; do not hand-edit.

## Reproduce
```bash
pip install -r requirements.txt
python code/clean_dataset.py      # writes data/space_pk_master_clean.csv + dataset_census.json
python code/pbpk_model.py         # Earth calibration + PBPK figures
python code/mission_sensitivity.py
streamlit run app.py
```

## Paper
Manuscript in preparation; preprint planned on bioRxiv prior to journal submission.

## Author
K. Prakash Raj, M.Pharm — Sri Balaji Vidyapeeth, Puducherry, India
GitHub: [@prakash023-hub](https://github.com/prakash023-hub)

## Citation
See `CITATION.cff`.

## Licence
MIT — see `LICENSE`.

## Acknowledgement
AI-assisted tools (Claude, Anthropic) were used for code structuring, data cleaning, and manuscript drafting. All scientific interpretation and conclusions are the author's own.
