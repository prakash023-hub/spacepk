# SpacePK Framework

**Integrated cheminformatics → PBPK → Bayesian PopPK for Earth vs spaceflight pharmacokinetics**

Physiologically-based pharmacokinetic modeling with mission-phase space physiology modifiers and Bayesian uncertainty quantification for astronaut dose adjustment.

## Novel contribution

SpacePK connects three layers that are usually published separately:

1. **RDKit molecular descriptors** → PBPK parameterization (BCS, P-gp, Lipinski)
2. **7-compartment mechanistic PBPK** with acute / adaptation / chronic mission-phase physiology
3. **Bayesian PopPK** with 95% credible intervals on space dose adjustment

Literature anchors: Gandia 2003, Kovachevich 2009, Polyakov 2021, Dello Russo 2022.

> **Dose normalization (important):** pooled anchors are normalized to a common
> 500 mg reference (`Cmax × 500/Dose`) because the source studies use different
> doses (Gandia 1 g, Kovachevich 500 mg, Polyakov 625 mg). Paracetamol PK is
> approximately dose-linear over 500–1000 mg. All PK code lives in a single
> source of truth, `code/core.py`; every figure, table, and the app import it.

## Structure

```
spacepk/
├── app.py              # Streamlit flight-surgeon decision support UI
├── run_app.sh          # Launch interactive app
├── run_pipeline.sh     # Run all layers + sensitivity
├── code/
│   ├── core.py                 # SINGLE SOURCE OF TRUTH — PBPK engine + dose logic
│   ├── drug_properties.py      # Layer 1 — RDKit + space modifiers
│   ├── pbpk_model.py           # Layer 2 — compatibility shim → core.py
│   ├── clean_dataset.py        # Data harmonization + dose normalization
│   ├── bayesian_popPK.py       # Layer 3 — PyMC Bayesian PopPK
│   └── mission_sensitivity.py  # Mission-day sensitivity analysis
├── data/               # PK master datasets (cleaned: 11 papers, 10 drugs, 190 points)
├── figures/            # Output plots
└── papers/             # Source PDFs
```

## Quick start

### Setup (Mac + Anaconda)

```bash
conda install -c conda-forge rdkit pymc numpy scipy matplotlib arviz pytensor
pip install "arviz>=0.17,<1.0" streamlit
cd ~/spacepk
chmod +x run_app.sh run_pipeline.sh
```

### Run full research pipeline

```bash
./run_pipeline.sh
```

Publication-quality Bayesian run (slower):

```bash
SPACEPK_FAST=0 ./run_pipeline.sh
```

### Launch interactive app

```bash
./run_app.sh
```

Or:

```bash
streamlit run app.py
```

## Layers

| Layer | Script | Output |
|-------|--------|--------|
| 1 | `drug_properties.py` | Molecular properties + space-adjusted PK inputs |
| 2 | `pbpk_model.py` | Earth vs space concentration curves + dose recommendation |
| 3 | `bayesian_popPK.py` | Posterior distributions + dose 95% CI |
| + | `mission_sensitivity.py` | Dose factor vs mission day (1–180) |

## Drugs modeled

The ISS / spaceflight medical catalog includes analgesics (Paracetamol, Ibuprofen, Tramadol), antibiotics (Ciprofloxacin, Amoxicillin, Doxycycline), antiemetics (Promethazine, Scopolamine, Ondansetron), and cardiovascular agents (Nifedipine, Verapamil, Propranolol, Furosemide) among others.

**10 drugs** have direct spaceflight / head-down-tilt (HDT) PK literature in the cleaned dataset (11 papers, 190 data points). Additional catalog compounds are modeled with literature-estimated PK parameters and are clearly labeled as predicted, not observed.

## GitHub

https://github.com/prakash023-hub/spacepk

## URLs for publication / submission

| URL | Use in paper? |
|-----|----------------|
| `http://localhost:8501` | **No** — only on your Mac |
| `http://192.168.x.x:8501` | **No** — home Wi‑Fi only |
| `http://117.254.x.x:8501` | **No** — your public IP; dies when laptop sleeps |
| **https://github.com/prakash023-hub/spacepk** | **Yes** — code & reproducibility |
| **https://prakash023-hub-spacepk-app-kua6sq.streamlit.app** | **Yes** — live demo (after deploy) |

### Deploy public demo (Streamlit Community Cloud)

1. Push repo to GitHub (already done)
2. Go to [share.streamlit.io](https://share.streamlit.io) → Sign in with GitHub
3. **New app** → Repository `prakash023-hub/spacepk` → Branch `main` → Main file `app.py`
4. Advanced settings → use `environment.yml` (includes RDKit via conda-forge)
5. Deploy → copy URL like `https://spacepk.streamlit.app` into your paper/Devpost

## Citation (draft)

> SpacePK: An Integrated Cheminformatics–PBPK–Bayesian Framework for Earth–Spaceflight Pharmacokinetic Comparison. CPT: Pharmacometrics & Systems Pharmacology (in preparation).

## Limitations

- Small N in real spaceflight PK studies; Bayesian priors carry uncertainty
- HDT is an Earth analog — not identical to microgravity
- No formulation stability (radiation, humidity) or drug–drug interactions yet
