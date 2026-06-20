# SpacePK Framework

Physiologically-based pharmacokinetic (PBPK) modeling for Earth vs spaceflight conditions.

## Structure

```
spacepk/
├── code/           # Layer 1–3 Python pipeline
├── data/           # PK master datasets (CSV)
├── figures/        # Output plots
├── papers/         # Source PDFs
└── outputs/        # Generated results
```

## Layers

| Layer | File | Purpose |
|-------|------|---------|
| 1 | `drug_properties.py` | RDKit molecular properties → PBPK inputs |
| 2 | `pbpk_model.py` | 7-compartment PBPK, Earth vs space |
| 3 | `bayesian_popPK.py` | Bayesian PopPK (PyMC) |

## Setup (Mac)

```bash
cd ~/spacepk
chmod +x setup.sh
./setup.sh
```

Or manually with Anaconda:

```bash
conda install -c conda-forge rdkit pymc numpy scipy matplotlib arviz pytensor
cd ~/spacepk/code
python drug_properties.py
python pbpk_model.py
python bayesian_popPK.py
```

## Data

- `data/space_pk_master_v2.csv` — 63 rows, 15 papers
- `data/paracetamol_pk_clean.csv` — paracetamol subset

## Git

```bash
git remote add origin https://github.com/YOUR_USERNAME/spacepk.git
git push -u origin main
```

**Do not commit:** large PDFs if repo size is an issue — use Git LFS or omit papers/
