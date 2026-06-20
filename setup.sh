#!/bin/bash
# SpacePK — install dependencies (Mac, Anaconda recommended)
set -e
cd "$(dirname "$0")"

echo "=== SpacePK Setup ==="

if command -v conda >/dev/null 2>&1; then
  echo "Installing via conda-forge (best for RDKit on Mac)..."
  conda install -y -c conda-forge rdkit pymc numpy scipy matplotlib arviz pytensor
else
  echo "No conda found. Create venv and install pip packages..."
  python3 -m venv .venv
  source .venv/bin/activate
  pip install -U pip
  pip install -r requirements.txt
  echo ""
  echo "WARNING: RDKit still needed. Run:"
  echo "  conda install -c conda-forge rdkit"
fi

echo ""
echo "=== Test ==="
cd code
python drug_properties.py
python pbpk_model.py
python bayesian_popPK.py
echo ""
echo "Done. Figures in ../figures/"
