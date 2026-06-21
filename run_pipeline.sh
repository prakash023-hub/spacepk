#!/usr/bin/env bash
set -e
cd "$(dirname "$0")/code"
echo "=== SpacePK Full Pipeline ==="
python drug_properties.py
python pbpk_model.py
python mission_sensitivity.py
if [ "${SPACEPK_FAST:-1}" = "1" ]; then
  python bayesian_popPK.py --fast
else
  python bayesian_popPK.py --full
fi
echo "=== Done — figures in ../figures/ ==="
