"""
Clean SpacePK master dataset — single source of truth.

Fixes:
  - Acetaminophen → Paracetamol
  - Drop 'General ADME review' pseudo-drug
  - Collapse DelloRusso2022_* → DelloRusso2022; Polyakov2021_SF → Polyakov2021
  - Harmonize Cmax/AUC to µg/mL and µg·h/mL
  - Flag / correct impossible mg/mL acetaminophen labels (values are µg/mL-scale)
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "data" / "space_pk_master_v2.csv"
OUT = ROOT / "data" / "space_pk_master_clean.csv"
META = ROOT / "data" / "dataset_census.json"

PAPER_MAP = {
    "DelloRusso2022_ADME": "DelloRusso2022",
    "DelloRusso2022_Table1a": "DelloRusso2022",
    "DelloRusso2022_Table1b": "DelloRusso2022",
    "DelloRusso2022_Table1c": "DelloRusso2022",
    "Polyakov2021_SF": "Polyakov2021",
}

DRUG_MAP = {
    "Acetaminophen": "Paracetamol",
    "Scopolamine/Dextroamphetamine": "Scopolamine",
}


def _to_float(x):
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return np.nan
    s = str(x).strip()
    if s.upper() in {"", "NA", "NAN", "NONE"}:
        return np.nan
    try:
        return float(s)
    except ValueError:
        return np.nan


def convert_cmax(value, unit, drug):
    v = _to_float(value)
    if np.isnan(v):
        return np.nan
    u = (unit or "").strip().lower().replace("μ", "µ")
    # Mis-labeled: Dello Russo table values ~10 look like µg/mL, not mg/mL
    if u in {"mg/ml"} and v < 100:
        return v  # treat as µg/mL
    if u in {"µg/ml", "ug/ml", "mg/l"}:
        return v
    if u in {"mg/ml"}:
        return v * 1000.0
    if u in {"ng/ml"}:
        return v * 0.001
    if u in {"µg/l", "ug/l"}:
        return v * 0.001
    if u in {"pg/ml"}:
        return v * 1e-6
    return v


def convert_auc(value, unit):
    v = _to_float(value)
    if np.isnan(v):
        return np.nan
    u = (unit or "").strip().lower().replace("μ", "µ").replace("*", "·")
    if u in {"µg·h/ml", "ug·h/ml", "µg*h/ml", "ug*h/ml", "mg/l*h", "mg·h/l", "mg/l·h"}:
        return v
    if u in {"ng/ml*h", "ng·h/ml", "ng/ml·h"}:
        return v * 0.001
    if u in {"µg/l*h", "ug/l*h", "µg·h/l"}:
        return v * 0.001
    if u in {"mg/l*min", "mg·min/l"}:
        return v / 60.0
    return v


def clean():
    df = pd.read_csv(SRC)
    df = df[df["Drug"].astype(str).str.strip() != "General ADME review"].copy()

    df["Paper_ID"] = df["Paper_ID"].replace(PAPER_MAP)
    df["Drug"] = df["Drug"].replace(DRUG_MAP)

    df["Cmax_ug_mL"] = [
        convert_cmax(v, u, d)
        for v, u, d in zip(df["Cmax_value"], df["Cmax_unit"], df["Drug"])
    ]
    df["AUC0t_ug_h_mL"] = [
        convert_auc(v, u) for v, u in zip(df["AUC0t_value"], df["AUC0t_unit"])
    ]
    df["AUC0inf_ug_h_mL"] = [
        convert_auc(v, u) for v, u in zip(df["AUC0inf_value"], df["AUC0inf_unit"])
    ]
    df["Tmax_h"] = pd.to_numeric(df["Tmax_h"], errors="coerce")
    # Convert Tmax in minutes if values look like minutes (>24 and Notes mention min)
    # Rumble used "35.6 min" stored as number with unit elsewhere — leave numeric as-is if < 24

    # Rumble rows: Tmax_h column actually minutes in source — fix if Dose 500 and Tmax > 24
    mask_rumble = (df["Paper_ID"] == "Rumble1991") & (df["Tmax_h"] > 24)
    df.loc[mask_rumble, "Tmax_h"] = df.loc[mask_rumble, "Tmax_h"] / 60.0

    # Half-life: Rumble stored minutes
    df["t_half_h"] = pd.to_numeric(df["t_half_h"], errors="coerce")
    mask_t12 = (df["Paper_ID"] == "Rumble1991") & (df["t_half_h"] > 24)
    df.loc[mask_t12, "t_half_h"] = df.loc[mask_t12, "t_half_h"] / 60.0

    df["Cmax_unit_harmonized"] = "µg/mL"
    df["AUC_unit_harmonized"] = "µg·h/mL"

    # ── DOSE NORMALISATION (critical) ─────────────────────────────────────────
    # Studies pool DIFFERENT doses (Gandia 1000 mg, Kovachevich 500 mg,
    # Polyakov 625 mg). Pooling raw Cmax/AUC is invalid. Paracetamol PK is
    # approximately dose-LINEAR over 500–1000 mg, so we normalise to a common
    # 500 mg reference before pooling / anchoring.
    dose = pd.to_numeric(df["Dose_mg"], errors="coerce")
    scale = 500.0 / dose.where(dose > 0, np.nan)
    df["Cmax_ugmL_per_500mg"] = df["Cmax_ug_mL"] * scale
    df["AUC_ughmL_per_500mg"] = df["AUC0t_ug_h_mL"].fillna(df["AUC0inf_ug_h_mL"]) * scale

    # Count data points (non-null PK fields)
    point_cols = ["Cmax_ug_mL", "Tmax_h", "AUC0t_ug_h_mL", "AUC0inf_ug_h_mL", "t_half_h"]
    n_points = int(df[point_cols].notna().sum().sum())

    census = {
        "source": str(SRC.name),
        "rows": int(len(df)),
        "unique_papers": int(df["Paper_ID"].nunique()),
        "papers": sorted(df["Paper_ID"].dropna().unique().tolist()),
        "unique_drugs": int(df["Drug"].nunique()),
        "drugs": sorted(df["Drug"].dropna().unique().tolist()),
        "data_points": n_points,
        "notes": [
            "Acetaminophen merged into Paracetamol",
            "DelloRusso2022_* collapsed to one paper",
            "Polyakov2021_SF collapsed into Polyakov2021",
            "General ADME review row removed",
            "All Cmax reported as µg/mL; mg/mL labels with Cmax<100 treated as unit typos",
        ],
    }

    # ── Recompute Earth vs Space anchors: RAW vs NORMALISED ───────────────────
    para = df[df["Drug"] == "Paracetamol"].copy()
    cond = para["Condition"].astype(str).str.lower()
    is_space = cond.str.contains("space|iss|mir|in-flight|hdt day|bed rest|anoh|flight day", regex=True) \
        & ~cond.str.contains("earth|ambulat|ground|pre-flight|before br|normal|1g|background", regex=True)
    is_earth = ~is_space

    def _mean(series, mask):
        v = pd.to_numeric(series[mask], errors="coerce").dropna()
        return round(float(v.mean()), 3) if len(v) else None

    anchors = {
        "paracetamol_Cmax_ugmL": {
            "raw_earth": _mean(para["Cmax_ug_mL"], is_earth),
            "raw_space": _mean(para["Cmax_ug_mL"], is_space),
            "norm500_earth": _mean(para["Cmax_ugmL_per_500mg"], is_earth),
            "norm500_space": _mean(para["Cmax_ugmL_per_500mg"], is_space),
        },
        "paracetamol_AUC_ughmL": {
            "raw_earth": _mean(para["AUC0t_ug_h_mL"], is_earth),
            "raw_space": _mean(para["AUC0t_ug_h_mL"], is_space),
            "norm500_earth": _mean(para["AUC_ughmL_per_500mg"], is_earth),
            "norm500_space": _mean(para["AUC_ughmL_per_500mg"], is_space),
        },
        "linearity_assumption": "Paracetamol PK approximately dose-linear 500–1000 mg; Cmax/AUC scaled by 500/Dose_mg.",
    }
    census["anchors"] = anchors

    df.to_csv(OUT, index=False)
    META.write_text(json.dumps(census, indent=2))
    print(json.dumps(census, indent=2))
    print("\n── Paracetamol anchors: RAW (invalid pooling) vs NORMALISED to 500 mg ──")
    print(json.dumps(anchors, indent=2))
    print(f"\nWrote {OUT}")
    return df, census


if __name__ == "__main__":
    clean()
