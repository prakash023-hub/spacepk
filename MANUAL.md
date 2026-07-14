# SpacePK — Plain-Language Manual (what has been done, and why)

_Last updated: 2026-07-15_

This document explains the whole project in simple terms: what SpacePK is, what
was broken, what was fixed, and exactly how to run and verify it. It is written
so that you (and any reviewer, or NASA/ISRO reader) can follow it without
digging through code.

---

## 1. What SpacePK is

SpacePK predicts how the human body handles a drug **in space** compared to
**on Earth**, and recommends a **space dose**. It does this in layers:

1. **Chemistry (RDKit):** from a drug's structure it reads molecular weight,
   logP, polar surface area, BCS class, P-gp, etc.
2. **7-compartment PBPK model:** a physiological model of the body
   (GI → portal → liver → venous → arterial → tissue → kidney). It solves
   differential equations to get drug concentration over time.
3. **Space physiology modifiers:** blood flows and volumes change by mission
   phase — acute (days 1–3), adaptation (4–14), chronic (15+).
4. **Bayesian PopPK (PyMC):** uncertainty on the key parameters (F, ka, ke, Vd)
   with credible intervals.
5. **Dose recommendation:** a single primary space dose (Cmax-matched), with an
   AUC-matched value reported as a secondary sensitivity check.

**Goal:** a defensible, reproducible tool for astronaut drug dosing —
useful for future NASA / ISRO (Gaganyaan) missions and a Q1 journal
(CPT: Pharmacometrics & Systems Pharmacology) + bioRxiv preprint.

---

## 2. The five bugs that were found and fixed

These are the problems that would have sunk the paper. All are now fixed.

### Bug 1 — The 7-compartment model didn't use Vd (underpredicted 4.5–11×)
**Before:** `Vd` (volume of distribution) and `ke` were read from the data but
never used in the equations. So the drug was effectively trapped in ~3 L of
blood, and predicted concentrations were far too low.
**Fix:** `calibrate_kp()` now sets the tissue:plasma partition coefficient (Kp)
so the model's steady-state volume equals the drug's real Vd. This is standard
PBPK practice (Rodgers & Rowland / Poulin & Theil). Vd is now **alive**.

### Bug 2 — First-pass metabolism was counted twice
**Before:** literature bioavailability `F` already includes first-pass liver
loss, but the model **also** ran the absorbed drug through the liver — double
counting, so concentrations dropped again.
**Fix:** absorbed drug (already multiplied by `F`) now enters the **venous**
blood directly. The liver only processes recirculating drug. No double penalty.

### Bug 3 — The dataset pooled different doses (the biggest hidden bug)
**Before:** the "pooled Earth Cmax = 9.41 µg/mL for 500 mg" was wrong. It
averaged **1000 mg** data (Gandia) with **500 mg** data (Kovachevich) and
**625 mg** data (Polyakov) without normalizing. Everything downstream — priors,
figures, headline numbers — inherited this error.
**Fix:** `clean_dataset.py` now adds dose-normalized columns
(`Cmax_ugmL_per_500mg`, `AUC_ughmL_per_500mg`) scaled by `500 / Dose_mg`, using
the stated assumption that paracetamol PK is approximately dose-linear over
500–1000 mg.

| Paracetamol (per 500 mg) | RAW (wrong) | NORMALIZED (correct) |
|---|---|---|
| Cmax Earth | 11.12 | **7.14 µg/mL** |
| Cmax Space | 10.53 | **6.00 µg/mL** |
| AUC Earth | 31.12 | **21.76 µg·h/mL** |
| AUC Space | 33.32 | **17.93 µg·h/mL** |

The honest finding: **space exposure is lower than Earth** for paracetamol.

### Bug 4 — No single source of truth (5 different doses for one drug)
**Before:** the app, the figure scripts, and the tables each computed the model
independently, so they disagreed (five different paracetamol space doses).
**Fix:** all PK code now lives in **`code/core.py`**. `pbpk_model.py` is a thin
shim that re-exports from it. The app and every figure script import from
`core`. They physically cannot disagree now.

### Bug 5 — Rounding destroyed low-concentration drugs; fake rescaling
**Before:** `round(x, 3)` inside the extraction collapsed drugs like
scopolamine to "no change"; and some plotting scripts rescaled the simulated
curve to the literature value (fake validation).
**Fix:** `extract_pk()` returns raw floats (rounding only happens at display).
No production code rescales a simulated curve to a literature number.

---

## 3. Data cleaning summary

`code/clean_dataset.py` produces `data/space_pk_master_clean.csv` and
`data/dataset_census.json`. It:

- Merges **Acetaminophen → Paracetamol** (same drug).
- Collapses `DelloRusso2022_*` into one paper; `Polyakov2021_SF` into `Polyakov2021`.
- Removes the invalid "General ADME review" row.
- Harmonizes all Cmax to µg/mL and AUC to µg·h/mL (fixes the 1000× mg/mL typos).
- Adds dose-normalized (per-500 mg) columns and recomputes Earth/Space anchors.

**Honest census:** 11 papers · 10 drugs (with observed data) · 190 data points.
(Previous claims of "15 papers / 41 drugs / 230 points" were wrong.)

---

## 4. Current validation (paracetamol, Earth, 500 mg)

| Parameter | Model | Literature | Status |
|---|---|---|---|
| Cmax | 4.65 µg/mL | 5.13 (Kovachevich) | ✅ −9.4% |
| t½ | 2.73 h | ~2.7 h | ✅ +1% |
| AUC | ~25.5 µg·h/mL | 21.8 (AUC0t) → higher AUC0inf | ✅ in range |
| Tmax | 1.16 h | ~1.7 h | ⚠️ runs fast (disclosed) |

**Disclosed limitation:** Tmax is early because absorption is a single
first-order process. This is stated as a limitation, not hidden. Cmax, t½, and
AUC validate acceptably for a first mechanistic space-PBPK.

---

## 5. How to run it

```bash
cd ~/spacepk

# 1. Clean + normalize the dataset (writes clean CSV + census + anchors)
python code/clean_dataset.py

# 2. Validate the model + regenerate PBPK figure
python code/pbpk_model.py

# 3. Mission-phase sensitivity figures
python code/mission_sensitivity.py

# 4. Launch the interactive app
streamlit run app.py
```

Single-source-of-truth check (must show matches only in `core.py`):

```bash
grep -rn "def pbpk_odes\|def extract_pk\|def recommend_dose" --include=*.py code/
```

---

## 6. What still remains (for the paper)

- Rerun the Bayesian layer with more draws/tuning (widen or fix `ka` prior;
  the old run showed a `ka` boundary spike = not identifiable).
- Regenerate manuscript tables from `core.py` outputs (never hand-typed).
- Write the manuscript around the honest finding (space exposure lower;
  Cmax-matched dosing) and the disclosed limitations.
- Submit to bioRxiv, then CPT:PSP.

---

## 7. Novelty (defensible claims)

1. First **multi-organ (7-compartment) PBPK** with mission-phase space
   physiology modifiers (a 1-compartment model literally cannot claim
   "multi-organ space physiology").
2. Vd-calibrated Kp so the mechanistic model reproduces literature Vd.
3. Bayesian uncertainty on the **space** dose adjustment, not just Earth.
4. An open, reproducible flight-surgeon decision-support app for
   Gaganyaan / ISS-relevant drugs.
