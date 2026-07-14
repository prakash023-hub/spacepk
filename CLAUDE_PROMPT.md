# Prompt to paste into Claude (project handoff)

Copy everything inside the code block below into a new Claude chat. It gives
Claude the full, accurate state of the project so it can help without
re-introducing old bugs.

```
You are helping me finish SpacePK, a computational pharmacokinetics project for
Q1 publication (CPT: Pharmacometrics & Systems Pharmacology) + a bioRxiv
preprint, and as a demo for future NASA / ISRO (Gaganyaan) drug-dosing support.

WHAT SPACEPK DOES
It predicts how a drug behaves in space vs on Earth and recommends a space dose.
Layers: (1) RDKit chemistry, (2) a 7-compartment PBPK ODE model
(GI→portal→liver→venous→arterial→tissue→kidney), (3) mission-phase space
physiology modifiers (acute/adaptation/chronic), (4) Bayesian PopPK (PyMC),
(5) dose recommendation (Cmax-matched primary, AUC-matched secondary).

SINGLE SOURCE OF TRUTH (do not violate this)
ALL pharmacokinetic code lives in code/core.py:
  get_physiology, calibrate_kp, pbpk_odes, extract_pk, recommend_dose,
  run_pbpk_analysis.
code/pbpk_model.py is only a compatibility shim that re-exports from core.py.
app.py and every figure/table script MUST import from core. Never re-implement
the model anywhere else. Verify with:
  grep -rn "def pbpk_odes\|def extract_pk\|def recommend_dose" --include=*.py code/
It must return matches ONLY inside code/core.py.

BUGS ALREADY FIXED — do not reintroduce them
1. Vd must stay ALIVE: distribution volume is enforced via calibrate_kp() so
   steady-state Vss ≈ literature Vd (Rodgers & Rowland practice). Do NOT
   restructure the model so Vd stops affecting concentrations.
2. NO double first-pass: literature F already includes first-pass loss, so
   absorbed drug (ka*A_gi*F) enters VENOUS blood directly. The liver only
   processes recirculating drug. Do not route absorption through the liver.
3. DOSE NORMALIZATION: never pool Cmax/AUC across studies at different doses.
   Use the per-500 mg columns (Cmax_ugmL_per_500mg = Cmax * 500/Dose_mg),
   assuming approximate dose-linearity for paracetamol over 500–1000 mg.
4. NEVER rescale a simulated curve to a literature value (e.g. 9.41/C.max()).
   If validation fails, fix the model — do not fake it.
5. Do NOT round inside extract_pk(); round only at the display layer
   (rounding there previously destroyed low-conc drugs like scopolamine).
6. recommend_dose() returns ONE primary number (default criterion='cmax');
   AUC-matched is secondary sensitivity only.

HONEST DATASET FACTS (do not inflate)
Cleaned dataset: 11 papers, 10 drugs with observed data, 190 data points.
(Old claims of 15 papers / 41 drugs / 230 points were WRONG.)
Dose-normalized paracetamol anchors (per 500 mg):
  Cmax Earth 7.14, Space 6.00 µg/mL ; AUC Earth 21.76, Space 17.93 µg·h/mL.
Space exposure is LOWER than Earth — that is the real finding.

CURRENT VALIDATION (paracetamol, Earth, 500 mg)
Cmax 4.65 vs 5.13 lit (−9.4%) ✓ ; t½ 2.73 vs ~2.7 ✓ ; AUC ~25.5 (in AUC0t→
AUC0inf range) ✓ ; Tmax 1.16 vs ~1.7 — runs fast, DISCLOSED limitation.

WHAT I NEED HELP WITH NEXT
- Improve the Bayesian layer (more draws/tune; fix or widen the ka prior — the
  old run had a ka boundary spike = not identifiable).
- Regenerate all manuscript tables programmatically from core.py (never typed).
- Draft the manuscript around the honest finding + disclosed limitations.
- Prep bioRxiv + CPT:PSP submission.

RULES FOR YOU
- Always read code/core.py and MANUAL.md before changing anything.
- Keep everything reproducible and honest; disclose limitations rather than
  hiding them. This is for a real journal and space agencies.
- Ask me before changing scope or deleting data.
```
