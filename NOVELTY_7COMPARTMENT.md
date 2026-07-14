# How SpacePK stays novel with a fixed 7-compartment PBPK

## Short answer

**We did NOT drop 7-compartment because 1-compartment is “more novel.”**  
We temporarily used 1-compartment only while the 7-compartment code was **broken** (`Vd` unused → Earth Cmax wrong by ~5–11×).

**Now the 7-compartment model is restored and fixed.** Novelty is preserved.

---

## What was wrong before (not a novelty choice)

| Issue | Effect |
|-------|--------|
| `Vd` assigned but not used for distribution | Systemic conc lived in ~3.3 L venous volume |
| Literature `F` + hepatic first-pass both applied | Double first-pass → extra underprediction |
| Arterial equation poorly closed | Weak mass balance |

That made **any** “7-compartment PBPK” claim scientifically false — reviewers would catch the Earth validation failure immediately.

---

## What “fixed 7-compartment” means now

1. **Vd-calibrated Kp**  
   `Kp = (Vd − V_blood) / V_perfused` so steady-state volume ≈ literature Vd.

2. **Honest oral bioavailability**  
   Literature `F` already includes first-pass → absorbed drug enters **systemic venous** blood; liver/kidney clear **recirculating** drug. No double first-pass.

3. **Real 7 compartments still simulated**  
   GI, portal, liver, arterial, venous, tissue, kidney — with space-modified flows/volumes by mission phase.

4. **Earth validation (paracetamol)**  
   - 500 mg: Cmax **4.65** vs Kovachevich ~5.13 (**−9.4%**)  
   - 1000 mg: Cmax **9.30** vs ~9.41 (**−1.2%**)

---

## How to write the novelty (paper language)

**Do say:**
> We present a **7-compartment PBPK framework for microgravity**, with tissue:plasma partition calibrated to literature Vd, mission-phase physiological modifiers (acute / adaptation / chronic), multi-drug ISS formulary coverage, and Bayesian PopPK synthesis of sparse astronaut/bedrest data for **Gaganyaan-relevant dosing**.

**Do not say:**
> “First ever PBPK in history” / “Neural ODE 99.8%” / “15 papers” if census says 11.

**Honest census (after CSV clean):**
- **11 unique papers**
- **10 drugs with extracted PK rows**
- ISS catalog in the app can still be **41 drugs** (model-predicted space PK) — label tiers clearly

---

## Why 7-compartment beats 1-compartment for novelty

| | 1-compartment | 7-compartment (fixed) |
|--|---------------|------------------------|
| Matches F,ka,ke,Vd | Yes | Yes (now) |
| Organ-level space physiology (hepatic/renal flow) | No | **Yes** |
| Tissue vs plasma under fluid shift | No | **Yes** |
| Reviewer expectation for “PBPK” title | Weak | **Required** |
| Gaganyaan story (where drug goes in body) | Thin | **Strong** |

**Bottom line:** Keep the **7-compartment title**. The novelty is *microgravity-modified multi-compartment physiology + mission-phase dosing + sparse-data Bayesian layer* — not the number of ODEs alone, and not a broken model with a fancy name.
