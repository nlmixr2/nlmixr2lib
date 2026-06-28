# Canonical compartment and metabolite-suffix names

This file is the authoritative register of compartment / state names and metabolite-suffix tokens used in nlmixr2lib models. Every compartment that appears in a model's `model()` block (and every metabolite / sibling-drug suffix used to derive a compartment / parameter / residual-SD name) is expected to match one of the canonical entries below — or to fit one of the regex constants documented in the header section. The register is seeded from the 2026-05-28 naming audit and extended whenever a new paper introduces a state that isn't yet registered.

## How to use this register

1. **Before adding a compartment / state to a new model**, search this file (by canonical name and by source alias) for the concept you need.
2. **If the canonical name exists**, use it exactly. Document the source-paper rename in a code comment if the paper used a different name.
3. **If the source paper uses an alias listed under an existing canonical name**, prefer the canonical name. Aliases are documented for cross-reference, not as a free pass to introduce the deprecated form in new models.
4. **If the state is not in this register at all**, propose a new entry with a canonical name, type, role, source aliases, and example models. Verify with the user before committing. The addition is part of the model's PR.
5. **Do not modify existing model files when you discover a missing entry**; simply register the canonical here. Retrofitting existing models is a separate effort.

## Entry schema

Each canonical entry is an H3 heading whose first whitespace-separated token (before the parenthetical) is the canonical name. Required fields:

```yaml
- name: <CANONICAL_NAME>
  type: compartment | metabolite-suffix
  role: <one-sentence description of mechanistic role>
  source_aliases:
    - <ALIAS_NAME> -- used in <model.R>
  example_models:
    - <model.R>
  notes: <free text>
```

The `Type:` field is the routing tag the runtime parser uses to assign the entry to the appropriate static vector:

- `compartment` → `compartments` (used by `.matchesCompartment()` and the compartment-name validator)
- `metabolite-suffix` → `registeredMetabolites` (used by `.endsWithMetabolite()` and the `<canonical>_<suffix>` compartment / parameter / residual-SD checks)

A single token can appear under both Types (e.g., `lzd` is both a bare drug-state PK compartment and a metabolite suffix for linezolid). In that case the file carries two H3 entries with the same canonical name but distinct `Type:` fields; the runtime parser routes each entry to its respective list.

## Regex constants (kept in R, not migrated)

The following pattern constants remain hard-coded in `R/conventions.R::.nlmixr2libConventionsStatic` because they are structural regular expressions rather than name lists:

- `compartmentRegex = "^(transit|effect|precursor|lat|depot)[0-9]+$"` — numbered-chain compartments: transit absorption chains (`transit1`, `transit2`, ...), effect-compartment chains (`effect1`, `effect2`, ...), precursor pools for delayed-feedback IDR (`precursor1`, `precursor2`, ...), latent chains (`lat1`, ...), parallel-absorption depots (`depot1`, `depot2`, ...). Numeric suffix is required (single-state members use the bare canonical `effect` / `depot`).
- `darCompartmentRegex = "^dar[0-9]+_(central|peripheral[0-9]?)$"` — DAR-numbered ADC isoform compartments (`dar0_central`, `dar4_peripheral1`, ...).
- `targetLocationRegex = "^(target|complex)_(csf|isf|peripheral[0-9]?)$"` — target species in physiologic / numbered-peripheral compartments (`target_csf`, `target_isf`, `target_peripheral`, `target_peripheral1`, `complex_peripheral`, ...).
- `pbpkSubCompartmentRegex = "^(bc|eu|eb|fr|is|int|mrna|luc)_(liver|lung|kidney|spleen|heart|muscle|skin|adipose|bone|brain|small_intestine|large_intestine|pancreas|thymus|portal|remainder|other|hepatic|fat|rapidly_perfused|slowly_perfused|venous|arterial|urine|gut)$` — membrane-limited PBPK sub-compartments: vascular blood cells (`bc_`), endosomal unbound (`eu_`), endosomal FcRn-bound (`eb_`), endosomal free FcRn (`fr_`), interstitial space (`is_`), intracellular (`int_`), mRNA pool (`mrna_`), luciferase reporter (`luc_`).
- `compartmentRegex` and the four extension patterns above are extended only when a new paper introduces a structurally new shape. Adding a new spelled-out organ to the `pbpkSubCompartmentRegex` is a routine extension; introducing a new chain prefix is a naming-audit decision.

---

## Standard drug-disposition compartments

### depot (**canonical absorption depot**)
- **Type:** compartment
- **Role:** Absorption / dosing depot upstream of the central compartment. Receives dosing events for oral, subcutaneous, intramuscular, or other extravascular routes.
- **Source aliases:** none.
- **Example models:** universal in oral / subcutaneous popPK extractions.
- **Notes:** Numbered variants `depot1`, `depot2`, ... are accepted via `compartmentRegex` for parallel-absorption models. Route-specific variants (`depot_im`, `depot_oral`, `depot_brain`) are registered separately.

### central (**canonical central compartment**)
- **Type:** compartment
- **Role:** Central compartment holding the drug pool in plasma / systemic circulation; the conventional output state for plasma concentration `Cc = central / vc`.
- **Source aliases:** none.
- **Example models:** universal in popPK extractions.

### peripheral1 (**canonical first peripheral compartment**)
- **Type:** compartment
- **Role:** First peripheral / tissue distribution compartment in 2- and 3-compartment models. Connected to central via inter-compartmental clearance `q` (and volume `vp`).
- **Source aliases:** none.
- **Example models:** any 2-compartment or higher popPK extraction.

### peripheral2 (**canonical second peripheral compartment**)
- **Type:** compartment
- **Role:** Second peripheral compartment in 3- and 4-compartment models. Connected to central via `q2` / `vp2`.
- **Source aliases:** none.
- **Example models:** 3-compartment popPK extractions.

### peripheral3 (**canonical third peripheral compartment**)
- **Type:** compartment
- **Role:** Third peripheral compartment in 4-compartment models.
- **Source aliases:** none.
- **Example models:** 4-compartment popPK extractions.

### effect (**canonical effect compartment**)
- **Type:** compartment
- **Role:** Generic effect compartment (Sheiner 1979) used to introduce a hysteresis between plasma concentration and PD response.
- **Source aliases:** none.
- **Example models:** PK-PD models with effect-compartment hysteresis.
- **Notes:** Numbered variants `effect1`, `effect2`, ... are accepted via `compartmentRegex` for multi-effect chains.

### target (**canonical target compartment**)
- **Type:** compartment
- **Role:** Free target species in TMDD / receptor-binding models (amount).
- **Source aliases:** none.
- **Example models:** TMDD popPK extractions.
- **Notes:** Location-suffixed variants (`target_csf`, `target_isf`, `target_peripheral`) are accepted via `targetLocationRegex`.

### complex (**canonical drug-target complex compartment**)
- **Type:** compartment
- **Role:** Drug-target complex in TMDD / receptor-binding models (amount).
- **Source aliases:** none.
- **Example models:** TMDD popPK extractions.
- **Notes:** Location-suffixed variants (`complex_csf`, `complex_isf`, `complex_peripheral`) are accepted via `targetLocationRegex`.

### total_target (**canonical total target compartment**)
- **Type:** compartment
- **Role:** Total target species (free + drug-bound complex) used by quasi-steady-state TMDD parameterisations.
- **Source aliases:** none.
- **Example models:** QSS-TMDD popPK extractions.

---

## Semi-physiological organ states

### liver (**canonical liver compartment**)
- **Type:** compartment
- **Role:** Liver organ state used by paper-specific extraction-ratio first-pass models and whole-organ PBPK extractions.
- **Source aliases:**
  - `liv` -- deprecated.
- **Example models:** `Xie_2019_agomelatine.R`, `Ayyar_2024_givosiran.R`, `Gilkey_2015_DiRnanoparticle.R`.
- **Notes:** Always use the full English name; never `liv`.

### kidney (**canonical kidney compartment**)
- **Type:** compartment
- **Role:** Kidney organ state used by paper-specific extraction-ratio first-pass models and whole-organ PBPK extractions.
- **Source aliases:**
  - `kid` -- deprecated.
- **Example models:** `Ayyar_2024_givosiran.R`, `Gilkey_2015_DiRnanoparticle.R`.

### cumhaz (**canonical cumulative-hazard state, single-hazard models**)
- **Type:** compartment
- **Role:** Cumulative-hazard state for time-to-event / dropout sub-models. Integrates instantaneous hazard so that `survival = exp(-cumhaz)`. Use the bare `cumhaz` name when a model has only one hazard; reserve `cumhaz_<type>` (e.g., `cumhaz_os`, `cumhaz_drop`) for multi-hazard models that need to disambiguate.
- **Source aliases:**
  - `cumHazard` -- prior canonical name for the generic single-hazard form used in `Zecchin_2016.R` (pre-2026-06-19 lowercase + drop-suffix standardization).
- **Example models:** `Girard_2012_pimasertib.R`, `Zecchin_2016.R`.
- **Notes:** Source NONMEM idiom is `$MODEL COMP=(CUMHAZ)` with `DADT(<cumhaz>) = HAZARD`. The pre-2026-06-19 register carried a separate `cumHazard` canonical for single-hazard models (used in `Zecchin_2016.R`); the 2026-06-19 audit collapsed it into this entry under the operator's rule that single-hazard models drop any suffix and that all cumulative-hazard compartment names are uniformly lowercase.

### renal_cortex (**canonical renal-cortex accumulation compartment**)
- **Type:** compartment
- **Role:** Renal-cortex accumulation compartment used by aminoglycoside nephrotoxicity models. Tracks drug amount sequestered in the renal cortex via saturable uptake from central plus first-order tubular reabsorption back out (Rougier 2003 / Croes 2011 mechanism).
- **Source aliases:** none.
- **Example models:** `Llanos-Paez_2017_gentamicin.R`.

### csf (**canonical cerebrospinal-fluid compartment**)
- **Type:** compartment
- **Role:** Cerebrospinal-fluid physiologic compartment used by mechanistic mAb / TMDD models with multiple body-fluid distribution volumes.
- **Source aliases:** none.
- **Example models:** `Perez-Ruixo_2025_posdinemab.R`.

### isf (**canonical interstitial-fluid compartment**)
- **Type:** compartment
- **Role:** Interstitial-fluid physiologic compartment used by mechanistic mAb / TMDD models with multiple body-fluid distribution volumes.
- **Source aliases:** none.
- **Example models:** `Perez-Ruixo_2025_posdinemab.R`.

### ecf (**canonical brain / tumor extracellular-fluid compartment**)
- **Type:** compartment
- **Role:** Brain / tumor extracellular-fluid compartment used by cerebral-microdialysis-based CNS-distribution popPK models. Carries unbound drug delivered via influx (`clin`) and efflux (`clef`) clearances.
- **Source aliases:** none.
- **Example models:** `Campagne_2019_cyclophosphamide_mouse.R`.

---

## Brain-region namespace

The `brain_<region>` namespace was adopted 2026-05-28 to disambiguate brain-anatomical compartments from same-named non-brain compartments (e.g., the bare `cortex` and `renal_cortex`). Each state holds the extracellular drug concentration in the named region; total brain concentration including residual plasma is derived as `Cbrain_<region>` in `model()`. Bare region names (`cerebellum`, `hippocampus`, `striatum`, `choroid_plexus`, `brain_ecf`) are deprecated in favour of the prefixed forms.

### brain_cerebellum (**canonical cerebellum compartment**)
- **Type:** compartment
- **Role:** Cerebellum extracellular drug compartment in mAb brain-distribution PK and brain-PBPK extractions.
- **Source aliases:**
  - `cerebellum` -- deprecated bare form.
- **Example models:** `Grimm_2023_gantenerumab.R`, `Grimm_2023_trontinemab.R`.

### brain_hippocampus (**canonical hippocampus compartment**)
- **Type:** compartment
- **Role:** Hippocampus extracellular drug compartment.
- **Source aliases:**
  - `hippocampus` -- deprecated bare form.
- **Example models:** `Grimm_2023_gantenerumab.R`, `Grimm_2023_trontinemab.R`.

### brain_striatum (**canonical striatum compartment**)
- **Type:** compartment
- **Role:** Striatum extracellular drug compartment.
- **Source aliases:**
  - `striatum` -- deprecated bare form.
- **Example models:** `Stevens_2012_remoxipride.R`, `Grimm_2023_gantenerumab.R`.

### brain_cortex (**canonical brain cortex compartment**)
- **Type:** compartment
- **Role:** Brain cortex extracellular drug compartment.
- **Source aliases:** none.
- **Example models:** `Grimm_2023_gantenerumab.R`, `Grimm_2023_trontinemab.R`.

### brain_choroid_plexus (**canonical choroid plexus compartment**)
- **Type:** compartment
- **Role:** Choroid plexus extracellular drug compartment.
- **Source aliases:**
  - `choroid_plexus` -- deprecated bare form.
- **Example models:** `Grimm_2023_gantenerumab.R`.

### brain_csf (**canonical brain CSF compartment**)
- **Type:** compartment
- **Role:** Brain cerebrospinal-fluid compartment for brain-PBPK and brain-distribution PK models. Use this name for the single lumped CSF compartment when a model does not need to resolve CSF anatomical subdivisions; use the `brain_csf_<location>` family below when a model distinguishes CSF subregions (LV / TFV / CM / SAS) mechanistically.
- **Source aliases:** none.
- **Example models:** `Grimm_2023_gantenerumab.R`, `Grimm_2023_trontinemab.R`, `Syvanen_2012_quinidine_rat.R` (where `brain_csf` carries the single brain extracellular sampling compartment in the absence of an ECF / CSF distinction).

### brain_ecf (**canonical brain parenchymal extracellular-fluid compartment**)
- **Type:** compartment
- **Role:** Brain parenchymal extracellular fluid (interstitial fluid surrounding brain cells) in systems-based PK (SBPK) and mechanistic brain-distribution models that distinguish brain parenchymal ECF from CSF as separate mass-balance compartments. Drug enters via passive and / or transporter-mediated clearance from plasma and flows toward the CSF ventricular system at the physiological brain-ECF flow rate `Q_ECF` (rat: ~0.2 uL/min).
- **Source aliases:** none.
- **Example models:** `Westerhout_2013_quinidine.R` (founding example; rat SBPK model jointly fit to parenchymal brain ECF, four CSF subregions, and total brain).
- **Notes:** Distinct from `brain_csf` (the cerebrospinal-fluid compartment) and from the bare `ecf` (the general brain / tumor extracellular-fluid compartment used by simpler microdialysis models without an explicit CSF compartment, e.g., `Campagne_2019_cyclophosphamide_mouse.R`). Pre-2026-06-24 `brain_ecf` was a deprecated alias of `brain_csf`; the operator restored the name 2026-06-24 (response to task 215 sidecar) under the role-based `brain_<region>` namespace so SBPK papers that distinguish parenchymal ECF from CSF can encode both. Existing pre-2026-06-24 models using `brain_csf` for the parenchymal extracellular compartment (e.g., Syvanen 2012) are NOT retrofitted.

### brain_csf_lv (**canonical CSF lateral-ventricle compartment**)
- **Type:** compartment
- **Role:** Cerebrospinal-fluid compartment in the lateral ventricle, the first CSF compartment downstream of choroid-plexus secretion in the rat ventricular system. Used by mechanistic intra-brain SBPK models that resolve CSF anatomical subregions; drug enters via blood-CSF barrier transport from plasma and via parenchymal brain-ECF flow, and exits to `brain_csf_tfv` at the CSF flow rate `Q_CSF` (rat: ~2.2 uL/min).
- **Source aliases:** none.
- **Example models:** `Westerhout_2013_quinidine.R` (founding example; physiological volume `V_LV` = 50 uL in the rat).

### brain_csf_tfv (**canonical CSF third + fourth ventricle compartment**)
- **Type:** compartment
- **Role:** Combined third and fourth ventricle CSF compartment in mechanistic intra-brain SBPK models. Receives CSF flow from `brain_csf_lv` and discharges into `brain_csf_cm` at `Q_CSF`. Drug enters via blood-CSF barrier transport from plasma.
- **Source aliases:** none.
- **Example models:** `Westerhout_2013_quinidine.R` (founding example; physiological volume `V_TFV` = 50 uL in the rat; in Westerhout 2013 the plasma-to-TFV transfer clearance was structurally assumed equal to the plasma-to-LV transfer clearance because no TFV microdialysis was performed).

### brain_csf_cm (**canonical CSF cisterna-magna compartment**)
- **Type:** compartment
- **Role:** Cisterna-magna CSF compartment in mechanistic intra-brain SBPK models. Receives CSF flow from `brain_csf_tfv` and discharges into `brain_csf_sas` at `Q_CSF`. Drug enters via blood-CSF barrier transport from plasma; commonly the second microdialysis sampling site in multi-probe rat SBPK designs.
- **Source aliases:** none.
- **Example models:** `Westerhout_2013_quinidine.R` (founding example; physiological volume `V_CM` = 17 uL in the rat).

### brain_csf_sas (**canonical CSF subarachnoid-space compartment**)
- **Type:** compartment
- **Role:** Subarachnoid-space CSF compartment (cranial + spinal SAS combined) in mechanistic intra-brain SBPK models. Terminal CSF compartment that returns drug to systemic plasma via arachnoid villi at `Q_CSF`.
- **Source aliases:** none.
- **Example models:** `Westerhout_2013_quinidine.R` (founding example; physiological volume `V_SAS` = 180 uL in the rat; this compartment closes the CSF loop by feeding back into the plasma mass balance).

### brain_deep (**canonical deep brain compartment**)
- **Type:** compartment
- **Role:** Deep brain tissue compartment in brain-PBPK extractions.
- **Source aliases:** none.
- **Example models:** `Grimm_2023_gantenerumab.R`.

### brain_vascular (**canonical brain vascular compartment**)
- **Type:** compartment
- **Role:** Drug in cerebral capillary blood (volume `Vbv`, fed by cerebral blood flow `CLbv` from systemic central) in hybrid physiology-based PK-PD models that resolve the blood-brain barrier transport as two coupled states.
- **Source aliases:** none.
- **Example models:** `Johnson_2011_olanzapine_rat.R`.

### brain_extravascular (**canonical brain extravascular compartment**)
- **Type:** compartment
- **Role:** Drug in brain tissue beyond the BBB (volume `Vbev`, fed via the BBB-clearance term `CLbev`). Paired with `brain_vascular`; transport between the two states is driven by the unbound concentration on each side via the fixed `fu_plasma` and `fu_brain` fractions.
- **Source aliases:** none.
- **Example models:** `Johnson_2011_olanzapine_rat.R`.

---

## Friberg myelosuppression chains

### circ (**canonical circulating-cell compartment**)
- **Type:** compartment
- **Role:** Friberg-style myelosuppression circulating-cell compartment (Friberg 2002 paclitaxel and derivatives). Terminal compartment of a `precursor1` ... `precursorN` maturation chain.
- **Source aliases:** none.
- **Example models:** Friberg-style myelosuppression PD extractions.
- **Notes:** Replaces a paper-naming `central` for circulating neutrophils / platelets / lymphocytes when the model is a maturation chain rather than a classical-PK central compartment. Suffix-form `circ_<celltype>` (e.g., `circ_anc`, `circ_plt`) is accepted for paired-output multi-cell-type models via the registered `_anc` / `_plt` / `_wbc` metabolite suffixes.

---

## Urinary excretion

### urine (**canonical urinary-excretion compartment**)
- **Type:** compartment
- **Role:** Urinary-excretion compartment for renally cleared small molecules. Tracks the cumulative excreted amount of the parent drug.
- **Source aliases:** none.
- **Example models:** renally cleared small-molecule popPK extractions.
- **Notes:** Per-metabolite `urine_<metab>` is accepted via the registered metabolite suffixes (e.g., `urine_apap`, `urine_morphine`). For PBPK extractions, the `a_urine` form on the PBPK organ-amount namespace is used instead.

---

## mPBPK exception (Cao 2013 family)

The Cao 2013 mAb mPBPK family uses paper-anatomical compartment names that are an explicit exception to the standard `central` / `peripheral1` / `peripheral2` convention. The physiological meaning of plasma / tight / leaky / lymph is load-bearing and would be lost under the generic `peripheralN` renaming. Codified 2026-05-28 per the naming audit.

### plasma (**canonical mPBPK plasma compartment**)
- **Type:** compartment
- **Role:** Plasma compartment in the Cao 2013 mAb mPBPK family.
- **Source aliases:** none.
- **Example models:** `Cao_2013_mAb_mPBPK.R` (and Cao_2013_* variants), `Yuan_2019_concizumab.R`.

### tight (**canonical mPBPK tight-tissue compartment**)
- **Type:** compartment
- **Role:** Tight-junction tissue compartment in the Cao 2013 mAb mPBPK family.
- **Source aliases:**
  - `sigma1` -- deprecated paper-mechanistic name (per the 2026-05-28 naming audit R9 rename).
- **Example models:** `Cao_2013_mAb_mPBPK.R`, `Yuan_2019_concizumab.R`.

### leaky (**canonical mPBPK leaky-tissue compartment**)
- **Type:** compartment
- **Role:** Leaky-tissue compartment in the Cao 2013 mAb mPBPK family.
- **Source aliases:**
  - `sigma2` -- deprecated paper-mechanistic name (per the 2026-05-28 naming audit R9 rename).
- **Example models:** `Cao_2013_mAb_mPBPK.R`, `Yuan_2019_concizumab.R`.

### lymph (**canonical mPBPK lymph compartment**)
- **Type:** compartment
- **Role:** Lymph compartment in the Cao 2013 mAb mPBPK family.
- **Source aliases:** none.
- **Example models:** `Cao_2013_mAb_mPBPK.R`, `Yuan_2019_concizumab.R`.

---

## Enterohepatic recirculation

### gallbladder (**canonical gallbladder compartment**)
- **Type:** compartment
- **Role:** Gallbladder / biliary recirculation compartment in enterohepatic-circulation (EHC) popPK models. Drug accumulates from central via biliary excretion (`k12`) and re-enters central after a delay (`k21` gated by gallbladder-emptying time `tg`), producing the characteristic second-peak phenomenon.
- **Source aliases:** none.
- **Example models:** `Ide_2009_pravastatin.R`.

---

## Soluble-receptor biomarkers

### svegfr2 (**canonical soluble VEGFR-2 plasma compartment**)
- **Type:** compartment
- **Role:** Soluble vascular endothelial growth factor receptor 2 plasma compartment used by indirect-response biomarker PD models for angiogenesis inhibitors.
- **Source aliases:** none.
- **Example models:** `Ait-Oudhia_2016_sunitinib.R`, `Hansson_2013a_sunitinib.R`.

### svegfr3 (**canonical soluble VEGFR-3 plasma compartment**)
- **Type:** compartment
- **Role:** Soluble VEGFR-3 turnover compartment used alongside `svegfr2` in the Hansson 2013 sunitinib biomarker panel.
- **Source aliases:** none.
- **Example models:** `Hansson_2013a_sunitinib.R`.
- **Notes:** Registered 2026-05-28 per the naming audit.

---

## Tumor growth inhibition (TGI) states

### tumor (**canonical tumor compartment**)
- **Type:** compartment
- **Role:** Tumor / tumour-size compartment in oncology TGI models.
- **Source aliases:** none.
- **Example models:** `Ait-Oudhia_2016_sunitinib.R`, `NA_NA_sunitinib.R`, `Schindler_2016_sunitinib.R`, `Wilbaux_2015_paclitaxel.R`.

### tumor_size (**canonical TGI tumor-size output state**)
- **Type:** compartment
- **Role:** Snake-case canonical output-state name for the TGI template family.
- **Source aliases:**
  - `Ts` -- legacy.
  - `ts` -- legacy.
- **Example models:** `tgi_no_sat_*.R`, `tgi_sat_*.R`, `Ouerdani_2015_pazopanib.R`, `Mazzocco_2015.R`, `Zecchin_2016.R`, `Wilson_2015_sunitinib_irinotecan_mouse.R`.
- **Notes:** Registered 2026-05-28 per the naming audit for the TGI template family.

### carrying_capacity (**canonical TGI saturable growth ceiling**)
- **Type:** compartment
- **Role:** Gompertz / generalised-logistic ceiling `K` used alongside `tumor_size` in saturable-growth TGI variants.
- **Source aliases:** none.
- **Example models:** `tgi_sat_*.R`.
- **Notes:** Registered 2026-05-28 per the naming audit.

### cycling_cells (**canonical proliferating-cell pool**)
- **Type:** compartment
- **Role:** Proliferating-cell pool in the Simeoni 2004 / Wilson 2015 oncology TGI cell-cycle decomposition. Drug-driven killing routes cells through three damaged-cell transit compartments before clearance.
- **Source aliases:** none.
- **Example models:** `Simeoni_2004.R`, `Wilson_2015_sunitinib_irinotecan_mouse.R`.
- **Notes:** Codified 2026-05-28 per the naming audit.

### damaged_cells1 (**canonical first damaged-cell transit compartment**)
- **Type:** compartment
- **Role:** First damaged-cell transit compartment in Simeoni 2004 cell-cycle decomposition.
- **Source aliases:** none.
- **Example models:** `Simeoni_2004.R`, `Wilson_2015_sunitinib_irinotecan_mouse.R`.

### damaged_cells2 (**canonical second damaged-cell transit compartment**)
- **Type:** compartment
- **Role:** Second damaged-cell transit compartment in Simeoni 2004 cell-cycle decomposition.
- **Source aliases:** none.
- **Example models:** `Simeoni_2004.R`, `Wilson_2015_sunitinib_irinotecan_mouse.R`.

### damaged_cells3 (**canonical third damaged-cell transit compartment**)
- **Type:** compartment
- **Role:** Third damaged-cell transit compartment in Simeoni 2004 cell-cycle decomposition.
- **Source aliases:** none.
- **Example models:** `Simeoni_2004.R`, `Wilson_2015_sunitinib_irinotecan_mouse.R`.

---

## Endogenous metabolic species

### glucose (**canonical plasma glucose**)
- **Type:** compartment
- **Role:** Endogenous plasma glucose used by glucose / lactate turnover sub-models with drug-stimulated production. State holds a concentration (mmol/L), mirroring the source paper's mass-balance parameterisation. Also used by integrated glucose-insulin homeostasis models (e.g., Silber 2007 framework, Hong 2013 HGC / MTT models) as the dynamic-state glucose amount or concentration; per-model `units` field documents which.
- **Source aliases:** none.
- **Example models:** `Oualha_2014_epinephrine.R`, `Hong_2013_glucose_insulin_HGC.R`, `Hong_2013_glucose_insulin_MTT.R`.

### insulin (**canonical plasma insulin compartment**)
- **Type:** compartment
- **Role:** Endogenous plasma insulin used by integrated glucose-insulin homeostasis models (Silber 2007 framework, Jauslin 2007 OGTT framework, Hong 2013 HGC / MTT models). State holds insulin amount (mU) or concentration (mU/L or pmol/L) consistent with the paper's mass-balance parameterisation; the per-model `units` field documents which. Distinct from the existing `INS` (time-varying plasma-insulin regressor covariate) and `INS_BL` (baseline plasma-insulin covariate) -- those are exogenous inputs that drive other models; `insulin` is the dynamic state when insulin is itself a modelled quantity with its own ODE (production / secretion plus elimination).
- **Source aliases:** none.
- **Example models:** `Hong_2013_glucose_insulin_HGC.R` (founding example; insulin amount mU, dynamic state with Gaussian first-phase plus linear second-phase secretion and first-order CLI/VI elimination), `Hong_2013_glucose_insulin_MTT.R` (insulin amount mU, dynamic state with power-function + Emax-incretin-stimulated secretion and first-order CLI/VI elimination).
- **Notes:** Companion canonical to the existing `glucose` (plasma glucose; Oualha 2014 epinephrine), `igf1` (IGF-1 plasma biomarker; somatropin / GH PK-PD), `prolactin`, `nefa`, and `lactate` plasma-biomarker compartments. Adding `insulin` as canonical (rather than declaring it via `paper_specific_compartments`) reflects its high generalisability -- any integrated glucose-insulin model will need a dynamic-state `insulin` compartment, in the same way that `glucose` is canonical rather than paper-specific. The Bizzotto 2016 glucose-kinetics model uses INS as a regressor (no `insulin` state) and so does not exercise this canonical. `NA_NA_paracetamol.R` uses INS as a regressor too but declares `effect_ins` as a paper-specific effect compartment for insulin's delayed action on glucose elimination -- distinct from the dynamic-state-as-`insulin` use here.

### lactate (**canonical plasma lactate**)
- **Type:** compartment
- **Role:** Endogenous plasma lactate produced at the rate of glucose elimination and first-order eliminated. State holds a concentration (mmol/L).
- **Source aliases:** none.
- **Example models:** `Oualha_2014_epinephrine.R`.

### nefa (**canonical plasma non-esterified fatty acids**)
- **Type:** compartment
- **Role:** Plasma non-esterified fatty acids (NEFA / free fatty acids) used by lipid-turnover PD models with feedback control. State holds a concentration (mmol/L).
- **Source aliases:** none.
- **Example models:** `Ahlstrom_2010_niacin.R`.
- **Notes:** NiAc inhibits hydrolysis of TG to NEFA; NEFA formation is also suppressed by a moderator transit chain (`precursor1` .. `precursor8`) representing insulin-like delayed feedback, with a NiAc-independent capillary release term setting the lower physiological limit.

---

## Purine metabolism (Hill-McManus 2017)

### xanthine (**canonical serum xanthine**)
- **Type:** compartment
- **Role:** Serum xanthine amount (mg) in xanthine-oxidase / uric-acid turnover models.
- **Source aliases:** none.
- **Example models:** `Hill-McManus_2017_uricLT.R` (doi:10.1111/bcp.13427).

### urate (**canonical serum urate**)
- **Type:** compartment
- **Role:** Serum urate amount (mg).
- **Source aliases:** none.
- **Example models:** `Hill-McManus_2017_uricLT.R`.

### xanthine_urine (**canonical urinary xanthine excretion**)
- **Type:** compartment
- **Role:** Cumulative urinary xanthine excretion (mg) integrated from `CLX` renal-clearance outflow for direct comparison with 24-h urinary collection data.
- **Source aliases:** none.
- **Example models:** `Hill-McManus_2017_uricLT.R`.

### urate_urine (**canonical urinary urate excretion**)
- **Type:** compartment
- **Role:** Cumulative urinary urate excretion (mg) integrated from `CLUA` renal-clearance outflow.
- **Source aliases:** none.
- **Example models:** `Hill-McManus_2017_uricLT.R`.

---

## Multistate Tuberculosis Pharmacometric (MTP) bacterial states

The MTP framework partitions the bacterial population into three states. The original Clewe series uses the `*bugs` form; later Svensson / Wicha rifampicin papers use the bare `fast` / `slow` / `nonm` form. Both are canonical under the MTP exception, registered 2026-05-28 per the naming audit.

### fast (**canonical fast-multiplying bacteria**)
- **Type:** compartment
- **Role:** Fast-multiplying bacterial subpopulation in MTP TB models.
- **Source aliases:** none.
- **Example models:** `Svensson_2016_rifampicin.R`, `Wicha_2018_rifampicin.R`.

### slow (**canonical slow-multiplying bacteria**)
- **Type:** compartment
- **Role:** Slow-multiplying bacterial subpopulation in MTP TB models.
- **Source aliases:** none.
- **Example models:** `Svensson_2016_rifampicin.R`, `Wicha_2018_rifampicin.R`.

### nonm (**canonical non-multiplying bacteria**)
- **Type:** compartment
- **Role:** Non-multiplying bacterial subpopulation in MTP TB models.
- **Source aliases:** none.
- **Example models:** `Svensson_2016_rifampicin.R`, `Wicha_2018_rifampicin.R`.

### fbugs (**canonical fast-multiplying bacteria (Clewe form)**)
- **Type:** compartment
- **Role:** Fast-multiplying bacterial subpopulation in the Clewe-series MTP form.
- **Source aliases:** none.
- **Example models:** `Clewe_2016_rifampicin.R`, `Chen_2017_TB_MTP_GPDI_mouse.R`, `Clewe_2018_TB_MTP_GPDI_in_vitro.R`.

### sbugs (**canonical slow-multiplying bacteria (Clewe form)**)
- **Type:** compartment
- **Role:** Slow-multiplying bacterial subpopulation in the Clewe-series MTP form.
- **Source aliases:** none.
- **Example models:** `Clewe_2016_rifampicin.R`, `Chen_2017_TB_MTP_GPDI_mouse.R`, `Clewe_2018_TB_MTP_GPDI_in_vitro.R`.

### nbugs (**canonical non-multiplying bacteria (Clewe form)**)
- **Type:** compartment
- **Role:** Non-multiplying bacterial subpopulation in the Clewe-series MTP form.
- **Source aliases:** none.
- **Example models:** `Clewe_2016_rifampicin.R`, `Chen_2017_TB_MTP_GPDI_mouse.R`, `Clewe_2018_TB_MTP_GPDI_in_vitro.R`.

---

## Enzyme-induction reservoirs

### enzyme (**canonical enzyme pool (Wicha form)**)
- **Type:** compartment
- **Role:** Bare `enzyme` compartment for the autoinduction mass-action term in Wicha 2018 / Svensson 2018 rifampicin autoinduction popPK.
- **Source aliases:** none.
- **Example models:** `Wicha_2018_rifampicin.R`, `Svensson_2018_rifampicin.R`.
- **Notes:** Registered 2026-05-28 per the naming audit.

### enz_pool (**canonical enzyme pool (Clewe form)**)
- **Type:** compartment
- **Role:** Central enzyme pool that drives time-varying clearance via an indirect-response loop in Clewe 2015 / Svensson 2016 rifampicin autoinduction.
- **Source aliases:** none.
- **Example models:** `Clewe_2015_rifampicin.R`, `Svensson_2016_rifampicin.R`.

### enzyme_2b6 (**canonical CYP2B6 enzyme pool**)
- **Type:** compartment
- **Role:** Parallel CYP2B6 enzyme pool for autoinduction popPK models in which a drug induces CYP2B6 alongside other isoenzymes; the relative enzyme amount drives time-varying CL of the CYP2B6-mediated arm. Extension of the single-`enzyme` Wicha form to a CYP-isoenzyme-resolved pair.
- **Source aliases:** none.
- **Example models:** `Heathman_2024_efavirenz.R`.
- **Notes:** Initial condition `enzyme_2b6(0) <- 1` (relative to baseline). Founding example: `Heathman_2024_efavirenz.R` (efavirenz autoinduction of both CYP2B6 and CYP2A6).

### enzyme_2a6 (**canonical CYP2A6 enzyme pool**)
- **Type:** compartment
- **Role:** Parallel CYP2A6 enzyme pool for autoinduction popPK models in which a drug induces CYP2A6 alongside other isoenzymes; the relative enzyme amount drives time-varying CL of the CYP2A6-mediated arm.
- **Source aliases:** none.
- **Example models:** `Heathman_2024_efavirenz.R`.
- **Notes:** Initial condition `enzyme_2a6(0) <- 1` (relative to baseline). Founding example: `Heathman_2024_efavirenz.R`.

---

## DAS28 disease-activity score

### das28 (**canonical DAS28 output compartment**)
- **Type:** compartment
- **Role:** DAS28 disease-activity score output compartment used by rheumatoid-arthritis PD models. Single PD output, paper-named.
- **Source aliases:** none.
- **Example models:** `Frey_2013_tocilizumab.R`, `Ma_2020_sarilumab_das28crp.R`.
- **Notes:** Registered 2026-05-28 per the naming audit.

---

## Body-weight PD output

### bw (**canonical body-weight PD output**)
- **Type:** compartment
- **Role:** Body-weight PD output compartment used by drug-induced weight-change models. The state is the kg / g body-weight value with first-order turnover driven by drug-modulated production.
- **Source aliases:** none.
- **Example models:** `Han_2015_sibutramine.R`, `Thorsted_2016_somatropin_rat.R`.
- **Notes:** Registered 2026-05-28 per the naming audit.

### igf1 (**canonical IGF-1 plasma biomarker**)
- **Type:** compartment
- **Role:** IGF-1 (insulin-like growth factor 1) plasma biomarker compartment used by somatropin / GH PK/PD models. Stimulated by central GH via an Emax function; drives downstream body-weight dynamics.
- **Source aliases:** none.
- **Example models:** `Thorsted_2016_somatropin_rat.R`, `Thorsted_2016_somatropin_human.R`.
- **Notes:** Registered 2026-05-28 per the naming audit.

### BMD (**canonical bone mineral density**)
- **Type:** compartment
- **Role:** Bone mineral density PD output, reported in g/cm^2 (femoral neck / lumbar spine / total hip per the source paper's DXA region). Used as the modelled endpoint in osteoporosis disease-progression models (e.g., piecewise-linear menopausal BMD trajectories driven by years since final menstrual period; semi-mechanistic bone-remodelling models linking RANKL / OPG / PTH / vitamin D dynamics to BMD turnover). Usually algebraic in `model()`; can be an ODE state in mechanistic bone-remodelling models.
- **Source aliases:** `bmd`. Capitalisation matches the convention for upper-case observation variables (`Cc`).
- **Example models:** `Plan_2012_bmd_fracture.R`.
- **Notes:** Founding example Plan 2012 PAGE poster; the BMD output is algebraic (BMD <- b * (1 + sum(slope * piece) / 100)) and is the LHS of the residual-error declaration `BMD ~ add(addSd)`. Distinct from the `vp_bone` PBPK bone-vascular compartment, which is a drug-distribution state, not an output endpoint.

---

## Gastric / GI transit

### stomach (**canonical stomach compartment**)
- **Type:** compartment
- **Role:** Gastric / stomach compartment used by gastric-emptying transit models where the gastric mass-balance is resolved as a distinct state ahead of the duodenal absorption depot.
- **Source aliases:** none.
- **Example models:** `Guiastrennec_2016_paracetamol.R`, `Back_2018_fenofibrate.R`.
- **Notes:** Registered 2026-05-28 per the naming audit.

### duodenum (**canonical duodenum compartment**)
- **Type:** compartment
- **Role:** Duodenum compartment in GI-segment paracetamol models. Canonical small-intestine subsegment downstream of the stomach.
- **Source aliases:** none.
- **Example models:** `NA_NA_paracetamol.R`, `Allegaert_2015_paracetamol.R`.

### jejunum (**canonical jejunum compartment**)
- **Type:** compartment
- **Role:** Jejunum compartment in GI-segment paracetamol models. Canonical small-intestine subsegment downstream of the duodenum.
- **Source aliases:** none.
- **Example models:** `NA_NA_paracetamol.R`, `Allegaert_2015_paracetamol.R`.

### ileum (**canonical ileum compartment**)
- **Type:** compartment
- **Role:** Ileum compartment in GI-segment paracetamol models. Canonical small-intestine subsegment downstream of the jejunum.
- **Source aliases:** none.
- **Example models:** `NA_NA_paracetamol.R`, `Allegaert_2015_paracetamol.R`.

---

## PBPK organ-amount compartments (a_<organ> namespace)

PBPK organ-amount compartments used by mass-balance whole-body PBPK extractions. Each state holds the drug amount in the named organ. Spelled-out organ vocabulary per the operator clarification on the 2026-05-28 naming audit; abbreviated organ names (e.g., `a_li`, `a_ki`) are deprecated.

### a_liver (**canonical PBPK liver-amount compartment**)
- **Type:** compartment
- **Role:** Liver organ compartment in mass-balance PBPK extractions. Mass balance: `Q_L * (CA - CVL)`.
- **Source aliases:** none.
- **Example models:** `Zurlinden_2016_paracetamol.R`.

### a_hepatic (**canonical PBPK hepatic metabolite pool**)
- **Type:** compartment
- **Role:** Hepatic intermediate metabolite pool. Distinct from `a_liver`; Zurlinden treats hepatic metabolite formation as a separate compartment from the liver organ for the APAP-sulfate and APAP-glucuronide species.
- **Source aliases:** none.
- **Example models:** `Zurlinden_2016_paracetamol.R`.

### a_fat (**canonical PBPK adipose-amount compartment**)
- **Type:** compartment
- **Role:** Adipose tissue compartment in mass-balance PBPK extractions.
- **Source aliases:** none.
- **Example models:** `Zurlinden_2016_paracetamol.R`.

### a_muscle (**canonical PBPK muscle-amount compartment**)
- **Type:** compartment
- **Role:** Muscle tissue compartment in mass-balance PBPK extractions.
- **Source aliases:** none.
- **Example models:** `Zurlinden_2016_paracetamol.R`.

### a_kidney (**canonical PBPK kidney-amount compartment**)
- **Type:** compartment
- **Role:** Kidney organ compartment in mass-balance PBPK extractions.
- **Source aliases:** none.
- **Example models:** `Zurlinden_2016_paracetamol.R`.

### a_rapidly_perfused (**canonical PBPK rapidly-perfused tissue compartment**)
- **Type:** compartment
- **Role:** Lumped rapidly perfused tissues compartment for highly perfused organs not modelled individually.
- **Source aliases:** none.
- **Example models:** `Zurlinden_2016_paracetamol.R`.

### a_slowly_perfused (**canonical PBPK slowly-perfused tissue compartment**)
- **Type:** compartment
- **Role:** Lumped slowly perfused tissues compartment.
- **Source aliases:** none.
- **Example models:** `Zurlinden_2016_paracetamol.R`.

### a_venous (**canonical PBPK venous-blood compartment**)
- **Type:** compartment
- **Role:** Venous blood compartment in mass-balance PBPK extractions.
- **Source aliases:** none.
- **Example models:** `Zurlinden_2016_paracetamol.R`.

### a_arterial (**canonical PBPK arterial-blood compartment**)
- **Type:** compartment
- **Role:** Arterial blood compartment in mass-balance PBPK extractions.
- **Source aliases:** none.
- **Example models:** `Zurlinden_2016_paracetamol.R`.

### a_urine (**canonical PBPK urinary excretion compartment**)
- **Type:** compartment
- **Role:** Urinary excretion compartment on the PBPK organ-amount namespace. Distinct from the bare `urine` compartment; PBPK models track urine on the `a_<organ>` namespace alongside the other amount-tracking compartments.
- **Source aliases:** none.
- **Example models:** `Zurlinden_2016_paracetamol.R`.
- **Notes:** Non-PBPK renal-clearance models use the bare `urine` form. Both are acceptable for their respective conventions.

### a_gut (**canonical PBPK gut absorption compartment**)
- **Type:** compartment
- **Role:** Gut absorption / intestinal compartment in mass-balance PBPK extractions.
- **Source aliases:** none.
- **Example models:** `Zurlinden_2016_paracetamol.R`.

### a_fast (**canonical bimodal-disease fast-progression arm**)
- **Type:** compartment
- **Role:** Bimodal disease-progression state for the fast-progression arm in Delor 2013 Alzheimer mixture-of-progression-rates PD model. Per-subject mixture weight selects between this fast arm and the slow arm.
- **Source aliases:** none.
- **Example models:** `Delor_2013_alzheimer.R`.
- **Notes:** Distinct from the PBPK perfusion compartment `a_rapidly_perfused` above (different mechanistic role).

### a_slow (**canonical bimodal-disease slow-progression arm**)
- **Type:** compartment
- **Role:** Bimodal disease-progression state for the slow-progression arm. Paired with `a_fast`.
- **Source aliases:** none.
- **Example models:** `Delor_2013_alzheimer.R`.

---

## PBPK vascular-concentration compartments (vp_<organ> namespace)

PBPK organ-vascular concentration compartments used by membrane-limited PBPK extractions where each organ vascular volume is a distinct state. Spelled-out organ vocabulary per the 2026-05-28 naming audit; abbreviated forms (`vp_li` / `vp_lu` / `vp_ki` / `vp_sp` / `vp_he` / `vp_ht` / `vp_mu` / `vp_sk` / `vp_ad` / `vp_bo` / `vp_br` / `vp_si` / `vp_lr` / `vp_pa` / `vp_th` / `vp_po` / `vp_re` / `vp_ot`) are deprecated.

### vp_liver (**canonical PBPK liver vascular concentration**)
- **Type:** compartment
- **Role:** Vascular concentration in the liver organ in membrane-limited PBPK.
- **Source aliases:**
  - `vp_li` -- deprecated abbreviated form.
- **Example models:** `Parhiz_2024_mRNA_LNP.R`, `Shah_2012_mAb_PBPK.R`.

### vp_lung (**canonical PBPK lung vascular concentration**)
- **Type:** compartment
- **Role:** Vascular concentration in the lung organ.
- **Source aliases:**
  - `vp_lu` -- deprecated.
- **Example models:** `Parhiz_2024_mRNA_LNP.R`, `Shah_2012_mAb_PBPK.R`.

### vp_kidney (**canonical PBPK kidney vascular concentration**)
- **Type:** compartment
- **Role:** Vascular concentration in the kidney organ.
- **Source aliases:**
  - `vp_ki` -- deprecated.
- **Example models:** `Parhiz_2024_mRNA_LNP.R`, `Shah_2012_mAb_PBPK.R`.

### vp_spleen (**canonical PBPK spleen vascular concentration**)
- **Type:** compartment
- **Role:** Vascular concentration in the spleen organ.
- **Source aliases:**
  - `vp_sp` -- deprecated.
- **Example models:** `Parhiz_2024_mRNA_LNP.R`, `Shah_2012_mAb_PBPK.R`.

### vp_heart (**canonical PBPK heart vascular concentration**)
- **Type:** compartment
- **Role:** Vascular concentration in the heart organ.
- **Source aliases:**
  - `vp_he` -- deprecated.
  - `vp_ht` -- deprecated.
- **Example models:** `Parhiz_2024_mRNA_LNP.R`, `Shah_2012_mAb_PBPK.R`.

### vp_muscle (**canonical PBPK muscle vascular concentration**)
- **Type:** compartment
- **Role:** Vascular concentration in muscle tissue.
- **Source aliases:**
  - `vp_mu` -- deprecated.
- **Example models:** `Parhiz_2024_mRNA_LNP.R`, `Shah_2012_mAb_PBPK.R`.

### vp_skin (**canonical PBPK skin vascular concentration**)
- **Type:** compartment
- **Role:** Vascular concentration in skin tissue.
- **Source aliases:**
  - `vp_sk` -- deprecated.
- **Example models:** `Parhiz_2024_mRNA_LNP.R`, `Shah_2012_mAb_PBPK.R`.

### vp_adipose (**canonical PBPK adipose vascular concentration**)
- **Type:** compartment
- **Role:** Vascular concentration in adipose tissue.
- **Source aliases:**
  - `vp_ad` -- deprecated.
- **Example models:** `Parhiz_2024_mRNA_LNP.R`, `Shah_2012_mAb_PBPK.R`.

### vp_bone (**canonical PBPK bone vascular concentration**)
- **Type:** compartment
- **Role:** Vascular concentration in bone tissue.
- **Source aliases:**
  - `vp_bo` -- deprecated.
- **Example models:** `Parhiz_2024_mRNA_LNP.R`, `Shah_2012_mAb_PBPK.R`.

### vp_brain (**canonical PBPK brain vascular concentration**)
- **Type:** compartment
- **Role:** Vascular concentration in the brain organ.
- **Source aliases:**
  - `vp_br` -- deprecated.
- **Example models:** `Parhiz_2024_mRNA_LNP.R`, `Shah_2012_mAb_PBPK.R`.

### vp_small_intestine (**canonical PBPK small-intestine vascular concentration**)
- **Type:** compartment
- **Role:** Vascular concentration in the small intestine.
- **Source aliases:**
  - `vp_si` -- deprecated.
- **Example models:** `Parhiz_2024_mRNA_LNP.R`, `Shah_2012_mAb_PBPK.R`.

### vp_large_intestine (**canonical PBPK large-intestine vascular concentration**)
- **Type:** compartment
- **Role:** Vascular concentration in the large intestine.
- **Source aliases:**
  - `vp_lr` -- deprecated.
- **Example models:** `Parhiz_2024_mRNA_LNP.R`, `Shah_2012_mAb_PBPK.R`.

### vp_pancreas (**canonical PBPK pancreas vascular concentration**)
- **Type:** compartment
- **Role:** Vascular concentration in the pancreas.
- **Source aliases:**
  - `vp_pa` -- deprecated.
- **Example models:** `Parhiz_2024_mRNA_LNP.R`, `Shah_2012_mAb_PBPK.R`.

### vp_thymus (**canonical PBPK thymus vascular concentration**)
- **Type:** compartment
- **Role:** Vascular concentration in the thymus.
- **Source aliases:**
  - `vp_th` -- deprecated.
- **Example models:** `Parhiz_2024_mRNA_LNP.R`, `Shah_2012_mAb_PBPK.R`.

### vp_portal (**canonical PBPK portal vascular concentration**)
- **Type:** compartment
- **Role:** Vascular concentration in portal circulation.
- **Source aliases:**
  - `vp_po` -- deprecated.
- **Example models:** `Parhiz_2024_mRNA_LNP.R`, `Shah_2012_mAb_PBPK.R`.

### vp_remainder (**canonical PBPK remainder vascular concentration**)
- **Type:** compartment
- **Role:** Vascular concentration in the lumped remainder compartment.
- **Source aliases:**
  - `vp_re` -- deprecated.
- **Example models:** `Parhiz_2024_mRNA_LNP.R`, `Shah_2012_mAb_PBPK.R`.

### vp_other (**canonical PBPK other vascular concentration**)
- **Type:** compartment
- **Role:** Vascular concentration in the "other" lumped compartment.
- **Source aliases:**
  - `vp_ot` -- deprecated.
- **Example models:** `Parhiz_2024_mRNA_LNP.R`, `Shah_2012_mAb_PBPK.R`.

---

## Whole-body blood / lymph compartments

### blood (**canonical whole-body blood compartment**)
- **Type:** compartment
- **Role:** Whole-body central blood compartment in membrane-limited PBPK extractions.
- **Source aliases:** none.
- **Example models:** `Parhiz_2024_mRNA_LNP.R`.

### bldeg (**canonical blood-pool LNP degradation reservoir**)
- **Type:** compartment
- **Role:** Blood-pool LNP degradation reservoir in mRNA-LNP PBPK.
- **Source aliases:** none.
- **Example models:** `Parhiz_2024_mRNA_LNP.R`.

### bcc (**canonical central blood cells compartment**)
- **Type:** compartment
- **Role:** Central blood cells compartment in mAb PBPK models.
- **Source aliases:** none.
- **Example models:** `Shah_2012_mAb_PBPK.R`.

### lnode (**canonical lymph-node return compartment**)
- **Type:** compartment
- **Role:** Lymph-node return compartment in mAb PBPK models.
- **Source aliases:** none.
- **Example models:** `Shah_2012_mAb_PBPK.R`.

---

## Standard clinical PD-output biomarkers

These are internationally standardised clinical abbreviations registered as canonical compartments / output-state names so single-output PD models using them pass the relaxed `Cc` rule. Registered 2026-05-28 per the naming audit (operator decision: spell out paper-mechanistic names but the standard clinical abbreviations are themselves canonical and need not be expanded).

### ANC (**canonical absolute neutrophil count**)
- **Type:** compartment
- **Role:** Absolute neutrophil count PD output.
- **Source aliases:** none.
- **Example models:** myelosuppression PD models with ANC output.

### PLT (**canonical platelet count**)
- **Type:** compartment
- **Role:** Platelet count PD output.
- **Source aliases:** none.
- **Example models:** thrombocytopenia PD models.

### WBC (**canonical white blood cell count**)
- **Type:** compartment
- **Role:** White blood cell count PD output.
- **Source aliases:** none.
- **Example models:** leukopenia PD models.

### RBC (**canonical red blood cell count**)
- **Type:** compartment
- **Role:** Red blood cell count PD output.
- **Source aliases:** none.
- **Example models:** anemia PD models.

### INR (**canonical international normalised ratio**)
- **Type:** compartment
- **Role:** Coagulation test (international normalised ratio) PD output.
- **Source aliases:** none.
- **Example models:** warfarin / anticoagulant PD models.

### PT (**canonical prothrombin time**)
- **Type:** compartment
- **Role:** Prothrombin time coagulation-test PD output.
- **Source aliases:** none.
- **Example models:** anticoagulant PD models.

### aPTT (**canonical activated partial thromboplastin time**)
- **Type:** compartment
- **Role:** Activated partial thromboplastin time coagulation-test PD output.
- **Source aliases:** none.
- **Example models:** anticoagulant PD models.

### ECT (**canonical ecarin clotting time**)
- **Type:** compartment
- **Role:** Ecarin clotting time coagulation-test PD output (selective ecarin-activated thrombin assay used to monitor direct thrombin inhibitors).
- **Source aliases:** none.
- **Example models:** direct thrombin inhibitor PD models (dabigatran, argatroban, hirudin / hirulog, melagatran).

### hb (**canonical hemoglobin**)
- **Type:** compartment
- **Role:** Hemoglobin PD output.
- **Source aliases:** none.
- **Example models:** anemia / EPO PD models.

### Hba1c (**canonical glycated hemoglobin (HbA1c) PD output**)
- **Type:** compartment
- **Role:** Glycated hemoglobin (HbA1c, %) PD output state computed as the glycated-fraction percent of the total RBC pool in IGRH-style integrated glucose-RBC-HbA1c models. Counterpart to the existing `HBA1C` covariate canonical in `inst/references/covariate-columns.md`; the covariate represents the baseline / observed clinical measurement, and `Hba1c` is the modelled dynamic output state.
- **Source aliases:**
  - `HBA` -- Lledo-Garcia 2013 / Kjellsson 2015 / Bosch 2025 supplement S2 NONMEM `$ERROR HBA = GLY/TOT*100`.
- **Example models:** `Bosch_2025_glp1ra_hba1c.R` (24-state IGRH transit-chain HbA1c sub-model; founding example).
- **Notes:** Camel-case `Hba1c` matches the convention used by other mixed-case clinical-abbreviation PD-output canonicals in this register (`P1NP`, `QTc`, `aPTT`). Single PD output for IGRH-based HbA1c models.

### PSA (**canonical prostate-specific antigen**)
- **Type:** compartment
- **Role:** Prostate-specific antigen PD output (oncology).
- **Source aliases:** none.
- **Example models:** prostate-cancer PD models.

### P1NP (**canonical procollagen type I N-terminal propeptide**)
- **Type:** compartment
- **Role:** Procollagen type I N-terminal propeptide bone-formation biomarker.
- **Source aliases:** none.
- **Example models:** osteoporosis PD models.

### OC (**canonical osteocalcin**)
- **Type:** compartment
- **Role:** Osteocalcin bone-turnover biomarker.
- **Source aliases:** none.
- **Example models:** osteoporosis PD models.

### TT (**canonical total testosterone / thrombin time**)
- **Type:** compartment
- **Role:** Total testosterone (endocrinology) or thrombin time (coagulation) — paper-dependent; both share the TT abbreviation in the contexts where it appears.
- **Source aliases:** none.
- **Example models:** endocrinology / coagulation PD models.

### QTc (**canonical heart-rate-corrected QT interval**)
- **Type:** compartment
- **Role:** Heart-rate-corrected QT interval (electrocardiographic PD endpoint), typically expressed in ms. Used as the observation variable in direct-effect Emax models of drug-induced QT prolongation (cardiac-safety / thorough-QT studies, e.g. quinidine, moxifloxacin, sotalol). The choice of correction formula (Bazett, Fridericia, individual-correction) is paper-dependent and recorded in `description`; the canonical name covers the corrected-interval output regardless of which correction was applied.
- **Source aliases:** `QTcB` (Bazett), `QTcF` (Fridericia), `QTcI` (individual correction) — translate to `QTc` and record the correction in the model file's description / vignette.
- **Example models:** `Shin_2006_quinidine_QT.R` (Bazett-corrected QT interval; founding example).

### serumK (**canonical serum potassium**)
- **Type:** compartment
- **Role:** Serum potassium concentration PD output / turnover-state, in mmol/L. Used as the observation variable in indirect-response / turnover models of drug-induced potassium shifts (mineralocorticoid-receptor antagonists, potassium-sparing diuretics, RAAS inhibitors). Standard clinical-laboratory biomarker reported on essentially every comprehensive metabolic panel; KDIGO and ESC thresholds for hyperkalemia are at 5.5 and 6.0 mmol/L. Distinct from any drug-PK central compartment because the modelled species is the endogenous electrolyte rather than the dosed drug.
- **Source aliases:** `K`, `K+`, `serum_K`, `POTAS` -- translate to `serumK` when assembling input data; document the source-paper symbol in the model file's `description`.
- **Example models:** `Goulooze_2022_finerenone.R` (FIDELIO-DKD Phase III PKPD turnover model for finerenone effect on serum K; founding example).

---

## Bacterial-count PD outputs

### cfu (**canonical colony-forming-unit count**)
- **Type:** compartment
- **Role:** Linear colony-forming-unit count PD output. In Clewe 2016 rifampicin: `cfu = fbugs + sbugs` with proportional residual error.
- **Source aliases:** none.
- **Example models:** `Clewe_2016_rifampicin.R`.

### log_cfu (**canonical log CFU PD output**)
- **Type:** compartment
- **Role:** Log-transformed sputum / culture CFU output. Universal TB-PK/PD endpoint.
- **Source aliases:** none.
- **Example models:** `Clewe_2018_TB_MTP_GPDI_in_vitro.R`, `Khan_2015_rifampicin.R`, `Mohamed_2016_colistin_meropenem.R`, `Sadouki_2025_meropenem.R`, `Svensson_2016_rifampicin.R`, `Wicha_2018_rifampicin.R`.
- **Notes:** Transform base (ln vs log10) is paper-dependent and documented in each source file. Registered 2026-05-28 per the naming audit.

### MBL (**canonical mean bacterial load**)
- **Type:** compartment
- **Role:** Mean bacterial load PD output (Svensson 2017 bedaquiline). Capitalised form; lowercase `mbl` is a registered bare alias.
- **Source aliases:** none.
- **Example models:** `Svensson_2017_bedaquiline.R`.

---

## Paper-specific PD-endpoint output states

Each entry below is a paper-mechanistic PD endpoint registered as a canonical compartment / output-state name so single-output PD models that use them pass the relaxed `Cc` rule. Registered 2026-05-28 per the naming audit.

### ADAS_cog (**canonical Alzheimer Disease Assessment Scale - cognitive subscale**)
- **Type:** compartment
- **Role:** ADAS-cog cognitive PD endpoint.
- **Source aliases:** none.
- **Example models:** Alzheimer's PD models.

### ADAS_NORM (**canonical ADAS normalised PD output**)
- **Type:** compartment
- **Role:** ADAS normalised (per-paper rescaling) PD endpoint.
- **Source aliases:** none.
- **Example models:** Alzheimer's PD models.

### cdr_mix (**canonical Clinical Dementia Rating mixture output**)
- **Type:** compartment
- **Role:** Clinical Dementia Rating mixture-of-progression-rates output.
- **Source aliases:** none.
- **Example models:** Alzheimer's PD models.

### deltaUPDRS (**canonical change-from-baseline total UPDRS score**)
- **Type:** compartment
- **Role:** Change-from-baseline in total Unified Parkinson's Disease Rating Scale (UPDRS) score PD output, used as the modelled endpoint in algebraic Parkinson's disease-progression models that combine a linear disease-progression slope with an asymptotic short-term symptomatic-effect component.
- **Source aliases:** `Delta UPDRS`, the Lee 2011 paper's `Delta_UPDRS_it` notation.
- **Example models:** `Lee_2011_parkinson_progression.R`.

### tumor_vol (**canonical TGI tumour-volume output state**)
- **Type:** compartment
- **Role:** Tumour volume output state in TGI models.
- **Source aliases:** none.
- **Example models:** `Lobo_2002.R`, `Simeoni_2004.R`.

### aescore (**canonical composite adverse-event score**)
- **Type:** compartment
- **Role:** Composite adverse-event score PD output.
- **Source aliases:** none.
- **Example models:** `Girard_2012_pimasertib.R`.

### bcva (**canonical best-corrected visual acuity**)
- **Type:** compartment
- **Role:** Best-corrected visual acuity PD output (ophthalmology).
- **Source aliases:** none.
- **Example models:** `Mulyukov_2018_ranibizumab.R`.

### score (**canonical generic pain score**)
- **Type:** compartment
- **Role:** Generic pain-score PD output.
- **Source aliases:** none.
- **Example models:** `Plan_2012.R`.

### vas_pred (**canonical visual-analog-scale prediction**)
- **Type:** compartment
- **Role:** Visual-analog-scale prediction PD output (Valitalo 2017 morphine).
- **Source aliases:** none.
- **Example models:** `Valitalo_2017_morphine.R`.

### fatigue_grade (**canonical fatigue grade**)
- **Type:** compartment
- **Role:** Fatigue grade PD output (Hansson 2013c sunitinib).
- **Source aliases:** none.
- **Example models:** `Hansson_2013c_sunitinib.R`.

### walkDist (**canonical 6-minute walk distance**)
- **Type:** compartment
- **Role:** 6-minute walk-test distance PD output (Hamuro 2017 DMD).
- **Source aliases:** none.
- **Example models:** `Hamuro_2017_DMD.R`.

### fev1pp (**canonical FEV1 percent predicted**)
- **Type:** compartment
- **Role:** FEV1 percent predicted pulmonary PD output (cystic fibrosis).
- **Source aliases:** none.
- **Example models:** `Harun_2019_cystic_fibrosis.R`.

### msHeadacheDays (**canonical monthly headache-day count**)
- **Type:** compartment
- **Role:** Monthly headache-day count PD output.
- **Source aliases:** none.
- **Example models:** `FiedlerKelly_2020_fremanezumab.R`.

### migraineDays (**canonical monthly migraine-day count**)
- **Type:** compartment
- **Role:** Monthly migraine-day count PD output.
- **Source aliases:** none.
- **Example models:** `FiedlerKelly_2020_fremanezumab.R`.

### viralLoad (**canonical viral load**)
- **Type:** compartment
- **Role:** Viral load PD output (virology).
- **Source aliases:** none.
- **Example models:** `Koloskoff_2025_ganciclovir.R`.

### prob_roc (**canonical ROC-style logistic PD probability output**)
- **Type:** compartment
- **Role:** Probability output for ROC-style logistic PD models.
- **Source aliases:** none.
- **Example models:** `Shin_2014_sevoflurane.R`.

### prolactin (**canonical serum prolactin output**)
- **Type:** compartment
- **Role:** Serum prolactin PD output (endocrinology).
- **Source aliases:** none.
- **Example models:** `Stevens_2012_remoxipride.R`.

### aaaSize (**canonical abdominal aortic aneurysm size**)
- **Type:** compartment
- **Role:** Abdominal aortic aneurysm size PD output.
- **Source aliases:** none.
- **Example models:** `Sherer_2012_AAA.R`.

### cel_count (**canonical MS lesion cell-count**)
- **Type:** compartment
- **Role:** Cell counts in multiple-sclerosis lesions PD output.
- **Source aliases:** none.
- **Example models:** `VelezdeMendizabal_2013_multiple_sclerosis.R`.

### G (**canonical endogenous glucose output**)
- **Type:** compartment
- **Role:** Endogenous glucose PD output (Bizzotto 2016 glucose).
- **Source aliases:** none.
- **Example models:** `Bizzotto_2016_glucose.R`.

---

## PBPK bare organ-amount compartments (Zhang 2011 family)

PBPK bare organ-amount compartments used by Zhang 2011 nutlin3a and similar full-body PBPK extractions that don't prefix the organ name with `a_` / `vp_`. New PBPK extractions should prefer the spelled-out `a_<organ>` namespace, but the bare forms remain canonical for paper-mechanistic models that already use them. Registered 2026-05-29 per the naming-audit compartment-warning cleanup.

### venous (**canonical bare venous-blood compartment**)
- **Type:** compartment
- **Role:** Bare venous blood compartment in Zhang 2011 PBPK and similar.
- **Source aliases:** none.
- **Example models:** `Zhang_2011_nutlin3a.R`.

### arterial (**canonical bare arterial-blood compartment**)
- **Type:** compartment
- **Role:** Bare arterial blood compartment in Zhang 2011 PBPK and similar.
- **Source aliases:** none.
- **Example models:** `Zhang_2011_nutlin3a.R`.

### adipose (**canonical bare adipose compartment**)
- **Type:** compartment
- **Role:** Bare adipose tissue compartment.
- **Source aliases:** none.
- **Example models:** `Zhang_2011_nutlin3a.R`.

### adrenal (**canonical bare adrenal compartment**)
- **Type:** compartment
- **Role:** Bare adrenal organ compartment.
- **Source aliases:** none.
- **Example models:** `Zhang_2011_nutlin3a.R`.

### bonemarrow (**canonical bare bone-marrow compartment**)
- **Type:** compartment
- **Role:** Bare bone-marrow compartment.
- **Source aliases:** none.
- **Example models:** `Zhang_2011_nutlin3a.R`.

### muscle (**canonical bare muscle compartment**)
- **Type:** compartment
- **Role:** Bare muscle tissue compartment.
- **Source aliases:** none.
- **Example models:** `Zhang_2011_nutlin3a.R`.

### spleen (**canonical bare spleen compartment**)
- **Type:** compartment
- **Role:** Bare spleen organ compartment.
- **Source aliases:** none.
- **Example models:** `Zhang_2011_nutlin3a.R`.

### intestine (**canonical bare intestine compartment**)
- **Type:** compartment
- **Role:** Bare intestine compartment.
- **Source aliases:** none.
- **Example models:** `Zhang_2011_nutlin3a.R`.

### retina (**canonical bare retina compartment**)
- **Type:** compartment
- **Role:** Bare retina compartment used in ocular PBPK extractions.
- **Source aliases:** none.
- **Example models:** ocular-PK models.

### vitreous (**canonical bare vitreous compartment**)
- **Type:** compartment
- **Role:** Bare vitreous-humor compartment used in ocular PBPK extractions.
- **Source aliases:** none.
- **Example models:** ocular-PK models.

### res_vasc (**canonical bare lumped-remainder vascular compartment**)
- **Type:** compartment
- **Role:** Lumped remainder vascular compartment in PBPK.
- **Source aliases:** none.
- **Example models:** `Zhang_2011_nutlin3a.R`.

### res_tis (**canonical bare lumped-remainder tissue compartment**)
- **Type:** compartment
- **Role:** Lumped remainder tissue compartment in PBPK.
- **Source aliases:** none.
- **Example models:** `Zhang_2011_nutlin3a.R`.

### lung (**canonical bare lung compartment**)
- **Type:** compartment
- **Role:** Bare lung organ compartment in full-body PBPK extractions.
- **Source aliases:** none.
- **Example models:** `Zhang_2011_nutlin3a.R`.

### brain (**canonical bare brain compartment**)
- **Type:** compartment
- **Role:** Bare brain organ compartment in full-body PBPK extractions.
- **Source aliases:** none.
- **Example models:** `Zhang_2011_nutlin3a.R`.

### heart (**canonical bare heart compartment**)
- **Type:** compartment
- **Role:** Bare heart organ compartment in full-body PBPK extractions. Total tissue (well-stirred) drug concentration in the myocardium; paired with `lung`, `liver`, `kidney`, `spleen`, `brain` etc. in whole-body PBPK extractions that resolve heart as a distinct organ. The token `heart` already appears in the `vp_heart` canonical entry and in the `pbpkSubCompartmentRegex` valid-organ list (alongside `lung`, `kidney`, `spleen`), so this entry registers the bare-organ form for parity with the surrounding canonicals.
- **Source aliases:** none.
- **Example models:** `An_2012_mitoxantrone_mouse_pbpk.R`, `An_2012_mitoxantrone_human_pbpk.R`.

### bone (**canonical bare bone compartment**)
- **Type:** compartment
- **Role:** Bare bone tissue compartment in full-body PBPK extractions.
- **Source aliases:** none.
- **Example models:** `Zhang_2011_nutlin3a.R`.

### other (**canonical bare "other" lumped compartment**)
- **Type:** compartment
- **Role:** Bare "other" lumped tissue compartment in full-body PBPK extractions.
- **Source aliases:** none.
- **Example models:** `Zhang_2011_nutlin3a.R`.

---

## Adaptive-resistance bacterial states

### aron (**canonical adaptive-resistance on state**)
- **Type:** compartment
- **Role:** Adaptive-resistance "on" state in time-kill TB / antibiotic-combination PD models. Drives the dynamic isoniazid EC50 adaptive-resistance switch.
- **Source aliases:** none.
- **Example models:** `Clewe_2018_TB_MTP_GPDI_in_vitro.R`.
- **Notes:** Drug-suffixed `aron_<drug>` forms are accepted via the registered drug-suffix metabolite tokens.

### aroff (**canonical adaptive-resistance off state**)
- **Type:** compartment
- **Role:** Adaptive-resistance "off" state. Paired with `aron`.
- **Source aliases:** none.
- **Example models:** `Clewe_2018_TB_MTP_GPDI_in_vitro.R`.

---

## Bare drug-state PK compartments (combination antibiotic / TB)

Central-compartment drug states named after the drug INN abbreviation. Sibling to (and accepted alongside) the `central_<drug>` canonical-with-metab-suffix form used by Chen 2017 TB MTP-GPDI mouse and similar.

### lzd (**canonical linezolid bare drug-state compartment**)
- **Type:** compartment
- **Role:** Linezolid bare drug-state PK compartment in combination-antibiotic time-kill PD models.
- **Source aliases:** none.
- **Example models:** `Wicha_2017_linezolid_meropenem_vancomycin.R`.

### mer (**canonical meropenem bare drug-state compartment**)
- **Type:** compartment
- **Role:** Meropenem bare drug-state PK compartment.
- **Source aliases:** none.
- **Example models:** `Wicha_2017_linezolid_meropenem_vancomycin.R`.

### mero (**canonical meropenem long-form bare drug-state**)
- **Type:** compartment
- **Role:** Meropenem (long-form INN) bare drug-state PK compartment.
- **Source aliases:** none.
- **Example models:** `Mohamed_2016_colistin_meropenem.R`.

### van (**canonical vancomycin bare drug-state compartment**)
- **Type:** compartment
- **Role:** Vancomycin bare drug-state PK compartment.
- **Source aliases:** none.
- **Example models:** `Wicha_2017_linezolid_meropenem_vancomycin.R`.

### col (**canonical colistin bare drug-state compartment**)
- **Type:** compartment
- **Role:** Colistin bare drug-state PK compartment.
- **Source aliases:** none.
- **Example models:** `Mohamed_2016_colistin_meropenem.R`.

### dap (**canonical daptomycin bare drug-state compartment**)
- **Type:** compartment
- **Role:** Daptomycin bare drug-state PK compartment.
- **Source aliases:** none.
- **Example models:** combination-antibiotic time-kill PD models.

### cmem (**canonical meropenem bath-concentration state**)
- **Type:** compartment
- **Role:** Meropenem bath / medium concentration state (mg/L) in hollow-fiber infection-model time-kill PD; dosed by the user and declines at the simulated HFIM elimination rate, driving the Hill-type meropenem killing term on each bacterial subpopulation.
- **Source aliases:** none.
- **Example models:** `Rees_2018_meropenem_ciprofloxacin.R`.

### ccip (**canonical ciprofloxacin bath-concentration state**)
- **Type:** compartment
- **Role:** Ciprofloxacin bath / medium concentration state (mg/L) in hollow-fiber infection-model time-kill PD; dosed by the user, drives the Emax ciprofloxacin killing term, and lowers the effective meropenem KC50 through the mechanistic-synergy term. Paired with `cmem`.
- **Source aliases:** none.
- **Example models:** `Rees_2018_meropenem_ciprofloxacin.R`.

---

## Bacterial subpopulation states

Lowercase / uppercase casing variants and mutation suffixes used by combination-antibiotic time-kill models.

### S (**canonical susceptible bacterial subpopulation**)
- **Type:** compartment
- **Role:** Susceptible bacterial population.
- **Source aliases:** none.
- **Example models:** `Mohamed_2016_colistin_meropenem.R`.
- **Notes:** Resistant mutant subpopulation tracked via `S_mut` (mutation suffix permitted by the metabolite-suffix register).

### R (**canonical resistant bacterial subpopulation**)
- **Type:** compartment
- **Role:** Resistant bacterial population. Paired with `S`.
- **Source aliases:** none.
- **Example models:** `Mohamed_2016_colistin_meropenem.R`.

### Fbugs (**canonical fast-multiplying bacteria (capitalised Clewe casing)**)
- **Type:** compartment
- **Role:** Fast-multiplying bacterial population (capitalised Clewe-series casing variant).
- **Source aliases:** none.
- **Example models:** `Clewe_2016_rifampicin.R` (variant casing).

### Sbugs (**canonical slow-multiplying bacteria (capitalised Clewe casing)**)
- **Type:** compartment
- **Role:** Slow-multiplying bacterial population (capitalised casing variant).
- **Source aliases:** none.
- **Example models:** `Clewe_2016_rifampicin.R` (variant casing).

### Nbugs (**canonical non-multiplying bacteria (capitalised Clewe casing)**)
- **Type:** compartment
- **Role:** Non-multiplying bacterial population (capitalised casing variant).
- **Source aliases:** none.
- **Example models:** `Clewe_2016_rifampicin.R` (variant casing).

---

## Bacterial subpopulation states (semi-mechanistic time-kill / HFIM PD)

Semi-mechanistic time-kill / hollow-fiber-infection-model (HFIM) PD models (Bulitta / Wicha / Landersdorfer life-cycle growth lineage) partition the bacterial population by antibiotic-resistance phenotype. The canonical scheme spells the phenotype out so the resistance status of each subpopulation is self-documenting in the state name, rather than using the terse single-letter `s` / `i` / `r` source labels:

- **Single-drug models** name each subpopulation `bact_<phenotype>`, where `<phenotype>` is one of the spelled-out resistance phenotypes `susceptible`, `intermediate`, or `resistant`.
- **Combination-therapy (two-drug) models** name each subpopulation by its joint per-drug status as a spelled-out compound `bact_<drug1pheno>_<drug2pheno>` (the two phenotype tokens are in the model's drug order; e.g., for a meropenem + ciprofloxacin model `bact_resistant_intermediate` is the meropenem-resistant / ciprofloxacin-intermediate subpopulation).
- An **optional trailing digit** indexes the Bulitta / Wicha two-state bacterial life cycle: `1` = vegetative / resting state, `2` = replicating state (the state-2 cells replicate back into state-1 daughter cells). Subpopulations without a life-cycle split carry no trailing digit.

These states are **not** registered as individual H3 entries; they are matched at runtime by the `bacterialSubpopRegex` constant in `R/conventions.R`:

```
^bact_(susceptible|intermediate|resistant)(_(susceptible|intermediate|resistant))?[0-9]*$
```

(The regex is a structural pattern and lives in R alongside the other `*Regex` compartment constants documented in the "Regex constants" header section, not as a name list in this file.)

- **Example models:**
  - `Garonzik_2016_daptomycin.R` — single-drug: `bact_susceptible1` / `bact_susceptible2`, `bact_intermediate1` / `bact_intermediate2`, `bact_resistant1` / `bact_resistant2` (three subpopulations of decreasing daptomycin susceptibility, each with the two-state life cycle).
  - `Rees_2018_meropenem_ciprofloxacin.R` — two-drug (meropenem + ciprofloxacin): `bact_susceptible_susceptible1` / `2`, `bact_resistant_intermediate1` / `2`, `bact_intermediate_resistant1` / `2`.
  - `Landersdorfer_2018_imipenem_tobramycin.R` — two-drug (imipenem + tobramycin): same `bact_<drug1pheno>_<drug2pheno>` compound scheme with the two-state life-cycle digit.

---

## Mann 2022 respiratory / cerebrovascular physiology states

Physiological state variables of the Magosso / Ursino respiratory and cerebrovascular control model with the Mann 2022 opioid-induced-ventilatory-depression and cardiac-arrest extensions (`Mann_2022_respiratory_physiology.R`). The CAR (fraction of mu receptors bound by an opioid agonist) input from the binding layer drives reductions in the ventilatory drives. Registered 2026-05-30 per the naming-warning resolution.

### palv_co2 (**canonical alveolar / arterial CO2 partial pressure**)
- **Type:** compartment
- **Role:** Alveolar / arterial CO2 partial pressure (mmHg) gas-exchange state; balances minute-ventilation washout against mixed-venous CO2 delivery (Magosso / Ursino gas-exchange ODE).
- **Source aliases:** none.
- **Example models:** `Mann_2022_respiratory_physiology.R`.

### palv_o2 (**canonical alveolar / arterial O2 partial pressure**)
- **Type:** compartment
- **Role:** Alveolar / arterial O2 partial pressure (mmHg) gas-exchange state; balances inspired-O2 uptake against mixed-venous O2 delivery and is the quantity compared with the cardiac-arrest threshold.
- **Source aliases:** none.
- **Example models:** `Mann_2022_respiratory_physiology.R`.

### cb_co2 (**canonical brain-blood CO2 content**)
- **Type:** compartment
- **Role:** Brain blood-gas CO2 content state (Spencer dissociation units); driven by cerebral blood flow times arterial-minus-venous CO2 difference plus brain CO2 metabolic production.
- **Source aliases:** none.
- **Example models:** `Mann_2022_respiratory_physiology.R`.

### cb_o2 (**canonical brain-blood O2 content**)
- **Type:** compartment
- **Role:** Brain blood-gas O2 content state (Spencer units); driven by cerebral blood flow times arterial-minus-venous O2 difference plus brain O2 metabolic consumption.
- **Source aliases:** none.
- **Example models:** `Mann_2022_respiratory_physiology.R`.

### ct_co2 (**canonical peripheral-tissue blood CO2 content**)
- **Type:** compartment
- **Role:** Peripheral (non-brain) tissue blood-gas CO2 content state; driven by tissue blood flow times arterial-minus-venous CO2 difference plus tissue CO2 metabolic production.
- **Source aliases:** none.
- **Example models:** `Mann_2022_respiratory_physiology.R`.

### ct_o2 (**canonical peripheral-tissue blood O2 content**)
- **Type:** compartment
- **Role:** Peripheral (non-brain) tissue blood-gas O2 content state; driven by tissue blood flow times arterial-minus-venous O2 difference plus tissue O2 metabolic consumption.
- **Source aliases:** none.
- **Example models:** `Mann_2022_respiratory_physiology.R`.

### yco2 (**canonical filtered peripheral CO2 chemoreflex signal**)
- **Type:** compartment
- **Role:** First-order-filtered peripheral-chemoreflex CO2 input signal (dimensionless); modulates cerebral and peripheral blood flow with time-constant tau_co2.
- **Source aliases:** none.
- **Example models:** `Mann_2022_respiratory_physiology.R`.

### yo2 (**canonical filtered peripheral O2 chemoreflex signal**)
- **Type:** compartment
- **Role:** First-order-filtered peripheral-chemoreflex O2 input signal (dimensionless); modulates blood flow with time-constant tau_o2. Paired with `yco2`.
- **Source aliases:** none.
- **Example models:** `Mann_2022_respiratory_physiology.R`.

### dp_state (**canonical peripheral chemoreflex ventilatory drive**)
- **Type:** compartment
- **Role:** Peripheral-chemoreflex ventilatory drive state (L/min); opioid-attenuated (factor 1 - CAR^P1) and filtered with time-constant tau_Dp.
- **Source aliases:** none.
- **Example models:** `Mann_2022_respiratory_physiology.R`.

### dc_state (**canonical central chemoreflex ventilatory drive**)
- **Type:** compartment
- **Role:** Central-chemoreflex ventilatory drive state (L/min); driven by the brain-CO2-minus-baseline error, opioid-attenuated, and filtered with time-constant tau_Dc. Paired with `dp_state`.
- **Source aliases:** none.
- **Example models:** `Mann_2022_respiratory_physiology.R`.

### alpha_h (**canonical central hypoxic ventilatory-depression factor**)
- **Type:** compartment
- **Role:** Central hypoxic ventilatory-depression factor (dimensionless), Mann 2022 alphaH; first-order relaxes toward a brain-O2-dependent target and multiplies the peripheral drive in the total-ventilation synthesis.
- **Source aliases:** none.
- **Example models:** `Mann_2022_respiratory_physiology.R`.

### t_pao2_below (**canonical sub-threshold-PaO2 dwell-time accumulator**)
- **Type:** compartment
- **Role:** Accumulator (min) of time arterial O2 partial pressure has spent below the critical cardiac-arrest threshold; grows while PaO2 is below threshold and slowly re-arms above it, gating the cardiovascular-collapse trigger.
- **Source aliases:** none.
- **Example models:** `Mann_2022_respiratory_physiology.R`.

### im_arrest (**canonical cardiovascular-collapse cardiac-output multiplier**)
- **Type:** compartment
- **Role:** Cardiovascular-collapse multiplier (dimensionless, 0..1) on cardiac output; once the sustained-hypoxemia trigger fires it decays toward zero, driving total cardiac output toward the Mann 2022 cardiac-arrest floor.
- **Source aliases:** none.
- **Example models:** `Mann_2022_respiratory_physiology.R`.

---

## Mann 2022 mu-opioid receptor occupancy states

Receptor-occupancy state variables of the Mann 2022 competitive mu-opioid receptor binding layer (`Mann_2022_mu_receptor_binding.R`), tracking simultaneous agonist and antagonist occupancy of a shared receptor pool. Registered 2026-05-30 per the naming-warning resolution.

### RL_op (**canonical mu-opioid agonist receptor-occupancy fraction**)
- **Type:** compartment
- **Role:** Fraction of the mu-opioid receptor pool bound by the opioid agonist (0..1); follows the multi-ligand competitive binding ODE and is the CAR output piped into the respiratory-physiology layer.
- **Source aliases:** none.
- **Example models:** `Mann_2022_mu_receptor_binding.R`.

### RL_antag (**canonical mu-opioid antagonist receptor-occupancy fraction**)
- **Type:** compartment
- **Role:** Fraction of the mu-opioid receptor pool bound by the opioid antagonist (0..1), competing with `RL_op` for the shared free-receptor fraction `R_free = 1 - RL_op - RL_antag`.
- **Source aliases:** none.
- **Example models:** `Mann_2022_mu_receptor_binding.R`.

---

## Inflammatory-mediator PD states

Indirect-response state variables of the Xiang 2018 baicalein anti-inflammatory cellular PD cascade (`Xiang_2018_baicalein.R`): a TNF-alpha -> {IL-6, iNOS -> NO} indirect-response chain in LPS-stimulated RAW264.7 macrophages. Registered 2026-05-30 per the naming-warning resolution.

### tnf (**canonical TNF-alpha indirect-response state**)
- **Type:** compartment
- **Role:** TNF-alpha indirect-response state (pg/mL) with LPS-stimulated zero-order production inhibited log-linearly by baicalein, and first-order elimination; drives the downstream IL-6 and iNOS responses via delay states.
- **Source aliases:** none.
- **Example models:** `Xiang_2018_baicalein.R`.

### il6 (**canonical IL-6 indirect-response state**)
- **Type:** compartment
- **Role:** IL-6 indirect-response state (pg/mL) produced at a rate proportional to the lag-delayed TNF-alpha signal, with first-order elimination.
- **Source aliases:** none.
- **Example models:** `Xiang_2018_baicalein.R`.
- **Notes:** Distinct namespace from the `IL6` covariate (an upstream interleukin-6 covariate column); the lowercase `il6` compartment is the modelled IL-6 PD state, not a covariate.

### inos (**canonical iNOS-expression indirect-response state**)
- **Type:** compartment
- **Role:** Inducible nitric-oxide-synthase (iNOS) expression state (relative to the t = 0 control) produced from the lag-delayed TNF-alpha signal; elimination held at zero per source to match the post-12.5 h plateau.
- **Source aliases:** none.
- **Example models:** `Xiang_2018_baicalein.R`.

### no (**canonical nitric-oxide indirect-response state**)
- **Type:** compartment
- **Role:** Nitric-oxide state (uM) produced from iNOS via an iNOS^delta amplification term, with elimination held at zero per source.
- **Source aliases:** none.
- **Example models:** `Xiang_2018_baicalein.R`.

---

## Radiation tumor-growth-inhibition states

State variables specific to the Cardilin 2018 combination radiation + radiosensitizer tumor-growth-inhibition model (`Cardilin_2018_radiation_radiosensitizer_mouse.R`), where linear-quadratic radiation kill routes proliferating cells into an irradiated-cell chain that divides at most once more before dying. Registered 2026-05-30 per the naming-warning resolution.

### irrad1 (**canonical first irradiated-cell pool**)
- **Type:** compartment
- **Role:** First irradiated-cell pool; receives proliferating cells killed by the linear-quadratic radiation hazard at each fraction and either dies or progresses to a final division.
- **Source aliases:** none.
- **Example models:** `Cardilin_2018_radiation_radiosensitizer_mouse.R`.

### irrad2 (**canonical second irradiated-cell pool**)
- **Type:** compartment
- **Role:** Second irradiated-cell pool fed by the post-division progression from `irrad1` (factor-2 source for the one final division before death); first-order natural death thereafter.
- **Source aliases:** none.
- **Example models:** `Cardilin_2018_radiation_radiosensitizer_mouse.R`.

### radDepot (**canonical radiation-timing trigger compartment**)
- **Type:** compartment
- **Role:** Radiation-timing trigger compartment; a unit bolus is dosed in at each irradiation time and decays fast so that the kill hazard integrates to the linear-quadratic lethal-lesion number per fraction (a Dirac-delta numerical device, not a fitted state).
- **Source aliases:** none.
- **Example models:** `Cardilin_2018_radiation_radiosensitizer_mouse.R`.

---

## Viral-dynamics states (Neumann target-cell model)

State variables of the Neumann-style three-state HCV target-cell viral-dynamics model used in the Wang 2018 daclatasvir / asunaprevir integrated PK / viral-dynamic model (`Wang_2018_daclatasvir_asunaprevir.R`). Registered 2026-05-30 per the naming-warning resolution.

### virus (**canonical free-virus / virion pool**)
- **Type:** compartment
- **Role:** Free-virus / virion pool (Neumann state V); produced by productively-infected cells (drug-inhibited via the combination antiviral effect) and cleared first-order; the log10 of this state is the viral-load output.
- **Source aliases:** none.
- **Example models:** `Wang_2018_daclatasvir_asunaprevir.R`.

### infected (**canonical productively-infected cells**)
- **Type:** compartment
- **Role:** Productively-infected hepatocytes (Neumann state I); produced by infection of target cells by free virus and lost first-order at rate delta.
- **Source aliases:** none.
- **Example models:** `Wang_2018_daclatasvir_asunaprevir.R`.
- **Notes:** The uninfected target-cell state reuses the existing canonical `target` compartment.

---

## Airway interstitial-fluid (ISF) mAb / target species

Airway interstitial-fluid (ISF) species of the Rymut 2023 mechanistic anti-tryptase mAb (MTPS9579A) PK/PD model (`Rymut_2023_anti_tryptase.R`), where free mAb and mAb-monomer complex are delivered from the systemic circulation to the airway ISF via lymph flow. Registered 2026-05-30 per the naming-warning resolution.

### mab_isf (**canonical free mAb in airway ISF**)
- **Type:** compartment
- **Role:** Free MTPS9579A (anti-tryptase mAb) concentration in the airway ISF (nM); enters via lymph influx and binds tetrameric and monomeric tryptase in the ISF.
- **Source aliases:** none.
- **Example models:** `Rymut_2023_anti_tryptase.R`.

### monomer_isf (**canonical free monomeric tryptase in airway ISF**)
- **Type:** compartment
- **Role:** Free inactive monomeric tryptase concentration in the airway ISF (nM); generated by spontaneous tetramer dissociation and by mAb-induced tetramer disruption, eliminated first-order, and bound by free mAb.
- **Source aliases:** none.
- **Example models:** `Rymut_2023_anti_tryptase.R`.

### complex_monomer_isf (**canonical mAb-monomer complex in airway ISF**)
- **Type:** compartment
- **Role:** MTPS9579A-monomer complex concentration in the airway ISF (nM); formed from free mAb plus monomeric tryptase (and from tetramer-complex disruption), and also receiving the systemic mAb-monomer complex via lymph influx.
- **Source aliases:** none.
- **Example models:** `Rymut_2023_anti_tryptase.R`.
- **Notes:** The active-tetramer species `target_isf` and the mAb-tetramer complex `complex_isf` in the same model already pass via `targetLocationRegex`.

---

## Body-composition / disease-risk PD outputs

Population body-composition / disease-risk PD output states from the Oniki 2018 elderly-Japanese health-screening companion models. Registered 2026-05-30 per the naming-warning resolution.

### bmi (**canonical body-mass-index PD output**)
- **Type:** compartment
- **Role:** Body-mass-index PD output (kg/m^2); the typical BMI as a power-of-age scalar with a female sex multiplier and a DsbA-L T/T additive shift, with log-normal between-subject variability. Sibling of `bw` / `weight`.
- **Source aliases:** none.
- **Example models:** `Oniki_2018_bmi.R`.
- **Notes:** Lowercase `bmi` is the modelled BMI PD output; distinct from the uppercase `BMI` covariate column that drives the companion NAFLD-risk model.

### p_nafld (**canonical NAFLD-probability PD output**)
- **Type:** compartment
- **Role:** Probability of nonalcoholic fatty liver disease (NAFLD) PD output (0..1); the expit of a baseline logit floor plus a sigmoidal-Emax function of (BMI - 17) with genotype / lab covariate effects. Sibling of `prob_roc`.
- **Source aliases:** none.
- **Example models:** `Oniki_2018_nafld_risk.R`.

---

## Bare drug-effect mechanistic states

### gro (**canonical growing-bacteria state**)
- **Type:** compartment
- **Role:** Growing bacterial state in combination-PD models.
- **Source aliases:** none.
- **Example models:** `Wicha_2017_linezolid_meropenem_vancomycin.R`.

### repl (**canonical replicating-bacteria state**)
- **Type:** compartment
- **Role:** Replicating bacterial state in combination-PD models.
- **Source aliases:** none.
- **Example models:** `Wicha_2017_linezolid_meropenem_vancomycin.R`.

### pers (**canonical persistent-bacteria state**)
- **Type:** compartment
- **Role:** Persistent (non-replicating) bacterial state in combination-PD models.
- **Source aliases:** none.
- **Example models:** `Wicha_2017_linezolid_meropenem_vancomycin.R`.

---

## Hormonal PD output

### glp1 (**canonical glucagon-like peptide 1 output**)
- **Type:** compartment
- **Role:** Glucagon-like peptide 1 hormone PD output.
- **Source aliases:** none.
- **Example models:** `NA_NA_paracetamol.R`.

### gip (**canonical glucose-dependent insulinotropic polypeptide output**)
- **Type:** compartment
- **Role:** Glucose-dependent insulinotropic polypeptide PD output.
- **Source aliases:** none.
- **Example models:** `NA_NA_paracetamol.R`.

---

## Lab values and endogenous biomarker compartments

Standard clinical-biomarker / endogenous-output compartments. Widely-recognised clinical lab values, endogenous biomarkers, immune-cell populations, and standard organ-anatomy compartments used as PD output states. Registered 2026-05-29 per the naming-audit long-tail compartment cleanup.

### crp (**canonical C-reactive protein PD output**)
- **Type:** compartment
- **Role:** C-reactive protein biomarker PD output.
- **Source aliases:** none.
- **Example models:** `Yang_2016_dilmapimod.R`, `AitOudhia_2012_IL1beta.R`.

### sdma (**canonical symmetric dimethylarginine PD output**)
- **Type:** compartment
- **Role:** Symmetric dimethylarginine biomarker PD output.
- **Source aliases:** none.
- **Example models:** `Guo_2022_PRMT5.R`.

### ldl (**canonical LDL-cholesterol PD output**)
- **Type:** compartment
- **Role:** LDL-cholesterol biomarker PD output.
- **Source aliases:** none.
- **Example models:** `Pu_2021_evinacumab.R`.

### cox2 (**canonical COX-2 enzyme activity PD output**)
- **Type:** compartment
- **Role:** COX-2 enzyme activity PD output.
- **Source aliases:** none.
- **Example models:** `VasquezBahena_2009_lumiracoxib.R`.

### ast (**canonical aspartate aminotransferase PD output**)
- **Type:** compartment
- **Role:** Aspartate aminotransferase liver-function biomarker.
- **Source aliases:** none.
- **Example models:** `Yang_2024_axatilimab.R`.

### cpk (**canonical creatine phosphokinase PD output**)
- **Type:** compartment
- **Role:** Creatine phosphokinase muscle biomarker.
- **Source aliases:** none.
- **Example models:** `Yang_2024_axatilimab.R`.

### csf1 (**canonical colony-stimulating factor 1 PD output**)
- **Type:** compartment
- **Role:** Colony-stimulating factor 1 biomarker PD output.
- **Source aliases:** none.
- **Example models:** `Yang_2024_axatilimab.R`.

### igg (**canonical IgG endogenous turnover compartment**)
- **Type:** compartment
- **Role:** IgG endogenous turnover compartment.
- **Source aliases:** none.
- **Example models:** `Kim_2006_igg_model.R`.

### total_igg (**canonical total serum IgG compartment**)
- **Type:** compartment
- **Role:** Total serum IgG PD output.
- **Source aliases:**
  - `total_IgG` -- prior canonical name (pre-2026-06-19 case standardization).
- **Example models:** `Valenzuela_2025_nipocalimab.R`.
- **Notes:** Renamed from `total_IgG` to `total_igg` on 2026-06-19 per the canonical-register standardization audit (operator decision: compartment names follow the all-lowercase nlmixr2 convention; the mixed-case `total_IgG` was an outlier).

### phe (**canonical phenylalanine PD output**)
- **Type:** compartment
- **Role:** Phenylalanine biomarker PD output (phenylketonuria models).
- **Source aliases:** none.
- **Example models:** `Charbonneau_2021_phenylalanine.R`.

### pth (**canonical parathyroid hormone PD output**)
- **Type:** compartment
- **Role:** Parathyroid hormone biomarker PD output.
- **Source aliases:** none.
- **Example models:** `Ahn_2014.R`.

### ca (**canonical serum calcium PD output**)
- **Type:** compartment
- **Role:** Serum calcium biomarker PD output.
- **Source aliases:** none.
- **Example models:** `Ahn_2014.R`.

### ca_unobs (**canonical unobserved calcium pool**)
- **Type:** compartment
- **Role:** Unobserved calcium pool used in calcium homeostasis models.
- **Source aliases:** none.
- **Example models:** `Ahn_2014.R`.

### thb (**canonical total hemoglobin PD output**)
- **Type:** compartment
- **Role:** Total hemoglobin biomarker PD output (erythropoiesis models).
- **Source aliases:** none.
- **Example models:** `Tetschke_2018_erythropoiesis.R`.

### iron (**canonical serum iron PD output**)
- **Type:** compartment
- **Role:** Serum iron biomarker PD output for iron-metabolism / iron-status turnover models. The state carries the iron concentration directly (umol/L in Angeli 2016) rather than an amount, so the synthesis and elimination rate constants `ksyn_iron` / `kout_iron` have concentration / time and 1 / time units respectively.
- **Source aliases:**
  - `Ir` -- Angeli 2016 paper notation; same orientation, no transformation.
- **Example models:** `Angeli_2016_iron_hepcidin.R`.
- **Notes:** Registered 2026-06-03 alongside the Angeli 2016 iron / hepcidin joint turnover extraction. Reserved for serum iron as a PD biomarker (concentration-state IDR / turnover models); paper-mechanistic intracellular iron pools or membrane-limited PBPK iron sub-compartments should use a distinct namespaced canonical when they arise.

### hep (**canonical serum hepcidin PD output**)
- **Type:** compartment
- **Role:** Serum hepcidin biomarker PD output for hepcidin-driven iron-regulation turnover models. The state carries the hepcidin concentration directly (nmol/L in Angeli 2016) and is paired with `iron` in the Angeli 2016 joint model via a multiplicative coupling on hepcidin synthesis.
- **Source aliases:**
  - `He` -- Angeli 2016 paper notation; same orientation, no transformation.
- **Example models:** `Angeli_2016_iron_hepcidin.R`.
- **Notes:** Registered 2026-06-03 alongside the Angeli 2016 iron / hepcidin joint turnover extraction. The short `hep` form was chosen to mirror the paper's `He` notation and to follow the existing biomarker-state precedent of short lowercase names (`phe`, `pth`, `thb`, `psa`).

### psa (**canonical prostate-specific antigen (lowercase form)**)
- **Type:** compartment
- **Role:** Prostate-specific antigen PD output (lowercase form alongside the canonical capitalised `PSA`).
- **Source aliases:** none.
- **Example models:** `Wilbaux_2015_PSA.R`.

### sld (**canonical sum of longest diameters PD output**)
- **Type:** compartment
- **Role:** Sum of longest diameters RECIST-style TGI endpoint.
- **Source aliases:** none.
- **Example models:** `Schindler_2016_sunitinib.R`.

### mbl (**canonical bare mean bacterial load**)
- **Type:** compartment
- **Role:** Bare-case alias of the registered observation `MBL` (Svensson 2017 bedaquiline).
- **Source aliases:** none.
- **Example models:** `Svensson_2017_bedaquiline.R`.

### aaa (**canonical bare AAA size alias**)
- **Type:** compartment
- **Role:** Bare-case alias of the registered observation `aaaSize` (Sherer 2012 AAA).
- **Source aliases:** none.
- **Example models:** `Sherer_2012_AAA.R`.

### serum (**canonical generic serum compartment**)
- **Type:** compartment
- **Role:** Generic serum compartment.
- **Source aliases:** none.
- **Example models:** `Aksenov_2018_uric_acid.R`.

---

## Cell populations and lymphoid tissues

### bcell (**canonical B-lymphocyte central pool**)
- **Type:** compartment
- **Role:** B-lymphocyte central pool.
- **Source aliases:** none.
- **Example models:** `Yu_2022_ofatumumab.R`.

### bcell_periph (**canonical B-lymphocyte peripheral pool**)
- **Type:** compartment
- **Role:** B-lymphocyte peripheral pool. Paired with `bcell`.
- **Source aliases:** none.
- **Example models:** `Yu_2022_ofatumumab.R`.

### pbmc (**canonical peripheral blood mononuclear cells**)
- **Type:** compartment
- **Role:** Peripheral blood mononuclear cells PD output.
- **Source aliases:** none.
- **Example models:** `Sampson_2014_azithromycin.R`.

### pmn (**canonical polymorphonuclear leukocytes**)
- **Type:** compartment
- **Role:** Polymorphonuclear leukocytes PD output.
- **Source aliases:** none.
- **Example models:** `Sampson_2014_azithromycin.R`.

### erythrocytes (**canonical erythrocyte pool**)
- **Type:** compartment
- **Role:** Red blood cell pool.
- **Source aliases:** none.
- **Example models:** `Dao_2020_sultiame.R`.

### cells (**canonical generic cell population**)
- **Type:** compartment
- **Role:** Generic cell population PD output.
- **Source aliases:** none.
- **Example models:** `Jager_2011_gemtuzumab.R`.

### lactotroph (**canonical anterior-pituitary lactotroph cells**)
- **Type:** compartment
- **Role:** Anterior-pituitary lactotroph cells.
- **Source aliases:** none.
- **Example models:** `Stevens_2012_remoxipride.R`.

### ncmc (**canonical non-classical monocytes**)
- **Type:** compartment
- **Role:** Non-classical monocytes PD output.
- **Source aliases:** none.
- **Example models:** `Yang_2024_axatilimab.R`.

---

## PD output and regulatory states

### urine_vol (**canonical urine-volume PD output**)
- **Type:** compartment
- **Role:** Urine volume PD output.
- **Source aliases:** none.
- **Example models:** `Heuberger_2018_salbutamol.R`.

### weight (**canonical bare body-weight PD output**)
- **Type:** compartment
- **Role:** Body weight PD output.
- **Source aliases:** none.
- **Example models:** `Choy_2016_T2DM.R`.

### fc (**canonical Fc-receptor pool**)
- **Type:** compartment
- **Role:** Fc-receptor / IgG-Fc pool.
- **Source aliases:** none.
- **Example models:** `Aguiar_2021_ustekinumab.R`.

### gut (**canonical bare gut compartment**)
- **Type:** compartment
- **Role:** Bare gut compartment. Bare alias of the PBPK `a_gut` canonical.
- **Source aliases:** none.
- **Example models:** `Charbonneau_2021_phenylalanine.R`.

### bacteria (**canonical generic bacterial pool**)
- **Type:** compartment
- **Role:** Generic bacterial pool. Bare alias of the registered `cfu` canonical.
- **Source aliases:** none.
- **Example models:** `Sadouki_2025_meropenem.R`.

---

## PD biomarker chain (Ait-Oudhia 2012 IL-1beta cascade)

The Ait-Oudhia 2012 canakinumab IL-1beta -> CRP transit cascade: `crp1` / `crp2` / `crp3` are the three CRP transit compartments feeding the acute-response level `acrl`.

### crp1 (**canonical CRP first transit compartment**)
- **Type:** compartment
- **Role:** First CRP transit compartment in Ait-Oudhia 2012 IL-1beta cascade.
- **Source aliases:** none.
- **Example models:** `AitOudhia_2012_canakinumab.R`.

### crp2 (**canonical CRP second transit compartment**)
- **Type:** compartment
- **Role:** Second CRP transit compartment.
- **Source aliases:** none.
- **Example models:** `AitOudhia_2012_canakinumab.R`.

### crp3 (**canonical CRP third transit compartment**)
- **Type:** compartment
- **Role:** Third CRP transit compartment.
- **Source aliases:** none.
- **Example models:** `AitOudhia_2012_canakinumab.R`.

### acrl (**canonical acute-response level**)
- **Type:** compartment
- **Role:** Acute-response level fed by the CRP transit cascade.
- **Source aliases:** none.
- **Example models:** `AitOudhia_2012_canakinumab.R`.

---

## Hansson 2013 sunitinib soluble biomarkers

### vegf (**canonical VEGF biomarker compartment**)
- **Type:** compartment
- **Role:** Soluble VEGF biomarker.
- **Source aliases:** none.
- **Example models:** `Hansson_2013_sunitinib.R`.

### skit (**canonical sKIT biomarker compartment**)
- **Type:** compartment
- **Role:** Soluble c-KIT biomarker.
- **Source aliases:** none.
- **Example models:** `Hansson_2013_sunitinib.R`.

### skit_drug (**canonical sKIT drug-arm output**)
- **Type:** compartment
- **Role:** sKIT drug-arm output (Hansson 2013b).
- **Source aliases:** none.
- **Example models:** `Hansson_2013b_sunitinib.R`.

### skit_pla (**canonical sKIT placebo-arm output**)
- **Type:** compartment
- **Role:** sKIT placebo-arm output (Hansson 2013b).
- **Source aliases:** none.
- **Example models:** `Hansson_2013b_sunitinib.R`.

---

## Keizer 2011 E7820 integrin biomarker

### integrin (**canonical alpha2-integrin biomarker compartment**)
- **Type:** compartment
- **Role:** Platelet alpha2-integrin expression turnover pool, the biomarker driven by an indirect-response model in the Keizer 2011 PK/PD analysis of E7820. Used both in the preclinical mouse model (units: % of reference platelet flow-cytometry signal) and in the clinical model (units: MESF, molecules of equivalent soluble fluorochrome) -- the same state with species-specific units.
- **Source aliases:** I (paper symbol).
- **Example models:** `Keizer_2011_E7820_mouse.R`, `Keizer_2011_E7820_human.R`.

---

## Survival / dropout cumulative-hazard compartments

### cumhaz_os (**canonical overall-survival cumulative-hazard**)
- **Type:** compartment
- **Role:** Overall-survival cumulative-hazard state in oncology TTE sub-models.
- **Source aliases:**
  - `cumHaz_os` -- prior canonical name (pre-2026-06-19 lowercase standardization).
- **Example models:** `Schindler_2016_sunitinib.R`.
- **Notes:** Renamed from `cumHaz_os` to `cumhaz_os` on 2026-06-19 per the canonical-register standardization audit (operator decision: compartment names follow the all-lowercase nlmixr2 convention; the cumulative-hazard family is now uniformly lowercase across `cumhaz`, `cumhaz_os`, `cumhaz_drop`).

### cumhaz_drop (**canonical dropout cumulative-hazard**)
- **Type:** compartment
- **Role:** Dropout cumulative-hazard state.
- **Source aliases:**
  - `cumHaz_drop` -- prior canonical name (pre-2026-06-19 lowercase standardization).
- **Example models:** `Schindler_2016_sunitinib.R`.
- **Notes:** Renamed from `cumHaz_drop` to `cumhaz_drop` on 2026-06-19 per the canonical-register standardization audit (lowercase compartment-name standardization).

---

## MBMA placebo / drug arm output compartments

The Li 2015 taspoglutide MBMA model maintains separate placebo and drug arms for each clinical endpoint. The placebo arm captures the background placebo response; the drug arm carries the drug-driven delta.

### fpg_placebo (**canonical fasting plasma glucose placebo arm**)
- **Type:** compartment
- **Role:** Fasting plasma glucose placebo-arm output state.
- **Source aliases:** none.
- **Example models:** `Li_2015_taspoglutide_MBMA.R`.

### fpg_drug (**canonical fasting plasma glucose drug arm**)
- **Type:** compartment
- **Role:** Fasting plasma glucose drug-arm output state.
- **Source aliases:** none.
- **Example models:** `Li_2015_taspoglutide_MBMA.R`.

### hba1c_placebo (**canonical HbA1c placebo arm**)
- **Type:** compartment
- **Role:** HbA1c placebo-arm output state.
- **Source aliases:** none.
- **Example models:** `Li_2015_taspoglutide_MBMA.R`.

### hba1c_drug (**canonical HbA1c drug arm**)
- **Type:** compartment
- **Role:** HbA1c drug-arm output state.
- **Source aliases:** none.
- **Example models:** `Li_2015_taspoglutide_MBMA.R`.

---

## Depot dosing-route compartments

The `depot_<route>` pattern distinguishes parallel dosing routes when a model carries more than one. Sibling of the canonical `depot` / numbered `depot1`, `depot2` forms.

### depot_im (**canonical intramuscular depot**)
- **Type:** compartment
- **Role:** Intramuscular depot used in parallel-route popPK models.
- **Source aliases:** none.
- **Example models:** `Dunn_2025_tranexamic_acid.R`.

### depot_oral (**canonical oral depot**)
- **Type:** compartment
- **Role:** Oral depot used in parallel-route popPK models alongside `depot_im`.
- **Source aliases:** none.
- **Example models:** `Dunn_2025_tranexamic_acid.R`.

### depot_brain (**canonical intranasal direct-to-brain depot**)
- **Type:** compartment
- **Role:** Intranasal direct-to-brain depot.
- **Source aliases:** none.
- **Example models:** `Stevens_2012_remoxipride.R`.

---

## K-PD virtual drug compartments

K-PD (kinetic-pharmacodynamic) models treat dose as entering a hypothetical body-amount compartment with first-order elimination and no measured drug concentration; the effect of the dose is driven by the amount in that compartment. The canonical K-PD virtual drug compartment is `depot_kpd`, with elimination rate constant `lkel` (log-transformed) / `kel` (bare) -- the same canonical pair used by single-rate-constant PK elimination. Combination K-PD models with two or more parallel K-PD compartments use the drug-suffixed form `depot_kpd_<drug>` (with paired `lkel_<drug>` / `kel_<drug>` rates), where `<drug>` is registered as a metabolite-suffix token below.

### depot_kpd (**canonical K-PD virtual drug compartment**)
- **Type:** compartment
- **Role:** Single hypothetical body-amount compartment used in K-PD models. Receives dose events; decays at first-order rate `kel`. The effect is driven directly by `depot_kpd` (or by `depot_kpd / vc` when a derived "drug-delivery rate concentration" is used).
- **Source aliases:**
  - `kpdConc` -- used in `Mazzocco_2015_temozolomide.R`.
  - `depot` (when the model has no extravascular absorption depot and the lone depot serves as the K-PD virtual drug compartment) -- used in `Shoji_2017_fosdagrocorat_oc.R`, `Shoji_2017_fosdagrocorat_p1np.R`, `vanHasselt_2015_eribulin.R`, `Xia_2024_warfarin.R`.
- **Example models:** `Mazzocco_2015_temozolomide.R`, `Shoji_2017_fosdagrocorat_oc.R`, `Shoji_2017_fosdagrocorat_p1np.R`, `vanHasselt_2015_eribulin.R`, `Xia_2024_warfarin.R`.
- **Notes:** Drug-suffixed variants `depot_kpd_<drug>` are accepted for combination K-PD models via the metabolite-suffix mechanism, where `<drug>` is a registered drug-name suffix below (e.g., `depot_kpd_sunitinib`, `depot_kpd_irinotecan` in Wilson 2015). Canonical `depot_kpd` adopted 2026-05-30 per the K-PD canonical-name retrofit (see `memory/kpd-model-canonical-standards.md`).

### sunitinib (**canonical sunitinib K-PD drug-name suffix**)
- **Type:** metabolite-suffix
- **Role:** Sunitinib drug-name suffix for combination K-PD compartments and rate constants (`depot_kpd_sunitinib`, `lkel_sunitinib`, `kel_sunitinib`).
- **Source aliases:** none.
- **Example models:** `Wilson_2015_sunitinib_irinotecan_mouse.R`.
- **Notes:** Full INN name (lowercase) rather than the TB-style 3-letter abbreviation because combination K-PD models are rare enough that semantic clarity wins; the abbreviated `sun` form is reserved for a future paper that needs it.

### irinotecan (**canonical irinotecan K-PD drug-name suffix**)
- **Type:** metabolite-suffix
- **Role:** Irinotecan drug-name suffix for combination K-PD compartments and rate constants (`depot_kpd_irinotecan`, `lkel_irinotecan`, `kel_irinotecan`).
- **Source aliases:** none.
- **Example models:** `Wilson_2015_sunitinib_irinotecan_mouse.R`.
- **Notes:** Full INN name (lowercase) for the same reason as `sunitinib`.

---

## Metabolite / sibling-drug / payload suffixes

These tokens may appear as a trailing `_<suffix>` on a canonical compartment, parameter, or residual-SD name to denote a non-parent species tracked alongside the parent. Examples: `central_mmae` (MMAE payload central compartment), `lcl_lesn` (lesinurad clearance in the Hill-McManus dual-urate-lowering-therapy model), `propSd_dxd` (Dxd payload residual proportional SD). The suffix matches via `endsWith(name, "_<suffix>")`; the prefix must be canonical under the relevant validator's check.

### mmae (**canonical MMAE payload suffix**)
- **Type:** metabolite-suffix
- **Role:** Monomethyl auristatin E (MMAE) ADC payload species suffix.
- **Source aliases:** none.
- **Example models:** ADC popPK extractions with MMAE payload.

### dxd (**canonical Dxd payload suffix**)
- **Type:** metabolite-suffix
- **Role:** DXd (exatecan derivative) ADC payload species suffix.
- **Source aliases:** none.
- **Example models:** ADC popPK extractions with DXd payload.

### sn38 (**canonical SN-38 payload suffix**)
- **Type:** metabolite-suffix
- **Role:** SN-38 (topoisomerase-I-inhibitor) ADC payload species suffix.
- **Source aliases:** none.
- **Example models:** ADC popPK extractions with SN-38 payload.

### dm4 (**canonical DM4 payload suffix**)
- **Type:** metabolite-suffix
- **Role:** DM4 (maytansinoid) ADC payload species suffix.
- **Source aliases:** none.
- **Example models:** ADC popPK extractions with DM4 payload.

### medm4 (**canonical Me-DM4 payload suffix**)
- **Type:** metabolite-suffix
- **Role:** Me-DM4 (S-methyl DM4) ADC payload metabolite species suffix.
- **Source aliases:** none.
- **Example models:** ADC popPK extractions with Me-DM4 metabolite.

### mcmmaf (**canonical MC-MMAF payload suffix**)
- **Type:** metabolite-suffix
- **Role:** MC-MMAF (maleimidocaproyl-monomethyl auristatin F) ADC payload species suffix.
- **Source aliases:** none.
- **Example models:** ADC popPK extractions with MC-MMAF payload.

### complex (**canonical complex suffix**)
- **Type:** metabolite-suffix
- **Role:** Generic drug-target complex species suffix used by TMDD models.
- **Source aliases:** none.
- **Example models:** TMDD popPK extractions.
- **Notes:** Same token as the bare `complex` compartment; both Types co-exist for the same canonical name.

### ige (**canonical IgE binding-partner suffix**)
- **Type:** metabolite-suffix
- **Role:** IgE binding-partner species suffix (omalizumab and similar).
- **Source aliases:** none.
- **Example models:** anti-IgE mAb popPK extractions.

### il1b (**canonical IL-1beta binding-partner suffix**)
- **Type:** metabolite-suffix
- **Role:** IL-1beta binding-partner species suffix (canakinumab).
- **Source aliases:** none.
- **Example models:** `AitOudhia_2012_canakinumab.R`.

### tab (**canonical total-antibody (TAb) suffix**)
- **Type:** metabolite-suffix
- **Role:** Total-antibody species suffix in ADC popPK (TAb = drug-loaded + unconjugated antibody).
- **Source aliases:** none.
- **Example models:** ADC popPK extractions tracking TAb.

### nab (**canonical neutralising-antibody (NAb) suffix**)
- **Type:** metabolite-suffix
- **Role:** Neutralising-antibody species suffix in immunogenicity-aware mAb popPK.
- **Source aliases:** none.
- **Example models:** mAb popPK extractions tracking NAb.

### dar0 (**canonical DAR-0 ADC isoform suffix**)
- **Type:** metabolite-suffix
- **Role:** Drug-to-antibody ratio 0 (unconjugated antibody) ADC isoform species suffix.
- **Source aliases:** none.
- **Example models:** DAR-resolved ADC popPK extractions.

### dar1 (**canonical DAR-1 ADC isoform suffix**)
- **Type:** metabolite-suffix
- **Role:** DAR-1 ADC isoform species suffix.
- **Source aliases:** none.
- **Example models:** DAR-resolved ADC popPK extractions.

### dar2 (**canonical DAR-2 ADC isoform suffix**)
- **Type:** metabolite-suffix
- **Role:** DAR-2 ADC isoform species suffix.
- **Source aliases:** none.
- **Example models:** DAR-resolved ADC popPK extractions.

### dar3 (**canonical DAR-3 ADC isoform suffix**)
- **Type:** metabolite-suffix
- **Role:** DAR-3 ADC isoform species suffix.
- **Source aliases:** none.
- **Example models:** DAR-resolved ADC popPK extractions.

### dar4 (**canonical DAR-4 ADC isoform suffix**)
- **Type:** metabolite-suffix
- **Role:** DAR-4 ADC isoform species suffix.
- **Source aliases:** none.
- **Example models:** DAR-resolved ADC popPK extractions.

### dar5 (**canonical DAR-5 ADC isoform suffix**)
- **Type:** metabolite-suffix
- **Role:** DAR-5 ADC isoform species suffix.
- **Source aliases:** none.
- **Example models:** DAR-resolved ADC popPK extractions.

### dar6 (**canonical DAR-6 ADC isoform suffix**)
- **Type:** metabolite-suffix
- **Role:** DAR-6 ADC isoform species suffix.
- **Source aliases:** none.
- **Example models:** DAR-resolved ADC popPK extractions.

### dar7 (**canonical DAR-7 ADC isoform suffix**)
- **Type:** metabolite-suffix
- **Role:** DAR-7 ADC isoform species suffix.
- **Source aliases:** none.
- **Example models:** DAR-resolved ADC popPK extractions.

### dar8 (**canonical DAR-8 ADC isoform suffix**)
- **Type:** metabolite-suffix
- **Role:** DAR-8 ADC isoform species suffix.
- **Source aliases:** none.
- **Example models:** DAR-resolved ADC popPK extractions.

---

## Small-molecule metabolite suffixes

### 3oh (**canonical 3-hydroxy-agomelatine suffix**)
- **Type:** metabolite-suffix
- **Role:** 3-hydroxy-agomelatine metabolite of agomelatine (Xie 2019).
- **Source aliases:** none.
- **Example models:** `Xie_2019_agomelatine.R`.
- **Notes:** Suffix starts with a digit; the convention check matches on `endsWith(name, "_<metab>")` rather than treating the metabolite name as an R identifier.

### 7dm (**canonical 7-desmethyl-agomelatine suffix**)
- **Type:** metabolite-suffix
- **Role:** 7-desmethyl-agomelatine metabolite of agomelatine (Xie 2019).
- **Source aliases:** none.
- **Example models:** `Xie_2019_agomelatine.R`.

### 8oh (**canonical 8-hydroxy-efavirenz suffix**)
- **Type:** metabolite-suffix
- **Role:** 8-hydroxy-efavirenz, the primary CYP2B6-formed metabolite of efavirenz; further metabolised by CYP2B6 and UGT2B7 (Heathman 2024).
- **Source aliases:** none.
- **Example models:** `Heathman_2024_efavirenz.R`.
- **Notes:** Suffix starts with a digit; the convention check matches on `endsWith(name, "_<metab>")` rather than treating the metabolite name as an R identifier. Founding example: `Heathman_2024_efavirenz.R` (full 3-analyte EFV + 8-OH + 7-OH popPK with CYP2B6 / CYP2A6 enzyme-turnover autoinduction).

### 7oh (**canonical 7-hydroxy-efavirenz suffix**)
- **Type:** metabolite-suffix
- **Role:** 7-hydroxy-efavirenz, the CYP2A6-formed metabolite of efavirenz (Heathman 2024).
- **Source aliases:** none.
- **Example models:** `Heathman_2024_efavirenz.R`.
- **Notes:** Suffix starts with a digit; the convention check matches on `endsWith(name, "_<metab>")` rather than treating the metabolite name as an R identifier. Founding example: `Heathman_2024_efavirenz.R`.

### m1 (**canonical paper-named M1 metabolite suffix**)
- **Type:** metabolite-suffix
- **Role:** Generic M1 metabolite suffix used by parent + metabolite popPK extractions in which the source paper names the major active metabolite simply "M1" without supplying a chemical name. Each model file's `description` / `reference` text identifies the specific compound; the `m1` suffix is the shared notational token. Disambiguation by drug context: in `Lehr_2010_tesofensine.R`, M1 is the major CYP3A4-formed metabolite of tesofensine.
- **Source aliases:** none.
- **Example models:** `Lehr_2010_tesofensine.R` (tesofensine M1, CYP3A4-formed, in vivo five-fold lower dopamine-reuptake potency than parent per Lehr 2010 Methods reference 17).

### m2 (**canonical N-desmethyl-bedaquiline (M2) suffix**)
- **Type:** metabolite-suffix
- **Role:** N-desmethyl-bedaquiline (M2) metabolite of bedaquiline.
- **Source aliases:** none.
- **Example models:** `Svensson_2013_bedaquiline.R`, `Svensson_2016_bedaquiline.R` (DDMODEL00000219).

### m3 (**canonical N,N-bis-desmethyl-bedaquiline (M3) suffix**)
- **Type:** metabolite-suffix
- **Role:** N,N-bis-desmethyl-bedaquiline (M3) metabolite of bedaquiline; the downstream demethylation product of M2 (responsible enzyme(s) not identified in vitro but suspected CYP3A4-mediated demethylation by analogy with the BDQ -> M2 step).
- **Source aliases:** none.
- **Example models:** `Svensson_2013_bedaquiline.R`.
- **Notes:** Distinct from `m3g` (morphine-3-glucuronide) — the suffix matcher uses `endsWith(name, "_m3")` vs `endsWith(name, "_m3g")` and these do not collide. Registered alongside the Svensson 2013 bedaquiline extraction (the first BDQ paper to model the M3 metabolite).

### m8 (**canonical hydroxy-tert-butylamide (M8) suffix**)
- **Type:** metabolite-suffix
- **Role:** Hydroxy-tert-butylamide (M8) active metabolite of nelfinavir; formed by CYP2C19-mediated hydroxylation of nelfinavir and eliminated by CYP3A4. Equipotent to the parent drug against HIV-1 protease.
- **Source aliases:** none.
- **Example models:** `Hirt_2006_nelfinavir.R`.

### endox (**canonical endoxifen suffix**)
- **Type:** metabolite-suffix
- **Role:** Endoxifen (4-hydroxy-N-desmethyltamoxifen), major active metabolite of tamoxifen.
- **Source aliases:**
  - `endx` -- prior canonical name (pre-2026-06-19 readability standardization).
- **Example models:** `TerHeine_2014_tamoxifen.R`.
- **Notes:** Renamed from `endx` to `endox` on 2026-06-19 per the canonical-register standardization audit (operator decision: keep the contracted form but include the "o" so the suffix is readable as "endox" without vowel-stripping; `endx` was an opaque consonant cluster).

### megx (**canonical MEGX lidocaine metabolite suffix**)
- **Type:** metabolite-suffix
- **Role:** Monoethylglycinexylidide (MEGX) lidocaine metabolite (LID -> MEGX via CYP1A2/3A4).
- **Source aliases:** none.
- **Example models:** `DDMODEL00000281.R`, `NA_NA_lidocaine.R`.

### gx (**canonical GX lidocaine metabolite suffix**)
- **Type:** metabolite-suffix
- **Role:** Glycinexylidide (GX) lidocaine metabolite (MEGX -> GX).
- **Source aliases:** none.
- **Example models:** `DDMODEL00000281.R`, `NA_NA_lidocaine.R`.

### xyl (**canonical 2,6-xylidide lidocaine metabolite suffix**)
- **Type:** metabolite-suffix
- **Role:** 2,6-Xylidide (LID -> 2,6-XYL minor pathway) lidocaine metabolite.
- **Source aliases:** none.
- **Example models:** `DDMODEL00000281.R`, `NA_NA_lidocaine.R`.

### m3g (**canonical morphine-3-glucuronide suffix**)
- **Type:** metabolite-suffix
- **Role:** Morphine-3-glucuronide, major glucuronide metabolite of morphine.
- **Source aliases:** none.
- **Example models:** `Knibbe_2009_morphine.R` (DDMODEL00000248).

### m6g (**canonical morphine-6-glucuronide suffix**)
- **Type:** metabolite-suffix
- **Role:** Morphine-6-glucuronide, major glucuronide metabolite of morphine.
- **Source aliases:** none.
- **Example models:** `Knibbe_2009_morphine.R` (DDMODEL00000248).

### gluc (**canonical glucuronide phase-II conjugate suffix**)
- **Type:** metabolite-suffix
- **Role:** Phase-II glucuronide conjugate suffix (Allegaert 2015 paracetamol-glucuronide).
- **Source aliases:** none.
- **Example models:** `Allegaert_2015_paracetamol.R` (DDMODEL00000267).

### sulf (**canonical sulphate phase-II conjugate suffix**)
- **Type:** metabolite-suffix
- **Role:** Phase-II sulphate conjugate suffix (Allegaert 2015 paracetamol-sulphate).
- **Source aliases:** none.
- **Example models:** `Allegaert_2015_paracetamol.R` (DDMODEL00000267).

### apapg (**canonical APAP-glucuronide suffix**)
- **Type:** metabolite-suffix
- **Role:** Paracetamol (APAP) phase-II glucuronide metabolite. Used across all paracetamol PBPK / popPK models for the glucuronide species, regardless of source-paper notation.
- **Source aliases:**
  - `ag` -- deprecated Zurlinden 2016 shorthand (collides with silver chemistry symbol); migrated to `apapg` on 2026-06-19.
- **Example models:** `Cook_2016_paracetamol.R` (DDMODEL00000271), `Zurlinden_2016_paracetamol.R` (migrated from `ag` to `apapg` on 2026-06-19).
- **Notes:** On 2026-06-19 the Zurlinden 2016 paracetamol PBPK model was migrated from its prior `ag` shorthand to this canonical, per the canonical-register standardization audit (chemistry-symbol collision fix: `Ag` = silver).

### apaps (**canonical APAP-sulphate suffix**)
- **Type:** metabolite-suffix
- **Role:** Paracetamol (APAP) phase-II sulphate metabolite. Used across all paracetamol PBPK / popPK models for the sulphate species, regardless of source-paper notation.
- **Source aliases:**
  - `as` -- deprecated Zurlinden 2016 shorthand (collides with R reserved word `as.numeric` etc.); migrated to `apaps` on 2026-06-19.
- **Example models:** `Cook_2016_paracetamol.R` (DDMODEL00000271), `Zurlinden_2016_paracetamol.R` (migrated from `as` to `apaps` on 2026-06-19).
- **Notes:** On 2026-06-19 the Zurlinden 2016 paracetamol PBPK model was migrated from its prior `as` shorthand to this canonical, per the canonical-register standardization audit (R-reserved-word collision fix: bare `as` clashes with `as.numeric`, `as.integer`, `as.character`, etc.).

### col (**canonical colistin metabolite suffix**)
- **Type:** metabolite-suffix
- **Role:** Colistin, the active polymyxin generated in vivo by hydrolysis of the prodrug colistimethate sodium (CMS).
- **Source aliases:** none.
- **Example models:** `LeuppiTaegtmeyer_2019_CMS.R` (DDMODEL00000295).
- **Notes:** Same token as the bare `col` drug-state compartment; both Types co-exist for the same canonical name.

### dihydroart (**canonical dihydroartemisinin suffix**)
- **Type:** metabolite-suffix
- **Role:** Dihydroartemisinin, active metabolite of artesunate.
- **Source aliases:**
  - `dha` -- prior canonical name (pre-2026-06-19 readability standardization).
- **Example models:** `Birgersson_2019_artesunate.R` (DDMODEL00000297).
- **Notes:** Renamed from `dha` to `dihydroart` on 2026-06-19 per the canonical-register standardization audit (operator decision: the three-letter abbreviation is widely used in the malaria literature but is ambiguous - `DHA` collides with docosahexaenoic acid in nutrition contexts; `dihydroart` is unambiguous while still contracted relative to the full "dihydroartemisinin" spelling).

### ohi (**canonical hydroxy-itraconazole suffix**)
- **Type:** metabolite-suffix
- **Role:** Hydroxy-itraconazole (OH-ITZ), major active metabolite of itraconazole produced by CYP3A4 hydroxylation.
- **Source aliases:** none.
- **Example models:** `Hennig_2006_itraconazole.R`, `Hennig_2007_itraconazole.R`.

### ohcla (**canonical 14-(R)-hydroxy-clarithromycin suffix**)
- **Type:** metabolite-suffix
- **Role:** 14-(R)-hydroxy-clarithromycin (14-OH-CLA), the principal active metabolite of clarithromycin formed predominantly by CYP3A4-mediated stereospecific hydroxylation at the 14 position. Used as the metabolite suffix in parent + metabolite simultaneous popPK models (compartments `central_ohcla`, parameters `lcl_ohcla` / `lvc_ohcla`, residuals `addSd_ohcla`). Founding example: `Abduljalil_2009_clarithromycin.R`.
- **Source aliases:** none.
- **Example models:** `Abduljalil_2009_clarithromycin.R`.

### doxol (**canonical doxorubicinol suffix**)
- **Type:** metabolite-suffix
- **Role:** Doxorubicinol, C-13 alcohol metabolite of doxorubicin.
- **Source aliases:** none.
- **Example models:** `Kunarajah_2017_doxorubicin.R`.

### dol (**canonical daunorubicinol suffix**)
- **Type:** metabolite-suffix
- **Role:** Daunorubicinol (DOL), C-13 alcohol metabolite of daunorubicin formed primarily by carbonyl reductase 1 (CBR1) in adult AML patients.
- **Source aliases:** none.
- **Example models:** `Varatharajan_2016_daunorubicin.R` (doi:10.1007/s00280-016-3166-8).

### desacetylrbn (**canonical 25-O-desacetyl rifabutin suffix**)
- **Type:** metabolite-suffix
- **Role:** 25-O-desacetyl rifabutin, primary active metabolite of rifabutin formed by arylacetamide deacetylase.
- **Source aliases:**
  - `desrbn` -- prior canonical name (pre-2026-06-19 readability standardization).
- **Example models:** `Hennig_2015_rifabutin.R` (doi:10.1128/AAC.01195-15).
- **Notes:** Renamed from `desrbn` to `desacetylrbn` on 2026-06-19 per the canonical-register standardization audit (operator decision: spell out the `desacetyl` chemical descriptor rather than the opaque `des` abbreviation; the parallel `desacetylrpt` form would similarly replace `desrpt` if a future audit extends the rule).

### desrpt (**canonical 25-O-desacetyl rifapentine suffix**)
- **Type:** metabolite-suffix
- **Role:** 25-O-desacetyl rifapentine (25-DRFP), primary active metabolite of rifapentine, formed by enzymatic deacetylation; microbiologically active against Mycobacterium tuberculosis. Parallel naming to the registered `desrbn` (25-O-desacetyl rifabutin) suffix.
- **Source aliases:** `25-DRFP` / `25DRFP` / `metabolite_M` (paper narrative in Zvada 2010 Methods and Figure 1 caption).
- **Example models:** `Zvada_2010_rifapentine.R` (doi:10.1128/AAC.00345-10).

### az5104 (**canonical AZ5104 osimertinib metabolite suffix**)
- **Type:** metabolite-suffix
- **Role:** AZ5104 (N-desmethyl osimertinib), active EGFR-inhibitor metabolite of osimertinib formed predominantly via CYP3A4/5.
- **Source aliases:** none.
- **Example models:** `Brown_2017_osimertinib.R` (doi:10.1111/bcp.13223).

### ndmsel (**canonical N-desmethyl-selumetinib suffix**)
- **Type:** metabolite-suffix
- **Role:** N-desmethyl-selumetinib, active selumetinib metabolite (~3-5-fold more potent for MEK1 inhibition than parent), formed by oxidative N-demethylation.
- **Source aliases:**
  - `ndsel` -- prior canonical name (pre-2026-06-19 readability standardization).
- **Example models:** `Patel_2017_selumetinib.R` (doi:10.1002/psp4.12175).
- **Notes:** Renamed from `ndsel` to `ndmsel` on 2026-06-19 per the canonical-register standardization audit (operator decision: insert `m` so the contracted form reads as `n-desmethyl-sel` rather than `n-des-sel`; matches the N-desmethyl-`<drug>` pattern used by other contracted suffixes).

### dfcr (**canonical 5'-DFCR capecitabine metabolite suffix**)
- **Type:** metabolite-suffix
- **Role:** 5'-deoxy-5-fluorocytidine (5'-DFCR), formed in the liver by carboxylesterase from capecitabine.
- **Source aliases:** none.
- **Example models:** `Urien_2005_capecitabine.R` (doi:10.1007/s10928-005-0018-2).

### dfur (**canonical 5'-DFUR capecitabine metabolite suffix**)
- **Type:** metabolite-suffix
- **Role:** 5'-deoxy-5-fluorouridine (5'-DFUR), formed from 5'-DFCR by cytidine deaminase in liver and tumour cells.
- **Source aliases:** none.
- **Example models:** `Urien_2005_capecitabine.R`.

### 5fu (**canonical 5-fluorouracil capecitabine metabolite suffix**)
- **Type:** metabolite-suffix
- **Role:** 5-fluorouracil (5-FU), formed from 5'-DFUR by thymidine phosphorylase preferentially in tumour tissue.
- **Source aliases:** none.
- **Example models:** `Urien_2005_capecitabine.R`, `Blesch_2003_capecitabine.R`.

### fbal (**canonical alpha-fluoro-beta-alanine capecitabine catabolite suffix**)
- **Type:** metabolite-suffix
- **Role:** alpha-fluoro-beta-alanine (FBAL), the terminal renally-excreted catabolite of 5-FU produced by the dihydropyrimidine-dehydrogenase / dihydropyrimidinase / beta-ureidopropionase cascade. Plasma FBAL is the most abundant capecitabine-derived species and tracks systemic catabolic capacity; its renal clearance is a function of creatinine clearance.
- **Source aliases:** none.
- **Example models:** `Blesch_2003_capecitabine.R`.

### asn1 (**canonical AS(N-1)3' siRNA truncated antisense suffix**)
- **Type:** metabolite-suffix
- **Role:** AS(N-1)3' truncated antisense strand of GalNAc-conjugated siRNAs (givosiran and similar). Treated as the active metabolite equipotent with the parent for RISC loading and target mRNA silencing.
- **Source aliases:** none.
- **Example models:** `Ayyar_2024_givosiran.R` (doi:10.1016/j.xphs.2023.10.026).

### dfdu (**canonical 2',2'-difluorodeoxyuridine gemcitabine metabolite suffix**)
- **Type:** metabolite-suffix
- **Role:** 2',2'-difluorodeoxyuridine (dFdU), the principal inactive plasma metabolite of gemcitabine produced by cytidine deaminase in liver, kidney, blood, and other tissues; about 99% of dFdU is renally excreted (Abbruzzese 1991). In NONMEM ADVAN6 parent-metabolite popPK models dFdU parameters are reported as apparent values (CL_dFdU/F, Q_dFdU/F, V_C,dFdU/F, V_P,dFdU/F) because the fraction of gemcitabine converted to dFdU is not separately identifiable.
- **Source aliases:** none.
- **Example models:** `Jiang_2008_gemcitabine.R` (doi:10.1111/j.1365-2125.2007.03040.x).

### mhd (**canonical 10-monohydroxy oxcarbazepine suffix**)
- **Type:** metabolite-suffix
- **Role:** 10-monohydroxy derivative (MHD, "10-hydroxy-carbazepine"), primary active metabolite of oxcarbazepine produced by cytosolic arylketone reductases.
- **Source aliases:** none.
- **Example models:** `Rodrigues_2017_oxcarbazepine.R` (doi:10.1111/bcp.13392).

### r_enant (**canonical R-enantiomer suffix**)
- **Type:** metabolite-suffix
- **Role:** R-enantiomer suffix for enantiomer-resolved popPK models in which both enantiomers are followed in plasma but no interconversion is modelled.
- **Source aliases:**
  - `r` -- prior canonical name (pre-2026-06-19 disambiguation standardization).
- **Example models:** `Valitalo_2017_ketorolac.R` (doi:10.1111/bcp.13311).
- **Notes:** Treated as a "non-parent analyte" suffix; neither enantiomer is the parent. Renamed from bare `r` to `r_enant` on 2026-06-19 per the canonical-register standardization audit (operator decision: a single-letter `r` is too easy to collide with an unrelated short identifier or with R-language conventions; the `_enant` suffix makes the enantiomer meaning explicit and matches the `_enant` form used in similar metabolite-suffix registries).

### s_enant (**canonical S-enantiomer suffix**)
- **Type:** metabolite-suffix
- **Role:** S-enantiomer suffix for enantiomer-resolved popPK models.
- **Source aliases:**
  - `s` -- prior canonical name (pre-2026-06-19 disambiguation standardization).
- **Example models:** `Valitalo_2017_ketorolac.R`.
- **Notes:** Renamed from bare `s` to `s_enant` on 2026-06-19 per the canonical-register standardization audit (operator decision: a single-letter `s` is too easy to collide with an unrelated short identifier; the `_enant` suffix makes the enantiomer meaning explicit and matches sibling `r_enant`).

### noxide (**canonical roflumilast N-oxide suffix**)
- **Type:** metabolite-suffix
- **Role:** Roflumilast N-oxide, active metabolite of roflumilast (contributes ~90% of total PDE4 inhibitory activity).
- **Source aliases:** none.
- **Example models:** `Lahu_2010_roflumilast.R` (doi:10.2165/11536600-000000000-00000).

### norcloz (**canonical norclozapine (N-desmethylclozapine) suffix**)
- **Type:** metabolite-suffix
- **Role:** Norclozapine (N-desmethylclozapine), the primary pharmacologically active metabolite of clozapine formed predominantly by CYP1A2 (with secondary contributions from CYP2C19, CYP3A4, CYP2C9, and CYP2D6). Norclozapine retains receptor affinity at multiple monoaminergic and muscarinic targets and is routinely measured alongside clozapine in therapeutic-drug-monitoring (TDM) practice; the parent-to-metabolite ratio is itself a clinical descriptor of CYP1A2 activity.
- **Source aliases:** none.
- **Example models:** `Li_2012_clozapine.R` (doi:10.1038/aps.2012.71).

### cysmer (**canonical APAP cysteine+mercapturate suffix**)
- **Type:** metabolite-suffix
- **Role:** Combined acetaminophen cysteine + mercapturate compartment used by CYP2E1-oxidation popPK models that lump the two oxidation metabolites (the species are in rapid equilibrium with overlapping disposition).
- **Source aliases:** none.
- **Example models:** `vanRongen_2016_paracetamol.R` (doi:10.1007/s40262-015-0357-0).

### cpg2 (**canonical glucarpidase (CPG2) suffix**)
- **Type:** metabolite-suffix
- **Role:** Glucarpidase (CPG2), bacterial carboxypeptidase G2 enzyme given as rescue therapy after high-dose methotrexate. Co-administered perpetrator, not a metabolite of MTX.
- **Source aliases:** none.
- **Example models:** `Kimura_2023_methotrexate.R` (doi:10.21873/anticanres.16351).

### acdap (**canonical 3-N-acetyl-3,4-diaminopyridine suffix**)
- **Type:** metabolite-suffix
- **Role:** 3-N-acetyl-3,4-diaminopyridine (3-Ac DAP), inactive N-acetyl metabolite of 3,4-diaminopyridine (amifampridine) free base produced by N-acetyltransferases.
- **Source aliases:** none.
- **Example models:** `Thakkar_2017_amifampridine.R` (doi:10.1002/psp4.12218).

### h4 (**canonical clopidogrel H4 active thiol suffix**)
- **Type:** metabolite-suffix
- **Role:** H4 (active thiol) metabolite of clopidogrel: pharmacologically active species responsible for P2Y12 receptor inhibition, formed via sequential CYP-mediated oxidation (CYP2C19 dominant for second oxidation step). "H4" refers to the H4 stereoisomer specifically.
- **Source aliases:** none.
- **Example models:** `Danielak_2017_clopidogrel.R` (doi:10.1007/s00228-017-2334-z).

### mpag (**canonical mycophenolic acid glucuronide suffix**)
- **Type:** metabolite-suffix
- **Role:** Mycophenolic acid glucuronide (MPAG, 7-O-glucuronide phase II metabolite of mycophenolic acid produced by UGT1A9 and UGT2B7). Major plasma metabolite of mycophenolic acid after MMF dosing in renal transplant recipients.
- **Source aliases:** none.
- **Example models:** `deWinter_2009_mycophenolic.R` (doi:10.1007/s10928-009-9136-6).

### mpa (**canonical mycophenolic acid sibling-drug suffix**)
- **Type:** metabolite-suffix
- **Role:** Mycophenolic acid (MPA, the active moiety of co-administered mycophenolate mofetil), a co-medication sibling tracked alongside the parent tacrolimus in the TAC-MMF drug-drug-interaction model. Drives `depot_mpa` / `central_mpa` / `peripheral1_mpa` and the `propSd_mpa` residual.
- **Source aliases:** none.
- **Example models:** `Kim_2018_tacrolimus.R` (doi:10.1038/s41598-018-20071-3).

### acmpag (**canonical MPA acyl-glucuronide suffix**)
- **Type:** metabolite-suffix
- **Role:** Mycophenolic acid acyl-glucuronide (AcMPAG, the minor ~15% acyl-glucuronide phase II metabolite of mycophenolic acid). Distinct from the 7-O-glucuronide `mpag` suffix. Drives `central_acmpag` and the `propSd_acmpag` residual.
- **Source aliases:** none.
- **Example models:** `Kim_2018_tacrolimus.R`.

### asv (**canonical asunaprevir sibling-drug suffix**)
- **Type:** metabolite-suffix
- **Role:** Asunaprevir (ASV, NS3/4A protease inhibitor), a sibling direct-acting antiviral (DAA) co-administered with daclatasvir; neither is the "parent". Drives `depot_asv` / `central_asv` / `peripheral1_asv` / `effect_asv` and the `propSd_asv` / `addSd_asv` residuals.
- **Source aliases:** none.
- **Example models:** `Wang_2018_daclatasvir_asunaprevir.R` (doi:10.1038/aps.2017.84).

### cdb4453 (**canonical CDB-4453 monodemethylated-telapristone metabolite suffix**)
- **Type:** metabolite-suffix
- **Role:** Active monodemethylated metabolite of telapristone (CDB-4124). Removal of one N-methyl group on the C-17 side chain produces CDB-4453, a more polar (smaller apparent volume of distribution) metabolite with possible equipotent antiprogestational activity in vivo (Morris 2011 Discussion). Drives `central_cdb4453` and the `propSd_cdb4453` residual; the parent-side parameters use the canonical unsuffixed names (`lcl_pop1` / `lcl_pop2` / `lvc` / `lvp` / `lq` / `lka`).
- **Source aliases:** none (the paper uses the compound code CDB-4453 throughout).
- **Example models:** `Morris_2011_telapristone.R` (doi:10.1208/s12248-011-9304-7).
- **Notes:** The Morris 2011 model fixes the metabolite apparent volume V3/F to 1 L for identifiability (Fmet not separately identifiable from V3), so the estimated `fmetest` (= Fmet / V3_metab, units 1/L) numerically equals Fmet under that constraint. The metabolite compartment `central_cdb4453` therefore numerically equals the metabolite concentration (nmol/L) when V3 = 1 L. Ratified canonically on 2026-06-09 alongside the Morris 2011 telapristone extraction.

---

## Cell-type suffixes (Friberg multi-cell-type chains)

Cell-type suffixes used with Friberg-style `circ_<celltype>` myelosuppression compartments and `precursor1_<celltype>` ... `precursorN_<celltype>` maturation chains for paired-output multi-cell models. Registered 2026-05-28 per the naming audit.

### anc (**canonical ANC cell-type suffix**)
- **Type:** metabolite-suffix
- **Role:** Absolute neutrophil count cell-type suffix.
- **Source aliases:** none.
- **Example models:** `Han_2015_decitabine.R`.

### plt (**canonical platelet cell-type suffix**)
- **Type:** metabolite-suffix
- **Role:** Platelet cell-type suffix.
- **Source aliases:** none.
- **Example models:** Multi-cell-type myelosuppression PD models.

### wbc (**canonical WBC cell-type suffix**)
- **Type:** metabolite-suffix
- **Role:** White blood cell cell-type suffix.
- **Source aliases:** none.
- **Example models:** Multi-cell-type myelosuppression PD models.

---

## Parent-drug suffixes (urine excretion)

Parent-drug suffixes used with `urine_<X>` (or `central_<X>`) excretion compartments in parent + metabolite renal-elimination models where the parent itself is also tracked in urine. Registered 2026-05-28 per the naming audit.

### apap (**canonical paracetamol parent-drug suffix**)
- **Type:** metabolite-suffix
- **Role:** Paracetamol (APAP) parent-drug suffix for `urine_apap` tracking.
- **Source aliases:** none.
- **Example models:** `Allegaert_2015_paracetamol.R`.

### morphine (**canonical morphine parent-drug suffix**)
- **Type:** metabolite-suffix
- **Role:** Morphine parent-drug suffix for `urine_morphine` tracking.
- **Source aliases:** none.
- **Example models:** `Pierre_2017_morphine.R`.

---

## Sibling-drug suffixes (urate-lowering therapy)

Sibling-drug suffixes for the Hill-McManus 2017 dual-urate-lowering-therapy PKPD model, where febuxostat and lesinurad are co-administered and neither is the "parent". Both PK subsystems use canonical compartment / PK-param names with the drug suffix.

### febx (**canonical febuxostat sibling-drug suffix**)
- **Type:** metabolite-suffix
- **Role:** Febuxostat (xanthine oxidase inhibitor) sibling-drug suffix.
- **Source aliases:** none.
- **Example models:** `Hill-McManus_2017_uricLT.R` (doi:10.1111/bcp.13427).

### lesn (**canonical lesinurad sibling-drug suffix**)
- **Type:** metabolite-suffix
- **Role:** Lesinurad (URAT1 uricosuric) sibling-drug suffix.
- **Source aliases:** none.
- **Example models:** `Hill-McManus_2017_uricLT.R`.

---

## Neuromuscular-blocking-agent / reversal-agent sibling-drug suffixes

Sibling-drug suffix for the Kleijn 2011 sugammadex-mediated reversal of rocuronium-induced neuromuscular blockade PK-PD model, where sugammadex (the reversal agent, mentioned first in the paper title) is the unsuffixed parent and rocuronium (the substrate aminosteroid NMBA) carries the `roc` suffix throughout. The model also uses the existing `complex` registered suffix for the sugammadex-rocuronium inclusion complex compartments (with PK set equal to free sugammadex), so all three species coexist as separate compartment chains.

### roc (**canonical rocuronium sibling-drug suffix**)
- **Type:** metabolite-suffix
- **Role:** Rocuronium (aminosteroid neuromuscular blocking agent) sibling-drug suffix, paired with sugammadex as the unsuffixed parent in the Kleijn 2011 reversal PK-PD model. Drives `central_roc` / `peripheral1_roc` / `effect_roc` compartments, `lcl_roc` / `lvc_roc` / `lq_roc` / `lvp_roc` PK parameters, and the `propSd_roc` residual on total rocuronium plasma concentration. The effect compartment carries the rocuronium concentration at the neuromuscular junction that drives the sigmoid-Emax NMB readout; sugammadex itself has no effect compartment because its NMB-reversal action enters the model as an extra elimination route on `effect_roc`.
- **Source aliases:** none.
- **Example models:** `Kleijn_2011_sugammadex_rocuronium.R` (doi:10.1111/j.1365-2125.2011.04000.x).

---

## Combination antimalarial / antibiotic sibling-drug suffixes

### pyra (**canonical pyrimethamine sibling-drug suffix**)
- **Type:** metabolite-suffix
- **Role:** Pyrimethamine sibling-drug suffix (paired with sulfadoxine in Odongo 2015 / deKock 2017 sulfadoxine-pyrimethamine models). Drives `depot_pyra` / `central_pyra` / `peripheral1_pyra` PK subsystem.
- **Source aliases:** none.
- **Example models:** `Odongo_2015_SDX_PYR.R`, `deKock_2017_SDX_PYR.R`.

### mer (**canonical meropenem sibling-drug suffix**)
- **Type:** metabolite-suffix
- **Role:** Meropenem sibling-drug suffix (paired with gentamicin / ciprofloxacin in Sadouki 2025 and with linezolid / vancomycin in Wicha 2017).
- **Source aliases:** none.
- **Example models:** `Sadouki_2025_meropenem.R`, `Wicha_2017_linezolid_meropenem_vancomycin.R`.
- **Notes:** Same token as the bare `mer` drug-state compartment; both Types co-exist for the same canonical name.

---

## Deprecated Zurlinden 2016 paracetamol PBPK metabolite suffixes (subsumed into Cook 2016 forms)

Zurlinden 2016 paracetamol PBPK shorthand suffixes `as` and `ag` were deprecated on 2026-06-19 because they collide with R reserved words (`as.numeric`, `as.integer`, `as.character`, etc.) and chemistry symbols (Ag = silver). The Zurlinden 2016 paracetamol PBPK model is migrated to use the existing Cook 2016 `apaps` / `apapg` canonicals for the same chemical species. See the `apaps` / `apapg` entries above for the canonical names.

### as (**DEPRECATED -- use `apaps` instead**)
- **Type:** metabolite-suffix (deprecated)
- **Role:** Formerly APAP-sulfate suffix (Zurlinden notation). Same chemical species as `apaps` (Cook 2016 notation). Deprecated on 2026-06-19 because the bare `as` token collides with R reserved words (`as.numeric`, `as.integer`, etc.).
- **Source aliases:**
  - `as` -- deprecated canonical name used in `Zurlinden_2016_paracetamol.R` prior to the 2026-06-19 audit.
- **Example models:** `Zurlinden_2016_paracetamol.R`.
- **Notes:** Do not use for new models. The Zurlinden 2016 paracetamol PBPK model was migrated to `apaps` on 2026-06-19 per the canonical-register standardization audit (R-reserved-word collision fix).

### ag (**DEPRECATED -- use `apapg` instead**)
- **Type:** metabolite-suffix (deprecated)
- **Role:** Formerly APAP-glucuronide suffix (Zurlinden notation). Same chemical species as `apapg` (Cook 2016 notation). Deprecated on 2026-06-19 because the bare `ag` token collides with the silver chemical symbol (Ag).
- **Source aliases:**
  - `ag` -- deprecated canonical name used in `Zurlinden_2016_paracetamol.R` prior to the 2026-06-19 audit.
- **Example models:** `Zurlinden_2016_paracetamol.R`.
- **Notes:** Do not use for new models. The Zurlinden 2016 paracetamol PBPK model was migrated to `apapg` on 2026-06-19 per the canonical-register standardization audit (chemistry-symbol collision fix).

---

## Per-paper additions (2026-05-28 naming audit)

Per-paper metabolite / sibling-drug suffix additions discovered during the 2026-05-28 naming-audit eta + residual-SD cleanup pass. Each is the canonical lowercase abbreviation used by its source paper for a non-parent species tracked alongside the parent.

### 1ohm (**canonical 1'-hydroxymidazolam suffix**)
- **Type:** metabolite-suffix
- **Role:** 1'-hydroxymidazolam metabolite of midazolam.
- **Source aliases:** none.
- **Example models:** `Brussee_2018_midazolam_PBPK.R`, `Franken_2017_midazolam.R`.

### 4ohctx (**canonical 4-hydroxycyclophosphamide suffix**)
- **Type:** metabolite-suffix
- **Role:** 4-hydroxycyclophosphamide active metabolite of cyclophosphamide.
- **Source aliases:** none.
- **Example models:** `Campagne_2019_cyclophosphamide_mouse.R`.

### cepm (**canonical carboxyethylphosphoramide mustard suffix**)
- **Type:** metabolite-suffix
- **Role:** Carboxyethylphosphoramide-mustard metabolite of cyclophosphamide.
- **Source aliases:** none.
- **Example models:** `Campagne_2019_cyclophosphamide_mouse.R`.

### ftc (**canonical emtricitabine suffix**)
- **Type:** metabolite-suffix
- **Role:** Emtricitabine (FTC) sibling-drug suffix.
- **Source aliases:** none.
- **Example models:** `Chen_2016_tenofovir_emtricitabine.R`.

### tfvdp (**canonical tenofovir diphosphate suffix**)
- **Type:** metabolite-suffix
- **Role:** Tenofovir diphosphate active intracellular metabolite suffix.
- **Source aliases:** none.
- **Example models:** `Chen_2016_tenofovir_emtricitabine.R`.

### ftctp (**canonical emtricitabine triphosphate suffix**)
- **Type:** metabolite-suffix
- **Role:** Emtricitabine triphosphate active intracellular metabolite suffix.
- **Source aliases:** none.
- **Example models:** `Chen_2016_tenofovir_emtricitabine.R`.

### snk (**canonical S-norketamine suffix**)
- **Type:** metabolite-suffix
- **Role:** S-norketamine metabolite of S-ketamine.
- **Source aliases:** none.
- **Example models:** `Flint_2017_Sketamine.R`.

### acid (**canonical simvastatin acid suffix**)
- **Type:** metabolite-suffix
- **Role:** Simvastatin acid active metabolite of simvastatin.
- **Source aliases:** none.
- **Example models:** `Jin_2014_simvastatin.R`.

### act (**canonical ACT-333679 selexipag metabolite suffix**)
- **Type:** metabolite-suffix
- **Role:** ACT-333679 selexipag active metabolite suffix.
- **Source aliases:** none.
- **Example models:** `Krause_2017_selexipag.R`.

### rtv (**canonical ritonavir sibling-drug suffix**)
- **Type:** metabolite-suffix
- **Role:** Ritonavir sibling-drug / booster suffix.
- **Source aliases:** none.
- **Example models:** `Schipani_2013_atazanavir_ritonavir.R`.

### 9oh (**canonical 9-hydroxyrisperidone suffix**)
- **Type:** metabolite-suffix
- **Role:** 9-hydroxyrisperidone (paliperidone) active metabolite of risperidone.
- **Source aliases:** none.
- **Example models:** `Sherwin_2012_risperidone.R`.

### p88 (**canonical iloperidone metabolite P-88 (M1) suffix**)
- **Type:** metabolite-suffix
- **Role:** P-88 (also termed M1), the active iloperidone metabolite produced via cytosolic / microsomal enzymes (not CYP2D6); contributes to the iloperidone therapeutic profile via D2 / 5-HT2A receptor binding affinity comparable to the parent.
- **Source aliases:** none.
- **Example models:** `Pei_2016_iloperidone.R`.

### p95 (**canonical iloperidone metabolite P-95 (M2) suffix**)
- **Type:** metabolite-suffix
- **Role:** P-95 (also termed M2), the iloperidone metabolite produced by CYP2D6-mediated hydroxylation; pharmacologically inactive on the 5-HT2A receptor.
- **Source aliases:** none.
- **Example models:** `Pei_2016_iloperidone.R`.

### 5oh (**canonical 5-hydroxyomeprazole suffix**)
- **Type:** metabolite-suffix
- **Role:** 5-hydroxyomeprazole metabolite of omeprazole.
- **Source aliases:** none.
- **Example models:** `Zhao_2018_omeprazole.R`.

### sfn (**canonical omeprazole sulfone suffix**)
- **Type:** metabolite-suffix
- **Role:** Omeprazole sulfone metabolite suffix.
- **Source aliases:** none.
- **Example models:** `Zhao_2018_omeprazole.R`.

### d3og (**canonical dapagliflozin 3-O-glucuronide suffix**)
- **Type:** metabolite-suffix
- **Role:** Dapagliflozin 3-O-glucuronide metabolite suffix.
- **Source aliases:** none.
- **Example models:** `vanderWalt_2013_dapagliflozin.R`.

### su12662 (**canonical SU12662 sunitinib metabolite suffix**)
- **Type:** metabolite-suffix
- **Role:** SU12662 sunitinib active metabolite suffix.
- **Source aliases:** none.
- **Example models:** `Ait-Oudhia_2016_sunitinib.R`.

### dact (**canonical desacetylcefotaxime metabolite suffix**)
- **Type:** metabolite-suffix
- **Role:** Desacetylcefotaxime (DACT), the major active metabolite of cefotaxime (CTX) formed by hepatic deacetylation. Used as a metabolite suffix in parent + metabolite simultaneous popPK models where CTX and DACT are fitted jointly with the conventional 1:1 (CTX-equivalent) mass-balance assumption FDACT/CTX = 1.
- **Source aliases:**
  - `DACT` -- used in `Ahsman_2010_cefotaxime.R` (paper notation).
- **Example models:** `Ahsman_2010_cefotaxime.R`.

### tam (**canonical tamoxifen tracking-species suffix**)
- **Type:** metabolite-suffix
- **Role:** Tamoxifen tracking species suffix (documented as the ticagrelor-paired tracking species in Almquist 2016).
- **Source aliases:** none.
- **Example models:** `Almquist_2016_ticagrelor.R`.

### vact (**canonical vascular-active lesion-state suffix**)
- **Type:** metabolite-suffix
- **Role:** Vascular-active lesion-state suffix in Schindler 2017 imatinib tumor lesion-specific states.
- **Source aliases:** none.
- **Example models:** `Schindler_2017_imatinib.R`.

### vell (**canonical vascular-extracellular-lesion suffix**)
- **Type:** metabolite-suffix
- **Role:** Vascular-extracellular-lesion state suffix in Schindler 2017 imatinib lesion modeling.
- **Source aliases:** none.
- **Example models:** `Schindler_2017_imatinib.R`.

### dens (**canonical lesion-density state suffix**)
- **Type:** metabolite-suffix
- **Role:** Lesion-density state suffix in Schindler 2017 imatinib lesion modeling.
- **Source aliases:** none.
- **Example models:** `Schindler_2017_imatinib.R`.

### m1trans (**canonical M1-trans rolofylline metabolite suffix**)
- **Type:** metabolite-suffix
- **Role:** M1-trans active hydroxyl metabolite of rolofylline. CYP3A4-mediated hydroxylation of the parent adenosine A1 receptor antagonist generates a diastereomeric pair of M1 metabolites; M1-trans is the (trans-) stereoisomer tracked alongside the parent and the M1-cis stereoisomer in the Stroh 2013 simultaneous PK model. Drives `central_m1trans` / `peripheral1_m1trans` and the `propSd_m1trans` / `addSd_m1trans` residuals.
- **Source aliases:** none.
- **Example models:** `Stroh_2013_rolofylline.R` (doi:10.1208/s12248-012-9443-5).

### m1cis (**canonical M1-cis rolofylline metabolite suffix**)
- **Type:** metabolite-suffix
- **Role:** M1-cis active hydroxyl metabolite of rolofylline. Co-eluting (cis-) stereoisomer of the CYP3A4 hydroxylation pair; formed both directly from the parent (fraction FM) and via unidirectional stereochemical interconversion from M1-trans. Drives `central_m1cis` and the `propSd_m1cis` / `addSd_m1cis` residuals.
- **Source aliases:** none.
- **Example models:** `Stroh_2013_rolofylline.R` (doi:10.1208/s12248-012-9443-5).

---

## TB-treatment drug suffixes (combination antibiotic)

TB-treatment drug suffixes used in combination-antibiotic `central_<drug>` / `depot_<drug>` / `peripheral1_<drug>` PK subsystems. Each suffix is the canonical drug INN lowercase abbreviation.

### rif (**canonical rifampicin drug suffix**)
- **Type:** metabolite-suffix
- **Role:** Rifampicin drug suffix in combination TB models.
- **Source aliases:** none.
- **Example models:** `Chen_2017_TB_MTP_GPDI_mouse.R`, `Clewe_2018_TB_MTP_GPDI_in_vitro.R`.

### inh (**canonical isoniazid drug suffix**)
- **Type:** metabolite-suffix
- **Role:** Isoniazid drug suffix in combination TB models.
- **Source aliases:** none.
- **Example models:** `Chen_2017_TB_MTP_GPDI_mouse.R`.

### emb (**canonical ethambutol drug suffix**)
- **Type:** metabolite-suffix
- **Role:** Ethambutol drug suffix in combination TB models.
- **Source aliases:** none.
- **Example models:** `Chen_2017_TB_MTP_GPDI_mouse.R`.

---

## Antibiotic combination drug suffixes

Antibiotic combination-PK drug suffixes (linezolid, vancomycin, meropenem long form, colistin, daptomycin) for bare drug-state and combination-stratified PD models.

### lzd (**canonical linezolid drug suffix**)
- **Type:** metabolite-suffix
- **Role:** Linezolid drug suffix in combination antibiotic models.
- **Source aliases:** none.
- **Example models:** `Wicha_2017_linezolid_meropenem_vancomycin.R`.
- **Notes:** Same token as the bare `lzd` drug-state compartment; both Types co-exist for the same canonical name.

### van (**canonical vancomycin drug suffix**)
- **Type:** metabolite-suffix
- **Role:** Vancomycin drug suffix in combination antibiotic models.
- **Source aliases:** none.
- **Example models:** `Wicha_2017_linezolid_meropenem_vancomycin.R`.
- **Notes:** Same token as the bare `van` drug-state compartment.

### mero (**canonical meropenem long-form drug suffix**)
- **Type:** metabolite-suffix
- **Role:** Meropenem long-form drug suffix in combination antibiotic models.
- **Source aliases:** none.
- **Example models:** `Mohamed_2016_colistin_meropenem.R`.
- **Notes:** Same token as the bare `mero` drug-state compartment.

### col (**canonical colistin drug suffix**)
- **Type:** metabolite-suffix
- **Role:** Colistin drug suffix in combination antibiotic models (sibling of the active-polymyxin-from-CMS use case above; same canonical token).
- **Source aliases:** none.
- **Example models:** `Mohamed_2016_colistin_meropenem.R`.
- **Notes:** Same token as the bare `col` drug-state compartment.

### dap (**canonical daptomycin drug suffix**)
- **Type:** metabolite-suffix
- **Role:** Daptomycin drug suffix in combination antibiotic models.
- **Source aliases:** none.
- **Example models:** Combination antibiotic PD models with daptomycin.
- **Notes:** Same token as the bare `dap` drug-state compartment.

---

## Miscellaneous metabolite suffixes

### glu (**canonical paracetamol glucuronide suffix (template)**)
- **Type:** metabolite-suffix
- **Role:** Glucuronide suffix used by paracetamol PBPK template / placeholder extraction. Sibling of the Allegaert 2015 `gluc`.
- **Source aliases:** none.
- **Example models:** `NA_NA_paracetamol.R`.

### metab (**canonical generic-metabolite template suffix**)
- **Type:** metabolite-suffix
- **Role:** Generic metabolite suffix used by template / placeholder models that track an unnamed metabolite.
- **Source aliases:** none.
- **Example models:** `NA_NA_sunitinib.R`.
- **Notes:** Documented as "the active metabolite of the parent drug" without naming a specific INN.

### udca (**canonical ursodeoxycholic acid suffix**)
- **Type:** metabolite-suffix
- **Role:** Ursodeoxycholic acid (UDCA) suffix used by the Zuo 2016 enterohepatic-recycling UDCA PBPK (`stomach_udca` / `intestine_udca` / `portal_udca` / `liver_udca` / `biliary_udca` / `blood_udca`).
- **Source aliases:** none.
- **Example models:** `Zuo_2016_UDCA.R`.

### gudca (**canonical glycine-conjugated UDCA suffix**)
- **Type:** metabolite-suffix
- **Role:** Glycine-conjugated UDCA (GUDCA) conjugate-pool suffix. Sibling suffix of `udca`.
- **Source aliases:** none.
- **Example models:** `Zuo_2016_UDCA.R`.

### tudca (**canonical taurine-conjugated UDCA suffix**)
- **Type:** metabolite-suffix
- **Role:** Taurine-conjugated UDCA (TUDCA) conjugate-pool suffix. Sibling suffix of `udca`.
- **Source aliases:** none.
- **Example models:** `Zuo_2016_UDCA.R`.

### pza (**canonical pyrazinamide drug suffix**)
- **Type:** metabolite-suffix
- **Role:** Pyrazinamide (PZA), fourth first-line anti-TB drug, used in `depot_pza` / `central_pza` / `peripheral1_pza` PK subsystem.
- **Source aliases:** none.
- **Example models:** `Chen_2017_TB_MTP_GPDI_mouse.R`.

### oc (**canonical oseltamivir carboxylate suffix**)
- **Type:** metabolite-suffix
- **Role:** Oseltamivir carboxylate (OC), the active neuraminidase-inhibitor metabolite formed from the oseltamivir prodrug primarily via human carboxylesterase 1 (HCE1) in the liver.
- **Source aliases:** none.
- **Example models:** `Standing_2012_oseltamivir.R`.

### ppf (**canonical propofol active-metabolite suffix**)
- **Type:** metabolite-suffix
- **Role:** Propofol (2,6-diisopropylphenol), the active sedative-hypnotic metabolite liberated from the water-soluble prodrug fospropofol (GPI 15715, AQUAVAN) via systemic alkaline-phosphatase hydrolysis. Used as the metabolite suffix on `central_ppf` / `peripheral1_ppf` compartments and `Cc_ppf` observation in joint parent-prodrug + active-drug popPK models.
- **Source aliases:** `PR` (Gibiansky 2005 poster table headings: Vc_PR, K10_PR, K12_PR, K21_PR, K_GPI-PR).
- **Example models:** `Gibiansky_2005_fospropofol.R`.

---

## PBPK organ sub-compartment suffixes (Ayyar 2024)

PBPK organ sub-compartment suffixes used by Ayyar 2024 givosiran whole-organ extractions. The `<organ>_endo` / `<organ>_deep` / `<organ>_vas` shape lets each organ carry endosomal, sequestered-deep, and vascular pools alongside the existing `vp_<organ>` membrane-limited form.

### endo (**canonical endosomal-pool suffix**)
- **Type:** metabolite-suffix
- **Role:** Endosomal pool inside the named organ.
- **Source aliases:** none.
- **Example models:** `Ayyar_2024_givosiran.R`.

### deep (**canonical deep-bound pool suffix**)
- **Type:** metabolite-suffix
- **Role:** Deep-bound / sequestered pool inside the named organ.
- **Source aliases:** none.
- **Example models:** `Ayyar_2024_givosiran.R`.

### vas (**canonical vascular-pool suffix**)
- **Type:** metabolite-suffix
- **Role:** Vascular pool inside the named organ, alongside the existing `vp_<organ>` membrane-limited form.
- **Source aliases:** none.
- **Example models:** `Ayyar_2024_givosiran.R`.

---

## Stereoisomer / enantiomer suffixes (Jansson 2008 eflornithine)

L- and D-enantiomer suffixes for stereoselective popPK models that simultaneously track both isomers of a racemic drug. The `_l` and `_d` suffixes attach to canonical compartments (`depot_l`, `central_d`, `peripheral1_l`, etc.), canonical parameters (`lcl_l`, `lvc_d`, `ltmax_abs_l`, etc.), and observation outputs (`Cc_l`, `Cc_d`). The `_rac` suffix denotes the racemic-sum observation output `Cc_rac = Cc_l + Cc_d` and its residual SD; it has no associated state compartment.

### l (**canonical L-enantiomer suffix**)
- **Type:** metabolite-suffix
- **Role:** L-isomer (levorotatory / S-isomer) of a racemic drug, tracked alongside the D-isomer in a stereoselective popPK model.
- **Source aliases:** none.
- **Example models:** `Jansson_2008_eflornithine_rat.R`.
- **Notes:** Although enantiomers are not metabolites, the registered metabolite-suffix machinery handles the same `<canonical>_<token>` shape. Paired with `d` for the D-isomer and `rac` for the racemic sum. Future papers using R/S, +/-, or E/Z stereodescriptors should register separate suffixes if their notation differs.

### d (**canonical D-enantiomer suffix**)
- **Type:** metabolite-suffix
- **Role:** D-isomer (dextrorotatory / R-isomer) of a racemic drug, tracked alongside the L-isomer in a stereoselective popPK model.
- **Source aliases:** none.
- **Example models:** `Jansson_2008_eflornithine_rat.R`.

### rac (**canonical racemic-sum output suffix**)
- **Type:** metabolite-suffix
- **Role:** Racemic-sum observation output in a stereoselective popPK model: `Cc_rac = Cc_l + Cc_d` (and `propSd_Cc_rac` for the racemic-output residual SD).
- **Source aliases:** none.
- **Example models:** `Jansson_2008_eflornithine_rat.R`.
- **Notes:** Algebraic sum only; no `central_rac`, `peripheral1_rac`, or `depot_rac` compartment exists.

### bibf (**canonical BIBF 1202 nintedanib-metabolite suffix**)
- **Type:** metabolite-suffix
- **Role:** BIBF 1202, the main hydrolytic metabolite of nintedanib (BIBF 1120) formed by cleavage of the methyl ester. Used as the metabolite suffix in parent + metabolite simultaneous popPK models (compartments `depot_bibf`, `central_bibf`; parameters `lka_bibf`, `lvc_bibf`, `lcl_bibf`, `lfdepot_bibf`, `ltlag_bibf`; residual `expSd_bibf`). Founding example: `Schmid_2017_nintedanib.R`.
- **Source aliases:** none.
- **Example models:** `Schmid_2017_nintedanib.R`.
