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

### cumhaz (**canonical cumulative-hazard state**)
- **Type:** compartment
- **Role:** Cumulative-hazard state for time-to-event / dropout sub-models. Integrates instantaneous hazard so that `survival = exp(-cumhaz)`.
- **Source aliases:** none.
- **Example models:** `Girard_2012_pimasertib.R`.
- **Notes:** Source NONMEM idiom is `$MODEL COMP=(CUMHAZ)` with `DADT(<cumhaz>) = HAZARD`.

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
- **Role:** Brain cerebrospinal-fluid compartment for brain-PBPK models. Replaces the older `brain_ecf` for the cerebrospinal-fluid compartment per the 2026-05-28 naming audit.
- **Source aliases:**
  - `brain_ecf` -- deprecated; replaced by `brain_csf`.
- **Example models:** `Grimm_2023_gantenerumab.R`, `Grimm_2023_trontinemab.R`.

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
- **Example models:** `Ide_2009_pravastatin.R` and similar EHC models.

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
- **Example models:** `tgi_sat_*.R` family.
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
- **Role:** Endogenous plasma glucose used by glucose / lactate turnover sub-models with drug-stimulated production. State holds a concentration (mmol/L), mirroring the source paper's mass-balance parameterisation.
- **Source aliases:** none.
- **Example models:** `Oualha_2014_epinephrine.R`.

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

### hb (**canonical hemoglobin**)
- **Type:** compartment
- **Role:** Hemoglobin PD output.
- **Source aliases:** none.
- **Example models:** anemia / EPO PD models.

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
- **Example models:** `NA_NA_paracetamol.R` reference extraction.

### gip (**canonical glucose-dependent insulinotropic polypeptide output**)
- **Type:** compartment
- **Role:** Glucose-dependent insulinotropic polypeptide PD output.
- **Source aliases:** none.
- **Example models:** `NA_NA_paracetamol.R` reference extraction.

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

### total_IgG (**canonical total serum IgG compartment**)
- **Type:** compartment
- **Role:** Total serum IgG PD output.
- **Source aliases:** none.
- **Example models:** `Valenzuela_2025_nipocalimab.R`.

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

## Survival / dropout cumulative-hazard compartments

### cumHaz_os (**canonical overall-survival cumulative-hazard**)
- **Type:** compartment
- **Role:** Overall-survival cumulative-hazard state in oncology TTE sub-models.
- **Source aliases:** none.
- **Example models:** `Schindler_2016_sunitinib.R`.

### cumHaz_drop (**canonical dropout cumulative-hazard**)
- **Type:** compartment
- **Role:** Dropout cumulative-hazard state.
- **Source aliases:** none.
- **Example models:** `Schindler_2016_sunitinib.R`.

### cumHazard (**canonical generic cumulative-hazard**)
- **Type:** compartment
- **Role:** Generic cumulative-hazard state.
- **Source aliases:** none.
- **Example models:** `Zecchin_2016.R`.

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

### m2 (**canonical N-desmethyl-bedaquiline (M2) suffix**)
- **Type:** metabolite-suffix
- **Role:** N-desmethyl-bedaquiline (M2) metabolite of bedaquiline.
- **Source aliases:** none.
- **Example models:** `Svensson_2016_bedaquiline.R` (DDMODEL00000219).

### endx (**canonical endoxifen suffix**)
- **Type:** metabolite-suffix
- **Role:** Endoxifen (4-hydroxy-N-desmethyltamoxifen), major active metabolite of tamoxifen.
- **Source aliases:** none.
- **Example models:** `TerHeine_2014_tamoxifen.R`.

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

### apapg (**canonical APAP-glucuronide suffix (Cook notation)**)
- **Type:** metabolite-suffix
- **Role:** Paracetamol (APAP) phase-II glucuronide metabolite (Cook notation).
- **Source aliases:** none.
- **Example models:** `Cook_2016_paracetamol.R` (DDMODEL00000271).

### apaps (**canonical APAP-sulphate suffix (Cook notation)**)
- **Type:** metabolite-suffix
- **Role:** Paracetamol (APAP) phase-II sulphate metabolite (Cook notation).
- **Source aliases:** none.
- **Example models:** `Cook_2016_paracetamol.R` (DDMODEL00000271).

### col (**canonical colistin metabolite suffix**)
- **Type:** metabolite-suffix
- **Role:** Colistin, the active polymyxin generated in vivo by hydrolysis of the prodrug colistimethate sodium (CMS).
- **Source aliases:** none.
- **Example models:** `LeuppiTaegtmeyer_2019_CMS.R` (DDMODEL00000295).
- **Notes:** Same token as the bare `col` drug-state compartment; both Types co-exist for the same canonical name.

### dha (**canonical dihydroartemisinin suffix**)
- **Type:** metabolite-suffix
- **Role:** Dihydroartemisinin, active metabolite of artesunate.
- **Source aliases:** none.
- **Example models:** `Birgersson_2019_artesunate.R` (DDMODEL00000297).

### ohi (**canonical hydroxy-itraconazole suffix**)
- **Type:** metabolite-suffix
- **Role:** Hydroxy-itraconazole (OH-ITZ), major active metabolite of itraconazole produced by CYP3A4 hydroxylation.
- **Source aliases:** none.
- **Example models:** `Hennig_2006_itraconazole.R`, `Hennig_2007_itraconazole.R`.

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

### desrbn (**canonical 25-O-desacetyl rifabutin suffix**)
- **Type:** metabolite-suffix
- **Role:** 25-O-desacetyl rifabutin, primary active metabolite of rifabutin formed by arylacetamide deacetylase.
- **Source aliases:** none.
- **Example models:** `Hennig_2015_rifabutin.R` (doi:10.1128/AAC.01195-15).

### az5104 (**canonical AZ5104 osimertinib metabolite suffix**)
- **Type:** metabolite-suffix
- **Role:** AZ5104 (N-desmethyl osimertinib), active EGFR-inhibitor metabolite of osimertinib formed predominantly via CYP3A4/5.
- **Source aliases:** none.
- **Example models:** `Brown_2017_osimertinib.R` (doi:10.1111/bcp.13223).

### ndsel (**canonical N-desmethyl-selumetinib suffix**)
- **Type:** metabolite-suffix
- **Role:** N-desmethyl-selumetinib, active selumetinib metabolite (~3-5-fold more potent for MEK1 inhibition than parent), formed by oxidative N-demethylation.
- **Source aliases:** none.
- **Example models:** `Patel_2017_selumetinib.R` (doi:10.1002/psp4.12175).

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
- **Example models:** `Urien_2005_capecitabine.R`.

### asn1 (**canonical AS(N-1)3' siRNA truncated antisense suffix**)
- **Type:** metabolite-suffix
- **Role:** AS(N-1)3' truncated antisense strand of GalNAc-conjugated siRNAs (givosiran and similar). Treated as the active metabolite equipotent with the parent for RISC loading and target mRNA silencing.
- **Source aliases:** none.
- **Example models:** `Ayyar_2024_givosiran.R` (doi:10.1016/j.xphs.2023.10.026).

### mhd (**canonical 10-monohydroxy oxcarbazepine suffix**)
- **Type:** metabolite-suffix
- **Role:** 10-monohydroxy derivative (MHD, "10-hydroxy-carbazepine"), primary active metabolite of oxcarbazepine produced by cytosolic arylketone reductases.
- **Source aliases:** none.
- **Example models:** `Rodrigues_2017_oxcarbazepine.R` (doi:10.1111/bcp.13392).

### r (**canonical R-enantiomer suffix**)
- **Type:** metabolite-suffix
- **Role:** R-enantiomer suffix for enantiomer-resolved popPK models in which both enantiomers are followed in plasma but no interconversion is modelled.
- **Source aliases:** none.
- **Example models:** `Valitalo_2017_ketorolac.R` (doi:10.1111/bcp.13311).
- **Notes:** Treated as a "non-parent analyte" suffix; neither enantiomer is the parent.

### s (**canonical S-enantiomer suffix**)
- **Type:** metabolite-suffix
- **Role:** S-enantiomer suffix for enantiomer-resolved popPK models.
- **Source aliases:** none.
- **Example models:** `Valitalo_2017_ketorolac.R`.

### noxide (**canonical roflumilast N-oxide suffix**)
- **Type:** metabolite-suffix
- **Role:** Roflumilast N-oxide, active metabolite of roflumilast (contributes ~90% of total PDE4 inhibitory activity).
- **Source aliases:** none.
- **Example models:** `Lahu_2010_roflumilast.R` (doi:10.2165/11536600-000000000-00000).

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

---

## Cell-type suffixes (Friberg multi-cell-type chains)

Cell-type suffixes used with Friberg-style `circ_<celltype>` myelosuppression compartments and `precursor1_<celltype>` ... `precursorN_<celltype>` maturation chains for paired-output multi-cell models. Registered 2026-05-28 per the naming audit.

### anc (**canonical ANC cell-type suffix**)
- **Type:** metabolite-suffix
- **Role:** Absolute neutrophil count cell-type suffix.
- **Source aliases:** none.
- **Example models:** `Han_2015_decitabine.R` and similar.

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

## Zurlinden 2016 paracetamol PBPK metabolite suffixes

Zurlinden 2016 paracetamol PBPK metabolite shorthand. Registered as separate suffixes (alongside the Cook 2016 `apapg` / `apaps`) so the Zurlinden `a_<organ>_as` / `a_<organ>_ag` compartment names pass the metabolite-suffix check without rewriting the source-paper notation.

### as (**canonical APAP-sulfate suffix (Zurlinden notation)**)
- **Type:** metabolite-suffix
- **Role:** APAP-sulfate metabolite suffix (Zurlinden notation; same chemical species as the Cook 2016 `apaps`).
- **Source aliases:** none.
- **Example models:** `Zurlinden_2016_paracetamol.R`.

### ag (**canonical APAP-glucuronide suffix (Zurlinden notation)**)
- **Type:** metabolite-suffix
- **Role:** APAP-glucuronide metabolite suffix (Zurlinden notation; same chemical species as the Cook 2016 `apapg`).
- **Source aliases:** none.
- **Example models:** `Zurlinden_2016_paracetamol.R`.

---

## Per-paper additions (2026-05-28 naming audit)

Per-paper metabolite / sibling-drug suffix additions discovered during the 2026-05-28 naming-audit eta + residual-SD cleanup pass. Each is the canonical lowercase abbreviation used by its source paper for a non-parent species tracked alongside the parent.

### 1ohm (**canonical 1'-hydroxymidazolam suffix**)
- **Type:** metabolite-suffix
- **Role:** 1'-hydroxymidazolam metabolite of midazolam.
- **Source aliases:** none.
- **Example models:** `Brussee_2018_midazolam_PBPK.R`.

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
