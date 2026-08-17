# Canonical compartment and metabolite-suffix names

This file is the authoritative register of compartment / state names and metabolite-suffix tokens used in nlmixr2lib models. Every compartment that appears in a model's `model()` block (and every metabolite / sibling-drug suffix used to derive a compartment / parameter / residual-SD name) is expected to match one of the canonical entries below -- or to fit one of the regex constants documented in the header section. The register is seeded from the 2026-05-28 naming audit and extended whenever a new paper introduces a state that isn't yet registered.

## How to use this register

1. **Before adding a compartment / state to a new model**, search this file (by canonical name and by source alias) for the concept you need.
2. **If the canonical name exists**, use it exactly. Document the source-paper rename in a code comment if the paper used a different name.
3. **If the source paper uses an alias listed under an existing canonical name**, prefer the canonical name. Aliases are documented for cross-reference, not as a free pass to introduce the deprecated form in new models.
4. **If the state is not in this register at all**, propose a new entry with a canonical name, type, role, source aliases, and example models. Verify with the user before committing. The addition is part of the model's PR.
5. **Do not modify existing model files when you discover a missing entry**; simply register the canonical here. Retrofitting existing models is a separate effort.
6. **Never add a second entry for a name that already has one at the same `Type`.** Extend the existing block instead -- add your source alias and example model to it. A name may appear twice only when the two entries carry *different* `Type` values (`col`, `complex`, `dap`, `lzd`, `mer`, `mero`, `plasma` and `van` are each both a bare `compartment` and a `metabolite-suffix`). This register is resolved in document order, last one wins, so a same-`Type` repeat silently discards the earlier block along with any alias or example recorded only there. `buildModelDb()` fails the build on one.

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

- `compartmentRegex = "^(transit|effect|precursor|lat|depot)[0-9]+$"` -- numbered-chain compartments: transit absorption chains (`transit1`, `transit2`, ...), effect-compartment chains (`effect1`, `effect2`, ...), precursor pools for delayed-feedback IDR (`precursor1`, `precursor2`, ...), latent chains (`lat1`, ...), parallel-absorption depots (`depot1`, `depot2`, ...). Numeric suffix is required (single-state members use the bare canonical `effect` / `depot`).
- `darCompartmentRegex = "^dar[0-9]+_(central|peripheral[0-9]?)$"` -- DAR-numbered ADC isoform compartments (`dar0_central`, `dar4_peripheral1`, ...).
- `targetLocationRegex = "^(target|complex)_(csf|isf|peripheral[0-9]?)$"` -- target species in physiologic / numbered-peripheral compartments (`target_csf`, `target_isf`, `target_peripheral`, `target_peripheral1`, `complex_peripheral`, ...).
- `pbpkSubCompartmentRegex = "^(bc|eu|eb|fr|is|int|mrna|luc)_(liver|lung|kidney|spleen|heart|muscle|skin|adipose|bone|brain|small_intestine|large_intestine|pancreas|thymus|portal|remainder|other|hepatic|fat|rapidly_perfused|slowly_perfused|venous|arterial|urine|gut|tumor)$` -- membrane-limited PBPK sub-compartments: vascular blood cells (`bc_`), endosomal unbound (`eu_`), endosomal FcRn-bound (`eb_`), endosomal free FcRn (`fr_`), interstitial space (`is_`), intracellular (`int_`), mRNA pool (`mrna_`), luciferase reporter (`luc_`).
- `rbcCompartmentRegex = "^rbc_[a-z0-9]+$"` -- intracellular drug / active-metabolite pools inside red blood cells, carried as ODE states in concentration units (`rbc_mtx`, `rbc_tgn`). Deliberately kept out of `registeredMetabolites` because the analyte is frequently the *parent* drug (methotrexate), and recording a parent drug in the metabolite register would mislead later readers of that list. See the "Intracellular red-cell analyte pools" section below for the naming rule and the per-analyte entries.
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

### brain_icf (**canonical brain intracellular-fluid compartment**)
- **Type:** compartment
- **Role:** Brain intracellular fluid (the aggregate cytosolic space of all brain cells) in whole-CNS PBPK models that resolve intracellular from extracellular brain drug distribution. Drug reaches it only through the brain cell membrane (`brain_cell_membrane`), never directly from `brain_ecf`, and exchanges onward with the acidic lysosomal space (`brain_lysosome`). Carried as a total unbound concentration of which the pH-dependent fraction `PHF_ICF` is uncharged; the mouse volume is 80 percent of total brain volume (288 uL).
- **Source aliases:** `brainICF`, `brain_ICF`, `ICF`.
- **Example models:** `Saleh_2023_quinidine_mouse_pbpk.R` (founding example; the LeiCNS-PK3.0 mouse family, 10 drugs).
- **Notes:** Distinct from the membrane-limited PBPK sub-compartment pattern `int_brain` (`pbpkSubCompartmentRegex`), which denotes a generic intracellular organ sub-compartment in whole-body PBPK models. `brain_icf` is the brain-specific role under the `brain_<region>` namespace and is paired with `brain_ecf`, `brain_cell_membrane` and `brain_lysosome` as one mechanistic intra-brain unit.

### brain_cell_membrane (**canonical brain cell-membrane phospholipid compartment**)
- **Type:** compartment
- **Role:** Non-specific phospholipid binding compartment representing the brain cell membranes, introduced in LeiCNS-PK3.0 to replace the instantaneous "binding factor" of LeiCNS-PK1.0. It sits between `brain_ecf` and `brain_icf` and exchanges with BOTH faces using the same water-to-oil / oil-to-water clearance pair (`CL_wo` = `P0` x `SA_BCM`, `CL_ow` = `CL_wo` / 10^logP), so at equilibrium the membrane-to-ECF concentration ratio equals the octanol/water partition coefficient. Its volume is the brain phospholipid volume fraction (0.05) times total brain volume.
- **Source aliases:** `BCM`, `brain_BCM`, `brainCellMembrane`, `phospholipid`.
- **Example models:** `Saleh_2023_quinidine_mouse_pbpk.R` (founding example; the LeiCNS-PK3.0 mouse family, 10 drugs).
- **Notes:** Named for the anatomical structure rather than the binding substance so it composes with the `brain_<region>` namespace. Because equilibration across the membrane is fast (half-life well under a minute for lipophilic drugs), this state is numerically near-instantaneous; it is nonetheless carried explicitly because it is what makes the intra-brain steady-state identity `P_oct/water = C_BCM / (C_ECF x PHF_ECF)` hold.

### brain_lysosome (**canonical brain lysosomal compartment**)
- **Type:** compartment
- **Role:** Aggregate acidic lysosomal space of the brain cells (pH 5.5) in whole-CNS PBPK models that account for lysosomal trapping of basic drugs. Exchanges only with `brain_icf` across the lysosomal membrane at `Q_LYSO` = `P0` x `SA_LYSO`, with the neutral fraction on each side (`PHF_ICF`, `PHF_LYS`) gating the transfer; the large pH gradient is what drives accumulation of bases. Mouse volume is 1.25 percent of the brain ICF volume (3.6 uL).
- **Source aliases:** `lysosome`, `LYSO`, `brain_LYS`.
- **Example models:** `Saleh_2023_quinidine_mouse_pbpk.R` (founding example; the LeiCNS-PK3.0 mouse family, 10 drugs).
- **Notes:** Prefixed under the `brain_<region>` namespace even though lysosomes are subcellular rather than anatomical regions, so that the whole intra-brain unit (`brain_ecf` / `brain_cell_membrane` / `brain_icf` / `brain_lysosome`) shares one prefix. A non-brain lysosomal compartment in a whole-body PBPK model would need its own registered name.

### brain_deep (**canonical deep brain compartment**)
- **Type:** compartment
- **Role:** Deep brain tissue compartment in brain-PBPK extractions.
- **Source aliases:** none.
- **Example models:** `Grimm_2023_gantenerumab.R`.

### brain_vascular (**canonical brain vascular compartment**)
- **Type:** compartment
- **Role:** Drug in cerebral capillary blood (volume `Vbv`, fed by cerebral blood flow `CLbv` from systemic central) in hybrid physiology-based PK-PD models that resolve the blood-brain barrier transport as two coupled states. In whole-CNS PBPK models this is also the brain microvasculature: the intravascular space that exchanges with systemic plasma at cerebral blood flow `Q_CBF` and is the donor compartment for BOTH the blood-brain barrier and the blood-CSF barrier.
- **Source aliases:** `brain_MV`, `brainMV`, `MV`, `microvasculature`.
- **Example models:** `Johnson_2011_olanzapine_rat.R`; `Saleh_2023_quinidine_mouse_pbpk.R` (LeiCNS-PK3.0 mouse family, where the compartment is the paper's `brainMV`, volume 5 uL).
- **Notes:** The LeiCNS-PK3.0 family's `brainMV` is this role, not a new one -- it holds the total drug concentration in cerebral capillary blood and the unbound driving concentration is `fu_plasma` times the total. A separate `brain_mv` canonical was proposed during extraction of the LeiCNS mouse family and rejected as a synonym of this entry.

### brain_extravascular (**canonical brain extravascular compartment**)
- **Type:** compartment
- **Role:** Drug in brain tissue beyond the BBB (volume `Vbev`, fed via the BBB-clearance term `CLbev`). Paired with `brain_vascular`; transport between the two states is driven by the unbound concentration on each side via the fixed `fu_plasma` and `fu_brain` fractions.
- **Source aliases:** none.
- **Example models:** `Johnson_2011_olanzapine_rat.R`.

---

## Lactation / mother-to-infant dyad namespace

Adopted 2026-08-06 with the first lactation model in the library (`Wattanakul_2024_primaquine.R` / `Wattanakul_2024_primaquine_motherinfant.R`; operator sidecar `oare_PMC11078975` request-001 / response-001, question q1, option A).

A lactation-transfer model has two structurally distinct sides. The maternal side is an ordinary popPK model plus one breast-milk compartment per analyte; that compartment is a physiologic matrix like `csf` or `isf`, so it gets the bare canonical `milk` (metabolite forms `milk_<metab>`). The infant side is a *second person's* complete PK model carried inside the same model object, so bare compartment names would collide with the mother's: `central` would be ambiguous between the two members of the dyad. The `infant_<canonical>` prefix resolves this the same way `brain_<region>` resolves the bare-`cortex` collision -- every canonical compartment name may be prefixed with `infant_` to denote the breastfed-partner copy of that state, and metabolite suffixes compose normally (`infant_central_<metab>`).

The corresponding derived observation variables are `Cmilk` / `Cmilk_<metab>` and `Cinfant` / `Cinfant_<metab>`. Infant-partner covariates are `WT_INFANT` and `AGE_INFANT` (see `covariate-columns.md`); the milk-transfer parameters are `lq_milk`, `pcmilk`, `cfcap`, and `kmilkinf` (see `parameter-names.md`).

### milk (**canonical breast-milk compartment**)
- **Type:** compartment
- **Role:** Breast-milk physiologic compartment of a lactating subject. Holds the amount of the parent drug in the milk that the infant ingests at one feed; its volume is normally derived from the infant's daily milk intake divided by the number of feeds per day rather than estimated. Metabolite forms take the standard suffix (`milk_cpq`). Matches the `milk` entry already present in `conventions$specimenVocabulary`.
- **Source aliases:** none.
- **Example models:** `Wattanakul_2024_primaquine.R`, `Wattanakul_2024_primaquine_motherinfant.R` (breast-milk compartments for primaquine and carboxyprimaquine, exchanging with their respective central compartments through a shared apparent inter-compartmental clearance and an analyte-specific milk:plasma partition coefficient, gated by a square-wave breastfeeding function).

### infant_depot (**canonical breastfed-infant dose compartment**)
- **Type:** compartment
- **Role:** Dose compartment of the breastfed infant in a mother-to-infant dyad model. Receives drug from the mother's `milk` compartment during a feed rather than from an external dose event.
- **Source aliases:** none.
- **Example models:** `Wattanakul_2024_primaquine_motherinfant.R`.

### infant_transit1, infant_transit2, infant_transit3, infant_transit4 (**canonical breastfed-infant transit-absorption chain**)
- **Type:** compartment
- **Role:** Transit-absorption chain of the breastfed infant in a mother-to-infant dyad model. Same semantics as the bare `transit<n>` chain, applied to the infant partner.
- **Source aliases:** none.
- **Example models:** `Wattanakul_2024_primaquine_motherinfant.R` (two transit compartments per analyte, mean transit time fixed to a paediatric literature value).

### infant_central (**canonical breastfed-infant central compartment**)
- **Type:** compartment
- **Role:** Central (plasma) compartment of the breastfed infant in a mother-to-infant dyad model. Metabolite forms take the standard suffix (`infant_central_cpq`).
- **Source aliases:** none.
- **Example models:** `Wattanakul_2024_primaquine_motherinfant.R`.

### infant_peripheral1, infant_peripheral2 (**canonical breastfed-infant peripheral compartments**)
- **Type:** compartment
- **Role:** Peripheral distribution compartments of the breastfed infant in a mother-to-infant dyad model. Registered alongside the rest of the namespace so a future dyad model with a multi-compartment infant disposition does not have to re-open the naming question; not exercised by the founding models, whose infant disposition is one-compartment per analyte.
- **Source aliases:** none.
- **Example models:** none yet (namespace member registered with `Wattanakul_2024_primaquine_motherinfant.R`).

---

## Friberg myelosuppression chains

### circ (**canonical circulating-cell compartment**)
- **Type:** compartment
- **Role:** Friberg-style myelosuppression circulating-cell compartment (Friberg 2002 paclitaxel and derivatives). Terminal compartment of a `precursor1` ... `precursorN` maturation chain.
- **Source aliases:** none.
- **Example models:** Friberg-style myelosuppression PD extractions.
- **Notes:** Replaces a paper-naming `central` for circulating neutrophils / platelets / lymphocytes when the model is a maturation chain rather than a classical-PK central compartment. Suffix-form `circ_<celltype>` (e.g., `circ_anc`, `circ_plt`) is accepted for paired-output multi-cell-type models via the registered `_anc` / `_plt` / `_wbc` metabolite suffixes.

### prol (**canonical Friberg proliferating-cell pool**)
- **Type:** compartment
- **Role:** Self-renewing proliferating-progenitor pool at the head of a Friberg-Karlsson myelosuppression chain (`prol` -> `transit1` ... `transitN` -> `circ`). The drug effect and the `(circ0 / circ)^gamma` rebound-feedback term act on the proliferation rate of this state.
- **Source aliases:** none.
- **Example models:** `Hansson_2013_sunitinib_myelosuppression.R`.
- **Notes:** Registered 2026-06-28. Names the proliferating head state explicitly when the source model separates proliferation from a numbered `transit1..N` maturation chain feeding the canonical `circ` circulating pool (some library models instead use the `precursor1..N` chain form documented under `circ`). Distinct from `cycling_cells` (the Simeoni 2004 oncology TGI proliferating pool); `prol` is the hematopoietic stem/progenitor pool of the neutrophil / platelet / leukocyte maturation cascade.

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
- **Notes:** `tumor` is also a registered organ of the `pbpkSubCompartmentRegex` header pattern, so a permeability-limited tumour / target-tissue PBPK model carries drug in `is_tumor` (the vascular + interstitial, i.e. extracellular, space) and `int_tumor` (the intracellular space) -- exactly as `is_liver` / `int_liver` do for the liver. Source aliases for those two sub-compartments include `C_tumourEx` / `C_tumourIC` (`Aoki_2024_intratarget_microdosing_pbpk.R`, whose Table 1 calls the same two spaces "volume of target" and "volume of inner-cellular space of the target"). Note the distinction from the bare `tumor` state above, which holds a tumour *size* in a TGI model rather than a drug concentration, and from `ecf`, the single lumped brain / tumour extracellular-fluid space used by microdialysis models that do not resolve an intracellular compartment.

### tumor_size, TS (**canonical TGI tumor-size output state**)
- **Type:** compartment
- **Role:** Snake-case canonical output-state name for the TGI template family (`tumor_size`); `TS` is the equivalent upper-case abbreviation (RECIST 1.1 sum of longest diameters of target lesions, mm) registered as a canonical sibling because Struemper 2025 uses `TS` directly as the model observation variable.
- **Source aliases:**
  - `Ts` -- deprecated legacy lower-case form.
  - `ts` -- deprecated legacy lower-case form.
- **Example models:** `tgi_no_sat_*.R`, `tgi_sat_*.R`, `Ouerdani_2015_pazopanib.R`, `Mazzocco_2015.R`, `Zecchin_2016.R`, `Wilson_2015_sunitinib_irinotecan_mouse.R`, `Struemper_2025_tumorsize_OS_nsclc.R` (as `TS`).
- **Notes:** Registered 2026-05-28 per the naming audit for the TGI template family. `TS` added 2026-06-28 as a canonical sibling name (upper-case RECIST sum-of-longest-diameters abbreviation) for Struemper 2025, which observes `TS = growth + shrink - TSb` from the Stein bi-exponential `growth` / `shrink` states; the related time-varying covariate column for an observed-TS data input is `TUM_SLD` in `covariate-columns.md`. The deprecated lower-case `Ts` / `ts` forms remain aliases (no active model uses them as the bare observation).

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

### growth (**canonical Stein-model tumor-growth sub-state**)
- **Type:** compartment
- **Role:** Exponentially growing (treatment-resistant) sub-population of the Stein bi-exponential tumor-size model: `d/dt(growth) = kge * growth` with `growth(0) = TSb`, so `growth(t) = TSb * exp(kge * t)`. Paired with `shrink`; the observed tumor size is `TS = growth + shrink - TSb`.
- **Source aliases:** none.
- **Example models:** `Struemper_2025_tumorsize_OS_nsclc.R`.
- **Notes:** Registered 2026-06-28. The Stein (2008) bi-exponential TS model decomposes tumor dynamics into a treatment-resistant growing fraction (`growth`, rate `kge`) and a treatment-sensitive shrinking fraction (`shrink`, rate `kse`). Distinct from the saturable-growth `tumor_size` / `carrying_capacity` family and from the Simeoni `cycling_cells` chain.

### shrink (**canonical Stein-model tumor-shrinkage sub-state**)
- **Type:** compartment
- **Role:** Exponentially shrinking (treatment-sensitive) sub-population of the Stein bi-exponential tumor-size model: `d/dt(shrink) = -kse * shrink` with `shrink(0) = TSb`, so `shrink(t) = TSb * exp(-kse * t)`. Paired with `growth`; the observed tumor size is `TS = growth + shrink - TSb`.
- **Source aliases:** none.
- **Example models:** `Struemper_2025_tumorsize_OS_nsclc.R`.
- **Notes:** Registered 2026-06-28. See `growth` for the Stein bi-exponential decomposition.

### growth_ctdna (**canonical Stein-model ctDNA-growth sub-state**)
- **Type:** compartment
- **Role:** Circulating-tumor-DNA counterpart of the `growth` Stein sub-state: `d/dt(growth_ctdna) = kge_ctdna * growth_ctdna` with `growth_ctdna(0) = rbase_ctdna`, so `growth_ctdna(t) = rbase_ctdna * exp(kge_ctdna * t)`. Paired with `shrink_ctdna`; the observed ctDNA level is `ctdna = growth_ctdna + shrink_ctdna - rbase_ctdna`. The `_ctdna` endpoint suffix disambiguates the ctDNA Stein pair from the tumor-size pair (`growth` / `shrink`) when a single model fits both endpoints jointly.
- **Source aliases:** none.
- **Example models:** `Ribba_2022_ctdna.R`, `Ribba_2022_ctdna_sld_joint.R`.
- **Notes:** Registered 2026-07-28 alongside the Ribba 2022 ctDNA extraction, the library's first ctDNA-modality model. `rbase_ctdna` is the model-side baseline, which in Ribba 2022 is `log10(CTDNA)` because the source paper fits the Stein form to base-10 log-transformed MMPM; a model fitting ctDNA on the natural scale would set `rbase_ctdna <- CTDNA` instead. See `growth` for the underlying Stein (2011) bi-exponential decomposition.

### shrink_ctdna (**canonical Stein-model ctDNA-decay sub-state**)
- **Type:** compartment
- **Role:** Circulating-tumor-DNA counterpart of the `shrink` Stein sub-state: `d/dt(shrink_ctdna) = -kse_ctdna * shrink_ctdna` with `shrink_ctdna(0) = rbase_ctdna`, so `shrink_ctdna(t) = rbase_ctdna * exp(-kse_ctdna * t)`. Paired with `growth_ctdna`.
- **Source aliases:** none.
- **Example models:** `Ribba_2022_ctdna.R`, `Ribba_2022_ctdna_sld_joint.R`.
- **Notes:** Registered 2026-07-28. In a joint ctDNA / tumor-size model the ctDNA decay rate may be tied to the tumor-size decay rate rather than estimated freely -- Ribba 2022 sets `kse_ctdna = zeta * kse` (see `zeta` in `parameter-names.md`). See `growth_ctdna` for the paired state.

### ctdna (**canonical ctDNA output state**)
- **Type:** compartment
- **Role:** Circulating-tumor-DNA observation variable for liquid-biopsy biomarker models, analogous to `TS` for RECIST tumor size. In Stein bi-exponential ctDNA models it is the algebraic combination `ctdna = growth_ctdna + shrink_ctdna - rbase_ctdna`.
- **Source aliases:**
  - `ctDNA` -- the source-paper mixed-case spelling; the canonical form is all-lower-case for consistency with the rest of the compartment register.
- **Example models:** `Ribba_2022_ctdna.R`, `Ribba_2022_ctdna_sld_joint.R`.
- **Notes:** Registered 2026-07-28. The scale on which `ctdna` is expressed is model-specific and must be recorded in the model's `units` metadata: Ribba 2022 works on log10(MMPM) because the source paper log-transformed the data before fitting, so the residual-error parameter `addSd` is in log10 units, not MMPM. The related baseline covariate column is `CTDNA` (untransformed MMPM) in `covariate-columns.md`.

---

## CD3 bispecific / trispecific binding states (tumor microenvironment)

The three canonicals below describe the drug-target-effector-cell mass-action binding equilibrium at the tumor microenvironment used by CD3 (and analogous immune-effector) bispecific / trispecific QSP models. Each compartment holds an accessible-void concentration (nM) in the tumor extracellular space; the "free" (unbound) target concentrations enter the mass-action rate law after being divided by the tumor void fraction `void_frac`. Companion parameters `kon_cd3`, `koff_cd3`, `kon_pcad`, `koff_pcad` and receptor densities `cd3_receptors`, `mpcad`, `tumor_cells_g` live in `parameter-names.md`.

### drug_cd3_tumor (**canonical drug-CD3 dimer in tumor**)
- **Type:** compartment
- **Role:** Drug-CD3 receptor dimer concentration in the tumor extracellular space (nM). Formed by mass-action binding of free drug to CD3 receptors on T cells present in the TME; can bind additional free P-cadherin (or the target antigen for other CD3 bispecifics) to yield the productive `trimer` complex, or dissociate back to free drug and CD3.
- **Source aliases:**
  - `DCD3t` -- Betts 2019 paper notation.
- **Example models:** `Betts_2019_pf_06671008_qsp.R`.
- **Notes:** Founding example Betts 2019 (Eq 14 dDCD3t/dt). CD3 concentration is derived from a T-cell density (cells/L) and per-cell receptor count, divided by Avogadro; see `parameter-names.md` § "Receptor densities". Sibling `drug_pcad_tumor` and `trimer` complete the CD3 x TAA binding triad.

### drug_pcad_tumor (**canonical drug-target-antigen dimer in tumor**)
- **Type:** compartment
- **Role:** Drug bound to the tumor-associated antigen (P-cadherin, HER2, CEA, etc.) receptor in the tumor extracellular space (nM). Formed by mass-action binding of free drug to antigen receptors on tumor cells; can bind additional free CD3 (or other CD3 arm) to yield the productive `trimer` complex, or dissociate back to free drug and antigen. Includes an internalization loss term (`kint`) representing receptor-mediated endocytosis of the drug-antigen complex.
- **Source aliases:**
  - `DPcadt` -- Betts 2019 paper notation for the drug-P-cadherin dimer in tumor.
- **Example models:** `Betts_2019_pf_06671008_qsp.R`.
- **Notes:** Founding example Betts 2019 (Eq 15 dDPcadt/dt). Canonical name is `drug_pcad_tumor` for the founding P-cadherin example; per-antigen suffixed forms (`drug_her2_tumor`, `drug_cea_tumor`, etc.) can be registered as new canonicals when future CD3 bispecifics targeting other antigens are extracted. The internalization rate `kint` distinguishes this compartment from `drug_cd3_tumor`, which is not internalized in the current model class.

### trimer (**canonical drug-CD3-antigen ternary complex in tumor**)
- **Type:** compartment
- **Role:** Productive drug-CD3-antigen ternary complex (trimer) in the tumor extracellular space (nM). Forms an immune synapse-mimetic bridge between a T cell and an antigen-expressing tumor cell; the paper's PD driver linking bispecific PK to tumor cell killing via `kkill = kmax * trimer / (kc50 + trimer)`. Reversible mass-action formation from either dimer (`drug_cd3_tumor` + free antigen, or `drug_pcad_tumor` + free CD3) with dissociation back to either dimer.
- **Source aliases:**
  - `Trimer` -- Betts 2019 paper notation (deprecated capitalized form).
- **Example models:** `Betts_2019_pf_06671008_qsp.R`.
- **Notes:** Founding example Betts 2019 (Eq 16 dTrimer/dt). The trimer concentration drives the tumor-killing rate `kkill` and is the mechanistic PD linker of the CD3-bispecific class. Trimer formation exhibits the bell-shaped concentration-response phenomenon (Betts 2019 Fig 1b): trimer decreases at very high drug concentrations because the equilibrium shifts toward drug-CD3 and drug-antigen dimers rather than the productive trimer.

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

### enzyme_1a2 (**canonical CYP1A2 enzyme pool**)
- **Type:** compartment
- **Role:** CYP1A2 enzyme pool in an isoenzyme-resolved enzyme-turnover model; the relative enzyme amount is the isoenzyme's activity as a fraction of its untreated baseline. Same `enzyme_<isoform>` family as `enzyme_2b6` / `enzyme_2a6`, but driven by an exogenous or endogenous modulator rather than by drug autoinduction.
- **Source aliases:** none.
- **Example models:** `Willemin_2024_interleukin6_cyp_pbpk.R`.
- **Notes:** Initial condition `enzyme_1a2(0) <- 1` (relative to baseline). Founding example: `Willemin_2024_interleukin6_cyp_pbpk.R` (interleukin-6 net *induction* of CYP1A2 during cytokine release syndrome, alongside suppression of four other isoenzymes).

### enzyme_2c9 (**canonical CYP2C9 enzyme pool**)
- **Type:** compartment
- **Role:** CYP2C9 enzyme pool in an isoenzyme-resolved enzyme-turnover model; the relative enzyme amount is the isoenzyme's activity as a fraction of its untreated baseline.
- **Source aliases:** none.
- **Example models:** `Willemin_2024_interleukin6_cyp_pbpk.R`.
- **Notes:** Initial condition `enzyme_2c9(0) <- 1` (relative to baseline). Founding example: `Willemin_2024_interleukin6_cyp_pbpk.R`.

### enzyme_2c19 (**canonical CYP2C19 enzyme pool**)
- **Type:** compartment
- **Role:** CYP2C19 enzyme pool in an isoenzyme-resolved enzyme-turnover model; the relative enzyme amount is the isoenzyme's activity as a fraction of its untreated baseline.
- **Source aliases:** none.
- **Example models:** `Willemin_2024_interleukin6_cyp_pbpk.R`.
- **Notes:** Initial condition `enzyme_2c19(0) <- 1` (relative to baseline). Founding example: `Willemin_2024_interleukin6_cyp_pbpk.R`.

### enzyme_3a4 (**canonical CYP3A4 enzyme pool**)
- **Type:** compartment
- **Role:** CYP3A4 enzyme pool in an isoenzyme-resolved enzyme-turnover model; the relative enzyme amount is the isoenzyme's activity as a fraction of its untreated baseline.
- **Source aliases:** none.
- **Example models:** `Willemin_2024_interleukin6_cyp_pbpk.R`.
- **Notes:** Initial condition `enzyme_3a4(0) <- 1` (relative to baseline). Founding example: `Willemin_2024_interleukin6_cyp_pbpk.R`. Use bare `enzyme_3a4` when the model resolves several isoenzymes within a single tissue; when the *same* isoenzyme is resolved across several organs, append the organ to the isoenzyme rather than replacing it -- `enzyme_3a4_liver` / `enzyme_3a4_gut`. An organ-only form (`enzyme_liver`) is **not** used: it does not say which isoenzyme, and it would collide the moment a second isoenzyme were resolved in the same organ. Ratified 2026-08-06 with `Chen_2024_interleukin6_cyp3a_pbpk.R`.

### enzyme_3a5 (**canonical CYP3A5 enzyme pool**)
- **Type:** compartment
- **Role:** CYP3A5 enzyme pool in an isoenzyme-resolved enzyme-turnover model; the relative enzyme amount is the isoenzyme's activity as a fraction of its untreated baseline.
- **Source aliases:** none.
- **Example models:** `Willemin_2024_interleukin6_cyp_pbpk.R`.
- **Notes:** Initial condition `enzyme_3a5(0) <- 1` (relative to baseline). Founding example: `Willemin_2024_interleukin6_cyp_pbpk.R`. Kept separate from `enzyme_3a4` even when a source assumes identical potencies for the two, because the two isoenzymes drive different victim drugs and a user may want to break the assumption.

### enzyme_3a4_liver (**canonical hepatic CYP3A4 enzyme pool**)
- **Type:** compartment
- **Role:** Hepatic CYP3A4 enzyme pool in a tissue-resolved enzyme-turnover model, carrying an *absolute* enzyme abundance (pmol/mg microsomal protein) rather than a fraction of baseline. Turns over as `d/dt(enzyme_3a4_liver) <- kdeg_liver * bl_enzyme_3a4_liver * fsupp - kdeg_liver * enzyme_3a4_liver`, where `fsupp` is the cytokine-driven suppression of synthesis.
- **Source aliases:**
  - `CYP3A4 in liver` -- Chen 2024 Table 1 notation.
- **Example models:** `Chen_2024_interleukin6_cyp3a_pbpk.R` (founding example; baseline 137 pmol/mg microsomal protein in healthy volunteers, 82.2 in rheumatoid arthritis).
- **Notes:** Registered 2026-08-06. Initial condition is the baseline abundance `enzyme_3a4_liver(0) <- bl_enzyme_3a4_liver`, **not** 1 -- this is the absolute-abundance form, distinct from the relative-to-baseline `enzyme_3a4`. The isoenzyme is named before the organ so the two naming axes compose: use `enzyme_3a4` when several isoenzymes are resolved in one tissue, and `enzyme_3a4_<organ>` when one isoenzyme is resolved across organs. An organ-only `enzyme_liver` is deliberately not used -- it does not state the isoenzyme. Paired with `enzyme_3a4_gut`.

### enzyme_3a4_gut (**canonical intestinal CYP3A4 enzyme pool**)
- **Type:** compartment
- **Role:** Intestinal (small-bowel) CYP3A4 enzyme pool in a tissue-resolved enzyme-turnover model, carrying an *absolute* enzyme abundance (nmol per small intestine) rather than a fraction of baseline. Turns over as `d/dt(enzyme_3a4_gut) <- kdeg_gut * bl_enzyme_3a4_gut * fsupp - kdeg_gut * enzyme_3a4_gut`.
- **Source aliases:**
  - `CYP3A4 in gut` -- Chen 2024 Table 1 notation.
- **Example models:** `Chen_2024_interleukin6_cyp3a_pbpk.R` (founding example; baseline 66.2 nmol/small intestine in healthy volunteers, 40.0 in rheumatoid arthritis).
- **Notes:** Registered 2026-08-06. Initial condition `enzyme_3a4_gut(0) <- bl_enzyme_3a4_gut`. Kept separate from `enzyme_3a4_liver` because the gut pool has both a different baseline unit (an amount per organ, not a per-mg-protein density) and a different degradation half-life, and because intestinal first-pass CYP3A4 is the arm that drives oral victim-drug interactions. See `enzyme_3a4_liver` for the naming rule.

---

## DAS28 disease-activity score

### das28 (**canonical DAS28 output compartment**)
- **Type:** compartment
- **Role:** DAS28 disease-activity score output compartment used by rheumatoid-arthritis PD models. Single PD output, paper-named.
- **Source aliases:** none.
- **Example models:** `Frey_2013_tocilizumab.R`, `Ma_2020_sarilumab_das28crp.R`.
- **Notes:** Registered 2026-05-28 per the naming audit.

### das28cfb (**canonical DAS28 change-from-baseline output compartment**)
- **Type:** compartment
- **Role:** DAS28 change-from-baseline (DAS28cfb) PD output used by rheumatoid-arthritis models that fit the paper-declared change score rather than the absolute DAS28 value. Companion to `das28` (which holds the absolute score); use `das28cfb` when the paper's own equation targets `DAS28cfb = f(t, Cij)`, with the change interpretation as a negative-going quantity for treatment improvement. Same canonical-lower-case-name convention as `das28` / `deltaUPDRS`.
- **Source aliases:** `DAS28cfb` -- Williams 2016 paper notation.
- **Example models:** `Williams_2016_rituximab_das28cfb.R`.
- **Notes:** Ratified canonically on 2026-07-09 alongside the Williams 2016 rituximab-biosimilar extraction (the paper fits the DAS28 change-from-baseline directly via a 1 - exp(fnon-C + fC) transformation, so the modelled endpoint is the change score rather than the absolute DAS28). Follows the same lowercase-name convention as `das28` / `deltaUPDRS`.

---

## Abuse-liability and opioid-withdrawal clinical scores

### druglikingvascfb (**canonical period-corrected drug liking VAS change-score output compartment**)
- **Type:** compartment
- **Role:** Period-corrected drug liking maximum-effect (Emax) visual analog scale (VAS) score output for opioid abuse-liability / blockade PD models. The endpoint is a CHANGE score: the pre-challenge VAS value is subtracted from the maximum value recorded after a fixed opioid challenge, so the state is a difference on the VAS scale rather than an absolute reading. Use `druglikingvascfb` whenever the paper's own equation targets the period-corrected (change) quantity; a future model fitting the ABSOLUTE drug liking VAS reading should register a companion `druglikingvas`, exactly as `das28` and `das28cfb` are separated.
- **Source aliases:** `VAS(Cp)`, `drug liking Emax VAS score` -- Walsh 2024 paper notation.
- **Example models:** `Walsh_2024_buprenorphine_drugLiking.R` (Imax blockade of hydromorphone-18-mg drug liking by subcutaneous depot buprenorphine CAM2038; logit-transformed onto -1 to 52).
- **Notes:** Ratified canonically on 2026-08-05 alongside the Walsh 2024 CAM2038 extraction. The `cfb` suffix follows the `das28` / `das28cfb` precedent for marking a change-from-baseline variant, and the lowercase run-together form follows `das28` / `deltaUPDRS`. Note the deliberate distinction from the metabolite-suffix token `vas` (canonical vascular-pool suffix, Ayyar 2024 PBPK) -- the full `druglikingvascfb` name avoids that collision, which is why a bare `vas` observation variable must never be used for a VAS score. Walsh 2024 bounds the logit transformation at -1 to 52 (not 0 to 100) because the period-corrected bipolar score is a change from a neutral 50, so 52 is just above the largest attainable increase.

### desiretousevas (**canonical desire-to-use VAS output compartment**)
- **Type:** compartment
- **Role:** Desire to use visual analog scale (VAS) score output for opioid craving / substance-use PD models. Absolute (not change-corrected) reading on a 0-100 mm unipolar scale where 0 = no effect. Companion endpoint to `druglikingvascfb` in opioid-blockade analyses.
- **Source aliases:** `desire to use VAS score` -- Walsh 2024 paper notation.
- **Example models:** `Walsh_2024_buprenorphine_desireToUse.R` (Imax suppression of pre-challenge desire to use by subcutaneous depot buprenorphine CAM2038, with an exponential onset-delay term on Imax; logit-transformed onto -1 to 101).
- **Notes:** Ratified canonically on 2026-08-05 alongside the Walsh 2024 CAM2038 extraction. No `cfb` suffix because the modelled quantity is the absolute pre-challenge score, in contrast to its sister endpoint `druglikingvascfb`. The logit bounds -1 to 101 pad the 0-100 unipolar scale by one unit either side so that predictions of exactly 0 and 100 remain attainable.

### cows (**canonical Clinical Opiate Withdrawal Scale total-score output compartment**)
- **Type:** compartment
- **Role:** Clinical Opiate Withdrawal Scale (COWS) total score output for opioid-withdrawal PD models. Eleven observer-rated withdrawal signs each scored 0-4, giving an integer total of 0-44 (45 ordered categories) in which 5-12 denotes mild symptoms. Because the endpoint is a bounded integer, models of this state typically live on the standard-normal quantile (probit) scale and map back via the normal CDF -- see `probitbase` / `probitbase_low` in `parameter-names.md`.
- **Source aliases:** `COWS score`, `COWS total score` -- Walsh 2024 paper notation.
- **Example models:** `Walsh_2024_buprenorphine_cows.R` (bounded-integer Imax model of COWS suppression by subcutaneous depot buprenorphine CAM2038, encoded as `probitNorm(addSd, 0, 45)`).
- **Notes:** Ratified canonically on 2026-08-05 alongside the Walsh 2024 CAM2038 extraction. Registered as the bare score acronym (no scale-type suffix) because COWS is a single unambiguous instrument with one total score, unlike the VAS family where the scale type must be named. The `0, 45` bounds on the probit transform are the 45 ordered categories, not the 0-44 score maximum.

---

## Depression-severity clinical scores

### madrsenh (**canonical MADRS enhancement-rate output compartment**)
- **Type:** compartment
- **Role:** Enhancement rate in depression severity: the percentage reduction in the Montgomery-Asberg Depression Rating Scale (MADRS) total score from the pre-treatment baseline, `100 * (MADRS_baseline - MADRS_t) / MADRS_baseline`. Used as the modelled endpoint by antidepressant PD models that fit the paper-declared relative improvement rather than the absolute score. **Sign convention: POSITIVE for clinical improvement**, because the underlying MADRS score falls as depression improves. This is the opposite orientation to `das28cfb`, which is a raw change score and therefore negative-going for improvement -- the two must not be treated as interchangeable shapes. A future model fitting the ABSOLUTE MADRS reading should register a companion `madrs`, and one fitting the raw (non-percentage) change should register `madrscfb`, exactly as `das28` and `das28cfb` are separated.
- **Source aliases:** `EFF`, `enhancement rate` -- Shigetome 2025 paper notation and NONMEM `$PRED` block.
- **Example models:** `Shigetome_2025_paroxetine_madrs.R` (Emax model in treatment duration, `EFF = Emax * Time / (ET50 + Time)`, with cumulative first-week paroxetine exposure on Emax and the week-1 MADRS score on ET50).
- **Notes:** Ratified canonically alongside the Shigetome 2025 paroxetine extraction, following the same lowercase run-together convention as `das28` / `das28cfb` / `deltaUPDRS` / `cows` / `druglikingvascfb`. The state is a percentage on a roughly 0-100 scale but is not bounded by the model: a negative value is attainable when the score worsens from baseline, and the Emax term itself can exceed 100 at high exposure, so this endpoint should not be given a bounded (logit / probit) transform unless the source paper declares one. Distinct from the covariate `SCORE_MADRS` in `covariate-columns.md`, which carries an absolute MADRS reading at a stated visit; the baseline reading that defines this endpoint's denominator is part of the endpoint, not a covariate.

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
- **Notes:** Founding example Plan 2012 PAGE poster; the BMD output is algebraic (BMD <- b * (1 + sum(slope * piece) / 100)) and is the LHS of the residual-error declaration `BMD ~ add(addSd)`. Distinct from the `vp_bone` PBPK bone-vascular compartment, which is a drug-distribution state, not an output endpoint. Anatomically-specific variants `BMD_LS` (lumbar spine) and `BMD_TH` (total hip) are registered as separate canonicals for models where the two DXA regions are distinct ODE states.

### BMD_LS (**canonical lumbar-spine bone mineral density**)
- **Type:** compartment
- **Role:** Lumbar-spine (L2-L4) bone mineral density ODE state, reported in g/cm^2. Used by mechanistic bone-remodelling models that treat the lumbar-spine DXA region as an ODE compartment distinct from the femoral-neck / total-hip region (Berkhout 2015 indirect-response BMD equation, driven by relative osteoblast (`osteoblast`) production coefficient D_AOB and relative osteoclast (`osteoclast`) degradation coefficient D_AOC). The state's initial condition is a BMI-corrected baseline `BMD_LS_0 * (1 + BMI_frac_LS * (BMI - BMI_ref))` with exponential IIV.
- **Source aliases:** `LS-BMD`, `LSBMD`, `B_LSBMD` (Berkhout 2015 NONMEM compartment name).
- **Example models:** `Berkhout_2015_osteoporosis_placebo_qsp.R`.
- **Notes:** Capitalisation matches the parent `BMD` canonical. Founding example Berkhout 2015 (doi:10.1002/psp4.12006). Distinct from the general `BMD` canonical (algebraic in Plan 2012; ODE state here) by anatomical region.

### BMD_TH (**canonical total-hip bone mineral density**)
- **Type:** compartment
- **Role:** Total-hip (femoral neck + trochanter + intertrochanteric area) bone mineral density ODE state, reported in g/cm^2. Sibling of `BMD_LS` for models with two co-fit DXA regions (Berkhout 2015 indirect-response BMD equation). Same functional form as `BMD_LS` with independent baseline (`BMD_TH_0`), BMI covariate coefficient (`BMI_frac_TH`), and turnover rate (`k_in_TH`); the two states share the D_AOB / D_AOC osteoblast / osteoclast coupling coefficients.
- **Source aliases:** `TH-BMD`, `THBMD`, `B_THBMD` (Berkhout 2015 NONMEM compartment name).
- **Example models:** `Berkhout_2015_osteoporosis_placebo_qsp.R`.
- **Notes:** Founding example Berkhout 2015 (doi:10.1002/psp4.12006).

### osteoblast (**canonical relative active-osteoblast state (bone-remodelling QSP)**)
- **Type:** compartment
- **Role:** Dimensionless relative active-osteoblast state (y = B / B_0) in reduced Lemaire / Post 2013 bone-remodelling QSP models. Starts at 1 at menopause onset (baseline). Driven by k_B * (piz1 - y) where piz1 = pi_z / pi_1 and pi_z = z / (z + z_s). Feeds the BSAP transducer (`BSAP = BSAP_0 * y^q_BSAP`) and the BMD indirect-response production term (`(1 + D_AOB * osteoblast)`).
- **Source aliases:** `y`, `AOB`, `BM_AOR` (Berkhout 2015 NONMEM compartment name).
- **Example models:** `Berkhout_2015_osteoporosis_placebo_qsp.R`.
- **Notes:** Registered 2026-07-24 alongside the Berkhout 2015 extraction. Founding example Berkhout 2015 (doi:10.1002/psp4.12006); the state originates in the mechanistic model of Lemaire et al. J Theor Biol 2004 as reduced by Schmidt et al. 2011 and applied to placebo-arm osteoporosis by Post et al. 2013 / Berkhout et al. 2015. Sibling state `osteoclast`.

### osteoclast (**canonical relative active-osteoclast state (bone-remodelling QSP)**)
- **Type:** compartment
- **Role:** Dimensionless relative active-osteoclast state (z = C / C_0) in reduced Lemaire / Post 2013 bone-remodelling QSP models. Starts at 1 at menopause onset (baseline). Driven by D_A * pi_1 * (fdbf - piz1 * z) where fdbf = (y * (1 + b_baseline) / (1 + f(t) * piz1^2)) * PCa carries the disease-progression `f(t) = exp(-k_estrogen * t)` and the placebo (calcium) modulation `PCa`. Feeds the NTX transducer (`NTX = NTX_0 * z^q_NTX`) and the BMD indirect-response degradation term (`(1 + D_AOC * osteoclast)`).
- **Source aliases:** `z`, `AOC`, `BN_AOB` (Berkhout 2015 NONMEM compartment name; note the paper's NONMEM $MODEL block has a typographic swap of AOR / AOB labels between compartments 2 and 3 but the equations DADT(2) / DADT(3) match this canonical assignment where osteoblast is `y` and osteoclast is `z`).
- **Example models:** `Berkhout_2015_osteoporosis_placebo_qsp.R`.
- **Notes:** Registered 2026-07-24 alongside the Berkhout 2015 extraction. Founding example Berkhout 2015 (doi:10.1002/psp4.12006). Sibling state `osteoblast`.

---

## Gastric / GI transit

### stomach (**canonical stomach compartment**)
- **Type:** compartment
- **Role:** Gastric / stomach compartment used by gastric-emptying transit models where the gastric mass-balance is resolved as a distinct state ahead of the duodenal absorption depot.
- **Source aliases:** none.
- **Example models:** `Guiastrennec_2016_paracetamol.R`, `Back_2018_fenofibrate.R`.
- **Notes:** Registered 2026-05-28 per the naming audit.

### gastric_remaining (**canonical algebraic gastric-emptying percent-of-meal-remaining observation**)
- **Type:** compartment
- **Role:** Algebraic observation variable holding the percentage of a test meal remaining in the stomach at time `t` (0 - 100). Used by pure-algebraic gastric-emptying meta-analyses that fit a % remaining vs. time curve directly (no drug, no dose, no mass-balance ODE). Distinct from the `stomach` compartment (which holds a mass-balance drug amount as an ODE state driving downstream absorption); `gastric_remaining` is the algebraic dependent variable of the meta-analysis and is not integrated as a state.
- **Source aliases:** paper narrative `GE` (Bonner 2015 Eq. 1), `% remaining`.
- **Example models:** `Bonner_2015_gastric_emptying.R` (double Weibull mixture % remaining = (100 - PR) * exp(-(t/gamma1)^beta1) + PR * exp(-(t/gamma2)^beta2)).
- **Notes:** Registered 2026-07-24 alongside the Bonner 2015 gastric-emptying meta-analysis extraction. Reserved for algebraic-observable use in gastric-emptying meta-analyses; do not overload for ODE-integrated drug-amount states (use `stomach` for those).

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

### a_skin (**canonical PBPK skin-amount compartment**)
- **Type:** compartment
- **Role:** Skin organ compartment in mass-balance PBPK extractions. Flow-limited: `Q_skin * (CA - C_skin/pc_skin)`.
- **Source aliases:** none.
- **Example models:** `Decrane_2023_oxyfluorfen_rat.R` (founding example), `Decrane_2023_oxyfluorfen_human.R`.
- **Notes:** The bare `skin` form is also canonical; PBPK organ-amount models use the `a_` prefix alongside `a_liver` / `a_kidney` / `a_muscle`.

### a_brain (**canonical PBPK brain-amount compartment**)
- **Type:** compartment
- **Role:** Brain organ compartment in mass-balance PBPK extractions, treated as a single flow-limited lump. Distinct from the `brain_<region>` namespace, which resolves anatomically separate CNS sub-regions for brain-distribution models.
- **Source aliases:** none.
- **Example models:** `Decrane_2023_oxyfluorfen_rat.R` (founding example), `Decrane_2023_oxyfluorfen_human.R`.

### a_blood (**canonical PBPK single well-mixed blood pool**)
- **Type:** compartment
- **Role:** Single well-mixed blood pool in a PBPK model that has no separate arterial / venous split and no lung compartment: it simultaneously supplies the arterial concentration `Cart` to every tissue and receives all venous return.
- **Source aliases:** none.
- **Example models:** `Decrane_2023_oxyfluorfen_rat.R` (founding example), `Decrane_2023_oxyfluorfen_human.R`.
- **Notes:** Use `a_arterial` + `a_venous` when the source model splits the two; use `a_blood` only when the source carries a single pool. Distinct from the bare `blood` form used by membrane-limited extractions such as `Parhiz_2024_mRNA_LNP.R`.

### a_thyroid_blood (**canonical PBPK thyroid vascular sub-compartment**)
- **Type:** compartment
- **Role:** Vascular sub-compartment of a diffusion-limited thyroid. Exchanges with the systemic blood pool by perfusion (`Q_thy`) and with the thyroid tissue sub-compartment by passive diffusion across a permeability-area product (`PA_thy`). Because it is blood, it carries no tissue:blood partition coefficient.
- **Source aliases:**
  - `Athyblood` / `Cthyblood` -- Decrane 2023 (paper p. 3).
- **Example models:** `Decrane_2023_oxyfluorfen_rat.R` (founding example), `Decrane_2023_oxyfluorfen_human.R`.
- **Notes:** Paired with `a_thyroid_tissue`; a thyroid split into these two states is the standard structure in the thyroid-disruptor / HPT-axis modelling family (Ekerot 2013, Willemin & Lumen 2017, Handa 2021).

### a_thyroid_tissue (**canonical PBPK thyroid tissue sub-compartment**)
- **Type:** compartment
- **Role:** Tissue sub-compartment of a diffusion-limited thyroid, where thyroid-hormone synthesis takes place and where a thyroid-disrupting chemical exerts its sodium-iodide-symporter (NIS) or thyroid-peroxidase (TPO) inhibition. Exchanges only with `a_thyroid_blood`, across the permeability-area product `PA_thy`.
- **Source aliases:**
  - `Athytissue` / `Cthytissue` -- Decrane 2023 (paper p. 3).
- **Example models:** `Decrane_2023_oxyfluorfen_rat.R` (founding example), `Decrane_2023_oxyfluorfen_human.R`.

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

### vp_subcutaneous (**canonical PBPK subcutaneous-injection-site vascular concentration**)
- **Type:** compartment
- **Role:** Vascular concentration in the subcutaneous injection-site compartment
  of a whole-body PBPK model that bifurcates skin into a small
  dose-receiving subcutaneous depot plus a "rest of skin" compartment. The
  companion sub-compartments follow the standard membrane-limited pattern
  (`bc_subcutaneous`, `eu_subcutaneous`, `eb_subcutaneous`,
  `fr_subcutaneous`, `is_subcutaneous`); the SC dose is administered into
  `is_subcutaneous` (the interstitial space), not into a `depot`.
  Distinct from the canonical `depot`, which is an empirical
  first-order absorption compartment with no physiological volume or flow.
- **Source aliases:**
  - `SC` -- Kumar 2024 Figure 1 / Table 1.
- **Example models:** `Kumar_2024_mAb_popPBPK_sc.R`.

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

### ANC, anc (**canonical absolute neutrophil count**)
- **Type:** compartment
- **Role:** Absolute neutrophil count PD output. `ANC` is the standard upper-case clinical abbreviation; `anc` is the lower-case sibling registered as canonical because some Friberg myelosuppression models assign the circulating-pool observation to a lower-case `anc <- circ` output variable (matching the all-lowercase ODE-state convention).
- **Source aliases:** none.
- **Example models:** myelosuppression PD models with ANC output; `Hansson_2013_sunitinib_myelosuppression.R` (as lower-case `anc`).
- **Notes:** `anc` added 2026-06-28 as a canonical lower-case sibling. For paired-output multi-cell-type Friberg models prefer the `circ_anc` suffix form (see `circ`); the bare lower-case `anc` is the single-output PD-observation form. Distinct in Type from the `anc` cell-type metabolite suffix (see the Cell-type suffixes section), which the runtime parser routes separately.

### PLT (**canonical platelet count**)
- **Type:** compartment
- **Role:** Platelet count PD output.
- **Source aliases:** none.
- **Example models:** thrombocytopenia PD models.

### PRU (**canonical P2Y12 reaction unit platelet-reactivity state**)
- **Type:** compartment
- **Role:** Platelet reactivity expressed in P2Y12 reaction units (PRU), the device-standardised readout of the Accumetrics / Werfen VerifyNow P2Y12 assay, carried as a turnover pool whose fractional loss rate is stimulated by an active antiplatelet metabolite. Reported normal baseline range 180-376 PRU; the therapeutic window targeted for P2Y12 inhibition is 70-150 PRU. Used both as the ODE state and as the observation output name.
- **Source aliases:**
  - `PRU` -- universal in the antiplatelet PK-PD literature.
  - `P2Y12 reaction unit` -- expanded form used in figure axes and table headers.
- **Example models:** `Jung_2024_clopidogrel.R` (doi:10.1002/psp4.13053; `d/dt(PRU) = kin - kout * (1 + emax * C_h4^hill / (ec50^hill + C_h4^hill)) * PRU`, `PRU(0) = kin / kout = 212.67`).
- **Notes:** Registered 2026-08-06 (sidecar `oare_PMC10787215` request-001 q2, option A). `pru` is the canonical LOWERCASE sibling for the ODE state itself, and `PRU` is the observation-output name -- the same split already registered for `ANC` / `anc`, and required because `checkModelConventions()` enforces lowercase ODE-state names while the biomarker observables are uppercase. A model therefore writes `d/dt(pru) <- ...; pru(0) <- kin / kout; PRU <- pru; PRU ~ add(addSd_PRU)`. Uppercase for the observable per the rule stated in the `OSTCALC` entry ("Uppercase to match the sibling biomarker observables") and to match the surrounding clinical-acronym biomarker family (`INR`, `PT`, `aPTT`, `TT`, `PLT`, `P1NP`, `WBC`). Deliberately distinct from the two existing platelet entries, which are different quantities: `PLT` is platelet COUNT (cells/volume) and `integrin` is platelet alpha2-integrin EXPRESSION; `PRU` is platelet REACTIVITY (an assay-defined aggregation-response unit) and is not convertible to either. Registered as a general canonical rather than a paper-specific state because the VerifyNow readout is the shared endpoint across the antiplatelet class, so prasugrel / ticagrelor / cangrelor PK-PD extractions should reuse it.

### pru (**canonical lowercase P2Y12-reaction-unit ODE state**)
- **Type:** compartment
- **Role:** Lowercase sibling of `PRU`, used for the platelet-reactivity turnover ODE state itself so that the state name follows the all-lowercase ODE-state convention while the observation output keeps the uppercase clinical acronym. See the `PRU` entry above for the quantity, the assay, and the reference ranges.
- **Source aliases:** none.
- **Example models:** `Jung_2024_clopidogrel.R` (`d/dt(pru)`, `pru(0) <- kin / kout`, observed as `PRU <- pru`).
- **Notes:** Registered 2026-08-06 alongside `PRU`. Exactly parallel to the registered `ANC` / `anc` pair: the uppercase form is the biomarker observable, the lowercase form is the state.

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

### OSTCALC (**canonical osteocalcin**)
- **Type:** compartment
- **Role:** Osteocalcin bone-turnover biomarker. Renamed from the former `OC` canonical on 2026-07-26: `OC` was overloaded across five unrelated concepts in this codebase (osteocalcin, oseltamivir carboxylate, oral contraceptive, ovarian cancer, and the omega-correlation column of a source table). Uppercase to match the sibling biomarker observables (`P1NP`, `PSA`, `PLT`, `WBC`, `MBL`) and, in particular, the sister model `Shoji_2017_fosdagrocorat_p1np`, which names its observable `P1NP` in the identical structural slot.
- **Source aliases:** none.
- **Example models:** osteoporosis PD models.

### TT (**canonical total testosterone / thrombin time**)
- **Type:** compartment
- **Role:** Total testosterone (endocrinology) or thrombin time (coagulation) -- paper-dependent; both share the TT abbreviation in the contexts where it appears.
- **Source aliases:** none.
- **Example models:** endocrinology / coagulation PD models.

### QTc, QTcF, QTcI, QTcS (**canonical heart-rate-corrected QT interval**)
- **Type:** compartment
- **Role:** Heart-rate-corrected QT interval (electrocardiographic PD endpoint), typically expressed in ms. Used as the observation variable in direct-effect / linear concentration-QTc models of drug-induced QT prolongation (cardiac-safety / thorough-QT studies, e.g. quinidine, moxifloxacin, sotalol, glasdegib). `QTc` is the generic canonical; `QTcF` (Fridericia correction) and `QTcS` (study-specific correction) are registered as canonical sibling names because the Fostvedt 2021 glasdegib models use the correction-specific name directly as the observation variable.
- **Source aliases:** `QTcB` (Bazett), `QTcI` (individual correction) -- translate to `QTc` and record the correction in the model file's description / vignette.
- **Example models:** `Shin_2006_quinidine_QT.R` (Bazett-corrected QT interval; founding example), `Fostvedt_2021_glasdegib_QTcF.R` (Fridericia, as `QTcF`), `Fostvedt_2021_glasdegib_QTcS.R` (study-specific correction, as `QTcS`).
- **Notes:** `QTcF` / `QTcS` promoted from translate-to-`QTc` aliases to canonical sibling names 2026-06-28 so single-output models that name the observation by its specific correction (rather than the generic `QTc`) pass the convention check. New models should still prefer the generic `QTc` where the correction is incidental; use the specific name only when the correction is the defining feature of the endpoint (as in the paired Fostvedt 2021 QTcF / QTcS glasdegib analyses).

### serumK (**canonical serum potassium**)
- **Type:** compartment
- **Role:** Serum potassium concentration PD output / turnover-state, in mmol/L. Used as the observation variable in indirect-response / turnover models of drug-induced potassium shifts (mineralocorticoid-receptor antagonists, potassium-sparing diuretics, RAAS inhibitors). Standard clinical-laboratory biomarker reported on essentially every comprehensive metabolic panel; KDIGO and ESC thresholds for hyperkalemia are at 5.5 and 6.0 mmol/L. Distinct from any drug-PK central compartment because the modelled species is the endogenous electrolyte rather than the dosed drug.
- **Source aliases:** `K`, `K+`, `serum_K`, `POTAS` -- translate to `serumK` when assembling input data; document the source-paper symbol in the model file's `description`.
- **Example models:** `Goulooze_2022_finerenone.R` (FIDELIO-DKD Phase III PKPD turnover model for finerenone effect on serum K; founding example).

### uacr, uacrObs (**canonical urine albumin-to-creatinine ratio**)
- **Type:** compartment
- **Role:** Urine albumin-to-creatinine ratio PD output / disease-progression state, in mg/g. Used in longitudinal albuminuria disease-progression models of chronic kidney disease, where UACR is integrated as a state whose fractional rate of change carries the progression rate and the drug effect (mineralocorticoid-receptor antagonists, SGLT2 inhibitors, RAAS inhibitors). Distinct from the covariate column `UACR`, which is the time-fixed observed BASELINE value used to scale parameters; this compartment is the modelled time course. Distinct from any drug-PK central compartment because the modelled species is an endogenous damage biomarker rather than the dosed drug. `uacrObs` is the canonical sibling name for the **observed / model-predicted** UACR when a model separates it from the integrated state -- registered for the same reason `QTcF` / `QTcS` are registered alongside `QTc`, so that a single-output model can name its observation variable directly. Use the sibling only when the two genuinely differ; when the state IS the observation, observe `uacr`.
- **Source aliases:** `UACR`, `A(1)` (Goulooze 2022 UACR control stream, `COMP=(PROGR)`), `A(3)` (Goulooze 2022 eGFR control stream, `COMP=(PROGRUACR)`) -- translate to `uacr`. `IPREDUACR` (Goulooze 2022 `$ERROR`) -- translate to `uacrObs`.
- **Example models:** `Goulooze_2022_finerenone_uacr.R` (founding example; FIDELIO-DKD Phase III UACR disease-progression model with a power-function finerenone effect through an effect compartment), `Goulooze_2022_finerenone_egfr.R` (same state embedded as the driver of chronic eGFR decline).
- **Notes:** The state / observation split is load-bearing in the founding example: Goulooze 2022 integrates the drug-free natural progression of UACR (`NATUACR = A(3)` in its own `$ERROR` block) and applies the finerenone and SGLT2i effects algebraically at the observation, which is what makes those effects instantaneously reversible.

### egfr, egfrObs (**canonical estimated glomerular filtration rate**)
- **Type:** compartment
- **Role:** Estimated glomerular filtration rate PD output / disease-progression state, in mL/min/1.73 m^2. Used in longitudinal renal-function-decline models where eGFR is integrated as a state carrying a chronic decline slope, an approach to a low-eGFR stabilisation point, and acute reversible haemodynamic drug effects. Distinct from the covariate column `CRCL`, which is the time-fixed observed BASELINE renal function used to scale parameters; this compartment is the modelled time course. `egfrObs` is the canonical sibling name for the **observed / model-predicted** eGFR when a model separates it from the integrated state; see the `uacr` / `uacrObs` entry for the rationale.
- **Source aliases:** `eGFR`, `EGFR`, `EGFREPI`, `A(1)` (Goulooze 2022 eGFR control stream, `COMP=(PROGR)`) -- translate to `egfr`. `IPRED` of the eGFR control stream -- translate to `egfrObs`.
- **Example models:** `Goulooze_2022_finerenone_egfr.R` (founding example; FIDELIO-DKD Phase III eGFR model with a constant chronic slope, an exponential stabilisation function toward 16.1 mL/min/1.73 m^2, a reversible acute finerenone effect via an effect compartment, and model-predicted UACR driving the chronic slope).
- **Notes:** In the founding example the acute finerenone effect and the acute SGLT2 inhibitor effect are applied to `egfr` at the observation rather than inside `d/dt(egfr)`, which is precisely what makes the acute eGFR decline fully reversible on discontinuation.

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

### hfs_grade (**canonical hand-foot syndrome grade**)
- **Type:** compartment
- **Role:** Hand-foot syndrome (palmar-plantar erythrodysesthesia) NCI-CTC ordinal grade (0 / 1 / 2 / 3+) PD output. In Hansson 2013 sunitinib the reported output is the expected HFS grade given the previous state, derived from a first-order-Markov + proportional-odds transition model.
- **Source aliases:** none.
- **Example models:** `Hansson_2013_sunitinib_hfs.R`.
- **Notes:** Registered 2026-06-28. Sibling endpoint to `fatigue_grade` (Hansson 2013c sunitinib), which shares the same Markov + proportional-odds structure.

### dbp (**canonical diastolic blood pressure PD state**)
- **Type:** compartment
- **Role:** Diastolic blood pressure (mmHg) indirect-response turnover state and PD output. In Hansson 2013 sunitinib, `dbp` turns over via a stimulated zero-order production `kin` and first-order loss `kout`, with sunitinib AUC linearly stimulating `kin`; the state both carries the ODE and is the single observation variable.
- **Source aliases:** none.
- **Example models:** `Hansson_2013_sunitinib_dbp.R`.
- **Notes:** Registered 2026-06-28. Holds a blood-pressure value (mmHg), not a drug concentration; the related drug-induced *relative* change covariate used downstream is `DBP_REL` in `covariate-columns.md`.

### sbp (**canonical systolic blood pressure PD state**)
- **Type:** compartment
- **Role:** Systolic blood pressure (mmHg) indirect-response turnover state and PD output; the systolic sibling of `dbp`. In Ibrahim 2023 ibrutinib, `sbp` turns over via a zero-order production `kin` and first-order loss `kout` fed by a single upstream `transit1` compartment, with the daily ibrutinib AUC(0-24) stimulating `kin` through an Emax function; the state both carries the ODE and is the single observation variable.
- **Source aliases:** none.
- **Example models:** `Ibrahim_2023_ibrutinib_sbp.R`.
- **Notes:** Registered 2026-07-30 alongside the Ibrahim 2023 ibrutinib extraction, completing the blood-pressure pair with the previously registered `dbp` (Hansson 2013 sunitinib). Holds a blood-pressure value (mmHg), not a drug concentration. Systolic and diastolic pressure are fitted as separate models in both founding papers, so keep them as two states in two files rather than collapsing them into one multi-output model; the related drug-induced *relative* change covariate used downstream is `DBP_REL` in `covariate-columns.md`.

### bm (**canonical delayed biomarker-signal effect state**)
- **Type:** compartment
- **Role:** Effect-compartment-smoothed (first-order `ke0`) delayed biomarker-signal state that drives a downstream PD endpoint. In Hansson 2013 sunitinib (HFS / fatigue models) `bm` is the delayed relative-change-from-baseline signal of soluble VEGFR-3: `d/dt(bm) = ke0 * (bm_input - bm)` with `bm_input = (svegfr3 - baseline) / baseline`, which shifts the proportional-odds baseline logits.
- **Source aliases:** none.
- **Example models:** `Hansson_2013_sunitinib_hfs.R`.
- **Notes:** Registered 2026-06-28. `bm` is the paper's notation for the delayed sVEGFR-3 "biomarker" signal; it is a transduction / effect-compartment state, not a clinical endpoint, and is scoped to the Hansson 2013 sunitinib adverse-effect models. Reuse only for the same delayed-biomarker-signal role, not as a generic short name.

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

### hae_attacks (**canonical per-4-week hereditary-angioedema attack-count PD output**)
- **Type:** compartment
- **Role:** Per-4-week normalized hereditary-angioedema (HAE) attack count / rate PD output used by count-outcome exposure-response models that relate on-treatment HAE attack rate to a plasma-prekallikrein exposure driver. The state holds the Poisson-rate mean at each observation row (algebraic sigmoidal Emax function of `PKK`); the observation model is `hae_attacks ~ pois(attack_rate)`.
- **Source aliases:** none.
- **Example models:** `Singh_2025_donidalorsen.R`.
- **Notes:** Ratified canonically alongside the Singh 2025 donidalorsen exposure-response extraction. Companion covariates are `PKK` (time-varying prekallikrein exposure driver), `PKK_BL` (per-subject baseline PKK), and `HAERATE_BL` (per-subject baseline attack rate) documented in `inst/references/covariate-columns.md`. Sibling count-outcome PD-output canonicals include `msHeadacheDays`, `migraineDays`, `cel_count`, and `score` (each paper-specific by disease domain).

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

### dopamine (**canonical extracellular dopamine output**)
- **Type:** compartment
- **Role:** Extracellular dopamine PD output (CNS neurotransmitter turnover); the state holds the dopamine concentration in brain interstitial fluid as sampled by intracerebral microdialysis.
- **Source aliases:** none.
- **Example models:** `Dias_2024_quetiapine_rat.R`.
- **Notes:** Neurotransmitter sibling of the existing `prolactin` endocrine-biomarker output, and registered as canonical for the same reason `insulin` and `glucose` are: an extracellular-neurotransmitter state generalises to any CNS turnover / precursor-pool model, so it should not be hidden behind `paper_specific_compartments`. Name the state for the neurotransmitter itself; a paired precursor pool uses the canonical `precursor1` chain rather than a `dopamine_pool` variant, and the interval-integral bookkeeping state that microdialysis models need to form a collection-interval mean stays paper-specific (`auc_dopamine` in the founding model).

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

### nows (**canonical neonatal opioid withdrawal severity-score PD state**)
- **Type:** compartment
- **Role:** Indirect-response turnover state and single-output PD endpoint for the MOTHER NAS (Modification of Finnegan / Modified Neonatal Abstinence Scoring) severity score in neonatal opioid withdrawal syndrome (NOWS). The state carries the ODE `d/dt(nows) = kin * (1 + nowst) - kout * nows * effect_drug` (Eudy-Byrne 2021 Results / Model development), where `nowst = nowsmax * exp(-nowsm * PNA_days)` describes the natural withdrawal-severity decay with postnatal age and `effect_drug = 1 + emax * C / (ec50 + C)` describes opioid-agonist stimulation of NAS-score elimination.
- **Source aliases:** `NOWS`, `MOTHER_NAS`, `MNAS`, `NAS` -- equivalent paper notation for the same clinical scoring endpoint. NOWS = Neonatal Opioid Withdrawal Syndrome; NAS = Neonatal Abstinence Syndrome (older synonymous term); MOTHER NAS = the specific 19-item Modified Neonatal Abstinence scoring instrument used in the MOTHER, BBORN, and BPHORE clinical trials.
- **Example models:** `EudyByrne_2021_buprenorphine.R`.
- **Notes:** Registered 2026-07-25 alongside the Eudy-Byrne 2021 buprenorphine PD extraction. Companion paper-specific parameters `nowsmax` (unitless score baseline at PNA = 0) and `nowsm` (natural NAS decay rate, 1/day) live in `parameter-names.md`. Holds a NAS severity score (unitless integer scale, typically 0-40), not a drug concentration; the paired drug driver is buprenorphine plasma concentration supplied via the covariate `Cbuprenorphine`.

---

### MCC (**canonical maximum cystometric capacity**)
- **Type:** compartment
- **Role:** Maximum cystometric capacity (MCC), the urodynamic bladder-volume PD endpoint measured by multichannel cystometry, in mL. Used in exposure-response models of overactive bladder (OAB) and neurogenic detrusor overactivity (NDO) antimuscarinics, where the maximum attainable MCC is commonly anchored to the age-based pediatric expected bladder capacity (EBC) rather than estimated. All-caps because MCC is the standard urodynamics abbreviation and is never spelled out in the source literature after first use, matching the `ANC` / `PLT` / `WBC` / `RBC` clinical-abbreviation precedent in this register.
- **Source aliases:** none.
- **Example models:** `Sano_2023_fesoterodine_mcc.R` (Emax exposure-response of 5-HMT average steady-state concentration on MCC in pediatric NDO; baseline MCC and the EBC ceiling both scale with age by the same `(AGE + 1)/13` factor, and residual error is combined proportional plus additive as `propSd_MCC` / `addSd_MCC`).

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
- **Role:** Bare lung organ compartment in full-body PBPK extractions. Also used for the *uninvolved* (non-lesion) lung tissue in site-of-action models that resolve a diseased `lesion` compartment separately -- see `lesion`.
- **Source aliases:**
  - `CUlung` -- uninvolved-lung state in `Mehta_2023_bedaquiline_mpbpk.R` (Mehta 2023 ESM S2).
- **Example models:** `Zhang_2011_nutlin3a.R`, `Mehta_2023_bedaquiline_mpbpk.R`.

### lesion (**canonical site-of-action disease-lesion compartment**)
- **Type:** compartment
- **Role:** Drug concentration within a pathological tissue lesion treated as a distinct site of action -- e.g. a tuberculosis cavitary / necrotic lung lesion, whose caseum is a survival niche for non-replicating bacteria and where exposure differs from both plasma and the surrounding healthy organ. The compartment equilibrates with systemic blood scaled by a unitless penetration ratio, `d/dt(lesion) = k_le * (Cbld * R_le - lesion)` with `k_le = Qc / V_lesion`, so the state holds a **concentration**, not an amount. Pair it with the bare organ compartment for the uninvolved tissue (`lung` for pulmonary TB). Metabolite variants follow the standard `<canonical>_<metab>` pattern (e.g. `lesion_m2`). Distinct from the oncology lesion-state *suffixes* `vact` / `vell` / `dens` (Schindler 2017) and from `cel_count` (a multiple-sclerosis lesion cell-count PD output), none of which is a drug-concentration compartment.
- **Source aliases:**
  - `Cles` -- lesion state in `Mehta_2023_bedaquiline_mpbpk.R` (Mehta 2023 ESM S2).
- **Example models:** `Mehta_2023_bedaquiline_mpbpk.R` (founding example; also `lesion_m2` for the M2 metabolite), `Mehta_2023_pretomanid_mpbpk.R`, `Mehta_2023_pyrazinamide_mpbpk.R`.

### pleura (**canonical pleural-space fluid compartment**)
- **Type:** compartment
- **Role:** The pleural space, a serous cavity modelled as a fluid sub-compartment of the lung. Pleural liquid is a microvascular filtrate that flows in through the parietal pleural capillaries and is removed via lymphatic stomata in the parietal pleura, so the state is fed by lung efflux at a pleural fluid flow `q_pleura` and drained at the same flow: `d/dt(pleura) = q_pleura * cv_lung - q_pleura * Cpleura`. Because it is a fluid space rather than a perfused tissue it carries **no** tissue:plasma partition coefficient -- its outflow uses the pleural concentration directly, not `Cpleura / kp_pleura`. Distinct from `isf` (the generic mAb interstitial-fluid space) and from `ecf` (brain / tumour extracellular fluid in microdialysis models): the pleural cavity is an anatomically separate serous cavity and a model could legitimately carry `pleura` alongside either of those. Pleural tuberculosis is the second most common form of extrapulmonary TB, and pleural effusion is also a site of interest for oncology and anti-infective penetration studies.
- **Source aliases:**
  - `Pl` / `C_Pl` -- pleural state in `Ramachandran_2023_*_pbpk.R` (Appendix S1 section 3).
- **Example models:** `Ramachandran_2023_rifampicin_pbpk.R` (founding example), `Ramachandran_2023_ethambutol_pbpk.R`, `Ramachandran_2023_isoniazid_pbpk.R`, `Ramachandran_2023_pyrazinamide_pbpk.R`.
- **Notes:** Pair with `lnode` when a model resolves both major extrapulmonary-TB sites. The volume and flow are per-kg quantities (0.3 mL/kg and 0.15 mL/kg/h in the founding example), not fractions of body weight or cardiac output.

### gut_lumen (**canonical gut-luminal drug reservoir**)
- **Type:** compartment
- **Role:** Non-absorbed drug held in the intestinal lumen, as an **amount** rather than a concentration. Receives biliary / hepatic output and drains by two competing first-order routes: reabsorption back into the perfused `gut` tissue at `kr` (enterohepatic circulation) and faecal transit out of the body at `kF`. Distinct from `gallbladder`, which models a storage organ with delayed, gated emptying producing a discrete secondary peak -- the gut lumen is a continuously-draining reservoir with no emptying delay. Also distinct from `gut` / `a_gut`, which are the perfused gut *tissue* with their own volume, blood flow, lymph flow, and partition coefficient; a model that carries enterohepatic recycling needs both states simultaneously.
- **Source aliases:**
  - `GL` / `A_GL` -- gut lumen state in `Ramachandran_2023_*_pbpk.R` (Appendix S1 section 1).
- **Example models:** `Ramachandran_2023_rifampicin_pbpk.R` (founding example; `kr = 0.17 /h` for rifampicin enterohepatic circulation), `Ramachandran_2023_ethambutol_pbpk.R`, `Ramachandran_2023_isoniazid_pbpk.R`, `Ramachandran_2023_pyrazinamide_pbpk.R` (all three with `kr = 0`, so the lumen is a terminal faecal-transit sink).

### brain (**canonical bare brain compartment**)
- **Type:** compartment
- **Role:** Bare brain organ compartment in full-body PBPK extractions.
- **Source aliases:** none.
- **Example models:** `Zhang_2011_nutlin3a.R`.

### trachea (**canonical bare trachea compartment**)
- **Type:** compartment
- **Role:** Bare trachea organ compartment in PBPK extractions that resolve the conducting airway as a distinct target tissue. Holds the total tissue (well-stirred) drug amount and is paired with `lung`, `blood` and `other` in respiratory-target minimal-PBPK models. Registered for parity with the surrounding bare-organ canonicals (`lung`, `brain`, `heart`, `bone`, `skin`, `other`); it is a plain anatomical organ name rather than a mechanism-specific state, so the metabolite-suffixed form `trachea_<metab>` (e.g. `trachea_dihydroart`) follows the usual parent + metabolite rule without a separate registration.
- **Source aliases:** none.
- **Example models:** `Kang_2023_pyronaridine_hamster_pbpk.R`, `Kang_2023_artesunate_hamster_pbpk.R`.

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

### skin (**canonical bare skin compartment**)
- **Type:** compartment
- **Role:** Bare skin tissue compartment in full-body PBPK extractions. Total tissue (well-stirred) drug concentration; paired with `lung`, `liver`, `kidney`, `spleen`, `brain`, `heart`, `muscle`, `adipose`, `bone`, `other` etc. in whole-body PBPK extractions that resolve skin as a distinct organ. The token `skin` already appears in the `vp_skin` canonical entry and in the `pbpkSubCompartmentRegex` valid-organ list, so this entry registers the bare-organ form for parity with the surrounding canonicals.
- **Source aliases:** none.
- **Example models:** `Gaohua_2012_pregnancy_pbpk_caffeine.R`, `Gaohua_2012_pregnancy_pbpk_metoprolol.R`, `Gaohua_2012_pregnancy_pbpk_midazolam.R`, `Levitt_2005_propofol_pbpk.R`.

### pancreas (**canonical bare pancreas compartment**)
- **Type:** compartment
- **Role:** Bare pancreas organ compartment in full-body PBPK extractions. Total tissue (well-stirred) drug concentration; one of the four splanchnic organs (with `stomach`, `gut`, `spleen`) whose venous outflow drains into `liver` via the portal vein rather than returning directly to blood. The token `pancreas` already appears in the `vp_pancreas` canonical entry and in the `pbpkSubCompartmentRegex` valid-organ list (alongside `lung`, `kidney`, `spleen`, `thymus`), so this entry registers the bare-organ form for parity with the surrounding canonicals -- exactly as `heart` and `skin` above.
- **Source aliases:**
  - `PAN` -- NONMEM `$MODEL` compartment label in Yau 2023 Appendix S1.
- **Example models:** `Yau_2023_diazepam_pbpk_kpu_human.R`, `Yau_2023_diazepam_pbpk_scalar_human.R`, `Yau_2023_diazepam_pbpk_kpu_rat.R`, `Yau_2023_diazepam_pbpk_scalar_rat.R`.

### tendon (**canonical bare tendon / connective-tissue compartment**)
- **Type:** compartment
- **Role:** Bare tendon (dense connective tissue) compartment in full-body PBPK extractions built on the Levitt PKQuest standard-human physiology, which resolves tendon as its own organ with a distinct mass and perfusion (3 kg, 0.6 L/h/kg at the 70 kg reference) separate from the `other` lumped remainder (5.56 kg, 1.2 L/h/kg). Total tissue (well-stirred) drug amount, paired with the surrounding bare organ compartments.
- **Source aliases:**
  - `Tendon` -- Table S1 row label in Pei 2023 (Pharmaceutics 15:2580) and Table 1 of Levitt 2002 (BMC Clin Pharmacol 2:5) / Levitt & Schnider 2005 (BMC Anesthesiol 5:4).
- **Example models:** `Pei_2023_tacrolimus_pbpk.R`.
- **Notes:** `Levitt_2005_propofol_pbpk.R` lumps the same source physiology's tendon into `other` (`v_other = 5.524 + 3` kg), which is exact for propofol because both organs share a partition coefficient AND the model needed only their combined return to venous blood. Lumping is NOT generally equivalent: tendon and `other` differ two-fold in perfusion (0.6 vs 1.2 L/h/kg), so a drug whose distribution into connective tissue is perfusion-limited equilibrates at different rates in the two organs even when their Kp values coincide. Extractions that keep the source's tendon row separate should use this canonical rather than re-lumping. Ratified 2026-08-05 alongside the Pei 2023 tacrolimus PBPK extraction (sidecar `oare_PMC10675244` q3 = A).

---

### skin_fat (**canonical bare lumped skin + fat compartment**)
- **Type:** compartment
- **Role:** Skin and subcutaneous fat lumped into a **single** well-stirred tissue compartment carrying one volume, one plasma flow and one tissue:plasma partition coefficient. Distinct from the separate bare `skin` and `adipose` canonicals: this is the combined tissue, used when a source model does not resolve the two organs independently. Use it only when the source genuinely lumps them -- a model that reports separate partition coefficients for skin and fat must use `skin` + `adipose`. The motivating case is veterinary residue-depletion PBPK: "skin and fat in natural proportions" is the statutory edible-tissue matrix against which both the Chinese (GB 31650) and European (Reg. 37/2010) maximum residue limits for poultry are set, so the combined tissue is the regulated unit and cannot be split without inventing per-organ partition coefficients the source does not report. Component volume fractions may still be documented separately (e.g. Yang 2023 Table 4 gives `Vcsk` 13.38% and `Vcfa` 13.40% of body weight, summing to the 26.78% of the combined compartment) while the ODE system carries exactly one state. Expected to recur across the poultry and swine residue-PBPK family. The associated partition coefficient follows the `lk_<organ>` pattern as `lk_skf`.
- **Source aliases:** `skin + fat`, `skin+fat`, `sf` (Yang 2023 subscript).
- **Example models:** `Yang_2023_diclazuril_chicken_pbpk.R` (founding example).

### salivary_gland (**canonical bare salivary-gland compartment**)
- **Type:** compartment
- **Role:** Bare salivary-gland tissue compartment; the lumped state representing all major salivary glands (parotid and submandibular) in semi-physiologic distribution / dosimetry models. Registered for parity with the surrounding bare-organ canonicals (`liver`, `kidney`, `spleen`, `pancreas`, `skin`, `heart`, `other`) rather than as a new mechanistic role. Snake-cased on the multi-word-organ pattern already used by `small_intestine` / `large_intestine` / `renal_cortex`, and singular (one lumped state, not one per gland) on the pattern used by `kidney` for both kidneys.
- **Source aliases:**
  - `compartment 2` / "salivary glands" -- Siebinga 2023 numbers the six states 1-6 (blood, salivary gland, kidney, liver, tumor, other) and refers to the lumped state as "salivary glands" (plural).
- **Example models:** `Siebinga_2023_lu177psma617.R`.
- **Notes:** The salivary glands are the dose-limiting organ for PSMA-targeted radioligand therapy (and a target of interest for any drug or radiotracer with salivary uptake), so the state generalises beyond one paper and warrants a canonical rather than a `paper_specific_compartments` declaration. In Siebinga 2023 this is the only compartment with saturable (capacity-limited) uptake, parameterised by a maximum binding capacity `bmax`; the canonical name carries no commitment to saturable vs first-order kinetics.

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

### ctob (**canonical tobramycin bath-concentration state**)
- **Type:** compartment
- **Role:** Tobramycin bath / medium concentration state (mg/L) in static-concentration time-kill and hollow-fiber infection-model time-kill PD; dosed by the user, drives the Hill-type tobramycin killing term on each bacterial subpopulation, and disrupts the bacterial outer membrane to lower the effective meropenem KC50 through the mechanistic-synergy term. Pairs with `cmem` in meropenem-plus-tobramycin combination MBMs.
- **Source aliases:** none.
- **Example models:** `Landersdorfer_2018_meropenem_tobramycin_PAO1.R`, `Landersdorfer_2018_meropenem_tobramycin_PAOmutS.R`.

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
  - `Garonzik_2016_daptomycin.R` -- single-drug: `bact_susceptible1` / `bact_susceptible2`, `bact_intermediate1` / `bact_intermediate2`, `bact_resistant1` / `bact_resistant2` (three subpopulations of decreasing daptomycin susceptibility, each with the two-state life cycle).
  - `Rees_2018_meropenem_ciprofloxacin.R` -- two-drug (meropenem + ciprofloxacin): `bact_susceptible_susceptible1` / `2`, `bact_resistant_intermediate1` / `2`, `bact_intermediate_resistant1` / `2`.
  - `Landersdorfer_2018_imipenem_tobramycin.R` -- two-drug (imipenem + tobramycin): same `bact_<drug1pheno>_<drug2pheno>` compound scheme with the two-state life-cycle digit.

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

### shbg (**canonical sex hormone-binding globulin biomarker**)
- **Type:** compartment
- **Role:** Sex hormone-binding globulin (SHBG) serum concentration (nmol/L) modelled as an indirect-response turnover state. In Jensen 2023 LNG-IUS 52 mg the SHBG state has a zero-order synthesis rate `kin_shbg = rbase_shbg * kout_shbg` linearly inhibited by a delayed unbound-LNG signal (`kin_shbg * (1 - inh)`) and a first-order elimination rate `kout_shbg`; the drug-free steady state equals `rbase_shbg`. SHBG also binds LNG with dissociation constant KDS = 1.82 nmol/L, coupling the SHBG turnover to the unbound-LNG PK via a closed-form free-fraction equation.
- **Source aliases:** none.
- **Example models:** `Jensen_2023_lngIus52mg.R`.
- **Notes:** Registered 2026-07-24 following the `igf1` / `pth` / `iron` pattern for named endogenous protein biomarker states in indirect-response turnover models. SHBG has generalisable scope beyond a single paper (any drug modulating sex-hormone binding — contraceptives, HRT, androgens — may reuse this canonical).

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

### antipla2r (**canonical anti-PLA2R autoantibody titer biomarker**)
- **Type:** compartment
- **Role:** Serum anti-phospholipase A2 receptor (anti-PLA2R) IgG AUTOANTIBODY titer, the pathogenic-antibody biomarker of primary membranous nephropathy and the surrogate the field uses to judge immunological remission. Assayed by ELISA and reported in U/mL. Used as a PD output state whose decline follows B-cell-depleting therapy.
- **Source aliases:**
  - `anti-PLA2R titer`, `anti-PLA2R antibody titer` -- paper prose forms.
- **Example models:** `Liang_2024_rituximab_pla2r.R` (founding example; mono-exponential decline after rituximab, Liang 2024 Eq. 5).
- **Notes:** Registered 2026-08-06 (sidecar `oare_PMC11002205` request-001 q1, option A). The name deliberately marks the AUTOANTIBODY rather than its antigen: PLA2R itself is the podocyte receptor, so the bare name `pla2r` is left free for a future model of the receptor or of the immune-complex deposition it mediates. The `anti` prefix is what carries that distinction; the alternative `pla2r_ab` was rejected because `_ab` is not an established compartment suffix in this register (`tab` / `nab` / `ige` exist only as metabolite / binding-partner suffixes) and would have introduced a second new convention for one model. Joins the same one-word-analyte biomarker family as `igg`, `total_igg`, `crp`, `vegf` and `shbg`.

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

### t4_thyroid (**canonical thyroxine pool in thyroid tissue**)
- **Type:** compartment
- **Role:** Thyroxine (T4) held in thyroid tissue, as an **amount**. Produced by a zero-order thyroidal synthesis rate that a thyroid-disrupting chemical scales down and that serum TSH scales up; drains by deiodination to T3 within the gland and by first-order transfer into serum.
- **Source aliases:**
  - `AT4thy` -- Decrane 2023 (paper p. 3).
- **Example models:** `Decrane_2023_oxyfluorfen_rat.R` (founding example), `Decrane_2023_oxyfluorfen_human.R`.
- **Notes:** The location suffix is load-bearing: T4 exists as two kinetically distinct pools (thyroid tissue and serum), so a bare `t4` would be ambiguous. Paired with `t4_serum`. This is the hypothalamic-pituitary-thyroid (HPT) axis naming family shared with `t3_thyroid`, `t3_serum` and `tsh_serum`.

### t3_thyroid (**canonical triiodothyronine pool in thyroid tissue**)
- **Type:** compartment
- **Role:** Triiodothyronine (T3) held in thyroid tissue, as an **amount**. Produced both by direct thyroidal synthesis and by intrathyroidal deiodination of T4; drains by intrathyroidal metabolism and by first-order transfer into serum.
- **Source aliases:**
  - `AT3thy` -- Decrane 2023 (paper p. 3).
- **Example models:** `Decrane_2023_oxyfluorfen_rat.R` (founding example), `Decrane_2023_oxyfluorfen_human.R`.

### t4_serum (**canonical thyroxine pool in serum**)
- **Type:** compartment
- **Role:** Circulating thyroxine (T4) in serum, as an **amount** distributed in a T4-specific volume of distribution. Fed by transfer out of the thyroid; drains by systemic deiodination to T3 and by non-deiodinative loss. Its deviation from the basal amount drives the HPT negative-feedback terms on TSH production and turnover.
- **Source aliases:**
  - `AT4srm` -- Decrane 2023 (paper p. 4).
- **Example models:** `Decrane_2023_oxyfluorfen_rat.R` (founding example), `Decrane_2023_oxyfluorfen_human.R`.

### t3_serum (**canonical triiodothyronine pool in serum**)
- **Type:** compartment
- **Role:** Circulating triiodothyronine (T3) in serum, as an **amount** distributed in a T3-specific volume of distribution. Fed by transfer out of the thyroid and by systemic deiodination of serum T4; drains by first-order systemic loss.
- **Source aliases:**
  - `AT3srm` -- Decrane 2023 (paper p. 4).
- **Example models:** `Decrane_2023_oxyfluorfen_rat.R` (founding example), `Decrane_2023_oxyfluorfen_human.R`.

### tsh_serum (**canonical thyrotropin pool in serum**)
- **Type:** compartment
- **Role:** Circulating thyroid-stimulating hormone (thyrotropin, TSH) in serum, as an **amount**. Zero-order production is up-regulated and first-order turnover is slowed when serum T4 falls below its basal level; the resulting TSH concentration stimulates thyroidal T4 synthesis. This state closes the HPT feedback loop.
- **Source aliases:**
  - `ATSH` -- Decrane 2023 (paper p. 4).
- **Example models:** `Decrane_2023_oxyfluorfen_rat.R` (founding example), `Decrane_2023_oxyfluorfen_human.R`.
- **Notes:** Distinct from `tsh` should a future model carry TSH as a directly observed concentration biomarker rather than an amount state; the `_serum` suffix marks the HPT-axis amount-state family.

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
- **Notes:** This is a red-cell **count** state (cells per volume). It is NOT the state to use for a drug or active metabolite accumulating *inside* erythrocytes -- that role belongs to the `rbc_<analyte>` family (see "Intracellular red-cell analyte pools"), whose states carry drug concentration (e.g. umol/L). The two are physically different quantities and must not be conflated.

### rbc_mtx (**canonical intracellular red-cell methotrexate pool**)
- **Type:** compartment
- **Role:** Methotrexate (predominantly as methotrexate polyglutamates) accumulated inside red blood cells, carried as an ODE state in concentration units (umol/L). Fed by an influx term driven by the plasma methotrexate concentration and drained by a first-order efflux / elimination rate constant.
- **Source aliases:**
  - `XEMTX` -- used in `Gebhard_2023_methotrexate.R` (paper symbol `X_E^MTX`; the observation is reported as "E-MTX").
- **Example models:** `Gebhard_2023_methotrexate.R`, `Gebhard_2023_mercaptopurine_anc.R`.
- **Notes:** Registered 2026-07-30 as the founding member of the `rbc_<analyte>` family (see below). Distinct from `erythrocytes`, which is a red-cell *count* pool.

### rbc_tgn (**canonical intracellular red-cell thioguanine-nucleotide pool**)
- **Type:** compartment
- **Role:** 6-thioguanine nucleotides (TGN), the active metabolites of 6-mercaptopurine, accumulated inside red blood cells and carried as an ODE state in concentration units (umol/L). Fed by a saturable (Michaelis-Menten) or linear influx term driven by the plasma 6-mercaptopurine concentration and drained by a first-order efflux / elimination rate constant.
- **Source aliases:**
  - `XE6MP` -- used in `Gebhard_2023_mercaptopurine.R` (paper symbol `X_E^6MP`; the observation is reported as "E-TGN").
- **Example models:** `Gebhard_2023_mercaptopurine.R`, `Gebhard_2023_mercaptopurine_anc.R`.
- **Notes:** Registered 2026-07-30 alongside `rbc_mtx`. In the Gebhard 2023 PKPD model this state is the driver of the Friberg myelosuppression effect function (`Edrug = slope * rbc_tgn`).

---

## Intracellular red-cell analyte pools (`rbc_<analyte>` namespace)

Naming rule: `rbc_<analyte>`, where `<analyte>` is the lowercase name of the species measured inside the erythrocyte (`mtx` for methotrexate / methotrexate polyglutamates, `tgn` for 6-thioguanine nucleotides). The state holds a **concentration** (typically umol/L), not an amount, because red-cell metabolite assays are reported per unit of red-cell volume (or per mmol hemoglobin, converted).

Validated by `rbcCompartmentRegex = "^rbc_[a-z0-9]+$"` in `R/conventions.R`, NOT by the `<canonical>_<metabolite-suffix>` machinery. The distinction is deliberate: for methotrexate the intracellular analyte is the *parent* drug, so listing `mtx` in `registeredMetabolites` would record a parent drug in the metabolite register and mislead anyone reading that list. The regex lets any new red-cell analyte follow without further register edits, while keeping `registeredMetabolites` semantically clean.

Use `rbc_<analyte>` when the model carries the red-cell pool as an ODE state with its own influx / efflux kinetics. When a paper instead models red-cell drug **algebraically** as a fixed proportion of plasma (e.g. `Crbc <- slprbc * Cc` in `Gastonguay_2005_efaproxiral.R`, or the saturable whole-blood binding in `Storset_2014_tacrolimus.R`), no compartment is needed and no `rbc_` state should be introduced.

The corresponding transport parameters are `lkinf_rbc` / `lkeff_rbc` (first-order influx / efflux rate constants) and `lvmax_rbc` / `lkm_rbc` (saturable influx); see `parameter-names.md`.

---

## Cell populations and lymphoid tissues (continued)

### cells (**canonical generic cell population**)
- **Type:** compartment
- **Role:** Generic cell population PD output.
- **Source aliases:** none.
- **Example models:** `Jager_2011_gemtuzumab.R`.

### viability (**canonical normalised cell-viability PD state**)
- **Type:** compartment
- **Role:** Percent cell viability of an in-vitro cytotoxicity assay, normalised to the untreated vehicle control, carried as a dynamic PD state (`d/dt(viability)`). Its baseline is the fixed or fitted `rbase` (typically 100%). Used for colorimetric / metabolic viability readouts (CCK-8, MTT, resazurin) where the measurement is an absorbance ratio and no absolute cell number is observed.
- **Source aliases:**
  - `R` -- Mody 2023 notation, `dR/dt = kg*R - K3*R`.
- **Example models:** `Mody_2023_doxorubicin_dexrazoxane_jimt1.R` (founding example), `Mody_2023_doxorubicin_dexrazoxane_mdamb468.R`, `Mody_2023_doxorubicin_dexrazoxane_clinical_jimt1.R`, `Mody_2023_doxorubicin_dexrazoxane_clinical_mdamb468.R`, `Mody_2023_doxorubicin_dexrazoxane.R`.
- **Notes:** Distinct from `cells`, which holds an absolute cell count or number. `cells` was explicitly considered and rejected for this role: a percentage normalised to a control is not a count, and conflating the two would make the units of a `cells` state unknowable from its name. Because the state IS the observed quantity, `viability` is also canonical as a single-output observation variable (`viability ~ add(addSd_viability)`). Ratified by operator sidecar on 2026-08-04 with the Mody 2023 breast-cancer extraction.

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

### cumhaz_pfs (**canonical progression-free-survival cumulative-hazard**)
- **Type:** compartment
- **Role:** Progression-free-survival cumulative-hazard state in oncology TTE sub-models that expose PFS alongside OS as two parallel time-to-event endpoints. Integrates the PFS hazard so the PFS survivor function is `S_pfs = exp(-cumhaz_pfs)`; paired with `cumhaz_os` for the OS sub-model when both endpoints are simulated jointly. Distinct from `cumhaz_os` (overall-survival endpoint), `cumhaz_drop` (dropout / trial withdrawal), and `cumhaz_cens` (informative censoring distribution).
- **Source aliases:** none.
- **Example models:** `Franzese_2026_pdl1_nsclc_mbma.R` (MBMA semi-parametric proportional hazards for PFS with monthly discrete baseline hazard from Franzese 2026 Table S2; five treatment-category-specific ORR-slopes, PD-(L)1-monotherapy hazard-intercept shift, and chemotherapy time-dependent baseline-hazard shift).
- **Notes:** Registered 2026-07-24 alongside the Franzese 2026 MBMA extraction. Member of the `cumhaz_<type>` multi-hazard family alongside `cumhaz`, `cumhaz_os`, `cumhaz_drop`, `cumhaz_cens`, `cumhaz_1st`, and `cumhaz_2nd`. Applies whenever an oncology model separately parameterises the PFS hazard alongside an OS hazard (both endpoints simulated in the same solve so the two survival curves are internally consistent given the same shared covariates and study-strata random effects).

### cumhaz_drop (**canonical dropout cumulative-hazard**)
- **Type:** compartment
- **Role:** Dropout cumulative-hazard state.
- **Source aliases:**
  - `cumHaz_drop` -- prior canonical name (pre-2026-06-19 lowercase standardization).
- **Example models:** `Schindler_2016_sunitinib.R`.
- **Notes:** Renamed from `cumHaz_drop` to `cumhaz_drop` on 2026-06-19 per the canonical-register standardization audit (lowercase compartment-name standardization).

### cumhaz_cens (**canonical censoring cumulative-hazard**)
- **Type:** compartment
- **Role:** Censoring-distribution cumulative-hazard state in time-to-event sub-models that describe informative dropout / censoring with a separate hazard alongside the event hazard. Integrates the censoring hazard so the censoring survival is `sur_cens = exp(-cumhaz_cens)`, used for simulation-based dropout.
- **Source aliases:** none.
- **Example models:** `Hansson_2013_sunitinib_os.R`.
- **Notes:** Registered 2026-06-28. Member of the `cumhaz_<type>` multi-hazard family alongside `cumhaz`, `cumhaz_os`, and `cumhaz_drop`; `_cens` denotes the censoring distribution specifically (in Hansson 2013 a separate Weibull `lambdacens` / `alphacens` censoring model exposed for forward dropout simulation).

### cumhaz_1st (**canonical first-event cumulative-hazard (repeated-time-to-event sub-model)**)
- **Type:** compartment
- **Role:** Cumulative-hazard state for the time-to-first-event Weibull sub-model in repeated-time-to-event models with distinct first-event and subsequent-event hazard sub-models (Abrantes 2010 / Lindauer 2017 parameterisation). Integrates the first-event hazard from t = 0 so the first-event survivor function is `sur = exp(-cumhaz_1st)` (the paper's `sur` for "time to first seizure"). Paired with `cumhaz_2nd` for the subsequent-event sub-model.
- **Source aliases:** none.
- **Example models:** `Lindauer_2017_lacosamide_seizure.R` (Weibull baseline hazard `hazard_1st = lam1_eff * p1 * (lam1_eff * (t + del))^(p1 - 1)`; the LHS survivor function `sur = exp(-cumhaz_1st)` is the observation variable for the first-seizure TTE sub-model).
- **Notes:** Registered 2026-07-03. Member of the `cumhaz_<type>` multi-hazard family alongside `cumhaz`, `cumhaz_os`, `cumhaz_drop`, and `cumhaz_cens`. The paired 2nd-event canonical is `cumhaz_2nd`.

### cumhaz_2nd (**canonical second-and-subsequent-event cumulative-hazard (repeated-time-to-event sub-model)**)
- **Type:** compartment
- **Role:** Cumulative-hazard state for the time-to-second-and-subsequent-event Weibull sub-model in repeated-time-to-event models with distinct first-event and subsequent-event hazard sub-models (Abrantes 2010 / Lindauer 2017 parameterisation). Integrates the subsequent-event hazard from t = 0; the derived subsequent-event survivor function `sur_2nd = exp(-cumhaz_2nd)` is exposed as a diagnostic output. Paired with `cumhaz_1st` for the first-event sub-model.
- **Source aliases:** none.
- **Example models:** `Lindauer_2017_lacosamide_seizure.R` (Weibull baseline hazard `hazard_2nd = lam2_eff * p2 * (lam2_eff * (t + del))^(p2 - 1)` with an IIV eta on ln(lam2), reflecting the paper's Table 3 finding of substantial between-subject variability in the subsequent-seizure hazard; SD 2.03 on ln(k2)).
- **Notes:** Registered 2026-07-03. Member of the `cumhaz_<type>` multi-hazard family. Applies whenever a model separately parameterises the first-event hazard (`cumhaz_1st`) and the subsequent-event hazard (`cumhaz_2nd`) with different Weibull scale / shape and, typically, different covariate coefficients on each sub-model, as first proposed by Abrantes et al. and adopted by Lindauer 2017 for lacosamide.

### sur (**canonical survival-probability output**)
- **Type:** compartment
- **Role:** Survival-probability output `S(t)` of a time-to-event sub-model, derived from the cumulative hazard (`sur = exp(-cumhaz)` for a proportional-hazard form, or `sur = 1 - Phi(z)` for an accelerated-failure-time log-normal form). Single PD output for forward-simulation TTE models.
- **Source aliases:** none.
- **Example models:** `Hansson_2013_sunitinib_os.R` (Weibull PH `sur = exp(-cumhaz)`; the observation variable), `Struemper_2025_tumorsize_OS_nsclc.R` (AFT log-normal `sur = 1 - pnorm(z_os)`, derived OS output).
- **Notes:** Registered 2026-06-28. Founding models expose `sur` with a small placeholder residual so the nlmixr2 likelihood machinery accepts the forward-simulation model.

### prob_scc (**canonical sputum-culture-conversion probability output**)
- **Type:** compartment
- **Role:** Marginal probability of occupying the sputum-culture-converted state at time `t` in a tuberculosis treatment-outcome multistate model, i.e. the probability that a patient has achieved sputum culture conversion (SCC) and has not since relapsed, dropped out or died. Sputum culture conversion is the standard early efficacy endpoint in TB drug development, so `prob_scc` is the natural single PD output to expose for a TB multistate model whose remaining state-occupancy probabilities (active TB, recurrent TB, dropout, death) are returned alongside it as ordinary `model()` outputs.
- **Source aliases:** none.
- **Example models:** `Lin_2024_TB_multistate.R` (five-state pharmacometric multistate model; `prob_scc <- s_converted` is the observation variable and carries the placeholder additive residual, while `prob_active_tb`, `prob_recurrent_tb`, `prob_dropout` and `prob_death` expose the other four state-occupancy probabilities).
- **Notes:** A probability output in `[0, 1]`, not a concentration or an amount. Distinct from `sur` (a survival probability derived from a cumulative hazard in a time-to-event sub-model) because `prob_scc` is a state-occupancy probability of a *transient, re-enterable* state: a patient can leave the converted state for recurrent TB and return to it, so `prob_scc` is not monotone in time and is not the complement of any cumulative hazard. Follows the `prob_<endpoint>` output-naming shape founded by `prob_roc`. Founding models expose `prob_scc` with a small placeholder residual so the nlmixr2 likelihood machinery accepts the forward-simulation model; the source analysis maximises an exact multistate event likelihood on the observed categorical state and has no observation-error model.

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

### dox (**canonical doxorubicin drug-name suffix**)
- **Type:** metabolite-suffix
- **Role:** Doxorubicin drug-name suffix for combination-therapy compartments, parameters and residual SDs (`central_dox`, `peripheral1_dox`, `peripheral2_dox`, `transit1_dox`, `conc_dox`, `lkmax_dox`, `lkc50_dox`, `lktr_dox`, `lkdeg_dox`, `addSd_Cc_dox`).
- **Source aliases:**
  - `DOX` -- Mody 2023 notation.
- **Example models:** `Mody_2023_doxorubicin_dexrazoxane_jimt1.R` (founding example), `Mody_2023_doxorubicin_dexrazoxane_mdamb468.R`, `Mody_2023_doxorubicin_dexrazoxane_clinical_jimt1.R`, `Mody_2023_doxorubicin_dexrazoxane_clinical_mdamb468.R`.
- **Notes:** The source paper's own abbreviation and unambiguous across the registry. Ratified by operator sidecar on 2026-08-04 together with `dexrazoxane`.

### dexrazoxane (**canonical dexrazoxane drug-name suffix**)
- **Type:** metabolite-suffix
- **Role:** Dexrazoxane drug-name suffix for combination-therapy compartments, parameters and residual SDs (`central_dexrazoxane`, `peripheral1_dexrazoxane`, `transit1_dexrazoxane`, `conc_dexrazoxane`, `lkmax_dexrazoxane`, `lkc50_dexrazoxane`, `lktr_dexrazoxane`, `lkdeg_dexrazoxane`, `addSd_Cc_dexrazoxane`).
- **Source aliases:**
  - `DEX` -- Mody 2023 notation.
- **Example models:** `Mody_2023_doxorubicin_dexrazoxane_jimt1.R` (founding example), `Mody_2023_doxorubicin_dexrazoxane_mdamb468.R`, `Mody_2023_doxorubicin_dexrazoxane_clinical_jimt1.R`, `Mody_2023_doxorubicin_dexrazoxane_clinical_mdamb468.R`.
- **Notes:** Deliberately spelled out in full rather than abbreviated to the paper's `DEX`, because `dex` would collide with dexamethasone and dexmedetomidine and permanently burn the token on the less common drug. The resulting asymmetry with the abbreviated `dox` is intentional; ratified by operator sidecar on 2026-08-04.

---

## In-vitro single-cell mRNA-translation states (Frohlich 2018 QSP)

The Frohlich 2018 multi-experiment NLME single-cell mRNA-transfection model (npj Syst Biol Appl 4:42) introduces three mechanistic states for an in-vitro gene-expression assay and one log-fluorescence observation. The bare names are short to match the source code (`transfection_ribo_syms.m` in the authors' Zenodo deposit doi:10.5281/zenodo.1228899) and the paper's symbol set (m, GFP, R, r-m complex).

### mrna (**canonical free cytoplasmic messenger RNA**)
- **Type:** compartment
- **Role:** Free (translatable) cytoplasmic mRNA copies per cell, normalised by the per-cell transfected mRNA amount `m0` so that the bolus event sets `mrna = 1`. Used by single-cell translation-kinetics models that explicitly model mRNA degradation, ribosome binding, and protein production. Distinct from the PBPK sub-compartment-regex prefix `mrna_<organ>` (membrane-limited intracellular pool in whole-body Parhiz 2024 mRNA-LNP / Sasaki 2022 BNT162b2 models), which is a different mechanistic concept (tissue-resolved bulk mRNA) and is matched by `pbpkSubCompartmentRegex`.
- **Source aliases:**
  - `MRNA` -- uppercase in `transfection_ribo_syms.m` (Frohlich 2018 Zenodo deposit).
  - `m` -- the paper's symbol for the normalised mRNA abundance.
  - `m_prime` (`m'`) -- the un-normalised mRNA abundance before the m0-normalisation transformation (paper Methods, Mathematical models for GFP translation).
- **Example models:** `Frohlich_2018_mRNA_translation.R`.
- **Notes:** State holds a dimensionless normalised copy-number-per-cell (m0 units). Bolus event encoded as a dose of amount 1 to `mrna` with `alag(mrna) = t0` to realise the per-cell transfection-onset time.

### gfp (**canonical green-fluorescent-protein reporter state**)
- **Type:** compartment
- **Role:** Concentration / amount per cell of the fluorescent reporter protein (eGFP or destabilized eGFP / d2eGFP) produced by ribosomal translation of transfected mRNA. Used by single-cell mRNA-transfection / reporter-translation assays where the observed quantity is the optical intensity from a fluorescent reporter. Includes the absorbed-into-state scale factor (paper-transformed `scale * GFP`); the linear observable is `gfp + offset` and the log observable is `log(gfp + offset)`.
- **Source aliases:**
  - `GFP` -- uppercase in `transfection_ribo_syms.m` (Frohlich 2018 Zenodo deposit).
  - `p` -- the paper's symbol for protein abundance.
- **Example models:** `Frohlich_2018_mRNA_translation.R`.
- **Notes:** Both eGFP and d2eGFP cohorts share the same `gfp` state; the cohort distinction is carried by the `STUDY_d2eGFP` covariate (it selects between two protein degradation rates `gamma_eGFP` and `gamma_d2eGFP`). Single-cell concentration / amount; no body-volume scaling.

### ribo (**canonical free ribosome state (single-cell mRNA-translation models)**)
- **Type:** compartment
- **Role:** Free (unbound, available-for-translation) ribosome concentration per cell in single-cell mRNA-transfection / reporter-translation assays with ribosomal rate-limited translation. Used by Frohlich 2018 model (ii) and (iv) (ribosomal-binding extensions of the standard two-stage gene-expression model). The bound mRNA-ribosome complex is the conservation residual `fracr0_m0 - ribo` (total ribosomes = free + bound = R0/m0 normalised), so no separate `complex` state is required.
- **Source aliases:**
  - `RIBO` -- uppercase in `transfection_ribo_syms.m` (Frohlich 2018 Zenodo deposit).
  - `R` -- the paper's symbol for ribosome abundance.
- **Example models:** `Frohlich_2018_mRNA_translation.R`.
- **Notes:** Dimensionless normalised ribosome count (m0 units). Initial condition is the per-cell free-ribosome / m0 ratio `fracr0_m0` (paper parameter `R0/m0`). Used jointly with the `mrna` and `gfp` states; the three-state system `(mrna, gfp, ribo)` describes one cohort. Distinct from `bact_<phenotype>` bacterial subpopulation states (different mechanism) and from `enzyme` / `enz_pool` (autoinduction reservoirs).

### logfluor (**canonical log-fluorescence-intensity observable (in-vitro single-cell)**)
- **Type:** compartment
- **Role:** Single-output observation variable `logfluor = log(gfp + offset)` for in-vitro single-cell mRNA-transfection / reporter-translation assays. The log-transform is applied to the linear `gfp + offset` predictor so that the residual error is additive on the log scale (the source paper uses this transform to convert multiplicative measurement noise on linear-fluorescence intensity into additive noise on log-fluorescence; Frohlich 2018 Methods, Data acquisition and quantitative image analysis). Used with `~ add(addSd_logfluor)` as the residual-error structure.
- **Source aliases:**
  - `y` -- the paper's symbol for the observable (Methods, "the output fluorescence y is assumed to be the sum of a signal proportional to the amount of eGFP (with scaling parameter scale) plus background fluorescence (offset)").
- **Example models:** `Frohlich_2018_mRNA_translation.R`.
- **Notes:** Not a state in the ODE system (`logfluor` is derived inside `model()` from the `gfp` state and the `offset` parameter). Registered as a canonical compartment so the single-output observation-variable check accepts `logfluor ~ add(addSd_logfluor)` for in-vitro / mechanistic models without forcing a rename to `Cc`.

---

## Muliaditan 2025 TfR-mediated brain mPBPK states

The Muliaditan 2025 translational minimal-PBPK (mPBPK) model for transferrin-receptor (TfR)-mediated brain delivery of monoclonal antibodies introduces a membrane-limited, FcRn-recycling brain-distribution scheme with its own `<location>_<subtype>` compartment namespace. Registered individually (rather than via a regex) per the 2026-06-28 naming-audit decision: the scheme is bespoke to this two-file model family (`Muliaditan_2025_mab_mpbpk_human.R`, `Muliaditan_2025_mab_mpbpk_nhp.R`), uses only three barrier locations (so it lacks the combinatorial organ explosion that motivates the `pbpkSubCompartmentRegex`), and each TfR-binding bookkeeping state has a distinct mechanistic role that warrants an explicit entry. All states hold an antibody **amount** (nmol); concentrations are derived in `model()` as `c_<state> = <state> / v_<volume>`. The two models also use the already-canonical `central` (plasma), `brain_vascular`, `csf`, and `lymph` states, which are not re-registered here.

The namespace has two families:

1. **Membrane-limited sub-compartments** `<location>_<subtype>` where location is `tissue` (lumped peripheral tissue), `bbb` (blood-brain-barrier endothelium), or `bcsfb` (blood-CSF-barrier / choroid-plexus epithelium), and subtype is `vasc` (vascular space), `endo_u` (endosomal unbound drug), `endo_b` (endosomal FcRn-bound drug), `isf` (interstitial fluid), or `fcrn` (free FcRn receptor pool); plus `brain_isf` (brain parenchymal ISF). This is the inverse-ordering analogue of the existing `<subtype>_<organ>` `pbpkSubCompartmentRegex` (Shah 2012 / Parhiz 2024).
2. **TfR-binding bookkeeping states** tracking the drug-TfR complex (`complex_<location>`) and the unbound-TfR pools at each barrier face (luminal `_lum` = blood side, abluminal `_abl` = brain / CSF side), plus the neuronal TfR sink. `delta_utfr_<barrier>` carries the deviation of luminal unbound TfR from its baseline so the luminal level is `utfr0 + delta`.

### tissue_vasc (**canonical mPBPK peripheral-tissue vascular space**)
- **Type:** compartment
- **Role:** Antibody amount in the lumped peripheral-tissue vascular space (volume `v_tissue_v`); exchanges with plasma by tissue plasma flow `qt` and lymph, and receives FcRn-recycled drug from `tissue_endo_b`.
- **Source aliases:** none.
- **Example models:** `Muliaditan_2025_mab_mpbpk_human.R`, `Muliaditan_2025_mab_mpbpk_nhp.R`.

### tissue_endo_u (**canonical mPBPK peripheral-tissue endosomal unbound drug**)
- **Type:** compartment
- **Role:** Unbound antibody in the peripheral-tissue endosomal space (volume `v_tissue_e`); taken up by pinocytosis (`clup_t`), then either binds FcRn (-> `tissue_endo_b`) or is degraded (`kdeg_endo`).
- **Source aliases:** none.
- **Example models:** `Muliaditan_2025_mab_mpbpk_human.R`, `Muliaditan_2025_mab_mpbpk_nhp.R`.

### tissue_endo_b (**canonical mPBPK peripheral-tissue endosomal FcRn-bound drug**)
- **Type:** compartment
- **Role:** FcRn-bound antibody in the peripheral-tissue endosome; recycled back to `tissue_vasc` (fraction `fr_tissue`) and `tissue_isf` (`1 - fr_tissue`) by `clup_t`.
- **Source aliases:** none.
- **Example models:** `Muliaditan_2025_mab_mpbpk_human.R`, `Muliaditan_2025_mab_mpbpk_nhp.R`.

### tissue_isf (**canonical mPBPK peripheral-tissue interstitial fluid**)
- **Type:** compartment
- **Role:** Antibody amount in the lumped peripheral-tissue interstitial fluid (volume `v_tissue_i`); fed by lymph-reflection convection from `tissue_vasc` and FcRn recycling, drained by lymph flow `lt`.
- **Source aliases:** none.
- **Example models:** `Muliaditan_2025_mab_mpbpk_human.R`, `Muliaditan_2025_mab_mpbpk_nhp.R`.

### tissue_fcrn (**canonical mPBPK peripheral-tissue free FcRn pool**)
- **Type:** compartment
- **Role:** Free (unoccupied) FcRn receptor pool in the peripheral-tissue endosome; consumed by binding `tissue_endo_u` and regenerated when `tissue_endo_b` is recycled out. Initial condition `fcrn_baseline * v_tissue_e`.
- **Source aliases:** none.
- **Example models:** `Muliaditan_2025_mab_mpbpk_human.R`, `Muliaditan_2025_mab_mpbpk_nhp.R`.

### bbb_endo_u (**canonical mPBPK blood-brain-barrier endosomal unbound drug**)
- **Type:** compartment
- **Role:** Unbound antibody in the blood-brain-barrier (BBB) endothelial endosome (volume `v_bbb_e`); taken up from `brain_vascular` and `brain_isf` (`clup_bbb`), then binds FcRn or is degraded.
- **Source aliases:** none.
- **Example models:** `Muliaditan_2025_mab_mpbpk_human.R`, `Muliaditan_2025_mab_mpbpk_nhp.R`.

### bbb_endo_b (**canonical mPBPK blood-brain-barrier endosomal FcRn-bound drug**)
- **Type:** compartment
- **Role:** FcRn-bound antibody in the BBB endothelial endosome; recycled to `brain_vascular` (fraction `fr_brain`) and `brain_isf` (`1 - fr_brain`) by `clup_bbb`. This trans-endothelial recycling is the mechanism of antibody brain entry.
- **Source aliases:** none.
- **Example models:** `Muliaditan_2025_mab_mpbpk_human.R`, `Muliaditan_2025_mab_mpbpk_nhp.R`.

### bbb_fcrn (**canonical mPBPK blood-brain-barrier free FcRn pool**)
- **Type:** compartment
- **Role:** Free FcRn receptor pool in the BBB endothelial endosome; binds `bbb_endo_u` and is regenerated on `bbb_endo_b` recycling. Initial condition `fcrn_baseline * v_bbb_e`.
- **Source aliases:** none.
- **Example models:** `Muliaditan_2025_mab_mpbpk_human.R`, `Muliaditan_2025_mab_mpbpk_nhp.R`.

### bcsfb_endo_u (**canonical mPBPK blood-CSF-barrier endosomal unbound drug**)
- **Type:** compartment
- **Role:** Unbound antibody in the blood-CSF-barrier (BCSFB, choroid-plexus) epithelial endosome (volume `v_bcsfb_e`); taken up from `brain_vascular` and `csf` (`clup_bcsfb`), then binds FcRn or is degraded.
- **Source aliases:** none.
- **Example models:** `Muliaditan_2025_mab_mpbpk_human.R`, `Muliaditan_2025_mab_mpbpk_nhp.R`.

### bcsfb_endo_b (**canonical mPBPK blood-CSF-barrier endosomal FcRn-bound drug**)
- **Type:** compartment
- **Role:** FcRn-bound antibody in the BCSFB epithelial endosome; recycled to `brain_vascular` (fraction `fr_brain`) and `csf` (`1 - fr_brain`) by `clup_bcsfb`.
- **Source aliases:** none.
- **Example models:** `Muliaditan_2025_mab_mpbpk_human.R`, `Muliaditan_2025_mab_mpbpk_nhp.R`.

### bcsfb_fcrn (**canonical mPBPK blood-CSF-barrier free FcRn pool**)
- **Type:** compartment
- **Role:** Free FcRn receptor pool in the BCSFB epithelial endosome; binds `bcsfb_endo_u` and is regenerated on `bcsfb_endo_b` recycling. Initial condition `fcrn_baseline * v_bcsfb_e`.
- **Source aliases:** none.
- **Example models:** `Muliaditan_2025_mab_mpbpk_human.R`, `Muliaditan_2025_mab_mpbpk_nhp.R`.

### brain_isf (**canonical mPBPK brain parenchymal interstitial fluid**)
- **Type:** compartment
- **Role:** Antibody amount in the brain parenchymal interstitial fluid (volume `v_brain_i`); fed by BBB ECF convection and BBB FcRn recycling, exchanges with `csf` at the brain-ECF flow, and is the site of neuronal-TfR binding. The whole-brain homogenate observation `Cbrain` is computed from this plus the BBB/BCSFB endosomal pools.
- **Source aliases:** none.
- **Example models:** `Muliaditan_2025_mab_mpbpk_human.R`, `Muliaditan_2025_mab_mpbpk_nhp.R`.

### complex_plasma (**canonical mPBPK drug-TfR complex, plasma**)
- **Type:** compartment
- **Role:** Drug-soluble-TfR complex in plasma (volume `vp`); formed by `kon_t` binding of plasma drug to plasma TfR (`tfrpt`), dissociating at `koff_t`, and internalized/eliminated at `kint`.
- **Source aliases:** none.
- **Example models:** `Muliaditan_2025_mab_mpbpk_human.R`, `Muliaditan_2025_mab_mpbpk_nhp.R`.

### complex_bbb_lum (**canonical mPBPK drug-TfR complex, BBB luminal**)
- **Type:** compartment
- **Role:** Drug-TfR complex on the luminal (blood) face of the BBB; formed from `brain_vascular` drug and luminal BBB TfR, transcytosed to the abluminal face at `ktrans`, and internalized at `kint`.
- **Source aliases:** none.
- **Example models:** `Muliaditan_2025_mab_mpbpk_human.R`, `Muliaditan_2025_mab_mpbpk_nhp.R`.

### complex_bbb_abl (**canonical mPBPK drug-TfR complex, BBB abluminal**)
- **Type:** compartment
- **Role:** Drug-TfR complex on the abluminal (brain) face of the BBB; arrives by transcytosis from `complex_bbb_lum`, dissociates to release drug into `brain_isf`, and is internalized at `kint`.
- **Source aliases:** none.
- **Example models:** `Muliaditan_2025_mab_mpbpk_human.R`, `Muliaditan_2025_mab_mpbpk_nhp.R`.

### complex_bcsfb_lum (**canonical mPBPK drug-TfR complex, BCSFB luminal**)
- **Type:** compartment
- **Role:** Drug-TfR complex on the luminal (blood) face of the BCSFB; formed from `brain_vascular` drug and luminal BCSFB TfR, transcytosed to the abluminal face at `ktrans`, internalized at `kint`.
- **Source aliases:** none.
- **Example models:** `Muliaditan_2025_mab_mpbpk_human.R`, `Muliaditan_2025_mab_mpbpk_nhp.R`.

### complex_bcsfb_abl (**canonical mPBPK drug-TfR complex, BCSFB abluminal**)
- **Type:** compartment
- **Role:** Drug-TfR complex on the abluminal (CSF) face of the BCSFB; arrives by transcytosis from `complex_bcsfb_lum`, dissociates to release drug into `csf`, internalized at `kint`.
- **Source aliases:** none.
- **Example models:** `Muliaditan_2025_mab_mpbpk_human.R`, `Muliaditan_2025_mab_mpbpk_nhp.R`.

### complex_neuron (**canonical mPBPK drug-TfR complex, neuronal**)
- **Type:** compartment
- **Role:** Drug bound to neuronal TfR (total `tfrtotn`) in the brain ISF; formed by `kon_t` binding of `brain_isf` drug to unbound neuronal TfR, dissociating at `koff_t`, internalized at `kint`. Represents the target-mediated neuronal sink.
- **Source aliases:** none.
- **Example models:** `Muliaditan_2025_mab_mpbpk_human.R`, `Muliaditan_2025_mab_mpbpk_nhp.R`.

### delta_utfr_bbb (**canonical mPBPK luminal unbound-TfR deviation, BBB**)
- **Type:** compartment
- **Role:** Deviation of the luminal unbound TfR concentration on the BBB from its baseline `utfr0_bbb`, so the luminal level is `utfr_bbb_lum = utfr0_bbb + c_delta_utfr_bbb`. Driven by TfR synthesis/degradation, drug binding, and abluminal-TfR recycling.
- **Source aliases:** none.
- **Example models:** `Muliaditan_2025_mab_mpbpk_human.R`, `Muliaditan_2025_mab_mpbpk_nhp.R`.

### delta_utfr_bcsfb (**canonical mPBPK luminal unbound-TfR deviation, BCSFB**)
- **Type:** compartment
- **Role:** Deviation of the luminal unbound TfR concentration on the BCSFB from its baseline `utfr0_bcsfb`, so `utfr_bcsfb_lum = utfr0_bcsfb + c_delta_utfr_bcsfb`. BCSFB counterpart of `delta_utfr_bbb`.
- **Source aliases:** none.
- **Example models:** `Muliaditan_2025_mab_mpbpk_human.R`, `Muliaditan_2025_mab_mpbpk_nhp.R`.

### utfr_bbb_abl (**canonical mPBPK abluminal unbound TfR, BBB**)
- **Type:** compartment
- **Role:** Unbound TfR pool on the abluminal (brain) face of the BBB; consumed by binding abluminal drug to form `complex_bbb_abl` and regenerated by recycling at `krec_utfr`.
- **Source aliases:** none.
- **Example models:** `Muliaditan_2025_mab_mpbpk_human.R`, `Muliaditan_2025_mab_mpbpk_nhp.R`.

### utfr_bcsfb_abl (**canonical mPBPK abluminal unbound TfR, BCSFB**)
- **Type:** compartment
- **Role:** Unbound TfR pool on the abluminal (CSF) face of the BCSFB; consumed by binding abluminal drug to form `complex_bcsfb_abl` and regenerated by recycling at `krec_utfr`. BCSFB counterpart of `utfr_bbb_abl`.
- **Source aliases:** none.
- **Example models:** `Muliaditan_2025_mab_mpbpk_human.R`, `Muliaditan_2025_mab_mpbpk_nhp.R`.

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

### np (**canonical nanoparticle-conjugated species suffix**)
- **Type:** metabolite-suffix
- **Role:** Nanoparticle- / carrier-conjugated form of a drug in dual-species nanoparticle biodistribution models, where the bare canonical compartment name holds the released (free) drug and the `_np` suffix holds the still-conjugated drug travelling with its carrier.
- **Source aliases:** none.
- **Example models:** `Vasalou_2023_dendriticNanoparticle_mouse.R`, `Vasalou_2023_dendriticNanoparticle_rat.R`, `Vasalou_2023_dendriticNanoparticle_dog.R`, `Vasalou_2023_dendriticNanoparticle_human.R`.
- **Notes:** Registered for the library's first dual-species (carrier-bound + released) nanoparticle model. Use when a model carries the conjugated and the released drug through the *same* set of anatomical compartments, giving pairs such as `blood` / `blood_np`, `liver` / `liver_np`, `spleen` / `spleen_np`, `other` / `other_np`, and matching observation variables `Cc` / `Cc_np`. The bare-name-holds-free-drug orientation follows the rest of the library, where an unsuffixed state is the unbound / active species and a suffix marks a bound or carrier-associated form. Distinct from `nab` (neutralising antibody) and from the ADC payload suffixes (`mmae`, `dxd`, `sn38`), which name a specific chemical entity rather than the conjugation state; prefer a named payload suffix whenever the released moiety has a published INN, and reserve `np` for carriers whose API is unnamed. Operator-ratified 2026-07-29.

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

### 7ohmtx (**canonical 7-hydroxy-methotrexate suffix**)
- **Type:** metabolite-suffix
- **Role:** 7-hydroxy-methotrexate (7-OH-MTX), the principal circulating metabolite of methotrexate formed by aldehyde-oxidase-mediated 7-hydroxylation in the liver. Used in parent + metabolite joint popPK extractions of high-dose intravenous MTX therapy.
- **Source aliases:** none.
- **Example models:** `Joerger_2006_methotrexate.R` (3-cmt MTX parent + 2-cmt 7-OH-MTX metabolite, joint NONMEM ADVAN5 fit; metabolic fraction fixed at 10 percent per Joerger 2006 Results page 75).
- **Notes:** Suffix starts with a digit; the convention check matches on `endsWith(name, "_<metab>")` rather than treating the metabolite name as an R identifier, following the `3oh` / `7dm` precedent.

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
- **Notes:** Distinct from `m3g` (morphine-3-glucuronide) -- the suffix matcher uses `endsWith(name, "_m3")` vs `endsWith(name, "_m3g")` and these do not collide. Registered alongside the Svensson 2013 bedaquiline extraction (the first BDQ paper to model the M3 metabolite).

### m4 (**canonical edoxaban M4 metabolite suffix**)
- **Type:** metabolite-suffix
- **Role:** M4, the main human metabolite of edoxaban, formed by carboxylesterase-1-mediated hydrolysis and pharmacologically active. Used as the metabolite suffix on the `central_m4` compartment, the `lcl_m4` / `lvc_m4` parameters and the `Cc_m4` observation in joint parent + metabolite popPK models. In Jonsson 2015 the fraction of edoxaban converted to M4 was not identifiable, so the model assumes that all edoxaban not excreted renally forms M4 (`fm,M4 = 1 - fe = CLNR/CL`); the resulting `cl_m4` and `vc_m4` are therefore APPARENT values, overestimated relative to the true ones by a factor `1/fm,M4`.
- **Source aliases:** `M-4` -- the hyphenated form used in some edoxaban literature and regulatory documents.
- **Example models:** `Jonsson_2015_edoxaban.R` (doi:10.1002/jcph.541; joint edoxaban + M4 renal-impairment popPK with a urine compartment).
- **Notes:** Fourth member of the generic numbered-metabolite suffix family alongside `m1`, `m2` and `m3`, following the same convention of lowercasing a source paper's `M<n>` designation. Distinct from `m3g` / `m6g` (morphine glucuronides) and from `m2` / `m3` (bedaquiline demethylation products) -- the suffix matcher compares with `endsWith(name, "_m4")` so there is no collision. Do not confuse with `dm4` / `medm4`, which are the maytansinoid ADC payload suffixes.

### m8 (**canonical hydroxy-tert-butylamide (M8) suffix**)
- **Type:** metabolite-suffix
- **Role:** Hydroxy-tert-butylamide (M8) active metabolite of nelfinavir; formed by CYP2C19-mediated hydroxylation of nelfinavir and eliminated by CYP3A4. Equipotent to the parent drug against HIV-1 protease.
- **Source aliases:** none.
- **Example models:** `Hirt_2006_nelfinavir.R`.

### ko739 (**canonical KO-739 ziftomenib active metabolite suffix**)
- **Type:** metabolite-suffix
- **Role:** KO-739, one of the two active minor metabolites of the menin inhibitor ziftomenib (Kura Oncology development-code compound designation). Formed via CYP3A4-mediated biotransformation of ziftomenib in parallel with KO-516; reported to be less than 10% of total drug-related exposure in plasma. Used in parent + metabolite joint popPK extractions of ziftomenib.
- **Source aliases:** none.
- **Example models:** `Mitra_2026_ziftomenib.R` (sequential two-stage popPK; KO-739 central + one peripheral compartment; 1:1 in-vitro-anchored KO-739:KO-516 metabolic split; Table 1).
- **Notes:** Follows the paper-named metabolite suffix convention established by `m1` / `m2` / `m3` / `m8`. Kura development-code compound designation retained instead of a chemical-name shorthand because the paper does not disclose the chemical identity of the metabolite.

### ko516 (**canonical KO-516 ziftomenib active metabolite suffix**)
- **Type:** metabolite-suffix
- **Role:** KO-516, one of the two active minor metabolites of the menin inhibitor ziftomenib (Kura Oncology development-code compound designation). Formed via CYP3A4-mediated biotransformation of ziftomenib in parallel with KO-739; reported to be less than 10% of total drug-related exposure in plasma. Used in parent + metabolite joint popPK extractions of ziftomenib.
- **Source aliases:** none.
- **Example models:** `Mitra_2026_ziftomenib.R` (sequential two-stage popPK; KO-516 central + one peripheral compartment; 1:1 in-vitro-anchored KO-739:KO-516 metabolic split; Table 1).
- **Notes:** Follows the paper-named metabolite suffix convention established by `m1` / `m2` / `m3` / `m8`. Kura development-code compound designation retained instead of a chemical-name shorthand because the paper does not disclose the chemical identity of the metabolite.

### endox (**canonical endoxifen suffix**)
- **Type:** metabolite-suffix
- **Role:** Endoxifen (4-hydroxy-N-desmethyltamoxifen), major active metabolite of tamoxifen.
- **Source aliases:**
  - `endx` -- prior canonical name (pre-2026-06-19 readability standardization).
- **Example models:** `TerHeine_2014_tamoxifen.R`.
- **Notes:** Renamed from `endx` to `endox` on 2026-06-19 per the canonical-register standardization audit (operator decision: keep the contracted form but include the "o" so the suffix is readable as "endox" without vowel-stripping; `endx` was an opaque consonant cluster).

### canrenone (**canonical canrenone suffix**)
- **Type:** metabolite-suffix
- **Role:** Canrenone, the principal active metabolite of spironolactone and the species that actually occupies the mineralocorticoid receptor. Spironolactone is extensively metabolised; the canrenone concentration -- not the parent -- drives the Emax receptor-inhibition term in mineralocorticoid-receptor-antagonist PK/PD and QSP models, so canrenone carries its own central and peripheral disposition compartments alongside the parent.
- **Source aliases:**
  - `canrenone` / `canrenone_C2` -- the source-listing state names in Meid 2024 Appendix C, mapped to `central_canrenone` and `peripheral1_canrenone`.
- **Example models:** `Meid_2024_spironolactone_qsp.R` (potassium-homeostasis QSP; two-compartment canrenone disposition formed from spironolactone at `Spiro_Fmetabolized` = 0.19311 of parent clearance, driving mineralocorticoid-receptor occupancy).
- **Notes:** Full metabolite name kept rather than a contraction: it is short, unambiguous, and already the standard name in the spironolactone literature. Follows the same parent + metabolite pattern as `m1` / `endox` / `megx`.

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

### col (**canonical colistin metabolite / sibling-drug suffix**)
- **Type:** metabolite-suffix
- **Role:** Colistin as a suffixed species. Two use cases share the suffix: the active polymyxin generated in vivo by hydrolysis of the prodrug colistimethate sodium (CMS), and colistin as one member of a combination-antibiotic model.
- **Source aliases:** none.
- **Example models:** `LeuppiTaegtmeyer_2019_CMS.R` (DDMODEL00000295; colistin formed from CMS), `Mohamed_2016_colistin_meropenem.R` (colistin as a sibling drug alongside meropenem).
- **Notes:** Same token as the bare `col` drug-state compartment; both Types co-exist for the same canonical name. Merged 2026-08-06 from two separate `metabolite-suffix` blocks that had been added independently for the two use cases. Because the register resolves document-order last-writes-win, having both meant the CMS-prodrug description was silently discarded; a repeated name is only safe when the Types differ.

### cpq (**canonical carboxyprimaquine suffix**)
- **Type:** metabolite-suffix
- **Role:** Carboxyprimaquine, the major circulating metabolite of primaquine formed by monoamine-oxidase-A-mediated oxidative deamination of the terminal amine of the aminopentyl side chain. Conversion is 1:1 molar, which is why joint primaquine + carboxyprimaquine models are parameterised on the molar scale. Used in parent + metabolite simultaneous popPK extractions.
- **Source aliases:** none.
- **Example models:** `Wattanakul_2024_primaquine.R`, `Wattanakul_2024_primaquine_motherinfant.R` (one-compartment carboxyprimaquine disposition fed both by first-pass metabolism from the last transit compartment and by the entirety of systemic primaquine clearance, plus its own breast-milk compartment).

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

### eddp (**canonical EDDP methadone metabolite suffix**)
- **Type:** metabolite-suffix
- **Role:** 2-Ethylidene-1,5-dimethyl-3,3-diphenylpyrrolidine (EDDP), the principal (largely inactive) N-demethylation cyclisation metabolite of methadone. Formed stereoselectively by CYP3A, CYP2B6, and CYP2C19 from R- and S-methadone. Used as the metabolite suffix in parent + metabolite simultaneous popPK models of methadone enantiomers (compartments `central_eddp`, parameters `lcl_eddp`, residuals `propSd_eddp`).
- **Source aliases:** none.
- **Example models:** `Aruldhas_2021_R_methadone.R`, `Aruldhas_2021_S_methadone.R`.

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

### ast5902 (**canonical AST5902 furmonertinib active metabolite suffix**)
- **Type:** metabolite-suffix
- **Role:** AST5902, the CYP3A4-derived active metabolite of furmonertinib (AST2818 / alflutinib), a third-generation irreversible EGFR TKI approved for NSCLC with EGFR-sensitising / T790M mutations. Used as a compartment / parameter suffix in joint parent-plus-metabolite popPK models where furmonertinib and AST5902 are followed simultaneously in plasma.
- **Source aliases:** none.
- **Example models:** `Zou_2022_furmonertinib.R` (doi:10.1038/s41401-021-00798-y).

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

### oxy (**canonical oxypurinol suffix**)
- **Type:** metabolite-suffix
- **Role:** Oxypurinol (1H-purine-2,6,8(3H)-trione), the principal long-lived active metabolite of allopurinol produced by xanthine oxidase and aldehyde oxidase catalysed 2-hydroxylation. Oxypurinol is itself a xanthine-oxidase inhibitor and is the main contributor to sustained urate-lowering effect during chronic allopurinol therapy (allopurinol has a short plasma half-life of ~1-2 h while oxypurinol persists for ~24 h with a renally-cleared elimination). Used as a metabolite suffix in parent-plus-metabolite joint or sequential popPK models for allopurinol.
- **Source aliases:** none.
- **Example models:** `Wright_2013_allopurinol.R` (Wright 2013 doi:10.1007/s00228-013-1478-8 -- sequential 2-cpt allopurinol + 1-cpt oxypurinol model; the metabolite is central_oxy with Cc_oxy observation and independent renal + non-renal clearance components sensitive to CLcr, fat-free mass, and concomitant diuretic use).

### napa (**canonical N-acetylprocainamide suffix**)
- **Type:** metabolite-suffix
- **Role:** N-acetylprocainamide (NAPA), the major active metabolite of procainamide formed by hepatic N-acetylation (NAT2). NAPA is a Vaughan-Williams Class III antiarrhythmic in its own right, is predominantly renally cleared, and accumulates alongside the parent in chronic-kidney-disease patients -- parent-plus-NAPA combined exposure is the therapeutic-window target for procainamide dosing in CKD.
- **Source aliases:** none.
- **Example models:** `Mohamed_2013_procainamide.R` (doi:10.1053/j.ajkd.2013.02.358).

### norcloz (**canonical norclozapine (N-desmethylclozapine) suffix**)
- **Type:** metabolite-suffix
- **Role:** Norclozapine (N-desmethylclozapine), the primary pharmacologically active metabolite of clozapine formed predominantly by CYP1A2 (with secondary contributions from CYP2C19, CYP3A4, CYP2C9, and CYP2D6). Norclozapine retains receptor affinity at multiple monoaminergic and muscarinic targets and is routinely measured alongside clozapine in therapeutic-drug-monitoring (TDM) practice; the parent-to-metabolite ratio is itself a clinical descriptor of CYP1A2 activity.
- **Source aliases:** none.
- **Example models:** `Li_2012_clozapine.R` (doi:10.1038/aps.2012.71).

### dnef (**canonical desmethyl-nefopam suffix**)
- **Type:** metabolite-suffix
- **Role:** Desmethyl-nefopam (nor-nefopam), N-demethyl metabolite of nefopam formed by CYP3A4-mediated N-demethylation. Followed in plasma alongside the parent in joint parent + metabolite popPK models of postoperative nefopam analgesia. Used as the metabolite suffix on `central_dnef` compartments, `lcl_dnef` / `lvc_dnef` / `lcl_form_dnef` parameters, and `addSd_dnef` / `propSd_dnef` residual SDs.
- **Source aliases:** none.
- **Example models:** `Djerada_2014_nefopam.R` (doi:10.1111/bcp.12291).

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

### cns7054 (**canonical CNS 7054 remimazolam metabolite suffix**)
- **Type:** metabolite-suffix
- **Role:** CNS 7054, the pharmacologically inactive carboxylic-acid metabolite of remimazolam formed by carboxylesterase-1 (CES-1) hydrolysis of the parent's methyl-propanoate ester. It is the sole quantified metabolite of remimazolam, is cleared predominantly renally (70-90% of its elimination), and clears roughly 19-fold more slowly than the parent, so it accumulates well above parent concentrations during infusion. Drives `central_cns7054` / `peripheral1_cns7054`, the `lcl_cns7054` / `lvc_cns7054` / `lq_cns7054` / `lvp_cns7054` parameters, and the `propSd_cns7054` / `addSd_cns7054` residuals; parent-side parameters keep the canonical unsuffixed names.
- **Source aliases:** `CNS7054` and `CNS 7054` (both spellings appear in the source papers).
- **Example models:** `Chen_2024_remimazolam.R` (doi:10.3390/pharmaceutics16091122).
- **Notes:** Chen 2024 states no stoichiometric conversion factor for the parent-to-metabolite transfer, so the amount transfer is 1:1 on a mass basis and the fitted CNS 7054 volumes and clearance are apparent values in remimazolam-mass equivalents. The parent (439 g/mol) and metabolite (425.1 g/mol) molecular weights differ by only 3.2 percent, which is within the reported RSEs of the metabolite disposition parameters (3-5 percent).

---

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

### sbt (**canonical sulbactam sibling-drug suffix**)
- **Type:** metabolite-suffix
- **Role:** Sulbactam (beta-lactamase-inhibitor beta-lactam) sibling-drug suffix, used when sulbactam is co-modelled with a partner beta-lactam or beta-lactamase inhibitor and is not the unsuffixed parent. Drives `central_sbt` / `peripheral1_sbt` compartments, `lcl_sbt` / `lvc_sbt` / `lq_sbt` / `lvp_sbt` PK parameters, `f_renal_sbt` and the `e_<cov>_<param>_sbt` covariate-effect forms, and the `propSd_sbt` / `addSd_sbt` residuals on the `Cc_sbt` sulbactam plasma concentration.
- **Source aliases:**
  - `SUL` -- abbreviation used in `Cammarata_2024_sulbactam_durlobactam.R` (supplement Abbreviations Listing).
  - `(2)` -- NONMEM L2 drug index for sulbactam in `Soto_2014_ampicillin_sulbactam.R`.
- **Example models:**
  - `Cammarata_2024_sulbactam_durlobactam.R` (doi:10.1128/aac.00485-24) where durlobactam is the unsuffixed parent (NONMEM compartments 1-2) and sulbactam carries the suffix (NONMEM compartments 3-4).
  - `Soto_2014_ampicillin_sulbactam.R` (doi:10.1111/bcp.12232), where ampicillin is the unsuffixed parent.
- **Notes:** Registered 2026-07-28 alongside the Cammarata 2024 sulbactam-durlobactam extraction. `Soto_2014_ampicillin_sulbactam.R` had been using this suffix since its own extraction without a register entry, so registering it also clears that model's pre-existing `central_sbt` / `peripheral1_sbt` / `propSd_sbt` convention warnings. Sulbactam is a sibling drug, not a metabolite: in both source papers the two analytes are dosed as a fixed-ratio combination and fitted simultaneously without interconversion.

---

## Anthracycline cardioprotection sibling-drug suffixes

Sibling-drug suffix for anthracycline + cardioprotectant PK/TD models (e.g., Mody 2023 doxorubicin + dexrazoxane), where doxorubicin (the anthracycline responsible for the antitumor and cardiotoxic effect) is the unsuffixed parent and dexrazoxane (the FDA-approved cardioprotectant, iron chelator via ADR-925 and Top2 catalytic inhibitor) carries the `dex` suffix.

### dex (**canonical dexrazoxane sibling-drug suffix**)
- **Type:** metabolite-suffix
- **Role:** Dexrazoxane (ICRF-187, bisdioxopiperazine cardioprotectant) sibling-drug suffix, paired with doxorubicin as the unsuffixed parent in anthracycline + cardioprotectant PK/TD models. Drives `central_dex` / `peripheral1_dex` PK compartments, `lcl_dex` / `lvc_dex` / `lq_dex` / `lvp_dex` PK parameters, `limax_dex` / `lic50_dex` PD-inhibition parameters (DEX's Hill inhibition on the DOX-induced stimulation of cardiomyocyte death), and the `addSd_dex` residual on plasma concentration.
- **Source aliases:** none (the paper uses DEX throughout; the FDA-approved brand names Zinecard / Cardioxane / Totect refer to the same molecule).
- **Example models:** `Mody_2023_doxorubicin_dexrazoxane.R` (doi:10.1038/s41598-023-29964-4).

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

### tfv (**canonical tenofovir suffix**)
- **Type:** metabolite-suffix
- **Role:** Tenofovir (TFV), the pharmacologically relevant plasma moiety formed by complete, irreversible intracellular hydrolysis of the ester prodrug tenofovir alafenamide (TAF). Used as the metabolite suffix on `central_tfv` compartments, `lcl_tfv` / `lvc_tfv` parameters, and the `Cc_tfv` observation in joint prodrug + active-moiety popPK models where tenofovir alafenamide itself is carried as a state and therefore keeps the canonical unsuffixed `central` / `Cc` names.
- **Source aliases:** `TFV` (the standard antiretroviral-literature abbreviation, used by Thoueille 2023 in `CL_TFV`, `V_TFV`, `sigma_addTFV`).
- **Example models:** `Thoueille_2023_tenofovir_alafenamide.R`.
- **Notes:** Operator-ratified sidecar 2026-07-29 (oare_PMC10232258 request-001 Q1, option A). Apply this suffix ONLY when the prodrug is itself a modelled compartment, per the standing "the parent / dosed species always wins canonical naming" rule; models that dose into a tenofovir compartment directly, without carrying a tenofovir-disoproxil-fumarate or tenofovir-alafenamide state, keep tenofovir as the bare canonical `central` / `Cc` (existing precedents: `Baheti_2011_tenofovir.R`, `Chen_2016_tenofovir_emtricitabine.R`, and the two tenofovir-alone siblings `Thoueille_2023_tenofovir_full.R` / `Thoueille_2023_tenofovir_reduced.R`). Consequently `central` denotes tenofovir alafenamide in `Thoueille_2023_tenofovir_alafenamide.R` but tenofovir in that paper's other two model files -- an intended consequence of the parent-wins rule, not an inconsistency. Distinct from [[tfvdp]] (tenofovir diphosphate, the intracellular active anabolite) and from `taf`, which is deliberately NOT registered. STRING-COLLISION NOTE: `tfv` also appears as the tail of the unrelated compartment `brain_csf_tfv`, where it abbreviates "third + fourth ventricle" in rat intra-brain SBPK models. That makes `brain_csf_tfv` parseable both as a canonical compartment and as `brain_csf` + `_tfv` metabolite suffix by the `endsWith(name, "_<metab>")` test in `checkModelConventions()`; the collision is harmless because `brain_csf` is itself canonical and so the compartment passes the check either way.

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

### val (**canonical valsartan sibling-drug suffix**)
- **Type:** metabolite-suffix
- **Role:** Valsartan (angiotensin-II receptor blocker) sibling-drug suffix, paired with amlodipine as the unsuffixed parent in the Heo 2016 amlodipine + valsartan combined PK/PD antihypertensive-interaction model. Drives `central_val` / `peripheral1_val` / `effect_val` compartments, `lcl_val` / `lvc_val` / `lq_val` / `lvp_val` / `ld1_val` PK parameters, and the `propSd_val` / `addSd_val` residuals on valsartan plasma concentration.
- **Source aliases:** none.
- **Example models:** `Heo_2016_amlodipine_valsartan.R` (doi:10.1111/bcp.13082).

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

### dap (**canonical daptomycin drug suffix**)
- **Type:** metabolite-suffix
- **Role:** Daptomycin drug suffix in combination antibiotic models.
- **Source aliases:** none.
- **Example models:** Combination antibiotic PD models with daptomycin.
- **Notes:** Same token as the bare `dap` drug-state compartment.

### cipro (**canonical ciprofloxacin metabolite / sibling-drug suffix**)
- **Type:** metabolite-suffix
- **Role:** Ciprofloxacin species suffix. Ciprofloxacin is both a fluoroquinolone antibacterial in its own right and the principal active metabolite of enrofloxacin (formed by hepatic N-dealkylation), so the suffix serves parent + metabolite extractions in which enrofloxacin is the parent (`central_cipro`, `lvc_cipro`, `lcl_cipro`, `propSd_cipro`, `Cc_cipro`) and would equally serve a combination model in which ciprofloxacin is a sibling drug.
- **Source aliases:** `Cipro`, `CIP` (Foster 2023 Figure S2 and Table 3 printed forms).
- **Example models:** `Foster_2023_enrofloxacin_ciprofloxacin_cat.R` (founding example; enrofloxacin -> ciprofloxacin formation clearance in cats with reduced kidney function).
- **Notes:** Spelled `cipro` rather than `cip` to stay unambiguous against the covariate register's `CONMED_CIP` and against the `ccip` bath-concentration compartment of the Rees 2018 hollow-fiber meropenem + ciprofloxacin model, which is a distinct state (a dosed medium concentration in a time-kill experiment), not a metabolite species suffix.

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

### oselcarb (**canonical oseltamivir carboxylate suffix**)
- **Type:** metabolite-suffix
- **Role:** Oseltamivir carboxylate (OC), the active neuraminidase-inhibitor metabolite formed from the oseltamivir prodrug primarily via human carboxylesterase 1 (HCE1) in the liver.
- **Source aliases:**
  - `OC` -- the abbreviation used throughout the oseltamivir literature (Chairat 2016, Kamal 2013, Standing 2012) and retained in prose / plot labels.
- **Example models:** `Standing_2012_oseltamivir.R`, `Chairat_2016_oseltamivir.R`, `Kamal_2013_oseltamivir.R`.
- **Notes:** Renamed from the former `oc` suffix on 2026-07-26 to remove the `OC`/`oc` case-only distinction from osteocalcin. Used as `central_oselcarb`, `cl_oselcarb`, `lcl_oselcarb`, `lvc_oselcarb`, `Cc_oselcarb`, following the `enaat` metabolite-suffix pattern.

### cloca (**canonical clopidogrel carboxylic acid suffix**)
- **Type:** metabolite-suffix
- **Role:** Clopidogrel carboxylic acid (CLO-CA), the pharmacologically inactive hydrolysis product formed from clopidogrel by carboxylesterase-1 (CES1). It is the major circulating clopidogrel species -- roughly 85-90% of the absorbed dose is routed to it, and plasma concentrations are about 2000-fold higher than those of the parent, so it is assayed in ug/mL where clopidogrel is assayed in ng/mL. Distinct from `h4`, the active thiol metabolite of the same parent drug, which sits on the competing CYP-mediated oxidation branch; a joint model can carry both suffixes.
- **Source aliases:** `CLO-CA`, `SR 26334` (the Sanofi development code for the carboxylic acid metabolite), `carbo` (Jung 2024's abbreviation, used in `fm_carbo` / `sigma_carbo`), `m2` (Jung 2024 Table 2 parameter subscripts `V_m2`, `V_p2`, `CL_m2`, `Q_m2`, a positional index rather than an analyte name).
- **Example models:** `Pejcic_2024_clopidogrel.R` (`central_cloca`, `peripheral1_cloca`, `Cc_cloca`, `lcl_cloca`, `lvc_cloca`, `lq_cloca`, `lvp_cloca`; doi:10.3390/pharmaceutics16050685), `Jung_2024_clopidogrel.R` (same suffix set, carried alongside `central_h4` in a joint parent + both-metabolites model; doi:10.1002/psp4.13053).
- **Notes:** The generic `acid` suffix is unavailable for this analyte -- it is already the canonical simvastatin acid suffix -- and a bare carboxylate abbreviation is a poor canonical anyway (the `oc` -> `oselcarb` rename records exactly that collision pressure), so the drug-prefixed `cloca` form is used.

### enaat (**canonical enalaprilat metabolite suffix**)
- **Type:** metabolite-suffix
- **Role:** Enalaprilat (ENAAT), the pharmacologically active diacid metabolite of enalapril formed by hepatic carboxylesterase 1 (CES1) hydrolysis of the inactive ester prodrug. Used as the metabolite suffix on `central_enaat` compartments, `lcl_enaat` / `lvc_enaat` parameters, and `Cc_enaat` observation in joint parent-prodrug + active-metabolite popPK models.
- **Source aliases:** `ENAAT` (NONMEM-style abbreviation used by Steichert 2025).
- **Example models:** `Steichert_2025_enalapril_enalaprilat_pediatric.R`, `Luo_2024_enalapril_pbpk.R`.

### benat (**canonical benazeprilat metabolite suffix**)
- **Type:** metabolite-suffix
- **Role:** Benazeprilat, the pharmacologically active diacid metabolite of the ester prodrug benazepril, formed by hepatic carboxylesterase 1 (CES1) hydrolysis and eliminated renally. Drives `central_benat` / `peripheral1_benat` / `liver_benat` / `kidney_benat` compartments, `lvc_benat` / `lcl_int_k_benat` / `kp_liver_benat` parameters and the `Cc_benat` plasma-concentration observation.
- **Source aliases:** none.
- **Example models:** `Luo_2024_benazepril_pbpk.R`.
- **Notes:** Registered 2026-08-03 alongside the Luo 2024 CES1 semi-PBPK family. Follows the `enaat` (enalaprilat) precedent: drug stem plus `at` for the active diacid of an ACE-inhibitor ester prodrug.

### cilat (**canonical cilazaprilat metabolite suffix**)
- **Type:** metabolite-suffix
- **Role:** Cilazaprilat, the pharmacologically active diacid metabolite of the ester prodrug cilazapril, formed by hepatic CES1 hydrolysis and eliminated renally.
- **Source aliases:** none.
- **Example models:** `Luo_2024_cilazapril_pbpk.R`.
- **Notes:** Registered 2026-08-03 alongside the Luo 2024 CES1 semi-PBPK family, per the `enaat` precedent.

### temat (**canonical temocaprilat metabolite suffix**)
- **Type:** metabolite-suffix
- **Role:** Temocaprilat, the pharmacologically active diacid metabolite of the ester prodrug temocapril, formed by hepatic CES1 hydrolysis. Unusually among the ACE-inhibitor diacids it is eliminated by BOTH renal and MRP2-mediated biliary routes, so it also drives a `lcl_int_bile_temat` biliary intrinsic-clearance parameter.
- **Source aliases:** none.
- **Example models:** `Luo_2024_temocapril_pbpk.R`.
- **Notes:** Registered 2026-08-03 alongside the Luo 2024 CES1 semi-PBPK family, per the `enaat` precedent.

### perat (**canonical perindoprilat metabolite suffix**)
- **Type:** metabolite-suffix
- **Role:** Perindoprilat, the pharmacologically active diacid metabolite of the ester prodrug perindopril, formed by hepatic CES1 hydrolysis and eliminated renally. Note that perindopril's parallel UGT glucuronidation route does NOT form perindoprilat, so a model carrying this suffix must keep the two hepatic pathways separate.
- **Source aliases:** none.
- **Example models:** `Luo_2024_perindopril_pbpk.R`.
- **Notes:** Registered 2026-08-03 alongside the Luo 2024 CES1 semi-PBPK family, per the `enaat` precedent.

### oselc (**canonical oseltamivir carboxylate metabolite suffix**)
- **Type:** metabolite-suffix
- **Role:** Oseltamivir carboxylate, the antivirally active carboxylate metabolite of the ester prodrug oseltamivir, formed by hepatic CES1 hydrolysis and eliminated renally by combined glomerular filtration and tubular secretion.
- **Source aliases:** `OC` -- the abbreviation used by Luo 2024 Table 9 and Section 3.1.6.
- **Example models:** `Luo_2024_oseltamivir_pbpk.R`.
- **Notes:** Registered 2026-08-03 alongside the Luo 2024 CES1 semi-PBPK family. Deliberately NOT named `oc`: the former `OC` canonical was retired on 2026-07-26 precisely because it was overloaded across five unrelated concepts, oseltamivir carboxylate among them (see the `osteocalcin` entry). `oselc` keeps the drug stem explicit so the collision cannot recur.

### ppf (**canonical propofol active-metabolite suffix**)
- **Type:** metabolite-suffix
- **Role:** Propofol (2,6-diisopropylphenol), the active sedative-hypnotic metabolite liberated from the water-soluble prodrug fospropofol (GPI 15715, AQUAVAN) via systemic alkaline-phosphatase hydrolysis. Used as the metabolite suffix on `central_ppf` / `peripheral1_ppf` compartments and `Cc_ppf` observation in joint parent-prodrug + active-drug popPK models.
- **Source aliases:** `PR` (Gibiansky 2005 poster table headings: Vc_PR, K10_PR, K12_PR, K21_PR, K_GPI-PR).
- **Example models:** `Gibiansky_2005_fospropofol.R`.

### ceftaroline (**canonical ceftaroline active-metabolite suffix**)
- **Type:** metabolite-suffix
- **Role:** Ceftaroline, the active anti-MRSA cephalosporin liberated from the water-soluble N-phosphono prodrug ceftaroline fosamil by systemic phosphatase hydrolysis. Used as the metabolite suffix on `central_ceftaroline` / `peripheral1_ceftaroline` / `peripheral2_ceftaroline` compartments, `lcl_ceftaroline` / `lvc_ceftaroline` / `lq_ceftaroline` / `lvp_ceftaroline` / `lq2_ceftaroline` / `lvp2_ceftaroline` parameters, and the `Cc_ceftaroline` observation in joint parent-prodrug + active-metabolite popPK models. The prodrug ceftaroline fosamil is the dosed species and is itself a modelled compartment, so per the standing "the parent / dosed species always wins canonical naming" rule it keeps the unsuffixed `central` / `peripheral1` / `Cc` names.
- **Source aliases:**
  - `C` -- the subscript Riccobene 2016 Supplemental Table 1 uses to distinguish ceftaroline (CLc, Vcc, Q1c, Vp1c, Q2c, Vp2c) from ceftaroline fosamil (CLcf, Vccf, Qcf, Vpcf), against the prodrug's `cf` subscript.
  - `CPT` -- the abbreviation used in the ceftaroline literature (ceftaroline fosamil is `CPT-F`); deliberately NOT adopted as the suffix token because `CPT` also abbreviates camptothecin.
- **Example models:** `Riccobene_2016_ceftaroline.R`.
- **Notes:** Spelled out rather than abbreviated, following the `sunitinib` / `irinotecan` precedent, because every short form collides with an unrelated agent (`cpt` with camptothecin, `cef`/`cft` with the rest of the cephalosporin class, of which nlmixr2lib already carries more than twenty). Models that dose a ceftaroline compartment directly without carrying a ceftaroline fosamil state would keep ceftaroline as the bare canonical `central` / `Cc`, exactly as documented for [[tfv]].

---

## PBPK organ extracellular / cellular split (Parmar 2023 mPBPK family)

Minimal-PBPK models that resolve each organ into a rapid-equilibrium pool (vascular + interstitial, taken to be in instantaneous equilibrium with blood) and a slow-equilibrium cellular pool use the `<organ>_extracellular` / `<organ>_cellular` suffix pair. The two suffixes are complements of each other and must be used together: an organ carrying one without the other is a naming error, because the pair's whole point is that the organ's total concentration is the amount-weighted sum of exactly these two pools divided by the organ volume.

The names are spelled out in full per the 2026-05-28 anti-abbreviation audit. In particular, do **not** shorten the rapid-equilibrium pool to `<organ>_rap`: that would read as a per-organ analogue of the existing `a_rapidly_perfused` compartment, which means something different (a lumped rapidly-perfused ORGAN GROUP, not an intra-organ sub-pool).

### extracellular (**canonical rapid-equilibrium vascular+interstitial pool suffix**)
- **Type:** metabolite-suffix
- **Role:** Rapid-equilibrium sub-compartment of the named organ, lumping the vascular and interstitial spaces which are assumed to be in instantaneous equilibrium with blood. Its effective volume is `V_vascular + V_interstitial / k(b/p)` so that a single blood-scale concentration describes both spaces.
- **Source aliases:** `rapid equilibrium sub-compartment`, `V+I` -- Parmar 2023 Table 4 / Table 7 row descriptions and supplement equations S3, S7, S10, S13, S16 (free-form prose labels, not data column names).
- **Example models:** `Parmar_2023_spectinamide_1599_mouse_pbpk.R`, `Parmar_2023_spectinamide_1599_rat_pbpk.R`, `Parmar_2023_spectinamide_1810_mouse_pbpk.R`, `Parmar_2023_spectinamide_1810_rat_pbpk.R`.
- **Notes:** Distinct from the membrane-limited `is_<organ>` interstitial-space prefix in `pbpkSubCompartmentRegex`, which carries interstitium alone with the vascular space held as a separate `bc_` / `vp_` state; `<organ>_extracellular` deliberately lumps the two because the source model assumes instantaneous vascular-interstitial equilibrium. Paired with `cellular`.

### cellular (**canonical slow-equilibrium cellular pool suffix**)
- **Type:** metabolite-suffix
- **Role:** Slow-equilibrium cellular sub-compartment of the named organ, coupled to the organ's `_extracellular` pool by a first-order influx acting on the unbound fraction (`K(I->C) x fu`) and a first-order back flux (`K(C->I)`) on the total cellular amount.
- **Source aliases:** `cellular sub-compartment` -- Parmar 2023 Table 4 / Table 7 row descriptions and supplement equations S4, S8, S11, S14, S17 (prose label, not a data column name).
- **Example models:** `Parmar_2023_spectinamide_1599_mouse_pbpk.R`, `Parmar_2023_spectinamide_1599_rat_pbpk.R`, `Parmar_2023_spectinamide_1810_mouse_pbpk.R`, `Parmar_2023_spectinamide_1810_rat_pbpk.R`.
- **Notes:** Distinct from the membrane-limited `int_<organ>` intracellular prefix in `pbpkSubCompartmentRegex`, which pairs with `is_<organ>` under a permeability-surface-area coupling rather than the pair of first-order rate constants used here. Paired with `extracellular`.

---

## Epithelial lining fluid

### elf (**canonical epithelial-lining-fluid compartment**)
- **Type:** compartment
- **Role:** Epithelial lining fluid of the lung airspace: the apical aqueous layer that inhaled or intratracheally instilled drug is deposited into before it distributes into the lung tissue and reaches the systemic circulation. Carries its own physiological volume and its own unbound fraction `fu_elf` (typically derived from plasma `fu` and the plasma/ELF albumin ratio), and exchanges with the lung's `lung_cellular` and `lung_extracellular` pools.
- **Source aliases:** `ELF` -- Parmar 2023 Figure 2, Table 2 (`V ELF`, `fu ELF`) and supplement equations S22-S26.
- **Example models:** `Parmar_2023_spectinamide_1599_mouse_pbpk.R`.
- **Notes:** A physiological airspace compartment, not an absorption depot: the inhaled or intratracheal dose lands in a separate `depot2` and reaches `elf` by first-order absorption with its own bioavailability, exactly as an oral dose reaches `central` from `depot`. Whole-lung concentration includes the ELF amount (Parmar 2023 supplement S26). Related but distinct from `isf` (generic interstitial fluid) and from the `brain_csf*` namespace; a future inhalation model that resolves regional airway ELF should register `elf_<region>` names rather than overload the bare `elf`.

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

---

## PBPK permeability-limited tissue subcompartment suffixes (Gaohua 2023)

Permeability-limited whole-body PBPK subcompartment suffixes. Each tissue carries four subcompartments -- residual blood cells, residual plasma, extracellular water and intracellular water -- with passive permeation between adjacent pairs and active uptake / efflux transporters on the cell membrane between the two water spaces. The same `_plasma` / `_bc` pair also describes the blood compartments (`venous_plasma`, `arterial_plasma`, `portal_plasma` and their `_bc` partners), which have no extracellular or intracellular water. This is the inverse-ordering analogue of the `<subtype>_<organ>` `pbpkSubCompartmentRegex` and extends the `<organ>_<suffix>` shape founded by Ayyar 2024, so the two families compose. Registered per the operator naming decision of 2026-08-05.

### plasma (**canonical residual-plasma subcompartment suffix**)
- **Type:** metabolite-suffix
- **Role:** Residual plasma within the named tissue (the paper's "tissue plasma", TP) -- the plasma fraction of the blood left in the tissue after bleeding. For a blood compartment it is that compartment's plasma space. Perfused at the compartment's plasma flow and exchanging with `<organ>_bc` and `<organ>_ew`.
- **Source aliases:** `TP`, `tp`, `plasma_<organ>`.
- **Example models:** `Gaohua_2023_permeabilityLimited_pbpk.R`.

### bc (**canonical residual-blood-cell subcompartment suffix**)
- **Type:** metabolite-suffix
- **Role:** Residual blood cells within the named tissue (the paper's "tissue blood cells", TC), or a blood compartment's blood-cell space. Perfused at the compartment's blood-cell flow and exchanging only with `<organ>_plasma`. Echoes the existing `bcc` (central blood cells) and the `bc_<organ>` prefix form in `pbpkSubCompartmentRegex`.
- **Source aliases:** `TC`, `tc`, `rbc_<organ>`, `RBC`.
- **Example models:** `Gaohua_2023_permeabilityLimited_pbpk.R`.

### ew (**canonical extracellular-water subcompartment suffix**)
- **Type:** metabolite-suffix
- **Role:** Extracellular water of the named tissue, as implemented in the Simcyp Simulator tissue-composition scheme. Sits between the residual plasma and the intracellular water; the cell membrane separating it from `<organ>_iw` is where active uptake and efflux transporters act. Distinct from `is_<organ>` (interstitial space) in `pbpkSubCompartmentRegex`, which is an antibody-distribution concept paired with endosomal FcRn recycling rather than a small-molecule tissue-composition water space.
- **Source aliases:** `EW`.
- **Example models:** `Gaohua_2023_permeabilityLimited_pbpk.R`.

### iw (**canonical intracellular-water subcompartment suffix**)
- **Type:** metabolite-suffix
- **Role:** Intracellular water of the named tissue. Reached from `<organ>_ew` by passive permeation and by active uptake, and is the usual site of intracellular metabolism. Distinct from `int_<organ>` in `pbpkSubCompartmentRegex`, which denotes a bulk intracellular pool in mRNA-LNP models rather than a tissue-composition water space.
- **Source aliases:** `IW`.
- **Example models:** `Gaohua_2023_permeabilityLimited_pbpk.R`.

### portal (**canonical bare portal-vein blood compartment**)
- **Type:** compartment
- **Role:** Portal vein blood compartment, collecting the venous outflow of the splanchnic organs (pancreas, spleen, gut) and delivering it to the liver alongside the hepatic artery. A blood compartment, so it carries only the `_plasma` and `_bc` subcompartments. Complements the already-canonical bare `venous` and `arterial` roots and the `vp_portal` vascular-concentration form; the bare root is registered so `portal_plasma` / `portal_bc` compose with the subcompartment suffixes above.
- **Source aliases:** `PV`, `pv`.
- **Example models:** `Gaohua_2023_permeabilityLimited_pbpk.R`.
