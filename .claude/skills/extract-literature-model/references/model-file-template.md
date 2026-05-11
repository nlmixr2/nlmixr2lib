# Model file template

Copy this skeleton into `inst/modeldb/<category>/<FirstAuthor>_<Year>_<drug>.R` and fill in. Placeholders are wrapped in `<>`. Inline comments that start with `#!` are instructions for the person filling the template and must be deleted before committing.

Related references:
- `naming-conventions.md` — parameter / compartment / IIV naming, plus the NONMEM → nlmixr2 syntax-translation section
- `covariate-columns.md` — required register of covariate column names (consult before inventing any new name)
- `verification-checklist.md` — run through after first pass

```r
<FirstAuthor>_<Year>_<drug> <- function() {
  description <- "<One-sentence summary, e.g., 'Two-compartment population PK model for <drug> in <population>'>"
  #! Non-human studies: PREPEND the species to description, e.g.,
  #!   description <- "Preclinical (rat, Sprague-Dawley). Two-compartment popPK model for ..."
  #!   description <- "In vitro (SKBR3 cell line). Mechanistic HER2 receptor trafficking model."
  #! so the species is visible in modellib() listings without inspecting population$species.
  #! For preclinical/in vitro models: filename should also carry a species suffix
  #!   <FirstAuthor>_<Year>_<drug>_<species>.R  (e.g., Geldof_2008_fluvoxamine_rat.R)
  #! to keep human and animal extractions of the same drug as separate, non-colliding files.
  #!
  #! PBPK / QSP / MBMA / mechanism-class models: PREPEND the class label, e.g.,
  #!   description <- "PBPK (whole-body, SimCYP V12.1). APAP disposition in pregnancy ..."
  #!   description <- "QSP. 5-lipoxygenase pathway with 33 ODEs, 113 parameters ..."
  #!   description <- "MBMA. Trastuzumab + T-DXd from 103 trials; between-study variability only ..."
  #! and use a class suffix on the filename when a popPK extraction of the same paper /
  #! drug already exists:
  #!   <FirstAuthor>_<Year>_<drug>_pbpk.R  /  _qsp.R  /  _mbma.R
  #!
  #! CRITICAL for PBPK / QSP / MBMA: every parameter in ini() must come from the paper
  #! text, an on-disk supplement, or an on-disk upstream paper. Never substitute from
  #! training-data class-typical values or platform defaults (SimCYP / GastroPlus / OSP
  #! built-in libraries). When a parameter is missing in all on-disk sources, use the
  #! Phase 4 missing-parameter sidecar. See SKILL.md Phase 1 step 3a.
  reference   <- "<Full citation with DOI>"
  vignette    <- "<FirstAuthor>_<Year>_<drug>"  #! basename of vignettes/articles/<FirstAuthor>_<Year>_<drug>.Rmd; used by list-of-models to link here
  units       <- list(time = "<day|hour>", dosing = "<mg|mg/kg>", concentration = "<ug/mL|ng/mL>")

  covariateData <- list(
    #! One entry per covariate. Canonical names come from inst/references/covariate-columns.md.
    #! If the source paper uses a different column name, set source_name and, if values must be
    #! transformed (e.g., SEXM -> SEXF), note the transformation in notes and confirm with user.
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying; used for allometric scaling with reference weight <ref_wt> kg.",
      source_name        = "WT"
    )
    #! Add further covariates here.
  )

  population <- list(
    #! Common fields, all optional except n_subjects. Units and wording are not fixed —
    #! match the source. Pick whichever age/weight units match the study population.
    #! Example A (adult study):
    #!   age_range     = "18-75 years"
    #!   age_median    = "52 years"
    #!   weight_range  = "50-120 kg"
    #!   weight_median = "78 kg"
    #!   disease_state = "moderate-to-severe atopic dermatitis"
    #!   dose_range    = "150-300 mg SC Q2W"
    #! Example B (pediatric study):
    #!   age_range     = "0-24 months"
    #!   age_median    = "3 months"
    #!   weight_range  = "2-12 kg"
    #!   weight_median = "5 kg"
    #!   disease_state = "healthy infants at risk for RSV"
    #!   dose_range    = "50-200 mg IM single dose"
    #! Additional keys welcome (ga_range, renal_function, hepatic_function, co_medication, ...).
    species        = "<required: 'human', 'rat (Sprague-Dawley)', 'mouse (HBCx-9 PDX)', 'beagle dog', 'in vitro (SKBR3 cell line)', 'human + rat', etc.>",
    n_subjects     = <integer>,
    n_studies      = <integer>,
    age_range      = "<adult: '18-75 years' | pediatric: '0-24 months'>",
    age_median     = "<adult: '52 years' | pediatric: '3 months'>",
    weight_range   = "<adult: '50-120 kg' | pediatric: '2-12 kg'>",
    weight_median  = "<adult: '78 kg' | pediatric: '5 kg'>",
    sex_female_pct = <numeric>,
    race_ethnicity = c(White = <pct>, Black = <pct>, Asian = <pct>, Other = <pct>),
    disease_state  = "<e.g., 'moderate-to-severe atopic dermatitis' or 'healthy infants at risk for RSV'>",
    dose_range     = "<e.g., '150-300 mg SC Q2W' or '50-200 mg IM single dose'>",
    regions        = "<e.g., North America, EU>",
    notes          = "<free text; cite the Table in the source that lists baseline demographics>"
  )

  ini({
    #! Every parameter has:
    #!   1. A log-transform prefix where applicable.
    #!   2. A label() with units and a plain-English description.
    #!   3. A trailing in-file comment pointing to the source location the value came from.
    #! Example: lcl <- log(0.0388); label("CL for 70 kg adult (L/day)")  # Clegg 2024 Table 3
    #!
    #! Wrap fixed values in fixed(): if the source paper holds the parameter
    #! constant rather than estimating it, the value belongs in fixed(...).
    #! See SKILL.md § "Fixed parameters in ini()" for the full source-signal
    #! checklist (NONMEM FIX flags, "fixed at <value>" prose, allometric
    #! exponents reported without uncertainty, F1=1 anchors, parameters
    #! inherited from upstream papers, etc.).

    # Structural parameters — reference values for <reference weight / age>
    lka     <- log(<value>); label("<description with units>")  # <source location>
    lcl     <- log(<value>); label("<description with units>")  # <source location>
    lvc     <- log(<value>); label("<description with units>")  # <source location>
    # Add lvp, lq, lfdepot, etc. as applicable
    # If F1 was fixed to 1: lfdepot <- fixed(log(1)); label("Bioavailability") # <source location>

    # Allometric / maturation parameters (if applicable)
    # Estimated:           allo_cl <- <value>; label("Allometric exponent on CL (unitless)")
    # Fixed at canonical:  allo_cl <- fixed(0.75); label("Allometric exponent on CL (unitless)")
    allo_cl <- <value>; label("Allometric exponent on CL (unitless)")  # <source location>

    # Covariate effects — one per covariate/parameter combination
    # Wrap in fixed() if the paper held this coefficient constant.
    e_<cov>_<param> <- <value>; label("<Effect of <cov> on <param>>")  # <source location>

    # IIV — eta + transformed parameter name. Use a block for correlated IIV.
    # Use ~ fixed(<var>) for IIVs the paper held constant from a prior publication.
    # Use fixed(0) inside a c(...) block for off-diagonals NONMEM fixed to zero.
    etalcl + etalvc ~ c(<var_cl>, <cov_cl_vc>, <var_vc>)  # <source location>
    etalka          ~ <var_ka>                              # <source location>

    # Residual error
    # Wrap in fixed() if NONMEM $SIGMA had FIX (LTBS pattern often does this).
    propSd <- <value>; label("Proportional residual error (fraction)")  # <source location>
    # addSd  <- <value>; label("Additive residual error (<units>)")     # uncomment if combined
  })

  model({
    #! Order inside model():
    #!   1. Derived covariate terms (maturation, race/season multipliers)
    #!   2. Individual PK parameters (exp(lX + etaX) * scaling)
    #!   3. Micro-constants (kel, k12, k21)
    #!   4. ODE system (d/dt(...))
    #!   5. Bioavailability / lag-time adjustments
    #!   6. Observation variable (Cc) and error model

    # 1. Derived terms
    # maturation_cl <- 1 - (1 - beta_cl) * exp(-age_term * log(2) / t50_cl)
    # race_cl       <- 1 + e_race_black_cl * RACE_BLACK

    # 2. Individual parameters
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) * (WT / <ref_wt>)^allo_cl  # * maturation_cl * race_cl ...
    vc <- exp(lvc + etalvc) * (WT / <ref_wt>)^allo_v

    # 3. Micro-constants (if using explicit ODEs)
    kel <- cl / vc
    # k12 <- q / vc
    # k21 <- q / vp

    # 4. ODE system
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central
    # d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # 5. Bioavailability
    # f(depot) <- exp(lfdepot)

    # 6. Observation and error
    Cc <- central / vc
    Cc ~ prop(propSd)
    # Cc ~ add(addSd) + prop(propSd)
  })
}
```

## Notes

- When the source paper uses `linCmt()`-compatible structure (simple 1/2/3-compartment with first-order or IV absorption), an equally valid body is:
  ```r
  model({
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) * (WT / <ref_wt>)^allo_cl
    vc <- exp(lvc + etalvc) * (WT / <ref_wt>)^allo_v
    vp <- exp(lvp) * (WT / <ref_wt>)^allo_v
    q  <- exp(lq)  * (WT / <ref_wt>)^allo_cl
    Cc <- linCmt()
    Cc ~ prop(propSd)
  })
  ```
  Prefer explicit ODEs if the source includes mechanism-specific structure (TMDD, transit chains, enterohepatic cycling).

- Helper functions (e.g., to derive one-hot race indicators from a raw race column) live in `R/` and are exported with roxygen. Name them `<modelname>_<action>()`, e.g., `Clegg_2024_nirsevimab_derive_race_indicators()`. Only add if the covariate encoding is non-trivial.

- Before committing: delete every `#!` instruction comment. Leave the source-location `# <source>` comments — they are the audit trail.

## Documenting non-paper-derived parameter values

When a parameter value did not come from the paper's text or tables — e.g. the
operator read it off a graphical figure, an author supplied it via email
correspondence, or it was carried from an upstream-task model file — record the
provenance inline on the parameter line so the source-trace is unambiguous:

```r
ini({
  # paper-derived (Table 4)
  lcl <- log(0.225) ; label("Clearance (L/day)")

  # author correspondence (J. Almquist email 2026-04-29);
  # see tracking/operator_followups.md F12
  lkdeg <- log(0.0231) ; label("Receptor degradation rate (1/day)")

  # operator-extracted from Figure 2 (digitised); ±10% uncertainty
  lvdxd <- log(0.038) ; label("DXd payload volume (L/kg) — figure-derived")
})
```

Mirror the call-out in the vignette's Assumptions and deviations section (see
`vignette-template.md`).

## `depends_on` lineage in the `reference` field

When a model's PK structure is fixed from another paper already in nlmixr2lib
(typically a `depends_on` task), record both citations so the dependency is
visible in `?modellib`, `checkModelConventions()` output, and `readModelDb()`
metadata — not only in the task YAML:

```r
reference <- paste(
  "Iwaki et al. (2021) J Clin Pharmacol 61(11):1488-1499.",
  "doi:10.1002/jcph.1953.",
  "PK structure adapted from Yao et al. (2018)",
  "J Clin Pharmacol 58(11):1488-1497; see",
  "modellib('Yao_2018_guselkumab').",
  sep = " "
)
```

Use this whenever the task YAML carries a `depends_on:` key referencing an
upstream nlmixr2lib model.
