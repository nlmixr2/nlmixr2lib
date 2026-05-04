Zhu_2017_lebrikizumab <- function() {
  description <- "Lebrikizumab population PK model (Zhu 2017): two-compartment model with first-order absorption after SC dosing in adults with moderate-to-severe asthma."
  reference <- "Zhu R, Zheng Y, Dirks NL, et al. Model-based clinical pharmacology profiling and exposure-response relationships of the efficacy and biomarker of lebrikizumab in patients with moderate-to-severe asthma. Pulmonary Pharmacology & Therapeutics. 2017;46:88-98. doi:10.1016/j.pupt.2017.08.010"
  vignette <- "Zhu_2017_lebrikizumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on CL, Vc, Vp, and Q; normalized as WT/70 per Table 3 footnote.",
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on CL; normalized as AGE/40 per Table 3 footnote.",
      source_name        = "AGE"
    ),
    SEXF = list(
      description        = "Biological sex indicator, 1 = female, 0 = male",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Source paper reports the covariate in the canonical SEXF encoding (1 = female).",
      source_name        = "SEXF"
    ),
    ADA_POS = list(
      description        = "Anti-drug antibody positivity",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ADA-negative)",
      notes              = "Source used column 'ADA' with the semantic of 'ever positive'. Renamed to the canonical ADA_POS per covariate-columns.md; confirm that the time-frame semantics of ADA_POS (current/observation-time) match the source 'ever-positive' definition before applying the effect.",
      source_name        = "ADA"
    ),
    RACE_BLACK = list(
      description        = "Black or African American race indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (White / reference race group)",
      notes              = "Power-style multiplicative effect on CL.",
      source_name        = "RACE_BLACK"
    ),
    RACE_ASIAN = list(
      description        = "Asian race indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (White / reference race group)",
      notes              = "Power-style multiplicative effect on CL.",
      source_name        = "RACE_ASIAN"
    ),
    RACE_OTHER = list(
      description        = "Race category 'Other' indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (White / reference race group)",
      notes              = "Power-style multiplicative effect on CL.",
      source_name        = "RACE_OTHER"
    ),
    FORM_NS0 = list(
      description        = "NS0 cell-line formulation indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (reference CHO formulation)",
      notes              = "Typically 0 in routine use; affects ka and bioavailability.",
      source_name        = "FORM_NS0"
    ),
    FORM_CHO_PHASE2 = list(
      description        = "CHO Phase 2 formulation indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (reference formulation)",
      notes              = "Typically 0 in routine use; affects ka and bioavailability.",
      source_name        = "FORM_CHO_PHASE2"
    )
  )

  population <- list(
    n_subjects     = 2148L,
    n_observations = 21917L,
    n_studies      = 6L,
    age_range      = "18-75 years (adults only)",
    age_median     = "48 years",
    weight_range   = "40-165 kg",
    weight_median  = "~77 kg",
    sex_female_pct = 58,
    race_ethnicity = "White majority; Black, Asian, and 'Other' categories also represented (each with a separate CL effect in the final model).",
    disease_state  = "Pooled analysis across 6 studies: 2 Phase I studies in healthy volunteers (n=114), 1 Phase II study in asthma, 1 Phase II study in atopic dermatitis, 1 Phase II study in idiopathic pulmonary fibrosis, and the Phase III MILLY program in moderate-to-severe asthma.",
    dose_range     = "37.5-250 mg SC (pooled analysis also included some IV data from the Phase I studies).",
    regions        = "Multi-regional (not reported in detail in Zhu 2017).",
    notes          = "Reference covariate values are WT = 70 kg and AGE = 40 years (Table 3 footnote). Three formulations were evaluated: the reference CHO formulation used in late development, an early-development NS0 formulation, and an interim CHO formulation used in Phase 2 ('CHO Phase 2'); indicator covariates FORM_NS0 and FORM_CHO_PHASE2 are both 0 for the reference formulation."
  )

  ini({
    lcl <- log(0.156); label("Clearance (L/day)")
    lvc <- log(4.10); label("Central volume of distribution (L)")
    lvp <- log(1.45); label("Peripheral volume of distribution (L)")
    lq <- log(0.284); label("Intercompartmental clearance (L/day)")
    lka <- log(0.239); label("Absorption rate (1/day)")
    lfdepot <- log(0.856); label("Subcutaneous bioavailability (fraction)")

    # Zhu 2017 Table 3 reports WT effect on CL as 1.00; ambiguous whether this was
    # fixed (theta locked at 1.00) or estimated to ~1.00. Kept as estimated; flag
    # for follow-up if the intended behavior is fixed allometry.
    e_wt_cl <- 1.00; label("Effect of body weight on clearance (unitless)")
    e_wt_vc <- 0.814; label("Effect of body weight on central volume (unitless)")
    e_wt_vp <- 0.692; label("Effect of body weight on peripheral volume (unitless)")
    e_wt_q <- 0.479; label("Effect of body weight on intercompartmental clearance (unitless)")
    e_age_cl <- 0.0241; label("Effect of age on clearance (unitless)")
    e_sexf_cl <- 1.06; label("Effect of sex on clearance (unitless)")
    e_race_black_cl <- 1.07; label("Effect of race (black or African American) on clearance (unitless)")
    e_race_asian_cl <- 1.09; label("Effect of race (Asian) on clearance (unitless)")
    e_race_other_cl <- 1.11; label("Effect of race (other) on clearance (unitless)")
    e_form_ns0_ka <- 0.981; label("Effect of NS0 formulation on absorption rate (unitless)")
    e_form_cho_phase2_ka <- 0.989; label("Effect of CHO formulation used during Phase 2 on absorption rate (unitless)")
    e_form_ns0_fdepot <- 1.00; label("Effect of NS0 formulation on bioavailability (unitless)")
    e_form_cho_phase2_fdepot <- 0.973; label("Effect of CHO formulation used during Phase 2 on bioavailability (unitless)")
    e_ada_pos_cl <- 1.04; label("Effect of anti-drug antibody (ADA) positivity on clearance (unitless)")

    # IIV variance-covariance matrix (omega^2 / cov) from Zhu 2017 Table 3.
    # Lower-triangular order is var(CL); cov(CL,Vc), var(Vc); cov(CL,ka), cov(Vc,ka), var(ka).
    # A prior version of this file stored sqrt(variance) in these slots; this is the fix.
    etalcl + etalvc + etalka ~
      c(
        0.105,
        0.0832, 0.124,
        0.00203, 0.00439, 0.154
      )

    CcpropSd <- 0.0490; label("Proportional residual error (fraction)")
    CcaddSd <- 0.00154; label("Additive residual error (ug/mL)")
  })
  model({
    # Normalized continuous covariate values based on footnote to Table 3
    WTNORM <- WT/70
    AGENORM <- AGE/40

    cl <-
      exp(lcl + etalcl) *
      WTNORM^e_wt_cl * AGENORM^e_age_cl * e_sexf_cl^SEXF *
      e_race_black_cl^RACE_BLACK * e_race_asian_cl^RACE_ASIAN * e_race_other_cl^RACE_OTHER *
      e_ada_pos_cl^ADA_POS
    vc <- exp(lvc + etalvc) * WTNORM^e_wt_vc
    vp <- exp(lvp) * WTNORM^e_wt_vp
    q <- exp(lq) * WTNORM^e_wt_q
    ka <- exp(lka + etalka) * e_form_ns0_ka^FORM_NS0 * e_form_cho_phase2_ka^FORM_CHO_PHASE2
    fdepot <- exp(lfdepot) * e_form_ns0_fdepot^FORM_NS0 * e_form_cho_phase2_fdepot^FORM_CHO_PHASE2
    Cc <- linCmt()
    f(depot) <- fdepot
    Cc ~ add(CcaddSd) + prop(CcpropSd)
  })
}
