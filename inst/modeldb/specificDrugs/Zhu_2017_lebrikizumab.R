Zhu_2017_lebrikizumab <- function() {
  description <- "Lebrikizumab population PK model (Zhu 2017): two-compartment model with first-order absorption after SC dosing in adults with moderate-to-severe asthma."
  reference <- "Zhu R, Zheng Y, Dirks NL, et al. Model-based clinical pharmacology profiling and exposure-response relationships of the efficacy and biomarker of lebrikizumab in patients with moderate-to-severe asthma. Pulmonary Pharmacology & Therapeutics. 2017;46:88-98. doi:10.1016/j.pupt.2017.08.010"
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
    n_subjects     = "TODO: from source paper",
    n_studies      = "TODO: from source paper",
    age_range      = "TODO: from source paper",
    age_median     = "TODO: from source paper",
    weight_range   = "TODO: from source paper",
    weight_median  = "TODO: from source paper",
    sex_female_pct = "TODO: from source paper",
    race_ethnicity = "TODO: from source paper",
    disease_state  = "Adults with moderate-to-severe asthma (primary indication in Zhu 2017); model also pooled data from atopic dermatitis program formulations.",
    dose_range     = "TODO: from source paper",
    regions        = "TODO: from source paper",
    notes          = "Reference covariate values are WT = 70 kg and AGE = 40 years (Table 3 footnote)."
  )

  ini({
    lcl <- log(0.156); label("Clearance (L/day)")
    lvc <- log(4.10); label("Central volume of distribution (L)")
    lvp <- log(1.45); label("Peripheral volume of distribution (L)")
    lq <- log(0.284); label("Intercompartmental clearance (L/day)")
    lka <- log(0.239); label("Absorption rate (1/day)")
    lfdepot <- log(0.856); label("Subcutaneous bioavailability (fraction)")

    e_cl_wt <- 1.00; label("Effect of body weight on clearance (unitless)")
    e_vc_wt <- 0.814; label("Effect of body weight on central volume (unitless)")
    e_vp_wt <- 0.692; label("Effect of body weight on peripheral volume (unitless)")
    e_q_wt <- 0.479; label("Effect of body weight on intercompartmentl clearance (unitless)")
    e_cl_age <- 0.0241; label("Effect of age on clearance (unitless)")
    e_cl_sexf <- 1.06; label("Effect of sex on clearance (unitless)")
    e_cl_race_black <- 1.07; label("Effect of race (black or African American) on clearance (unitless)")
    e_cl_race_asian <- 1.09; label("Effect of race (Asian) on clearance (unitless)")
    e_cl_race_other <- 1.11; label("Effect of race (other) on clearance (unitless)")
    e_ka_form_nso <- 0.981; label("Effect of NSO formulation on absorption rate (unitless)")
    e_ka_form_cho_phase2 <- 0.989; label("Effect of CHO formulation used during Phase 2 on absorption rate (unitless)")
    e_f_form_nso <- 1.00; label("Effect of NSO formulation on bioavailability (unitless)")
    e_f_form_cho_phase2 <- 0.973; label("Effect of CHO formulation used during Phase 2 on bioavailability (unitless)")
    e_cl_ada_positive<- 1.04; label("Effect of anti-drug antibody (ADA) positivity on clearance (unitless)")

    # converted from covariance matrix reported in Table 3
    etalcl + etalvc + etalka ~
      c(
        0.32403703,
        0.28844410, 0.35213634,
        0.04505552, 0.06625708, 0.39242834
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
      WTNORM^e_cl_wt * AGENORM^e_cl_age * e_cl_sexf^SEXF *
      e_cl_race_black^RACE_BLACK * e_cl_race_asian^RACE_ASIAN * e_cl_race_other^RACE_OTHER *
      e_cl_ada_positive^ADA_POS
    vc <- exp(lvc + etalvc) * WTNORM^e_vc_wt
    vp <- exp(lvp) * WTNORM^e_vp_wt
    q <- exp(lq) * WTNORM^e_q_wt
    ka <- exp(lka + etalka) * e_ka_form_nso^FORM_NS0 * e_ka_form_cho_phase2^FORM_CHO_PHASE2
    fdepot <- exp(lfdepot) * e_f_form_nso^FORM_NS0 * e_f_form_cho_phase2^FORM_CHO_PHASE2
    Cc <- linCmt()
    f(depot) <- fdepot
    Cc ~ add(CcaddSd) + prop(CcpropSd)
  })
}
