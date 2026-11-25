Baverel_2019_tremelimumab <- function() {
  description <- "Two-compartment linear population PK model for tremelimumab (anti-CTLA-4 IgG2) in adults with unresectable pleural or peritoneal malignant mesothelioma from the DETERMINE trial (Baverel 2019); final model without covariates, used to derive steady-state exposure for an overall-survival exposure-response analysis"
  reference <- "Baverel P, Roskos L, Tatipalli M, Lee N, Stockman P, Taboada M, Vicini P, Horgan K, Narwal R. Exposure-Response Analysis of Overall Survival for Tremelimumab in Unresectable Malignant Mesothelioma: The Confounding Effect of Disease Status. Clin Transl Sci. 2019;12(5):450-458. doi:10.1111/cts.12633"
  vignette <- "Baverel_2019_tremelimumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  compartmentData <- list(
    central = list(analyte = "tremelimumab", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "tremelimumab", units = "mg", specimen = "plasma", verified = TRUE)
  )

  # No covariates were retained in the final model (Supplementary Figure S1
  # caption: 'none were kept in the final PK model of tremelimumab'). The
  # screened covariates and the effect sizes printed in Table 2 are recorded
  # below for provenance only.
  covariateData <- list()

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Female sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Table 2 intermediate model 1: 'CL % reduction for female patients vs. male' = -0.19 (dOFV = -14, P = 0.00018). Statistically significant but changed AUCss by less than the 30% clinical-relevance threshold (Supplementary Figure S1), so not retained in the final model.",
      source_name = "SEX"
    ),
    CRP = list(
      description = "Baseline C-reactive protein (standard assay)",
      units = "mg/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Table 2 intermediate model 3: proportional linear effect on CL from the 33 mg/L median, coefficient 0.0034 per mg/L (dOFV = -44, P < 1e-10). Not retained (effect on AUCss < 30%, Supplementary Figure S1).",
      source_name = "CRP"
    ),
    TUMSZ = list(
      description = "Baseline tumor size (target-lesion size)",
      units = "mm",
      type = "continuous",
      reference_category = NULL,
      notes = "Table 2 'Final model' row 7: proportional linear effect on CL from the 97 mm median, coefficient 0.00058 per mm (dOFV = -4.3, P = 0.039). Not retained in the final PK model used for exposure derivation (effect on AUCss < 30%, Supplementary Figure S1).",
      source_name = "TUMSIZE"
    ),
    ALB = list(
      description = "Baseline serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Table 2 intermediate model 2: proportional linear effect on CL from the 31 g/L median, coefficient 0.00094 (dOFV = -1, P = 0.32); not significant.",
      source_name = "ALB"
    ),
    ECOG_GE1 = list(
      description = "ECOG performance status >= 1 indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (ECOG 0)",
      notes = "Table 2 intermediate model 4: ECOG effect on CL = 0.019 (dOFV = -1.4, P = 0.24); not significant.",
      source_name = "ECOG"
    ),
    LDH = list(
      description = "Baseline serum lactate dehydrogenase",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Table 2 intermediate model 6: LDH effect on CL = -2.1e-04 (dOFV = -1.3, P = 0.26); not significant. The centring value for LDH is not printed.",
      source_name = "LDH"
    ),
    ADA_POS = list(
      description = "Post-baseline anti-drug-antibody positive indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (ADA negative)",
      notes = "Results, PK modeling: 15 of 377 ADA-evaluable tremelimumab-arm patients (4%) were ADA-positive post-baseline; 'ADA had no apparent effect on PK'.",
      source_name = "ADA"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 376L,
    n_studies = 1L,
    n_observations = 1328L,
    disease_state = "Unresectable pleural or peritoneal malignant mesothelioma, second- or third-line after one or two prior systemic treatments including a platinum-based regimen",
    dose_range = "Tremelimumab 10 mg/kg IV every 4 weeks for seven doses, then every 12 weeks",
    regions = "Multinational (DETERMINE, NCT01843374)",
    crp_median = "33 mg/L",
    albumin_median = "31 g/L",
    tumor_size_median = "97 mm",
    ada_positive_pct = 4,
    notes = "Phase IIb randomized, double-blind, placebo-controlled DETERMINE trial (tremelimumab n = 382, placebo n = 189; 2:1). The PK data set had 376 patients with 1,328 evaluable concentrations (six patients had no evaluable PK). Baseline demographics (age, weight, sex, race) are not tabulated in Baverel 2019; they are reported in the primary trial publication (Maio et al., Lancet Oncol 2017). Covariate medians come from the Table 2 footnote. Further screened covariates not representable as register columns: EORTC prognostic status (Table 2 intermediate model 5, effect 0.068 on CL, P = 0.16) and epithelioid histology (intermediate model 7, 8.4% increase on CL, P = 0.19); neither retained."
  )

  ini({
    # Final population PK estimates, Baverel 2019 Table 1. Table 1 header
    # 'BSV (Omega) (%)' with footnote 'random effect ... variance Omega^2'
    # makes the printed percent 100 * Omega (the SD of eta), so the variance
    # is (BSV/100)^2.
    lcl <- log(0.310); label("Clearance CL (L/day)") # Table 1: CL = 0.310 L/day
    lvc <- log(3.85); label("Central volume of distribution V1 (L)") # Table 1: V1 = 3.85 L
    lvp <- log(1.72); label("Peripheral volume of distribution V2 (L)") # Table 1: V2 = 1.72 L
    lq <- log(0.273); label("Intercompartmental clearance Q (L/day)") # Table 1: Q = 0.273 L/day

    etalcl ~ 0.380^2 # Table 1: BSV CL = 38.0%
    etalvc ~ 0.325^2 # Table 1: BSV V1 = 32.5%
    etalvp ~ 0.252^2 # Table 1: BSV V2 = 25.2%

    propSd <- 0.377; label("Proportional residual error (fraction)") # Table 1: CV% (proportional epsilon) = 37.7%
  })
  model({
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp + etalvp)
    q <- exp(lq)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Dose in mg and volumes in L give mg/L = ug/mL.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
