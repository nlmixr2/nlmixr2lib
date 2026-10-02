Li_2020_bevacizumab <- function() {
  description <- paste(
    "Two-compartment population PK model with zero-order IV infusion and first-order elimination",
    "for bevacizumab (PF-06439535 biosimilar and EU-sourced reference Avastin, pooled) in adults with",
    "advanced non-squamous non-small cell lung cancer (Li 2020). Baseline body weight enters CL and",
    "V1 as power terms normalised to 71 kg; male sex increases CL and V1 as fractional shifts; the",
    "drug-product (PF-06439535 vs bevacizumab-EU) multiplier on CL and V1 was retained by the authors",
    "for the similarity assessment despite not being statistically significant."
  )
  reference <- paste(
    "Li CSW, Sweeney K, Cronenberger C.",
    "Population pharmacokinetic modeling of PF-06439535 (a bevacizumab biosimilar) and reference",
    "bevacizumab (Avastin) in patients with advanced non-squamous non-small cell lung cancer.",
    "Cancer Chemother Pharmacol. 2020;85(3):487-499. doi:10.1007/s00280-019-03946-8"
  )
  vignette <- "Li_2020_bevacizumab"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters CL and V1 as (WT/71)^theta (Li 2020 final-model equations, Results 'Final PK model').",
        "71 kg is the median baseline weight of the PK population (Table 1)."
      ),
      source_name = "BWT"
    ),
    SEXF = list(
      description = "Biological sex, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (female; the typical patient in Li 2020 is a 71-kg female)",
      notes = paste(
        "Li 2020 codes sex as Male = 1, Female = 2 (Figure 2 axis labels) with female as the reference",
        "level of the fractional-change Eq. 4. The model derives male = 1 - SEXF and applies",
        "(1 + theta * male) to CL and V1."
      ),
      source_name = "Sex (Male=1, Female=2)"
    ),
    TRT_PF06439535 = list(
      description = "PF-06439535 bevacizumab-biosimilar treatment-arm indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (reference bevacizumab sourced from the European Union, Avastin)",
      notes = paste(
        "1 = PF-06439535 (bevacizumab biosimilar), 0 = bevacizumab-EU. Li 2020 Eq. 5 parameterises the",
        "drug-product effect as TVP = theta^COV; the final-model equations print '* (1.02 in",
        "PF-06439535)' on CL and '* (1.07 in PF-06439535)' on V1. The effect was not statistically",
        "significant in the stepwise covariate search (Table 2 footnote b) and its 95% bootstrap CIs",
        "include unity; it was retained to quantify drug-product similarity."
      ),
      source_name = "DP (PF-06439535=1, Bevacizumab-EU=2)"
    )
  )

  covariatesDataExcluded <- list(
    ALB = list(
      description = "Baseline serum albumin",
      units = "g/dL",
      type = "continuous",
      notes = "Tested in the stepwise covariate search on CL and V1 but not retained (Results; Online Resource Table S2)."
    ),
    ALT = list(
      description = "Baseline alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Tested in the stepwise covariate search on CL and V1 but not retained (Results; Online Resource Table S2)."
    ),
    ALP = list(
      description = "Baseline alkaline phosphatase",
      units = "U/L",
      type = "continuous",
      notes = "Tested in the stepwise covariate search on CL and V1 but not retained (Results; Online Resource Table S2)."
    )
  )

  compartmentData <- list(
    central = list(
      analyte = "bevacizumab",
      units = "mg",
      specimen = "serum",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "bevacizumab",
      units = "mg",
      specimen = "serum",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 705L,
    n_studies = 1L,
    weight_range = "28.0-135 kg",
    weight_median = "71.0 kg",
    sex_female_pct = 35.2,
    race_ethnicity = c(White = 88.7, Asian = 10.6, Black = 0.567, Other = 0.142),
    disease_state = "Advanced non-squamous non-small cell lung cancer, first-line treatment",
    dose_range = paste(
      "15 mg/kg IV every 21 days (90-min first infusion, 60-min second, 30-min thereafter) with",
      "paclitaxel + carboplatin for 4-6 cycles, then bevacizumab monotherapy for up to ~1 year"
    ),
    regions = "Multinational (Study B7391003, NCT02364999)",
    notes = paste(
      "Demographics from Li 2020 Table 1 (PK population N = 705: 351 PF-06439535, 354 bevacizumab-EU).",
      "8632 serum concentrations analysed (LLOQ 250 ng/mL; post-dose BLQ excluded by the M1 method;",
      "27 influential |CWRES| or |IWRES| > 6 outliers omitted). Sparse sampling: pre-dose troughs",
      "every cycle and 1-h post-infusion peaks on Cycle 1 Day 1 and Cycle 5 Day 1. 64.8% male;",
      "ECOG 0/1 = 28.8/71.2%. NONMEM 7.2, FOCE-I."
    )
  )

  ini({
    lcl <- log(0.0113)
    label("Clearance for a 71-kg female on bevacizumab-EU (L/h)") # Table 2 'CL (L/h)' = 0.0113
    lvc <- log(2.99)
    label("Central volume for a 71-kg female on bevacizumab-EU (L)") # Table 2 'V1 (L)' = 2.99
    lq <- log(0.269)
    label("Intercompartmental clearance (L/h)") # Table 2 'Q (L/h)' = 0.269
    lvp <- log(6.09)
    label("Peripheral volume (L)") # Table 2 'V2 (L)' = 6.09

    e_wt_cl <- 0.354
    label("Power exponent of (WT/71) on CL (unitless)") # Table 2 'BWT effect on CL' = 0.354
    e_wt_vc <- 0.468
    label("Power exponent of (WT/71) on V1 (unitless)") # Table 2 'BWT effect on V1' = 0.468
    e_male_cl <- 0.262
    label("Fractional increase in CL for males (unitless)") # Table 2 'Sex effect on CL' = 0.262; final-model equation '1.262 in male'
    e_male_vc <- 0.247
    label("Fractional increase in V1 for males (unitless)") # Table 2 'Sex effect on V1' = 0.247; final-model equation '1.247 in male'
    e_pf06439535_cl <- 1.02
    label("Multiplier on CL for PF-06439535 vs bevacizumab-EU (unitless)") # Table 2 'Drug product effect on CL' = 1.02 (Eq. 5, theta^COV)
    e_pf06439535_vc <- 1.07
    label("Multiplier on V1 for PF-06439535 vs bevacizumab-EU (unitless)") # Table 2 'Drug product effect on V1' = 1.07 (Eq. 5, theta^COV)

    etalcl ~ 0.0871 # Table 2 'omega2 CL' = 0.0871 (Discussion: 29.5% CV)
    etalvc ~ 0.117 # Table 2 'omega2 V1' = 0.117 (Discussion: 34.2% CV)

    expSd <- 0.284
    label("Additive residual error on log-transformed concentrations (SD, log scale)") # Table 2 'Residual additive error' = 0.284; Eq. 2 ln(Y) = ln(F) + W*eps with var(eps) = 1
  })
  model({
    male <- 1 - SEXF

    cl <- exp(lcl + etalcl) *
      (WT / 71)^e_wt_cl *
      (1 + e_male_cl * male) *
      e_pf06439535_cl^TRT_PF06439535
    vc <- exp(lvc + etalvc) *
      (WT / 71)^e_wt_vc *
      (1 + e_male_vc * male) *
      e_pf06439535_vc^TRT_PF06439535
    q <- exp(lq)
    vp <- exp(lvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d / dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d / dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Dose in mg, volume in L -> Cc in mg/L (= ug/mL)
    Cc <- central / vc
    Cc ~ lnorm(expSd)
  })
}
