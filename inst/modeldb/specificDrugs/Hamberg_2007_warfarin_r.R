Hamberg_2007_warfarin_r <- function() {
  description <- "R-warfarin population PK (1-compartment, first-order absorption) with age as the only structural covariate on CL_R (Hamberg 2007). R-warfarin was not found to contribute (additive or competitive) to the INR PD; the companion file Hamberg_2007_warfarin_s carries the S-warfarin PK and the INR PD model."
  reference <- paste(
    "Hamberg A-K, Dahl M-L, Barban M, Scordo MG, Wadelius M, Pengo V, Padrini R, Jonsson EN.",
    "A PK-PD Model for Predicting the Impact of Age, CYP2C9, and VKORC1 Genotype",
    "on Individualization of Warfarin Therapy.",
    "Clin Pharmacol Ther. 2007;81(4):529-538. doi:10.1038/sj.clpt.6100084.",
    "PMID: 17301738.",
    "R-warfarin PK parameters and age effect from Table 3; structural equations from the Appendix."
  )
  vignette <- "Hamberg_2007_warfarin_pkpd_pgx"

  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    AGE = list(
      description        = "Subject age at baseline",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed at baseline. Linear scale on CL_R with reference age 71 years (cohort median); Hamberg 2007 Appendix Eq 7. CL_R decreases approximately 0.98% per year of age above 71. The only structural covariate retained on CL_R; sex and weight were tested and not retained.",
      source_name        = "AGE"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 150L,
    n_studies      = 2L,
    age_range      = "22-87 years",
    age_median     = "71 years",
    weight_range   = "45-120 kg",
    weight_median  = "80 kg",
    sex_female_pct = 34,
    race_ethnicity = "Italian (not stratified by race in source)",
    disease_state  = "Adults on long-term warfarin anticoagulant therapy for thromboembolic prophylaxis",
    dose_range     = "Single 10 mg racemic dose (Study I) and 6.25-78.75 mg/week maintenance (median 29.375 mg/week)",
    regions        = "Italy (Thrombosis Center, University of Padova)",
    n_pk_records   = "171 R-warfarin observations after single-dose + 150 after stable maintenance dosing (Study I + Study II)",
    notes          = paste(
      "Pooled across Italian Studies I (n=57) and II (n=93); same patient cohort as the companion S-warfarin model.",
      "See Hamberg 2007 Table 1 for demographics; Table 3 for R-warfarin parameter estimates."
    )
  )

  ini({
    # ============================================================
    # R-warfarin 1-cmt PK -- Hamberg 2007 Table 3 + Appendix
    # ============================================================
    lcl       <- log(0.139)        ; label("Apparent oral CL_R typical (L/h), age 71 y")    # Table 3 (estimated, RSE 3.13%)
    lvc       <- log(12.8)         ; label("Apparent volume of distribution V (L)")         # Table 3 (estimated, RSE 4.08%)
    lka       <- fixed(log(2))     ; label("R-warfarin absorption rate Ka_R (1/h)")          # Table 3 (fixed, sensitivity confirmed across 0.5-5; Methods)

    # Age effect: CL_R decreases ~0.98% per year of age above 71; Table 3 reports magnitude only,
    # prose ("CLR was reduced with increasing age, decreasing by approximately 10% per decade")
    # establishes the negative sign.
    e_age_cl  <- -0.00980          ; label("Age effect on CL_R (fractional change per year, ref 71 y)")  # Table 3 magnitude 0.98 (% change/year), RSE 25.4% (sign from prose)

    # ============================================================
    # IIV block: omega(CL_R) 25.5%, omega(V) 26.8%, with covariance term between CL_R and V.
    #   omega values from Hamberg 2007 Table 3 (CL_R RSE 17.6%, V RSE 39.8%); variances on log scale
    #   via log(1 + CV^2). The off-diagonal covariance is not numerically reported in Table 3; the
    #   starting value below preserves the structural form and will be re-estimated by the user.
    # ============================================================
    etalcl + etalvc ~ c(0.0630, 0.020, 0.0696)

    # ============================================================
    # Residual error: Hamberg 2007 Table 3 reports two log-additive (~ constant-CV) residuals --
    # single-dose (0.0865, RSE 17.3%) and steady-state (0.289, RSE 8.10%). Encoded as the larger
    # (steady-state) value; see vignette Errata.
    # ============================================================
    propSd <- 0.289                ; label("R-warfarin proportional residual SD (fraction; Hamberg 2007 steady-state value 0.289; single-dose was 0.0865)")  # Table 3 (estimated, RSE 8.10%)
  })

  model({
    # ---- Individual structural parameters ----
    cl <- exp(lcl + etalcl) * (1 + e_age_cl * (AGE - 71))
    vc <- exp(lvc + etalvc)
    ka <- exp(lka)

    # ---- Micro-constants ----
    kel <- cl / vc

    # ---- R-warfarin 1-compartment PK ----
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # ---- Observation ----
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
