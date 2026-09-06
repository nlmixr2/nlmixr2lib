Hughes_2022_brepocitinib <- function() {
  description <- "One-compartment population pharmacokinetic model with first-order absorption and a tablet-only absorption lag for oral brepocitinib (PF-06700841, a dual TYK2/JAK1 inhibitor) in healthy participants and patients with plaque psoriasis or alopecia areata. Apparent clearance and apparent central volume carry allometric body-weight scaling with exponents fixed to literature values at a 70 kg reference, and apparent clearance is 24.3 percent lower in Asian participants. A high-fat meal slows absorption by 69.9 percent and lowers relative bioavailability by 28.3 percent, while doses of 175 mg and above raise relative bioavailability by 35.1 percent. Interindividual variability on CL/F and Vc/F is correlated, and the log-scale residual error differs between the phase 1 healthy-volunteer studies and the phase 2 patient studies."
  reference <- paste(
    "Hughes JH, Qiu R, Banfield C, Dowty ME, Nicholas T.",
    "Population pharmacokinetics of oral brepocitinib in healthy volunteers",
    "and patients.",
    "Clinical Pharmacology in Drug Development. 2022;11(12):1447-1456.",
    "doi:10.1002/cpdd.1163",
    sep = " "
  )
  vignette <- "Hughes_2022_brepocitinib"
  paper_specific_residual_sds <- c("expSdHv", "expSdPat")
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed baseline body weight entered as an allometric power",
        "function referenced to 70 kg, with exponents fixed to literature",
        "values of 0.75 on CL/Frel and 1.0 on Vc/Frel (Hughes 2022 Table 3",
        "rows 'Effect of weight on CL/F (70 kg reference)' and 'Effect of",
        "weight on Vc/F (70 kg reference)', both marked fixed; Equations 9",
        "and 10). The Full Covariate Model section states the exponents were",
        "fixed because allometric scaling on CL/F could not be estimated with",
        "adequate precision during full-model estimation. The 70 kg reference",
        "is a rounded literature value, not the observed median: Table 2",
        "reports a median of 83 kg (mean 86.6, SD 22.0, range 45-204 kg).",
        "The Discussion notes the Asian subgroup averaged 78.5 kg versus",
        "87.4 kg in the remainder, so part of the observed Asian exposure",
        "difference is already absorbed by this allometric term."
      ),
      source_name        = "BWT"
    ),
    RACE_ASIAN = list(
      description        = "Asian race indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Asian; the pooled White, African American and Other race groups)",
      notes              = paste(
        "1 = Asian, 0 = otherwise. Enters CL/Frel as a fractional shift,",
        "CL/Frel = 18.7 * (BWT/70)^0.75 * (1 - 0.243 * Asian) (Hughes 2022",
        "Equation 9; Table 3 row 'Effect of Asian subjects on CL/F' =",
        "-0.243, 95% CI -0.363 to -0.123), i.e. 24.3% lower apparent",
        "clearance. 30 of 379 subjects (8%) were Asian (Table 2). The",
        "Discussion flags this as resting on a limited number of Asian",
        "subjects with dense sampling and states that the mechanism is",
        "unclear given that CYP3A4/5 and CYP1A2 polymorphisms are not",
        "primary drivers of brepocitinib PK variability. The reference",
        "population for CL/Frel stated in the Figure 1 caption is a 70 kg",
        "White male patient with psoriasis aged 42 years taking 30 mg once",
        "daily."
      ),
      source_name        = "Asian"
    ),
    FED_HIGHFAT = list(
      description        = "High-fat-meal-at-dosing indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (dosed without regard to food, the modal condition in the pooled phase 2 studies)",
      notes              = paste(
        "1 = the dose was administered with a high-fat meal (800-1000",
        "calories, 500-600 of them from fat, per the NCT02310750 food-effect",
        "arm described in Methods / Study Data), 0 = any other meal",
        "condition. Acts on two parameters: ka = 3.46 * (1 - 0.699 * Fed)",
        "(Hughes 2022 Equation 7) and Frel = 1 * (1 - 0.283 * Fed) *",
        "(1 + 0.351 * Dose>=175mg) (Equation 8). The Discussion notes that",
        "the model-based 28.3% reduction in relative bioavailability is",
        "somewhat larger than the ~18% AUC reduction reported by",
        "noncompartmental analysis, but that the difference is not",
        "significant given the covariate 95% CI of 16.9%-39.1%. The paper",
        "reports that no covariate effect could be estimated to separate",
        "explicitly fasted subjects from subjects dosed without regard for",
        "food, so the reference level pools those two conditions."
      ),
      source_name        = "Fed"
    ),
    FORM_TABLET = list(
      description        = "Tablet-versus-oral-suspension formulation indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (oral suspension, used in the first-in-human study NCT02310750)",
      notes              = paste(
        "1 = oral tablet, 0 = oral suspension. Hughes 2022 Equation 6 is",
        "written on the complementary suspension indicator, Alag = 0.24 *",
        "(1 - 1 * Suspension), with the coefficient fixed at -1.00 (Table 3",
        "row 'Effect of suspension formulation on Alag'), so the absorption",
        "lag collapses to exactly zero for the suspension and equals the",
        "estimated 0.240 h for the tablet. This file stores the canonical",
        "tablet-oriented column and derives Suspension = 1 - FORM_TABLET",
        "inside model() so the printed equation is reproduced verbatim;",
        "the same column and orientation are used by the companion",
        "Maleki_2024_brepocitinib.R model of the same programme. 313 of 379",
        "subjects (83%) received tablets and 66 (17%) received the",
        "suspension (Table 2). Formulation was also statistically",
        "significant on ka during full-model estimation but was not carried",
        "into the final model because the resulting ka estimates were not",
        "meaningfully different from the population value."
      ),
      source_name        = "Suspension"
    ),
    DOSE_HIGH = list(
      description        = "Indicator for administered brepocitinib doses of 175 mg and above",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (administered dose below 175 mg; the 100 mg cohort sits in the reference group)",
      notes              = paste(
        "1 = administered brepocitinib dose of 175 mg or more, which in this",
        "analysis is the 175 mg once-daily multiple-dose arm and the 200 mg",
        "single-dose arms of NCT02310750 and NCT03656952; 0 otherwise. Gates",
        "a 35.1% uplift in relative bioavailability, Frel = 1 * (1 - 0.283 *",
        "Fed) * (1 + 0.351 * Dose>=175mg) (Hughes 2022 Equation 8; Table 3",
        "row 'Effect of dose on Frel (dose >= 175 mg)' = 0.351, 95% CI 0.194",
        "to 0.508). The threshold is a step function only because no doses",
        "between 100 and 175 mg were studied: the Results state that a",
        "continuous dose-Frel relationship was explored but gave minimal",
        "improvement in objective function value and predictive performance,",
        "and the Discussion argues the step estimate is likely the maximum",
        "effect with the true relationship rising gradually between 100 and",
        "175 mg. NOTE the threshold differs from the companion",
        "Maleki_2024_brepocitinib.R model, which sets DOSE_HIGH = 1 above",
        "100 mg/day; both papers flag the same 175 and 200 mg arms, so the",
        "two encodings select the same subjects in practice."
      ),
      source_name        = "Dose>=175mg"
    ),
    DIS_HEALTHY = list(
      description        = "Healthy-participant cohort indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (patient; the pooled plaque psoriasis and alopecia areata cohorts)",
      notes              = paste(
        "1 = healthy participant enrolled in one of the three phase 1",
        "studies (NCT02310750, NCT03236493, NCT03656952), 0 = patient",
        "enrolled in one of the two phase 2 studies (NCT02969018 plaque",
        "psoriasis, NCT02974868 alopecia areata). Used ONLY to select the",
        "residual-error magnitude: Hughes 2022 Base Model states 'a",
        "proportional error model was used, with separate error parameters",
        "estimated for phase 1 and phase 2 studies', and Table 3 reports",
        "0.527 for phase 1 and 0.875 for phase 2. The Abstract restates the",
        "same two numbers as applying to 'healthy volunteers' and",
        "'patients' respectively, which is what licenses mapping the",
        "phase 1 / phase 2 study split onto this health-status column: every",
        "phase 1 subject in the analysis was a healthy participant and every",
        "phase 2 subject was a patient. 92 of 379 subjects (24%) were",
        "healthy, 210 (55%) had plaque psoriasis and 77 (20%) had alopecia",
        "areata (Table 2). Patient status was ALSO screened as a covariate",
        "on CL/F and on Frel during full-model estimation and was retained",
        "in neither: the Discussion reports 'no statistically or clinically",
        "significant difference between healthy participants and patient",
        "populations with plaque psoriasis or alopecia areata' in apparent",
        "clearance."
      ),
      source_name        = "Patient type"
    )
  )

  # Covariates that Hughes 2022 Table 1 lists as evaluated during model
  # development but that the final model does not retain. Documentation only;
  # not referenced in model().
  covariatesDataExcluded <- list(
    AGE = list(
      description        = "Age at baseline",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened on CL/F and on Vc/F (Hughes 2022 Table 1) and not retained",
        "in the full or final model. Table 2 reports mean 41.9 years (SD",
        "12.9), median 42, range 18-75. The Figure 1 reference populations",
        "are quoted at 42 years."
      ),
      source_name        = "Age"
    ),
    SEXF = list(
      description        = "Female sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = paste(
        "Screened on CL/F and on Vc/F (Hughes 2022 Table 1) and not retained",
        "in the full or final model. Table 2 reports 263 male (69%) and 116",
        "female (31%) subjects. The Figure 1 reference populations are male."
      ),
      source_name        = "Sex"
    ),
    CRCL = list(
      description        = "Creatinine clearance at baseline",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Reported in the Hughes 2022 Table 2 demographic summary (mean",
        "132.7 mL/min, SD 38.7; median 124, range 58-279) but not listed",
        "among the covariates evaluated in Table 1 and not present in the",
        "final model. Renal elimination of unchanged brepocitinib is low",
        "(urinary recovery <16%, Introduction), so a renal covariate would",
        "not be expected to be informative."
      ),
      source_name        = "Creatinine clearance"
    )
  )

  compartmentData <- list(
    depot   = list(analyte = "brepocitinib", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "brepocitinib", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 379,
    n_studies      = 5,
    n_observations = 5541,
    age_range      = "18-75 years",
    age_median     = "42 years",
    weight_range   = "45-204 kg",
    weight_median  = "83 kg",
    sex_female_pct = 31,
    race_ethnicity = c(White = 82, `African American` = 5, Asian = 8, Other = 6),
    disease_state  = paste(
      "Healthy participants (24%) pooled with patients with plaque psoriasis",
      "(55%) and alopecia areata (20%)"
    ),
    renal_function = "Creatinine clearance median 124 mL/min (range 58-279); no renal-impairment cohort",
    dose_range     = paste(
      "1-200 mg oral. Phase 1: single doses of 1, 3, 10, 30, 100 and 200 mg,",
      "multiple doses of 10, 30, 100 and 175 mg once daily or 50 mg twice",
      "daily for 10 days, and 100 mg once daily for 10 days in Japanese",
      "participants. Phase 2: 30 or 60 mg once-daily induction followed by 10",
      "or 30 mg once-daily maintenance. Tablet (83%) or oral suspension (17%)"
    ),
    regions        = "Not reported by region; includes a dedicated phase 1 study in healthy Japanese participants (NCT03236493)",
    notes          = paste(
      "Baseline demographics from Hughes 2022 Table 2; study descriptions",
      "from Methods / Study Data. Three phase 1 studies (NCT02310750,",
      "NCT03236493, NCT03656952) and two phase 2 studies (NCT02969018 plaque",
      "psoriasis, NCT02974868 alopecia areata). Of the 5541 plasma",
      "observations, 10.8% were missing or below the 0.2 ng/mL limit of",
      "quantification; BLQ values were treated as missing with no imputation,",
      "and an M3 sensitivity analysis changed only the CL/F interindividual",
      "variability. Maintenance-period samples from the phase 2 patients",
      "assigned to 100 mg once-weekly maintenance were excluded because ~68%",
      "of the 511 samples in that period were missing, leaving dose-timing",
      "uncertain; their induction-period samples were retained. Estimation",
      "was by FOCE with interaction in NONMEM 7.4.1. This is the earlier of",
      "two published brepocitinib population PK analyses from the same",
      "programme; the later Maleki 2024 analysis (Maleki_2024_brepocitinib.R)",
      "pools nine studies and 775 subjects into a two-compartment model."
    )
  )

  ini({
    # ==================================================================
    # Structural parameters. The reference individual for these typical
    # values is a 70 kg non-Asian subject taking a tablet, dosed below
    # 175 mg and without a high-fat meal (Hughes 2022 Equations 6-10 with
    # every indicator set to 0).
    #
    # Naming note: Hughes 2022 writes CL/Frel and Vc/Frel for the apparent
    # clearance and apparent central volume BEFORE correction by relative
    # bioavailability (the sentence following Equation 10). Relative
    # bioavailability enters separately through f(depot), so lcl and lvc
    # below are the paper's CL/Frel and Vc/Frel exactly as tabulated.
    # ==================================================================
    lka     <- log(3.46)     ; label("Absorption rate constant (ka, 1/h)")                              # Table 3 'First-order absorption rate constant (ka)' = 3.46 (95% CI 3.04 to 3.88); Equation 7. Table 3 and the Abstract print the unit as 'h'; a first-order rate constant is 1/h, and 3.46 1/h is what reproduces the reported 5.04 h half-life population and the sub-1 h NCA tmax
    lcl     <- log(18.7)     ; label("Apparent oral clearance before Frel correction (CL/Frel, L/h)")   # Table 3 'Apparent clearance (CL/F)' = 18.7 (16.9 to 20.5); Equation 9
    lvc     <- log(136)      ; label("Apparent central volume before Frel correction (Vc/Frel, L)")     # Table 3 'Apparent volume of distribution (Vc/F)' = 136 (124 to 148); Equation 10
    ltlag   <- log(0.240)    ; label("Absorption lag time of the tablet formulation (h)")               # Table 3 'Absorption lag (Alag)' = 0.240 (0.234 to 0.246); Equation 6
    lfdepot <- fixed(log(1)) ; label("Relative bioavailability of the reference condition (unitless)")  # Equation 8: Frel,i = 1 * (1 - 0.283 * Fed) * (1 + 0.351 * Dose>=175mg). The leading 1 is a structural anchor with no Table 3 row, not an estimate

    # ==================================================================
    # Covariate effects. Hughes 2022 Equation 4 gives the categorical form
    # as P_i = theta_P * (1 + theta_PCOV) when the indicator is 1, so every
    # coefficient below is a FRACTIONAL change and is used as
    # (1 + coefficient * indicator). Equations 6-10 print the same
    # relationships with the sign moved outside the coefficient, e.g.
    # (1 - 0.699 * Fed); the Table 3 values carry the sign, so the two
    # printings agree.
    #
    # Equation 5 gives the continuous form as a power function referenced
    # to the covariate median; for body weight the reference is the
    # literature 70 kg rather than the observed 83 kg median (Full
    # Covariate Model section).
    # ==================================================================
    e_wt_cl                <- fixed(0.750)  ; label("Allometric body-weight exponent on CL/Frel, 70 kg reference (unitless)")     # Table 3 'Effect of weight on CL/F (70 kg reference)' = 0.750 (fixed); Equation 9
    e_wt_vc                <- fixed(1.00)   ; label("Allometric body-weight exponent on Vc/Frel, 70 kg reference (unitless)")     # Table 3 'Effect of weight on Vc/F (70 kg reference)' = 1.00 (fixed); Equation 10
    e_form_suspension_tlag <- fixed(-1.00)  ; label("Fractional change in the absorption lag for the oral suspension (unitless)") # Table 3 'Effect of suspension formulation on Alag' = -1.00 (fixed); Equation 6: Alag = 0.24 * (1 - 1 * Suspension), i.e. exactly no lag for the suspension
    e_fed_highfat_ka       <- -0.699        ; label("Fractional change in ka with a high-fat meal (unitless)")                    # Table 3 'Effect of high-fat meal on ka' = -0.699 (-0.899 to -0.499); Equation 7
    e_race_asian_cl        <- -0.243        ; label("Fractional change in CL/Frel in Asian participants (unitless)")              # Table 3 'Effect of Asian subjects on CL/F' = -0.243 (-0.363 to -0.123); Equation 9
    e_fed_highfat_fdepot   <- -0.283        ; label("Fractional change in relative bioavailability with a high-fat meal (unitless)")        # Table 3 'Effect of high-fat meal on Frel' = -0.283 (-0.389 to -0.177); Equation 8
    e_dose_high_fdepot     <- 0.351         ; label("Fractional change in relative bioavailability at doses of 175 mg and above (unitless)") # Table 3 'Effect of dose on Frel (dose >= 175 mg)' = 0.351 (0.194 to 0.508); Equation 8

    # ==================================================================
    # Interindividual variability: log-normal on CL/Frel and Vc/Frel with a
    # full 2x2 variance-covariance matrix (Base Model section: 'IIV was
    # included on CL/F and Vc/F with a full variance-covariance matrix').
    # Equations 9 and 10 apply the etas as exp(eta), and Equation 1 defines
    # P_i = theta_P * exp(eta_i) with variance omega_P^2.
    #
    # SCALE OF THE REPORTED %CV -- the load-bearing reading of Table 3.
    # The two IIV rows are labelled 'omega CL/F (% CV)' = 78.0 and
    # 'omega Vc/F (%CV)' = 60.5, i.e. the tabulated quantity IS omega,
    # printed as a percentage, NOT the log-normal coefficient of variation
    # sqrt(exp(omega^2) - 1). Two independent confirmations in the same
    # table: (1) the residual-error rows are labelled the same way --
    # 'Proportional RUV (phase 1; CV)' = 0.527 -- and Equation 2 defines
    # that quantity as theta_pro in ln(DV) = ln(IPRED) + theta_pro * eps
    # with Var(eps) fixed to 1, so it is unambiguously a log-scale SD that
    # the paper calls a CV; (2) the reported 95% CIs are exactly symmetric
    # about the point estimates on this scale (78.0 with 57.5 to 98.5, and
    # 60.5 with 27.3 to 93.5), as an SD-scale interval would be. So:
    #   var(etalcl) = 0.780^2 = 0.6084
    #   var(etalvc) = 0.605^2 = 0.366025
    #   cov         = rho * omega_CL * omega_Vc
    #               = 0.760 * 0.780 * 0.605 = 0.358644
    # For reference, the equivalent true log-normal CVs are
    # sqrt(exp(0.6084) - 1) = 91.5% and sqrt(exp(0.366025) - 1) = 66.3%.
    #
    # Element-by-element provenance of the lower-triangular block below.
    # These comments sit ABOVE the block rather than inline: a comment
    # inside a multi-line `~ c(...)` breaks rxode2's comment-to-label()
    # reparse with "unexpected ';'", and buildModelDb() does not catch it.
    #   0.6084   var(etalcl)          Table 3 'omega CL/F (% CV)' = 78.0
    #                                 (57.5 to 98.5) -> 0.780^2
    #   0.358644 cov(etalcl, etalvc)  Table 3 'rho CL/F-Vc/F' = 0.760
    #                                 (0.403 to 1.12) -> 0.760*0.780*0.605
    #   0.366025 var(etalvc)          Table 3 'omega Vc/F (%CV)' = 60.5
    #                                 (27.3 to 93.5) -> 0.605^2
    # ==================================================================
    etalcl + etalvc ~ c(
      0.6084,
      0.358644, 0.366025
    )

    # ==================================================================
    # Residual error. Equation 2 is a log-transform-both-sides proportional
    # model, ln(DV_ij) = ln(IPRED_ij) + theta_pro * eps_ij with eps ~ N(0, 1)
    # fixed, which is the nlmixr2 lnorm() endpoint with theta_pro as the
    # log-scale SD. The Table 3 row 'eps_res 1.00 (fixed)' is the NONMEM
    # SIGMA that Equation 2 fixes to 1 and has no nlmixr2 counterpart.
    #
    # Two magnitudes are estimated, one per study phase (Base Model
    # section); the stratum suffixes follow the pattern registered in
    # references/parameter-names.md and used by Friberg_2012_voriconazole.R.
    # ==================================================================
    expSdHv  <- 0.527 ; label("Log-scale residual SD in the phase 1 healthy-volunteer studies")  # Table 3 'Proportional RUV (phase 1;CV)' = 0.527 (0.444 to 0.610); Abstract: '52.7% CV in healthy volunteers'
    expSdPat <- 0.875 ; label("Log-scale residual SD in the phase 2 patient studies")            # Table 3 'Proportional RUV (phase 2;CV)' = 0.875 (0.812 to 0.938); Abstract: '87.5% CV in patients'
  })

  model({
    # ---- 1. Covariate-derived multipliers ----
    # Equation 5 power form, referenced to the literature 70 kg rather than
    # the observed 83 kg median.
    allomCl <- (WT / 70.0)^e_wt_cl
    allomVc <- (WT / 70.0)^e_wt_vc

    # Hughes 2022 Equation 6 is written on a suspension indicator. The
    # canonical register column is tablet-oriented, so the complement is
    # formed here and the printed equation is reproduced verbatim below.
    formSuspension <- 1.0 - FORM_TABLET

    # ---- 2. Individual parameters ----
    # Equation 9: CL/Frel,i = 18.7 * (BWT/70)^0.75 * (1 - 0.243*Asian) * exp(eta)
    cl <- exp(lcl + etalcl) * allomCl * (1.0 + e_race_asian_cl * RACE_ASIAN)

    # Equation 10: Vc/Frel,i = 136 * (BWT/70)^1 * exp(eta). No covariate
    # beyond allometry was retained on the volume.
    vc <- exp(lvc + etalvc) * allomVc

    # Equation 7: ka,i = 3.46 * (1 - 0.699*Fed). No IIV was estimated on ka.
    ka <- exp(lka) * (1.0 + e_fed_highfat_ka * FED_HIGHFAT)

    # Equation 6: Alag,i = 0.24 * (1 - 1*Suspension), so the lag is 0.240 h
    # for the tablet and exactly 0 for the suspension.
    tlag <- exp(ltlag) * (1.0 + e_form_suspension_tlag * formSuspension)

    # Equation 8: Frel,i = 1 * (1 - 0.283*Fed) * (1 + 0.351*Dose>=175mg)
    fdepot <- exp(lfdepot) *
      (1.0 + e_fed_highfat_fdepot * FED_HIGHFAT) *
      (1.0 + e_dose_high_fdepot * DOSE_HIGH)

    # ---- 3. Micro-constants ----
    kel <- cl / vc

    # ---- 4. ODE system: one compartment, first-order absorption ----
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # ---- 5. Dose input ----
    f(depot)    <- fdepot
    alag(depot) <- tlag

    # ---- 6. Observation and residual error ----
    # Amounts are mg and volumes are L, so central / vc is mg/L; the factor
    # of 1000 converts to ng/mL, the unit of the bioanalytical assay
    # (LLOQ 0.2 ng/mL, Methods / Study Data). The residual error is
    # scale-free, so this conversion does not disturb it.
    Cc <- 1000.0 * central / vc

    expSdCohort <- expSdHv * DIS_HEALTHY + expSdPat * (1.0 - DIS_HEALTHY)

    Cc ~ lnorm(expSdCohort)
  })
}
