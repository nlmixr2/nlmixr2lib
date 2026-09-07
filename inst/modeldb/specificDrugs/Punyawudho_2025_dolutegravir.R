Punyawudho_2025_dolutegravir <- function() {
  description <- "One-compartment population PK model for dolutegravir co-administered with rifampicin in Thai people living with HIV and tuberculosis, with lagged first-order absorption, allometric weight scaling, a negative exponential total-bilirubin effect on apparent clearance (-25.7% per mg/dL), and between-occasion variability on bioavailability and the absorption rate constant"
  reference <- paste(
    "Punyawudho B, Chanruang A, Ueaphongsukkit T, Gatechompol S, Ubolyam S,",
    "Cho YS, Shin JG, Avihingsanon A. The population pharmacokinetics of",
    "dolutegravir co-administered with rifampicin in Thai people living with",
    "HIV: Assessment of alternative dosing regimens. CPT Pharmacometrics Syst",
    "Pharmacol. 2025;14(1):95-104. doi:10.1002/psp4.13244.",
    sep = " "
  )
  vignette <- "Punyawudho_2025_dolutegravir"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric scaling on CL/F and V/F with a reference weight of 60 kg, the value the paper itself uses when quoting the typical clearance ('the estimated CL/F of DTG was 2.82 L/h among PLWH weighing 60 kg'); the cohort medians were 59.3 kg (once-daily arm) and 60.2 kg (twice-daily arm), Punyawudho 2025 Table 1. The exponents are fixed, not estimated: Punyawudho 2025 Methods, 'The allometric exponents were fixed to the values of 0.75 and 1 for CL/F and V/F, respectively'. Total body weight was selected over fat-free mass, which fit worse (Results, 'Allometric scaling with body weight outperformed fat-free mass').",
      source_name        = "body weight"
    ),
    TBILI = list(
      description        = "Total serum bilirubin",
      units              = "umol/L (SI canonical); the source paper reports mg/dL and the model converts inline via TBILI / 17.1",
      type               = "continuous",
      reference_category = NULL,
      notes              = "The only covariate retained in the final model (Punyawudho 2025 Results). Enters CL/F as an exponential effect centred on the cohort median of 0.38 mg/dL: exp(-0.297 * (TBILI_mgdL - 0.38)), so a 1 mg/dL rise in total bilirubin lowers CL/F by 1 - exp(-0.297) = 25.7%, the figure quoted in the Results. Median-centred continuous covariates with linear, power and exponential forms were screened (Methods); the exponential form is the one printed in the paper's CL/F equation. Mechanism: dolutegravir and bilirubin are both cleared by UGT1A1, so competition for the enzyme raises dolutegravir exposure in hyperbilirubinaemia (Discussion). The register's canonical unit is SI umol/L, so pass umol/L in the data; 1 mg/dL = 17.1 umol/L and the median 0.38 mg/dL is 6.50 umol/L. The paper's simulations stratify by the DAIDS hyperbilirubinaemia grades: normal 0.1-1.29, grade 1 1.3-1.89, grade 2 1.9-3.09, grade 3 3.1-6.09 and grade 4 >6.1 mg/dL (Methods, 'Simulations for evaluating optimal dosage regimens').",
      source_name        = "total bilirubin"
    ),
    OCC = list(
      description        = "Integer-valued occasion indicator for between-occasion-variability multiplexing",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = "Punyawudho 2025 Methods defines an occasion explicitly and states the count: 'each occasion was defined as a dosing event with at least one blood sample. As a result, there were two occasions: the pre-dose occasion and the post-dose occasion.' Two occasions are therefore encoded, with occasion 1 carrying the estimated variance and occasion 2 fixing it to the same value (the NONMEM $OMEGA BLOCK(1) SAME idiom used throughout this register -- see Kawuma_2023_dolutegravir.R, Chen_2023_nemonoxacin.R, Bihorel_2023_molnupiravir.R). Between-occasion variability was carried on the absorption parameters only, and it replaced rather than supplemented between-subject variability there (Methods, 'IOV was tested either following the inclusion of the IIV or substituting the IIV on absorption parameters (F and absorption rate constant; Ka)'; Results, 'The addition of the IOV on the absorption parameters (F and Ka) significantly improved the fit'). For single-occasion records pass OCC = 1; for a multi-dose simulation alternate OCC between 1 and 2 across successive doses so that consecutive administrations draw independent absorption behaviour.",
      source_name        = "occasion"
    )
  )

  compartmentData <- list(
    depot   = list(analyte = "dolutegravir", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "dolutegravir", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 40L,
    n_studies      = 1L,
    n_observations = 332L,
    age_range      = "21.6-60.5 years (arm means 37.5 and 35.6 years)",
    weight_range   = "41.1-86.0 kg (arm means 59.3 and 60.2 kg)",
    sex_female_pct = 12.5,
    disease_state  = "HIV/tuberculosis co-infection; treatment-naive people living with HIV newly diagnosed with tuberculosis, all receiving rifampicin-based anti-tuberculosis therapy",
    dose_range     = "dolutegravir 50 mg once daily with food (n = 20) or 50 mg twice daily without food (n = 20), each with rifampicin 450 mg daily (35-49 kg) or 600 mg daily (>=50 kg)",
    regions        = "Thailand (HIV-NAT, Thai Red Cross AIDS Research Centre, Bangkok)",
    notes          = "Cross-sectional analysis nested in NCT03731559. Intensive sampling at week 4: pre-dose and 1, 2, 4, 6, 8, 10 and 12 h post-dose in both arms, plus a 24 h sample in the once-daily arm. 332 dolutegravir concentrations from 40 participants entered the analysis; one participant whose concentrations were consistently below the 0.1 mg/L LLOQ was removed except for the pre-dose sample, and single below-LLOQ values were imputed at LLOQ/2 (Punyawudho 2025 Methods and Results). Serum creatinine medians were 0.891 and 0.895 mg/dL and total bilirubin medians 0.380 and 0.350 mg/dL (Table 1). EVERY participant received rifampicin, so the parameter values below are the co-administered (induced) values -- the model cannot separate the rifampicin effect from the baseline, and the authors say so explicitly (Discussion limitations: 'the pharmacokinetics of DTG in the absence of rifampicin and the effect of rifampicin on the pharmacokinetics of DTG cannot be determined'). Do not use this model for dolutegravir given without rifampicin."
  )

  ini({
    # =========================================================================
    # Structural parameters. All are apparent (oral) values for a typical
    # 60 kg individual with a total bilirubin of 0.38 mg/dL, taken co-
    # administered with rifampicin, from Punyawudho 2025 Table 2 column
    # 'NONMEM Point estimate'. The 95% CIs quoted below are the asymptotic
    # NONMEM intervals in the same table; a 1000-sample non-parametric
    # bootstrap (97.3% successful minimisations) agreed closely.
    # =========================================================================
    lcl <- log(2.82)
    label("Apparent oral clearance CL/F with rifampicin co-administration (L/h)")
    # Punyawudho 2025 Table 2 row 'CL/F (L/h)' = 2.82, %RSE 7.2%, 95% CI
    # 2.42-3.22 (bootstrap median 2.84, 95% CI 2.43-3.20). Also stated in the
    # Results: 'the estimated CL/F of DTG was 2.82 L/h among PLWH weighing
    # 60 kg and having a total bilirubin of 0.38 mg/dL'.

    lvc <- log(19.8)
    label("Apparent volume of distribution V/F (L)")
    # Punyawudho 2025 Table 2 row 'V/F (L)' = 19.8, %RSE 8.2%, 95% CI
    # 16.6-22.9 (bootstrap median 19.9, 95% CI 17.0-23.5).

    lka <- log(1.41)
    label("First-order absorption rate constant ka (1/h)")
    # Punyawudho 2025 Table 2 row 'Ka (h-1)' = 1.41, %RSE 38.6%, 95% CI
    # 0.344-2.48 (bootstrap median 1.47, 95% CI 0.344-2.47).

    ltlag <- log(0.562)
    label("Absorption lag time (h)")
    # Punyawudho 2025 Table 2 row 'Lag time (h)' = 0.562, %RSE 33.6%, 95% CI
    # 0.192-0.932 (bootstrap median 0.563, 95% CI 0.189-0.854). Adding the lag
    # time improved the fit by dOFV = -15.9 (Results).

    lfdepot <- fixed(log(1))
    label("Relative bioavailability F (unitless)")
    # Punyawudho 2025 Table 2 row 'F' = '1 (fixed)'; also Results, 'The
    # bioavailability (F) was fixed to 1'. Held at 1 because the study has no
    # intravenous reference arm; the between-occasion variability below is
    # nevertheless estimated around it.

    # =========================================================================
    # Covariate effects.
    # =========================================================================
    e_wt_cl <- fixed(0.75)
    label("Allometric exponent on apparent clearance (unitless)")
    # Punyawudho 2025 Methods: 'The allometric exponents were fixed to the
    # values of 0.75 and 1 for CL/F and V/F, respectively.' The reference
    # weight of 60 kg is the one printed in the paper's CL/F equation
    # (Results): CL/F = 2.82 x (body weight / 60)^0.75 x
    # exp(-0.297 x (total bilirubin - 0.38)).

    e_wt_vc <- fixed(1)
    label("Allometric exponent on apparent volume of distribution (unitless)")
    # Punyawudho 2025 Methods (fixed at 1 for V/F). The paper prints only the
    # CL/F equation, but the Methods statement covers both parameters and the
    # 60 kg reference is shared.

    e_tbili_cl <- -0.297
    label("Exponential total-bilirubin coefficient on CL/F (per mg/dL)")
    # Punyawudho 2025 Table 2 row 'CL-Bilirubin' = -0.297, %RSE 9.7%, 95% CI
    # -0.347 to -0.240 (bootstrap median -0.295, 95% CI -0.512 to -0.118), and
    # the Results equation exp(-0.297 x (total bilirubin - 0.38)). Self-check:
    # 1 - exp(-0.297) = 0.2570, matching the Results sentence 'An elevation of
    # 1 mg/dL in total bilirubin decreased the CL/F of DTG by 25.7%'.

    # =========================================================================
    # Between-subject variability. Punyawudho 2025 Table 2 section
    # 'Inter-individual/Inter-occasion variability (%CV)'. The random effects
    # are log-normal (Methods: 'The inter-individual variability (IIV) and
    # inter-occasion variability (IOV) were assumed to be log-normally
    # distributed'). The paper does not print the formula behind its '%CV'
    # column, so omega^2 = log(1 + (CV/100)^2) is used here -- the exact
    # log-normal coefficient of variation CV = sqrt(exp(omega^2) - 1) rather
    # than the sqrt(omega^2) approximation. The arbitration is the vignette
    # section 'Which %CV convention?', which re-solves one cohort under both
    # readings with common random numbers and scores each against the paper's
    # own Table 3 target-attainment percentages: the exact reading is closer
    # on both target columns (mean absolute deviation over the 30 cells of
    # 0.67 vs 0.70 percentage points on the IC90 column and 2.16 vs 2.19 on
    # the EC90 column). The margin is small -- BOTH readings reproduce
    # Table 3 well, and the choice only moves the bioavailability IOV variance
    # from 0.599 to 0.821 -- so this is the one convention question worth
    # re-checking against the NONMEM control stream, should it become
    # available. It does not affect propSd: for a proportional error model the
    # tabulated %CV is the standard deviation directly.
    #
    # BSV was retained on clearance only. The Results state that 'The IIV of
    # V/F and lag time could not be precisely estimated, thus the V/F and lag
    # time were estimated without its IIV', and Table 2 prints no IIV row for
    # ka or F because the between-occasion terms replaced them there.
    # =========================================================================
    etalcl ~ 0.0358313
    # Punyawudho 2025 Table 2 row 'IIV-CL' = 19.1 %CV, %RSE 14.58%, 95% CI
    # 12.8-23.7 (bootstrap median 18.9, 95% CI 12.2-24.1);
    # omega^2 = log(1 + 0.191^2) = 0.0358313.

    # =========================================================================
    # Between-occasion variability on the two absorption parameters. One
    # variance per parameter is reported and is shared by both occasions, so
    # occasion 1 carries the estimated value and occasion 2 fixes it to the
    # same number (NONMEM $OMEGA BLOCK(1) SAME). Two occasions, not more:
    # Punyawudho 2025 Methods, 'there were two occasions: the pre-dose
    # occasion and the post-dose occasion'.
    # =========================================================================
    etaiov_fdepot_1 ~ 0.5992957
    # Punyawudho 2025 Table 2 row 'IOV-F1' = 90.6 %CV, %RSE 20.3%, 95% CI
    # 86.9-94.1 (bootstrap median 93.3, 95% CI 70.0-134);
    # omega^2 = log(1 + 0.906^2) = 0.5992957 (estimated, occasion 1).
    etaiov_fdepot_2 ~ fixed(0.5992957)  # SAME-equivalent: equal to the occasion-1 variance

    etaiov_ka_1 ~ 0.2393318
    # Punyawudho 2025 Table 2 row 'IOV-Ka' = 52.0 %CV, %RSE 18.1%, 95% CI
    # 28.1-67.9 (bootstrap median 48.9, 95% CI 24.9-72.5);
    # omega^2 = log(1 + 0.520^2) = 0.2393318 (estimated, occasion 1).
    etaiov_ka_2 ~ fixed(0.2393318)  # SAME-equivalent: equal to the occasion-1 variance

    # =========================================================================
    # Residual unexplained variability. Proportional only (Punyawudho 2025
    # Methods: 'The residual unexplained variability (RUV) was characterized by
    # proportional error model'); no additive component is reported.
    # =========================================================================
    propSd <- 0.156
    label("Proportional residual error (fraction)")
    # Punyawudho 2025 Table 2 row 'RUV prop' = 15.6 %CV, %RSE 16.4%, 95% CI
    # 9.27-20.0 (bootstrap median 15.3, 95% CI 9.84-20.3). For a proportional
    # error model the tabulated %CV is the standard deviation of the
    # proportional term directly, independently of the log-normal convention
    # question above.
  })

  model({
    # 1. Decompose the integer occasion column into binary indicators and
    #    multiplex the between-occasion etas onto bioavailability and the
    #    absorption rate constant. For single-occasion records pass OCC = 1.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)

    iov_fdepot <- oc1 * etaiov_fdepot_1 + oc2 * etaiov_fdepot_2
    iov_ka <- oc1 * etaiov_ka_1 + oc2 * etaiov_ka_2

    # 2. The canonical TBILI column is SI umol/L; the source paper's covariate
    #    equation, its median of 0.38 and its coefficient of -0.297 are all on
    #    the US-convention mg/dL scale, so convert inline (1 mg/dL = 17.1
    #    umol/L). The median 0.38 mg/dL is 6.50 umol/L.
    tbili_mgdL <- TBILI / 17.1

    # 3. Individual parameters. Punyawudho 2025 Results:
    #      CL/F = 2.82 x (body weight / 60)^0.75
    #                  x exp(-0.297 x (total bilirubin - 0.38))
    #    with the volume scaled allometrically on the same 60 kg reference at
    #    an exponent of 1 (Methods).
    cl <- exp(lcl + etalcl) * (WT / 60)^e_wt_cl *
      exp(e_tbili_cl * (tbili_mgdL - 0.38))
    vc <- exp(lvc) * (WT / 60)^e_wt_vc
    ka <- exp(lka + iov_ka)
    tlag <- exp(ltlag)
    fdepot <- exp(lfdepot + iov_fdepot)

    # 4. Micro-constant.
    kel <- cl / vc

    # 5. One-compartment disposition with lagged first-order absorption from an
    #    oral depot and first-order elimination from the central compartment
    #    (Punyawudho 2025 Results: 'a one-compartment model with first-order
    #    absorption and elimination'; 'The addition of lag time for describing
    #    the delayed absorption improved the model fit'). A second compartment
    #    was tested and did not improve the fit.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # 6. Relative bioavailability and absorption lag time on the depot.
    f(depot) <- fdepot
    alag(depot) <- tlag

    # 7. Observation and residual error. Dose in mg with vc in L gives Cc in
    #    mg/L, the unit the source reports concentrations and its 0.064 mg/L
    #    in vitro protein-adjusted IC90 and 0.3 mg/L in vivo EC90 targets in.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
