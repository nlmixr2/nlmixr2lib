Hennig_2007_itraconazole <- function() {
  description <- "Two-compartment population PK model for oral itraconazole and its one-compartment hydroxy-itraconazole metabolite in adult cystic fibrosis patients (Hennig 2007), with first-order absorption from a depot, formulation-specific absorption rate constants and bioavailability for capsule vs. oral solution selected by the binary FORM_CAPSULE covariate, and a single absorption lag-time shared across both formulations. The fraction of itraconazole metabolised to hydroxy-itraconazole is fixed to 1; metabolite parameters are reported as CL_m/(F*f_m) and V_m/(F*f_m)."
  reference <- paste(
    "Hennig S, Waterhouse TH, Bell SC, France M, Wainwright CE, Miller H,",
    "Charles BG, Duffull SB. (2007).",
    "A D-optimal designed population pharmacokinetic study of oral",
    "itraconazole in adult cystic fibrosis patients.",
    "Br J Clin Pharmacol 63(4):438-450.",
    "doi:10.1111/j.1365-2125.2006.02778.x.",
    sep = " "
  )
  vignette <- "Hennig_2007_itraconazole"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    FORM_CAPSULE = list(
      description        = "Itraconazole oral formulation indicator (1 = Sporanox capsule, 0 = Sporanox oral solution)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (oral solution; F_rel fixed to 1 in the source model)",
      notes              = paste(
        "Per-dose-record indicator (cross-over: each subject received both",
        "formulations on two separate occasions 72 h apart).",
        "Selects between two formulation-specific typical-value absorption",
        "rate constants (k_a cap = 0.032 1/h, k_a sol = 0.125 1/h;",
        "Hennig 2007 Table 3) and applies the relative-bioavailability",
        "factor F_rel = 0.817 to the capsule arm only, with F = 1 fixed",
        "for the oral solution (Hennig 2007 Methods: 'Frel was determined",
        "by fixing Frel to 1 for administration of the oral solution and",
        "then estimating Frel for the capsule').",
        "The single absorption lag time t_lag = 19.3 min is shared across",
        "both formulations (Hennig 2007 Table 3, t_lag row).",
        "Source NONMEM column is `PREP` with PREP = 1 = capsule and",
        "PREP = 0 = oral solution; the canonical FORM_CAPSULE column has the",
        "same orientation."
      ),
      source_name        = "PREP"
    )
  )

  population <- list(
    n_subjects     = 30L,
    n_studies      = 1L,
    age_range      = "16-61 years",
    age_median     = "25 years",
    weight_range   = "46-86 kg",
    weight_median  = "57 kg",
    height_range   = "149-186 cm",
    height_median  = "169 cm",
    lbw_range      = "36-64 kg (Cheymol formula)",
    lbw_median     = "46 kg",
    sex_female_pct = 40,
    race_ethnicity = NULL,
    disease_state  = paste(
      "Adult cystic fibrosis patients in hospital for management of a chest",
      "exacerbation, recruited from the Department of Thoracic Medicine,",
      "The Prince Charles Hospital, Brisbane. Patients were not taking",
      "itraconazole for clinical indications at the time of the study and",
      "had no impaired hepatic function.",
      "Pancreatic insufficiency (typical of CF) was managed by enzyme",
      "supplementation with meals on study days; subjects on antireflux",
      "medications (e.g. proton pump inhibitors) timed those medications",
      "at least 2 h after the itraconazole dose to ensure capsule",
      "absorption."
    ),
    dose_range     = paste(
      "Single 200 mg oral itraconazole on each of two occasions 72 h apart",
      "(cross-over): two 100 mg Sporanox capsules and 20 mL of 10 mg/mL",
      "Sporanox oral solution. The simulation chapter of the paper extends",
      "to multiple-dose regimens (200, 400-600 mg twice daily, and",
      "300-1000 mg twice daily over 7 days) for target-attainment",
      "predictions, but the parameter estimates were obtained from the",
      "single-dose cross-over data."
    ),
    regions        = "Brisbane, Queensland, Australia (single-site).",
    notes          = paste(
      "Hennig 2007 Table 2 baseline demographics. Median (range) age 25",
      "(16-61) years, weight 57 (46-86) kg, height 169 (149-186) cm,",
      "lean body weight 46 (36-64) kg, 13 (7-21) co-medications per",
      "patient. 241 plasma samples obtained over 72 h post-dose; 46.0%",
      "of itraconazole and 27.8% of hydroxy-itraconazole concentrations",
      "were below the limit of detection (0.04 mg/L; assay LOQ 0.075",
      "mg/L). The race / ethnicity composition is not reported in the",
      "paper. Eight blood samples per subject were drawn within optimal",
      "sampling windows determined by D-optimal design."
    )
  )

  ini({
    # Structural absorption: separate first-order rate constants for capsule
    # and oral solution. The FORM_CAPSULE covariate selects the relevant typical
    # value (and IIV) inside model(); both etas exist for every subject in
    # the cross-over design but only one is "active" for any given dose
    # record. Source $PK block:
    #   IF (PREP.EQ.1) KA = THETA(3) * EXP(ETA(4))   ; capsule
    #   IF (PREP.EQ.0) KA = THETA(4) * EXP(ETA(5))   ; oral solution
    lka_cap <- log(0.032) ; label("Absorption rate Ka for capsule (1/h)")               # Hennig 2007 Table 3, k_a cap = 0.032 (RSE 46.7%)
    lka_sol <- log(0.125) ; label("Absorption rate Ka for oral solution (1/h)")         # Hennig 2007 Table 3, k_a sol = 0.125 (RSE 44.2%)

    # Apparent itraconazole disposition: 2-cmt with all elimination from
    # central described as the hydroxy-itraconazole formation pathway
    # (k20 = 0; fm = 1 fixed). Apparent values are CL_p/F and V_c/F where
    # F is the relative bioavailability of the formulation (F = 1 for the
    # oral solution, F = F_rel for the capsule).
    lcl <- log(31.5) ; label("Apparent itraconazole clearance CL_p/F (L/h)")            # Hennig 2007 Table 3, Cl_p = 31.5 (RSE 14.0%)
    lvc <- log(56.7) ; label("Apparent itraconazole central volume V_c/F (L)")          # Hennig 2007 Table 3, V_c = 56.7 (RSE 33.9%)
    lq  <- log(71.3) ; label("Apparent itraconazole inter-compartmental clearance Q/F (L/h)")  # Hennig 2007 Table 3, Q = 71.3 (RSE 35.3%)
    lvp <- log(2090) ; label("Apparent itraconazole peripheral volume V_per/F (L)")     # Hennig 2007 Table 3, V_per = 2090 (RSE 35.0%)

    # Apparent hydroxy-itraconazole disposition: 1-cmt with linear
    # elimination. CL_m and V_m are reported as CL_m/(F*f_m) and
    # V_m/(F*f_m); the fraction metabolised f_m is fixed to 1
    # (Hennig 2007 Methods: 'The fraction (fm) of itraconazole metabolised
    # to hydroxy-itraconazole has not been reported in the literature and
    # was therefore fixed to 1').
    lcl_ohi <- log(18.3) ; label("Apparent hydroxy-itraconazole clearance CL_m/(F*f_m) (L/h)")  # Hennig 2007 Table 3, CL_m = 18.3 (RSE 12.9%)
    lvc_ohi <- log(2.67) ; label("Apparent hydroxy-itraconazole volume V_m/(F*f_m) (L)")        # Hennig 2007 Table 3, V_m = 2.67 (RSE 49.8%)

    # Relative oral bioavailability of the capsule vs. oral solution
    # (F_solution fixed to 1) and a single absorption lag-time shared
    # across both formulations.
    lfdepot <- log(0.817)     ; label("Relative bioavailability F_rel of capsule (oral solution F = 1 fixed)")  # Hennig 2007 Table 3, F_rel = 0.817 (RSE 23.5%)
    llag    <- log(19.3 / 60) ; label("Absorption lag time (h; converted from t_lag = 19.3 min)")               # Hennig 2007 Table 3, t_lag = 19.3 min (RSE 1.68%)

    # Inter-individual variability. Final-model BSV CV% values from
    # Hennig 2007 Table 3 converted to log-normal omega^2 via
    # omega^2 = log(1 + CV^2). The example BQL .ctl bundled with the
    # paper used a $OMEGA BLOCK(2) on (etalcl, etalvc) but the
    # off-diagonal final estimate is not reported in the paper text or
    # Table 3, so the two etas are modelled here as independent (see the
    # validation vignette's Assumptions and deviations section).
    etalcl     ~ 0.04770   # log(1 + 0.221^2); Hennig 2007 Table 3 BSV Cl_p = 22.1 % CV
    etalvc     ~ 0.46838   # log(1 + 0.773^2); Hennig 2007 Table 3 BSV V_c = 77.3 % CV
    etalka_cap ~ 0.61222   # log(1 + 0.919^2); Hennig 2007 Table 3 BSV k_a cap = 91.9 % CV
    etalka_sol ~ 0.75607   # log(1 + 1.063^2); Hennig 2007 Table 3 BSV k_a sol = 106.3 % CV
    etalfdepot ~ 0.32787   # log(1 + 0.623^2); Hennig 2007 Table 3 BSV F_rel = 62.3 % CV

    # Residual error: proportional for both itraconazole and hydroxy-
    # itraconazole (Hennig 2007 Methods: 'A proportional error model was
    # shown to be sufficient to characterize the residual unexplained
    # variability for both the parent drug and the metabolite').
    propSd     <- 0.408 ; label("Proportional residual error on itraconazole (fraction)")          # Hennig 2007 Table 3, residual itraconazole 40.8 % CV
    propSd_ohi <- 0.479 ; label("Proportional residual error on hydroxy-itraconazole (fraction)")  # Hennig 2007 Table 3, residual hydroxy-itraconazole 47.9 % CV
  })

  model({
    # Formulation-specific first-order absorption: FORM_CAPSULE selects between
    # the capsule and oral-solution typical values and IIV deviations.
    ka_cap <- exp(lka_cap + etalka_cap)
    ka_sol <- exp(lka_sol + etalka_sol)
    ka     <- FORM_CAPSULE * ka_cap + (1 - FORM_CAPSULE) * ka_sol

    # Apparent itraconazole and hydroxy-itraconazole PK parameters.
    cl     <- exp(lcl + etalcl)
    vc     <- exp(lvc + etalvc)
    q      <- exp(lq)
    vp     <- exp(lvp)
    cl_ohi <- exp(lcl_ohi)
    vc_ohi <- exp(lvc_ohi)

    # Relative bioavailability and absorption lag-time. F = 1 (fixed) for
    # the oral solution; F = F_rel * exp(etalfdepot) for the capsule.
    fdepot     <- exp(lfdepot + etalfdepot)
    f(depot)   <- (1 - FORM_CAPSULE) + FORM_CAPSULE * fdepot
    lag(depot) <- exp(llag)

    # Two-compartment parent disposition with first-order metabolism to
    # hydroxy-itraconazole; one-compartment metabolite with linear
    # elimination. The cl/vc term out of central is the formation
    # clearance because f_m is fixed to 1, so all parent leaving central
    # is converted to metabolite (k20 = 0 in the source $PK block).
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - (cl / vc) * central -
                          (q / vc) * central + (q / vp) * peripheral1
    d/dt(peripheral1) <-  (q / vc) * central - (q / vp) * peripheral1
    d/dt(central_ohi) <-  (cl / vc) * central - (cl_ohi / vc_ohi) * central_ohi

    # Observation equations. Cc and Cc_ohi are plasma concentrations of
    # itraconazole and hydroxy-itraconazole in mg/L; both carry
    # independent proportional residual errors.
    Cc     <- central     / vc
    Cc_ohi <- central_ohi / vc_ohi

    Cc     ~ prop(propSd)
    Cc_ohi ~ prop(propSd_ohi)
  })
}
