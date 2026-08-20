Franck_2021_ganciclovir <- function() {
  description <- paste(
    "Two-compartment population PK model for ganciclovir after IV ganciclovir and",
    "oral valganciclovir in pediatric solid-organ and stem-cell transplant",
    "recipients (Franck 2021), with first-order absorption after a lag time,",
    "estimated bioavailability, allometric body-weight scaling, a creatinine-clearance",
    "power effect on clearance, and a purely additive residual error.",
    "Parameters transcribed from the Yang 2023 ganciclovir / valganciclovir",
    "population-PK model repository (Table 3), not from the primary publication;",
    "re-verify against Franck 2021 when the primary is obtained.",
    sep = " "
  )
  reference <- paste(
    "Franck B, Autmizguine J, Marquet P, Ovetchkine P, Woillard JB. Population",
    "pharmacokinetics of ganciclovir and valganciclovir in paediatric solid organ",
    "and stem cell transplant recipients. Br J Clin Pharmacol. 2021.",
    "Cited as reference 9 of the Yang 2023 model repository.",
    "Parameters transcribed from Yang W, Mak W, Gwee A, Gu M, Wu Y, Shi Y, He Q,",
    "Xiang X, Han B, Zhu X. Establishment and Evaluation of a Parametric Population",
    "Pharmacokinetic Model Repository for Ganciclovir and Valganciclovir.",
    "Pharmaceutics. 2023;15(7):1801. doi:10.3390/pharmaceutics15071801 (Table 3).",
    sep = " "
  )
  vignette <- "Yang_2023_ganciclovir_model_repository"
  units    <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Derived mechanically; verified = FALSE means it has
  # NOT been checked against the source paper.
  compartmentData <- list(
    depot       = list(analyte = "ganciclovir", units = "mg", specimen = "administration site", verified = FALSE),
    central     = list(analyte = "ganciclovir", units = "mg", specimen = "plasma", verified = FALSE),
    peripheral1 = list(analyte = "ganciclovir", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometric scaling referenced to the cohort median 26.7 kg: exponent 0.75",
        "on CL, exponent 1 (linear) on Vc and Vp. Q is not weight-scaled in this",
        "model. Cohort median 26.7 kg, range 5.96-87 kg (Yang 2023 Table 2).",
        sep = " "
      ),
      source_name        = "BW"
    ),
    CRCL = list(
      description        = "Creatinine clearance, BSA-normalized",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Yang 2023 Table 3 footnote defines CrCL as creatinine clearance in",
        "mL/min/1.73 m^2 (BSA-normalized), distinct from the raw-mL/min `CLcr`",
        "used by other studies in the same repository. Power effect on CL",
        "referenced to 149.8 mL/min/1.73 m^2. This is the covariate the Yang 2023",
        "MAP-Bayesian AUC calculator requires alongside body weight",
        "(Yang 2023 Section 3.4.2).",
        sep = " "
      ),
      source_name        = "CrCL"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 50L,
    n_studies      = 1L,
    n_observations = 580L,
    age_median     = "7.5 years (range 0.5-17.4)",
    weight_median  = "26.7 kg (range 5.96-87)",
    sex_female_pct = 40,
    race_ethnicity = "Not reported; ethnicity was tested as a covariate and not retained.",
    disease_state  = "Pediatric solid-organ transplant (SOT) and stem-cell transplant (SCT) recipients.",
    dose_range     = paste(
      "Pre-emptive approach for prevention of CMV disease: IV ganciclovir",
      "5 mg/kg/12 h or oral valganciclovir 10 mg/kg/12 h.",
      sep = " "
    ),
    regions        = "Canada (retrospective study).",
    bioassay       = "HPLC, LLOQ 0.039 mg/L.",
    notes          = paste(
      "Demographics and dosing from Yang 2023 Table 2. Intensive sampling for both",
      "the ganciclovir and valganciclovir arms. Covariates tested: weight, BSA, sex,",
      "age, ethnicity, transplant type, formulation, serum creatinine, urea and",
      "CrCL; weight and CrCL were retained (Yang 2023 Table 4). This study did NOT",
      "convert valganciclovir doses to ganciclovir equivalents, so oral doses are",
      "supplied as valganciclovir milligrams and the estimated bioavailability",
      "F = 0.43 absorbs the molecular-weight ratio.",
      "Yang 2023 used this model as the basis of its MAP-Bayesian AUC calculator.",
      sep = " "
    )
  )

  ini({
    # Structural PK -- Yang 2023 Table 3, Franck et al. (2021) row. Reference
    # subject: WT = 26.7 kg, CrCL = 149.8 mL/min/1.73 m^2. Clearances in L/h,
    # volumes in L, ka in 1/h, lag time in h.
    lcl   <- log(6.9) ; label("Clearance at WT = 26.7 kg and CrCL = 149.8 mL/min/1.73 m^2 (CL, L/h)") # Yang 2023 Table 3 (Franck 2021): CL = 6.9 * (BW/26.7)^0.75 * (CrCL/149.8)^0.88
    lvc   <- log(9.7) ; label("Central volume of distribution at WT = 26.7 kg (Vc, L)")               # Yang 2023 Table 3 (Franck 2021): Vc = 9.7 * (BW/26.7)
    lq    <- log(10.9); label("Inter-compartmental clearance (Q, L/h; not weight-scaled)")            # Yang 2023 Table 3 (Franck 2021): Q = 10.9
    lvp   <- log(7.6) ; label("Peripheral volume of distribution at WT = 26.7 kg (Vp, L)")            # Yang 2023 Table 3 (Franck 2021): Vp = 7.6 * (BW/26.7)
    lka   <- log(0.73); label("First-order oral absorption rate constant (ka, 1/h)")                  # Yang 2023 Table 3 (Franck 2021): Ka = 0.73
    ltlag <- log(0.33); label("Absorption lag time (Tlag, h)")                                        # Yang 2023 Table 3 (Franck 2021): Tlag = 0.33

    # Oral bioavailability of ganciclovir from a valganciclovir milligram dose
    # (no molecular-weight conversion was applied by this study).
    lfdepot <- log(0.43); label("Oral bioavailability of ganciclovir from valganciclovir (F, fraction)") # Yang 2023 Table 3 (Franck 2021): F = 0.43

    # Covariate effects. The CrCL exponent 0.88 is a non-canonical estimated
    # value; the body-weight exponents 0.75 (CL) and 1 (Vc, Vp) are the canonical
    # allometric values imposed by the authors and are encoded as fixed.
    e_crcl_cl  <- 0.88      ; label("Power exponent of CrCL on CL (unitless; reference 149.8 mL/min/1.73 m^2)") # Yang 2023 Table 3 (Franck 2021): (CrCL/149.8)^0.88
    e_wt_cl    <- fixed(0.75); label("Allometric exponent of WT on CL (unitless; reference 26.7 kg)")           # Yang 2023 Table 3 (Franck 2021): (BW/26.7)^0.75
    e_wt_vc_vp <- fixed(1)   ; label("Allometric exponent of WT on Vc and Vp (unitless; linear, reference 26.7 kg)") # Yang 2023 Table 3 (Franck 2021): (BW/26.7) on Vc and Vp

    # Between-subject variability. Yang 2023 Methods: %CV = sqrt(omega^2) * 100%,
    # so variance = (BSV% / 100)^2. Q, Vp and Tlag carry no BSV in the source table.
    etalcl     ~ 0.439569  # Yang 2023 Table 3 (Franck 2021): BSV CL = 66.3% -> 0.663^2
    etalvc     ~ 0.589824  # Yang 2023 Table 3 (Franck 2021): BSV Vc = 76.8% -> 0.768^2
    etalka     ~ 0.700569  # Yang 2023 Table 3 (Franck 2021): BSV ka = 83.7% -> 0.837^2
    etalfdepot ~ 0.310249  # Yang 2023 Table 3 (Franck 2021): BSV F  = 55.7% -> 0.557^2

    # Residual unexplained variability: additive only (the single study in this
    # repository with a pure additive residual error).
    addSd <- 0.98; label("Additive residual error (mg/L)")  # Yang 2023 Table 3 (Franck 2021): 0.98 mg/L additive error
  })

  model({
    cl <- exp(lcl + etalcl) * (WT / 26.7)^e_wt_cl * (CRCL / 149.8)^e_crcl_cl
    vc <- exp(lvc + etalvc) * (WT / 26.7)^e_wt_vc_vp
    q  <- exp(lq)
    vp <- exp(lvp)           * (WT / 26.7)^e_wt_vc_vp
    ka <- exp(lka + etalka)
    tlag   <- exp(ltlag)
    fdepot <- exp(lfdepot + etalfdepot)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Oral valganciclovir doses enter `depot` as valganciclovir milligrams;
    # IV ganciclovir doses go directly to `central`.
    f(depot)    <- fdepot
    alag(depot) <- tlag

    # Dose in mg, volume in L -> concentration in mg/L.
    Cc <- central / vc
    Cc ~ add(addSd)
  })
}
