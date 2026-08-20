Lalagkas_2023_ganciclovir <- function() {
  description <- paste(
    "Two-compartment population PK model for ganciclovir after IV ganciclovir",
    "and oral valganciclovir in Caucasian solid-organ transplant recipients with",
    "established CMV infection (Lalagkas 2023), with first-order absorption after",
    "a lag time, estimated bioavailability, allometric body-weight scaling, and a",
    "CKD-EPI eGFR power effect on clearance.",
    "Parameters transcribed from the Yang 2023 ganciclovir / valganciclovir",
    "population-PK model repository (Table 3), not from the primary publication;",
    "re-verify against Lalagkas 2023 when the primary is obtained.",
    sep = " "
  )
  reference <- paste(
    "Lalagkas PN, Iliou J, Rigo R, Miarons M, Fernandez-Alarcon B, Bestard O,",
    "Cruzado JM, Melilli E, Torras J, Grinyo JM, et al. Comparison of Three Renal",
    "Function Formulas for Ganciclovir/Valganciclovir Dose Individualization in",
    "CMV-Infected Solid Organ Transplantation Patients Using a Population Approach.",
    "Clin Pharmacokinet. 2023;62(6):861-880. doi:10.1007/s40262-023-01245-3.",
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
        "Allometric scaling referenced to 70 kg: exponent 0.75 on CL and Q,",
        "exponent 1 (linear) on Vc and Vp. Cohort median 68 kg, range 43-131 kg",
        "(Yang 2023 Table 2).",
        sep = " "
      ),
      source_name        = "BW"
    ),
    CRCL = list(
      description        = "CKD-EPI-estimated glomerular filtration rate",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Renal function estimated with the Chronic Kidney Disease Epidemiology",
        "Collaboration (CKD-EPI) equation, which returns mL/min/1.73 m^2. Power",
        "effect on CL referenced to 55 mL/min/1.73 m^2. The primary publication",
        "compared three renal-function formulas (CKD-EPI, and two others); the",
        "Yang 2023 repository transcribes the CKD-EPI variant, which is the model",
        "encoded here.",
        sep = " "
      ),
      source_name        = "CKD-EPI"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 60L,
    n_studies      = 1L,
    n_observations = 640L,
    age_median     = "57 years (range 22-78)",
    weight_median  = "68 kg (range 43-131)",
    sex_female_pct = 35,
    race_ethnicity = "Caucasian (all subjects).",
    disease_state  = paste(
      "Caucasian patients with established CMV infection undergoing allogeneic",
      "solid-organ transplantation (kidney, liver and heart).",
      sep = " "
    ),
    dose_range     = paste(
      "IV ganciclovir for CMV infection at 2.5 mg/kg/12 h, 2.5 mg/kg/24 h, or",
      "1.25 mg/kg/24 h. Oral valganciclovir for CMV infection at 450 mg every 12,",
      "24 or 48 h; for prophylaxis 450 mg every 24, 48 or 84 h.",
      sep = " "
    ),
    regions        = "Spain (prospective study).",
    bioassay       = "HPLC, LLOQ 0.5 mg/L.",
    notes          = paste(
      "Demographics and dosing from Yang 2023 Table 2 (the ganciclovir /",
      "valganciclovir model-repository review). Sampling was both intensive",
      "(0 pre-dose, 0.5, 1, 1.5, 2, 3, 4, 6, 8, 10 and 12 h post dose) and sparse",
      "(0.5-1.5, 4-5 and 6-8 h post dose). Externally validated in 22 subjects.",
      "This study converted valganciclovir doses to their ganciclovir-equivalent",
      "content by multiplying by 0.72 (the ratio of the ganciclovir and",
      "valganciclovir molecular weights); oral doses supplied to this model must",
      "therefore already be expressed as ganciclovir equivalents.",
      sep = " "
    )
  )

  ini({
    # Structural PK -- Yang 2023 Table 3, Lalagkas et al. (2023) row. Reference
    # subject: WT = 70 kg, CKD-EPI eGFR = 55 mL/min/1.73 m^2. Clearances in L/h,
    # volumes in L, ka in 1/h, lag time in h. All values are typical (population)
    # values of the final model; the review reports no standard errors.
    lcl   <- log(6.93) ; label("Clearance at WT = 70 kg and CKD-EPI = 55 mL/min/1.73 m^2 (CL, L/h)")  # Yang 2023 Table 3 (Lalagkas 2023): CL = 6.93 * (CKD-EPI/55)^0.817 * (BW/70)^0.75
    lvc   <- log(43.1) ; label("Central volume of distribution at WT = 70 kg (Vc, L)")                 # Yang 2023 Table 3 (Lalagkas 2023): Vc = 43.1 * (BW/70)
    lq    <- log(9.23) ; label("Inter-compartmental clearance at WT = 70 kg (Q, L/h)")                 # Yang 2023 Table 3 (Lalagkas 2023): Q = 9.23 * (BW/70)^0.75
    # Vp = 219 L gives Vss ~262 L at 70 kg -- by far the largest volume among the
    # 16 models in the Yang 2023 repository, and well above the ~0.7 L/kg
    # literature volume of distribution for ganciclovir. Transcribed verbatim from
    # Table 3; it produces a long terminal half-life (~34 h in the validation
    # vignette) while contributing little to AUC over a 24 h interval. Flagged in
    # the validation vignette for re-verification against the primary.
    lvp   <- log(219)  ; label("Peripheral volume of distribution at WT = 70 kg (Vp, L)")              # Yang 2023 Table 3 (Lalagkas 2023): Vp = 219 * (BW/70)
    lka   <- log(0.766); label("First-order oral absorption rate constant (ka, 1/h)")                  # Yang 2023 Table 3 (Lalagkas 2023): Ka = 0.766
    ltlag <- log(0.331); label("Absorption lag time (Tlag, h)")                                        # Yang 2023 Table 3 (Lalagkas 2023): Tlag = 0.331

    # Oral bioavailability of ganciclovir from valganciclovir. The review reports
    # F = 0.699 on the linear scale with 16.6% between-subject variability and
    # states that BSV was exponential in every included study, so F carries an
    # exponential eta here (see the vignette Assumptions and deviations section
    # for the F > 1 caveat this parameterisation permits).
    lfdepot <- log(0.699); label("Oral bioavailability of ganciclovir from valganciclovir (F, fraction)")  # Yang 2023 Table 3 (Lalagkas 2023): F = 0.699

    # Covariate effects. The CKD-EPI exponent 0.817 is a non-canonical value and
    # is therefore an estimated effect; the body-weight exponents 0.75 (CL, Q) and
    # 1 (Vc, Vp) are the canonical allometric values imposed by the authors and
    # are encoded as fixed.
    e_crcl_cl <- 0.817      ; label("Power exponent of CKD-EPI eGFR on CL (unitless; reference 55 mL/min/1.73 m^2)") # Yang 2023 Table 3 (Lalagkas 2023): (CKD-EPI/55)^0.817
    e_wt_cl_q <- fixed(0.75); label("Allometric exponent of WT on CL and Q (unitless; reference 70 kg)")             # Yang 2023 Table 3 (Lalagkas 2023): (BW/70)^0.75 on CL and Q
    e_wt_vc_vp <- fixed(1)  ; label("Allometric exponent of WT on Vc and Vp (unitless; linear, reference 70 kg)")    # Yang 2023 Table 3 (Lalagkas 2023): (BW/70) on Vc and Vp

    # Between-subject variability. Yang 2023 Methods define the reported
    # percentage as %CV = sqrt(omega^2) * 100%, i.e. omega equals the tabulated
    # percentage divided by 100, so the variance is (BSV% / 100)^2. Q and Tlag
    # carry no BSV in the source table.
    etalcl     ~ 0.089401  # Yang 2023 Table 3 (Lalagkas 2023): BSV CL = 29.9%  -> 0.299^2
    etalvc     ~ 0.130321  # Yang 2023 Table 3 (Lalagkas 2023): BSV Vc = 36.1%  -> 0.361^2
    etalvp     ~ 1.069156  # Yang 2023 Table 3 (Lalagkas 2023): BSV Vp = 103.4% -> 1.034^2
    etalka     ~ 0.208849  # Yang 2023 Table 3 (Lalagkas 2023): BSV ka = 45.7%  -> 0.457^2
    etalfdepot ~ 0.027556  # Yang 2023 Table 3 (Lalagkas 2023): BSV F  = 16.6%  -> 0.166^2

    # Residual unexplained variability: combined proportional + additive.
    propSd <- 0.282; label("Proportional residual error (fraction)")  # Yang 2023 Table 3 (Lalagkas 2023): 28.2% proportional error
    addSd  <- 0.237; label("Additive residual error (mg/L)")          # Yang 2023 Table 3 (Lalagkas 2023): 0.237 mg/L additive error
  })

  model({
    # Individual parameters: allometric body-weight scaling on all four
    # disposition parameters plus a CKD-EPI eGFR power effect on CL.
    cl <- exp(lcl + etalcl) * (CRCL / 55)^e_crcl_cl * (WT / 70)^e_wt_cl_q
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc_vp
    q  <- exp(lq)           * (WT / 70)^e_wt_cl_q
    vp <- exp(lvp + etalvp) * (WT / 70)^e_wt_vc_vp
    ka <- exp(lka + etalka)
    tlag   <- exp(ltlag)
    fdepot <- exp(lfdepot + etalfdepot)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Oral valganciclovir doses enter `depot` (already expressed as
    # ganciclovir-equivalent milligrams, i.e. valganciclovir mg * 0.72);
    # IV ganciclovir doses go directly to `central`.
    f(depot)    <- fdepot
    alag(depot) <- tlag

    # Dose in mg, volume in L -> concentration in mg/L.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
