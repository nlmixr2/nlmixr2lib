Quintairos_2021_tacrolimus <- function() {
  description <- paste(
    "Two-compartment population PK model for whole-blood tacrolimus after oral",
    "twice-daily dosing in adult de novo kidney transplant recipients during the",
    "first 6 months post-transplant (Quintairos 2021). First-order absorption with",
    "an absorption lag time and first-order elimination; fixed-exponent allometric",
    "scaling on body weight (0.75 on CL/F and Q/F; 1 on Vc/F and Vp/F); no other",
    "covariates retained. Between-subject variability on CL/F, Q/F and Vc/F, small",
    "fixed-variance random effects on ka, Vp/F and the lag time, and an additive",
    "residual error on log-transformed concentrations (log-normal)."
  )
  reference <- paste(
    "Quintairos L, Colom H, Millan O, Fortuna V, Espinosa C, Guirado L, Budde K,",
    "Sommerer C, Lizana A, Lopez-Pua Y, Brunet M. Early prognostic performance of",
    "miR155-5p monitoring for the risk of rejection: Logistic regression with a",
    "population pharmacokinetic approach in adult kidney transplant patients.",
    "PLoS ONE. 2021;16(1):e0245880. doi:10.1371/journal.pone.0245880.",
    "Parameter values from Table 4; random-effect and residual-error structure from",
    "the S1 Appendix NONMEM control stream ('Tacrolimus PK model').",
    sep = " "
  )
  vignette <- "Quintairos_2021_tacrolimus_mycophenolic_acid"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-varying (the deposited dataset S1 Table updates WGT at each visit).",
        "Reference weight 70 kg (Table 4 units 'L/h/70 kg', 'L/70 kg'). Exponent",
        "0.75 on CL/F and Q/F and 1 on Vc/F and Vp/F, both fixed a priori (Methods",
        "'Pharmacokinetic models'). ka and the lag time are not scaled. Cohort",
        "median 73 kg (IQR 62.9-86.8; Table 1)."
      ),
      source_name = "WGT"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "tacrolimus", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "tacrolimus", units = "mg", specimen = "whole blood", verified = TRUE),
    peripheral1 = list(analyte = "tacrolimus", units = "mg", specimen = "whole blood", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 58L,
    n_studies = 1L,
    n_observations = "1102 tacrolimus whole-blood concentrations",
    age_range = "Median 48 years (IQR 38-58); recipients older than 70 years were excluded",
    age_median = "48 years",
    weight_range = "Median 73 kg (IQR 62.9-86.8)",
    weight_median = "73 kg",
    sex_female_pct = 34.5,
    race_ethnicity = "Caucasian only (Discussion)",
    disease_state = paste(
      "Adult de novo kidney transplant recipients (living donor 30, deceased donor 28)",
      "in the first 6 months post-transplant; median GFR 44 mL/min (IQR 15-55).",
      "8 of 58 patients (14%) had biopsy-proven acute rejection."
    ),
    dose_range = paste(
      "Oral tacrolimus (Prograf) twice daily, therapeutic-drug-monitoring adjusted;",
      "mean dose 14.6 mg at week 1 falling to 5.29 mg at month 6 (Table 2)."
    ),
    regions = "Germany (Charite Berlin, Heidelberg) and Spain (Fundacio Puigvert, Barcelona)",
    co_medication = "Mycophenolate mofetil (Myfenax), methylprednisolone, basiliximab induction (2 x 20 mg)",
    notes = paste(
      "Pharmacokinetic sub-study (58 of 80 patients) of a European multicentre",
      "prospective observational study (EudraCT 2013-001817-33). Sampling at 0, 0.5,",
      "1, 1.5, 2, 3, 4, 6, 8 and 12 h post-dose at week 1 and at 0, 1.5, 2 and 4 h at",
      "months 1, 2, 3 and 6. Baseline demographics in Table 1."
    )
  )

  ini({
    # Structural parameters: Quintairos 2021 Table 4, tacrolimus column
    # (reference 70 kg).
    lcl <- log(16.5); label("Apparent clearance CL/F for a 70 kg patient (L/h)")  # Table 4 tacrolimus CL = 16.5 L/h/70 kg (RSE 10%)
    lvc <- log(311); label("Apparent central volume Vc/F for a 70 kg patient (L)")  # Table 4 tacrolimus VC = 311 L/70 kg (RSE 9%)
    lq <- log(20.5); label("Apparent intercompartmental clearance Q/F for a 70 kg patient (L/h)")  # Table 4 tacrolimus Q = 20.5 L/h/70 kg (RSE 12%)
    lvp <- log(56300); label("Apparent peripheral volume Vp/F for a 70 kg patient (L)")  # Table 4 tacrolimus VP = 56300 L/70 kg (RSE 8%)
    lka <- log(3.08); label("Absorption rate constant (1/h)")  # Table 4 tacrolimus KA = 3.08 1/h (RSE 39%)
    ltlag <- log(0.295); label("Absorption lag time (h)")  # Table 4 tacrolimus tLAG = 0.295 h (RSE 22%)

    # Allometric exponents, fixed a priori (Methods 'Pharmacokinetic models':
    # 'Fixed allometric exponents of 0.75 and 1 were applied to flow parameters
    # and distribution volumes, respectively').
    e_wt_cl <- fixed(0.75); label("Allometric exponent on CL/F and Q/F (unitless)")  # Methods, Pharmacokinetic models
    e_wt_vc <- fixed(1); label("Allometric exponent on Vc/F and Vp/F (unitless)")  # Methods, Pharmacokinetic models

    # Between-subject variability. Table 4 reports BSV as sqrt(omega) x 100:
    # CL 57.6% -> 0.576^2 = 0.3318; Q 68.9% -> 0.689^2 = 0.4747;
    # Vc 55.6% -> 0.556^2 = 0.3091. The control stream codes these three as an
    # $OMEGA BLOCK(3), but the final covariances are not reported, so the
    # off-diagonals are zero here (see the vignette).
    etalcl + etalq + etalvc ~ c(0.3318, 0, 0.4747, 0, 0, 0.3091)  # Table 4 tacrolimus BSV CL 57.6%, Q 68.9%, VC 55.6% (sqrt(omega) x 100)
    # Fixed small random effects kept as coded in the control stream
    # ('$OMEGA 0.01 FIX' on KA, V3 and tLAG; Table 4 lists them as not estimated).
    etalka ~ fixed(0.01)  # control stream $OMEGA 0.01 FIX ; KA
    etalvp ~ fixed(0.01)  # control stream $OMEGA 0.01 FIX ; V3
    etaltlag ~ fixed(0.01)  # control stream $OMEGA 0.01 FIX ; tLAG

    # Residual error: additive on log-transformed concentrations (control stream
    # Y = LOG(F) + EPS(1)), i.e. log-normal error in linear space. Table 4 prints
    # it as 'Proportional 36.6%' = sqrt(sigma) x 100.
    expSd <- 0.366; label("Residual error SD on the log-concentration scale")  # Table 4 tacrolimus residual error 36.6% (RSE 12%)
  })

  model({
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    vp <- exp(lvp + etalvp) * (WT / 70)^e_wt_vc
    q <- exp(lq + etalq) * (WT / 70)^e_wt_cl
    ka <- exp(lka + etalka)
    tlag <- exp(ltlag + etaltlag)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    alag(depot) <- tlag

    # Dose in mg, volume in L -> mg/L; x 1000 gives ng/mL (control stream
    # S2 = V2/1000).
    Cc <- central / vc * 1000
    Cc ~ lnorm(expSd)
  })
}
