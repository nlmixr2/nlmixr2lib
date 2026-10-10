Quintairos_2021_mycophenolic_acid <- function() {
  description <- paste(
    "Two-compartment population PK model for mycophenolic acid (MPA) after oral",
    "mycophenolate mofetil (doses entered as MPA molar equivalents) in adult de novo",
    "kidney transplant recipients during the first 6 months post-transplant",
    "(Quintairos 2021). First-order absorption with an absorption lag time and",
    "first-order elimination; fixed-exponent allometric scaling on body weight",
    "(0.75 on CL/F and Q/F; 1 on Vc/F and Vp/F; and, as in the published control",
    "stream, a linear weight scaling on ka). Vp/F is fixed at 800 L. Correlated",
    "between-subject variability on CL/F, Vc/F and Vp/F, small fixed-variance",
    "random effects on Q/F, ka and the lag time, and proportional residual error."
  )
  reference <- paste(
    "Quintairos L, Colom H, Millan O, Fortuna V, Espinosa C, Guirado L, Budde K,",
    "Sommerer C, Lizana A, Lopez-Pua Y, Brunet M. Early prognostic performance of",
    "miR155-5p monitoring for the risk of rejection: Logistic regression with a",
    "population pharmacokinetic approach in adult kidney transplant patients.",
    "PLoS ONE. 2021;16(1):e0245880. doi:10.1371/journal.pone.0245880.",
    "Parameter values from Table 4; model structure from the S1 Appendix NONMEM",
    "control stream ('MPA PK model').",
    sep = " "
  )
  vignette <- "Quintairos_2021_tacrolimus_mycophenolic_acid"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-varying (the deposited dataset S2 Table updates WGT at each visit).",
        "Reference weight 70 kg (control stream MCOV = 70). Exponent 0.75 on CL/F",
        "and Q/F, 1 on Vc/F and Vp/F, and 1 (linear) on ka; all exponents fixed.",
        "Cohort median 73 kg (IQR 62.9-86.8; Table 1)."
      ),
      source_name = "WGT"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "mycophenolic acid", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "mycophenolic acid", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "mycophenolic acid", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 58L,
    n_studies = 1L,
    n_observations = "1071 MPA plasma concentrations",
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
      "Oral mycophenolate mofetil (Myfenax) twice daily, mean 1547 mg per dose",
      "(IQR 1250-2000; Table 2), converted to MPA molar equivalents for modelling",
      "(1000 mg MMF = 738.96 mg MPA in the deposited dataset)."
    ),
    regions = "Germany (Charite Berlin, Heidelberg) and Spain (Fundacio Puigvert, Barcelona)",
    co_medication = "Tacrolimus (Prograf), methylprednisolone, basiliximab induction (2 x 20 mg)",
    notes = paste(
      "Pharmacokinetic sub-study (58 of 80 patients) of a European multicentre",
      "prospective observational study (EudraCT 2013-001817-33). Sampling at 0, 0.5,",
      "1, 1.5, 2, 3, 4, 6, 8 and 12 h post-dose at week 1 and at 0, 1.5, 2 and 4 h at",
      "months 1, 2, 3 and 6. Baseline demographics in Table 1."
    )
  )

  ini({
    # Structural parameters: Quintairos 2021 Table 4, MPA column (reference 70 kg).
    lcl <- log(11.8); label("Apparent clearance CL/F for a 70 kg patient (L/h)")  # Table 4 MPA CL = 11.8 L/h/70 kg (RSE 5%)
    lvc <- log(106); label("Apparent central volume Vc/F for a 70 kg patient (L)")  # Table 4 MPA VC = 106 L/70 kg (RSE 22%)
    lq <- log(37.1); label("Apparent intercompartmental clearance Q/F for a 70 kg patient (L/h)")  # Table 4 MPA Q = 37.1 L/h/70 kg (RSE 9%)
    lvp <- fixed(log(800)); label("Apparent peripheral volume Vp/F for a 70 kg patient (L)")  # Table 4 MPA VP = 800 FIX; control stream THETA(3) '800 FIX'
    lka <- log(1.79); label("Absorption rate constant for a 70 kg patient (1/h)")  # Table 4 MPA KA = 1.79 1/h (RSE 16%)
    ltlag <- log(0.243); label("Absorption lag time (h)")  # Table 4 MPA tLAG = 0.243 h (RSE 29%)

    # Allometric exponents, fixed a priori (Methods 'Pharmacokinetic models';
    # control stream ECOV = (WGT/70)**0.75, DCOV = WGT/70).
    e_wt_cl <- fixed(0.75); label("Allometric exponent on CL/F and Q/F (unitless)")  # Methods; control stream ECOV
    e_wt_vc <- fixed(1); label("Allometric exponent on Vc/F and Vp/F (unitless)")  # Methods; control stream DCOV
    # The ka scaling is coded in the control stream but not described in the
    # article text (Methods scales 'disposition' parameters; Table 4 prints ka
    # in 1/h). Kept as executed; see the vignette.
    e_wt_ka <- fixed(1); label("Weight exponent on ka, as coded in the control stream (unitless)")  # control stream KA = EXP(MU_5+ETA(5)) * DCOV

    # Between-subject variability. Table 4 reports BSV as sqrt(omega) x 100:
    # CL 34.9% -> 0.349^2 = 0.1218; Vc 133.8% -> 1.338^2 = 1.790;
    # Vp 164.6% -> 1.646^2 = 2.709. The control stream codes these three as an
    # $OMEGA BLOCK(3), but the final covariances are not reported, so the
    # off-diagonals are zero here (see the vignette).
    etalcl + etalvc + etalvp ~ c(0.1218, 0, 1.790, 0, 0, 2.709)  # Table 4 MPA BSV CL 34.9%, VC 133.8%, VP 164.6% (sqrt(omega) x 100)
    # Fixed small random effects kept as coded in the control stream
    # ('$OMEGA 0.01 FIX' on Q, KA and tLAG; Table 4 lists them as not estimated).
    etalq ~ fixed(0.01)  # control stream $OMEGA 0.01 FIX ; Q
    etalka ~ fixed(0.01)  # control stream $OMEGA 0.01 FIX ; KA
    etaltlag ~ fixed(0.01)  # control stream $OMEGA 0.01 FIX ; tLAG

    # Residual error: proportional (control stream Y = F + F*EPS(1)).
    propSd <- 0.553; label("Proportional residual error (fraction)")  # Table 4 MPA residual error proportional 55.3% (RSE 6%)
  })

  model({
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    vp <- exp(lvp + etalvp) * (WT / 70)^e_wt_vc
    q <- exp(lq + etalq) * (WT / 70)^e_wt_cl
    ka <- exp(lka + etalka) * (WT / 70)^e_wt_ka
    tlag <- exp(ltlag + etaltlag)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    alag(depot) <- tlag

    # Dose in mg of MPA, volume in L -> mg/L (control stream S2 = V2).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
