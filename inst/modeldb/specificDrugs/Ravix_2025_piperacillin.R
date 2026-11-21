Ravix_2025_piperacillin <- function() {
  description <- "Two-compartment IV population PK model for piperacillin in hospitalised children aged 1 month to 18 years, including patients with augmented renal clearance (Ravix 2025). All four disposition parameters scale allometrically with body weight (fixed exponents 0.75 on CL and Q, 1 on V1 and V2, reference 27.6 kg); clearance is additionally a sigmoid (Hill) function of bedside-Schwartz creatinine clearance that saturates at 9.3 L/h (EC50 62 mL/min/1.73 m^2, Hill 3.5). Exponential inter-individual variability on CL only; proportional residual error with separate magnitudes for samples drawn under intermittent and continuous infusion."
  reference <- "Ravix A, Di Paolo ER, Dao K, Perez MH, Crisinel PA, Buclin T, Andre P, Marzolini C, Csajka C, Guidi M. Piperacillin dosing in pediatric patients and augmented renal clearance: are current guidelines sufficient to achieve therapeutic concentration coverage? J Antimicrob Chemother. 2025;80(12):3419-3430. doi:10.1093/jac/dkaf387"
  vignette <- "Ravix_2025_piperacillin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying (recorded at the time of each sample; Table 2 footnote b). Allometric reference 27.6 kg, the cohort median (Table 3 footnote 'BW ... with a median of 27.6 kg'; File S1 'MBW = 27.6'). Cohort range 5-73 kg (Table 2). Missing values were linearly interpolated when a measurement existed within 10 days (Methods, Data).",
      source_name = "BW"
    ),
    CRCL = list(
      description = "Creatinine clearance estimated by the bedside Schwartz equation, BSA-normalised",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Bedside Schwartz equation CrCL = 36.5 x height (cm) / serum creatinine (umol/L) (Table 2 footnote; Methods, Data). Time-varying. Enters CL through the uncentred Hill function CRCL^hill / (crcl50^hill + CRCL^hill) (File S1 '$PK': FCLCR = CLCR**HILL/(TM**HILL + CLCR**HILL)), so CL approaches its maximum 9.3 L/h (at 27.6 kg) as CRCL rises. Cohort median 146, range 23-540 mL/min/1.73 m^2; 22% of samples at >= 200 and 38% at 130-199 mL/min/1.73 m^2 (Table 2). Only 3 samples (1%) were below 30 mL/min/1.73 m^2, so the source did not simulate that range.",
      source_name = "CLCR"
    ),
    REGIMEN_CONTINF = list(
      description = "Sample drawn under continuous infusion (1) versus intermittent infusion (0)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Selects the residual-error magnitude only (File S1 '$ERROR': IF (CONTINUE.EQ.0) Y = IPRED + W*EPS(1); IF (CONTINUE.EQ.1) Y = IPRED + W*EPS(2)); no structural parameter depends on it. 0 = intermittent infusion (median 1 h, range 1 min to 4 h); 1 = continuous infusion (median 24 h, range 2-28 h) (Table 2). 334 of the 386 concentrations were drawn under intermittent and 52 under continuous infusion (Results). Supply on the observation rows; the value on dose rows has no effect.",
      source_name = "CONTINUE"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      notes = "Screened in the univariate step (Methods, Model building); not significant (Results: 'no significant effect of sex, age, or urea (dOFV > -5, P > 0.03)'). Table S1 lists the tested model."
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened in the univariate step; not significant (Results). The source notes that a maturation function could not be estimated because only 3% of patients were younger than 1 year (Discussion)."
    ),
    BUN = list(
      description = "Blood urea",
      units = "mmol/L",
      type = "continuous",
      notes = "Screened in the univariate step ('urea'); not significant (Results). Cohort median 3.1, range 0.4-40 mmol/L (Table 2)."
    )
  )

  compartmentData <- list(
    central = list(analyte = "piperacillin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "piperacillin", units = "mg", specimen = "not applicable", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 104L,
    n_studies = 1L,
    n_concentrations = 386L,
    age_range = "1 month to 18 years at inclusion; 0.3-18.1 years at the time of sampling",
    age_median = "9.2 years",
    weight_range = "5-73 kg",
    weight_median = "27.6 kg",
    sex_female_pct = 38.5,
    disease_state = "Hospitalised children receiving intravenous piperacillin/tazobactam with routine therapeutic drug monitoring: oncology (36% of concentrations), intensive care (18%), intermediate care (15%), surgery (11%) and general paediatrics (20%). 72% of infection events had no bacteria identified. Patients on extracorporeal membrane oxygenation or renal replacement therapy were excluded.",
    dose_range = "Intermittent infusion median 97 (8-180) mg/kg per dose, median 3 doses/day, infusion median 1 h (1 min to 4 h); continuous infusion median 321 (64-540) mg/kg/day over a median 24 h (Table 2).",
    regions = "Switzerland (Lausanne University Hospital)",
    renal_function = "Bedside-Schwartz CrCL median 146 (23-540) mL/min/1.73 m^2; 22% of samples >= 200, 38% 130-199, 24% 90-129, 10% 60-89, 5% 30-59 and 1% < 30 mL/min/1.73 m^2 (Table 2).",
    notes = "Demographics from Ravix 2025 Table 2 (median (min, max) at the time of sampling). Retrospective TDM data collected 2012-2022; LC-MS/MS assay, LLOQ 0.08 mg/L. Table 2 prints 'Female/Male 40 (31%)/64 (62%)'; 40 of 104 is 38.5%, which is the value recorded here. NONMEM 7.4.3 FOCE-I (File S1). Bootstrap with 2000 replicates; pcVPC Figure 1."
  )

  ini({
    # Clearance (Table 3 'Final equations' and File S1 $PK):
    #   CL = thetaCL * (BW/27.6)^0.75 * CrCL^Hill / (EC50^Hill + CrCL^Hill) * exp(eta)
    # thetaCL is therefore the asymptotic maximum at saturating CrCL, for a
    # 27.6 kg child.
    lclmax <- log(9.30)
    label("Maximum (asymptotic) clearance at saturating CrCL for a 27.6 kg child (L/h)") # Table 3, 'thetaCL (L/h)' = 9.30 (RSE 8%)
    lcrcl50 <- log(62.00)
    label("Creatinine clearance giving half-maximal clearance (mL/min/1.73 m^2)") # Table 3, 'EC50 (mL/min/1.73 m2)' = 62.00 (RSE 5%)
    lhill <- log(3.50)
    label("Hill coefficient of the clearance / creatinine clearance relationship (unitless)") # Table 3, 'Hill' = 3.50 (RSE 33%)

    lvc <- log(8.40)
    label("Central volume of distribution for a 27.6 kg child (L)") # Table 3, 'thetaV1 (L)' = 8.40 (RSE 12%)
    lq <- log(0.63)
    label("Intercompartmental clearance for a 27.6 kg child (L/h)") # Table 3, 'thetaQ (L/h)' = 0.63 (RSE 23%)
    lvp <- log(2.00)
    label("Peripheral volume of distribution for a 27.6 kg child (L)") # Table 3, 'thetaV2 (L)' = 2.00 (RSE 21%)

    # Allometric exponents, all fixed (Methods 'allometric exponent fixed at
    # 0.75 for clearances and 1 for volumes'; File S1 $THETA 5-8 FIX).
    e_wt_cl <- fixed(0.75)
    label("Allometric exponent of body weight on CL (unitless)") # File S1 $THETA '0.75 FIX ; BW sur CL'
    e_wt_vc <- fixed(1)
    label("Allometric exponent of body weight on V1 (unitless)") # File S1 $THETA '1 FIX ; BW sur V1'
    e_wt_q <- fixed(0.75)
    label("Allometric exponent of body weight on Q (unitless)") # File S1 $THETA '0.75 FIX ; BW sur Q'
    e_wt_vp <- fixed(1)
    label("Allometric exponent of body weight on V2 (unitless)") # File S1 $THETA '1 FIX ; BW sur V2'

    # IIV on CL only (IIV on V1 was significant but dropped for stability).
    # Table 3 prints 31% CV. The CV is taken as sqrt(omega^2): the File S1
    # $OMEGA value 0.0955 gives sqrt = 30.9%, whereas sqrt(exp(0.0955) - 1)
    # would print as 32%, and both $SIGMA values (0.393, 0.19) convert to the
    # printed 63% and 44% by the same square root.
    etalclmax ~ 0.31^2 # Table 3, 'IIVCL (CV%)' = 31.00 (RSE 11%, shrinkage 8%); CV = sqrt(omega^2)

    # Proportional residual error, one magnitude per infusion type (File S1
    # $ERROR: Y = IPRED + IPRED*EPS(1) when CONTINUE = 0, EPS(2) when 1).
    propSd_intermittent <- 0.63
    label("Proportional residual error, samples under intermittent infusion (fraction)") # Table 3, 'sigmaIntermittent (CV%)' = 63.00 (RSE 5%)
    propSd_continuous <- 0.44
    label("Proportional residual error, samples under continuous infusion (fraction)") # Table 3, 'sigmacontinuous (CV%)' = 44.00 (RSE 15%)
  })

  model({
    # Covariate terms (Table 3 'Final equations'; File S1 $PK).
    wt_ratio <- WT / 27.6
    crcl50 <- exp(lcrcl50)
    hill <- exp(lhill)
    f_crcl <- CRCL^hill / (crcl50^hill + CRCL^hill)

    # Individual parameters.
    cl <- exp(lclmax + etalclmax) * wt_ratio^e_wt_cl * f_crcl
    vc <- exp(lvc) * wt_ratio^e_wt_vc
    q <- exp(lq) * wt_ratio^e_wt_q
    vp <- exp(lvp) * wt_ratio^e_wt_vp

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Dose in mg, vc in L -> mg/L.
    Cc <- central / vc

    # Residual magnitude selected by the infusion type of the sample.
    propSd_inf <- propSd_intermittent * (1 - REGIMEN_CONTINF) +
      propSd_continuous * REGIMEN_CONTINF
    Cc ~ prop(propSd_inf)
  })
}
