Chauzy_2022_ceftaroline <- function() {
  description <- "Joint population PK model for the prodrug ceftaroline fosamil, its active moiety ceftaroline and the inactive open-ring metabolite ceftaroline M-1 in 18 mechanically ventilated ICU adults with early-onset pneumonia and augmented renal clearance (measured urinary creatinine clearance 83-309 mL/min) receiving ceftaroline fosamil 600 mg every 12 h as a 1 h IV infusion. One-compartment ceftaroline fosamil with complete conversion to ceftaroline; two-compartment ceftaroline whose clearance forms M-1 (fraction fm, unidentified); two-compartment M-1 with apparent (/fm) parameters. Molar-mass ratios convert the amounts between the three analytes. Power effects of creatinine clearance (centred at 180 mL/min) on ceftaroline clearance and M-1 apparent clearance; full-block exponential IIV on ceftaroline clearance, ceftaroline peripheral volume and M-1 apparent clearance; additive error for ceftaroline fosamil and proportional errors for ceftaroline and M-1. All concentrations are total plasma concentrations."
  reference <- "Chauzy A, Gregoire N, Ferrandiere M, Lasocki S, Ashenoune K, Seguin P, Boisson M, Couet W, Marchand S, Mimoz O, Dahyot-Fizelier C. Population pharmacokinetic/pharmacodynamic study suggests continuous infusion of ceftaroline daily dose in ventilated critical care patients with early-onset pneumonia and augmented renal clearance. J Antimicrob Chemother. 2022;77(11):3173-3179. doi:10.1093/jac/dkac299. Model structure is Supplementary Figure S2 and the Supplementary 'Population pharmacokinetic analysis' methods (molecular weights, complete fosamil-to-ceftaroline conversion, apparent /fm M-1 parameters); final parameter estimates and simulated secondary parameters are Supplementary Table S2; the IIV and residual variance-covariance matrix is Supplementary Table S3; per-patient covariates are Supplementary Table S1; the CLCR-clearance power relationship is main-text Equation 2 and Figure 2."
  vignette <- "Chauzy_2022_ceftaroline"
  units <- list(
    time = "h",
    dosing = "mg (ceftaroline fosamil, the administered prodrug)",
    concentration = "mg/L (total plasma ceftaroline fosamil for Cc, ceftaroline for Cc_ceftaroline, ceftaroline M-1 for Cc_m1)"
  )

  # What each ODE state holds. The M-1 states hold amount/fm because the
  # fraction of ceftaroline clearance forming M-1 was not identifiable and the
  # M-1 disposition parameters were estimated as apparent (/fm) values.
  compartmentData <- list(
    central = list(analyte = "ceftaroline fosamil", units = "mg", specimen = "plasma", verified = TRUE),
    central_ceftaroline = list(analyte = "ceftaroline", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_ceftaroline = list(analyte = "ceftaroline", units = "mg", specimen = "tissue", verified = TRUE),
    central_m1 = list(
      analyte = "ceftaroline M-1",
      units = "mg (apparent, i.e. amount/fm)",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1_m1 = list(
      analyte = "ceftaroline M-1",
      units = "mg (apparent, i.e. amount/fm)",
      specimen = "tissue",
      verified = TRUE
    )
  )

  covariateData <- list(
    CRCL = list(
      description = "Measured urinary creatinine clearance (24 h urine collection), NOT BSA-normalised",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Computed from plasma and urine creatinine concentrations and urine flow over 24 h on each PK sampling day",
        "(Study population; Table S1 footnote 1). Measured on two occasions (PK1 = first dose, PK2 = fifth to ninth dose)",
        "and used as a time-varying covariate: CLCR varied by about +/- 20% between occasions (Figure S1).",
        "Enters as power terms (CRCL/180)^0.328 on ceftaroline clearance and (CRCL/180)^0.419 on M-1 apparent",
        "clearance; 180 mL/min is the cohort median. Observed range 83-309 mL/min (all patients had augmented or",
        "normal renal function; inclusion required MDRD eGFR > 80 mL/min/1.73 m^2)."
      ),
      source_name = "CLCR"
    )
  )

  # Screened in the stepwise covariate model (Supplementary 'Covariate
  # analysis') but not retained in the final model; no estimates reported.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened (SCM, forward p < 0.01 / backward p < 0.001); not retained."
    ),
    WT = list(description = "Body weight", units = "kg", type = "continuous", notes = "Screened (SCM); not retained."),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      notes = "Screened (SCM); not retained."
    ),
    TPRO = list(
      description = "Total serum protein",
      units = "g/L",
      type = "continuous",
      notes = "Screened (SCM); not retained."
    ),
    SAPS_II = list(
      description = "Simplified Acute Physiology Score II",
      units = "points",
      type = "continuous",
      notes = "Screened (SCM); not retained."
    ),
    SOFA = list(
      description = "Sequential Organ Failure Assessment score",
      units = "points",
      type = "continuous",
      notes = "Screened (SCM); not retained."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 18L,
    n_studies = 1L,
    age_range = "21-77 years (mean 46, SD 17)",
    weight_range = "53-95.6 kg (mean 73, SD 12)",
    height_range = "1.50-1.85 m (mean 1.71)",
    bmi_range = "19-38.2 kg/m^2 (mean 25.2)",
    sex_female_pct = 27.8,
    race_ethnicity = c(Caucasian = 100),
    disease_state = "Mechanically ventilated ICU adults with early-onset (within 7 days of hospital admission) pneumonia and creatinine clearance above 80 mL/min/1.73 m^2 (augmented renal clearance in most patients); SAPS II 20-60 (mean 45)",
    renal_function = "Measured urinary CLCR 83-267 mL/min on PK1 (mean 181.5) and 100-309 mL/min on PK2 (mean 180.6); median 180 mL/min",
    dose_range = "Ceftaroline fosamil 600 mg every 12 h as a 1 h IV infusion for at least 3 days",
    regions = "France (five university-hospital ICUs: Poitiers, Tours, Angers, Nantes, Rennes)",
    notes = paste(
      "Prospective open-label PK study, NCT03025841, February 2017 - May 2018 (Study design; demographics in Table 1 and per-patient Table S1).",
      "Seven total plasma samples (pre-dose, 1, 2, 4, 6, 9 and 12 h) on PK1 (first dose) and PK2 (fifth to ninth dose); 15 of 18 patients had a PK2 occasion.",
      "Ceftaroline fosamil, ceftaroline and ceftaroline M-1 were assayed by LC-MS/MS (LOQ 0.050 mg/L each) and fit simultaneously in NONMEM 7.4 with FOCE-I; BLQ data handled by the M3 method.",
      "PK/PD simulations used an unbound fraction of 0.8 (20% protein binding, literature value) applied to total ceftaroline concentrations; fu is not part of model()."
    )
  )

  ini({
    # Ceftaroline fosamil (prodrug), one compartment. Doses are mg of
    # ceftaroline fosamil infused into central.
    lcl <- log(668);  label("Ceftaroline fosamil clearance, entirely conversion to ceftaroline (L/h)")  # Table S2 CLfosamil = 668 L/h (RSE 14%)
    lvc <- log(44.9); label("Ceftaroline fosamil volume of distribution (V1, L)")                     # Table S2 V1 = 44.9 L (RSE 29%)

    # Ceftaroline, two compartments.
    lcl_ceftaroline <- log(10.6); label("Ceftaroline clearance at CLCR = 180 mL/min (L/h)")          # Table S2 CLceftaroline,pop = 10.6 L/h (RSE 4.4%)
    e_crcl_cl_ceftaroline <- 0.328; label("Exponent of (CRCL/180) on ceftaroline clearance (unitless)") # Table S2 CLCR,cov1 = 0.328 (RSE 47%)
    lvc_ceftaroline <- log(13.2); label("Ceftaroline central volume (V2, L)")                         # Table S2 V2 = 13.2 L (RSE 10%)
    lq_ceftaroline <- log(6.79);  label("Ceftaroline intercompartmental clearance (Q1, L/h)")         # Table S2 Q1 = 6.79 L/h (RSE 14%)
    lvp_ceftaroline <- log(12);   label("Ceftaroline peripheral volume (V3, L)")                      # Table S2 V3 = 12 L (RSE 9.7%)

    # Ceftaroline M-1, two compartments, apparent parameters (divided by fm).
    lcl_m1 <- log(65.1); label("Ceftaroline M-1 apparent clearance CL/fm at CLCR = 180 mL/min (L/h)") # Table S2 CLM-1,pop = 65.1 L/h (RSE 10%)
    e_crcl_cl_m1 <- 0.419; label("Exponent of (CRCL/180) on M-1 apparent clearance (unitless)")       # Table S2 CLCR,cov2 = 0.419 (RSE 59%)
    lvc_m1 <- log(18);   label("Ceftaroline M-1 apparent central volume V4/fm (L)")                   # Table S2 V4/fm = 18 L (RSE 56%)
    lq_m1 <- log(299);   label("Ceftaroline M-1 apparent intercompartmental clearance Q2/fm (L/h)")   # Table S2 Q2/fm = 299 L/h (RSE 12%)
    lvp_m1 <- log(211);  label("Ceftaroline M-1 apparent peripheral volume V5/fm (L)")                # Table S2 V5/fm = 211 L (RSE 6.2%)

    # Full-block IIV, variances and covariances from Table S3. The diagonal
    # matches the Table S2 %CV column as sqrt(variance): 0.150, 0.321, 0.309.
    etalcl_ceftaroline + etalvp_ceftaroline + etalcl_m1 ~ c(
      0.0226,
      -0.0274, 0.103,
      0.0294, -0.0359, 0.0953
    ) # Table S3 omega block (CLceftaroline, V3, CLM-1/fm)

    # Residual error. Table S3 SIGMA for ceftaroline and M-1 proportional errors
    # (0.052 and 0.0493) are variances; their square roots are the Table S2
    # values 22.8% and 22.2%. The two were correlated (covariance 0.0347,
    # correlation 0.685, NONMEM L2 item); cross-endpoint residual correlation
    # cannot be expressed in nlmixr2, so the errors are independent here.
    addSd <- 0.271;              label("Ceftaroline fosamil additive residual error SD (mg/L)")   # Table S2 additive error = 0.271 mg/L (RSE 22%)
    propSd_ceftaroline <- 0.228; label("Ceftaroline proportional residual error SD (fraction)")   # Table S2 proportional error = 22.8% (RSE 14%)
    propSd_m1 <- 0.222;          label("Ceftaroline M-1 proportional residual error SD (fraction)") # Table S2 proportional error = 22.2% (RSE 15%)
  })

  model({
    # Molecular weights (Supplementary 'Population pharmacokinetic analysis'):
    # ceftaroline fosamil 684.7, ceftaroline 604.7, ceftaroline M-1 622.7 g/mol.
    mw_fosamil <- 684.7
    mw_ceftaroline <- 604.7
    mw_m1 <- 622.7

    # Covariates centred at the cohort median CLCR of 180 mL/min (Equation 2).
    ref_crcl <- 180

    cl <- exp(lcl)
    vc <- exp(lvc)

    cl_ceftaroline <- exp(lcl_ceftaroline + etalcl_ceftaroline) * (CRCL / ref_crcl)^e_crcl_cl_ceftaroline
    vc_ceftaroline <- exp(lvc_ceftaroline)
    q_ceftaroline <- exp(lq_ceftaroline)
    vp_ceftaroline <- exp(lvp_ceftaroline + etalvp_ceftaroline)

    cl_m1 <- exp(lcl_m1 + etalcl_m1) * (CRCL / ref_crcl)^e_crcl_cl_m1
    vc_m1 <- exp(lvc_m1)
    q_m1 <- exp(lq_m1)
    vp_m1 <- exp(lvp_m1)

    kel <- cl / vc
    kel_ceftaroline <- cl_ceftaroline / vc_ceftaroline
    k12_ceftaroline <- q_ceftaroline / vc_ceftaroline
    k21_ceftaroline <- q_ceftaroline / vp_ceftaroline
    kel_m1 <- cl_m1 / vc_m1
    k12_m1 <- q_m1 / vc_m1
    k21_m1 <- q_m1 / vp_m1

    # Ceftaroline fosamil is converted completely to ceftaroline (fraction 1),
    # so its whole clearance flux is the ceftaroline formation flux, scaled by
    # the molar-mass ratio. The fraction fm of ceftaroline clearance forming
    # M-1 is unknown; the M-1 states hold amount/fm, so the full ceftaroline
    # clearance flux (scaled by the molar-mass ratio) enters central_m1 and the
    # M-1 disposition uses the apparent /fm parameters (Figure S2).
    d/dt(central) <- -kel * central
    d/dt(central_ceftaroline) <- kel * central * mw_ceftaroline / mw_fosamil -
      kel_ceftaroline * central_ceftaroline -
      k12_ceftaroline * central_ceftaroline + k21_ceftaroline * peripheral1_ceftaroline
    d/dt(peripheral1_ceftaroline) <- k12_ceftaroline * central_ceftaroline - k21_ceftaroline * peripheral1_ceftaroline
    d/dt(central_m1) <- kel_ceftaroline * central_ceftaroline * mw_m1 / mw_ceftaroline -
      kel_m1 * central_m1 -
      k12_m1 * central_m1 + k21_m1 * peripheral1_m1
    d/dt(peripheral1_m1) <- k12_m1 * central_m1 - k21_m1 * peripheral1_m1

    # Total plasma concentrations (mg/L). Unbound ceftaroline for fT>MIC is
    # 0.8 * Cc_ceftaroline (assumed fu, Methods 'PTA and CFR').
    Cc <- central / vc
    Cc_ceftaroline <- central_ceftaroline / vc_ceftaroline
    Cc_m1 <- central_m1 / vc_m1

    Cc ~ add(addSd)
    Cc_ceftaroline ~ prop(propSd_ceftaroline)
    Cc_m1 ~ prop(propSd_m1)
  })
}
