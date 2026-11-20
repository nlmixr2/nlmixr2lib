Liu_2026_meropenem <- function() {
  description <- paste(
    "One-compartment intravenous population PK model for continuous-infusion meropenem",
    "linked to an indirect-response (turnover) PD model for plasma C-reactive protein",
    "(C-RP) in 141 adult onco-haematological patients with febrile neutropenia (Liu 2026).",
    "Meropenem clearance rises exponentially with CKD-EPI estimated glomerular filtration",
    "rate; the effect is UNCENTERED, so exp(lcl) is the clearance extrapolated to eGFR = 0.",
    "The volume of distribution is fixed to 28.5 L from Contejean 2020, with an estimated",
    "random effect. Total plasma meropenem fully inhibits C-RP production (no Imax term,",
    "IC50 41.5 mg/L), and the production rate kin itself decays exponentially with time",
    "since the start of therapy, kin * exp(-k * t), which lets C-RP first rise and then",
    "fall. The C-RP baseline R0 is estimated separately from kin and kout, so the C-RP",
    "state does NOT start at the kin/kout steady state.",
    sep = " "
  )
  reference <- paste(
    "Liu C, Cojutti PG, van Hasselt JGC, Bartoletti M, Giannella M, Zinzani P, Viale P,",
    "Pea F. Relationship between continuous infusion meropenem PK/PD target attainment",
    "and C-reactive protein dynamics in onco-haematologic patients with febrile",
    "neutropenia. Br J Clin Pharmacol. 2026;92(10):3646-3655.",
    "doi:10.1002/bcp.70667",
    sep = " "
  )
  vignette <- "Liu_2026_meropenem"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "meropenem", units = "mg", specimen = "plasma", verified = TRUE),
    # The C-RP state holds a CONCENTRATION (mg/dL): the paper's equation is
    # written directly on R, 'C-RP concentration in plasma', with no volume.
    crp = list(analyte = "C-reactive protein", units = "mg/dL", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = "CKD-EPI estimated glomerular filtration rate, BSA-normalised",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Source name eGFR. Methods: 'Estimated glomerular filtration rate was calculated by",
        "means of the Chronic Kidney Disease Epidemiology (CKD-EPI) formula', expressed in",
        "mL/min/1.73 m2; serum creatinine was collected at each TDM assessment, so the column",
        "may be time-varying. Table 1 median 97.0 (IQR 80.0-109.0), Results range 17.0-195.6.",
        "The effect is applied UNCENTERED and EXPONENTIALLY, cl = exp(lcl) * exp(e_crcl_cl *",
        "CRCL), which is the Monolix parameterisation for an untransformed continuous",
        "covariate and the form this group printed verbatim for meropenem in critically ill",
        "patients (modellib('Troisi_2024_meropenem')) and for dalbavancin",
        "(modellib('Cojutti_2024_dalbavancin')). Only this form reproduces the magnitude of",
        "the reported median clearance (13.34 L/h): at the cohort median eGFR it gives",
        "4.99 * exp(0.009 * 97) = 11.9 L/h, whereas a centred form would give 4.99 L/h."
      ),
      source_name = "eGFR"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Tested as a potential PK/PD covariate (Methods) but not retained. Table 1 median 70 kg (IQR 62-80); Results range 42-118 kg."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Tested as a potential PK/PD covariate (Methods) but not retained."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      notes = "Gender tested as a potential PK/PD covariate (Methods) but not retained. Table 1: 57 female / 84 male (40.4% / 59.6%)."
    ),
    CONMED_ANTIMRSA = list(
      description = "Concomitant anti-Gram-positive antibacterial indicator",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Concomitant treatment with anti-Gram-positive agents was tested as a potential",
        "PK/PD covariate (Methods) but not retained; 81/141 (57.4%) received one. In the",
        "sibling critically-ill analysis (Troisi 2024) the same indicator WAS retained on",
        "IC50; here it was not."
      )
    )
  )

  population <- list(
    species = "human",
    n_subjects = 141L,
    n_studies = 1L,
    n_observations = "283 meropenem steady-state plasma concentrations and 1373 plasma C-RP concentrations",
    age_range = "19-81 years",
    age_median = "59 years (IQR 52-64)",
    weight_range = "42-118 kg",
    weight_median = "70 kg (IQR 62-80)",
    sex_female_pct = 40.4,
    race_ethnicity = "Not reported; single-centre Italian cohort.",
    disease_state = paste(
      "Adult onco-haematological patients with febrile neutropenia treated with",
      "TDM-guided continuous-infusion meropenem. Underlying disease (Table 1): acute",
      "myeloid leukaemia 91 (64.5%), non-Hodgkin lymphoma 23 (16.3%), acute lymphoid",
      "leukaemia 9 (6.4%), multiple myeloma 6 (4.3%), chronic myeloid leukaemia 3,",
      "chronic lymphoid leukaemia 2, other 7. Microbiological isolates in 46/141 (32.6%).",
      "Co-treatment with antifungals in 134/141 (95.0%) and with an anti-Gram-positive",
      "agent in 81/141 (57.4%)."
    ),
    renal_function = "CKD-EPI eGFR median 97.0 mL/min/1.73 m^2 (IQR 80.0-109.0; range 17.0-195.6).",
    baseline_crp = "Median 9.8 mg/dL (IQR 4.9-17.4); median 12 (range 1-32) subsequent C-RP assessments per patient.",
    dose_range = paste(
      "Loading dose 2 g over 30 min, then 1 g q8h infused over 8 h (eGFR >= 60",
      "mL/min/1.73 m2) or 0.5 g q6h infused over 6 h (eGFR < 60), i.e. continuous infusion,",
      "TDM-adjusted to a steady-state concentration of 8-16 mg/L. Table 1 median",
      "maintenance dose 1 g q6h (IQR 1 g q8h to 1 g q6h); median treatment duration 10",
      "days (range 2-32). Observed Css median 11.4 mg/L (range 3.0-39.5)."
    ),
    regions = "Italy (IRCCS Azienda Ospedaliero-Universitaria di Bologna)",
    notes = paste(
      "Retrospective single-centre study, January 2021 to July 2023. Meropenem by",
      "LC-MS/MS (LLOQ 0.3 mg/L); C-RP LLOQ 0.1 mg/dL; no observation was below either",
      "LLOQ. Fitted in Monolix 2023R with a sequential approach: the PK model first, then",
      "the PD model with the PK population values and random-effect SDs fixed. Covariates",
      "also tested but not retained: type of onco-haematological disease, concomitant",
      "antifungals and treatment duration (with weight, BMI, gender and anti-Gram-positive",
      "co-therapy; see covariatesDataExcluded)."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # PK -- Liu 2026 Table 2, 'PK/PD model Value' column. The PK block was
    # estimated in the first (PK-only) step and then held at these values
    # while the PD parameters were estimated (Methods, 'Population PK/PD
    # analysis'); they are estimates, not literature constants, so only V
    # is wrapped in fixed().
    # ------------------------------------------------------------------
    lcl <- log(4.99); label("Meropenem clearance intercept at CRCL = 0 (L/h)")                         # Table 2: CL = 4.99 L/h (RSE 15.4%; bootstrap median 4.86, 5th-95th 3.95-7.23)
    lvc <- fixed(log(28.5)); label("Meropenem volume of distribution (L)")                             # Table 2: V = 28.5 L, no RSE; Methods: 'Volume of distribution was fixed to 28.5 L, based on Contejean et al.'
    e_crcl_cl <- 0.009; label("Exponential coefficient of CRCL on CL (per mL/min/1.73 m^2; uncentered)") # Table 2: beta_eGFR = 0.009 (RSE 17.7%; bootstrap median 0.0093, 5th-95th 0.0054-0.0115)

    # ------------------------------------------------------------------
    # PD -- indirect turnover model with inhibition of C-RP production,
    # Methods equation: dR/dt = kin * (1 - Cc / (Cc + IC50)) * e^(-k*t) - kout * R.
    # R0, kin and kout are all estimated, so R0 != kin / kout.
    # ------------------------------------------------------------------
    lrbase <- log(5.82); label("Baseline plasma C-RP, R0 (mg/dL)")                                        # Table 2: R0 = 5.82 mg/dL (RSE 13.2%; bootstrap median 7.08, 5th-95th 3.85-9.28)
    lkout <- log(0.02); label("First-order C-RP degradation rate constant (1/h)")                        # Table 2: kout = 0.02 1/h (RSE 4.86%; bootstrap median 0.02, 5th-95th 0.02-0.02)
    lic50 <- log(41.5); label("Meropenem concentration for half-maximal inhibition of C-RP production (mg/L)") # Table 2: IC50 = 41.5 mg/L (RSE 23.6%; bootstrap median 79, 5th-95th 23.72-405.58)
    lkin <- log(0.49); label("Zero-order C-RP production rate constant at t = 0 (mg/dL/h)")              # Table 2: kin = 0.49 mg/dL*1/h (RSE 6.09%; bootstrap median 0.42, 5th-95th 0.30-0.56)
    lkin_exp_kdes <- log(0.0037); label("First-order rate constant of the exponential decline of kin with time, k (1/h)") # Table 2: k = 0.0037 1/h (RSE 15.7%; bootstrap median 0.0035, 5th-95th 0.0019-0.0050)

    # ------------------------------------------------------------------
    # Inter-individual variability. Table 2 'Random effects' rows are
    # Monolix omegas, the SD of the normal random effect on the log scale
    # (Methods: 'Exponential random effects ... A log-normal distribution of
    # population parameters was applied'), so each variance is omega^2. No
    # correlation rows are printed, so the matrix is diagonal.
    # ------------------------------------------------------------------
    etalcl ~ 0.2025           # 0.45^2; Table 2: omega CL = 0.45 (RSE 7.34%)
    etalvc ~ 2.6896           # 1.64^2; Table 2: omega V = 1.64 (RSE 12.7%), estimated although the typical V is held constant
    etalrbase ~ 1.5625        # 1.25^2; Table 2: omega R0 = 1.25 (RSE 8.35%)
    etalkout ~ 0.0441         # 0.21^2; Table 2: omega kout = 0.21 (RSE 22.4%)
    etalic50 ~ 2.8224         # 1.68^2; Table 2: omega IC50 = 1.68 (RSE 11.2%)
    etalkin ~ 0.2209          # 0.47^2; Table 2: omega kin = 0.47 (RSE 10.2%)
    etalkin_exp_kdes ~ 2.2801 # 1.51^2; Table 2: omega k = 1.51 (RSE 7.59%)

    # ------------------------------------------------------------------
    # Residual error. Table 2 footnote: b1 proportional error of the PK
    # model; a2 constant and b2 proportional error of the PD model.
    # ------------------------------------------------------------------
    propSd <- 0.27; label("Proportional residual error on meropenem plasma concentration (fraction)") # Table 2: b1 = 0.27 (RSE 5.99%)
    addSd_crp <- 0.11; label("Additive residual error on plasma C-RP (mg/dL)")                         # Table 2: a2 = 0.11 (RSE 24.6%)
    propSd_crp <- 0.31; label("Proportional residual error on plasma C-RP (fraction)")                 # Table 2: b2 = 0.31 (RSE 2.88%)
  })

  model({
    # Individual PK parameters; the eGFR effect is exponential and uncentered.
    cl <- exp(lcl + etalcl) * exp(e_crcl_cl * CRCL)
    vc <- exp(lvc + etalvc)
    kel <- cl / vc

    # Individual PD parameters.
    rbase <- exp(lrbase + etalrbase)
    kout <- exp(lkout + etalkout)
    ic50 <- exp(lic50 + etalic50)
    kin <- exp(lkin + etalkin)
    kin_exp_kdes <- exp(lkin_exp_kdes + etalkin_exp_kdes)

    # Total plasma meropenem (mg/L), the driver of C-RP production inhibition.
    Cc <- central / vc

    # Methods equation verbatim. t is time since the start of meropenem
    # (the first dose is at t = 0). Inhibition is FULL (no Imax term), and
    # the production rate decays as exp(-k * t), the paper's 'k ... introduced
    # to discriminate between patients with an upward or downward trend'.
    d/dt(central) <- -kel * central
    d/dt(crp) <- kin * (1 - Cc / (Cc + ic50)) * exp(-kin_exp_kdes * t) - kout * crp
    crp(0) <- rbase

    Cc ~ prop(propSd)
    # Monolix 'combined' error with its default combined1 form, sd = a + b * f.
    crp ~ add(addSd_crp) + prop(propSd_crp) + combined1()
  })
}
