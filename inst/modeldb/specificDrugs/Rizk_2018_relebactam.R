Rizk_2018_relebactam <- function() {
  description <- paste(
    "Two-compartment intravenous population PK model for relebactam in healthy adult volunteers",
    "(Rizk 2018), fitted by a naive pooled approach to total plasma concentrations and to",
    "urea-corrected epithelial lining fluid (ELF) concentrations from bronchoalveolar lavage after",
    "five doses of relebactam 250 mg with imipenem-cilastatin 500 mg every 6 h. ELF is not a",
    "kinetic compartment: the ELF concentration is a time-invariant partition coefficient (0.553)",
    "times the predicted unbound plasma concentration, with the unbound fraction in plasma fixed",
    "at 0.8. No covariates and no between-subject variability (naive pooled fit).",
    sep = " "
  )
  reference <- paste(
    "Rizk ML, Rhee EG, Jumes PA, Gotfried MH, Zhao T, Mangin E, Bi S, Chavez-Eng CM, Zhang Z,",
    "Butterton JR (2018). Intrapulmonary pharmacokinetics of relebactam, a novel beta-lactamase",
    "inhibitor, dosed in combination with imipenem-cilastatin in healthy subjects.",
    "Antimicrob Agents Chemother. 2018;62(3):e01411-17. doi:10.1128/AAC.01411-17",
    sep = " "
  )
  vignette <- "Rizk_2018_relebactam"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    # ELF is not an ODE state: Rizk 2018 Eq. 3 is an algebraic partition of the
    # unbound central concentration, so only the two disposition states exist.
    central = list(analyte = "relebactam", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "relebactam", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species = "human",
    n_subjects = 16,
    n_studies = 1,
    age_range = "24-42 years",
    sex_female_pct = 17.6,
    disease_state = "healthy volunteers",
    dose_range = paste(
      "Relebactam 250 mg with imipenem-cilastatin 500 mg, IV infusion over 30 min every 6 h for",
      "5 doses"
    ),
    regions = "United States (single site, Phoenix AZ; MK-7655 protocol 007)",
    sampling = paste(
      "Plasma before the first and fifth doses and at bronchoscopy; one bronchoscopy/BAL per subject",
      "at 0.5, 1, 1.5 or 3 h after the fifth dose (4 subjects per time point)."
    ),
    notes = paste(
      "17 subjects entered (14 male, 3 female); 16 completed and form the per-protocol PK population",
      "(one discontinued during dose 2 for an adverse event). Body weight was not reported;",
      "eligibility required BMI <= 32 kg/m^2 and creatinine clearance >= 80 mL/min.",
      "Rizk 2018 Results (first paragraph) and Materials and Methods."
    )
  )

  ini({
    # Structural parameters (Table 4). Naive pooled fit, so these are the only
    # source of variability-free typical values; no IIV was estimated.
    lcl <- log(9.17); label("Clearance (L/h)")                                  # Table 4: thetaCL = 9.17 L/h (RSE 5.04%)
    lvc <- log(15.3); label("Central volume of distribution (L)")               # Table 4: thetaVc = 15.3 L (RSE 15.9%)
    lvp <- log(10.6); label("Peripheral volume of distribution (L)")            # Table 4: thetaVp = 10.6 L (RSE 384.9%)
    lq <- log(2.64); label("Intercompartmental clearance (L/h)")                # Table 4: thetaQ2 = 2.64 L/h (RSE 107.2%)

    # ELF partition coefficient on the unbound plasma concentration (Eq. 3)
    lfelf <- log(0.553); label("ELF:unbound-plasma partition coefficient (unitless)") # Table 4: thetakELF = 0.553 (RSE 9.06%)

    # Unbound fraction in plasma, fixed from protein-binding data (Methods: fu,REL ~80%)
    fu <- fixed(0.8); label("Fraction of relebactam unbound in plasma (unitless)")    # Methods after Eq. 3: 'free fraction of relebactam in plasma (~80%)'; Results: 'relebactam is 80% unbound in plasma'

    # Residual error - plasma (combined) and ELF (proportional)
    addSd <- fixed(0.01); label("Additive residual SD for plasma Cc (mg/L)")           # Table 4: sigma addi = 0.01 mg/liter (fixed)
    propSd <- 0.28; label("Proportional residual SD for plasma Cc (fraction)")         # Table 4: sigma prop = 28.0 percent CV (RSE 49.5%)
    propSd_Celf <- 0.392; label("Proportional residual SD for ELF Celf (fraction)")    # Table 4: sigma prop,kELF = 39.2 percent CV (RSE 44.2%)
  })

  model({
    cl <- exp(lcl)
    vc <- exp(lvc)
    vp <- exp(lvp)
    q <- exp(lq)
    felf <- exp(lfelf)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Rizk 2018 Eq. 1-2 (A1 = central, A2 = peripheral1)
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Total plasma concentration (mg/L)
    Cc <- central / vc

    # Unbound plasma concentration (mg/L); the ELF free fraction is assumed 100%
    Cu <- fu * Cc

    # Rizk 2018 Eq. 3: C_ELF = kELF * fu,REL * A1 / V1
    Celf <- felf * Cu

    Cc ~ add(addSd) + prop(propSd)
    Celf ~ prop(propSd_Celf)
  })
}
