Tanigawara_2012_garenoxacin_interim <- function() {
  description <- paste(
    "One-compartment population PK model with first-order absorption for oral",
    "garenoxacin in Japanese healthy volunteers and patients with respiratory",
    "tract infection: the interim end-of-phase-II model (Tanigawara 2012",
    "Model 2, phase I-II data, NONMEM V first-order method) that drove the",
    "paper's Monte Carlo dose selection. Apparent clearance is",
    "linear-additive in Cockcroft-Gault creatinine clearance; apparent volume",
    "and absorption rate have no covariates. The final phase I-III model is",
    "Tanigawara_2012_garenoxacin."
  )
  reference <- paste(
    "Tanigawara Y, Nozawa K, Tsuda H. Optimal dose finding of garenoxacin",
    "based on population pharmacokinetics/pharmacodynamics and Monte Carlo",
    "simulation. Eur J Clin Pharmacol. 2012;68(1):39-53.",
    "doi:10.1007/s00228-011-1095-3"
  )
  vignette <- "Tanigawara_2012_garenoxacin"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  compartmentData <- list(
    depot = list(analyte = "garenoxacin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "garenoxacin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = "Creatinine clearance by the Cockcroft-Gault equation on total body weight, raw mL/min (NOT BSA-normalized)",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Tanigawara 2012 Table 2 footnote b (Cockcroft-Gault, 0.85 factor for",
        "women). Linear-additive effect on CL/F centred at 102 mL/min, the",
        "equation printed in Table 4 for Model 2 and the arithmetic mean of",
        "the 104 phase I/II subjects (Table 2: 102.4 mL/min, range",
        "34.4-178.7). Table 3 prints the same Model 2 row with the final",
        "model's centring of 85; see the vignette for why 102 is used."
      ),
      source_name = "CLCR"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 104L,
    n_studies = 4L,
    n_observations = 1396L,
    age_range = "20-74 years",
    weight_range = "34-83 kg",
    sex_female_pct = 26.0,
    race_ethnicity = c(Japanese = 100),
    disease_state = paste(
      "49 healthy volunteers (three phase I studies) and 55 patients with",
      "community-acquired pneumonia or secondary infection of chronic",
      "respiratory disease (phase II)."
    ),
    dose_range = paste(
      "Phase I 100-600 mg single dose, 200 or 400 mg once daily for 7-14 days,",
      "and a 400 mg fed/fasted study; phase II 200 or 400 mg once daily for",
      "7-14 days. All oral."
    ),
    regions = "Japan",
    renal_function = "Cockcroft-Gault CLCR 34.4-178.7 mL/min (mean 102.4); phase II patients alone 84.6 +/- 24.8 mL/min",
    notes = paste(
      "Tanigawara 2012 Tables 1-2 and Results. Men 77 of 104. Estimated in",
      "2002 with NONMEM V, first-order (FO) method. Used for the Monte Carlo",
      "target-attainment (Figure 3) and trough-above-MPC (Figure 4)",
      "simulations."
    )
  )

  ini({
    # Tanigawara 2012 Table 4, 'End of phase II stage' column (Model 2, NONMEM
    # V FO). CL/F = theta1 + theta2 x (CLCR - 102) in mL/min; V/F = theta4.
    lka <- log(4.72); label("First-order absorption rate constant (1/h)") # Table 4 Model 2: ka = 4.72 +/- 0.54 1/h
    lcl <- log(66.7); label("Apparent clearance CL/F at CLCR 102 mL/min (mL/min)") # Table 4 Model 2: theta1 = 66.7 +/- 1.7 mL/min
    lvc <- log(61.8); label("Apparent volume V/F (L)") # Table 4 Model 2: theta4 = 61.8 +/- 1.7 L

    e_crcl_cl <- 0.167; label("Additive change in CL/F per mL/min of (CRCL - 102) (mL/min per mL/min)") # Table 4 Model 2: theta2 = 0.167 +/- 0.036

    # Table 4 Model 2 omega^2 on the log scale: CL/F 0.0424 (CV 20.6%), V/F
    # 0.0443 (CV 21.0%), covariance 0.0330 (r = 0.761), ka 2.63 (CV 162%).
    etalcl + etalvc ~ c(0.0424, 0.0330, 0.0443)
    etalka ~ 2.63

    # Eq. 1 constant-CV error; Table 4 sigma^2 = 0.0645 (CV 25.4%), so the
    # proportional SD is sqrt(0.0645) = 0.254.
    propSd <- 0.254; label("Proportional residual error (SD, fraction)") # Table 4 Model 2: sigma^2 = 0.0645 +/- 0.0057, CV 25.4%
  })
  model({
    # Linear-additive CL/F (mL/min); x 60 / 1000 converts mL/min to L/h.
    cl <- (exp(lcl) + e_crcl_cl * (CRCL - 102)) * exp(etalcl) * 60 / 1000
    vc <- exp(lvc + etalvc)
    ka <- exp(lka + etalka)

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
