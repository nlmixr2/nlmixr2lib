Jullien_2005_abacavir <- function() {
  description <- "Two-compartment population PK model for abacavir in HIV-infected adults (Jullien 2005); apparent clearance scales with body weight via an estimated power exponent, Q/F is fixed when BW is added to the model"
  reference <- "Jullien V, Treluyer J-M, Chappuy H, Dimet J, Rey E, Dupin N, Salmon D, Pons G, Urien S. Weight related differences in the pharmacokinetics of abacavir in HIV-infected patients. Br J Clin Pharmacol. 2005;59(2):183-188. doi:10.1111/j.1365-2125.2004.02259.x"
  vignette <- "Jullien_2005_abacavir"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline body weight. Estimated power-form effect on apparent clearance (CL/F): CL/F = TV(CL/F) * (WT/65)^e_wt_cl. Reference is the population median (65 kg).",
      source_name        = "BW"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 188,
    n_studies      = 1,
    age_range      = "16.1-72.8 years",
    age_median     = "40 years",
    weight_range   = "36-102 kg",
    weight_median  = "65 kg",
    sex_female_pct = 19.7,
    disease_state  = "HIV-infected adults on highly active antiretroviral therapy (HAART) receiving abacavir 300 mg BID",
    dose_range     = "300 mg oral BID",
    regions        = "France",
    n_observations = 344,
    notes          = "Retrospective routine therapeutic-drug-monitoring data; 344 plasma samples (30 BLQ set to half the LOQ of 0.01 mg/L), median 1 sample per patient (range 1-7); 151 men and 37 women (19.7% female). Demographics from Jullien 2005 Table 1."
  )

  ini({
    # Structural parameters (apparent oral PK; F unknown).
    # Reference weight for the CL covariate term is the population median (65 kg).
    lka <- log(1.8);       label("Absorption rate constant (Ka, 1/h)")           # Table 3, TV(Ka) = 1.8 1/h
    lcl <- log(47.5);      label("Apparent clearance at 65 kg (CL/F, L/h)")      # Table 3, TV(CL/F) = 47.5 L/h
    lvc <- log(75);        label("Apparent central volume (Vc/F, L)")            # Table 3, TV(Vc/F) = 75 L
    lvp <- log(24);        label("Apparent peripheral volume (Vp/F, L)")         # Table 3, TV(Vp/F) = 24 L
    lq  <- fixed(log(10)); label("Apparent intercompartmental clearance (Q/F, L/h); fixed at the basic-model typical value when BW was added to CL/F")  # Table 3, TV(Q/F) = 10 L/h (fixed)

    # Estimated power-form covariate effect on CL/F.
    # CL/F = TV(CL/F) * (WT/65)^e_wt_cl
    e_wt_cl <- 0.80; label("Power exponent for body weight on CL/F (unitless)")  # Table 3, theta_BW = 0.80

    # IIV (exponential model: P_i = TV * exp(eta_i)); variances on the log scale.
    # IIV on Q/F and Vp/F could not be estimated (Results, p.185).
    # Correlated CL/F and Vc/F: omega^2_CL = 0.18, omega^2_Vc = 0.51, COV(CL,V) = 0.27
    # (Pearson r = 0.27 / sqrt(0.18 * 0.51) ~= 0.89).
    etalcl + etalvc ~ c(0.18,
                        0.27, 0.51)  # Table 3
    etalka          ~ 1.13           # Table 3, omega^2_Ka = 1.13

    # Residual error -- additive on the linear concentration scale.
    # sigma^2 = 0.039 (mg/L)^2 -> SD = sqrt(0.039) ~= 0.1975 mg/L.
    addSd <- 0.1975; label("Additive residual error (mg/L)")  # Table 3, sigma^2 = 0.039 (mg/L)^2
  })

  model({
    # Individual PK parameters; power-form weight effect on CL/F.
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) * (WT / 65)^e_wt_cl
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp)
    q  <- exp(lq)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                                k12 * central - k21 * peripheral1

    # Concentration: dose in mg, volumes in L -> mg/L (= ug/mL).
    Cc <- central / vc
    Cc ~ add(addSd)
  })
}
