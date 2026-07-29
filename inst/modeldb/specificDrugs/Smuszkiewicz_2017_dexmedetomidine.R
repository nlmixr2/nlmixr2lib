Smuszkiewicz_2017_dexmedetomidine <- function() {
  description <- "Two-compartment population PK model for intravenous dexmedetomidine continuous infusion in adult ICU patients undergoing analgosedation (Smuszkiewicz 2017). 27 medical and surgical ICU patients (17 male, 10 female; median age 59.5 y, median weight 75 kg) received continuous infusions of 0.1-1.5 ug/kg/h for 23.7-102 h. Age, sex, body weight, infusion duration, pretreatment SOFA score, and inotrope use were screened as covariates but none reached statistical significance, so the final model contains no covariate effects. IIVs on Vc, CL, Vp, and Q are diagonal (no clear correlations). Proportional residual error."
  reference   <- paste(
    "Smuszkiewicz P, Wiczling P, Ber J, Warzybok J, Malkiewicz T, Matysiak J,",
    "Klupczynska A, Trojanowska I, Kokot Z, Grzeskowiak E, Krzyzanski W, Bienert A.",
    "Pharmacokinetics of dexmedetomidine during analgosedation in ICU patients.",
    "J Pharmacokinet Pharmacodyn. 2018;45(2):277-284.",
    "doi:10.1007/s10928-017-9564-7.",
    sep = " "
  )
  vignette <- "Smuszkiewicz_2017_dexmedetomidine"
  units <- list(
    time          = "h",
    dosing        = "ug",
    concentration = "ng/mL"
  )

  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = 27L,
    n_studies      = 1L,
    age_range      = "19-84 years (median 59.5)",
    age_median     = "59.5 years",
    weight_range   = "45-100 kg (median 75)",
    weight_median  = "75 kg",
    sex_female_pct = 37.0,
    disease_state  = "Adults staying in a mixed medical and surgical intensive care unit; respiratory insufficiency requiring analgosedation and mechanical ventilation, or hyperactive delirium / agitation refractory to haloperidol in intubated and/or extubated ICU patients. Pretreatment SOFA score median 12 (range 5-16); inotrope use in 21 of 27.",
    dose_range     = "Continuous intravenous infusion of dexmedetomidine (Dexdor) without a loading dose. Initial infusion rate 0.8-1 ug/kg/h followed by maintenance infusion 0.4-1.5 ug/kg/h, titrated to a modified Ramsay sedation score of 2-3 and stopped on significant hemodynamic instability, after extubation, or at physician discretion. Median infusion duration 42.8 h (23.7-102), median total dose 1.55 mg (0.29-6.67), median rate 0.51 ug/kg/h (0.1-1.5).",
    regions        = "Poland (Poznan University of Medical Sciences)",
    notes          = "Single-center observational study; bioethics permissions 213/13 and 572/16. 368 dexmedetomidine concentrations from 27 patients (Fig. 1, Table 1). 22 patients sampled at 0, 1, 4, 8, 12, 16, 20 h after infusion start and 0, 5, 10, 20, 60 min and 2, 4, 6 h after infusion cessation; 5 patients followed a longer protocol (0, 2, 8, 24, 32, 48, 56, 72, 80 h on infusion and 0, 5, 10, 15, 30, 60 min and 2, 4, 6, 12 h post-infusion). Plasma assay: LC-MS/MS, calibration range 0.05-20 ng/mL, within-day CV < 10%, no measurements below the LLOQ. Model estimation: NONMEM 7.3.0 with gfortran 9.0, FOCE-style discrimination via likelihood-ratio test and AIC, 1000-replicate nonparametric bootstrap, prediction-corrected VPC. Covariates screened: WT, AGE, SEXF, infusion duration, pretreatment SOFA, inotrope use; none retained."
  )

  ini({
    # Structural parameters from Smuszkiewicz 2017 Table 2 'theta, estimate'
    # column (population typical values). Two-compartment IV-infusion model;
    # dexmedetomidine was administered as a continuous infusion into the
    # central compartment without a loading dose.
    lvc <- log(27.0);  label("Volume of central compartment (Vc, L)")             # Table 2 (theta_VC = 27.0 L)
    lcl <- log(38.5);  label("Systemic clearance (CL, L/h)")                      # Table 2 (theta_CL = 38.5 L/h)
    lvp <- log(87.6);  label("Volume of peripheral compartment (Vp, L)")          # Table 2 (theta_VT = 87.6 L)
    lq  <- log(46.4);  label("Inter-compartmental clearance (Q, L/h)")            # Table 2 (theta_Q  = 46.4 L/h)

    # Inter-individual variability. Table 2 reports each omega^2 as %CV; for a
    # log-normal IIV, omega^2 = log(CV^2 + 1) where CV is the decimal fraction.
    # Methods: "Inter-individual variability (IIV) for all PK parameters was
    # modeled assuming a lognormal distribution"; "There was no clear
    # correlations between the inter-individual random effect estimates" so
    # the four etas are uncorrelated (diagonal omega).
    etalvc ~ 0.9312  # Table 2 omega^2_VC = 124 %CV -> log(1 + 1.24^2)  = 0.9312
    etalcl ~ 0.3361  # Table 2 omega^2_CL = 63.2 %CV -> log(1 + 0.632^2) = 0.3361
    etalvp ~ 0.5832  # Table 2 omega^2_VT = 89.0 %CV -> log(1 + 0.89^2)  = 0.5832
    etalq  ~ 0.5034  # Table 2 omega^2_Q  = 80.9 %CV -> log(1 + 0.809^2) = 0.5034

    # Residual error: proportional only. Methods: "the additive or combined
    # additive and proportional error models did not lead to model
    # improvements based on AIC criterion." Table 2 reports sigma^2 as 24
    # %CV, which is the SD of the proportional residual on the linear scale
    # (sqrt(0.24^2) = 0.24); use the SD directly.
    propSd <- 0.24; label("Proportional residual error SD (fraction)")            # Table 2 (sigma^2 = 24 %CV -> propSd = 0.24)
  })

  model({
    # Individual PK parameters (log-normal IIV; no covariate effects retained
    # in the final model per Results: "None of the covariates were found to
    # be statistically significant in this study").
    vc <- exp(lvc + etalvc)
    cl <- exp(lcl + etalcl)
    vp <- exp(lvp + etalvp)
    q  <- exp(lq  + etalq)

    # Micro-rate constants for the explicit two-compartment ODE system.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Two-compartment open model with IV infusion into the central
    # compartment (cmt = central with a rate column on the dose row;
    # rxode2 adds the infusion rate to d/dt(central) over the dose
    # duration). Smuszkiewicz 2017 PK model:
    #   dCp/dt = Infusion(t)/Vc - (CL/Vc)*Cp - (Q/Vc)*(Cp - CT)
    #   dCT/dt = (Q/Vp)*(Cp - CT)
    # written here in amounts.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # Observation and proportional residual error (Methods: "e_prop,C
    # represents the proportional residual random error"; Cp = central / Vc).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
