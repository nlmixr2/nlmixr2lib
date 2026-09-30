Umpierrez_2021_cyclosporine <- function() {
  description <- "Two-compartment oral population PK model with lagged first-order absorption for whole-blood cyclosporine in Uruguayan transplant and autoimmune-disease patients monitored at steady state, with a power effect of Cockcroft-Gault creatinine clearance on CL/F, between-subject variability on CL/F and Q/F, and correlated inter-occasion variability on Ka and CL/F plus inter-occasion variability on the lag time (Umpierrez 2021)"
  reference <- paste(
    "Umpierrez M, Guevara N, Ibarra M, Fagiolino P, Vazquez M, Maldonado C.",
    "Development of a population pharmacokinetic model for cyclosporine from",
    "therapeutic drug monitoring data. Biomed Res Int. 2021;2021:3108749.",
    "doi:10.1155/2021/3108749.",
    "Parameter estimates are Umpierrez 2021 Table 2; the variability and",
    "residual-error forms are Methods equations (1)-(3) and the creatinine",
    "clearance relationship is Results equation (8).",
    sep = " "
  )
  vignette <- "Umpierrez_2021_cyclosporine"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    CRCL = list(
      description = "Creatinine clearance estimated from serum creatinine with the Cockcroft-Gault formula, raw mL/min (NOT BSA-normalized)",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on CL/F: CL = CLpop * (CRCL / 98.62)^(-0.204) (Umpierrez 2021 Results equation (8) and the Table 2 row 'beta CL-CLCr', which prints the same relationship on the log scale). The reference 98.62 mL/min is the Group A (model-building) cohort MEAN, not a median (Umpierrez 2021 Table 1 'Cl Cr (mean, mL/min) 98.62 (417.92-13.79)', restated in the Results as 'the mean value for creatinine clearance (98.62 mL/min)'). The exponent is negative, so CL/F rises as renal function falls; the paper attributes this to a redistribution of cardiac output toward the splanchnic (metabolising) region when renal blood flow is reduced. The Cockcroft-Gault value is not normalized to body surface area in the source (Methods section 2.3.1).",
      source_name = "ClCrea"
    ),
    OCC = list(
      description = "Integer-valued occasion indicator (one steady-state monitoring visit = one occasion) for inter-occasion variability on Ka, CL/F and the lag time",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = "Occasions 1 to 4 each carry their own inter-occasion random effects, all sharing the magnitudes of Umpierrez 2021 Table 2 (occasions 2-4 are fixed to the occasion-1 values). The paper does not print the number of occasions per patient; four slots cover the model-building data (621 samples in 37 patients at up to five samples per profile is about 3.4 profiles per patient) and the three occasions of the prospective evaluation (Umpierrez 2021 Figure 2). For a longer monitoring history, re-use the four values cyclically: the occasion slots are exchangeable. OCC = 0 switches inter-occasion variability off.",
      source_name = "OCC"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened in the stepwise covariate analysis but not retained (Umpierrez 2021 Methods section 2.3.1 and Results)."
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened in the stepwise covariate analysis but not retained (Umpierrez 2021 Methods section 2.3.1 and Results)."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      notes = "Screened in the stepwise covariate analysis but not retained (Umpierrez 2021 Methods section 2.3.1 and Results)."
    )
  )

  compartmentData <- list(
    depot = list(analyte = "cyclosporine", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "cyclosporine", units = "mg", specimen = "whole blood", verified = TRUE),
    peripheral1 = list(analyte = "cyclosporine", units = "mg", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 37L,
    n_studies = 1L,
    age_range = "mean 34.4 years (SD 15.85)",
    weight_range = "mean 64.3 kg (SD 11.0)",
    sex_female_pct = 56.8,
    race_ethnicity = "Uruguayan patients; race not reported",
    disease_state = "Patients receiving oral cyclosporine at steady state for kidney transplantation (14), kidney autoimmune disease (21), liver autoimmune disease (1) or bone-marrow transplantation (1)",
    dose_range = "Individualised oral cyclosporine maintenance regimens from routine therapeutic drug monitoring; doses not reported",
    regions = "Uruguay (Montevideo)",
    renal_function = "Cockcroft-Gault creatinine clearance mean 98.62 mL/min (range 13.79-417.92); serum creatinine mean 1.11 mg/dL (range 0.2-5.9)",
    notes = "Model-building cohort (Group A) of Umpierrez 2021 Table 1: 37 patients (16 male / 21 female) with at least one four-sample steady-state profile, 621 whole-blood concentrations sampled at 0, 1, 2, 3 and 4 h postdose and assayed by chemiluminescent microparticle immunoassay (LLOQ 12.5 ng/mL). A separate prospective cohort (Group B: 16 patients, 81 concentrations) was used only to evaluate Bayesian-forecasting performance. Estimation used SAEM in Monolix 2019R1."
  )

  ini({
    # Structural parameters, Umpierrez 2021 Table 2. All are apparent (divided
    # by the unknown oral bioavailability F). Reference subject: CRCL = 98.62 mL/min.
    ltlag <- log(0.512); label("Absorption lag time (h)")                              # Umpierrez 2021 Table 2: Tlag = 0.512 h (RSE 8.48%)
    lka   <- log(0.523); label("First-order absorption rate constant Ka (1/h)")        # Umpierrez 2021 Table 2: Ka = 0.523 1/h (RSE 8.54%)
    lcl   <- log(30.3);  label("Apparent clearance CL/F at CRCL = 98.62 mL/min (L/h)") # Umpierrez 2021 Table 2: Cl/F = 30.3 L/h (RSE 8.25%)
    lq    <- log(17.0);  label("Apparent intercompartmental clearance Q/F (L/h)")      # Umpierrez 2021 Table 2: Q/F = 17.0 L/h (RSE 12.1%)
    lvc   <- log(17.9);  label("Apparent central volume V1/F (L)")                     # Umpierrez 2021 Table 2: V1 = 17.9 L (RSE 17.6%)
    lvp   <- log(400);   label("Apparent peripheral volume V2/F (L)")                  # Umpierrez 2021 Table 2: V2 = 400 L (RSE 45.6%)

    # Creatinine clearance on CL/F: log(CLi) = log(CLpop) + beta * log(CRCL / 98.62) + eta.
    e_crcl_cl <- -0.204; label("Power exponent of (CRCL / 98.62) on CL/F (unitless)") # Umpierrez 2021 Table 2: beta CL-CLCr = -0.204 (RSE 43.7%); Results equation (8)

    # Between-subject variability. Table 2 prints CV%, derived from the
    # log-normal variance by Methods equation (2), CV = 100 * sqrt(exp(omega^2) - 1),
    # so omega^2 = log(1 + CV^2).
    etalcl ~ 0.14704 # Umpierrez 2021 Table 2: IIV Cl = 39.8% (RSE 16.6%); omega^2 = log(1 + 0.398^2)
    etalq ~ 0.24922  # Umpierrez 2021 Table 2: IIV Q = 53.2% (RSE 18.7%); omega^2 = log(1 + 0.532^2)

    # Inter-occasion variability, same equation (2) conversion. Ka carries no
    # between-subject variability, so the Table 2 Ka-Cl correlation (-0.551) can
    # only sit at the occasion level (Monolix correlates random effects within
    # one variability level): cov = -0.551 * sqrt(0.24426 * 0.13488) = -0.10001.
    # Occasions 2-4 share the occasion-1 magnitudes (fixed).
    etaiov_lka_1 + etaiov_lcl_1 ~ c(
      0.24426,
      -0.10001, 0.13488
    ) # Umpierrez 2021 Table 2: IOV ka = 52.6% (RSE 9.85%), IOV Cl = 38.0% (RSE 8.53%), Ka-Cl correlation = -0.551 (RSE 20.2%)
    etaiov_lka_2 + etaiov_lcl_2 ~ c(
      fixed(0.24426),
      fixed(-0.10001), fixed(0.13488)
    ) # occasion 2 shares the occasion-1 covariance matrix
    etaiov_lka_3 + etaiov_lcl_3 ~ c(
      fixed(0.24426),
      fixed(-0.10001), fixed(0.13488)
    ) # occasion 3 shares the occasion-1 covariance matrix
    etaiov_lka_4 + etaiov_lcl_4 ~ c(
      fixed(0.24426),
      fixed(-0.10001), fixed(0.13488)
    ) # occasion 4 shares the occasion-1 covariance matrix
    etaiov_ltlag_1 ~ 0.25672        # Umpierrez 2021 Table 2: IOV tlag = 54.1% (RSE 12.4%); omega^2 = log(1 + 0.541^2)
    etaiov_ltlag_2 ~ fixed(0.25672) # occasion 2 shares the occasion-1 variance
    etaiov_ltlag_3 ~ fixed(0.25672) # occasion 3 shares the occasion-1 variance
    etaiov_ltlag_4 ~ fixed(0.25672) # occasion 4 shares the occasion-1 variance

    # Combined residual error, Methods equation (3): Y = C * (1 + e_prop) + e_add
    # with independent normal e_prop and e_add (variances summed).
    propSd <- 0.228; label("Proportional residual error (fraction)") # Umpierrez 2021 Table 2: Prop = 0.228 (RSE 8.86%); the row header reads '(%)' but the value is the fraction 0.228 (22.8%)
    addSd  <- 7.52;  label("Additive residual error (ng/mL)")        # Umpierrez 2021 Table 2: Add = 7.52 ng/mL (RSE 45.0%)
  })

  model({
    # Inter-occasion variability, multiplexed by the occasion indicator. OCC = 0
    # zeroes every indicator and switches inter-occasion variability off.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    iov_ka <- oc1 * etaiov_lka_1 + oc2 * etaiov_lka_2 +
      oc3 * etaiov_lka_3 + oc4 * etaiov_lka_4
    iov_cl <- oc1 * etaiov_lcl_1 + oc2 * etaiov_lcl_2 +
      oc3 * etaiov_lcl_3 + oc4 * etaiov_lcl_4
    iov_tlag <- oc1 * etaiov_ltlag_1 + oc2 * etaiov_ltlag_2 +
      oc3 * etaiov_ltlag_3 + oc4 * etaiov_ltlag_4

    # Individual parameters (Methods equation (1), exponential random effects;
    # Results equation (8) for the creatinine clearance effect).
    tlag <- exp(ltlag + iov_tlag)
    ka <- exp(lka + iov_ka)
    cl <- exp(lcl + etalcl + iov_cl) * (CRCL / 98.62)^e_crcl_cl
    q <- exp(lq + etalq)
    vc <- exp(lvc)
    vp <- exp(lvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    alag(depot) <- tlag

    # Dose in mg, volume in L: mg/L = ug/mL; x 1000 gives ng/mL.
    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
