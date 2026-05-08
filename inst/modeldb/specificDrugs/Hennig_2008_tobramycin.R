Hennig_2008_tobramycin <- function() {
  description <- "Two-compartment population PK model for once-daily IV tobramycin in paediatric cystic fibrosis patients (Hennig 2008), with allometric weight scaling on CL, Q, Vc, and Vper (reference 70 kg, exponent 3/4 for clearances and 1 for volumes), full-block correlated between-subject variability on CL/Vc/Vper, a fixed 30 min infusion duration into the central compartment, and an estimated lag time between infusion start and drug entry into the patient's vein."
  reference <- paste(
    "Hennig S, Norris R, Kirkpatrick CMJ. (2008).",
    "Target concentration intervention is needed for tobramycin dosing in paediatric patients",
    "with cystic fibrosis - a population pharmacokinetic study.",
    "Br J Clin Pharmacol 65(4):502-510. doi:10.1111/j.1365-2125.2007.03045.x"
  )
  vignette <- "Hennig_2008_tobramycin"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric scaling on CL, Q (exponent 3/4) and Vc, Vper (exponent 1) with reference weight 70 kg, per Hennig 2008 Table 2 footnote.",
      source_name        = "WT"
    )
  )

  population <- list(
    n_subjects     = 35L,
    n_studies      = 1L,
    age_range      = "0.5-17.9 years",
    age_mean       = "9.5 years",
    weight_range   = "6.0-72.6 kg",
    weight_mean    = "34.0 kg",
    height_range   = "60.0-178.0 cm",
    sex_female_pct = 60.0,
    disease_state  = "Paediatric cystic fibrosis patients receiving once-daily intravenous tobramycin for pulmonary Pseudomonas aeruginosa infection.",
    dose_range     = "70-560 mg/day (mean 9.6 mg/kg, range 6.9-15.2 mg/kg) once-daily IV infusion over 30 min.",
    n_observations = 318L,
    obs_per_subject = "9.1 (range 2-26)",
    occasions_per_subject = "4.6 (range 1-14)",
    serum_creatinine_range = "20.0-73.0 umol/L (mean 44.0)",
    crcl_range     = "23.0-194.7 mL/min (Cockcroft-Gault, mean 105.7); only one patient had CrCL <50 mL/min.",
    regions        = "Single-centre Australian cohort: Mater Health Services, Brisbane (data collected July 2005 - September 2006).",
    notes          = "Baseline demographics from Hennig 2008 Table 1. Renal impairment was rare in this cohort, so the model is most informative for paediatric CF patients with normal serum creatinine."
  )

  ini({
    # Structural parameters -- typical values referenced to a 70 kg patient.
    # All four PK parameters in Hennig 2008 Table 2 (Covariate model column) are
    # the typical values per 70 kg; allometric scaling per the table footnote
    # (f = 3/4 for clearance, f = 1 for volume) is applied to CL, Q, Vc, Vper.
    lcl    <- log(6.37);   label("Clearance for a 70 kg patient (CL, L/h)")           # Hennig 2008 Table 2 (Covariate model)
    lvc    <- log(18.70);  label("Central volume of distribution for a 70 kg patient (Vc, L)")  # Hennig 2008 Table 2 (Covariate model)
    lq     <- log(0.39);   label("Intercompartmental clearance for a 70 kg patient (Q, L/h)")   # Hennig 2008 Table 2 (Covariate model)
    lvp    <- log(1.32);   label("Peripheral volume of distribution for a 70 kg patient (Vper, L)")  # Hennig 2008 Table 2 (Covariate model)

    # Allometric exponents (Hennig 2008 Table 2 footnote: fixed)
    e_wt_cl_q  <- fixed(0.75); label("Allometric (WT) exponent on CL and Q (unitless)")  # Hennig 2008 Table 2 footnote
    e_wt_vc_vp <- fixed(1.00); label("Allometric (WT) exponent on Vc and Vper (unitless)")  # Hennig 2008 Table 2 footnote

    # Dose-record adjustments
    llag <- log(0.40);   label("Lag between infusion hang time and drug entry into vein (h)")  # Hennig 2008 Table 2 (tlag, Covariate model)
    ldur <- fixed(log(0.5)); label("Fixed infusion duration into central compartment (h, hospital protocol)")  # Hennig 2008 Table 2 (D2, fixed); Methods page 503

    # Inter-individual variability (full block on CL, Vc, Vper).
    # Conversion of CV% to log-scale variance: omega^2 = log(1 + CV^2).
    #   BSV CL  : 11.70% CV -> log(1 + 0.117^2)  = 0.013596
    #   BSV Vc  : 11.66% CV -> log(1 + 0.1166^2) = 0.013505
    #   BSV Vper: 41.95% CV -> log(1 + 0.4195^2) = 0.162130
    # Covariances from reported correlations (R) and the omega standard
    # deviations omega_i = sqrt(omega_i^2):
    #   omega_CL = 0.11660, omega_Vc = 0.11621, omega_Vper = 0.40265
    #   cov(CL, Vc)   = 0.73 * 0.11660 * 0.11621 = 0.009891
    #   cov(CL, Vper) = 0.49 * 0.11660 * 0.40265 = 0.023000
    #   cov(Vc, Vper) = 0.27 * 0.11621 * 0.40265 = 0.012632
    # nlmixr2 expects the lower triangle in row-major order:
    #   var_CL,
    #   cov(CL, Vc),  var_Vc,
    #   cov(CL, Vper), cov(Vc, Vper), var_Vper
    etalcl + etalvc + etalvp ~ c(0.013596,
                                 0.009891, 0.013505,
                                 0.023000, 0.012632, 0.162130)  # Hennig 2008 Table 2 (BSV CL/Vc/Vper, Covariate model + correlations)

    # Residual error -- exponential model in NONMEM (additive on log-scale)
    # maps to proportional in nlmixr2 linear space; CV = 19.0%.
    propSd <- 0.190; label("Proportional residual error (fraction)")  # Hennig 2008 Table 2 (Residual variability, Covariate model)
  })
  model({
    # Allometric scaling on every PK parameter (reference 70 kg);
    # WT enters as a continuous, optionally time-varying covariate.
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl_q
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc_vp
    q  <- exp(lq)           * (WT / 70)^e_wt_cl_q
    vp <- exp(lvp + etalvp) * (WT / 70)^e_wt_vc_vp

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Fixed 30 min nominal infusion duration and estimated infusion lag
    dur(central) <- exp(ldur)
    lag(central) <- exp(llag)

    # Concentration in central compartment: dose in mg, volume in L -> mg/L
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
