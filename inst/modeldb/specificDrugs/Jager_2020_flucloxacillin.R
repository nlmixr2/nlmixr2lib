Jager_2020_flucloxacillin <- function() {
  description <- "Two-compartment joint total/unbound population PK model for intravenous flucloxacillin in critically ill adults (Jager 2020). The disposition is carried on the UNBOUND concentration Cu = central/V1 with linear unbound clearance and linear intercompartmental exchange; the measured TOTAL plasma concentration is reconstructed as Cc = Cu + Bmax * Cu / (Kd + Cu), a single saturable binding site whose capacity Bmax rises with serum albumin. Total and unbound concentrations are BOTH observed endpoints, each with its own additive error on the log scale. Unbound clearance carries a power term on CKD-EPI eGFR. Inter-individual variability on CL and Bmax."
  reference <- "Jager NGL, van Hest RM, Xie J, Wong G, Ulldemolins M, Bruggemann RJM, Lipman J, Roberts JA. Optimization of flucloxacillin dosing regimens in critically ill patients using population pharmacokinetic modelling of total and unbound concentrations. J Antimicrob Chemother. 2020;75(9):2641-2649. doi:10.1093/jac/dkaa187. PMCID PMC7443729. Binding equation from main-text Equation 3 (Supplementary Appendix 2 Equation 1); covariate equations from main-text Equations 4-5; parameter estimates from main-text Table 2 (final model column); unit conversion (flucloxacillin MW 453.9 g/mol) and residual-error form from Supplementary Appendix 2."
  vignette <- "Jager_2020_flucloxacillin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(
      analyte = "flucloxacillin (unbound)",
      units = "mg",
      specimen = "plasma",
      verified = TRUE,
      notes = "Holds the UNBOUND flucloxacillin amount in the central compartment; dividing by V1 gives the unbound plasma concentration Cu. Supplementary Appendix 2 states that the PK parameters were 'estimated for unbound flucloxacillin' and that total concentrations were obtained 'by combining the PK of unbound flucloxacillin with protein binding parameters', and main-text Figure 3b labels CL as the clearance 'of unbound flucloxacillin'. The main-text Results confirm the consequence of this structure: in the Monte Carlo simulations 'serum albumin concentrations did not affect unbound concentrations' (Figure 4b and d), which holds only when the ODE runs on the unbound amount and albumin enters solely through the bound term of the total-concentration observation."
    ),
    peripheral1 = list(
      analyte = "flucloxacillin (unbound)",
      units = "mg",
      specimen = "tissue",
      verified = TRUE,
      notes = "Unbound peripheral amount; exchanges with central at intercompartmental clearance Q (Table 2)."
    )
  )

  covariateData <- list(
    ALB = list(
      description = "Serum albumin concentration",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Reference 20 g/L, printed as the literal denominator of main-text Equation 4, Bmax (mmol/L) = 0.469 * (Albumin/20)^1.51. The equation is in g/L: the cohort median is 21 g/L (Table 1, IQR 15-34) and the Monte Carlo simulations used the 10th/90th percentiles 15 and 30 g/L. Enters on the saturable binding capacity Bmax only, so albumin changes TOTAL but not unbound concentrations. Every patient in the model-building and both external cohorts had albumin below 35 g/L (Discussion, limitations); the authors caution against extrapolating above that.",
      source_name = "Albumin"
    ),
    CRCL = list(
      description = "Estimated glomerular filtration rate by the CKD-EPI creatinine equation",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Reference 90, printed as the literal denominator of main-text Equation 5, CL (L/h) = 55.4 * (eGFR/90)^0.809. The paper labels eGFR in mL/min throughout (Table 1, Figure 3b, Methods), but the covariate is stated to be the CKD-EPI creatinine equation (Abstract; Supplementary Appendix 2), whose native output is BSA-normalized mL/min/1.73 m^2; no de-normalization step is described, so the column is stored under the BSA-normalized canonical. Cohort median 96 (IQR 26-166). Four model-building patients (11%) were on renal replacement therapy; RRT was tested and not retained (Discussion), so the model carries no dialysis term.",
      source_name = "eGFR"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened (Supplementary Appendix 2, covariate model) and not significant. Cohort median 95 kg (IQR 73-120)."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Significant on V1 in the univariate step (Supplementary Appendix 2) but dropped in backward elimination; not in the final model."
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened (Supplementary Appendix 2) and not significant."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = male",
      notes = "Screened as a categorical covariate (Supplementary Appendix 2) and not significant."
    ),
    RRT = list(
      description = "Renal replacement therapy indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = no RRT",
      notes = "Significant on CL in the univariate step (Supplementary Appendix 2) but not retained once eGFR was in the model; the Discussion attributes this to only four RRT patients."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 35,
    n_studies = 2,
    age_range = "median 52 years (IQR 43-67)",
    weight_range = "median 95 kg (IQR 73-120); BMI median 31 kg/m^2 (IQR 25-35)",
    sex_female_pct = 34,
    race_ethnicity = "Not reported (single-centre Australian ICU)",
    disease_state = "Critically ill adults in a tertiary-referral ICU treated with intermittent intravenous flucloxacillin for (suspected) infection, mostly bloodstream or respiratory; median SOFA score 8 (IQR 5-13).",
    renal_function = "Median eGFR (CKD-EPI) 96 (IQR 26-166); 4 patients (11%) on renal replacement therapy.",
    albumin = "Median serum albumin 21 g/L (IQR 15-34); all patients below 35 g/L.",
    dose_range = "Intermittent 30-minute infusions from 1 g q6h to 2 g q2h, at the discretion of the treating intensivist",
    regions = "Australia (Royal Brisbane and Women's Hospital ICU); pooled from a prospective PK study (2009, n = 10 hypoalbuminaemic patients) and a beta-lactam therapeutic-drug-monitoring programme (2012-2014, n = 25)",
    notes = "79 total and 104 unbound concentrations (total 1.0-202 mg/L, unbound 0.1-30 mg/L) sampled at a median 30 h after the start of treatment; both total and unbound measured in 16 patients, unbound only in 19. Observed protein binding 63.4-97.2% (median 89.4%). Fitted in NONMEM 7.1.2 (FOCE-I, ADVAN6) on log-transformed data; 1000-sample bootstrap. External validation in a Brisbane (n = 20, unbound only) and a Nijmegen, Netherlands (n = 14, total and unbound) ICU cohort."
  )

  ini({
    # ----------------------------------------------------------------
    # Unbound disposition (main-text Table 2, "Final model" column;
    # estimate [RSE%]). All four are UNBOUND-referenced -- see
    # compartmentData$central.
    # ----------------------------------------------------------------
    lcl <- log(55.4); label("Unbound clearance CL at eGFR 90 (L/h)")               # Table 2 final model CL = 55.4 L/h [11.4%]; bootstrap 55.2 (42.8-68.1); Equation 5
    lvc <- log(52.7); label("Unbound central volume of distribution V1 (L)")       # Table 2 final model V1 = 52.7 L [12.0%]; bootstrap 53.3 (36.9-68.6)
    lvp <- log(56.8); label("Unbound peripheral volume of distribution V2 (L)")    # Table 2 final model V2 = 56.8 L [11.8%]; bootstrap 57.3 (41.6-72.0)
    lq  <- log(67.2); label("Unbound intercompartmental clearance Q (L/h)")        # Table 2 final model Q = 67.2 L/h [26.0%]; bootstrap 65.6 (32.6-102)

    # ----------------------------------------------------------------
    # Saturable binding (main-text Equation 3). Bmax and Kd are in
    # mmol/L as printed; model() converts with MW 453.9 g/mol
    # (Supplementary Appendix 2, structural model).
    # ----------------------------------------------------------------
    lbmax_pb <- log(0.469);  label("Maximum binding capacity Bmax at albumin 20 g/L (mmol/L)")  # Table 2 final model Bmax = 0.469 mmol/L [14.1%]; bootstrap 0.478 (0.316-0.622); Equation 4
    lkd_pb   <- log(0.0441); label("Equilibrium dissociation constant Kd (mmol/L)")             # Table 2 final model Kd = 0.0441 mmol/L [16.6%]; bootstrap 0.0450 (0.0260-0.0621)

    # ----------------------------------------------------------------
    # Covariate exponents (main-text Equations 4-5; Table 2 'Covariates')
    # ----------------------------------------------------------------
    e_alb_bmax_pb <- 1.51;  label("Exponent on (ALB / 20) for Bmax")    # Table 2 covariate 'albumin' = 1.51 [28.6%]; bootstrap 1.52 (0.521-2.50); Equation 4
    e_crcl_cl     <- 0.809; label("Exponent on (CRCL / 90) for unbound CL") # Table 2 covariate 'eGFR' = 0.809 [24.2%]; bootstrap 0.809 (0.365-1.02); Equation 5

    # ----------------------------------------------------------------
    # Between-patient variability (Table 2 'BPV', %CV; exponential per
    # Supplementary Appendix 2). Converted as omega^2 = log(CV^2 + 1):
    #   Bmax 30.4% -> log(1 + 0.304^2) = 0.088384
    #   CL   71.6% -> log(1 + 0.716^2) = 0.413876
    # BPV on Kd and V1 was tested and not supported (Appendix 2).
    # ----------------------------------------------------------------
    etalbmax_pb ~ 0.088384 # Table 2 final model BPV Bmax = 30.4 %CV [19.2%]
    etalcl      ~ 0.413876 # Table 2 final model BPV CL = 71.6 %CV [15.8%]

    # ----------------------------------------------------------------
    # Residual error. Supplementary Appendix 2: 'additive error model for
    # logarithmically transformed data ... 0.16 for total flucloxacillin
    # and 0.22 for unbound flucloxacillin', i.e. log-scale SDs (the
    # appendix's '17-18%' and '24-25%' equal exp(SD) - 1). Table 2 heads
    # the same numbers 'proportional error', the linear-scale reading of
    # a log-additive error. Cc (total) takes the suffix-free name.
    # ----------------------------------------------------------------
    expSd    <- 0.160; label("Log-scale additive residual SD on total concentration Cc") # Table 2 final model 'proportional error, total flucloxacillin' = 0.160 [11.6%]
    expSd_Cu <- 0.222; label("Log-scale additive residual SD on unbound concentration Cu") # Table 2 final model 'proportional error, unbound flucloxacillin' = 0.222 [11.2%]
  })

  model({
    # 1. Individual unbound disposition parameters. CL carries the
    #    CKD-EPI eGFR power term (Equation 5); V1, V2 and Q are typical
    #    values only.
    cl <- exp(lcl + etalcl) * (CRCL / 90)^e_crcl_cl
    vc <- exp(lvc)
    vp <- exp(lvp)
    q  <- exp(lq)

    # 2. Binding capacity with the albumin power term (Equation 4) and
    #    the dissociation constant, both in mmol/L.
    bmax_pb <- exp(lbmax_pb + etalbmax_pb) * (ALB / 20)^e_alb_bmax_pb
    kd_pb   <- exp(lkd_pb)

    # 3. Micro-constants
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 4. Two-compartment ODE system on the UNBOUND amount. Flucloxacillin
    #    is given as an intravenous infusion straight into central.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # 5. Observations. Equation 3, Ctotal = Cu + Bmax * Cu / (Kd + Cu),
    #    evaluated in mmol/L as fitted (Supplementary Appendix 2
    #    converted all concentrations with MW 453.9 g/mol) and reported
    #    back in mg/L.
    mw_fluclox <- 453.9
    Cu    <- central / vc
    Cu_mm <- Cu / mw_fluclox
    Cc    <- (Cu_mm + bmax_pb * Cu_mm / (kd_pb + Cu_mm)) * mw_fluclox

    Cc ~ lnorm(expSd)
    Cu ~ lnorm(expSd_Cu)
  })
}
