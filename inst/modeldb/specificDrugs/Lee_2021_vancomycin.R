Lee_2021_vancomycin <- function() {
  description <- "One-compartment intravenous population PK model for vancomycin in Korean neonates treated in a neonatal intensive care unit (Lee 2021), developed from routine peak and trough therapeutic-drug-monitoring concentrations. Clearance scales allometrically with body weight (fixed exponent 0.75, 70 kg reference) and as power functions of postmenstrual age (reference 31.7 weeks) and Schwartz creatinine clearance (reference 50.3 mL/min/1.73 m^2); volume of distribution scales linearly with body weight (fixed exponent 1, 70 kg reference). Correlated log-normal between-subject variability on clearance and volume and a combined additive-plus-proportional residual error."
  reference <- "Lee SM, Yang S, Kang S, Chang MJ. Population pharmacokinetics and dose optimization of vancomycin in neonates. Sci Rep. 2021;11:6168. doi:10.1038/s41598-021-85529-3"
  vignette <- "Lee_2021_vancomycin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Vancomycin is given as an intravenous infusion (typically 1 h; Lee 2021
  # Results 'Patients' and Methods 'Data collection'), so the dose enters
  # `central` directly. Lee 2021 Methods 'Data collection' states that serum
  # vancomycin concentrations were measured by chemiluminescent microparticle
  # immunoassay (Abbott ARCHITECT i).
  compartmentData <- list(
    central = list(analyte = "vancomycin", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight at the time of the vancomycin record",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying in the deposited NONMEM dataset (Supplementary MOESM4). Enters CL as (WT/70)^0.75 and V as (WT/70)^1 (Lee 2021 Equations 1-2 and 3-4; Table 2 'Final model'); both exponents are fixed, not estimated (Abstract: 'using fixed powers (0.75 and 1 ...)'). The 70 kg reference is a standardised-adult scaling convention; the cohort itself spans 0.5-5.9 kg (Table 1, median 1.8 kg).",
      source_name = "WT"
    ),
    PAGE = list(
      description = "Postmenstrual age (gestational age plus postnatal age)",
      units = "weeks",
      type = "continuous",
      reference_category = NULL,
      notes = "Lee 2021 writes PMA in WEEKS (Table 1; Equation 1 reference 31.7), so this model declares PAGE in weeks rather than the register default of months (see the PAGE register entry). Enters CL as (PMA/31.7)^0.795. The reference 31.7 weeks is the median PMA across the records of the deposited dataset (31.69; Supplementary MOESM4), not the per-subject Table 1 median of 35.6 weeks.",
      source_name = "PMA"
    ),
    CRCL = list(
      description = "Creatinine clearance estimated with the Schwartz equation, BSA-normalised",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Lee 2021 Methods 'Data collection': CLcr (mL/min/1.73 m^2) = length (cm) x k / SCr (mg/dL), with k = 0.45 for infants 1 to 52 weeks old (Schwartz). Enters CL as (CLcr/50.3)^0.741; 50.3 is the Table 1 median and also the median across records of the deposited dataset (50.29). Serum creatinine by kinetic Jaffe (compensated, IDMS-traceable).",
      source_name = "CLCR"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 207L,
    n_studies = 1L,
    n_sites = 1L,
    n_observations = 900L,
    age_range = "Gestational age 23.3-41.5 weeks (median 31.5); postnatal age at examination 0-16.4 weeks (median 2.3); postmenstrual age at examination 24.0-48.4 weeks (median 35.6)",
    weight_range = "Body weight at examination 0.5-5.9 kg (median 1.8); birth weight 0.5-5.4 kg (median 1.5)",
    sex_female_pct = 50,
    race_ethnicity = "Korean",
    disease_state = "Neonates in a neonatal intensive care unit treated with vancomycin for more than 24 h (suspected or confirmed gram-positive sepsis), with at least one steady-state concentration. Neonates with acute kidney injury before vancomycin initiation (urine output < 1 mL/kg/h or serum creatinine >= 1.5 mg/dL) were excluded.",
    dose_range = "Mostly 10 mg/kg as a 1 h intravenous infusion every 8 or 12 h per Neofax (interval by PMA and postnatal age); lower or higher doses in individual patients",
    regions = "South Korea (Gangnam Severance Hospital, Seoul)",
    renal_function = "Schwartz creatinine clearance median 50.3 (range 6.8-140.3) mL/min/1.73 m^2; serum creatinine median 0.5 (0.2-2.6) mg/dL",
    notes = "Retrospective chart review of neonates admitted January 2008 to April 2017 (Lee 2021 Methods, Table 1). Only peak (1 h after the end of infusion) and trough (0.5 h before the next dose) concentrations were sampled. Lower limit of quantification 3.0 mg/L; below-quantification values were kept as observed. Estimated in NONMEM 7.3 with FOCE-I; the final control stream and the full analysis dataset are deposited as Supplementary MOESM3 and MOESM4."
  )

  ini({
    # ===== Structural PK: Lee 2021 Table 2 and Equations (1)-(2) =====
    #   TVCL (L/h) = theta1 * (WT/70)^0.75 * (PMA/31.7)^theta2 * (CLcr/50.3)^theta3
    #   TVV  (L)   = theta4 * (WT/70)^1
    # Units L/h and L are printed in the Table 2 footnote ('TVCL typical value
    # of clearance (L/h), TVV typical value of volume of distribution (L)').
    lcl <- log(2.09); label("Clearance at WT = 70 kg, PMA = 31.7 weeks, CLcr = 50.3 mL/min/1.73 m^2 (L/h)") # Table 2, theta1 = 2.09 (RSE 3%); Equation (1)
    lvc <- log(45.6); label("Volume of distribution at WT = 70 kg (L)") # Table 2, theta4 = 45.6 (RSE 5%); Equation (2)

    # Allometric (CL) and isometric (V) body-weight exponents, fixed a priori
    # (Abstract: 'using fixed powers (0.75 and 1, respectively, for clearance
    # and volume)'; control stream MOESM3 hard-codes both).
    e_wt_cl <- fixed(0.75); label("Allometric exponent of (WT / 70 kg) on CL (unitless)") # Equation (1) / Equation (3): exponent 0.75, fixed
    e_wt_vc <- fixed(1); label("Isometric exponent of (WT / 70 kg) on V (unitless)") # Equation (2) / Equation (4): exponent 1, fixed

    e_page_cl <- 0.795; label("Power exponent of (PMA / 31.7 weeks) on CL (unitless)") # Table 2, theta2 = 0.795 (RSE 26%); Equation (1)
    e_crcl_cl <- 0.741; label("Power exponent of (CLcr / 50.3 mL/min/1.73 m^2) on CL (unitless)") # Table 2, theta3 = 0.741 (RSE 6%); Equation (1)

    # ===== Between-subject variability =====
    # Exponential (log-normal) etas on CL and V with an OMEGA BLOCK(2)
    # (Results 'Pharmacokinetic modeling': correlation between CL and V,
    # dOFV = -47.9; control stream MOESM3 '$OMEGA BLOCK(2)'). Table 2 prints
    # the two diagonal elements as variances (see the vignette for the check
    # against the supplementary PTA tables that rules out the SD reading) but
    # does not print the covariance. The covariance below is
    # 0.48 * sqrt(0.123 * 0.260) = 0.0858. The correlation 0.48 is the value
    # that minimises the FOCEi objective function on the deposited dataset
    # (MOESM4) with every other Table 2 value held fixed. A free re-estimation
    # of the deposited control stream (MOESM3) gives 0.51. See the vignette.
    etalcl + etalvc ~ c(0.123, 0.0858, 0.260) # Table 2 'omega CL' 0.123 and 'omega V' 0.260 (variances); covariance not printed (correlation 0.48 profiled on the deposited data)

    # ===== Residual error =====
    # Control stream MOESM3: W = SQRT(THETA(3)**2 + THETA(4)**2 * IPRED**2),
    # Y = IPRED + W * EPS(1), SIGMA 1 FIX -- a combined error whose additive
    # and proportional SDs are thetas and add in variance (nlmixr2 combined2,
    # the default for prop() + add()). Kept as printed, although the deposited
    # dataset supports a smaller proportional SD (about 0.39); see the vignette.
    propSd <- 0.583; label("Proportional residual error (fraction)") # Table 2, 'sigma proportional (%CV)' = 58.3%
    addSd <- 2.015; label("Additive residual error (mg/L)") # Table 2, 'sigma additive (mg/L)' = 2.015
  })
  model({
    # ----- Individual PK parameters (Lee 2021 Equations 1-2) -----
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl * (PAGE / 31.7)^e_page_cl * (CRCL / 50.3)^e_crcl_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc

    # ----- Micro-constant and ODE (NONMEM ADVAN1 TRANS2) -----
    kel <- cl / vc
    d/dt(central) <- -kel * central

    # ----- Output -----
    # Dose in mg and vc in L give mg/L (= ug/mL, the units of the dataset DV).
    Cc <- central / vc
    Cc ~ prop(propSd) + add(addSd)
  })
}
