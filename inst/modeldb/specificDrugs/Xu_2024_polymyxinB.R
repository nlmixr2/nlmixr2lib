Xu_2024_polymyxinB <- function() {
  description <- "Two-compartment IV infusion population PK model for polymyxin B (summed PMB1 + PMB2) in 11 critically ill adults on venoarterial or venovenous ECMO with carbapenem-resistant Gram-negative infections (Xu 2024). Clearance is 3.75 L/h at an ECMO blood flow of 3.48 L/min and falls linearly by 32.1% per L/min of additional flow, CL = 3.75 * (1 - 0.321 * (Q_ECMO - 3.48)). The central volume is 20.41 L in a patient receiving CRRT and 52.5% smaller (9.69 L) without CRRT. Peripheral volume 9.86 L and intercompartmental clearance 3.82 L/h carry no covariates. Exponential IIV on CL (omega^2 = 0.038) and on the central volume (omega^2 = 0.062); additive residual error 0.361 mg/L."
  reference <- "Xu M, Chen N, Yu YW, Pan XY, Li T. Pharmacokinetic Changes and Influencing Factors of Polymyxin B in Different ECMO Modes. Infect Drug Resist. 2024;17:5815-5825. doi:10.2147/IDR.S486169"
  vignette <- "Xu_2024_polymyxinB"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Total polymyxin B (PMB1 + PMB2 summed) in plasma by LC-MS/MS (Xu 2024
  # Drug Assay section); the peripheral compartment is the model's tissue
  # distribution space.
  compartmentData <- list(
    central = list(analyte = "polymyxinB", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "polymyxinB", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    Q_ECMO = list(
      description = "ECMO circuit blood flow rate",
      units = "L/min",
      type = "continuous",
      reference_category = NULL,
      notes = "Source covariate 'speadECMO' in the final-model CL equation, defined in the Table 2 footnote as 'spead of ECMO (L/min)'; the L/min unit identifies it as the delivered blood flow rather than the pump speed (Table 1 reports pump speed separately in RPM, median 3500). Table 1 median (IQR) 3.60 (3.00, 3.65) L/min. Enters CL as the linear centred form 1 + dCLdspeedECMO * (Q_ECMO - 3.48), with 3.48 L/min the centring value printed in the equation (not the Table 1 median of 3.60). Treated as a per-subject value at the steady-state sampling occasion; the paper does not describe flow as time-varying. Because the effect is linear with a negative slope (-0.321 per L/min), predicted CL reaches zero at Q_ECMO = 3.48 + 1/0.321 = 6.60 L/min; do not extrapolate far outside the observed 3.0-3.65 L/min interquartile range.",
      source_name = "speadECMO"
    ),
    RRT_CRRT_STATUS = list(
      description = "Subject-level binary indicator for continuous renal replacement therapy during the PK sampling period",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Source covariate CRRT. Table 1: 8/11 (72.73%) received CRRT, all of it continuous venovenous haemodialysis (CVVHD; Ultraflux AV1000S) for anuria, connected in the ECMO circuit (outflow at the post-oxygenator, inflow at the pre-oxygenator line). The paper estimates the effect in the Phoenix 'dVdCRRT0' parameterization, V = tvV * (1 + dVdCRRT0 * (if no CRRT)), so the typical central volume of 20.41 L is the CRRT value and the -0.525 coefficient applies to patients WITHOUT CRRT. Encoded here as (1 + e_nocrrt_vc * (1 - RRT_CRRT_STATUS)). Direction: CRRT is associated with a larger central volume, as the Discussion states ('patients receiving CRRT-ECMO had significantly increased the central compartment').",
      source_name = "CRRT"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Subject age",
      units = "years",
      type = "continuous",
      notes = "The Results text says 'CRRT and age were identified as the covariates for V', and the Table 2 abbreviation list defines 'dVdage, fixed parameter coefficient of age (year) to V', but neither the printed final-model V equation nor the Table 2 parameter rows contain an age term or estimate. The printed equation is taken as the final model, so age is not carried. Table 1 median (IQR) 56 (48, 64) years; the Results prose gives 60 (53.50-62.50)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 11L,
    n_studies = 1L,
    age_range = "Adults (>= 18 years); Table 1 median (IQR) 56 (48, 64) years (the Results prose reports 60, 53.50-62.50)",
    weight_range = "Median (IQR) 60 (57.05, 66.00) kg; BSA 1.69 (1.68, 1.79) m^2 (Table 1)",
    sex_female_pct = 27.27,
    race_ethnicity = "Not reported (single-centre Chinese ICU cohort)",
    disease_state = "Heart or respiratory failure requiring ECMO (8 VA-ECMO, 3 VV-ECMO; SOFA >= 2) with infection by carbapenem-resistant A. baumannii (7), K. pneumoniae (6) or P. aeruginosa (3), polymyxin B MIC 0.5 mg/L. Pneumonia 9/11, bloodstream infection 4/11. Median APACHE II 14, SOFA 10. ECMO duration median 16.18 days.",
    dose_range = "Loading dose 200-300 mg on day 1, then 100-150 mg (2.5-3 mg/kg/day) at the treating physician's discretion, each dose in 50 mL saline as a 1-h infusion. Sampling at steady state (after more than four doses) at 1, 1.5, 2, 4, 8 and 12 h after the start of infusion, pre- and post-oxygenator.",
    regions = "China (First Affiliated Hospital of Zhejiang University School of Medicine; single centre)",
    renal_function = "8/11 on CVVHD for anuria. Table 1 CrCL median (IQR) 79.70 (33.33, 94.89) mL/min; serum creatinine 137.27 (71.00, 184.50), labelled mg/dL but clearly umol/L.",
    ecmo = "Rotaflow centrifugal pump with polymethylpentene oxygenator (MAQUET), PVC tubing, 600 mL saline prime. Pump speed median 3500 (3000, 3770) RPM; blood flow 3.60 (3.00, 3.65) L/min.",
    notes = "Prospective study, January-December 2021. 132 plasma samples (66 pre- and 66 post-oxygenator; no significant difference between sites, Figure 1). Estimation: Phoenix NLME 8.1, FOCE; pcVPC with 1000 replicates and 1000-sample bootstrap."
  )

  ini({
    # Structural parameters: Xu 2024 Table 2, 'Final Model Estimate, Mean'
    # column. The reference patient receives CRRT and has an ECMO blood flow
    # of 3.48 L/min (final-model CL and V equations).
    lcl <- log(3.75); label("Clearance at Q_ECMO = 3.48 L/min (L/h)") # Table 2: tvCL = 3.75 L/h (CV 6.98%; bootstrap median 3.74, 95% CI 1.59-4.34)
    lvc <- log(20.41); label("Central volume of distribution with CRRT (L)") # Table 2: tvV = 20.41 L (CV 17.1%; bootstrap median 20.47, 95% CI 12.82-27.89)
    lvp <- log(9.86); label("Peripheral volume of distribution (L)") # Table 2: tvVp = 9.86 L (CV 53.4%; bootstrap median 10.34, 95% CI 4.06-24.34)
    lq <- log(3.82); label("Intercompartmental clearance (L/h)") # Table 2: tvQ = 3.82 L/h (CV 42.1%; bootstrap median 4.19, 95% CI 1.83-10.49)

    # Covariate effects (final-model equations printed in the Results)
    e_q_ecmo_cl <- -0.321; label("Fractional change in CL per L/min of ECMO blood flow above 3.48 L/min (1/(L/min))") # Table 2: dCLdspeedECMO = -0.321 (bootstrap median -0.315, 95% CI -0.710 to 0.00120); CL equation: tvCL * (1 + (speadECMO - 3.48 L/min) * dCLdspeadECMO)
    e_nocrrt_vc <- -0.525; label("Fractional change in central volume without CRRT (unitless)") # Table 2: dVdCRRT0 = -0.525 (bootstrap median -0.553, 95% CI -0.769 to -0.207); V equation: tvV * (1 + dVdCRRTO * (if no CRRT))

    # Between-subject variability: Table 2 rows 'omega2 V' and 'omega2 CL',
    # footnoted as variances; exponential IIV per Methods. No covariance
    # term is reported.
    etalcl ~ 0.038 # Table 2: omega2 CL = 0.038 (bootstrap median 0.049)
    etalvc ~ 0.062 # Table 2: omega2 V = 0.062 (bootstrap median 0.058)

    # Residual variability: Table 2 'stdev0' = 0.361. The Results select 'the
    # two-compartment model with an additive option', and the Figure 3 pcVPC
    # lower prediction band dips below zero, consistent with an additive
    # (not proportional) error in mg/L.
    addSd <- 0.361; label("Additive residual error (mg/L)") # Table 2: stdev0 = 0.361 (CV 12.0%; bootstrap median 0.351, 95% CI 0.256-0.447)
  })
  model({
    # Individual parameters (Xu 2024 final-model equations)
    cl <- exp(lcl + etalcl) * (1 + e_q_ecmo_cl * (Q_ECMO - 3.48))
    vc <- exp(lvc + etalvc) * (1 + e_nocrrt_vc * (1 - RRT_CRRT_STATUS))
    vp <- exp(lvp)
    q <- exp(lq)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Two-compartment disposition; polymyxin B given as a 1-h IV infusion
    # into the central compartment. Dose in mg and volume in L give mg/L.
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    Cc <- central / vc
    Cc ~ add(addSd)
  })
}
