Lee_2024_topiramate <- function() {
  description <- "One-compartment population PK model with first-order absorption and elimination for topiramate (TPM) in Korean adults with epilepsy undergoing routine therapeutic drug monitoring (Lee 2024). Apparent clearance CL/F carries an additive enzyme-inducing antiseizure-medication term (separate coefficients in L/h for phenytoin, carbamazepine, oxcarbazepine and phenobarbital), a power effect of creatinine clearance normalized to 90 mL/min, and a power effect of the topiramate daily dose normalized to 100 mg/day; apparent volume Vd/F scales allometrically with body weight (exponent fixed at 1) centred at 62 kg. Absorption rate constant fixed at 2 /h. Inter-individual variability on CL/F only (31.0 %CV) with a proportional residual error of 27.8%. The structural model was carried from Bae 2016 and all parameters were re-estimated on the present therapeutic-drug-monitoring dataset."
  reference <- paste(
    "Lee S, Kim HC, Jang Y, Lee HS, Ahn SJ, Lee ST, Jung KH, Park KI, Jung KY,",
    "Oh J, Lee S, Yu KS, Jang IJ, Lee S, Chu K, Lee SK.",
    "Topiramate dosage optimization for effective antiseizure management via",
    "population pharmacokinetic modeling.",
    "Ann Clin Transl Neurol. 2024;11(2):424-435. doi:10.1002/acn3.51962.",
    "Structural model carried from Bae EK, Lee J, Shin JW, et al. Factors",
    "influencing topiramate clearance in adult patients with epilepsy: a",
    "population pharmacokinetic analysis. Seizure. 2016;37:8-12;",
    "all parameters re-estimated on the Lee 2024 dataset (Lee 2024 Table S1).",
    sep = " "
  )
  vignette <- "Lee_2024_topiramate"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric scaling on Vd/F with the exponent fixed at 1, centred at 62 kg (Lee 2024 Table S1 equation Vd/F = theta2 * (WT/62); Results 'allometric scaling was applied with the exponent of 1'). Body weight was documented for only 237 of the 555 samples (Lee 2024 Table 1 footnote a); cohort mean 64.6 +/- 17.7 kg.",
      source_name        = "WT"
    ),
    CRCL = list(
      description        = "Creatinine clearance",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on CL/F normalized to 90 mL/min (Lee 2024 Table S1 theta8 = 0.277; equation (CLcr/90)^theta8). Lee 2024 does not state the creatinine-clearance assay or whether the value was body-surface-area normalized, so the raw mL/min form is assumed here (as in Delattre_2010_amikacin.R, Georges_2009_ceftazidime.R and Chen_2023_nemonoxacin.R); 90 mL/min is the conventional normal-renal-function anchor.",
      source_name        = "CLcr"
    ),
    CONMED_PHT = list(
      description        = "Concomitant phenytoin coadministration indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant phenytoin)",
      notes              = "Additive effect on CL/F of +1.02 L/h (Lee 2024 Table S1 theta6), i.e. +70% relative to the 1.45 L/h monotherapy intercept (Lee 2024 Discussion). Time-fixed here because the analysis required a stable antiseizure-medication regimen for at least one month before sampling (Lee 2024 Methods).",
      source_name        = "PHT"
    ),
    CONMED_CBZ = list(
      description        = "Concomitant carbamazepine coadministration indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant carbamazepine)",
      notes              = "Additive effect on CL/F of +0.703 L/h (Lee 2024 Table S1 theta4), i.e. +48% relative to the 1.45 L/h monotherapy intercept (Lee 2024 Discussion). Time-fixed here because the analysis required a stable antiseizure-medication regimen for at least one month before sampling (Lee 2024 Methods).",
      source_name        = "CBZ"
    ),
    CONMED_OXC = list(
      description        = "Concomitant oxcarbazepine coadministration indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant oxcarbazepine)",
      notes              = "Additive effect on CL/F of +0.419 L/h (Lee 2024 Table S1 theta5), i.e. +29% relative to the 1.45 L/h monotherapy intercept (Lee 2024 Discussion). Time-fixed here because the analysis required a stable antiseizure-medication regimen for at least one month before sampling (Lee 2024 Methods).",
      source_name        = "OXC"
    ),
    CONMED_PB = list(
      description        = "Concomitant phenobarbital coadministration indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant phenobarbital)",
      notes              = "Additive effect on CL/F of +0.376 L/h (Lee 2024 Table S1 theta7), i.e. +26% relative to the 1.45 L/h monotherapy intercept (Lee 2024 Discussion). Lee 2024 does not state whether primidone was pooled with phenobarbital.",
      source_name        = "PB"
    ),
    DOSE_TPM_MGD = list(
      description        = "Total daily topiramate dose",
      units              = "mg/day",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on CL/F normalized to 100 mg/day (Lee 2024 Table S1 theta9 = 0.193; equation (DOSE/100)^theta9), encoding the mild dose-dependence of topiramate apparent clearance. Per-record covariate carrying the current total daily dose level, updated when the prescriber changes the dose; distinct from the rxode2 event-table amt column, which carries the amount of each individual administration. For a twice-daily regimen the value is the sum across the day (89.2 mg b.i.d. -> DOSE_TPM_MGD = 178.4). Cohort mean 178.4 +/- 117.9 mg/day (Lee 2024 Table 1).",
      source_name        = "DOSE"
    )
  )

  compartmentData <- list(
    depot   = list(analyte = "topiramate", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "topiramate", units = "mg", specimen = "serum", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 389L,
    n_studies      = 1L,
    n_observations = 555L,
    age_mean       = "46.4 +/- 10.9 years",
    weight_mean    = "64.6 +/- 17.7 kg (documented for 237 of 555 samples)",
    sex_female_pct = 49.5,
    race_ethnicity = c(Korean = 100),
    disease_state  = "Adults with epilepsy on chronic topiramate therapy (53 patients on topiramate monotherapy, 340 on antiseizure-medication polytherapy; 55.3% of samples drawn while co-prescribed an enzyme-inducing antiseizure medication - carbamazepine, oxcarbazepine, phenytoin or phenobarbital). 5.6% of samples were taken during an insufficient response (< 50% seizure-frequency reduction), 94.4% during a sufficient response of which 57.7% were seizure-free.",
    dose_range     = "Mean total daily dose 178.4 +/- 117.9 mg/day (121.4 +/- 59.4 mg/day monotherapy; 187.5 +/- 122.3 mg/day polytherapy). Lee 2024 does not state the dosing frequency.",
    regions        = "South Korea (Seoul National University Hospital)",
    notes          = "Retrospective therapeutic-drug-monitoring cohort sampled between January 2017 and January 2022 (Lee 2024 Table 1). Regimens were stable for at least one month before sampling so steady state could be assumed, and the sampling time relative to the last dose was recorded for every sample. Inpatient measurements were excluded. Samples below the 0.25 mg/L lower limit of quantification were discarded (M1 method). Serum topiramate was measured by UPLC-MS/MS with a 0.25-50 mg/L calibration range."
  )

  ini({
    # ------------------------------------------------------------------
    # Structural parameters (Lee 2024 Table S1, "Estimation results").
    #   CL/F = (theta1 + theta4*CBZ + theta5*OXC + theta6*PHT + theta7*PB)
    #          * (CLcr/90)^theta8 * (DOSE/100)^theta9
    #   Vd/F = theta2 * (WT/62)
    # The published equation is confirmed by the paper's own Discussion,
    # which reports the comedication effects as +70% (PHT), +48% (CBZ),
    # +29% (OXC) and +26% (PB) of the 1.45 L/h intercept -- reproduced
    # exactly by theta6/theta1, theta4/theta1, theta5/theta1, theta7/theta1.
    # CL/F and Vd/F are apparent (oral, F not identifiable), so no
    # separate bioavailability parameter is estimated.
    # ------------------------------------------------------------------
    lcl <- log(1.45); label("Apparent clearance CL/F intercept, topiramate monotherapy at CLcr 90 mL/min and 100 mg/day (L/h)")  # Lee 2024 Table S1 theta1 = 1.45 (2.80% RSE)
    lvc <- log(78.1); label("Apparent volume of distribution Vd/F at the 62 kg reference weight (L)")                            # Lee 2024 Table S1 theta2 = 78.1 (22.0% RSE)
    lka <- fixed(log(2)); label("First-order absorption rate constant Ka (1/h)")                                                 # Lee 2024 Table S1 "Ka ... Fixed at 2"

    # ------------------------------------------------------------------
    # Covariate effects on CL/F. The four enzyme-inducing antiseizure
    # medications enter ADDITIVELY in L/h (not multiplicatively), so each
    # coefficient is an absolute clearance increment on the 1.45 L/h
    # monotherapy intercept before the CLcr and DOSE power terms apply.
    # ------------------------------------------------------------------
    e_pht_cl <- 1.02;  label("Additive effect of concomitant phenytoin on CL/F (L/h)")       # Lee 2024 Table S1 theta6 = 1.02 (24.8% RSE)
    e_cbz_cl <- 0.703; label("Additive effect of concomitant carbamazepine on CL/F (L/h)")   # Lee 2024 Table S1 theta4 = 0.703 (15.4% RSE)
    e_oxc_cl <- 0.419; label("Additive effect of concomitant oxcarbazepine on CL/F (L/h)")   # Lee 2024 Table S1 theta5 = 0.419 (21.9% RSE)
    e_pb_cl  <- 0.376; label("Additive effect of concomitant phenobarbital on CL/F (L/h)")   # Lee 2024 Table S1 theta7 = 0.376 (34.8% RSE)

    # Power exponents on CL/F.
    e_crcl_cl     <- 0.277; label("Power exponent of creatinine clearance on CL/F, normalized to 90 mL/min (unitless)")  # Lee 2024 Table S1 theta8 = 0.277 (42.6% RSE)
    e_dose_tpm_cl <- 0.193; label("Power exponent of total daily topiramate dose on CL/F, normalized to 100 mg/day (unitless)")  # Lee 2024 Table S1 theta9 = 0.193 (16.7% RSE)

    # Allometric exponent on Vd/F. Reported without uncertainty and stated
    # in the Results as applied with "the exponent of 1", i.e. held fixed.
    e_wt_vc <- fixed(1); label("Allometric exponent of body weight on Vd/F (unitless)")  # Lee 2024 Table S1 equation Vd/F = theta2*(WT/62); Results "allometric scaling was applied with the exponent of 1"

    # ------------------------------------------------------------------
    # Inter-individual variability. Lee 2024 estimated IIV on CL/F only
    # (no IIV on Vd/F or Ka was reported), applied to the fully
    # covariate-adjusted CL/F. Reported as 31.0 %CV; converted to the
    # internal log-scale variance as omega^2 = log(1 + CV^2)
    # = log(1 + 0.310^2) = 0.0917584.
    # ------------------------------------------------------------------
    etalcl ~ 0.0917584  # CV = 31.0% -> log(1 + 0.310^2); Lee 2024 Table S1 omega_CL/F = 31.0 (14.1% RSE)

    # ------------------------------------------------------------------
    # Residual error: proportional only, reported as 27.8%.
    # ------------------------------------------------------------------
    propSd <- 0.278; label("Proportional residual error on serum topiramate (fraction)")  # Lee 2024 Table S1 sigma_prop = 27.8% (6.20% RSE)
  })

  model({
    # 1. Apparent clearance. The comedication increments are additive in
    # L/h inside the parentheses; the CLcr and daily-dose power terms
    # multiply the whole intercept-plus-comedication sum; exponential IIV
    # multiplies the fully covariate-adjusted value (Lee 2024 Table S1).
    cl <-
      (exp(lcl) +
         e_pht_cl * CONMED_PHT +
         e_cbz_cl * CONMED_CBZ +
         e_oxc_cl * CONMED_OXC +
         e_pb_cl  * CONMED_PB) *
      (CRCL / 90)^e_crcl_cl *
      (DOSE_TPM_MGD / 100)^e_dose_tpm_cl *
      exp(etalcl)

    # 2. Apparent volume of distribution, allometric on body weight.
    vc <- exp(lvc) * (WT / 62)^e_wt_vc

    ka <- exp(lka)

    # 3. Micro-constant
    kel <- cl / vc

    # 4. ODE system. Oral doses enter depot (mg); central carries mg.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # 5. Serum topiramate concentration in mg/L (dose in mg, volume in L)
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
