MoralesJunior_2025_cefepime <- function() {
  description <- paste(
    "Two-compartment population PK model with first-order elimination for",
    "intravenously infused cefepime in critically ill children and young",
    "adults (1 month to 30 years) admitted to a pediatric intensive care",
    "unit (Morales Junior 2025; 100 patients, 510 opportunistically",
    "scavenged plasma concentrations). Allometric body-weight scaling is",
    "applied with exponents fixed at 0.75 on the clearance parameters and 1",
    "on the volumes, standardised to 70 kg. Clearance is 6.38 L/h at 70 kg",
    "and eGFR 147.6 mL/min/1.73 m^2 and scales with a power of 0.66 on",
    "BSA-normalized eGFR; the central volume is 15 L at 70 kg and expands",
    "exponentially with the cumulative percentage of fluid balance",
    "(exp(0.026 * CUM_FLUID_BAL_PCT)). Intercompartmental clearance is 3.65",
    "L/h and the peripheral volume 8.91 L at 70 kg, both without covariates",
    "and without interindividual variability (the sparse opportunistic",
    "sampling did not support random effects on them). Proportional",
    "residual error only. Patients receiving renal replacement therapy or",
    "ECMO were excluded and no neonates were studied. Externally validated",
    "against an independent 41-patient PICU cohort.",
    sep = " "
  )
  reference <- paste(
    "Morales Junior R, Hambrick HR, Mizuno T, Pavia KE, Paice KM, Tang P,",
    "Schuler E, Krallman KA, Johnson L, Collins M, Gibson A, Curry C,",
    "Kaplan J, Goldstein S, Tang Girdwood S (2025). Population",
    "Pharmacokinetics of Cefepime in Critically Ill Children and Young",
    "Adults: Model Development and External Validation for Monte Carlo",
    "Simulations and Model-Informed Precision Dosing. Clinical",
    "Pharmacokinetics 64(4):553-564. doi:10.1007/s40262-025-01485-5.",
    "PMID 39987410; PMC12041147. Final parameter estimates and the printed",
    "model equations are from Table 2.",
    sep = " "
  )
  vignette <- "MoralesJunior_2025_cefepime"
  units    <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Cefepime was given as an intravenous infusion and
  # TOTAL (not unbound) plasma cefepime was quantified by HPLC
  # (Morales Junior 2025 Sect. 2.4); the free fraction of 0.80 used for the
  # fT>MIC targets is applied outside the PK model (see
  # population$protein_binding).
  compartmentData <- list(
    central     = list(analyte = "cefepime", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "cefepime", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometric scaling with exponents FIXED at 0.75 for the clearance",
        "parameters (CL, Q) and 1 for the volume parameters (V1, V2),",
        "standardised to a typical 70 kg adult (Morales Junior 2025",
        "Sect. 2.6: 'employing a power function with a fixed exponent of",
        "0.75 for clearance parameters and 1 for volume parameters, and",
        "scaled to a typical adult weighing 70 kg'). The same reference",
        "weight is printed in every Table 2 equation. Cohort median 24.8 kg",
        "(IQR 11.9-53.6, Table 1); the reference weight is therefore far",
        "outside the studied range and is a normalisation constant, not a",
        "typical patient. Body weight was collected daily (ESM 1) so the",
        "covariate can be carried time-varying; note that the",
        "CUM_FLUID_BAL_PCT denominator is the PICU ADMISSION body weight,",
        "which is a different quantity from a time-varying WT.",
        sep = " "
      ),
      source_name        = "WT"
    ),
    CRCL = list(
      description        = "BSA-normalized estimated glomerular filtration rate (bedside Schwartz equation for age < 18 years, CKD-EPI for age >= 18 years)",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "BSA-normalized eGFR, NOT a raw Cockcroft-Gault creatinine",
        "clearance. Two different estimating equations are used within the",
        "one column, split by age: the bedside Schwartz equation for",
        "patients aged < 18 years and the CKD-EPI equation for those aged",
        ">= 18 years (Morales Junior 2025 Sect. 2.5 and the Table 1",
        "footnote). Enters clearance as the power term",
        "(CRCL / 147.6)^0.66 (Table 2 fixed-effects equation). The 147.6",
        "mL/min/1.73 m^2 normalising constant is printed only inside that",
        "equation; the paper does not state how it was chosen and it is not",
        "the cohort median (Table 1 median 128.3, IQR 91-171.6). Renal",
        "function spans kidney impairment through augmented renal",
        "clearance: at cefepime initiation 39% of patients had normal renal",
        "function, 41% augmented renal clearance and 20% kidney impairment",
        "(Sect. 3.1). Missing values on a study day were imputed by last",
        "observation carried backward, or by the population median when",
        "absent for the whole stay (Sect. 2.5), so the column is naturally",
        "time-varying. Must be strictly positive -- it enters a power term.",
        "Patients on any renal replacement modality (intermittent dialysis,",
        "CRRT, peritoneal dialysis) or ECMO were EXCLUDED (Sect. 2.2), so",
        "this covariate is never paired with a renal-replacement indicator",
        "in this model and the model must not be used for those patients.",
        sep = " "
      ),
      source_name        = "eGFR"
    ),
    CUM_FLUID_BAL_PCT = list(
      description        = "Cumulative fluid balance expressed as a percentage of PICU admission body weight",
      units              = "%",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Defined in Morales Junior 2025 Sect. 2.5 as 'the sum of all the",
        "previous days' and present day's percentage of fluid balance (sum",
        "of each day's net fluid balance / PICU admission body weight x",
        "100)', where the daily net fluid balance is the daily difference",
        "between all recorded intakes and all recorded outputs. It is",
        "therefore a running SUM over study days, normalised by the",
        "ADMISSION body weight (a fixed denominator), and it is naturally",
        "time-varying within subject and can be negative when a patient is",
        "net-negative. Enters the central volume as the EXPONENTIAL",
        "multiplier exp(0.026 * CUM_FLUID_BAL_PCT) (Table 2 V1 equation:",
        "'V1 = V1pop x (WT/70) x e^(beta x Cum%FB) x e^(etaV1)'), so a +10%",
        "cumulative fluid balance expands V1 by 30%. Cohort day-by-day",
        "medians run 3.3, 4.5, 5.5, 6.2, 4.9, 5.7 and 4.0% on study days",
        "1-7 (Table 1), with IQRs reaching -0.1 to 15.1%. The Monte Carlo",
        "simulations in Sect. 2.8 set this covariate to 0% for every",
        "simulated patient, and the Discussion states that 'if fluid",
        "balance data are unavailable, the model can still be used by",
        "assuming a fluid balance of 0, although this could reduce the",
        "accuracy of the volume of distribution estimation'.",
        sep = " "
      ),
      source_name        = "Cum%FB"
    )
  )

  # Screened covariates that were NOT retained in the final model. These are
  # documentation only -- none is referenced in model().
  covariatesDataExcluded <- list(
    SBP = list(
      description        = "Systolic blood pressure",
      units              = "mmHg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Significant on clearance during forward inclusion (dOFV > 3.84) but removed during backward elimination (Morales Junior 2025 Sect. 3.2). Collected daily as the lowest and highest recorded systolic pressure (ESM 1 Supplementary Table 1). No point estimate is published."
    ),
    HR = list(
      description        = "Heart rate",
      units              = "beats/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Significant on clearance during forward inclusion but removed during backward elimination (Morales Junior 2025 Sect. 3.2). Collected daily as the lowest and highest recorded heart rate (ESM 1 Supplementary Table 1). No point estimate is published."
    ),
    PRISM3 = list(
      description        = "Pediatric Risk of Mortality III (PRISM-III) severity-of-illness score",
      units              = "(score)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Significant on clearance during forward inclusion but removed during backward elimination (Morales Junior 2025 Sect. 3.2). ESM 1 Supplementary Table 1 records PRISM III alongside PIM2 and PIM3 as the pediatric mortality-risk scores collected. No point estimate is published."
    ),
    PMA = list(
      description        = "Postmenstrual age (gestational age at birth plus postnatal age)",
      units              = "weeks",
      type               = "continuous",
      reference_category = NULL,
      notes              = "A Hill-function maturation factor on renal clearance was tested (Morales Junior 2025 Sect. 2.6, citing Rhodin 2009) and is absent from the final Table 2 model. The Discussion attributes this to the absence of neonates from the cohort ('The absence of data from neonates restricted the applicability of the model to this age group, who have distinct maturation-related pharmacokinetic profiles'). No maturation parameters are published."
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 100L,
    n_studies        = 1L,
    n_sites          = 1L,
    n_concentrations = 510L,
    age_range        = "1 month to 30 years (protocol inclusion range)",
    age_median       = "7.6 years (IQR 1.6-16; Table 1)",
    weight_median    = "24.8 kg (IQR 11.9-53.6; Table 1)",
    sex_female_pct   = 45,
    race_ethnicity   = "Not reported. Table 1 tabulates only age band, sex, body weight, serum creatinine, eGFR, serum albumin, mechanical ventilation, vasopressor treatment and daily cumulative percentage of fluid balance.",
    disease_state    = paste(
      "Critically ill children and young adults admitted to the pediatric",
      "intensive care unit at Cincinnati Children's Hospital Medical Center",
      "who received at least one dose of cefepime and had at least one total",
      "cefepime concentration measured. Age bands: 27 infants (1 month to",
      "< 2 years, 27%), 32 children (2 to < 12 years, 32%), 23 adolescents",
      "(12 to < 18 years, 23%) and 18 young adults (18 to < 30 years, 18%).",
      "Illness severity is reflected in 45% receiving mechanical ventilation",
      "and 41% vasopressor treatment at some point during follow-up.",
      "Baseline (study day 1) serum creatinine 0.38 mg/dL (IQR 0.23-0.63)",
      "and serum albumin 2.9 g/dL (IQR 2.4-3.3) (Table 1). Patients",
      "receiving intermittent dialysis, continuous renal replacement",
      "therapy, peritoneal dialysis or ECMO were EXCLUDED.",
      sep = " "
    ),
    renal_function   = "At cefepime initiation 39 patients (39%) had normal renal function, 41 (41%) augmented renal clearance and 20 (20%) kidney impairment (Sect. 3.1). Cohort eGFR median 128.3 mL/min/1.73 m^2 (IQR 91-171.6). Renal replacement therapy of any modality was an exclusion criterion, so the model carries no dialysis clearance arm.",
    fluid_status     = "Cumulative percentage of fluid balance, medians by study day: 3.3, 4.5, 5.5, 6.2, 4.9, 5.7, 4.0% on days 1-7 (Table 1).",
    dose_range       = "Real-world clinician-directed dosing. The institutional standard is 50 mg/kg per dose (maximum 2000 mg per dose) every 8 h as a 30-min infusion in patients with preserved kidney function (Sect. 2.3). Monte Carlo simulations explored 50 mg/kg (maximum 2000 mg) every 6, 8, 12 or 24 h over 30 min or 3 h, and continuous infusions of 30-180 mg/kg per day (maximum 6000 mg/day).",
    protein_binding  = "Not fitted. TOTAL plasma cefepime was quantified; the free concentrations used for the fT>MIC targets were obtained by assuming a fixed 20% protein binding, i.e. an unbound fraction of 0.80 (Sect. 2.8). The model therefore outputs TOTAL plasma cefepime as Cc; multiply by 0.80 for free drug.",
    sampling         = "Scavenged opportunistic residual plasma sampling during the first 7 days of beta-lactam therapy; samples drawn during a cefepime infusion were excluded. Median 4 samples per patient (range 1-22). Assay linear 0.5-200 ug/mL by HPLC, within- and between-day CV below 15%. Three of 510 samples (0.6%) were below the limit of quantification and were imputed at half the LLOQ (Beal M6).",
    regions          = "United States (single centre, Cincinnati, Ohio)",
    notes            = paste(
      "Model-building data were collected October 2018 - November 2021",
      "under IRB #2018-3245. The model was externally validated against a",
      "separate prospective PICU cohort (November 2022 - August 2024; 41",
      "patients, 234 concentrations, median age 11.9 years, median weight",
      "28.8 kg, median eGFR 62.6 mL/min/1.73 m^2), giving a population-level",
      "median prediction error (MDPE) of 6.9% with median absolute",
      "prediction error (MDAPE) 34.5%, and individual-level MDPE -1.8% with",
      "MDAPE 18.3% (Sect. 3.4). Note a reporting inconsistency in the",
      "source: Sect. 2.6 states the fit was done in NONMEM 7.5 with FOCE-I,",
      "while the Table 2 column header reads 'Stochastic approximation',",
      "the residual-error parameter is named 'b' and Sect. 2.8 uses Simulx",
      "-- all Monolix/SAEM conventions. The discrepancy is in the software",
      "description only; it does not affect any reported value or its",
      "interpretation (see the vignette Assumptions and deviations).",
      sep = " "
    )
  )

  ini({
    # -----------------------------------------------------------------
    # Structural parameters. All values are the "Estimate" column of
    # Morales Junior 2025 Table 2; the bootstrap columns (median and
    # 2.5th / 97.5th percentiles of 1000 replicates) are parameter
    # precision and are deliberately not carried into ini().
    # -----------------------------------------------------------------

    lcl <- log(6.38)
    label("Clearance at WT = 70 kg and eGFR = 147.6 mL/min/1.73 m^2 (L/h)")
    # Morales Junior 2025 Table 2, "CL (L/h/70kg^0.75)" = 6.38 (RSE 5%;
    # bootstrap median 6.41, 95% CI 5.85-7.04). Repeated in Sect. 3.2:
    # "Clearance was 6.38 L/h/70 kg^0.75 with an IIV of 38.1%".

    e_wt_cl_q <- fixed(0.75)
    label("Allometric exponent on (WT / 70) shared by CL and Q (unitless)")
    # Morales Junior 2025 Sect. 2.6: "employing a power function with a
    # fixed exponent of 0.75 for clearance parameters and 1 for volume
    # parameters, and scaled to a typical adult weighing 70 kg".
    # Reported without RSE in Table 2 because it was not estimated; the
    # exponent is printed inside both the CL and the Q equations.

    e_crcl_cl <- 0.66
    label("Power exponent on (CRCL / 147.6) for CL (unitless)")
    # Morales Junior 2025 Table 2, "beta_eGFR" = 0.66 (RSE 11%; bootstrap
    # median 0.66, 95% CI 0.51-0.79). Enters as the printed Table 2
    # equation CL = CLpop * (WT/70)^0.75 * (eGFR/147.6)^beta * e^(etaCl).
    # Estimated, not fixed: Sect. 3.2 records that removing eGFR from
    # clearance in backward elimination raised the OFV by 132.3.

    lvc <- log(15)
    label("Central volume of distribution at WT = 70 kg and zero cumulative fluid balance (L)")
    # Morales Junior 2025 Table 2, "V1 (L/70kg)" = 15 (RSE 28%; bootstrap
    # median 16.41, 95% CI 10.74-20.39). Repeated in Sect. 3.2 and in the
    # Discussion ("the central volume of distribution [was] ... 15 L/70 kg").

    e_wt_vc_vp <- fixed(1)
    label("Allometric exponent on (WT / 70) shared by V1 and V2 (unitless)")
    # Morales Junior 2025 Sect. 2.6 (same sentence as e_wt_cl_q): the
    # volume exponent was fixed at 1. Table 2 prints the volume equations
    # with a bare (WT/70) ratio and no exponent, which is the exponent-1
    # form.

    e_cum_fluid_bal_pct_vc <- 0.026
    label("Exponential coefficient of cumulative percentage of fluid balance on V1 (per %)")
    # Morales Junior 2025 Table 2, "beta_Cum%FB" = 0.026 (RSE 44%;
    # bootstrap median 0.022, 95% CI 0.002-0.041). Enters EXPONENTIALLY,
    # not proportionally, per the printed Table 2 equation
    # V1 = V1pop * (WT/70) * e^(beta * Cum%FB) * e^(etaV1). Estimated, not
    # fixed: Sect. 3.2 records that removing it in backward elimination
    # raised the OFV by 14.1. The RSE of 44% exceeds the paper's stated
    # 40% retention criterion (Sect. 2.6); the authors nonetheless retained
    # it on the dOFV criterion, and the bootstrap 95% CI (0.002-0.041)
    # excludes zero.

    lq <- log(3.65)
    label("Intercompartmental clearance at WT = 70 kg (L/h)")
    # Morales Junior 2025 Table 2, "Q (L/h/70kg^0.75)" = 3.65 (RSE 54%;
    # bootstrap median 3.09, 95% CI 1.12-6.70). Repeated in Sect. 3.2.
    # No covariate other than allometric weight (Table 2 equation
    # Q = Qpop * (WT/70)^0.75).

    lvp <- log(8.91)
    label("Peripheral volume of distribution at WT = 70 kg (L)")
    # Morales Junior 2025 Table 2, "V2 (L/70kg)" = 8.91 (RSE 21%;
    # bootstrap median 8.92, 95% CI 6.55-12.98). Repeated in Sect. 3.2.
    # No covariate other than allometric weight (Table 2 equation
    # V2 = V2pop * (WT/70)).

    # -----------------------------------------------------------------
    # Interindividual variability. Table 2 reports IIV as a coefficient
    # of variation, and its footnote gives the transform explicitly:
    # "The IIV is expressed as coefficient of variation (%) calculated as
    # sqrt(e^(omega^2) - 1) x 100, where omega^2 corresponds to the
    # variance of the random effects". Inverting,
    # omega^2 = log(1 + (CV/100)^2), which is the log-normal variance
    # rxode2 wants:
    #   CL: log(1 + 0.381^2) = 0.1355452
    #   V1: log(1 + 0.149^2) = 0.0219582
    # Random effects on Q and V2 are ABSENT by design, not unreported:
    # Sect. 3.2 states "The data did not support the inclusion of random
    # effects on intercompartmental clearance and peripheral volume of
    # distribution, likely because of the sparse sampling approach". No
    # eta correlations are reported, so the two etas are independent.
    # -----------------------------------------------------------------

    etalcl ~ 0.1355452
    # Morales Junior 2025 Table 2, "IIV CL" = 38.1% CV (RSE 10%; eta
    # shrinkage 4.2%; bootstrap median 37.9%, 95% CI 29.6-45.3%).

    etalvc ~ 0.0219582
    # Morales Junior 2025 Table 2, "IIV V1" = 14.9% CV (RSE 55%; eta
    # shrinkage 74.5%; bootstrap median 16.4%, 95% CI 3.9-28.8%). The
    # Discussion warns that the 74.5% shrinkage means the V1 individual
    # estimates and their diagnostics should be interpreted with caution.

    propSd <- 0.319
    label("Proportional residual error (fraction)")
    # Morales Junior 2025 Table 2, "Error model parameter: proportional
    # only", row "b" = 31.9% (RSE 2%; bootstrap median 31.4%, 95% CI
    # 27.4-35.7%). Sect. 3.2: "The residual variability was best described
    # by the proportional error model". A combined and an additive error
    # model were both tested and rejected (Sect. 2.6), so there is no
    # additive term to encode.
  })

  model({
    # Normalising constants, exactly as printed in the Morales Junior 2025
    # Table 2 fixed-effects equations.
    wt_ref   <- 70    # kg
    crcl_ref <- 147.6 # mL/min/1.73 m^2

    # Table 2: CL = CLpop * (WT/70)^0.75 * (eGFR/147.6)^beta * e^(etaCl)
    cl <- exp(lcl + etalcl) * (WT / wt_ref)^e_wt_cl_q * (CRCL / crcl_ref)^e_crcl_cl

    # Table 2: V1 = V1pop * (WT/70) * e^(beta * Cum%FB) * e^(etaV1)
    # The fluid-balance term is EXPONENTIAL, not proportional.
    vc <-
      exp(lvc + etalvc) *
      (WT / wt_ref)^e_wt_vc_vp *
      exp(e_cum_fluid_bal_pct_vc * CUM_FLUID_BAL_PCT)

    # Table 2: Q = Qpop * (WT/70)^0.75 and V2 = V2pop * (WT/70).
    # Both share the weight exponents with CL and V1 respectively and
    # carry no other covariate and no random effect.
    q  <- exp(lq)  * (WT / wt_ref)^e_wt_cl_q
    vp <- exp(lvp) * (WT / wt_ref)^e_wt_vc_vp

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(central)     <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <-          k12 * central - k21 * peripheral1

    # TOTAL plasma cefepime, which is what the assay measured. The paper's
    # fT>MIC targets are evaluated on free drug, obtained by multiplying by
    # the assumed unbound fraction 0.80 (20% protein binding); that is a
    # literature constant applied outside the PK model and is recorded in
    # population$protein_binding rather than in ini().
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
