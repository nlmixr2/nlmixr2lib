Reeder_2025_tazobactam <- function() {
  description <- paste(
    "One-compartment population PK model with linear elimination for the",
    "tazobactam component of intravenously infused piperacillin-tazobactam",
    "in critically ill adults in the ICU (Reeder 2025; 108 patients, 229",
    "plasma concentrations from opportunistic sampling). Clearance is the sum",
    "of a renal arm proportional to Cockcroft-Gault creatinine clearance on",
    "lean body weight and a constant non-renal arm; in patients on continuous",
    "renal replacement therapy (CRRT) an estimated CRRT clearance replaces",
    "the renal arm while the non-renal arm is retained. Uncorrelated",
    "interindividual variability on total clearance and volume and a",
    "combined proportional plus additive residual error. The companion",
    "piperacillin model is Reeder_2025_piperacillin.",
    sep = " "
  )
  reference <- paste(
    "Reeder JA, Creech CB, Nation RL, Gu K, Nalbant D, Wu N,",
    "Jimenez-Truque N, Fissell W, Rolsma S, Fishbane N, Kirkpatrick C, Patel",
    "P, Watanabe A, Landersdorfer C, Winokur P, An G (2025). Utilizing an",
    "Opportunistic Clinical Study and Population-Based Pharmacokinetic Models",
    "to Identify Rational Empiric Dosing Regimens for Piperacillin-Tazobactam",
    "in Critically Ill Patients. The Journal of Clinical Pharmacology",
    "65(4):452-465. doi:10.1002/jcph.6161.",
    sep = " "
  )
  vignette <- "Reeder_2025_piperacillin_tazobactam"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "tazobactam", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = "Cockcroft-Gault creatinine clearance computed on lean body weight (raw, not BSA-normalized)",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Source name CLCR,LBW. Methods 'Population PK Model Development':",
        "'estimated creatinine clearance (CLcr) (both CLcr,TBW and CLcr,LBW",
        "calculated using the Cockcroft-Gault formula)'; Table 1 footnote:",
        "'CLcr,LBW, estimated creatinine clearance using lean body weight'.",
        "The value is a raw mL/min clearance and is NOT BSA-normalized,",
        "matching the raw-mL/min variant of the canonical CRCL column used by",
        "the sibling Rolsma_2025_cefepime.R model from the same research",
        "group. Enters the renal clearance arm as the power term",
        "(CRCL / 51)^1, 51 mL/min being the Table 1 cohort median (range 8 to",
        "173). The CLcr,TBW variant (median 82, range 10 to 289 mL/min) was",
        "screened and not retained. The renal arm is switched off when",
        "RRT_CRRT_STATUS = 1, so the CRCL value of a CRRT patient does not",
        "affect the prediction."
      ),
      source_name = "CLCR,LBW"
    ),
    RRT_CRRT_STATUS = list(
      description = "Continuous renal replacement therapy (1 = on CRRT, 0 = not)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no CRRT)",
      notes = paste(
        "Source name CRRT. 10 of 112 patients (8.9%) were on CRRT, 8 on",
        "continuous veno-venous hemodiafiltration and 2 on continuous",
        "veno-venous hemofiltration (Results 'Patient Characteristics').",
        "Table 2 footnote b: total clearance in CRRT patients is '(CLR x (1 -",
        "CRRT) + CLNR + (CLCRRT x CRRT)) = CLNR + CLCRRT', so CRRT replaces",
        "the renal arm only and the non-renal arm is retained. Subject-level",
        "indicator, treated as time-fixed."
      ),
      source_name = "CRRT"
    )
  )

  covariatesDataExcluded <- list(
    LBM = list(
      description = "Lean body weight (formula not stated by the source)",
      units = "kg",
      type = "continuous",
      notes = "Screened (Methods covariate list, 'lean body weight (LBW)') and retained for piperacillin non-renal clearance but NOT for tazobactam (Results: 'For tazobactam, two significant covariates were identified: CLCR,LBW on CLR and CRRT on total clearance'). Cohort median 57 kg, range 30 to 89 (Table 1)."
    ),
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened (Methods covariate list, 'total body weight (TBW)') and not retained. Cohort median 92 kg, range 41 to 160 (Table 1)."
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened (Methods covariate list) and not retained. Cohort median 64 years, range 21 to 88 (Table 1)."
    ),
    SEXF = list(
      description = "Sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      notes = "Screened (Methods covariate list) and not retained. 50 of 112 patients (45%) female (Table 1)."
    ),
    SOFA = list(
      description = "Sequential Organ Failure Assessment score",
      units = "(score)",
      type = "continuous",
      notes = "Screened as the illness-severity covariate and not retained. Printed as 5 (3-7) under the Table 1 'Median (Range)' header."
    ),
    DIS_ARF = list(
      description = "Acute kidney injury (acute renal failure) indicator",
      units = "(binary)",
      type = "binary",
      notes = "Screened and not retained (Discussion: 'these covariates were not significant, as they did not provide any additional significant reduction in IIV'). 29 of 112 patients (26%) (Table 1)."
    ),
    DIS_SEPSIS = list(
      description = "Sepsis indicator",
      units = "(binary)",
      type = "binary",
      notes = "Screened and not retained. 37 of 112 patients (33%) (Table 1)."
    ),
    DIS_SEPTIC_SHOCK = list(
      description = "Septic shock indicator",
      units = "(binary)",
      type = "binary",
      notes = "Screened and not retained. 30 of 112 patients (27%) (Table 1)."
    ),
    MECH_VENT = list(
      description = "Mechanical ventilation indicator",
      units = "(binary)",
      type = "binary",
      notes = "Screened and not retained. 44 of 112 patients (39%) (Table 1)."
    ),
    HEPATIC_FUNCTION = list(
      description = "Hepatic function (normal / not normal)",
      units = "(categorical)",
      type = "categorical",
      notes = "Screened and not retained. 27 normal, 30 not normal and 55 missing (Table 1); the paper does not define the classification."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 108L,
    n_studies = 1L,
    n_observations = 229L,
    age_range = "21 to 88 years",
    age_median = "64 years",
    weight_range = "41 to 160 kg total body weight; 30 to 89 kg lean body weight",
    weight_median = "92 kg total body weight; 57 kg lean body weight",
    sex_female_pct = 45,
    race_ethnicity = c(
      White = 90,
      `Black or African American` = 6,
      Asian = 1,
      Unknown = 3,
      `Hispanic (ethnicity)` = 4
    ),
    disease_state = paste(
      "Critically ill adults admitted to the ICU and treated with",
      "piperacillin-tazobactam as standard of care. Indications (patients",
      "could have more than one): pneumonia 32, sepsis 37, septic shock 30,",
      "urinary tract infection 8, respiratory tract infection 7, abdominal",
      "infection 7, bacteremia 5. Acute kidney injury 26%, mechanical",
      "ventilation 39%, CRRT 8.9% (10 patients). Median SOFA 5."
    ),
    renal_function = "Cockcroft-Gault CLcr on lean body weight median 51 mL/min (range 8 to 173); on total body weight median 82 mL/min (range 10 to 289).",
    dose_range = paste(
      "Piperacillin-tazobactam 3375 mg (3000 mg piperacillin / 375 mg",
      "tazobactam) as 4-hour intravenous infusions in all but one patient,",
      "who received 4500 mg (4000 / 500 mg). Dosing intervals were mostly",
      "every 8 hours (TID or Q8H), with occasional Q12H, QD and QID periods",
      "(Supplementary Table S1). 764 doses in total."
    ),
    regions = "United States (University of Iowa; Vanderbilt University Medical Center).",
    notes = paste(
      "Prospective opportunistic-sampling study; 1 to 8 samples per patient",
      "(mean 2.2). Four of the 112 enrolled patients were excluded from the",
      "tazobactam analysis as outliers (base-model WRES > 5). LC-MS/MS",
      "assay, calibration range 0.25 to 150 mg/L for tazobactam. NONMEM",
      "7.4.3, FOCE with interaction, ADVAN2. Demographics from Table 1;",
      "final estimates from Table 2. The target-attainment simulations in",
      "the paper set the tazobactam free fraction to 0.70; that",
      "factor is a post-processing step and is not part of this model."
    )
  )

  ini({
    # Structural parameters -- Reeder 2025 Table 2 (tazobactam), 'Model
    # Estimate (RSE%)' column. At the reference covariate (CLcr,LBW 51
    # mL/min, no CRRT) total clearance is 2.38 + 1.63 = 4.01 L/h, the Table 2
    # 'CLT, nonCRRT' row; on CRRT it is 1.63 + 0.957 = 2.59 L/h, the 'CLT,
    # CRRT' row.
    lcl_renal <- log(2.38); label("Renal clearance at CLcr,LBW = 51 mL/min (L/h)") # Table 2 row CLR = 2.38 L/h (RSE 23%)
    lcl_nonren <- log(1.63); label("Non-renal clearance (L/h)") # Table 2 row CLNR = 1.63 L/h (RSE 38%)
    lcl_crrt <- log(0.957); label("Clearance by CRRT, replacing renal clearance (L/h)") # Table 2 row CLCRRT = 0.957 L/h (RSE 70%)
    lvc <- log(25.3); label("Volume of distribution (L)") # Table 2 row V = 25.3 L (RSE 14%)

    # Covariate exponent. Table 2 prints it with no RSE and no bootstrap
    # interval (dashes in both columns), as a round theory value, so it is
    # encoded as held constant rather than estimated.
    e_crcl_cl_renal <- fixed(1.00); label("Power exponent of CLcr,LBW / 51 on renal clearance (unitless)") # Table 2 row 'CLR * (CLCR,LBW/51 mL/min)^theta' = 1.00

    # Interindividual variability, exponential model (Methods 'IIV of CL and
    # V were estimated using an exponential model'), so the tabulated CV%
    # converts as omega^2 = log(CV^2 + 1). Table 2 footnote c: the IIV on CL
    # applies to total clearance. No CL-V covariance is reported for the final
    # model, so the two etas are uncorrelated.
    etalcl ~ 0.148420 # Table 2 row 'IIV CL' = 40% CV -> log(0.40^2+1) = 0.148420
    etalvc ~ 0.097490 # Table 2 row 'IIV V' = 32% CV -> log(0.32^2+1) = 0.097490

    # Residual error: combined proportional plus additive (Results 'both
    # models were best explained by a combined proportional and additive
    # error model').
    propSd <- 0.40; label("Proportional residual error (fraction)") # Table 2 row CVCP = 40%
    addSd <- 0.255; label("Additive residual error (mg/L)") # Table 2 row SDCP = 0.255 mg/L
  })

  model({
    # 1. Clearance arms (Table 2 and footnote b).
    cl_renal <- exp(lcl_renal) * (CRCL / 51)^e_crcl_cl_renal
    cl_nonren <- exp(lcl_nonren)
    cl_crrt <- exp(lcl_crrt)

    # 2. Individual PK parameters. Footnote b: CL = CLR x (1 - CRRT) + CLNR +
    # CLCRRT x CRRT, so CRRT replaces only the renal arm. Footnote c: the
    # exponential IIV acts on this total clearance.
    cl <- (cl_renal * (1 - RRT_CRRT_STATUS) + cl_nonren + cl_crrt * RRT_CRRT_STATUS) * exp(etalcl)
    vc <- exp(lvc + etalvc)

    # 3. Micro-constant.
    kel <- cl / vc

    # 4. One-compartment disposition with zero-order (infusion) input and
    # first-order elimination (Results 'one-compartment model with
    # zero-order input and first-order elimination'). Doses enter central as
    # infusions; dose amounts are of tazobactam, not of the
    # piperacillin-tazobactam combination product.
    d/dt(central) <- -kel * central

    # 5. Observation: total plasma tazobactam (mg/L).
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
