Komatsu_2015_digoxin <- function() {
  description <- paste0(
    "Steady-state trough population PK model for oral digoxin in 192 adult ",
    "Japanese cardiology patients (Komatsu 2015). The source analysis used ",
    "routine therapeutic-drug-monitoring trough samples drawn only after ",
    "steady state had been reached, so it fits no absorption or distribution ",
    "process at all: the structural model is the algebraic steady-state ",
    "identity Css = D / (CL * tau), i.e. the trough concentration is the ",
    "daily dose divided by the apparent oral clearance. Apparent oral ",
    "clearance is an ADDITIVE linear function of raw Cockcroft-Gault ",
    "creatinine clearance, multiplied by a fractional factor for the ABSENCE ",
    "of concomitant amiodarone (Komatsu 2015 Results final model and ",
    "Table 3). There are therefore no ODE states, no volume of distribution ",
    "and no absorption rate constant -- the paper's own Discussion states ",
    "'our population model didn't consider volume of distribution or ",
    "absorption phase'. Exponential between-subject variability on CL/F and ",
    "a proportional residual error. Companion digoxin models with full ",
    "disposition structure: Zhou_2010_digoxin, Jelliffe_2014_digoxin, ",
    "Gu_2025_digoxin_pbpk."
  )
  reference <- paste(
    "Komatsu T, Morita M, Miyaji F, Inomata T, Ako J, Atsuda K.",
    "Population pharmacokinetics and optimization of the dosing regimen of",
    "digoxin in adult patients.",
    "J Pharm Health Care Sci. 2015;1:25.",
    "doi:10.1186/s40780-015-0023-6."
  )
  vignette <- "Komatsu_2015_digoxin"
  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "ng/mL"
  )

  covariateData <- list(
    CRCL = list(
      description = paste0(
        "Creatinine clearance estimated from serum creatinine by the ",
        "Cockcroft-Gault method (Komatsu 2015 Methods, Data source), in RAW ",
        "mL/min and NOT BSA-normalized. Time-fixed per subject in this ",
        "analysis: the source is a cross-sectional therapeutic-drug-",
        "monitoring dataset and the paper reports one CLcr per patient."
      ),
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "NOT centred and NOT normalized -- CRCL enters the ADDITIVE linear ",
        "clearance term `exp(lcl) + e_crcl_cl * CRCL` directly, so ",
        "`exp(lcl)` = 1.21 L/h is the clearance intercept extrapolated to ",
        "CRCL = 0 (anuria) rather than a cohort typical value. Cohort mean ",
        "56.17 +/- 33.76 mL/min (Komatsu 2015 Table 1); the paper's dosing ",
        "simulations span 5-130 mL/min, which is the range over which the ",
        "model should be exercised. Patients on dialysis or with rapidly ",
        "deteriorating renal function were excluded, so the model is not ",
        "informed at the extreme low end."
      ),
      source_name = "CLcr"
    ),
    CONMED_AMIO = list(
      description = paste0(
        "1 = subject is coadministered amiodarone, 0 = no concomitant ",
        "amiodarone. 15 of 192 patients (7.8%) were on amiodarone ",
        "(Komatsu 2015 Table 1, Combination medication)."
      ),
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant amiodarone)",
      notes = paste0(
        "ORIENTATION IS INVERTED RELATIVE TO THE SOURCE COLUMN. Komatsu ",
        "2015 defines its own indicator as 'AMD is 0 in the case of ",
        "concomitant administration of amiodarone and 1 otherwise' ",
        "(Abstract, Results), which is the opposite polarity to the ",
        "canonical CONMED_AMIO. The model therefore reconstructs the ",
        "paper's indicator as `amd <- 1 - CONMED_AMIO` before applying the ",
        "published coefficient, so the encoded numbers are the paper's own ",
        "1.21 / 0.0532 / 0.787 with no re-parameterisation. Net effect: ",
        "amiodarone LOWERS apparent digoxin clearance by a factor of 1.787 ",
        "(equivalently, raises the steady-state trough 1.787-fold), which ",
        "agrees with the paper's own statement that 'amiodarone increased ",
        "the trough level of digoxin concentration by approximately ",
        "two-fold' and with its mechanism (inhibition of renal tubular ",
        "secretion and of P-glycoprotein). See the vignette Source trace ",
        "for the five-way dosing-nomogram arithmetic that confirms this ",
        "polarity."
      ),
      source_name = "AMD"
    ),
    DOSE_DIGOXIN_MGD = list(
      description = paste0(
        "The patient's own total daily maintenance dose of oral digoxin ",
        "(mg/day), i.e. the administered amount D divided by the dosing ",
        "interval tau. Entered as a data column rather than as an rxode2 ",
        "dosing event because this model has no compartment to dose into: ",
        "the paper's structural model is the algebraic steady-state ",
        "identity Css = D / (CL * tau) and the dose therefore appears as a ",
        "model INPUT, not as an event record."
      ),
      units = "mg/d",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Komatsu 2015 Table 1 records five regimens, all of which reduce to ",
        "a daily dose: 0.125 mg every 3 days (0.04167 mg/d, n = 7), ",
        "0.125 mg every 2 days (0.0625 mg/d, n = 17), 0.0625 mg/day ",
        "(n = 14), 0.125 mg/day (n = 200) and 0.25 mg/day (n = 49); counts ",
        "are observations, not patients. Must be strictly positive. The ",
        "paper's dosing nomogram (Fig. 4 / Abstract) additionally uses ",
        "0.1875 mg/day. Collapsing D and tau into a single daily-dose ",
        "column is algebraically identical to carrying them separately, ",
        "because they only ever appear in the model as the ratio D/tau."
      ),
      source_name = "Dosage"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 192L,
    n_studies = 1L,
    n_observations = paste0(
      "287 steady-state trough serum digoxin concentrations (Komatsu 2015 ",
      "Table 1). Samples were drawn before the morning dose at least one ",
      "week after digoxin was started."
    ),
    age_range = "Adults; mean 71 +/- 12 years (Komatsu 2015 Table 1)",
    age_median = NA_character_,
    weight_range = "Mean 55.47 +/- 11.94 kg (Komatsu 2015 Table 1)",
    weight_median = NA_character_,
    sex_female_pct = 37.0,
    race_ethnicity = c(Asian = 100),
    disease_state = paste0(
      "Japanese cardiology inpatients and outpatients receiving digoxin for ",
      "congestive heart failure or atrial fibrillation. Left-ventricular ",
      "ejection fraction >= 40% in 156 patients and < 40% in 36 (Komatsu ",
      "2015 Table 1). Patients with major hepatic or gastrointestinal ",
      "disorders, on dialysis, or with rapidly deteriorating renal function ",
      "were excluded."
    ),
    dose_range = paste0(
      "0.0625-0.25 mg/day oral maintenance dosing, including 0.125 mg every ",
      "2 days and every 3 days (Komatsu 2015 Table 1)"
    ),
    regions = "Japan (single centre; Kitasato University Hospital)",
    renal_function = paste0(
      "Cockcroft-Gault CLcr 56.17 +/- 33.76 mL/min (Komatsu 2015 Table 1); ",
      "the paper's dosing simulations span 5-130 mL/min"
    ),
    co_medication = paste0(
      "Amiodarone 15, amlodipine 21, atorvastatin 14, azelnidipine 13, ",
      "bisoprolol 28, carvedilol 53, nifedipine 13, spironolactone 35, ",
      "tolvaptan 8, class I antiarrhythmics 12, class IV antiarrhythmics 31 ",
      "(Komatsu 2015 Table 1). Amlodipine, atorvastatin, bisoprolol, ",
      "carvedilol and tolvaptan were each significant on CL in univariate ",
      "screening (Table 2) but only amiodarone survived backward elimination."
    ),
    notes = paste0(
      "Observed serum digoxin concentration 0.90 +/- 0.56 ng/mL (Komatsu ",
      "2015 Table 1); assay was a cloned enzyme immunoassay with a 0.2 ng/mL ",
      "lower limit of detection and intra- and inter-assay CV below 10%. ",
      "Fitted in NONMEM VI with FOCE; the final model was validated by 200 ",
      "bootstrap replicates (180 converged, Table 4)."
    )
  )

  ini({
    # Structural parameters. Komatsu 2015 final model (Results, and Table 3):
    #   CL/F (L/h) = (1.21 + 0.0532 * CLcr [mL/min]) * (1 + 0.787 * AMD)
    # with AMD = 0 when amiodarone is coadministered and 1 otherwise. The
    # clearance is ADDITIVE and linear in raw creatinine clearance, so
    # exp(lcl) is the intercept at CLcr = 0 rather than a cohort typical
    # value. Table 4 gives the same three thetas with standard errors from
    # the original dataset, alongside their bootstrap means.
    lcl <- log(1.21)
    label("Apparent oral clearance intercept of the additive linear CL/F ~ CRCL model, at CRCL = 0 and with amiodarone present (L/h)")  # Komatsu 2015 Table 4 theta 1 = 1.21 +/- 0.21; same value in the Results final-model equation and Table 3

    e_crcl_cl <- 0.0532
    label("CRCL slope on CL/F (L/h per mL/min)")  # Komatsu 2015 Table 4 theta 2 = 0.0532 +/- 0.0068; same value in the Results final-model equation and Table 3

    e_amio_cl <- 0.787
    label("Fractional increase in CL/F when amiodarone is ABSENT, relative to the amiodarone-present reference (unitless)")  # Komatsu 2015 Table 4 theta 3 = 0.787 +/- 0.187; applied to the paper's AMD indicator, which is 1 when amiodarone is absent

    # Between-subject variability. Exponential IIV on CL/F (Komatsu 2015
    # Methods, Pharmacokinetics model: Pi = TV(Pi) * exp(eta_i)). Table 4
    # reports the VARIANCE from the original dataset as 0.104 +/- 0.017;
    # Table 3 and the Results text report the same quantity on the SD scale
    # as 'Interindividual variance omega(CL) = 32.2 %', and
    # sqrt(0.104) = 0.3225, so the two printings agree and the value below is
    # the variance taken directly from Table 4 with no back-transformation.
    # The Table 4 bootstrap column (0.324) is on the SD scale, not the
    # variance scale -- do not read it as a second variance estimate.
    etalcl ~ 0.104  # Komatsu 2015 Table 4 'omega CL' final-model column = 0.104 +/- 0.017 (variance); cross-checks against the 32.2 % of Table 3

    # Residual error. Proportional (Komatsu 2015 Methods, Pharmacokinetics
    # model: Cobs = Cpred * (1 + eps)). Table 4 reports the VARIANCE as
    # 0.065 +/- 0.010 and Table 3 the SD as 'Intraindividual variance
    # sigma = 25.5 %'; sqrt(0.065) = 0.2550, so the two printings agree.
    # nlmixr2's prop() takes the SD, hence sqrt(0.065).
    propSd <- sqrt(0.065)
    label("Proportional residual error (fraction)")  # Komatsu 2015 Table 4 'sigma' final-model column = 0.065 +/- 0.010 (variance); cross-checks against the 25.5 % of Table 3
  })

  model({
    # Reconstruct the paper's own amiodarone indicator from the canonical
    # column. Komatsu 2015 defines AMD = 0 when amiodarone is coadministered
    # and 1 otherwise, which is the inverse of CONMED_AMIO's canonical
    # 1 = coadministered orientation.
    amd <- 1 - CONMED_AMIO

    # Apparent oral clearance, Komatsu 2015 Results final model:
    #   CL/F (L/h) = (1.21 + 0.0532 * CLcr) * (1 + 0.787 * AMD)
    # Additive linear in raw Cockcroft-Gault CLcr (mL/min), multiplied by the
    # amiodarone-absence factor, with exponential between-subject variability
    # per Methods (Pi = TV(Pi) * exp(eta_i)).
    cl <- (exp(lcl) + e_crcl_cl * CRCL) * (1 + e_amio_cl * amd) * exp(etalcl)

    # Steady-state trough concentration, Komatsu 2015 Methods
    # (Pharmacokinetics model): Css_ij = D_ij / (CL_ij * tau_ij). D/tau is
    # carried as the single daily-dose column DOSE_DIGOXIN_MGD, so
    #   Css = DOSE_DIGOXIN_MGD / CL.
    # Unit reconciliation: DOSE_DIGOXIN_MGD is mg/day and cl is L/h, so
    # cl * 24 converts clearance to L/day; (mg/day)/(L/day) = mg/L = ug/mL,
    # and the factor 1000 converts ug/mL to the paper's reported ng/mL.
    # There is no ODE and no dosing event -- the model is algebraic, exactly
    # as published, and every predicted value is a steady-state trough.
    Cc <- 1000 * DOSE_DIGOXIN_MGD / (24 * cl)
    Cc ~ prop(propSd)
  })
}
