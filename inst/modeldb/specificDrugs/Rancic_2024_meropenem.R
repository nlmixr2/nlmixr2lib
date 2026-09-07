Rancic_2024_meropenem <- function() {
  description <- "One-compartment intravenous population PK model for meropenem in critically ill adults in intensive care (Rancic 2024). Clearance combines a multiplicative power term in serum creatinine and white blood cell count with additive shifts for hypertension and for concomitant vancomycin or colistimethate; the central volume carries no covariates. NOTE: the published central volume (2.05 L) is roughly ten-fold smaller than the value implied by the paper's own reported peak concentrations, so simulated peaks are correspondingly high - see the validation vignette."
  reference <- paste(
    "Rancic A, Milosavljevic MN, Rosic N, Milovanovic D, Folic M,",
    "Ruzic Zecevic D, Petrovic N, Milojevic Corbic M, Dabanovic V,",
    "Jankovic SM.",
    "Population pharmacokinetics of meropenem in critically ill patients.",
    "Open Med (Wars). 2024;19(1):20241004.",
    "doi:10.1515/med-2024-1004. PMCID PMC11278387.",
    sep = " "
  )
  vignette <- "Rancic_2024_meropenem"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "meropenem", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CREAT = list(
      description        = "Serum creatinine concentration.",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters clearance as the RAW (uncentred, unnormalised) power term CREAT^0.000001.",
        "The exponent is estimated at essentially zero (Table 3 reports the estimate, its",
        "standard error and both confidence limits as 0.0000; the Results equation prints",
        "0.000001), so the term evaluates to 1.0000 for every physiologically plausible",
        "creatinine and the covariate is numerically inert despite having survived backward",
        "deletion at p < 0.01. It is retained here because the paper retained it.",
        "Cohort mean 94.47 +/- 63.68 umol/L, range 40 - 452 (Table 1).",
        "Note that this model's creatinine effect on clearance is POSITIVE in direction",
        "(Discussion: 'higher serum concentrations of creatinine were associated with larger",
        "meropenem clearance'), which the authors attribute to augmented renal clearance in",
        "critically ill patients rather than to renal impairment.",
        sep = " "
      ),
      source_name        = "CRE"
    ),
    WBC = list(
      description        = "Total white blood cell (leukocyte) count.",
      units              = "10^9 cells/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters clearance as the RAW (uncentred, unnormalised) power term WBC^(-0.165).",
        "The paper does not state the units of WBC and does not report its distribution in",
        "Table 1, so the unit was settled by internal consistency: the base model's typical",
        "clearance is 3.80 L/h, and the final model's covariate-adjusted clearance can only",
        "reproduce that magnitude if WBC is on the 10^9/L scale. Solving",
        "5.29 * WBC^(-0.165) + 0.23 = 3.80 (0.23 being the cohort-average contribution of the",
        "additive vancomycin and colistimethate terms at their Table 1 prevalences) gives",
        "WBC = 10.8 x 10^9/L, a typical leukocytosis for this population. On a cells/uL scale",
        "the same arithmetic would require WBC around 10,800 and give a typical clearance of",
        "1.13 L/h, a third of the base-model value, so cells/uL is excluded.",
        "Direction: the ESTIMATE IS NEGATIVE (-0.165, 95% CI -0.219 to -0.110; Table 3), so",
        "clearance DECREASES as the leukocyte count rises. The Discussion and Conclusion",
        "state the opposite ('associated with greater clearance of meropenem'). The signed",
        "estimate and its confidence interval are the primary record and are what is encoded;",
        "see the vignette Errata.",
        sep = " "
      ),
      source_name        = "WBCs"
    ),
    DIS_HYPERT = list(
      description        = "Hypertension recorded as a comorbidity at study entry.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no hypertension)",
      notes              = paste(
        "Enters clearance as an ADDITIVE shift (+0.000001 L/h when present), i.e. it is",
        "numerically inert: Table 3 reports the estimate, its standard error and both",
        "confidence limits as 0.0000, and the Results equation prints 0.000001. Retained",
        "here because the paper retained it in the final model (backward deletion at",
        "p < 0.01). 33 of 101 patients (32.7%) were hypertensive (Table 1).",
        sep = " "
      ),
      source_name        = "HTA"
    ),
    CONMED_VANCOMYCIN = list(
      description        = "Concomitant intravenous vancomycin during the meropenem course.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant vancomycin)",
      notes              = paste(
        "Enters clearance as an ADDITIVE shift of +0.825 L/h (95% CI 0.770 - 0.879;",
        "Table 3), not as a multiplicative factor. 17 of 101 patients (16.8%) received",
        "vancomycin and 3 (2.97%) received both vancomycin and colistin, so the two",
        "comedication indicators are not mutually exclusive and their additive shifts sum.",
        "The authors read the effect as a marker of augmented renal clearance in the",
        "sicker patients who need combination therapy rather than as a mechanistic",
        "drug-drug interaction (Discussion).",
        sep = " "
      ),
      source_name        = "VAN"
    ),
    CONMED_COLISTIMETHATE = list(
      description        = "Concomitant intravenous colistimethate sodium during the meropenem course.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant colistimethate)",
      notes              = paste(
        "Enters clearance as an ADDITIVE shift of +1.28 L/h (95% CI 1.23 - 1.33;",
        "Table 3). 7 of 101 patients (6.9%) received colistin, 3 of whom (2.97%) also",
        "received vancomycin. Same interpretive caveat as CONMED_VANCOMYCIN: the paper",
        "reads it as a severity / augmented-renal-clearance marker (Discussion), and the",
        "estimate rests on 7 patients, which the authors name as the study's main",
        "limitation.",
        sep = " "
      ),
      source_name        = "COL"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age at study entry.",
      units       = "years",
      type        = "continuous",
      notes       = paste(
        "Screened as (AGE/50)^theta3 and significant univariately (MOF 1277.958 vs 1283.053,",
        "difference 5.095, p < 0.05; Table 2) and retained in the full model, but dropped at",
        "backward deletion, which required a MOF increase above 6.64. Cohort mean",
        "62.37 +/- 14.89 years, range 21 - 86 (Table 1)."
      )
    ),
    WT = list(
      description = "Total body weight.",
      units       = "kg",
      type        = "continuous",
      notes       = paste(
        "Screened as (TBW/70)^theta4 and significant univariately (MOF 1266.640, difference",
        "16.413, p < 0.01; Table 2, where the row label is misprinted 'TWB') and retained in",
        "the full model, but dropped at backward deletion. The final model therefore carries",
        "NO body-size scaling at all, on either clearance or volume. Cohort mean",
        "78.97 +/- 13.76 kg, range 48 - 130 (Table 1)."
      )
    ),
    SEXF = list(
      description = "Female sex indicator.",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Listed among the demographic data collected for every patient (Methods 2.1) and",
        "among the 24 screened covariates, but it is not one of the 13 covariates carried",
        "into the full model (Table 2) and no estimate is reported. 39 of 101 patients",
        "(38.6%) were female (Table 1)."
      )
    ),
    DIS_CANCER = list(
      description = "Malignant disease (neoplasm) recorded as a comorbidity at study entry.",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened as an additive shift (theta10 x NEO) and NOT significant univariately",
        "(MOF 1280.925, difference 2.128, p > 0.05; Table 2). 18 of 101 patients (17.8%)",
        "had a neoplasm (Table 1)."
      )
    ),
    AST = list(
      description = "Aspartate aminotransferase activity.",
      units       = "U/L",
      type        = "continuous",
      notes       = "Collected as one of the basic biochemical parameters screened as covariates (Methods 2.1); not among the 13 covariates carried into the full model and no estimate or cohort summary is reported."
    ),
    ALT = list(
      description = "Alanine aminotransferase activity.",
      units       = "U/L",
      type        = "continuous",
      notes       = "Collected as one of the basic biochemical parameters screened as covariates (Methods 2.1); not among the 13 covariates carried into the full model and no estimate or cohort summary is reported."
    ),
    RBC = list(
      description = "Red blood cell count.",
      units       = "10^12 cells/L",
      type        = "continuous",
      notes       = "Collected as one of the basic biochemical parameters screened as covariates (Methods 2.1); not among the 13 covariates carried into the full model and no estimate or cohort summary is reported."
    ),
    HGB = list(
      description = "Hemoglobin concentration.",
      units       = "g/L",
      type        = "continuous",
      notes       = "Collected as one of the basic biochemical parameters screened as covariates (Methods 2.1); not among the 13 covariates carried into the full model and no estimate or cohort summary is reported."
    ),
    PLT = list(
      description = "Platelet count.",
      units       = "10^9 cells/L",
      type        = "continuous",
      notes       = "Collected as one of the basic biochemical parameters screened as covariates (Methods 2.1); not among the 13 covariates carried into the full model and no estimate or cohort summary is reported."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 101L,
    n_studies      = 1L,
    n_observations = 202L,
    age_range      = "21 - 86 years (mean 62.37 +/- 14.89)",
    weight_range   = "48 - 130 kg (mean 78.97 +/- 13.76)",
    sex_female_pct = 38.6,
    disease_state  = paste(
      "Critically ill adults in the intensive care unit with severe infection (meningitis,",
      "pneumonia, sepsis, septic shock or febrile neutropenia) caused by multi-resistant",
      "Gram-negative bacteria. Comorbidities recorded: hypertension 32.7%, chronic renal",
      "failure 11.9%, neoplasm 17.8%, cerebral infarction 20.8%, pneumonia 31.7%, urinary",
      "tract infection 35.6% (Table 1). Serum creatinine mean 94.47 +/- 63.68 umol/L,",
      "range 40 - 452, so most patients were not renally impaired; the authors emphasise",
      "augmented renal clearance rather than renal failure as the dominant phenomenon.",
      sep = " "
    ),
    dose_range     = paste(
      "1,000 - 2,000 mg meropenem every 8 or 12 h by intermittent intravenous infusion",
      "(total daily dose mean 3,000 +/- 692.82 mg/day, range 2,000 - 6,000). Patients were",
      "sampled only after at least 3 days of continuous therapy, i.e. at steady state.",
      sep = " "
    ),
    regions        = "Serbia (Intensive Care Unit, University Clinical Centre Kragujevac).",
    notes          = paste(
      "Prospective observational case-series. Two plasma samples per patient: the first",
      "5 - 30 min after the end of the infusion (mean 40.69 +/- 16.67 mg/L, range",
      "13.07 - 88.95) and the second 3 - 4 h after the end of the infusion (mean",
      "12.55 +/- 7.62 mg/L, range 2.06 - 36.28); assayed by HPLC-DAD at 300 nm.",
      "NONMEM (version 5, level 1.1) ADVAN1, i.e. one compartment without absorption.",
      "24 covariates were screened univariately; 13 entered the full model (MOF 1207.153)",
      "and 5 survived backward deletion at p < 0.01 into the final model (MOF 1224.035",
      "against a base model MOF of 1283.053). Validated by 100 bootstrap replicates",
      "(Table 5), whose point estimates differ substantially from the final model's -",
      "see the vignette Errata. Comorbidity flags for chronic renal failure, cerebral",
      "infarction, pneumonia and urinary tract infection were also screened (Table 2) but",
      "are not recorded in covariatesDataExcluded because they carry no effect in the",
      "final model and their mapping onto the existing DIS_ register entries is not",
      "settled; see the vignette for the full screen.",
      sep = " "
    )
  )

  ini({
    # =========================================================================
    # Structural parameters, Table 3 'Parameter estimates for the final model'.
    #
    # theta1 is the INTERCEPT of the multiplicative power term, not the typical
    # clearance itself: the final-model equation (Results, and reproduced in the
    # Abstract) is
    #
    #   CL (L/h) = 5.29 * CRE^0.000001 * WBCs^(-0.165)
    #              + 0.000001 * HTA + 0.825 * VAN + 1.28 * COL
    #
    # so the covariate-adjusted typical clearance for a patient with a leukocyte
    # count near the cohort norm and no comedication is about 3.6 - 3.8 L/h, in
    # line with the base model's 3.80 L/h (Results paragraph 2).
    # =========================================================================
    lcl <- log(5.29); label("Clearance intercept theta1 (CL, L/h)")             # Table 3 row 'Clearance (L/h) (theta 1)'
    lvc <- log(2.05); label("Central volume of distribution (V, L)")            # Table 3 row 'Volume of distribution (L) (theta 2)'

    # =========================================================================
    # Covariate effects on clearance, Table 3. Two of the five are estimated at
    # (essentially) zero and are numerically inert; they are encoded at the
    # values the paper's Results equation prints because the paper retained them
    # in the final model. The two power terms take RAW covariate values - the
    # final-model equation has no centring or normalising constant, unlike the
    # univariate screens in Table 2, which used (AGE/50) and (TBW/70).
    # =========================================================================
    e_creat_cl                 <- 0.000001; label("Power exponent on raw serum creatinine for CL (unitless)")   # Results final-model equation; Table 3 row 'Effect of CRE' (0.0000, SE 0.0000)
    e_wbc_cl                   <- -0.165;   label("Power exponent on raw leukocyte count for CL (unitless)")    # Table 3 row 'Effect of WBCs' (-0.165, 95% CI -0.219 to -0.110)
    e_dis_hypert_cl            <- 0.000001; label("Additive shift on CL for hypertension (L/h)")                # Results final-model equation; Table 3 row 'Effect of HTA' (0.0000, SE 0.0000)
    e_conmed_vancomycin_cl     <- 0.825;    label("Additive shift on CL for concomitant vancomycin (L/h)")      # Table 3 row 'Effect of VAN' (0.825, 95% CI 0.770 - 0.879)
    e_conmed_colistimethate_cl <- 1.28;     label("Additive shift on CL for concomitant colistimethate (L/h)")  # Table 3 row 'Effect of COL' (1.28, 95% CI 1.23 - 1.33)

    # =========================================================================
    # Between-subject variability. Exponential on clearance, multiplying the
    # whole covariate-adjusted expression (Table 2 base and full model rows, both
    # of the form '... x Exp[ETA(1)]'). Table 3 reports omega^2 directly as
    # 'Inter-individual variance of CL (omega^2 CL) = 0.0215', so no CV-to-
    # variance conversion is needed; that variance corresponds to a 14.7% CV.
    #
    # The Results prose restates the same number as '2.15%' ('Inter- and
    # intra-individual variability were 2.15 and 44%, respectively'), i.e. it
    # multiplies the VARIANCE by 100 and calls the product a percentage. The
    # table's explicit omega^2 / sigma^2 notation, standard errors and confidence
    # limits are the primary record and are what is used here.
    #
    # No variability is reported on the volume of distribution: the Table 3
    # footnote defines an ETA(2) on V, but no estimate for it appears in any
    # table, so V is a typical value only.
    # =========================================================================
    etalcl ~ 0.0215  # Table 3 row 'Inter-individual variance of CL (omega^2 CL)' = 0.0215 (SE 0.0146, 95% CI 0.0000 - 0.0501)

    # =========================================================================
    # Residual error. Table 3 reports 'Residual error variance (sigma^2 CL) =
    # 0.44' (SE 0.066, 95% CI 0.310 - 0.570), so the residual SD is
    # sqrt(0.44) = 0.6633.
    #
    # Form: Methods 2.3 says 'Exponential and additive error models were used to
    # evaluate the estimation of intra-individual variability in clearance and
    # residual error in concentration' - i.e. exponential for the clearance ETA
    # (confirmed by the Exp[ETA(1)] rows of Table 2) and additive for the
    # residual. In NONMEM an 'additive' residual on log-transformed observations
    # is exactly a proportional / log-normal residual once back-transformed to
    # the linear concentration scale, and that is the reading encoded here,
    # because an additive residual of 0.66 mg/L on the LINEAR scale is
    # arithmetically impossible for this dataset: it is 1.6% of the mean
    # observed peak of 40.69 mg/L, whereas the paper's own final model reports
    # RMSPE 15.86 (Table 4) and observed between-patient CVs of 41% on the peak
    # and 61% on the 3 - 4 h sample (Table 1), none of which a 0.66 mg/L additive
    # term can generate. A proportional SD of 0.6633 can.
    #
    # The Results prose restates sigma^2 = 0.44 as '44%', by the same
    # variance-times-100 convention it applies to omega^2. See vignette Errata.
    # =========================================================================
    propSd <- 0.6633; label("Proportional residual SD on Cc (fraction)")  # sqrt(0.44); Table 3 row 'Residual error variance (sigma^2 CL)'
  })

  model({
    # -----------------------------------------------------------------------
    # 1. Individual clearance. The exponential ETA multiplies the ENTIRE
    #    covariate-adjusted bracket, additive terms included, exactly as
    #    written in the Table 2 full-model row:
    #
    #      CL = [theta1 * (...)^... * (...)^... + theta8*HTA + ... ] * Exp[ETA(1)]
    #
    #    Covariates enter RAW: creatinine in umol/L and the leukocyte count in
    #    10^9/L, with no centring constant anywhere in the final-model equation.
    # -----------------------------------------------------------------------
    cl <- (exp(lcl) * CREAT^e_creat_cl * WBC^e_wbc_cl +
             e_dis_hypert_cl * DIS_HYPERT +
             e_conmed_vancomycin_cl * CONMED_VANCOMYCIN +
             e_conmed_colistimethate_cl * CONMED_COLISTIMETHATE) * exp(etalcl)

    # -----------------------------------------------------------------------
    # 2. Central volume - a typical value with no covariates and no reported
    #    between-subject variability (Table 3).
    # -----------------------------------------------------------------------
    vc <- exp(lvc)

    # -----------------------------------------------------------------------
    # 3. One-compartment intravenous disposition (NONMEM ADVAN1; Methods 2.3
    #    'It implies the application of a single-compartment model without
    #    absorption'). The state carries drug amount in mg and the volume is in
    #    L, so Cc is in mg/L, which is the mg/L equivalent of the ug/mL the
    #    paper reports.
    # -----------------------------------------------------------------------
    kel <- cl / vc

    d/dt(central) <- -kel * central

    # -----------------------------------------------------------------------
    # 4. Observation and residual error.
    # -----------------------------------------------------------------------
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
