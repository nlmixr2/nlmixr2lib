Komatsu_2016_penicillinG <- function() {
  description <- paste0(
    "One-compartment intravenous population PK model for penicillin G ",
    "(benzylpenicillin, given as the potassium salt) in 25 Japanese adults ",
    "treated for suspected or documented infective endocarditis (Komatsu ",
    "2016, 46 serum samples, NONMEM VI ADVAN1 TRANS2). Clearance is a ",
    "THROUGH-ORIGIN linear function of Cockcroft-Gault creatinine ",
    "clearance with no intercept term, CL (L/h) = 0.21 x CLcr (mL/min), so ",
    "a patient with no residual renal function is predicted to have no ",
    "penicillin G clearance at all; the volume of distribution 28.9 L ",
    "carries between-subject variability but no covariate. Body weight, ",
    "serum creatinine, ALT, sex and age were all individually significant ",
    "on CL in forward inclusion (Table 3) but only CLcr survived backward ",
    "elimination, and body weight was not significant on Vd. Because the ",
    "samples were drawn only at the trough and 2 or 3 h after a dose, the ",
    "authors chose one compartment over the three-compartment structure ",
    "published elsewhere for penicillin G. The companion static ",
    "exposure-response model for clinical outcome is ",
    "Komatsu_2016_penicillinG_clinical_success."
  )
  reference <- paste(
    "Komatsu T, Inomata T, Watanabe I, Kobayashi M, Kokubun H, Ako J, Atsuda K.",
    "Population pharmacokinetic analysis and dosing regimen optimization of",
    "penicillin G in patients with infective endocarditis.",
    "J Pharm Health Care Sci. 2016 Apr 5;2:9.",
    "doi:10.1186/s40780-016-0043-x. PMCID PMC4820900.",
    sep = " "
  )
  vignette <- "Komatsu_2016_penicillinG"
  # Doses in mg and a volume in L give Cc in mg/L, which is numerically
  # identical to the ug/mL that Komatsu 2016 prints throughout; mg/L is
  # declared here because it is the dimensionally consistent statement
  # against units$dosing.
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    CRCL = list(
      description        = "Creatinine clearance estimated from serum creatinine by the Cockcroft-Gault equation. RAW mL/min, NOT BSA-normalized.",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Komatsu 2016 Methods, Data source: 'CLcr was estimated from the ",
        "serum creatinine level by the Cockcroft-Gault method'. Cohort mean ",
        "82.52 mL/min (SD 33.29, range 11-144; Table 1). Enters clearance as ",
        "a bare product with NO centering and NO intercept: the Table 5 ",
        "final model is CL (L/h) = theta1 x CLcr (mL/min) with theta1 = ",
        "0.21, so the covariate is not a multiplier on a typical clearance ",
        "but the whole of it, and CL is exactly proportional to renal ",
        "function. This through-origin form is the unusual feature of this ",
        "model and the reason the parameter label carries the compound unit ",
        "L/h per mL/min rather than L/h. The authors extrapolate it in their ",
        "own Monte Carlo dosing simulations down to CLcr = 5 mL/min, below ",
        "the observed minimum of 11 mL/min, and their nomogram (Figure 5) ",
        "recommends a regimen for CLcr below 20 mL/min; at CLcr = 0 the ",
        "model predicts no elimination at all, which is not physiologic ",
        "(penicillin G retains a non-renal clearance component), so do not ",
        "use this model at or near zero renal function."
      ),
      source_name        = "CLcr"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened on both CL and Vd (Komatsu 2016 Table 3). On CL the linear form theta1 + theta2 x BW dropped the objective function by 32.9 (P < 0.001) in forward inclusion but did not survive backward elimination once CLcr was in the model; on Vd the form theta3 + theta4 x BW gave a -2 log likelihood change of exactly 0 and is reported as not significant. No retained point estimate exists."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units       = "mg/dL",
      type        = "continuous",
      notes       = "Screened on CL as theta1 + theta2 x (1/sCr) (Komatsu 2016 Table 3, -2 log likelihood change 75.919, P < 0.001) but eliminated in favour of CLcr, which is itself derived from serum creatinine by Cockcroft-Gault and is the stronger predictor (-2 log likelihood change 117.303). Cohort mean 0.92 mg/dL (SD 0.56, range 0.51-3.26; Table 1). No retained point estimate exists."
    ),
    ALT = list(
      description = "Alanine aminotransferase, screened as a binary split at 40 IU/L",
      units       = "IU/L",
      type        = "binary",
      notes       = "Screened on CL as the multiplicative binary theta1 x theta2^ALT with ALT > 40 coded 0 and ALT < 40 coded 1 (Komatsu 2016 Table 3, -2 log likelihood change 19.27, P < 0.001), but not retained after backward elimination. Cohort mean 26.80 IU/L (SD 26.45, range 4-106; Table 1). No retained point estimate exists."
    ),
    SEXF = list(
      description = "Sex",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened on CL as the multiplicative binary theta1 x theta2^sex with male coded 1 and female coded 0 (Komatsu 2016 Table 3, -2 log likelihood change 31.334, P < 0.001), but not retained. Note the source codes MALE as 1; the canonical SEXF codes FEMALE as 1, so any future use of this screen would need SEXF = 1 - sex. Cohort 16 male : 9 female (Table 1). No retained point estimate exists."
    ),
    AGE = list(
      description = "Age, screened as a binary split at 65 years",
      units       = "years",
      type        = "binary",
      notes       = "Screened on CL as the multiplicative binary theta1 x theta2^age with age > 65 coded 1 and age < 64 coded 0 (Komatsu 2016 Table 3, -2 log likelihood change 8.608, P < 0.005 -- the weakest of the six screened covariates), but not retained. Cohort mean 54 years (SD 17, range 21-83; Table 1). No retained point estimate exists."
    )
  )

  compartmentData <- list(
    central = list(analyte = "penicillin G", units = "mg", specimen = "serum", verified = TRUE)
  )

  population <- list(
    species         = "human",
    n_subjects      = 25L,
    n_studies       = 1L,
    n_observations  = "46 serum penicillin G concentrations (Komatsu 2016 Table 1)",
    age_range       = "21-83 years",
    age_median      = "mean 54 years (SD 17); no median reported",
    weight_range    = "33-86.9 kg",
    weight_median   = "mean 55.35 kg; no median reported (see notes on the printed SD)",
    sex_female_pct  = 36,
    race_ethnicity  = "Not reported; single-centre Japanese cohort (Kitasato University Hospital, Sagamihara)",
    disease_state   = "Suspected or documented infective endocarditis. Viridans group streptococci were isolated in 21 of the 25 patients; 15 of those 21 responded to penicillin G and 6 failed (Table 2). Treatment failure was defined as persistence of fever and/or bacteremia requiring a change of antibiotic, or infection-related mortality within 30 days.",
    renal_function  = "Creatinine clearance (Cockcroft-Gault) mean 82.52 mL/min (SD 33.29, range 11-144); serum creatinine mean 0.92 mg/dL (SD 0.56, range 0.51-3.26)",
    hepatic_function = "Alanine aminotransferase mean 26.80 IU/L (SD 26.45, range 4-106)",
    dose_range      = "Penicillin G potassium (Meiji Seika Pharma) intravenously; the individual clinical regimens are not tabulated. The paper's Monte Carlo dosing simulations span 0.5 million IU every 6 h to 4 million IU every 4 h, plus 1 million IU/h by continuous infusion (24 million IU/day).",
    regions         = "Japan (Kitasato University Hospital, Sagamihara, Kanagawa); patients treated between January 1997 and April 2013",
    notes           = paste0(
      "Baseline demographics are Komatsu 2016 Table 1. Observed serum ",
      "penicillin G concentrations spanned 0.5-212.3 ug/mL (mean 33.2, SD ",
      "45.2); the lower limit of detection of the HPLC assay was 0.5 ug/mL. ",
      "Samples were drawn immediately before a dose and 2 or 3 h after it, ",
      "which is why the authors could not resolve a distribution phase. ",
      "The Table 1 body-weight SD is printed as 33.31 kg against a range of ",
      "33-86.9 kg and a mean of 55.35 kg, which is not attainable; it ",
      "appears to be a transcription of the adjacent CLcr SD of 33.29 and is ",
      "recorded here as printed. Weight is not used by this model."
    )
  )

  ini({
    # ==================================================================
    # Structural parameters: Komatsu 2016 Table 5 'Final population
    # pharmacokinetic parameters of penicillin G', whose point estimates
    # agree with the final-model column of Table 4 and with the
    # final-model equation printed in the Results text ('The final model
    # was: CL (L/h) = 0.21 x CLcr (mL/min), Vd (L) = 28.9').
    #
    # theta1 is NOT a clearance. It is the slope of clearance on
    # creatinine clearance in a through-origin line, so its unit is
    # (L/h) per (mL/min) and the model has no clearance intercept. At
    # the cohort mean CLcr of 82.52 mL/min it gives CL = 17.3 L/h.
    #
    # Falsifier walked before committing these values: the Discussion
    # states a serum half-life of 0.79 h 'calculated by keeping the CLcr
    # fixed at 120 mL/min'. 0.693 * 28.9 / (0.21 * 120) = 0.795 h, which
    # reproduces the printed figure and pins both parameters, the
    # through-origin form, and the mL/min unit of the covariate
    # simultaneously.
    # ==================================================================
    lcl <- log(0.21);  label("Slope of clearance on Cockcroft-Gault creatinine clearance, CL = lcl_exp x CRCL (L/h per mL/min)")  # Komatsu 2016 Table 5 row 'CL (L/h) = theta 1 x CLcr (mL/min)', theta1 = 0.21 (RSE 8.81%, bootstrap 95% CI 0.171-0.249); Table 4 final-model column gives the same 0.21 (95% CI 0.173-0.246)
    lvc <- log(28.9); label("Volume of distribution (L)")                                                                        # Komatsu 2016 Table 5 row 'Vd (L) = theta 2', theta2 = 28.9 (RSE 8.58%, bootstrap 95% CI 23.4-34.7); Table 4 final-model column gives the same 28.9 (95% CI 24.0-33.7)

    # ==================================================================
    # Between-subject variability. Komatsu 2016 Methods, Pharmacokinetic
    # calculations: 'Inter-individual variability of the parameters was
    # best explained by an exponential error model (Pi = TV(Pi) x
    # exp(eta i))' with 'variance of omega 2'.
    #
    # Omega SCALE is settled, not assumed: the Results text states 'The
    # coefficients of variation of the inter-individual variability
    # (omega 2) of CL, Vd, and the residual variability (sigma 2) were
    # 28.8, 32.4, and 17.4 %, respectively'. sqrt(0.0835) = 28.9 %,
    # sqrt(0.104) = 32.2 %, sqrt(0.0304) = 17.4 % -- all three tabulated
    # numbers are VARIANCES whose square roots are the printed CVs, so
    # the Table 5 column may be used directly as an rxode2 omega.
    # ==================================================================
    etalcl ~ 0.0835  # Komatsu 2016 Table 5 row 'eta CL' = 0.0835 (a variance; RSE 33.1%, bootstrap 95% CI 0.0171-0.1397); Table 4 final-model column 0.0835, 95% CI 0.0292-0.137
    etalvc ~ 0.104   # Komatsu 2016 Table 5 row 'eta vd' = 0.104 (a variance; RSE 5.22%, bootstrap 95% CI 0.0087-0.1203); Table 4 final-model column 0.104, 95% CI 0.094-0.115

    # ==================================================================
    # Residual variability. Komatsu 2016 Methods: 'The residual
    # (intra-individual) variability of the parameters was also
    # explained by a proportional error model (Cobs,ij = Cpred,ij x (1 +
    # eps ij))' with 'variance sigma 2'. Table 5 reports sigma^2 =
    # 0.0304; rxode2's prop() takes the SD, so sqrt(0.0304) = 0.174356.
    # ==================================================================
    propSd <- 0.174356; label("Proportional residual error SD (unitless fraction of the predicted concentration)")  # sqrt of Komatsu 2016 Table 5 row 'epsilon' = 0.0304 (a variance; RSE 5.85%, bootstrap 95% CI 0.0265-0.0342); the Results text prints the corresponding CV as 17.4%
  })

  model({
    # Clearance is the covariate term alone -- there is no typical-value
    # intercept to multiply. Writing it as exp(lcl + etalcl) * CRCL keeps
    # the log-normal between-subject distribution the authors used while
    # reproducing CL = 0.21 x CLcr at the typical value.
    cl <- exp(lcl + etalcl) * CRCL
    vc <- exp(lvc + etalvc)

    kel <- cl / vc

    # One compartment with first-order elimination (NONMEM ADVAN1
    # TRANS2). Penicillin G was given intravenously, so there is no
    # absorption compartment and no bioavailability term; supply an
    # infusion through the event table's rate or duration column, or
    # dose as a bolus.
    d/dt(central) <- -kel * central

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
