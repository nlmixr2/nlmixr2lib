Launay_2024_ceftazidime <- function() {
  description <- "One-compartment population PK model for ceftazidime in septic critically ill adults receiving continuous intravenous infusion, with CKD-EPI estimated GFR as a power covariate on clearance. Developed to size the loading dose needed before continuous infusion; the ICU volume of distribution (88 L) is roughly 2.4-fold the median of previously published ICU models (37.2 L), which is why a 2 g loading dose is shown to be insufficient."
  reference <- paste(
    "Launay M, Ollier E, Kably B, Le Louedec F, Thiery G, Lanoiselee J,",
    "Perinel-Ragey S. Loading Dose of Ceftazidime Needs to Be Increased in",
    "Critically Ill Patients: A Retrospective Study to Evaluate Recommended",
    "Loading Dose with Pharmacokinetic Modelling.",
    "Antibiotics. 2024;13(8):756. doi:10.3390/antibiotics13080756",
    sep = " "
  )
  vignette <- "Launay_2024_ceftazidime"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    CRCL = list(
      description        = "Glomerular filtration rate estimated with the CKD-EPI equation, BSA-normalized",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Creatinine-based CKD-EPI eGFR, BSA-normalized to mL/min/1.73 m^2",
        "(Launay 2024 Methods section 4.4 and Table 2 footnote). Enters as a",
        "power term (CRCL / 73.90)^0.9 on clearance; 73.90 mL/min/1.73 m^2 is",
        "the population median stated in the Table 2 footnote (NOT the mean of",
        "87.6 given in Table 1 -- the footnote is explicit that GFRmedian is a",
        "median). Cohort mean (SD) was 87.6 (74.2), so the distribution is",
        "strongly right-skewed and spans augmented renal clearance through",
        "renal failure. Patients on continuous renal replacement therapy were",
        "excluded from the validation cohort.",
        sep = " "
      ),
      source_name        = "GFR"
    )
  )

  # Covariates that Launay 2024 screened in the stepwise covariate search
  # (Methods section 4.4) but did NOT retain in the final model. Recorded for
  # provenance only; none is referenced in model().
  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = paste(
        "Screened but not retained. Mean (SD) 91.3 (25.3) kg (Table 1). The",
        "Discussion explicitly notes that weight and BMI were not identified as",
        "significant covariates on Vd, attributed to ~48% of the cohort being",
        "obese and therefore to a high mean Vd with comparatively little spread.",
        sep = " "
      )
    ),
    BMI = list(
      description = "Body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = paste(
        "Screened but not retained. Mean (SD) 31.2 (9.5) kg/m^2, 47.7% obese",
        "(Table 1). Available only in the model-development dataset.",
        sep = " "
      )
    ),
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened but not retained. Mean (SD) 64.5 (11.9) years (Table 1)."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened but not retained. 67 of 86 (77.9%) were male, i.e. 22.1% female (Table 1)."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units       = "not reported",
      type        = "continuous",
      notes       = paste(
        "Screened but not retained as a covariate in its own right (Methods",
        "section 4.4 lists serum creatinine at baseline and at sampling). It",
        "enters the model only indirectly, through the CKD-EPI eGFR carried in",
        "CRCL. Neither the unit nor summary statistics are reported.",
        sep = " "
      )
    ),
    TPRO = list(
      description = "Serum total protein concentration",
      units       = "not reported",
      type        = "continuous",
      notes       = "Screened but not retained (Methods section 4.4). No summary statistics reported."
    ),
    DIS_COVID19 = list(
      description = "Ongoing COVID-19 infection indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened but not retained (Methods section 4.4 lists 'status of ongoing COVID infection'). Prevalence not reported."
    )
  )

  compartmentData <- list(
    central = list(analyte = "ceftazidime", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 86,
    n_studies      = 6,
    n_observations = 223,
    age_mean       = "64.5 years (SD 11.9)",
    weight_mean    = "91.3 kg (SD 25.3)",
    bmi_mean       = "31.2 kg/m^2 (SD 9.5); 47.7% obese",
    sex_female_pct = 22.1,
    renal_function = "CKD-EPI eGFR mean 87.6 mL/min/1.73 m^2 (SD 74.2), population median 73.90 mL/min/1.73 m^2",
    disease_state  = "septic critically ill adults in intensive care",
    dose_range     = "2 g intravenous loading dose followed by continuous infusion, median 6 g/day",
    regions        = "France (six ICUs in the Saint-Etienne area)",
    notes          = paste(
      "Retrospective therapeutic drug monitoring cohort, 1 November 2019 to",
      "31 October 2021; 86 patients contributing 223 samples (1 to 9 per",
      "patient), all drawn at least 6 h after the start of the continuous",
      "infusion (Launay 2024 Table 1 and Methods sections 4.2-4.3). Assay",
      "calibration range 8-150 mg/L; values below the limit of quantification",
      "were set to 4 mg/L. An independent external validation cohort of 32",
      "patients (32 samples, one Paris ICU, mean eGFR 86.4 mL/min/1.73 m^2,",
      "mean weight 77.4 kg) was used for evaluation only and did not",
      "contribute to parameter estimation; patients on continuous renal",
      "replacement therapy were excluded from it.",
      sep = " "
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Structural parameters. Launay 2024 Table 2, "Population" block. The
    # model is one-compartment with intravenous input; the paper reports
    # the disposition parameters directly as CL and Vd (no bioavailability
    # term, since every dose is intravenous).
    # ------------------------------------------------------------------
    lcl <- log(4.45); label("Clearance (L/h)")                 # Launay 2024 Table 2, Population CL = 4.45 L/h (RSE 7.3%)
    lvc <- log(88.0); label("Volume of distribution (L)")      # Launay 2024 Table 2, Population Vd = 88.0 L (RSE 18.3%)

    # ------------------------------------------------------------------
    # Covariate effect. Launay 2024 Table 2 prints the full equation:
    #     CL(i) = CL(pop) * (GFRi / GFRmedian)^(beta^CL_GFR) * exp(eta_CL)
    # with GFRmedian = 73.90 mL/min/1.73 m^2 given in the table footnote.
    # The Results section calls this "eGFR as a covariate on CL using
    # allometric scaling", i.e. a power model, matching the printed form.
    # ------------------------------------------------------------------
    e_crcl_cl <- 0.9; label("Power exponent on (CRCL / 73.90) for clearance (unitless)")  # Launay 2024 Table 2, Covariate effect beta^CL_GFR = 0.9 (RSE 15.5%)

    # ------------------------------------------------------------------
    # Interindividual variability. Table 2 labels this block
    # "Interindividual variability (standard deviation)", so the printed
    # 0.46 and 0.57 are SDs of the random effects and are squared here to
    # give variances. The exponential form exp(eta) is printed explicitly
    # in the Table 2 equation for CL.
    #
    # The Vd SD is independently confirmed by the paper's own loading-dose
    # arithmetic: Methods section 4.4 gives LD = TC * Vd * exp(kappa *
    # omega_V) with TC = 60 mg/L and kappa = 0.84 (the 80th-percentile
    # normal quantile), and Table 4 reports 8.5 g for this model.
    # 60 * 88 * exp(0.84 * 0.57) / 1000 = 8.52 g, which rounds to the
    # printed 8.5 g. That reproduces only if omega_Vd is an SD of 0.57 on
    # the log scale, so the "(standard deviation)" label is correct and the
    # value is not a variance.
    #
    # No off-diagonal covariance is reported, so the etas are independent.
    # Eta shrinkage was 38.9% (CL) and 82.5% (Vd); the very high Vd
    # shrinkage follows from the sampling design -- Methods section 4.4 and
    # the Discussion note that no sample earlier than 12 h after the
    # loading dose was available, so individual Vd is only weakly informed.
    # ------------------------------------------------------------------
    etalcl ~ 0.46^2  # Launay 2024 Table 2, IIV (standard deviation) CL = 0.46 (RSE 11.5%, shrinkage 38.9%)
    etalvc ~ 0.57^2  # Launay 2024 Table 2, IIV (standard deviation) Vd = 0.57 (RSE 25.6%, shrinkage 82.5%)

    # ------------------------------------------------------------------
    # Residual unexplained variability. Launay 2024 Table 2 prints this as
    # "Error model / Additive (mg/L)  0.39  (RSE 6.2%)" and the Results
    # text says "one-compartment with additive error model".
    #
    # The printed VALUE (0.39) is used unchanged, but the printed UNIT
    # ("mg/L", i.e. additive in linear concentration space) is falsified by
    # the paper's own Figure 1B, and is therefore treated as a units typo:
    # 0.39 is a standard deviation on the LOG concentration scale, i.e.
    # additive on the log scale, which is exactly a log-normal residual in
    # linear space.
    #
    # Evidence, all internal to the paper:
    #   1. Magnitude. The assay calibration range is 8-150 mg/L and the
    #      target band is 35-80 mg/L (Methods section 4.3). An additive SD
    #      of 0.39 mg/L is ~1% of a typical observation -- roughly two
    #      orders of magnitude tighter than any TDM popPK residual.
    #   2. Figure 1B (observed vs individual predictions) shows scatter of
    #      tens of mg/L at every concentration. Under a 0.39 mg/L additive
    #      residual every point would sit within +/-0.8 mg/L of the identity
    #      line, i.e. visually on it. Digitising panel B gives an SD of
    #      log(observed / individual-predicted) of about 0.33, and the
    #      ratio quantiles are flat from 20 to 120 mg/L (q05/q95 ratios
    #      0.64/1.78 in the 20-30 mg/L bin vs 0.30/1.79 in the 80-120 mg/L
    #      bin). A constant RATIO spread across a six-fold concentration
    #      range is a relative error, not an additive one. The measured
    #      0.33 sits just below the reported 0.39 exactly as expected,
    #      because residuals taken against shrunken individual predictions
    #      understate sigma.
    #   3. The reported eta shrinkages (38.9% on CL, 82.5% on Vd) require a
    #      substantial residual; a ~1% residual on steady-state samples
    #      would make individual CL nearly exactly identified and drive CL
    #      shrinkage toward zero.
    #   4. The paper's own printed accuracy metric, which is the decisive
    #      evidence and needs no digitising. Methods section 4.5 defines
    #      PE (%) = (Cpred - Cobs) / Cobs -- a RELATIVE error -- and calls
    #      the model acceptable if MDAPE (its median absolute value) is
    #      <= 30%. Results report "MPDE and MDAPE for individual
    #      predictions were between -3.3 and 16.9%". An additive 0.39 mg/L
    #      residual against observations of tens of mg/L implies an MDAPE
    #      of ~1%, which cannot be reconciled with a reported 16.9%. A 0.39
    #      SD on the log scale implies a median absolute relative deviation
    #      of exp(0.674 * 0.39) - 1, about 30%, falling into the reported
    #      range once individual predictions absorb part of the residual.
    #
    # Encoded as `lnorm(expSd)` rather than `prop(propSd)` because the
    # paper's own word is "additive" -- additive on the log scale is
    # log-normal, and MONOLIX (used here, Methods section 4.4) names that
    # error model "constant", which authors routinely report as "additive".
    # The two forms differ materially at this magnitude: a 0.39 log-scale
    # SD is a 40.5% CV, not 39%. Same precedent and reasoning as
    # Sano_2023_fesoterodine.R.
    # ------------------------------------------------------------------
    expSd <- 0.39; label("Additive residual SD on the log-transformed concentration scale (log-normal)")  # Launay 2024 Table 2, Error model = 0.39 (RSE 6.2%); printed unit "mg/L" is a typo, see the note above
  })

  model({
    # Reference (median) eGFR for the clearance covariate, Launay 2024
    # Table 2 footnote: "GFRmedian (i.e., 73.90 mL/min/1.73 m^2)".
    crclRef <- 73.90

    # Individual parameters, Launay 2024 Table 2 equation:
    #   CL(i) = CL(pop) * (GFRi / GFRmedian)^(beta^CL_GFR) * exp(eta_CL)
    # Vd carries exponential IIV and no covariate (the Discussion states
    # that weight and BMI were screened but were not significant on Vd).
    cl <- exp(lcl + etalcl) * (CRCL / crclRef)^e_crcl_cl
    vc <- exp(lvc + etalvc)

    kel <- cl / vc

    # One-compartment intravenous disposition. Doses (loading dose and
    # continuous infusion alike) go to `central`.
    d/dt(central) <- -kel * central

    Cc <- central / vc
    Cc ~ lnorm(expSd)
  })
}
