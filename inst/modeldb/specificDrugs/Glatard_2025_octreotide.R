Glatard_2025_octreotide <- function() {
  description <- paste(
    "Joint two-compartment population PK model of octreotide after",
    "subcutaneous administration of the sustained-release depot CAM2029",
    "or immediate-release (IR) octreotide, in healthy participants and",
    "patients with acromegaly. CAM2029 release is described by two",
    "simultaneous first-order absorption processes -- a fraction f_fast",
    "of the dose via the faster process (mean absorption time MAT_fast)",
    "and the remainder via the slow process (MAT_slow) -- whereas",
    "octreotide IR uses a single, much faster first-order process.",
    "Disposition is two-compartmental with first-order elimination.",
    "Body weight is scaled allometrically (exponent 0.75 on CL and Q,",
    "1 on V and Vp, reference 75 kg) and injection in the thigh or",
    "buttock shortens MAT_fast relative to the abdomen. Bioavailability",
    "was fixed at 1 for both formulations. Separate log-scale residual",
    "error terms, each with its own inter-individual variability, apply",
    "to CAM2029 and to octreotide IR observations."
  )
  reference <- paste(
    "Glatard A, Friberg-Hietala S, Keutzer L, Hansson A, Johnsson M,",
    "Tiberg F.",
    "Population Pharmacokinetic Analysis of an Octreotide Depot",
    "(CAM2029) in the Treatment of Acromegaly.",
    "Clin Pharmacokinet 2025;64(7):1079-1092.",
    "doi:10.1007/s40262-025-01522-3.",
    sep = " "
  )
  vignette <- "Glatard_2025_octreotide"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # `logitffast` / `etalogitffast` follow the registered logit transform
  # prefix; `expSdCam2029` / `expSdIr` are formulation-stratified residual
  # SDs (Friberg_2012_voriconazole precedent) with an exponential IIV on
  # each (Chandasana 2024 / Yamamoto 2023 eta-on-epsilon precedent).
  paper_specific_etas <- c(
    "etalogitffast", "etaexpSdCam2029", "etaexpSdIr"
  )
  paper_specific_residual_sds <- c("expSdCam2029", "expSdIr")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometric scaling with reference weight 75 kg; exponent fixed at",
        "0.75 on CL and Q and at 1 on V and Vp (Table 3; Eqs 7-10)."
      ),
      source_name        = "WT"
    ),
    INJSITE_THIGH = list(
      description        = "Subcutaneous injection into the thigh",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (abdomen)",
      notes              = paste(
        "Per-dose-record indicator. Fractional change of -0.351 on MAT_fast",
        "for CAM2029 relative to abdominal injection (Table 3; Eq 11).",
        "Mutually exclusive with INJSITE_BUTTOCK; both 0 = abdomen.",
        "Applies to CAM2029 only -- octreotide IR absorption carries no",
        "injection-site effect in the published model."
      ),
      source_name        = "Injection site"
    ),
    INJSITE_BUTTOCK = list(
      description        = "Subcutaneous injection into the buttock",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (abdomen)",
      notes              = paste(
        "Per-dose-record indicator. Fractional change of -0.527 on MAT_fast",
        "for CAM2029 relative to abdominal injection (Table 3; Eq 11).",
        "Supported by only 1.5% of observations (from 5.1% of participants);",
        "the authors advise interpreting this effect with caution.",
        "Mutually exclusive with INJSITE_THIGH; both 0 = abdomen."
      ),
      source_name        = "Injection site"
    ),
    FORM_OCTREOTIDE_IR = list(
      description        = "Immediate-release octreotide formulation indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (CAM2029 sustained-release depot)",
      notes              = paste(
        "Per-observation-record indicator selecting which of the two",
        "formulation-specific log-scale residual error terms applies",
        "(Table 3: additive RUV on the log scale 0.393 for CAM2029 vs",
        "0.204 for octreotide IR, each with its own exponential IIV).",
        "It does NOT route the dose -- dosing compartment does that:",
        "CAM2029 doses go to `depot` and `depot2`, octreotide IR to",
        "`depot3`. In trial HS-19-664 the same participants received",
        "octreotide IR and, after a washout, CAM2029, so the indicator is",
        "time-varying within a participant."
      ),
      source_name        = "Treatment"
    ),
    PRIOR_OCTREOTIDE = list(
      description        = paste(
        "Record is a pre-first-dose baseline sample in a participant",
        "already receiving octreotide (LAR) or lanreotide therapy"
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (octreotide-naive, or any record at/after the first study dose)",
      notes              = paste(
        "Per-record indicator, NOT a per-subject flag. Gates the estimated",
        "'concentration due to pre-treatment' (`rbase`, 0.433 ng/mL) that",
        "the authors added as a fudge factor for participants who entered",
        "the phase 3 trials on stable octreotide LAR / lanreotide, or who",
        "rolled over from HS-18-633 already on CAM2029 (Methods 2.3.1).",
        "Online Resource 2 shows the 'Pre-treatment' group contributed",
        "exactly 1.0 observation per participant (22 of 46 in HS-18-633,",
        "49 of 95 in HS-19-647), i.e. the single pre-dose baseline sample.",
        "Set to 0 for ordinary simulation; the model prediction from",
        "study dosing is identically 0 at those baseline records, so",
        "adding `rbase` to the prediction and replacing the prediction",
        "with `rbase` are the same thing there."
      ),
      source_name        = "Pre-treatment"
    )
  )

  # Screened in the covariate analysis (Table 2 / Online Resource 4) but NOT
  # retained in the final model, so no point estimate is available. Recorded
  # here to preserve the provenance of the paper's covariate screen.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age", units = "years", type = "continuous",
      notes = "Exploratory covariate on CL and on the absorption parameters; not retained (Sect. 3.3)."
    ),
    SEXF = list(
      description = "Female sex", units = "(binary)", type = "binary",
      notes = "Exploratory covariate on CL and on the absorption parameters; not retained (Sect. 3.3)."
    ),
    BMI = list(
      description = "Body mass index", units = "kg/m^2", type = "continuous",
      notes = paste(
        "Tested on the absorption parameters rather than WT because that was",
        "judged more appropriate for the subcutaneous route; not retained",
        "(Sect. 2.3.2, 3.3)."
      )
    ),
    CRCL = list(
      description = "Creatinine clearance", units = "mL/min", type = "continuous",
      notes = "Exploratory covariate on CL; not retained (Sect. 3.3)."
    ),
    AST = list(
      description = "Aspartate aminotransferase", units = "U/L", type = "continuous",
      notes = "Exploratory covariate on CL; not retained (Sect. 3.3)."
    ),
    TBILI = list(
      description = "Total bilirubin", units = "umol/L", type = "continuous",
      notes = paste(
        "Reached significance as an exponential effect on CL in the SCM",
        "forward step (Online Resource 4, dOFV -15.88) but was deliberately",
        "dropped from the final model because the authors judged the",
        "causality more likely to run the other way -- octreotide affecting",
        "bilirubinaemia (Table 2 footnote c; Sect. 3.2)."
      )
    ),
    DEVICE_AI = list(
      description = "Pre-filled pen (autoinjector) vs pre-filled syringe", units = "(binary)",
      type = "binary",
      notes = paste(
        "Drug delivery system / device design; no evidence of an effect on",
        "PK, consistent with CAM2029 being a true solution (Sect. 3.3, 4.1)."
      )
    )
  )

  compartmentData <- list(
    depot       = list(analyte = "octreotide", units = "mg", specimen = "administration site", verified = TRUE),
    depot2      = list(analyte = "octreotide", units = "mg", specimen = "administration site", verified = TRUE),
    depot3      = list(analyte = "octreotide", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "octreotide", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "octreotide", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 216L,
    n_studies      = 3L,
    n_observations = 4098L,
    age_range      = "18-83 years",
    age_median     = "47 years",
    weight_range   = "50.5-144 kg",
    weight_median  = "79.8 kg",
    bmi_range      = "19.2-50.7 kg/m^2",
    sex_female_pct = 58,
    race_ethnicity = c(White = 96, Asian = 2.3, Black = 0.46, Other = 0.93),
    disease_state  = paste(
      "75 healthy participants (phase 1 trial HS-19-664) and 141 patients",
      "with acromegaly (phase 3 trials HS-18-633 and HS-19-647). Patients",
      "entering the phase 3 trials were on a stable monthly dose of",
      "octreotide LAR or lanreotide autogel for at least 3 months. No",
      "marked PK difference between healthy participants and patients was",
      "identified."
    ),
    dose_range     = paste(
      "CAM2029 10 or 20 mg SC every 4 weeks (single dose and repeated",
      "dosing); octreotide IR 0.25 mg SC every 8 h (4 doses) in the",
      "phase 1 trial."
    ),
    renal_function = "Creatinine clearance 54.3-246 mL/min (median 126); adequate renal function was an eligibility criterion",
    hepatic_function = paste(
      "AST 4.70-119 U/L (median 20.0); total bilirubin 1.37-34.9 umol/L",
      "(median 7.35); adequate hepatic function was an eligibility",
      "criterion, so severe renal or hepatic impairment could not be",
      "evaluated."
    ),
    regions        = "Not reported (EudraCT 2020-002643-35; NCT04076462; NCT04125836)",
    notes          = paste(
      "Baseline characteristics from Table 1 and Online Resource 3.",
      "Injection sites across the 4098 observations: abdomen 88.0%,",
      "thigh 10.3%, buttock 1.5%, missing 0.2%. Assay LLOQ 0.0286 ng/mL",
      "(UPLC-MS/MS). Estimation by FOCEI in NONMEM 7.5 on natural-log",
      "transformed concentrations. 30 participants rolled over from",
      "HS-18-633 into HS-19-647 and were treated as the same individual."
    )
  )

  ini({
    # ---------------------------------------------------------------
    # Allometric exponents -- Table 3 rows "WT on CL, Q" and
    # "WT on V, Vp", both reported (FIX); Eqs 7-10 print the reference
    # weight of 75 kg explicitly.
    # ---------------------------------------------------------------
    e_wt_cl_q  <- fixed(0.750); label("Allometric exponent of WT on CL and Q (unitless)")   # Table 3: 0.750 (FIX); Eqs 7-8
    e_wt_vc_vp <- fixed(1.00);  label("Allometric exponent of WT on V and Vp (unitless)")   # Table 3: 1.00 (FIX); Eqs 9-10

    # ---------------------------------------------------------------
    # Disposition -- Table 3. Values are for the 75 kg reference
    # individual. Total V = 5.98 + 16.1 = 22.08 L, which the paper
    # quotes in Sect. 4.1.
    # ---------------------------------------------------------------
    lcl <- log(9.59); label("Clearance CL (L/h)")                            # Table 3: 9.59 (RSE 2.12%); Eq 7
    lvc <- log(5.98); label("Central volume of distribution V (L)")          # Table 3: 5.98 (RSE 4.40%); Eq 9
    lq  <- log(3.59); label("Inter-compartmental clearance Q (L/h)")         # Table 3: 3.59 (RSE 4.37%); Eq 8
    lvp <- log(16.1); label("Peripheral volume of distribution Vp (L)")      # Table 3: 16.1 (RSE 6.23%); Eq 10

    # ---------------------------------------------------------------
    # Absorption -- Table 3. The paper parameterises each first-order
    # process by its mean absorption time (MAT = 1 / ka). Sect. 4.1
    # corroborates the slow process: the pop-PK slow absorption rate
    # constant is quoted as 0.0022 1/h, and 1 / 459 h = 0.002179 1/h.
    # ---------------------------------------------------------------
    lmatfast <- log(89.2); label("Mean absorption time of the fast CAM2029 release process, MAT_fast (h)") # Table 3: 89.2 (RSE 4.52%); Eq 11
    lmatslow <- log(459);  label("Mean absorption time of the slow CAM2029 release process, MAT_slow (h)") # Table 3: 459 (RSE 3.46%)
    lmatir   <- log(1.28); label("Mean absorption time of octreotide IR, MAT_octreotide_IR (h)")           # Table 3: 1.28 (RSE 3.85%)

    # f_fast is a fraction bounded on (0, 1) and is carried on the logit
    # scale so that individual values stay in range. logit(0.307) =
    # log(0.307 / (1 - 0.307)) = -0.8137.
    logitffast <- log(0.307 / (1 - 0.307)); label("Logit of the CAM2029 dose fraction absorbed by the fast process, f_fast (logit units)") # Table 3: f_fast = 0.307 (RSE 4.05%)

    # Bioavailability was fixed at 1 for BOTH CAM2029 and octreotide IR
    # to let the IIV in the extent of absorption be estimated
    # (Sect. 2.3.1). Sect. 4.2 confirms the near-complete bioavailability.
    lfdepot <- fixed(log(1)); label("Bioavailability F of CAM2029 and octreotide IR (unitless)") # Sect. 2.3.1: F fixed to 1 for both formulations

    # ---------------------------------------------------------------
    # Injection-site effect on MAT_fast -- Table 3 / Eq 11, entered as a
    # fractional change relative to the abdomen (Eq 5).
    # ---------------------------------------------------------------
    e_injsite_thigh_matfast   <- -0.351; label("Fractional change in MAT_fast for injection in the thigh vs abdomen (unitless)")   # Table 3: -0.351 (RSE 14.6%); Eq 11
    e_injsite_buttock_matfast <- -0.527; label("Fractional change in MAT_fast for injection in the buttock vs abdomen (unitless)") # Table 3: -0.527 (RSE 18.6%); Eq 11

    # ---------------------------------------------------------------
    # Pre-treatment baseline concentration -- the "fudge factor" of
    # Sect. 2.3.1, describing the octreotide concentration measured
    # before the first study dose in participants already on octreotide
    # LAR / lanreotide. Gated by PRIOR_OCTREOTIDE.
    # ---------------------------------------------------------------
    lrbase <- log(0.433); label("Plasma octreotide concentration due to pre-treatment (ng/mL)") # Table 3: 0.433 (RSE 12.6%)

    # ---------------------------------------------------------------
    # Inter-individual variability -- Table 3. Every IIV row is reported
    # in the "CV" unit column, so omega = CV on the exponential (log)
    # scale and omega^2 is the variance below. The single exception is
    # IIV f_fast, whose unit cell is blank -- consistent with an omega on
    # the logit scale, which has no CV interpretation (see the vignette
    # Assumptions section for the Table 4 check that supports this).
    # ---------------------------------------------------------------
    etalcl         ~ 0.222^2; # Table 3: IIV CL, CV 0.222 (RSE 6.05%, shrinkage 12.1%)
    etalmatfast    ~ 0.194^2; # Table 3: IIV MAT_fast CAM2029, CV 0.194 (RSE 23.0%, shrinkage 47.8%)
    etalmatslow    ~ 0.482^2; # Table 3: IIV MAT_slow CAM2029, CV 0.482 (RSE 6.79%, shrinkage 9.15%)
    etalogitffast  ~ 0.373^2; # Table 3: IIV f_fast CAM2029, 0.373 (RSE 16.3%, shrinkage 43.3%) -- unit cell blank, taken as the logit-scale SD
    etalrbase      ~ 0.720^2; # Table 3: IIV concentration due to pre-treatment, CV 0.720 (RSE 12.2%, shrinkage 12.4%)

    # IIV on the residual error itself (the eta on epsilon of Eq 2),
    # included by the authors to absorb between-individual differences in
    # RUV and outliers.
    etaexpSdCam2029 ~ 0.410^2; # Table 3: IIV RUV CAM2029, CV 0.410 (RSE 7.00%, shrinkage 3.75%)
    etaexpSdIr      ~ 0.256^2; # Table 3: IIV RUV octreotide IR, CV 0.256 (RSE 15.5%, shrinkage 18.3%)

    # ---------------------------------------------------------------
    # Residual unexplained variability -- Eq 2 is additive on the
    # natural-log scale, i.e. log-normal residual error on the linear
    # scale, hence `lnorm()`. Two separate terms, one per formulation.
    # ---------------------------------------------------------------
    expSdCam2029 <- 0.393; label("Log-scale residual SD for CAM2029 observations (typical value, before the exp(eta) modifier)")       # Table 3: additive RUV log scale CAM2029 0.393 (RSE 3.6%, shrinkage 1.42%)
    expSdIr      <- 0.204; label("Log-scale residual SD for octreotide IR observations (typical value, before the exp(eta) modifier)") # Table 3: additive RUV log scale octreotide IR 0.204 (RSE 4.3%, shrinkage 1.20%)
  })

  model({
    # ---- 1. Covariate terms -------------------------------------------------
    # Eqs 7-10: allometric scaling on a reference weight of 75 kg.
    wt_cl_q  <- (WT / 75)^e_wt_cl_q
    wt_vc_vp <- (WT / 75)^e_wt_vc_vp

    # Eq 11: fractional change on MAT_fast, abdomen as reference. The
    # indicators are mutually exclusive, so this additive form reproduces
    # the paper's piecewise selection exactly.
    site_matfast <- 1 +
      e_injsite_thigh_matfast   * INJSITE_THIGH +
      e_injsite_buttock_matfast * INJSITE_BUTTOCK

    # ---- 2. Individual parameters -------------------------------------------
    cl <- exp(lcl + etalcl) * wt_cl_q
    vc <- exp(lvc) * wt_vc_vp
    q  <- exp(lq)  * wt_cl_q
    vp <- exp(lvp) * wt_vc_vp

    matfast <- exp(lmatfast + etalmatfast) * site_matfast
    matslow <- exp(lmatslow + etalmatslow)
    matir   <- exp(lmatir)
    ffast   <- expit(logitffast + etalogitffast)
    rbase   <- exp(lrbase + etalrbase)
    fdepot  <- exp(lfdepot)

    # First-order rate constants are the reciprocals of the mean
    # absorption times.
    kafast <- 1 / matfast
    kaslow <- 1 / matslow
    kair   <- 1 / matir

    # ---- 3. Micro-constants -------------------------------------------------
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ---- 4. ODE system ------------------------------------------------------
    # `depot` and `depot2` are the fast and slow CAM2029 release routes
    # (Fig. 1, left); `depot3` is the octreotide IR route (Fig. 1, top).
    d/dt(depot)  <- -kafast * depot
    d/dt(depot2) <- -kaslow * depot2
    d/dt(depot3) <- -kair   * depot3
    d/dt(central) <- kafast * depot + kaslow * depot2 + kair * depot3 -
      kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # ---- 5. Bioavailability -------------------------------------------------
    # A CAM2029 injection is entered as the SAME dose amount into both
    # `depot` and `depot2`; f() then partitions it f_fast / (1 - f_fast)
    # so the total absorbed mass equals the administered dose (F = 1).
    # An octreotide IR injection is entered into `depot3` alone.
    f(depot)  <- ffast * fdepot
    f(depot2) <- (1 - ffast) * fdepot
    f(depot3) <- fdepot

    # ---- 6. Observation and error model -------------------------------------
    # central is in mg and vc in L, so central / vc is mg/L; 1 mg/L =
    # 1000 ng/mL, the unit Table 3 and Table 4 report.
    Cc <- central / vc * 1000 + rbase * PRIOR_OCTREOTIDE

    # Eq 2 with the formulation-specific sigma of Table 3.
    expSd_i <- expSdCam2029 * exp(etaexpSdCam2029) * (1 - FORM_OCTREOTIDE_IR) +
      expSdIr * exp(etaexpSdIr) * FORM_OCTREOTIDE_IR
    Cc ~ lnorm(expSd_i)
  })
}
