Jang_2023_zaltoprofen <- function() {
  description <- paste(
    "Two-compartment population PK model with first-order absorption and no",
    "lag time for oral zaltoprofen in healthy Korean male adults (Jang 2023;",
    "n = 26; single 80 mg oral dose). Covariate effects on apparent clearance:",
    "Cockcroft-Gault creatinine clearance as a power function centred on the",
    "cohort median 104.38 mL/min (exponent 0.48) and serum albumin as a power",
    "function centred on the cohort median 4.90 g/dL (exponent -1.83).",
    "Inter-individual variability was estimated on Ka, CL/F and V2/F only;",
    "residual error is proportional. CAUTION: the published parameter set",
    "reproduced here does NOT reproduce the source paper's own Figures 1-2 or",
    "its quoted steady-state concentrations - it predicts roughly 7-fold lower",
    "exposure and a roughly 5-fold earlier Tmax than the observed data plotted",
    "in the same paper. The values are shipped exactly as printed in",
    "supplementary Table S3; see the validation vignette's Errata section for",
    "the full quantitative demonstration before using this model."
  )
  reference <- paste(
    "Jang JH, Jeong SH, Lee YB (2023).",
    "Population Pharmacokinetic Modeling of Zaltoprofen in Healthy Adults:",
    "Exploring the Dosage Regimen.",
    "Pharmaceuticals 16(2):161.",
    "doi:10.3390/ph16020161.",
    "All parameter values are from supplementary Table S3;",
    "the structural equations are main-text Eqs. (1)-(5).",
    sep = " "
  )
  vignette <- "Jang_2023_zaltoprofen"
  units    <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    # Oral tablet; the source study administered a single 80 mg zaltoprofen
    # tablet (Methods 4.2 / ref [19] bioequivalence study).
    depot       = list(analyte = "zaltoprofen", units = "mg", specimen = "administration site", verified = TRUE),
    # Plasma zaltoprofen was the measured analyte (Methods 4.2 and Table S6:
    # "plasma zaltoprofen concentration-time curves").
    central     = list(analyte = "zaltoprofen", units = "mg", specimen = "plasma", verified = TRUE),
    # V2/F is named only as the peripheral-compartment volume (main text after
    # Eq. (5): "V and CL2 indicate volume and clearance in the central and
    # peripheral compartments"). No biological matrix is stated for it, so the
    # specimen assignment is conventional rather than source-verified.
    peripheral1 = list(analyte = "zaltoprofen", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    CRCL = list(
      description        = paste(
        "Creatinine clearance calculated with the Cockcroft-Gault equation,",
        "reported as raw mL/min and NOT normalised to 1.73 m^2 body surface",
        "area."
      ),
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject (calculated from the pre-dose 0 h blank serum",
        "sample; Supplementary Materials Subsection 1). Power effect on CL/F:",
        "CL/F = tvCL/F * (CrCL / mCrCL)^dCL/FdCrCL with dCL/FdCrCL = 0.48 per",
        "main-text Eq. (2) and supplementary Table S3. The reference mCrCL is",
        "the cohort MEDIAN 104.38 mL/min, stated repeatedly in the main text",
        "(Results 2.4 and the Figure 2-5 captions: 'the group with median",
        "values of CrCL and albumin of 104.38 mL/min and 4.90 g/dL').",
        "Table S5 reports the cohort MEAN as 107.53 +/- 17.28 mL/min, which is",
        "a different statistic and is NOT the centring constant.",
        "The Cockcroft-Gault form used is given in Supplementary Materials",
        "Subsection 1: CrCL = (140 - age) * body weight (kg) /",
        "(serum creatinine (mg/dL) * 72); the cohort was entirely male so no",
        "0.85 female multiplier applies. Note the model was fitted only over a",
        "narrow healthy range (mean 107.53, SD 17.28), while the paper's own",
        "dosing simulations (Results 2.4) extrapolate the term to 80 and",
        "130 mL/min. MDRD-formula GFR (101.34 +/- 11.53 mL/min, Table S5) was",
        "measured on the same subjects and screened separately (Figure S3D)",
        "but was not the retained covariate; it is not encoded here."
      ),
      source_name        = "CrCL"
    ),
    ALB = list(
      description        = "Serum albumin measured in the pre-dose (0 h) blank serum sample by reflectance spectrophotometry.",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject. Source reports albumin in US-convention g/dL",
        "(Table S5: 4.92 +/- 0.18 g/dL; centring median 4.90 g/dL), so",
        "model() applies the inline conversion alb_gdL <- ALB * 0.1 required",
        "by the ALB register entry, whose canonical unit is SI g/L. The",
        "reference 4.90 g/dL therefore equals 49.0 g/L. Power effect on CL/F:",
        "CL/F = tvCL/F * (Albumin / mAlbumin)^dCL/FdAlbumin with",
        "dCL/FdAlbumin = -1.83 per main-text Eq. (2) and supplementary",
        "Table S3. The negative exponent means HIGHER albumin gives LOWER",
        "apparent clearance, which the paper attributes to zaltoprofen's high",
        "plasma protein binding restricting the unbound fraction available for",
        "elimination. mAlbumin = 4.90 g/dL is the cohort median stated in the",
        "main text (Results 2.4 and the Figure 2-5 captions); the Table S5",
        "mean of 4.92 g/dL is a different statistic and is NOT the centring",
        "constant. The paper's own simulations extrapolate the term to 3.5 and",
        "5.5 g/dL (i.e. 35 and 55 g/L)."
      ),
      source_name        = "Albumin"
    )
  )

  # Covariates the source screened but did NOT retain in the final model. These
  # are documentation only -- they are deliberately absent from model(). Keeping
  # them here preserves the provenance of the paper's covariate screen without
  # raising a "declared but not referenced" convention warning.
  covariatesDataExcluded <- list(
    CYP2C9_S3_COUNT = list(
      description        = "Count of CYP2C9*3 reduced-function alleles (0 for *1/*1, 1 for *1/*3).",
      units              = "(count)",
      type               = "count",
      reference_category = "0 (CYP2C9*1/*1, wild type)",
      notes              = paste(
        "Only the *1 and *3 alleles were detected in this cohort",
        "(Supplementary Materials Subsection 2), so the genotype covariate is",
        "the two-level contrast *1/*1 versus *1/*3. Tested on Ka, V2/F and",
        "CL/F (Table S2 models 2-4, dOFV -1.544, -3.722 and -3.845) and again",
        "alongside the retained covariates (models 9 and 11). It met the",
        "forward-selection criterion on CL/F but failed backward elimination",
        "at p < 0.01 and produced no visual improvement in the goodness-of-fit",
        "plots (Results 2.1), so it is absent from the final model. Of the",
        "seven NCA parameters compared between genotypes in Figure S1, only",
        "Tmax differed significantly (p = 0.024); Cmax, AUC0-t, AUCinf, V/F,",
        "CL/F and T1/2 did not. No coefficient is reported for it anywhere in",
        "the paper, so it cannot be encoded even optionally."
      ),
      source_name        = "CYP2C9 genotype"
    ),
    BSA = list(
      description        = "Body surface area computed with the Mosteller formula.",
      units              = "m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Cohort 1.76 +/- 0.12 m^2 (Table S5); computed as",
        "sqrt(height (cm) * weight (kg) / 3600) (Supplementary Materials",
        "Subsection 1). Passed the 30% linear-regression screen against",
        "individual CL/F (Figure S3E) and was tested formally on CL/F",
        "(Table S2 model 7, dOFV -2.076), but did not reach the",
        "forward-selection threshold and was not retained. No coefficient is",
        "reported."
      ),
      source_name        = "BSA"
    ),
    BMI = list(
      description        = "Body mass index (Kaup index).",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Cohort 21.74 +/- 2.60 kg/m^2 (Table S5). Screened against individual",
        "CL/F and V/F (Figures S4D and S4F) with a linear correlation below",
        "the 30% cutoff, so it was excluded before the formal stepwise",
        "covariate analysis and never entered Table S2. No coefficient is",
        "reported."
      ),
      source_name        = "BMI"
    ),
    ALT = list(
      description        = "Alanine aminotransferase.",
      units              = "IU/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Cohort 18.23 +/- 8.28 IU/L (Table S5). Screened against individual",
        "CL/F (Figure S4A) with a linear correlation below the 30% cutoff and",
        "excluded before the stepwise covariate analysis. No coefficient is",
        "reported."
      ),
      source_name        = "ALT"
    ),
    ALP = list(
      description        = "Alkaline phosphatase.",
      units              = "IU/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Cohort 73.23 +/- 15.30 IU/L (Table S5). Screened against individual",
        "CL/F (Figure S4B) with a linear correlation below the 30% cutoff and",
        "excluded before the stepwise covariate analysis. No coefficient is",
        "reported."
      ),
      source_name        = "ALP"
    ),
    BUN = list(
      description        = "Blood urea nitrogen.",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Cohort 12.93 +/- 2.29 mg/dL (Table S5). Screened against individual",
        "CL/F (Figure S4C) with a linear correlation below the 30% cutoff and",
        "excluded before the stepwise covariate analysis. No coefficient is",
        "reported."
      ),
      source_name        = "BUN"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 26L,
    n_studies      = 1L,
    age_mean       = "23.19 +/- 2.26 years (mean +/- SD; Table S5)",
    weight_mean    = "64.73 +/- 8.08 kg (mean +/- SD; Table S5)",
    height_mean    = "172.64 +/- 5.95 cm (mean +/- SD; Table S5)",
    sex_female_pct = 0,
    race_ethnicity = c(Asian = 100),
    disease_state  = paste(
      "Healthy adult volunteers. The modelling dataset is the",
      "pharmacokinetic arm of a previously reported bioequivalence study",
      "(Methods 4.2, citing ref [19]), approved by the Institutional Review",
      "Board of the Institute of Bioequivalence and Bridging Study, Chonnam",
      "National University (Approval No. 060118, 2005-11-25). All",
      "biochemical parameters in Table S5 are within normal limits."
    ),
    dose_range     = paste(
      "Single 80 mg oral dose of zaltoprofen. No multiple-dose data were",
      "used to fit the model; the paper's 80 mg q24h and q8h multiple-dose",
      "profiles (Figures 3-5) are simulations from the fitted model, not",
      "observations."
    ),
    renal_function = paste(
      "Cockcroft-Gault CrCL 107.53 +/- 17.28 mL/min and MDRD GFR",
      "101.34 +/- 11.53 mL/min (Table S5); serum creatinine",
      "0.98 +/- 0.09 mg/dL. The centring median used by the model is",
      "104.38 mL/min. No subject had impaired renal function."
    ),
    genotype       = paste(
      "CYP2C9 genotyped by PCR-RFLP on exon 7 (A1075->C); only the *1 and *3",
      "alleles were detected, so subjects are either *1/*1 or *1/*3",
      "(Supplementary Materials Subsection 2). Genotype concordance against",
      "direct sequencing was 100%. Genotype was not retained as a covariate."
    ),
    regions        = "Republic of Korea (Gwangju; Chonnam National University).",
    notes          = paste(
      "Baseline demographics from supplementary Table S5. The model was fitted",
      "in Phoenix NLME 8.3 by first-order conditional estimation with extended",
      "least squares and eta-epsilon interaction (Methods 4.2), and evaluated",
      "with goodness-of-fit plots, a 1000-replicate non-parametric bootstrap",
      "(Table S4), a 100-simulation VPC (Figures S7-S8) and NPDE (Figure S6).",
      "External validation used published mean profiles from three prior",
      "studies (Table S6: Kang 2006 80 mg n = 26; Lee 2006 160 mg n = 24;",
      "Li 2011 80/160/240 mg single dose and 80 mg q8h multiple dose n = 12).",
      "Those external data are observed means only and did not inform the fit."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # All values are supplementary Table S3, "Estimate" column. Typical
    # values refer to the covariate reference subject, i.e. the cohort
    # medians CrCL = 104.38 mL/min and albumin = 4.90 g/dL, because both
    # covariate terms in Eq. (2) are median-normalised ratios that equal 1
    # at those values.
    #
    # SEE THE VIGNETTE ERRATA: this printed parameter set does not reproduce
    # the paper's own Figures 1-2. It is shipped verbatim as printed, per the
    # standing rule that printed values have authority over figures.
    # ------------------------------------------------------------------

    lka <- log(1.73);  label("Absorption rate constant (Ka, 1/h)")                       # Table S3, tv Ka (SE 0.15, RSE 8.56%)
    lvc <- log(4.88);  label("Apparent central volume of distribution (V/F, L)")         # Table S3, tv V/F (SE 1.28, RSE 26.32%)
    lcl <- log(43.70); label("Apparent clearance (CL/F, L/h) at CrCL 104.38 mL/min and albumin 4.90 g/dL")  # Table S3, tv CL/F (SE 2.12, RSE 4.85%)
    lvp <- log(40.56); label("Apparent peripheral volume of distribution (V2/F, L)")     # Table S3, tv V2/F (SE 7.02, RSE 17.32%)
    lq  <- log(5.61);  label("Apparent intercompartmental clearance (CL2/F, L/h)")       # Table S3, tv CL2/F (SE 0.63, RSE 11.16%)

    # Covariate effects on CL/F. Both are power exponents on a
    # median-normalised ratio (main-text Eq. (2)).
    e_crcl_cl <- 0.48;  label("Power exponent of CrCL on CL/F (unitless)")               # Table S3, dCL/FdCrCL (SE 0.17, RSE 35.09%)
    e_alb_cl  <- -1.83; label("Power exponent of serum albumin on CL/F (unitless)")      # Table S3, dCL/FdAlbumin (SE 0.64, RSE 35.18%)

    # Inter-individual variability. Eqs. (1)-(5) place exponential IIV on Ka,
    # CL/F and V2/F only; V/F and CL2/F carry no eta, which matches Table S1
    # step 02-04-06 ("Remove IIV V/F, CL2/F", the selected IIV model).
    #
    # Table S3's "omega^2" column IS a variance, and its "IIV (%)" column
    # carries two more significant figures of the same quantity. The two
    # columns are related by omega = IIV(%) / 100 (the NONMEM/Phoenix
    # approximate-CV convention), which reproduces every rounded omega^2 in
    # the table:
    #   Ka   : 0.4063^2 = 0.16508 -> prints as 0.17
    #   V2/F : 0.4855^2 = 0.23571 -> prints as 0.24
    #   CL/F : 0.0980^2 = 0.00960 -> prints as 0.01
    # The alternative log-normal reading omega = sqrt(log(1 + CV^2)) does NOT
    # reproduce them (it gives 0.15, 0.21 and 0.01), so it is ruled out. The
    # higher-precision back-calculated variances are used below; the vignette
    # re-derives this in the "iiv-precision-check" chunk.
    etalka ~ 0.16508   # Table S3, omega^2 Ka   = 0.17 (RSE 33.95%, shrinkage 0.09, IIV 40.63%)
    etalcl ~ 0.009604  # Table S3, omega^2 CL/F = 0.01 (RSE 46.60%, shrinkage 0.34, IIV  9.80%)
    etalvp ~ 0.23571   # Table S3, omega^2 V2/F = 0.24 (RSE 32.31%, shrinkage 0.21, IIV 48.55%)

    # Residual error: proportional only. Table S1 selected model 02-04
    # ("Proportional", d-2LL = -474.97 against the additive model 02); the
    # mixed model 02-06 gained nothing (d-2LL = -474.95 on one more
    # parameter).
    propSd <- 0.38; label("Proportional residual error (fraction)")                      # Table S3, epsilon (SE 0.02, RSE 6.03%)
  })

  model({
    # 1. Derived covariate terms.
    # The ALB register stores serum albumin in SI g/L; Jang 2023 calibrated
    # the exponent against US-convention g/dL, so convert before the ratio.
    alb_gdL <- ALB * 0.1  # SI g/L -> US-convention g/dL (49.0 g/L == 4.90 g/dL, the reference)

    # 2. Individual parameters, exactly as main-text Eqs. (1)-(5).
    ka <- exp(lka + etalka)                                                    # Eq. (5)
    vc <- exp(lvc)                                                             # Eq. (1); no IIV
    cl <- exp(lcl + etalcl) * (CRCL / 104.38)^e_crcl_cl * (alb_gdL / 4.90)^e_alb_cl  # Eq. (2)
    vp <- exp(lvp + etalvp)                                                    # Eq. (3)
    q  <- exp(lq)                                                              # Eq. (4); no IIV

    # 3. ODE system. Absorption is first order with no lag time (Table S1
    # selected model 02 over 02-01 "Add lag-time", d-2LL = +70.58).
    #
    # The distribution terms are written in MACRO form (q / vc, q / vp) rather
    # than as k12 / k21 micro-constants. rxSolve()'s default useLinCmt = TRUE
    # detector matches on canonical macro names; a two-compartment body written
    # with k12 / k21 is silently collapsed to a one-compartment closed form
    # (peripheral1 dropped, terminal half-life wrong) while total AUC still
    # equals Dose / CL, so the usual sanity check does not catch it. The macro
    # form keeps the default solve correct for every downstream user.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - (cl / vc) * central - (q / vc) * central + (q / vp) * peripheral1
    d/dt(peripheral1) <-  (q / vc) * central - (q / vp) * peripheral1

    # 4. Observation. Volumes and clearances are apparent (X/F), so no separate
    # bioavailability term is estimable and none is reported. Dose in mg over
    # volume in L gives mg/L, which is numerically ug/mL as reported by the
    # paper (e.g. the 0.43 ug/mL steady-state mean in Results 2.4).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
