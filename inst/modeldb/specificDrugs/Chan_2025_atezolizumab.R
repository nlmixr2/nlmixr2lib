Chan_2025_atezolizumab <- function() {
  description <- paste0(
    "Two-compartment population PK model for atezolizumab (anti-PD-L1 ",
    "IgG1) describing BOTH intravenous and subcutaneous administration in ",
    "adults with locally advanced or metastatic non-small cell lung ",
    "cancer (Chan 2025, n = 435, 3100 serum concentrations, the phase Ib ",
    "+ phase III study IMscin001, NCT03735121). This is an EXTENSION ",
    "model: the authors held every disposition parameter (CL, Vc, Vp, Q), ",
    "every disposition covariate effect and the entire CL/Vc/Vp ",
    "between-subject covariance block FIXED at the values of the ",
    "historical intravenous atezolizumab population PK model of Stroh ",
    "2017, on the argument that disposition is intrinsic to the molecule ",
    "and does not change with route, and estimated ONLY the subcutaneous ",
    "absorption layer: a first-order absorption rate constant KA = 0.304 ",
    "1/day and a bioavailability F1 = 71.8 percent, plus their ",
    "between-subject variances, two absorption covariate effects and the ",
    "residual error. Bioavailability is carried on the LOGIT scale and ",
    "the hemoglobin covariate MULTIPLIES that logit rather than adding ",
    "to it, exactly as written in the supplementary NONMEM control ",
    "stream. Subcutaneous doses enter the depot compartment (F1 applies); ",
    "intravenous doses enter the central compartment directly (F1 does ",
    "not apply). Six companion exposure-response models in the ",
    "Chan_2025_atezolizumab_* family."
  )
  reference <- paste(
    "Chan P, Liu SN, Gosselin N, Sauve Z, Marchand M, Lin A,",
    "Herraez-Baranda L, Zanghi J, Shearer-Kang E, Liu X, Wu B, Chanu P.",
    "Population pharmacokinetics and exposure-response of subcutaneous",
    "atezolizumab in patients with non-small cell lung cancer.",
    "CPT Pharmacometrics Syst Pharmacol. 2025;14(4):726-737.",
    "doi:10.1002/psp4.13310.",
    "Structural values are transcribed from Table 2 and from the NONMEM",
    "control stream supplied as Supporting Information file s002.CTL.",
    "The FIXED disposition parameters, their covariate effects and the",
    "CL/Vc/Vp covariance block originate in the historical intravenous",
    "model of Stroh M, Winter H, Marchand M, Claret L, Eppler S, Ruppel J,",
    "et al. Clinical pharmacokinetics and pharmacodynamics of atezolizumab",
    "in metastatic urothelial carcinoma. Clin Pharmacol Ther.",
    "2017;102(2):305-312. doi:10.1002/cpt.587;",
    "Chan 2025 reprints every one of those values in its own Table 2 and",
    "control stream, so nothing here is taken from an off-disk source.",
    sep = " "
  )
  vignette <- "Chan_2025_atezolizumab_sc_nsclc"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power effects on CL (exponent 0.808) and Vc (exponent 0.559),",
        "both normalized to a reference weight of 77 kg. Both are FIXED",
        "from the historical intravenous model (Chan 2025 Table 2 rows",
        "'Bodyweight on CL (BWT/77 in kg)' and 'Bodyweight on Vc (BWT/77",
        "in kg)', both marked FIX). The control stream imputes a missing",
        "weight to 68.1 kg -- note that 68.1 kg is the imputation value,",
        "NOT the normalization reference, which is 77 kg; do not conflate",
        "them. Body weight was also screened on KA (control-stream block",
        "KABWT, normalized to 68.1 kg) but its exponent is $THETA 0 FIX,",
        "so it contributes a factor of exactly 1 and is deliberately",
        "absent from model(). Cohort medians span 65.4-73.2 kg across the",
        "five IMscin001 cohorts (Chan 2025 Table 1); the phase III",
        "subcutaneous arm median is 67.8 kg [30.0, 117]."
      ),
      source_name        = "BWT"
    ),
    ALB = list(
      description        = "Baseline serum albumin",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "The only covariate acting on all three of CL, Vc and KA. Power",
        "effects normalized to 40 g/L throughout: CL exponent -1.12",
        "(FIXED), Vc exponent -0.350 (FIXED), KA exponent +0.795",
        "(ESTIMATED in Chan 2025, RSE 33.1 percent). Higher albumin",
        "therefore lowers clearance and central volume but speeds",
        "subcutaneous absorption, shortening Tmax. Reported in g/L (SI):",
        "Chan 2025 Table 1 gives cohort medians of 39.0-41.2 g/L and the",
        "covariate reference of 40 confirms the g/L scale. Albumin is the",
        "single largest covariate effect on Cycle-1 AUC0-21d: Chan 2025",
        "Results reports every covariate effect falling inside 80-125",
        "percent EXCEPT the 5th-percentile albumin value, which gives",
        "77.9 percent at 26 g/L. The control stream sets ALB to a missing",
        "sentinel when the recorded value is below 5 g/L (a unit-error",
        "guard) and then imputes the reference 40 g/L for CL and Vc, but",
        "leaves the KA factor at 1 for a missing albumin."
      ),
      source_name        = "ALBU (recoded to ALB in $PK)"
    ),
    TUMSZ = list(
      description        = "Baseline tumor burden (sum of longest diameters of target lesions)",
      units              = "mm",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power effect on CL only, exponent 0.125, normalized to 63 mm,",
        "FIXED from the historical intravenous model (Chan 2025 Table 2",
        "row 'Tumor burden on CL (Tumor burden/63mm)', marked FIX;",
        "control-stream block CLBSLD, source column BSLD). The control",
        "stream imputes a missing tumor burden to 74 mm, which is NOT the",
        "63 mm normalization reference. Cohort medians span 56.0-87.0 mm",
        "(Chan 2025 Table 1); the phase III subcutaneous arm median is",
        "79.5 mm [10.0, 319]."
      ),
      source_name        = "BSLD"
    ),
    ADA_POS = list(
      description        = "Treatment-emergent anti-drug-antibody status (1 = ADA-positive, 0 = ADA-negative)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ADA-negative)",
      notes              = paste(
        "FRACTIONAL, not exponential: the effect enters as",
        "cl *= (1 + 0.159 * ADA_POS), i.e. a 15.9 percent higher typical",
        "clearance in ADA-positive patients. Chan 2025 Table 2 labels the",
        "row 'ADA status on CL (additive effect for positive ADA)', which",
        "is ambiguous between an additive-on-log-scale and an",
        "additive-on-fractional-scale reading; the control stream settles",
        "it -- 'IF(ATAG.EQ.1) CLATAG=( 1 + THETA(6))' is unambiguously the",
        "(1 + theta) fractional form, so exp(0.159) = 1.172 is NOT the",
        "right factor. FIXED from the historical intravenous model. ADA",
        "was deliberately NOT re-evaluated on the subcutaneous absorption",
        "parameters (Chan 2025 Discussion) because ADA incidence was",
        "similar between the intravenous and subcutaneous arms; when it",
        "was tested in the sensitivity analysis it was not significant.",
        "A missing ADA status is IMPUTED by the source analysis rather",
        "than carried as its own level ('ATAG=ATAGIM ; ATAG imputed for",
        "missing' in $PK), so no ADA_MISSING indicator is needed here",
        "even though 10.6 percent of the phase III subcutaneous arm has a",
        "missing status (Chan 2025 Table 1)."
      ),
      source_name        = "ATAG (imputed from ATAGIM)"
    ),
    SEXF = list(
      description        = "Biological sex indicator (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = paste(
        "FRACTIONAL effects on both volumes: vc *= (1 - 0.129 * SEXF) and",
        "vp *= (1 - 0.272 * SEXF), i.e. females have a 12.9 percent lower",
        "central volume and a 27.2 percent lower peripheral volume. Both",
        "FIXED from the historical intravenous model. As with ADA_POS,",
        "Chan 2025 Table 2 calls these 'additive effect for female' and",
        "the control stream disambiguates: 'IF(SEX.EQ.2) V2SEX=( 1 +",
        "THETA(11))' with THETA(11) = -0.129. The source SEX column codes",
        "1 = male (the reference, flagged 'MOST COMMON') and 2 = female,",
        "so SEXF = SEX - 1 under that encoding and the sign of the",
        "coefficient is preserved. Sex has NO effect on CL in this model.",
        "The phase III subcutaneous arm is 29.3 percent female (Chan 2025",
        "Table 1)."
      ),
      source_name        = "SEX (1 = male, 2 = female)"
    ),
    HGB = list(
      description        = "Baseline hemoglobin",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "The only covariate on bioavailability, and the one whose",
        "functional form most needs the control stream to read correctly.",
        "It acts on the LOGIT of F1 and MULTIPLIES that logit rather than",
        "adding to it:",
        "F1 = expit(logit(0.718) * (HGB/123)^1.76 + eta). Chan 2025",
        "Table 2 records this only as 'Hemoglobin on F1 On logit scale",
        "(HGB/123 in g/L) 1.76', which does not by itself reveal that",
        "1.76 is a power exponent on a MULTIPLIER of the logit; the",
        "control-stream lines 'LF1 = LOG(THETA(14)/(1-THETA(14)))' and",
        "'F1 = EXP(LF1*F1HGB+...)/(1+EXP(LF1*F1HGB+...))' with",
        "'F1HGB=((HGB/123)**THETA(19))' establish it. At the reference",
        "123 g/L the multiplier is exactly 1 and F1 returns the typical",
        "0.718. Because logit(0.718) is POSITIVE, higher hemoglobin",
        "raises the logit and hence F1, matching the Chan 2025 Results",
        "statement that 'increasing hemoglobin was associated with a",
        "higher F1'. Reported in g/L: cohort medians 119-124 g/L (Chan",
        "2025 Table 1). The control stream sets HGB to a missing sentinel",
        "when the recorded value exceeds 1000 g/L (a unit-error guard)",
        "and then leaves the multiplier at 1. Hemoglobin was also",
        "screened on KA (control-stream block KAHGB) but its exponent is",
        "$THETA 0 FIX, so it is deliberately absent from model()."
      ),
      source_name        = "HGBU (recoded to HGB in $PK)"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = paste(
        "Screened on KA in the source control stream (block KAAGE,",
        "'KAAGE=((AGEN/64)**THETA(15))', reference 64 years) but its",
        "exponent is $THETA 0 FIX in the final model, so it contributes a",
        "factor of exactly 1. The control stream's second header line",
        "reads ';; 2. Description: no AGE KA', confirming that dropping",
        "the age-on-KA term is the deliberate final structure. Documented",
        "for provenance only; deliberately absent from model()."
      )
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened on the logit of F1 in the source control stream (block",
        "F1RACE, 'IF(RACEN.EQ.2) F1RACE=(THETA(21))') but THETA(21) is",
        "$THETA 0 FIX in the final model, so the term vanishes.",
        "Documented for provenance only; deliberately absent from",
        "model()."
      )
    ),
    RACE_HISPANIC = list(
      description = "Hispanic / Latino ethnicity indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened on the logit of F1 in the source control stream (block",
        "F1ETHN, 'IF(ETHN.EQ.1) F1ETHN=(THETA(20))') but THETA(20) is",
        "$THETA 0 FIX in the final model, so the term vanishes. The",
        "control stream's first header line reads ';; 1. Based on:",
        "noETHNF1', naming the parent run in which the ethnicity-on-F1",
        "term was already removed. Note that ethnicity DOES survive in a",
        "companion exposure-response model",
        "(Chan_2025_atezolizumab_isr). Documented for provenance only;",
        "deliberately absent from model()."
      )
    ),
    FORM_ATEZOLIZUMAB_COFORMULATED = list(
      description = "Subcutaneous formulation (co-formulated with rHuPH20 vs co-mixed at the site)",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Tested as a candidate covariate on the absorption parameters and",
        "NOT retained: Chan 2025 Results states 'the SC formulation was",
        "not identified as a statistically significant covariate on",
        "absorption parameters'. The phase Ib cohorts received",
        "atezolizumab co-mixed with rHuPH20 and the phase III",
        "subcutaneous cohort received a ready-to-use co-formulation, so",
        "the contrast was estimable. Documented for provenance only; no",
        "point estimate is reported anywhere in the paper, so it cannot",
        "be encoded even as a null."
      )
    ),
    INJSITE_THIGH = list(
      description = "Subcutaneous injection site (thigh vs abdomen)",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Tested as a candidate covariate on the absorption parameters and",
        "NOT retained: Chan 2025 Results states 'The site of SC",
        "administration and the SC formulation were not identified as",
        "statistically significant covariates on absorption parameters',",
        "and the Discussion notes this is consistent with most",
        "therapeutic proteins while cautioning that the abdomen subset is",
        "small (N = 39). This is a NEGATIVE result of real interest,",
        "because the earlier phase-Ib-only model (Chanu 2023, reference",
        "14 of the paper) had estimated site-specific bioavailabilities",
        "of 82.9 percent for thigh and 71.1 percent for abdomen. The",
        "single patient dosed in the upper arm was excluded from the",
        "analysis ('IGNORE(PTNM.EQ.20168) ; wrong site of injection').",
        "Documented for provenance only; no point estimate survives in",
        "the final model, so it cannot be encoded even as a null."
      )
    )
  )

  compartmentData <- list(
    depot       = list(analyte = "atezolizumab", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "atezolizumab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "atezolizumab", units = "mg", specimen = "serum", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 435,
    n_studies      = 1,
    n_observations = 3100,
    disease_state  = "locally advanced or metastatic non-small cell lung cancer",
    weight_range   = "30.0-117 kg (phase III subcutaneous arm; cohort medians 65.4-73.2 kg across all five cohorts)",
    sex_female_pct = 29.3,
    dose_range     = paste(
      "phase Ib: a single 1800 mg subcutaneous dose in the thigh",
      "(cohort 1, N = 13), 1200 mg subcutaneous every 2 weeks for three",
      "cycles in the thigh (cohort 2, N = 15), or 1800 mg subcutaneous",
      "every 3 weeks for three cycles with cycle 1 in the abdomen and",
      "cycles 2-3 in the thigh (cohort 3, N = 39), each co-mixed with",
      "recombinant human hyaluronidase PH20 and each followed by",
      "intravenous 1200 mg every 3 weeks thereafter; phase III:",
      "randomized 2:1 to subcutaneous 1875 mg every 3 weeks in the thigh",
      "as a ready-to-use rHuPH20 co-formulation (cohort 5, N = 246) or",
      "intravenous 1200 mg every 3 weeks (cohort 4, N = 122)"
    ),
    notes          = paste(
      "IMscin001 (NCT03735121), a two-part open-label study; part 1 is",
      "the phase Ib dose-finding portion and part 2 the phase III",
      "dose-confirmation portion. 435 of the 436 patients with PK data",
      "(99.8 percent) were analyzable. Baseline medians in the phase III",
      "subcutaneous arm (Chan 2025 Table 1): body weight 67.8 kg, tumor",
      "burden 79.5 mm, albumin 40.0 g/L, hemoglobin 123 g/L; 72.0 percent",
      "treatment-emergent-ADA-negative, 17.5 percent positive, 10.6",
      "percent missing. Each phase Ib patient received BOTH subcutaneous",
      "and intravenous atezolizumab and so served as their own control,",
      "which the Discussion identifies as what makes the bioavailability",
      "estimable; the phase III portion alone has no intensive",
      "absorption-phase sampling."
    )
  )

  ini({
    # ==================================================================
    # DISPOSITION -- every value below is FIXED.
    #
    # Chan 2025 Methods: "the typical values for the systemic parameters
    # that are associated with the two-compartment disposition model and
    # its covariates effects were fixed to those from the historical IV
    # popPK model, as theoretically the disposition parameters are
    # intrinsic to a molecule and should not change with different routes
    # of administration." Every one of these rows carries "FIX" in the
    # RSE column of Chan 2025 Table 2 and "FIX" in the control stream's
    # $THETA block, so all are encoded with fixed().
    # ==================================================================
    lcl <- fixed(log(0.200)); label("Clearance at reference covariates (L/day)")                                        # Chan 2025 Table 2 "CL (L/d) 0.200 FIX"; control stream $THETA "0.2 FIX ; CL"
    lvc <- fixed(log(3.28)); label("Central volume of distribution at reference covariates (L)")                        # Chan 2025 Table 2 "Vc (L) 3.28 FIX"; control stream $THETA "3.28 FIX ; V2"
    lvp <- fixed(log(3.63)); label("Peripheral volume of distribution at reference covariates (L)")                     # Chan 2025 Table 2 "Vp (L) 3.63 FIX"; control stream $THETA "3.63 FIX ; V3"
    lq <- fixed(log(0.546)); label("Intercompartmental clearance (L/day)")                                              # Chan 2025 Table 2 "Q (L/d) 0.546 FIX"; control stream $THETA "0.546 FIX ; Q"

    # ---- Disposition covariate effects (all FIXED) -------------------
    # The two categorical effects are FRACTIONAL (1 + theta), not
    # exponential. Chan 2025 Table 2 calls them "additive effect", which
    # is ambiguous; the control stream is explicit -- see the ADA_POS and
    # SEXF covariateData notes.
    e_alb_cl <- fixed(-1.12); label("Power exponent of baseline albumin on CL, normalized to 40 g/L (unitless)")        # Chan 2025 Table 2 "Albumin on CL (ALB/40 in g/L) -1.12 FIX"; control stream CLALBU=((ALB/40.00)**THETA(5)), $THETA "-1.12 FIX"
    e_ada_cl <- fixed(0.159); label("Fractional increase in CL for ADA-positive patients (unitless)")                   # Chan 2025 Table 2 "ADA status on CL (additive effect for positive ADA) 0.159 FIX"; control stream IF(ATAG.EQ.1) CLATAG=( 1 + THETA(6)), $THETA "0.159 FIX"
    e_tumsz_cl <- fixed(0.125); label("Power exponent of baseline tumor burden on CL, normalized to 63 mm (unitless)")  # Chan 2025 Table 2 "Tumor burden on CL (Tumor burden/63mm) 0.125 FIX"; control stream CLBSLD=((BSLD/63.00)**THETA(7)), $THETA "0.125 FIX"
    e_wt_cl <- fixed(0.808); label("Power exponent of baseline body weight on CL, normalized to 77 kg (unitless)")      # Chan 2025 Table 2 "Bodyweight on CL (BWT/77 in kg) 0.808 FIX"; control stream CLBWT=((BWT/77)**THETA(8)), $THETA "0.808 FIX"
    e_alb_vc <- fixed(-0.350); label("Power exponent of baseline albumin on Vc, normalized to 40 g/L (unitless)")       # Chan 2025 Table 2 "Albumin on Vc (ALB/40 in g/L) -0.350 FIX"; control stream V2ALBU=((ALB/40.00)**THETA(9)), $THETA "-0.35 FIX"
    e_wt_vc <- fixed(0.559); label("Power exponent of baseline body weight on Vc, normalized to 77 kg (unitless)")      # Chan 2025 Table 2 "Bodyweight on Vc (BWT/77 in kg) 0.559 FIX"; control stream V2BWT=((BWT/77)**THETA(10)), $THETA "0.559 FIX"
    e_sexf_vc <- fixed(-0.129); label("Fractional change in Vc for female sex (unitless)")                              # Chan 2025 Table 2 "Sex on Vc (additive effect for female) -0.129 FIX"; control stream IF(SEX.EQ.2) V2SEX=( 1 + THETA(11)), $THETA "-0.129 FIX"
    e_sexf_vp <- fixed(-0.272); label("Fractional change in Vp for female sex (unitless)")                              # Chan 2025 Table 2 "Sex on Vp (additive effect for female) -0.272 FIX"; control stream IF(SEX.EQ.2) V3SEX=( 1 + THETA(12)), $THETA "-0.272 FIX"

    # ==================================================================
    # ABSORPTION -- the layer Chan 2025 actually estimated.
    #
    # These four are the only structural parameters with a real RSE in
    # Chan 2025 Table 2, and the only $THETA entries in the control
    # stream without a FIX flag. NOTE that the control-stream $THETA
    # values (0.3 for KA, 0.7 for F1, 0.1 for both covariate exponents)
    # are INITIAL estimates; the final estimates below come from Table 2.
    # ==================================================================
    lka <- log(0.304); label("First-order subcutaneous absorption rate constant at reference albumin (1/day)")          # Chan 2025 Table 2 "KA (1/d) 0.304, RSE 3.0"; control stream $THETA (0,0.3) ; KA thigh/abdomen (initial estimate)

    # Bioavailability is carried on the LOGIT scale because the source
    # parameterises it that way:
    #   LF1 = LOG(THETA(14)/(1-THETA(14)))
    #   F1  = EXP(LF1*F1HGB + ETA(6)) / (1 + EXP(LF1*F1HGB + ETA(6)))
    # so the estimated quantity THETA(14) = 0.718 is F1 on the natural
    # 0-1 scale and logit() of it is what the covariate multiplies. The
    # eta ADDS to the covariate-scaled logit; it is not itself scaled.
    logitfdepot <- log(0.718 / (1 - 0.718)); label("Logit of subcutaneous bioavailability at reference hemoglobin (unitless logit; expit gives F1 = 0.718)")  # Chan 2025 Table 2 "F1 0.718, RSE 1.8"; control stream $THETA (0,0.7,1) ; F1 thigh/abdomen (initial estimate), transformed by LF1 = LOG(THETA(14)/(1-THETA(14)))

    e_alb_ka <- 0.795; label("Power exponent of baseline albumin on KA, normalized to 40 g/L (unitless)")               # Chan 2025 Table 2 "Albumin on KA (ALB/40 in g/L) 0.795, RSE 33.1"; control stream KAALB=((ALB/40.00)**THETA(16))
    e_hgb_logitfdepot <- 1.76; label("Power exponent of baseline hemoglobin on the MULTIPLIER of logit(F1), normalized to 123 g/L (unitless)")  # Chan 2025 Table 2 "Hemoglobin on F1 On logit scale (HGB/123 in g/L) 1.76, RSE 37.8"; control stream F1HGB=((HGB/123)**THETA(19))

    # ==================================================================
    # BETWEEN-SUBJECT VARIABILITY
    #
    # Chan 2025 Results: "Variability terms of the systemic parameters
    # were initially re-estimated, but this approach resulted in high
    # uncertainty on estimates of KA and F1. Therefore, the variability
    # terms of the systemic parameters were also fixed to the values
    # estimated in the historical IV model." Hence the CL/Vc/Vp block is
    # fixed() and only the two absorption etas are estimated.
    #
    # The block is taken verbatim from the control stream's
    # "$OMEGA BLOCK(3) FIX", which holds VARIANCES and COVARIANCES
    # directly (not correlations):
    #   var(CL)          = 0.0867   -> SD 29.4 percent
    #   cov(CL, Vc)      = 0.0182
    #   var(Vc)          = 0.0328   -> SD 18.1 percent
    #   cov(CL, Vp)      = -0.0234
    #   cov(Vc, Vp)      = 0.0265
    #   var(Vp)          = 0.114    -> SD 33.8 percent
    # Those three square roots reproduce the 29.4 / 18.1 / 33.8 percent
    # printed in the Chan 2025 Table 2 "Omega (standard deviation
    # scale)" rows, which confirms that Table 2 reports SDs while the
    # control stream reports variances. The CL-Vp covariance is
    # NEGATIVE; that sign is load-bearing and is not a transcription
    # slip.
    #
    # There is no eta on Q: the control stream carries "$OMEGA 0 FIX ;
    # POPIIV_Q" and Chan 2025 Table 2 lists no POPIIV Q row.
    # ==================================================================
    etalcl + etalvc + etalvp ~ fixed(c(
      0.0867,
      0.0182, 0.0328,
      -0.0234, 0.0265, 0.114
    ))

    # Absorption etas, ESTIMATED by Chan 2025. Table 2 prints them on the
    # standard-deviation scale as percentages, so the variances below are
    # the squares: 0.346^2 = 0.119716 and 0.830^2 = 0.6889.
    etalka ~ 0.119716                                                                                                  # Chan 2025 Table 2 'POPIIV KA 34.6%, RSE 23.7, shrinkage 38.1'
    etalogitfdepot ~ 0.6889                                                                                            # Chan 2025 Table 2 'POPIIV F1 83.0%, RSE 18.7, shrinkage 37.4'; this eta lives on the logit scale (control stream adds ETA(6) inside the expit)

    # ==================================================================
    # RESIDUAL UNEXPLAINED VARIABILITY
    #
    # Control stream $ERROR: "Y=IPRED*(1+EPS(1))+EPS(2)", i.e. combined
    # proportional plus additive on the linear concentration scale.
    # Chan 2025 Table 2 reports the two terms already on the SD scale
    # ("Proportional error (%) 19.0" and "Additive error (ug/mL) 15.4"),
    # so no square root is applied here. Cross-check: the control-stream
    # $SIGMA INITIAL values 0.033 and 166 are variances whose square
    # roots, 18.2 percent and 12.9 ug/mL, sit right beside the final
    # estimates -- consistent with Table 2 being SDs.
    # ==================================================================
    propSd <- 0.190; label("Proportional residual error (fraction)")                                                    # Chan 2025 Table 2 "Proportional error (%) 19.0, RSE 9.6, shrinkage 13.2"
    addSd <- 15.4; label("Additive residual error (ug/mL)")                                                             # Chan 2025 Table 2 "Additive error (ug/mL) 15.4, RSE 15.1"
  })

  model({
    # ------------------------------------------------------------------
    # Individual disposition parameters. Covariate forms transcribed from
    # the control-stream $PK blocks; the normalization constants 40 g/L,
    # 63 mm, 77 kg and 123 g/L are the historical-model reference values
    # printed inside those expressions and repeated in Chan 2025 Table 2.
    # ------------------------------------------------------------------
    cl <- exp(lcl + etalcl) *
      (ALB / 40)^e_alb_cl *
      (1 + e_ada_cl * ADA_POS) *
      (TUMSZ / 63)^e_tumsz_cl *
      (WT / 77)^e_wt_cl
    vc <- exp(lvc + etalvc) *
      (ALB / 40)^e_alb_vc *
      (WT / 77)^e_wt_vc *
      (1 + e_sexf_vc * SEXF)
    vp <- exp(lvp + etalvp) *
      (1 + e_sexf_vp * SEXF)
    q <- exp(lq)

    # ------------------------------------------------------------------
    # Absorption. KA carries a plain power covariate; F1 is built on the
    # logit scale with hemoglobin MULTIPLYING the logit, per the control
    # stream. At HGB = 123 g/L the multiplier is exactly 1, so a typical
    # patient returns expit(logit(0.718)) = 0.718.
    # ------------------------------------------------------------------
    ka <- exp(lka + etalka) * (ALB / 40)^e_alb_ka
    fdepot <- expit(logitfdepot * (HGB / 123)^e_hgb_logitfdepot + etalogitfdepot)

    # ------------------------------------------------------------------
    # Micro-constants (control stream $SUBROUTINE ADVAN4 TRANS4, which
    # parameterises on CL, V2, Q, V3, KA).
    # ------------------------------------------------------------------
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ------------------------------------------------------------------
    # ODE system. Compartment 1 of ADVAN4 is the subcutaneous depot,
    # compartment 2 the central (serum) compartment and compartment 3 the
    # peripheral compartment. Intravenous doses in IMscin001 are
    # administered directly into the central compartment and are NOT
    # subject to F1; only depot doses are.
    # ------------------------------------------------------------------
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    f(depot) <- fdepot

    # ------------------------------------------------------------------
    # Serum concentration. The control stream sets S2 = V2, so the
    # predicted concentration is the central amount divided by Vc, in
    # mg/L = ug/mL for a dose in mg and a volume in L.
    # ------------------------------------------------------------------
    Cc <- central / vc
    Cc ~ prop(propSd) + add(addSd)
  })
}
