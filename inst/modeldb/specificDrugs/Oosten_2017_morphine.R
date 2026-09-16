Oosten_2017_morphine <- function() {
  description <- paste(
    "Joint parent-metabolite population PK model for morphine and its two",
    "glucuronide metabolites (M3G, M6G) in 49 adult cancer patients treated",
    "for nociceptive cancer pain (Oosten 2017). Morphine: one-compartment",
    "disposition with three parallel first-order absorption routes",
    "(subcutaneous, oral immediate release, oral extended release with a",
    "lag time); the subcutaneous route is assumed completely bioavailable",
    "and the oral routes carry an extensive first pass. The oral first pass",
    "is a four-way partition of the absorbed dose - unchanged morphine,",
    "M3G, M6G and a residual lost fraction - parameterised as odds against",
    "the unchanged route, which is what makes oral bioavailability (0.372)",
    "and the first-pass metabolite fractions (0.355 for M3G, 0.0631 for",
    "M6G) sum to at most one by construction. M3G and M6G are each",
    "one-compartment models fed both by this oral first pass and by",
    "fixed-fraction transformation of systemic morphine clearance",
    "(0.57 to M3G, 0.10 to M6G, both fixed from literature); the two",
    "metabolites share a single estimated clearance and volume and differ",
    "only in how much of them is formed. Allometric body weight with",
    "theory-based exponents is applied a priori to every disposition",
    "parameter of every entity, and metabolite clearance rises linearly",
    "with estimated glomerular filtration rate up to a plateau at 90",
    "mL/min/1.73 m^2. All amounts are nmol of the free base and all",
    "concentrations nmol/L, so the 1:1 molar conversion of morphine to",
    "each glucuronide needs no molecular-weight factor."
  )
  reference <- paste(
    "Oosten AW, Abrantes JA, Jonsson S, Matic M, van Schaik RHN,",
    "de Bruijn P, van der Rijt CCD, Mathijssen RHJ.",
    "A Prospective Population Pharmacokinetic Study on Morphine Metabolism",
    "in Cancer Patients.",
    "Clin Pharmacokinet. 2017;56(6):649-659",
    "(published online 5 November 2016).",
    "doi:10.1007/s40262-016-0471-7.",
    sep = " "
  )
  vignette <- "Oosten_2017_morphine"
  units <- list(
    time = "h",
    dosing = "nmol",
    concentration = "nmol/L"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Oosten 2017 Fig. 1 (model diagram)
  # and Methods 2.4 ('Concentration data and doses of morphine were
  # expressed as free base in molar units (nmol/L and nmol, respectively)').
  compartmentData <- list(
    depot = list(analyte = "morphine", units = "nmol", specimen = "administration site", verified = TRUE),
    depot2 = list(analyte = "morphine", units = "nmol", specimen = "administration site", verified = TRUE),
    depot3 = list(analyte = "morphine", units = "nmol", specimen = "administration site", verified = TRUE),
    central = list(analyte = "morphine", units = "nmol", specimen = "plasma", verified = TRUE),
    central_m3g = list(analyte = "M3G", units = "nmol", specimen = "plasma", verified = TRUE),
    central_m6g = list(analyte = "M6G", units = "nmol", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Total body weight at study entry.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Allometric size descriptor applied a priori with theory-based",
        "exponents to every disposition parameter of every entity:",
        "(WT / 70)^0.75 on morphine clearance and on the common metabolite",
        "clearance, (WT / 70)^1 on the morphine volume and on the common",
        "metabolite volume (Table 3 footnotes c and d; Results 3.2",
        "'Allometric body weight with theory-based exponents was included",
        "a priori on all disposition parameters of all entities').",
        "Reference 70 kg. Study median 83 kg (range 53-140, Table 2)."
      ),
      source_name = "weight (kg)"
    ),
    CRCL = list(
      description = paste(
        "Estimated glomerular filtration rate from the Modification of Diet",
        "in Renal Disease (MDRD) formula, truncated above 90 mL/min/1.73 m^2."
      ),
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Linear (not power) covariate on the common metabolite clearance:",
        "CLmet = 4.71 * (WT/70)^0.75 * (1 + 0.0128 * (eGFR - 81)), Table 3",
        "footnote d. Methods: 'values > 90 mL/min/1.73 m^2 were truncated',",
        "and Results 3.2 states clearance is constant above that value, so",
        "the model applies min(CRCL, 90) before the linear term. Reference",
        "81 mL/min/1.73 m^2 = the study median (Table 2). The linear form",
        "goes negative below eGFR = 81 - 1/0.0128 = 3.1 mL/min/1.73 m^2,",
        "far below the observed range (33 to >90, Table 2); do not",
        "extrapolate the model there. Canonical register name CRCL covers",
        "MDRD-derived eGFR (same mapping as Franken_2015_morphine.R).",
        "Creatinine was missing for 4 of 49 patients and was imputed by the",
        "authors from a linear regression of eGFR on age and gender."
      ),
      source_name = "eGFR (MDRD, mL/min/1.73 m^2)"
    )
  )

  # Screened by the authors but NOT retained in the final model, so these
  # carry no parameter and must not appear in model().
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age at study entry.",
      units = "years",
      type = "continuous",
      notes = paste(
        "Tested on the pharmacokinetic profiles; 'Age did not statistically",
        "significantly improve the model fit (p > 0.01)' (Results 3.2).",
        "Median 60 years, range 38-80 (Table 2)."
      )
    ),
    SEXF = list(
      description = "Female sex indicator.",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Tested on morphine clearance; 'The inclusion of an effect of gender",
        "on clearance of morphine resulted in a decrease of 17.0% for females",
        "(p > 0.01) but was not retained in the model' (Results 3.2). The",
        "-17% point estimate is reported only in prose and has no reported",
        "uncertainty, so it is documented here rather than encoded.",
        "45% female (Table 2)."
      )
    ),
    SNP_UGT2B7_RS7438135 = list(
      description = "UGT2B7 -900G>A genotype (rs7438135).",
      units = "(categorical)",
      type = "categorical",
      notes = paste(
        "Tested on total morphine clearance and on the morphine metabolic",
        "clearances to M3G and M6G; no effect identified (p > 0.01,",
        "Results 3.3). GG 29%, GA 53%, AA 20% (Table 2)."
      )
    ),
    SNP_SLC22A1_RS12208357 = list(
      description = "SLC22A1 (OCT1) haplotype, built from rs72552763, rs12208357, rs34130495 and rs34059508.",
      units = "(categorical)",
      type = "categorical",
      notes = paste(
        "Tested on total morphine clearance and on the morphine metabolic",
        "clearances; no effect identified (p > 0.01, Results 3.3). Two",
        "active alleles 53%, one active/one inactive 37%, two inactive 10%",
        "(Table 2). Registered here under the rs12208357 member of the",
        "four-variant haplotype the authors constructed."
      )
    ),
    SNP_ABCC3_RS4793665 = list(
      description = "ABCC3 -211C>T genotype (rs4793665).",
      units = "(categorical)",
      type = "categorical",
      notes = paste(
        "Tested on total morphine clearance, on the morphine metabolic",
        "clearances, and on the clearance of the metabolites themselves;",
        "no effect identified (p > 0.01, Results 3.3). CC 10%, CT 57%,",
        "TT 33% (Table 2)."
      )
    )
  )

  population <- list(
    species = "human",
    n_subjects = 49,
    n_studies = 1,
    age_range = "38-80 years",
    age_median = "60 years",
    weight_range = "53-140 kg",
    weight_median = "83 kg",
    sex_female_pct = 45,
    race_ethnicity = c(Caucasian = 90, LatinAmerican = 2, UnknownOther = 8),
    disease_state = paste(
      "Adults admitted to the Erasmus MC Cancer Institute with moderate to",
      "severe nociceptive cancer-related pain, treated with morphine.",
      "89% had distant metastases; median WHO performance status 2",
      "(range 0-3). Primary tumour sites: breast 22%, colorectal 14%,",
      "prostate 12%, sarcoma 8%, other 43%. Both opioid-naive patients and",
      "patients rotating to morphine from fentanyl or oxycodone were",
      "eligible. Median baseline creatinine 72 umol/L (range 25-190),",
      "median eGFR 81 mL/min/1.73 m^2 (range 33 to >90, truncated at 90);",
      "7 patients had eGFR 33-57. Median serum albumin 40 g/L."
    ),
    dose_range = paste(
      "Median 2 mg/h for continuous subcutaneous morphine (range 0.6-14",
      "mg/h) and for subcutaneous bolus (range 0.6-10 mg); 40 mg twice",
      "daily for oral extended release (range 10-150 mg); 10 mg for oral",
      "immediate release (range 5-60 mg). Doses titrated to clinical",
      "effect. Routes during sampling: subcutaneous only 57%, oral",
      "extended plus immediate release 24%, oral immediate release only",
      "12%, both oral and subcutaneous consecutively 6%."
    ),
    regions = "The Netherlands (Erasmus MC Cancer Institute, Rotterdam)",
    notes = paste(
      "Demographics from Oosten 2017 Table 2. 410 plasma samples from 49",
      "patients, collected February 2010 to March 2014 (Dutch Trial",
      "Register NTR4369). NONMEM 7.3 with FOCE-I and eta-epsilon",
      "interaction; PsN 4.2.0, Xpose 4.4.1. Parameter precision from",
      "sampling-importance-resampling (SIR, five iterations) rather than",
      "bootstrap. Concentrations below the lower limit of quantitation",
      "(7.6% of morphine, 0.7% of M3G and 0.9% of M6G observations) were",
      "discarded. Doses and concentrations were handled as free base in",
      "molar units; the salts administered were morphine",
      "hydrochloride-3-water (375.84 mg/mmol, parenteral) and morphine",
      "5-sulphate-water (758.83 mg/mmol per 2 mmol morphine, oral).",
      "Assay ranges: morphine 1.00-100 ng/mL, M3G 10.0-1000 ng/mL,",
      "M6G 2.00-200 ng/mL. Model evaluated by prediction-corrected VPCs",
      "of concentrations and of metabolite:parent concentration ratios."
    )
  )

  ini({
    # ================================================================
    # Morphine absorption. Three parallel first-order routes; the
    # subcutaneous and oral immediate-release constants were fixed to
    # literature values (Results 3.2: 'Parameters describing the
    # absorption phases for subcutaneous and oral IR morphine were
    # fixed to literature values [27, 28]'), the extended-release one
    # was estimated ('p < 0.001 when compared with a fixed value of
    # 0.8 h-1').
    # ================================================================
    ltlag_er <- fixed(log(0.25))
    label("Absorption lag time for oral extended release (h)") # Table 3: t_lag,ER = 0.25 h, no RSE/CI reported
    lka_sc <- fixed(log(3.96))
    label("Absorption rate constant, subcutaneous (1/h)") # Table 3: ka,SC = 3.96 /h, literature value (Upton et al., ref 27)
    lka_oral_ir <- fixed(log(6.00))
    label("Absorption rate constant, oral immediate release (1/h)") # Table 3: ka,IR = 6.00 /h, literature value (refs 27, 28)
    lka_oral_er <- log(0.221)
    label("Absorption rate constant, oral extended release (1/h)") # Table 3: ka,ER = 0.221 /h (SIR RSE 17.7%, 95% CI 0.155-0.306)

    # ================================================================
    # Oral first-pass partition. Table 3 footnote b defines
    #   F_oral    = 1  / (1 + theta1 + theta2 + theta3)
    #   F1p,M3G   = t2 / (1 + theta1 + theta2 + theta3)
    #   F1p,M6G   = t1 / (1 + theta1 + theta2 + theta3)
    # so each theta is the odds of that first-pass route against the
    # route delivering unchanged morphine to the systemic circulation,
    # and the four shares sum to exactly one by construction. This is
    # the same odds-style partition that Padavia_2024_paracetamol.R
    # carries as lclrat_<metab>, moved from the systemic elimination
    # arms to the oral first pass, so the names take the first-pass
    # root lfprat_ rather than lclrat_ (which the register defines as
    # a ratio against the parent's UNCHANGED ELIMINATION clearance).
    # The thetas must stay primary: the etas below sit on theta1 and
    # theta3 and propagate non-linearly to all four shares through the
    # shared denominator, so an eta on a derived fraction cannot
    # reproduce the published random-effects structure.
    # ================================================================
    lfprat_m3g <- log(0.953)
    label("Odds of first-pass M3G formation against unchanged oral absorption (unitless)") # Table 3: theta2 = 0.953 (SIR RSE 8.95%, 95% CI 0.796-1.14)
    lfprat_m6g <- log(0.170)
    label("Odds of first-pass M6G formation against unchanged oral absorption (unitless)") # Table 3: theta1 = 0.170 (SIR RSE 10.4%, 95% CI 0.136-0.206)
    lfprat_other <- log(0.565)
    label("Odds of first-pass loss by other routes against unchanged oral absorption (unitless)") # Table 3: theta3 = 0.565 (SIR RSE 30.3%, 95% CI 0.310-1.01)

    # ================================================================
    # Morphine disposition, one compartment. Values are for a 70-kg
    # patient (Table 3 footnote c).
    # ================================================================
    lcl <- log(91.9)
    label("Morphine clearance at 70 kg (L/h)") # Table 3: CL,70kg = 91.9 L/h (SIR RSE 3.91%, 95% CI 85.8-99.9). The Discussion quotes 92.9 L/h for the same quantity; 91.9 is the tabulated estimate and is the value that reproduces the paper's own AUC ratios (see vignette Errata)
    lvc <- log(278)
    label("Morphine volume of distribution at 70 kg (L)") # Table 3: V,70kg = 278 L (SIR RSE 12.3%, 95% CI 221-351)
    e_wt_cl <- fixed(0.75)
    label("Allometric exponent of body weight on morphine clearance (unitless)") # Table 3 footnote c: theory-based exponent 0.75 for clearances
    e_wt_vc <- fixed(1)
    label("Allometric exponent of body weight on morphine volume (unitless)") # Table 3 footnote c: theory-based exponent 1 for volumes

    # ================================================================
    # Systemic fractions of morphine clearance forming each metabolite.
    # Results 3.2: 'The fractions of total morphine clearance forming
    # M3G and M6G were fixed to 0.57 and 0.10, respectively [4-6].'
    # Table 3 prints 0.573 and 0.104 for the same two rows; 0.57 and
    # 0.10 are used here because they are the values that reproduce
    # ALL SIX of the paper's own published AUC ratios exactly, while
    # 0.573 / 0.104 reproduce none of the subcutaneous ones (the
    # subcutaneous ratios have no first-pass term, so they pin these
    # fractions directly: M6G/morphine = Fm,M6G * CL / CLmet =
    # 0.10 * 91.9 / 4.71 = 1.95 as published, whereas 0.104 gives
    # 2.03). See the vignette's Errata section.
    # ================================================================
    fm_m3g <- fixed(0.57)
    label("Fraction of systemic morphine clearance forming M3G (unitless)") # Results 3.2 (literature refs 4-6); Table 3 prints 0.573
    fm_m6g <- fixed(0.10)
    label("Fraction of systemic morphine clearance forming M6G (unitless)") # Results 3.2 (literature refs 4-6); Table 3 prints 0.104

    # ================================================================
    # Metabolite disposition. Results 3.2: 'The metabolite disposition
    # parameters were estimated to common values, and the estimation of
    # separate clearance and volume parameters for each metabolite was
    # not found to be statistically significant.' The metabolite-
    # suffixed names below therefore carry the SAME single estimate
    # twice, to keep the source trace explicit at the point of use
    # (the convention Franken_2015_morphine.R uses for its shared
    # covariate exponents). The metabolites differ only in how much of
    # them is formed, and in their independent random effects below.
    # Reference subject: 70 kg, eGFR 81 mL/min/1.73 m^2.
    # ================================================================
    lcl_m3g <- log(4.71)
    label("Common metabolite clearance, M3G arm, at 70 kg and eGFR 81 (L/h)") # Table 3: CLmet,70kg = 4.71 L/h (SIR RSE 5.24%, 95% CI 4.24-5.20); single estimate common to M3G and M6G
    lcl_m6g <- log(4.71)
    label("Common metabolite clearance, M6G arm, at 70 kg and eGFR 81 (L/h)") # Table 3: CLmet,70kg = 4.71 L/h; the same single estimate as lcl_m3g
    lvc_m3g <- log(25.8)
    label("Common metabolite volume of distribution, M3G arm, at 70 kg (L)") # Table 3: Vmet,70kg = 25.8 L (SIR RSE 6.12%, 95% CI 22.8-29.0); single estimate common to M3G and M6G
    lvc_m6g <- log(25.8)
    label("Common metabolite volume of distribution, M6G arm, at 70 kg (L)") # Table 3: Vmet,70kg = 25.8 L; the same single estimate as lvc_m3g
    e_wt_cl_m3g <- fixed(0.75)
    label("Allometric exponent of body weight on M3G clearance (unitless)") # Table 3 footnote d: CLmet = 4.71 * (weight/70)^0.75 * (1 + 0.0128 * (eGFR - 81))
    e_wt_cl_m6g <- fixed(0.75)
    label("Allometric exponent of body weight on M6G clearance (unitless)") # Table 3 footnote d; the same theory-based exponent
    e_wt_vc_m3g <- fixed(1)
    label("Allometric exponent of body weight on M3G volume (unitless)") # Table 3 footnote c: theory-based exponent 1 for volumes
    e_wt_vc_m6g <- fixed(1)
    label("Allometric exponent of body weight on M6G volume (unitless)") # Table 3 footnote c; the same theory-based exponent

    # eGFR on metabolite clearance. LINEAR, not power: a fractional
    # change per mL/min/1.73 m^2 relative to the value at the median
    # eGFR of 81 (Table 3 row 'eGFR on CLmet,70kg').
    e_crcl_cl_m3g <- 0.0128
    label("Fractional change in M3G clearance per mL/min/1.73 m^2 of eGFR above 81 (unitless)") # Table 3: 0.0128 (SIR RSE 12.9%, 95% CI 0.00924-0.0156); single estimate common to both metabolites. Abstract cross-check: 4.71 * 0.0128 * 10 = 0.603 L/h per 10 mL/min/1.73 m^2, as published
    e_crcl_cl_m6g <- 0.0128
    label("Fractional change in M6G clearance per mL/min/1.73 m^2 of eGFR above 81 (unitless)") # Table 3: 0.0128; the same single estimate as e_crcl_cl_m3g

    # ================================================================
    # Interindividual variability. Table 3 reports each as a per cent
    # coefficient of variation with footnote a: 'For interindividual
    # and residual variability, %RSE is reported on the standard
    # deviation scale'. The printed per cent is 100 * omega (the
    # log-normal SD), NOT 100 * sqrt(exp(omega^2) - 1): the printed
    # SIR confidence intervals are the point estimate scaled by
    # 1 +/- 1.96 * RSE in that scale. The largest omega settles it -
    # for theta3, 98.2% with RSE 25.1% gives 49.9-146.5, matching the
    # published 50.1-146, whereas the exp form would give 43.6-187.
    # So variance = (CV/100)^2.
    # ================================================================
    etalka ~ 0.710^2
    # Table 3 'ka,all': 71.0% CV (eta-shrinkage 25.3%). One eta shared by all three route-specific ka (hence 'all'). Reported with no SIR RSE and no CI; every other NA row in Table 3 is a fixed or a derived quantity, and the absorption parameters were fixed from literature, so this variance is most likely fixed too - see the vignette Errata, where it is encoded as estimated because the paper does not say so explicitly.
    etalcl ~ 0.222^2
    # Table 3 'CL,70kg': 22.2% CV (SIR RSE 12.9%, 95% CI 16.9-27.8, eta-shrinkage 17.9%)
    etalvc ~ 0.747^2
    # Table 3 'V,70kg': 74.7% CV (SIR RSE 9.71%, 95% CI 60.4-88.6, eta-shrinkage 21.4%)

    # M3G and M6G clearance random effects are strongly correlated
    # (Table 3 'Correlation CL M3G - CL M6G' = 0.912, SIR RSE 11.8%,
    # 95% CI 0.864-0.952). Footnote e defines the reported correlation
    # as cov / sqrt(var1 * var2), so cov = rho * omega1 * omega2.
    etalcl_m3g + etalcl_m6g ~ c(
      0.362^2,
      0.912 * 0.362 * 0.368,
      0.368^2
    )
    # Table 3: CL M3G 36.2% CV (SIR RSE 10.6%, 95% CI 29.5-44.4, shrinkage 6.28%); CL M6G 36.8% CV (SIR RSE 11.8%, 95% CI 29.6-46.2, shrinkage 7.00%)

    etalvc_m3g ~ 0.247^2
    # Table 3 'V M3G': 24.7% CV (SIR RSE 18.4%, 95% CI 17.0-34.0, eta-shrinkage 30.7%). No correlation with V M6G is reported
    etalvc_m6g ~ 0.243^2
    # Table 3 'V M6G': 24.3% CV (SIR RSE 20.5%, 95% CI 16.3-34.6, eta-shrinkage 39.0%)

    # Random effects on the first-pass odds. Table 3 reports IIV on
    # theta1 and theta3 only; theta2 carries none.
    etalfprat_m6g ~ 0.150^2
    # Table 3 'theta1': 15.0% CV (SIR RSE 18.7%, 95% CI 8.91-19.9, eta-shrinkage 62.7%)
    etalfprat_other ~ 0.982^2
    # Table 3 'theta3': 98.2% CV (SIR RSE 25.1%, 95% CI 50.1-146, eta-shrinkage 57.1%)

    # ================================================================
    # Residual variability. Proportional for all three entities
    # (Table 3 'Residual variability'), reported as per cent CV on the
    # standard deviation scale per footnote a. Oosten 2017 additionally
    # estimated correlations BETWEEN the three residuals of one sample
    # using the NONMEM L2 data item (morphine-M3G 0.420, morphine-M6G
    # 0.386, M3G-M6G 0.918); nlmixr2 has no analogue, so those three
    # off-diagonal sigma terms are not represented here. They affect
    # only the within-sample scatter of a simulated concentration
    # RATIO, never a concentration, an AUC or a typical-value profile.
    # ================================================================
    propSd <- 0.286
    label("Proportional residual SD for morphine (fraction)") # Table 3: 28.6% CV (SIR RSE 4.31%, 95% CI 26.5-31.2, eps-shrinkage 8.45%)
    propSd_m3g <- 0.200
    label("Proportional residual SD for M3G (fraction)") # Table 3: 20.0% CV (SIR RSE 4.14%, 95% CI 18.5-21.7, eps-shrinkage 8.00%)
    propSd_m6g <- 0.239
    label("Proportional residual SD for M6G (fraction)") # Table 3: 23.9% CV (SIR RSE 4.04%, 95% CI 22.2-26.0, eps-shrinkage 8.00%)
  })

  model({
    # ------------------------------------------------------------
    # Reference covariate values (Table 2 medians / Table 3 footnotes).
    # ------------------------------------------------------------
    wt_ref <- 70 # Table 3 footnote c: parameters reported for a 70-kg patient
    egfr_ref <- 81 # Table 3 footnote d: reference eGFR 81 mL/min/1.73 m^2 (Table 2 median)
    egfr_max <- 90 # Methods: 'values > 90 mL/min/1.73 m^2 were truncated'; Results 3.2: clearance constant above 90

    egfr_use <- min(CRCL, egfr_max)

    # ------------------------------------------------------------
    # Route-specific absorption. One shared random effect ('ka, all').
    # ------------------------------------------------------------
    ka_sc <- exp(lka_sc + etalka)
    ka_oral_ir <- exp(lka_oral_ir + etalka)
    ka_oral_er <- exp(lka_oral_er + etalka)
    tlag_er <- exp(ltlag_er)

    # ------------------------------------------------------------
    # Oral first-pass partition (Table 3 footnote b). The absorbed
    # oral dose splits four ways at a single point: unchanged morphine
    # into central, M3G and M6G straight into their own compartments,
    # and a residual fraction lost. The partition is applied to the
    # ABSORPTION FLUX out of each oral depot rather than to the dose
    # record, which is the reading Fig. 1 draws (every arrow, solid
    # and dashed, leaves the oral dose box) and the one the odds
    # parameterisation implies (odds of competing first-order routes
    # out of a shared compartment). The two readings differ only in
    # the early metabolite profile after an extended-release dose;
    # every AUC and every steady state is identical. See the vignette
    # Errata.
    # ------------------------------------------------------------
    fprat_m3g <- exp(lfprat_m3g)
    fprat_m6g <- exp(lfprat_m6g + etalfprat_m6g)
    fprat_other <- exp(lfprat_other + etalfprat_other)
    fp_denom <- 1 + fprat_m6g + fprat_m3g + fprat_other
    fdepot <- 1 / fp_denom
    f1p_m3g <- fprat_m3g / fp_denom
    f1p_m6g <- fprat_m6g / fp_denom

    # ------------------------------------------------------------
    # Morphine disposition with allometric body weight.
    # ------------------------------------------------------------
    cl <- exp(lcl + etalcl) * (WT / wt_ref)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / wt_ref)^e_wt_vc

    # ------------------------------------------------------------
    # Metabolite disposition: one common typical value per parameter,
    # separate random effects, allometric weight, and a LINEAR eGFR
    # effect on clearance (Table 3 footnote d).
    # ------------------------------------------------------------
    egfr_factor_m3g <- 1 + e_crcl_cl_m3g * (egfr_use - egfr_ref)
    egfr_factor_m6g <- 1 + e_crcl_cl_m6g * (egfr_use - egfr_ref)
    cl_m3g <- exp(lcl_m3g + etalcl_m3g) * (WT / wt_ref)^e_wt_cl_m3g * egfr_factor_m3g
    cl_m6g <- exp(lcl_m6g + etalcl_m6g) * (WT / wt_ref)^e_wt_cl_m6g * egfr_factor_m6g
    vc_m3g <- exp(lvc_m3g + etalvc_m3g) * (WT / wt_ref)^e_wt_vc_m3g
    vc_m6g <- exp(lvc_m6g + etalvc_m6g) * (WT / wt_ref)^e_wt_vc_m6g

    # ------------------------------------------------------------
    # Fluxes. Amounts are nmol of free base throughout, so the 1:1
    # molar conversion of morphine to each glucuronide needs no
    # molecular-weight factor (Methods 2.4).
    # ------------------------------------------------------------
    abs_oral_ir <- ka_oral_ir * depot2
    abs_oral_er <- ka_oral_er * depot3
    abs_oral <- abs_oral_ir + abs_oral_er
    elim_morphine <- cl * central / vc

    # ------------------------------------------------------------
    # ODE system (Fig. 1). depot = subcutaneous (assumed completely
    # bioavailable, Methods 2.4), depot2 = oral immediate release,
    # depot3 = oral extended release (lag time tlag_er).
    # ------------------------------------------------------------
    d/dt(depot) <- -ka_sc * depot
    d/dt(depot2) <- -abs_oral_ir
    d/dt(depot3) <- -abs_oral_er
    d/dt(central) <- ka_sc * depot + fdepot * abs_oral - elim_morphine
    d/dt(central_m3g) <-
      fm_m3g * elim_morphine + f1p_m3g * abs_oral - cl_m3g * central_m3g / vc_m3g
    d/dt(central_m6g) <-
      fm_m6g * elim_morphine + f1p_m6g * abs_oral - cl_m6g * central_m6g / vc_m6g

    alag(depot3) <- tlag_er

    # ------------------------------------------------------------
    # Observations, nmol/L (Methods 2.4).
    # ------------------------------------------------------------
    Cc <- central / vc
    Cc_m3g <- central_m3g / vc_m3g
    Cc_m6g <- central_m6g / vc_m6g

    Cc ~ prop(propSd)
    Cc_m3g ~ prop(propSd_m3g)
    Cc_m6g ~ prop(propSd_m6g)
  })
}
