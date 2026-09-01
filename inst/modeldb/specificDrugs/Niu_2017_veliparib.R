Niu_2017_veliparib <- function() {
  description <- "Joint parent plus metabolite population PK model for the oral PARP inhibitor veliparib (ABT-888) and its primary active metabolite M8, in patients with BRCA 1/2-mutated cancer or PARP-sensitive tumor types. Veliparib disposition is two-compartment with first-order absorption and an absorption lag time; total apparent clearance CL/F is split into a renal arm CLR/F = CL/F * frenal, which is a power function of Cockcroft-Gault creatinine clearance, and a non-renal arm CLNR/F = CL/F * fm, which is the formation clearance that is assumed to produce all of the M8. M8 is described by a two-compartment model with first-order formation. Lean body mass enters as a power function on the veliparib central volume of distribution. The renal and metabolised fractions were not identifiable from these data and were held at frenal = 0.7 and fm = 1 - frenal = 0.3 from a published mass-balance study."
  reference <- paste(
    "Niu J., Scheuerell C., Mehrotra S., Karan S., Puhalla S., Kiesel B. F.,",
    "Ji J., Chu E., Gopalakrishnan M., Ivaturi V., Gobburu J., Beumer J. H.",
    "(2017). Parent-metabolite pharmacokinetic modeling and pharmacodynamics",
    "of veliparib (ABT-888), a PARP inhibitor, in patients with BRCA 1/2-mutated",
    "cancer or PARP-sensitive tumor types.",
    "The Journal of Clinical Pharmacology 57(8):977-987.",
    "doi:10.1002/jcph.892.",
    sep = " "
  )
  vignette <- "Niu_2017_veliparib"

  # Niu 2017 Methods "Software and Estimation Methods": "The molecular weights
  # of veliparib and M8 are similar, with a ratio of M8 to veliparib of 1.057.
  # Hence, no corrections for the concentration units were performed, and
  # nanograms per milliliter was used for analysis purposes." The parent ->
  # metabolite flux in model() below is therefore carried 1:1 in MASS units
  # with no molecular-weight ratio, exactly as the authors fitted it.
  #
  # Doses are in mg and volumes in L, so an ODE state divided by its volume is
  # mg/L == ug/mL. Both observables are multiplied by 1000 ng/ug so that Cc and
  # Cc_m8 come out in the ng/mL the paper reports (and so that the additive
  # residual SDs below stay verbatim in ng/mL).
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Niu 2017 Figure 1 (schematic of the
  # final parent-metabolite model) and Methods "Sampling and Assay for
  # Veliparib and M8 Serum Concentration" (plasma sampling).
  compartmentData <- list(
    depot           = list(analyte = "veliparib", units = "mg", specimen = "administration site", verified = TRUE),
    central         = list(analyte = "veliparib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1     = list(analyte = "veliparib", units = "mg", specimen = "plasma", verified = TRUE),
    central_m8      = list(analyte = "M8", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_m8  = list(analyte = "M8", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description        = "Creatinine clearance estimated by the Cockcroft-Gault formula, NOT body-surface-area normalised",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Raw (non-BSA-normalised) Cockcroft-Gault creatinine clearance in mL/min, the",
        "only renal-function covariate in the model. Niu 2017 Table 1 footnote: 'CLCR",
        "was estimated using the Cockcroft-Gault formula. Estimated CLCR values higher",
        "than 120 were assigned a value of 120 mL/min as a reasonable upper limit for",
        "the analysis.' That 120 mL/min CAP IS PART OF THE FITTED MODEL, not a display",
        "convention -- Niu 2017 Results states it was applied to 20 of the 71 patients",
        "before estimation, so a downstream user must cap CRCL at 120 mL/min before",
        "driving this model or the renal clearance arm will be extrapolated beyond the",
        "range the exponent was estimated over. Enters as a median-centred power term",
        "on the renal clearance arm: cl_renal = cl * frenal * (CRCL/95)^0.903",
        "(Niu 2017 Eq. 5). The 95 mL/min centring value matches the cohort mean CLCR",
        "of 95.2 mL/min (Table 1, Total column) and the Abstract's typical subject",
        "('CLCR, 95 mL/min'). Observed: mean 95.2 mL/min, range 48.5-120 (Table 1).",
        "The units are plain mL/min rather than the register's default",
        "mL/min/1.73 m^2; the same raw Cockcroft-Gault usage is already registered for",
        "CRCL via the CLCR source alias (precedents: Delattre_2010_amikacin.R,",
        "Jonsson_2015_edoxaban.R)."
      ),
      source_name        = "CLCR"
    ),
    LBM = list(
      description        = "Lean body mass, the size descriptor on the veliparib central volume of distribution",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters as a median-centred power term on the veliparib central volume:",
        "vc = vc * (LBM/48)^1.21 (Niu 2017 Eq. 7). The 48 kg centring value matches",
        "the cohort mean LBM of 47.6 kg (Table 1, Total column) and the Abstract's",
        "typical subject ('LBM, 48 kg'). Observed: mean 47.6 kg, range 36.6-63.7 kg",
        "(Table 1). BODY-COMPOSITION FORMULA NOT RESTATED: Niu 2017 Methods cites",
        "reference 16 (Cheymol G. Effects of obesity on pharmacokinetics: implications",
        "for drug therapy. Clin Pharmacokinet. 2000;39(3):215-231) for LBM but does not",
        "print the equation it used, so a downstream user reproducing LBM from height",
        "and weight must choose a formula (James, Boer or Hume) and should record which",
        "-- the choice materially changes LBM at a given weight and height. Niu 2017",
        "Discussion motivates LBM over total body weight: veliparib is hydrophilic and",
        "nearly 60% of this cohort had a BMI of 25 or larger. Total body weight was",
        "screened and REJECTED in favour of LBM (Supplemental Table S1 run 105, an",
        "allometric weight model with exponents fixed at 1 and 0.75, was rejected)."
      ),
      source_name        = "LBM"
    )
  )

  # Covariates that Niu 2017 SCREENED in the forward-addition / backward-
  # elimination step but did NOT retain in the final model (Methods "Covariate
  # Model Development": "The covariates that were screened included age, total
  # body weight (WT), lean body mass (LBM), body surface area (BSA), creatinine
  # clearance (CLCR), sex, and liver function (ALT, AST, and total bilirubin)").
  # Only CLCR and LBM survived. Niu 2017 Discussion: "After adding CLCR on
  # CLR/F, none of the other patient demographic measures, such as body weight,
  # LBM, or age, were identified as significant covariates because of possible
  # collinearity with CLCR." Recorded here so the paper's covariate screen is
  # preserved without declaring covariates that model() never references.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age", units = "years", type = "continuous",
      notes = "Screened, not retained. Cohort mean 53.5 years, range 28.0-84.0 (Niu 2017 Table 1)."
    ),
    WT = list(
      description = "Total body weight", units = "kg", type = "continuous",
      notes = paste(
        "Screened, not retained; LBM was preferred as the size descriptor on Vc/F.",
        "Supplemental Table S1 run 105 (allometric WT with exponents fixed at 1 for",
        "volumes and 0.75 for clearances) was Rejected: -2LL moved from -3389.09 to",
        "-3388.90, i.e. no improvement. Cohort mean 74.9 kg, range 45.0-119 (Table 1)."
      )
    ),
    BSA = list(
      description = "Body surface area, estimated by the Mosteller formula", units = "m^2",
      type = "continuous",
      notes = "Screened, not retained. Niu 2017 Methods cites the Mosteller formula (reference 18); BSA values are not tabulated in the paper."
    ),
    SEXF = list(
      description = "Biological sex indicator, 1 = female, 0 = male", units = "(binary)",
      type = "binary", reference_category = "0 (male)",
      notes = paste(
        "Screened, not retained. Niu 2017 coded sex as 0 for females and 1 for males",
        "(Eq. 4), i.e. the SEXM orientation; the canonical SEXF is the inverse",
        "(SEXF = 1 - SEXM). The cohort was 69/71 female (97%) with only 2 male",
        "patients (Table 1), so a sex effect was effectively not estimable here."
      )
    ),
    ALT = list(
      description = "Alanine aminotransferase", units = "U/L", type = "continuous",
      notes = "Screened as part of liver function, not retained. Cohort mean 28.3 U/L, range 7.00-135 (Niu 2017 Table 1)."
    ),
    AST = list(
      description = "Aspartate aminotransferase", units = "U/L", type = "continuous",
      notes = "Screened as part of liver function, not retained. Cohort mean 33.7 U/L, range 11.0-117 (Niu 2017 Table 1)."
    ),
    TBILI = list(
      description = "Total bilirubin", units = "mg/dL", type = "continuous",
      notes = "Screened as part of liver function, not retained. Cohort mean 0.51 mg/dL, range 0.10-1.10 (Niu 2017 Table 1); reported in mg/dL, not the register's canonical umol/L."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 67L,
    n_studies      = 1L,
    age_range      = "28.0-84.0 years",
    age_median     = "53.5 years (cohort mean; the paper reports mean and range, not median)",
    weight_range   = "45.0-119 kg",
    weight_median  = "74.9 kg (cohort mean; the paper reports mean and range, not median)",
    sex_female_pct = 97.2,
    race_ethnicity = "Not reported",
    disease_state  = paste(
      "Adults with BRCA 1/2-mutated cancer or PARP-sensitive tumor types, enrolled in a",
      "phase 1 multicenter, randomized, open-label dose-escalation study of chronically",
      "dosed single-agent veliparib (ClinicalTrials.gov NCT00892736). ECOG performance",
      "status 0 in 43 (61%), 1 in 22 (31%) and 2 in 6 (8%) of the 71 enrolled patients."
    ),
    dose_range     = paste(
      "Oral veliparib twice daily without regard to meals at 9 AM/PM dose levels:",
      "50/50, 100/50, 100/100, 150/100, 150/150, 200/200, 300/300, 400/400 and",
      "500/500 mg. Only the morning dose was given on day 1; twice-daily dosing started",
      "on day 2 and continued for at least one 28-day cycle."
    ),
    regions        = "Not reported (multicenter, United States)",
    renal_function = paste(
      "Cockcroft-Gault CLCR mean 95.2 mL/min, range 48.5-120 (Niu 2017 Table 1). For 20",
      "of the 71 patients the Cockcroft-Gault estimate exceeded 120 mL/min and was",
      "assigned 120 mL/min before estimation, so the fitted CLCR range is effectively",
      "48.5-120 mL/min and the model should not be extrapolated above that cap. No",
      "patients with clinically significant renal impairment were enrolled; a dedicated",
      "renal-dysfunction study (NCT01366144) was ongoing at publication."
    ),
    notes          = paste(
      "Baseline demographics from Niu 2017 Table 1 (Total column, n = 71 enrolled).",
      "THE PK MODEL WAS FIT TO 67 PATIENTS, not 71: Table 1 footnote states 'Of the 71",
      "patients, 67 had pharmacokinetic samples for veliparib, and 38 of the 67 had M8",
      "pharmacokinetic samples collected. Of the 71 patients, 41 had PAR levels in PBMC",
      "collected. Four patients did not have pharmacokinetic results but had PAR levels",
      "in PBMC.' n_subjects above is the 67 who contributed to the PK fit; the",
      "demographic summary statistics quoted in this metadata are the n = 71 Table 1",
      "values, which are the only ones the paper reports. Observations: 1214 veliparib",
      "and 656 M8 plasma concentrations, plus 295 PAR concentrations from 41 patients",
      "for the exploratory PD analysis. Other baseline values (Table 1, Total): lean",
      "body mass mean 47.6 kg (36.6-63.7), height mean 164 cm (152-179), BMI mean 27.8",
      "kg/m^2 (17.6-44.1), serum creatinine mean 0.78 mg/dL (0.4-1.37), albumin mean",
      "3.77 (2.00-4.50) reported under a 'g/L' header that is almost certainly g/dL.",
      "PK sampling on days 1 and 15 of cycle 1 predose and 0.5, 1, 1.5, 2, 3, 4, 6, 8",
      "and 24 h after the morning dose, plus a single predose day-4 draw at 200/200 mg",
      "and above. Estimation was FOCE-ELS in Phoenix NLME 1.4; precision from a",
      "250-sample nonparametric bootstrap. Final model OFV -3608."
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # All values are the FINAL population estimates from Niu 2017 Table 2,
    # column "Estimate". The %RSE column of that table was computed from the
    # 250-sample bootstrap, and the bootstrap medians agree with every
    # structural estimate to three digits, so Table 2 is internally
    # consistent and is used verbatim.
    #
    # KNOWN ABSTRACT / TABLE DISCREPANCY: the Abstract reports "central and
    # peripheral volume of distribution for veliparib were estimated as
    # ... 98.7 L, and 48.3 L", whereas Table 2 reports Vc/F = 99.2 L and
    # Vp/F = 47.8 L. Both pairs sum to the 147 L that the Discussion quotes
    # for Vd/F, so the split differs but the total does not. Table 2 is used
    # here because it is the final-parameter table, it carries %RSE and
    # bootstrap CIs, and its bootstrap medians (99.1 and 47.5 L) corroborate
    # it. See the vignette's Assumptions and deviations section.
    # ---------------------------------------------------------------------

    # -- Veliparib absorption. First-order with a lag time; no bioavailability
    #    term was estimated, so every volume and clearance below is APPARENT
    #    (X/F). Niu 2017 Results: "A 2-compartment model with first-order
    #    absorption and lag time adequately described the concentration-time
    #    profile of veliparib."
    lka   <- log(2.02)  ; label("Veliparib first-order absorption rate constant Ka (1/h, log scale)")    # Niu 2017 Table 2 'Ka (h-1) = 2.02' (%RSE 16.5; bootstrap median 2.02, 95%CI 1.59-2.93)
    ltlag <- log(0.272) ; label("Veliparib absorption lag time tlag (h, log scale)")                     # Niu 2017 Table 2 'tlag (h) = 0.272' (%RSE 9.53; bootstrap median 0.270, 95%CI 0.241-0.318)

    # -- Veliparib disposition, two-compartment.
    lvc   <- log(99.2)  ; label("Veliparib apparent central volume Vc/F at LBM = 48 kg (L, log scale)")  # Niu 2017 Table 2 'Vc/F (L) = 99.2' (%RSE 6.91; bootstrap median 99.1, 95%CI 89.5-119.6)
    lvp   <- log(47.8)  ; label("Veliparib apparent peripheral volume Vp/F (L, log scale)")              # Niu 2017 Table 2 'Vp/F (L) = 47.8' (%RSE 10.2; bootstrap median 47.5, 95%CI 31.2-56.7)
    lq    <- log(17.9)  ; label("Veliparib apparent intercompartmental clearance Q/F (L/h, log scale)")  # Niu 2017 Table 2 'Q/F (L h-1) = 17.9' (%RSE 18.1; bootstrap median 17.9, 95%CI 11.4-23.6)
    lcl   <- log(17.3)  ; label("Veliparib total apparent clearance CL/F at CLCR = 95 mL/min (L/h, log scale)")  # Niu 2017 Table 2 'CL/F (L h-1) = 17.3' (%RSE 3.85; bootstrap median 17.3, 95%CI 15.9-19.1)

    # -- Split of total apparent clearance into the renal and the M8-forming
    #    arm. Niu 2017 Methods "Base Model": "Using current data, the fraction
    #    of veliparib metabolized to M8 (fm) was not identifiable; however,
    #    from the literature, it is known that an average of 70% oral
    #    veliparib is renally excreted (frenal = 70%). Thus, veliparib not
    #    cleared by renal excretion was assumed to be metabolized into M8
    #    (fm = 1 - frenal = 30%)." Only fm is carried as a parameter here;
    #    frenal is recovered inside model() as 1 - fm, which is the paper's
    #    own identity (Figure 1 legend: "fm = 1 - frenal = 0.3") and so keeps
    #    the two from ever drifting apart.
    #
    #    CONSEQUENCE FOR M8: because the model routes ALL non-renally-cleared
    #    veliparib into M8, the M8 parameters below are apparent values in the
    #    same sense as the parent's. Niu 2017 Discussion notes the fraction
    #    actually reaching M8 is probably nearer 22% (30% non-renal x 72%
    #    CYP2D6 share of veliparib metabolism), which would scale Vc_met,
    #    Vp_met, Qmet and CLmet by 0.3/0.22 if it were used instead.
    fm    <- fixed(0.3) ; label("Fraction of veliparib apparent clearance forming M8 (unitless)")        # Niu 2017 Table 2 'fm = 0.3 (fix)' and Figure 1 legend '(apparent) fraction of veliparib converted to M8 (fm = 1-frenal = 0.3)'

    # -- M8 disposition, two-compartment with first-order formation. Niu 2017
    #    Results: "The M8 PK was adequately described with a 2-compartment
    #    model."
    lvc_m8 <- log(23.6) ; label("M8 apparent central volume Vc_met (L, log scale)")                      # Niu 2017 Table 2 'Vc_met (L) = 23.6' (%RSE 16.6; bootstrap median 23.5, 95%CI 16.3-33.8)
    lvp_m8 <- log(51.4) ; label("M8 apparent peripheral volume Vp_met (L, log scale)")                   # Niu 2017 Table 2 'Vp_met (L) = 51.4' (%RSE 15.4; bootstrap median 51.3, 95%CI 33.0-70.2)
    lq_m8  <- log(29.0) ; label("M8 apparent intercompartmental clearance Qmet (L/h, log scale)")        # Niu 2017 Table 2 'Qmet (L h-1) = 29.0' (%RSE 18.2; bootstrap median 29.0, 95%CI 21.0-43.6)
    lcl_m8 <- log(22.8) ; label("M8 apparent clearance CLmet (L/h, log scale)")                          # Niu 2017 Table 2 'CLmet (L h-1) = 22.8' (%RSE 7.74; bootstrap median 22.7, 95%CI 19.2-28.8)

    # -- Covariate effects. Both are median-centred power terms of the general
    #    form of Niu 2017 Eq. 3, theta_i = theta_TV * (cov_i/cov_med)^theta_cov.
    e_lbm_vc        <- 1.21  ; label("Power exponent on (LBM/48) for veliparib Vc/F (unitless)")         # Niu 2017 Table 2 'LBM effect on Vc/F, power_LBM = 1.21' (%RSE 20.8; bootstrap median 1.21, 95%CI 0.72-1.56); Eq. 7
    e_crcl_cl_renal <- 0.903 ; label("Power exponent on (CRCL/95) for veliparib CLR/F (unitless)")       # Niu 2017 Table 2 'CLCR effect on CLrenal/F, power_CLcr = 0.903' (%RSE 20.1; bootstrap median 0.901, 95%CI 0.525-1.33); Eq. 5

    # ---------------------------------------------------------------------
    # Inter-individual variability. Niu 2017 Results: "Between-subject
    # variability was estimated on Ka, tlag, Vc/F, Vp/F, CLR/F, CLNR/F,
    # Vc_met, CLmet, and Vp_met, whereas that on Q and Qmet was not estimated
    # in the final model." Supplemental Table S1 confirms the two removals
    # were the last two model-building steps (run 112 "remove ETA on Q",
    # Accepted; run 113 "remove ETA on Qmet", Final Model), so lq and lq_m8
    # deliberately carry no eta below.
    #
    # NO IIV CORRELATION BLOCK. Niu 2017 Results: "Correlation between random
    # effects of CLNR/F, CLR/F, and Vc/F decreased the OFV by 64 units but
    # worsened the parameter precision and increased the condition number to
    # 10^8. Hence, it was dropped from further analysis." The final model is
    # diagonal.
    #
    # BACK-TRANSFORM. Table 2 labels these rows "omega for Ka (CV%)",
    # "omega Vc/F (CV%)", etc. -- the symbol named is omega, and Niu 2017
    # Eq. 1 defines omega as the square root of the variance of eta
    # ("eta_i ... is assumed to follow normal distribution with a mean of
    # zero and a variance of omega^2"). The tabulated percentage is therefore
    # omega itself on the APPROXIMATE standard-deviation scale, i.e.
    # %CV = 100 * sqrt(omega^2), NOT the exact log-normal
    # 100 * sqrt(exp(omega^2) - 1). Every variance below is (%CV / 100)^2.
    # Corroboration: the %RSE column, which Table 2 says "was calculated from
    # bootstrap results", is roughly twice the RSE implied by the tabulated
    # bootstrap CI and point estimate (e.g. Ka: CI 64.0-114.4 about a
    # 78.1 estimate gives RSE ~16.5%, and Table 2 prints 33.7%), which is the
    # signature of an RSE reported on the VARIANCE scale next to an estimate
    # reported on the SD scale -- so the estimate column is not itself a
    # variance. Same convention and same reasoning as the sibling
    # parent-metabolite extraction Jonsson_2015_edoxaban.R. See the
    # vignette's Assumptions and deviations section.
    # ---------------------------------------------------------------------
    etalka        ~ 0.609961  # var = 0.781^2 ; Niu 2017 Table 2 'omega for Ka (CV%) = 78.1%' (%RSE 33.7; bootstrap median 78.1%, 95%CI 64.0-114.4%)
    etaltlag      ~ 0.208849  # var = 0.457^2 ; Niu 2017 Table 2 'omega for tlag (CV%) = 45.7%' (%RSE 22.3; bootstrap median 45.6%, 95%CI 24.6-45.9%)
    etalvc        ~ 0.071289  # var = 0.267^2 ; Niu 2017 Table 2 'omega Vc/F (CV%) = 26.7%' (%RSE 29.4; bootstrap median 26.7%, 95%CI 18.7-35.3%)
    etalvp        ~ 0.146689  # var = 0.383^2 ; Niu 2017 Table 2 'omega Vp/F (CV%) = 38.3%' (%RSE 47.6; bootstrap median 38.3%, 95%CI 29.0-47.3%)
    etalcl_renal  ~ 0.120409  # var = 0.347^2 ; Niu 2017 Table 2 'omega CLR/F (CV%) = 34.7%' (%RSE 29.5; bootstrap median 34.7%, 95%CI 27.2-45.2%)
    etalcl_nonren ~ 0.221841  # var = 0.471^2 ; Niu 2017 Table 2 'omega CLNR/F (CV%) = 47.1%' (%RSE 27.9; bootstrap median 47.1%, 95%CI 34.6-57.4%)
    etalvc_m8     ~ 0.383161  # var = 0.619^2 ; Niu 2017 Table 2 'omega Vc_met (CV%) = 61.9%' (%RSE 35.5; bootstrap median 62.0%, 95%CI 47.2-76.7%)
    etalvp_m8     ~ 0.439569  # var = 0.663^2 ; Niu 2017 Table 2 'omega Vp_met (CV%) = 66.3%' (%RSE 31.7; bootstrap median 66.3%, 95%CI 42.4-85.4%)
    etalcl_m8     ~ 0.107584  # var = 0.328^2 ; Niu 2017 Table 2 'omega CLmet (CV%) = 32.8%' (%RSE 22.2; bootstrap median 32.8%, 95%CI 21.6-39.0%)

    # ---------------------------------------------------------------------
    # Residual variability: proportional plus additive on both analytes.
    # Niu 2017 Results: "The proportional plus additive residual error model
    # best accounted for the unexplained variability of the observed
    # concentrations for veliparib and M8."
    #
    # M8 PROPORTIONAL ERROR IS A PRODUCT, not a directly tabulated number.
    # Table 2 leaves the "Proportional error for M8 (%CV), MultStdev" row
    # blank and gives the recipe in its own notes: "The proportional residual
    # error for M8 is expanded as: CMultStdev * Ratio * MultStdev, where
    # MultStdev is the coefficient variation for proportional error for M8
    # and the value (%CV) was fixed as 100%." With CMultStdev = 25.1%,
    # Ratio = 0.808 and MultStdev = 1.00, the M8 proportional SD is
    # 0.251 * 0.808 * 1.00 = 0.202808. This is Phoenix NLME's shared
    # multiplicative-stdev construct for a multi-observation model.
    #
    # NOT ENCODED: the correlation between the veliparib and M8 residual
    # errors. Niu 2017 Methods: "Correlation was considered between residual
    # error terms for veliparib and M8 because their quantitation by LC-MS
    # used a single shared internal standard", and Results: "The expected
    # correlation between residual errors of parent and metabolite because of
    # sampling from the same matrix was accounted for using a fixed-effect
    # correlation term". nlmixr2 has no idiomatic encoding for cross-endpoint
    # residual correlation, so the two error models are independent here; the
    # MAGNITUDES above are preserved exactly, only the coupling is lost.
    # Precedents that drop the same construct: Jonsson_2015_edoxaban.R
    # (edoxaban/M4 residual correlation 0.232), Svensson_2014_bedaquiline.R
    # (BDQ/M2 residual correlation 55%). See the vignette's Assumptions and
    # deviations section.
    # ---------------------------------------------------------------------
    propSd    <- 0.251    ; label("Veliparib proportional residual error (SD, fraction)")  # Niu 2017 Table 2 'Proportional error for veliparib (%CV), CMultStdev = 25.1%' (bootstrap median 26.0%, 95%CI 21.4-30.5%)
    addSd     <- 0.607    ; label("Veliparib additive residual error (ng/mL)")             # Niu 2017 Table 2 'Additive error for veliparib (ng/mL) = 0.607' (bootstrap median 0.183, 95%CI 0-0.607 -- poorly determined, see vignette Errata)
    propSd_m8 <- 0.202808 ; label("M8 proportional residual error (SD, fraction)")         # Niu 2017 Table 2 notes: CMultStdev * Ratio * MultStdev = 0.251 * 0.808 * 1.00; 'Ratio between veliparib and M8 proportional error, ratio = 0.808' (bootstrap median 0.806, 95%CI 0.639-0.926)
    addSd_m8  <- 3.37     ; label("M8 additive residual error (ng/mL)")                    # Niu 2017 Table 2 'Additive error for M8 (ng/mL) = 3.37' (bootstrap median 3.37, 95%CI 1.60-5.60)
  })

  model({
    # 1. Absorption. Niu 2017 Figure 1: dose enters the parent central
    #    compartment through Ka after a lag of tlag.
    ka   <- exp(lka + etalka)
    tlag <- exp(ltlag + etaltlag)

    # 2. Veliparib disposition.
    #    Eq. 7: Typical Vc/F = Vc/F * (LBM/48)^power_LBM, LBM in kg.
    vc <- exp(lvc + etalvc) * (LBM / 48)^e_lbm_vc
    vp <- exp(lvp + etalvp)
    q  <- exp(lq)  # no IIV; eta on Q removed in Supplemental Table S1 run 112

    # 3. Split of total apparent clearance. frenal is not an independent
    #    parameter -- Niu 2017 Figure 1 legend defines fm = 1 - frenal -- so
    #    it is recomputed here from fm rather than carried separately.
    #    Eq. 5: Typical CLR/F  = CL/F * frenal * (CLCR/95)^power_CLcr, CLCR in mL/min.
    #    Eq. 6: Typical CLNR/F = CL/F * fm.
    #    The two arms share the single estimated CL/F but carry SEPARATE etas,
    #    which is what Table 2's distinct "omega CLR/F" and "omega CLNR/F"
    #    rows report.
    frenal    <- 1 - fm
    cl_renal  <- exp(lcl + etalcl_renal)  * frenal * (CRCL / 95)^e_crcl_cl_renal
    cl_nonren <- exp(lcl + etalcl_nonren) * fm

    # 4. M8 disposition.
    vc_m8 <- exp(lvc_m8 + etalvc_m8)
    vp_m8 <- exp(lvp_m8 + etalvp_m8)
    q_m8  <- exp(lq_m8)  # no IIV; eta on Qmet removed in Supplemental Table S1 run 113
    cl_m8 <- exp(lcl_m8 + etalcl_m8)

    # 5. ODE system, mirroring Niu 2017 Figure 1. The non-renal (formation)
    #    clearance flux leaves veliparib central and enters M8 central 1:1 in
    #    MASS units: the molecular weights differ by only 5.7%, and Methods
    #    states no concentration-unit correction was applied.
    d/dt(depot)          <- -ka * depot
    d/dt(central)        <-  ka * depot -
                             cl_renal  * central / vc -
                             cl_nonren * central / vc -
                             q * central / vc + q * peripheral1 / vp
    d/dt(peripheral1)    <-  q * central / vc - q * peripheral1 / vp
    d/dt(central_m8)     <-  cl_nonren * central / vc -
                             cl_m8 * central_m8 / vc_m8 -
                             q_m8 * central_m8 / vc_m8 + q_m8 * peripheral1_m8 / vp_m8
    d/dt(peripheral1_m8) <-  q_m8 * central_m8 / vc_m8 - q_m8 * peripheral1_m8 / vp_m8

    # 6. Absorption lag time.
    alag(depot) <- tlag

    # 7. Observations. Both are declared BEFORE any residual-error line so the
    #    ODE states they reference are not parsed into the error block. The
    #    1000 is the ng/ug conversion described in the units comment above:
    #    state (mg) / volume (L) is mg/L == ug/mL, and the paper works in
    #    ng/mL.
    Cc    <- 1000 * central    / vc
    Cc_m8 <- 1000 * central_m8 / vc_m8

    Cc    ~ add(addSd)    + prop(propSd)
    Cc_m8 ~ add(addSd_m8) + prop(propSd_m8)
  })
}
