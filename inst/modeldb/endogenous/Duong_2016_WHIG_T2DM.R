Duong_2016_WHIG_T2DM <- function() {
  description <- paste(
    "Semi-mechanistic weight-HbA1c-insulin-glucose (WHIG) disease-progression",
    "model of placebo-treated adults with type 2 diabetes mellitus (T2DM)",
    "(Duong 2016 BJCP). Structural model is the Choy 2016 WHIG framework:",
    "a body-weight turnover compartment, a steady-state closed-form solution",
    "for the fasting-serum-insulin (FSI) / fasting-plasma-glucose (FPG)",
    "homeostatic feedback pair, and three transit compartments for glycated",
    "haemoglobin (HbA1c). Insulin sensitivity (IS) is driven by change in",
    "body weight (dWGT); ScaleEFS = 0.0458 per kg scales the IS-vs-weight",
    "sensitivity. baseline b-cell function (b0, logit) is shifted between",
    "Study 1 (newly diagnosed obese T2DM, b0 = -0.298 -> BF0 ~ 57 pct) and",
    "Studies 2 and 3 (advanced T2DM, b0 = 0.677 -> BF0 ~ 34 pct) via the",
    "STUDY_1 covariate. Placebo treatment effects are step functions active",
    "at t > 0 (EFDE, diet-and-exercise + placebo run-in) and at t >",
    "placebo-run-in phase (EFPL, Study 1 only; the placebo effect at the",
    "treatment phase was removed for Studies 2 and 3 because it was not",
    "significant). rB (rate of b-cell function loss) and MTT (HbA1c mean",
    "transit time) were fixed to the Choy 2016 upstream WHIG estimates",
    "(0.209 per year and 38.9 days) because the shorter trial duration of",
    "Studies 2 and 3 rendered them poorly identifiable in the joint fit;",
    "these fixed values are re-estimated for the pooled 3-study population.",
    "The FSI / FPG closed-form solution is the quadratic root of the linked",
    "steady-state ODEs (Choy 2016 Eqs 11 and 12): C * FSI^2 + 3.5 * A * C *",
    "FSI - 35.1 * A = 0 where A = EFB * B * 7.8 (K_in,FSI / K_out,FSI = 7.8",
    "is the healthy FSI_ss anchor) and C = EFS * IS0 (35.1 =",
    "K_in,FPG / K_out,FPG is the healthy FPG_ss * FSI_ss anchor). The",
    "3.5 mmol/L subtracted from FPG is the physiological FPG floor for",
    "insulin secretion (Choy 2016 Methods)."
  )
  reference <- paste(
    "Duong JK, de Winter W, Choy S, Plock N, Naik H, Krauwinkel W,",
    "Visser SAG, Verhamme KMC, Sturkenboom MCJM, Stricker BH, Danhof M",
    "(2017). The variability in beta-cell function in placebo-treated",
    "subjects with type 2 diabetes: application of the",
    "weight-HbA1c-insulin-glucose (WHIG) model.",
    "Br J Clin Pharmacol 83(3):487-497.",
    "doi:10.1111/bcp.13144.",
    "Structural WHIG framework (weight turnover, IS-from-dWGT, b-cell",
    "logistic decay, FSI-FPG steady-state homeostasis, 3-transit HbA1c)",
    "and the fixed values rB = 0.209 per year and MTT = 38.9 days are",
    "adapted from Choy S, Kjellsson MC, Karlsson MO, de Winter W (2016).",
    "Weight-HbA1c-insulin-glucose model for describing disease",
    "progression of type 2 diabetes.",
    "CPT Pharmacometrics Syst Pharmacol 5(1):11-19.",
    "doi:10.1002/psp4.12058.",
    "PMID 26844011; PMCID PMC4728293.",
    sep = " "
  )
  vignette <- "Duong_2016_WHIG_T2DM"
  paper_specific_compartments <- c("hba1c_1", "hba1c_2", "hba1c_3")

  units <- list(
    time          = "day",
    dosing        = "none (placebo-only; no drug administration)",
    concentration = paste(
      "Weight in kg; FSI in microU/mL (equivalent to mIU/L on the SI",
      "convention used in Duong 2017 Table 2); FPG in mmol/L; HbA1c in %"
    )
  )

  covariateData <- list(
    STUDY_1 = list(
      description        = paste(
        "Binary indicator for enrolment in Study 1 (NCT00236600) of the",
        "Duong 2017 pooled analysis (newly diagnosed, treatment-naive,",
        "obese T2DM subjects with a 6-week placebo run-in and 60-week",
        "placebo treatment phase, weight-loss diet-and-exercise arm)."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (Study 2 NCT01071850 or Study 3 NCT01117584 -- advanced T2DM cohorts)",
      notes              = paste(
        "Switches two parameters in the model:",
        "(a) baseline b-cell-function logit b0 (Duong 2017 Table 4:",
        "-0.298 for Study 1 vs 0.677 for Studies 2 and 3); and",
        "(b) placebo treatment-phase effect on weight EFPL (Duong 2017",
        "Methods: EFPL was retained only for Study 1 because the",
        "additional placebo effect during the treatment phase for",
        "Studies 2 and 3 was not significant). Time-fixed per subject."
      ),
      source_name        = "STUDY"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description        = "Age at study entry (screened, not retained in final model)",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened against b0, s0, rB, EFDE, EFPL, EFB via stepwise",
        "inclusion at P < 0.01 (Duong 2017 Methods (iii) Covariate",
        "modelling). Not significant; not in the final model.",
        "Baseline demographics (Duong 2017 Table 2): Study 1 median",
        "54 years (IQR 48-60), Study 2 median 55 (48-60),",
        "Study 3 median 57 (51-63)."
      ),
      source_name        = "AGE"
    ),
    SEXF = list(
      description        = "Biological sex indicator (1 = female, 0 = male; screened, not retained in final model)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = paste(
        "Screened against baseline and treatment-effect parameters.",
        "Not significant; not in the final model. Baseline",
        "demographics (Duong 2017 Table 2): Study 1 63 pct female,",
        "Study 2 56 pct, Study 3 45 pct. Paper uses 'gender'."
      ),
      source_name        = "gender"
    ),
    T_DIAG_DIAB = list(
      description        = "Time since T2DM diagnosis at study entry (screened, not retained in final model)",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened against b-cell parameters b0, s0, rB. Not",
        "significant despite a Discussion-noted trend towards",
        "lower b-cell function with longer T2DM duration.",
        "Baseline (Duong 2017 Table 2): Study 1 not known",
        "(assumed 1 year for the covariate test), Study 2 median",
        "2.7 years (IQR 1.3-4.6), Study 3 median 5.8 years",
        "(2.9-8.8). Confounder discussed in Duong 2017 Results:",
        "diagnosis delay of up to 6 years can obscure the",
        "biological effect."
      ),
      source_name        = "DIABDUR"
    ),
    ADHERENCE_PLACEBO_PCT = list(
      description        = "Placebo pill count adherence at end of study (screened, not retained in final model)",
      units              = "%",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened against EFDE (placebo run-in effect on weight).",
        "Not significant. Baseline (Duong 2017 Results):",
        "Study 2 median 99.8 pct (range 61.3-100), Study 3",
        "median 100 pct (range 87.1-100). Screening was to",
        "verify that low adherence did not confound the placebo",
        "response."
      ),
      source_name        = "COMPL"
    ),
    RACE = list(
      description        = "Race / ethnicity category (screened, not retained in final model)",
      units              = "(categorical)",
      type               = "categorical",
      reference_category = NULL,
      notes              = paste(
        "Screened against baseline parameters b0, s0. Not",
        "significant. Underlying race distribution not tabulated",
        "in the paper; the covariate was tested at the coarse",
        "categorical level Duong 2017 Methods used (paper term",
        "'ethnicity')."
      ),
      source_name        = "ETH"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 305L,
    n_studies      = 3L,
    age_range      = paste(
      "Median 54-57 years across studies (Duong 2017 Table 2:",
      "Study 1 median 54, IQR 48-60; Study 2 median 55, IQR 48-60;",
      "Study 3 median 57, IQR 51-63)"
    ),
    weight_range   = paste(
      "Study 1 median 104.2 kg (IQR 94.3-115.4); Study 2 median",
      "79.3 kg (68.6-87.9); Study 3 median 89.0 kg (81.2-97.6)."
    ),
    sex_female_pct = 55,
    disease_state  = paste(
      "Type 2 diabetes mellitus (T2DM). Study 1: newly diagnosed",
      "obese (BMI 27-50 kg / m^2), treatment-naive, HbA1c < 10.5 pct.",
      "Study 2: mixed treatment-naive (n = 28) and prior-treatment",
      "(n = 31), BMI 20-45 kg / m^2, HbA1c 6.8-9.5 pct. Study 3:",
      "inadequately controlled on 1.5 g metformin daily for > 6 weeks,",
      "HbA1c 7.0-9.5 pct. All arms modelled are placebo-only.",
      "Distinguishing feature captured by the STUDY_1 covariate is",
      "newly-diagnosed obese (Study 1) vs advanced T2DM (Studies 2",
      "and 3)."
    ),
    dose_range     = "None (placebo-only disease-progression model, no drug).",
    regions        = "Not tabulated in Duong 2017; the underlying three trials were multicentre.",
    notes          = paste(
      "Cohort demographics from Duong 2017 Table 2. Study 1 sample",
      "size is n = 181 (66-week trial, weight-loss counselling arm);",
      "Study 2 n = 59-66 (14-week trial, stable diet-and-exercise);",
      "Study 3 n = 64-65 (14-week trial, stable diet-and-exercise",
      "with prior metformin). Data volumes: 8587 observations",
      "(Study 1), 1526 (Study 2), 1554 (Study 3). ClinicalTrials.gov",
      "identifiers: NCT00236600 (Study 1), NCT01071850 (Study 2),",
      "NCT01117584 (Study 3)."
    )
  )

  ini({
    # -------------------------------------------------------------------
    # Weight turnover (Duong 2017 Eq 2, Table 4 'Weight' block).
    # -------------------------------------------------------------------
    lwgt_bl     <- log(102)          ; label("Baseline body weight WGT (kg; population median across the 3-study pooled cohort)")            # Duong 2017 Table 4: Baseline WGT = 102 (RSE 1 pct)
    lwgt_thalf  <- log(73.9)         ; label("Body weight half-life WGT_thalf (days); WGT_kout = ln(2) / WGT_thalf, WGT_kin = WGT_kout * baseline WGT")  # Duong 2017 Table 4: WGT_thalf = 73.9 (RSE 5 pct)

    # -------------------------------------------------------------------
    # Placebo treatment-effect step magnitudes (Duong 2017 Eq 1,
    # Table 4 'Treatment effects' block). EFDE and EFPL are reported
    # in percentage-point units and enter Eq 1 as (100 - EFDE - EFPL)
    # / 100. EFLOSS is the annual counter-effect (positive = weight
    # gain back over time).
    # -------------------------------------------------------------------
    efde        <- 3.0               ; label("Placebo, diet and exercise effect at placebo run-in EFDE (pct; step at t > 0, all subjects)")            # Duong 2017 Table 4: EFDE = 3.0 (RSE 16 pct)
    efpl        <- 3.46              ; label("Placebo effect at treatment phase EFPL (pct; step at t > run-in end, Study 1 only)")                     # Duong 2017 Table 4: EFPL = 3.46 (RSE 16 pct)
    efloss      <- 3.76              ; label("Rate of loss of placebo treatment effect EFLOSS (pct / year)")                                            # Duong 2017 Table 4: EFLOSS = 3.76 (RSE 24 pct)

    # -------------------------------------------------------------------
    # Insulin-sensitivity structural parameters (Duong 2017 Eqs 3-4,
    # Table 4 'Insulin sensitivity' block).
    # -------------------------------------------------------------------
    lscale_efs  <- log(0.0458)       ; label("ScaleEFS scaling factor for weight change on IS (per kg; log-normal IIV per paper)")                     # Duong 2017 Table 4: Scale EFS = 0.0458 (RSE 9 pct)
    s0          <- 0.963             ; label("Baseline insulin-sensitivity logit s0; IS_0 = 1 / (1 + exp(s0))")                                        # Duong 2017 Table 4: s0 = 0.963 (RSE 5 pct)

    # Duong 2017 Table 4 also reports a Box-Cox shape parameter
    # theta_shape = -0.476 (RSE 15 pct) on the etas0 IIV distribution
    # (Petersson 2009 semiparametric distribution). It is NOT
    # implemented here because the paper does not print the closed
    # form of the transformation and the Duong 2017 supplement (the
    # NONMEM control stream) is not on disk; the transformation
    # affects only the tails of the etas0 distribution and does not
    # shift the population median. See vignette Errata.

    # -------------------------------------------------------------------
    # b-cell function structural parameters (Duong 2017 Eqs 5-7,
    # Table 4 'b-cell function' block). Two population values for b0
    # by study (Duong 2017 modelled ISV on b0 as a fixed effect
    # partitioning Study 1 vs Studies 2 and 3).
    # -------------------------------------------------------------------
    b0_s1       <- -0.298            ; label("Baseline b-cell function logit b0 for Study 1 (newly diagnosed obese T2DM); BF0 = 1/(1+exp(b0))")        # Duong 2017 Table 4: b0 Study 1 = -0.298 (RSE 31 pct)
    b0_s23      <- 0.677             ; label("Baseline b-cell function logit b0 for Studies 2 and 3 (advanced T2DM)")                                  # Duong 2017 Table 4: b0 Studies 2 and 3 = 0.677 (RSE 17 pct)
    rb          <- fixed(0.209)      ; label("Rate of b-cell function loss per year rB (logistic-function units; FIXED to the Choy 2016 upstream WHIG value)")  # Duong 2017 Table 4: rB = 0.209 (fixed); Choy 2016 Table 1
    efbt        <- 0.0781            ; label("Treatment effect on b-cell function EFBT (unitless; step at t > 0, all subjects)")                       # Duong 2017 Table 4: EFBT = 0.0781 (RSE 31 pct)

    # -------------------------------------------------------------------
    # HbA1c structural parameters (Duong 2017 Eq 8 + Choy 2016 Eqs
    # 14-17 for compartments 2 and 3, Duong 2017 Table 4 'HbA1c'
    # block).
    # -------------------------------------------------------------------
    hba1c_kin   <- 0.0152            ; label("HbA1c compartment-1 production rate constant HbA1c_kin (pct per day per (mmol / L FPG))")                # Duong 2017 Table 4: HbA1c kin = 0.0152 (RSE 3 pct)
    mtt         <- fixed(38.9)       ; label("HbA1c mean transit time MTT (days; FIXED to the Choy 2016 upstream WHIG value); HbA1c_kout = 3 / MTT")   # Duong 2017 Table 4: MTT = 38.9 (fixed); Choy 2016 Table 1
    lppg        <- log(0.057)        ; label("FPG-independent HbA1c production residual PPG (pct per day; log-normal IIV per paper)")                  # Duong 2017 Table 4: PPG = 0.057 (RSE 7 pct)
    scale_ppg   <- 0.967             ; label("PPG scaling factor ScalePPG active at t > 0 (post-diet-and-exercise; unitless)")                         # Duong 2017 Table 4: Scale PPG = 0.967 (RSE 1 pct)

    # -------------------------------------------------------------------
    # IIVs (Duong 2017 Table 4 'IIV' column; log-normal for baseline
    # WGT, ScaleEFS, PPG per Methods 'residual variability was
    # described with proportional error models... IIV for PPG,
    # ScaleEFs, baseline weight and residual error models were
    # log-normally distributed, while all other parameters were
    # assumed to be normally distributed'). Correlations between IIVs
    # (n = 10 full omega block, Duong 2017 Table S3) are NOT encoded
    # here because Table S3 is only in the Duong 2017 supplement,
    # which was not obtainable at extraction time -- see vignette
    # Errata for the Table 4 diagonal-only footprint.
    # -------------------------------------------------------------------
    etalwgt_bl     ~ log(0.161^2 + 1)    # log-normal, CV 16.1 pct  (Duong 2017 Table 4: baseline WGT IIV 16.1 pct)
    etalscale_efs  ~ log(0.757^2 + 1)    # log-normal, CV 75.7 pct  (Duong 2017 Table 4: ScaleEFS IIV 75.7 pct)
    etalppg        ~ log(0.256^2 + 1)    # log-normal, CV 25.6 pct  (Duong 2017 Table 4: PPG IIV 25.6 pct)
    etab0          ~ 1.13^2              # normal, absolute SD 1.13 (Duong 2017 Table 4: b0 IIV 1.13 -- footnote a reported-as-absolute)
    etarb          ~ 0.408^2             # normal, absolute SD 0.408
    etaefbt        ~ 0.053^2             # normal, absolute SD 0.053
    etas0          ~ 0.485^2             # normal, absolute SD 0.485 (Box-Cox NOT applied here -- see vignette Errata)
    etaefde        ~ 21.3^2              # normal, absolute SD 21.3 (percentage-point units)
    etaefpl        ~ 28.9^2              # normal, absolute SD 28.9
    etaefloss      ~ 70.4^2              # normal, absolute SD 70.4

    # -------------------------------------------------------------------
    # Residual error (Duong 2017 Table 4 'Residual errors' block:
    # all four biomarkers use proportional residual error models per
    # Methods 'All observations were log-transformed prior to the
    # analysis and residual variability was described with
    # proportional error models for weight, FSI, FPG and HbA1c').
    # The subject-level IIV on residual error variance (FSI 41.5 pct,
    # FPG 29 pct, HbA1c 22 pct in Duong 2017 Table 4) is NOT encoded
    # here -- see vignette Errata.
    # -------------------------------------------------------------------
    propSd_WGT   <- 0.0096           ; label("Proportional residual SD on body weight (fraction)")   # Duong 2017 Table 4: sigma_WGT = 0.0096 (RSE 1 pct)
    propSd_FSI   <- 0.265            ; label("Proportional residual SD on FSI (fraction)")            # Duong 2017 Table 4: sigma_FSI = 0.265 (RSE 4 pct)
    propSd_FPG   <- 0.0841           ; label("Proportional residual SD on FPG (fraction)")            # Duong 2017 Table 4: sigma_FPG = 0.0841 (RSE 4 pct)
    propSd_HbA1c <- 0.0254           ; label("Proportional residual SD on HbA1c (fraction)")          # Duong 2017 Table 4: sigma_HbA1c = 0.0254 (RSE 3 pct)
  })

  model({
    # -------------------------------------------------------------------
    # 1. Individual structural parameters.
    # -------------------------------------------------------------------
    wgt_baseline <- exp(lwgt_bl     + etalwgt_bl)
    wgt_thalf    <- exp(lwgt_thalf)
    wgt_kout     <- log(2) / wgt_thalf
    wgt_kin_i    <- wgt_kout * wgt_baseline

    efde_i       <- efde   + etaefde
    efpl_i       <- efpl   + etaefpl
    efloss_i     <- efloss + etaefloss

    scale_efs_i  <- exp(lscale_efs + etalscale_efs)
    s0_i         <- s0     + etas0

    # Study-1 vs Studies-2-and-3 baseline b0 shift. Parens around each
    # THETA operand block are the mu-ref-parser workaround: without
    # them the parser reads STUDY_1 * b0_s1 + (1 - STUDY_1) * b0_s23
    # as a THETA + THETA / mu-reference candidate and errors out.
    b0_tv        <- ((STUDY_1) * (b0_s1)) + (((1.0) - (STUDY_1)) * (b0_s23))
    b0_i         <- b0_tv  + etab0
    rb_i         <- rb     + etarb
    efbt_i       <- efbt   + etaefbt

    hba1c_kout   <- 3.0 / mtt
    ppg_i        <- exp(lppg + etalppg)

    # -------------------------------------------------------------------
    # 2. Time-based placebo-treatment switches (Duong 2017 Eq 1, Eq 6,
    # Methods).
    #   OC1 = 1 for t > 0 (placebo-run-in step on b-cell function; all
    #   studies).
    #   EFDE-active step at t > 0 (all studies).
    #   EFPL-active step at t > 42 days AND Study 1 (6-week placebo
    #   run-in in Study 1; EFPL removed for Studies 2 and 3 in Duong
    #   2017 Results because the additional treatment-phase placebo
    #   effect was not significant).
    # -------------------------------------------------------------------
    oc1          <- (t > 0.0)
    is_after_runin <- (t > 42.0) * STUDY_1

    # -------------------------------------------------------------------
    # 3. Weight dynamics (Duong 2017 Eqs 1-2). EFW is the fractional
    # weight-input multiplier; steady-state at EFW = 1 gives WGT_ss =
    # wgt_baseline.
    # -------------------------------------------------------------------
    efw_treat   <- oc1 * efde_i + is_after_runin * efpl_i
    efw_loss    <- efloss_i * t / 365.0
    efw         <- ((100.0 - efw_treat) / 100.0) * ((100.0 + efw_loss) / 100.0)

    d/dt(weight) <- wgt_kin_i * efw - wgt_kout * weight
    weight(0)    <- wgt_baseline

    # -------------------------------------------------------------------
    # 4. Insulin sensitivity (Duong 2017 Eqs 3-4). At t = 0, dwgt = 0,
    # EFS = 1, so IS = IS_0.
    # -------------------------------------------------------------------
    dwgt         <- weight - wgt_baseline
    efs          <- 1.0 + scale_efs_i * dwgt
    is_0         <- 1.0 / (1.0 + exp(s0_i))

    # -------------------------------------------------------------------
    # 5. b-cell function (Duong 2017 Eqs 5-7). rB is a rate per year;
    # divide by 365 for a per-day time base.
    # -------------------------------------------------------------------
    bf           <- 1.0 / (1.0 + exp(b0_i + rb_i * t / 365.0))
    efb          <- 1.0 + efbt_i * oc1
    bfunc        <- bf * efb

    # -------------------------------------------------------------------
    # 6. FSI-FPG steady-state closed-form solution (Choy 2016 Eqs 11
    # and 12 with dFSI/dt = dFPG/dt = 0; Duong 2017 Methods refers to
    # 'a closed-form solution for FSI and FPG with steady-state
    # assumptions'). Physiological anchors: K_in,FSI / K_out,FSI =
    # 7.8 (healthy FSI_ss microU / mL); K_in,FPG / K_out,FPG = 35.1
    # (product of healthy FPG_ss = 4.5 mmol / L and healthy FSI_ss =
    # 7.8 microU / mL); FPG floor 3.5 mmol / L for insulin secretion
    # (Choy 2016 Methods citing Matthews 1985 HOMA and HOMA2). Let A =
    # EFB * B * 7.8 and C = EFS * IS0; then FSI_ss satisfies the
    # quadratic C * FSI^2 + 3.5 * A * C * FSI - 35.1 * A = 0. The
    # positive root is the biologically-meaningful solution.
    # -------------------------------------------------------------------
    aa           <- efb * bf * 7.8
    cc           <- efs * is_0
    disc         <- (3.5 * aa * cc) * (3.5 * aa * cc) + 4.0 * cc * 35.1 * aa
    fsi          <- (-3.5 * aa * cc + sqrt(disc)) / (2.0 * cc)
    fpg          <- 35.1 / (cc * fsi)

    # -------------------------------------------------------------------
    # 7. HbA1c three-transit-compartment chain (Choy 2016 Eqs 14-16).
    # ScalePPG is active at t > 0 (post-diet-and-exercise reduction in
    # PPG contribution). Steady-state initial conditions at t = 0
    # use the unscaled PPG input: HbA1c_(1..3)(0) = (PPG + HbA1c_kin *
    # FPG(0)) / HbA1c_kout so the sum-of-three-compartments equals
    # the observed baseline HbA1c (~ 6.7 pct for the Study-1 median
    # baseline FPG per Duong 2017 Table 2). We evaluate FPG(0) as the
    # steady-state expression above at t = 0 (efs = 1, efb = 1,
    # b = BF0).
    # -------------------------------------------------------------------
    ppg_effective <- ppg_i * (1.0 + (scale_ppg - 1.0) * oc1)

    d/dt(hba1c_1) <- ppg_effective + hba1c_kin * fpg - hba1c_kout * hba1c_1
    d/dt(hba1c_2) <- hba1c_kout * hba1c_1 - hba1c_kout * hba1c_2
    d/dt(hba1c_3) <- hba1c_kout * hba1c_2 - hba1c_kout * hba1c_3

    # Steady-state initial conditions. bf_ss / efb_ss / efs_ss all
    # equal at t = 0 (efs = 1, efb = 1, bf = BF0). aa_ss, cc_ss,
    # disc_ss are the t = 0 values of aa, cc, disc for the
    # closed-form FSI / FPG at baseline.
    bf_ss     <- 1.0 / (1.0 + exp(b0_i))
    aa_ss     <- bf_ss * 7.8
    cc_ss     <- is_0
    disc_ss   <- (3.5 * aa_ss * cc_ss) * (3.5 * aa_ss * cc_ss) + 4.0 * cc_ss * 35.1 * aa_ss
    fsi_ss    <- (-3.5 * aa_ss * cc_ss + sqrt(disc_ss)) / (2.0 * cc_ss)
    fpg_ss    <- 35.1 / (cc_ss * fsi_ss)
    hba1c_1_ss <- (ppg_i + hba1c_kin * fpg_ss) / hba1c_kout

    hba1c_1(0) <- hba1c_1_ss
    hba1c_2(0) <- hba1c_1_ss
    hba1c_3(0) <- hba1c_1_ss

    # -------------------------------------------------------------------
    # 8. Observations. Total HbA1c is the sum of the three transit
    # compartments (Choy 2016 Eq 13). FSI is in microU / mL (Duong 2017
    # Table 2 uses microU / mL; the equivalent mIU / L SI convention
    # is numerically identical).
    # -------------------------------------------------------------------
    WGT   <- weight
    FSI   <- fsi
    FPG   <- fpg
    HbA1c <- hba1c_1 + hba1c_2 + hba1c_3

    WGT   ~ prop(propSd_WGT)
    FSI   ~ prop(propSd_FSI)
    FPG   ~ prop(propSd_FPG)
    HbA1c ~ prop(propSd_HbA1c)
  })
}
