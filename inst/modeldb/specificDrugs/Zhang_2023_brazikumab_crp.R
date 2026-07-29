Zhang_2023_brazikumab_crp <- function() {
  description <- "Indirect-response PK/PD model of Crohn's Disease Activity Index (CDAI) for brazikumab with a constant inhibitory placebo effect and a baseline-C-reactive-protein-dependent sigmoid drug effect on the CDAI input rate, in adults with moderately to severely active Crohn's disease (Zhang 2023)"
  reference <- "Zhang N, Chan ML, Li J, Brohawn PZ, Sun B, Vainshtein I, Roskos LK, Faggioni R, Savic RM. Combining pharmacometric models with predictive and prognostic biomarkers for precision therapy in Crohn's disease: A case study of brazikumab. CPT Pharmacometrics Syst Pharmacol. 2023;12(12):1945-1959. doi:10.1002/psp4.13044"
  vignette <- "Zhang_2023_brazikumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  # The CDAI turnover state is the paper's efficacy endpoint. It is the direct
  # analogue of the registered `das28` disease-activity-score output compartment
  # but is not itself a registered canonical, so it is declared here rather than
  # silently introduced. Candidate for promotion to a canonical sibling of
  # `das28` in a future naming audit.
  paper_specific_compartments <- c("cdai")

  covariateData <- list(
    CRP = list(
      description        = "Baseline serum C-reactive protein concentration",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "BASELINE ONLY (pre-dose, day 0); not time-varying. Drives the sigmoid drug effect on",
        "the CDAI input rate via Equation 5:",
        "idrug = imax * CRP^hill / (ec50^hill + CRP^hill), with the paper's IB50 mapped to the",
        "canonical `ec50` = 8.03 mg/L and the paper's gamma mapped to `hill` = 2.07 (fixed).",
        "Unlike the IL-22 variant (hill = 20, effectively a step), hill = 2.07 gives a genuinely",
        "graded relationship, so the total inhibition varies continuously between 17.8% and",
        "42.4% with the individual BCRP level. Phase IIa median 18.2 mg/L (range 0.32-212.8) in",
        "the treatment arm and 9.55 (0.44-110.4) in the placebo arm (Table 1); pooled median",
        "15.7 mg/L. Measured by the Roche Integra immunoturbidimetric high-sensitivity assay,",
        "upper limit of normal 3.00 mg/L (Appendix S1). Note that the model-derived 8.03 mg/L",
        "cutoff is well BELOW the 15.7 mg/L cohort median, so it assigns more patients (38 vs.",
        "23) to the biomarker-high group than the median cutoff does (Discussion).",
        sep = " "
      ),
      source_name        = "BCRP"
    ),
    ON_TREATMENT = list(
      description        = "Active brazikumab treatment-arm indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (placebo arm)",
      notes              = paste(
        "Gates the drug effect entirely: the source control stream",
        "PSP-2023-0072-T-s05.mod encodes 'EFF = 0' followed by",
        "'IF(ARM.EQ.1) EFF = EMAX*BCRP**GAM/(EC50**GAM+BCRP**GAM)', so idrug is exactly zero in",
        "the placebo arm and the placebo arm's trajectory is driven by `iplac` alone.",
        "ON_TREATMENT rather than an exposure-driven effect is the correct canonical here",
        "because Zhang 2023 explicitly could not estimate a concentration-dependent drug",
        "effect (Results: only one dose level, 700 mg IV, was studied in phase IIa and the",
        "5th-95th percentile exposure range was a narrow 43.7-106.1 ug/mL, so 'the",
        "concentration dependent-drug effect for the dose level of 700 mg i.v. was considered",
        "as constant and cannot be estimated in the presence of other confounding factors').",
        "Source column ARM is coded 1 = treatment, 2 = placebo, so ON_TREATMENT = 2 - ARM;",
        "the dataset's own PLACEBO column is the exact complement (ON_TREATMENT = 1 - PLACEBO),",
        "verified against the supplement analysis dataset.",
        sep = " "
      ),
      source_name        = "ARM"
    ),
    ALB = list(
      description        = "Baseline serum albumin concentration",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form effect on CL centred at 39 g/L, carried unchanged from the final population PK model; see modellib('Zhang_2023_brazikumab') for the units caveat on the paper's 'mg/dL' label. Affects the PK layer only -- the CDAI drug effect is not exposure-driven.",
      source_name        = "BALBU"
    ),
    DIS_HEALTHY = list(
      description        = "Healthy-participant cohort indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (patient with Crohn's disease)",
      notes              = "Linear fractional effect on CL, carried unchanged from the final population PK model. Every subject in the phase IIa efficacy dataset is a patient with CD, so DIS_HEALTHY = 0 throughout that cohort; the term is retained so the PK layer stays identical to modellib('Zhang_2023_brazikumab').",
      source_name        = "GRP"
    ),
    SEXF = list(
      description        = "Female sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 (female; the reference category of the published Vc typical value)",
      notes              = "Linear fractional effect on Vc applied to the male category as (1 + 0.214 * (1 - SEXF)), carried unchanged from the final population PK model. Source column GNDR is male = 1 / female = 0, so SEXF = 1 - GNDR.",
      source_name        = "GNDR"
    )
  )

  covariatesDataExcluded <- list(
    SCORE_CDAI = list(
      description = "Baseline Crohn's Disease Activity Index (BCDAI)",
      units       = "(score, 0-600)",
      type        = "continuous",
      notes       = paste(
        "Zhang 2023's prognostic biomarker, but NOT a covariate in this model. The paper does",
        "not put a BCDAI covariate effect on any structural parameter; instead the baseline",
        "CDAI is an estimated model parameter (`lrbase`, typical value 318) carrying its own",
        "log-normal random effect, and the BCDAI/placebo relationship is expressed as the",
        "-100% CORRELATION between the baseline and placebo-effect random effects (Table 2,",
        "'Correlation between IIV of BCDAI score and placebo effect'). Figure S2 reports the",
        "post-hoc Pearson correlation between individual BCDAI and the individual",
        "model-estimated placebo effect for this BCRP variant, which is a diagnostic of that",
        "omega block rather than a covariate model. Listed here so the provenance of the",
        "paper's prognostic-biomarker claim is preserved; see the etalrbase/etaiplac block in",
        "ini().",
        sep = " "
      )
    ),
    IL22 = list(
      description = "Baseline serum interleukin-22 concentration",
      units       = "pg/mL",
      type        = "continuous",
      notes       = "The alternative predictive biomarker. Zhang 2023 fit BIL22- and BCRP-dependent drug effects as two SEPARATE models rather than jointly (Discussion: adding both 'did not result in a more statistically significant relationship ... possibly due to the small sample size of the treatment arm and that BIL22 and BCRP are 60% correlated (p = 5.332e-13)'). The BIL22 variant is modellib('Zhang_2023_brazikumab_il22')."
    ),
    FCP = list(
      description = "Baseline faecal calprotectin",
      units       = "ug/g",
      type        = "continuous",
      notes       = "Screened as a candidate predictive/prognostic biomarker in the efficacy model (Methods) and reported in Table 1 (median 639.0 ug/g, placebo arm); not statistically significant and not retained. Only BIL22, BCRP, and BCDAI met the screening criteria."
    ),
    IL17 = list(
      description = "Baseline serum interleukin-17",
      units       = "ng/dL",
      type        = "continuous",
      notes       = "Screened as a candidate biomarker (Methods 'BIL17'); Table 1 median 0.48 ng/dL. Not retained in the final model."
    ),
    LCN2 = list(
      description = "Baseline serum lipocalin-2",
      units       = "ng/mL",
      type        = "continuous",
      notes       = "Screened as a candidate biomarker (Methods 'BLCN2'); Table 1 median 214.2 ng/mL. Not retained in the final model."
    ),
    MIP3A = list(
      description = "Baseline macrophage inflammatory protein-3 alpha (CCL20)",
      units       = "pg/mL",
      type        = "continuous",
      notes       = "Screened as a candidate biomarker (Methods 'BMIP3A'); Table 1 median 22.7 pg/mL. Not retained in the final model."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 119L,
    n_studies      = 1L,
    age_range      = "18-61 years",
    age_median     = "35 years (34 treatment arm, 36 placebo arm)",
    weight_range   = "44.0-158.8 kg",
    weight_median  = "66.9 kg (65.8 treatment arm, 69.3 placebo arm)",
    sex_female_pct = 62.2,
    race_ethnicity = c(White = 93.3, Black = 5.0, Other = 1.7),
    disease_state  = "Moderately to severely active Crohn's disease in adults who had failed or were intolerant to anti-TNF-alpha therapy",
    dose_range     = "700 mg brazikumab IV over at least 60 min on day 1 and day 29, or matching placebo, during the 12-week double-blind induction period",
    regions        = "Multinational (NCT01714726)",
    trials         = "NCT01714726 (phase IIa)",
    baseline_covariates = list(
      cdai_median_all       = 317,
      cdai_median_treatment = 330,
      cdai_median_placebo   = 304,
      crp_median_all_mg_L   = 15.7,
      crp_median_treatment  = 18.2,
      crp_median_placebo    = 9.55,
      il22_median_all_pg_mL = 15.6,
      albumin_median_g_L    = 39
    ),
    notes          = paste(
      "Efficacy analysis population: the 119 patients of the phase IIa double-blind",
      "placebo-controlled induction period only (59 brazikumab, 60 placebo), with a median of 5",
      "CDAI observations per subject (range 1-5 treatment, 3-5 placebo; Table 1). Data from the",
      "100-week open-label period were EXCLUDED by the authors 'due to potential influence of",
      "unblinding on the clinical readout'. The source control stream",
      "PSP-2023-0072-T-s05.mod additionally restricts to TIME <= 90 days and drops one flagged",
      "outlier record, and fits the CDAI observations only (IGNORE(CMT.LE.3) removes the PK",
      "records), holding every PK parameter fixed. Of the 59 treated patients, 52 had a CDAI",
      "measurement around week 8 (Figure 5). Baseline demographics from Zhang 2023 Table 1;",
      "medians verified against the supplement analysis dataset PSP-2023-0072-T-s06.csv.",
      sep = " "
    )
  )

  ini({
    # ---- PK layer -------------------------------------------------------------
    # Carried from the final population PK model (Zhang 2023 Table 2, PK block);
    # see modellib('Zhang_2023_brazikumab'). Every PK parameter was held FIXED
    # while the efficacy parameters were estimated (PSP-2023-0072-T-s05.mod
    # $THETA rows 1-13 all carry FIX), so they are wrapped in fixed() here.
    # NOTE: the drug effect below is NOT exposure-driven, so the PK layer does not
    # feed the CDAI prediction; it is retained so this file simulates both the
    # concentration and the CDAI endpoint from one model.
    lcl     <- fixed(log(0.26));    label("Clearance (CL, L/day)")                                # Table 2: CL in female patients with CD = 0.26 L/d
    lvc     <- fixed(log(3.27));    label("Central volume of distribution (Vc, L)")               # Table 2: Vc in female subjects = 3.27 L
    lvp     <- fixed(log(2.64));    label("Peripheral volume of distribution (Vp, L)")            # Table 2: Vp = 2.64 L
    lq      <- fixed(log(0.412));   label("Intercompartmental clearance (Q, L/day)")              # Table 2: Q = 0.412 L/d
    lka     <- fixed(log(0.286));   label("First-order SC absorption rate constant (ka, 1/day)")  # Table 2: ka = 0.286 1/d
    ltlag   <- fixed(log(0.0296));  label("SC absorption lag time (Tlag, day)")                   # Table 2: Tlag = 0.0296 d
    lfdepot <- fixed(log(0.88));    label("SC bioavailability (F, fraction)")                     # Table 2: F = 0.88

    e_alb_cl         <- fixed(-1.32);  label("Power exponent of baseline albumin on CL (centred at 39 g/L)")  # Table 2 footnote a
    e_dis_healthy_cl <- fixed(-0.362); label("Healthy-participant fractional change in CL")                   # Table 2 footnote a
    e_male_vc        <- fixed(0.214);  label("Male fractional change in Vc vs. the female reference")         # Table 2 footnote b

    # ---- PD layer: indirect response on CDAI ----------------------------------
    lrbase    <- log(318);  label("Baseline CDAI score (BCDAI, score units)")                     # Table 2 CRP column: baseline CDAI = 318 (RSE 2%)
    lthalfrec <- log(11.7); label("CDAI remission half-life (T1/2, day)")                         # Table 2 CRP column: half-life HL = 11.7 d (RSE 8%)
    iplac     <- 0.178;     label("Constant inhibitory placebo effect on the CDAI input rate (fraction)")  # Table 2 CRP column: inhibitory placebo effect = 0.178 (RSE 16%)
    limax     <- log(0.246); label("Maximum CRP-dependent drug inhibition of the CDAI input rate (Imax, fraction)")  # Table 2 CRP column: Imax = 0.246 (RSE 10%)
    lec50     <- log(8.03); label("Baseline CRP achieving 50% of Imax (the paper's IB50, mg/L)")   # Table 2 CRP column: IB50 = 8.03 mg/L (RSE 10%)

    # Hill coefficient (the paper's gamma). Fixed at 2.07 because although "the
    # value of gamma in the initial estimation is moderate, it still caused large
    # uncertainty in the estimation of the steepness of the biomarker/drug effect
    # relationship. As such, the value of gamma was fixed as the initial estimate
    # of 2.07" (Discussion). Kept on the bare linear scale (the registered
    # canonical form) so the fixed value is exact rather than round-tripped
    # through log/exp.
    hill <- fixed(2.07); label("Hill coefficient of the CRP drug-effect sigmoid (gamma, unitless)")  # Table 2 CRP column: Hill coefficient gamma = 2.07 FIX

    # Inter-individual variability on the PD parameters. Zhang 2023 Methods:
    # BCDAI is log-normal (etalrbase enters as exp(lrbase + etalrbase)), whereas
    # "The placebo effect ... is modeled to have an additive variability (i.e.,
    # PLAC_i = PLAC_TV + eta_i)", so etaiplac enters additively.
    #
    # Table 2 reports the correlation between the two random effects as exactly
    # -100%, which makes the 2x2 block SINGULAR (determinant 0) and breaks
    # rxode2's Cholesky sampler with "chol(): decomposition failed". The
    # off-diagonal is therefore scaled by 0.99 (correlation -1.00 -> -0.99),
    # which is the smallest change that restores positive definiteness; both
    # published variances are kept exactly.
    #   full-correlation covariance = -sqrt(0.00997 * 0.0519) = -0.0227474
    etalrbase + etaiplac ~ c(
      0.00997,
      0.99 * -0.0227474, 0.0519
    )  # Table 2 CRP column: IIV of BCDAI score = 0.00997 (RSE 35%), IIV of placebo effect = 0.0519 (RSE 4%), correlation -100% (encoded as -99%; see the note above)

    # Residual error. PK: proportional only (see modellib('Zhang_2023_brazikumab')).
    # CDAI: combined proportional + additive, matching the Methods expression
    # Aobs = Apred + sqrt((Apred * eps1)^2 + eps2^2).
    propSd      <- fixed(0.249); label("Proportional residual error on brazikumab concentration (fraction)")  # Table 2 PK block: proportional error = 24.9%
    propSd_cdai <- 0.117;        label("Proportional residual error on CDAI (fraction)")                       # Table 2 CRP column: proportional error, efficacy = 0.117 (RSE 9%)
    addSd_cdai  <- 38.9;         label("Additive residual error on CDAI (score units)")                        # Table 2 CRP column: additive error, efficacy = 38.9 (RSE 8%)
  })

  model({
    # 1. Covariate-effect terms (PK layer).
    alb_cl <- (ALB / 39)^e_alb_cl
    dis_cl <- 1 + e_dis_healthy_cl * DIS_HEALTHY
    sex_vc <- 1 + e_male_vc * (1 - SEXF)

    # 2. Individual PK parameters and micro-constants.
    cl  <- exp(lcl) * alb_cl * dis_cl
    vc  <- exp(lvc) * sex_vc
    vp  <- exp(lvp)
    q   <- exp(lq)
    ka  <- exp(lka)
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 3. Individual PD parameters (Equations 1-4).
    rbase    <- exp(lrbase + etalrbase)
    thalfrec <- exp(lthalfrec)
    kout     <- log(2) / thalfrec
    kin      <- rbase * kout
    iplac_i  <- iplac + etaiplac
    imax     <- exp(limax)
    ec50     <- exp(lec50)

    # 4. Total inhibition of the CDAI input rate = placebo effect + drug effect.
    # The drug term is the Equation 5 sigmoid in the BASELINE CRP level (not in
    # concentration), and it is switched off entirely in the placebo arm.
    idrug  <- ON_TREATMENT * imax * CRP^hill / (ec50^hill + CRP^hill)
    itotal <- iplac_i + idrug

    # 5. ODE system. PK: Appendix S1 Equations 1-6. PD: indirect response with
    # inhibition of the zero-order input rate (Equations 1-4); the CDAI state
    # starts at its individual baseline, which is the steady state when itotal = 0.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(cdai)        <-  kin * (1 - itotal) - kout * cdai
    cdai(0)           <-  rbase

    # 6. Bioavailability and lag time apply to the SC depot only.
    f(depot)    <- exp(lfdepot)
    alag(depot) <- exp(ltlag)

    # 7. Observations. Dose in mg / volume in L gives mg/L = ug/mL.
    Cc <- central / vc
    Cc ~ prop(propSd)
    cdai ~ prop(propSd_cdai) + add(addSd_cdai)
  })
}
