Zhang_2023_brazikumab_il22 <- function() {
  description <- "Indirect-response PK/PD model of Crohn's Disease Activity Index (CDAI) for brazikumab with a constant inhibitory placebo effect and a baseline-serum-IL-22-dependent sigmoid drug effect on the CDAI input rate, in adults with moderately to severely active Crohn's disease (Zhang 2023)"
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
    IL22 = list(
      description        = "Baseline serum interleukin-22 concentration",
      units              = "pg/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "BASELINE ONLY (pre-dose, day 0); not time-varying. Drives the sigmoid drug effect on",
        "the CDAI input rate via Equation 5:",
        "idrug = imax * IL22^hill / (ec50^hill + IL22^hill), with the paper's IB50 mapped to the",
        "canonical `ec50` = 22.8 pg/mL and the paper's gamma mapped to `hill` = 20 (fixed).",
        "Because hill is effectively a step (see the `hill` note in ini()), IL22 acts as a",
        "near-dichotomous responder classifier at the 22.8 pg/mL cutoff. Phase IIa median 15.9",
        "pg/mL (range 1.00-711) in the treatment arm and 14.1 (1.00-170) in the placebo arm",
        "(Table 1). Measured by the R&D Systems Human IL-22 Quantikine ELISA (D2200),",
        "quantifiable range 10-800 pg/mL in 100% serum (Appendix S1).",
        sep = " "
      ),
      source_name        = "BIL22"
    ),
    ON_TREATMENT = list(
      description        = "Active brazikumab treatment-arm indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (placebo arm)",
      notes              = paste(
        "Gates the drug effect entirely: the source control stream",
        "PSP-2023-0072-T-s04.mod encodes 'EFF = 0' followed by",
        "'IF(ARM.EQ.1) EFF = EMAX*BIL22**GAM/(EC50**GAM+BIL22**GAM)', so idrug is exactly zero",
        "in the placebo arm and the placebo arm's trajectory is driven by `iplac` alone.",
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
        "'Correlation between IIV of BCDAI score and placebo effect'). Figure 3 reports a",
        "post-hoc Pearson correlation between individual BCDAI and the individual",
        "model-estimated placebo effect (p < 0.001), which is a diagnostic of that omega",
        "block rather than a covariate model. Listed here so the provenance of the paper's",
        "prognostic-biomarker claim is preserved; see the etalrbase/etaiplac block in ini().",
        sep = " "
      )
    ),
    CRP = list(
      description = "Baseline C-reactive protein concentration",
      units       = "mg/L",
      type        = "continuous",
      notes       = "The alternative predictive biomarker. Zhang 2023 fit BIL22- and BCRP-dependent drug effects as two SEPARATE models rather than jointly (Discussion: adding both 'did not result in a more statistically significant relationship ... possibly due to the small sample size of the treatment arm and that BIL22 and BCRP are 60% correlated (p = 5.332e-13)'). The BCRP variant is modellib('Zhang_2023_brazikumab_crp')."
    ),
    FCP = list(
      description = "Baseline faecal calprotectin",
      units       = "ug/g",
      type        = "continuous",
      notes       = "Screened as a candidate predictive/prognostic biomarker in the efficacy model (Methods) and reported in Table 1 (median 628.5 ug/g, treatment arm); not statistically significant and not retained. Only BIL22, BCRP, and BCDAI met the screening criteria."
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
      notes       = "Screened as a candidate biomarker (Methods 'BLCN2'); Table 1 median 215.4 ng/mL. Not retained in the final model."
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
      il22_median_all_pg_mL  = 15.6,
      il22_median_treatment  = 15.9,
      il22_median_placebo    = 14.1,
      crp_median_all_mg_L    = 15.7,
      albumin_median_g_L     = 39
    ),
    notes          = paste(
      "Efficacy analysis population: the 119 patients of the phase IIa double-blind",
      "placebo-controlled induction period only (59 brazikumab, 60 placebo), with a median of 5",
      "CDAI observations per subject (range 1-5 treatment, 3-5 placebo; Table 1). Data from the",
      "100-week open-label period were EXCLUDED by the authors 'due to potential influence of",
      "unblinding on the clinical readout'. The source control stream",
      "PSP-2023-0072-T-s04.mod additionally restricts to TIME <= 90 days and drops one flagged",
      "outlier record, and fits the CDAI observations only (IGNORE(CMT.LE.3) removes the PK",
      "records), holding every PK parameter fixed. Of the 59 treated patients, 52 had a CDAI",
      "measurement around week 8 (Figure 5). Baseline demographics from Zhang 2023 Table 1;",
      "medians verified against the supplement analysis dataset PSP-2023-0072-T-s03.csv.",
      sep = " "
    )
  )

  ini({
    # ---- PK layer -------------------------------------------------------------
    # Carried from the final population PK model (Zhang 2023 Table 2, PK block);
    # see modellib('Zhang_2023_brazikumab'). Every PK parameter was held FIXED
    # while the efficacy parameters were estimated (PSP-2023-0072-T-s04.mod
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
    lrbase    <- log(318);  label("Baseline CDAI score (BCDAI, score units)")                     # Table 2 IL-22 column: baseline CDAI = 318 (RSE 3%)
    lthalfrec <- log(11.6); label("CDAI remission half-life (T1/2, day)")                         # Table 2 IL-22 column: half-life HL = 11.6 d (RSE 33%)
    iplac     <- 0.209;     label("Constant inhibitory placebo effect on the CDAI input rate (fraction)")  # Table 2 IL-22 column: inhibitory placebo effect = 0.209 (RSE 7%)
    limax     <- log(0.297); label("Maximum IL-22-dependent drug inhibition of the CDAI input rate (Imax, fraction)")  # Table 2 IL-22 column: Imax = 0.297 (RSE 30%)
    lec50     <- log(22.8); label("Baseline IL-22 achieving 50% of Imax (the paper's IB50, pg/mL)")  # Table 2 IL-22 column: IB50 = 22.8 pg/mL (RSE 10%)

    # Hill coefficient (the paper's gamma). Fixed at 20 because it could not be
    # estimated with acceptable precision: "a value of gamma around 20 can be
    # easily affected in the estimation step and thus results in a large
    # uncertainty in the initial estimation, and it was fixed as 20 in the final
    # biomarker/efficacy model" (Discussion). gamma = 20 makes the sigmoid
    # effectively a step at ec50 -- the paper reads it exactly that way ("This
    # large gamma value indicates a steep relationship, which approximates a
    # dichotomous function"). Kept on the bare linear scale (the registered
    # canonical form) so the fixed value is exact rather than round-tripped
    # through log/exp.
    hill <- fixed(20); label("Hill coefficient of the IL-22 drug-effect sigmoid (gamma, unitless)")  # Table 2 IL-22 column: Hill coefficient gamma = 20 FIX

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
    #   full-correlation covariance = -sqrt(0.0106 * 0.0509) = -0.0232280
    etalrbase + etaiplac ~ c(
      0.0106,
      0.99 * -0.0232280, 0.0509
    )  # Table 2 IL-22 column: IIV of BCDAI score = 0.0106 (RSE 74%), IIV of placebo effect = 0.0509 (RSE 26%), correlation -100% (encoded as -99%; see the note above)

    # Residual error. PK: proportional only (see modellib('Zhang_2023_brazikumab')).
    # CDAI: combined proportional + additive, matching the Methods expression
    # Aobs = Apred + sqrt((Apred * eps1)^2 + eps2^2).
    propSd      <- fixed(0.249); label("Proportional residual error on brazikumab concentration (fraction)")  # Table 2 PK block: proportional error = 24.9%
    propSd_cdai <- 0.116;        label("Proportional residual error on CDAI (fraction)")                       # Table 2 IL-22 column: proportional error, efficacy = 0.116 (RSE 42%)
    addSd_cdai  <- 38.9;         label("Additive residual error on CDAI (score units)")                        # Table 2 IL-22 column: additive error, efficacy = 38.9 (RSE 28%)
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
    # The drug term is the Equation 5 sigmoid in the BASELINE IL-22 level (not in
    # concentration), and it is switched off entirely in the placebo arm.
    idrug  <- ON_TREATMENT * imax * IL22^hill / (ec50^hill + IL22^hill)
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
