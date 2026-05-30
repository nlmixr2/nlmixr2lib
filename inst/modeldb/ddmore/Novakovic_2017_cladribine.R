Novakovic_2017_cladribine <- function() {
  description <- "Item Response Theory (IRT) model of EDSS disability progression in patients with multiple sclerosis treated with cladribine (Novakovic 2017). Eight EDSS functional system subscores (Pyramidal, Cerebellar, Brainstem, Sensory, Bowel/Bladder, Visual, Mental, Ambulation) are linked to a latent disability variable that follows a power-law disease progression in time. The model embeds an exposure-dependent symptomatic drug effect (Emax on cumulative cladribine dose adjusted for creatinine clearance) and an exposure-independent fractional protective effect on disease progression, plus full Random Effects on covariates (FREM) for Age, months since diagnosis (MSD), and exacerbation rate baseline (EXNB)."
  reference <- paste(
    "Novakovic AM, Krekels EHJ, Munafo A, Ueckert S, Karlsson MO. (2017).",
    "Application of item response theory to modeling of expanded disability status scale",
    "in clinical trials. AAPS J 19(1):172-179.",
    "doi:10.1208/s12248-016-9977-z.",
    "DDMORE Foundation Model Repository: DDMODEL00000223.",
    sep = " "
  )
  vignette <- "Novakovic_2017_cladribine"
  paper_specific_etas <- c("etap1", "etap2", "etap3", "etap4", "etap5")

  units <- list(
    time          = "day",
    dosing        = "mg",
    concentration = "mg/(no concentration output - IRT model; primary outputs are unitless EDSS subscore probabilities and a latent disability variable rather than drug concentration)",
    response      = "EDSS subscore (unitless ordered-categorical; 0-5 for Pyramidal/Cerebellar/Bowel-Bladder, 0-4 for Brainstem/Mental, 0-6 for Sensory/Visual, 0-9 for Ambulation), latent disability (unitless IRT scale)"
  )
  ddmore_id    <- "DDMODEL00000223"
  replicate_of <- NULL

  covariateData <- list(
    TRT = list(
      description        = "Treatment-cohort indicator. 0 = placebo, 1 = cladribine 3.5 mg/kg cumulative dose cohort, 2 = cladribine 5.25 mg/kg cumulative dose cohort.",
      units              = "(categorical)",
      type               = "categorical",
      reference_category = 0,
      notes              = "Gates the symptomatic and protective drug-effect terms via `TRT >= 1` (active treatment) in conjunction with the time > 0 condition. The categorical level (1 vs 2) is informational; the dose-response is driven by the time-varying cumulative-dose covariate CD, not by TRT itself.",
      source_name        = "TRT"
    ),
    CD = list(
      description        = "Time-varying cumulative cladribine dose administered to date (mg total dose, not body-weight-normalized).",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Substitute for cladribine PK exposure. Rises stepwise across the dosing schedule (CLARITY-style 3.5 or 5.25 mg/kg cumulative-dose regimens delivered as short oral pulses across two years) and stays constant between dosing weeks. Combined with capped creatinine clearance (CRL = min(CRCL, 150)) into the exposure surrogate `EXPS = CD * 104.5 / CRL` that drives the symptomatic Emax effect on disease progression.",
      source_name        = "CD"
    ),
    CRCL = list(
      description        = "Creatinine clearance.",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Used to scale the cumulative-dose-based exposure surrogate via a capped form: CRL = min(CRCL, 150). The constant 104.5 mL/min in the EXPS calculation is the typical-population CRCL anchor. The Novakovic 2017 study did not BSA-normalize CRCL, so this entry stores raw mL/min rather than the BSA-normalized canonical form.",
      source_name        = "CRCL"
    )
  )

  population <- list(
    n_subjects     = NA_integer_,
    n_studies      = NA_integer_,
    age_range      = "FREM-modeled, mean = 38.6 years",
    weight_range   = NA_character_,
    sex_female_pct = NA_real_,
    disease_state  = "Adults with relapsing-remitting multiple sclerosis (RRMS). The model was developed using pooled data from the CLARITY (cladribine) phase III trial program.",
    dose_range     = "Two cladribine cumulative-dose cohorts: 3.5 mg/kg and 5.25 mg/kg over two years, plus a placebo arm. Dose is delivered as short oral pulses; the model treats cumulative dose (CD) as the exposure surrogate rather than modelling cladribine PK explicitly.",
    regions        = NA_character_,
    notes          = paste(
      "Per-trial demographic detail (n_subjects, weight, sex distribution, geographic mix) is",
      "not reproduced in the DDMORE Foundation Model Repository bundle for DDMODEL00000223,",
      "and the Novakovic 2017 publication is not on disk in this worktree. The model includes",
      "Full Random Effects on covariates (FREM) for Age, months since diagnosis (MSD), and the",
      "annualised exacerbation-rate baseline (EXNB) - the FREM means are exposed as parameters",
      "(`age_mean`, `msd_mean`, `exnb_mean`) and the FREM continuous-covariate observations are",
      "exposed as outputs (`age_pred`, `msd_pred`, `exnb_pred`) for completeness. See the",
      "validation vignette's Assumptions and deviations section for the full caveat list."
    )
  )

  ini({
    # ----------------------------------------------------------------------------------------------
    # All values are FINAL ESTIMATES from the
    # `Output_real_Novakovic_2016_multiplesclerosis_cladribine_irt.lst` FINAL PARAMETER ESTIMATE
    # block (TH 1-75 line at .lst lines 1061-1067, OMEGA at 1077-1092, SIGMA at 1102). See the
    # vignette's Assumptions and deviations section: the .lst reports MINIMIZATION TERMINATED DUE
    # TO ROUNDING ERRORS (sig digits 0.6), so these values are taken as the deposited final
    # estimates and the convergence status is documented in the vignette rather than treated as
    # a defect of the bundle. The 8 EDSS functional-system subscores follow the standard EDSS
    # Kurtzke (1983) categorisation.
    # ----------------------------------------------------------------------------------------------

    # Disease-progression structural parameters. Latent disability is dimensionless; time enters
    # as (t/365)^prog_power, so prog_slope is in latent-units/(year^prog_power).
    prog_slope <- 0.0870 ; label("Disease-progression slope (per year^prog_power)")  # TH 1, .lst line 1061
    prog_power <- 0.707  ; label("Disease-progression time-power exponent (unitless)")  # TH 2, .lst line 1061

    # Drug-effect parameters. The symptomatic effect is Emax-style on the cumulative-dose /
    # capped-CRCL exposure surrogate; the protective effect is exposure-independent (a constant
    # fractional reduction of the disease-progression term while on treatment). Emax is
    # log-transformed (paired with the multiplicative eta on Emax via `exp(lemax + etalemax)`)
    # so that nlmixr2's checkModelConventions() sees the canonical eta<param>/<param> pairing.
    lemax    <- log(0.171) ; label("log(Emax) of the symptomatic drug effect (latent disability units)")  # TH 18 (back-transformed value)
    ec50     <- 407        ; label("EC50 of the symptomatic drug effect on the exposure surrogate (mg per mL/min, scaled)")  # TH 19
    prot_eff <- 0.228      ; label("Fractional protective effect on disease progression while on treatment (0-1)")  # TH 20

    # FREM continuous-covariate population means. The .mod treats Age / MSD / EXNB rows
    # (RTYPE 1/2/3) as additional observations whose predicted value is `mean + P_i + EPS(1)`
    # with EPS(1) variance fixed at 1e-5 - i.e., essentially deterministic FREM observations.
    age_mean  <- 38.6  ; label("FREM population-mean baseline age (years)")  # TH 21
    msd_mean  <- 8.74  ; label("FREM population-mean months since multiple-sclerosis diagnosis")  # TH 22
    exnb_mean <- 1.35  ; label("FREM population-mean annualised exacerbation rate at baseline")  # TH 23

    # Per-item IRT thresholds and slopes. Each EDSS functional-system subscore is treated as
    # an ordered-categorical observation with category-specific cumulative-logit thresholds
    # b_<item>_1..b_<item>_K and a single discrimination slope a_<item>. The first threshold
    # (b_<item>_1) can be negative (indicates the latent value at which P(score >= 1) = 0.5);
    # subsequent thresholds are reported as positive *increments* in the source .mod.

    # Item 1: Pyramidal (EDSS scale 0-5). Source THETA(24)..THETA(29).
    b_pyr_1 <- -1.57 ; label("Pyramidal item: cumulative-logit threshold for score 1 (latent units)")  # TH 24
    b_pyr_2 <-  1.27 ; label("Pyramidal item: incremental threshold from 1 to 2 (latent units)")        # TH 25
    b_pyr_3 <-  0.823; label("Pyramidal item: incremental threshold from 2 to 3 (latent units)")        # TH 26
    b_pyr_4 <-  1.43 ; label("Pyramidal item: incremental threshold from 3 to 4 (latent units)")        # TH 27
    b_pyr_5 <-  1.31 ; label("Pyramidal item: incremental threshold from 4 to 5 (latent units)")        # TH 28
    a_pyr   <-  3.16 ; label("Pyramidal item: discrimination slope (1/latent unit)")                    # TH 29

    # Item 2: Cerebellar (EDSS scale 0-5). Source THETA(30)..THETA(35).
    b_cer_1 <- -0.929; label("Cerebellar item: cumulative-logit threshold for score 1 (latent units)")  # TH 30
    b_cer_2 <-  0.979; label("Cerebellar item: incremental threshold from 1 to 2 (latent units)")        # TH 31
    b_cer_3 <-  1.03 ; label("Cerebellar item: incremental threshold from 2 to 3 (latent units)")        # TH 32
    b_cer_4 <-  1.64 ; label("Cerebellar item: incremental threshold from 3 to 4 (latent units)")        # TH 33
    b_cer_5 <-  1.11 ; label("Cerebellar item: incremental threshold from 4 to 5 (latent units)")        # TH 34
    a_cer   <-  2.81 ; label("Cerebellar item: discrimination slope (1/latent unit)")                    # TH 35

    # Item 3: Brainstem (EDSS scale 0-4). Source THETA(36)..THETA(40).
    b_bs_1  <- -0.121; label("Brainstem item: cumulative-logit threshold for score 1 (latent units)")  # TH 36
    b_bs_2  <-  1.74 ; label("Brainstem item: incremental threshold from 1 to 2 (latent units)")        # TH 37
    b_bs_3  <-  2.02 ; label("Brainstem item: incremental threshold from 2 to 3 (latent units)")        # TH 38
    b_bs_4  <-  2.81 ; label("Brainstem item: incremental threshold from 3 to 4 (latent units)")        # TH 39
    a_bs    <-  1.02 ; label("Brainstem item: discrimination slope (1/latent unit)")                    # TH 40

    # Item 4: Sensory (EDSS scale 0-6). Source THETA(41)..THETA(47).
    b_sen_1 <- -0.799; label("Sensory item: cumulative-logit threshold for score 1 (latent units)")  # TH 41
    b_sen_2 <-  1.20 ; label("Sensory item: incremental threshold from 1 to 2 (latent units)")        # TH 42
    b_sen_3 <-  1.97 ; label("Sensory item: incremental threshold from 2 to 3 (latent units)")        # TH 43
    b_sen_4 <-  2.39 ; label("Sensory item: incremental threshold from 3 to 4 (latent units)")        # TH 44
    b_sen_5 <-  3.02 ; label("Sensory item: incremental threshold from 4 to 5 (latent units)")        # TH 45
    b_sen_6 <-  2.35 ; label("Sensory item: incremental threshold from 5 to 6 (latent units)")        # TH 46
    a_sen   <-  0.980; label("Sensory item: discrimination slope (1/latent unit)")                    # TH 47

    # Item 5: Bowel/Bladder (EDSS scale 0-5). Source THETA(48)..THETA(53).
    b_bb_1  <- -0.151; label("Bowel/Bladder item: cumulative-logit threshold for score 1 (latent units)")  # TH 48
    b_bb_2  <-  1.90 ; label("Bowel/Bladder item: incremental threshold from 1 to 2 (latent units)")        # TH 49
    b_bb_3  <-  1.54 ; label("Bowel/Bladder item: incremental threshold from 2 to 3 (latent units)")        # TH 50
    b_bb_4  <-  3.03 ; label("Bowel/Bladder item: incremental threshold from 3 to 4 (latent units)")        # TH 51
    b_bb_5  <-  0.962; label("Bowel/Bladder item: incremental threshold from 4 to 5 (latent units)")        # TH 52
    a_bb    <-  1.25 ; label("Bowel/Bladder item: discrimination slope (1/latent unit)")                    # TH 53

    # Item 6: Visual (EDSS scale 0-6). Source THETA(54)..THETA(60).
    b_vis_1 <-  0.0351; label("Visual item: cumulative-logit threshold for score 1 (latent units)")  # TH 54
    b_vis_2 <-  3.87  ; label("Visual item: incremental threshold from 1 to 2 (latent units)")        # TH 55
    b_vis_3 <-  2.76  ; label("Visual item: incremental threshold from 2 to 3 (latent units)")        # TH 56
    b_vis_4 <-  1.76  ; label("Visual item: incremental threshold from 3 to 4 (latent units)")        # TH 57
    b_vis_5 <-  1.77  ; label("Visual item: incremental threshold from 4 to 5 (latent units)")        # TH 58
    b_vis_6 <-  1.75  ; label("Visual item: incremental threshold from 5 to 6 (latent units)")        # TH 59
    a_vis   <-  0.424 ; label("Visual item: discrimination slope (1/latent unit)")                    # TH 60

    # Item 7: Mental (EDSS scale 0-4). Source THETA(61)..THETA(65).
    b_men_1 <-  0.400 ; label("Mental item: cumulative-logit threshold for score 1 (latent units)")  # TH 61
    b_men_2 <-  1.12  ; label("Mental item: incremental threshold from 1 to 2 (latent units)")        # TH 62
    b_men_3 <-  4.17  ; label("Mental item: incremental threshold from 2 to 3 (latent units)")        # TH 63
    b_men_4 <-  2.90  ; label("Mental item: incremental threshold from 3 to 4 (latent units)")        # TH 64
    a_men   <-  0.907 ; label("Mental item: discrimination slope (1/latent unit)")                    # TH 65

    # Item 8: Ambulation (EDSS scale 0-9). Source THETA(66)..THETA(75).
    b_amb_1 <-  1.14  ; label("Ambulation item: cumulative-logit threshold for score 1 (latent units)")  # TH 66
    b_amb_2 <-  0.243 ; label("Ambulation item: incremental threshold from 1 to 2 (latent units)")        # TH 67
    b_amb_3 <-  0.225 ; label("Ambulation item: incremental threshold from 2 to 3 (latent units)")        # TH 68
    b_amb_4 <-  0.423 ; label("Ambulation item: incremental threshold from 3 to 4 (latent units)")        # TH 69
    b_amb_5 <-  0.439 ; label("Ambulation item: incremental threshold from 4 to 5 (latent units)")        # TH 70
    b_amb_6 <-  0.477 ; label("Ambulation item: incremental threshold from 5 to 6 (latent units)")        # TH 71
    b_amb_7 <-  0.139 ; label("Ambulation item: incremental threshold from 6 to 7 (latent units)")        # TH 72
    b_amb_8 <-  0.257 ; label("Ambulation item: incremental threshold from 7 to 8 (latent units)")        # TH 73
    b_amb_9 <-  0.349 ; label("Ambulation item: incremental threshold from 8 to 9 (latent units)")        # TH 74
    a_amb   <-  3.71  ; label("Ambulation item: discrimination slope (1/latent unit)")                    # TH 75

    # Correlated 5-dimensional latent random-effect vector. The source .mod implements this via
    # a Cholesky decomposition (.mod $PRED block lines 35-65) of a 5x5 covariance matrix whose
    # diagonals are var_p1..var_p5 (THETA 3-7) and whose 10 off-diagonals are
    # cov_ij = sqrt(var_i * var_j) * cor_ij with cor_ij from THETA 8-17. Translated here as a
    # direct correlated-IIV BLOCK so that nlmixr2's mu-reference parser sees one eta per latent
    # dimension (etap1..etap5) rather than a non-mu-referenced linear combination.
    #
    # var_p1   = 1.000 FIX                      # TH 3 (FIX, identifiability anchor)
    # var_p2   = 0.199                          # TH 4
    # var_p3   = 99.5                           # TH 5
    # var_p4   = 54.2                           # TH 6
    # var_p5   = 0.371                          # TH 7
    # cov_ij   = sqrt(var_i * var_j) * cor_ij with cor_ij from TH 8-17:
    #   cov_12 = sqrt(1.000 *  0.199) *  0.114  = 0.05085      # TH 8 (Dis-Slope)
    #   cov_13 = sqrt(1.000 * 99.500) *  0.265  = 2.6434       # TH 9 (Dis-Age)
    #   cov_23 = sqrt(0.199 * 99.500) *  0.118  = 0.5250       # TH 10 (Slope-Age)
    #   cov_14 = sqrt(1.000 * 54.200) *  0.272  = 2.0024       # TH 11 (Dis-MSD)
    #   cov_24 = sqrt(0.199 * 54.200) *  0.0894 = 0.2936       # TH 12 (Slope-MSD)
    #   cov_34 = sqrt(99.500 * 54.200) *  0.456  = 33.4800     # TH 13 (Age-MSD)
    #   cov_15 = sqrt(1.000 *  0.371) *  0.0454 = 0.02765      # TH 14 (Dis-EXNB)
    #   cov_25 = sqrt(0.199 *  0.371) *  0.0685 = 0.01862      # TH 15 (Slope-EXNB)
    #   cov_35 = sqrt(99.500 *  0.371) * -0.0911 = -0.5534     # TH 16 (Age-EXNB)
    #   cov_45 = sqrt(54.200 *  0.371) * -0.120  = -0.5383     # TH 17 (MSD-EXNB)
    # Lower-triangle row-major: var1, cov12, var2, cov13, cov23, var3, cov14, cov24, cov34, var4,
    # cov15, cov25, cov35, cov45, var5.
    etap1 + etap2 + etap3 + etap4 + etap5 ~ c(
      1.000,
      0.05085,  0.199,
      2.6434,   0.5250,   99.500,
      2.0024,   0.2936,   33.4800,  54.200,
      0.02765,  0.01862,  -0.5534,  -0.5383,  0.371
    )
    # Sixth random effect on the symptomatic Emax (multiplicative log-normal). Source OMEGA(6,6)
    # = 2.20 (.lst line 1092). Mu-references emax via `efsm_ind <- emax * exp(etalemax)`.
    etalemax  ~ 2.20

    # Residual-error parameters for the FREM continuous-covariate observations (Age, MSD, EXNB).
    # The source uses a single shared EPS(1) with SIGMA fixed at 1e-5, i.e., essentially
    # deterministic FREM observations. Translated as additive error with SD = sqrt(1e-5) on each
    # of the three FREM outputs, kept FIXED so re-fit attempts do not re-estimate them. Names
    # follow the multi-output residual-error convention `addSd_<output>`.
    addSd_age_pred  <- fixed(0.00316) ; label("Additive residual SD on FREM age observation (years)")
    addSd_msd_pred  <- fixed(0.00316) ; label("Additive residual SD on FREM MSD observation (months)")
    addSd_exnb_pred <- fixed(0.00316) ; label("Additive residual SD on FREM EXNB observation")
  })

  model({
    # ----------------------------------------------------------------------------------------------
    # 1. Symptomatic and protective drug-effect terms.
    # The symptomatic effect is Emax-style on the cumulative-dose / capped-CRCL exposure
    # surrogate; both effects are zero for placebo subjects (TRT == 0) and for any record at
    # baseline (t == 0). The CRL = min(CRCL, 150) cap mirrors the .mod IF(CRCL.GT.150) CRL=150.
    # ----------------------------------------------------------------------------------------------
    efsm_ind <- exp(lemax + etalemax)
    crl      <- CRCL - (CRCL - 150) * (CRCL > 150)
    exps     <- CD * 104.5 / crl
    on_trt   <- (TRT >= 1) * (t > 0)
    efss     <- on_trt * efsm_ind * exps / (exps + ec50)
    efpp     <- on_trt * prot_eff

    # 2. Latent disease-progression trajectory. The (t/365)^prog_power term scales day-units
    # to year-units before the power; matches the source PD = P1 + (slope+P2)*(TIME/365)^power.
    # Split across simple lines so nlmixr2's mu-reference parser sees at most one population
    # parameter per expression that contains an eta-derived variable.
    pd_baseline    <- etap1
    prog_slope_ind <- prog_slope + etap2
    prog_time      <- (t / 365)^prog_power
    pd_prog_term   <- prog_slope_ind * prog_time * (1 - efpp)
    pd             <- pd_baseline + pd_prog_term - efss

    # 3. Per-item IRT logits and category probabilities. Each item's cumulative-logit at
    # category boundary k is `a_<item> * (pd - sum(b_<item>_1..b_<item>_k))`. The per-category
    # probability follows the standard ordered-categorical decomposition
    # P(score = k) = P(score >= k) - P(score >= k+1) with P(score >= 0) = 1 and
    # P(score >= maxscore+1) = 0. The expected score reduces to the sum of survival
    # probabilities sum(P(score >= k)) for k = 1..maxscore.

    # Item 1: Pyramidal (0-5)
    pyr_cum1  <- b_pyr_1
    pyr_cum2  <- pyr_cum1 + b_pyr_2
    pyr_cum3  <- pyr_cum2 + b_pyr_3
    pyr_cum4  <- pyr_cum3 + b_pyr_4
    pyr_cum5  <- pyr_cum4 + b_pyr_5
    pge_pyr_1 <- 1 / (1 + exp(-a_pyr * (pd - pyr_cum1)))
    pge_pyr_2 <- 1 / (1 + exp(-a_pyr * (pd - pyr_cum2)))
    pge_pyr_3 <- 1 / (1 + exp(-a_pyr * (pd - pyr_cum3)))
    pge_pyr_4 <- 1 / (1 + exp(-a_pyr * (pd - pyr_cum4)))
    pge_pyr_5 <- 1 / (1 + exp(-a_pyr * (pd - pyr_cum5)))
    pyramidal <- pge_pyr_1 + pge_pyr_2 + pge_pyr_3 + pge_pyr_4 + pge_pyr_5

    # Item 2: Cerebellar (0-5)
    cer_cum1  <- b_cer_1
    cer_cum2  <- cer_cum1 + b_cer_2
    cer_cum3  <- cer_cum2 + b_cer_3
    cer_cum4  <- cer_cum3 + b_cer_4
    cer_cum5  <- cer_cum4 + b_cer_5
    pge_cer_1 <- 1 / (1 + exp(-a_cer * (pd - cer_cum1)))
    pge_cer_2 <- 1 / (1 + exp(-a_cer * (pd - cer_cum2)))
    pge_cer_3 <- 1 / (1 + exp(-a_cer * (pd - cer_cum3)))
    pge_cer_4 <- 1 / (1 + exp(-a_cer * (pd - cer_cum4)))
    pge_cer_5 <- 1 / (1 + exp(-a_cer * (pd - cer_cum5)))
    cerebellar <- pge_cer_1 + pge_cer_2 + pge_cer_3 + pge_cer_4 + pge_cer_5

    # Item 3: Brainstem (0-4)
    bs_cum1  <- b_bs_1
    bs_cum2  <- bs_cum1 + b_bs_2
    bs_cum3  <- bs_cum2 + b_bs_3
    bs_cum4  <- bs_cum3 + b_bs_4
    pge_bs_1 <- 1 / (1 + exp(-a_bs * (pd - bs_cum1)))
    pge_bs_2 <- 1 / (1 + exp(-a_bs * (pd - bs_cum2)))
    pge_bs_3 <- 1 / (1 + exp(-a_bs * (pd - bs_cum3)))
    pge_bs_4 <- 1 / (1 + exp(-a_bs * (pd - bs_cum4)))
    brainstem <- pge_bs_1 + pge_bs_2 + pge_bs_3 + pge_bs_4

    # Item 4: Sensory (0-6)
    sen_cum1  <- b_sen_1
    sen_cum2  <- sen_cum1 + b_sen_2
    sen_cum3  <- sen_cum2 + b_sen_3
    sen_cum4  <- sen_cum3 + b_sen_4
    sen_cum5  <- sen_cum4 + b_sen_5
    sen_cum6  <- sen_cum5 + b_sen_6
    pge_sen_1 <- 1 / (1 + exp(-a_sen * (pd - sen_cum1)))
    pge_sen_2 <- 1 / (1 + exp(-a_sen * (pd - sen_cum2)))
    pge_sen_3 <- 1 / (1 + exp(-a_sen * (pd - sen_cum3)))
    pge_sen_4 <- 1 / (1 + exp(-a_sen * (pd - sen_cum4)))
    pge_sen_5 <- 1 / (1 + exp(-a_sen * (pd - sen_cum5)))
    pge_sen_6 <- 1 / (1 + exp(-a_sen * (pd - sen_cum6)))
    sensory <- pge_sen_1 + pge_sen_2 + pge_sen_3 + pge_sen_4 + pge_sen_5 + pge_sen_6

    # Item 5: Bowel/Bladder (0-5)
    bb_cum1  <- b_bb_1
    bb_cum2  <- bb_cum1 + b_bb_2
    bb_cum3  <- bb_cum2 + b_bb_3
    bb_cum4  <- bb_cum3 + b_bb_4
    bb_cum5  <- bb_cum4 + b_bb_5
    pge_bb_1 <- 1 / (1 + exp(-a_bb * (pd - bb_cum1)))
    pge_bb_2 <- 1 / (1 + exp(-a_bb * (pd - bb_cum2)))
    pge_bb_3 <- 1 / (1 + exp(-a_bb * (pd - bb_cum3)))
    pge_bb_4 <- 1 / (1 + exp(-a_bb * (pd - bb_cum4)))
    pge_bb_5 <- 1 / (1 + exp(-a_bb * (pd - bb_cum5)))
    bowel_bladder <- pge_bb_1 + pge_bb_2 + pge_bb_3 + pge_bb_4 + pge_bb_5

    # Item 6: Visual (0-6)
    vis_cum1  <- b_vis_1
    vis_cum2  <- vis_cum1 + b_vis_2
    vis_cum3  <- vis_cum2 + b_vis_3
    vis_cum4  <- vis_cum3 + b_vis_4
    vis_cum5  <- vis_cum4 + b_vis_5
    vis_cum6  <- vis_cum5 + b_vis_6
    pge_vis_1 <- 1 / (1 + exp(-a_vis * (pd - vis_cum1)))
    pge_vis_2 <- 1 / (1 + exp(-a_vis * (pd - vis_cum2)))
    pge_vis_3 <- 1 / (1 + exp(-a_vis * (pd - vis_cum3)))
    pge_vis_4 <- 1 / (1 + exp(-a_vis * (pd - vis_cum4)))
    pge_vis_5 <- 1 / (1 + exp(-a_vis * (pd - vis_cum5)))
    pge_vis_6 <- 1 / (1 + exp(-a_vis * (pd - vis_cum6)))
    visual <- pge_vis_1 + pge_vis_2 + pge_vis_3 + pge_vis_4 + pge_vis_5 + pge_vis_6

    # Item 7: Mental (0-4)
    men_cum1  <- b_men_1
    men_cum2  <- men_cum1 + b_men_2
    men_cum3  <- men_cum2 + b_men_3
    men_cum4  <- men_cum3 + b_men_4
    pge_men_1 <- 1 / (1 + exp(-a_men * (pd - men_cum1)))
    pge_men_2 <- 1 / (1 + exp(-a_men * (pd - men_cum2)))
    pge_men_3 <- 1 / (1 + exp(-a_men * (pd - men_cum3)))
    pge_men_4 <- 1 / (1 + exp(-a_men * (pd - men_cum4)))
    mental <- pge_men_1 + pge_men_2 + pge_men_3 + pge_men_4

    # Item 8: Ambulation (0-9)
    amb_cum1  <- b_amb_1
    amb_cum2  <- amb_cum1 + b_amb_2
    amb_cum3  <- amb_cum2 + b_amb_3
    amb_cum4  <- amb_cum3 + b_amb_4
    amb_cum5  <- amb_cum4 + b_amb_5
    amb_cum6  <- amb_cum5 + b_amb_6
    amb_cum7  <- amb_cum6 + b_amb_7
    amb_cum8  <- amb_cum7 + b_amb_8
    amb_cum9  <- amb_cum8 + b_amb_9
    pge_amb_1 <- 1 / (1 + exp(-a_amb * (pd - amb_cum1)))
    pge_amb_2 <- 1 / (1 + exp(-a_amb * (pd - amb_cum2)))
    pge_amb_3 <- 1 / (1 + exp(-a_amb * (pd - amb_cum3)))
    pge_amb_4 <- 1 / (1 + exp(-a_amb * (pd - amb_cum4)))
    pge_amb_5 <- 1 / (1 + exp(-a_amb * (pd - amb_cum5)))
    pge_amb_6 <- 1 / (1 + exp(-a_amb * (pd - amb_cum6)))
    pge_amb_7 <- 1 / (1 + exp(-a_amb * (pd - amb_cum7)))
    pge_amb_8 <- 1 / (1 + exp(-a_amb * (pd - amb_cum8)))
    pge_amb_9 <- 1 / (1 + exp(-a_amb * (pd - amb_cum9)))
    ambulation <- pge_amb_1 + pge_amb_2 + pge_amb_3 + pge_amb_4 + pge_amb_5 + pge_amb_6 + pge_amb_7 + pge_amb_8 + pge_amb_9

    # 4. FREM continuous-covariate observations. Each is the typical-population mean shifted by
    # the corresponding latent random-effect dimension; the residual error is essentially zero
    # (sigma fixed at 1e-5 in the source). These outputs let the FREM-style covariate
    # observations be simulated alongside the IRT trajectory.
    age_pred  <- age_mean  + etap3
    msd_pred  <- msd_mean  + etap4
    exnb_pred <- exnb_mean + etap5

    age_pred  ~ add(addSd_age_pred)
    msd_pred  ~ add(addSd_msd_pred)
    exnb_pred ~ add(addSd_exnb_pred)
  })
}
