Girard_2012_pimasertib <- function() {
  description <- "Joint K-PD / cumulative-logit Markov / Weibull-TTE-dropout model for ocular adverse events and treatment discontinuation in advanced solid-tumour and hematological-malignancy patients dosed with the MEK inhibitor pimasertib in two phase I dose-escalation studies (Girard 2012; DDMODEL00000215)"
  reference <- paste(
    "Girard P, Brockhaus B, Massimini G, Asiatiani E, Rejeb N, Rajeswaran RA,",
    "Lupfert C, von Richter O, Munafo A. (2012).",
    "Simultaneous ocular adverse event and treatment discontinuation model of pimasertib.",
    "PAGE 21 (2012) Abstr 2458 [www.page-meeting.org/?abstract=2458].",
    "DDMORE Foundation Model Repository: DDMODEL00000215.",
    sep = " "
  )
  vignette <- "Girard_2012_pimasertib"
  paper_specific_etas <- c("etalogit")

  units <- list(
    time          = "week",
    dosing        = "(K-PD AUC; AMT column carries the per-week pimasertib AUC, ng*h/mL or paper unit)",
    concentration = "(K-PD exposure proxy; same AUC units as the AMT dose)"
  )
  ddmore_id    <- "DDMODEL00000215"
  replicate_of <- NULL

  covariateData <- list(
    HYPERT = list(
      description        = "History of hypertension comorbidity (binary; 1 = prior or current hypertension, 0 = none)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no hypertension)",
      notes              = "Source column MHHY (medical history of hypertension). Time-fixed at study entry. Acts as additive shift on logit(P(AE >= k)) for both cumulative-logit thresholds.",
      source_name        = "MHHY"
    ),
    REGI_BID = list(
      description        = "Twice-daily-vs-once-daily dosing-regimen indicator (binary; 1 = BID, 0 = QD)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (QD or other non-BID schedule)",
      notes              = "Source column BID. Time-fixed per subject. Acts as additive shift on logit(P(AE >= k)) for both cumulative-logit thresholds. Distinct from the more-granular source REGI column (numeric regimen code) which Girard 2012 reduces to the binary BID indicator before testing it as a covariate.",
      source_name        = "BID"
    ),
    CMAX_M1 = list(
      description        = "Empirical-Bayes maximum plasma pimasertib concentration during the first month of dosing",
      units              = "ng/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source column CMAXM1. Per-subject summary derived from the upstream Girard 2012 pimasertib popPK model (not in nlmixr2lib at extraction time). Acts as additive shift on logit(P(AE >= k)) centred at 0 ng/mL: source MED17 = 0 sets the reference, so the linear term is theta * (CMAX_M1 - 0) = theta * CMAX_M1 directly.",
      source_name        = "CMAXM1"
    ),
    DOSE = list(
      description        = "Per-subject assigned daily dose of pimasertib",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source column DOSE. Time-fixed per subject (dose-escalation study; one dose level per cohort). Enters the dropout-hazard Weibull as a multiplicative exponential effect: hazard(t) ~ lambda * alpha * t^(alpha-1) * exp(beta * DOSE). Distinct from the K-PD AMT column, which carries the per-week pimasertib AUC and dosed into the central compartment.",
      source_name        = "DOSE"
    ),
    PREV_AE_SCORE = list(
      description        = "Previous-observation ocular-AE CTCAE grade (Markov-state covariate; integer 0-3)",
      units              = "(CTCAE ocular-toxicity grade, 0..3)",
      type               = "count",
      reference_category = NULL,
      notes              = "Source column PREVSCOR -- the lagged-DV column carried via NONMEM's IF (TIME.EQ.0) PREVSCOR=0 / PREVSCOR=DV idiom in $PK / $ERROR. Set to 0 at the first observation per subject; updated each subsequent observation to the previous-step sampled CTCAE grade. The cumulative-logit thresholds (b01/b11/b21 for AE >= 1, b02/b12/b22 for AE >= 2) and the Emax term are conditioned on the FPS group: FPS0 = (PREV == 0), FPS1 = (PREV in {1,2}), FPS2 = (PREV >= 3). Values 1 and 2 collapse into the FPS1 stratum per the source Markov grouping; in the simulated dataset shipped with DDMODEL00000215 the observed DVID == 2 values are 0/1/2/3.",
      source_name        = "PREVSCOR"
    )
  )

  population <- list(
    n_subjects     = 199,
    n_studies      = 2,
    age_range      = "(adult oncology cohort; specific range not extracted -- the linked PAGE 21 (2012) Abstr 2458 publication is conference-abstract-only and was not on disk for cross-check at extraction time)",
    weight_range   = "(not extracted)",
    sex_female_pct = "(not extracted)",
    disease_state  = "Advanced solid tumours and hematological malignancies (two phase I dose-escalation studies)",
    dose_range     = "1-255 mg/day pimasertib (orally; QD or BID schedules pooled across the two phase I studies; observed daily-dose values in the bundled simulated dataset: 1, 1.5, 2, 2.5, 3.5, 5, 7, 14, 16, 28, 30, 45, 46, 60, 68, 84, 90, 94, 120, 150, 195, 255 mg)",
    notes          = "n_subjects = 199 and n_observations = 3655 (DVID == 2 ocular-AE-grade rows after IGNORE(DVID.EQ.3)) read from Output_real_Pimasertib_AeDropout.lst run header; the listing reports ESTIMATION-EVALUATION (MAXEVALS=0) on the original-data fit, so the THETA / OMEGA values in the .lst FINAL PARAMETER ESTIMATE block equal the .mod $THETA / $OMEGA initials and are the publication's reported final estimates. Demographic detail (age, weight, sex split, race) is not derivable from the DDMORE bundle; the linked PAGE 21 (2012) abstract (URL www.page-meeting.org/?abstract=2458) is a conference abstract not available as a downloadable PDF and was not on disk for cross-check."
  )

  ini({
    # All values are taken from Output_real_Pimasertib_AeDropout.lst FINAL PARAMETER
    # ESTIMATE block (THETA TH 1..18; OMEGA ETA1 / ETA2). The .lst reports an
    # ESTIMATION-EVALUATION run (MAXEVALS=0): NONMEM evaluates the objective at
    # the supplied THETA / OMEGA without re-estimating, so the .mod $THETA / $OMEGA
    # blocks carry the publication's final estimates as their "initial" values
    # and the .lst FINAL block reproduces them. See ddmore-source.md.

    # Cumulative-logit thresholds for P(AE >= 1) (Markov-conditioned on the
    # previous score grouping FPS0/FPS1/FPS2 = (PREV == 0) / (PREV in {1,2}) /
    # (PREV >= 3)).
    b01 <- -6.12
    label("Logit intercept for P(ocular AE >= grade 1) given PREV_AE_SCORE = 0")
    # .lst TH 1 = -6.12E+00 (B01).

    b11 <- 1.67
    label("Increment to logit intercept for P(AE >= 1) given PREV_AE_SCORE in {1,2}")
    # .lst TH 2 = 1.67E+00 (B11).

    b21 <- 1.59
    label("Increment to logit intercept for P(AE >= 1) given PREV_AE_SCORE >= 3")
    # .lst TH 3 = 1.59E+00 (B21).

    # Cumulative-logit increment from logit(P(AE >= 1)) to logit(P(AE >= 2))
    # (also Markov-conditioned). Negative values keep PC2 < PC1 (monotone CDF).
    b02 <- -3.20
    label("Logit gap from P(AE >= 1) to P(AE >= 2) given PREV_AE_SCORE = 0")
    # .lst TH 4 = -3.20E+00 (B02).

    b12 <- -7.65
    label("Logit gap from P(AE >= 1) to P(AE >= 2) given PREV_AE_SCORE in {1,2}")
    # .lst TH 5 = -7.65E+00 (B12).

    b22 <- -0.214
    label("Logit gap from P(AE >= 1) to P(AE >= 2) given PREV_AE_SCORE >= 3")
    # .lst TH 6 = -2.14E-01 (B22).

    lkel <- log(2.33)
    label("Log of K-PD elimination rate constant for the exposure proxy (1/week)")
    # .lst TH 7 = 2.33E+00 (TVK).

    emax0 <- 4.04
    label("Maximum sigmoidal Emax effect on AE-score logit given PREV_AE_SCORE = 0")
    # .lst TH 8 = 4.04E+00 (EMAX0).

    emax1 <- -0.483
    label("Maximum sigmoidal Emax effect on AE-score logit given PREV_AE_SCORE >= 1")
    # .lst TH 9 = -4.83E-01 (EMAX1).

    # emax2 -- .lst TH10 = 0 FIX (EMAX2) -- is intentionally omitted from
    # ini() because nlmixr2 rejects ini-only parameters that are never read
    # by model(). The corresponding IF(PREVSCOR.GE.3) EMAX = EMAX2 branch
    # in the source $PK is commented out in the executable, so emax2 is
    # never actually consumed (the live emax-selection is emax0 vs emax1
    # only). The fixed-at-0 declaration is preserved here as a comment so
    # the .lst -> model parameter mapping stays exhaustive.

    led50 <- 7.69
    label("Log of ED50 for the sigmoidal Emax effect (units of K-PD exposure proxy)")
    # .lst TH11 = 7.69E+00 (LNED50). Source: ED50 = EXP(THETA(11)) * (CL_IND / 39.4)^THETA(12),
    # but THETA(12) is fixed at 0 (.lst TH12 = 0 FIX), so the CL_IND term collapses to 1
    # and ED50 reduces to exp(led50).

    e_hypert_logit <- 0.539
    label("Additive shift on the AE-score cumulative logit when HYPERT = 1 (logit units)")
    # .lst TH13 = 5.39E-01 (Cov11_MHHY). Source: ALL2 += COV11 * TCOV11 with COV11 = MHHY.

    e_regi_bid_logit <- -0.399
    label("Additive shift on the AE-score cumulative logit when REGI_BID = 1 (logit units)")
    # .lst TH14 = -3.99E-01 (Cov13_BID). Source: ALL2 += COV13 * TCOV13 with COV13 = BID.

    e_cmax_m1_logit <- 0.000902
    label("Additive shift on the AE-score cumulative logit per ng/mL of CMAX_M1 (logit units / (ng/mL))")
    # .lst TH15 = 9.02E-04 (COV17_CMAXM1). Source: ALL1 = (COV17 - MED17) * TCOV17 with COV17 = CMAXM1, MED17 = 0.

    llambda <- -3.32
    label("Log of Weibull baseline-hazard scale parameter (1/week^alpha at TIME = 1)")
    # .lst TH16 = -3.32E+00 (LNLAMBDA). Source: LAMBDA = EXP(THETA(16)).

    lalpha <- 0.232
    label("Log of Weibull baseline-hazard shape parameter (unitless)")
    # .lst TH17 = 2.32E-01 (LNALPHA). Source: ALPHA = EXP(THETA(17)).

    e_dose_haz <- 0.00416
    label("Linear coefficient on DOSE in the dropout hazard log-rate (1/mg)")
    # .lst TH18 = 4.16E-03 (BETA1). Source: HAZARD = LAMBDA * ALPHA * TIME^(ALPHA-1) * EXP(BETA1 * DOSE).

    # IIV: a single shared eta on the cumulative-logit shift (NONMEM IIV =
    # ETA(1) added to AA1 and AA2 alike). ETA(2) on K is fixed at 0 in the
    # source ($OMEGA(2) = 0 FIX), so no IIV on the K-PD elimination rate.
    etalogit ~ 0.786
    # .lst OMEGA ETA1 diag = 7.86E-01.
  })

  model({
    # Time scale: source uses TIME = WEEK (NONMEM column rename TIME=WEEK in
    # $INPUT). units$time = "week", so `time` here is weeks directly.

    # 1. Markov-state grouping: collapse PREV_AE_SCORE to the three FPS strata
    #    used in the source $PK block (FPS0 = no prior AE, FPS1 = prior grade
    #    1 or 2, FPS2 = prior grade >= 3).
    fps0 <- (PREV_AE_SCORE < 0.5)
    fps1 <- (PREV_AE_SCORE >= 0.5) * (PREV_AE_SCORE < 2.5)
    fps2 <- (PREV_AE_SCORE >= 2.5)

    # 2. Per-stratum cumulative-logit intercepts for AE >= 1 and AE >= 2.
    a1 <- fps0 * b01 + fps1 * b11 + fps2 * b21
    a2 <- a1 + (fps0 * b02 + fps1 * b12 + fps2 * b22)

    # 3. Per-stratum Emax (the source comments out the IF(PREVSCOR.GE.3) branch
    #    that would have selected emax2; we follow the executable convention
    #    and pick emax0 for FPS0 and emax1 for FPS1+FPS2).
    emax_iv <- fps0 * emax0 + (1 - fps0) * emax1

    # 4. K-PD exposure compartment: dosed with AMT = AUC, decays at rate kel.
    kel <- exp(lkel)
    d/dt(central) <- -kel * central
    expo <- central

    # 5. Sigmoidal Emax effect on the cumulative logit.
    ed50 <- exp(led50)
    eff  <- emax_iv * expo / (expo + ed50)

    # 6. Continuous-covariate logit shift (centred at 0 ng/mL for CMAX_M1 per
    #    source MED17 = 0).
    cov_logit <- e_hypert_logit * HYPERT + e_regi_bid_logit * REGI_BID + e_cmax_m1_logit * CMAX_M1

    # 7. Cumulative-logit linear predictors (with shared eta on both thresholds).
    aa1 <- a1 + etalogit + eff + cov_logit
    aa2 <- a2 + etalogit + eff + cov_logit

    # 8. Cumulative and per-grade probabilities (CTCAE 0 / 1-2 / >=3 strata
    #    on the cumulative-logit scale; in the source, DV in {1, 2} both map
    #    to P1 and DV >= 3 maps to P2).
    pc0 <- 1
    pc1 <- expit(aa1)
    pc2 <- expit(aa2)
    p2  <- pc2
    p1  <- pc1 - pc2
    p0  <- pc0 - pc1

    # 9. Weibull dropout sub-model (state cumhaz integrates the instantaneous
    #    hazard; the source DVID == 4 / DVID == 5 likelihood expressions
    #    Y = exp(-cumhaz)*hazard / Y = exp(-cumhaz) consume cumhaz).
    lambda <- exp(llambda)
    alpha  <- exp(lalpha)
    hazard <- lambda * alpha * (time + 1e-12)^(alpha - 1) * exp(e_dose_haz * DOSE)
    d/dt(cumhaz) <- hazard
    survival <- exp(-cumhaz)

    # 10. Observation. The source NONMEM model uses LIKE LAPLACE with a
    #     conditional categorical-likelihood Y per DVID (DVID == 2 returns
    #     P0 / P1 / P2; DVID == 4 returns dropout-hazard density;
    #     DVID == 5 returns survival probability). nlmixr2 / rxode2 do not
    #     natively express this multi-DVID joint likelihood, so the formal
    #     observation here is the typical-value expected ordinal AE score
    #     (0 * P0 + 1 * P1 + 2 * P2 in {0..2}) modelled as Poisson -- purely
    #     a placeholder that lets the model parse and that drives F.3
    #     mechanistic-sanity simulation in the validation vignette. The full
    #     Markov / cumulative-logit / Weibull-TTE structure is encoded
    #     above and is exercised by the vignette via direct rxSolve()
    #     evaluation of pc1 / pc2 / p0 / p1 / p2 / hazard / survival; the
    #     placeholder Poisson is NOT a valid likelihood for fitting against
    #     the original Girard 2012 dataset and the deviation is documented
    #     in the vignette's Assumptions and deviations section.
    expected_aescore <- p1 + 2 * p2
    aescore <- expected_aescore
    aescore ~ pois(expected_aescore)
  })
}
