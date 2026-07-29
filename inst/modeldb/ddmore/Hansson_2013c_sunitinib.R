Hansson_2013c_sunitinib <- function() {
  description <- "Population PD model of fatigue (NCI-CTC grades 0 / 1 / 2 / 3+) in adults with imatinib-resistant gastrointestinal stromal tumours (GIST) on sunitinib. A first-order Markov + proportional-odds (PO) likelihood describes the fatigue-grade transition probabilities at each scheduled visit, conditional on the previous fatigue grade. The cumulative-logit baselines are shifted per starting state by a placebo coefficient on the relative change in plasma soluble VEGFR-3 (sVEGFR-3) from baseline; sVEGFR-3 itself follows an indirect-response turnover driven by the per-cycle drug-exposure summary AUC = DOSE / CLI. The PD model has no PK ODE and consumes individual posthoc upstream-PD parameters (BAS_SVEGFR3, MRT_SVEGFR3, EC50_SVEGFR3) and posthoc upstream-PK clearance (CLI) as data covariates. Random effects are diagonal across the four per-state baseline logits."
  reference <- paste(
    "Hansson EK, Ma G, Amantea MA, French J, Milligan PA, Friberg LE,",
    "Karlsson MO.",
    "PKPD modeling of predictors for adverse effects and overall survival",
    "in sunitinib-treated patients with GIST.",
    "CPT Pharmacometrics Syst Pharmacol. 2013;2(11):e85.",
    "doi:10.1038/psp.2013.62.",
    "DDMORE Foundation Model Repository: DDMODEL00000222.",
    "Upstream sVEGFR-3 biomarker dynamics adapted from",
    "Hansson 2013 (CPT Pharmacometrics Syst Pharmacol 2013;2:e84,",
    "doi:10.1038/psp.2013.61, DDMODEL00000197); see",
    "modellib('Hansson_2013a_sunitinib').",
    sep = " "
  )
  vignette <- "Hansson_2013c_sunitinib"
  units <- list(
    time = "hour",
    dosing = "mg",
    concentration = "(NCI-CTC fatigue grade 0-3+, ordinal)"
  )
  ddmore_id <- "DDMODEL00000222"
  replicate_of <- NULL

  covariateData <- list(
    DOSE = list(
      description        = "Current administered sunitinib daily dose (mg) carried as a time-varying data column. Set to 0 during off-cycles (4 weeks on / 2 weeks off in the Hansson 2013 GIST cohort) or for placebo subjects so the derived AUC = DOSE / CLI becomes 0.",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "DDMORE-bundle simulated dataset reports DOSE = 50 (treated, on-cycle) or 0 (off-cycle / placebo) at every record. The .mod feeds DOSE into AUC = DOSE / CLI in $PK at every event call, producing a per-cycle daily-AUC equivalent (mg*h/L). For typical-cohort vignette simulations the value is held at 50 mg during the 4-week on-cycles.",
      source_name        = "DOSE"
    ),
    CLI = list(
      description        = "Individual posthoc total plasma clearance (L/h) of sunitinib from the paper's upstream 2-compartment popPK fit. Per-subject, time-fixed.",
      units              = "L/h",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Required input. The DDMORE-bundle simulated dataset carries CLI values 30-43 L/h across its three subjects; this is consistent with the Houk et al. 2010 typical sunitinib CL and with the typical-value reference (32.819 L/h) used by the upstream Hansson 2013a biomarker PD model. For a re-fit or new-population simulation the user must supply individual CL drawn from a sunitinib popPK model (the Hansson 2013 paper text describes the upstream PK as a 'previously developed 2-compartment model'; that popPK is not extracted into nlmixr2lib).",
      source_name        = "CL"
    ),
    BAS_SVEGFR3 = list(
      description        = "Individual posthoc baseline sVEGFR-3 (pg/mL) from the upstream Hansson 2013a biomarker indirect-response PD fit (DDMODEL00000197). Per-subject, time-fixed; used both as the initial condition for the in-model svegfr3 state and as the denominator in the relative-change driver bm = (svegfr3 - BAS_SVEGFR3) / BAS_SVEGFR3.",
      units              = "pg/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Required input. The DDMORE-bundle simulated dataset carries BAS_SVEGFR3 values 42554-57365 pg/mL across its three subjects, consistent with the Hansson 2013a typical-value sVEGFR-3 baseline of 63900 pg/mL (the 011 vignette's typical-value figure shows steady-state at the BAS_SVEGFR3 set point, with on-cycle drug effect depleting it by ~50%). For new-population simulations either (a) simulate from `Hansson_2013a_sunitinib` to obtain individual posthoc baselines, or (b) set every subject to the typical 63900 pg/mL.",
      source_name        = "BAS3"
    ),
    MRT_SVEGFR3 = list(
      description        = "Individual posthoc mean residence time of sVEGFR-3 (h) from the upstream Hansson 2013a biomarker indirect-response PD fit (DDMODEL00000197). Per-subject, time-fixed; used as kout3 = 1 / MRT_SVEGFR3 inside model().",
      units              = "h",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Required input. The DDMORE-bundle simulated dataset carries MRT_SVEGFR3 values 313-408 h across its three subjects, consistent with the Hansson 2013a typical sVEGFR-3 MRT of 401 h. Same population strategy as BAS_SVEGFR3.",
      source_name        = "MRT3"
    ),
    EC50_SVEGFR3 = list(
      description        = "Individual posthoc EC50 of the simple-Imax drug effect on sVEGFR-3 (mg*h/L AUC) from the upstream Hansson 2013a biomarker indirect-response PD fit (DDMODEL00000197). Per-subject, time-fixed; appears in the drug-effect term eff3 = auc / (EC50_SVEGFR3 + auc).",
      units              = "mg*h/L (matches the auc = DOSE / CLI exposure summary)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Required input. The DDMORE-bundle simulated dataset carries EC50_SVEGFR3 values 1.0-2.8 mg*h/L across its three subjects, consistent with the Hansson 2013a typical sVEGFR-3 IC50 of 1.0 mg*h/L. Same population strategy as BAS_SVEGFR3.",
      source_name        = "EC53"
    )
  )

  population <- list(
    n_subjects     = 303L,
    n_studies      = 1L,
    age_range      = "adults with imatinib-resistant GIST (paper not on disk; Hansson 2013c CPT Pharmacometrics Syst Pharmacol 2013;2:e85 baseline-demographics table not available in the bundle)",
    weight_range   = "not reported in the DDMORE bundle",
    sex_female_pct = NA_real_,
    race_ethnicity = NULL,
    disease_state  = "Imatinib-resistant gastrointestinal stromal tumours (GIST). The Hansson 2013 multinational Phase III sunitinib trial (sunitinib 50 mg PO QD on a 4-weeks-on / 2-weeks-off schedule + placebo run-in) supplied the biomarker / AE / overall-survival dataset for the e84 / e85 / e86 paper trio; this fatigue Markov + proportional-odds model fits to the per-visit fatigue-grade data from that cohort.",
    dose_range     = "Sunitinib 50 mg PO QD on a 4-weeks-on / 2-weeks-off schedule (standard GIST regimen at the time of the source study). Placebo arm: no sunitinib.",
    regions        = "Phase III multinational trial; specific regions not reported in the DDMORE bundle.",
    biomarkers     = "Fatigue grade per NCI-CTC v3 (ordinal: 0 = none, 1 = mild, 2 = moderate, 3+ = severe or worse, with the .mod's PDV.GT.2 branch pooling grade 3 and grade 4). The upstream sVEGFR-3 biomarker dynamics are consumed as data covariates (BAS_SVEGFR3 + MRT_SVEGFR3 + EC50_SVEGFR3) rather than re-fitted.",
    notes          = "n_subjects = 303 carried over from the upstream Hansson 2013a biomarker model (same Phase III trial cohort); the AE-OS paper's fatigue analysis-set count is not derivable from the DDMORE bundle. Detailed baseline-demographics (age, weight, sex, race, prior-imatinib-duration distributions) are in the linked publication (CPT Pharmacometrics Syst Pharmacol 2013;2:e85, doi:10.1038/psp.2013.62), which was not on disk in /home/bill/github/mab_human_consensus/literature at extraction time. The bundle's simulated dataset is intentionally minimal (three subjects, daily observations over 27 days) and is not representative of the published cohort."
  )

  ini({
    # ----------------------------------------------------------------------
    # Final estimates from Output_real_Fatigue_GIST.lst_Fatigue_PSP_2014
    # FINAL PARAMETER ESTIMATE block (post-MINIMIZATION SUCCESSFUL).
    #
    # The bundle .mod's $THETA / $OMEGA are NOT the published values for the
    # parameters that are estimated (TH1-TH12 transition betas, TH13-TH16
    # placebo coefficients, OMEGA1-OMEGA4) -- the .mod is a MAXEVAL=0
    # evaluation-run scaffold whose initial values were left at rough
    # placeholders (e.g. TH4 = -2.78, OMEGA(i,i) = 0.01). The actual
    # published values come from the .lst's FINAL PARAMETER ESTIMATE block.
    #
    # Mapping between the .lst's run (which had 17 thetas: TH13 = KEO FIXED
    # 0.422, an unused effect-compartment relic; TH14-TH17 = PX0..PX3
    # placebo coefficients) and the bundle .mod's 16 thetas (TH13-TH16 =
    # PX0..PX3 placebo coefficients) is .mod TH(13+i) <-> .lst TH(14+i) for
    # i in 0..3. The structural sVEGFR-3 ODE is identical between the two
    # (the .lst's KEO-driven effect compartment was never fed into the EFF3
    # expression -- the DRG0/DRG1/DRG2/DRG3 lines that consumed A(2) are
    # commented out in the .mod header -- so EFF3 = AUC/(EC53+AUC) in both
    # cases with AUC = DOSE/CL; the .lst's CP rename is cosmetic).
    #
    # Reparameterization note: the .mod parameterizes the proportional-odds
    # cumulative logits as B1 + B2 + B3 with B2 <= 0 and B3 <= 0 (so that
    # P(grade>=1) >= P(grade>=2) >= P(grade>=3+) by construction). nlmixr2's
    # mu-reference parser does not allow more than one population parameter
    # in the additive expression that carries a single eta, so this model
    # uses the equivalent **cumulative-logit** parameterization: per state,
    # three separate baseline-logit population parameters
    # `clge1_px<i>` / `clge2_px<i>` / `clge3_px<i>` carrying the cumulative
    # sums TH(b1), TH(b1)+TH(b2), TH(b1)+TH(b2)+TH(b3). The math is
    # identical; the cost is that the B2 <= 0 / B3 <= 0 monotonicity
    # constraint is not enforced at the parameter level (it is enforced at
    # the typical-value level by initialising the cumulative sums in
    # decreasing order, and a user re-fitting the model should retain the
    # corresponding ordered-logit constraint via fit-time bounds).
    # ----------------------------------------------------------------------

    # --- PX0 cumulative baseline logits (drug-free, IIV-free) ---
    # TH1 = -5.85; TH2 = -1.14; TH3 = -1.60.
    clge1_px0 <- -5.85;  label("PX0 typical cumulative logit for grade>=1 given previous = 0 (= TH1)")               # .lst TH 1
    clge2_px0 <- -6.99;  label("PX0 typical cumulative logit for grade>=2 (= TH1+TH2 = -5.85 + -1.14)")              # .lst TH 1 + TH 2
    clge3_px0 <- -8.59;  label("PX0 typical cumulative logit for grade>=3+ (= TH1+TH2+TH3 = -5.85 + -1.14 + -1.60)") # .lst TH 1 + TH 2 + TH 3

    # --- PX1 cumulative baseline logits ---
    # TH4 = 2.63; TH5 = -10.7; TH6 = -1.77.
    clge1_px1 <-  2.63;   label("PX1 typical cumulative logit for grade>=1 given previous = 1 (= TH4)")                # .lst TH 4
    clge2_px1 <- -8.07;   label("PX1 typical cumulative logit for grade>=2 (= TH4+TH5 = 2.63 + -10.7)")                # .lst TH 4 + TH 5
    clge3_px1 <- -9.84;   label("PX1 typical cumulative logit for grade>=3+ (= TH4+TH5+TH6 = 2.63 + -10.7 + -1.77)")   # .lst TH 4 + TH 5 + TH 6

    # --- PX2 cumulative baseline logits ---
    # TH7 = 2.86; TH8 = -0.427; TH9 = -11.6.
    clge1_px2 <-  2.86;   label("PX2 typical cumulative logit for grade>=1 given previous = 2 (= TH7)")                 # .lst TH 7
    clge2_px2 <-  2.433;  label("PX2 typical cumulative logit for grade>=2 (= TH7+TH8 = 2.86 + -0.427)")                # .lst TH 7 + TH 8
    clge3_px2 <- -9.167;  label("PX2 typical cumulative logit for grade>=3+ (= TH7+TH8+TH9 = 2.86 + -0.427 + -11.6)")    # .lst TH 7 + TH 8 + TH 9

    # --- PX3 cumulative baseline logits ---
    # TH10 = 3.06; TH11 = -0.0903; TH12 = -0.636.
    clge1_px3 <-  3.06;     label("PX3 typical cumulative logit for grade>=1 given previous = 3+ (= TH10)")              # .lst TH10
    clge2_px3 <-  2.9697;   label("PX3 typical cumulative logit for grade>=2 (= TH10+TH11 = 3.06 + -0.0903)")            # .lst TH10 + TH11
    clge3_px3 <-  2.3337;   label("PX3 typical cumulative logit for grade>=3+ (= TH10+TH11+TH12 = 3.06 + -0.0903 + -0.636)") # .lst TH10 + TH11 + TH12

    # --- Placebo (sVEGFR-3 relative change) coefficients on the four B1 logits ---
    # The .mod reads VE3i = THETA(13+i) * BM where BM = (svegfr3 - BAS_SVEGFR3)/BAS_SVEGFR3.
    # Each VE3i shifts the corresponding starting-state baseline logit additively;
    # negative coefficients with negative BM (drug depletes sVEGFR-3) push the
    # baseline logit positive => higher fatigue probability under drug.
    e_bm_px0 <- -1.93;  label("Slope of sVEGFR-3 relative-change BM on PX0 baseline logit (log-odds units per unit BM)")  # .lst TH14 -> bundle .mod TH13
    e_bm_px1 <- -4.62;  label("Slope of sVEGFR-3 relative-change BM on PX1 baseline logit")                              # .lst TH15 -> bundle .mod TH14
    e_bm_px2 <- -4.64;  label("Slope of sVEGFR-3 relative-change BM on PX2 baseline logit")                              # .lst TH16 -> bundle .mod TH15
    e_bm_px3 <- -3.32;  label("Slope of sVEGFR-3 relative-change BM on PX3 baseline logit")                              # .lst TH17 -> bundle .mod TH16

    # --- IIV on the per-state baseline logits (diagonal $OMEGA in the .mod) ---
    # The .mod adds ETA(i) only to the B1 (grade>=1) logit per state; we attach
    # the same eta to all three cumulative logits per state inside model() so the
    # subject-level shift is shared across the three thresholds (matching the
    # .mod's per-state shared-eta structure).
    etaclge_px0 ~ 1.12;   label("IIV variance on PX0 cumulative-logit baselines (logit-scale variance, shared across grade>=1/2/3+ thresholds)")  # .lst OMEGA(1,1)
    etaclge_px1 ~ 1.57;   label("IIV variance on PX1 cumulative-logit baselines")                                                                  # .lst OMEGA(2,2)
    etaclge_px2 ~ 1.68;   label("IIV variance on PX2 cumulative-logit baselines")                                                                  # .lst OMEGA(3,3)
    etaclge_px3 ~ 0.708;  label("IIV variance on PX3 cumulative-logit baselines")                                                                  # .lst OMEGA(4,4)

    # --- Residual-error placeholder ---
    # The published likelihood is the joint Markov + proportional-odds form
    # the .mod implements via its conditional Y = P_xy assignment under
    # $ESTIMATION ... LAPLACE LIKE. nlmixr2 / rxode2 do not natively express
    # "previous DV"-conditioned discrete likelihoods, so this model file
    # provides the typical-value transition-probability outputs as
    # deterministic functions of time; for nlmixr2 fitting compatibility the
    # observation is taken to be the "expected fatigue grade given previous
    # state = 0" output (a continuous 0-3 quantity), with a placeholder
    # additive residual error the user can tune at fit time. See the
    # vignette's Assumptions and deviations section for the full deviation
    # rationale.
    addSd_fatigue_grade <- 0.5
    label("Placeholder additive residual error on the typical-value expected-fatigue output (units: ordinal grade); not from the source -- see vignette Assumptions")
  })

  model({
    # 1. Per-cycle drug-exposure summary (mg*h/L AUC).
    # mg / (L/h) = mg*h/L; matches the .mod's `AUC = DOSE / CL` line.
    auc <- DOSE / CLI

    # 2. Simple-Imax drug effect on sVEGFR-3 (Hill = 1, EMAX = 1; fixed in
    #    the .mod via EMAX = 1, no Hill exponentiation in the EFF3 line).
    eff3 <- auc / (EC50_SVEGFR3 + auc)

    # 3. sVEGFR-3 indirect-response turnover with Kin inhibition.
    kout3 <- 1 / MRT_SVEGFR3
    kin3  <- BAS_SVEGFR3 * kout3
    svegfr3(0)     <- BAS_SVEGFR3
    d/dt(svegfr3)  <- kin3 * (1 - eff3) - kout3 * svegfr3

    # 4. Biomarker driver: relative change in sVEGFR-3 from individual baseline.
    bm <- (svegfr3 - BAS_SVEGFR3) / BAS_SVEGFR3

    # 5. Cumulative logits per starting state. Each adds the per-state placebo
    #    shift (e_bm_px<i> * bm) and the per-state shared eta. The mu-ref
    #    expression for each cumulative logit is therefore (single typical-value
    #    pop param) + (covariate-effect param * state-derived bm) + (eta).
    pge10 <- expit(clge1_px0 + e_bm_px0 * bm + etaclge_px0)
    pge20 <- expit(clge2_px0 + e_bm_px0 * bm + etaclge_px0)
    pge30 <- expit(clge3_px0 + e_bm_px0 * bm + etaclge_px0)

    pge11 <- expit(clge1_px1 + e_bm_px1 * bm + etaclge_px1)
    pge21 <- expit(clge2_px1 + e_bm_px1 * bm + etaclge_px1)
    pge31 <- expit(clge3_px1 + e_bm_px1 * bm + etaclge_px1)

    pge12 <- expit(clge1_px2 + e_bm_px2 * bm + etaclge_px2)
    pge22 <- expit(clge2_px2 + e_bm_px2 * bm + etaclge_px2)
    pge32 <- expit(clge3_px2 + e_bm_px2 * bm + etaclge_px2)

    pge13 <- expit(clge1_px3 + e_bm_px3 * bm + etaclge_px3)
    pge23 <- expit(clge2_px3 + e_bm_px3 * bm + etaclge_px3)
    pge33 <- expit(clge3_px3 + e_bm_px3 * bm + etaclge_px3)

    # 6. Transition probabilities (each row = P(next = 0/1/2/3+ | previous = i)).
    # The four rows each sum to 1 by construction; monotonicity of the cumulative
    # logits (clge1 >= clge2 >= clge3 at the typical values) keeps the per-row
    # probabilities non-negative for the published parameter values.
    p00 <- 1     - pge10
    p10 <- pge10 - pge20
    p20 <- pge20 - pge30
    p30 <- pge30

    p01 <- 1     - pge11
    p11 <- pge11 - pge21
    p21 <- pge21 - pge31
    p31 <- pge31

    p02 <- 1     - pge12
    p12 <- pge12 - pge22
    p22 <- pge22 - pge32
    p32 <- pge32

    p03 <- 1     - pge13
    p13 <- pge13 - pge23
    p23 <- pge23 - pge33
    p33 <- pge33

    # 7. Expected fatigue grade given previous = 0 (the natural-history /
    #    placebo-arm reference output for typical-value visualization).
    expected_grade_from0 <- 0 * p00 + 1 * p10 + 2 * p20 + 3 * p30
    fatigue_grade        <- expected_grade_from0

    fatigue_grade ~ add(addSd_fatigue_grade)
  })
}
