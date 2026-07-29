Hansson_2013_sunitinib_hfs <- function() {
  description <- "Population PD model of hand-foot syndrome (HFS; NCI-CTC grades 0 / 1 / 2 / 3+) in adults with imatinib-resistant gastrointestinal stromal tumours (GIST) on sunitinib. A first-order Markov + proportional-odds (PO) likelihood describes the per-visit HFS-grade transition probabilities conditional on the previous grade. The cumulative-logit baselines per starting state are shifted by a starting-state-specific slope on a delayed sVEGFR-3 relative-change signal (effect-compartment-smoothed by ke0 = 0.347/h per Hansson 2013 e85). sVEGFR-3 itself is simulated in-model as a one-compartment indirect-response turnover driven by the per-cycle drug-exposure summary AUC = DOSE / CLI. The PD model has no PK ODE and consumes individual posthoc upstream-PD parameters (BAS_SVEGFR3, MRT_SVEGFR3, EC50_SVEGFR3) and posthoc upstream-PK clearance (CLI) as data covariates. Random effects are diagonal across the four per-state baseline logits; the >=3 state IIV is fixed to zero (NE in Hansson 2013 Table 3)."
  reference <- paste(
    "Hansson EK, Ma G, Amantea MA, French J, Milligan PA, Friberg LE,",
    "Karlsson MO.",
    "PKPD modeling of predictors for adverse effects and overall survival",
    "in sunitinib-treated patients with GIST.",
    "CPT Pharmacometrics Syst Pharmacol. 2013;2(11):e85.",
    "doi:10.1038/psp.2013.62.",
    "Sister model files from the same paper:",
    "modellib('Hansson_2013_sunitinib_myelosuppression'),",
    "modellib('Hansson_2013_sunitinib_dbp'),",
    "modellib('Hansson_2013c_sunitinib') [fatigue, structurally parallel to this HFS model],",
    "modellib('Hansson_2013_sunitinib_os').",
    "Upstream sVEGFR-3 biomarker dynamics adapted from",
    "Hansson 2013 (CPT Pharmacometrics Syst Pharmacol 2013;2:e84,",
    "doi:10.1038/psp.2013.61, DDMODEL00000197); see",
    "modellib('Hansson_2013a_sunitinib').",
    "Mixed-effects Markov + proportional-odds form adapted from",
    "Henin E et al. Clin Pharmacol Ther 2009;85(4):418-425,",
    "doi:10.1038/clpt.2008.220 and Zingmark PH, Kagedal M, Karlsson MO,",
    "J Pharmacokinet Pharmacodyn 2005;32(2):261-281,",
    "doi:10.1007/s10928-005-0034-2.",
    sep = " "
  )
  vignette <- "Hansson_2013_sunitinib_hfs"
  units <- list(
    time = "hour",
    dosing = "mg",
    concentration = "(NCI-CTC HFS grade 0-3+, ordinal)"
  )

  covariateData <- list(
    DOSE = list(
      description        = "Current administered sunitinib daily dose (mg) carried as a time-varying data column. Set to 0 during off-cycles or for placebo subjects so the derived AUC = DOSE / CLI becomes 0.",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Hansson 2013 e85 Methods describes sunitinib at 25-75 mg PO QD on 4/2, 2/2, 2/1, or continuous schedules. For typical-cohort vignette simulations the value is held at 50 mg during on-cycles of a 4/2 schedule.",
      source_name        = "DOSE"
    ),
    CLI = list(
      description        = "Individual posthoc total plasma clearance (L/h) of sunitinib from the paper's upstream 2-compartment popPK fit. Per-subject, time-fixed.",
      units              = "L/h",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Required input. The companion Hansson_2013a / Hansson_2013c sunitinib model files use a typical-value reference of 32.819 L/h.",
      source_name        = "CL"
    ),
    BAS_SVEGFR3 = list(
      description        = "Individual posthoc baseline sVEGFR-3 (pg/mL) from the upstream Hansson 2013a biomarker indirect-response PD fit (DDMODEL00000197). Per-subject, time-fixed; used as the initial condition for the in-model svegfr3 state and as the denominator in the relative-change driver bm_input = (svegfr3 - BAS_SVEGFR3) / BAS_SVEGFR3.",
      units              = "pg/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Required input. Hansson 2013 e84 reports a typical sVEGFR-3 baseline of 63900 pg/mL.",
      source_name        = "BAS3"
    ),
    MRT_SVEGFR3 = list(
      description        = "Individual posthoc mean residence time of sVEGFR-3 (h) from the upstream Hansson 2013a biomarker indirect-response PD fit (DDMODEL00000197). Per-subject, time-fixed; used as kout3 = 1 / MRT_SVEGFR3 inside model().",
      units              = "h",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Required input. The Hansson 2013a typical sVEGFR-3 MRT is 401 h.",
      source_name        = "MRT3"
    ),
    EC50_SVEGFR3 = list(
      description        = "Individual posthoc EC50 of the simple-Imax drug effect on sVEGFR-3 (mg*h/L AUC) from the upstream Hansson 2013a biomarker indirect-response PD fit (DDMODEL00000197). Per-subject, time-fixed; appears in the drug-effect term eff3 = auc / (EC50_SVEGFR3 + auc).",
      units              = "mg*h/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Required input. The Hansson 2013a typical IC50 is 1.0 mg*h/L.",
      source_name        = "EC53"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 303L,
    n_studies      = 4L,
    age_range      = "adults with imatinib-resistant GIST (paper text reports n=303 pooled across phases I-III; per-cohort baseline-demographics table not in the trimmed PDF section)",
    weight_range   = "not reported in the on-disk trimmed paper text",
    sex_female_pct = NA_real_,
    race_ethnicity = NULL,
    disease_state  = "imatinib-resistant gastrointestinal stromal tumours (GIST). Pooled four sunitinib studies (Demetri 2006, George 2009, Shirao 2010, Maki 2005).",
    dose_range     = "sunitinib 25-75 mg PO QD on a 4/2, 2/2, 2/1 or continuous schedule (Table 1).",
    regions        = "multinational (study 1004) and Japanese (study 1045).",
    biomarkers     = "HFS grade per NCI-CTC v3 (ordinal: 0 = none, 1, 2, 3+ = grouped with grade 4 due to rare grade 4 events per paper Methods). Observed grade distribution per study: study 1004 (0: 83%, 1: 5.0%, 2: 6.9%, 3: 5.4%); study 1047 (0: 100%); study 1045 (0: 14%, 1: 20%, 2: 34%, 3: 31%); study 013 (NA). The upstream sVEGFR-3 biomarker dynamics are consumed as data covariates (BAS_SVEGFR3 + MRT_SVEGFR3 + EC50_SVEGFR3) rather than re-fitted.",
    notes          = "n_subjects = 303 reported in Hansson 2013 e85 Methods. Per-cohort baseline demographics are not transcribed in the trimmed paper text. HFS grade 4 was reported in 0% of patients (Methods) so was grouped with grade 3 into a single category."
  )

  ini({
    # ----------------------------------------------------------------------
    # Final estimates from Hansson 2013 e85 Table 3 'HFS model' column.
    # The .mod parameterizes the proportional-odds cumulative logits as
    # B1 + B2 + B3 with B2 <= 0 and B3 <= 0 (so that P(grade>=1) >=
    # P(grade>=2) >= P(grade>=3+) by construction). nlmixr2's mu-reference
    # parser does not allow more than one population parameter in the
    # additive expression that carries a single eta, so this model uses
    # the equivalent **cumulative-logit** parameterization: per state,
    # three separate baseline-logit population parameters
    # `clge1_px<i>` / `clge2_px<i>` / `clge3_px<i>` carrying the cumulative
    # sums TH(b1), TH(b1)+TH(b2), TH(b1)+TH(b2)+TH(b3). This matches the
    # parameterization used by Hansson_2013c_sunitinib (fatigue) for the
    # same paper.
    # ----------------------------------------------------------------------

    # --- PX0 cumulative baseline logits (drug-free, IIV-free typical values) ---
    # Table 3 HFS: B1|0 = -10.4; B2|0 = -0.974; B>=3|0 = -1.59
    clge1_px0 <- -10.4;   label("PX0 typical cumulative logit for HFS grade>=1 given previous = 0 (= Table 3 HFS B1|0)")           # Table 3 HFS B1|0 = -10.4 (RSE 11%)
    clge2_px0 <- -11.374; label("PX0 typical cumulative logit for HFS grade>=2 (= B1|0 + B2|0 = -10.4 + -0.974)")                 # Table 3 HFS B1|0 + B2|0 = -10.4 + -0.974 (RSE 13%)
    clge3_px0 <- -12.964; label("PX0 typical cumulative logit for HFS grade>=3+ (= B1|0 + B2|0 + B>=3|0 = -10.4 + -0.974 + -1.59)") # Table 3 HFS B>=3|0 = -1.59 (RSE 19%)

    # --- PX1 cumulative baseline logits ---
    # Table 3 HFS: B1|1 = 2.29; B2|1 = -9.53; B>=3|1 = -1.33
    clge1_px1 <-  2.29;   label("PX1 typical cumulative logit for HFS grade>=1 given previous = 1 (= Table 3 HFS B1|1)")           # Table 3 HFS B1|1 = 2.29 (RSE 14%)
    clge2_px1 <- -7.24;   label("PX1 typical cumulative logit for HFS grade>=2 (= B1|1 + B2|1 = 2.29 + -9.53)")                    # Table 3 HFS B2|1 = -9.53 (RSE 5.0%)
    clge3_px1 <- -8.57;   label("PX1 typical cumulative logit for HFS grade>=3+ (= B1|1 + B2|1 + B>=3|1 = 2.29 + -9.53 + -1.33)")  # Table 3 HFS B>=3|1 = -1.33 (RSE 24%)

    # --- PX2 cumulative baseline logits ---
    # Table 3 HFS: B1|2 = 3.04; B2|2 = -0.747; B>=3|2 = -9.09
    clge1_px2 <-  3.04;   label("PX2 typical cumulative logit for HFS grade>=1 given previous = 2 (= Table 3 HFS B1|2)")            # Table 3 HFS B1|2 = 3.04 (RSE 15%)
    clge2_px2 <-  2.293;  label("PX2 typical cumulative logit for HFS grade>=2 (= B1|2 + B2|2 = 3.04 + -0.747)")                    # Table 3 HFS B2|2 = -0.747 (RSE 14%)
    clge3_px2 <- -6.797;  label("PX2 typical cumulative logit for HFS grade>=3+ (= B1|2 + B2|2 + B>=3|2 = 3.04 + -0.747 + -9.09)")  # Table 3 HFS B>=3|2 = -9.09 (RSE 5.1%)

    # --- PX3+ cumulative baseline logits ---
    # Table 3 HFS: B1|>3 = 3.4; B2|>3 = -1.65; B>=3|>=3 = -0.453
    clge1_px3 <-  3.4;    label("PX3+ typical cumulative logit for HFS grade>=1 given previous >= 3+ (= Table 3 HFS B1|>3)")        # Table 3 HFS B1|>3 = 3.4 (RSE 21%)
    clge2_px3 <-  1.75;   label("PX3+ typical cumulative logit for HFS grade>=2 (= B1|>3 + B2|>3 = 3.4 + -1.65)")                   # Table 3 HFS B2|>3 = -1.65 (RSE 23%)
    clge3_px3 <-  1.297;  label("PX3+ typical cumulative logit for HFS grade>=3+ (= B1|>3 + B2|>3 + B>=3|>=3 = 3.4 + -1.65 + -0.453)") # Table 3 HFS B>=3|>=3 = -0.453 (RSE 37%)

    # --- Placebo (sVEGFR-3 relative-change) slope coefficients on the four baseline logits ---
    # Table 3 HFS Slope x|<i> values. Negative coefficients with negative bm
    # (drug depletes sVEGFR-3 -> bm < 0) push the baseline logit positive
    # -> higher HFS probability under drug. See vignette for sign-convention
    # mapping to the Methods text 'Increasing sVEGFR-3 REL ... associated
    # with increased probability and severity of HFS'.
    e_bm_px0 <- -8.00; label("Slope of effect-compartment sVEGFR-3 relative-change (bm) on PX0 HFS baseline logit (log-odds units per unit bm; Table 3 HFS row 'Slope x|0')") # Table 3 HFS Slope x|0 (RSE 14%)
    e_bm_px1 <- -6.00; label("Slope of bm on PX1 HFS baseline logit (Table 3 HFS row 'Slope x|1')") # Table 3 HFS Slope x|1 (RSE 18%)
    e_bm_px2 <- -3.23; label("Slope of bm on PX2 HFS baseline logit (Table 3 HFS row 'Slope x|2')") # Table 3 HFS Slope x|2 (RSE 43%)
    e_bm_px3 <- -4.75; label("Slope of bm on PX3+ HFS baseline logit (Table 3 HFS row 'Slope x|>=3')") # Table 3 HFS Slope x|>=3 (RSE 32%)

    # --- Effect-compartment rate constant (fixed in source per Methods text) ---
    # Hansson 2013 e85 Results: 'Incorporation of an effect compartment
    # into the model [ke0 = 0.424/hour (fatigue) and 0.347/hour (HFS)]
    # significantly improved both the fatigue and HFS models'. The value
    # is reported point-estimated to three significant figures with no
    # uncertainty -> wrap in fixed().
    lke0 <- fixed(log(0.347)); label("Effect-compartment first-order rate constant ke0 (1/h; Hansson 2013 e85 Results for HFS effect compartment)") # Hansson 2013 e85 Results: ke0 = 0.347/h fixed

    # ----------------------------------------------------------------------
    # IIV on the per-state baseline logits. Table 3 HFS reports omega values
    # (SDs); convert to variances. The PX3+ state IIV is 'NE' (Not Estimated)
    # in Table 3 -> fixed at 0.
    # ----------------------------------------------------------------------
    etaclge_px0 ~ 9.4249  # Table 3 HFS omega x|0 = 3.07 -> var = 3.07^2 = 9.4249 (RSE 67%)
    etaclge_px1 ~ 0.8136  # Table 3 HFS omega x|1 = 0.902 -> var = 0.902^2 = 0.8136 (RSE 54%)
    etaclge_px2 ~ 0.0729  # Table 3 HFS omega x|2 = 0.270 -> var = 0.270^2 = 0.0729 (RSE 118%)
    etaclge_px3 ~ fixed(0.0001) # Table 3 HFS omega x|>=3 = NE (Not Estimated); fixed at tiny non-zero value (rxode2 requires non-zero variance)

    # ----------------------------------------------------------------------
    # Residual-error placeholder
    # The published likelihood is the joint Markov + proportional-odds form
    # implemented via the conditional Y = P_xy assignment under
    # $ESTIMATION ... LAPLACE LIKE. nlmixr2 / rxode2 do not natively express
    # "previous DV"-conditioned discrete likelihoods, so this model file
    # provides the typical-value transition-probability outputs as
    # deterministic functions of time; for nlmixr2 fitting compatibility
    # the observation is taken to be the "expected HFS grade given previous
    # state = 0" output (a continuous 0-3 quantity), with a placeholder
    # additive residual error the user can tune at fit time. See the
    # vignette's Assumptions and deviations section. Mirrors the
    # Hansson_2013c_sunitinib (fatigue) deviation exactly.
    # ----------------------------------------------------------------------
    addSd_hfs_grade <- 0.5
    label("Placeholder additive residual error on the typical-value expected-HFS-grade output (units: ordinal grade); not from the source -- see vignette Assumptions")
  })

  model({
    # 1. Per-cycle drug-exposure summary (mg*h/L AUC).
    auc <- DOSE / CLI

    # 2. Simple-Imax drug effect on sVEGFR-3.
    eff3 <- auc / (EC50_SVEGFR3 + auc)

    # 3. sVEGFR-3 indirect-response turnover (Kin inhibition).
    kout3 <- 1 / MRT_SVEGFR3
    kin3  <- BAS_SVEGFR3 * kout3
    svegfr3(0)     <- BAS_SVEGFR3
    d/dt(svegfr3)  <- kin3 * (1 - eff3) - kout3 * svegfr3

    # 4. Effect-compartment-smoothed sVEGFR-3 relative-change signal.
    # Eq. 2 of the paper uses a delayed bm to drive the proportional-odds
    # logits. The raw bm = (svegfr3 - BAS_SVEGFR3) / BAS_SVEGFR3 is
    # negative under drug; the delayed effect-compartment output (`bm`)
    # tracks the raw signal with first-order delay ke0 = 0.347/h.
    bm_input <- (svegfr3 - BAS_SVEGFR3) / BAS_SVEGFR3
    ke0      <- exp(lke0)
    d/dt(bm) <- ke0 * (bm_input - bm)
    bm(0)    <- 0

    # 5. Cumulative logits per starting state. Each adds the per-state
    # slope shift (e_bm_px<i> * bm) and the per-state shared eta.
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

    # 6. Transition probabilities (each row sums to 1).
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

    # 7. Expected HFS grade given previous = 0 (natural-history /
    # placebo-arm reference output for typical-value visualisation).
    expected_grade_from0 <- 0 * p00 + 1 * p10 + 2 * p20 + 3 * p30
    hfs_grade            <- expected_grade_from0

    hfs_grade ~ add(addSd_hfs_grade)
  })
}
