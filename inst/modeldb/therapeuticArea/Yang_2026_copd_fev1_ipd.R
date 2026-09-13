Yang_2026_copd_fev1_ipd <- function() {
  description <- "Longitudinal individual-patient-data (IPD) disease-progression and drug-effect model of morning trough forced expiratory volume in 1 second (FEV1) in chronic obstructive pulmonary disease (COPD), fit to 2,241 patients in two 24-week randomized fluticasone furoate / vilanterol trials (NCT01053988 and NCT01054885) pooled by Yang 2026. FEV1 = baseline - linear disease progression + Emax dose-response for vilanterol and for fluticasone furoate; no placebo effect was supported. Baseline carries exponential IIV and effects of age, GOLD stage (hockey-stick with a knee at stage 3), current-smoking status and sex; the disease-progression slope and the vilanterol reference efficacy carry additive IIV and GOLD-stage effects. Residual error is a power function of the prediction with its own log-normal IIV. This is the individual-level half of the paper's combined aggregated-data + IPD (ADIPD) model -- see Yang_2026_copd_fev1_adipd_mbma for the 23-compound combined model. There is no PK layer: the drug effect is driven by the per-arm total daily dose supplied as covariate columns."

  reference <- paste(
    "Yang L, Llanos-Paez C, Yang S, Ambery C, Berges A, Kjellsson MC,",
    "Karlsson MO. A Combined Model-Based Meta-Analysis of Aggregated and",
    "Individual FEV1 Data From Randomized COPD Trials.",
    "CPT Pharmacometrics Syst Pharmacol. 2026;15(1):e70059.",
    "doi:10.1002/psp4.70059.",
    "Final parameter estimates are in Supporting Information Table S2;",
    "the model equations are in the Supporting Information section",
    "'NONMEM control stream for the IPD model'.",
    sep = " "
  )

  vignette <- "Yang_2026_copd_fev1"

  # The observation is absolute FEV1 in litres -- the canonical `FEV1`
  # compartment (registered with the uppercase paper spelling), also used by
  # Zhang_2025_dupilumab_fev1 (L) and Jin_2025_benralizumab_fev1 (mL). Distinct
  # from `fev1pp`, the percent-predicted surface used by
  # Harun_2019_cysticFibrosis.

  units <- list(
    time          = "week (weeks since randomization; the disease-progression slope is reported per year and is divided by 52 inside model())",
    dosing        = "ug/day (per-arm TOTAL DAILY dose supplied through the CONMED_<drug>_DOSE covariate columns, NOT as rxode2 dose events; this model has no PK layer)",
    concentration = "L (FEV1 absolute volume, observation FEV1)"
  )

  covariateData <- list(
    AGE = list(
      description        = "Subject age at randomization.",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Centred at 62 years in this model (Supporting Information IPD control stream: 'BAGE = ( 1 + THETA(10)*(AGE - 62))'), which is the pooled median and mean age of the two IPD studies (Table S1: mean 62, median 62, range 40-85). Age enters BOTH the baseline (e_age_base) and the vilanterol reference efficacy (e_age_effref_vi). NOTE the centring differs from the sibling combined model Yang_2026_copd_fev1_adipd_mbma, which centres age at 63.4 years because it pools the 298 aggregated-data studies as well; the two centrings are NOT interchangeable.",
      source_name        = "AGE"
    ),
    SEXF = list(
      description        = "Female-sex indicator. 1 = female, 0 = male.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male). The source parameterised this as SEX with 1 = male (the most common level, 69.7% of the pooled IPD cohort per Table S1) as the reference, so the published coefficient is carried on the female side of the contrast and applies unchanged to SEXF.",
      notes              = "Source column SEX is coded 1 = male; the canonical SEXF is its complement (SEXF = 1 - SEX). The source control stream writes 'IF(SEX.EQ.1) BSEX = 1' and 'IF(SEX.EQ.0) BSEX = (1 + THETA(14))', i.e. the estimated coefficient THETA(14) = -0.248 multiplies the FEMALE level with male as the reference. Re-expressed on SEXF the identical algebra is '1 + e_sexf_base * SEXF', so the published value and its sign are carried over unchanged and no reference-category flip is involved. Direction confirmed twice in the paper text ('female ... related to lower baseline', Abstract and Section 3.3) and against the sibling combined model, where the equivalent centred form (1 + 0.276*(SEX - 0.671)) gives the same female/male baseline ratio of 0.75.",
      source_name        = "SEX (1 = male)"
    ),
    SMOKE = list(
      description        = "Current-smoker indicator at baseline. 1 = current smoker, 0 = non-current (former) smoker.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 (current smoker) is the source's reference level because it is the most common (54.2% of the pooled IPD cohort, Table S1); the estimated coefficient is carried on the non-current-smoker side.",
      notes              = "Source column FL_SMOK, coded 1 = current smoker, which matches the canonical SMOKE coding exactly (no transformation). The control stream writes 'IF(FL_SMOK.EQ.1) BFL_SMOK = 1' and 'IF(FL_SMOK.EQ.0) BFL_SMOK = (1 + THETA(13))', so the coefficient applies to the NON-smoker level; model() encodes this as '1 + e_nonsmoke_base * (1 - SMOKE)'. The paper cautions that this association is not causal -- patients tend to stop smoking at more severe COPD stages (Section 4, refs 26-28).",
      source_name        = "FL_SMOK"
    ),
    DIS_COPD_GOLD = list(
      description        = "GOLD spirometric severity stage at screening, as an ordinal 1-4 category (1 = mild, 2 = moderate, 3 = severe, 4 = very severe).",
      units              = "(ordinal stage 1-4)",
      type               = "ordinal",
      reference_category = "3 (severe), the pooled cohort median, used as the centring constant and as the knee of the hockey-stick relationships.",
      notes              = "Source column FL_COPD; Table S1 footnote 1 defines it as the '% predicted GOLD Stage Category at the screening phase'. The pooled IPD distribution is stage 1: 0.09%, stage 2: 46.5%, stage 3: 44.1%, stage 4: 8.66%, with 0.625% missing (imputed to the median, 3). Carried as a SINGLE ordinal column rather than pre-binned indicators because the source fits a piecewise-LINEAR (hockey-stick) effect in the raw stage number, with separate slopes below and above the median stage 3; indicator columns could not reproduce that form. This follows the single-ordinal-column rationale already recorded for SMOKE_TTFC_SCORE. The stage enters the baseline (two slopes), the vilanterol reference efficacy (two slopes) and the disease-progression slope (one slope).",
      source_name        = "FL_COPD"
    ),
    CONMED_VILANTEROL_DOSE = list(
      description        = "Per-arm total daily vilanterol dose.",
      units              = "ug/day",
      type               = "continuous",
      reference_category = "0 (no vilanterol in this arm, which collapses the vilanterol Emax term to zero).",
      notes              = "Source column avdostotN for the slot whose drgNoN equals 25 (vilanterol). The reference dose against which the reported reference efficacy is quoted is 25 ug/day (control stream 'REFDVIL = 25'), i.e. vilanterol 25 ug q.d. Dose levels present in the two IPD studies are 0 and 25 ug/day. The unit must match the ED50 unit -- ED50_vilanterol = exp(0.709) = 2.03 ug/day.",
      source_name        = "avdostot (drgNo = 25)"
    ),
    CONMED_FLUTICASONEFUROATE_DOSE = list(
      description        = "Per-arm total daily fluticasone furoate dose.",
      units              = "ug/day",
      type               = "continuous",
      reference_category = "0 (no fluticasone furoate in this arm, which collapses the fluticasone furoate Emax term to zero).",
      notes              = "Source column avdostotN for the slot whose drgNoN equals 32 (fluticasone furoate). The reference dose against which the reported reference efficacy is quoted is 100 ug/day (control stream 'REFDFF = 100'), i.e. fluticasone furoate 100 ug q.d. Dose levels present in the two IPD studies are 0, 50, 100 and 200 ug/day (Kerwin 2013 gave 50/25 and 100/25 ug FF/VI; Martinez 2013 gave 100/25 and 200/25 ug FF/VI). The unit must match the ED50 unit -- ED50_fluticasone furoate = exp(2.43) = 11.4 ug/day.",
      source_name        = "avdostot (drgNo = 32)"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 2241L,
    n_studies       = 2L,
    age_range       = "40-85 years (pooled; mean 62, SD 8.8, median 62)",
    age_median      = "62 years",
    weight_range    = "not reported in the source (body weight was not a covariate in this model)",
    sex_female_pct  = 30.3,
    disease_state   = "Moderate-to-very-severe chronic obstructive pulmonary disease. GOLD spirometric stage at screening: stage 1 0.09%, stage 2 46.5%, stage 3 44.1%, stage 4 8.66% (0.625% missing). 54.2% were current smokers at baseline. Mean on-treatment FEV1 1.372 L (SD 0.509, range 0.300-4.140).",
    dose_range      = "Fluticasone furoate 0, 50, 100 or 200 ug q.d. and vilanterol 0 or 25 ug q.d., alone or in combination, over 24 weeks; placebo-controlled.",
    regions         = "Multicentre (both trials were international multicentre studies).",
    trials_included = "NCT01053988 (Kerwin EM et al., Respir Med 2013;107:560-569; n = 1025) and NCT01054885 (Martinez FJ et al., Respir Med 2013;107:550-559; n = 1216). Both are 24-week randomized, smoking-status-stratified, placebo-controlled, double-blind, parallel-group, multicentre studies.",
    notes           = "Demographics are from Supporting Information Table S1. sex_female_pct computed as 1 - 1563/2241 = 30.3% from the reported male counts (681 + 882 = 1563 of 2241). The two studies are also among the 298 studies of the aggregated-data set, but were REMOVED from the aggregated data when building the combined ADIPD model so that no observation is used twice (Section 2.1). The drug-interaction term between fluticasone furoate and vilanterol was tested and NOT retained (dOFV < 0.1; Section 3.3), so the two drug effects are purely additive here."
  )

  ini({
    # ==================================================================
    # Structural model (Supporting Information IPD control stream $PRED;
    # the same skeleton as paper Equation 10 without the placebo,
    # post-bronchodilator and background-therapy terms):
    #
    #   FEV1 = B - DP + P + drug effects
    #   B    = base * BCOV * exp(etabase)
    #   DP   = (dps * DPSCOV + etadps) * t / 52
    #   P    = 0                       (placebo Emax fixed to zero)
    #
    # All VALUES below are the final estimates from Supporting Information
    # Table S2. They are deliberately NOT taken from the control stream's
    # $THETA / $OMEGA records, which carry that run's INITIAL estimates
    # (e.g. 1.18 vs the final 1.17 for baseline).
    # ==================================================================

    base <- 1.17    ; label("Typical baseline FEV1 at the reference covariate values (L)")           # Table S2 'Typical baseline FEV1 (L)' = 1.17 (RSE 1.10%)
    dps  <- 0.0241  ; label("Typical disease-progression slope, a DECLINE in FEV1 (L/year)")         # Table S2 'Disease progression slope (L/year)' = 0.0241 (RSE 35.60%); enters model() with a minus sign, so a positive value is a decline

    # Placebo effect: the IPD model retained NO placebo effect (Section 3.3,
    # 'the model structure of no placebo effect ... were applied'); the source
    # control stream carries it as '(0) FIX ; Placebo.Emax (3)'.
    pmx  <- fixed(0) ; label("Immediate placebo effect (L)") # IPD control stream $THETA (3): '(0) FIX ; Placebo.Emax'; no placebo effect was supported, Section 3.3

    # ---- Drug effects: Emax in total daily dose ------------------------
    # The source parameterises each drug by its efficacy AT A REFERENCE DOSE
    # rather than by Emax directly:
    #   Emax = EffRef / RefDose * (ED50 + RefDose)
    # so that EffRef is the effect at RefDose. RefDose is 25 ug/day for
    # vilanterol and 100 ug/day for fluticasone furoate.
    effref_vi <- 0.117  ; label("Vilanterol efficacy at the reference dose of 25 ug/day (L)")           # Table S2 'Reference efficacy of vilanterol 25 ug q.d. (L)' = 0.117 (RSE 6%)
    effref_ff <- 0.0304 ; label("Fluticasone furoate efficacy at the reference dose of 100 ug/day (L)") # Table S2 'Reference efficacy of fluticasone furoate 100 ug q.d. (L)' = 0.0304 (RSE 13.30%)

    # Both ED50s were FIXED in the IPD model (Table S2 footnote 2: '*means the
    # parameter was fixed'); they are estimated in the sibling combined model.
    led50_vi <- fixed(0.709) ; label("log ED50 for vilanterol (log ug/day)")           # Table S2 'Log of ED50 for vilanterol' = 0.709* (fixed); ED50 = exp(0.709) = 2.03 ug/day
    led50_ff <- fixed(2.43)  ; label("log ED50 for fluticasone furoate (log ug/day)")  # Table S2 'Log of ED50 for fluticasone furoate q.d.' = 2.43* (fixed); ED50 = exp(2.43) = 11.4 ug/day

    # ---- Covariate effects on baseline ---------------------------------
    # All are LINEAR-PROPORTIONAL deviations, i.e. (1 + theta * (cov - ref)).
    e_age_base      <- -0.0144 ; label("Fractional change in baseline FEV1 per year of age above 62")                        # Table S2 'Covariate effect of age on baseline' = -0.0144 (RSE 4.10%); control stream 'BAGE = (1 + THETA(10)*(AGE - 62))'
    e_gold_lo_base  <- -0.47   ; label("Fractional change in baseline FEV1 per GOLD stage above 3, for stages <= 3")         # Table S2 'Covariate effect of disease severity on baseline (first slope of hockey-stick model, cutoff value is the median of disease severity =3)' = -0.47 (RSE 3.30%)
    e_gold_hi_base  <- -0.368  ; label("Fractional change in baseline FEV1 per GOLD stage above 3, for stages > 3")          # Table S2 'Covariate effect of disease severity on baseline (second slope of hockey-stick model ...)' = -0.368 (RSE 3.40%)
    e_nonsmoke_base <- -0.0419 ; label("Fractional change in baseline FEV1 for a non-current smoker vs a current smoker")    # Table S2 'Covariate effect of non-smoke on baseline' = -0.0419 (RSE 24.60%)
    e_sexf_base     <- -0.248  ; label("Fractional change in baseline FEV1 for a female vs a male subject")                  # Table S2 'Covariate effect of female on baseline' = -0.248 (RSE 3.30%)

    # ---- Covariate effects on the vilanterol reference efficacy --------
    e_age_effref_vi     <- -0.0119 ; label("Fractional change in vilanterol efficacy per year of age above 62")              # Table S2 'Covariate effect of age on efficacy of vilanterol' = -0.0119 (RSE 33.60%)
    e_gold_lo_effref_vi <- -0.163  ; label("Fractional change in vilanterol efficacy per GOLD stage above 3, stages <= 3")   # Table S2 'Covariate effect of disease severity on efficacy of vilanterol (first slope)' = -0.163 (RSE 56.40%)
    e_gold_hi_effref_vi <- -0.533  ; label("Fractional change in vilanterol efficacy per GOLD stage above 3, stages > 3")    # Table S2 'Covariate effect of disease severity on efficacy of vilanterol (second slope)' = -0.533 (RSE 26.60%)

    # ---- Covariate effect on the disease-progression slope -------------
    e_gold_dps <- -1.89 ; label("Fractional change in the disease-progression slope per GOLD stage above 3")                 # Table S2 'Covariate effect of disease severity on disease progression slope' = -1.89 (RSE 51%); a single slope, no hockey stick

    # ---- Residual error: SD = powSd_FEV1 * FEV1^powExp_FEV1 * exp(etapowSd_FEV1) --------
    powExp_FEV1 <- 0.645 ; label("Power of the prediction in the residual-error model (unitless)")  # Table S2 'Power on prediction of residual error' = 0.645 (RSE 4.50%)
    powSd_FEV1  <- sqrt(0.00603) ; label("Residual-error scale (SD of EPS1; L^(1-powExp_FEV1))")         # Table S2 'Variance of residual error' = 0.00603 (RSE 2.50%); NONMEM $SIGMA is a VARIANCE, so the SD nlmixr2 wants is sqrt(0.00603) = 0.0777

    # ==================================================================
    # Random effects. Table S2 reports every one of these on the VARIANCE
    # scale ('IIV variance of ...'), which is the scale nlmixr2's `~` wants,
    # so no back-transformation is applied. etabase is EXPONENTIAL on the
    # baseline; etadps and etaeffref_vi are ADDITIVE on their parameters
    # (Section 3.3: 'exponential IIV on baseline, and additive IIV on
    # disease progression/vilanterol efficacy'); etapowSd_FEV1 is exponential on
    # the residual-error magnitude.
    # ==================================================================

    etabase       ~ 0.0534   # Table S2 'IIV variance of baseline' = 0.0534 (RSE 3.50%)
    etadps        ~ 0.0605   # Table S2 'IIV variance of disease progression slope' = 0.0605 (RSE 8.60%)
    etaeffref_vi  ~ 0.00698  # Table S2 'IIV variance of the efficacy of vilanterol' = 0.00698 (RSE 19.20%)
    etapowSd_FEV1        ~ 0.128    # Table S2 'IIV variance of residual error' = 0.128 (RSE 6.10%)
  })

  model({
    # ---- 1. Covariate multipliers on the baseline ----------------------
    # Hockey-stick in the GOLD stage with the knee at the cohort median
    # stage 3: a separate slope below and above the knee. Both branches
    # pass through 1 at GOLD stage 3, so the function is continuous.
    gold_dev  <- DIS_COPD_GOLD - 3
    b_gold    <- 1 + (e_gold_lo_base * (DIS_COPD_GOLD <= 3) +
                        e_gold_hi_base * (DIS_COPD_GOLD > 3)) * gold_dev
    b_age     <- 1 + e_age_base * (AGE - 62)
    # The source codes the smoking coefficient on the NON-current-smoker
    # level with current smoker as the reference, hence (1 - SMOKE).
    b_smoke   <- 1 + e_nonsmoke_base * (1 - SMOKE)
    # The source codes the sex coefficient on the FEMALE level with male as
    # the reference; SEXF is 1 for a female subject.
    b_sexf    <- 1 + e_sexf_base * SEXF

    # ---- 2. Baseline FEV1 with exponential IIV -------------------------
    base_i <- base * b_age * b_gold * b_smoke * b_sexf * exp(etabase)

    # ---- 3. Disease progression (a DECLINE; subtracted below) ----------
    # A single GOLD-stage slope, no hockey stick, and ADDITIVE IIV.
    dps_i <- dps * (1 + e_gold_dps * gold_dev) + etadps
    dp    <- dps_i * t / 52

    # ---- 4. Placebo effect: a constant for t > 0 (fixed to zero here) --
    plac <- pmx * (t > 0)

    # ---- 5. Drug effects ------------------------------------------------
    # Emax is reconstructed from the efficacy at the reference dose:
    #   Emax = EffRef / RefDose * (ED50 + RefDose).
    # Age and GOLD stage modify the vilanterol efficacy with the same
    # linear-proportional / hockey-stick forms used on the baseline; the
    # fluticasone furoate efficacy carries no covariates and no IIV
    # (its IIV was fixed to zero in the source $OMEGA).
    ed50_vi <- exp(led50_vi)
    ed50_ff <- exp(led50_ff)

    d_gold_vi <- 1 + (e_gold_lo_effref_vi * (DIS_COPD_GOLD <= 3) +
                        e_gold_hi_effref_vi * (DIS_COPD_GOLD > 3)) * gold_dev
    d_age_vi  <- 1 + e_age_effref_vi * (AGE - 62)

    emax_vi <- (effref_vi + etaeffref_vi) / 25 * (ed50_vi + 25) * d_age_vi * d_gold_vi
    emax_ff <- effref_ff / 100 * (ed50_ff + 100)

    # Hyperbolic Emax in the per-arm total daily dose. A zero dose column
    # collapses its term to zero, so a placebo arm contributes nothing and
    # the drug-flag bookkeeping of the source control stream is redundant.
    de_vi <- emax_vi * CONMED_VILANTEROL_DOSE / (CONMED_VILANTEROL_DOSE + ed50_vi)
    de_ff <- emax_ff * CONMED_FLUTICASONEFUROATE_DOSE / (CONMED_FLUTICASONEFUROATE_DOSE + ed50_ff)

    # No drug effect is expressed at or before randomization (source
    # control stream: 'IF (TRTNO.GT.1.AND.TIME.GT.0) TREAT = 1').
    # The fluticasone furoate / vilanterol interaction term was tested and
    # NOT retained, so the two effects are additive.
    drug_all <- (t > 0) * (de_vi + de_ff)

    # ---- 6. Prediction and residual error -------------------------------
    FEV1 <- base_i - dp + plac + drug_all

    # Power residual error whose magnitude carries its own log-normal IIV
    # (source control stream: 'W = F**THETA(9)*EXP(ETA(4))', 'Y = F + EPS(1)*W').
    ruvScale <- powSd_FEV1 * exp(etapowSd_FEV1)
    FEV1 ~ pow(ruvScale, powExp_FEV1)
  })
}
