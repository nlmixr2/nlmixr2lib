Hansson_2013b_sunitinib <- function() {
  description <- "Population PD tumor growth inhibition model for sunitinib in adults with imatinib-resistant gastrointestinal stromal tumours (GIST). The longitudinal sum of longest tumor diameters (SLD) is modelled as exponential growth (KG) with three additive shrinkage drivers -- exposure-driven (KDRUG * AUC), and the model-predicted relative-from-baseline changes in soluble KIT (sKIT) and soluble VEGFR-3 (sVEGFR-3) acting through the rate constants KSKIT and KVEGFR3 -- modulated by an exponential time-dependent resistance term (LAMBDA). The two soluble-biomarker time-courses are simulated in-model as 3 indirect-response compartments (treated sKIT, placebo / untreated sKIT, sVEGFR-3) driven by simple-Imax inhibition of Kin with the per-cycle exposure summary AUC = DOSE / CLI; the placebo-arm sKIT compartment carries a linear disease-progression term. The PD model has no PK ODE and consumes individual posthoc upstream-PD parameters (BAS_SKIT, MRT_SKIT, EC50_SKIT, SLOPE_SKIT, BAS_SVEGFR3, MRT_SVEGFR3, EC50_SVEGFR3) plus posthoc upstream-PK clearance (CLI) and observed baseline tumor size (TUMSZ, mm) as data covariates. IIV on tumor growth (KG), drug effect (KDRUG), and the sKIT-driven shrinkage rate (KSKIT); LAMBDA's IIV is held at zero in the source; the IPP-style baseline-residual eta is fixed-variance 1 with proportional residual scaling (Dansirikul / Silber / Karlsson 2008)."
  reference <- paste(
    "Hansson EK, Amantea MA, Westwood P, Milligan PA, Houk BE,",
    "French J, Karlsson MO, Friberg LE.",
    "PKPD modeling of VEGF, sVEGFR-2, sVEGFR-3, and sKIT as predictors of",
    "tumor dynamics and overall survival following sunitinib treatment in GIST.",
    "CPT Pharmacometrics Syst Pharmacol. 2013;2(11):e84.",
    "doi:10.1038/psp.2013.61.",
    "DDMORE Foundation Model Repository: DDMODEL00000198.",
    "Companion biomarker indirect-response model:",
    "modellib('Hansson_2013a_sunitinib') (DDMODEL00000197).",
    sep = " "
  )
  vignette <- "Hansson_2013b_sunitinib"
  units <- list(time = "hour", dosing = "mg", concentration = "mm (tumor SLD)")
  ddmore_id <- "DDMODEL00000198"
  replicate_of <- NULL

  covariateData <- list(
    DOSE = list(
      description        = "Current administered sunitinib daily dose (mg) carried as a time-varying data column. Set to 0 during off-cycles (4 weeks on / 2 weeks off in the Hansson 2013 GIST cohort) or for placebo subjects so the derived AUC = DOSE / CLI becomes 0.",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "DDMORE-bundle simulated dataset reports DOSE = 50 (treated) or 0 (placebo) at every record. The .mod feeds DOSE into AUC = DOS/CL in $PK at every event call, producing a per-cycle daily-AUC equivalent (mg*h/L). For typical-cohort vignette simulations the value is held at 50 mg during 4-week on-cycles and 0 mg during 2-week off-cycles.",
      source_name        = "DOS"
    ),
    CLI = list(
      description        = "Individual posthoc total plasma clearance (L/h) of sunitinib from the paper's upstream 2-compartment popPK fit. Per-subject, time-fixed.",
      units              = "L/h",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Required input. The DDMORE-bundle simulated dataset carries CLI = 32.819 L/h for subject 1; this value (broadly consistent with the typical sunitinib CL reported by Houk et al. 2010) is used as the typical-value reference for the validation vignette's virtual cohort. For a re-fit or new-population simulation the user must supply individual CL drawn from a sunitinib popPK model.",
      source_name        = "CL"
    ),
    TUMSZ = list(
      description        = "Observed baseline tumor size (sum of longest diameters of target lesions, mm) at study entry; per-subject, time-fixed. Used both as the deterministic component of the tumor ODE initial condition (`tumor(0) = TUMSZ * (1 + etaibase * propSd_tumor)`, the IPP-style proportional baseline-residual construction of Dansirikul 2008) and as the typical-value scale for the proportional residual error.",
      units              = "mm",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Required input. The .mod reads OBASE dynamically from DV at TIME=0, FLAG=4 (`IF(TIME.EQ.0.AND.FLAG.EQ.4) THEN OBASE = DV`). nlmixr2 / rxode2 cannot replicate the in-record assignment idiom, so the observed baseline is supplied as a per-subject covariate column instead. Median baseline SLD in the cohort: 194 mm (study 1004), 108 mm (study 1047), 166 mm (study 1045), 255 mm (study 013) per Hansson 2013 Table 1.",
      source_name        = "(read from DV at TIME=0/FLAG=4 via .mod IF block)"
    ),
    BAS_SKIT = list(
      description        = "Individual posthoc baseline sKIT (pg/mL) from the upstream Hansson 2013a biomarker indirect-response PD fit (DDMODEL00000197); per-subject, time-fixed; used as the initial condition for both the treated and placebo sKIT compartments.",
      units              = "pg/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Required input. The Hansson 2013 e84 paper reports a typical sKIT baseline of 39200 pg/mL (Table 2). For a typical-value simulation set every subject to that value; for an IIV simulation either (a) simulate from `Hansson_2013a_sunitinib` and take the per-subject posthoc baseline, or (b) draw from a log-normal centred at 39200 pg/mL with the upstream IIV (~50% CV).",
      source_name        = "SBAS"
    ),
    MRT_SKIT = list(
      description        = "Individual posthoc mean residence time of sKIT (h) from the upstream Hansson 2013a biomarker indirect-response PD fit (DDMODEL00000197); per-subject, time-fixed; appears as kout_skit = 1 / MRT_SKIT inside model().",
      units              = "h",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Required input. The Hansson 2013 e84 paper reports a typical sKIT MRT of 101 days = 2424 h (Table 2; ~2430 h matching Hansson_2013a_sunitinib's typical value). Same population-input strategy as BAS_SKIT.",
      source_name        = "SMRT"
    ),
    EC50_SKIT = list(
      description        = "Individual posthoc EC50 of the simple-Imax drug effect on sKIT (mg*h/L AUC) from the upstream Hansson 2013a biomarker indirect-response PD fit (DDMODEL00000197); per-subject, time-fixed; appears in the drug-effect term eff_skit = auc / (EC50_SKIT + auc).",
      units              = "mg*h/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Required input. The Hansson 2013 e84 paper reports a typical (common across the four biomarkers) IC50 of 1.0 mg*h/L (Table 2). Same population-input strategy as BAS_SKIT.",
      source_name        = "SEC5"
    ),
    SLOPE_SKIT = list(
      description        = "Individual posthoc linear disease-progression slope on the placebo / untreated sKIT compartment (1/h) from the upstream Hansson 2013a biomarker indirect-response PD fit (DDMODEL00000197); per-subject, time-fixed; appears in the placebo-arm Kin expression dps = BAS_SKIT * (1 + SLOPE_SKIT * t).",
      units              = "1/h",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Required input. The Hansson 2013 e84 paper reports a typical disease-progression slope of 0.0261 / month shared between VEGF and sKIT (Table 2); converted to 1/h that is approximately 3.5e-5 / h, matching Hansson_2013a_sunitinib's typical value. Same population-input strategy as BAS_SKIT.",
      source_name        = "SLO"
    ),
    BAS_SVEGFR3 = list(
      description        = "Individual posthoc baseline sVEGFR-3 (pg/mL) from the upstream Hansson 2013a biomarker indirect-response PD fit (DDMODEL00000197); per-subject, time-fixed; used as the initial condition for the svegfr3 state and as the denominator in the relative-change driver bm_svegfr3 = (svegfr3 - BAS_SVEGFR3) / BAS_SVEGFR3.",
      units              = "pg/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Required input. The Hansson 2013 e84 paper reports a typical sVEGFR-3 baseline of 63900 pg/mL (Table 2). Same population-input strategy as BAS_SKIT.",
      source_name        = "BAS3"
    ),
    MRT_SVEGFR3 = list(
      description        = "Individual posthoc mean residence time of sVEGFR-3 (h) from the upstream Hansson 2013a biomarker indirect-response PD fit (DDMODEL00000197); per-subject, time-fixed; appears as kout_svegfr3 = 1 / MRT_SVEGFR3 inside model().",
      units              = "h",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Required input. The Hansson 2013 e84 paper reports a typical sVEGFR-3 MRT of 16.7 days = 401 h (Table 2). Same population-input strategy as BAS_SKIT.",
      source_name        = "MRT3"
    ),
    EC50_SVEGFR3 = list(
      description        = "Individual posthoc EC50 of the simple-Imax drug effect on sVEGFR-3 (mg*h/L AUC) from the upstream Hansson 2013a biomarker indirect-response PD fit (DDMODEL00000197); per-subject, time-fixed; appears in the drug-effect term eff_svegfr3 = auc / (EC50_SVEGFR3 + auc).",
      units              = "mg*h/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Required input. The Hansson 2013 e84 paper reports a typical (common across the four biomarkers) IC50 of 1.0 mg*h/L (Table 2). Same population-input strategy as BAS_SKIT.",
      source_name        = "EC53"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 303L,
    n_studies      = 4L,
    age_range      = "adults with imatinib-resistant GIST (Hansson 2013 e84 Table 1 lists baseline tumor size by study but does not break out age / weight / sex / race in the on-disk Methods section of the trimmed PDF)",
    weight_range   = "not reported in the on-disk paper trimmed text",
    sex_female_pct = NA_real_,
    race_ethnicity = NULL,
    disease_state  = "Imatinib-resistant gastrointestinal stromal tumours (GIST). Pooled four sunitinib studies: Demetri 2006 (NCT00075218 / study 1004; placebo-controlled phase III; 202 active + 47 placebo), George 2009 (NCT00428220 / study 1047; phase II continuous-dosing 37.5 mg QD; n=13 in this analysis subset), Shirao 2010 (study 1045; Japanese phase I/II; 25-75 mg QD on a 4/2 schedule; n=36), Maki 2005 (study 013; phase I/II 25-75 mg QD on a 2/1 or 2/2 schedule; n=52).",
    dose_range     = "Sunitinib 25-75 mg PO QD on a 4/2, 2/2, 2/1 (weeks on / weeks off) or continuous treatment schedule. The largest cohort (study 1004) used 50 mg QD on a 4/2 schedule. Placebo arm: no sunitinib (study 1004 only).",
    regions        = "Phase III multinational (study 1004); Japanese phase I/II (study 1045); other studies regions not stated in the trimmed paper text.",
    biomarkers     = "Tumor size endpoint: sum of longest tumor diameters (SLD, mm). Baseline (median, range): 194 (35-822) mm in study 1004; 108 (29-191) mm in study 1047; 166 (31-644) mm in study 1045; 255 (55-687) mm in study 013 (Hansson 2013 e84 Table 1). Two soluble-biomarker time-courses driven from the upstream Hansson 2013a fit: sKIT and sVEGFR-3.",
    notes          = "n_subjects = 303 reported in Hansson 2013 e84 Methods (`pooled four clinical studies, comprising a total of 303 patients with imatinib-resistant GIST`) and confirmed by the .lst header (`TOT. NO. OF INDIVIDUALS: 303`). The .lst also reports `TOT. NO. OF OBS RECS: 973` (tumor SLD observations only; FLAG=4 records). Detailed baseline demographics (age, weight, sex, race) at the cohort level are not in the trimmed paper text; populating those keys requires reading the original paper's untrimmed PDF or its supplement."
  )

  ini({
    # ----------------------------------------------------------------------
    # Final estimates from Output_real_TGI_GIST.lst FINAL PARAMETER ESTIMATE
    # block (post-MINIMIZATION SUCCESSFUL). Cross-checked against Hansson 2013
    # CPT Pharmacometrics Syst Pharmacol 2013;2:e84 Table 3 (Final tumor
    # growth inhibition model parameter estimates). The .lst's final estimates
    # equal the .mod's $THETA initial values to four significant figures
    # because the run started near-converged; both equal the published Table 3
    # values.
    #
    # Source paper THETA convention: rate constants reported in 1/week. The
    # .mod and the model() block convert to 1/h via the standard /24/7
    # multiplier so all rates carry units of `units$time = "hour"`.
    # ----------------------------------------------------------------------

    # Tumor growth rate (paper Table 3: KG = 0.0118/week = 7.024e-5 /h).
    lkg     <- log(0.0118 / 24 / 7); label("Tumor growth rate constant KG (1/h; paper Table 3 = 0.0118/week)")  # .lst TH 1 = 1.18E-02 /week
    # sKIT-driven tumor shrinkage rate (paper Table 3: KsKIT = -0.00282/week
    # / AUC; the negative sign in the paper indicates effect direction. The
    # .mod stores the magnitude as a positive THETA(2) and inverts the sign
    # in the tumor ODE via `(-SKIT)` -- see model() block). Magnitude
    # 0.00282/week = 1.679e-5 /h.
    lksk    <- log(0.00282 / 24 / 7); label("sKIT-driven tumor shrinkage rate KSKIT magnitude (1/h; paper Table 3 = 0.00282/week)") # .lst TH 2 = 2.82E-03 /week
    # Resistance / regrowth rate (paper Table 3: lambda = 0.0217/week =
    # 1.292e-4 /h).
    llam    <- log(0.0217 / 24 / 7); label("Resistance appearance / tumor regrowth rate LAMBDA (1/h; paper Table 3 = 0.0217/week)") # .lst TH 3 = 2.17E-02 /week
    # Proportional residual error (paper Table 3: 12.5%); also the scale
    # multiplier on the IPP-style baseline-residual eta (`W1 = THETA(4)*OBASE`
    # in the .mod; see model() block).
    propSd  <- 0.125; label("Proportional residual error on tumor SLD (fraction); also IPP baseline-residual scale (paper Table 3 = 12.5%)") # .lst TH 4 = 1.25E-01
    # Exposure-driven shrinkage rate (paper Table 3: KDRUG = 0.0050/week per
    # AUC unit = 2.994e-5 /h per (mg*h/L)).
    lkdrug  <- log(0.00503 / 24 / 7); label("Exposure-driven tumor shrinkage rate KDRUG (1/h per (mg*h/L) AUC; paper Table 3 = 0.0050/week)") # .lst TH 5 = 5.03E-03 /week
    # sVEGFR-3-driven tumor shrinkage rate (paper Table 3: KsVEGFR3 =
    # -0.0371/week; same effect-direction convention as KSKIT, encoded as
    # positive magnitude with sign-inversion in the tumor ODE).
    lkv3    <- log(0.0371 / 24 / 7); label("sVEGFR-3-driven tumor shrinkage rate KVEGFR3 magnitude (1/h; paper Table 3 = 0.0371/week)") # .lst TH 6 = 3.71E-02 /week

    # ----------------------------------------------------------------------
    # Inter-individual variability -- Output_real_TGI_GIST.lst FINAL OMEGA
    # block. omega^2 values on the internal log-normal scale; Hansson 2013
    # Table 3 reports CV% as `sqrt(omega2) * 100` (the small-omega
    # approximation), so CV% values back-compute as: sqrt(0.290)=54%,
    # sqrt(5.91)=243%, sqrt(1.42)=119% -- all match Table 3.
    #
    # OMEGA(3) (LAMBDA) and OMEGA(5) (IBASE) are FIX in the source. The
    # zero-variance LAMBDA eta is omitted from ini() entirely (a fixed-zero
    # IIV is mathematically the absence of IIV); the unit-variance IBASE eta
    # is retained as `etaibase ~ fixed(1)` because it is consumed inside the
    # baseline initial-condition expression in model() to draw the baseline
    # around the observed TUMSZ with a residual std of `propSd * TUMSZ`
    # (the IPP construction of Dansirikul / Silber / Karlsson 2008).
    # ----------------------------------------------------------------------
    etalkg    ~ 0.290           # .lst FINAL OMEGA(1,1) = 2.90E-01 (paper Table 3 IIV CV = 54%)
    etalksk   ~ 5.91            # .lst FINAL OMEGA(2,2) = 5.91E+00 (paper Table 3 IIV CV = 243%)
    etalkdrug ~ 1.42            # .lst FINAL OMEGA(4,4) = 1.42E+00 (paper Table 3 IIV CV = 119%)
    etaibase  ~ fixed(1)        # .lst FINAL OMEGA(5,5) = 1.00E+00 FIX (IPP baseline-residual eta; std = 1)
  })

  model({
    # 1. Per-cycle drug-exposure summary (mg*h/L AUC). mg / (L/h) = mg*h/L;
    #    matches the .mod's `AUC = DOS/CL` line.
    auc <- DOSE / CLI

    # 2. Individual structural parameters (lognormal IIV; paper Table 3).
    kg     <- exp(lkg     + etalkg)
    ksk    <- exp(lksk    + etalksk)
    lam    <- exp(llam)              # OMEGA(3) FIX 0 in source -- no IIV on lambda
    kdrug  <- exp(lkdrug  + etalkdrug)
    kv3    <- exp(lkv3)              # the .mod assigns no eta to KVEG3 (KVEG3 = TVKV3, no EXP(ETA))

    # 3. Biomarker indirect-response sub-model (sKIT treated, sKIT placebo,
    #    sVEGFR-3) with simple-Imax drug effect on Kin. All three Kin / Kout
    #    rates and the linear disease-progression slope come from the
    #    individual posthoc upstream-PD covariates BAS_SKIT / MRT_SKIT /
    #    EC50_SKIT / SLOPE_SKIT and BAS_SVEGFR3 / MRT_SVEGFR3 / EC50_SVEGFR3
    #    -- see covariateData notes for population-strategy options.
    kout_skit    <- 1 / MRT_SKIT
    kout_svegfr3 <- 1 / MRT_SVEGFR3

    eff_skit    <- auc / (EC50_SKIT    + auc)   # IMAX1 = 1 in .mod
    eff_svegfr3 <- auc / (EC50_SVEGFR3 + auc)

    # Placebo / untreated sKIT carries the linear DP term; treated sKIT
    # has the DP-modulated Kin scaled by (1 - eff_skit).
    dps_skit  <- BAS_SKIT * (1 + SLOPE_SKIT * t)
    kin_skit  <- dps_skit * kout_skit
    kin_svegfr3 <- BAS_SVEGFR3 * kout_svegfr3

    # Initial conditions. Treated and placebo sKIT both start at the
    # upstream-fit per-subject baseline; sVEGFR-3 likewise. Tumor IC uses the
    # IPP-style baseline-residual construction:
    #   IBASE = OBASE + ETA(5)*W1, W1 = THETA(4)*OBASE
    #   = TUMSZ * (1 + etaibase * propSd),  with etaibase ~ N(0, 1)
    # propSd is the same proportional residual SD that drives the tumor
    # observation residual error -- the paper / .mod intentionally tie
    # the baseline-residual scale to the observation-residual scale.
    skit_drug(0) <- BAS_SKIT
    skit_pla(0)  <- BAS_SKIT
    svegfr3(0)   <- BAS_SVEGFR3
    tumor(0)     <- TUMSZ * (1 + etaibase * propSd)

    # 4. Biomarker ODEs. Drug inhibits production of treated sKIT and of
    #    sVEGFR-3; placebo sKIT has no drug term.
    d/dt(skit_drug) <- kin_skit * (1 - eff_skit) - kout_skit * skit_drug
    d/dt(skit_pla)  <- kin_skit                  - kout_skit * skit_pla
    d/dt(svegfr3)   <- kin_svegfr3 * (1 - eff_svegfr3) - kout_svegfr3 * svegfr3

    # 5. Biomarker-driven tumor-shrinkage contributions. The .mod defines:
    #      SKIT  = ((A(1) - A(2)) / A(2)) * KSKIT   (treated - placebo) / placebo
    #      VEG3  = ((A(3) - BASE3) / BASE3) * KVEG3 (sVEGFR3 - baseline) / baseline
    #    and folds them into the tumor ODE as `(AUC1 + (-SKIT) + (-VEG3))`,
    #    i.e. with sign-inversion so the contributions are POSITIVE under
    #    drug (when the biomarkers are depressed below their reference). The
    #    same end result is obtained here by writing the relative-change with
    #    the numerator's sign already inverted (`reference - state`), so the
    #    shrinkage contributions stay positive without further sign
    #    manipulation. Reference for sKIT is the contemporaneous placebo /
    #    untreated sKIT compartment (skit_pla), NOT the time-fixed baseline
    #    BAS_SKIT -- skit_pla drifts upward over time via the linear
    #    disease-progression term, so the relative-change driver tracks the
    #    drug-induced deviation from the simulated untreated trajectory.
    skit_eff    <- ((skit_pla     - skit_drug) / skit_pla)    * ksk   # = ((A(2) - A(1)) / A(2)) * KSKIT
    svegfr3_eff <- ((BAS_SVEGFR3  - svegfr3)   / BAS_SVEGFR3) * kv3   # = ((BASE3 - A(3)) / BASE3) * KVEG3
    auc_eff     <- auc * kdrug

    # 6. Tumor SLD ODE. Hansson 2013 e84 Eq. 4:
    #      dY/dt = KG * Y - (AUC*KDRUG + RC_sKIT*KsKIT + RC_sVEGFR3*KsVEGFR3)
    #              * exp(-LAMBDA * t) * Y
    # where RC_X is the relative change of biomarker X from baseline (so
    # that the parenthesised group is positive under drug, when biomarkers
    # have been suppressed) and the exp(-LAMBDA*t) attenuates the drug /
    # biomarker effects over time, capturing acquired resistance.
    d/dt(tumor) <- kg * tumor - (auc_eff + skit_eff + svegfr3_eff) * exp(-lam * t) * tumor

    # 7. Observation: tumor SLD with proportional residual error (paper
    #    Table 3: 12.5%; .mod $ERROR `Y = IPRED + IPRED*THETA(4)*EPS(1)`).
    tumor ~ prop(propSd)
  })
}
