Hansson_2013_sunitinib_myelosuppression <- function() {
  description <- "Semi-physiological Friberg-Karlsson myelosuppression model for sunitinib in adults with imatinib-resistant gastrointestinal stromal tumours (GIST). Absolute neutrophil count (ANC) is described by a self-renewing proliferating progenitor pool, three transit compartments reflecting cell maturation, and a circulating-neutrophil pool, with an Emax drug-effect function driven by the relative change in soluble VEGFR-3 from baseline (sVEGFR-3 REL) inhibiting proliferation and a (ANC0/circ)^gamma feedback term. sVEGFR-3 itself is simulated in-model as a one-compartment indirect-response turnover with simple-Imax inhibition of Kin by the per-cycle drug-exposure summary AUC = DOSE / CLI. The PD model has no PK ODE and consumes individual posthoc upstream-PD parameters (BAS_SVEGFR3, MRT_SVEGFR3, EC50_SVEGFR3) plus posthoc upstream-PK clearance (CLI) and a Japanese-cohort indicator (RACE_JAPANESE) as data covariates. The Japanese-cohort indicator switches the typical baseline ANC0 between 4.94 (non-Japanese) and 3.69 (10^9/L) per Hansson 2013 Table 2."
  reference <- paste(
    "Hansson EK, Ma G, Amantea MA, French J, Milligan PA, Friberg LE,",
    "Karlsson MO.",
    "PKPD modeling of predictors for adverse effects and overall survival",
    "in sunitinib-treated patients with GIST.",
    "CPT Pharmacometrics Syst Pharmacol. 2013;2(11):e85.",
    "doi:10.1038/psp.2013.62.",
    "Sister model files from the same paper:",
    "modellib('Hansson_2013_sunitinib_dbp'),",
    "modellib('Hansson_2013c_sunitinib') [fatigue],",
    "modellib('Hansson_2013_sunitinib_hfs'),",
    "modellib('Hansson_2013_sunitinib_os').",
    "Upstream sVEGFR-3 biomarker dynamics adapted from",
    "Hansson 2013 (CPT Pharmacometrics Syst Pharmacol 2013;2:e84,",
    "doi:10.1038/psp.2013.61, DDMODEL00000197); see",
    "modellib('Hansson_2013a_sunitinib').",
    "Friberg myelosuppression backbone:",
    "Friberg LE et al. J Clin Oncol 2002;20(24):4713-4721,",
    "doi:10.1200/JCO.2002.02.140.",
    sep = " "
  )
  vignette <- "Hansson_2013_sunitinib_myelosuppression"
  units <- list(time = "hour", dosing = "mg", concentration = "10^9 cells/L (ANC)")

  covariateData <- list(
    DOSE = list(
      description        = "Current administered sunitinib daily dose (mg) carried as a time-varying data column. Set to 0 during off-cycles or for placebo subjects so the derived AUC = DOSE / CLI becomes 0.",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "The paper Methods describes sunitinib administered at 25-75 mg PO QD on 4/2, 2/2, 2/1, or continuous schedules across studies 1004, 1047, 1045, and 013 (Table 1). For typical-cohort vignette simulations the value is held at 50 mg during on-cycles of a 4/2 schedule, 0 mg during off-cycles, matching the largest cohort (Demetri 2006 / Study 1004).",
      source_name        = "DOSE"
    ),
    CLI = list(
      description        = "Individual posthoc total plasma clearance (L/h) of sunitinib from the paper's upstream 2-compartment popPK fit. Per-subject, time-fixed.",
      units              = "L/h",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Required input. The Hansson 2013 e85 Methods describes the upstream popPK as a previously developed 2-compartment model (Houk et al. 2009 Clin Cancer Res 15:2497-2506; that popPK is not packaged in nlmixr2lib at extraction time). The companion Hansson_2013a / Hansson_2013c sunitinib model files use a typical-value reference of 32.819 L/h for typical-cohort vignettes; the same value is used here so the in-model svegfr3 dynamics match the upstream Hansson 2013a typical-value figure.",
      source_name        = "CL"
    ),
    BAS_SVEGFR3 = list(
      description        = "Individual posthoc baseline sVEGFR-3 (pg/mL) from the upstream Hansson 2013a biomarker indirect-response PD fit (DDMODEL00000197). Per-subject, time-fixed; used as the initial condition for the in-model svegfr3 state and as the denominator in the sVEGFR-3 REL driver svegfr3_rel = (BAS_SVEGFR3 - svegfr3) / BAS_SVEGFR3.",
      units              = "pg/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Required input. Hansson 2013 e84 (companion paper) reports a typical sVEGFR-3 baseline of 63900 pg/mL. For new-population simulations either (a) simulate from `Hansson_2013a_sunitinib` to obtain individual posthoc baselines, or (b) set every subject to the typical 63900 pg/mL.",
      source_name        = "BAS3"
    ),
    MRT_SVEGFR3 = list(
      description        = "Individual posthoc mean residence time of sVEGFR-3 (h) from the upstream Hansson 2013a biomarker indirect-response PD fit (DDMODEL00000197). Per-subject, time-fixed; used as kout3 = 1 / MRT_SVEGFR3 inside model().",
      units              = "h",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Required input. The Hansson 2013a typical sVEGFR-3 MRT is 401 h (Hansson 2013 e84 Table 2). Same population strategy as BAS_SVEGFR3.",
      source_name        = "MRT3"
    ),
    EC50_SVEGFR3 = list(
      description        = "Individual posthoc EC50 of the simple-Imax drug effect on sVEGFR-3 (mg*h/L AUC) from the upstream Hansson 2013a biomarker indirect-response PD fit (DDMODEL00000197). Per-subject, time-fixed; appears in the drug-effect term eff_svegfr3 = auc / (EC50_SVEGFR3 + auc).",
      units              = "mg*h/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Required input. The Hansson 2013a typical (shared across the four biomarkers) IC50 is 1.0 mg*h/L (Hansson 2013 e84 Table 2). Same population strategy as BAS_SVEGFR3.",
      source_name        = "EC53"
    ),
    RACE_JAPANESE = list(
      description        = "Binary indicator: 1 = Japanese subject (Study 1045 from Shirao 2010 phase I/II), 0 = non-Japanese (Studies 1004, 1047, 013). Switches the typical baseline ANC0 between 4.94 (non-Japanese) and 3.69 (Japanese) 10^9/L per Hansson 2013 Table 2 row 'ANC0: Study 45'.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Japanese)",
      notes              = "Hansson 2013 estimated a separate baseline parameter for Study 1045 to account for lower ANC levels in Japanese patients (Results section: 'A separate baseline parameter (ANC0) was estimated to account for lower ANC levels in Study 45, which was conducted in Japanese patients.'). For typical-cohort vignette simulations set every subject to 0 (non-Japanese, ANC0 = 4.94).",
      source_name        = "(derived from study identifier; Study 1045 in the paper = Japanese cohort)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 303L,
    n_studies      = 4L,
    age_range      = "adults with imatinib-resistant GIST (paper text reports n=303 pooled across phases I-III; per-cohort baseline-demographics table not in the trimmed PDF section that includes Methods + Results + Tables)",
    weight_range   = "not reported in the on-disk trimmed paper text",
    sex_female_pct = NA_real_,
    race_ethnicity = "majority non-Japanese (Studies 1004, 1047, 013); Japanese subgroup is Study 1045 (n=36) per Shirao 2010",
    disease_state  = "imatinib-resistant gastrointestinal stromal tumours (GIST). Pooled four sunitinib studies: Demetri 2006 (study 1004; placebo-controlled phase III; 202 active + 47 placebo), George 2009 (study 1047; phase II continuous-dosing 37.5 mg QD; n=13 in this analysis subset), Shirao 2010 (study 1045; Japanese phase I/II; 25-75 mg QD on a 4/2 schedule; n=36), Maki 2005 (study 013; phase I/II 25-75 mg QD on a 2/1 or 2/2 schedule; n=52).",
    dose_range     = "sunitinib 25-75 mg PO QD on a 4/2, 2/2, 2/1 (weeks on / weeks off) or continuous treatment schedule (Table 1). The largest cohort (study 1004) used 50 mg QD on a 4/2 schedule. Placebo arm: no sunitinib.",
    regions        = "phase III multinational (study 1004); Japanese phase I/II (study 1045); other studies' regions not stated in the trimmed paper text.",
    biomarkers     = "Absolute neutrophil count (ANC) measured serially across treatment cycles. Median (range) observed ANC during treatment: 3.1 (0.080-20) in study 1004, 1.8 (0.010-7.5) in study 1047, 2.1 (0.28-12) in study 1045, 2.6 (0.16-15) in study 013 (Hansson 2013 Table 1).",
    notes          = "n_subjects = 303 reported in Hansson 2013 e85 Methods ('analyzed data were from four clinical trials in phases I-III, which comprised patients with imatinib-resistant malignant GIST treated with sunitinib... totaling 303 patients'). Per-cohort baseline demographics (age, weight, sex, race) are not transcribed in the trimmed paper text; the cohort breakdown by study is from Table 1."
  )

  ini({
    # ------------------------------------------------------------------
    # Structural PD parameters from Hansson 2013 e85 Table 2 'Myelosuppression
    # model' block. The paper Table 2 reports point estimates with relative
    # standard error (RSE %) for the typical population; IIV is reported as
    # CV% which back-transforms to omega^2 = log(CV^2 + 1).
    # ------------------------------------------------------------------

    # Baseline ANC for non-Japanese subjects (Studies 1004, 1047, 013).
    lanc0      <- log(4.94);  label("Baseline ANC for non-Japanese subjects (10^9/L; Hansson 2013 Table 2 row 'ANC0')") # Table 2 ANC0 = 4.94 (RSE 2.8%)

    # Multiplicative shift for Japanese subjects (Study 1045): ANC0_japanese / ANC0_other = 3.69 / 4.94 = 0.747.
    # Encoded as an exponential effect on log(ANC0): log(3.69/4.94) = -0.292.
    e_japanese_anc0 <- log(3.69 / 4.94); label("Multiplicative shift on log(ANC0) for Japanese subjects (Study 1045) yielding ANC0_japanese = 3.69 10^9/L (Hansson 2013 Table 2 row 'ANC0: Study 45')") # Table 2 ANC0 Study 45 = 3.69 (RSE 6.9%)

    # Mean transit time through non-sensitive compartments (proliferation + 3 transits).
    lmtt       <- log(248);   label("Mean transit time MTT through the proliferation -> transit chain (h; Hansson 2013 Table 2)") # Table 2 MTT = 248 h (RSE 3.6%)

    # Maximum effect of sVEGFR-3 REL on cell-loss / proliferation.
    lanc_emax  <- log(0.520); label("Maximum drug effect Emax on neutrophil proliferation (unitless fraction; Hansson 2013 Table 2 row 'ANC Emax')") # Table 2 ANC Emax = 0.520 (RSE 9.1%)

    # EC50 for sVEGFR-3 REL effect. Paper Table 2 reports units 'pg.hour/l' in
    # the row label but the driver of the Emax function is the unitless
    # relative change in sVEGFR-3 from baseline (Methods: 'longitudinal
    # model-predicted relative change in sVEGFR-3 from baseline (sVEGFR-3 REL)
    # was the better descriptor'). The 0.552 value is interpreted as the
    # unitless sVEGFR-3 REL value at which the drug effect is half-maximal
    # (i.e., when sVEGFR-3 has been depressed by 55.2% below baseline). The
    # 'pg.hour/l' units in the Table 2 row label are an editing artifact
    # carried over from competing AUC-driven model rows -- see vignette
    # Assumptions and deviations.
    lanc_ec50  <- log(0.552); label("EC50 for sVEGFR-3 REL drug effect (unitless fractional reduction in sVEGFR-3; Hansson 2013 Table 2 row 'ANC EC50')") # Table 2 ANC EC50 = 0.552 (RSE 17%)

    # Feedback factor on (ANC0 / circ).
    gamma      <- 0.362;      label("Feedback exponent gamma on (ANC0 / circ) (unitless; Hansson 2013 Table 2 row 'gamma')") # Table 2 gamma = 0.362 (RSE 7.4%)

    # ------------------------------------------------------------------
    # Inter-individual variability. Table 2 reports IIV as CV%; conversion
    # to log-normal omega^2 = log((CV/100)^2 + 1).
    # ANC0 IIV = 42%  -> omega2 = log(0.42^2 + 1) = 0.1633
    # MTT  IIV = 17%  -> omega2 = log(0.17^2 + 1) = 0.02852
    # Emax IIV = 13%  -> omega2 = log(0.13^2 + 1) = 0.01678
    # EC50 IIV = 46%  -> omega2 = log(0.46^2 + 1) = 0.1929
    # ANC0 and Emax are reported as correlated 90% (Results section).
    # cov(ANC0, Emax) = corr * sqrt(var_ANC0 * var_Emax)
    #                 = 0.90 * sqrt(0.1633 * 0.01678) = 0.04707
    # ------------------------------------------------------------------
    etalanc0 + etalanc_emax ~ c(0.1633,
                                0.04707, 0.01678) # Table 2 ANC0 CV=42%, ANC Emax CV=13%, corr=0.90
    etalmtt       ~ 0.02852  # Table 2 MTT CV=17%
    etalanc_ec50  ~ 0.1929   # Table 2 ANC EC50 CV=46%

    # ------------------------------------------------------------------
    # Residual error. Hansson 2013 Methods: 'Residual variability was
    # described by an additive (on Box-Cox scale) error model' with the
    # Box-Cox lambda fixed at 0.2. Table 2 reports residual error 0.406
    # on the Box-Cox-transformed scale. nlmixr2 / rxode2 do not yet
    # natively express the Box-Cox-residual transformation as a one-line
    # error block, so this model implements a simple additive residual on
    # the linear scale as a forward-simulation surrogate; the published
    # 0.406 value is preserved as the additive SD with a vignette
    # Assumptions note. See model() for the alternative linear-scale
    # tolerance interpretation.
    # ------------------------------------------------------------------
    addSd_anc <- 0.406; label("Residual error magnitude on the source Box-Cox-transformed ANC scale (lambda = 0.2; Hansson 2013 Table 2 row 'Residual error'). Used as a linear-scale additive SD in this nlmixr2 forward-simulation port; see vignette Assumptions for the deviation.") # Table 2 residual error = 0.406 (RSE 4.3%)
  })

  model({
    # ---- 1. Per-cycle drug-exposure summary (mg*h/L AUC) ----
    # mg / (L/h) = mg*h/L; matches the Hansson_2013c_sunitinib convention.
    auc <- DOSE / CLI

    # ---- 2. Upstream sVEGFR-3 turnover (simple-Imax inhibition of Kin) ----
    # Reproduces the Hansson 2013a typical-value sVEGFR-3 dynamics using the
    # per-subject upstream-PD covariates BAS_SVEGFR3, MRT_SVEGFR3, EC50_SVEGFR3.
    kout_svegfr3 <- 1 / MRT_SVEGFR3
    kin_svegfr3  <- BAS_SVEGFR3 * kout_svegfr3
    eff_svegfr3  <- auc / (EC50_SVEGFR3 + auc)
    svegfr3(0)     <- BAS_SVEGFR3
    d/dt(svegfr3)  <- kin_svegfr3 * (1 - eff_svegfr3) - kout_svegfr3 * svegfr3

    # ---- 3. Relative change in sVEGFR-3 from baseline (positive when drug-treated) ----
    # Hansson 2013 e85 Methods: 'Increasing sVEGFR-3 REL, i.e., a more
    # pronounced reduction in sVEGFR-3, was associated with increased
    # probability and severity ...'. The driver is therefore the magnitude
    # of the reduction below baseline; encoded with the sign such that
    # drug-driven sVEGFR-3 depletion gives positive svegfr3_rel.
    svegfr3_rel <- (BAS_SVEGFR3 - svegfr3) / BAS_SVEGFR3

    # ---- 4. Individual ANC structural parameters ----
    # Baseline ANC switches by RACE_JAPANESE (multiplicative log shift).
    anc0     <- exp(lanc0 + e_japanese_anc0 * RACE_JAPANESE + etalanc0)
    mtt      <- exp(lmtt  + etalmtt)
    anc_emax <- exp(lanc_emax + etalanc_emax)
    anc_ec50 <- exp(lanc_ec50 + etalanc_ec50)

    # ---- 5. Drug-effect Emax on the proliferation rate ----
    edrug <- anc_emax * svegfr3_rel / (anc_ec50 + svegfr3_rel)
    feed  <- (anc0 / circ)^gamma

    # ---- 6. Friberg myelosuppression chain ----
    # ktr = (n_trans + 1) / MTT with n_trans = 3 (paper Methods: 'three
    # transit compartments reflecting cell maturation'); matches the
    # Friberg 2002 / Wahlby 2004 / Friberg_2002_paclitaxel.R precedent.
    # Hansson 2013 fixed the half-life of circulating neutrophils to 7 h
    # 'to enhance the physiological interpretation of the model'; in this
    # forward-simulation port the circulating-pool elimination shares ktr
    # with the rest of the chain, which is the standard Friberg form. The
    # MTT = 248 h value is the estimated maturation time and is preserved
    # as the source-published value. See vignette Assumptions and
    # deviations for the implication for circulating-pool half-life.
    ktr <- 4 / mtt

    d/dt(prol)     <- ktr * prol * (1 - edrug) * feed - ktr * prol
    d/dt(transit1) <- ktr * prol     - ktr * transit1
    d/dt(transit2) <- ktr * transit1 - ktr * transit2
    d/dt(transit3) <- ktr * transit2 - ktr * transit3
    d/dt(circ)     <- ktr * transit3 - ktr * circ

    prol(0)     <- anc0
    transit1(0) <- anc0
    transit2(0) <- anc0
    transit3(0) <- anc0
    circ(0)     <- anc0

    # ---- 7. Observation ----
    # Observation is ANC (= circulating-neutrophil pool) in 10^9/L. See
    # vignette Assumptions for the Box-Cox-residual deviation.
    anc <- circ
    anc ~ add(addSd_anc)
  })
}
