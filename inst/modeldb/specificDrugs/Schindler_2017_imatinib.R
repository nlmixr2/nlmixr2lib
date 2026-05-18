Schindler_2017_imatinib <- function() {
  description <- paste(
    "Joint tumor-dynamics PD model for imatinib-treated GIST liver metastases",
    "(Schindler 2017). Three size metrics (maximum transaxial diameter MTD in",
    "mm, software-segmented actual volume Vactual in mL, calculated ellipsoidal",
    "volume Vellipsoid in mL) follow a logistic tumor-growth model with a",
    "linear DOSE-dependent shrinkage term and a mono-exponential drug-effect",
    "washout (resistance development). Tumor density (Hounsfield units) follows",
    "an indirect-response model in which imatinib linearly stimulates the loss",
    "rate. Each subject can carry up to two liver lesions (lesion 1 has the",
    "larger baseline by convention); the binary covariate MIX_LARGE_BASE",
    "selects between a mixture subpopulation with larger lesion baselines",
    "(MIX_LARGE_BASE = 1, P = 0.348) and a smaller-baseline subpopulation",
    "(MIX_LARGE_BASE = 0). Drug exposure enters via the daily dose normalized",
    "to the median 400 mg, so DOSE is supplied as a per-record time-varying",
    "covariate (in mg/day). The OS and PFS time-to-event arms of the source",
    "publication are not encoded as ODE compartments here (see vignette",
    "Assumptions and deviations)."
  )
  reference <- paste(
    "Schindler E, Krishnan SM, Mathijssen RHJ, Ruggiero A, Schiavon G,",
    "Friberg LE. (2017). Pharmacometric modeling of liver metastases'",
    "diameter, volume, and density and their relation to clinical outcome",
    "in imatinib-treated patients with gastrointestinal stromal tumors.",
    "CPT Pharmacometrics Syst Pharmacol 6(7):449-457.",
    "doi:10.1002/psp4.12195."
  )
  vignette <- "Schindler_2017_imatinib"
  units <- list(
    time          = "week",
    dosing        = "mg/day",
    concentration = "n/a (non-PK outputs: MTD in mm, Vactual / Vellipsoid in mL, density in HU)"
  )

  covariateData <- list(
    DOSE = list(
      description        = paste(
        "Current daily oral imatinib dose level (in mg/day) supplied as a",
        "per-record time-varying covariate. Enters the size and density",
        "models through DOSE / 400 (the median starting dose), so the",
        "drug-effect term Kdrug,S * (DOSE / 400) * exp(-k * t) scales",
        "linearly with dose intensity."
      ),
      units              = "mg/day",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Schindler 2017 Methods 'Maximum transaxial diameter ...' subsection",
        "and Results paragraph after Eq. 2 specify the dose normalisation:",
        "'a linear effect with slope Kdrug,S driven by imatinib daily dose,",
        "normalized by the median dose (400 mg).' Reference daily dose is",
        "400 mg (74 of 77 patients started at 400 mg; 3 started at 800 mg).",
        "30 of 77 (39 %) underwent dose escalation to a median of 800 mg",
        "(range 600-1200) and 8 (10 %) had dose reductions (median 300 mg,",
        "range 200-300) per Table 1. Set DOSE to 0 during any off-treatment",
        "/ drug-holiday record so the drug-effect term vanishes; for a",
        "typical patient on 400 mg/day continuously, DOSE = 400 throughout."
      ),
      source_name        = "Dose"
    ),
    MIX_LARGE_BASE = list(
      description        = paste(
        "Per-subject binary mixture-model class indicator. 1 = subject",
        "classified to the subpopulation with larger baseline lesions",
        "(typically 76.6 mm MTD / 161 mL Vactual / 187 mL Vellipsoid for",
        "lesion 1, per Table 2); 0 = subject classified to the",
        "smaller-baseline subpopulation (20.9 mm / 3.45 mL / 3.93 mL for",
        "lesion 1)."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (smaller-baseline subpopulation)",
      notes              = paste(
        "Not a measured patient covariate. Per-subject latent class",
        "assignment from a NONMEM $MIXTURE block (Schindler 2017 Methods,",
        "Maximum transaxial diameter / actual volumes / ellipsoidal volume",
        "subsection: 'Semiparametric distributions and mixture models were",
        "investigated to describe the observed bimodal distribution in",
        "baseline MTD, Vactual, and Vellipsoid'). Population probability of",
        "MIX_LARGE_BASE = 1 is the estimated mixture fraction Ppop1 = 0.348",
        "(Table 2; common to MTD, Vactual, and Vellipsoid models; 95 % CI",
        "not reported). Of patients with one target lesion 22 % were in",
        "subpopulation 1 vs 39 % of those with two target lesions. For",
        "typical-value simulation set MIX_LARGE_BASE = 1 (larger-baseline",
        "phenotype; mean simulations dominated by the larger lesion volume)",
        "or 0 (smaller-baseline phenotype). For population simulation, draw",
        "MIX_LARGE_BASE ~ Bernoulli(0.348) per subject. The reference",
        "category (= 0 = smaller baseline) is chosen so the binary",
        "numerically matches the paper's mixture indicator orientation."
      ),
      source_name        = "MIXTURE (NONMEM $MIXTURE assignment)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 77L,
    n_studies      = 2L,
    age_range      = "34-83 years (median 62)",
    weight_range   = NA_character_,
    sex_female_pct = 39.0,
    race_ethnicity = NA_character_,
    disease_state  = paste(
      "Adults with gastrointestinal stromal tumor (GIST) and at least one",
      "liver metastasis assessed by computed tomography. 77 patients total",
      "with up to two liver target lesions per patient (60 of 77 with two",
      "lesions, 17 with one). Lesions were numbered so that lesion 1 had the",
      "largest baseline maximum transaxial diameter (MTD). 35 / 77 patients",
      "(45 %) received second-line therapy and beyond; 30 / 77 (39 %) had",
      "dose escalation and 8 / 77 (10 %) had dose reductions during",
      "follow-up. Tumor observations: 502 measurements per size metric and",
      "496 density measurements collected from 136 lesions, with median",
      "follow-up of 360 days (range 82-495). Overall survival followed for",
      "a median of 4.5 years (range 0.79-13) and progression-free survival",
      "for a median of 3.4 years (range 0.25-13); 43 / 77 deaths and",
      "50 / 77 progression events were recorded as of 17 February 2015."
    ),
    dose_range     = paste(
      "Oral imatinib first-line therapy at starting dose 400 mg / day",
      "(n = 74) or 800 mg / day (n = 3); dose escalation to a median 800",
      "(range 600-1200) mg / day in 30 patients (39 %); dose reductions to",
      "a median 300 (range 200-300) mg / day in 8 patients (10 %)."
    ),
    regions        = NA_character_,
    notes          = paste(
      "Demographics from Schindler 2017 Table 1. The cohort pools two",
      "retrospective, non-interventional datasets previously published by",
      "Schiavon et al. (2014, Eur J Cancer 50:972-980) and used as part of",
      "Eechoute et al. (2012, Clin Cancer Res 18:5780-5787); ethical",
      "standards followed as described in Schiavon 2014. Information on",
      "subsequent systemic therapy (sunitinib, regorafenib) was recorded",
      "but is not part of the modelled exposure-response."
    )
  )

  ini({
    # =====================================================================
    # MTD (maximum transaxial diameter, mm) -- Schindler 2017 Table 2
    # =====================================================================
    # Mixture subpopulation 1 (larger baseline; Ppop1 = 0.348) typical S0
    lS0_mtd_pop1_l1 <- log(76.6);   label("MTD baseline S0 -- lesion 1, mixture subpop 1 (larger, mm)")  # Table 2 row 'S0, pop1' lesion 1
    lS0_mtd_pop1_l2 <- log(41.9);   label("MTD baseline S0 -- lesion 2, mixture subpop 1 (larger, mm)")  # Table 2 row 'S0, pop1' lesion 2

    # Mixture subpopulation 2 (smaller baseline) typical S0
    lS0_mtd_pop2_l1 <- log(20.9);   label("MTD baseline S0 -- lesion 1, mixture subpop 2 (smaller, mm)") # Table 2 row 'S0, pop2' lesion 1
    lS0_mtd_pop2_l2 <- log(14.2);   label("MTD baseline S0 -- lesion 2, mixture subpop 2 (smaller, mm)") # Table 2 row 'S0, pop2' lesion 2

    # Per-lesion carrying capacity (shared across subpopulations)
    lSmax_mtd_l1    <- log(171);    label("MTD carrying capacity Smax -- lesion 1 (mm)")                  # Table 2 row 'Smax' lesion 1
    lSmax_mtd_l2    <- log(125);    label("MTD carrying capacity Smax -- lesion 2 (mm)")                  # Table 2 row 'Smax' lesion 2

    # Tumor-growth, resistance and drug-effect rates
    lKG_mtd         <- log(0.00176); label("MTD logistic tumor growth rate KG (1/week)")                  # Table 2 row 'KG'
    lk_mtd          <- log(0.0475);  label("MTD drug-effect resistance decay rate k (1/week)")            # Table 2 row 'k'
    lKdrug_mtd      <- log(0.0124);  label("MTD drug-effect slope Kdrug,S per (DOSE / 400) (1/week)")     # Table 2 row 'Kdrug,S'

    # =====================================================================
    # Vactual (software-segmented actual volume, mL) -- Schindler 2017 Table 2
    # =====================================================================
    lS0_vact_pop1_l1 <- log(161);     label("Vactual S0 -- lesion 1, mixture subpop 1 (mL)")              # Table 2 row 'S0, pop1' Vactual lesion 1
    lS0_vact_pop1_l2 <- log(29.7);    label("Vactual S0 -- lesion 2, mixture subpop 1 (mL)")              # Table 2 row 'S0, pop1' Vactual lesion 2
    lS0_vact_pop2_l1 <- log(3.45);    label("Vactual S0 -- lesion 1, mixture subpop 2 (mL)")              # Table 2 row 'S0, pop2' Vactual lesion 1
    lS0_vact_pop2_l2 <- log(1.21);    label("Vactual S0 -- lesion 2, mixture subpop 2 (mL)")              # Table 2 row 'S0, pop2' Vactual lesion 2
    lSmax_vact_l1    <- log(1190);    label("Vactual carrying capacity Smax -- lesion 1 (mL)")            # Table 2 row 'Smax' Vactual lesion 1
    lSmax_vact_l2    <- log(540);     label("Vactual carrying capacity Smax -- lesion 2 (mL)")            # Table 2 row 'Smax' Vactual lesion 2
    lKG_vact         <- log(0.00861); label("Vactual logistic tumor growth rate KG (1/week)")             # Table 2 row 'KG' Vactual
    lk_vact          <- log(0.0469);  label("Vactual drug-effect resistance decay rate k (1/week)")       # Table 2 row 'k' Vactual
    lKdrug_vact      <- log(0.0547);  label("Vactual drug-effect slope Kdrug,S per (DOSE / 400) (1/week)") # Table 2 row 'Kdrug,S' Vactual

    # =====================================================================
    # Vellipsoid (calculated ellipsoidal volume, mL) -- Schindler 2017 Table 2
    # =====================================================================
    # Per Table 2 footnote d: Smax lesion 1 95 % LL profile CI = (685, 3260) mL.
    lS0_vell_pop1_l1 <- log(187);     label("Vellipsoid S0 -- lesion 1, mixture subpop 1 (mL)")            # Table 2 row 'S0, pop1' Vellipsoid lesion 1
    lS0_vell_pop1_l2 <- log(33.4);    label("Vellipsoid S0 -- lesion 2, mixture subpop 1 (mL)")            # Table 2 row 'S0, pop1' Vellipsoid lesion 2
    lS0_vell_pop2_l1 <- log(3.93);    label("Vellipsoid S0 -- lesion 1, mixture subpop 2 (mL)")            # Table 2 row 'S0, pop2' Vellipsoid lesion 1
    lS0_vell_pop2_l2 <- log(1.27);    label("Vellipsoid S0 -- lesion 2, mixture subpop 2 (mL)")            # Table 2 row 'S0, pop2' Vellipsoid lesion 2
    lSmax_vell_l1    <- log(1230);    label("Vellipsoid carrying capacity Smax -- lesion 1 (mL)")          # Table 2 row 'Smax' Vellipsoid lesion 1 (95 % LL profile CI 685-3260)
    lSmax_vell_l2    <- log(588);     label("Vellipsoid carrying capacity Smax -- lesion 2 (mL)")          # Table 2 row 'Smax' Vellipsoid lesion 2
    lKG_vell         <- log(0.00882); label("Vellipsoid logistic tumor growth rate KG (1/week)")           # Table 2 row 'KG' Vellipsoid
    lk_vell          <- log(0.0508);  label("Vellipsoid drug-effect resistance decay rate k (1/week)")     # Table 2 row 'k' Vellipsoid
    lKdrug_vell      <- log(0.0610);  label("Vellipsoid drug-effect slope Kdrug,S per (DOSE / 400) (1/week)") # Table 2 row 'Kdrug,S' Vellipsoid

    # Mixture probability Ppop1 = 0.348 (common to MTD, Vactual, Vellipsoid;
    # Table 2 row 'Ppop1' with 18 % RSE) is recorded in
    # covariateData[[MIX_LARGE_BASE]]$notes rather than as a model parameter
    # because the mixture class is supplied as the per-subject covariate
    # MIX_LARGE_BASE; users draw it externally via
    # MIX_LARGE_BASE ~ Bernoulli(0.348) per Schindler 2017 Table 2.

    # =====================================================================
    # Density (Hounsfield units) -- Schindler 2017 Table 3
    # =====================================================================
    # Box-Cox shape parameter on the combined IIV + ILV random effects of D0
    # is reported as -1.06 with 95 % LL profile CI (-2.01, -0.397) in
    # Table 3 footnote a. Not represented natively here -- the per-subject
    # and per-lesion deviations on D0 are encoded as standard log-normal
    # etas because nlmixr2 does not have a native Box-Cox transform on
    # random-effect realizations. See the vignette's Assumptions and
    # deviations section for the deviation note.
    lD0              <- log(59.0);    label("Density baseline D0 (HU)")                                    # Table 3 row 'D0'
    lkout_dens       <- log(0.0935);  label("Density loss rate kout (1/week); MRT = 1 / kout ~ 75 days")   # Table 3 row 'kout' (MRT = 75 d confirmed in Results)
    lKdrug_dens      <- log(0.154);   label("Density drug-effect slope Kdrug,D per (DOSE / 400) (unitless multiplier on kout)") # Table 3 row 'Kdrug,D'

    # =====================================================================
    # Residual unexplained variability (LTBS in source -> proportional)
    # =====================================================================
    # NONMEM "additive on log-scale" residual model maps to nlmixr2's prop()
    # (proportional in linear space) per Methods 'Residual unexplained
    # variability was described by additive models on a logarithmic scale.'
    # RUV % from Tables 2 and 3 are CV %; entered directly as SD on the
    # fractional proportional scale. The source reports one RUV per size
    # metric (shared across the two lesions); nlmixr2 requires one
    # residual-error parameter per observed compartment, so the per-lesion
    # propSd values are declared separately with identical numeric values.
    propSd_mtd_l1   <- 0.140; label("MTD lesion 1 proportional residual SD (fraction)")        # Table 2 row 'RUV (%)' MTD = 14.0 % (shared across lesions)
    propSd_mtd_l2   <- 0.140; label("MTD lesion 2 proportional residual SD (fraction)")        # Table 2 row 'RUV (%)' MTD = 14.0 % (shared across lesions)
    propSd_vact_l1  <- 0.368; label("Vactual lesion 1 proportional residual SD (fraction)")    # Table 2 row 'RUV (%)' Vactual = 36.8 % (shared across lesions)
    propSd_vact_l2  <- 0.368; label("Vactual lesion 2 proportional residual SD (fraction)")    # Table 2 row 'RUV (%)' Vactual = 36.8 % (shared across lesions)
    propSd_vell_l1  <- 0.433; label("Vellipsoid lesion 1 proportional residual SD (fraction)") # Table 2 row 'RUV (%)' Vellipsoid = 43.3 % (shared across lesions)
    propSd_vell_l2  <- 0.433; label("Vellipsoid lesion 2 proportional residual SD (fraction)") # Table 2 row 'RUV (%)' Vellipsoid = 43.3 % (shared across lesions)
    propSd_dens_l1  <- 0.206; label("Density lesion 1 proportional residual SD (fraction)")    # Table 3 row 'RUV (%)' = 20.6 % (shared across lesions)
    propSd_dens_l2  <- 0.206; label("Density lesion 2 proportional residual SD (fraction)")    # Table 3 row 'RUV (%)' = 20.6 % (shared across lesions)

    # =====================================================================
    # Inter-individual variability (IIV)
    # =====================================================================
    # Variance on the log-eta scale is computed from the CV % entries of
    # Tables 2 and 3 via omega^2 = log(1 + (CV/100)^2). Per Table 2 footnote
    # b, the S0 IIV variance is shared between lesions within each
    # subpopulation and model, so a single eta per (metric x subpop) is
    # multiplied onto both lesion 1 and lesion 2 inside model().

    # --- MTD ---
    etalS0_mtd_pop1 ~ 0.1996         # Table 2 MTD S0 IIV pop1 = 47 % CV -> log(1 + 0.47^2)
    etalS0_mtd_pop2 ~ 0.0807         # Table 2 MTD S0 IIV pop2 = 29 % CV -> log(1 + 0.29^2)
    etalKG_mtd      ~ 1.3581         # Table 2 MTD KG IIV       = 170 % CV -> log(1 + 1.70^2)
    etalKdrug_mtd   ~ 0.4655         # Table 2 MTD Kdrug,S IIV  = 77 % CV  -> log(1 + 0.77^2)

    # --- Vactual ---
    etalS0_vact_pop1 ~ 1.0852        # Table 2 Vactual S0 IIV pop1 = 140 % CV -> log(1 + 1.40^2)
    etalS0_vact_pop2 ~ 0.4559        # Table 2 Vactual S0 IIV pop2 = 76 % CV  -> log(1 + 0.76^2)
    etalKG_vact      ~ 1.0376        # Table 2 Vactual KG IIV       = 135 % CV -> log(1 + 1.35^2)
    etalKdrug_vact   ~ 0.2728        # Table 2 Vactual Kdrug,S IIV  = 56 % CV  -> log(1 + 0.56^2)

    # --- Vellipsoid ---
    etalS0_vell_pop1 ~ 1.0852        # Table 2 Vellipsoid S0 IIV pop1 = 140 % CV -> log(1 + 1.40^2)
    etalS0_vell_pop2 ~ 0.4753        # Table 2 Vellipsoid S0 IIV pop2 = 78 % CV  -> log(1 + 0.78^2)
    etalKG_vell      ~ 1.0852        # Table 2 Vellipsoid KG IIV      = 140 % CV -> log(1 + 1.40^2)
    etalKdrug_vell   ~ 0.1625        # Table 2 Vellipsoid Kdrug,S IIV = 42 % CV  -> log(1 + 0.42^2)

    # --- Density ---
    etalD0           ~ 0.0862        # Table 3 D0 IIV       = 30 % CV  -> log(1 + 0.30^2)
    etalKdrug_dens   ~ 0.8920        # Table 3 Kdrug,D IIV  = 120 % CV -> log(1 + 1.20^2)

    # =====================================================================
    # Inter-lesion variability (ILV) on density parameters
    # =====================================================================
    # The two ILV etas per parameter share a common variance (Methods Eq. 1
    # and Schindler 2017 Results 'IIV magnitude ... larger than the ILV
    # magnitude (120 vs 53 % CV)'). Each lesion draws its own per-lesion eta
    # from the same variance; no cross-lesion covariance.
    etalD0_les1      ~ fixed(0.03189)  # Table 3 D0 ILV       = 18 % CV  -> log(1 + 0.18^2)
    etalD0_les2      ~ fixed(0.03189)  # Table 3 D0 ILV (lesion 2, shared variance)
    etalKdrug_dens_les1 ~ fixed(0.2474) # Table 3 Kdrug,D ILV  = 53 % CV  -> log(1 + 0.53^2)
    etalKdrug_dens_les2 ~ fixed(0.2474) # Table 3 Kdrug,D ILV (lesion 2, shared variance)
  })

  model({
    # ============================================================
    # 1. Derived mixture indicators
    # ============================================================
    # MIX_LARGE_BASE is the per-subject Bernoulli draw selecting between
    # mixture subpopulations; mix1 / mix2 are the two-class indicators that
    # gate the lesion-baseline typical values and the subpop-specific IIV
    # etas below. Per Methods, this latent class is independent of size
    # metric (one assignment per subject), so the same mix1 / mix2 are
    # reused across MTD, Vactual, and Vellipsoid.
    mix1 <- MIX_LARGE_BASE
    mix2 <- 1 - mix1

    # Normalized daily-dose driver (Methods 'a linear effect ... driven by
    # imatinib daily dose, normalized by the median dose (400 mg)').
    ddose <- DOSE / 400

    # ============================================================
    # 2. Per-subject IIV etas selected by subpopulation
    # ============================================================
    # The subpopulation-specific IIV etas combine into a single per-subject
    # additive log-deviation; only the eta whose variance matches the
    # subject's mixture class contributes. This reproduces the source
    # model's per-class OMEGA structure for S0 while keeping a single
    # composite eta per (metric x subject) downstream.
    eta_S0_mtd  <- etalS0_mtd_pop1  * mix1 + etalS0_mtd_pop2  * mix2
    eta_S0_vact <- etalS0_vact_pop1 * mix1 + etalS0_vact_pop2 * mix2
    eta_S0_vell <- etalS0_vell_pop1 * mix1 + etalS0_vell_pop2 * mix2

    # ============================================================
    # 3. Typical baselines per metric per lesion (mixture-selected)
    # ============================================================
    S0_mtd_l1_typ <- exp(lS0_mtd_pop1_l1) * mix1 + exp(lS0_mtd_pop2_l1) * mix2
    S0_mtd_l2_typ <- exp(lS0_mtd_pop1_l2) * mix1 + exp(lS0_mtd_pop2_l2) * mix2

    S0_vact_l1_typ <- exp(lS0_vact_pop1_l1) * mix1 + exp(lS0_vact_pop2_l1) * mix2
    S0_vact_l2_typ <- exp(lS0_vact_pop1_l2) * mix1 + exp(lS0_vact_pop2_l2) * mix2

    S0_vell_l1_typ <- exp(lS0_vell_pop1_l1) * mix1 + exp(lS0_vell_pop2_l1) * mix2
    S0_vell_l2_typ <- exp(lS0_vell_pop1_l2) * mix1 + exp(lS0_vell_pop2_l2) * mix2

    # ============================================================
    # 4. Individual parameters (with IIV applied)
    # ============================================================
    # MTD
    S0_mtd_l1   <- S0_mtd_l1_typ * exp(eta_S0_mtd)
    S0_mtd_l2   <- S0_mtd_l2_typ * exp(eta_S0_mtd)
    Smax_mtd_l1 <- exp(lSmax_mtd_l1)
    Smax_mtd_l2 <- exp(lSmax_mtd_l2)
    KG_mtd      <- exp(lKG_mtd + etalKG_mtd)
    k_mtd       <- exp(lk_mtd)
    Kdrug_mtd   <- exp(lKdrug_mtd + etalKdrug_mtd)

    # Vactual
    S0_vact_l1   <- S0_vact_l1_typ * exp(eta_S0_vact)
    S0_vact_l2   <- S0_vact_l2_typ * exp(eta_S0_vact)
    Smax_vact_l1 <- exp(lSmax_vact_l1)
    Smax_vact_l2 <- exp(lSmax_vact_l2)
    KG_vact      <- exp(lKG_vact + etalKG_vact)
    k_vact       <- exp(lk_vact)
    Kdrug_vact   <- exp(lKdrug_vact + etalKdrug_vact)

    # Vellipsoid
    S0_vell_l1   <- S0_vell_l1_typ * exp(eta_S0_vell)
    S0_vell_l2   <- S0_vell_l2_typ * exp(eta_S0_vell)
    Smax_vell_l1 <- exp(lSmax_vell_l1)
    Smax_vell_l2 <- exp(lSmax_vell_l2)
    KG_vell      <- exp(lKG_vell + etalKG_vell)
    k_vell       <- exp(lk_vell)
    Kdrug_vell   <- exp(lKdrug_vell + etalKdrug_vell)

    # Density (indirect response). The source applies a Box-Cox transform
    # of shape -1.06 to the combined IIV + ILV deviation on D0 to handle a
    # skewed random-effects distribution; here the deviation is encoded as
    # a standard log-normal because nlmixr2 does not natively support
    # Box-Cox transformed random effects (see vignette Assumptions and
    # deviations). The composite (theta + eta_IIV + eta_ILV) form cannot
    # be expressed on a single mu-referenced line in nlmixr2; the
    # subject-typical IIV value is built first, then each lesion
    # multiplies by its own per-lesion eta draw.
    D0_typ        <- exp(lD0 + etalD0)
    Kdrug_dens_typ <- exp(lKdrug_dens + etalKdrug_dens)
    D0_l1        <- D0_typ * exp(etalD0_les1)
    D0_l2        <- D0_typ * exp(etalD0_les2)
    kout_dens    <- exp(lkout_dens)
    Kdrug_dens_l1 <- Kdrug_dens_typ * exp(etalKdrug_dens_les1)
    Kdrug_dens_l2 <- Kdrug_dens_typ * exp(etalKdrug_dens_les2)

    # ============================================================
    # 5. Initial conditions
    # ============================================================
    mtd_l1(0)      <- S0_mtd_l1
    mtd_l2(0)      <- S0_mtd_l2
    vactual_l1(0)  <- S0_vact_l1
    vactual_l2(0)  <- S0_vact_l2
    vellipsoid_l1(0) <- S0_vell_l1
    vellipsoid_l2(0) <- S0_vell_l2
    density_l1(0)  <- D0_l1
    density_l2(0)  <- D0_l2

    # ============================================================
    # 6. ODE system
    # ============================================================
    # Logistic tumor growth with linear drug-effect shrinkage and
    # exponential resistance development (Schindler 2017 Eq. 2). Each
    # lesion of each size metric follows the same equation with its own
    # carrying capacity and per-subject scaled baseline.
    d/dt(mtd_l1) <- KG_mtd  * mtd_l1 * (1 - mtd_l1 / Smax_mtd_l1) -
                    Kdrug_mtd * ddose * exp(-k_mtd * t) * mtd_l1
    d/dt(mtd_l2) <- KG_mtd  * mtd_l2 * (1 - mtd_l2 / Smax_mtd_l2) -
                    Kdrug_mtd * ddose * exp(-k_mtd * t) * mtd_l2

    d/dt(vactual_l1) <- KG_vact * vactual_l1 * (1 - vactual_l1 / Smax_vact_l1) -
                        Kdrug_vact * ddose * exp(-k_vact * t) * vactual_l1
    d/dt(vactual_l2) <- KG_vact * vactual_l2 * (1 - vactual_l2 / Smax_vact_l2) -
                        Kdrug_vact * ddose * exp(-k_vact * t) * vactual_l2

    d/dt(vellipsoid_l1) <- KG_vell * vellipsoid_l1 * (1 - vellipsoid_l1 / Smax_vell_l1) -
                           Kdrug_vell * ddose * exp(-k_vell * t) * vellipsoid_l1
    d/dt(vellipsoid_l2) <- KG_vell * vellipsoid_l2 * (1 - vellipsoid_l2 / Smax_vell_l2) -
                           Kdrug_vell * ddose * exp(-k_vell * t) * vellipsoid_l2

    # Indirect response model with stimulation of output (Schindler 2017
    # Eq. 3). Rin is parameterized as kout * D0 so the state is at
    # steady-state D0 with no drug (Methods 'The rate constant for
    # production of response was parameterized as Rin = kout * D0').
    d/dt(density_l1) <- kout_dens * D0_l1 -
                        kout_dens * (1 + Kdrug_dens_l1 * ddose) * density_l1
    d/dt(density_l2) <- kout_dens * D0_l2 -
                        kout_dens * (1 + Kdrug_dens_l2 * ddose) * density_l2

    # ============================================================
    # 7. Observation models (proportional residual error per metric;
    #    shared across the subject's two lesions, per the source)
    # ============================================================
    mtd_l1        ~ prop(propSd_mtd_l1)
    mtd_l2        ~ prop(propSd_mtd_l2)
    vactual_l1    ~ prop(propSd_vact_l1)
    vactual_l2    ~ prop(propSd_vact_l2)
    vellipsoid_l1 ~ prop(propSd_vell_l1)
    vellipsoid_l2 ~ prop(propSd_vell_l2)
    density_l1    ~ prop(propSd_dens_l1)
    density_l2    ~ prop(propSd_dens_l2)
  })
}
