Mann_2022_fentanyl_iv <- function() {
  description <- paste(
    "Three-compartment IV fentanyl population PK with a first-order",
    "biophase (effect-site) equilibrium compartment, used as the agonist",
    "input layer of the Mann 2022 translational opioid-overdose model.",
    "Parameter values are the Algera 2021 popPK fit re-tabulated in Mann",
    "2022 Supplement 1 Table S1 (intravenous fentanyl, healthy opioid-",
    "naive and chronic opioid-user volunteers pooled, n = 30). Allometric",
    "scaling: CL and inter-compartmental clearances on (WT/70)^0.75,",
    "volumes on (WT/70). Outputs plasma concentration Cc in ng/mL and",
    "effect-site Ce in both ng/mL and pM for downstream consumption by",
    "the Mann 2022 mu-opioid receptor binding model. Intended for use as",
    "the IV-fentanyl agonist input in a simulated overdose-rescue chain;",
    "no residual error is reported in the source supplement."
  )
  reference <- paste(
    "Mann J, Samieegohar M, Chaturbedi A, Zirkle J, Han X, Ahmadi SF,",
    "Eshleman A, Janowsky A, Wolfrum K, Swanson T, Bloom S, Dahan A,",
    "Olofsen E, Florian J, Strauss DG, Li Z.",
    "Development of a Translational Model to Assess the Impact of Opioid",
    "Overdose and Naloxone Dosing on Respiratory Depression and",
    "Cardiac Arrest. Clin Pharmacol Ther. 2022;112(5):1020-1032.",
    "doi:10.1002/cpt.2696. PMID 35766413.",
    "Upstream popPK source for fentanyl IV: Algera MH et al.",
    "Tolerance to opioid-induced respiratory depression in chronic",
    "high-dose opioid users: a model-based comparison with opioid-naive",
    "individuals. Clin Pharmacol Ther. 2021;109(3):637-645.",
    "doi:10.1002/cpt.1972 (re-tabulated in Mann 2022 Table S1).",
    "FDA code repository for the integrated model:",
    "https://github.com/FDA/Mechanistic-PK-PD-Model-to-Rescue-Opioid-Overdose"
  )
  vignette <- "Laffont_2025_opioid_overdose_reversal_simulation"
  units <- list(
    time          = "min",
    dosing        = "mg",
    concentration = "ng/mL"
  )

  covariateData <- list(
    WT = list(
      description        = "Total body weight (Mann 2022 simulations fix WT = 70 kg; the allometric scaling block keeps WT exposed as a covariate so subject-level body weight can vary in downstream composition)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometric a priori scaling per Mann 2022 Supplement 1",
        "Pharmacokinetic Component section: CL = CL_TV * (WT/70)^0.75",
        "and V = V_TV * (WT/70). The inter-compartmental clearances Q2",
        "and Q3 follow the CL exponent (0.75); the peripheral volumes",
        "V2 and V3 follow the V exponent (1). Reference weight 70 kg.",
        "The Mann 2022 virtual-population simulations all assume",
        "WT = 70 kg (Supplement 1 Pharmacokinetic Component, last",
        "paragraph)."
      ),
      source_name        = "weight"
    ),
    Q_TOTAL_LPM = list(
      description        = "Total cardiac output Qb + Qt feeding the FDA delaymymod.c lines 358-368 shock-state Q_Scale feedback that concentrates the opioid in the central / biophase compartment as hyperperfusion drives Q above baseline.",
      units              = "L/min",
      type               = "continuous",
      reference_category = "4.87 (baseline, gives Q_Scale ~ 1, no amplification)",
      notes              = paste(
        "When this PK layer runs standalone (no downstream physiology),",
        "leave Q_TOTAL_LPM at the baseline 4.87 L/min reference so",
        "Q_Scale = 1 + 1 / (1 + exp((1.6 - Q_TOTAL_LPM / 4.87) / 0.05))",
        "evaluates to ~1 and the effective central volume is unchanged.",
        "When this PK layer is composed with",
        "Mann_2022_respiratory_physiology, pass that model's Q_total",
        "(qb + qt) as Q_TOTAL_LPM so that the FDA shock-state feedback",
        "engages during overdose-induced chemoreflex hyperperfusion -",
        "during which Q rises to roughly 1.7-2 x baseline and Q_Scale",
        "saturates at the 2 x cap, effectively halving the apparent",
        "central volume and doubling the effect-site concentration.",
        "Source: FDA delaymymod.c Q_0 = 4.87/60 L/s baseline (line 353)",
        "and Q_Scale sigmoid (lines 358-368)."
      ),
      source_name        = "(none; new canonical covariate registered for this model)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 30L,
    n_studies      = 1L,
    age_range      = "Adult (Algera 2021 cohort: 18-55 years)",
    weight_range   = "Algera 2021 cohort weight range; Mann 2022 virtual-population simulations assume 70 kg",
    sex_female_pct = NA_real_,
    race_ethnicity = NULL,
    disease_state  = paste(
      "Pooled healthy opioid-naive volunteers (n = 14) and chronic",
      "high-dose opioid-user volunteers (n = 16) from Algera 2021;",
      "Mann 2022 uses the pooled-population fit for the IV-fentanyl",
      "PK component."
    ),
    dose_range     = paste(
      "Algera 2021: IV fentanyl 75 or 150 ug/70 kg (opioid-naive)",
      "and 250 or 350 ug/70 kg (chronic opioid users), administered",
      "as 90-second infusions. Mann 2022 simulates IV bolus overdose",
      "scenarios at 1.625 mg and 2.965 mg in chronic users."
    ),
    regions        = NA_character_,
    notes          = paste(
      "PK component imported from Algera 2021 nonlinear mixed-effects",
      "fit without re-estimation in Mann 2022. omega^2 values in",
      "Table S1 are reported on the log scale (variance of the",
      "log-transformed parameter); CV% checks: sqrt(exp(0.087)-1) =",
      "0.30 (30% CV) for CL etc. The biophase equilibrium rate",
      "constant k1 has omega^2 = 0 in Table S1 (no IIV reported).",
      "Effect-site concentration Ce_pM is the downstream input to",
      "Mann_2022_mu_receptor_binding (L_op slot)."
    ),
    upstream_model = "Algera MH et al. Clin Pharmacol Ther. 2021;109(3):637-645 (doi:10.1002/cpt.1972)."
  )

  ini({
    # Structural PK parameters (Mann 2022 Supplement 1 Table S1, intravenous
    # fentanyl row -- values converted from per-second to per-minute by
    # multiplying rate constants by 60).
    lcl  <- log(1.26)
    label("Clearance (L/min)")                                              # Table S1: CL = 0.021 L/s = 1.26 L/min
    lvc  <- log(10.5)
    label("Central volume V1 (L)")                                          # Table S1: V1 = 10.5 L
    lvp  <- log(14.4)
    label("First peripheral volume V2 (L)")                                 # Table S1: V2 = 14.4 L
    lvp2 <- log(166)
    label("Second peripheral volume V3 (L)")                                # Table S1: V3 = 166 L
    lq   <- log(2.04)
    label("Inter-compartmental clearance Q2, central <-> peripheral1 (L/min)") # Table S1: Q2 = 0.034 L/s = 2.04 L/min
    lq2  <- log(2.28)
    label("Inter-compartmental clearance Q3, central <-> peripheral2 (L/min)") # Table S1: Q3 = 0.038 L/s = 2.28 L/min
    lk1  <- fixed(log(0.18))
    label("Biophase equilibration rate k1 (1/min, FIXED, no IIV reported)") # Table S1: k1 = 0.003 1/s = 0.18/min; reported with omega^2 = 0 (no IIV)

    # Fixed allometric exponents (a priori, Algera 2021 / Mann 2022)
    e_wt_cl <- fixed(0.75)
    label("Body-weight allometric exponent on CL and Q (unitless, FIXED)")  # Supplement 1 Pharmacokinetic Component: CL = CL_TV * (WT/70)^0.75
    e_wt_vc  <- fixed(1)
    label("Body-weight allometric exponent on V1, V2, V3 (unitless, FIXED)") # Supplement 1 Pharmacokinetic Component: V = V_TV * (WT/70)

    # IIV -- Table S1 reports log-scale variances omega^2 directly. The
    # parenthesised %CV values are derived from sqrt(exp(omega^2) - 1).
    etalcl  ~ 0.087
    # Table S1: omega^2 CL = 0.087 (30% CV)
    etalvc  ~ 0.43
    # Table S1: omega^2 V1 = 0.43 (73% CV)
    etalvp  ~ 0.52
    # Table S1: omega^2 V2 = 0.52 (82% CV)
    etalvp2 ~ 0.060
    # Table S1: omega^2 V3 = 0.060 (25% CV)
    etalq   ~ 0.39
    # Table S1: omega^2 Q2 = 0.39 (69% CV)
    etalq2  ~ 0.098
    # Table S1: omega^2 Q3 = 0.098 (32% CV)
  })

  model({
    # Allometric weight scaling factor (70 kg reference)
    wt_ratio <- WT / 70

    # Individual structural parameters
    cl  <- exp(lcl  + etalcl)  * wt_ratio^e_wt_cl
    vc  <- exp(lvc  + etalvc)  * wt_ratio^e_wt_vc
    vp  <- exp(lvp  + etalvp)  * wt_ratio^e_wt_vc
    vp2 <- exp(lvp2 + etalvp2) * wt_ratio^e_wt_vc
    q   <- exp(lq   + etalq)   * wt_ratio^e_wt_cl
    q2  <- exp(lq2  + etalq2)  * wt_ratio^e_wt_cl
    k1  <- exp(lk1)

    # Molecular weight constant for plasma -> molar conversion. Fentanyl
    # free base 336.47 g/mol; the FDA reference implementation uses 336.4
    # (delaymymod.c parm Mmass) -- match that to keep the binding-layer
    # input numerically consistent.
    fentanyl_mw_g_mol <- 336.4

    # FDA delaymymod.c shock-state PK amplification (lines 358-368).
    # Q_TOTAL_LPM is supplied as a covariate; standalone use leaves it
    # at baseline 4.87 L/min so q_scale collapses to ~1 (no
    # amplification). When composed with Mann_2022_respiratory_physiology
    # the hyperperfusion-driven Q rises above baseline during overdose,
    # q_scale saturates toward 2, vc_eff = vc / q_scale halves, and the
    # effect-site concentration doubles.
    q_scale_raw <- 1 + 1 / (1 + exp((1.6 - Q_TOTAL_LPM / 4.87) / 0.05))
    q_scale <- (q_scale_raw < 1) * 1 + (q_scale_raw > 2) * 2 +
               (q_scale_raw >= 1) * (q_scale_raw <= 2) * q_scale_raw
    vc_eff <- vc / q_scale

    # Three-compartment IV disposition with biophase effect compartment.
    # The disposition (central + peripherals + clearance) keeps the
    # nominal vc; only the effect-site equation uses vc_eff per FDA
    # delaymymod.c line 594.
    d/dt(central)     <- -(cl + q + q2) / vc * central +
                          q  / vp  * peripheral1 +
                          q2 / vp2 * peripheral2
    d/dt(peripheral1) <-  q  / vc * central - q  / vp  * peripheral1
    d/dt(peripheral2) <-  q2 / vc * central - q2 / vp2 * peripheral2
    d/dt(effect)      <-  k1 * (central / vc_eff - effect)

    # Outputs
    Cc    <- (central / vc) * 1000                                          # plasma fentanyl (ng/mL) = (mg/L) * 1000
    Ce    <- effect * 1000                                                  # effect-site fentanyl (ng/mL)
    Ce_pM <- effect * 1e9 / fentanyl_mw_g_mol                               # effect-site fentanyl (pM); feeds Mann_2022_mu_receptor_binding L_op slot
  })
}
