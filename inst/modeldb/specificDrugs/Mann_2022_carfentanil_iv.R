Mann_2022_carfentanil_iv <- function() {
  description <- paste(
    "Three-compartment IV carfentanil population PK with a first-order",
    "biophase (effect-site) equilibrium compartment, used as the agonist",
    "input layer of the Mann 2022 translational opioid-overdose model.",
    "Mann 2022 had no full carfentanil clinical PK study to fit -- only a",
    "single microdose case report (Minkowski 2012) reporting a roughly",
    "45-minute plasma half-life. They therefore took the fentanyl",
    "Algera 2021 PK micro-constants and applied a fixed set of rate-",
    "constant modifications (k_el and k13 divided by 10; k21 and k31",
    "multiplied by 10; k1 increased to 10/min) to reproduce the longer",
    "plasma persistence and faster effect-site equilibration of",
    "carfentanil. The resulting macro-constants encoded in ini() below",
    "yield exactly those micro-constants when combined with",
    "(WT/70)^0.75 / (WT/70) allometric scaling. Inter-subject",
    "variability is assumed equal to the fentanyl model (same omega^2",
    "values, no carfentanil-specific IIV estimated). Outputs plasma",
    "concentration Cc in ng/mL and effect-site Ce in pM for downstream",
    "consumption by Mann_2022_mu_receptor_binding."
  )
  reference <- paste(
    "Mann J, Samieegohar M, Chaturbedi A, Zirkle J, Han X, Ahmadi SF,",
    "Eshleman A, Janowsky A, Wolfrum K, Swanson T, Bloom S, Dahan A,",
    "Olofsen E, Florian J, Strauss DG, Li Z.",
    "Development of a Translational Model to Assess the Impact of Opioid",
    "Overdose and Naloxone Dosing on Respiratory Depression and",
    "Cardiac Arrest. Clin Pharmacol Ther. 2022;112(5):1020-1032.",
    "doi:10.1002/cpt.2696. PMID 35766413.",
    "Half-life anchor: Minkowski CP, Epstein D, Frost JJ, Gorelick DA.",
    "Differential response to IV carfentanil in chronic cocaine users",
    "and healthy controls. Addict Biol. 2012;17(1):149-155.",
    "doi:10.1111/j.1369-1600.2010.00280.x.",
    "Modified-rate-constant prescription:",
    "https://github.com/FDA/Mechanistic-PK-PD-Model-to-Rescue-Opioid-Overdose,",
    "Figure_7/simulateToGetOD_IM.R lines 169-183",
    "(useCarfentanilLikePK = useFasterBiophaseEquilibriumForCarfentanil = yes)."
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
        "Allometric scaling inherited from the fentanyl PK model that",
        "carfentanil was derived from: CL = CL_TV * (WT/70)^0.75 and",
        "V = V_TV * (WT/70). Reference weight 70 kg."
      ),
      source_name        = "weight"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 1L,
    n_studies      = 1L,
    age_range      = "Adult (Minkowski 2012 microdose single-subject anchor)",
    weight_range   = "Mann 2022 virtual-population simulations assume 70 kg",
    sex_female_pct = NA_real_,
    race_ethnicity = NULL,
    disease_state  = paste(
      "No carfentanil-specific popPK clinical study exists. The Mann",
      "2022 carfentanil PK model is a structural derivative of the",
      "Algera 2021 fentanyl PK with hand-set rate-constant modifications",
      "calibrated to reproduce the ~45 minute plasma half-life observed",
      "in the single Minkowski 2012 microdose case report."
    ),
    dose_range     = paste(
      "Minkowski 2012 microdose anchor; Mann 2022 simulates IV bolus",
      "overdose scenarios at 0.012 mg (medium) and 0.022 mg (high) in",
      "chronic users (dose equivalence derived from the carfentanil-to-",
      "fentanyl ratio of minimum cardiac-arrest-inducing dose)."
    ),
    regions        = NA_character_,
    notes          = paste(
      "PK micro-constants from FDA simulateToGetOD_IM.R lines 169-183:",
      "k_el (kout) = fentanyl k_el / 10; k13 = fentanyl k13 / 10;",
      "k21 = fentanyl k21 * 10; k31 = fentanyl k31 * 10; k1 = 10/min",
      "(faster biophase equilibration, replacing fentanyl k1 = 0.18/min).",
      "The macro-constant decomposition below holds V1 unchanged at",
      "10.5 L (same as fentanyl), keeps k12 (Q2/V1) unchanged so Q2",
      "stays 2.04 L/min, and absorbs the other rate-constant changes",
      "into CL, V2, V3, and Q3. Mass-balance check (FDA delaymymod.c",
      "and Mann 2022 Figure S3B): the resulting plasma profile reaches",
      "~50% of microdose concentration at 45 minutes. IIV is assumed",
      "equal to fentanyl (same omega^2 values, no carfentanil-specific",
      "IIV reported)."
    )
  )

  ini({
    # Carfentanil PK macro-constants (Mann 2022 / FDA simulateToGetOD_IM.R
    # lines 169-183, applied to Algera 2021 fentanyl PK micro-constants).
    # All values in nlmixr2-canonical per-minute units; derivation in the
    # in-file source-trace comments below.

    lcl  <- log(0.126)
    label("Clearance (L/min)")                                              # FDA: kout (k_el = CL/V1) / 10 from fentanyl. fentanyl CL = 1.26 L/min -> 0.126 L/min.
    lvc  <- log(10.5)
    label("Central volume V1 (L)")                                          # FDA: V1 (= VP) unchanged from fentanyl = 10.5 L.
    lvp  <- log(1.44)
    label("First peripheral volume V2 (L)")                                 # FDA: k21 (= Q2/V2) * 10; with Q2 held at fentanyl value, V2 = fentanyl V2 / 10 = 14.4 / 10 = 1.44 L.
    lvp2 <- log(1.66)
    label("Second peripheral volume V3 (L)")                                # FDA: k13 (= Q3/V1) / 10 and k31 (= Q3/V3) * 10; combined Q3_carf = Q3_fent / 10 and V3_carf = V3_fent / 100 = 166 / 100 = 1.66 L. Verified: V3 = Q3_carf / k31_carf = 0.228 / 0.137 = 1.66 L
    lq   <- log(2.04)
    label("Inter-compartmental clearance Q2 (L/min)")                       # FDA: k12 (= Q2/V1) unchanged; Q2 = fentanyl Q2 = 2.04 L/min.
    lq2  <- log(0.228)
    label("Inter-compartmental clearance Q3 (L/min)")                       # FDA: k13 / 10 implies Q3 / 10; Q3 = fentanyl Q3 / 10 = 2.28 / 10 = 0.228 L/min.
    lk1  <- fixed(log(10))
    label("Biophase equilibration rate k1 (1/min, FIXED, faster than fentanyl)") # FDA: useFasterBiophaseEquilibriumForCarfentanil = yes -> k1 = 10/min (replaces fentanyl 0.18/min).

    # Fixed allometric exponents (inherited from Algera 2021 fentanyl PK)
    e_wt_cl <- fixed(0.75)
    label("Body-weight allometric exponent on CL and Q (unitless, FIXED)")  # Inherited from fentanyl PK; Mann 2022 carfentanil uses same WT scaling
    e_wt_v  <- fixed(1)
    label("Body-weight allometric exponent on V1, V2, V3 (unitless, FIXED)") # Inherited from fentanyl PK

    # IIV -- same omega^2 values as fentanyl (Mann 2022 Supplement 1 Table
    # S1 footnote: "Carfentanil PK parameters were based on these fentanyl
    # parameters with the same %CV but modified rate constants").
    etalcl  ~ 0.087
    # Same as fentanyl Table S1: omega^2 CL = 0.087 (30% CV)
    etalvc  ~ 0.43
    # Same as fentanyl Table S1: omega^2 V1 = 0.43 (73% CV)
    etalvp  ~ 0.52
    # Same as fentanyl Table S1: omega^2 V2 = 0.52 (82% CV)
    etalvp2 ~ 0.060
    # Same as fentanyl Table S1: omega^2 V3 = 0.060 (25% CV)
    etalq   ~ 0.39
    # Same as fentanyl Table S1: omega^2 Q2 = 0.39 (69% CV)
    etalq2  ~ 0.098
    # Same as fentanyl Table S1: omega^2 Q3 = 0.098 (32% CV)
  })

  model({
    # Allometric weight scaling factor (70 kg reference)
    wt_ratio <- WT / 70

    # Individual structural parameters
    cl  <- exp(lcl  + etalcl)  * wt_ratio^e_wt_cl
    vc  <- exp(lvc  + etalvc)  * wt_ratio^e_wt_v
    vp  <- exp(lvp  + etalvp)  * wt_ratio^e_wt_v
    vp2 <- exp(lvp2 + etalvp2) * wt_ratio^e_wt_v
    q   <- exp(lq   + etalq)   * wt_ratio^e_wt_cl
    q2  <- exp(lq2  + etalq2)  * wt_ratio^e_wt_cl
    k1  <- exp(lk1)

    # Carfentanil free base 394.5 g/mol per FDA delaymymod.c (Mmass
    # parameter when simulateToGetOD_IM.R sets opioid = "carfentanil").
    carfentanil_mw_g_mol <- 394.5

    # Three-compartment IV disposition with biophase effect compartment
    d/dt(central)     <- -(cl + q + q2) / vc * central +
                          q  / vp  * peripheral1 +
                          q2 / vp2 * peripheral2
    d/dt(peripheral1) <-  q  / vc * central - q  / vp  * peripheral1
    d/dt(peripheral2) <-  q2 / vc * central - q2 / vp2 * peripheral2
    d/dt(effect)      <-  k1 * (central / vc - effect)

    # Outputs
    Cc    <- (central / vc) * 1000                                          # plasma carfentanil (ng/mL) = (mg/L) * 1000
    Ce    <- effect * 1000                                                  # effect-site carfentanil (ng/mL)
    Ce_pM <- effect * 1e9 / carfentanil_mw_g_mol                            # effect-site carfentanil (pM); feeds Mann_2022_mu_receptor_binding L_op slot
  })
}
