Goryanin_2026_semaglutide_aging_qsp <- function() {
  description <- paste(
    "QSP. Deterministic systems-pharmacology model of aging under GLP-1",
    "receptor agonist (semaglutide) combination therapy in adults with",
    "obesity and/or type 2 diabetes, time in weeks, no between-subject or",
    "residual variability. Encoded verbatim from the authors' deposited SBML",
    "file (model IQANOVA_GLP1_Ageing_QSP, Supporting Information S2), whose",
    "values equal the paper's Table 2. Layer 1: the once-weekly semaglutide",
    "dose D enters a saturating exposure E = D / (1 + D / ED50); E drives",
    "logistic approach-to-maximum weight loss and HbA1c reduction, and",
    "nausea prevalence induced by E and gated by a tolerance state. Layer 2:",
    "four aging states (damage, repair capacity, frailty, biological age",
    "gap) driven by two time-averaged scenario scalars, drug_d (damage",
    "suppression) and drug_r (repair induction), that represent the drug",
    "combination: 0.40 / 0.10 for semaglutide + SGLT2 inhibitor + metformin",
    "and 0.50 / 0.60 (the SBML defaults) for semaglutide + SGLT2 inhibitor +",
    "rapamycin. Layer 3: the cystatin C and leptin biomarker states carry the",
    "SBML kinetic law k_prod * pot_semaglutide * (-11) - k_deg * B verbatim;",
    "the GDF15, adiponectin and eGDR reactions have empty kinetic laws in the",
    "SBML and are not encoded. KNOWN DEFECTS OF THE SOURCE: the deposited",
    "model does not reproduce the paper's own calibration Figure 2 (68-week",
    "weight loss 6.8% / 9.4% at 1.0 / 2.4 mg versus 10.4% / 14.9% in the",
    "calibration data; nausea decays monotonically from 2% instead of",
    "peaking near 21-25% at week 12) or most bars of Figure 6, and the",
    "cystatin C and leptin states settle at impossible changes below -100%",
    "(-156% and -282%). See the vignette before using any output."
  )
  reference <- paste(
    "Goryanin I, Damms B, Goryanin I.",
    "A Systems Pharmacology Model of Aging Identifies Optimal Combination",
    "Therapies With Secondary Benefits on Weight Loss and Metabolic Health.",
    "CPT Pharmacometrics Syst Pharmacol. 2026;15:e70322.",
    "doi:10.1002/psp4.70322.",
    "Model code: Supporting Information S2 (SBML L3V2,",
    "IQANOVA_GLP1_Ageing_QSP); identifiability report S6; MIDD regulatory",
    "report S7.",
    sep = " "
  )
  vignette <- "Goryanin_2026_semaglutide_aging_qsp"
  units <- list(
    time = "week",
    dosing = "mg (semaglutide once-weekly subcutaneous dose, supplied as the covariate DOSE_SEMAGLUTIDE_INJ_MG; no dose events)",
    concentration = paste(
      "bwloss: % body-weight loss from baseline (positive = loss);",
      "hba1cred: absolute HbA1c reduction (% points, positive = reduction);",
      "nausea: nausea prevalence (%); moderator1: tolerance (unitless, 0-1);",
      "damage, repair, frailty, bioagegap: arbitrary units;",
      "cysc, leptin: % change from baseline"
    )
  )

  covariateData <- list(
    DOSE_SEMAGLUTIDE_INJ_MG = list(
      description = "Once-weekly subcutaneous semaglutide maintenance dose",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Constant over the simulation (the model has no titration and no",
        "PK). Enters the saturating exposure E = D / (1 + D / ED50); the",
        "calibration arms are 1.0 mg (STEP 2) and 2.4 mg (STEP 1), and the",
        "SBML default is dose_sema = 2.4 mg. Only the weight, HbA1c and",
        "nausea states depend on it."
      ),
      source_name = "dose_sema (SBML); D (paper Layer Coupling equations)"
    )
  )

  compartmentData <- list(
    bwloss = list(
      analyte = "body weight (% loss from baseline)",
      units = "%",
      specimen = "not applicable",
      verified = TRUE
    ),
    hba1cred = list(
      analyte = "HbA1c (absolute reduction from baseline)",
      units = "%",
      specimen = "not applicable",
      verified = TRUE
    ),
    nausea = list(analyte = "nausea prevalence", units = "%", specimen = "not applicable", verified = TRUE),
    moderator1 = list(analyte = "nausea tolerance", units = "unitless", specimen = "not applicable", verified = TRUE),
    damage = list(
      analyte = "accumulated damage",
      units = "arbitrary units",
      specimen = "not applicable",
      verified = TRUE
    ),
    repair = list(analyte = "repair capacity", units = "arbitrary units", specimen = "not applicable", verified = TRUE),
    frailty = list(analyte = "frailty index", units = "arbitrary units", specimen = "not applicable", verified = TRUE),
    bioagegap = list(
      analyte = "biological age gap",
      units = "arbitrary units",
      specimen = "not applicable",
      verified = TRUE
    ),
    cysc = list(
      analyte = "cystatin C (% change from baseline)",
      units = "%",
      specimen = "serum",
      verified = TRUE
    ),
    leptin = list(
      analyte = "leptin (% change from baseline)",
      units = "%",
      specimen = "serum",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    n_studies = 2L,
    disease_state = "adults with obesity and/or type 2 diabetes",
    dose_range = "semaglutide 1.0 or 2.4 mg SC once weekly for 68 weeks",
    notes = paste(
      "Not fitted to individual data. The semaglutide layer was calibrated",
      "to arm-level mean trajectories digitised from the STEP 1 (2.4 mg) and",
      "STEP 2 (1.0 and 2.4 mg) trials at weeks 0, 4, 12, 24, 36 or 40, 52",
      "and 68 (Supporting Information S3, all_synthetic_data.csv). The aging",
      "and biomarker layers were anchored to literature-derived population",
      "targets (frailty, rapalogue repair induction, metformin damage",
      "modulation, GDF15 / kidney / eGDR studies). The model represents a",
      "single typical adult; there is no virtual population."
    )
  )

  ini({
    # Layer 1 -- semaglutide exposure, weight, HbA1c and nausea. Values are
    # the paper's Table 2 and the SBML listOfParameters (identical).
    ed50 <- fixed(2.0); label("Semaglutide dose giving half-maximal exposure E (mg)") # Table 2 'EC50' 2.0 mg, 'Fixed from literature in scenario simulations'; SBML EC50 = 2; S7 R5 'ASSUMED'
    kbwloss <- 0.20; label("Rate constant of weight loss per unit exposure (1/week)") # Table 2 'k_W' 0.20 wk-1; SBML k_weight = 0.2
    emax_bwloss <- 15.0; label("Maximum body-weight loss (% of baseline)") # Table 2 'Emax_W' 15.0 %; SBML Emax_weight = 15
    khba1cred <- 0.15; label("Rate constant of HbA1c reduction per unit exposure (1/week)") # Table 2 'k_A' 0.15 wk-1; SBML k_HbA1c = 0.15
    emax_hba1cred <- 1.80; label("Maximum absolute HbA1c reduction (% points)") # Table 2 'Emax_A' 1.80 %; SBML Emax_HbA1c = 1.8
    knausea <- 0.35; label("Rate constant of nausea induction per unit exposure (1/week)") # Table 2 'k_N' 0.35 wk-1; SBML k_nausea = 0.35
    kout_nausea <- fixed(0.15); label("First-order loss rate of nausea prevalence (1/week)") # Methods 'dN/dt = k_N E (1 - T) - 0.15 N'; SBML R03 literal 0.15
    ktol <- 0.04; label("Rate constant of nausea-tolerance build-up (1/(week * % prevalence))") # Table 2 'k_T' 0.04 wk-1; SBML k_tolerance = 0.04

    # Layer 2 -- aging states.
    kdam <- 0.10; label("Damage accumulation rate (a.u./week)") # Table 2 'k_dam' 0.10 wk-1; SBML k_damage = 0.1
    krep_clear <- 0.80; label("Repair-mediated damage clearance rate constant (1/(week * a.u.))") # Table 2 'k_rep' 0.80 wk-1; SBML k_repair_clear = 0.8
    krep_deg <- 0.25; label("Repair-capacity decay rate constant (1/week)") # Table 2 'k_Rdeg' 0.25 wk-1; SBML k_repair_decay = 0.25
    krapa <- 0.45; label("Repair-capacity induction rate per unit drug_r (a.u./week)") # Table 2 'k_Rap' 0.45 wk-1; SBML k_rapamycin = 0.45
    kfra <- 0.30; label("Frailty induction rate constant by damage (1/week)") # S6 aging-layer RSE table 'k_F' 0.3000; S7 R5 'k_F' 0.30 wk-1; SBML k_frailty = 0.3
    kfra_rep <- 0.50; label("Frailty suppression rate constant by repair capacity (1/week)") # S6 aging-layer RSE table 'k_Frep' 0.5000; SBML k_frailty_rep = 0.5
    kbag <- 0.40; label("Biological-age-gap rate constant (1/week)") # S6 aging-layer RSE table 'k_B' 0.4000; SBML k_BioAgeGap = 0.4

    # Scenario scalars for the drug combination (time-averaged). Defaults are
    # the SBML's, which wires drug_D_ageing / drug_R_ageing into R05 / R06.
    drug_d <- fixed(0.50); label("Fractional suppression of damage accumulation by the combination (fraction)") # Methods 'drug_D = 0.50' (aging combination; 0.40 metabolic); SBML drug_D_ageing = 0.5
    drug_r <- fixed(0.60); label("Repair-induction scalar of the combination (unitless)") # Methods 'drug_R = 0.60' (aging combination; 0.10 metabolic); SBML drug_R_ageing = 0.6

    # Layer 3 -- biomarker production / degradation constants (SBML only).
    kin_cysc <- 6.8; label("Cystatin C production constant (% change/week)") # SBML k_CysC_prod = 6.8
    kout_cysc <- 0.48; label("Cystatin C degradation rate constant (1/week)") # SBML k_CysC_deg = 0.48
    kin_leptin <- 14.1; label("Leptin production constant (% change/week)") # SBML k_Leptin_prod = 14.1
    kout_leptin <- 0.55; label("Leptin degradation rate constant (1/week)") # SBML k_Leptin_deg = 0.55
  })

  model({
    # Layer 1 exposure. SBML assignment rule E_sema and report S7 Eq. 1 use
    # D / (1 + D / EC50); the main text prints D / (EC50 + D), which is
    # EC50-fold smaller. The deposited code is followed.
    expo <- DOSE_SEMAGLUTIDE_INJ_MG / (1 + DOSE_SEMAGLUTIDE_INJ_MG / ed50)

    # Layer 1 states (SBML reactions R01-R04; S7 Eqs. 2-5).
    d/dt(bwloss) <- kbwloss * expo * (1 - bwloss / emax_bwloss)
    d/dt(hba1cred) <- khba1cred * expo * (1 - hba1cred / emax_hba1cred)
    d/dt(nausea) <- knausea * expo * (1 - moderator1) - kout_nausea * nausea
    d/dt(moderator1) <- ktol * nausea * (1 - moderator1)

    # Layer 2 aging states (SBML reactions R05-R08; S7 Eqs. 6-9). The SBML
    # and S7 write the age-gap rate as k_B * (DAM - 0.5 * REP); the main text
    # prints k_B * DAM - 0.5 * REP. The deposited code is followed.
    d/dt(damage) <- kdam * (1 - drug_d) - krep_clear * repair * damage
    d/dt(repair) <- krapa * drug_r - krep_deg * repair
    d/dt(frailty) <- kfra * damage - kfra_rep * repair
    d/dt(bioagegap) <- kbag * (damage - 0.5 * repair)

    # Layer 3 biomarkers, SBML reactions R09_CysC and R09_LEP verbatim:
    # k_prod * pot_semaglutide * (-11) - k_deg * B with pot_semaglutide = 1.
    # The paper's sigma_i / gamma_i are not reported; -11 is the literal in
    # the deposited kinetic law.
    d/dt(cysc) <- kin_cysc * 1 * (-11) - kout_cysc * cysc
    d/dt(leptin) <- kin_leptin * 1 * (-11) - kout_leptin * leptin

    # SBML initialAmount: N = 2 %, REP = 0.1 a.u.; all other states 0.
    nausea(0) <- 2
    repair(0) <- 0.1
  })
}
