AbdullahKoolmees_2024_posaconazole_pbpk <- function() {
  description <- paste(
    "PBPK (14-compartment whole-body, perfusion-limited, R/deSolve).",
    "Posaconazole disposition in a typical 73-kg adult after oral dosing,",
    "built by Abdullah-Koolmees et al. (2024) as the comparator arm to",
    "their voriconazole model, to test whether flucloxacillin lowers",
    "posaconazole exposure the way it lowers voriconazole exposure. The",
    "structure, physiology and blood flows are identical to the",
    "voriconazole model - twelve perfusion-limited tissues plus venous and",
    "arterial blood pools, with gut and spleen draining into the liver -",
    "but the drug is far more highly bound (fraction unbound 0.02) and its",
    "elimination is linear: posaconazole is 83% excreted unchanged and only",
    "minimally metabolised, by UGT1A4 glucuronidation rather than by",
    "cytochrome P450, so a single population hepatic clearance replaces",
    "the saturable multi-CYP liver model. Because the pregnane X receptor",
    "upregulation that flucloxacillin causes acts on CYP enzymes, the",
    "authors assumed posaconazole metabolism is unaffected, which the",
    "model reproduces as unchanged plasma concentrations across every",
    "interaction and inflammation scenario. The companion voriconazole",
    "model is modellib('AbdullahKoolmees_2024_voriconazole_pbpk')."
  )
  reference <- paste(
    "Abdullah-Koolmees H, van den Nieuwendijk JF, ten Hoope SMK, de Leeuw",
    "DC, Franken LGW, Said MM, Seefat MR, Swart EL, Hendrikse NH, Bartelink",
    "IH. Whole Body Physiologically Based Pharmacokinetic Model to Explain",
    "A Patient With Drug-Drug Interaction Between Voriconazole and",
    "Flucloxacillin. Eur J Drug Metab Pharmacokinet. 2024;49(6):689-699.",
    "doi:10.1007/s13318-024-00916-1.",
    "Model equations from Electronic Supplementary Material Appendix I",
    "equations 8-14 and 17 and the accompanying full R/deSolve source",
    "listing (function PBPK_Model_Posaconazole); drug-specific parameters",
    "from Table 1 and from the supplement's posaconazole parameter block;",
    "tissue Kpu, volumes and blood flows from Table 3.",
    sep = " "
  )
  vignette <- "AbdullahKoolmees_2024_voriconazole_flucloxacillin_ddi"
  units    <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    FU = list(
      description        = "Fraction of posaconazole unbound in plasma",
      units              = "fraction (unitless)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Scenario-level, not per-subject-measured, and constant at 0.02",
        "across every simulated scenario. Table 1 gives fraction unbound",
        "0.02 and the Results state explicitly that 'the fraction unbound",
        "posaconazole remained 2%' under both the albumin changes and the",
        "elevated-CRP scenarios, which is the paper's mechanistic reason",
        "for posaconazole being spared. Kept as a data column, matching",
        "the companion voriconazole model, so the two can be driven from",
        "one scenario table. FU enters twice: it converts every tissue Kpu",
        "to a tissue:plasma Kp (Appendix I equation 7, Kp = Kpu * Fu) and",
        "it selects the unbound drug presented to hepatic clearance."
      ),
      source_name        = "FuP"
    ),
    CONMED_FLUCLOXACILLIN = list(
      description        = "Concomitant flucloxacillin administration indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no flucloxacillin)",
      notes              = paste(
        "1 = flucloxacillin co-administered. Declared and carried at its",
        "reference value only: this covariate has NO structural effect in",
        "the posaconazole model, and that is the paper's finding rather",
        "than an omission. Flucloxacillin's perpetrator mechanism is",
        "pregnane-X-receptor-mediated upregulation of CYP3A4, CYP2C9 and",
        "CYP2C19; posaconazole is 83% excreted unchanged and its minor",
        "hepatic route is UGT1A4 glucuronidation (17%), so Methods 2.3",
        "states the assumption that increased UGT1A4 expression through",
        "PXR upregulation 'did not significantly affect posaconazole",
        "metabolism'. Table 6 accordingly reports posaconazole plasma",
        "concentrations as 'Unchanged' in scenarios 3 to 6. The supplement",
        "posaconazole ODE function contains no induction switch. The",
        "column is retained so one scenario table drives both azole",
        "models; see the vignette Errata for the case reports the",
        "Discussion cites of posaconazole levels falling on flucloxacillin,",
        "which this model deliberately does not reproduce."
      ),
      source_name        = "flucloxacillin"
    )
  )

  compartmentData <- list(
    depot    = list(analyte = "posaconazole", units = "mg", specimen = "administration site", verified = TRUE),
    gut      = list(analyte = "posaconazole", units = "mg", specimen = "tissue", verified = TRUE),
    liver    = list(analyte = "posaconazole", units = "mg", specimen = "tissue", verified = TRUE),
    kidney   = list(analyte = "posaconazole", units = "mg", specimen = "tissue", verified = TRUE),
    lung     = list(analyte = "posaconazole", units = "mg", specimen = "tissue", verified = TRUE),
    adipose  = list(analyte = "posaconazole", units = "mg", specimen = "tissue", verified = TRUE),
    heart    = list(analyte = "posaconazole", units = "mg", specimen = "tissue", verified = TRUE),
    brain    = list(analyte = "posaconazole", units = "mg", specimen = "tissue", verified = TRUE),
    bone     = list(analyte = "posaconazole", units = "mg", specimen = "tissue", verified = TRUE),
    muscle   = list(analyte = "posaconazole", units = "mg", specimen = "tissue", verified = TRUE),
    spleen   = list(analyte = "posaconazole", units = "mg", specimen = "tissue", verified = TRUE),
    other    = list(analyte = "posaconazole", units = "mg", specimen = "tissue", verified = TRUE),
    venous   = list(analyte = "posaconazole", units = "mg", specimen = "whole blood", verified = TRUE),
    arterial = list(analyte = "posaconazole", units = "mg", specimen = "whole blood", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 1L,
    n_studies      = 1L,
    age_range      = "48 years (the index case patient)",
    weight_range   = "50-70 kg for the case patient; 73 kg reference adult used for the physiology",
    sex_female_pct = 100,
    disease_state  = paste(
      "Relapsed acute myeloid leukaemia with Staphylococcus aureus",
      "cellulitis and probable pulmonary aspergillosis, treated in the",
      "intensive care unit; posaconazole was substituted for voriconazole",
      "after an undetectable voriconazole trough"
    ),
    dose_range     = paste(
      "Oral posaconazole 300 mg twice daily in the case patient (Fig. 4B);",
      "the tissue-comparison simulation used 400 mg once daily"
    ),
    regions        = "Netherlands (Amsterdam UMC)",
    notes          = paste(
      "Deterministic single-patient simulation model, not a population fit;",
      "no IIV and no residual-error estimate are reported. Posaconazole",
      "trough concentrations measured in the index case were 2.35 mg/L",
      "(2 days after switching), 5.21 mg/L (7 days later) and 3.35 mg/L (at",
      "readmission), all above the 1.5 mg/L target, while flucloxacillin",
      "therapy continued unchanged. Model qualification used healthy-",
      "volunteer plasma concentrations from the Chen 2020 posaconazole",
      "review and post-mortem lung and liver tissue concentrations from",
      "five deceased patients in Blennow 2014, against a 2-fold acceptance",
      "criterion (Table 6). The Discussion records a 2.7-fold",
      "under-prediction of steady-state plasma concentrations after oral",
      "dosing and attributes it to posaconazole's highly variable",
      "bioavailability (reported 8-50%), of which the model deliberately",
      "adopted the 8% worst case."
    )
  )

  ini({
    # ---- Physiology: 73-kg reference adult -------------------------------
    # Identical to the companion voriconazole model; supplement Appendix I
    # R source "System Parameters" block, tabulated in Table 3.
    v_adipose  <- fixed(18.2);  label("Adipose tissue volume (L)")                # Table 3; supplement Vad
    v_bone     <- fixed(10.5);  label("Bone tissue volume (L)")                   # Table 3; supplement Vbo
    v_brain    <- fixed(1.45);  label("Brain tissue volume (L)")                  # Table 3; supplement Vbr
    v_gut      <- fixed(0.65);  label("Gut wall volume (L)")                      # Table 3 '0.65 (wall)'; supplement VguWall
    v_heart    <- fixed(0.33);  label("Heart tissue volume (L)")                  # Table 3; supplement Vhe
    v_kidney   <- fixed(0.31);  label("Kidney tissue volume (L)")                 # Table 3; supplement Vki
    v_liver    <- fixed(1.8);   label("Liver tissue volume (L)")                  # Table 3; supplement Vli
    v_lung     <- fixed(0.5);   label("Lung tissue volume (L)")                   # Table 3; supplement Vlu
    v_muscle   <- fixed(29);    label("Muscle tissue volume (L)")                 # Table 3; supplement Vmu
    v_spleen   <- fixed(0.15);  label("Spleen tissue volume (L)")                 # Table 3; supplement Vsp
    v_blood    <- fixed(5.6);   label("Total blood volume (L)")                   # Table 3; supplement Vbl

    wt_ref     <- fixed(73);    label("Reference body weight used to close the volume balance (kg)")  # supplement `Weight = 73`
    f_venous   <- fixed(0.705); label("Fraction of blood volume that is venous (unitless)")           # supplement `Vve = 0.705*Vbl`
    f_arterial <- fixed(0.295); label("Fraction of blood volume that is arterial (unitless)")         # supplement `Var = 0.295*Vbl`

    q_co        <- fixed(390);   label("Cardiac output (L/h)")                             # supplement `CA = 6.5 * 60`
    fq_adipose  <- fixed(0.05);  label("Fraction of cardiac output to adipose (unitless)") # supplement `Qad = 0.05 * CA` -> 19.5 L/h (Table 3)
    fq_bone     <- fixed(0.05);  label("Fraction of cardiac output to bone (unitless)")    # supplement `Qbo = 0.05 * CA` -> 19.5 L/h (Table 3)
    fq_brain    <- fixed(0.12);  label("Fraction of cardiac output to brain (unitless)")   # supplement `Qbr = 0.12 * CA` -> 46.8 L/h (Table 3)
    fq_gut      <- fixed(0.15);  label("Fraction of cardiac output to gut (unitless)")     # supplement `Qgu = 0.15 * CA` -> 58.5 L/h (Table 3 prints 58.8)
    fq_heart    <- fixed(0.04);  label("Fraction of cardiac output to heart (unitless)")   # supplement `Qhe = 0.04 * CA` -> 15.6 L/h
    fq_kidney   <- fixed(0.19);  label("Fraction of cardiac output to kidney (unitless)")  # supplement `Qki = 0.19 * CA` -> 74.1 L/h (Table 3)
    fq_muscle   <- fixed(0.17);  label("Fraction of cardiac output to muscle (unitless)")  # supplement `Qmu = 0.17 * CA` -> 66.3 L/h (Table 3)
    fq_spleen   <- fixed(0.03);  label("Fraction of cardiac output to spleen (unitless)")  # supplement `Qsp = 0.03 * CA` -> 11.7 L/h (Table 3)
    fq_ha       <- fixed(0.065); label("Fraction of cardiac output to the hepatic artery (unitless)") # supplement `Qha = 0.065 * CA` -> 25.35 L/h

    # ---- Posaconazole tissue-to-unbound-plasma partition coefficients ----
    # Table 3 column "Kpu Posaconazole", predicted with the Rodgers (2005,
    # 2006) weak-base model (Appendix I equations 1-6). Posaconazole is far
    # more lipophilic than voriconazole (log P 5.41 vs 1.8), which is why
    # every Kpu is several-fold larger.
    # Adipose is the one Kpu NOT taken from Table 3, for the same reason as in
    # the companion voriconazole model: Appendix I prescribes the vegetable
    # oil:water coefficient for adipose, Table 1 tabulates it as "(adipose)",
    # and the supplement assigns `P_voP <- 1.115*PP - 1.35` = 4.682 (Table 1
    # prints 4.68). Appendix I equations 1-4 with that value give 6.5501;
    # Table 3 prints 7.17, which those equations give only with the octanol
    # logP. Posaconazole is nearly insensitive to the choice (fraction unbound
    # 0.02 makes Kp_adipose tiny either way - the simulated steady-state mean
    # moves by under 0.01%), but the same rule is applied for consistency with
    # the voriconazole model, where it matters. See the vignette Errata.
    kpu_adipose <- fixed(6.5501); label("Kpu, adipose (unitless)")          # Appendix I eqs 1-4 with the Table 1 Pv (adipose) row; Table 3 prints 7.17 - see vignette Errata
    kpu_bone    <- fixed(5.44);  label("Kpu, bone (unitless)")              # Table 3
    kpu_brain   <- fixed(3.33);  label("Kpu, brain (unitless)")             # Table 3
    kpu_gut     <- fixed(8.73);  label("Kpu, gut (unitless)")               # Table 3
    kpu_heart   <- fixed(8.57);  label("Kpu, heart (unitless)")             # Table 3
    kpu_kidney  <- fixed(7.35);  label("Kpu, kidney (unitless)")            # Table 3
    kpu_liver   <- fixed(5.07);  label("Kpu, liver (unitless)")             # Table 3
    kpu_lung    <- fixed(11.21); label("Kpu, lung (unitless)")              # Table 3
    kpu_muscle  <- fixed(3.95);  label("Kpu, muscle (unitless)")            # Table 3
    kpu_spleen  <- fixed(5.64);  label("Kpu, spleen (unitless)")            # Table 3
    kpu_other   <- fixed(7.09);  label("Kpu, rest of the body (unitless)")  # Table 3 'Rest of the body'

    # ---- Posaconazole drug-specific parameters ---------------------------
    bp      <- fixed(1);           label("Blood:plasma concentration ratio (unitless)")   # Table 1 'Blood:plasma ratio (BP)'
    lka     <- fixed(log(0.795));  label("First-order absorption rate constant (1/h)")    # supplement `KaP <- 0.795` (Pena-Lorenzo tablet popPK); Table 1 prints 0.96 - see vignette Errata
    lfdepot <- fixed(log(0.08));   label("Oral bioavailability factor on the absorption rate (unitless)")  # supplement `FP <- 0.08`; Discussion 'we tested the worst-case scenario of only 8% availability'; Table 1 prints 25% - see vignette Errata

    # ---- Hepatic elimination --------------------------------------------
    # Posaconazole is minimally metabolised (17% UGT1A4 glucuronidation)
    # and 83% excreted unchanged, so a single linear population hepatic
    # clearance replaces the voriconazole saturable-CYP liver model
    # (Appendix I "Drug metabolism and elimination", CLliver).
    lclpop <- fixed(log(8.02)); label("Population plasma clearance before bioavailability correction (L/h)")  # Table 1 'CL (L/h) 8.02'; supplement `Pos_CL_pop = 8.02/FP`

    # ---- Residual error --------------------------------------------------
    # The source performed one deterministic simulation per scenario and
    # reports no residual-error model; fixed at zero rather than invented.
    propSd <- fixed(0); label("Proportional residual error (fraction; ZERO - not reported by the source)")  # Methods 2.4: single-patient simulation, no variability estimated
  })

  model({
    # ================= Physiology ========================================
    v_venous   <- f_venous * v_blood
    v_arterial <- f_arterial * v_blood
    v_other    <- wt_ref - (v_liver + v_kidney + v_spleen + v_heart + v_lung +
                            v_bone + v_brain + v_muscle + v_adipose +
                            v_gut + v_blood)

    q_adipose <- fq_adipose * q_co
    q_bone    <- fq_bone    * q_co
    q_brain   <- fq_brain   * q_co
    q_gut     <- fq_gut     * q_co
    q_heart   <- fq_heart   * q_co
    q_kidney  <- fq_kidney  * q_co
    q_muscle  <- fq_muscle  * q_co
    q_spleen  <- fq_spleen  * q_co
    q_ha      <- fq_ha      * q_co
    q_lung    <- q_co
    q_liver   <- q_gut + q_spleen + q_ha
    q_other   <- q_lung - (q_liver + q_kidney + q_bone + q_heart +
                           q_muscle + q_adipose + q_brain)

    # ================= Individual parameters =============================
    ka <- exp(lka)
    fa <- exp(lfdepot)

    # The source divides the published plasma clearance by the
    # bioavailability factor (`Pos_CL_pop = 8.02/FP`), because the same
    # factor scales the absorption rate rather than the absorbed fraction,
    # so the whole oral dose still reaches the systemic circulation. See
    # the vignette Errata.
    clh <- exp(lclpop) / fa

    # Tissue:plasma partition coefficients, Appendix I equation 7.
    kp_adipose <- kpu_adipose * FU
    kp_bone    <- kpu_bone    * FU
    kp_brain   <- kpu_brain   * FU
    kp_gut     <- kpu_gut     * FU
    kp_heart   <- kpu_heart   * FU
    kp_kidney  <- kpu_kidney  * FU
    kp_liver   <- kpu_liver   * FU
    kp_lung    <- kpu_lung    * FU
    kp_muscle  <- kpu_muscle  * FU
    kp_spleen  <- kpu_spleen  * FU
    kp_other   <- kpu_other   * FU

    # Flucloxacillin is carried for scenario-table compatibility with the
    # companion voriconazole model but has no structural effect here: the
    # authors assumed PXR-driven UGT1A4 upregulation does not meaningfully
    # alter posaconazole metabolism (Methods 2.3), and Table 6 reports
    # posaconazole as "Unchanged" in every interaction scenario.
    flucloxacillin_effect_clh <- 1 + 0 * CONMED_FLUCLOXACILLIN
    clh_eff <- clh * flucloxacillin_effect_clh

    # ================= Concentrations ====================================
    c_venous   <- venous   / v_venous
    c_arterial <- arterial / v_arterial
    c_liver    <- liver    / v_liver

    # ================= ODEs (Appendix I equations 8-14 and 17) ===========
    d/dt(depot) <- -ka * fa * depot

    # Non-eliminating perfusion-limited tissues (equation 8).
    d/dt(adipose) <- q_adipose * c_arterial - q_adipose * ((adipose / v_adipose) / (kp_adipose / bp))
    d/dt(bone)    <- q_bone    * c_arterial - q_bone    * ((bone    / v_bone)    / (kp_bone    / bp))
    d/dt(brain)   <- q_brain   * c_arterial - q_brain   * ((brain   / v_brain)   / (kp_brain   / bp))
    d/dt(heart)   <- q_heart   * c_arterial - q_heart   * ((heart   / v_heart)   / (kp_heart   / bp))
    d/dt(muscle)  <- q_muscle  * c_arterial - q_muscle  * ((muscle  / v_muscle)  / (kp_muscle  / bp))
    d/dt(spleen)  <- q_spleen  * c_arterial - q_spleen  * ((spleen  / v_spleen)  / (kp_spleen  / bp))
    d/dt(other)   <- q_other   * c_arterial - q_other   * ((other   / v_other)   / (kp_other   / bp))

    # Kidney is a non-eliminating compartment in the posaconazole ODE
    # system: the supplement's posaconazole function carries no renal
    # clearance term even though Table S1 attributes 13% of elimination to
    # the renal route. See the vignette Errata.
    d/dt(kidney) <- q_kidney * c_arterial - q_kidney * ((kidney / v_kidney) / (kp_kidney / bp))

    # Lung sits between the venous pool and the arterial pool (equation 9).
    d/dt(lung) <- q_lung * c_venous - q_lung * ((lung / v_lung) / (kp_lung / bp))

    # Gut wall, receiving absorbed drug from the lumen (equations 12-13).
    d/dt(gut) <- q_gut * c_arterial - q_gut * ((gut / v_gut) / (kp_gut / bp)) +
                 ka * fa * depot

    # Liver: hepatic-artery inflow plus portal inflow from gut and spleen,
    # minus hepatic venous outflow and linear clearance of unbound drug.
    d/dt(liver) <- q_gut    * ((gut    / v_gut)    / (kp_gut    / bp)) +
                   q_spleen * ((spleen / v_spleen) / (kp_spleen / bp)) +
                   q_ha     * c_arterial -
                   q_liver  * ((liver  / v_liver)  / (kp_liver  / bp)) -
                   clh_eff * ((c_liver * FU) / (kp_liver / bp))

    # Venous pool (equation 10).
    d/dt(venous) <- -q_lung * c_venous +
                     q_adipose * ((adipose / v_adipose) / (kp_adipose / bp)) +
                     q_heart   * ((heart   / v_heart)   / (kp_heart   / bp)) +
                     q_brain   * ((brain   / v_brain)   / (kp_brain   / bp)) +
                     q_bone    * ((bone    / v_bone)    / (kp_bone    / bp)) +
                     q_muscle  * ((muscle  / v_muscle)  / (kp_muscle  / bp)) +
                     q_kidney  * ((kidney  / v_kidney)  / (kp_kidney  / bp)) +
                     q_liver   * ((liver   / v_liver)   / (kp_liver   / bp)) +
                     q_other   * ((other   / v_other)   / (kp_other   / bp))

    # Arterial pool (equation 11).
    d/dt(arterial) <- q_lung * ((lung / v_lung) / (kp_lung / bp)) -
                      (q_adipose + q_heart + q_brain + q_bone + q_muscle +
                       q_kidney + q_gut + q_spleen + q_ha + q_other) * c_arterial

    # ================= Observation =======================================
    # Total venous plasma concentration (supplement `capture_CP =
    # Cvenous / bpP`). Tissue concentrations are exposed for the Table 6
    # lung and liver comparisons.
    Cc     <- c_venous / bp
    Clung  <- lung  / v_lung
    Cliver <- c_liver

    Cc ~ prop(propSd)
  })
}
