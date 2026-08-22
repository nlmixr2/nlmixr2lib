AbdullahKoolmees_2024_voriconazole_pbpk <- function() {
  description <- paste(
    "PBPK (14-compartment whole-body, perfusion-limited, R/deSolve).",
    "Voriconazole disposition in a typical 73-kg adult after oral dosing,",
    "built by Abdullah-Koolmees et al. (2024) to explain the clinically",
    "observed drug-drug interaction with flucloxacillin. Twelve",
    "perfusion-limited tissues plus venous and arterial blood pools are",
    "driven by Rodgers-type tissue-to-unbound-plasma partition",
    "coefficients (Kpu). Gut and spleen drain into the liver through the",
    "portal circulation; the liver is the site of saturable",
    "Michaelis-Menten metabolism by CYP2C19, CYP3A4 and CYP2C9 acting on",
    "unbound drug, scaled to a whole-organ clearance through the",
    "microsomal-protein-per-gram-liver route, and the kidney carries a",
    "small linear clearance. Two disease/interaction perturbations are",
    "carried as covariates: flucloxacillin co-administration",
    "(CONMED_FLUCLOXACILLIN) lowers the CYP Michaelis constants from a",
    "user-settable onset time, representing PXR-mediated CYP induction,",
    "and an active severe infection (DIS_INFECT_ACTIVE) lowers the CYP",
    "Vmax values. Plasma protein binding enters as the fraction-unbound",
    "covariate FU, which the source varies between the healthy/DDI and",
    "the ICU scenarios. The companion posaconazole model is",
    "modellib('AbdullahKoolmees_2024_posaconazole_pbpk')."
  )
  reference <- paste(
    "Abdullah-Koolmees H, van den Nieuwendijk JF, ten Hoope SMK, de Leeuw",
    "DC, Franken LGW, Said MM, Seefat MR, Swart EL, Hendrikse NH, Bartelink",
    "IH. Whole Body Physiologically Based Pharmacokinetic Model to Explain",
    "A Patient With Drug-Drug Interaction Between Voriconazole and",
    "Flucloxacillin. Eur J Drug Metab Pharmacokinet. 2024;49(6):689-699.",
    "doi:10.1007/s13318-024-00916-1.",
    "Model equations from Electronic Supplementary Material Appendix I",
    "equations 8-18 and the accompanying full R/deSolve source listing;",
    "drug-specific parameters from Table 1; tissue Kpu, volumes and blood",
    "flows from Table 3; CYP Michaelis constants, their flucloxacillin-",
    "induced values and the infection-related Vmax fold decreases from",
    "Table 5.",
    sep = " "
  )
  vignette <- "AbdullahKoolmees_2024_voriconazole_flucloxacillin_ddi"
  units    <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    FU = list(
      description        = "Fraction of voriconazole unbound in plasma",
      units              = "fraction (unitless)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Scenario-level, not per-subject-measured. Table 4 assigns 0.42 to",
        "the healthy volunteer, to the case patient and to every",
        "flucloxacillin-DDI scenario (plasma albumin median 40 g/L), and",
        "0.49 to the ICU scenarios (albumin median 30 g/L). The supplement",
        "R source carries the same three constants as fup = 0.42,",
        "fup_ddi = 0.42 and fup_icu = 0.49. FU enters twice: it converts",
        "every tissue Kpu to a tissue:plasma Kp (Appendix I equation 7,",
        "Kp = Kpu * Fu) and it selects the unbound drug presented to",
        "hepatic and renal clearance. Supply it as a data column so the",
        "published scenarios are reproducible without editing the model."
      ),
      source_name        = "fup"
    ),
    CONMED_FLUCLOXACILLIN = list(
      description        = "Concomitant flucloxacillin administration indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no flucloxacillin)",
      notes              = paste(
        "1 = flucloxacillin co-administered. Flucloxacillin binds the",
        "pregnane X receptor, which upregulates CYP3A4, CYP2C9 and CYP2C19",
        "(Introduction and Discussion). The source implements this as a",
        "step reduction of all three CYP Michaelis constants once",
        "simulation time exceeds tind (Table 5 column 'Km upregulation",
        "(t > 24)'; supplement source `if (t < 24) ... else ...`). The",
        "covariate gates that switch; tind sets when it fires."
      ),
      source_name        = "flucloxacillin"
    ),
    DIS_INFECT_ACTIVE = list(
      description        = "Active severe infection / inflammation episode indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no active severe infection)",
      notes              = paste(
        "1 = the record falls inside a severe infectious-disease episode.",
        "Clinical criterion in this source is the C-reactive-protein band",
        "used to define the simulated patients in Table 4 (bacterial",
        "infected 40-200 mg/L; ICU severely infected > 200 mg/L; healthy",
        "and case patient < 40 mg/L). Divides each CYP Vmax by the",
        "corresponding Table 5 fold decrease. Setting it to 0 reproduces",
        "the supplement's R source exactly - that listing implements only",
        "the flucloxacillin Km switch and leaves Vmax unmodified, so the",
        "Table 5 Vmax column is carried here as an explicitly gated,",
        "off-by-default extension. See the vignette Errata."
      ),
      source_name        = "CRP band"
    )
  )

  compartmentData <- list(
    depot    = list(analyte = "voriconazole", units = "mg", specimen = "administration site", verified = TRUE),
    gut      = list(analyte = "voriconazole", units = "mg", specimen = "tissue", verified = TRUE),
    liver    = list(analyte = "voriconazole", units = "mg", specimen = "tissue", verified = TRUE),
    kidney   = list(analyte = "voriconazole", units = "mg", specimen = "tissue", verified = TRUE),
    lung     = list(analyte = "voriconazole", units = "mg", specimen = "tissue", verified = TRUE),
    adipose  = list(analyte = "voriconazole", units = "mg", specimen = "tissue", verified = TRUE),
    heart    = list(analyte = "voriconazole", units = "mg", specimen = "tissue", verified = TRUE),
    brain    = list(analyte = "voriconazole", units = "mg", specimen = "tissue", verified = TRUE),
    bone     = list(analyte = "voriconazole", units = "mg", specimen = "tissue", verified = TRUE),
    muscle   = list(analyte = "voriconazole", units = "mg", specimen = "tissue", verified = TRUE),
    spleen   = list(analyte = "voriconazole", units = "mg", specimen = "tissue", verified = TRUE),
    other    = list(analyte = "voriconazole", units = "mg", specimen = "tissue", verified = TRUE),
    venous   = list(analyte = "voriconazole", units = "mg", specimen = "whole blood", verified = TRUE),
    arterial = list(analyte = "voriconazole", units = "mg", specimen = "whole blood", verified = TRUE)
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
      "intensive care unit"
    ),
    dose_range     = paste(
      "Oral voriconazole 400 mg twice daily as a 24-h loading dose",
      "followed by 200 mg twice daily maintenance, with intravenous",
      "flucloxacillin 6 g daily throughout"
    ),
    regions        = "Netherlands (Amsterdam UMC)",
    notes          = paste(
      "Deterministic single-patient simulation model, not a population fit.",
      "The authors state explicitly that no pharmacokinetic variability or",
      "covariate variance could be included because of the retrospective",
      "single-case design, so one simulation per scenario was run and the",
      "model carries neither IIV nor a residual-error estimate. Six typical",
      "patients were simulated (Table 4): the case patient (ICU + DDI, low",
      "CRP), a healthy patient, a bacterial infected patient, a DDI patient",
      "with high CRP, an ICU non-infected patient and an ICU severely",
      "infected patient. Model qualification used external literature",
      "rather than a fitted dataset: healthy-volunteer plasma profiles from",
      "Purkins 2002, steady-state trough concentrations across inflammation",
      "stages from van Wanrooy 2014, and post-mortem lung and liver tissue",
      "concentrations from eight deceased patients in Weiler 2011. The",
      "acceptance criterion was agreement within 2-fold (Table 6). The ex",
      "vivo albumin-competition experiment used residual therapeutic-drug-",
      "monitoring samples from 5 patients."
    )
  )

  ini({
    # ---- Physiology: 73-kg reference adult -------------------------------
    # Supplement Appendix I R source, "System Parameters" block. Volumes in
    # L; the tabulated values are reproduced in the paper's Table 3.
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

    # Cardiac output and the tissue blood-flow fractions of cardiac output.
    # Supplement gives CA = 6.5 L/min * 60 = 390 L/h and each Q as a
    # fraction of CA; the resulting L/h values appear in Table 3.
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

    # ---- Voriconazole tissue-to-unbound-plasma partition coefficients ----
    # Table 3 column "Kpu Voriconazole". Predicted with the Rodgers (2005,
    # 2006) weak-base tissue-composition model (Appendix I equations 1-6);
    # the values below are the published predictions, not re-derived here.
    # Adipose is the one Kpu NOT taken from Table 3. Appendix I states
    # "Vegetable oil:water partition coefficient (Pv) is used for adipose
    # tissue", Table 1 tabulates that Pv row explicitly as "(adipose)", and the
    # supplement source assigns it (`P_voV <- 1.115*PV - 1.35` = 0.657, matching
    # the Table 1 value 0.66). Evaluating Appendix I equations 1-4 with Pv gives
    # 0.7815; Table 3 instead prints 1.76, which is what the same equations give
    # if the ordinary octanol logP is used. Table 6's predicted concentrations
    # are reproduced to within 4% by 0.7815 and are 1.45-2.7 fold off with 1.76,
    # so 1.76 is a Table 3 transcription defect. See the vignette Errata.
    kpu_adipose <- fixed(0.7815); label("Kpu, adipose (unitless)")         # Appendix I eqs 1-4 with the Table 1 Pv (adipose) row; Table 3 prints 1.76 - see vignette Errata
    kpu_bone    <- fixed(0.62); label("Kpu, bone (unitless)")              # Table 3
    kpu_brain   <- fixed(0.91); label("Kpu, brain (unitless)")             # Table 3
    kpu_gut     <- fixed(1.06); label("Kpu, gut (unitless)")               # Table 3
    kpu_heart   <- fixed(1.03); label("Kpu, heart (unitless)")             # Table 3
    kpu_kidney  <- fixed(1.00); label("Kpu, kidney (unitless)")            # Table 3
    kpu_liver   <- fixed(0.89); label("Kpu, liver (unitless)")             # Table 3
    kpu_lung    <- fixed(1.08); label("Kpu, lung (unitless)")              # Table 3
    kpu_muscle  <- fixed(0.86); label("Kpu, muscle (unitless)")            # Table 3
    kpu_spleen  <- fixed(0.93); label("Kpu, spleen (unitless)")            # Table 3
    kpu_other   <- fixed(1.03); label("Kpu, rest of the body (unitless)")  # Table 3 'Rest of the body'

    # ---- Voriconazole drug-specific parameters ---------------------------
    bp     <- fixed(1);     label("Blood:plasma concentration ratio (unitless)")   # Table 1 'Blood:plasma ratio (BP)'
    lka    <- fixed(log(0.849)); label("First-order absorption rate constant (1/h)")  # Table 1 'Ka'; supplement KaV
    lfdepot <- fixed(log(0.96)); label("Oral bioavailability factor on the absorption rate (unitless)")  # Table 1 'F (%) 96'; supplement FV

    # ---- Hepatic metabolism ---------------------------------------------
    # Saturable CYP metabolism of unbound drug in the liver, scaled to a
    # whole-organ clearance (Appendix I equations 15-16).
    fumic <- fixed(0.711); label("Fraction unbound in liver microsomes (unitless)")             # supplement `fumic = 0.711`
    mppgl <- fixed(30.3);  label("Microsomal protein per gram of liver (mg/g)")                 # supplement `MPPGL = 30.3`

    vmax_cyp2c19 <- fixed(40);    label("CYP2C19 Vmax (pmol/min/pmol enzyme)")   # Table 1 'Vmax (pmole/min/pmole) CYP2C19'; supplement Vmax_2C19
    vmax_cyp3a4  <- fixed(32.2);  label("CYP3A4 Vmax (pmol/min/pmol enzyme)")    # Table 1 'CYP3A4'; supplement Vmax_3A4
    vmax_cyp2c9  <- fixed(0.056); label("CYP2C9 Vmax (pmol/min/pmol enzyme)")    # Table 1 'CYP3C9' (sic); supplement Vmax_2C9

    km_cyp2c19 <- fixed(9.3);   label("CYP2C19 Michaelis constant, baseline (uM)")  # Table 1 and Table 5 'Km value at baseline'; supplement Km_2C19
    km_cyp3a4  <- fixed(834.7); label("CYP3A4 Michaelis constant, baseline (uM)")   # Table 1 and Table 5; supplement Km_3A4
    km_cyp2c9  <- fixed(20);    label("CYP2C9 Michaelis constant, baseline (uM)")   # Table 1 and Table 5; supplement Km_2C9

    # Flucloxacillin-induced (PXR-upregulated) Michaelis constants.
    km_cyp2c19_ind <- fixed(3.72);   label("CYP2C19 Michaelis constant under flucloxacillin induction (uM)")  # Table 5 'Km upregulation (t > 24)'; supplement Km_2C19_mRNA
    km_cyp3a4_ind  <- fixed(181.45); label("CYP3A4 Michaelis constant under flucloxacillin induction (uM)")   # Table 5; supplement Km_3A4_mRNA
    km_cyp2c9_ind  <- fixed(13.33);  label("CYP2C9 Michaelis constant under flucloxacillin induction (uM)")   # Table 5; supplement Km_2C9_mRNA

    tind <- fixed(24); label("Time after simulation start at which flucloxacillin CYP induction begins (h)")  # supplement `if (t < 24) ... else ...`; Table 5 header 't > 24'

    # Fold decrease in CYP Vmax during a severe infectious-disease episode.
    finf_cyp2c19 <- fixed(1.79); label("Fold decrease in CYP2C19 Vmax during severe infection (unitless)")  # Table 5 'Decrease in Vmax activity (fold)'
    finf_cyp3a4  <- fixed(4.6);  label("Fold decrease in CYP3A4 Vmax during severe infection (unitless)")   # Table 5
    finf_cyp2c9  <- fixed(1.5);  label("Fold decrease in CYP2C9 Vmax during severe infection (unitless)")   # Table 5

    # ---- Renal elimination ----------------------------------------------
    lclr <- fixed(log(0.096)); label("Renal clearance of unbound voriconazole (L/h)")  # supplement `CL_ki = 0.096`; Table 1 renal excretion < 2%

    # ---- Residual error --------------------------------------------------
    # The source performed one deterministic simulation per scenario and
    # reports no residual-error model; fixed at zero rather than invented.
    propSd <- fixed(0); label("Proportional residual error (fraction; ZERO - not reported by the source)")  # Methods 2.4: single-patient simulation, no variability estimated
  })

  model({
    # ================= Physiology ========================================
    # Blood sub-volumes and the rest-of-body volume that closes the
    # 73-kg balance (supplement `Vre = Weight - (sum of tissue volumes)`).
    v_venous   <- f_venous * v_blood
    v_arterial <- f_arterial * v_blood
    v_other    <- wt_ref - (v_liver + v_kidney + v_spleen + v_heart + v_lung +
                            v_bone + v_brain + v_muscle + v_adipose +
                            v_gut + v_blood)

    # Tissue blood flows as fractions of cardiac output. Total hepatic flow
    # is portal (gut + spleen) plus hepatic artery, and the rest-of-body
    # flow is the residual that closes the balance against cardiac output.
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
    ka  <- exp(lka)
    fa  <- exp(lfdepot)
    clr <- exp(lclr)

    # Tissue:plasma partition coefficients, Appendix I equation 7
    # (Kp = Kpu * Funbound). FU is the scenario fraction unbound.
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

    # ================= Disease and interaction perturbations =============
    # Flucloxacillin binds PXR, which upregulates CYP3A4 / CYP2C9 /
    # CYP2C19; the source implements the induction as a step change of the
    # Michaelis constants once time exceeds tind.
    induced <- CONMED_FLUCLOXACILLIN * (t >= tind)
    km_cyp2c19_eff <- km_cyp2c19 + induced * (km_cyp2c19_ind - km_cyp2c19)
    km_cyp3a4_eff  <- km_cyp3a4  + induced * (km_cyp3a4_ind  - km_cyp3a4)
    km_cyp2c9_eff  <- km_cyp2c9  + induced * (km_cyp2c9_ind  - km_cyp2c9)

    # Severe infection lowers CYP Vmax by the Table 5 fold decreases.
    # DIS_INFECT_ACTIVE = 0 leaves Vmax at its baseline value and
    # reproduces the supplement's R source exactly.
    vmax_cyp2c19_eff <- vmax_cyp2c19 / (1 + DIS_INFECT_ACTIVE * (finf_cyp2c19 - 1))
    vmax_cyp3a4_eff  <- vmax_cyp3a4  / (1 + DIS_INFECT_ACTIVE * (finf_cyp3a4  - 1))
    vmax_cyp2c9_eff  <- vmax_cyp2c9  / (1 + DIS_INFECT_ACTIVE * (finf_cyp2c9  - 1))

    # ================= Concentrations ====================================
    c_venous   <- venous   / v_venous
    c_arterial <- arterial / v_arterial
    c_liver    <- liver    / v_liver

    # ================= Hepatic clearance =================================
    # Unbound liver concentration drives the three saturable CYP arms
    # (Appendix I equation 15), summed without contribution weighting
    # exactly as the supplement source does, then scaled to a whole-liver
    # clearance in L/h through MPPGL and fumic (equation 16).
    cu_liver <- c_liver * FU
    clint <- vmax_cyp2c19_eff / (km_cyp2c19_eff + cu_liver) +
             vmax_cyp3a4_eff  / (km_cyp3a4_eff  + cu_liver) +
             vmax_cyp2c9_eff  / (km_cyp2c9_eff  + cu_liver)
    clh <- clint * mppgl * v_liver * 1000 * 60 * 1e-6 / fumic

    # ================= ODEs (Appendix I equations 8-18) ==================
    # Oral dose enters the gut lumen and is delivered to the gut wall.
    # Note that the source multiplies the absorption rate constant by the
    # bioavailability factor rather than splitting the dose, so the whole
    # lumen amount eventually reaches the gut; see the vignette Errata.
    d/dt(depot) <- -ka * fa * depot

    # Non-eliminating perfusion-limited tissues (equation 8).
    d/dt(adipose) <- q_adipose * c_arterial - q_adipose * ((adipose / v_adipose) / (kp_adipose / bp))
    d/dt(bone)    <- q_bone    * c_arterial - q_bone    * ((bone    / v_bone)    / (kp_bone    / bp))
    d/dt(brain)   <- q_brain   * c_arterial - q_brain   * ((brain   / v_brain)   / (kp_brain   / bp))
    d/dt(heart)   <- q_heart   * c_arterial - q_heart   * ((heart   / v_heart)   / (kp_heart   / bp))
    d/dt(muscle)  <- q_muscle  * c_arterial - q_muscle  * ((muscle  / v_muscle)  / (kp_muscle  / bp))
    d/dt(spleen)  <- q_spleen  * c_arterial - q_spleen  * ((spleen  / v_spleen)  / (kp_spleen  / bp))
    d/dt(other)   <- q_other   * c_arterial - q_other   * ((other   / v_other)   / (kp_other   / bp))

    # Lung sits between the venous pool and the arterial pool (equation 9).
    d/dt(lung) <- q_lung * c_venous - q_lung * ((lung / v_lung) / (kp_lung / bp))

    # Gut wall, receiving absorbed drug from the lumen (equations 12-13).
    d/dt(gut) <- q_gut * c_arterial - q_gut * ((gut / v_gut) / (kp_gut / bp)) +
                 ka * fa * depot

    # Liver: hepatic-artery inflow plus portal inflow from gut and spleen,
    # minus hepatic venous outflow and saturable metabolism of unbound
    # drug (equation 17).
    d/dt(liver) <- q_gut    * ((gut    / v_gut)    / (kp_gut    / bp)) +
                   q_spleen * ((spleen / v_spleen) / (kp_spleen / bp)) +
                   q_ha     * c_arterial -
                   q_liver  * ((liver  / v_liver)  / (kp_liver  / bp)) -
                   clh * ((c_liver * FU) / (kp_liver / bp))

    # Kidney with linear clearance of unbound drug (equation 18).
    d/dt(kidney) <- q_kidney * c_arterial -
                    q_kidney * ((kidney / v_kidney) / (kp_kidney / bp)) -
                    clr * (((kidney / v_kidney) * FU) / (kp_kidney / bp))

    # Venous pool: receives every tissue that drains directly to blood
    # (gut and spleen drain to the liver instead) and feeds the lung
    # (equation 10).
    d/dt(venous) <- -q_lung * c_venous +
                     q_adipose * ((adipose / v_adipose) / (kp_adipose / bp)) +
                     q_heart   * ((heart   / v_heart)   / (kp_heart   / bp)) +
                     q_brain   * ((brain   / v_brain)   / (kp_brain   / bp)) +
                     q_bone    * ((bone    / v_bone)    / (kp_bone    / bp)) +
                     q_muscle  * ((muscle  / v_muscle)  / (kp_muscle  / bp)) +
                     q_kidney  * ((kidney  / v_kidney)  / (kp_kidney  / bp)) +
                     q_liver   * ((liver   / v_liver)   / (kp_liver   / bp)) +
                     q_other   * ((other   / v_other)   / (kp_other   / bp))

    # Arterial pool: fed by the lung, drained by every arterial branch
    # (equation 11).
    d/dt(arterial) <- q_lung * ((lung / v_lung) / (kp_lung / bp)) -
                      (q_adipose + q_heart + q_brain + q_bone + q_muscle +
                       q_kidney + q_gut + q_spleen + q_ha + q_other) * c_arterial

    # ================= Observation =======================================
    # Total venous plasma concentration (supplement `capture_CP =
    # Cvenous / bp`). Tissue concentrations are exposed for the Table 6
    # lung and liver comparisons.
    Cc     <- c_venous / bp
    Clung  <- lung  / v_lung
    Cliver <- c_liver

    Cc ~ prop(propSd)
  })
}
