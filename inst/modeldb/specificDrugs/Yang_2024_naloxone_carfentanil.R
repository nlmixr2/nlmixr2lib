Yang_2024_naloxone_carfentanil <- function() {
  description <- paste(
    "QSP. Mechanistic PK-PD model of carfentanil-induced respiratory",
    "depression and its reversal by intramuscular naloxone auto-injector",
    "(NAI). Couples a two-compartment carfentanil PK model (the parent,",
    "unsuffixed) to the Yang 2024 naloxone PK model (depot plus three",
    "transit absorption compartments and a two-compartment disposition,",
    "all _naloxone-suffixed). Both drugs equilibrate with their own",
    "biophase compartment, expressed in pM rather than ng/mL because the",
    "mu-opioid receptor binding constants are on a molar scale; a single",
    "receptor-occupancy state RL_op carries the fraction of receptors",
    "bound by carfentanil, with naloxone competing through a",
    "quasi-steady-state Hill term. Ventilation is then recovered from",
    "occupancy through the FDA linear ventilation-versus-CO2",
    "steady-state relationship rather than a linear transduction",
    "function. No human carfentanil PK exists, so CL and V were obtained",
    "by interspecies allometric scaling; ini() carries Scenario 1",
    "(regression through the mouse data point, half-life about 5 h), the",
    "case published as control stream (D). For Scenario 2 (regression",
    "through the rabbit data point, half-life about 1 h) override",
    "lcl = log(1.165), lvc = log(14) and lvp = log(100) -- Q is",
    "unchanged. All IIV in this model was assumed at 15 percent."
  )
  reference <- paste(
    "Yang TE, Del Bene F, Lavezzi SM, Iavarone L, Zhang J, Kim J,",
    "Gjurich B, Kessler C. Mechanistic pharmacokinetic-pharmacodynamic",
    "modeling and simulations of naloxone auto-injector 10 mg reversal",
    "of opioid-induced respiratory depression.",
    "CPT Pharmacometrics Syst Pharmacol. 2024;13(10):1722-1733.",
    "doi:10.1002/psp4.13215. PMCID PMC11494827.",
    "Parameter source: Appendix S1 Table S8 (Carfentanil Population",
    "PK-PD Model Structure and Parameter Estimates) and Appendix S1",
    "example NONMEM control stream (D) Carfentanil, whose $THETA /",
    "$OMEGA blocks and $DES / $ERROR equations this file reproduces",
    "one-for-one.",
    "The mu-opioid receptor binding constants (kon, koff, n for both",
    "carfentanil and naloxone) and the ventilation-CO2 constants G,",
    "Bmax, P1 and P2 were taken by Yang 2024 from the FDA repository",
    "https://github.com/FDA/Mechanistic-PK-PD-Model-to-Rescue-Opioid-Overdose",
    "(Appendix S1 reference 8), the same source published as Mann J,",
    "Samieegohar M, Chaturbedi A, et al. Development of a Translational",
    "Model to Assess the Impact of Opioid Overdose and Naloxone Dosing",
    "on Respiratory Depression and Cardiac Arrest.",
    "Clin Pharmacol Ther. 2022;112(5):1020-1032. doi:10.1002/cpt.2696.",
    "The carfentanil allometric scaling used mouse data from Smith LC,",
    "Bremer PT, Hwang CS, et al. J Am Chem Soc. 2019;141(26):10489-10503",
    "and rabbit data from Rigg JR, Wong TY, Horsewood P, Hewson JR.",
    "Br J Anaesth. 1981;53(12):1337-1345 together with Feasel M,",
    "Randall WR, dissertation 2017, http://hdl.handle.net/10713/6764.",
    "Naloxone PK parameters are Yang 2024 Appendix S1 Table S4."
  )
  vignette <- "Yang_2024_naloxone_opioid_reversal"
  units <- list(
    time          = "min",
    dosing        = "ug",
    concentration = "ng/mL"
  )
  dosing <- c("central", "depot_naloxone")

  compartmentData <- list(
    central              = list(analyte = "carfentanil", units = "ug", specimen = "plasma", verified = TRUE),
    peripheral1          = list(analyte = "carfentanil", units = "ug", specimen = "plasma", verified = TRUE),
    effect               = list(analyte = "carfentanil biophase (effect-site) concentration", units = "pM", specimen = "not applicable", verified = TRUE),
    depot_naloxone       = list(analyte = "naloxone", units = "ug", specimen = "administration site", verified = TRUE),
    transit1_naloxone    = list(analyte = "naloxone", units = "ug", specimen = "administration site", verified = TRUE),
    transit2_naloxone    = list(analyte = "naloxone", units = "ug", specimen = "administration site", verified = TRUE),
    transit3_naloxone    = list(analyte = "naloxone", units = "ug", specimen = "administration site", verified = TRUE),
    central_naloxone     = list(analyte = "naloxone", units = "ug", specimen = "plasma", verified = TRUE),
    peripheral1_naloxone = list(analyte = "naloxone", units = "ug", specimen = "plasma", verified = TRUE),
    effect_naloxone      = list(analyte = "naloxone biophase (effect-site) concentration", units = "pM", specimen = "not applicable", verified = TRUE),
    RL_op                = list(analyte = "fraction of mu-opioid receptors bound by carfentanil", units = "fraction", specimen = "not applicable", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Acts only on naloxone apparent clearance:",
        "CL/F = 3.26 * (WT/70)^0.538 (Yang 2024 Table S4). The reference",
        "weight of 70 kg is explicit in control stream (D)",
        "('WT= 70 ; kg'). The allometrically scaled carfentanil PK",
        "parameters are point predictions for a 70 kg human and carry",
        "no further weight scaling."
      ),
      source_name        = "WT"
    )
  )

  population <- list(
    species        = "human (carfentanil PK extrapolated from mouse and rabbit by interspecies allometric scaling)",
    n_subjects     = NA_integer_,
    n_studies      = 2L,
    age_range      = "Healthy adults (naloxone layer 23-54 years)",
    weight_range   = "Simulations fix WT = 70 kg; naloxone PK data spanned 57.2-100.2 kg",
    sex_female_pct = NA_real_,
    race_ethnicity = NULL,
    disease_state  = paste(
      "Healthy adult virtual population. No human carfentanil clinical",
      "PK or PD study exists; the carfentanil layer is entirely",
      "constructed, and Yang 2024 states explicitly that validation of",
      "this model could not be conducted for want of human data. The",
      "naloxone PK layer is Yang 2024's own auto-injector population PK",
      "in 48 healthy adults and IS data-derived."
    ),
    dose_range     = paste(
      "Yang 2024 simulated IV carfentanil at 0.2, 2.2 and 4.2 ug/kg.",
      "The lowest is the human equivalent of the half-maximal effective",
      "dose (bradypnea and/or loss of posture) in non-human primates",
      "given subcutaneously; the highest is 21-fold higher and",
      "represents a lethal dose. Each was paired with naloxone",
      "auto-injector 2 mg or 10 mg IM/SC given once ventilation had",
      "fallen 60 percent from baseline, and a prophylactic arm gave",
      "NAI 10 mg at 5, 15, 30 or 60 min BEFORE the 4.2 ug/kg dose."
    ),
    regions        = NA_character_,
    pk_scenarios   = paste(
      "Two allometric scenarios are reported throughout Yang 2024.",
      "Scenario 1 (encoded in ini() here, and the one published as",
      "control stream (D)): regression through the MOUSE data point,",
      "CL 0.705 L/min, V1 36 L, V2 272 L, terminal half-life about 5 h,",
      "consistent with a human carfentanil case report. Scenario 2:",
      "regression through the RABBIT data point, CL 1.165 L/min,",
      "V1 14 L, V2 100 L, half-life about 1 h, consistent with a report",
      "of human PK to 90 min. Q2 is 3.51 L/min in both, assumed equal",
      "to fentanyl. Yang 2024 notes Q has limited impact because Cmax",
      "depends on V1 only. Scenario 2's smaller V1 produces a sharper",
      "initial ventilation drop but a faster reversal; Scenario 1's",
      "long half-life produces renarcotisation after about 2 h."
    ),
    notes          = paste(
      "Carfentanil PK was constructed by assuming carfentanil and",
      "fentanyl scale with the same allometric exponent b, fitting the",
      "coefficient log(a) through a single animal species, and then",
      "splitting the scaled total volume V into V1 and V2 using the",
      "fentanyl proportion V1 ~ 12 percent of V (V = 308 L in Scenario",
      "1, 114 L in Scenario 2).",
      "Baseline ventilation VB was carried as a per-subject data column",
      "(E0) in the control stream; Appendix S1 states VB and its IIV",
      "were assumed to be 24 (printed as mL/min, evidently a typo for",
      "L/min -- it matches the buprenorphine model's V0 of 24.0 L/min,",
      "and the ventilation equation returns VB exactly when occupancy",
      "is zero) and 21.8 percent.",
      "Every IIV term in the carfentanil layer was assumed at 15",
      "percent because no data support an estimate."
    )
  )

  ini({
    # Every value below is held FIX in Appendix S1 control stream (D).
    # The carfentanil PK block is an allometric point prediction, the
    # receptor-binding and ventilation-CO2 blocks come from the FDA
    # repository, and the naloxone PK block is Yang 2024 Table S4; none
    # was fitted to carfentanil data, of which none exists in humans.

    # --- Naloxone PK (Yang 2024 Table S4; control stream (D) THETA 1-6)
    lktr_naloxone <- fixed(log(0.696))
    label("Naloxone transit absorption rate constant (KTR, 1/min)")            # (D) THETA 1 = 0.696; Table S4
    lcl_naloxone <- fixed(log(3.26))
    label("Naloxone apparent clearance at 70 kg (CL/F, L/min)")                # (D) THETA 2 = 3.26; Table S4
    e_wt_cl_naloxone <- fixed(0.538)
    label("Allometric exponent of body weight on naloxone CL/F (unitless)")    # (D) THETA 3 = 0.538; Table S4
    lvc_naloxone <- fixed(log(404))
    label("Naloxone apparent central volume (V2/F, L)")                        # (D) THETA 4 = 404; Table S4
    lq_naloxone <- fixed(log(0.0847))
    label("Naloxone apparent intercompartmental clearance (Q/F, L/min)")       # (D) THETA 5 = 0.0847; Table S4
    lvp_naloxone <- fixed(log(81.8))
    label("Naloxone apparent peripheral volume (V3/F, L)")                     # (D) THETA 6 = 81.8; Table S4

    # --- Carfentanil PK, Scenario 1 (Table S8; control stream (D)
    # THETA 7-10). For Scenario 2 use log(1.165), log(14), log(100) for
    # lcl, lvc and lvp respectively; lq is identical in both scenarios.
    lcl <- fixed(log(0.705))
    label("Carfentanil clearance, allometric Scenario 1 (CL, L/min)")          # (D) THETA 7 = 0.705; Table S8 Scenario 1 (Scenario 2: 1.165)
    lvc <- fixed(log(36))
    label("Carfentanil central volume, allometric Scenario 1 (V1, L)")         # (D) THETA 8 = 36; Table S8 Scenario 1 (Scenario 2: 14)
    lq <- fixed(log(3.51))
    label("Carfentanil intercompartmental clearance, assumed equal to fentanyl (Q2, L/min)")  # (D) THETA 9 = 3.51; Table S8, both scenarios
    lvp <- fixed(log(272))
    label("Carfentanil peripheral volume, allometric Scenario 1 (V2, L)")      # (D) THETA 10 = 272; Table S8 Scenario 1 (Scenario 2: 100)

    # --- Biophase equilibration (equation 1)
    lke0_naloxone <- fixed(log(0.106))
    label("Naloxone biophase equilibration rate constant (ke0, 1/min)")        # (D) THETA 11 = 0.106; Table S8
    lke0 <- fixed(log(0.0422))
    label("Carfentanil biophase equilibration rate constant, assumed equal to fentanyl (ke0, 1/min)")  # (D) THETA 12 = 0.0422; Table S8

    # --- MOP receptor association / dissociation (equation 13).
    # Concentrations are in pM, so kon has units pM^-n min^-1. These are
    # the FDA repository values; they equal the Mann 2022 Supplement
    # Table S2 per-second constants multiplied by 60 (carfentanil
    # 9.95e-06*60 = 5.97e-04; naloxone 1.67e-04*60 = 1.002e-02;
    # 2.47e-04*60 = 1.482e-02; 3.96e-02*60 = 2.376), an independent
    # confirmation that the pM scale is the intended one.
    lkon <- fixed(log(0.00059724))
    label("Carfentanil receptor association rate constant (kon, pM^-n/min)")   # (D) THETA 13 = 0.00059724; Table S8
    lkoff <- fixed(log(0.014802))
    label("Carfentanil receptor dissociation rate constant (koff, 1/min)")     # (D) THETA 14 = 0.014802; Table S8
    lhill <- fixed(log(1.025))
    label("Carfentanil binding slope of the dose-effect relationship (n, unitless)")  # (D) THETA 15 = 1.025; Table S8
    lkon_naloxone <- fixed(log(0.010014))
    label("Naloxone receptor association rate constant (kon, pM^-n/min)")      # (D) THETA 16 = 0.010014; Table S8
    lkoff_naloxone <- fixed(log(2.3754))
    label("Naloxone receptor dissociation rate constant (koff, 1/min)")        # (D) THETA 17 = 2.3754; Table S8
    lhill_naloxone <- fixed(log(0.858))
    label("Naloxone binding slope of the dose-effect relationship (n, unitless)")  # (D) THETA 18 = 0.858; Table S8

    # --- Ventilation-CO2 steady-state relationship (equation 14).
    # Control stream (D) $ERROR hardcodes these four FDA constants
    # rather than declaring them as THETAs; they are lifted to ini()
    # here so a downstream user can see and vary them.
    g0 <- fixed(0.42)
    label("Baseline slope of the ventilation-PeCO2 curve without opioid (G, L/min/mmHg)")  # (D) $ERROR G0 = 0.42; Appendix S1 Carfentanil section
    bmax_co2 <- fixed(29.65)
    label("Maximum opioid-induced shift of the ventilation-PeCO2 curve (Bmax, mmHg)")      # (D) $ERROR BMAX = 29.65; Appendix S1 Carfentanil section
    p1 <- fixed(5.2)
    label("Occupancy exponent on the ventilation-PeCO2 slope term (P1, unitless)")         # (D) $ERROR P1 = 5.2; Appendix S1 Carfentanil section
    p2 <- fixed(1.629)
    label("Occupancy exponent on the ventilation-PeCO2 offset term (P2, unitless)")        # (D) $ERROR P2 = 1.629; Appendix S1 Carfentanil section
    le0 <- fixed(log(24))
    label("Baseline minute ventilation (VB, L/min)")                           # Appendix S1 Carfentanil section: VB assumed 24 (printed mL/min, a typo for L/min) with 21.8% IIV

    # --- IIV. Control stream (D) $OMEGA entries 1-12, in order. The
    # three naloxone-PK terms are the Table S4 estimates; every
    # carfentanil-layer term is 0.0225 = 0.15^2, each annotated
    # "; assumed 15 % IIV" in the control stream because no carfentanil
    # data exist to estimate them from.
    etalktr_naloxone ~ fixed(0.111)
    etalcl_naloxone ~ fixed(0.0129)
    etalvc_naloxone ~ fixed(0.0658)
    etalcl ~ fixed(0.0225)
    etalvc ~ fixed(0.0225)
    etalke0 ~ fixed(0.0225)
    etalkon ~ fixed(0.0225)
    etalkoff ~ fixed(0.0225)
    etalhill ~ fixed(0.0225)
    etalkon_naloxone ~ fixed(0.0225)
    etalkoff_naloxone ~ fixed(0.0225)
    etalhill_naloxone ~ fixed(0.0225)
    # Baseline ventilation IIV, Appendix S1 Carfentanil section: 21.8%,
    # giving the variance 0.218^2 = 0.047524. Control stream (D) has no
    # $OMEGA entry for it because E0 arrived as a per-subject data
    # column already carrying the variability.
    etale0 ~ fixed(0.047524)

    # Control stream (D) sets the PD additive error THETA(19) to
    # "(0 FIX)" and Table S8 reports no residual-error row for either
    # endpoint, so both residuals are explicit zeros rather than
    # invented values. See the vignette Assumptions and deviations.
    propSd <- fixed(0)
    label("Carfentanil proportional residual error SD (fraction; not reported)")  # Table S8 reports no PK residual-error row
    addSd <- fixed(0)
    label("Ventilation additive residual error SD (L/min)")                    # (D) THETA 19 = 0 FIX
  })

  model({
    # Molecular weights used to convert plasma concentrations from
    # ng/mL to pM, because the receptor binding constants above are on
    # a molar scale. Control stream (D) $DES writes the conversion as
    # CP*(10**6)/MW with NMW = 327.4 and CMW = 394.512: a concentration
    # of 1 ng/mL = 1 ug/L is 1e-6 g/L, i.e. 1e-6/MW mol/L = 1e6/MW pM.
    mw_naloxone <- 327.4
    mw_carfentanil <- 394.512
    # Control stream (D) adds this floor inside the receptor ODE ("DL =
    # 0.000000001 ; to avoid computation difficulty LOG(A(10)) or
    # LOG(A(7)) when A(10) or A(7) are zero"), because it raises the
    # concentrations to a non-integer power via EXP(n*LOG(.)).
    dl <- 1e-9

    # ---- Naloxone individual parameters
    ktr_naloxone <- exp(lktr_naloxone + etalktr_naloxone)
    cl_naloxone <- exp(lcl_naloxone + etalcl_naloxone) * (WT / 70)^e_wt_cl_naloxone
    vc_naloxone <- exp(lvc_naloxone + etalvc_naloxone)
    q_naloxone <- exp(lq_naloxone)
    vp_naloxone <- exp(lvp_naloxone)
    ke0_naloxone <- exp(lke0_naloxone)
    kon_naloxone <- exp(lkon_naloxone + etalkon_naloxone)
    koff_naloxone <- exp(lkoff_naloxone + etalkoff_naloxone)
    hill_naloxone <- exp(lhill_naloxone + etalhill_naloxone)

    # ---- Carfentanil individual parameters
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q <- exp(lq)
    vp <- exp(lvp)
    ke0 <- exp(lke0 + etalke0)
    kon <- exp(lkon + etalkon)
    koff <- exp(lkoff + etalkoff)
    hill <- exp(lhill + etalhill)
    e0 <- exp(le0 + etale0)

    # ---- Naloxone PK, control stream (D) $DES DADT(1)-DADT(6).
    d/dt(depot_naloxone) <- -ktr_naloxone * depot_naloxone
    d/dt(transit1_naloxone) <- ktr_naloxone * depot_naloxone - ktr_naloxone * transit1_naloxone
    d/dt(transit2_naloxone) <- ktr_naloxone * transit1_naloxone - ktr_naloxone * transit2_naloxone
    d/dt(transit3_naloxone) <- ktr_naloxone * transit2_naloxone - ktr_naloxone * transit3_naloxone
    d/dt(central_naloxone) <- ktr_naloxone * transit3_naloxone -
      (cl_naloxone / vc_naloxone) * central_naloxone -
      (q_naloxone / vc_naloxone) * central_naloxone +
      (q_naloxone / vp_naloxone) * peripheral1_naloxone
    d/dt(peripheral1_naloxone) <- (q_naloxone / vc_naloxone) * central_naloxone -
      (q_naloxone / vp_naloxone) * peripheral1_naloxone
    # Naloxone biophase, expressed in pM. Control stream (D):
    # DADT(7) = KE0N*(NCP*(10**6)/NMW - A(7)).
    d/dt(effect_naloxone) <- ke0_naloxone *
      (central_naloxone / vc_naloxone * 1e6 / mw_naloxone - effect_naloxone)

    # ---- Carfentanil PK, control stream (D) $DES DADT(8)-DADT(9).
    d/dt(central) <- -(cl / vc) * central -
      (q / vc) * central + (q / vp) * peripheral1
    d/dt(peripheral1) <- (q / vc) * central - (q / vp) * peripheral1
    # Carfentanil biophase, expressed in pM. Control stream (D):
    # DADT(10) = KE0C*(CCP*(10**6)/CMW - A(10)).
    d/dt(effect) <- ke0 * (central / vc * 1e6 / mw_carfentanil - effect)

    # ---- MOP receptor kinetics with competitive naloxone inhibition.
    # Control stream (D) DADT(11), with EXP(n*LOG(x + DL)) written here
    # as the equivalent (x + DL)^n:
    #   dRL/dt = kon*(Ce_carf+dl)^n1*(1 - RL)
    #            *(1 - (Ce_nal+dl)^n2/((koff_nal/kon_nal) + (Ce_nal+dl)^n2))
    #            - koff*RL
    # The naloxone term is the quasi-steady-state occupancy: only the
    # RATIO koff_naloxone/kon_naloxone enters, which is legitimate here
    # because naloxone unbinds with a half-life of log(2)/2.3754 = 0.29
    # min, far faster than anything else in the system. Note that the
    # free-receptor factor is (1 - RL_op) alone: naloxone occupancy is
    # folded into the bracketed term rather than carried as a second
    # receptor state.
    d/dt(RL_op) <- kon * (effect + dl)^hill * (1 - RL_op) *
      (1 - (effect_naloxone + dl)^hill_naloxone /
        (koff_naloxone / kon_naloxone + (effect_naloxone + dl)^hill_naloxone)) -
      koff * RL_op

    # Control stream (D) initialises both biophase states at 1e-7
    # rather than 0 for the same log-domain reason as dl.
    effect_naloxone(0) <- 1e-7
    effect(0) <- 1e-7

    # ---- Outputs
    Cc <- central / vc
    Cnal <- central_naloxone / vc_naloxone
    # Linear ventilation-versus-CO2 relationship at steady state,
    # equation 14 / control stream (D) $ERROR:
    #   E = (G - G*RL^P1) * (E0/G - Bmax*RL^P2)
    # At zero occupancy this returns G * (E0/G) = E0 exactly, and at
    # full occupancy the first factor vanishes and ventilation is zero.
    VE <- (g0 - g0 * RL_op^p1) * (e0 / g0 - bmax_co2 * RL_op^p2)
    ERATIO <- VE / e0

    Cc ~ prop(propSd)
    VE ~ add(addSd)
  })
}
