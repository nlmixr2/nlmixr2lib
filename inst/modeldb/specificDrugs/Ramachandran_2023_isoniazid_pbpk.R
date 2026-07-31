Ramachandran_2023_isoniazid_pbpk <- function() {
  description <- paste(
    "PBPK (whole-body, 18 ODEs, MATLAB ode15s). Oral isoniazid disposition",
    "at extrapulmonary tuberculosis (EPTB) sites in a 70-kg reference adult",
    "male (Ramachandran and Gadgil 2023, CPT Pharmacometrics Syst Pharmacol).",
    "Seventeen perfusion-limited well-stirred tissue compartments connected by",
    "blood and lymph flow, plus a gut-lumen state for enterohepatic recycling.",
    "The model's two novel states are a pleural-fluid compartment (a filtrate",
    "of the lung that drains via lymphatics) and a consolidated lymph-node",
    "compartment that collects afferent lymph from every organ except bone and",
    "spleen and returns it to venous blood - the two most common EPTB sites.",
    "Physiological volumes, blood flows, and lymph flows are fixed literature",
    "fractions of body weight, cardiac output (5200 mL/min), and afferent lymph",
    "flow (8 L/day) respectively (Appendix S1 Tables S2 and S3). Tissue:plasma",
    "partition coefficients were computed by the Rodgers and Rowland method",
    "from the drug physicochemistry in Table S4 and are tabulated in Table S5.",
    "Only the first-order oral absorption rate ka and total systemic clearance",
    "CL were estimated, by weighted least squares against reported plasma",
    "concentrations after a 300 mg oral dose (Appendix S1 Table S7). Slow and",
    "fast acetylators are modelled as two discrete cases selected by the",
    "NAT2_SLOW covariate, each with its own ka, CL, and fractional renal",
    "clearance fR; the partition coefficients are identical for the two",
    "groups. Clearance is split into a renal component fR * CL driven by",
    "arterial concentration and a hepatic component (1 - fR) * CL driven by",
    "the hepatic inlet concentration; the hepatic output enters the gut",
    "lumen, which for isoniazid is a terminal sink cleared only by faecal",
    "transit at kF = 0.252 /h (the paper assigns a non-zero gut reabsorption",
    "rate only to rifampicin). The model is a deterministic typical-value",
    "simulation: the paper reports no interindividual variability and no",
    "residual error model, so propSd is fixed at 0."
  )
  reference <- paste(
    "Ramachandran A, Gadgil CJ. A physiologically-based pharmacokinetic model",
    "for tuberculosis drug disposition at extrapulmonary sites.",
    "CPT Pharmacometrics Syst Pharmacol. 2023;12(9):1274-1284.",
    "doi:10.1002/psp4.13008. Model equations, physiological constants, and all",
    "drug-specific parameters are in Appendix S1 (Supporting Information).",
    "Whole-body topology adapted by the authors from Lyons MA, Reisfeld B,",
    "Yang RSH, Lenaerts AJ. A physiologically based pharmacokinetic model of",
    "rifampin in mice. Antimicrob Agents Chemother. 2013;57(4):1763-1771."
  )
  vignette <- "Ramachandran_2023_tuberculosis_eptb_pbpk"
  units    <- list(time = "hour", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    NAT2_SLOW = list(
      description        = "NAT2 slow-acetylator phenotype indicator (1 = slow acetylator, 0 = fast acetylator)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (fast acetylator)",
      notes              = paste(
        "The paper models slow and fast acetylators as two discrete cases,",
        "not as a continuous covariate effect: 'Slow and rapid metabolizers",
        "of isoniazid are considered as discrete cases, each with different",
        "absorption rates and clearance (CL)' (Introduction). Appendix S1",
        "Table S7 reports a separate (ka, CL) pair for each group and",
        "Table S4 reports a separate fractional renal clearance fR (0.07",
        "fast, 0.29 slow). Calibration used studies that specifically",
        "measured plasma concentrations in individuals already phenotyped as",
        "slow or fast acetylators (Gallicano et al. 1994; Appendix S1",
        "Table S6). Tissue:plasma partition coefficients are identical for",
        "the two groups. Set NAT2_SLOW = 1 to simulate a slow acetylator."
      ),
      source_name        = "acetylator status (SA / FA)"
    )
  )

  # Body weight is the scalar that converts every Table S3 volume fraction
  # into a litre volume, and cardiac output / afferent lymph flow are
  # likewise fixed at the Table S2 reference-male values. The paper does not
  # fit or vary any of these, so no covariate enters model(); the reference
  # physiology is hardcoded. Documented here to preserve the provenance
  # while passing the checkModelConventions() unused-covariate check.
  covariatesDataExcluded <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Fixed at the 70-kg 'assumed male individual' of Appendix S1",
        "Table S2. Tissue volumes are fractions of body weight and the",
        "pleural volume / pleural flow are per-kg quantities (0.3 mL/kg,",
        "0.15 mL/kg/h, Table S3), so a different body weight would rescale",
        "the whole physiology. The paper simulates only the 70-kg",
        "reference subject; users wanting other body sizes must rescale",
        "the physiological constants in model()."
      ),
      source_name        = "BW"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = NA_integer_,
    n_studies      = NA_integer_,
    age_range      = "adult",
    weight_range   = "70-kg reference male (Appendix S1 Table S2, 'Physiological parameters for the assumed male individual')",
    sex_female_pct = 0,
    disease_state  = paste(
      "Drug-susceptible tuberculosis is assumed. The paper simulates a",
      "healthy 70-kg reference adult male physiology and states that the",
      "pharmacokinetics of the drugs are taken to be similar in people with",
      "and without TB; age-dependent PK properties are not considered."
    ),
    dose_range     = paste(
      "Oral isoniazid 300 mg single dose. The same dose is used for model",
      "calibration (Figure 2), for plasma / pleura / lung validation",
      "(Figures 3-5), and for the Day 7 EPTB site simulations (Figure 6);",
      "only the source data set differs between calibration and validation."
    ),
    regions        = NA_character_,
    notes          = paste(
      "No individual-level data were fit. ka and CL were estimated by",
      "weighted least squares (fitnlm, MATLAB 2020) against digitised mean",
      "plasma concentration-time data from previously published studies",
      "(Gallicano et al. 1994, separately for fast and slow acetylators;",
      "Appendix S1 Table S6 lists the per-point calibration weights).",
      "Because several data sets were available for each acetylator type,",
      "separate sets of studies were used for calibration and validation.",
      "The model is therefore a typical-value forward simulation tool, not",
      "an individual-level popPK model."
    )
  )

  ini({
    # The ONLY estimated quantities in the model. Everything else in
    # model() is a fixed physiological or physicochemical constant taken
    # from Appendix S1 Tables S2-S5, with the source cited inline.
    # Separate (ka, CL) pairs for the two acetylator phenotypes, selected in
    # model() by the NAT2_SLOW indicator.
    lka_fast <- log(2.86);  label("Absorption rate constant ka, fast acetylators (1/h)")   # Appendix S1 Table S7: isoniazid fast ka = 2.86 /h
    lka_slow <- log(4.11);  label("Absorption rate constant ka, slow acetylators (1/h)")   # Appendix S1 Table S7: isoniazid slow ka = 4.11 /h
    lcl_fast <- log(24.56); label("Total systemic clearance CL, fast acetylators (L/h)")   # Appendix S1 Table S7: isoniazid fast CL = 24.56 L/h
    lcl_slow <- log(9.16);  label("Total systemic clearance CL, slow acetylators (L/h)")   # Appendix S1 Table S7: isoniazid slow CL = 9.16 L/h

    # The paper fits by weighted least squares and reports no residual-error
    # model and no interindividual variability (goodness of fit is reported
    # only as a correlation coefficient r, Figures S1 and S2). Fixing propSd
    # at 0 keeps the model a purely deterministic forward simulation rather
    # than inventing a variance the source does not report. See the vignette
    # Assumptions and deviations section.
    propSd <- fixed(0); label("Proportional residual error (fraction; not reported by the source)")
  })

  model({
    # =================================================================
    # Physiological constants -- Appendix S1 Table S2 (70-kg male)
    # =================================================================
    bw    <- 70          # body weight (kg)
    qc    <- 5200 * 60 / 1000    # cardiac output: 5200 mL/min = 312 L/h
    lymph <- 8 / 24              # afferent lymph flow: 8 L/day = 0.3333 L/h
    kf    <- 0.252               # gut lumen transit rate (1/h)

    # =================================================================
    # Tissue volumes (L) -- Table S3 fractions of body weight.
    # Arterial blood, venous blood, and lymph node are given as absolute
    # litres in Table S3 (1.8 L, 3.6 L, 0.274 L) via a /70 fraction.
    # =================================================================
    v_lung     <- 0.0076  * bw   # 0.532 L
    v_brain    <- 0.02    * bw   # 1.4 L
    v_adipose  <- 0.2142  * bw   # 14.994 L (density 0.916 g/cm^3, Table S3 note a)
    v_heart    <- 0.0047  * bw   # 0.329 L
    v_muscle   <- 0.4     * bw   # 28 L
    v_bone     <- 0.1429  * bw   # 10.003 L (density 1.92 g/cm^3, Table S3 note a)
    v_skin     <- 0.0371  * bw   # 2.597 L
    v_kidney   <- 0.0044  * bw   # 0.308 L
    v_spleen   <- 0.0026  * bw   # 0.182 L
    v_gut      <- 0.0171  * bw   # 1.197 L
    v_liver    <- 0.0257  * bw   # 1.799 L
    v_other    <- 0.04264 * bw   # 2.9848 L (1 - sum of all other volume fractions, Table S3 note d)
    v_lnode    <- 0.274          # combined volume of lymph nodes (Table S3 note b, Shah and Betts 2012)
    v_arterial <- 1.8            # arterial blood (Table S3 note c, Igari et al.)
    v_venous   <- 3.6            # venous blood (Table S3 note c, Igari et al.)
    v_pleura   <- 0.3 / 1000 * bw  # 0.3 mL/kg = 0.021 L (Table S3, not a fraction)

    # =================================================================
    # Blood flow rates (L/h) -- Table S3 fractions of cardiac output.
    # Spleen and gut are reported as absolute mL/min over the 5200 mL/min
    # cardiac output (77/5200 and 1100/5200, Table S3 note e). Liver total
    # flow is the sum of hepatic artery, gut, and spleen flow.
    # =================================================================
    q_brain   <- 0.12    * qc          # 37.44 L/h
    q_adipose <- 0.05    * qc          # 15.6 L/h
    q_heart   <- 0.04    * qc          # 12.48 L/h
    q_muscle  <- 0.17    * qc          # 53.04 L/h
    q_bone    <- 0.05    * qc          # 15.6 L/h
    q_skin    <- 0.05    * qc          # 15.6 L/h
    q_kidney  <- 0.19    * qc          # 59.28 L/h
    q_spleen  <- 77 / 5200   * qc      # 4.62 L/h
    q_gut     <- 1100 / 5200 * qc      # 66 L/h
    q_ha      <- 0.06    * qc          # hepatic artery, 18.72 L/h
    q_other   <- 0.04365 * qc          # 13.6188 L/h (1 - sum of the rest, Table S3 note f)
    q_liver   <- q_ha + q_gut + q_spleen   # 89.34 L/h (Table S3: Q_Li = Q_LA + Q_Gu + Q_Sp)
    q_pleura  <- 0.15 / 1000 * bw      # pleural fluid flow 0.15 mL/kg/h = 0.0105 L/h (Table S3)

    # Total arterial outflow, i.e. the sum in the arterial-blood ODE. Runs
    # over every tissue perfused directly from arterial blood (the lung is
    # excluded; the liver enters as its hepatic-artery branch q_ha, since
    # gut and spleen flow reach the liver via the portal route).
    q_arterial_out <- q_brain + q_adipose + q_heart + q_muscle + q_bone +
      q_skin + q_kidney + q_spleen + q_gut + q_ha + q_other   # = 311.999 L/h ~ qc

    # =================================================================
    # Lymph flow rates (L/h) -- Table S3 fractions of afferent lymph flow.
    # Bone and spleen contribute no lymph (paper Methods).
    # =================================================================
    l_lung    <- 0.03   * lymph
    l_brain   <- 0.0105 * lymph
    l_adipose <- 0.128  * lymph
    l_heart   <- 0.01   * lymph
    l_muscle  <- 0.16   * lymph
    l_bone    <- 0                  # Table S3: bone lymph flow fraction = 0
    l_skin    <- 0.0703 * lymph
    l_kidney  <- 0.085  * lymph
    l_spleen  <- 0                  # Table S3: spleen lymph flow fraction = 0
    l_gut     <- 0.12   * lymph
    l_liver   <- 0.33   * lymph
    l_other   <- 0.0562 * lymph     # 1 - sum of the rest (Table S3 note g)
    # Efferent lymph leaving the node equals the total afferent lymph
    # (the Table S3 fractions sum to exactly 1.000).
    l_lnode   <- lymph

    # =================================================================
    # Drug-specific constants -- rifampicin
    # =================================================================
    # Tissue:plasma partition coefficients, Rodgers and Rowland method
    # (Appendix S1 Table S5, isoniazid column). "Others" is the median of
    # all the other values (Table S5 note a). The pleura is a fluid space
    # and carries NO partition coefficient (Appendix S1 section 3).
    kp_lung     <- 0.7662
    kp_brain    <- 0.7537
    kp_adipose  <- 0.1543
    kp_heart    <- 0.7551
    kp_muscle   <- 0.7208
    kp_bone     <- 0.4330
    kp_skin     <- 0.6675
    kp_kidney   <- 0.7441
    kp_spleen   <- 0.7605
    kp_gut      <- 0.7429
    kp_liver    <- 0.7212
    kp_lnode    <- 0.7556
    kp_other    <- 0.7435

    # Fractional renal clearance (Appendix S1 Table S4, isoniazid fR):
    # 0.07 for fast acetylators, 0.29 for slow acetylators.
    f_renal <- 0.07 * (1 - NAT2_SLOW) + 0.29 * NAT2_SLOW
    # Enterohepatic recycling is a rifampicin-only feature of this model:
    # Methods assigns a non-zero gut reabsorption rate to rifampicin alone
    # ("The gut reabsorption rate for rifampicin was set to 0.17/h based on
    # the reported value"). For isoniazid the gut lumen is therefore a
    # terminal sink, cleared only by faecal transit at kF.
    kr <- 0

    # =================================================================
    # Estimated parameters -- acetylator-phenotype selection
    # =================================================================
    lka_i <- lka_fast * (1 - NAT2_SLOW) + lka_slow * NAT2_SLOW
    lcl_i <- lcl_fast * (1 - NAT2_SLOW) + lcl_slow * NAT2_SLOW
    ka <- exp(lka_i)
    cl <- exp(lcl_i)

    # =================================================================
    # Tissue concentrations (ug/mL == mg/L) and the venous-equilibrium
    # concentrations leaving each tissue, CvT = CT / PT.
    # =================================================================
    Cvenous   <- venous   / v_venous     # plasma (venous blood) concentration
    Carterial <- arterial / v_arterial

    Clung     <- lung     / v_lung
    Cbrain    <- brain    / v_brain
    Cadipose  <- adipose  / v_adipose
    Cheart    <- heart    / v_heart
    Cmuscle   <- muscle   / v_muscle
    Cbone     <- bone     / v_bone
    Cskin     <- skin     / v_skin
    Ckidney   <- kidney   / v_kidney
    Cspleen   <- spleen   / v_spleen
    Cgut      <- gut      / v_gut
    Cliver    <- liver    / v_liver
    Cother    <- other    / v_other
    Clnode    <- lnode    / v_lnode
    Cpleura   <- pleura   / v_pleura     # no partition coefficient: fluid space

    cv_lung    <- Clung    / kp_lung
    cv_brain   <- Cbrain   / kp_brain
    cv_adipose <- Cadipose / kp_adipose
    cv_heart   <- Cheart   / kp_heart
    cv_muscle  <- Cmuscle  / kp_muscle
    cv_bone    <- Cbone    / kp_bone
    cv_skin    <- Cskin    / kp_skin
    cv_kidney  <- Ckidney  / kp_kidney
    cv_spleen  <- Cspleen  / kp_spleen
    cv_gut     <- Cgut     / kp_gut
    cv_liver   <- Cliver   / kp_liver
    cv_other   <- Cother   / kp_other
    cv_lnode   <- Clnode   / kp_lnode

    # =================================================================
    # Clearance terms
    # =================================================================
    # Renal clearance fR * CL, driven by the ARTERIAL concentration, is
    # removed inside the kidney ODE (Appendix S1 section 1, Kidney).
    cl_renal <- f_renal * cl * Carterial

    # Hepatic clearance (1 - fR) * CL, driven by the flow-weighted hepatic
    # inlet concentration. Appendix S1 PRINTS this driver with the tissue
    # concentrations C_Sp and C_Gu, but the venous-equilibrium forms C_VSp
    # and C_VGu (which appear in the inflow terms of the same equations) are
    # used here. The literal printed form does not reproduce the authors'
    # own published simulations: at the 400 mg calibration dose the paper's
    # Figure 2 ethambutol curve peaks at ~0.69 ug/mL, which the
    # venous-equilibrium form reproduces (0.675, -2%) while the literal
    # printed form gives 0.393 (-43%). Ethambutol is the discriminating
    # case because it has by far the largest gut partition coefficient
    # (Kp_gut = 3.21 vs 1.08 rifampicin, 0.74 isoniazid, 0.71
    # pyrazinamide), so only for ethambutol do the two readings diverge
    # materially; the other three drugs match Figure 2 under either form.
    # Treated as a typographical error in Appendix S1. See the vignette
    # Assumptions and deviations section for the full four-drug comparison.
    c_inlet     <- (q_ha * Carterial + q_spleen * cv_spleen + q_gut * cv_gut) / q_liver
    cl_hepatic  <- (1 - f_renal) * cl * c_inlet

    # =================================================================
    # ODE system -- Appendix S1 section 1, 18 equations.
    # States hold AMOUNTS (mg); the published equations are written as
    # V_T * dC_T/dt, which is identical to dA_T/dt.
    # =================================================================

    # Oral dose input.
    d/dt(depot) <- -ka * depot

    # Lungs. The three outflow terms are written as published; they sum to
    # qc * cv_lung, splitting the lung efflux into arterial return,
    # direct lymphatic drainage, and pleural filtrate.
    d/dt(lung) <- qc * Cvenous -
      (qc - l_lung) * cv_lung -
      (l_lung - q_pleura) * cv_lung -
      q_pleura * cv_lung

    # Pleura: a filtrate of the lung, drained by lymphatics. Outflow uses
    # the pleural concentration directly (fluid space, no partitioning).
    d/dt(pleura) <- q_pleura * cv_lung - q_pleura * Cpleura

    # Non-eliminating tissues WITH afferent lymph.
    d/dt(brain)   <- q_brain   * Carterial - (q_brain   - l_brain)   * cv_brain   - l_brain   * cv_brain
    d/dt(heart)   <- q_heart   * Carterial - (q_heart   - l_heart)   * cv_heart   - l_heart   * cv_heart
    d/dt(adipose) <- q_adipose * Carterial - (q_adipose - l_adipose) * cv_adipose - l_adipose * cv_adipose
    d/dt(muscle)  <- q_muscle  * Carterial - (q_muscle  - l_muscle)  * cv_muscle  - l_muscle  * cv_muscle
    d/dt(skin)    <- q_skin    * Carterial - (q_skin    - l_skin)    * cv_skin    - l_skin    * cv_skin
    d/dt(other)   <- q_other   * Carterial - (q_other   - l_other)   * cv_other   - l_other   * cv_other

    # Non-eliminating tissues WITHOUT afferent lymph (bone, spleen).
    d/dt(bone)   <- q_bone   * Carterial - q_bone   * cv_bone
    d/dt(spleen) <- q_spleen * Carterial - q_spleen * cv_spleen

    # Kidney -- renal elimination.
    d/dt(kidney) <- q_kidney * Carterial -
      (q_kidney - l_kidney) * cv_kidney -
      l_kidney * cv_kidney -
      cl_renal

    # Gut tissue -- arterial inflow, portal outflow to the liver, oral
    # absorption from the depot, and reabsorption from the gut lumen.
    d/dt(gut) <- q_gut * Carterial -
      (q_gut - l_gut) * cv_gut -
      l_gut * cv_gut +
      ka * depot +
      kr * gut_lumen

    # Liver -- hepatic artery + portal (spleen, gut) inflow, hepatic
    # elimination.
    d/dt(liver) <- q_ha * Carterial +
      q_spleen * cv_spleen +
      (q_gut - l_gut) * cv_gut -
      (q_liver - l_liver) * cv_liver -
      l_liver * cv_liver -
      cl_hepatic

    # Gut lumen -- receives the hepatic (biliary) output, exits by
    # reabsorption into the gut tissue (kr) and by faecal transit (kF).
    d/dt(gut_lumen) <- cl_hepatic - kr * gut_lumen - kf * gut_lumen

    # Lymph node -- collects afferent lymph from every tissue with a
    # non-zero lymph flow and drains to venous blood.
    d/dt(lnode) <- l_lung    * cv_lung +
      l_brain   * cv_brain +
      l_adipose * cv_adipose +
      l_heart   * cv_heart +
      l_muscle  * cv_muscle +
      l_skin    * cv_skin +
      l_kidney  * cv_kidney +
      l_gut     * cv_gut +
      l_liver   * cv_liver +
      l_other   * cv_other -
      l_lnode   * cv_lnode

    # Arterial blood -- fed by the lung, drained by every perfused tissue.
    d/dt(arterial) <- (qc - l_lung) * cv_lung - q_arterial_out * Carterial

    # Venous blood -- collects the post-lymph venous return of every
    # tissue that drains directly to the vena cava, plus the liver and the
    # lymph node, and delivers cardiac output to the lung.
    d/dt(venous) <- (q_brain   - l_brain)   * cv_brain +
      (q_adipose - l_adipose) * cv_adipose +
      (q_heart   - l_heart)   * cv_heart +
      (q_muscle  - l_muscle)  * cv_muscle +
      (q_bone    - l_bone)    * cv_bone +
      (q_skin    - l_skin)    * cv_skin +
      (q_kidney  - l_kidney)  * cv_kidney +
      (q_other   - l_other)   * cv_other +
      (q_liver   - l_liver)   * cv_liver +
      l_lnode * cv_lnode -
      qc * Cvenous

    Cc <- Cvenous
    Cc ~ prop(propSd)
  })
}
