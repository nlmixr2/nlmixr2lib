Yau_2023_diazepam_rat_pbpk <- function() {
  description <- paste(
    "PBPK (simplified whole-body, 14 compartments). Preclinical (rat).",
    "Diazepam disposition in rat described by the Yau 2023 simplified",
    "physiologically based pharmacokinetic model with four common tissue",
    "Kpu scaling factors (Model 3D: k-means clustering on rat tissue",
    "composition). Thirteen perfusion-limited tissues (adipose, bone,",
    "brain, gut, heart, kidney, liver, muscle, pancreas, skin, spleen,",
    "stomach and a lumped rest of body) drain to a single central state",
    "that lumps arterial blood, venous blood and lung at quasi-steady",
    "state. Unbound tissue:plasma partition coefficients are predicted",
    "with the Rodgers and Rowland equation for weak bases and then",
    "multiplied by four estimated group scalars, so only the scalars and",
    "the residual error were estimated; total blood clearance was fixed",
    "to the observed 0.915 L/h. Model 3D was one of the two best-fitting",
    "models for diazepam in rat and gave a Vss,b of 1.11 L against an",
    "observed 0.91 L. Fitted to five digitised literature rat studies.",
    sep = " "
  )
  reference <- paste(
    "Yau E, Gertz M, Ogungbenro K, Aarons L, Olivares-Morales A.",
    "A 'middle-out approach' for the prediction of human drug",
    "disposition from preclinical data using simplified",
    "physiologically based pharmacokinetic (PBPK) models. CPT",
    "Pharmacometrics Syst Pharmacol. 2023;12(3):346-359.",
    "doi:10.1002/psp4.12915.",
    "Model structure, the Rodgers and Rowland Kpu implementation and",
    "the Vss,b equation are defined in the companion paper cited",
    "throughout as reference 20: Yau E, Olivares-Morales A,",
    "Ogungbenro K, Aarons L, Gertz M. Investigation of simplified",
    "physiologically-based pharmacokinetic models in rat and human.",
    "CPT Pharmacometrics Syst Pharmacol. 2023;12(3):333-345.",
    "doi:10.1002/psp4.12911.",
    sep = " "
  )
  vignette <- "Yau_2023_middle_out_pbpk_translation"
  units    <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list()

  population <- list(
    species        = "rat (Wistar and Sprague-Dawley)",
    n_subjects     = NA_integer_,
    n_studies      = 5L,
    age_range      = "adult; one study contributed middle-aged and old Wistar rats",
    weight_range   = "250 g reference rat for the physiology (Table S6)",
    sex_female_pct = 0,
    disease_state  = "healthy rats",
    dose_range     = "1 mg to 5 mg/kg intravenous, bolus or 5 min infusion (Table S2)",
    regions        = NA_character_,
    notes          = paste(
      "Five published rat intravenous studies (Yau 2023 Table S2),",
      "contributing a mixture of average and individual plasma or blood",
      "concentration-time profiles digitised with WebPlotDigitizer. All",
      "reported animals were male. Observed rat Vss,b was 0.91 L and the",
      "observed blood clearance 0.915 L/h; clearance was fixed to that",
      "value during fitting."
    )
  )

  ini({
    # Total blood clearance, fixed to the observed value so that all
    # investigated distribution models were compared on an equal
    # footing (paper Methods: "the CL was fixed to the observed CL
    # value").
    lcl <- fixed(log(0.915)); label("Total blood clearance CL (L/h)")
    # ^ Yau 2023 Results, Rat: "Using a blood CL of 0.915 L/h in rat"

    # Tissue Kpu scaling factors -- Model 3D (four scalars,
    # k-means clustering). Estimated on the log scale in NONMEM.
    lsf1 <- log(3.2); label("Kpu scalar group 1 (bone, brain, muscle, pancreas, rest of body)")
    # ^ Yau 2023 Table 2, Model 3D
    lsf2 <- log(23.5); label("Kpu scalar group 2 (kidney, spleen, liver)")
    # ^ Yau 2023 Table 2, Model 3D
    lsf3 <- log(2.02); label("Kpu scalar group 3 (skin, lung, gut, stomach, heart)")
    # ^ Yau 2023 Table 2, Model 3D
    lsf4 <- log(0.265); label("Kpu scalar group 4 (adipose)")
    # ^ Yau 2023 Table 2, Model 3D

    # Between-subject variability on clearance, reported as %CV;
    # omega^2 = log(CV^2 + 1) per the Table S7 footnote.
    etalcl ~ 0.065413   # 26% CV

    propSd <- 0.196; label("Proportional residual error (fraction)")
    # ^ Yau 2023 Table 2, Model 3D proportional error 19.6% (RSE 8.8%)
  })

  model({
    # ===== Drug physicochemistry (Yau 2023 Table S1) =====
    logp <- 2.82   # octanol-water partition coefficient
    pka  <- 3.4   # acid dissociation constant (basic pKa < 7)
    fup  <- 0.1   # fraction unbound in plasma, rat
    bp   <- 0.836   # blood:plasma ratio, rat
    fe   <- 0.009   # fraction excreted unchanged in urine, rat
    fub  <- fup / bp   # fraction unbound in blood

    cl <- exp(lcl + etalcl)

    # ===== 250 g reference rat physiology =====
    # Source: companion paper Table S1 (three decimal places; the lead
    # paper's Table S6 prints the same reference animal to two).
    # Blood flows tabulated in mL/min; converted to L/h (x 0.06).
    q_co       <- 83.085 * 0.06   # cardiac output = arterial = venous = lung flow
    q_adipose   <- 5.816 * 0.06
    q_bone      <- 10.136 * 0.06
    q_brain     <- 1.662 * 0.06
    q_gut       <- 9.139 * 0.06
    q_heart     <- 4.071 * 0.06
    q_kidney    <- 11.715 * 0.06
    q_muscle    <- 23.098 * 0.06
    q_other     <- 7.228 * 0.06
    q_pancreas  <- 1.496 * 0.06
    q_skin      <- 4.819 * 0.06
    q_spleen    <- 0.831 * 0.06
    q_stomach   <- 1.08 * 0.06

    # Tissue volumes tabulated in mL; converted to L (x 0.001).
    v_adipose   <- 19.792 * 0.001
    v_arterial  <- 6.538 * 0.001
    v_bone      <- 5.401 * 0.001
    v_brain     <- 1.37 * 0.001
    v_gut       <- 5.385 * 0.001
    v_heart     <- 0.793 * 0.001
    v_kidney    <- 1.755 * 0.001
    v_liver     <- 8.472 * 0.001
    v_lung      <- 1.202 * 0.001
    v_muscle    <- 97.188 * 0.001
    v_other     <- 27.938 * 0.001
    v_pancreas  <- 0.769 * 0.001
    v_skin      <- 45.745 * 0.001
    v_spleen    <- 0.481 * 0.001
    v_stomach   <- 1.106 * 0.001
    v_venous    <- 13.077 * 0.001

    # Portal-vein pool and hepatic-artery flow. The tabulated
    # "Liver" flow is the portal-vein pool, not total hepatic flow:
    # in rat (12.546 = 9.139 + 1.080 + 0.831 + 1.496 mL/min) and in
    # monkey (123.4 = 90.4 + 11 + 14.9 + 7.1 mL/min) that row equals
    # the splanchnic sum exactly. q_pv is therefore rebuilt from the
    # splanchnic rows and the hepatic artery is recovered as the
    # cardiac-output balance, which makes the ODE system
    # mass-conservative by construction and reproduces the companion
    # paper's own rat control-stream fractions exactly
    # (FCO_PV = 0.151, FCO_HA = 0.024, summing to 0.175 of cardiac
    # output = the 17.50% recovered here). Total hepatic flow then
    # sits at 17.5% of cardiac output in rat, 20.7% in monkey and
    # 20.5% in human. Ratified by the operator (sidecar request-002
    # q2), which also keeps the twice-tabulated rest-of-body flow.
    q_pv <- q_gut + q_stomach + q_spleen + q_pancreas
    q_ha <- q_co - (q_adipose + q_bone + q_brain + q_gut + q_heart +
                    q_kidney + q_muscle + q_other + q_pancreas +
                    q_skin + q_spleen + q_stomach)
    q_hv <- q_ha + q_pv
    # For this reference individual q_pv = 0.75276 L/h, q_ha = 0.11964 L/h and
    # total hepatic blood flow q_hv = 0.8724 L/h.

    # ===== Tissue composition (Yau 2023 Table S3) =====
    # Only the Rodgers & Rowland weak-base inputs are carried:
    # fEW / fIW (extra- and intracellular water), vNL / vNP
    # (neutral lipid and neutral phospholipid) and AR (albumin
    # ratio). Plasma vNL / vNP feed the KAPR affinity constant.
    plas_vnl <- 0.0023
    plas_vnp <- 0.0013

    few_adipose  <- 0.135
    fiw_adipose  <- 0.017
    vnl_adipose  <- 0.853
    vnp_adipose  <- 0.0016
    ar_adipose   <- 0.049

    few_bone     <- 0.1
    fiw_bone     <- 0.346
    vnl_bone     <- 0.0174
    vnp_bone     <- 0.0016
    ar_bone      <- 0.1

    few_brain    <- 0.162
    fiw_brain    <- 0.62
    vnl_brain    <- 0.0391
    vnp_brain    <- 0.0015
    ar_brain     <- 0.048

    few_gut      <- 0.282
    fiw_gut      <- 0.475
    vnl_gut      <- 0.0375
    vnp_gut      <- 0.0124
    ar_gut       <- 0.158

    few_heart    <- 0.32
    fiw_heart    <- 0.456
    vnl_heart    <- 0.0135
    vnp_heart    <- 0.0106
    ar_heart     <- 0.157

    few_kidney   <- 0.273
    fiw_kidney   <- 0.483
    vnl_kidney   <- 0.0121
    vnp_kidney   <- 0.024
    ar_kidney    <- 0.13

    few_liver    <- 0.161
    fiw_liver    <- 0.573
    vnl_liver    <- 0.0135
    vnp_liver    <- 0.0238
    ar_liver     <- 0.086

    few_lung     <- 0.336
    fiw_lung     <- 0.446
    vnl_lung     <- 0.0215
    vnp_lung     <- 0.0123
    ar_lung      <- 0.212

    few_muscle   <- 0.118
    fiw_muscle   <- 0.63
    vnl_muscle   <- 0.01
    vnp_muscle   <- 0.0072
    ar_muscle    <- 0.064

    few_pancreas <- 0.12
    fiw_pancreas <- 0.664
    vnl_pancreas <- 0.0403
    vnp_pancreas <- 0.009
    ar_pancreas  <- 0.06

    few_skin     <- 0.382
    fiw_skin     <- 0.291
    vnl_skin     <- 0.0603
    vnp_skin     <- 0.0044
    ar_skin      <- 0.277

    few_spleen   <- 0.207
    fiw_spleen   <- 0.579
    vnl_spleen   <- 0.0071
    vnp_spleen   <- 0.0107
    ar_spleen    <- 0.097

    # ===== Rodgers & Rowland unbound tissue:plasma partitioning =====
    # Diazepam, midazolam and basmisanil are all weak bases whose
    # only basic pKa is below 7, so the Rodgers & Rowland (2006)
    # equation for acids / very weak bases / neutrals applies, with
    # albumin (AR) as the extracellular binding protein. Encoded
    # exactly as the CLASS == 2 and PKA < 7 branch of the companion
    # paper's NONMEM $PK block (Yau 2023 doi:10.1002/psp4.12911
    # Appendix S1, "Example NONMEM model code for 14 compartmental
    # PBPK model with 3 common scalars").
    ph_p  <- 7.4   # plasma pH
    ph_iw <- 7     # intracellular-water pH
    fi1 <- 10^(pka - ph_iw)   # fraction ionised, intracellular water
    fi2 <- 10^(pka - ph_p)    # fraction ionised, plasma
    pp1 <- 10^logp                    # octanol:water, all tissues but adipose
    kvo <- 10^(1.115 * logp - 1.34)   # vegetable-oil:water, adipose only
    # Affinity for extracellular albumin (weak bases / acids).
    kapr <- 1 / fup - 1 -
      (pp1 * plas_vnl + ((0.3 * pp1) + 0.7) * plas_vnp) / (1 + fi2)

    kpu_adipose <- (1 + fi1) * fiw_adipose / (1 + fi2) + few_adipose +
      kapr * ar_adipose +
      (kvo * vnl_adipose + ((0.3 * kvo) + 0.7) * vnp_adipose) / (1 + fi2)
    kpu_bone <- (1 + fi1) * fiw_bone / (1 + fi2) + few_bone +
      kapr * ar_bone +
      (pp1 * vnl_bone + ((0.3 * pp1) + 0.7) * vnp_bone) / (1 + fi2)
    kpu_brain <- (1 + fi1) * fiw_brain / (1 + fi2) + few_brain +
      kapr * ar_brain +
      (pp1 * vnl_brain + ((0.3 * pp1) + 0.7) * vnp_brain) / (1 + fi2)
    kpu_gut <- (1 + fi1) * fiw_gut / (1 + fi2) + few_gut +
      kapr * ar_gut +
      (pp1 * vnl_gut + ((0.3 * pp1) + 0.7) * vnp_gut) / (1 + fi2)
    kpu_heart <- (1 + fi1) * fiw_heart / (1 + fi2) + few_heart +
      kapr * ar_heart +
      (pp1 * vnl_heart + ((0.3 * pp1) + 0.7) * vnp_heart) / (1 + fi2)
    kpu_kidney <- (1 + fi1) * fiw_kidney / (1 + fi2) + few_kidney +
      kapr * ar_kidney +
      (pp1 * vnl_kidney + ((0.3 * pp1) + 0.7) * vnp_kidney) / (1 + fi2)
    kpu_liver <- (1 + fi1) * fiw_liver / (1 + fi2) + few_liver +
      kapr * ar_liver +
      (pp1 * vnl_liver + ((0.3 * pp1) + 0.7) * vnp_liver) / (1 + fi2)
    kpu_lung <- (1 + fi1) * fiw_lung / (1 + fi2) + few_lung +
      kapr * ar_lung +
      (pp1 * vnl_lung + ((0.3 * pp1) + 0.7) * vnp_lung) / (1 + fi2)
    kpu_muscle <- (1 + fi1) * fiw_muscle / (1 + fi2) + few_muscle +
      kapr * ar_muscle +
      (pp1 * vnl_muscle + ((0.3 * pp1) + 0.7) * vnp_muscle) / (1 + fi2)
    kpu_pancreas <- (1 + fi1) * fiw_pancreas / (1 + fi2) + few_pancreas +
      kapr * ar_pancreas +
      (pp1 * vnl_pancreas + ((0.3 * pp1) + 0.7) * vnp_pancreas) / (1 + fi2)
    kpu_skin <- (1 + fi1) * fiw_skin / (1 + fi2) + few_skin +
      kapr * ar_skin +
      (pp1 * vnl_skin + ((0.3 * pp1) + 0.7) * vnp_skin) / (1 + fi2)
    kpu_spleen <- (1 + fi1) * fiw_spleen / (1 + fi2) + few_spleen +
      kapr * ar_spleen +
      (pp1 * vnl_spleen + ((0.3 * pp1) + 0.7) * vnp_spleen) / (1 + fi2)

    # Tables S3-S5 carry no stomach or rest-of-body composition row.
    # The source NONMEM code drives the stomach from the gut Kpu and
    # the rest of body from the muscle Kpu; Table S11 confirms this
    # (its rest-of-body and muscle Kpus are identical, and its gut
    # and stomach Kpus are identical for Model 3C).
    kpu_stomach <- kpu_gut
    kpu_other   <- kpu_muscle

    # Model 3D tissue -> scalar mapping (Yau 2023 Table 1: four
    # scalars, k-means clustering on rat tissue-composition data).
    sf1 <- exp(lsf1)   # bone, brain, muscle, pancreas, rest of body
    sf2 <- exp(lsf2)   # kidney, spleen, liver
    sf3 <- exp(lsf3)   # skin, lung, gut, stomach, heart
    sf4 <- exp(lsf4)   # adipose

    # Kb = Kpu * SF * fu_b (companion Appendix S1 Eq S7 with Eq 10).
    kb_bone      <- kpu_bone      * sf1 * fub
    kb_brain     <- kpu_brain     * sf1 * fub
    kb_muscle    <- kpu_muscle    * sf1 * fub
    kb_pancreas  <- kpu_pancreas  * sf1 * fub
    kb_other     <- kpu_other     * sf1 * fub
    kb_kidney    <- kpu_kidney    * sf2 * fub
    kb_spleen    <- kpu_spleen    * sf2 * fub
    kb_liver     <- kpu_liver     * sf2 * fub
    kb_skin      <- kpu_skin      * sf3 * fub
    kb_lung      <- kpu_lung      * sf3 * fub
    kb_gut       <- kpu_gut       * sf3 * fub
    kb_stomach   <- kpu_stomach   * sf3 * fub
    kb_heart     <- kpu_heart     * sf3 * fub
    kb_adipose   <- kpu_adipose   * sf4 * fub

    # ===== Clearance partitioning (companion Appendix S1) =====
    # Total blood clearance is fixed to the observed value; the
    # hepatic share is CL * (1 - fe) and is converted to an intrinsic
    # clearance through the well-stirred model (Eq S6). When the
    # observed blood clearance exceeds total hepatic blood flow -- the
    # case the source paper notes for rat -- the source code caps
    # hepatic extraction at 0.99 and assigns the remainder to the
    # renal route so that total clearance still equals the observed
    # value.
    clh_target <- cl * (1 - fe)
    if (clh_target >= q_hv) {
      clint <- 0.99 * q_hv * q_hv / ((q_hv - 0.99 * q_hv) * fub)
      clr   <- cl - (clint * fub * q_hv / (clint * fub + q_hv))
    } else {
      clint <- clh_target * q_hv / ((q_hv - clh_target) * fub)
      clr   <- cl * fe
    }
    clh <- clint * fub * q_hv / (clint * fub + q_hv)

    # ===== Central compartment =====
    # Arterial blood, venous blood and lung are assumed to reach
    # quasi-steady state instantaneously and are lumped into one
    # central state (companion Appendix S1 Eq S17). The amount-to-
    # blood-concentration divisor is V_central * Kb_central =
    # V_art + V_ven + V_lung * Kb_lung.
    v_conv  <- v_arterial + v_venous + v_lung * kb_lung
    c_blood <- central / v_conv

    # ===== Perfusion-limited tissue ODEs (Eq S18) =====
    d/dt(adipose ) <- q_adipose   * (c_blood - adipose  / (v_adipose   * kb_adipose))
    d/dt(bone    ) <- q_bone      * (c_blood - bone     / (v_bone      * kb_bone))
    d/dt(brain   ) <- q_brain     * (c_blood - brain    / (v_brain     * kb_brain))
    d/dt(gut     ) <- q_gut       * (c_blood - gut      / (v_gut       * kb_gut))
    d/dt(heart   ) <- q_heart     * (c_blood - heart    / (v_heart     * kb_heart))
    d/dt(kidney  ) <- q_kidney    * (c_blood - kidney   / (v_kidney    * kb_kidney))
    d/dt(muscle  ) <- q_muscle    * (c_blood - muscle   / (v_muscle    * kb_muscle))
    d/dt(other   ) <- q_other     * (c_blood - other    / (v_other     * kb_other))
    d/dt(pancreas) <- q_pancreas  * (c_blood - pancreas / (v_pancreas  * kb_pancreas))
    d/dt(skin    ) <- q_skin      * (c_blood - skin     / (v_skin      * kb_skin))
    d/dt(spleen  ) <- q_spleen    * (c_blood - spleen   / (v_spleen    * kb_spleen))
    d/dt(stomach ) <- q_stomach   * (c_blood - stomach  / (v_stomach   * kb_stomach))

    # Liver receives hepatic-artery inflow plus the drained
    # splanchnic organs and eliminates drug (Eq S19).
    d/dt(liver) <- q_ha * c_blood +
      q_gut      * gut      / (v_gut      * kb_gut) +
      q_stomach  * stomach  / (v_stomach  * kb_stomach) +
      q_spleen   * spleen   / (v_spleen   * kb_spleen) +
      q_pancreas * pancreas / (v_pancreas * kb_pancreas) -
      (q_hv + clint * fub) * liver / (v_liver * kb_liver)

    # Venous return from the non-splanchnic tissues plus the hepatic
    # vein, less cardiac output and renal elimination (Eq S17).
    d/dt(central) <-
      q_adipose * adipose / (v_adipose * kb_adipose) +
      q_bone    * bone    / (v_bone    * kb_bone) +
      q_brain   * brain   / (v_brain   * kb_brain) +
      q_heart   * heart   / (v_heart   * kb_heart) +
      q_kidney  * kidney  / (v_kidney  * kb_kidney) +
      q_hv      * liver   / (v_liver   * kb_liver) +
      q_muscle  * muscle  / (v_muscle  * kb_muscle) +
      q_other   * other   / (v_other   * kb_other) +
      q_skin    * skin    / (v_skin    * kb_skin) -
      (q_co + clr) * c_blood

    # ===== Steady-state volume of distribution on whole blood =====
    # Companion Appendix S1 Eq S8; only the liver carries an
    # extraction-ratio correction. Reported as an output so the
    # vignette can check it against the published value.
    vssb <- (kpu_lung     * sf3 * v_lung +
             kpu_kidney   * sf2 * v_kidney +
             kpu_heart    * sf3 * v_heart +
             kpu_spleen   * sf2 * v_spleen +
             kpu_liver    * sf2 * v_liver * (1 - clh / q_hv) +
             kpu_stomach  * sf3 * v_stomach +
             kpu_gut      * sf3 * v_gut +
             kpu_pancreas * sf1 * v_pancreas +
             kpu_muscle   * sf1 * v_muscle +
             kpu_bone     * sf1 * v_bone +
             kpu_skin     * sf3 * v_skin +
             kpu_brain    * sf1 * v_brain +
             kpu_other    * sf1 * v_other +
             kpu_adipose  * sf4 * v_adipose) * fup / bp +
            v_arterial + v_venous
    # Published Vss,b for this fit: 1.11 L.

    # The source fitted plasma concentrations for most studies
    # (NONMEM $ERROR: CPL = CBL / BP), with a proportional error.
    Cc <- c_blood / bp
    Cc ~ prop(propSd)
  })
}
