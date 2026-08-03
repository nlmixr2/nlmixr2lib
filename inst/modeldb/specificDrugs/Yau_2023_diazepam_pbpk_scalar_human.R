Yau_2023_diazepam_pbpk_scalar_human <- function() {
  description <- paste(
    "PBPK (14-compartment whole-body, four common Kpu scalars). Diazepam",
    "disposition in healthy adults after IV administration. Variant 2.2 of",
    "the Yau et al. (2023) simplification strategies: the full kinetic",
    "structure and every individual tissue blood flow and volume are",
    "retained, tissue Kpu values are first predicted bottom-up with the",
    "Rodgers and Rowland (R&R) equations from tissue composition and",
    "diazepam physicochemistry, and then corrected by a small number of",
    "estimated multiplicative scaling factors shared within",
    "composition-based clusters (Eq 10: Kpu_i = Kpu_predRR,i * SF_i). The",
    "scalar therefore quantifies the systematic bias of the bottom-up R&R",
    "prediction rather than replacing it, which is why the paper judged this",
    "parameterisation the most promising for interspecies translation - the",
    "bias is assumed to be similar across species. The R&R weak-base branch",
    "is used because diazepam's basic pKa is 3.4 (< 7), so tissue binding is",
    "driven by extracellular albumin rather than acidic phospholipids. The",
    "lung, arterial and venous blood are treated as a single",
    "quasi-steady-state 'central' compartment; the remaining 13 tissues are",
    "perfusion-limited well-stirred compartments, with stomach, gut, spleen",
    "and pancreas draining into the liver via the portal vein. The four",
    "clusters are (1) bone, brain, muscle, pancreas and rest of body,",
    "(2) lung, gut, stomach, kidney, heart, spleen and liver, (3) skin and",
    "(4) adipose. Together with the four-common-Kpu variant this was one of",
    "the two models the paper found best described the human diazepam data."
  )
  reference <- paste(
    "Yau E, Olivares-Morales A, Ogungbenro K, Aarons L, Gertz M.",
    "Investigation of simplified physiologically-based pharmacokinetic",
    "models in rat and human. CPT Pharmacometrics Syst Pharmacol.",
    "2023;12(3):333-345. doi:10.1002/psp4.12911.",
    "Parameter estimates from Table 2 column '4 scalars'; cluster",
    "membership from the Table 2 footnote; scalar definition from Eq 10;",
    "Rodgers and Rowland equations and tissue composition from the",
    "Appendix S1 example NONMEM control stream and Table S2; model",
    "equations from Appendix S1 Eq S17-S19; reference-man physiology from",
    "Table S1.",
    sep = " "
  )
  vignette <- "Yau_2023_diazepam_pbpk"
  units    <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = 55L,
    n_studies      = 7L,
    age_range      = "adults; study 2 contrasted young vs elderly subjects",
    weight_range   = "70-kg reference man used for the physiological constants (Table S1)",
    sex_female_pct = NA_real_,
    disease_state  = "healthy volunteers / normal subjects",
    dose_range     = paste(
      "IV diazepam 0.1-0.15 mg/kg or 10 mg, given as bolus or as short",
      "infusions of 0.375-2 min across the seven contributing studies",
      "(Table S4)"
    ),
    regions        = NA_character_,
    notes          = paste(
      "35 human IV concentration-time profiles pooled from seven published",
      "studies (Table S4), digitised with WebPlotDigitizer 4.2. Of the 35",
      "profiles, 28 are individual profiles and 7 are study-average",
      "profiles; n_subjects = 55 is the sum of the reported subject counts",
      "across the seven studies. Because most records are arithmetic-mean",
      "profiles from different protocols and populations, the estimated",
      "'IIV' on clearance is confounded with inter-study variability.",
      "Diazepam fu_p = 0.009 and blood:plasma ratio = 0.559 were measured",
      "internally by equilibrium dialysis and LC-MS/MS (Appendix S1);",
      "logP = 2.82 and pKa = 3.4 are from the literature (paper ref 29).",
      "Urinary excretion fraction fe = 0.0005 is from the literature",
      "(paper ref 30)."
    )
  )

  ini({
    # Table 2, column "4 scalars". Values are listed in the paper together
    # with their coefficient of variation in parentheses.
    lcl <- log(3.56); label("Total blood clearance (L/h)")                                              # Yau 2023 Table 2, 4 scalars, CLb = 3.56 L/h (CV 5%)

    # Multiplicative correction factors on the Rodgers and Rowland predicted
    # Kpu, one per tissue-composition cluster (Eq 10; Table 2 footnote,
    # four-cluster hierarchical solution). A value of 1 would mean the
    # bottom-up R&R prediction needs no correction.
    lsf1 <- log(0.206); label("Kpu scaling factor, cluster 1 - bone, brain, muscle, pancreas, rest of body (unitless)")  # Yau 2023 Table 2 (0.206, CV 42%)
    lsf2 <- log(5.70);  label("Kpu scaling factor, cluster 2 - lung, gut, stomach, kidney, heart, spleen, liver (unitless)")  # Yau 2023 Table 2 (5.70, CV 11%)
    lsf3 <- log(6.42);  label("Kpu scaling factor, cluster 3 - skin (unitless)")                        # Yau 2023 Table 2 (6.42, CV 14%)
    lsf4 <- log(8.67);  label("Kpu scaling factor, cluster 4 - adipose (unitless)")                     # Yau 2023 Table 2 (8.67, CV 6%)

    # Exponential IIV on clearance only (paper Data analysis). Reported as
    # 35.2% CV; back-transformed with omega^2 = log(1 + CV^2).
    etalcl ~ 0.1168171; label("IIV on total blood clearance")                                           # Yau 2023 Table 2, IIV CLb = 35.2% (CV 20%)

    # Proportional residual error (supplement $ERROR: W = IPRED*RUV2).
    propSd <- 0.386; label("Proportional residual error (fraction)")                                     # Yau 2023 Table 2, residual error 38.6% (CV 8%)
  })

  model({
    # ================= Drug-specific constants =================
    fup  <- 0.009      # fraction unbound in human plasma (Appendix S1, measured)
    bp    <- 0.559     # blood:plasma concentration ratio, human (Appendix S1, measured)
    fexc  <- 0.0005    # fraction of dose excreted unchanged in urine (paper ref 30)
    fub   <- fup / bp  # fraction unbound in blood

    # ================= Rodgers and Rowland predicted Kpu =================
    # Appendix S1 example NONMEM control stream, "Rowdgers and Rowland model
    # for calculating tissue Kpu". Diazepam is a base with pKa 3.4 < 7, so
    # the weak-base branch applies: the extracellular-protein term uses the
    # albumin ratio (AR) and there is no acidic-phospholipid term.
    logp  <- 2.82      # paper "Diazepam data" (ref 29)
    pka   <- 3.4       # paper "Diazepam data" (ref 29); basic pKa
    ph_p  <- 7.4       # plasma pH
    ph_iw <- 7.0       # intracellular water pH
    fi1 <- 10^(pka - ph_iw)          # fraction ionised in intracellular water
    fi2 <- 10^(pka - ph_p)           # fraction ionised in plasma
    pow <- 10^logp                   # octanol:water partition coefficient
    kvo <- 10^(1.115 * logp - 1.34)  # vegetable-oil:water; adipose only

    # Affinity for extracellular albumin, from plasma neutral lipid and
    # neutral phospholipid fractions (Table S2 human plasma row: vNL =
    # 0.0032, vNP = 0.0021).
    kapr <- (1 / fup) - 1 -
      ((pow * 0.0032 + ((0.3 * pow) + 0.7) * 0.0021) / (1 + fi2))

    # Per-tissue predicted Kpu. Numbers in each expression are, in order,
    # fIW, fEW, AR, vNL and vNP from Table S2 (human tissue composition).
    kpurr_adipose  <- ((1 + fi1) * 0.017) / (1 + fi2) + 0.135 + kapr * 0.049 +
      ((kvo * 0.79 + ((0.3 * kvo) + 0.7) * 0.002) / (1 + fi2))
    kpurr_bone     <- ((1 + fi1) * 0.346) / (1 + fi2) + 0.100 + kapr * 0.100 +
      ((pow * 0.074 + ((0.3 * pow) + 0.7) * 0.0011) / (1 + fi2))
    kpurr_brain    <- ((1 + fi1) * 0.620) / (1 + fi2) + 0.162 + kapr * 0.048 +
      ((pow * 0.051 + ((0.3 * pow) + 0.7) * 0.0565) / (1 + fi2))
    kpurr_gut      <- ((1 + fi1) * 0.475) / (1 + fi2) + 0.282 + kapr * 0.158 +
      ((pow * 0.0487 + ((0.3 * pow) + 0.7) * 0.0163) / (1 + fi2))
    kpurr_heart    <- ((1 + fi1) * 0.456) / (1 + fi2) + 0.320 + kapr * 0.157 +
      ((pow * 0.0115 + ((0.3 * pow) + 0.7) * 0.0166) / (1 + fi2))
    kpurr_kidney   <- ((1 + fi1) * 0.483) / (1 + fi2) + 0.273 + kapr * 0.130 +
      ((pow * 0.0207 + ((0.3 * pow) + 0.7) * 0.0162) / (1 + fi2))
    kpurr_liver    <- ((1 + fi1) * 0.573) / (1 + fi2) + 0.161 + kapr * 0.086 +
      ((pow * 0.0348 + ((0.3 * pow) + 0.7) * 0.0252) / (1 + fi2))
    kpurr_lung     <- ((1 + fi1) * 0.446) / (1 + fi2) + 0.336 + kapr * 0.212 +
      ((pow * 0.003 + ((0.3 * pow) + 0.7) * 0.009) / (1 + fi2))
    kpurr_muscle   <- ((1 + fi1) * 0.630) / (1 + fi2) + 0.079 + kapr * 0.064 +
      ((pow * 0.022 + ((0.3 * pow) + 0.7) * 0.0072) / (1 + fi2))
    kpurr_pancreas <- ((1 + fi1) * 0.664) / (1 + fi2) + 0.120 + kapr * 0.060 +
      ((pow * 0.0403 + ((0.3 * pow) + 0.7) * 0.009) / (1 + fi2))
    kpurr_skin     <- ((1 + fi1) * 0.291) / (1 + fi2) + 0.382 + kapr * 0.277 +
      ((pow * 0.0284 + ((0.3 * pow) + 0.7) * 0.0111) / (1 + fi2))
    kpurr_spleen   <- ((1 + fi1) * 0.579) / (1 + fi2) + 0.207 + kapr * 0.097 +
      ((pow * 0.0201 + ((0.3 * pow) + 0.7) * 0.0198) / (1 + fi2))
    # Table S2 has no stomach row; the control stream uses the gut Kpu for
    # stomach (DADT(11) uses KPU_GU) and the muscle Kpu for the rest of body
    # ("Assume same Kpu for Rest of body as for Muscle": KPU_RO = KPU_MU).
    kpurr_stomach  <- kpurr_gut
    kpurr_other    <- kpurr_muscle

    # ================= Physiology: 70-kg reference man (Table S1) =========
    q_adipose  <- 0.292 * 60
    q_bone     <- 0.292 * 60
    q_brain    <- 0.701 * 60
    q_gut      <- 0.526 * 60
    q_heart    <- 0.234 * 60
    q_kidney   <- 1.109 * 60
    q_liver    <- 1.489 * 60   # total hepatic flow (hepatic artery + portal vein)
    q_muscle   <- 0.993 * 60
    q_pancreas <- 0.058 * 60
    q_skin     <- 0.292 * 60
    q_spleen   <- 0.175 * 60
    q_stomach  <- 0.058 * 60
    q_co       <- 5.839 * 60   # cardiac output

    v_adipose   <- 15.619
    v_bone      <- 5.210
    v_brain     <- 1.346
    v_gut       <- 1.010
    v_heart     <- 0.316
    v_kidney    <- 0.296
    v_liver     <- 1.666
    v_lung      <- 0.512
    v_muscle    <- 26.923
    v_pancreas  <- 0.094
    v_skin      <- 2.497
    v_spleen    <- 0.175
    v_stomach   <- 0.141
    v_arterial  <- 1.329
    v_venous    <- 3.988
    v_other     <- 2.914

    q_pv <- q_gut + q_stomach + q_spleen + q_pancreas
    q_ha <- q_liver - q_pv
    q_other <- q_co - (q_adipose + q_bone + q_brain + q_gut + q_heart +
                       q_kidney + q_ha + q_muscle + q_pancreas +
                       q_skin + q_spleen + q_stomach)

    # ================= Individual parameters =================
    cl  <- exp(lcl + etalcl)
    sf1 <- exp(lsf1)
    sf2 <- exp(lsf2)
    sf3 <- exp(lsf3)
    sf4 <- exp(lsf4)

    # ---- Eq 10: Kpu_i = Kpu_predRR,i * SF_cluster(i) ----
    # Cluster membership per the Table 2 footnote (4 scalars).
    kpu_bone     <- kpurr_bone     * sf1
    kpu_brain    <- kpurr_brain    * sf1
    kpu_muscle   <- kpurr_muscle   * sf1
    kpu_pancreas <- kpurr_pancreas * sf1
    kpu_other    <- kpurr_other    * sf1
    kpu_lung     <- kpurr_lung     * sf2
    kpu_gut      <- kpurr_gut      * sf2
    kpu_stomach  <- kpurr_stomach  * sf2
    kpu_kidney   <- kpurr_kidney   * sf2
    kpu_heart    <- kpurr_heart    * sf2
    kpu_spleen   <- kpurr_spleen   * sf2
    kpu_liver    <- kpurr_liver    * sf2
    kpu_skin     <- kpurr_skin     * sf3
    kpu_adipose  <- kpurr_adipose  * sf4

    # ---- Tissue:blood partition coefficients, Eq S7 ----
    kb_adipose  <- kpu_adipose  * fub
    kb_bone     <- kpu_bone     * fub
    kb_brain    <- kpu_brain    * fub
    kb_gut      <- kpu_gut      * fub
    kb_heart    <- kpu_heart    * fub
    kb_kidney   <- kpu_kidney   * fub
    kb_liver    <- kpu_liver    * fub
    kb_lung     <- kpu_lung     * fub
    kb_muscle   <- kpu_muscle   * fub
    kb_pancreas <- kpu_pancreas * fub
    kb_skin     <- kpu_skin     * fub
    kb_spleen   <- kpu_spleen   * fub
    kb_stomach  <- kpu_stomach  * fub
    kb_other    <- kpu_other    * fub

    # ================= Clearance partitioning (Appendix S1 $PK) ===========
    clh   <- min(cl * (1 - fexc), 0.99 * q_liver)
    clr   <- cl - clh
    clint <- clh * q_liver / ((q_liver - clh) * fub)   # Eq S6 rearranged
    erh   <- clh / q_liver

    # ================= Central compartment volume term ===================
    v_conv <- v_arterial + v_venous + v_lung * kb_lung
    cblood <- central / v_conv

    # ================= ODEs (Eq S17-S19) =================================
    d/dt(central) <-
        q_adipose * (adipose  / v_adipose  / kb_adipose) +
        q_bone    * (bone     / v_bone     / kb_bone) +
        q_brain   * (brain    / v_brain    / kb_brain) +
        q_heart   * (heart    / v_heart    / kb_heart) +
        q_kidney  * (kidney   / v_kidney   / kb_kidney) +
        q_liver   * (liver    / v_liver    / kb_liver) +
        q_muscle  * (muscle   / v_muscle   / kb_muscle) +
        q_skin    * (skin     / v_skin     / kb_skin) +
        q_other   * (other    / v_other    / kb_other) -
        q_co * cblood - clr * cblood

    d/dt(adipose)  <- q_adipose  * cblood - q_adipose  * (adipose  / v_adipose  / kb_adipose)
    d/dt(bone)     <- q_bone     * cblood - q_bone     * (bone     / v_bone     / kb_bone)
    d/dt(brain)    <- q_brain    * cblood - q_brain    * (brain    / v_brain    / kb_brain)
    d/dt(heart)    <- q_heart    * cblood - q_heart    * (heart    / v_heart    / kb_heart)
    d/dt(kidney)   <- q_kidney   * cblood - q_kidney   * (kidney   / v_kidney   / kb_kidney)
    d/dt(muscle)   <- q_muscle   * cblood - q_muscle   * (muscle   / v_muscle   / kb_muscle)
    d/dt(skin)     <- q_skin     * cblood - q_skin     * (skin     / v_skin     / kb_skin)
    d/dt(other)    <- q_other    * cblood - q_other    * (other    / v_other    / kb_other)

    d/dt(gut)      <- q_gut      * cblood - q_gut      * (gut      / v_gut      / kb_gut)
    d/dt(stomach)  <- q_stomach  * cblood - q_stomach  * (stomach  / v_stomach  / kb_stomach)
    d/dt(spleen)   <- q_spleen   * cblood - q_spleen   * (spleen   / v_spleen   / kb_spleen)
    d/dt(pancreas) <- q_pancreas * cblood - q_pancreas * (pancreas / v_pancreas / kb_pancreas)

    d/dt(liver) <-
        q_ha * cblood +
        q_gut      * (gut      / v_gut      / kb_gut) +
        q_stomach  * (stomach  / v_stomach  / kb_stomach) +
        q_spleen   * (spleen   / v_spleen   / kb_spleen) +
        q_pancreas * (pancreas / v_pancreas / kb_pancreas) -
        (liver / v_liver) * ((q_liver + clint * fub) / kb_liver)

    # ================= Derived quantities and observation ================
    # Eq S8; reported as "Vss median" in Table 2 (159 L).
    vssb <- (kpu_lung     * v_lung +
             kpu_kidney   * v_kidney +
             kpu_heart    * v_heart +
             kpu_spleen   * v_spleen +
             kpu_liver    * v_liver * (1 - erh) +
             kpu_stomach  * v_stomach +
             kpu_gut      * v_gut +
             kpu_pancreas * v_pancreas +
             kpu_muscle   * v_muscle +
             kpu_bone     * v_bone +
             kpu_skin     * v_skin +
             kpu_brain    * v_brain +
             kpu_other    * v_other +
             kpu_adipose  * v_adipose) * fub + v_arterial + v_venous

    Cc <- cblood / bp
    Cc ~ prop(propSd)
  })
}
