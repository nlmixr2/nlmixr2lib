Yau_2023_diazepam_pbpk_scalar_rat <- function() {
  description <- paste(
    "Preclinical (rat). PBPK (14-compartment whole-body, four common Kpu",
    "scalars, hierarchical clustering). Diazepam disposition in rats after IV",
    "administration - the rat counterpart of",
    "Yau_2023_diazepam_pbpk_scalar_human. The full kinetic structure and",
    "every individual tissue blood flow and volume are retained; tissue Kpu",
    "values are first predicted bottom-up with the Rodgers and Rowland (R&R)",
    "equations from rat tissue composition and diazepam physicochemistry, and",
    "then corrected by a small number of estimated multiplicative scaling",
    "factors shared within composition-based clusters",
    "(Eq 10: Kpu_i = Kpu_predRR,i * SF_i). The scalar therefore quantifies",
    "the systematic bias of the bottom-up R&R prediction rather than",
    "replacing it, which is why the paper judged this parameterisation the",
    "most promising for interspecies translation. The R&R weak-base branch is",
    "used because diazepam's basic pKa is 3.4 (< 7). The lung, arterial and",
    "venous blood are treated as a single quasi-steady-state 'central'",
    "compartment; the remaining 13 tissues are perfusion-limited well-stirred",
    "compartments, with stomach, gut, spleen and pancreas draining into the",
    "liver via the portal vein. Clustering here is the hierarchical solution,",
    "which the paper selected as one of its three primary rat candidates: the",
    "four clusters are (1) bone, brain, muscle, pancreas and rest of body,",
    "(2) lung, gut, stomach, kidney, heart, spleen and liver, (3) skin and",
    "(4) adipose. Rat physiology is scaled from body weight using the",
    "authors' own control-stream construction, which reproduces the 250-g",
    "reference-rat values of Table S1 exactly. Hepatic clearance is capped at",
    "0.99 * Q_liver as the control stream and the Figure S4 footnote specify."
  )
  reference <- paste(
    "Yau E, Olivares-Morales A, Ogungbenro K, Aarons L, Gertz M.",
    "Investigation of simplified physiologically-based pharmacokinetic",
    "models in rat and human. CPT Pharmacometrics Syst Pharmacol.",
    "2023;12(3):333-345. doi:10.1002/psp4.12911.",
    "Parameter estimates from Table S6 column '4 scalars (H)'; cluster",
    "membership from the Table S6 footnote; scalar definition from Eq 10;",
    "Rodgers and Rowland equations and rat tissue composition from the",
    "Appendix S1 example NONMEM control stream and Table S3; model",
    "equations from Appendix S1 Eq S17-S19; rat physiological scaling and",
    "the hepatic-clearance cap from the Appendix S1 example NONMEM",
    "control stream.",
    sep = " "
  )
  vignette <- "Yau_2023_diazepam_pbpk"
  units    <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Drives every rat physiological constant. Cardiac output is",
        "0.235 * WT^0.75 L/min (Brown, Arms and Travis 1988, cited in the",
        "Appendix S1 control stream); tissue blood flows are fixed",
        "fractions of cardiac output and tissue volumes are fixed",
        "fractional body weights divided by tissue densities (Brown et al.",
        "1997 Table 5). At WT = 0.25 kg this reproduces the 250-g",
        "reference-rat blood flows and volumes of Table S1 to three",
        "decimal places. Rat doses in the source studies are expressed in",
        "mg/kg (Table S5), so WT is also needed to convert to the mg dose",
        "amount used by this model."
      ),
      source_name        = "BW"
    )
  )

  population <- list(
    species        = "rat (male Wistar and Sprague-Dawley)",
    n_subjects     = 51L,
    n_studies      = 5L,
    age_range      = paste(
      "adult; Table S5 study 4 contrasted middle-aged vs old rats",
      "(clearance 1.1 vs 3.3 L/h)"
    ),
    weight_range   = "250-g reference rat used for the physiological constants (Table S1)",
    sex_female_pct = 0,
    disease_state  = "healthy (normal) rats",
    dose_range     = paste(
      "IV diazepam 1.2-5 mg/kg as a bolus, or 1 mg as a 5-min infusion",
      "across the five contributing studies (Table S5)"
    ),
    regions        = NA_character_,
    notes          = paste(
      "6 rat IV concentration-time profiles pooled from five published",
      "studies (Table S5), all study-average profiles, digitised with",
      "WebPlotDigitizer 4.2; n_subjects = 51 is the sum of the reported",
      "animal counts across the five studies. Strains are male Wistar",
      "(studies 1, 2, 4, 5) and male Sprague-Dawley (study 3). The rat data",
      "are described by the paper as 'more limited and variable' than the",
      "human data; IIV on clearance is confounded with inter-study",
      "variability and IIV on the scaling factors was deliberately not",
      "estimated so the model variants could be compared. Experimental rat",
      "Kpu values for diazepam and the corresponding R&R predictions are",
      "available in Table S7 for comparison. Diazepam fu_p = 0.1 and",
      "blood:plasma ratio = 0.836 were measured internally (Appendix S1);",
      "logP = 2.82 and pKa = 3.4 are from the literature (paper ref 29).",
      "Urinary excretion fraction fe = 0.009 is from the literature",
      "(paper ref 31)."
    )
  )

  ini({
    # Table S6, column "4 scalars (H)". Values are listed in the supplement
    # together with their coefficient of variation in parentheses.
    lcl <- log(1.19); label("Total blood clearance (L/h)")                                               # Yau 2023 Table S6, 4 scalars (H), CLb = 1.19 L/h (CV 31%)

    # Multiplicative correction factors on the Rodgers and Rowland predicted
    # Kpu, one per tissue-composition cluster (Eq 10; Table S6 footnote,
    # four-cluster hierarchical solution).
    lsf1 <- log(3.25);  label("Kpu scaling factor, cluster 1 - bone, brain, muscle, pancreas, rest of body (unitless)")  # Yau 2023 Table S6 (3.25, CV 10%)
    lsf2 <- log(21.5);  label("Kpu scaling factor, cluster 2 - lung, gut, stomach, kidney, heart, spleen, liver (unitless)")  # Yau 2023 Table S6 (21.5, CV 7%)
    lsf3 <- log(2.07);  label("Kpu scaling factor, cluster 3 - skin (unitless)")                         # Yau 2023 Table S6 (2.07, CV 20%)
    lsf4 <- log(0.204); label("Kpu scaling factor, cluster 4 - adipose (unitless)")                      # Yau 2023 Table S6 (0.204, CV 6%)

    # Exponential IIV on clearance only (paper Data analysis). Reported as
    # 15.1% CV; back-transformed with omega^2 = log(1 + CV^2).
    etalcl ~ 0.0225448; label("IIV on total blood clearance")                                            # Yau 2023 Table S6, IIV CLp = 15.1% (CV 31%)

    # Proportional residual error (supplement $ERROR: W = IPRED*RUV2).
    propSd <- 0.219; label("Proportional residual error (fraction)")                                      # Yau 2023 Table S6, residual error 21.9% (CV 5%)
  })

  model({
    # ================= Drug-specific constants =================
    fup  <- 0.1        # fraction unbound in rat plasma (Appendix S1, measured)
    bp    <- 0.836     # blood:plasma concentration ratio, rat (Appendix S1, measured)
    fexc  <- 0.009     # fraction of dose excreted unchanged in urine (paper ref 31)
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

    # Affinity for extracellular albumin, from rat plasma neutral lipid and
    # neutral phospholipid fractions (Table S3 rat plasma row: vNL = 0.0023,
    # vNP = 0.0013).
    kapr <- (1 / fup) - 1 -
      ((pow * 0.0023 + ((0.3 * pow) + 0.7) * 0.0013) / (1 + fi2))

    # Per-tissue predicted Kpu. Numbers in each expression are, in order,
    # fIW, fEW, AR, vNL and vNP from Table S3 (rat tissue composition).
    kpurr_adipose  <- ((1 + fi1) * 0.017) / (1 + fi2) + 0.135 + kapr * 0.049 +
      ((kvo * 0.853 + ((0.3 * kvo) + 0.7) * 0.0016) / (1 + fi2))
    kpurr_bone     <- ((1 + fi1) * 0.346) / (1 + fi2) + 0.100 + kapr * 0.100 +
      ((pow * 0.0174 + ((0.3 * pow) + 0.7) * 0.0016) / (1 + fi2))
    kpurr_brain    <- ((1 + fi1) * 0.620) / (1 + fi2) + 0.162 + kapr * 0.048 +
      ((pow * 0.0391 + ((0.3 * pow) + 0.7) * 0.0015) / (1 + fi2))
    kpurr_gut      <- ((1 + fi1) * 0.475) / (1 + fi2) + 0.282 + kapr * 0.158 +
      ((pow * 0.0375 + ((0.3 * pow) + 0.7) * 0.0124) / (1 + fi2))
    kpurr_heart    <- ((1 + fi1) * 0.456) / (1 + fi2) + 0.320 + kapr * 0.157 +
      ((pow * 0.0135 + ((0.3 * pow) + 0.7) * 0.0106) / (1 + fi2))
    kpurr_kidney   <- ((1 + fi1) * 0.483) / (1 + fi2) + 0.273 + kapr * 0.130 +
      ((pow * 0.0121 + ((0.3 * pow) + 0.7) * 0.0240) / (1 + fi2))
    kpurr_liver    <- ((1 + fi1) * 0.573) / (1 + fi2) + 0.161 + kapr * 0.086 +
      ((pow * 0.0135 + ((0.3 * pow) + 0.7) * 0.0238) / (1 + fi2))
    kpurr_lung     <- ((1 + fi1) * 0.446) / (1 + fi2) + 0.336 + kapr * 0.212 +
      ((pow * 0.0215 + ((0.3 * pow) + 0.7) * 0.0123) / (1 + fi2))
    kpurr_muscle   <- ((1 + fi1) * 0.630) / (1 + fi2) + 0.118 + kapr * 0.064 +
      ((pow * 0.01 + ((0.3 * pow) + 0.7) * 0.0072) / (1 + fi2))
    kpurr_pancreas <- ((1 + fi1) * 0.664) / (1 + fi2) + 0.120 + kapr * 0.060 +
      ((pow * 0.0403 + ((0.3 * pow) + 0.7) * 0.009) / (1 + fi2))
    kpurr_skin     <- ((1 + fi1) * 0.291) / (1 + fi2) + 0.382 + kapr * 0.277 +
      ((pow * 0.0603 + ((0.3 * pow) + 0.7) * 0.0044) / (1 + fi2))
    kpurr_spleen   <- ((1 + fi1) * 0.579) / (1 + fi2) + 0.207 + kapr * 0.097 +
      ((pow * 0.0071 + ((0.3 * pow) + 0.7) * 0.0107) / (1 + fi2))
    # Table S3 has no stomach row; the control stream uses the gut Kpu for
    # stomach (DADT(11) uses KPU_GU) and the muscle Kpu for the rest of body
    # ("Assume same Kpu for Rest of body as for Muscle": KPU_RO = KPU_MU).
    kpurr_stomach  <- kpurr_gut
    kpurr_other    <- kpurr_muscle

    # ================= Physiology: body-weight-scaled rat =================
    co <- 0.235 * WT^0.75 * 60

    fco_adipose  <- 0.07
    fco_bone     <- 0.122
    fco_brain    <- 0.02
    fco_gut      <- 0.110
    fco_stomach  <- 0.013
    fco_heart    <- 0.049
    fco_kidney   <- 0.141
    fco_ha       <- 0.024
    fco_pv       <- 0.151
    fco_muscle   <- 0.278
    fco_pancreas <- 0.018
    fco_skin     <- 0.058
    fco_spleen   <- 0.010
    fco_other <- 1 - (fco_adipose + fco_bone + fco_brain + fco_heart +
                      fco_kidney + fco_muscle + fco_skin + fco_gut +
                      fco_stomach + fco_spleen + fco_pancreas + fco_ha)

    q_adipose  <- fco_adipose  * co
    q_bone     <- fco_bone     * co
    q_brain    <- fco_brain    * co
    q_gut      <- fco_gut      * co
    q_stomach  <- fco_stomach  * co
    q_heart    <- fco_heart    * co
    q_kidney   <- fco_kidney   * co
    q_muscle   <- fco_muscle   * co
    q_pancreas <- fco_pancreas * co
    q_skin     <- fco_skin     * co
    q_spleen   <- fco_spleen   * co
    q_other    <- fco_other    * co
    q_ha       <- fco_ha       * co
    q_liver    <- (fco_ha + fco_pv) * co
    q_co       <- co

    den_tis <- 1.040
    v_adipose  <- 0.076   * WT / 0.916
    v_bone     <- 0.04148 * WT / 1.990
    v_brain    <- 0.0057  * WT / 1.036
    v_gut      <- 0.0224  * WT / 1.043
    v_stomach  <- 0.0046  * WT / 1.050
    v_heart    <- 0.0033  * WT / 1.030
    v_kidney   <- 0.0073  * WT / 1.050
    v_liver    <- 0.0366  * WT / 1.080
    v_lung     <- 0.005   * WT / 1.051
    v_muscle   <- 0.4043  * WT / 1.041
    v_pancreas <- 0.0032  * WT / 1.045
    v_skin     <- 0.1903  * WT / 1.116
    v_spleen   <- 0.002   * WT / 1.054
    v_arterial <- 0.0272  * WT / den_tis
    v_venous   <- 0.0544  * WT / den_tis
    fw_other <- 1 - (0.076 + 0.04148 + 0.0057 + 0.0224 + 0.0046 + 0.0033 +
                     0.0073 + 0.0366 + 0.005 + 0.4043 + 0.0032 + 0.1903 +
                     0.002 + 0.0272 + 0.0544)
    v_other <- fw_other * WT / den_tis

    # ================= Individual parameters =================
    cl  <- exp(lcl + etalcl)
    sf1 <- exp(lsf1)
    sf2 <- exp(lsf2)
    sf3 <- exp(lsf3)
    sf4 <- exp(lsf4)

    # ---- Eq 10: Kpu_i = Kpu_predRR,i * SF_cluster(i) ----
    # Cluster membership per the Table S6 footnote (4 scalars, hierarchical).
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
    # Eq S8; reported as "Vss median" in Table S6 (1.68 L).
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
