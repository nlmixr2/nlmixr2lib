Yau_2023_diazepam_pbpk_kpu_rat <- function() {
  description <- paste(
    "Preclinical (rat). PBPK (14-compartment whole-body, four common Kpu",
    "clusters, k-means clustering). Diazepam disposition in rats after IV",
    "administration - the rat counterpart of Yau_2023_diazepam_pbpk_kpu_human.",
    "The full kinetic structure and every individual tissue blood flow and",
    "volume are retained; only the number of unknown distribution parameters",
    "is reduced, by grouping tissues a priori on similarity of tissue",
    "composition and assigning one common unbound tissue:plasma partition",
    "coefficient (Kpu) per group. The lung, arterial and venous blood are",
    "treated as a single quasi-steady-state 'central' compartment; the",
    "remaining 13 tissues are perfusion-limited well-stirred compartments,",
    "with stomach, gut, spleen and pancreas draining into the liver via the",
    "portal vein. Clustering here is the k-means solution, which the paper",
    "selected as one of its three primary rat candidates: the four clusters",
    "are (1) bone, brain, muscle, pancreas and rest of body, (2) kidney,",
    "spleen and liver, (3) skin, lung, gut, stomach and heart, and",
    "(4) adipose. Note this differs from the hierarchical grouping used for",
    "the human model, where lung and the splanchnic organs cluster with",
    "kidney and liver. Rat physiology is scaled from body weight using the",
    "authors' own control-stream construction, which reproduces the 250-g",
    "reference-rat values of Table S1 exactly. Because the rat hepatic",
    "clearance implied by CLb exceeds hepatic blood flow, hepatic clearance",
    "is capped at 0.99 * Q_liver as the authors' control stream and the",
    "Figure S4 footnote specify; the balance is carried as renal clearance."
  )
  reference <- paste(
    "Yau E, Olivares-Morales A, Ogungbenro K, Aarons L, Gertz M.",
    "Investigation of simplified physiologically-based pharmacokinetic",
    "models in rat and human. CPT Pharmacometrics Syst Pharmacol.",
    "2023;12(3):333-345. doi:10.1002/psp4.12911.",
    "Parameter estimates from Table S6 column '4 common Kpus (Km)';",
    "cluster membership from the Table S6 footnote; model equations from",
    "Appendix S1 Eq S17-S19 and the example NONMEM control stream; rat",
    "physiological scaling and the hepatic-clearance cap from the",
    "Appendix S1 example NONMEM control stream.",
    sep = " "
  )
  vignette <- "Yau_2023_diazepam_pbpk"
  units    <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. analyte/specimen proposed by a local model from the
  # model description; units derived from the units block. verified = FALSE
  # means NOT checked against the source paper.
  compartmentData <- list(
    central  = list(analyte = "diazepam", units = "mg", specimen = "plasma", verified = FALSE),
    adipose  = list(analyte = "diazepam", units = "mg", specimen = "tissue", verified = FALSE),
    bone     = list(analyte = "diazepam", units = "mg", specimen = "tissue", verified = FALSE),
    brain    = list(analyte = "diazepam", units = "mg", specimen = "tissue", verified = FALSE),
    heart    = list(analyte = "diazepam", units = "mg", specimen = "tissue", verified = FALSE),
    kidney   = list(analyte = "diazepam", units = "mg", specimen = "tissue", verified = FALSE),
    muscle   = list(analyte = "diazepam", units = "mg", specimen = "tissue", verified = FALSE),
    skin     = list(analyte = "diazepam", units = "mg", specimen = "tissue", verified = FALSE),
    other    = list(analyte = "diazepam", units = "mg", specimen = "tissue", verified = FALSE),
    gut      = list(analyte = "diazepam", units = "mg", specimen = "tissue", verified = FALSE),
    stomach  = list(analyte = "diazepam", units = "mg", specimen = "tissue", verified = FALSE),
    spleen   = list(analyte = "diazepam", units = "mg", specimen = "tissue", verified = FALSE),
    pancreas = list(analyte = "diazepam", units = "mg", specimen = "tissue", verified = FALSE),
    liver    = list(analyte = "diazepam", units = "mg", specimen = "tissue", verified = FALSE)
  )

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
      "variability and IIV on the Kpu parameters was deliberately not",
      "estimated so the model variants could be compared. Experimental rat",
      "Kpu values for diazepam are available in Table S7 for comparison",
      "against these estimates. Diazepam fu_p = 0.1 and blood:plasma ratio",
      "= 0.836 were measured internally (Appendix S1). Urinary excretion",
      "fraction fe = 0.009 is from the literature (paper ref 31)."
    )
  )

  ini({
    # Table S6, column "4 common Kpus (Km)". Values are listed in the
    # supplement together with their coefficient of variation in parentheses.
    lcl <- log(1.14); label("Total blood clearance (L/h)")                                               # Yau 2023 Table S6, 4 common Kpus (Km), CLb = 1.14 L/h (CV 33%)

    # Common unbound tissue:plasma partition coefficients, one per
    # tissue-composition cluster (Table S6 footnote, four-cluster k-means).
    lkpu1 <- log(29.4); label("Common Kpu, cluster 1 - bone, brain, muscle, pancreas, rest of body (unitless)")  # Yau 2023 Table S6 (29.4, CV 4%)
    lkpu2 <- log(446);  label("Common Kpu, cluster 2 - kidney, spleen, liver (unitless)")                # Yau 2023 Table S6 (446, CV 2%)
    lkpu3 <- log(89.1); label("Common Kpu, cluster 3 - skin, lung, gut, stomach, heart (unitless)")       # Yau 2023 Table S6 (89.1, CV 3%)
    lkpu4 <- log(18.2); label("Common Kpu, cluster 4 - adipose (unitless)")                              # Yau 2023 Table S6 (18.2, CV 5%)

    # Exponential IIV on clearance only (paper Data analysis). Reported as
    # 14.9% CV; back-transformed with omega^2 = log(1 + CV^2).
    etalcl ~ 0.0219575; label("IIV on total blood clearance")                                            # Yau 2023 Table S6, IIV CLp = 14.9% (CV 33%)

    # Proportional residual error (supplement $ERROR: W = IPRED*RUV2).
    propSd <- 0.200; label("Proportional residual error (fraction)")                                      # Yau 2023 Table S6, residual error 20.0% (CV 4%)
  })

  model({
    # ================= Drug-specific constants =================
    fup  <- 0.1        # fraction unbound in rat plasma (Appendix S1, measured)
    bp    <- 0.836     # blood:plasma concentration ratio, rat (Appendix S1, measured)
    fexc  <- 0.009     # fraction of dose excreted unchanged in urine (paper ref 31)
    fub   <- fup / bp  # fraction unbound in blood

    # ================= Physiology: body-weight-scaled rat =================
    # Appendix S1 example NONMEM control stream, $PK block. Cardiac output
    # in L/min (Brown, Arms and Travis 1988); converted to L/h here.
    co <- 0.235 * WT^0.75 * 60

    # Fractions of cardiac output (Brown 1997, Kuwahira 1994). HA = hepatic
    # artery; PV = portal vein draining stomach + gut + spleen + pancreas.
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
    q_liver    <- (fco_ha + fco_pv) * co   # total hepatic flow (Q_HV)
    q_co       <- co

    # Tissue densities (kg/L) and fractional body weights (Brown et al. 1997
    # Table 5). Volumes in L for WT in kg.
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
    cl   <- exp(lcl + etalcl)
    kpu1 <- exp(lkpu1)
    kpu2 <- exp(lkpu2)
    kpu3 <- exp(lkpu3)
    kpu4 <- exp(lkpu4)

    # ---- Cluster membership (Table S6 footnote, 4 common Kpus, k-means) ----
    kpu_bone     <- kpu1
    kpu_brain    <- kpu1
    kpu_muscle   <- kpu1
    kpu_pancreas <- kpu1
    kpu_other    <- kpu1
    kpu_kidney   <- kpu2
    kpu_spleen   <- kpu2
    kpu_liver    <- kpu2
    kpu_skin     <- kpu3
    kpu_lung     <- kpu3
    kpu_gut      <- kpu3
    kpu_stomach  <- kpu3
    kpu_heart    <- kpu3
    kpu_adipose  <- kpu4

    # ---- Tissue:blood partition coefficients, Eq S7: Kb = Kpu * fu_p / BP ----
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
    # For the rat diazepam estimates the implied hepatic clearance exceeds
    # hepatic blood flow, so the cap is active and a substantial part of
    # total clearance is carried as the "renal" term - a known limitation of
    # the rat fit, stated in the Figure S4 footnote.
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
    # Eq S8; reported as "Vss median" in Table S6 (1.19 L).
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
