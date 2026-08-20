Yau_2023_diazepam_pbpk_kpu_human <- function() {
  description <- paste(
    "PBPK (14-compartment whole-body, four common Kpu clusters). Diazepam",
    "disposition in healthy adults after IV administration. The second of",
    "two whole-body PBPK simplification strategies investigated by Yau et",
    "al. (2023): rather than lumping compartments, the full kinetic",
    "structure and every individual tissue blood flow and volume are",
    "retained, and only the number of unknown distribution parameters is",
    "reduced. Tissues are grouped a priori by similarity of tissue",
    "composition using hierarchical clustering (drug-independent), and each",
    "cluster is assigned one common unbound tissue:plasma partition",
    "coefficient (Kpu). This reduces the distribution parameters from 13 to",
    "4 while leaving the kinetic structure of the whole-body model intact,",
    "so individual tissue concentration-time profiles remain available. The",
    "lung, arterial and venous blood are treated as a single",
    "quasi-steady-state 'central' compartment (the only kinetic assumption);",
    "the remaining 13 tissues are perfusion-limited well-stirred",
    "compartments. Stomach, gut, spleen and pancreas drain into the liver",
    "via the portal vein; hepatic elimination uses the well-stirred model",
    "and renal elimination is placed on the central compartment. The four",
    "clusters are (1) bone, brain, muscle, pancreas and rest of body,",
    "(2) lung, gut, stomach, kidney, heart, spleen and liver, (3) skin and",
    "(4) adipose. Together with the four-scalar variant this was one of the",
    "two models the paper found best described the human diazepam data."
  )
  reference <- paste(
    "Yau E, Olivares-Morales A, Ogungbenro K, Aarons L, Gertz M.",
    "Investigation of simplified physiologically-based pharmacokinetic",
    "models in rat and human. CPT Pharmacometrics Syst Pharmacol.",
    "2023;12(3):333-345. doi:10.1002/psp4.12911.",
    "Parameter estimates from Table 2 column '4 common Kpus'; cluster",
    "membership from the Table 2 footnote; model equations from",
    "Appendix S1 Eq S17-S19 and the example NONMEM control stream;",
    "reference-man physiology from Table S1.",
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
      "internally by equilibrium dialysis and LC-MS/MS (Appendix S1).",
      "Urinary excretion fraction fe = 0.0005 is from the literature",
      "(paper ref 30)."
    )
  )

  ini({
    # Table 2, column "4 common Kpus". Values are listed in the paper
    # together with their coefficient of variation in parentheses.
    lcl <- log(3.67); label("Total blood clearance (L/h)")                                              # Yau 2023 Table 2, 4 common Kpus, CLb = 3.67 L/h (CV 4%)

    # Common unbound tissue:plasma partition coefficients, one per
    # tissue-composition cluster (Table 2 footnote, four-cluster
    # hierarchical solution).
    lkpu1 <- log(32.4); label("Common Kpu, cluster 1 - bone, brain, muscle, pancreas, rest of body (unitless)")  # Yau 2023 Table 2 (32.4, CV 6%)
    lkpu2 <- log(89.1); label("Common Kpu, cluster 2 - lung, gut, stomach, kidney, heart, spleen, liver (unitless)")  # Yau 2023 Table 2 (89.1, CV 10%)
    lkpu3 <- log(324);  label("Common Kpu, cluster 3 - skin (unitless)")                                # Yau 2023 Table 2 (324, CV 5%)
    lkpu4 <- log(483);  label("Common Kpu, cluster 4 - adipose (unitless)")                             # Yau 2023 Table 2 (483, CV 8%)

    # Exponential IIV on clearance only (paper Data analysis). Reported as
    # 32.1% CV; back-transformed with omega^2 = log(1 + CV^2).
    etalcl ~ 0.0980834; label("IIV on total blood clearance")                                           # Yau 2023 Table 2, IIV CLb = 32.1% (CV 19%)

    # Proportional residual error (supplement $ERROR: W = IPRED*RUV2).
    propSd <- 0.386; label("Proportional residual error (fraction)")                                     # Yau 2023 Table 2, residual error 38.6% (CV 8%)
  })

  model({
    # ================= Drug-specific constants =================
    fup  <- 0.009      # fraction unbound in human plasma (Appendix S1, measured)
    bp    <- 0.559     # blood:plasma concentration ratio, human (Appendix S1, measured)
    fexc  <- 0.0005    # fraction of dose excreted unchanged in urine (paper ref 30)
    fub   <- fup / bp  # fraction unbound in blood

    # ================= Physiology: 70-kg reference man (Table S1) =========
    # Blood flows reported in L/min, converted to L/h; volumes in L.
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
    q_co       <- 5.839 * 60   # cardiac output = Table S1 lung/arterial/venous flow

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

    # Portal vein pool (stomach + gut + spleen + pancreas) and hepatic
    # artery (Appendix S1 $PK: Q_HV = Q_HA + Q_PV).
    q_pv <- q_gut + q_stomach + q_spleen + q_pancreas
    q_ha <- q_liver - q_pv

    # Rest-of-body flow taken as the residual that closes the mass balance
    # against cardiac output; see the vignette "Assumptions and deviations".
    q_other <- q_co - (q_adipose + q_bone + q_brain + q_gut + q_heart +
                       q_kidney + q_ha + q_muscle + q_pancreas +
                       q_skin + q_spleen + q_stomach)

    # ================= Individual parameters =================
    cl   <- exp(lcl + etalcl)
    kpu1 <- exp(lkpu1)
    kpu2 <- exp(lkpu2)
    kpu3 <- exp(lkpu3)
    kpu4 <- exp(lkpu4)

    # ---- Cluster membership (Table 2 footnote, 4 common Kpus) ----
    kpu_bone     <- kpu1
    kpu_brain    <- kpu1
    kpu_muscle   <- kpu1
    kpu_pancreas <- kpu1
    kpu_other    <- kpu1
    kpu_lung     <- kpu2
    kpu_gut      <- kpu2
    kpu_stomach  <- kpu2
    kpu_kidney   <- kpu2
    kpu_heart    <- kpu2
    kpu_spleen   <- kpu2
    kpu_liver    <- kpu2
    kpu_skin     <- kpu3
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
    # Hepatic clearance is capped at 0.99 * hepatic blood flow so the
    # well-stirred inversion for intrinsic clearance stays finite; the
    # balance is carried as renal clearance. For the human diazepam
    # estimates the cap is inactive (CL_H is ~4% of hepatic flow).
    clh   <- min(cl * (1 - fexc), 0.99 * q_liver)
    clr   <- cl - clh
    clint <- clh * q_liver / ((q_liver - clh) * fub)   # Eq S6 rearranged
    erh   <- clh / q_liver                             # hepatic extraction ratio

    # ================= Central compartment volume term ===================
    # V_CONV = V_arterial + V_venous + V_lung * Kb_lung, i.e. the product of
    # the central volume and its composite partition coefficient
    # (Kb_central = (V_art + V_ven + V_lung*Kb_lung)/(V_art + V_ven + V_lung)).
    v_conv <- v_arterial + v_venous + v_lung * kb_lung
    cblood <- central / v_conv

    # ================= ODEs (Eq S17-S19) =================================
    # Central (Eq S17): receives venous return from every tissue that drains
    # directly to blood - i.e. all except stomach, gut, spleen and pancreas,
    # which drain into the liver - and loses drug to cardiac output and to
    # renal clearance.
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

    # Non-eliminating perfusion-limited tissues (Eq S18)
    d/dt(adipose)  <- q_adipose  * cblood - q_adipose  * (adipose  / v_adipose  / kb_adipose)
    d/dt(bone)     <- q_bone     * cblood - q_bone     * (bone     / v_bone     / kb_bone)
    d/dt(brain)    <- q_brain    * cblood - q_brain    * (brain    / v_brain    / kb_brain)
    d/dt(heart)    <- q_heart    * cblood - q_heart    * (heart    / v_heart    / kb_heart)
    d/dt(kidney)   <- q_kidney   * cblood - q_kidney   * (kidney   / v_kidney   / kb_kidney)
    d/dt(muscle)   <- q_muscle   * cblood - q_muscle   * (muscle   / v_muscle   / kb_muscle)
    d/dt(skin)     <- q_skin     * cblood - q_skin     * (skin     / v_skin     / kb_skin)
    d/dt(other)    <- q_other    * cblood - q_other    * (other    / v_other    / kb_other)

    # Splanchnic tissues draining into the portal vein (Eq S18)
    d/dt(gut)      <- q_gut      * cblood - q_gut      * (gut      / v_gut      / kb_gut)
    d/dt(stomach)  <- q_stomach  * cblood - q_stomach  * (stomach  / v_stomach  / kb_stomach)
    d/dt(spleen)   <- q_spleen   * cblood - q_spleen   * (spleen   / v_spleen   / kb_spleen)
    d/dt(pancreas) <- q_pancreas * cblood - q_pancreas * (pancreas / v_pancreas / kb_pancreas)

    # Liver (Eq S19): hepatic-artery inflow plus portal inflow from the four
    # splanchnic organs, minus hepatic venous outflow and intrinsic
    # metabolic clearance acting on unbound drug.
    d/dt(liver) <-
        q_ha * cblood +
        q_gut      * (gut      / v_gut      / kb_gut) +
        q_stomach  * (stomach  / v_stomach  / kb_stomach) +
        q_spleen   * (spleen   / v_spleen   / kb_spleen) +
        q_pancreas * (pancreas / v_pancreas / kb_pancreas) -
        (liver / v_liver) * ((q_liver + clint * fub) / kb_liver)

    # ================= Derived quantities and observation ================
    # Volume of distribution at steady state on a whole-blood basis (Eq S8);
    # only the eliminating liver is weighted by (1 - ER). Reported as
    # "Vss median" in Table 2 (159 L).
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

    Cc <- cblood / bp   # plasma concentration (supplement $ERROR: CPL = CBL/BP)
    Cc ~ prop(propSd)
  })
}
