Yau_2023_diazepam_pbpk_lumped_human <- function() {
  description <- paste(
    "PBPK (kinetically lumped 3-compartment, physiologically parameterised).",
    "Diazepam disposition in healthy adults after IV administration. One of",
    "two whole-body PBPK simplification strategies investigated by Yau et al.",
    "(2023): Nestorov kinetic lumping collapses the 16-state whole-body PBPK",
    "model into three compartments on the basis of tissue time constants",
    "(paper Table 1). The lung, arterial and venous blood are lumped into a",
    "quasi-steady-state 'central' compartment; all remaining tissues except",
    "adipose form the moderately equilibrating 'peripheral1' compartment;",
    "adipose alone forms the slowly equilibrating 'peripheral2' compartment.",
    "Micro-constants are not free parameters - they are computed from the",
    "70-kg reference-man blood flows and volumes (Table S1) and the",
    "estimated unbound tissue:plasma partition coefficients via supplement",
    "Eq S12-S16, so the model retains physiological interpretability while",
    "having the same number of free parameters as an empirical 3-compartment",
    "model. Only total blood clearance and three lumped Kpu values are",
    "estimated. Kpu1 is the lung Kpu (blood + lung central pool), Kpu2 the",
    "common Kpu of the peripheral1 tissues and Kpu3 the adipose Kpu.",
    "Elimination is split into hepatic and renal components using the",
    "well-stirred model and the reported urinary excretion fraction",
    "fe = 0.0005; because only 0.05% of the dose is excreted unchanged the",
    "renal term is negligible and clearance is effectively all hepatic."
  )
  reference <- paste(
    "Yau E, Olivares-Morales A, Ogungbenro K, Aarons L, Gertz M.",
    "Investigation of simplified physiologically-based pharmacokinetic",
    "models in rat and human. CPT Pharmacometrics Syst Pharmacol.",
    "2023;12(3):333-345. doi:10.1002/psp4.12911.",
    "Parameter estimates from Table 2 column 'Lumped 3-compartment';",
    "lumping equations from Appendix S1 Eq S9-S16; reference-man",
    "physiology from Table S1.",
    sep = " "
  )
  vignette <- "Yau_2023_diazepam_pbpk"
  units    <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Derived mechanically; verified = FALSE means it has
  # NOT been checked against the source paper.
  compartmentData <- list(
    central     = list(analyte = "diazepam pbpk lumped", units = "mg", specimen = "plasma", verified = FALSE),
    peripheral1 = list(analyte = "diazepam pbpk lumped", units = "mg", specimen = "plasma", verified = FALSE),
    peripheral2 = list(analyte = "diazepam pbpk lumped", units = "mg", specimen = "plasma", verified = FALSE)
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
      "profiles, 28 are individual profiles (4 from Table S4 study 6, 23",
      "from study 7, 1 from study 4) and 7 are study-average profiles;",
      "n_subjects = 55 is the sum of the reported subject counts across the",
      "seven studies. Because most records are arithmetic-mean profiles",
      "from different protocols and populations, the estimated 'IIV' on",
      "clearance is confounded with inter-study variability - the paper",
      "notes this explicitly and recommends individual or geometric-mean",
      "data for future analyses. Diazepam fu_p = 0.009 and blood:plasma",
      "ratio = 0.559 were measured internally by equilibrium dialysis and",
      "LC-MS/MS (Appendix S1). Urinary excretion fraction fe = 0.0005 is",
      "from the literature (paper ref 30)."
    )
  )

  ini({
    # Table 2, column "Lumped 3-compartment". Values are listed in the
    # paper together with their coefficient of variation in parentheses.
    lcl <- log(3.63); label("Total blood clearance (L/h)")                                                   # Yau 2023 Table 2, lumped 3-compartment, CLb = 3.63 L/h (CV 4%)

    # Lumped unbound tissue:plasma partition coefficients. Group membership
    # is given in the Table 2 footnote and matches the Table 1 "Compartment"
    # assignment for human and the Appendix S1 V1/V2/V3 definitions.
    lkpu1 <- log(1150); label("Kpu of the central lump - blood + lung (unitless)")                           # Yau 2023 Table 2 (1150, CV 5%); footnote "Kpu1 corresponds to the lumped compartment (blood and lungs)"
    lkpu2 <- log(24.3); label("Kpu of the peripheral1 lump - all tissues except adipose (unitless)")          # Yau 2023 Table 2 (24.3, CV 10%); footnote lists muscle, bone, skin, brain, rest of body, kidneys, heart, spleen, liver, pancreas, gut, stomach
    lkpu3 <- log(483);  label("Kpu of the peripheral2 lump - adipose (unitless)")                             # Yau 2023 Table 2 (483, CV 2%); footnote "Kpu3 to (adipose)"

    # Exponential IIV on clearance only (paper Data analysis: "An
    # exponential model was used to account for interindividual
    # variability"). Reported as 33.1% CV; back-transformed with the
    # log-normal relation omega^2 = log(1 + CV^2) = log(1 + 0.331^2).
    etalcl ~ 0.104259; label("IIV on total blood clearance")                                                 # Yau 2023 Table 2, IIV CLb = 33.1% (CV 17%)

    # Proportional residual error (supplement $ERROR: W = IPRED*RUV2).
    propSd <- 0.393; label("Proportional residual error (fraction)")                                          # Yau 2023 Table 2, residual error 39.3% (CV 8%)
  })

  model({
    # ================= Drug-specific constants =================
    # Measured internally by equilibrium dialysis / LC-MS/MS for human
    # plasma and blood (paper "Diazepam data" and Appendix S1
    # "Determination of plasma protein binding and blood to plasma
    # partitioning ratio").
    fup  <- 0.009      # fraction unbound in human plasma
    bp    <- 0.559     # blood:plasma concentration ratio, human
    fexc  <- 0.0005    # fraction of dose excreted unchanged in urine (paper ref 30)
    fub   <- fup / bp  # fraction unbound in blood

    # ================= Physiology: 70-kg reference man =================
    # Table S1 gives blood flows in L/min and volumes in L; flows are
    # converted to L/h to match units$time = "h".
    q_adipose  <- 0.292 * 60
    q_bone     <- 0.292 * 60
    q_brain    <- 0.701 * 60
    q_gut      <- 0.526 * 60
    q_heart    <- 0.234 * 60
    q_kidney   <- 1.109 * 60
    q_liver    <- 1.489 * 60   # total hepatic (hepatic artery + portal vein)
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

    # Portal vein drains stomach, gut, spleen and pancreas into the liver;
    # the hepatic artery carries the balance of the total hepatic flow
    # (Appendix S1 $PK: Q_HV = Q_HA + Q_PV, with the portal pool defined as
    # stomach + gut + spleen + pancreas).
    q_pv <- q_gut + q_stomach + q_spleen + q_pancreas
    q_ha <- q_liver - q_pv

    # Rest-of-body blood flow. Table S1 reports 0.730 L/min, but the sum of
    # all reported arterial outflows then exceeds the reported cardiac
    # output by 5%, which would make the model non-mass-conservative.
    # Q_other is therefore taken as the residual that closes the balance,
    # exactly as the authors' own rat control stream does
    # (FCO_ROB = 1 - sum(all other fractions), Appendix S1 $PK). Only the
    # hepatic-artery share of the hepatic flow is subtracted here, because
    # the splanchnic organs are perfused directly from arterial blood and
    # are already counted individually. This yields 0.437 L/min and is
    # independently confirmed by the Appendix S1 lumped-model constraint
    # Q1 = Q2 + Q3 = cardiac output, which only holds for 0.437. See the
    # vignette "Assumptions and deviations".
    q_other <- q_co - (q_adipose + q_bone + q_brain + q_gut + q_heart +
                       q_kidney + q_ha + q_muscle + q_pancreas +
                       q_skin + q_spleen + q_stomach)
    v_other <- 2.914

    # ================= Individual parameters =================
    cl   <- exp(lcl + etalcl)
    kpu1 <- exp(lkpu1)
    kpu2 <- exp(lkpu2)
    kpu3 <- exp(lkpu3)

    # Blood tissue:plasma partition coefficients, Eq S7: Kb = Kpu * fu_p / BP
    kb_central <- kpu1 * fub
    kb_p1      <- kpu2 * fub
    kb_p2      <- kpu3 * fub

    # ================= Lumped compartment flows and volumes =================
    # Appendix S1, "For the lumped model in man":
    #   Q2 = Q_muscle + Q_bone + Q_skin + Q_brain + Q_RoB + Q_liver +
    #        Q_heart + Q_kidney       (the supplement writes the
    #        hepato-splanchnic term as "Q_portal vein"; only the total
    #        hepatic flow Q_liver makes Q1 = Q2 + Q3 equal cardiac output)
    #   Q3 = Q_adipose
    #   V1 = V_lung + V_arterial + V_venous
    #   V2 = V_muscle + V_bone + V_skin + V_brain + V_RoB + V_liver +
    #        V_stomach + V_gut + V_pancreas + V_spleen + V_heart + V_kidney
    #   V3 = V_adipose
    q_p1 <- q_muscle + q_bone + q_skin + q_brain + q_other + q_liver +
            q_heart + q_kidney
    q_p2 <- q_adipose
    v_p1 <- v_muscle + v_bone + v_skin + v_brain + v_other + v_liver +
            v_stomach + v_gut + v_pancreas + v_spleen + v_heart + v_kidney
    v_p2 <- v_adipose

    # ================= Clearance partitioning =================
    clh <- cl * (1 - fexc)   # hepatic blood clearance
    clr <- cl * fexc         # renal blood clearance

    # ================= Effective (Kb-weighted) volumes =================
    # V_CONV is the central-compartment volume term of Eq S12/S15/S16:
    # (V1 - V_art - V_ven) * Kb1 + V_art + V_ven = V_lung * Kb1 + V_art + V_ven.
    v_conv <- v_arterial + v_venous + v_lung * kb_central

    # Eq S13 denominator: the peripheral1 volume is reduced by the
    # extracted fraction of the eliminating tissues, i.e. weighted by
    # sum(Vi * (1 - Ei)) / sum(Vi).
    v_p1_eff <- kb_p1 * (v_p1 - v_liver * clh / q_liver -
                                v_kidney * clr / q_kidney)
    v_p2_eff <- kb_p2 * v_p2

    # ================= Micro-constants (Eq S12-S16) =================
    k12 <- q_p1 / v_conv
    k21 <- q_p1 / v_p1_eff
    k13 <- q_p2 / v_conv
    k31 <- q_p2 / v_p2_eff
    kel <- cl   / v_conv

    # ================= ODEs =================
    d/dt(central)     <- k21 * peripheral1 + k31 * peripheral2 -
                         (k12 + k13 + kel) * central
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2

    # ================= Observation =================
    # The central state is a blood pool; plasma concentration follows by
    # dividing by the blood:plasma ratio (supplement $ERROR: CPL = CBL/BP).
    cblood <- central / v_conv

    # Volume of distribution at steady state on a whole-blood basis
    # (Eq S8 / Eq 9), reported as "Vss median" in Table 2 (154 L).
    vssb <- v_arterial + v_venous + kb_central * v_lung + v_p1_eff + v_p2_eff

    Cc <- cblood / bp
    Cc ~ prop(propSd)
  })
}
