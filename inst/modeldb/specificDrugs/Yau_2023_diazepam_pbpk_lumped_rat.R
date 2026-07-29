Yau_2023_diazepam_pbpk_lumped_rat <- function() {
  description <- paste(
    "Preclinical (rat). PBPK (kinetically lumped 3-compartment,",
    "physiologically parameterised). Diazepam disposition in rats after IV",
    "administration - the rat counterpart of",
    "Yau_2023_diazepam_pbpk_lumped_human. Nestorov kinetic lumping collapses",
    "the 16-state whole-body PBPK model into three compartments on the basis",
    "of tissue time constants (paper Table 1). The lung, arterial and venous",
    "blood are lumped into a quasi-steady-state 'central' compartment; all",
    "remaining tissues except skin form the moderately equilibrating",
    "'peripheral1' compartment; skin alone forms the slowly equilibrating",
    "'peripheral2' compartment. Note the lump differs from the human model:",
    "skin is the slow compartment in rat (skin is a much larger fraction of",
    "body weight in rat) whereas adipose is the slow compartment in human,",
    "which is precisely the cross-species scaling limitation the paper",
    "discusses. Micro-constants are computed from body-weight-scaled blood",
    "flows and volumes via supplement Eq S12-S16, so only total blood",
    "clearance and three lumped Kpu values are estimated. Rat physiology is",
    "scaled from body weight using the authors' own control-stream",
    "construction (cardiac output = 0.235 * BW^0.75, fractional cardiac",
    "outputs and fractional tissue weights), which reproduces the 250-g",
    "reference-rat values of Table S1 exactly. Because the rat hepatic",
    "clearance implied by CLb exceeds hepatic blood flow, hepatic clearance",
    "is capped at 0.99 * Q_liver exactly as the authors' control stream and",
    "Figure S4 footnote specify; the balance is carried as renal clearance."
  )
  reference <- paste(
    "Yau E, Olivares-Morales A, Ogungbenro K, Aarons L, Gertz M.",
    "Investigation of simplified physiologically-based pharmacokinetic",
    "models in rat and human. CPT Pharmacometrics Syst Pharmacol.",
    "2023;12(3):333-345. doi:10.1002/psp4.12911.",
    "Parameter estimates from Table S6 column 'Lumped 3-comp'; lumping",
    "equations from Appendix S1 Eq S9-S16; rat physiological scaling and",
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
      "(studies 1, 2, 4, 5) and male Sprague-Dawley (study 3). The rat",
      "data are described by the paper as 'more limited and variable' than",
      "the human data, which rendered the analysis less stable and precise;",
      "IIV on clearance is confounded with inter-study variability and IIV",
      "on the Kpu parameters was deliberately not estimated so that the",
      "model variants could be compared. Diazepam fu_p = 0.1 and",
      "blood:plasma ratio = 0.836 were measured internally by equilibrium",
      "dialysis and LC-MS/MS (Appendix S1). Urinary excretion fraction",
      "fe = 0.009 is from the literature (paper ref 31)."
    )
  )

  ini({
    # Table S6, column "Lumped 3-comp". Values are listed in the supplement
    # together with their coefficient of variation in parentheses.
    lcl <- log(1.14); label("Total blood clearance (L/h)")                                        # Yau 2023 Table S6, Lumped 3-comp, CLb = 1.14 L/h (CV 47%)

    # Lumped unbound tissue:plasma partition coefficients. Group membership
    # follows the Appendix S1 rat V1/V2/V3 definitions - see the vignette
    # "Assumptions and deviations" for the conflict with the Table S6
    # footnote and Table 1 rank column, and how it was resolved.
    lkpu1 <- log(73.0); label("Kpu of the central lump - blood + lung (unitless)")                # Yau 2023 Table S6 (73.0, CV 54%)
    lkpu2 <- log(30.6); label("Kpu of the peripheral1 lump - all tissues except skin (unitless)")  # Yau 2023 Table S6 (30.6, CV 0.3%)
    lkpu3 <- log(90.0); label("Kpu of the peripheral2 lump - skin (unitless)")                     # Yau 2023 Table S6 (90.0, CV 3%)

    # Exponential IIV on clearance only (paper Data analysis). Reported as
    # 15.1% CV; back-transformed with omega^2 = log(1 + CV^2).
    etalcl ~ 0.0225448; label("IIV on total blood clearance")                                     # Yau 2023 Table S6, IIV CLp = 15.1 (CV 33%)

    # Proportional residual error (supplement $ERROR: W = IPRED*RUV2).
    propSd <- 0.202; label("Proportional residual error (fraction)")                               # Yau 2023 Table S6, residual error 20.2% (CV 9%)
  })

  model({
    # ================= Drug-specific constants =================
    # Measured internally by equilibrium dialysis / LC-MS/MS for rat plasma
    # and blood (Appendix S1 "Determination of plasma protein binding and
    # blood to plasma partitioning ratio").
    fup  <- 0.1        # fraction unbound in rat plasma
    bp    <- 0.836     # blood:plasma concentration ratio, rat
    fexc  <- 0.009     # fraction of dose excreted unchanged in urine (paper ref 31)
    fub   <- fup / bp  # fraction unbound in blood

    # ================= Physiology: body-weight-scaled rat =================
    # Appendix S1 example NONMEM control stream, $PK block. Cardiac output
    # in L/min (Brown, Arms and Travis 1988); converted to L/h here.
    co <- 0.235 * WT^0.75 * 60

    # Fractions of cardiac output supplying each tissue (Brown 1997,
    # Kuwahira 1994). HA = hepatic artery; PV = portal vein draining
    # stomach + gut + spleen + pancreas.
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
    # Rest of body takes the residual so the flows sum to cardiac output.
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
    q_liver    <- (fco_ha + fco_pv) * co   # total hepatic flow (Q_HV)
    q_co       <- co

    # Tissue densities (kg/L) and fractional body weights (Brown et al.
    # 1997 Table 5). Volumes in L for WT in kg.
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
    # Rest of body takes the residual fractional body weight.
    fw_other <- 1 - (0.076 + 0.04148 + 0.0057 + 0.0224 + 0.0046 + 0.0033 +
                     0.0073 + 0.0366 + 0.005 + 0.4043 + 0.0032 + 0.1903 +
                     0.002 + 0.0272 + 0.0544)
    v_other <- fw_other * WT / den_tis

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
    # Appendix S1, "For the lumped model in rat":
    #   Q2 = Q_muscle + Q_bone + Q_adipose + Q_brain + Q_RoB + Q_liver +
    #        Q_heart + Q_kidney      (the supplement writes the
    #        hepato-splanchnic term as "Q_portal vein"; only the total
    #        hepatic flow Q_liver makes Q1 = Q2 + Q3 equal cardiac output)
    #   Q3 = Q_skin
    #   V1 = V_lung + V_arterial + V_venous
    #   V2 = V_muscle + V_bone + V_adipose + V_brain + V_RoB + V_liver +
    #        V_stomach + V_gut + V_pancreas + V_spleen + V_heart + V_kidney
    #   V3 = V_skin
    q_p1 <- q_muscle + q_bone + q_adipose + q_brain + q_other + q_liver +
            q_heart + q_kidney
    q_p2 <- q_skin
    v_p1 <- v_muscle + v_bone + v_adipose + v_brain + v_other + v_liver +
            v_stomach + v_gut + v_pancreas + v_spleen + v_heart + v_kidney
    v_p2 <- v_skin

    # ================= Clearance partitioning =================
    # Appendix S1 $PK: when the hepatic clearance implied by total blood
    # clearance meets or exceeds hepatic blood flow, hepatic clearance is
    # held at 0.99 * Q_liver and the balance is carried as renal clearance
    # (also stated in the Figure S4 footnote). For the rat diazepam
    # estimates this branch is active.
    clh <- min(cl * (1 - fexc), 0.99 * q_liver)
    clr <- cl - clh

    # ================= Effective (Kb-weighted) volumes =================
    v_conv <- v_arterial + v_venous + v_lung * kb_central
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
    cblood <- central / v_conv

    # Volume of distribution at steady state on a whole-blood basis
    # (Eq S8 / Eq 9), reported as "Vss median" in Table S6 (1.13 L).
    vssb <- v_arterial + v_venous + kb_central * v_lung + v_p1_eff + v_p2_eff

    Cc <- cblood / bp
    Cc ~ prop(propSd)
  })
}
