Groenendaal_2007_morphine_rat_eeg <- function() {
  description <- paste(
    "Preclinical (rat, male Wistar). Population PK-PD model for the effect of",
    "morphine on the amplitude of the delta frequency band (0.5-4.5 Hz) of the",
    "rat EEG, with an extended catenary biophase distribution model and a",
    "P-glycoprotein (Pgp) interaction at the blood-brain barrier (Groenendaal",
    "2007, Br J Pharmacol 151(6):713-720). Blood disposition is a",
    "three-compartment model that serves purely as the input function; its",
    "parameters were not re-estimated here and are fixed from the companion",
    "paper (Groenendaal 2007, Br J Pharmacol 151(5):701-712, Table 2), which",
    "carries linear body-weight covariate effects on CL and V2. Biophase",
    "distribution is the 'tank-in-series' extended catenary model of Upton",
    "2000: two sequential biophase compartments holding CONCENTRATIONS, a",
    "transfer compartment (effect1, paper symbol Cet) fed from blood at rate",
    "k1e and drained at the same rate k1e, followed by the effect compartment",
    "(effect2, paper symbol Ce) drained at rate keo. The asymmetric form",
    "(k1e != keo) was retained over the symmetric form on objective function",
    "(24671 vs 24936); the simpler one-compartment effect-site model could be",
    "fit in neither its symmetric nor its asymmetric form. Co-infusion of the",
    "Pgp inhibitor GF120918 (elacridar) reduces keo by 64 percent and leaves",
    "k1e unchanged. Effect-compartment concentrations drive EEG amplitude",
    "through a sigmoidal Emax model. Inter-animal variability is exponential",
    "on keo and PROPORTIONAL on E0 and Emax; all other variances were fixed to",
    "zero by the authors. A covariate contrasting the EEG and EEG-microdialysis",
    "experimental methods was screened on E0 and Emax and not retained.",
    sep = " "
  )
  reference <- paste(
    "Groenendaal D, Freijer J, de Mik D, Bouw MR, Danhof M, de Lange ECM.",
    "Influence of biophase distribution and P-glycoprotein interaction on",
    "pharmacokinetic-pharmacodynamic modelling of the effects of morphine on",
    "the EEG.",
    "Br J Pharmacol. 2007;151(6):713-720.",
    "doi:10.1038/sj.bjp.0707258.",
    "Blood pharmacokinetic parameters fixed from the companion paper:",
    "Groenendaal D, Freijer J, de Mik D, Bouw MR, Danhof M, de Lange ECM.",
    "Population pharmacokinetic modelling of non-linear brain distribution of",
    "morphine: influence of active saturable influx and P-glycoprotein",
    "mediated efflux.",
    "Br J Pharmacol. 2007;151(5):701-712;",
    "doi:10.1038/sj.bjp.0707257, Table 2.",
    sep = " "
  )
  vignette <- "Groenendaal_2007_morphine_rat_eeg"

  units <- list(
    time          = "min",
    dosing        = "ng",
    concentration = "ng/mL"
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Per-animal, time-fixed. Enters the blood-PK layer only, as the linear",
        "centred form of the companion paper's Equation 7,",
        "P_i = theta1 * (1 + theta2 * (BW_i - median BW)), on CL and on V2",
        "(peripheral1 volume). The companion paper never prints the value of",
        "'median BW'; 0.300 kg is used here, the rounded N-weighted mean of the",
        "seven per-group means in companion Table 1 (0.260-0.306 kg, N-weighted",
        "mean 0.294 kg). See the vignette Assumptions section.",
        sep = " "
      ),
      source_name        = "BW"
    ),
    CONMED_ELACRIDAR = list(
      description        = "Continuous co-infusion of the P-glycoprotein inhibitor GF120918 (elacridar)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = vehicle co-infusion",
      notes              = paste(
        "Experimentally administered Pgp blockade, not a clinical",
        "co-medication. Regimen (companion paper, Experimental procedures):",
        "6 mg/kg 1-min intravenous bolus in dimethyl sulphoxide followed by a",
        "continuous 25 ng/min infusion, started 120 min before morphine and",
        "maintained for the whole experiment; mean steady-state blood",
        "concentration 214 ng/mL. Set to 1 for the GF120918 arm and 0 for the",
        "vehicle arm; per-animal and time-fixed within an experiment. Acts on",
        "keo only (lead paper Table 1); k1e was unaffected.",
        sep = " "
      ),
      source_name        = "GF120918"
    )
  )

  # Screened by the authors but NOT retained in the final model, so it is
  # documented rather than encoded (an unreferenced covariateData entry would
  # raise a "declared but not referenced" convention warning).
  covariatesDataExcluded <- list(
    STUDY_EEG_MD = list(
      description        = "Experimental sub-study: EEG-microdialysis animals versus EEG-only animals",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = EEG-only experiment",
      notes              = paste(
        "Lead paper Equation 6, P_i = theta1 * (1 - METHOD_i) + theta2 * METHOD_i,",
        "with METHOD = 1 for EEG-MD and 0 for EEG. Screened on E0 and Emax to",
        "check whether removing three EEG electrodes and implanting a",
        "microdialysis probe shifted the baseline EEG. Results, 'PK-PD analysis",
        "of the EEG effect': 'Since no differences were observed in E0 and Emax",
        "values between the EEG and the EEG-MD group, a single parameter value",
        "was estimated.' Not retained; no point estimate is reported.",
        sep = " "
      ),
      source_name        = "METHOD"
    )
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. The two biophase states are hypothetical compartments of
  # negligible volume that hold CONCENTRATIONS (lead paper Equations 3 and 4
  # are written in concentration units), so their units are ng/mL rather than
  # an amount and their specimen is "not applicable".
  compartmentData <- list(
    central     = list(analyte = "morphine", units = "ng",    specimen = "whole blood",    verified = TRUE),
    peripheral1 = list(analyte = "morphine", units = "ng",    specimen = "whole blood",    verified = TRUE),
    peripheral2 = list(analyte = "morphine", units = "ng",    specimen = "whole blood",    verified = TRUE),
    effect1     = list(analyte = "morphine", units = "ng/mL", specimen = "not applicable", verified = TRUE),
    effect2     = list(analyte = "morphine", units = "ng/mL", specimen = "not applicable", verified = TRUE)
  )

  population <- list(
    species       = "rat (male Wistar)",
    n_subjects    = 68,
    n_studies     = 1,
    weight_range  = "0.25-0.35 kg",
    weight_median = "0.30 kg",
    sex_female_pct = 0,
    disease_state = "healthy, chronically instrumented conscious rats",
    dose_range    = "4, 10 or 40 mg/kg morphine hydrochloride as a 10-min intravenous infusion",
    regions       = "Leiden, The Netherlands (single laboratory)",
    co_medication = paste(
      "All animals received a continuous midazolam infusion (5.5 mg/kg/h,",
      "Wagner loading scheme) to prevent opioid-induced seizure activity, giving",
      "a mean steady-state blood concentration of 937 ng/mL. Animals in the Pgp",
      "inhibition arm additionally received GF120918 (elacridar) 6 mg/kg bolus",
      "plus 25 ng/min continuous infusion. Animals receiving 40 mg/kg morphine",
      "were artificially ventilated and given vecuronium bromide (0.15 mg bolus,",
      "0.10 mg as needed) for muscle rigidity.",
      sep = " "
    ),
    notes = paste(
      "Groups are the EEG (experiment 1) and EEG-microdialysis (experiment 2)",
      "arms of companion Table 1: EEG 4 mg/kg vehicle N = 7 (0.294 kg),",
      "EEG 10 mg/kg vehicle N = 7 (0.260 kg), EEG 40 mg/kg vehicle N = 5",
      "(0.273 kg), EEG/MD 4 mg/kg vehicle N = 14 (0.297 kg), EEG/MD 4 mg/kg +",
      "GF120918 N = 20 (0.300 kg), EEG/MD 40 mg/kg vehicle N = 15 (0.306 kg).",
      "The microdialysis-only arm (N = 3) of the companion paper contributed no",
      "EEG data and is excluded. EEG endpoint is the delta-band (0.5-4.5 Hz)",
      "amplitude in uV, averaged over 3-10 min bins out to 360 min.",
      sep = " "
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Blood pharmacokinetics: three-compartment input function.
    #
    # The lead paper does not re-fit the blood PK. Figure 1 legend: "The blood
    # PKs were described with a three-compartment model and used as input
    # function for morphine in the brain"; Data analysis: "Individual PK
    # parameter estimates were used as input for the PD models." Every value in
    # this block is therefore fixed() and traced to Table 2 of the companion
    # paper (doi:10.1038/sj.bjp.0707257), which fit ADVAN11 TRANS4 with
    # parameters CL, Q2, Q3, V1, V2, V3. NONMEM V1/V2/V3 map to the canonical
    # vc/vp/vp2 and Q2/Q3 to q/q2.
    #
    # Note on companion Table 2: the LLCI-ULCI column is transposed between the
    # two "Slope factor" rows. Recomputing each row's 95% CI as estimate +/-
    # 1.96 * CV% * estimate reproduces the printed interval exactly for V1, Q2,
    # V2-intercept, Q3 and V3; for the slope factors it matches only after
    # swapping the two intervals (CL slope 5.35 at CV 25.2% gives 2.70-7.99, the
    # interval printed on the V2 row, and V2 slope 8.50 at CV 17.1% gives
    # 5.66-11.3, the interval printed on the CL row). Only the Estimate and CV%
    # columns are used here, and those are internally consistent as printed.
    lcl <- fixed(log(20.0))
    label("Morphine blood clearance at median body weight (CL, mL/min)")                    # companion Table 2, Cl Intercept = 20.0 (CV 5.6%)
    e_wt_cl <- fixed(5.35)
    label("Linear body-weight slope on CL, centred on median body weight (per kg)")         # companion Table 2, Cl Slope factor = 5.35 (CV 25.2%)
    lvc <- fixed(log(68.1))
    label("Central volume of distribution (V1, mL)")                                        # companion Table 2, V1 = 68.1 (CV 16.7%)
    lq <- fixed(log(15.5))
    label("Intercompartmental clearance to peripheral1 (Q2, mL/min)")                       # companion Table 2, Q2 = 15.5 (CV 11.3%)
    lvp <- fixed(log(739))
    label("First peripheral volume at median body weight (V2, mL)")                         # companion Table 2, V2 Intercept = 739 (CV 7.6%)
    e_wt_vp <- fixed(8.50)
    label("Linear body-weight slope on V2, centred on median body weight (per kg)")         # companion Table 2, V2 Slope factor = 8.50 (CV 17.1%)
    lq2 <- fixed(log(17.8))
    label("Intercompartmental clearance to peripheral2 (Q3, mL/min)")                       # companion Table 2, Q3 = 17.8 (CV 18.4%)
    lvp2 <- fixed(log(133))
    label("Second peripheral volume (V3, mL)")                                              # companion Table 2, V3 = 133 (CV 15.9%)

    # ------------------------------------------------------------------
    # Biophase distribution: extended catenary model, estimated in this paper.
    # Lead paper Equations 3 and 4:
    #   dCet/dt = k1e * Cb  - k1e * Cet
    #   dCe/dt  = k1e * Cet - keo * Ce
    lk1e <- log(0.0378)
    label("Rate constant for transport through the transfer compartment (k1e, 1/min)")      # lead Table 1, k1e = 0.0378 (CV 8.4%; 95% CI 0.0315-0.0441)
    lke0 <- log(0.0426)
    label("Rate constant for loss from the effect compartment, vehicle arm (keo, 1/min)")   # lead Table 1, keo -GF120918 = 0.0426 (CV 10.0%; 95% CI 0.0342-0.0510)

    # Lead paper Equation 7, P_i = theta3 * (1 + theta4 * GF120918_i). Table 1
    # gives theta4 = -0.644 and reports keo(+GF120918) = 0.0152; the identity
    # 0.0426 * (1 - 0.644) = 0.01517 reproduces it. The Results text describes
    # this as a 64 percent decrease (the abstract rounds it to 60 percent).
    e_conmed_elacridar_ke0 <- -0.644
    label("Fractional change in keo under GF120918 (elacridar) co-infusion (unitless)")     # lead Table 1, Pgp inhibition factor = -0.644 (CV 7.3%; 95% CI -0.736 to -0.552)

    # ------------------------------------------------------------------
    # Sigmoidal Emax pharmacodynamics, lead paper Equation 5:
    #   E = E0 + Emax * Ce^nH / (EC50^nH + Ce^nH)
    # (Figure 1 renders the denominator of panel (a) as "EC50 - Ce^nH"; that is
    # a PDF sign artefact, Equation 5 in the text is the authoritative form.)
    le0 <- log(44.6)
    label("Baseline no-drug EEG delta-band amplitude (E0, uV)")                             # lead Table 2, E0 = 44.6 (CV 2.3%; 95% CI 42.6-46.6)
    lemax <- log(44.5)
    label("Maximum drug-induced increase in EEG delta-band amplitude (Emax, uV)")           # lead Table 2, Emax = 44.5 (CV 8.0%; 95% CI 37.5-51.5)
    lec50 <- log(451)
    label("Effect-compartment concentration giving half-maximal effect (EC50, ng/mL)")      # lead Table 2, EC50 = 451 (CV 17.3%; 95% CI 298-604)
    lhill <- log(2.32)
    label("Hill slope factor of the sigmoidal Emax model (nH, unitless)")                   # lead Table 2, nH = 2.32 (CV 10.4%; 95% CI 1.85-2.79)

    # ------------------------------------------------------------------
    # Inter-animal variability. The published tables label these rows "omega^2",
    # so each value is a VARIANCE and is used unchanged.
    #
    # Blood-PK variances are inherited with the rest of the blood-PK block and
    # are therefore fixed().
    etalcl ~ fixed(0.129)                                                                   # companion Table 2, omega^2 Cl = 0.129 (CV 17.2%)
    etalvp ~ fixed(0.099)                                                                   # companion Table 2, omega^2 V2 = 0.099 (CV 24.7%)

    # keo carries EXPONENTIAL variability, lead paper Equation 10,
    # P_i = P_typ * exp(eta_i).
    etalke0 ~ 0.237                                                                         # lead Table 1, omega^2 keo = 0.237 (CV 20.2%)

    # E0 and Emax carry PROPORTIONAL variability, lead paper Equation 8,
    # P_i = P_typ * (1 + eta_i). The model() block applies them in that form;
    # the eta names keep the canonical eta + transformed-parameter spelling.
    etale0 ~ 0.034                                                                          # lead Table 2, omega^2 E0 = 0.034 (CV 17.8%)
    etalemax ~ 0.121                                                                        # lead Table 2, omega^2 Emax = 0.121 (CV 24.1%)

    # ------------------------------------------------------------------
    # Residual error. Both papers report a proportional error model,
    # C_obs = C_pred * (1 + eps) with eps ~ N(0, sigma^2), and tabulate
    # sigma^2; propSd is the standard deviation, so sqrt() is taken here.
    propSd <- fixed(sqrt(0.074))
    label("Proportional residual error on blood morphine concentration (fraction)")         # companion Table 2, proportional error sigma^2 = 0.074 (CV 10.2%)
    propSd_eeg <- sqrt(0.027)
    label("Proportional residual error on EEG delta-band amplitude (fraction)")             # lead Table 2, proportional error sigma^2 = 0.027 (CV 7.6%)
  })

  model({
    # ---- 1. Reference values -----------------------------------------
    # Centring weight for the companion paper's Equation 7. The paper prints
    # the equation but never the value of "median BW"; see covariateData$WT.
    wt_median <- 0.300 # kg

    # ---- 2. Individual parameters ------------------------------------
    cl <- exp(lcl + etalcl) * (1 + e_wt_cl * (WT - wt_median))
    vc <- exp(lvc)
    q <- exp(lq)
    vp <- exp(lvp + etalvp) * (1 + e_wt_vp * (WT - wt_median))
    q2 <- exp(lq2)
    vp2 <- exp(lvp2)

    k1e <- exp(lk1e)
    ke0 <- exp(lke0 + etalke0) * (1 + e_conmed_elacridar_ke0 * CONMED_ELACRIDAR)

    # Proportional inter-animal variability (Equation 8), not exponential.
    e0 <- exp(le0) * (1 + etale0)
    emax <- exp(lemax) * (1 + etalemax)
    ec50 <- exp(lec50)
    hill <- exp(lhill)

    # ---- 3. Micro-constants ------------------------------------------
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # ---- 4. Blood disposition (three-compartment input function) ------
    d/dt(central) <- -(kel + k12 + k13) * central + k21 * peripheral1 + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2

    Cc <- central / vc

    # ---- 5. Extended catenary biophase (Equations 3 and 4) ------------
    # Both states hold concentrations (ng/mL), so blood concentration Cc -- not
    # the central amount -- is the driving input. The transfer compartment is
    # symmetric in k1e by construction; only the terminal effect compartment
    # carries the distinct loss constant keo.
    d/dt(effect1) <- k1e * Cc - k1e * effect1
    d/dt(effect2) <- k1e * effect1 - ke0 * effect2

    # ---- 6. Observations ---------------------------------------------
    # ce is clamped at zero so that the fractional Hill power cannot be applied
    # to a slightly negative solver value (both states start at exactly zero).
    ce <- max(effect2, 0.0)
    eeg <- e0 + emax * ce^hill / (ec50^hill + ce^hill)

    Cc ~ prop(propSd)
    eeg ~ prop(propSd_eeg)
  })
}
