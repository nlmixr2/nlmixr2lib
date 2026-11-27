Gomeni_2019_methylphenidate_ivivc <- function() {
  description <- paste(
    "Convolution-based level A in vitro/in vivo correlation (IVIVC) for",
    "methylphenidate from three Ritalin LA (bimodal extended-release capsule)",
    "formulations with slow, medium and fast in vitro dissolution, in healthy",
    "adults (Gomeni 2019). The fraction of the dose still unreleased in vitro",
    "is a double Weibull function of time, fitted per formulation and held",
    "fixed. The in vivo release is that in vitro function evaluated at a",
    "time-scaled clock tt = b1 + b2 * t^b3 (Levy-type time scaling with a",
    "linear intercept a1 and slope a2), common to all formulations. The in",
    "vivo input rate is the central-difference derivative of the scaled",
    "release (step 0.01 h, as in the authors' NONMEM code), multiplied by the",
    "relative bioavailability F and the dose, and feeds a one-compartment",
    "disposition whose elimination rate and volume are fixed from the Ritalin",
    "IR fit (Gomeni_2019_methylphenidate_ir). The formulation is selected by",
    "the FORM_METHYLPHENIDATE_MEDIUM and FORM_METHYLPHENIDATE_FAST indicators",
    "(slow is the reference). Fitted to mean (digitised) data, so the model",
    "has no between-subject variability; additive residual error. The input",
    "is computed from the most recent dose only, so the model is exact for a",
    "single dose and for once-daily dosing once the previous dose has been",
    "fully released."
  )
  reference <- paste(
    "Gomeni R, Fang L, Bressolle-Gomeni F, Spencer TJ, Faraone SV, Babiskin A.",
    "A General Framework for Assessing In vitro/In vivo Correlation as a Tool",
    "for Maximizing the Benefit-Risk Ratio of a Treatment Using a",
    "Convolution-Based Modeling Approach. CPT Pharmacometrics Syst Pharmacol.",
    "2019;8(2):97-106. doi:10.1002/psp4.12378. Dissolution parameters from",
    "Table S2 and the NONMEM control stream (Supplementary model code); IVIVC",
    "estimates from Table S3 Scenario 2 (common F, the final model);",
    "disposition fixed from Table S1 and the control stream."
  )
  vignette <- "Gomeni_2019_methylphenidate"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # The dose is recorded in depot, which tracks the amount not yet released
  # in vivo. The input rate itself is computed from podo(depot) and the time
  # since that dose, as in the authors' NONMEM $DES block (DOSE * KAB).
  dosing <- "depot"

  compartmentData <- list(
    depot = list(analyte = "methylphenidate", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "methylphenidate", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    FORM_METHYLPHENIDATE_MEDIUM = list(
      description = "1 = the dose is the medium-release Ritalin LA formulation; 0 = slow- or fast-release formulation.",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "Selects the medium formulation's fixed in vitro dissolution set",
        "(Table S2; control stream ID 2). Slow release is the reference",
        "(both indicators 0). Setting both indicators to 1 is not a valid",
        "formulation."
      ),
      source_name = "ID (control stream ID.EQ.2)"
    ),
    FORM_METHYLPHENIDATE_FAST = list(
      description = "1 = the dose is the fast-release Ritalin LA formulation; 0 = slow- or medium-release formulation.",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "Selects the fast formulation's fixed in vitro dissolution set",
        "(Table S2; control stream ID 3). Slow release is the reference",
        "(both indicators 0)."
      ),
      source_name = "ID (control stream ID.EQ.3)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 16L,
    n_studies = 1L,
    studies = paste(
      "Wang Y et al. Biopharm Drug Dispos 2004;25:91-98: four-period,",
      "four-treatment, single-dose randomised crossover in healthy volunteers",
      "after at least a 10-hour fast; three Ritalin LA 40 mg capsule",
      "formulations (slow, medium and fast in vitro release) and Ritalin IR",
      "(two 20 mg tablets 4 h apart). 17 enrolled, 16 in the PK analysis."
    ),
    disease_state = "Healthy volunteers",
    dose_range = "40 mg single oral dose of each Ritalin LA formulation; 10-40 mg simulated",
    regions = "Not reported",
    notes = paste(
      "The analysis used the published MEAN concentration-time profiles,",
      "digitised from Wang 2004 with Plot Digitizer 2.0, and the mean in",
      "vitro dissolution profiles. Age, weight and sex are not reported in",
      "Gomeni 2019. The SKAMP exposure-response layer (Eqs. 13-15) that the",
      "paper drives with this model takes its parameter values from Gomeni",
      "2017 (Clin Pharmacol Ther 102:951, doi:10.1002/cpt.684), which are not",
      "printed in Gomeni 2019, and is not part of this model."
    )
  )

  ini({
    # In vitro double Weibull, Eq. 9: unreleased fraction
    # ff * exp(-(t/td)^ss) + (1 - ff) * exp(-(t/td1)^ss1). The canonical rate
    # scalers are reciprocals of the printed times, ra = 1/td and ra2 = 1/td1.
    # Table S2 prints td to two decimals; the control stream carries the
    # values the IVIVC fit used, which are taken here.
    lra_slow <- fixed(log(1 / 0.186)); label("Slow formulation: first-process in vitro Weibull rate scaler 1/td (1/h)") # Control stream ID 1: TD = 0.186 h (Table S2 slow: td = 0.19 h)
    lgam1_slow <- fixed(log(1.7)); label("Slow formulation: first-process in vitro Weibull shape ss (unitless)") # Control stream ID 1: SS = 1.7 (Table S2 slow: ss = 1.70)
    lra2_slow <- fixed(log(1 / 9.91)); label("Slow formulation: second-process in vitro Weibull rate scaler 1/td1 (1/h)") # Control stream ID 1: TD1 = 9.91 h (Table S2 slow: td1 = 9.91 h)
    lgam2_slow <- fixed(log(3.59)); label("Slow formulation: second-process in vitro Weibull shape ss1 (unitless)") # Control stream ID 1: SS1 = 3.59 (Table S2 slow: ss1 = 3.59)
    logitfrel_slow <- fixed(log(0.5 / (1 - 0.5))); label("Slow formulation: fraction of the dose in the first release process ff (logit of fraction)") # Control stream ID 1: FF = 0.5 (Table S2 slow: ff = 0.50)
    lra_medium <- fixed(log(1 / 0.159)); label("Medium formulation: first-process in vitro Weibull rate scaler 1/td (1/h)") # Control stream ID 2: TD = 0.159 h (Table S2 medium: td = 0.16 h)
    lgam1_medium <- fixed(log(6.55)); label("Medium formulation: first-process in vitro Weibull shape ss (unitless)") # Control stream ID 2: SS = 6.55 (Table S2 medium: ss = 6.55)
    lra2_medium <- fixed(log(1 / 7.55)); label("Medium formulation: second-process in vitro Weibull rate scaler 1/td1 (1/h)") # Control stream ID 2: TD1 = 7.55 h (Table S2 medium: td1 = 7.55 h)
    lgam2_medium <- fixed(log(2.63)); label("Medium formulation: second-process in vitro Weibull shape ss1 (unitless)") # Control stream ID 2: SS1 = 2.63 (Table S2 medium: ss1 = 2.63)
    logitfrel_medium <- fixed(log(0.502 / (1 - 0.502))); label("Medium formulation: fraction of the dose in the first release process ff (logit of fraction)") # Control stream ID 2: FF = 0.502 (Table S2 medium: ff = 0.50)
    lra_fast <- fixed(log(1 / 0.17)); label("Fast formulation: first-process in vitro Weibull rate scaler 1/td (1/h)") # Control stream ID 3: TD = 0.17 h (Table S2 fast: td = 0.17 h)
    lgam1_fast <- fixed(log(1.79)); label("Fast formulation: first-process in vitro Weibull shape ss (unitless)") # Control stream ID 3: SS = 1.79 (Table S2 fast: ss = 1.79)
    lra2_fast <- fixed(log(1 / 5.42)); label("Fast formulation: second-process in vitro Weibull rate scaler 1/td1 (1/h)") # Control stream ID 3: TD1 = 5.42 h (Table S2 fast: td1 = 5.42 h)
    lgam2_fast <- fixed(log(1.75)); label("Fast formulation: second-process in vitro Weibull shape ss1 (unitless)") # Control stream ID 3: SS1 = 1.75 (Table S2 fast: ss1 = 1.75)
    logitfrel_fast <- fixed(log(0.471 / (1 - 0.471))); label("Fast formulation: fraction of the dose in the first release process ff (logit of fraction)") # Control stream ID 3: FF = 0.471 (Table S2 fast: ff = 0.47)

    # IVIVC time scaling, Eq. 3 / Eq. 10: rvivo(t) = a1 + a2 * rvitro(tt),
    # tt = b1 + b2 * t^b3. a1, a2 and b1 are fixed at the identity.
    ivivc_int <- fixed(0); label("IVIVC linear intercept a1 (fraction)") # Table S3 Scenario 2: A1 = 0 FIX
    ivivc_slope <- fixed(1); label("IVIVC linear slope a2 (unitless)") # Table S3 Scenario 2: A2 = 1 FIX
    ivivc_tshift <- fixed(0); label("IVIVC time shift b1 (h of in vitro time)") # Table S3 Scenario 2: B1 = 0 FIX
    livivc_tscale <- log(0.19); label("IVIVC time-scaling factor b2 (in vitro h per in vivo h^b3)") # Table S3 Scenario 2: B2 = 0.19 (SE 0.002, RSE 1.2%)
    livivc_tshape <- log(2.24); label("IVIVC time-shaping exponent b3 (unitless)") # Table S3 Scenario 2: B3 = 2.24 (SE 0.053, RSE 2.4%)
    lfdepot <- log(0.94); label("Bioavailability of Ritalin LA relative to Ritalin IR, common to the three formulations (fraction)") # Table S3 Scenario 2: F = 0.94 (SE 0.018, RSE 1.9%)

    # Disposition fixed from the Ritalin IR fit (step 1 of the analysis).
    lkel <- fixed(log(0.254)); label("Elimination rate constant kel (1/h)") # Control stream: KEL = 0.254 (Table S1: kel = 0.25 1/h)
    lvc <- fixed(log(1200)); label("Volume of distribution V (L)") # Control stream: S1 = 1.2 with dose in mg and concentration in ng/mL, i.e. 1.2 kL = 1200 L (Table S1: V = 1.20)

    addSd <- 1.13; label("Additive residual error (ng/mL)") # Table S3 Scenario 2: err_add = 1.13 (SE 0.099, RSE 8.8%)
  })

  model({
    # Formulation selection: slow is the reference.
    form_slow <- 1 - FORM_METHYLPHENIDATE_MEDIUM - FORM_METHYLPHENIDATE_FAST
    ra <- exp(lra_slow) * form_slow + exp(lra_medium) * FORM_METHYLPHENIDATE_MEDIUM + exp(lra_fast) * FORM_METHYLPHENIDATE_FAST
    gam1 <- exp(lgam1_slow) * form_slow + exp(lgam1_medium) * FORM_METHYLPHENIDATE_MEDIUM + exp(lgam1_fast) * FORM_METHYLPHENIDATE_FAST
    ra2 <- exp(lra2_slow) * form_slow + exp(lra2_medium) * FORM_METHYLPHENIDATE_MEDIUM + exp(lra2_fast) * FORM_METHYLPHENIDATE_FAST
    gam2 <- exp(lgam2_slow) * form_slow + exp(lgam2_medium) * FORM_METHYLPHENIDATE_MEDIUM + exp(lgam2_fast) * FORM_METHYLPHENIDATE_FAST
    frel <- expit(logitfrel_slow) * form_slow + expit(logitfrel_medium) * FORM_METHYLPHENIDATE_MEDIUM + expit(logitfrel_fast) * FORM_METHYLPHENIDATE_FAST

    ivivc_tscale <- exp(livivc_tscale)
    ivivc_tshape <- exp(livivc_tshape)
    fdepot <- exp(lfdepot)
    kel <- exp(lkel)
    vc <- exp(lvc)

    # Central-difference derivative of the scaled release (Eq. 11), with the
    # control stream's step DELT = 0.01 h and its floors at zero. The code
    # applies a1 and a2 to the UNRELEASED fraction and clamps the result to
    # [0, 1]; at the identity values this is the same as Eq. 3.
    tdose <- tad0(depot)
    t_lo <- max(tdose - 0.01, 0)
    t_hi <- tdose + 0.01
    tt_lo <- max(ivivc_tshift + ivivc_tscale * t_lo^ivivc_tshape, 0)
    tt_hi <- max(ivivc_tshift + ivivc_tscale * t_hi^ivivc_tshape, 0)
    unrel_lo <- min(max(ivivc_int + ivivc_slope * (frel * exp(-(ra * tt_lo)^gam1) + (1 - frel) * exp(-(ra2 * tt_lo)^gam2)), 0), 1)
    unrel_hi <- min(max(ivivc_int + ivivc_slope * (frel * exp(-(ra * tt_hi)^gam1) + (1 - frel) * exp(-(ra2 * tt_hi)^gam2)), 0), 1)

    dose_last <- podo(depot)
    if (is.na(dose_last)) {
      dose_last <- 0
    }
    rel_rate <- dose_last * (unrel_lo - unrel_hi) / (t_hi - t_lo)

    d/dt(depot) <- -rel_rate
    d/dt(central) <- fdepot * rel_rate - kel * central

    # Dose in mg and V in L give mg/L; x 1000 for ng/mL.
    Cc <- 1000 * central / vc
    Cc ~ add(addSd)
  })
}
