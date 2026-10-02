# Population PK model for buprenorphine after sublingual run-in dosing and
# monthly subcutaneous BUP-XR (RBP-6000, SUBLOCADE) depot injections in adults
# with opioid use disorder (Jones 2021, Clin Pharmacokinet 60(4):527-540;
# doi:10.1007/s40262-020-00957-0).

Jones_2021_buprenorphine <- function() {
  description <- paste(
    "Two-compartment population pharmacokinetic model for buprenorphine in",
    "adults with opioid use disorder receiving sublingual (SL) buprenorphine",
    "during a run-in period followed by up to 12 monthly subcutaneous (SC)",
    "injections of BUP-XR (RBP-6000, SUBLOCADE; buprenorphine in the ATRIGEL",
    "delivery system, 50-300 mg). SL buprenorphine is absorbed first-order",
    "(ka) from a single SL depot with a bioavailability fdepot relative to",
    "BUP-XR, modified by the SL formulation (buprenorphine/naloxone film vs",
    "tablet, on both ka and fdepot) and reduced for SL doses of 16 mg or",
    "more. The BUP-XR dose is split by a logit-normal fraction frel into a",
    "fast first-order pathway (ka_fast, the early 'initial burst' peak at",
    "about 24 h) and a slow pathway of one SC depot (ka_slow) and one",
    "transit compartment (ktr) that mimics slow release from the solidified",
    "depot. All clearances and volumes are apparent (relative to BUP-XR) and",
    "allometrically scaled by body weight (exponents 0.75 and 1, reference",
    "70 kg); body mass index is a power covariate on CL/F and ka_fast",
    "(reference 24.8 kg/m^2). Combined proportional and additive residual",
    "error."
  )
  reference <- paste(
    "Jones AK, Ngaimisi E, Gopalakrishnan M, Young MA, Laffont CM (2021).",
    "Population Pharmacokinetics of a Monthly Buprenorphine Depot Injection",
    "for the Treatment of Opioid Use Disorder: A Combined Analysis of Phase",
    "II and Phase III Trials. Clinical Pharmacokinetics 60(4):527-540.",
    "doi:10.1007/s40262-020-00957-0."
  )
  vignette <- "Jones_2021_buprenorphine"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Declared explicitly: buildModelDb()'s auto-detection only recognises
  # compartments literally named `depot` / `central`. `depot` is the SL
  # absorption compartment (Figure 2, compartment 1). A BUP-XR injection is
  # dosed into BOTH `depot_fast1` (Figure 2, compartment 2) and `depot_slow1`
  # (compartment 3) with the same amount; f() in model() splits it between
  # the fast and slow pathways (F2 and F3 = 1 - F2 in Figure 2).
  dosing <- c("depot", "depot_fast1", "depot_slow1")

  # Every ODE state holds an amount of buprenorphine in mg. Only plasma
  # buprenorphine derived from `central` was assayed (Methods 2.2, LC-MS/MS,
  # calibration range 0.050-25.0 ng/mL).
  compartmentData <- list(
    depot = list(analyte = "buprenorphine", units = "mg", specimen = "administration site", verified = TRUE),
    depot_fast1 = list(analyte = "buprenorphine", units = "mg", specimen = "administration site", verified = TRUE),
    depot_slow1 = list(analyte = "buprenorphine", units = "mg", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "buprenorphine", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "buprenorphine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "buprenorphine", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Allometric scaling of all clearances (CL/F, Q/F; exponent 0.75) and",
        "volumes (V4/F, V5/F; exponent 1) with reference weight 70 kg",
        "(Methods 2.3.1; Table 3 footnote). Baseline value used (Methods",
        "2.3.2: 'baseline values were used as no systematic trend was",
        "observed over time')."
      ),
      source_name = "Weight"
    ),
    BMI = list(
      description = "Baseline body mass index",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power covariate on CL/F (exponent -0.362) and on the fast SC",
        "absorption rate constant k24 (exponent -1.32), normalised to the",
        "median BMI of 24.8 kg/m^2 (Table 3 footnote). Baseline value."
      ),
      source_name = "BMI"
    ),
    FORM_BPN_FILM = list(
      description = "Sublingual buprenorphine formulation: buprenorphine/naloxone SL film (1) vs buprenorphine SL tablet (0)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (buprenorphine SL tablet, used in the phase IIa Study 1 run-in)",
      notes = paste(
        "1 = buprenorphine/naloxone SL film (the phase III Studies 2 and 3",
        "run-in), 0 = buprenorphine SL tablet (Study 1). Multiplies the SL",
        "absorption rate constant by FRK14 = 0.636 and the SL relative",
        "bioavailability by FRF1 = 1.47 (Table 3; Methods 2.3.2). Affects",
        "only the SL `depot`; no effect on BUP-XR."
      ),
      source_name = "SL film vs tablet (FRK14, FRF1)"
    ),
    DOSE_BPN_SL_MG = list(
      description = "Administered sublingual buprenorphine dose",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Daily SL buprenorphine dose in mg on the SL dose records. SL",
        "relative bioavailability F1 is multiplied by F1DOSE = 0.765 when",
        "the SL dose is 16 mg or more (Section 3.1; Table 3). Unused for",
        "BUP-XR records (the value only scales the SL `depot`)."
      ),
      source_name = "SL dose (F1DOSE threshold 16 mg)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 570L,
    n_studies = 3L,
    n_observations = 19686L,
    age_range = "19-64 years",
    age_mean = "38.8 years (SD 11.5)",
    weight_range = "46.1-132.0 kg",
    weight_mean = "76.5 kg (SD 15.5)",
    bmi_range = "18.0-35.0 kg/m^2 (median 24.8)",
    sex_female_pct = 32.1,
    race_ethnicity = c(White = 69.5, `Black/African American` = 28.2, Other = 2.3),
    disease_state = "Treatment-seeking adults with opioid use disorder (DSM-IV-TR opioid dependence in Study 1; moderate or severe DSM-5 OUD in Studies 2 and 3)",
    dose_range = paste(
      "SL buprenorphine run-in 8-24 mg/day (tablet in Study 1,",
      "buprenorphine/naloxone film in Studies 2-3), then BUP-XR 50-300 mg SC",
      "every 28 days for up to 12 injections (phase III regimens 300/100 mg",
      "and 300/300 mg)"
    ),
    regions = "USA",
    notes = paste(
      "Phase IIa multiple-ascending-dose Study 1 (NCT01738503, 103",
      "subjects), phase III double-blind efficacy Study 2 (NCT02357901, 434",
      "subjects incl. 16 placebo run-in subjects), and phase III open-label",
      "long-term safety Study 3 (NCT02510014, 287 subjects); Tables 1-2."
    )
  )

  ini({
    # Disposition (apparent, relative to BUP-XR) at 70 kg and BMI 24.8 kg/m^2
    lcl <- log(52.2); label("Apparent clearance CL/F (L/h)") # Table 3 'CL/F' = 52.2 (1.5% RSE)
    lvc <- log(432); label("Apparent central volume V4/F (L)") # Table 3 'V4/F' = 432 (6.1% RSE)
    lq <- fixed(log(79.5)); label("Apparent intercompartmental clearance Q/F (L/h)") # Table 3 'Q/F' = 79.5 (fixed to Study 1 estimate, Section 3.2)
    lvp <- fixed(log(1110)); label("Apparent peripheral volume V5/F (L)") # Table 3 'V5/F' = 1110 (fixed to Study 1 estimate, Section 3.2)

    # SL absorption
    lka <- fixed(log(1.17)); label("SL absorption rate constant k14, tablet (1/h)") # Table 3 'k14' = 1.17 (fixed to Study 1 estimate)
    lfdepot <- fixed(log(0.185)); label("SL tablet bioavailability relative to BUP-XR F1, dose < 16 mg (fraction)") # Table 3 'F1' = 0.185 (fixed to Study 1 estimate)
    e_form_film_ka <- 0.636; label("Multiplier on k14 for SL film vs tablet FRK14 (unitless)") # Table 3 'FRK14' = 0.636 (11% RSE)
    e_form_film_fdepot <- 1.47; label("Multiplier on F1 for SL film vs tablet FRF1 (unitless)") # Table 3 'FRF1' = 1.47 (3.5% RSE)
    e_dose_high_fdepot <- fixed(0.765); label("Multiplier on F1 for SL dose >= 16 mg vs < 16 mg F1DOSE (unitless)") # Table 3 'F1DOSE' = 0.765 (fixed; Table S2 footnote b)

    # BUP-XR (SC) dual absorption
    lka_fast <- log(0.0277); label("Fast absorption rate constant from SC depot k24 (1/h)") # Table 3 'k24' = 0.0277 (5.0% RSE)
    lka_slow <- log(0.00392); label("Slow absorption rate constant from SC depot to transit k36 (1/h)") # Table 3 'k36' = 0.00392 (7.5% RSE)
    lktr <- log(0.000507); label("Rate constant from transit to central compartment k64 (1/h)") # Table 3 'k64' = 0.000507 (3.5% RSE)
    logitfrel <- log(0.0680 / (1 - 0.0680)); label("Fraction of SC dose absorbed by the fast process F2 (logit)") # Table 3 'F2' = 0.0680 (2.1% RSE), logit-normal

    # Covariate effects
    e_wt_cl_q <- fixed(0.75); label("Allometric exponent of body weight on CL/F and Q/F (unitless)") # Methods 2.3.1 (fixed 0.75); Table 3 footnote
    e_wt_vc_vp <- fixed(1); label("Allometric exponent of body weight on V4/F and V5/F (unitless)") # Methods 2.3.1 (fixed 1)
    e_bmi_cl <- -0.362; label("Power exponent of BMI/24.8 on CL/F (unitless)") # Table 3 'theta BMI (CL)' = -0.362 (21% RSE)
    e_bmi_ka_fast <- -1.32; label("Power exponent of BMI/24.8 on k24 (unitless)") # Table 3 'theta BMI (k24)' = -1.32 (14% RSE)

    # IIV (variances; diagonal -- off-diagonal elements not reported)
    etalcl ~ 0.0909 # Table 3 CL/F variance 0.0909 (30.9% CV)
    etalvc ~ 0.704 # Table 3 V4/F variance 0.704 (101% CV)
    etalq ~ fixed(0.334) # Table 3 Q/F variance 0.334 (Study 1 estimate; 62.9% CV)
    etalvp ~ fixed(0.941) # Table 3 V5/F variance 0.941 (Study 1 estimate; 125% CV)
    etalka ~ fixed(0.190) # Table 3 k14 variance 0.190 (Study 1 estimate; 45.7% CV)
    etalka_fast ~ 0.643 # Table 3 k24 variance 0.643 (95.0% CV)
    etalka_slow ~ 1.69 # Table 3 k36 variance 1.69 (210% CV)
    etalktr ~ 0.384 # Table 3 k64 variance 0.384 (68.4% CV)
    etalfdepot ~ fixed(0.195) # Table 3 F1 variance 0.195 (Study 1 estimate; 46.4% CV)
    etalogitfrel ~ 0.194 # Table 3 F2 variance 0.194 (logit-normal)

    # Residual error
    propSd <- 0.190; label("Proportional residual error (fraction)") # Table 3 'PROP' = 0.190 (0.66% RSE)
    addSd <- 0.0378; label("Additive residual error (ng/mL)") # Table 3 'ADD' = 0.0378 ng/mL (13% RSE)
  })
  model({
    # Individual disposition parameters (allometry on WT/70; BMI/24.8 power)
    cl <- exp(lcl + etalcl) * (BMI / 24.8)^e_bmi_cl * (WT / 70)^e_wt_cl_q
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc_vp
    q <- exp(lq + etalq) * (WT / 70)^e_wt_cl_q
    vp <- exp(lvp + etalvp) * (WT / 70)^e_wt_vc_vp

    # SL absorption: film vs tablet on k14 and F1; F1 reduced for doses >= 16 mg
    ka <- exp(lka + etalka) * e_form_film_ka^FORM_BPN_FILM
    dose_high <- 0
    if (DOSE_BPN_SL_MG >= 16) dose_high <- 1
    fdepot <- exp(lfdepot + etalfdepot) *
      e_form_film_fdepot^FORM_BPN_FILM *
      e_dose_high_fdepot^dose_high

    # BUP-XR dual absorption
    ka_fast <- exp(lka_fast + etalka_fast) * (BMI / 24.8)^e_bmi_ka_fast
    ka_slow <- exp(lka_slow + etalka_slow)
    ktr <- exp(lktr + etalktr)
    logitfrel_ind <- logitfrel + etalogitfrel
    frel <- expit(logitfrel_ind)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(depot_fast1) <- -ka_fast * depot_fast1
    d/dt(depot_slow1) <- -ka_slow * depot_slow1
    d/dt(transit1) <- ka_slow * depot_slow1 - ktr * transit1
    d/dt(central) <- ka * depot + ka_fast * depot_fast1 + ktr * transit1 -
      kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    f(depot) <- fdepot
    f(depot_fast1) <- frel
    f(depot_slow1) <- 1 - frel

    # mg / L -> ng/mL
    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
