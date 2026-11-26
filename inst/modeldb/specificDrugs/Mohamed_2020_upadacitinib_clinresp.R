Mohamed_2020_upadacitinib_clinresp <- function() {
  description <- paste(
    "Continuous-time Markov exposure-response model for clinical response",
    "(>= 30% reduction from baseline in very soft/liquid stool frequency",
    "and/or abdominal pain score, neither worse than baseline) in adults with",
    "moderately-to-severely active Crohn's disease during the 16-week CELEST",
    "induction period on upadacitinib (Mohamed 2020). Three states (no",
    "response, response, absorbing dropout) with first-order transitions; the",
    "time-varying upadacitinib plasma concentration multiplies the",
    "no-response -> response rate by 1 + Emax * C / (EC50 + C), and that rate",
    "is further multiplied by 1 + 3.60 during the first two weeks. The",
    "concentration driver is the Crohn's disease arm of the paper's population",
    "PK model (Mohamed_2020_upadacitinib), including the extended-release",
    "absorption branch the paper used to simulate once-daily ER regimens."
  )
  reference <- paste(
    "Mohamed M-EF, Klunder B, Lacerda AP, Othman AA. Exposure-Response",
    "Analyses for Upadacitinib Efficacy and Safety in the Crohn's Disease",
    "CELEST Study and Bridging to the Extended-Release Formulation. Clin",
    "Pharmacol Ther. 2020;107(3):639-649. doi:10.1002/cpt.1668"
  )
  vignette <- "Mohamed_2020_upadacitinib"

  # Markov state-occupancy probabilities (s_<state> namespace), not drug
  # compartments.
  paper_specific_compartments <- c("s_noresp", "s_resp", "s_dropout")
  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "ng/mL (Cc); prob_clinresp is a probability (0-1)"
  )

  compartmentData <- list(
    depot = list(analyte = "upadacitinib", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "upadacitinib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "upadacitinib", units = "mg", specimen = "plasma", verified = TRUE),
    s_noresp = list(
      analyte = "probability of being a clinical non-responder and on study",
      units = NA_character_,
      specimen = "not applicable",
      verified = TRUE
    ),
    s_resp = list(
      analyte = "probability of being a clinical responder and on study",
      units = NA_character_,
      specimen = "not applicable",
      verified = TRUE
    ),
    s_dropout = list(
      analyte = "probability of having dropped out",
      units = NA_character_,
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description = "Baseline total body weight.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "PK driver only: power effect on Vc/F, exponent 0.490, reference 74 kg (Table S3).",
      source_name = "WTKG"
    ),
    SEXF = list(
      description = "Biological sex indicator (1 = female, 0 = male).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male), the PK reference.",
      notes = "PK driver only: female CL/F ratio 0.850 and Vc/F ratio 0.748 (Table S3). No covariates were tested on the exposure-response relationship (Discussion).",
      source_name = "SEX"
    ),
    CRCL = list(
      description = "Baseline creatinine clearance by the Cockcroft-Gault equation (raw mL/min).",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = "PK driver only: power effect on CL/F, exponent 0.313, reference 106.9 mL/min (Table S7 control stream).",
      source_name = "CRCL"
    ),
    FORM_UPA_ER = list(
      description = "Upadacitinib extended-release formulation indicator (1 = extended-release tablet, 0 = immediate-release reference).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (immediate release, the formulation used in CELEST).",
      notes = "PK driver only. Switches to the Table S4 ER absorption branch for simulating ER regimens (Methods: the ER absorption parameters were 'incorporated into the exposure-response models'). Each ER administration needs two dose records: one to 'depot' and one to 'central' with rate = -2.",
      source_name = "FORM (1 = IR, 2 = ER)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 216L,
    n_studies = 1L,
    age_range = "mean 40.8 years (SD 12.9)",
    weight_range = "mean 75.7 kg (SD 20.1)",
    sex_female_pct = 57,
    race_ethnicity = c(White = 86, Black = 6, Hispanic = 7, Asian = 1, Other = 0.5),
    disease_state = "Moderately-to-severely active Crohn's disease (CDAI 220-450, SES-CD >= 6 or >= 4 for ileal-only disease) with inadequate response or intolerance to an immunomodulator or anti-TNF therapy; 96% anti-TNF experienced. Baseline CDAI 300 (SD 59.9), SES-CD 14.9 (SD 8.4).",
    dose_range = "Placebo or immediate-release upadacitinib 3, 6, 12 or 24 mg twice daily or 24 mg once daily for 16 weeks (induction).",
    regions = "Multinational",
    notes = "CELEST (NCT02365649) induction period, Table 1. Symptoms assessed at baseline and weeks 2, 4, 8, 12 and 16. Markov models fitted in NONMEM (Laplace) with individual empirical-Bayes PK parameters as input."
  )

  ini({
    # ---- PK driver: Mohamed_2020_upadacitinib, Crohn's disease arm ---------
    # The Markov fit took each patient's PK parameters as given, so the PK
    # values are fixed here. Male, CrCL 106.9 mL/min, 74 kg reference.
    lcl <- fixed(log(39.6))
    label("Apparent clearance CL/F, CD male at CrCL 106.9 mL/min (L/h)") # Table S3 'CL/F (L/h)' = 39.6; no CD effect on typical CL/F
    lvc <- fixed(log(141))
    label("Apparent central volume Vc/F, male at 74 kg (L)") # Table S3 'Vc/F (L)' = 141
    lvp <- fixed(log(66.3))
    label("Apparent peripheral volume Vp/F (L)") # Table S3 'Vp/F (L)' = 66.3
    lq <- fixed(log(3.39))
    label("Apparent intercompartmental clearance Q/F (L/h)") # Table S3 'Q/F (L/h)' = 3.39
    lka <- fixed(log(3.60 * 1.63))
    label("Absorption rate constant Ka, IR formulation, CD patient (1/h)") # Table S3 'KA (1/h)' = 3.60 x 'Ratio of KA in RA or CD Patients Compared to healthy subjects' = 1.63
    ltlag <- fixed(log(0.469))
    label("Absorption lag time, IR formulation (h)") # Table S3 'ALAG1 (h)' = 0.4690
    e_sexf_cl <- fixed(log(0.850))
    label("Log ratio of CL/F, female vs male (unitless)") # Table S3 = 0.850
    e_crcl_cl <- fixed(0.313)
    label("Power exponent of creatinine clearance on CL/F (unitless)") # Table S3 = 0.313
    e_sexf_vc <- fixed(log(0.748))
    label("Log ratio of Vc/F, female vs male (unitless)") # Table S3 = 0.748
    e_wt_vc <- fixed(0.490)
    label("Power exponent of body weight on Vc/F (unitless)") # Table S3 = 0.490
    lka_er <- fixed(log(0.100))
    label("First-order absorption rate constant Ka, ER formulation (1/h)") # Table S4 'KA ER (1/h)' = 0.100
    ltlag_er <- fixed(log(0.225))
    label("Absorption lag time of both ER input routes (h)") # Table S4 'ALAG1 ER (h)' = 0.225
    logitffo <- fixed(log(0.446 / (1 - 0.446)))
    label("Logit of the fraction of the ER dose absorbed by the first-order route (unitless)") # Table S4 first-order fraction (%) = 44.6
    ld2 <- fixed(log(6.11))
    label("Duration of the zero-order ER input into the central compartment (h)") # Table S4 'Zero-Order Infusion duration D2 (h)' = 6.11
    lfrel_er <- fixed(log(0.742))
    label("Bioavailability of the ER formulation relative to IR (fraction)") # Table S4 relative bioavailability (%) = 74.2

    # ---- Markov transition rates (Table S5, clinical response, 1/h) ----------
    # Placebo K01 and K10 were fixed from the placebo-only fit (Supplement
    # Methods). Named lambda<from><to> because k12 is the PK micro-constant.
    llambda01 <- fixed(log(0.000123))
    label("Placebo rate no response -> response (1/h)") # Table S5 'K01 [1/hour]' clinical response = 0.000123 fix
    llambda10 <- fixed(log(0.000158))
    label("Placebo rate response -> no response (1/h)") # Table S5 'K10 [1/hour]' clinical response = 0.000158 fix
    llambda0d <- log(0.0000912)
    label("Rate no response -> dropout (1/h)") # Table S5 'K0D [1/hour]' clinical response = 0.0000912
    llambda1d <- log(0.0000380)
    label("Rate response -> dropout (1/h)") # Table S5 'K1D [1/hour]' clinical response = 0.0000380
    lscale_lambda <- log(1.98)
    label("Multiplicative factor on both K01 and K10 at all times (unitless)") # Table S5 'Multiplicative factor for K01 and K10 for all times' = 1.98 (footnote f)

    # ---- Drug effect on K01 (Supplement Methods Equation 1) -------------------
    lemax <- log(2.92)
    label("Maximum fractional increase of K01 by upadacitinib (unitless)") # Table S5 'Emax' clinical response = 2.92
    lec50 <- log(15.6)
    label("Upadacitinib plasma concentration giving half of Emax (ng/mL)") # Table S5 'EC50 (ng/mL)' clinical response = 15.6
    lhill <- fixed(log(1))
    label("Hill coefficient of the concentration-effect relationship (unitless)") # Table S7 example control stream (clinical remission) 'HILL=THETA(9)' with THETA(9) = 1 FIX; same template assumed
    e_early_lambda01 <- 3.60
    label("Fractional increase of K01 during the first two weeks (unitless)") # Table S5 'Increase parameter for K01 in first 2 weeks' = 3.60; footnote c: K01 x (1 + theta1) for t <= 2 weeks

    # ---- PK interindividual variability (Table S3, CD values, fixed) -----------
    etalcl ~ fixed(0.169744) # Table S3 'ISV CL in subjects with CD (%)' = 41.2; 0.412^2
    etalvc ~ fixed(0.081796) # Table S3 'ISV Vc in subjects with RA and CD (%)' = 28.6; 0.286^2
    etalka ~ fixed(1.2769) # Table S3 'ISV KA (%)' = 113; 1.13^2
    etalka_er ~ fixed(1.2321) # Table S4 'KA ER' %ISV = 111; 1.11^2
    etaltlag_er ~ fixed(0.218089) # Table S4 'ALAG1 ER' %ISV = 46.7; 0.467^2
    etalogitffo ~ fixed(1.49) # Table S4 variance of the ISV on the first-order fraction = 1.49
    etald2 ~ fixed(1.4641) # Table S4 'D2' %ISV = 121; 1.21^2
    etalfrel_er ~ fixed(0.0289) # Table S4 relative bioavailability %ISV = 17.0; 0.170^2

    # Placeholder residual: the Markov likelihood is categorical (observed state
    # at each visit) and has no residual error.
    addSd_prob_clinresp <- fixed(0.001)
    label("Placeholder additive residual SD on prob_clinresp (not from source)") # not from source; see vignette Assumptions and deviations
  })

  model({
    # ---- PK driver ------------------------------------------------------------
    cl <- exp(lcl + etalcl + e_sexf_cl * SEXF + e_crcl_cl * log(CRCL / 106.9))
    vc <- exp(lvc + etalvc + e_sexf_vc * SEXF + e_wt_vc * log(WT / 74))
    vp <- exp(lvp)
    q <- exp(lq)

    ka_ir <- exp(lka + etalka)
    ka_er <- exp(lka_er + etalka_er)
    tlag_er <- exp(ltlag_er + etaltlag_er)
    ffo_er <- expit(logitffo + etalogitffo)
    frel_er <- exp(lfrel_er + etalfrel_er)

    ka <- ka_er * FORM_UPA_ER + ka_ir * (1 - FORM_UPA_ER)
    tlag <- tlag_er * FORM_UPA_ER + exp(ltlag) * (1 - FORM_UPA_ER)
    frel <- frel_er * FORM_UPA_ER + (1 - FORM_UPA_ER)
    ffo <- ffo_er * FORM_UPA_ER + (1 - FORM_UPA_ER)
    d2 <- exp(ld2 + etald2)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - (kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    f(depot) <- frel * ffo
    f(central) <- frel * (1 - ffo)
    dur(central) <- d2
    alag(depot) <- tlag
    alag(central) <- tlag

    # mg / L x 1000 = ng/mL, the unit of EC50.
    Cc <- 1000 * central / vc

    # ---- Drug effect ------------------------------------------------------------
    # Same form as the Table S7 example control stream (clinical remission):
    # EFF = 1 + EMAX*(C+minc)**HILL / ((IC50+minc)**HILL + (C+minc)**HILL)
    # with minc = 0.001 ng/mL; here C is the time-varying concentration.
    minc <- 0.001
    emax <- exp(lemax)
    ec50 <- exp(lec50)
    hill <- exp(lhill)
    eff <- 1 + emax * (Cc + minc)^hill / ((ec50 + minc)^hill + (Cc + minc)^hill)

    # K01 is multiplied by (1 + theta1) up to and including 336 h (Table S5
    # footnote c; control stream IF(T.LE.336)), with time 0 the first dose.
    early <- 0
    if (t <= 336) early <- 1
    scale_lambda <- exp(lscale_lambda)
    lambda01 <- exp(llambda01) * eff * (1 + e_early_lambda01 * early) * scale_lambda
    lambda10 <- exp(llambda10) * scale_lambda
    lambda0d <- exp(llambda0d)
    lambda1d <- exp(llambda1d)

    # ---- State-occupancy probabilities (Figure S3) -------------------------------
    # Everyone starts as a non-responder; dropout is absorbing.
    d/dt(s_noresp) <- -(lambda01 + lambda0d) * s_noresp + lambda10 * s_resp
    d/dt(s_resp) <- lambda01 * s_noresp - (lambda10 + lambda1d) * s_resp
    d/dt(s_dropout) <- lambda0d * s_noresp + lambda1d * s_resp
    s_noresp(0) <- 1

    # Responders who have not dropped out: the nonresponder-imputation (NRI)
    # proportion the paper simulates (Figure 4).
    prob_clinresp <- s_resp
    prob_clinresp ~ add(addSd_prob_clinresp)
  })
}
