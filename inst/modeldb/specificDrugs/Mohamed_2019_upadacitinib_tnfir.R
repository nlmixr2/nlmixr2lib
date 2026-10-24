Mohamed_2019_upadacitinib_tnfir <- function() {
  description <- paste(
    "Continuous-time Markov exposure-response model for ACR20 / ACR50 /",
    "ACR70 response and dropout in adults with rheumatoid arthritis and an",
    "inadequate response to anti-TNF therapy (anti-TNF-IR, phase IIb",
    "BALANCE I) treated with immediate-release upadacitinib on background",
    "methotrexate (Mohamed 2019). Five states (no response, ACR20, ACR50,",
    "ACR70, dropout) with first-order transitions; the upadacitinib plasma",
    "concentration multiplies the three upward transition rates by",
    "1 + Emax * Cp^hill / (EC50^hill + Cp^hill), with Emax scaled by a power",
    "of baseline hsCRP and raised by (1 + tvtshift) during the first two",
    "weeks of treatment. The concentration driver is the Klunder 2017",
    "two-compartment population PK model for patients with rheumatoid",
    "arthritis, which the paper used to generate each patient's",
    "concentration-time profile."
  )
  reference <- paste(
    "Mohamed M-EF, Klunder B, Camp HS, Othman AA. Exposure-Response Analyses",
    "of Upadacitinib Efficacy in Phase II Trials in Rheumatoid Arthritis and",
    "Basis for Phase III Dose Selection. Clin Pharmacol Ther.",
    "2019;106(6):1319-1327. doi:10.1002/cpt.1543. PK driver: Klunder B,",
    "Mohamed M-EF, Othman AA. Population pharmacokinetics of upadacitinib in",
    "healthy subjects and subjects with rheumatoid arthritis: analyses of",
    "phase I and II clinical trials. Clin Pharmacokinet. 2018;57(8):977-988.",
    "doi:10.1007/s40262-017-0605-6"
  )
  vignette <- "Mohamed_2019_upadacitinib"
  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "ng/mL (Cc); prob_acr20 / prob_acr50 / prob_acr70 are probabilities (0-1)"
  )

  compartmentData <- list(
    depot = list(analyte = "upadacitinib", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "upadacitinib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "upadacitinib", units = "mg", specimen = "plasma", verified = TRUE),
    s_noresp = list(
      analyte = "probability of being an ACR non-responder (below ACR20) and on study",
      units = NA_character_,
      specimen = "not applicable",
      verified = TRUE
    ),
    s_acr20 = list(
      analyte = "probability of being an ACR20 (but not ACR50) responder",
      units = NA_character_,
      specimen = "not applicable",
      verified = TRUE
    ),
    s_acr50 = list(
      analyte = "probability of being an ACR50 (but not ACR70) responder",
      units = NA_character_,
      specimen = "not applicable",
      verified = TRUE
    ),
    s_acr70 = list(
      analyte = "probability of being an ACR70 responder",
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
      notes = "PK driver only: power effect on Vc/F (exponent 0.50, reference 74 kg; Klunder 2017 Table 3 footnote d). BALANCE I mean 78 kg (range 42-134; Mohamed 2019 Table 1).",
      source_name = "WT"
    ),
    SEXF = list(
      description = "Biological sex indicator (1 = female, 0 = male).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male), the Klunder 2017 PK reference.",
      notes = "PK driver only: women have 14% lower CL/F and 25% lower Vc/F (Klunder 2017 Table 3). Sex was not retained in the BALANCE I Markov model (Mohamed 2019 Table S1).",
      source_name = "SEX"
    ),
    CRCL = list(
      description = "Baseline creatinine clearance by Cockcroft-Gault (raw mL/min, not BSA-normalized).",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = "PK driver only: power effect on CL/F (exponent 0.32, reference 107 mL/min; Klunder 2017 Table 3 footnote c). Same raw Cockcroft-Gault column as Klunder_2017_upadacitinib.",
      source_name = "CrCL"
    ),
    CRP = list(
      description = "Baseline high-sensitivity C-reactive protein (hsCRP), time-fixed.",
      units = "mg/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on Emax with exponent 0.276 (Table S1): patients with higher baseline hsCRP have a larger maximal drug effect. The paper does not print the centring value. It is taken as 13.4 mg/L, the BALANCE I baseline mean (Mohamed 2019 Table 1; SD 18.6, range 0.08-135). See the vignette for the sensitivity of the week-12 response to this choice.",
      source_name = "baseline hsCRP"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 276L,
    n_studies = 1L,
    age_range = "26-88 years (mean 57)",
    weight_range = "42-134 kg (mean 78)",
    sex_female_pct = 80,
    race_ethnicity = c(White = 89, Black = 8, Asian = 1, Other = 1),
    disease_state = "Moderate-to-severe active rheumatoid arthritis with inadequate response to at least one anti-TNF biologic, on stable background methotrexate (BALANCE I, NCT01960855)",
    dose_range = "Placebo or immediate-release upadacitinib 3, 6, 12 or 18 mg twice daily for 12 weeks",
    regions = "Multinational",
    baseline_das28crp = "5.77 +/- 0.92 (range 3.58-7.97)",
    baseline_crp = "13.4 +/- 18.6 mg/L (range 0.08-135)",
    anti_ccp_positive_pct = 83,
    notes = "Mohamed 2019 Table 1 (BALANCE I column). ACR visits at weeks 2, 4, 6, 8 and 12. Exposure-response fit in NONMEM 7.3 with individual Klunder 2017 PK profiles as input."
  )

  ini({
    # ---- PK driver: Klunder 2017 Table 3, rheumatoid-arthritis cohort -------
    # Same values as Klunder_2017_upadacitinib; only the RA arm is carried
    # because every BALANCE I patient had RA. Reference: male, WT 74 kg,
    # CrCL 107 mL/min.
    lka <- log(12.3); label("First-order absorption rate constant ka (1/h)") # Klunder 2017 Table 3: Ka = 12.3 1/h
    lcl <- log(39.7 * 0.76); label("Apparent oral clearance CL/F in a male RA patient at CrCL 107 mL/min (L/h)") # Klunder 2017 Table 3: CL/F = 39.7 L/h x theta1 (RA) = 0.76
    lvc <- log(146); label("Apparent central volume Vc/F in a male at WT 74 kg (L)") # Klunder 2017 Table 3: Vc/F = 146 L
    lvp <- log(64.3); label("Apparent peripheral volume Vp/F (L)") # Klunder 2017 Table 3: Vp/F = 64.3 L
    lq <- log(3.23); label("Apparent intercompartmental clearance Q/F (L/h)") # Klunder 2017 Table 3: Q/F = 3.23 L/h
    ltlag <- log(0.48); label("Absorption lag time (h)") # Klunder 2017 Table 3: ALAG1 = 0.48 h
    e_sexf_cl <- log(0.86); label("Log ratio of CL/F, female vs male") # Klunder 2017 Table 3: theta2 = 0.86
    e_crcl_cl <- 0.32; label("Power exponent of CrCL on CL/F") # Klunder 2017 Table 3: theta3 = 0.32
    e_sexf_vc <- log(0.75); label("Log ratio of Vc/F, female vs male") # Klunder 2017 Table 3: theta4 = 0.75
    e_wt_vc <- 0.50; label("Power exponent of body weight on Vc/F") # Klunder 2017 Table 3: theta5 = 0.50

    # ---- Markov transition rates (Mohamed 2019 Table S1, 1/day) --------------
    # The paper's K01 ... K3D. Named lambdaXY here because k12 is the PK
    # micro-constant; model() converts from 1/day to 1/h.
    llambda01 <- log(0.0209); label("Rate no response -> ACR20, placebo (1/day)") # Table S1: K01 = 0.0209 1/day
    llambda10 <- log(0.0339); label("Rate ACR20 -> no response (1/day)") # Table S1: K10 = 0.0339 1/day
    llambda12 <- log(0.0309); label("Rate ACR20 -> ACR50, placebo (1/day)") # Table S1: K12 = 0.0309 1/day
    llambda21 <- log(0.0589); label("Rate ACR50 -> ACR20 (1/day)") # Table S1: K21 = 0.0589 1/day
    llambda23 <- log(0.0269); label("Rate ACR50 -> ACR70, placebo (1/day)") # Table S1: K23 = 0.0269 1/day
    llambda32 <- log(0.0263); label("Rate ACR70 -> ACR50 (1/day)") # Table S1: K32 = 0.0263 1/day
    llambda0d <- log(0.00269); label("Rate no response -> dropout (1/day)") # Table S1: K0D = 0.00269 1/day
    llambda1d <- log(0.000246); label("Rate ACR20 -> dropout (1/day)") # Table S1: K1D = 0.000246 1/day
    llambda2d <- log(0.000871); label("Rate ACR50 -> dropout (1/day)") # Table S1: K2D = 0.000871 1/day
    llambda3d <- fixed(log(1e-10)); label("Rate ACR70 -> dropout (1/day)") # Table S1: K3D = 1E-10 (fixed)

    # ---- Drug effect on the upward transitions (Methods equation) -----------
    lemax <- log(0.766); label("Maximum fractional increase of K01, K12 and K23 at the reference hsCRP, after week 2 (unitless)") # Table S1: Emax = 0.766
    lec50 <- log(12.2); label("Upadacitinib plasma concentration giving half of Emax (ng/mL)") # Table S1: EC50 = 12.2 ng/mL
    lhill <- log(3.46); label("Hill coefficient of the concentration-effect relationship (unitless)") # Table S1: Hill factor = 3.46
    e_early_emax <- 2.15; label("Fractional increase of Emax during the first two weeks, tvtshift (unitless)") # Table S1: tvtshift = 2.15, applied as Emax x (1 + tvtshift) (footnote b)
    e_crp_emax <- 0.276; label("Power exponent of baseline hsCRP on Emax") # Table S1: 'Exponent for the Effect of Baseline hs-CRP on Emax' = 0.276

    # ---- Between-subject variability: PK driver only -------------------------
    # Klunder 2017 Table 3 RA-cohort %ISV, variance = (%ISV/100)^2. The Markov
    # model has no random effects (Table S1 lists none).
    etalka ~ 2.25 # Klunder 2017 Table 3: ISV Ka = 150%
    etalcl ~ 0.0676 # Klunder 2017 Table 3: ISV CL/F, RA = 26%
    etalvc ~ 0.0729 # Klunder 2017 Table 3: ISV Vc/F, RA = 27%

    # Placeholder residuals: the Markov likelihood is categorical (observed
    # state at each visit) and has no residual error.
    addSd_prob_acr20 <- fixed(0.001); label("Placeholder additive residual SD on prob_acr20 (not from source)") # not from source; see vignette Assumptions and deviations
    addSd_prob_acr50 <- fixed(0.001); label("Placeholder additive residual SD on prob_acr50 (not from source)") # not from source; see vignette Assumptions and deviations
    addSd_prob_acr70 <- fixed(0.001); label("Placeholder additive residual SD on prob_acr70 (not from source)") # not from source; see vignette Assumptions and deviations
  })

  model({
    # ---- PK driver (Klunder 2017) -------------------------------------------
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl + e_sexf_cl * SEXF + e_crcl_cl * log(CRCL / 107))
    vc <- exp(lvc + etalvc + e_sexf_vc * SEXF + e_wt_vc * log(WT / 74))
    vp <- exp(lvp)
    q <- exp(lq)
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    alag(depot) <- exp(ltlag)

    # mg / L = ug/mL; x 1000 gives ng/mL, the unit of EC50.
    Cc <- 1000 * central / vc

    # ---- Drug effect --------------------------------------------------------
    # Emax rises by (1 + tvtshift) for the first two weeks (336 h) after the
    # start of treatment at time 0 (Table S1 footnote b), and scales with
    # baseline hsCRP relative to 13.4 mg/L.
    early <- 0
    if (t < 336) early <- 1
    emax <- exp(lemax) * (CRP / 13.4)^e_crp_emax * (1 + e_early_emax * early)
    ec50 <- exp(lec50)
    hill <- exp(lhill)
    drug <- 1 + emax * Cc^hill / (ec50^hill + Cc^hill)

    # ---- Transition rates (1/day in the source, / 24 for 1/h) --------------
    lambda01 <- exp(llambda01) / 24 * drug
    lambda10 <- exp(llambda10) / 24
    lambda12 <- exp(llambda12) / 24 * drug
    lambda21 <- exp(llambda21) / 24
    lambda23 <- exp(llambda23) / 24 * drug
    lambda32 <- exp(llambda32) / 24
    lambda0d <- exp(llambda0d) / 24
    lambda1d <- exp(llambda1d) / 24
    lambda2d <- exp(llambda2d) / 24
    lambda3d <- exp(llambda3d) / 24

    # ---- State-occupancy probabilities (Figure 1) --------------------------
    # Everyone starts as a non-responder; dropout is absorbing.
    d/dt(s_noresp) <- -(lambda01 + lambda0d) * s_noresp + lambda10 * s_acr20
    d/dt(s_acr20) <- lambda01 * s_noresp - (lambda10 + lambda12 + lambda1d) * s_acr20 +
      lambda21 * s_acr50
    d/dt(s_acr50) <- lambda12 * s_acr20 - (lambda21 + lambda23 + lambda2d) * s_acr50 +
      lambda32 * s_acr70
    d/dt(s_acr70) <- lambda23 * s_acr50 - (lambda32 + lambda3d) * s_acr70
    d/dt(s_dropout) <- lambda0d * s_noresp + lambda1d * s_acr20 + lambda2d * s_acr50 +
      lambda3d * s_acr70
    s_noresp(0) <- 1

    # Responder probabilities with dropouts counted as non-responders (NRI):
    # an ACR50 or ACR70 responder is also an ACR20 responder.
    prob_acr20 <- s_acr20 + s_acr50 + s_acr70
    prob_acr50 <- s_acr50 + s_acr70
    prob_acr70 <- s_acr70

    prob_acr20 ~ add(addSd_prob_acr20)
    prob_acr50 ~ add(addSd_prob_acr50)
    prob_acr70 ~ add(addSd_prob_acr70)
  })
}
