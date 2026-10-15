Mohamed_2019_upadacitinib_mtxir <- function() {
  description <- paste(
    "Continuous-time Markov exposure-response model for ACR20 / ACR50 /",
    "ACR70 response and dropout in adults with rheumatoid arthritis and an",
    "inadequate response to methotrexate (MTX-IR, phase IIb BALANCE II)",
    "treated with immediate-release upadacitinib (Mohamed 2019). Five",
    "states (no response, ACR20, ACR50, ACR70, dropout) with first-order",
    "transitions; the upadacitinib plasma concentration multiplies the three",
    "upward transition rates by 1 + Emax * Cp / (EC50 + Cp). Placebo-rate",
    "covariates: sex on the no-response -> ACR20 rate, and baseline anti-CCP",
    "status and baseline DAS28-CRP on the ACR20 -> no-response rate. The",
    "concentration driver is the Klunder 2017 two-compartment population PK",
    "model for patients with rheumatoid arthritis, which the paper used to",
    "generate each patient's concentration-time profile."
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
      notes = "PK driver only: power effect on Vc/F (exponent 0.50, reference 74 kg; Klunder 2017 Table 3 footnote d). BALANCE II mean 75 kg (range 44-122; Mohamed 2019 Table 1).",
      source_name = "WT"
    ),
    SEXF = list(
      description = "Biological sex indicator (1 = female, 0 = male).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male). Mohamed 2019 Table S2 reports the 'Ratio of K01 in females relative to males' (1.7), so the K01 estimate is the male rate; the Klunder 2017 PK driver also uses males as reference.",
      notes = "Two roles: (1) Markov placebo rate, K01 (no response -> ACR20) is 1.7-fold higher in women (Table S2); (2) PK driver, women have 14% lower CL/F and 25% lower Vc/F (Klunder 2017 Table 3). The male-reference reading is confirmed by the paper's own placebo prediction: a BALANCE II-like cohort (79% female, 84% anti-CCP positive) reproduces the printed week-12 placebo NRI response of 47/18/8 percent only under this orientation (see vignette).",
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
    SCORE_DAS28CRP = list(
      description = "Baseline 28-joint Disease Activity Score using C-reactive protein (DAS28-CRP).",
      units = "(score)",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on the Markov K10 rate (ACR20 -> no response) with exponent -2.26 (Table S2), so patients with higher baseline disease activity lose an ACR20 response more slowly. The paper does not print the centring value; it is taken as 5.66, the BALANCE II baseline mean (Mohamed 2019 Table 1; range 3-8.25).",
      source_name = "baseline DAS28-CRP"
    ),
    ANTICCP_POS = list(
      description = "Baseline anti-cyclic citrullinated peptide (anti-CCP) antibody status (1 = positive, 0 = negative).",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (anti-CCP negative). Table S2 reports the 'Ratio of K10 in anti-CCP positive relative to negative status subject' (0.394), so the K10 estimate is the anti-CCP-negative rate.",
      notes = "Anti-CCP-positive patients lose an ACR20 response 61% more slowly (K10 multiplied by 0.394; Table S2). 251 of 298 BALANCE II patients (84%) were anti-CCP positive (Table 1).",
      source_name = "baseline anti-CCP status"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 298L,
    n_studies = 1L,
    age_range = "19-82 years (mean 55)",
    weight_range = "44-122 kg (mean 75)",
    sex_female_pct = 79,
    race_ethnicity = c(White = 98, Black = 1, Asian = 0, Other = 1),
    disease_state = "Moderate-to-severe active rheumatoid arthritis with inadequate response to methotrexate, biologic-naive, on stable background methotrexate (BALANCE II, NCT02066389)",
    dose_range = "Placebo or immediate-release upadacitinib 3, 6, 12 or 18 mg twice daily or 24 mg once daily for 12 weeks",
    regions = "Multinational",
    baseline_das28crp = "5.66 +/- 0.99 (range 3-8.25)",
    baseline_crp = "13.3 +/- 18.2 mg/L",
    anti_ccp_positive_pct = 84,
    notes = "Mohamed 2019 Table 1 (BALANCE II column). 300 patients enrolled, 298 analysed. ACR visits at weeks 2, 4, 6, 8 and 12. Exposure-response fit in NONMEM 7.3 with individual Klunder 2017 PK profiles as input."
  )

  ini({
    # ---- PK driver: Klunder 2017 Table 3, rheumatoid-arthritis cohort -------
    # Same values as Klunder_2017_upadacitinib; only the RA arm is carried
    # because every BALANCE II patient had RA. Reference: male, WT 74 kg,
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

    # ---- Markov transition rates (Mohamed 2019 Table S2, 1/day) --------------
    # The paper's K01 ... K3D. Named lambdaXY here because k12 is the PK
    # micro-constant; model() converts from 1/day to 1/h.
    llambda01 <- log(0.0132); label("Rate no response -> ACR20, male, placebo (1/day)") # Table S2: K01 = 0.0132 1/day
    llambda10 <- log(0.0550); label("Rate ACR20 -> no response, anti-CCP negative, DAS28-CRP 5.66 (1/day)") # Table S2: K10 = 0.0550 1/day
    llambda12 <- log(0.0204); label("Rate ACR20 -> ACR50, placebo (1/day)") # Table S2: K12 = 0.0204 1/day
    llambda21 <- log(0.0407); label("Rate ACR50 -> ACR20 (1/day)") # Table S2: K21 = 0.0407 1/day
    llambda23 <- log(0.0234); label("Rate ACR50 -> ACR70, placebo (1/day)") # Table S2: K23 = 0.0234 1/day
    llambda32 <- log(0.0240); label("Rate ACR70 -> ACR50 (1/day)") # Table S2: K32 = 0.0240 1/day
    llambda0d <- log(0.00120); label("Rate no response -> dropout (1/day)") # Table S2: K0D = 0.00120 1/day
    llambda1d <- log(0.000891); label("Rate ACR20 -> dropout (1/day)") # Table S2: K1D = 0.000891 1/day
    llambda2d <- log(0.000617); label("Rate ACR50 -> dropout (1/day)") # Table S2: K2D = 0.000617 1/day
    llambda3d <- log(0.000468); label("Rate ACR70 -> dropout (1/day)") # Table S2: K3D = 0.000468 1/day

    # ---- Drug effect on the upward transitions (Methods equation) -----------
    lemax <- log(1.19); label("Maximum fractional increase of K01, K12 and K23 (unitless)") # Table S2: Emax = 1.19
    lec50 <- log(2.97); label("Upadacitinib plasma concentration giving half of Emax (ng/mL)") # Table S2: EC50 = 2.97 ng/mL

    # ---- Placebo-rate covariate effects (Table S2) ---------------------------
    e_das28_lambda10 <- -2.26; label("Power exponent of baseline DAS28-CRP on K10") # Table S2: 'Exponent for the Effect of Baseline DAS28-CRP on K10' = -2.26
    e_sexf_lambda01 <- 1.7; label("Ratio of K01 in females relative to males") # Table S2: 'Ratio of K01 in females relative to males' = 1.7
    e_anticcp_lambda10 <- 0.394; label("Ratio of K10 in anti-CCP positive relative to negative") # Table S2: 'Ratio of K10 in anti-CCP positive relative to negative status subject' = 0.394

    # ---- Between-subject variability: PK driver only -------------------------
    # Klunder 2017 Table 3 RA-cohort %ISV, variance = (%ISV/100)^2. The Markov
    # model has no random effects (Table S2 lists none).
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

    # ---- Transition rates (1/day in the source, / 24 for 1/h) --------------
    # Placebo covariates (Table S2): sex on K01, DAS28-CRP and anti-CCP on K10.
    drug <- 1 + exp(lemax) * Cc / (exp(lec50) + Cc)

    lambda01 <- exp(llambda01) / 24 * e_sexf_lambda01^SEXF * drug
    lambda10 <- exp(llambda10) / 24 * (SCORE_DAS28CRP / 5.66)^e_das28_lambda10 *
      e_anticcp_lambda10^ANTICCP_POS
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
