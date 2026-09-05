Yao_2023_canagliflozin_mbma <- function() {
  description <- paste0(
    "MBMA. Population pharmacokinetic model-based meta-analysis of ",
    "canagliflozin: FOUR transit-absorption compartments feeding a ",
    "two-compartment disposition model with first-order elimination (twice ",
    "the transit chain length used for dapagliflozin and empagliflozin in ",
    "the same paper). Fitted to pooled summary-level plasma concentrations ",
    "digitised from published canagliflozin clinical studies (158 subjects, ",
    "76 of them healthy, doses 25-400 mg). No food effect was retained. ",
    "Variability is BETWEEN-STUDY (inter-study), encoded as study-level ",
    "etas, so the model simulates study-arm summary profiles and is NOT ",
    "suitable for individual-subject simulation. The companion ",
    "PK/PD/endpoint model is modellib('Yao_2023_sglt2_endpoints_mbma'), ",
    "which consumes this model's steady-state AUC(0-24 h) as the AUC_CANA ",
    "covariate."
  )

  reference <- paste(
    "Yao X, Zhou J, Song L, Ren Y, Hu P, Liu D.",
    "A model-based meta analysis study of sodium glucose co-transporter-2",
    "inhibitors. CPT Pharmacometrics Syst Pharmacol. 2023;12(4):487-499.",
    "doi:10.1002/psp4.12934.",
    sep = " "
  )
  vignette <- "Yao_2023_sglt2_inhibitors"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot       = list(analyte = "canagliflozin", units = "mg", specimen = "administration site", verified = TRUE),
    transit1    = list(analyte = "canagliflozin", units = "mg", specimen = "administration site", verified = TRUE),
    transit2    = list(analyte = "canagliflozin", units = "mg", specimen = "administration site", verified = TRUE),
    transit3    = list(analyte = "canagliflozin", units = "mg", specimen = "administration site", verified = TRUE),
    transit4    = list(analyte = "canagliflozin", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "canagliflozin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "canagliflozin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 158L,
    n_studies      = 6L,
    age_range      = "mean 46.82 years (SD 11.15)",
    weight_range   = "mean 75.58 kg (SD 6.78)",
    sex_female_pct = 27.15,
    disease_state  = paste0(
      "Pooled healthy subjects and patients with type 2 diabetes mellitus ",
      "(76 of 158 healthy). Studies in patients with moderate or severe ",
      "renal impairment or hepatic insufficiency were excluded; the PK/PD ",
      "data come from subjects with glomerular filtration rate above ",
      "60 mL/min/1.73 m2."
    ),
    dose_range     = "25, 50, 100, 188, 200, 300, 400 mg oral",
    regions        = "International (published clinical trials indexed on PubMed to July 2016)",
    notes          = paste0(
      "Model-based meta-analysis: the unit of observation is a published ",
      "summary-level (study-arm mean) plasma concentration, not an ",
      "individual measurement. Demographics from Yao 2023 Table S2 (PK row ",
      "for canagliflozin); study list and dose ranges from Table S1 ",
      "(citations S61-S66). BMI mean 26.44 kg/m2 (SD 3.02). Male 72.85% ",
      "(SD 19.18), so the female percentage is a pooled study-level ",
      "complement rather than a subject count."
    )
  )

  # ==========================================================================
  # ini(): Yao 2023 Table 1 (PK model of canagliflozin).
  #
  # Structure: Appendix S1 Equations S6-S12 -- FOUR transit compartments,
  # selected because they "best describe the absorption delay of
  # canagliflozin" (Appendix S1 section 1.1; Yao 2023 Discussion: "Two transit
  # compartments were applied for dapagliflozin and empagliflozin and four
  # transit compartments for canagliflozin"). Kt is both the dosing-compartment
  # emptying rate and the transit rate constant.
  #
  # Variability is BETWEEN-STUDY, not between-subject (Appendix S1 section 1.2,
  # Equation S13 P_ij = theta_i * exp(eta_ij)). Table 1's "IIV (%)" column is
  # read as an inter-study CV% and converted with omega^2 = log(CV^2 + 1).
  #
  # CL, Vc, CLD and VT are APPARENT (CL/F, Vc/F, ...): the paper fits only
  # oral data and reports no bioavailability term, so F is folded into them.
  # ==========================================================================
  ini({
    lcl <- log(12.0)
    label("Apparent clearance from the central compartment (CL/F, L/h)")  # Yao 2023 Table 1 (RSE 6.20%)

    lvc <- log(85.5)
    label("Apparent central volume of distribution (Vc/F, L)")  # Yao 2023 Table 1 (RSE 4.60%)

    lq <- log(9.77)
    label("Apparent distribution clearance between central and peripheral compartments (CLD/F, L/h)")  # Yao 2023 Table 1 (RSE 10.8%)

    lvp <- log(108)
    label("Apparent peripheral volume of distribution (VT/F, L)")  # Yao 2023 Table 1 (RSE 6.40%)

    lktr <- log(6.38)
    label("Transit / absorption rate constant Kt (1/h)")  # Yao 2023 Table 1 (RSE 4.60%)

    # ----- Between-study (inter-study) variability; exponential, diagonal ----
    eta_study_lcl  ~ 0.00347496   # Yao 2023 Table 1; log(0.0590^2 + 1)
    eta_study_lvc  ~ 0.000590316  # Yao 2023 Table 1; log(0.0243^2 + 1)
    eta_study_lq   ~ 0.0162512    # Yao 2023 Table 1; log(0.128^2 + 1)
    eta_study_lvp  ~ 0.000985474  # Yao 2023 Table 1; log(0.0314^2 + 1)
    eta_study_lktr ~ 0.00100439   # Yao 2023 Table 1; log(0.0317^2 + 1)

    # ----- Residual error (Appendix S1 Equation S14) ------------------------
    # Y_obs = Y_pred + Y_pred * eps1 * W + eps2 * W, W = sqrt(study sample
    # size). Table 1 reports a proportional variance only for canagliflozin
    # (no additive row), so the model carries a proportional term alone. The
    # value below is at W = 1; see the vignette Assumptions and deviations
    # section for the weighting.
    propSd <- 0.354965
    label("Proportional residual SD at unit study weight (fraction)")  # Yao 2023 Table 1; sqrt(0.126)
  })

  model({
    # ---- Study-level (inter-study-perturbed) parameters --------------------
    cl  <- exp(lcl + eta_study_lcl)
    vc  <- exp(lvc + eta_study_lvc)
    q   <- exp(lq  + eta_study_lq)
    vp  <- exp(lvp + eta_study_lvp)
    ktr <- exp(lktr + eta_study_lktr)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ---- Appendix S1 Equations S6-S12 --------------------------------------
    d/dt(depot)       <- -ktr * depot
    d/dt(transit1)    <-  ktr * depot    - ktr * transit1
    d/dt(transit2)    <-  ktr * transit1 - ktr * transit2
    d/dt(transit3)    <-  ktr * transit2 - ktr * transit3
    d/dt(transit4)    <-  ktr * transit3 - ktr * transit4
    d/dt(central)     <-  ktr * transit4 - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central  - k21 * peripheral1

    # Amounts are in mg and volumes in L, so central/vc is mg/L; x 1000 gives
    # ng/mL, the scale on which EC50 and the residual error are reported.
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
