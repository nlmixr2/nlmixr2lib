Yao_2023_empagliflozin_mbma <- function() {
  description <- paste0(
    "MBMA. Population pharmacokinetic model-based meta-analysis of ",
    "empagliflozin: two transit-absorption compartments feeding a ",
    "two-compartment disposition model with first-order elimination. Fitted ",
    "to pooled summary-level plasma concentrations digitised from published ",
    "empagliflozin clinical studies (204 subjects, 122 of them healthy, ",
    "doses 1-100 mg). No food effect was retained. Variability is ",
    "BETWEEN-STUDY (inter-study), encoded as study-level etas, so the model ",
    "simulates study-arm summary profiles and is NOT suitable for ",
    "individual-subject simulation. The companion PK/PD/endpoint model is ",
    "modellib('Yao_2023_sglt2_endpoints_mbma'), which consumes this model's ",
    "steady-state AUC(0-24 h) as the AUC_EMPA covariate."
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
    depot       = list(analyte = "empagliflozin", units = "mg", specimen = "administration site", verified = TRUE),
    transit1    = list(analyte = "empagliflozin", units = "mg", specimen = "administration site", verified = TRUE),
    transit2    = list(analyte = "empagliflozin", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "empagliflozin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "empagliflozin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 204L,
    n_studies      = 7L,
    age_range      = "mean 43.52 years (SD 15.0)",
    weight_range   = "mean 73.3 kg (SD 9.18)",
    sex_female_pct = 36.84,
    disease_state  = paste0(
      "Pooled healthy subjects and patients with type 2 diabetes mellitus ",
      "(122 of 204 healthy). Studies in patients with moderate or severe ",
      "renal impairment or hepatic insufficiency were excluded; the PK/PD ",
      "data come from subjects with glomerular filtration rate above ",
      "60 mL/min/1.73 m2."
    ),
    dose_range     = "1, 5, 10, 25, 50, 100 mg oral",
    regions        = "International (published clinical trials indexed on PubMed to July 2016)",
    notes          = paste0(
      "Model-based meta-analysis: the unit of observation is a published ",
      "summary-level (study-arm mean) plasma concentration, not an ",
      "individual measurement. Demographics from Yao 2023 Table S2 (PK row ",
      "for empagliflozin); study list and dose ranges from Table S1 ",
      "(citations S83-S89). BMI mean 24.84 kg/m2 (SD 3.01). Male 63.16% ",
      "(SD 27.53), so the female percentage is a pooled study-level ",
      "complement rather than a subject count."
    )
  )

  # ==========================================================================
  # ini(): Yao 2023 Table 1 (PK model of empagliflozin).
  #
  # Structure: Appendix S1 Equations S1-S5 (two transit compartments feeding a
  # two-compartment disposition model), the same chain length as
  # dapagliflozin. Kt is both the dosing-compartment emptying rate and the
  # transit rate constant.
  #
  # Variability is BETWEEN-STUDY, not between-subject (Appendix S1 section 1.2,
  # Equation S13 P_ij = theta_i * exp(eta_ij)). Table 1's "IIV (%)" column is
  # read as an inter-study CV% and converted with omega^2 = log(CV^2 + 1).
  #
  # CL, Vc, CLD and VT are APPARENT (CL/F, Vc/F, ...): the paper fits only
  # oral data and reports no bioavailability term, so F is folded into them.
  # ==========================================================================
  ini({
    lcl <- log(4.25)
    label("Apparent clearance from the central compartment (CL/F, L/h)")  # Yao 2023 Table 1 (RSE 6.10%)

    lvc <- log(30.6)
    label("Apparent central volume of distribution (Vc/F, L)")  # Yao 2023 Table 1 (RSE 9.50%)

    lq <- log(1.37)
    label("Apparent distribution clearance between central and peripheral compartments (CLD/F, L/h)")  # Yao 2023 Table 1 (RSE 13.5%)

    lvp <- log(28.3)
    label("Apparent peripheral volume of distribution (VT/F, L)")  # Yao 2023 Table 1 (RSE 26.5%)

    lktr <- log(4.13)
    label("Transit / absorption rate constant Kt (1/h)")  # Yao 2023 Table 1 (RSE 6.90%)

    # ----- Between-study (inter-study) variability; exponential, diagonal ----
    eta_study_lcl  ~ 0.00323238  # Yao 2023 Table 1; log(0.0569^2 + 1)
    eta_study_lvc  ~ 0.00316468  # Yao 2023 Table 1; log(0.0563^2 + 1)
    eta_study_lvp  ~ 0.0472652   # Yao 2023 Table 1; log(0.220^2 + 1)
    eta_study_lktr ~ 0.00320973  # Yao 2023 Table 1; log(0.0567^2 + 1)
    # No inter-study variability is reported for CLD (Table 1 IIV column "-").

    # ----- Residual error (Appendix S1 Equation S14) ------------------------
    # Y_obs = Y_pred + Y_pred * eps1 * W + eps2 * W, W = sqrt(study sample
    # size). Table 1 reports a proportional variance only for empagliflozin
    # (no additive row), so the model carries a proportional term alone. The
    # value below is at W = 1; see the vignette Assumptions and deviations
    # section for the weighting.
    propSd <- 0.721803
    label("Proportional residual SD at unit study weight (fraction)")  # Yao 2023 Table 1; sqrt(0.521)
  })

  model({
    # ---- Study-level (inter-study-perturbed) parameters --------------------
    cl  <- exp(lcl + eta_study_lcl)
    vc  <- exp(lvc + eta_study_lvc)
    q   <- exp(lq)
    vp  <- exp(lvp + eta_study_lvp)
    ktr <- exp(lktr + eta_study_lktr)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ---- Appendix S1 Equations S1-S5 ---------------------------------------
    d/dt(depot)       <- -ktr * depot
    d/dt(transit1)    <-  ktr * depot    - ktr * transit1
    d/dt(transit2)    <-  ktr * transit1 - ktr * transit2
    d/dt(central)     <-  ktr * transit2 - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central  - k21 * peripheral1

    # Amounts are in mg and volumes in L, so central/vc is mg/L; x 1000 gives
    # ng/mL, the scale on which EC50 and the residual error are reported.
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
