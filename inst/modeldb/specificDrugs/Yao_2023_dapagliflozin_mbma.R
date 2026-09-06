Yao_2023_dapagliflozin_mbma <- function() {
  description <- paste0(
    "MBMA. Population pharmacokinetic model-based meta-analysis of ",
    "dapagliflozin: two transit-absorption compartments feeding a ",
    "two-compartment disposition model with first-order elimination. Fitted ",
    "to 880 pooled summary-level plasma concentrations digitised from 34 ",
    "published dapagliflozin, canagliflozin and empagliflozin clinical ",
    "studies (this file holds the dapagliflozin arm; 201 subjects, 177 of ",
    "them healthy, doses 2.5-250 mg). Food slows absorption: the fed state ",
    "multiplies the transit rate constant Kt by 0.254. Variability is ",
    "BETWEEN-STUDY (inter-study), encoded as study-level etas, so the model ",
    "simulates study-arm summary profiles and is NOT suitable for ",
    "individual-subject simulation. The companion PK/PD/endpoint model is ",
    "modellib('Yao_2023_sglt2_endpoints_mbma'), which consumes this model's ",
    "steady-state AUC(0-24 h) as the AUC_DAPA covariate."
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

  covariateData <- list(
    FED = list(
      description        = "Fed-versus-fasted state at the time of dosing.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (fasted)",
      notes              = paste0(
        "Study-arm-level indicator in this MBMA (a pooled study arm was ",
        "dosed either fed or fasted). Multiplicative effect on the transit ",
        "rate constant, kt * 0.254^FED, i.e. food cuts Kt to about a ",
        "quarter of its fasted value and prolongs Tmax. Yao 2023 Table 1 ",
        "reports the coefficient as a bare 'Fed' covariate value of 0.254 ",
        "(RSE 29.8%) without printing the covariate equation; the ",
        "multiplicative reading is the only one consistent with the ",
        "Discussion statement that food 'could significantly influence the ",
        "Kt parameter (p < 0.005) in dapagliflozin, which is consistent ",
        "with a longer time to maximum concentration of dapagliflozin ",
        "under fed condition'. See the vignette Assumptions and deviations ",
        "section."
      ),
      source_name        = "Fed"
    )
  )

  compartmentData <- list(
    depot       = list(analyte = "dapagliflozin", units = "mg", specimen = "administration site", verified = TRUE),
    transit1    = list(analyte = "dapagliflozin", units = "mg", specimen = "administration site", verified = TRUE),
    transit2    = list(analyte = "dapagliflozin", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "dapagliflozin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "dapagliflozin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 201L,
    n_studies      = 8L,
    age_range      = "mean 38.1 years (SD 9.78)",
    weight_range   = "mean 76.38 kg (SD 14.58)",
    sex_female_pct = 25.68,
    disease_state  = paste0(
      "Pooled healthy subjects and patients with type 2 diabetes mellitus ",
      "(177 of 201 healthy). Studies in patients with moderate or severe ",
      "renal impairment or hepatic insufficiency were excluded; the PK/PD ",
      "data come from subjects with glomerular filtration rate above ",
      "60 mL/min/1.73 m2."
    ),
    dose_range     = "2.5, 5, 10, 20, 50, 100, 250 mg oral",
    regions        = "International (published clinical trials indexed on PubMed to July 2016)",
    notes          = paste0(
      "Model-based meta-analysis: the unit of observation is a published ",
      "summary-level (study-arm mean) plasma concentration, not an ",
      "individual measurement. Demographics from Yao 2023 Table S2 (PK ",
      "row for dapagliflozin); study list and dose ranges from Table S1 ",
      "(citations S28-S35). BMI mean 27.02 kg/m2 (SD 3.3). Male 74.32% ",
      "(SD 17.15), so the female percentage is a pooled study-level ",
      "complement rather than a subject count."
    )
  )

  # ==========================================================================
  # ini(): Yao 2023 Table 1 (PK model of dapagliflozin).
  #
  # Structure: Appendix S1 Equations S1-S5 (two transit compartments feeding a
  # two-compartment disposition model). Kt is both the dosing-compartment
  # emptying rate and the transit rate constant, so the whole absorption chain
  # runs at one rate.
  #
  # Variability is BETWEEN-STUDY, not between-subject (Appendix S1 section 1.2:
  # "Interstudy variability of all PK parameters was assumed to be
  # log-normally distributed and was described using exponential model",
  # Equation S13 P_ij = theta_i * exp(eta_ij)). Table 1's "IIV (%)" column is
  # therefore read as an inter-study CV% and converted with
  # omega^2 = log(CV^2 + 1). Etas are named eta_study_* so they cannot be
  # mistaken for the popPK between-subject convention.
  #
  # CL, Vc, CLD and VT are APPARENT (CL/F, Vc/F, ...): the paper fits only
  # oral data and reports no bioavailability term, so F is folded into them.
  # ==========================================================================
  ini({
    lcl <- log(19.5)
    label("Apparent clearance from the central compartment (CL/F, L/h)")  # Yao 2023 Table 1 (RSE 4.00%)

    lvc <- log(82.0)
    label("Apparent central volume of distribution (Vc/F, L)")  # Yao 2023 Table 1 (RSE 6.90%)

    lq <- log(10.3)
    label("Apparent distribution clearance between central and peripheral compartments (CLD/F, L/h)")  # Yao 2023 Table 1 (RSE 7.10%)

    lvp <- log(122)
    label("Apparent peripheral volume of distribution (VT/F, L)")  # Yao 2023 Table 1 (RSE 9.60%)

    lktr <- log(6.50)
    label("Transit / absorption rate constant Kt (1/h)")  # Yao 2023 Table 1 (RSE 12.4%)

    e_fed_ktr <- 0.254
    label("Multiplicative effect of the fed state on Kt (unitless)")  # Yao 2023 Table 1 "Fed" (RSE 29.8%)

    # ----- Between-study (inter-study) variability; exponential, diagonal ----
    eta_study_lcl  ~ 0.000624805  # Yao 2023 Table 1; log(0.0250^2 + 1)
    eta_study_lvc  ~ 0.00887189   # Yao 2023 Table 1; log(0.0944^2 + 1)
    eta_study_lvp  ~ 0.00339313   # Yao 2023 Table 1; log(0.0583^2 + 1)
    eta_study_lktr ~ 0.0515483    # Yao 2023 Table 1; log(0.230^2 + 1)
    # No inter-study variability is reported for CLD (Table 1 IIV column "-").

    # ----- Residual error (Appendix S1 Equation S14) ------------------------
    # Y_obs = Y_pred + Y_pred * eps1 * W + eps2 * W, W = sqrt(study sample
    # size). The values below are the reported sigma^2 back-transformed to SD
    # at W = 1; a downstream user simulating an N-subject study arm must apply
    # the paper's weight. The DIRECTION of the weight is not resolvable from
    # on-disk sources -- see the vignette Assumptions and deviations section.
    propSd <- 0.676018
    label("Proportional residual SD at unit study weight (fraction)")  # Yao 2023 Table 1; sqrt(0.457)

    addSd <- 0.679706
    label("Additive residual SD at unit study weight (ng/mL)")  # Yao 2023 Table 1; sqrt(0.462)
  })

  model({
    # ---- Study-level (inter-study-perturbed) parameters --------------------
    cl <- exp(lcl + eta_study_lcl)
    vc <- exp(lvc + eta_study_lvc)
    q  <- exp(lq)
    vp <- exp(lvp + eta_study_lvp)

    # Food slows absorption: Kt is multiplied by 0.254 in the fed state.
    ktr <- exp(lktr + eta_study_lktr) * e_fed_ktr^FED

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ---- Appendix S1 Equations S1-S5 ---------------------------------------
    # The dosing compartment empties at Kt, each transit compartment passes on
    # at Kt, and the second transit compartment feeds the central compartment
    # at Kt.
    d/dt(depot)       <- -ktr * depot
    d/dt(transit1)    <-  ktr * depot    - ktr * transit1
    d/dt(transit2)    <-  ktr * transit1 - ktr * transit2
    d/dt(central)     <-  ktr * transit2 - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central  - k21 * peripheral1

    # Amounts are in mg and volumes in L, so central/vc is mg/L; x 1000 gives
    # ng/mL, the scale on which EC50 and the residual error are reported.
    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
