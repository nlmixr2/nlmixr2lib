Saeheng_2024_atractylodesLancea_group2day14 <- function() {
  description <- paste0(
    "One-compartment population PK model with zero-order absorption (no lag) and linear ",
    "clearance for total Atractylodes lancea (Thunb.) DC. bioactivity in patients with ",
    "advanced-stage intrahepatic cholangiocarcinoma (Saeheng 2024, Table 2). This file holds ",
    "the Group 2 day-14 fit: the last day of the 1,000 mg once-daily lead-in period of the ",
    "dose-escalating Group 2 regimen (capsule formulation of standardized AL extract, CMC-AL). ",
    "Saeheng 2024 reported three separate fits of the same structural model to three different ",
    "occasions -- Group 1 day 1, Group 2 day 14 and Group 2 day 28 -- carried here as three ",
    "model files sharing one vignette; see also ",
    "modellib('Saeheng_2024_atractylodesLancea_group1day1') and ",
    "modellib('Saeheng_2024_atractylodesLancea_group2day28'). Volume and clearance are apparent ",
    "(V/F, CL/F): CMC-AL is dosed orally and no bioavailability term was estimated. The analyte ",
    "is total AL bioactivity, a bioassay-derived aggregate of the extract's active constituents, ",
    "not a single molecular entity. None of the screened covariates (sex, age, weight, height) ",
    "improved the model. Between-subject variability is lognormal on Tk0, V/F and CL/F; residual ",
    "error is combined additive plus proportional."
  )
  reference <- paste0(
    "Saeheng T, Karbwang J, Na-Bangchang K. Population-pharmacokinetic/pharmacodynamic model of ",
    "atractylodes lancea (Thunb.) DC. administration in patients with advanced-stage intrahepatic ",
    "cholangiocarcinoma: a dosage prediction. BMC Complement Med Ther. 2024;24(1):384. ",
    "doi:10.1186/s12906-024-04618-8."
  )
  vignette <- "Saeheng_2024_atractylodesLancea"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(
      analyte  = "Atractylodes lancea (Thunb.) DC. total bioactivity",
      units    = "mg",
      specimen = "plasma",
      verified = TRUE
    )
  )

  # Saeheng 2024 Results: "None of the investigated covariates (sex, age, weight,
  # and height) improved the model." They are recorded here for provenance of the
  # covariate screen; none is referenced in model().
  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Female sex indicator (1 = female, 0 = male).",
      units       = "(binary)",
      type        = "binary",
      notes       = paste0(
        "Screened as a candidate covariate in Saeheng 2024 but not retained: 'None of the ",
        "investigated covariates (sex, age, weight, and height) improved the model.' The ",
        "cohort's sex distribution is not tabulated in this paper; it is reported in the parent ",
        "phase 2A trial publication (Saeheng 2024 reference [8])."
      )
    ),
    AGE = list(
      description = "Patient age at study entry.",
      units       = "years",
      type        = "continuous",
      notes       = paste0(
        "Screened but not retained (Saeheng 2024 Results). No age summary is tabulated in this ",
        "paper; see the parent phase 2A trial publication (reference [8])."
      )
    ),
    WT = list(
      description = "Body weight.",
      units       = "kg",
      type        = "continuous",
      notes       = paste0(
        "Screened but not retained (Saeheng 2024 Results). Because weight was not retained, V/F ",
        "and CL/F carry no allometric scaling and are absolute (L, L/h), not per-kg."
      )
    ),
    HT = list(
      description = "Body height.",
      units       = "cm",
      type        = "continuous",
      notes       = "Screened but not retained (Saeheng 2024 Results)."
    )
  )

  population <- list(
    species       = "human",
    n_subjects    = 16L,
    n_studies     = 1L,
    disease_state = "advanced-stage intrahepatic cholangiocarcinoma (iCCA)",
    dose_range    = paste0(
      "Group 2 dose-escalating regimen: 1,000 mg once daily (9 capsules) for 14 days, then ",
      "1,500 mg once daily (14 capsules) for 14 days, then 2,000 mg once daily (18 capsules) for ",
      "62 days, plus standard supportive care; this model is the day-14 fit, i.e. the last day of ",
      "the 1,000 mg once-daily lead-in"
    ),
    regions = "Thailand (Sakhon Na-Kon Hospital, Sakhon Na-Kon)",
    notes   = paste0(
      "Group 2 arm of a single-centre, open-label, randomized, controlled phase 2A trial ",
      "(TCTR20210129007). Sixteen patients were randomized to Group 2 (Saeheng 2024 Results). ",
      "Plasma sampling on day 14 at 0, 0.25, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 6 and 8 h post-dose ",
      "(Saeheng 2024 Methods, 'Pharmacokinetic study'). Baseline demographics (age, weight, sex, ",
      "height) are not tabulated in this paper; they are reported in the parent phase 2A trial ",
      "publication (reference [8]). Fitted in MonolixSuite 2023R1 by SAEM."
    )
  )

  ini({
    # ==================================================================
    # Structural parameters -- Saeheng 2024 Table 2 ("Pharmacokinetic
    # parameters of total AL bioactivity in Group 2 (day 14)"), column
    # "Value". Monolix population parameters, reported on the natural
    # scale and log-transformed here because the paper states "The
    # distribution of pharmacokinetic parameters was normal in a
    # logarithmic scale" (Methods, pop-PK modelling).
    #
    # Table 2 labels the absorption parameter "Tk_pop" in the fixed-
    # effects block and "Omega_Tk0" in the random-effects block; the
    # parameter is Monolix's Tk0, the duration of the zero-order
    # absorption process. Absorption is zero-order "without delay"
    # (Results), so there is no lag time.
    # ==================================================================
    ld1 <- log(1.11);   label("Duration of zero-order absorption into central (h)")     # Table 2 Tk_pop (= Tk0) = 1.11 h, S.E. 0.12, RSE 11.2%
    lvc <- log(32.57);  label("Apparent central volume of distribution V/F (L)")        # Table 2 V_pop = 32.57 L, S.E. 2.9, RSE 8.89%
    lcl <- log(17.2);   label("Apparent clearance CL/F (L/h)")                          # Table 2 Cl_pop = 17.2 L/h, S.E. 2.03, RSE 11.8%

    # ==================================================================
    # Between-subject variability -- Saeheng 2024 Table 2, "Standard
    # deviation of Random effects". Monolix reports omega as the SD of
    # the random effect on the log scale, so the nlmixr2 variance is
    # omega^2 and the reported %CV column is sqrt(exp(omega^2) - 1)*100.
    #
    # The printed Omega_V value in Table 2 (8.89) cannot be a log-scale
    # SD: it would imply a CV of ~1e19%, not the 27.3% printed beside
    # it. 8.89 L is the SD of V on the LINEAR scale (8.89 / 32.57 =
    # 27.29%, matching the printed 27.3% CV; its printed S.E. 2.44
    # likewise gives 2.44 / 8.89 = 27.4%, matching the printed %RSE).
    # 8.89 is also the value printed as V_pop's %RSE two rows above,
    # consistent with a transcription slip. The Table 1 and Table 2
    # Omega_V value/S.E. cells are on the linear scale while every other
    # omega row -- including Omega_V in Table 3 -- is on the log scale.
    # The paper's Methods statement that all parameters are lognormal
    # governs, so the log-scale omega is recovered from the %CV column,
    # which is self-consistent across all three tables:
    # omega = sqrt(log(1 + CV^2)). See the vignette Errata.
    # ==================================================================
    etald1 ~ 0.122500  # Table 2 Omega_Tk0 = 0.35 (RSE 26.3%) -> variance 0.35^2; implied CV sqrt(exp(0.35^2)-1) = 36.1%, printed 36.57%
    etalvc ~ 0.071888  # Table 2 Omega_V printed on the linear scale (8.89 L); recovered from the printed %CV 27.3% as log(1 + 0.273^2) (omega 0.2681)
    etalcl ~ 0.193600  # Table 2 Omega_Cl = 0.44 (RSE 19.5%) -> variance 0.44^2; implied CV sqrt(exp(0.44^2)-1) = 46.2%, printed 46.6%

    # ==================================================================
    # Residual error -- Saeheng 2024 Table 2, "Error model parameters".
    # Monolix reports a combined error model as a (additive term, in the
    # concentration units of the data) and b (proportional coefficient).
    # The paper does not state whether the Monolix combined1
    # (sd = a + b*f) or combined2 (sd = sqrt(a^2 + (b*f)^2)) variance
    # form was used; rxode2's default combined2 form is used here, in
    # line with the existing library precedent (Chen_2024_febuxostat).
    # See the vignette Errata.
    # ==================================================================
    addSd  <- 1.39;  label("Additive residual SD (mg/L)")            # Table 2 a = 1.39, S.E. 0.33, RSE 24.0%
    propSd <- 0.17;  label("Proportional residual SD (fraction)")    # Table 2 b = 0.17, S.E. 0.087, RSE 19.3%
  })

  model({
    # 1. Individual parameters (lognormal, per Saeheng 2024 Methods).
    d1 <- exp(ld1 + etald1)
    vc <- exp(lvc + etalvc)
    cl <- exp(lcl + etalcl)

    # 2. Micro-constant.
    kel <- cl / vc

    # 3. ODE system -- one compartment with zero-order input directly
    #    into central. Dose events target cmt = "central" with rate = -2
    #    so that rxode2 uses the modelled dur(central) below.
    d/dt(central) <- -kel * central

    dur(central) <- d1

    # 4. Observation. The dose is in mg of CMC-AL, central is in mg and
    #    V/F is in L, so central / vc is in mg/L -- the units of the
    #    Cmax cut-offs Saeheng 2024 uses for dose selection (32.39 mg/L
    #    for tumour-progression inhibition, 21.42 mg/L for death
    #    prevention; Methods, "Criteria for optimal dosage regimens").
    Cc <- central / vc

    Cc ~ add(addSd) + prop(propSd)
  })
}
