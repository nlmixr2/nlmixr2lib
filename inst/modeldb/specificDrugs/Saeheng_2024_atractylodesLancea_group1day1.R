Saeheng_2024_atractylodesLancea_group1day1 <- function() {
  description <- paste0(
    "One-compartment population PK model with zero-order absorption (no lag) and linear ",
    "clearance for total Atractylodes lancea (Thunb.) DC. bioactivity in patients with ",
    "advanced-stage intrahepatic cholangiocarcinoma (Saeheng 2024, Table 1). This file holds ",
    "the Group 1 day-1 fit: the first once-daily 1,000 mg dose of the capsule formulation of ",
    "standardized AL extract (CMC-AL). Saeheng 2024 reported three separate fits of the same ",
    "structural model to three different occasions -- Group 1 day 1, Group 2 day 14 and Group 2 ",
    "day 28 -- carried here as three model files sharing one vignette; see also ",
    "modellib('Saeheng_2024_atractylodesLancea_group2day14') and ",
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
    n_subjects    = 15L,
    n_studies     = 1L,
    disease_state = "advanced-stage intrahepatic cholangiocarcinoma (iCCA)",
    dose_range    = paste0(
      "1,000 mg once daily (9 capsules) of the capsule formulation of standardized Atractylodes ",
      "lancea extract (CMC-AL) for 90 days, plus standard supportive care; this model is the ",
      "day-1 (first-dose) fit"
    ),
    regions = "Thailand (Sakhon Na-Kon Hospital, Sakhon Na-Kon)",
    notes   = paste0(
      "Group 1 arm of a single-centre, open-label, randomized, controlled phase 2A trial ",
      "(TCTR20210129007). Fifteen patients were randomized to Group 1 (Saeheng 2024 Results). ",
      "Plasma sampling on day 1 at 0, 0.25, 0.5, 1, 1.5, 2, 2.5, 3, 4, 5, 6 and 8 h post-dose ",
      "(Saeheng 2024 Methods, 'Pharmacokinetic study'). Baseline demographics (age, weight, sex, ",
      "height) are not tabulated in this paper; they are reported in the parent phase 2A trial ",
      "publication (reference [8]). Fitted in MonolixSuite 2023R1 by SAEM."
    )
  )

  ini({
    # ==================================================================
    # Structural parameters -- Saeheng 2024 Table 1 ("Pharmacokinetic
    # parameters of total AL bioactivity in Group 1 (day 1)"), column
    # "Value". Monolix population parameters, reported on the natural
    # scale and log-transformed here because the paper states "The
    # distribution of pharmacokinetic parameters was normal in a
    # logarithmic scale" (Methods, pop-PK modelling).
    #
    # Table 1 labels the absorption parameter "Tk_pop" in the fixed-
    # effects block and "Omega_Tk0" in the random-effects block; the
    # parameter is Monolix's Tk0, the duration of the zero-order
    # absorption process. Absorption is zero-order "without delay"
    # (Results), so there is no lag time.
    # ==================================================================
    ld1 <- log(0.85);   label("Duration of zero-order absorption into central (h)")     # Table 1 Tk_pop (= Tk0) = 0.85 h, S.E. 0.15, RSE 17.3%
    lvc <- log(42.82);  label("Apparent central volume of distribution V/F (L)")        # Table 1 V_pop = 42.82 L, S.E. 4.86, RSE 11.3%
    lcl <- log(20.9);   label("Apparent clearance CL/F (L/h)")                          # Table 1 Cl_pop = 20.9 L/h, S.E. 3.36, RSE 16.1%

    # ==================================================================
    # Between-subject variability -- Saeheng 2024 Table 1, "Standard
    # deviation of Random effects". Monolix reports omega as the SD of
    # the random effect on the log scale, so the nlmixr2 variance is
    # omega^2 and the reported %CV column is sqrt(exp(omega^2) - 1)*100.
    #
    # The printed Omega_V value in Table 1 (13.2) cannot be a log-scale
    # SD: it would imply a CV of ~1e40%, not the 30.82% printed beside
    # it. 13.2 L is the SD of V on the LINEAR scale (13.2 / 42.82 =
    # 30.83%, matching the printed 30.82% CV; its printed S.E. 4.25
    # likewise gives 4.25 / 13.2 = 32.2%, matching the printed %RSE).
    # The Table 1 and Table 2 Omega_V value/S.E. cells are therefore on
    # the linear scale while every other omega row -- including Omega_V
    # in Table 3 -- is on the log scale. The paper's Methods statement
    # that all parameters are lognormal governs, so the log-scale omega
    # is recovered from the %CV column, which is self-consistent across
    # all three tables: omega = sqrt(log(1 + CV^2)). See the vignette
    # Errata.
    # ==================================================================
    etald1 ~ 0.260100  # Table 1 Omega_Tk0 = 0.51 (RSE 27.2%) -> variance 0.51^2; implied CV sqrt(exp(0.51^2)-1) = 54.5%, printed 54.88%
    etalvc ~ 0.090736  # Table 1 Omega_V printed on the linear scale (13.2 L); recovered from the printed %CV 30.82% as log(1 + 0.3082^2) (omega 0.3012)
    etalcl ~ 0.291600  # Table 1 Omega_Cl = 0.54 (RSE 21.8%) -> variance 0.54^2; implied CV sqrt(exp(0.54^2)-1) = 58.2%, printed 57.91%

    # ==================================================================
    # Residual error -- Saeheng 2024 Table 1, "Error model parameters".
    # Monolix reports a combined error model as a (additive term, in the
    # concentration units of the data) and b (proportional coefficient).
    # The paper does not state whether the Monolix combined1
    # (sd = a + b*f) or combined2 (sd = sqrt(a^2 + (b*f)^2)) variance
    # form was used; rxode2's default combined2 form is used here, in
    # line with the existing library precedent (Chen_2024_febuxostat).
    # See the vignette Errata.
    # ==================================================================
    addSd  <- 0.36;  label("Additive residual SD (mg/L)")            # Table 1 a = 0.36, S.E. 0.1, RSE 28.4%
    propSd <- 0.32;  label("Proportional residual SD (fraction)")    # Table 1 b = 0.32, S.E. 0.038, RSE 11.8%
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
