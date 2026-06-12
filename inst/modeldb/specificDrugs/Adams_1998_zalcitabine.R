Adams_1998_zalcitabine <- function() {
  description <- "One-compartment first-order-absorption population PK model for oral zalcitabine (ddC) in HIV-infected adults (Adams 1998). Apparent clearance (CL/F = 14.8 L/h) and apparent volume of distribution (V/F = 87.6 L) were estimated from sparse-sampling clinic data; the absorption rate constant was not estimable in Adams 1998 (paper Discussion p. 412) and is fixed in this model to ka = 2.5 /h from primary single-dose PK data (Klecker 1988). No baseline covariates (age, sex, total body weight, calculated creatinine clearance, food, concomitant zidovudine) improved the basic fit and none were retained in the final model."
  reference <- paste(
    "Adams JM, Shelton MJ, Hewitt RG, DeRemer M, DiFrancesco R, Grasela TH, Morse GD.",
    "Zalcitabine population pharmacokinetics: application of radioimmunoassay.",
    "Antimicrob Agents Chemother. 1998 Feb;42(2):409-413.",
    "doi:10.1128/aac.42.2.409.",
    "The absorption rate constant ka was not estimable in Adams 1998",
    "(paper Discussion p. 412) and is fixed in this model to ka = 2.5 /h",
    "per primary single-dose ddC PK data from",
    "Klecker RW Jr, Collins JM, Yarchoan R, Thomas R, McAtee N, Broder S, Myers CE.",
    "Pharmacokinetics of 2',3'-dideoxycytidine in patients with AIDS and related disorders.",
    "J Clin Pharmacol. 1988 Sep;28(9):837-842. doi:10.1002/j.1552-4604.1988.tb03225.x.",
    "See model file ini() comments and validation-vignette Assumptions and deviations for the substitution rationale.",
    sep = " "
  )
  vignette <- "Adams_1998_zalcitabine"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list()

  covariatesDataExcluded <- list(
    WT = list(
      description = "Total body weight.",
      units       = "kg",
      type        = "continuous",
      notes       = "Reported in Adams 1998 Table 2 (79.1 +/- 15.0 kg; range 46.5-123) and screened on V/F and CL/F, but did not improve the basic model fit and was not retained (Results p. 412 col 2; Discussion p. 412 col 2). The paper attributes the null finding partly to limited heterogeneity in the cohort."
    ),
    AGE = list(
      description = "Subject age at enrolment.",
      units       = "years",
      type        = "continuous",
      notes       = "Reported in Adams 1998 Table 2 (38.6 +/- 7.13 years; range 27-56) and screened on CL/F, but not retained."
    ),
    SEXF = list(
      description = "Sex (1 = female, 0 = male).",
      units       = "(binary)",
      type        = "binary",
      notes       = "Cohort 39 male / 5 female (Adams 1998 Table 2). Screened on CL/F, not retained. Heavy male predominance limits inference."
    ),
    CRCL = list(
      description = "Calculated creatinine clearance (Cockcroft-Gault style, raw mL/min, NOT BSA-normalized).",
      units       = "mL/min",
      type        = "continuous",
      notes       = "Reported in Adams 1998 Table 2 (89.1 +/- 21.5 mL/min; range 53.6-146). Screened on CL/F, not retained. The paper notes that despite the 54-146 mL/min range, many patients clustered around a similar CLCR (Discussion p. 412 col 2). Zalcitabine is primarily renally excreted and the absence of a retained CLCR-on-CL/F effect is attributed to limited heterogeneity in this cohort, not to a biological argument that the drug is renal-clearance-independent."
    ),
    FOOD_CONCOMITANT = list(
      description = "Administration of zalcitabine with food (1 = with food, 0 = fasted) at the index dose, captured by clinic questionnaire.",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened on CL/F, not retained (Adams 1998 Methods p. 410; Results p. 412 col 2)."
    ),
    ZDV_CONCOMITANT = list(
      description = "Concomitant zidovudine (AZT) co-administration (1 = co-administered, 0 = not).",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened on CL/F, not retained (Adams 1998 Methods p. 410; Results p. 412 col 2). Investigated because zidovudine + zalcitabine was a commonly-used combination antiretroviral regimen at the time."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 44L,
    n_studies      = 1L,
    n_observations = 81L,
    age_range      = "27-56 years (mean 38.6 +/- 7.13 SD)",
    weight_range   = "46.5-123 kg (mean 79.1 +/- 15.0 SD)",
    crcl_range     = "53.6-146 mL/min (mean 89.1 +/- 21.5 SD; Cockcroft-Gault calculated)",
    sex_female_pct = 11.4,
    race_ethnicity = c(Caucasian = 73, AfricanAmerican = 9, Hispanic = 16, AmericanIndian = 2),
    disease_state  = "HIV-infected adults in routine clinic follow-up at the Erie County Medical Center Immunodeficiency Clinic (Buffalo, NY).",
    dose_range     = "0.375 mg or 0.75 mg oral every 8 hours (zalcitabine).",
    regions        = "United States (Buffalo, NY and Rochester, NY).",
    notes          = "Demographics from Adams 1998 Table 2. Sampling: 1.84 +/- 1.24 samples per patient (range 1-6), randomly timed during routine clinic visits; dosing-history and meal information obtained from clinic questionnaire. Risk-factor distribution: 64% homosexual activity, 16% intravenous drug use, 20% other or unknown. Observed zalcitabine concentrations ranged 2.01-8.57 ng/mL (assay LOQ = 2 ng/mL). HIV-infected adults; immunologic / virologic disease severity not stratified in the paper."
  )

  ini({
    # Structural parameters from Adams 1998 Table 3 ('Pharmacostatistical model and parameter estimates'; final model).
    # The paper reports CL/F and V/F (apparent oral parameters relative to bioavailability F).
    # Absorption rate constant ka was NOT estimable in Adams 1998 ('The absorption rate
    # constant could not be modeled because of the paucity of blood samples collected
    # early in a dosing interval.' -- Results p. 411 col 2; reiterated in Discussion
    # p. 412 col 1). Per operator guidance (2026-06-07 sidecar response), ka is fixed
    # in this model to a primary-PK literature value; see vignette Assumptions and
    # deviations for the substitution rationale.

    # Non-paper provenance: ka fixed at 2.5 /h from Klecker et al. 1988
    # (Klecker RW Jr, Collins JM, Yarchoan R, Thomas R, McAtee N, Broder S, Myers CE.
    #  Pharmacokinetics of 2',3'-dideoxycytidine in patients with AIDS and related
    #  disorders. J Clin Pharmacol 1988;28(9):837-842,
    #  doi:10.1002/j.1552-4604.1988.tb03225.x), an early dedicated single-dose PK
    #  study of oral zalcitabine in AIDS patients. Operator-named candidate primary
    #  source. The substitution is documented in the validation vignette.
    lka <- fixed(log(2.5)); label("First-order absorption rate constant (1/h, FIXED; non-paper provenance per Klecker 1988)")

    lcl <- log(14.8); label("Apparent clearance CL/F (L/h)")                     # Adams 1998 Table 3 / Results p. 412 col 2: CL/F = 14.8 L/h (0.19 L/h/kg; 95% CI 0.18-0.21 L/h/kg)
    lvc <- log(87.6); label("Apparent central volume of distribution V/F (L)")  # Adams 1998 Table 3 / Results p. 412 col 2: V/F = 87.6 L (1.18 L/kg; 95% CI 1.07-1.30 L/kg)

    # Inter-individual variability (Adams 1998 Table 3: CL = TVCL * exp(eta1);
    # V = TVV * exp(eta2); no off-diagonal covariance reported). The paper reports
    # the IIVs as %CV values; converted to log-scale variance via
    # omega^2 = log(1 + CV^2).
    etalcl ~ 0.0551  # %CV = 23.8 per Table 3 -> log(1 + 0.238^2) = 0.05511
    etalvc ~ 0.2559  # %CV = 54.0 per Table 3 -> log(1 + 0.540^2) = 0.25589

    # Residual error (Adams 1998 Table 3 'Error': y = F + F * eps1; i.e., proportional).
    # The paper reports 'Residual variability = 20.6%'.
    propSd <- 0.206; label("Proportional residual error (fraction)")  # Adams 1998 Table 3 / Results p. 412 col 2
  })

  model({
    # Individual PK parameters. No baseline covariates were retained in Adams 1998 final
    # model; no covariate scaling is applied.
    ka <- exp(lka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)

    # Micro-constant
    kel <- cl / vc

    # ODE system: one-compartment with first-order absorption from depot.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Observation in ng/mL (matching the radioimmunoassay calibration range reported
    # in Adams 1998 Table 1). Dose is in mg and V is in L, so central / vc has units
    # of mg/L = ug/mL; multiply by 1000 to express the simulation in ng/mL.
    Cc <- central / vc * 1000
    Cc ~ prop(propSd)
  })
}
