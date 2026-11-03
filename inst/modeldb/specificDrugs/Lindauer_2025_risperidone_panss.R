Lindauer_2025_risperidone_panss <- function() {
  description <- paste(
    "Exposure-response model for PANSS total score in adults with an acute",
    "exacerbation of schizophrenia treated with monthly intramuscular",
    "Risperidone ISM (75 or 100 mg) or placebo in the Phase III PRISMA-3",
    "trial (Lindauer 2025). Weibull placebo time course, an Emax drug effect",
    "driven by the individual average steady-state active-moiety",
    "(risperidone + 9-hydroxyrisperidone) concentration from the upstream",
    "Risperidone ISM population PK model with a first-order onset delay, and",
    "a linear drift term in actively treated patients. Includes the",
    "constant-baseline-hazard dropout model the authors used in simulation,",
    "driven by the predicted PANSS and the enrollment country."
  )
  reference <- paste(
    "Lindauer A, Snoeck E, Laveille C, Ayani I, Ochoa Diaz de Monasterioguren L,",
    "Almendros M, Martinez-Gonzalez J, Anta L, Gutierro I. Exposure-Efficacy",
    "Analysis and Dopamine D2 Receptor Occupancy in Adults with Schizophrenia",
    "after Treatment with the Monthly Intramuscular Injectable Risperidone ISM.",
    "J Clin Pharmacol. 2025;65(3):350-360. doi:10.1002/jcph.6152"
  )
  vignette <- "Lindauer_2025_risperidone"
  units <- list(time = "day", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    CAV = list(
      description = "Individual average steady-state active-moiety plasma concentration (Cav,ss) of risperidone + 9-hydroxyrisperidone",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Computed per subject and per PANSS / PK sampling occasion as",
        "Cav,ss = Dose / (CL/F * dosing interval), with CL/F the individual",
        "post hoc apparent clearance from the Risperidone ISM population PK",
        "model of Laveille 2024 (Br J Clin Pharmacol 90:2256; doi:10.1111/bcp.16115)",
        "and a 28-day dosing interval (Methods, 'Datafile Creation'). Set to 0",
        "for placebo subjects. The $INPUT column that the supplemental control",
        "stream feeds to the Emax term is named CONC ('Predicted concentration",
        "of the active moiety from the pop PK model'); the Methods define the",
        "driver as Cav,ss. Also gates the DRIFT term: drift applies only when",
        "CAV > 0 (the control stream's IF(DOSE.NE.0) TVDRIFT = THETA(8))."
      ),
      source_name = "CONC (Cav,ss)"
    ),
    PANSS_BL = list(
      description = "Observed baseline PANSS total score",
      units = "PANSS points (30-210)",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-fixed per subject. Enters as a centred linear effect on the",
        "placebo time constant TD, TD = TVTD * (1 + 0.0566 * (PANSS_BL - 95))",
        "(control stream TDBLPANSS block; Table 2 'BL~TD'). Distinct from the",
        "model's estimated baseline parameter BL, which carries IIV."
      ),
      source_name = "BLPANSS"
    ),
    REGION_USA = list(
      description = "Enrolled at a United States study site (1) vs Ukraine (0)",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (United States) for the baseline PANSS effect; the PRISMA-3 cohort was enrolled only in the USA and Ukraine",
      notes = paste(
        "PRISMA-3 enrolled in the USA (COUNTRY = 840, 56.4%) and Ukraine",
        "(COUNTRY = 804, 43.6%) only (Table 1), so the Ukraine indicator is",
        "1 - REGION_USA. Ukraine raises the typical baseline PANSS by 7.27%",
        "(Table 2 'COUNTRY~BL'). REGION_USA = 1 also raises the dropout",
        "hazard (Appendix S1 'beta_USA')."
      ),
      source_name = "COUNTRY"
    ),
    INJSITE_ARM = list(
      description = "Intramuscular injection into the deltoid (1) vs the gluteus (0)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (gluteal intramuscular injection)",
      notes = paste(
        "Source column SINJ (1 = gluteus, 2 = deltoid), so INJSITE_ARM =",
        "SINJ - 1. A deltoid injection shortens the half-life of the",
        "drug-effect onset KTHL by 69.1% (Table 2 'Inj.Site~KTHL'). Placebo",
        "subjects also carry an injection site; it has no effect because the",
        "drug effect is zero when CAV = 0."
      ),
      source_name = "SINJ"
    )
  )

  covariatesDataExcluded <- list(
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Tested in the forward-inclusion / backward-deletion covariate search on every parameter with IIV; not retained (Methods, Results)."
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Tested in the stepwise covariate search; not retained."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      notes = "Tested in the stepwise covariate search (source SEX, 1 = male, 2 = female); not retained."
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Tested in the stepwise covariate search; not retained."
    ),
    RACE_BLACK = list(
      description = "Black or African-American race indicator",
      units = "(binary)",
      type = "binary",
      notes = "Race (source RACE, 1 = Caucasian, 2 = African-American, 3 = Asian, 4 = Other) was tested in the stepwise covariate search; not retained."
    )
  )

  compartmentData <- list(
    cumhaz = list(
      analyte = "dropout cumulative hazard",
      units = "(unitless)",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 383L,
    n_studies = 1L,
    age_range = "18-64 years",
    age_median = "41 years",
    weight_range = "46.8-135.0 kg",
    weight_median = "81.7 kg",
    bmi_median = "27.7 kg/m^2 (range 18.4-40.1)",
    sex_female_pct = 34.1,
    disease_state = "Adults with schizophrenia in an acute exacerbation (baseline PANSS 80-120, median 95)",
    dose_range = "Risperidone ISM 75 or 100 mg, or placebo, intramuscular (gluteal or deltoid) on days 1, 29 and 57",
    regions = "USA (56.4%), Ukraine (43.6%)",
    notes = paste(
      "PRISMA-3 (NCT03160521) 12-week double-blind phase. The PANSS",
      "analysis dataset had 383 patients (131 placebo, 127 at 75 mg, 125 at",
      "100 mg) and 2403 PANSS measurements; Table 1 summarises the 390",
      "patients of the summary population. PANSS was measured on days 1, 4,",
      "8, 15, 29, 57 and 85. NONMEM 7.3, FOCE."
    )
  )

  ini({
    # Placebo sub-model; time in days
    lbl <- log(92.9); label("Typical baseline PANSS total score in US patients (PANSS points)") # Table 2 'BL' = 92.9
    ltd <- log(8.32); label("Time constant of the Weibull placebo response TD at PANSS_BL = 95 (day)") # Table 2 'TD' = 8.32
    pmax <- 0.117; label("Maximal placebo effect as a fractional reduction from baseline (fraction)") # Table 2 'PMAX' = 0.117
    pow <- 1.39; label("Weibull shape parameter of the placebo response (unitless)") # Table 2 'POW' = 1.39

    # Drug sub-model
    emax <- 0.0661; label("Maximal drug effect as an additional fractional reduction from baseline (fraction)") # Table 2 'EMAX' = 0.0661
    lec50 <- fixed(log(4.6)); label("Active-moiety Cav,ss giving half the maximal drug effect (ng/mL)") # Table 2 'EC50' = 4.6 (fixed to Ivaturi 2017)
    drift <- -0.0929; label("Linear PANSS drift in actively treated patients (PANSS points/day)") # Table 2 'DRIFT' = -0.0929
    lkthl <- log(8.96); label("Half-life of the drug-effect onset delay KTHL, gluteal injection (day)") # Table 2 'KTHL' = 8.96

    # Covariate effects
    e_ukr_bl <- 0.0727; label("Fractional increase in baseline PANSS for patients enrolled in Ukraine vs the USA (fraction)") # Table 2 'COUNTRY~BL' = 7.27%
    e_injsite_arm_kthl <- -0.691; label("Fractional change in KTHL for a deltoid vs gluteal injection (fraction)") # Table 2 'Inj.Site~KTHL' = -69.1%
    e_panss_bl_td <- 0.0566; label("Fractional change in TD per PANSS point of observed baseline above 95 (1/PANSS point)") # Table 2 'BL~TD' = 5.66%

    # IIV: exponential on TD, BL and KTHL (variance = log(1 + CV^2));
    # additive on PMAX, EMAX and DRIFT (variance = SD^2). Covariances are
    # r_ij * sqrt(var_i * var_j) from the Table 2 correlation coefficients.
    # Order TD, BL, PMAX, EMAX, DRIFT, KTHL as in the control stream BLOCK(6).
    etaltd + etalbl + etapmax + etaemax + etadrift + etalkthl ~ c(
      0.292587,
      -0.00928652, 0.00422905,
      -0.0190834, 0.00113935, 0.0144,
      0.0457071, -0.001132, -0.00563784, 0.028561,
      -0.0073932, -0.00109799, 0.000850632, 0.00878665, 0.017956,
      -0.413039, 0.000630476, 0.0503031, -0.0959665, -0.0102569, 0.852541
    ) # Table 2: IIV_TD 58.3% CV, IIV_BL 6.51% CV, IIV_PMAX SD 0.12, IIV_EMAX SD 0.169, IIV_DRIFT SD 0.134, IIV_KTHL 116% CV, and the 15 COR_ rows

    addSd <- 4.78; label("Additive residual error on PANSS total score (PANSS points)") # Table 2 'addRUV' = 4.78

    # Dropout model (simulation only): h = BHAZ * exp(beta_PANSS * PANSS + beta_USA * USA)
    llambda <- log(4.572e-5); label("Baseline dropout hazard BHAZ (1/day)") # Appendix S1 table 'BHAZ (10e-5/day)' = 4.572
    e_panss_haz <- 0.049; label("Log-hazard coefficient of dropout per PANSS point (1/PANSS point)") # Appendix S1 table 'beta_PANSS' = 0.049
    e_region_usa_haz <- 0.996; label("Log-hazard coefficient of dropout for US enrollment (unitless)") # Appendix S1 table 'beta_USA' = 0.996
  })

  model({
    # Typical values with covariates
    ukr <- 1 - REGION_USA
    active <- CAV > 0

    bl <- exp(lbl + etalbl) * (1 + e_ukr_bl * ukr)
    td <- exp(ltd + etaltd) * (1 + e_panss_bl_td * (PANSS_BL - 95))
    pmax_i <- pmax + etapmax
    emax_i <- emax + etaemax
    drift_i <- drift * active + etadrift
    kthl <- exp(lkthl + etalkthl) * (1 + e_injsite_arm_kthl * INJSITE_ARM)
    ec50 <- exp(lec50)
    kt <- log(2) / kthl

    # Placebo (Weibull) and drug (Emax with first-order onset) effects; t is
    # time since the first injection in days
    placebo_eff <- pmax_i * (1 - exp(-(t / td)^pow))
    drug_eff <- emax_i * CAV / (ec50 + CAV) * (1 - exp(-kt * t))

    panss <- bl * (1 - placebo_eff - drug_eff) + drift_i * t

    # Dropout hazard on the predicted PANSS (constant baseline hazard)
    hazard <- exp(llambda + e_panss_haz * panss + e_region_usa_haz * REGION_USA)
    d/dt(cumhaz) <- hazard
    sur <- exp(-cumhaz)

    panss ~ add(addSd)
  })
}
