Dhondt_2017_mavacoxib_cockatiel_iv <- function() {
  description <- paste(
    "Preclinical (cockatiel).",
    "One-compartment population PK model with first-order elimination for",
    "mavacoxib after a single 4 mg/kg intravenous bolus of an analytical",
    "standard solution (STD) to cockatiels (Nymphicus hollandicus).",
    "Fitted in Phoenix NLME (FOCE-ELS); Table 2, 'IV STD' block. The",
    "combination of a large volume of distribution (10.99 L/kg) and a very",
    "low clearance (0.036 L/h/kg) gives a terminal half-life of about",
    "212 h, more than 150-fold longer than celecoxib in the same birds and",
    "the reason the authors propose less frequent dosing. Volumes and",
    "clearances are per kilogram, so the dosed amount is ug/kg and",
    "central/vc lands directly in ng/mL, the assay units of Figure 3a.",
    "Body weight and sex were screened as covariates and neither was",
    "retained. See Dhondt_2017_mavacoxib_cockatiel_oral_std and",
    "Dhondt_2017_mavacoxib_cockatiel_oral_cf for the separately fitted",
    "oral arms.",
    sep = " "
  )
  reference <- paste(
    "Dhondt L, Devreese M, Croubels S, De Baere S, Haesendonck R,",
    "Goessens T, Gehring R, De Backer P, Antonissen G.",
    "Comparative population pharmacokinetics and absolute oral",
    "bioavailability of COX-2 selective inhibitors celecoxib, mavacoxib and",
    "meloxicam in cockatiels (Nymphicus hollandicus).",
    "Sci Rep. 2017;7(1):12043.",
    "doi:10.1038/s41598-017-12159-z.",
    sep = " "
  )
  vignette <- "Dhondt_2017_cox2inhibitors_cockatiel"
  units <- list(
    time = "h",
    dosing = "ug",
    concentration = "ng/mL"
  )

  compartmentData <- list(
    central = list(analyte = "mavacoxib", units = "ug", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list()

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Screened as a continuous covariate by stepwise forward-backward",
        "selection against the -2LL criterion; not significant for any of",
        "the drugs and therefore not retained. Methods 'Pharmacokinetic",
        "analysis'; Results 'Pharmacokinetic analysis'."
      ),
      source_name = "BW"
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = male",
      notes = paste(
        "The source screened 'gender' as a categorical covariate and did",
        "not retain it; the source does not state which sex was the",
        "reference category. This cohort was balanced 20 male / 20 female."
      ),
      source_name = "gender"
    )
  )

  population <- list(
    species = "cockatiel (Nymphicus hollandicus)",
    n_subjects = 40L,
    n_studies = 1L,
    age_range = "6-12 months",
    weight_median = "101 g",
    weight_range = "101 +/- 8 g (mean +/- SD)",
    sex_female_pct = 50,
    disease_state = "Healthy (no disease model)",
    dose_range = paste(
      "Single 4 mg/kg body weight intravenous bolus into the vena cutanea",
      "ulnaris (wing vein) of a mavacoxib analytical standard solution",
      "(5 mg/mL in polyethylene glycol 400:physiological saline, 75:25 v/v)."
    ),
    regions = "Belgium (Ghent University, Merelbeke)",
    notes = paste(
      "The 40 birds (20 male / 20 female) received mavacoxib STD both",
      "intravenously and orally in a two-way crossover with a three-month",
      "washout; the two arms were fitted as separate models and the oral",
      "arm is Dhondt_2017_mavacoxib_cockatiel_oral_std. A sparse sampling",
      "protocol was used because of the limited blood volume of cockatiels:",
      "sampling times were randomly allocated across birds with a maximum",
      "of three samples per bird, drawn before dosing and at 5, 15, 30 and",
      "45 min and 1, 2, 4, 6, 8, 12, 24, 48, 72, 96, 120, 168, 336, 672 and",
      "1008 h -- a 42-day window needed to characterise the long terminal",
      "phase. Mavacoxib was quantified in plasma by LC-MS/MS over",
      "5-5000 ng/mL (LOQ 5 ng/mL, LOD 0.25 ng/mL); values below the LOQ",
      "were excluded before fitting. Plasma protein binding measured",
      "separately was 97.02 +/- 0.32% and is not a model parameter. See",
      "Dhondt 2017 Methods 'Animals and experimental procedure' and",
      "'Mavacoxib PK study'."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Structural parameters -- Dhondt 2017 Table 2, 'IV STD' block.
    # Methods equation (2): C(t) = C0 * exp(-Cl/Vd * t), a one-compartment
    # model with first-order elimination; Methods 'Population
    # pharmacokinetics mavacoxib' states the mavacoxib IV structural model
    # is the same one-compartment first-order form used for celecoxib.
    #
    # Both values are confirmed by the table's own computed secondary
    # parameters: C0 = D/Vd = 4000/10.99 = 364 ng/mL (Table 2 reports
    # 363.75), AUC(0-inf) = D/Cl = 4000/0.036 = 111111 ng.h/mL (Table 2
    # reports 111238), Ke = Cl/Vd = 0.00328 /h (Table 2 reports 0.0033) and
    # T1/2el = 211 h (Table 2 reports 211.97).
    # ------------------------------------------------------------------
    lvc <- log(10.99); label("Volume of distribution Vd (log L/kg)") # Table 2, IV STD: Vd = 10.99 L/kg (RSE 5.51%)
    lcl <- log(0.036); label("Total body clearance Cl (log L/h/kg)") # Table 2, IV STD: Cl = 0.036 L/h.kg (RSE 12.16%)

    # ------------------------------------------------------------------
    # IIV -- Dhondt 2017 Table 2, 'omega' column, read as the SD of eta on
    # the log scale and squared here because nlmixr2 omega entries are
    # variances. Methods equation (4): P_i = theta_P * exp(eta_Pi), with
    # variance omega^2 and 'Interindividual variability is reported as
    # omega'. See Dhondt_2017_celecoxib_cockatiel_iv and the vignette
    # 'Assumptions and deviations' section for why the equations are
    # preferred over the caption's word 'variance'.
    #   Vd: 0.075 -> 0.005625   (about 7.5% CV)
    #   Cl: 0.524 -> 0.274576   (about 56% CV)
    # ------------------------------------------------------------------
    etalvc ~ 0.005625 # Table 2, IV STD omega for Vd = 0.075
    etalcl ~ 0.274576 # Table 2, IV STD omega for Cl = 0.524

    # ------------------------------------------------------------------
    # Residual error -- Methods 'Population pharmacokinetics mavacoxib':
    # 'Inter- and intra-individual variability were expressed according to
    # the exponential error model and the multiplicative residual error
    # model, respectively.' Equation (5) is C_obs = C_pred * (1 + epsilon),
    # which is nlmixr2's prop().
    # ------------------------------------------------------------------
    propSd <- 0.26; label("Proportional residual error (fraction)") # Table 2, IV STD: Res. Error = 0.26 (RSE 13.74%)
  })

  model({
    # 1. Individual parameters (exponential IIV, Methods equation (4)).
    vc <- exp(lvc + etalvc)
    cl <- exp(lcl + etalcl)

    # 2. Micro-constant.
    kel <- cl / vc

    # 3. ODE system -- Methods equation (2) written as a differential
    # equation. The intravenous dose enters central directly.
    d/dt(central) <- -kel * central

    # 4. Observation. Amount in ug/kg over a volume in L/kg gives ug/L,
    # which is ng/mL -- the assay units of Figure 3a.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
