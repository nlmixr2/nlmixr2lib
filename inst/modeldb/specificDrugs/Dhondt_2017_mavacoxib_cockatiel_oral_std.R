Dhondt_2017_mavacoxib_cockatiel_oral_std <- function() {
  description <- paste(
    "Preclinical (cockatiel).",
    "One-compartment population PK model with first-order absorption and",
    "first-order elimination for mavacoxib after a single 4 mg/kg oral",
    "intra-crop bolus of an analytical standard solution (STD) to",
    "cockatiels (Nymphicus hollandicus). Fitted in Phoenix NLME",
    "(FOCE-ELS); Table 2, 'PO STD' block. No absorption lag time was",
    "retained for any mavacoxib arm. The disposition parameters are",
    "apparent (Vd/F, Cl/F): no bioavailability term was estimated, so the",
    "full 4 mg/kg dose enters the depot. Absolute oral bioavailability was",
    "computed separately from the AUC ratio against the intravenous arm as",
    "F = 113%, i.e. essentially complete absorption. Slow absorption",
    "(Ka = 0.20 /h) against a 172 h terminal half-life gives a Tmax of",
    "about 20 h. Volumes and clearances are per kilogram, so the dosed",
    "amount is ug/kg and central/vc lands directly in ng/mL, the assay",
    "units of Figure 3b. Body weight and sex were screened as covariates",
    "and neither was retained. See Dhondt_2017_mavacoxib_cockatiel_iv and",
    "Dhondt_2017_mavacoxib_cockatiel_oral_cf for the other two arms.",
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
    depot = list(analyte = "mavacoxib", units = "ug", specimen = "administration site", verified = TRUE),
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
      "Single 4 mg/kg body weight oral intra-crop bolus of a mavacoxib",
      "analytical standard solution (5 mg/mL in polyethylene glycol",
      "400:physiological saline, 75:25 v/v), given by a curved ball-tipped",
      "feeding needle."
    ),
    regions = "Belgium (Ghent University, Merelbeke)",
    notes = paste(
      "The same 40 birds (20 male / 20 female) received mavacoxib STD both",
      "orally and intravenously in a two-way crossover with a three-month",
      "washout; the two arms were fitted as separate models and the",
      "intravenous arm is Dhondt_2017_mavacoxib_cockatiel_iv. A sparse",
      "sampling protocol was used because of the limited blood volume of",
      "cockatiels: sampling times were randomly allocated across birds with",
      "a maximum of three samples per bird, drawn before dosing and at 15,",
      "30 and 45 min and 1, 2, 4, 6, 8, 12, 24, 48, 72, 96, 120, 168, 336,",
      "672 and 1008 h. Mavacoxib was quantified in plasma by LC-MS/MS over",
      "5-5000 ng/mL (LOQ 5 ng/mL, LOD 0.25 ng/mL); values below the LOQ",
      "were excluded before fitting. Feed was withheld from 8 h before",
      "until 4 h after dosing, but each bird received a 2 mL intra-crop",
      "feed bolus immediately after dosing. Plasma protein binding measured",
      "separately was 97.02 +/- 0.32% and is not a model parameter. See",
      "Dhondt 2017 Methods 'Animals and experimental procedure' and",
      "'Mavacoxib PK study'."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Structural parameters -- Dhondt 2017 Table 2, 'PO STD' block.
    # Methods equation (3):
    #   C(t) = F*D*ka / (Vd*(ka - ke)) * (exp(-ke*t) - exp(-ka*t))
    # Table 2 reports Vd/F and Cl/F -- apparent values uncorrected for
    # bioavailability -- so F*D/Vd is D/(Vd/F) and the full dose is used.
    # No Tlag row appears for any mavacoxib arm.
    #
    # All three values are confirmed by the table's own computed secondary
    # parameters: Ke = (Cl/F)/(Vd/F) = 0.031/7.85 = 0.00395 /h (Table 2
    # reports 0.0040), Tmax = ln(ka/ke)/(ka - ke) = 20.0 h (Table 2 reports
    # 19.98), Cmax = 471 ng/mL (Table 2 reports 469.69) and AUC(0-inf) =
    # D/(Cl/F) = 4000/0.031 = 129032 ng.h/mL (Table 2 reports 126108, which
    # matches an unrounded Cl/F of 0.03172).
    # ------------------------------------------------------------------
    lka <- log(0.20); label("First-order absorption rate constant Ka (log 1/h)") # Table 2, PO STD: Ka = 0.20 /h (RSE 14.80%)
    lvc <- log(7.85); label("Apparent volume of distribution Vd/F (log L/kg)") # Table 2, PO STD: Vd/F = 7.85 L/kg (RSE 4.26%)
    lcl <- log(0.031); label("Apparent total body clearance Cl/F (log L/h/kg)") # Table 2, PO STD: Cl/F = 0.031 L/h.kg (RSE 10.80%)

    # ------------------------------------------------------------------
    # IIV -- Dhondt 2017 Table 2, 'omega' column, read as the SD of eta on
    # the log scale and squared here because nlmixr2 omega entries are
    # variances. See Dhondt_2017_celecoxib_cockatiel_iv and the vignette
    # 'Assumptions and deviations' section for why the equations are
    # preferred over the caption's word 'variance'.
    #   Ka:   0.171 -> 0.029241  (about 17% CV)
    #   Vd/F: 0.006 -> 3.6e-05   (about 0.6% CV)
    #   Cl/F: 0.260 -> 0.0676    (about 26% CV)
    # ------------------------------------------------------------------
    etalka ~ 0.029241 # Table 2, PO STD omega for Ka = 0.171
    etalvc ~ 3.6e-05 # Table 2, PO STD omega for Vd/F = 0.006
    etalcl ~ 0.0676 # Table 2, PO STD omega for Cl/F = 0.260

    # ------------------------------------------------------------------
    # Residual error -- Methods 'Population pharmacokinetics mavacoxib':
    # 'Inter- and intra-individual variability were expressed according to
    # the exponential error model and the multiplicative residual error
    # model, respectively.' Equation (5) is C_obs = C_pred * (1 + epsilon),
    # which is nlmixr2's prop().
    # ------------------------------------------------------------------
    propSd <- 0.29; label("Proportional residual error (fraction)") # Table 2, PO STD: Res. Error = 0.29 (RSE 14.40%)
  })

  model({
    # 1. Individual parameters (exponential IIV, Methods equation (4)).
    ka <- exp(lka + etalka)
    vc <- exp(lvc + etalvc)
    cl <- exp(lcl + etalcl)

    # 2. Micro-constant.
    kel <- cl / vc

    # 3. ODE system -- Methods equation (3) written as differential
    # equations. The oral dose enters the depot in full; bioavailability
    # is absorbed into the apparent Vd/F and Cl/F.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # 4. Observation. Amount in ug/kg over an apparent volume in L/kg
    # gives ug/L, which is ng/mL -- the assay units of Figure 3b.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
