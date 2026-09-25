Dhondt_2017_celecoxib_cockatiel_oral_cf <- function() {
  description <- paste(
    "Preclinical (cockatiel).",
    "One-compartment population PK model with a lagged first-order",
    "absorption and first-order elimination for celecoxib after a single",
    "10 mg/kg oral intra-crop bolus of the ground commercial tablet",
    "formulation (Celebrex, CF) to cockatiels (Nymphicus hollandicus).",
    "Fitted in Phoenix NLME (FOCE-ELS); Table 1, 'PO CF' block. The",
    "disposition parameters are apparent (Vd/F, Cl/F): no bioavailability",
    "term was estimated, so the full 10 mg/kg dose enters the depot.",
    "Absolute oral bioavailability was computed separately from the AUC",
    "ratio against the intravenous arm as F = 56%, notably lower than the",
    "110% of the standard solution. Volumes and clearances are per",
    "kilogram, so the dosed amount is ug/kg and central/vc lands directly",
    "in ng/mL, the assay units of Figure 2c. This is the only arm of the",
    "paper whose residual error is the combined additive-plus-proportional",
    "form. Body weight and sex were screened as covariates and neither was",
    "retained. See Dhondt_2017_celecoxib_cockatiel_iv and",
    "Dhondt_2017_celecoxib_cockatiel_oral_std for the other two arms.",
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
    depot = list(analyte = "celecoxib", units = "ug", specimen = "administration site", verified = TRUE),
    central = list(analyte = "celecoxib", units = "ug", specimen = "plasma", verified = TRUE)
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
        "reference category. This cohort was balanced 11 male / 11 female."
      ),
      source_name = "gender"
    )
  )

  population <- list(
    species = "cockatiel (Nymphicus hollandicus)",
    n_subjects = 22L,
    n_studies = 1L,
    age_range = "6-12 months",
    weight_median = "93 g",
    weight_range = "93 +/- 10 g (mean +/- SD)",
    sex_female_pct = 50,
    disease_state = "Healthy (no disease model)",
    dose_range = paste(
      "Single 10 mg/kg body weight oral intra-crop bolus of the commercial",
      "tablet formulation (Celebrex 100 mg, Pfizer), ground and blended",
      "with lactose to 10 mg celecoxib/g and suspended in physiological",
      "saline at 5 mg/mL immediately before dosing."
    ),
    regions = "Belgium (Ghent University, Merelbeke)",
    notes = paste(
      "A separate group of 22 birds (11 male / 11 female) from the 34 that",
      "received the standard solution. A sparse sampling protocol was used",
      "because of the limited blood volume of cockatiels: sampling times",
      "were randomly allocated across birds with a maximum of two samples",
      "per bird, drawn before dosing and at 15, 30 and 45 min and 1, 2, 4,",
      "6, 8, 12 and 24 h. Celecoxib was quantified in plasma by LC-MS/MS",
      "over 5-5000 ng/mL (LOQ 5 ng/mL, LOD 0.22 ng/mL); values below the",
      "LOQ were excluded before fitting. Feed was withheld from 8 h before",
      "until 4 h after dosing, but each bird received a 2 mL intra-crop",
      "feed bolus immediately after dosing. Plasma protein binding measured",
      "separately was 98.98 +/- 0.07% and is not a model parameter. See",
      "Dhondt 2017 Methods 'Animals and experimental procedure' and",
      "'Celecoxib PK study'."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Structural parameters -- Dhondt 2017 Table 1, 'PO CF' block.
    # Methods equation (3):
    #   C(t) = F*D*ka / (Vd*(ka - ke)) * (exp(-ke*t) - exp(-ka*t))
    # with an absorption lag time, which Results reports was retained for
    # the oral commercial formulations of celecoxib and meloxicam because
    # it significantly improved the fit.
    #
    # All four values are confirmed by the table's own computed secondary
    # parameters: Ke = (Cl/F)/(Vd/F) = 4.32/5.49 = 0.787 /h (Table 1
    # reports 0.79), T1/2el = 0.881 h (Table 1 reports 0.88),
    # Tmax = Tlag + ln(ka/ke)/(ka - ke) = 0.33 + 1.77 = 2.10 h (Table 1
    # reports 2.09), Cmax = 453 ng/mL (Table 1 reports 454.97) and
    # AUC(0-inf) = D/(Cl/F) = 10000/4.32 = 2315 ng.h/mL (Table 1 reports
    # 2312.24).
    # ------------------------------------------------------------------
    lka <- log(0.39); label("First-order absorption rate constant Ka (log 1/h)") # Table 1, PO CF: Ka = 0.39 /h (RSE 17.37%)
    ltlag <- log(0.33); label("Absorption lag time Tlag (log h)") # Table 1, PO CF: Tlag = 0.33 h (RSE 39.55%)
    lvc <- log(5.49); label("Apparent volume of distribution Vd/F (log L/kg)") # Table 1, PO CF: Vd/F = 5.49 L/kg (RSE 25.02%)
    lcl <- log(4.32); label("Apparent total body clearance Cl/F (log L/h/kg)") # Table 1, PO CF: Cl/F = 4.32 L/h.kg (RSE 11.67%)

    # ------------------------------------------------------------------
    # IIV -- Dhondt 2017 Table 1, 'omega' column, read as the SD of eta on
    # the log scale and squared here because nlmixr2 omega entries are
    # variances. See Dhondt_2017_celecoxib_cockatiel_iv and the vignette
    # 'Assumptions and deviations' section for why the equations are
    # preferred over the caption's word 'variance'.
    #   Ka:   < 0.001 -> encoded at the printed upper bound, 1e-06
    #   Tlag: 0.207   -> 0.042849  (about 21% CV)
    #   Vd/F: < 0.001 -> encoded at the printed upper bound, 1e-06
    #   Cl/F: 0.055   -> 0.003025  (about 5.5% CV)
    # ------------------------------------------------------------------
    etalka ~ 1e-06 # Table 1, PO CF omega for Ka, reported as '< 0.001'
    etaltlag ~ 0.042849 # Table 1, PO CF omega for Tlag = 0.207
    etalvc ~ 1e-06 # Table 1, PO CF omega for Vd/F, reported as '< 0.001'
    etalcl ~ 0.003025 # Table 1, PO CF omega for Cl/F = 0.055

    # ------------------------------------------------------------------
    # Residual error -- Methods: 'for the PO CF administration, the best
    # residual error model was log-additive (equation (7))'. Equation (7)
    # as printed is
    #   C_obs = C_pred + epsilon * sqrt(1 + C_pred^(2*(sigma_mult/sigma_add)^2))
    # The exponent is a typesetting error for a squared factor; the
    # intended form is Phoenix NLME's mix-ratio observation model
    #   C_obs = C_pred + epsilon * sqrt(1 + C_pred^2 * (sigma_mult/sigma_add)^2)
    # whose residual SD is sqrt(sigma_add^2 + (sigma_mult*C_pred)^2). As
    # printed the exponent would raise a concentration to a dimensionless
    # power, which is not a usable error model. The intended form is
    # exactly nlmixr2's default combined error, add() + prop(), so
    # sigma_add is addSd and sigma_mult is propSd.
    #
    # Table 1 lists sigma_mult explicitly in the estimated block and the
    # remaining 'Res. Error' row is sigma_add; this is the only arm in the
    # paper that carries two residual rows, which is what a two-parameter
    # error model requires. On the assay scale sigma_add = 15.18 ng/mL
    # sits sensibly above the 5 ng/mL LOQ.
    # ------------------------------------------------------------------
    addSd <- 15.18; label("Additive residual error SD (ng/mL)") # Table 1, PO CF: Res. Error = 15.18 (RSE 35.72%)
    propSd <- 0.41; label("Proportional residual error (fraction)") # Table 1, PO CF: sigma_mult = 0.41 (RSE 18.17%)
  })

  model({
    # 1. Individual parameters (exponential IIV, Methods equation (4)).
    ka <- exp(lka + etalka)
    tlag <- exp(ltlag + etaltlag)
    vc <- exp(lvc + etalvc)
    cl <- exp(lcl + etalcl)

    # 2. Micro-constant.
    kel <- cl / vc

    # 3. ODE system -- Methods equation (3) written as differential
    # equations. The oral dose enters the depot in full; bioavailability
    # is absorbed into the apparent Vd/F and Cl/F.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # 4. Absorption lag on the depot (Table 1, Tlag).
    alag(depot) <- tlag

    # 5. Observation. Amount in ug/kg over an apparent volume in L/kg
    # gives ug/L, which is ng/mL -- the assay units of Figure 2c.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
