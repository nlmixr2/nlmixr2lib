Dhondt_2017_meloxicam_cockatiel_oral <- function() {
  description <- paste(
    "Preclinical (cockatiel).",
    "One-compartment population PK model with a lagged first-order",
    "absorption and first-order elimination for meloxicam after a single",
    "1 mg/kg oral intra-crop bolus of the commercial oral suspension",
    "(Metacam oral suspension) to cockatiels (Nymphicus hollandicus).",
    "Fitted in Phoenix NLME (FOCE-ELS); Table 3, 'PO CF 1' block, the",
    "final model without enterohepatic recycling. The disposition",
    "parameters are apparent (Vd/F, Cl/F): no bioavailability term was",
    "estimated, so the full 1 mg/kg dose enters the depot. Absolute oral",
    "bioavailability was computed separately from the AUC ratio against",
    "the intravenous arm as F = 11%, strikingly low against the 38-74%",
    "reported for the same oral suspension in other avian species and the",
    "85-106% reported in mammals. Volumes and clearances are per kilogram,",
    "so the dosed amount is ug/kg and central/vc lands directly in ng/mL,",
    "the assay units of Figure 4b. Body weight and sex were screened as",
    "covariates and neither was retained. See",
    "Dhondt_2017_meloxicam_cockatiel_iv for the intravenous arm.",
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
    depot = list(analyte = "meloxicam", units = "ug", specimen = "administration site", verified = TRUE),
    central = list(analyte = "meloxicam", units = "ug", specimen = "plasma", verified = TRUE)
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
        "reference category. This cohort was balanced 12 male / 12 female."
      ),
      source_name = "gender"
    )
  )

  population <- list(
    species = "cockatiel (Nymphicus hollandicus)",
    n_subjects = 24L,
    n_studies = 1L,
    age_range = "6-12 months",
    weight_median = "101 g",
    weight_range = "101 +/- 12 g (mean +/- SD)",
    sex_female_pct = 50,
    disease_state = "Healthy (no disease model)",
    dose_range = paste(
      "Single 1 mg/kg body weight oral intra-crop bolus of the commercial",
      "oral suspension (Metacam oral suspension, 0.5 mg/mL, Boehringer",
      "Ingelheim), given by a curved ball-tipped feeding needle."
    ),
    regions = "Belgium (Ghent University, Merelbeke)",
    notes = paste(
      "The same 24 birds (12 male / 12 female) received meloxicam both",
      "orally and intravenously in a two-way crossover with a one-month",
      "washout; the two arms were fitted as separate models and the",
      "intravenous arm is Dhondt_2017_meloxicam_cockatiel_iv. A sparse",
      "sampling protocol was used because of the limited blood volume of",
      "cockatiels: sampling times were randomly allocated across birds with",
      "a maximum of two samples per bird, drawn before dosing and at 15, 30",
      "and 45 min and 1, 2, 4, 6, 8, 12 and 24 h. Meloxicam was quantified",
      "in plasma by LC-MS/MS over 10-5000 ng/mL (LOQ 10 ng/mL, LOD",
      "0.18 ng/mL); values below the LOQ were excluded before fitting. With",
      "a Cmax of only about 102 ng/mL the oral profile sits close to the",
      "LOQ, which is consistent with the large residual error. The mean",
      "profile shows a secondary rise about 4 h after the oral dose that",
      "the authors attribute to enterohepatic recycling, but the recycling",
      "model did not improve the fit and the structural model carried here",
      "does not reproduce that second peak. Feed was withheld from 8 h",
      "before until 4 h after dosing, but each bird received a 2 mL",
      "intra-crop feed bolus immediately after dosing. Plasma protein",
      "binding measured separately was 95.02 +/- 0.01% and is not a model",
      "parameter. See Dhondt 2017 Methods 'Animals and experimental",
      "procedure', 'Meloxicam PK study' and 'Population pharmacokinetics",
      "meloxicam'."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Structural parameters -- Dhondt 2017 Table 3, 'PO CF 1' block, i.e.
    # the model WITHOUT enterohepatic recycling. Methods equation (3):
    #   C(t) = F*D*ka / (Vd*(ka - ke)) * (exp(-ke*t) - exp(-ka*t))
    # with an absorption lag time, which Results reports was retained for
    # the oral commercial formulations of celecoxib and meloxicam because
    # it significantly improved the fit. Table 3 reports Vd/F and Cl/F --
    # apparent values uncorrected for bioavailability -- so F*D/Vd is
    # D/(Vd/F) and the full dose is used.
    #
    # All four values are confirmed by the table's own computed secondary
    # parameters: Ke = (Cl/F)/(Vd/F) = 3.38/4.40 = 0.768 /h (Table 3
    # reports 0.77), T1/2el = 0.902 h (Table 3 reports 0.90),
    # Tmax = Tlag + ln(ka/ke)/(ka - ke) = 0.23 + 1.04 = 1.27 h (Table 3
    # reports 1.27), Cmax = 102.5 ng/mL (Table 3 reports 102.31) and
    # AUC(0-inf) = D/(Cl/F) = 1000/3.38 = 296 ng.h/mL (Table 3 reports
    # 295.26).
    # ------------------------------------------------------------------
    lka <- log(1.19); label("First-order absorption rate constant Ka (log 1/h)") # Table 3, PO CF 1: Ka = 1.19 /h (RSE 55.88%)
    ltlag <- log(0.23); label("Absorption lag time Tlag (log h)") # Table 3, PO CF 1: Tlag = 0.23 h (RSE 8.60%)
    lvc <- log(4.40); label("Apparent volume of distribution Vd/F (log L/kg)") # Table 3, PO CF 1: Vd/F = 4.40 L/kg (RSE 27.97%)
    lcl <- log(3.38); label("Apparent total body clearance Cl/F (log L/h/kg)") # Table 3, PO CF 1: Cl/F = 3.38 L/h.kg (RSE 18.82%)

    # ------------------------------------------------------------------
    # IIV -- Dhondt 2017 Table 3, 'omega' column, read as the SD of eta on
    # the log scale and squared here because nlmixr2 omega entries are
    # variances. See Dhondt_2017_celecoxib_cockatiel_iv and the vignette
    # 'Assumptions and deviations' section for why the equations are
    # preferred over the caption's word 'variance'.
    #   Ka:   1.075   -> 1.155625  (about 148% CV)
    #   Tlag: < 0.001 -> encoded at the printed upper bound, 1e-06
    #   Vd/F: < 0.001 -> encoded at the printed upper bound, 1e-06
    #   Cl/F: 0.122   -> 0.014884  (about 12% CV)
    # ------------------------------------------------------------------
    etalka ~ 1.155625 # Table 3, PO CF 1 omega for Ka = 1.075
    etaltlag ~ 1e-06 # Table 3, PO CF 1 omega for Tlag, reported as '< 0.001'
    etalvc ~ 1e-06 # Table 3, PO CF 1 omega for Vd/F, reported as '< 0.001'
    etalcl ~ 0.014884 # Table 3, PO CF 1 omega for Cl/F = 0.122

    # ------------------------------------------------------------------
    # Residual error -- Methods 'Population pharmacokinetics meloxicam':
    # 'Inter- and intra-individual variability were expressed according to
    # the exponential error model and the multiplicative residual error
    # model, respectively.' Equation (5) is C_obs = C_pred * (1 + epsilon),
    # which is nlmixr2's prop(). The value exceeds 100%, consistent with an
    # oral profile whose Cmax of about 102 ng/mL is only ten times the
    # 10 ng/mL LOQ.
    # ------------------------------------------------------------------
    propSd <- 1.15; label("Proportional residual error (fraction)") # Table 3, PO CF 1: Res. Error = 1.15 (RSE 14.19%)
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

    # 4. Absorption lag on the depot (Table 3, Tlag).
    alag(depot) <- tlag

    # 5. Observation. Amount in ug/kg over an apparent volume in L/kg
    # gives ug/L, which is ng/mL -- the assay units of Figure 4b.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
