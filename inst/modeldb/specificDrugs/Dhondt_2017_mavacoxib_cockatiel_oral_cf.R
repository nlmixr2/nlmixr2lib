Dhondt_2017_mavacoxib_cockatiel_oral_cf <- function() {
  description <- paste(
    "Preclinical (cockatiel).",
    "One-compartment population PK model with first-order absorption and",
    "first-order elimination for mavacoxib after a single 4 mg/kg oral",
    "intra-crop bolus of the ground commercial tablet formulation",
    "(Trocoxil, CF) to cockatiels (Nymphicus hollandicus). Fitted in",
    "Phoenix NLME (FOCE-ELS); Table 2, 'PO CF' block. No absorption lag",
    "time was retained for any mavacoxib arm. The disposition parameters",
    "are apparent (Vd/F, Cl/F): no bioavailability term was estimated, so",
    "the full 4 mg/kg dose enters the depot. Absolute oral bioavailability",
    "was computed separately from the AUC ratio against the intravenous",
    "arm as F = 111%; unlike celecoxib, mavacoxib absorption was complete",
    "from both the standard solution and the commercial formulation.",
    "Volumes and clearances are per kilogram, so the dosed amount is ug/kg",
    "and central/vc lands directly in ng/mL, the assay units of Figure 3c.",
    "Body weight and sex were screened as covariates and neither was",
    "retained. See Dhondt_2017_mavacoxib_cockatiel_iv and",
    "Dhondt_2017_mavacoxib_cockatiel_oral_std for the other two arms.",
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
        "reference category. This cohort was balanced 13 male / 13 female."
      ),
      source_name = "gender"
    )
  )

  population <- list(
    species = "cockatiel (Nymphicus hollandicus)",
    n_subjects = 26L,
    n_studies = 1L,
    age_range = "6-12 months",
    weight_median = "93 g",
    weight_range = "93 +/- 9 g (mean +/- SD)",
    sex_female_pct = 50,
    disease_state = "Healthy (no disease model)",
    dose_range = paste(
      "Single 4 mg/kg body weight oral intra-crop bolus of the commercial",
      "tablet formulation (Trocoxil 20 mg, Zoetis), ground and blended with",
      "lactose to 4 mg mavacoxib/g and suspended in physiological saline at",
      "2 mg/mL immediately before dosing."
    ),
    regions = "Belgium (Ghent University, Merelbeke)",
    notes = paste(
      "A separate group of 26 birds (13 male / 13 female) from the 40 that",
      "received the standard solution. A sparse sampling protocol was used",
      "because of the limited blood volume of cockatiels: sampling times",
      "were randomly allocated across birds with a maximum of three samples",
      "per bird, drawn before dosing and at 15, 30 and 45 min and 1, 2, 4,",
      "6, 8, 12, 24, 48, 72, 96, 120, 168, 336, 672 and 1008 h. Mavacoxib",
      "was quantified in plasma by LC-MS/MS over 5-5000 ng/mL (LOQ",
      "5 ng/mL, LOD 0.25 ng/mL); values below the LOQ were excluded before",
      "fitting. Feed was withheld from 8 h before until 4 h after dosing,",
      "but each bird received a 2 mL intra-crop feed bolus immediately",
      "after dosing. Plasma protein binding measured separately was",
      "97.02 +/- 0.32% and is not a model parameter. See Dhondt 2017",
      "Methods 'Animals and experimental procedure' and 'Mavacoxib PK",
      "study'."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Structural parameters -- Dhondt 2017 Table 2, 'PO CF' block.
    # Methods equation (3):
    #   C(t) = F*D*ka / (Vd*(ka - ke)) * (exp(-ke*t) - exp(-ka*t))
    # Table 2 reports Vd/F and Cl/F -- apparent values uncorrected for
    # bioavailability -- so F*D/Vd is D/(Vd/F) and the full dose is used.
    #
    # All three values are confirmed by the table's own computed secondary
    # parameters: Ke = (Cl/F)/(Vd/F) = 0.033/6.35 = 0.0052 /h (Table 2
    # reports 0.0051), T1/2el = 133 h (Table 2 reports 135.41),
    # Tmax = ln(ka/ke)/(ka - ke) = 14.5 h (Table 2 reports 14.42),
    # Cmax = 584 ng/mL (Table 2 reports 584.66) and AUC(0-inf) =
    # D/(Cl/F) = 4000/0.033 = 121212 ng.h/mL (Table 2 reports 122962).
    # ------------------------------------------------------------------
    lka <- log(0.28); label("First-order absorption rate constant Ka (log 1/h)") # Table 2, PO CF: Ka = 0.28 /h (RSE 9.97%)
    lvc <- log(6.35); label("Apparent volume of distribution Vd/F (log L/kg)") # Table 2, PO CF: Vd/F = 6.35 L/kg (RSE 0.19%)
    lcl <- log(0.033); label("Apparent total body clearance Cl/F (log L/h/kg)") # Table 2, PO CF: Cl/F = 0.033 L/h.kg (RSE 1.35%)

    # ------------------------------------------------------------------
    # IIV -- Dhondt 2017 Table 2, 'omega' column, read as the SD of eta on
    # the log scale and squared here because nlmixr2 omega entries are
    # variances. See Dhondt_2017_celecoxib_cockatiel_iv and the vignette
    # 'Assumptions and deviations' section for why the equations are
    # preferred over the caption's word 'variance'.
    #   Ka:   0.992 -> 0.984064  (about 130% CV)
    #   Vd/F: 0.090 -> 0.0081    (about 9% CV)
    #   Cl/F: 0.252 -> 0.063504  (about 26% CV)
    # ------------------------------------------------------------------
    etalka ~ 0.984064 # Table 2, PO CF omega for Ka = 0.992
    etalvc ~ 0.0081 # Table 2, PO CF omega for Vd/F = 0.090
    etalcl ~ 0.063504 # Table 2, PO CF omega for Cl/F = 0.252

    # ------------------------------------------------------------------
    # Residual error -- Methods 'Population pharmacokinetics mavacoxib':
    # 'Inter- and intra-individual variability were expressed according to
    # the exponential error model and the multiplicative residual error
    # model, respectively.' Equation (5) is C_obs = C_pred * (1 + epsilon),
    # which is nlmixr2's prop().
    # ------------------------------------------------------------------
    propSd <- 0.28; label("Proportional residual error (fraction)") # Table 2, PO CF: Res. Error = 0.28 (RSE 39.64%)
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
    # gives ug/L, which is ng/mL -- the assay units of Figure 3c.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
