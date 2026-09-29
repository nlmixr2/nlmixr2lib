Dhondt_2017_celecoxib_cockatiel_oral_std <- function() {
  description <- paste(
    "Preclinical (cockatiel).",
    "One-compartment population PK model with first-order absorption and",
    "first-order elimination for celecoxib after a single 10 mg/kg oral",
    "intra-crop bolus of an analytical standard solution (STD) to",
    "cockatiels (Nymphicus hollandicus). Fitted in Phoenix NLME",
    "(FOCE-ELS); Table 1, 'PO STD' block. No absorption lag time was",
    "retained for this arm. The disposition parameters are apparent",
    "(Vd/F, Cl/F): no bioavailability term was estimated, so the full",
    "10 mg/kg dose enters the depot. Absolute oral bioavailability was",
    "computed separately from the AUC ratio against the intravenous arm",
    "as F = 110%. Volumes and clearances are per kilogram, so the dosed",
    "amount is ug/kg and central/vc lands directly in ng/mL, the assay",
    "units of Figure 2b. Body weight and sex were screened as covariates",
    "and neither was retained. See Dhondt_2017_celecoxib_cockatiel_iv for",
    "the intravenous arm and Dhondt_2017_celecoxib_cockatiel_oral_cf for",
    "the commercial-formulation oral arm.",
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
        "reference category. Cohorts were balanced 17 male / 17 female."
      ),
      source_name = "gender"
    )
  )

  population <- list(
    species = "cockatiel (Nymphicus hollandicus)",
    n_subjects = 34L,
    n_studies = 1L,
    age_range = "6-12 months",
    weight_median = "91 g",
    weight_range = "91 +/- 10 g (mean +/- SD)",
    sex_female_pct = 50,
    disease_state = "Healthy (no disease model)",
    dose_range = paste(
      "Single 10 mg/kg body weight oral intra-crop bolus of a celecoxib",
      "analytical standard solution (5 mg/mL in polyethylene glycol",
      "400:physiological saline, 75:25 v/v), given by a curved ball-tipped",
      "feeding needle."
    ),
    regions = "Belgium (Ghent University, Merelbeke)",
    notes = paste(
      "The same 34 birds (17 male / 17 female) received celecoxib STD both",
      "orally and intravenously in a two-way crossover with a one-month",
      "washout; the two arms were fitted as separate models and the",
      "intravenous arm is Dhondt_2017_celecoxib_cockatiel_iv. A sparse",
      "sampling protocol was used because of the limited blood volume of",
      "cockatiels: sampling times were randomly allocated across birds with",
      "a maximum of two samples per bird, drawn before dosing and at 15, 30",
      "and 45 min and 1, 2, 4, 6, 8, 12 and 24 h. Celecoxib was quantified",
      "in plasma by LC-MS/MS over 5-5000 ng/mL (LOQ 5 ng/mL, LOD",
      "0.22 ng/mL); values below the LOQ were excluded before fitting.",
      "Feed was withheld from 8 h before until 4 h after dosing, but each",
      "bird received a 2 mL intra-crop feed bolus immediately after dosing.",
      "Plasma protein binding measured separately was 98.98 +/- 0.07% and",
      "is not a model parameter. See Dhondt 2017 Methods 'Animals and",
      "experimental procedure' and 'Celecoxib PK study'."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Structural parameters -- Dhondt 2017 Table 1, 'PO STD' block.
    # Methods equation (3):
    #   C(t) = F*D*ka / (Vd*(ka - ke)) * (exp(-ke*t) - exp(-ka*t))
    # i.e. one compartment with first-order absorption and elimination.
    # Table 1 reports Vd/F and Cl/F -- apparent values uncorrected for
    # bioavailability -- so F*D/Vd is D/(Vd/F) and the full dose is used.
    #
    # All three values are confirmed by the table's own computed secondary
    # parameters: Ke = (Cl/F)/(Vd/F) = 2.19/5.87 = 0.373 /h (Table 1
    # reports 0.37), Tmax = ln(ka/ke)/(ka - ke) = 3.14 h (Table 1 reports
    # 3.11), Cmax = 528 ng/mL (Table 1 reports 535.51) and AUC(0-inf) =
    # D/(Cl/F) = 10000/2.19 = 4566 ng.h/mL (Table 1 reports 4573.68).
    #
    # No Tlag row appears in the 'PO STD' block: Results states that a lag
    # time was included only 'for the oral CF of celecoxib and meloxicam'.
    # ------------------------------------------------------------------
    lka <- log(0.27); label("First-order absorption rate constant Ka (log 1/h)") # Table 1, PO STD: Ka = 0.27 /h (RSE 8.98%)
    lvc <- log(5.87); label("Apparent volume of distribution Vd/F (log L/kg)") # Table 1, PO STD: Vd/F = 5.87 L/kg (RSE 17.43%)
    lcl <- log(2.19); label("Apparent total body clearance Cl/F (log L/h/kg)") # Table 1, PO STD: Cl/F = 2.19 L/h.kg (RSE 17.86%)

    # ------------------------------------------------------------------
    # IIV -- Dhondt 2017 Table 1, 'omega' column, read as the SD of eta on
    # the log scale and squared here because nlmixr2 omega entries are
    # variances. Methods equation (4): P_i = theta_P * exp(eta_Pi), with
    # variance omega^2 and 'Interindividual variability is reported as
    # omega'. The table captions instead call omega a variance; the
    # equations win, for the reason set out in
    # Dhondt_2017_celecoxib_cockatiel_iv and in the vignette 'Assumptions
    # and deviations' section.
    #   Vd/F: 0.948 -> 0.898704   (about 121% CV)
    #   Cl/F: 0.707 -> 0.499849   (about 79% CV)
    #   Ka:   0.063 -> 0.003969   (about 6.3% CV)
    # ------------------------------------------------------------------
    etalka ~ 0.003969 # Table 1, PO STD omega for Ka = 0.063
    etalvc ~ 0.898704 # Table 1, PO STD omega for Vd/F = 0.948
    etalcl ~ 0.499849 # Table 1, PO STD omega for Cl/F = 0.707

    # ------------------------------------------------------------------
    # Residual error -- Methods: 'For the IV and PO STD dosing, a
    # multiplicative (equation (5)) and additive (equation (6)) error model
    # was used ... respectively'. Equation (6) is C_obs = C_pred + epsilon,
    # which is nlmixr2's add(). The value is on the assay scale (ng/mL),
    # a little above the 5 ng/mL LOQ.
    # ------------------------------------------------------------------
    addSd <- 8.02; label("Additive residual error SD (ng/mL)") # Table 1, PO STD: Res. Error = 8.02 (RSE 34.93%)
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
    # gives ug/L, which is ng/mL -- the assay units of Figure 2b.
    Cc <- central / vc
    Cc ~ add(addSd)
  })
}
