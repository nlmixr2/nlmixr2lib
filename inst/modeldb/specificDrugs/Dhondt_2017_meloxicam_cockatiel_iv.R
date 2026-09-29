Dhondt_2017_meloxicam_cockatiel_iv <- function() {
  description <- paste(
    "Preclinical (cockatiel).",
    "One-compartment population PK model with first-order elimination for",
    "meloxicam after a single 1 mg/kg intravenous bolus of the commercial",
    "injectable formulation (Metacam suspension for injection) to",
    "cockatiels (Nymphicus hollandicus). Fitted in Phoenix NLME",
    "(FOCE-ELS); Table 3, 'IV CF 1' block, the final model without",
    "enterohepatic recycling. A variant that added enterohepatic recycling",
    "was fitted and rejected: Table 4 gives AIC 542.13 and BIC 560.45 with",
    "recycling against 531.72 and 539.36 without, and the authors",
    "attribute the lack of improvement to the small number of blood",
    "collection points at the secondary peak. The small volume of",
    "distribution (0.173 L/kg) matches the high plasma protein binding,",
    "and the high clearance gives a 0.31 h half-life, far shorter than in",
    "other psittacine species. Volumes and clearances are per kilogram, so",
    "the dosed amount is ug/kg and central/vc lands directly in ng/mL, the",
    "assay units of Figure 4a. Body weight and sex were screened as",
    "covariates and neither was retained. See",
    "Dhondt_2017_meloxicam_cockatiel_oral for the separately fitted oral",
    "arm.",
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
      "Single 1 mg/kg body weight intravenous bolus into the vena cutanea",
      "ulnaris (wing vein) of the commercial injectable formulation",
      "(Metacam suspension for injection, 5 mg/mL, Boehringer Ingelheim).",
      "No analytical standard solution was prepared for meloxicam because",
      "an injectable commercial formulation was available."
    ),
    regions = "Belgium (Ghent University, Merelbeke)",
    notes = paste(
      "The 24 birds (12 male / 12 female) received meloxicam both",
      "intravenously and orally in a two-way crossover with a one-month",
      "washout; the two arms were fitted as separate models and the oral",
      "arm is Dhondt_2017_meloxicam_cockatiel_oral. A sparse sampling",
      "protocol was used because of the limited blood volume of cockatiels:",
      "sampling times were randomly allocated across birds with a maximum",
      "of two samples per bird, drawn before dosing and at 5, 15, 30 and",
      "45 min and 1, 2, 4, 6, 8 and 12 h. Meloxicam was quantified in",
      "plasma by LC-MS/MS over 10-5000 ng/mL (LOQ 10 ng/mL, LOD",
      "0.18 ng/mL); values below the LOQ were excluded before fitting. The",
      "mean profile shows a secondary concentration rise 1-2 h after the",
      "intravenous dose that the authors attribute to enterohepatic",
      "recycling, but the recycling model did not improve the fit (Table 4)",
      "and the structural model carried here does not reproduce that second",
      "peak. Plasma protein binding measured separately was",
      "95.02 +/- 0.01% and is not a model parameter. See Dhondt 2017",
      "Methods 'Animals and experimental procedure', 'Meloxicam PK study'",
      "and 'Population pharmacokinetics meloxicam'."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Structural parameters -- Dhondt 2017 Table 3, 'IV CF 1' block, i.e.
    # the model WITHOUT enterohepatic recycling, which Table 4 selects as
    # final on AIC and BIC. Methods equation (2):
    #   C(t) = C0 * exp(-Cl/Vd * t)
    # a one-compartment model with first-order elimination.
    #
    # Both values are confirmed by the table's own computed secondary
    # parameters: C0 = D/Vd = 1000/0.173 = 5780 ng/mL (Table 3 reports
    # 5775.52), AUC(0-inf) = D/Cl = 1000/0.388 = 2577 ng.h/mL (Table 3
    # reports 2575.66), Ke = Cl/Vd = 2.243 /h (Table 3 reports 2.24) and
    # T1/2el = 0.309 h (Table 3 reports 0.31). The same identities do NOT
    # reproduce the 'IV CF 2' (recycling) block's printed C0 of 385.515 or
    # its AUC of 149.05, which is one more reason that rejected variant is
    # not carried here; see the vignette 'Assumptions and deviations'
    # section.
    # ------------------------------------------------------------------
    lvc <- log(0.173); label("Volume of distribution Vd (log L/kg)") # Table 3, IV CF 1: Vd = 0.173 L/kg (RSE 11.74%)
    lcl <- log(0.388); label("Total body clearance Cl (log L/h/kg)") # Table 3, IV CF 1: Cl = 0.388 L/h.kg (RSE 10.19%)

    # ------------------------------------------------------------------
    # IIV -- Dhondt 2017 Table 3, 'omega' column, read as the SD of eta on
    # the log scale and squared here because nlmixr2 omega entries are
    # variances. Methods equation (4): P_i = theta_P * exp(eta_Pi), with
    # variance omega^2 and 'Interindividual variability is reported as
    # omega'. See Dhondt_2017_celecoxib_cockatiel_iv and the vignette
    # 'Assumptions and deviations' section for why the equations are
    # preferred over the caption's word 'variance'.
    #   Vd: < 0.001 -> encoded at the printed upper bound, 1e-06, which
    #       keeps OMEGA positive definite.
    #   Cl: 0.089   -> 0.007921  (about 8.9% CV)
    # ------------------------------------------------------------------
    etalvc ~ 1e-06 # Table 3, IV CF 1 omega for Vd, reported as '< 0.001'
    etalcl ~ 0.007921 # Table 3, IV CF 1 omega for Cl = 0.089

    # ------------------------------------------------------------------
    # Residual error -- Methods 'Population pharmacokinetics meloxicam':
    # 'Inter- and intra-individual variability were expressed according to
    # the exponential error model and the multiplicative residual error
    # model, respectively.' Equation (5) is C_obs = C_pred * (1 + epsilon),
    # which is nlmixr2's prop().
    # ------------------------------------------------------------------
    propSd <- 0.39; label("Proportional residual error (fraction)") # Table 3, IV CF 1: Res. Error = 0.39 (RSE 26.43%)
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
    # which is ng/mL -- the assay units of Figure 4a.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
