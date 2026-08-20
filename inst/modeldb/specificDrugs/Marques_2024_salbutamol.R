Marques_2024_salbutamol <- function() {
  description <- paste(
    "One-compartment population PK model for inhaled salbutamol (600 ug",
    "single dose, unit-dose dry-powder inhaler) in healthy adults, with",
    "first-order absorption, no lag time, and linear elimination. Encodes",
    "the EXTERNAL-VALIDATION fit to real clinical data (GSK study",
    "NCT01984086, n = 30; Marques 2024 Table 4) -- NOT the paper's headline",
    "virtual-patient model (Table 3), whose central volume of 0.02e-9 L is",
    "not solvable at any ODE tolerance and overshoots the paper's own",
    "observed Cmax by four orders of magnitude. The authors declare a",
    "two-compartment structure, but their clinical-data fit returned",
    "Q = 1.30e-7 L/h and V2 = 0 exactly, both with relative standard errors",
    "above 1e6 percent (the authors' own footnote calls this",
    "overparametrization); the peripheral compartment therefore carries no",
    "flux and k21 = Q/V2 is undefined, so the model is encoded honestly as",
    "one compartment. No covariate effects are encoded: the paper reports",
    "which covariates were significant but publishes no coefficients,",
    "functional forms, or centering values anywhere in the article or its",
    "supplement (see covariatesDataExcluded and the vignette Errata)."
  )
  reference <- paste(
    "Marques L, Vale N. Toward Personalized Salbutamol Therapy: Validating",
    "Virtual Patient-Derived Population Pharmacokinetic Model with",
    "Real-World Data. Pharmaceutics. 2024;16(7):881.",
    "doi:10.3390/pharmaceutics16070881. Parameters from Table 4 (final",
    "model refitted to the GlaxoSmithKline clinical dataset, NCT01984086",
    "Part A Treatment A, accessed via Vivli)."
  )
  vignette <- "Marques_2024_salbutamol"
  # Dose is in ug (600 ug nominal) and volumes are in L, so central / vc is
  # ug/L = ng/mL. The paper tabulates concentrations in ug/mL; 1 ug/mL =
  # 1000 ng/mL. No values are converted here -- only the declared spelling
  # of the unit the model actually produces.
  units <- list(time = "h", dosing = "ug", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(
      analyte = "salbutamol", units = "ug",
      specimen = "administration site", verified = TRUE
    ),
    central = list(
      analyte = "salbutamol", units = "ug",
      specimen = "plasma", verified = TRUE
    )
  )

  # No covariate enters model(): Table 4 (the fit encoded here) has no
  # covariate rows at all. The covariates below were screened by the authors
  # -- on the SEPARATE virtual-patient model (Table 3) -- and reported as
  # statistically significant in Section 3.4, but the paper and its
  # supplement publish no coefficient, functional form, reference value or
  # centering constant for any of them, so none can be encoded. They are
  # documented here so the paper's headline claim is not lost.
  covariatesDataExcluded <- list(
    AGE = list(
      description        = "Age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Section 3.4 reports a significant effect of age on Cl and Q in the",
        "virtual-patient model; no coefficient is published. Table 5",
        "stratifies individual EBE parameters at 5-22 vs 23-65 years, but",
        "those are descriptive subgroup summaries, not covariate effects,",
        "and are in units irreconcilable with Table 3."
      ),
      source_name        = "Age"
    ),
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Section 3.4 reports a significant effect of weight on Cl, Q and V2",
        "in the virtual-patient model; no coefficient is published. Table 5",
        "stratifies at 17.77-75.00 vs 75.01-105.00 kg (descriptive only)."
      ),
      source_name        = "Weight"
    ),
    HT = list(
      description        = "Height",
      units              = "cm",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened as a continuous covariate (Section 2.4) and correlated",
        "with AUC, Cmax and clearance in the data exploration (Section 3.1),",
        "but not retained as a reported covariate effect on any parameter."
      ),
      source_name        = "Height"
    ),
    SEXF = list(
      description        = "Female sex indicator (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = male",
      notes              = paste(
        "Section 3.4 reports a significant effect of gender on ka and Cl in",
        "the virtual-patient model; no coefficient is published. Table 5",
        "reports both parameters as lower in females (descriptive only).",
        "The source reports gender as a two-level 'Female' / 'Male'",
        "categorical; SEXF is the canonical encoding."
      ),
      source_name        = "Gender"
    ),
    RACE_ASIAN_NORTHEAST = list(
      description        = paste(
        "North East Asian heritage indicator (1 = Japanese or Chinese,",
        "0 = American Indian / Alaskan Native)"
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 = American Indian / Alaskan Native",
      notes              = paste(
        "Section 3.4 reports a significant effect of race on Cl in the",
        "virtual-patient model; no coefficient is published. The virtual",
        "cohort has exactly two race levels (Table 1: American Indian or",
        "Alaskan Native 38%, East Asian [Japanese + Chinese] 62%), so this",
        "single indicator spans the whole covariate. Note that the CLINICAL",
        "cohort encoded by this model contains no East Asian subjects at all",
        "(Table 1: 90% White, 7% Mixed, 3% American Indian / Alaskan",
        "Native), which is a further reason the covariate layer cannot be",
        "carried over to the Table 4 fit."
      ),
      source_name        = "Race"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 30,
    n_studies      = 1,
    age_range      = "18-65 years (eligibility criterion)",
    age_median     = "26.8 years (mean, SD 4.8; Table 1)",
    weight_range   = ">= 50 kg (eligibility criterion); BMI 19.0-34.0 kg/m^2",
    weight_median  = "78.9 kg (mean, SD 17.6; Table 1)",
    sex_female_pct = 26.7,
    race_ethnicity = c(
      White = 90, Multi = 7, AmericanIndianOrAlaskanNative = 3, Asian = 0
    ),
    disease_state  = "Healthy nonsmoking volunteers (no asthma)",
    dose_range     = paste(
      "600 ug salbutamol (as sulfate) as a single dose: 3 x 200 ug blisters",
      "inhaled via unit-dose dry-powder inhaler (UD-DPI)"
    ),
    regions        = NA_character_,
    notes          = paste(
      "GlaxoSmithKline study NCT01984086, Part A Treatment A, accessed",
      "through Vivli. Open-label, randomized, crossover, two-cohort,",
      "single-dose study in healthy volunteers. Baseline demographics are",
      "the 'Clinical Study Patients (n = 30)' column of Marques 2024",
      "Table 1; height 177.8 cm (SD 9.5), BMI 24.7 kg/m^2 (SD 3.8). PK",
      "sampling at 0, 0.08, 0.17, 0.33, 0.50, 0.75, 1, 1.5, 2, 4, 6, 8, 10",
      "and 12 h. The virtual-patient cohort (n = 32) that produced Table 3",
      "is a DIFFERENT population and is not the population of this model."
    )
  )

  ini({
    # Structural parameters -- Marques 2024 Table 4 ("Estimates of the popPK
    # parameters of the final model for salbutamol DPI formulation (using
    # GSK-derived data)"), the external-validation fit. Bioavailability was
    # not estimated (Monolix F = 1 against the full 600 ug nominal dose), so
    # Cl and V1 are apparent (CL/F, Vc/F).
    lka <- log(13.55);  label("Absorption rate constant ka (1/h)")            # Table 4 Fixed Effects: ka = 13.55 1/h (RSE 18.0%)
    lcl <- log(34.93);  label("Apparent clearance CL/F (L/h)")                # Table 4 Fixed Effects: Cl = 34.93 L/h (RSE 16.6%)
    lvc <- log(162.93); label("Apparent central volume of distribution Vc/F (L)") # Table 4 Fixed Effects: V1 = 162.93 L (RSE 22.6%)

    # Table 4 also lists Q = 1.30e-7 L/h (RSE 2.83e6 %) and V2 = 0 (RSE
    # 1.35e7 %), flagged by the authors' own footnote as "very large standard
    # errors, potentially suggesting an overparametrization of the model".
    # V2 = 0 exactly makes k21 = Q/V2 undefined, and Q is eight orders of
    # magnitude below Cl, so the peripheral compartment carries no flux.
    # There is no honest two-compartment encoding of Table 4; lvp / lq and
    # their IIVs (IIV(Q) 1.17, IIV(V2) 0.30, corr(ka,Q) -0.83, all with RSE
    # > 1e6 %) are deliberately omitted. See the vignette Errata.

    # IIV. Monolix reports random effects as omega, the STANDARD DEVIATION of
    # the (log-normally distributed) random effect; nlmixr2's ini() takes
    # VARIANCES, so each value below is the published omega squared.
    etalka ~ 0.7921                                  # Table 4 Random Effects: IIV(ka) = 0.89 (RSE 15.4%); 0.89^2 = 0.7921
    etalcl + etalvc ~ c(0.2601,
                        0.192372, 0.1681)            # Table 4: IIV(Cl) = 0.51 (RSE 13.3%) -> 0.51^2 = 0.2601; IIV(V1) = 0.41 (RSE 13.4%) -> 0.41^2 = 0.1681; Correlation V1 and Cl = 0.92 (RSE 3.98) -> 0.92 * 0.51 * 0.41 = 0.192372

    # Residual error. Section 2.3 states the Monolix "combined 1" error model
    # (additive + proportional) was used for both fits, but neither Table 3
    # nor Table 4 nor Supplementary Table S2 publishes the a and b
    # coefficients. Both terms are encoded so the declared structure is
    # preserved, and both are fixed at zero rather than invented -- the
    # standing policy for unreported residual variability. A user refitting
    # on real data must relax these.
    addSd  <- fixed(0); label("Additive residual standard deviation (ng/mL; ZERO - not reported in source)")     # Section 2.3 declares "combined 1"; coefficient never published
    propSd <- fixed(0); label("Proportional residual standard deviation (fraction; ZERO - not reported in source)") # Section 2.3 declares "combined 1"; coefficient never published
  })

  model({
    # 1. Individual PK parameters (log-normal, per Section 2.3 "log-normal
    #    distributions for PK parameters").
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)

    # 2. Micro-constants.
    kel <- cl / vc

    # 3. ODE system: first-order absorption, no lag time, linear elimination
    #    (Section 3.2 / Section 3.3 structural description, reduced to one
    #    compartment as described in ini()).
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # 4. Observation. Dose in ug / volume in L gives ug/L = ng/mL.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
