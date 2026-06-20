Foo_2016_droperidol <- function() {
  description <- paste(
    "Two-compartment population PK model with first-order absorption",
    "for intramuscular droperidol in 41 acutely agitated adults",
    "presenting to the emergency department (Foo 2016). Absorption rate",
    "constant ka and its IIV are fixed (ka = 10 1/h, omega_ka^2 = 1)",
    "because the available samples did not characterise absorption. A",
    "single shared random effect drives both CL and Vc (Table 2 footnote",
    "a: 'The same random effect was used for both Vc and CL'); Q and Vp",
    "have no IIV. No covariates were retained -- coingestion of alcohol",
    "was screened but not associated with CL or Vc, and patient weight",
    "was not available (Methods)."
  )
  reference <- paste(
    "Foo LK, Duffull SB, Calver L, Schneider J, Isbister GK.",
    "Population pharmacokinetics of intramuscular droperidol in",
    "acutely agitated patients.",
    "Br J Clin Pharmacol. 2016;82(6):1550-1556.",
    "doi:10.1111/bcp.13093."
  )
  vignette <- "Foo_2016_droperidol"
  units <- list(
    time          = "h",
    dosing        = "mg",
    concentration = "ug/L"
  )

  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = 41L,
    n_studies      = 1L,
    age_range      = "16-62 years",
    age_median     = "33 years",
    weight_range   = "not recorded",
    weight_median  = "not recorded",
    sex_female_pct = 65.9,
    race_ethnicity = c(Unreported = 100),
    disease_state  = paste(
      "Adults (one 16-year-old recruited inadvertently) presenting to",
      "the emergency department with acute behavioural disturbance",
      "requiring physical restraint and parenteral sedation. Primary",
      "reasons for presentation: threatened or deliberate self-harm",
      "44%, alcohol intoxication 41%, drug-induced delirium 10%,",
      "psychosis 5%. Urine drug screen positive in 22%; breath alcohol",
      "levels recorded in 32 of 41 (median 0.24, range 0.01-0.39",
      "mg/dL)."
    ),
    dose_range     = paste(
      "Single intramuscular dose of 5 mg (n = 17) or 10 mg (n = 24);",
      "5 of 41 patients received additional droperidol after the",
      "index dose."
    ),
    regions        = "Australia (Calvary Mater Newcastle emergency department, New South Wales)",
    notes          = paste(
      "Subgroup of the randomised controlled trial ACTRN12607000527460",
      "comparing droperidol vs midazolam for sedation of acute",
      "behavioural disturbance. 128 plasma samples in total (median 3",
      "per subject, range 1-7) drawn under an empirical geometrically",
      "spaced design (5, 10, 30 min, 1, 2, 4, 8 h) initially, then a",
      "POPT-optimised design (5, 25, 40, 70 min, 2, 4, 10 h). Assay:",
      "HPLC-UV with LLOQ 5 ug/L; 4 observations below LLOQ in 4",
      "different subjects handled by M6 (first BLQ in each series set",
      "to LLOQ/2, subsequent commented out). Patient weight was not",
      "available (Methods); allometric / linear weight scaling could",
      "not be tested. Coingestion of alcohol was tested visually",
      "against CL and Vc and showed no association (Results)."
    )
  )

  ini({
    # ================================================================
    # Absorption -- ka and its IIV are FIXED.
    # Results: "an estimate of the first order rate constant of
    # absorption when fixed to 10 h-1 provided a stable model and
    # lowest objective function ... a value of between subject
    # variance of 1 (approximately equivalent to a between subject
    # CV% of 100%) ... was stable and provided the lowest value of
    # the objective function."
    # ================================================================
    lka <- fixed(log(10))
    label("Absorption rate constant (1/h, FIXED)")                                  # Table 2: ka = 10 h-1 (F = fixed)

    # ================================================================
    # Structural disposition -- Foo 2016 Table 2 final model.
    # Table 2 has a typographic duplication of 'Vp' in the parameter
    # name column; the abstract disambiguates: 'clearance of 41.9
    # l h-1 and volume of distribution of the central compartment of,
    # 73.6 l'. The first '73.6' row is therefore Vc, the second '79.8'
    # row is Vp.
    # ================================================================
    lcl <- log(41.9)
    label("Clearance CL (L/h)")                                                     # Table 2: CL = 41.9 (95% CI 34.8-49.0) L/h
    lvc <- log(73.6)
    label("Central volume of distribution Vc (L)")                                  # Table 2 (Vc, mislabelled 'Vp' in Table 2; confirmed by abstract): Vc = 73.6 (95% CI 51.1-96.1) L
    lq  <- log(71.5)
    label("Intercompartmental clearance Q (L/h)")                                   # Table 2: Q = 71.5 (95% CI 42.3-100.7) L/h
    lvp <- log(79.8)
    label("Peripheral volume of distribution Vp (L)")                               # Table 2: Vp = 79.8 (95% CI 58.8-100.8) L

    # ================================================================
    # Inter-individual variability.
    # ka: variance fixed at 1 (= log(1 + ~100% CV)^2 approximately;
    # the paper states 'a value of between subject variance of 1
    # (approximately equivalent to a between subject CV% of 100%)').
    # CL and Vc: a single shared random effect drives both parameters
    # (Table 2 footnote a). Reported as CV% = 51% for CL; Vc carries
    # the same value, footnote a. Encoded as a single eta (etalcl)
    # referenced from both cl and vc in model(); omega^2 =
    # log(1 + 0.51^2) = 0.231.
    # Q and Vp: '-' in Table 2 (no random effect estimated).
    # ================================================================
    etalka ~ fixed(1)                                                               # Table 2: omega_ka^2 = 1, FIXED (~ 100% CV)
    etalcl ~ log(1 + 0.51^2)                                                        # Table 2: CV_CL = 51% (95% CI 31.2-64.4%); same eta also drives Vc per footnote a -> omega^2 = log(1 + 0.51^2) = 0.231

    # ================================================================
    # Residual error -- combined proportional + additive (Methods:
    # 'two-compartment first-order input with first-order output
    # model with combined error model'). The additive component was
    # fixed at 0.0001 ug/L for numerical stability.
    # ================================================================
    propSd <- 0.22
    label("Proportional residual error (fraction)")                                  # Table 2: sigma (CV%) = 22% (95% CI 8.5-30.3%)
    addSd  <- fixed(0.0001)
    label("Additive residual error (ug/L, FIXED)")                                   # Table 2: sigma_add = 0.0001 ug/L (F = fixed)
  })

  model({
    # ----------------------------------------------------------------
    # Individual structural parameters.
    # CL and Vc share the single random effect etalcl per Table 2
    # footnote a. ka has its own (fixed-variance) eta. Q and Vp have
    # no IIV.
    # ----------------------------------------------------------------
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalcl)
    q  <- exp(lq)
    vp <- exp(lvp)

    # ----------------------------------------------------------------
    # Two-compartment model with first-order absorption from the IM
    # depot. The internal amounts are in mg (matching the dose unit);
    # the * 1000 converts mg/L to ug/L for the reported observation.
    # ----------------------------------------------------------------
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka  * depot
    d/dt(central)     <-  ka  * depot      - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central    - k21 * peripheral1

    Cc <- central / vc * 1000
    Cc ~ prop(propSd) + add(addSd)
  })
}
