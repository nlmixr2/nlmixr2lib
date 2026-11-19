Stepniewska_2009_amodiaquine <- function() {
  description <- paste(
    "Two-compartment population PK model of desethylamodiaquine (DeAQ), the",
    "principal active metabolite of amodiaquine, after a 3-day oral",
    "artesunate-amodiaquine course (fixed-dose combination or loose tablets) in",
    "African children aged 6 months to 5 years with uncomplicated falciparum",
    "malaria (Stepniewska 2009). Amodiaquine (AQ) dose enters a depot and appears",
    "as DeAQ with a fixed first-order rate constant; all clearances and volumes",
    "are apparent and per kilogram of body weight; IIV on CL/F only; residual",
    "error additive on the log scale."
  )
  reference <- paste(
    "Stepniewska K, Taylor W, Sirima SB, Ouedraogo EB, Ouedraogo A, Gansane A,",
    "Simpson JA, Morgan CC, White NJ, Kiechel JR.",
    "Population pharmacokinetics of artesunate and amodiaquine in African children.",
    "Malar J. 2009;8(Suppl 1):200. doi:10.1186/1475-2875-8-200.",
    sep = " "
  )
  vignette <- "Stepniewska_2009_artesunate_amodiaquine"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # The depot receives the administered amodiaquine base dose (mg). The paper
  # modelled DeAQ concentrations against the AQ dose with no molecular-weight
  # correction: the published AUC normalised to a 34 mg/kg dose (median 57.76
  # and 54.75 ug*h/mL by arm) equals 34 / CL/F with CL/F = 0.61-0.63 L/kg/h,
  # whereas a molar correction (327.8 / 355.9) would give about 51. Central
  # therefore holds DeAQ in AQ-dose mass units.
  compartmentData <- list(
    depot = list(analyte = "amodiaquine", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "desethylamodiaquine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "desethylamodiaquine", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "All structural parameters are reported per kilogram (Table 3: CL/F and Q in",
        "L/kg/h, Vcentral/F and Vperipheral/F in L/kg), so each is multiplied by WT",
        "to give absolute values for a dose in mg. No further weight effect was",
        "retained ('Inclusion of any of the covariates ... did not reveal any",
        "statistically significant associations with clearance', Results)."
      ),
      source_name = "weight"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened on CL/F (paper records age in months); not retained (Results, 'Desethylamodiaquine')."
    ),
    BODYTEMP = list(
      description = "Body temperature at enrolment",
      units = "degC",
      type = "continuous",
      notes = "Screened; not retained."
    ),
    SEXF = list(
      description = "Sex (1 = female)",
      units = "(binary)",
      type = "binary",
      notes = "Screened; not retained."
    ),
    LNPC = list(
      description = "Natural log of enrolment parasitaemia",
      units = "log(parasites/uL)",
      type = "continuous",
      notes = "Screened; not retained."
    ),
    RESP_RATE = list(
      description = "Respiratory rate at enrolment",
      units = "breaths/min",
      type = "continuous",
      notes = "Screened; not retained."
    ),
    GAMETOCYTE_POS = list(
      description = "Gametocytaemia present at enrolment (1 = yes)",
      units = "(binary)",
      type = "binary",
      notes = "Screened; not retained."
    ),
    FORM_FDC = list(
      description = "Fixed-dose co-formulation (1) versus loose tablets (0)",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Screened ('Inclusion of any of the covariates, including drug formulation,",
        "... did not reveal any statistically significant associations with",
        "clearance'); the loose-to-fixed AUC ratio was 1.043 (95% CI 0.956-1.138)."
      )
    )
  )

  population <- list(
    species = "human",
    n_subjects = 70L,
    n_subjects_pk = "70 children in the amodiaquine PK group (35 per arm); 3 excluded for vomiting and 2 for pre-dose amodiaquine, leaving 65 in the DeAQ analysis.",
    n_studies = 1L,
    age_range = "7-59 months (amodiaquine PK group, Table 1)",
    age_median = "26 months (loose tablets) and 29 months (fixed-dose combination)",
    weight_range = "6-31 kg",
    weight_median = "12 kg in both arms",
    sex_female_pct = 48.6,
    race_ethnicity = "Not reported by race; Burkinabe children",
    disease_state = "Acute uncomplicated Plasmodium falciparum malaria (>1,000 parasites/uL, axillary temperature >= 37.5 degC)",
    dose_range = paste(
      "Amodiaquine base once daily for 3 days, target 10 mg/kg/day (therapeutic",
      "range 7.5-15): fixed-dose combination 67.5 mg AQ per tablet, 1 tablet (<12",
      "months) or 2 tablets (12-60 months); loose Flavoquine 153 mg tablets, half or",
      "one tablet by the same age bands. Median total AQ dose about 34 mg/kg."
    ),
    sampling = "Pre-dose 1, 4 h after dose 3, Day 7 or 14, and Day 21 or 28.",
    regions = "Burkina Faso (Pissy health district, Ouagadougou)",
    notes = paste(
      "Demographics from Stepniewska 2009 Table 1 (amodiaquine PK group);",
      "sex_female_pct derived from 16/35 and 20/35 males. Analysis in S-PLUS 6.0."
    )
  )

  ini({
    # Fixed because each child contributed at most three detectable samples
    # (Methods, 'Pharmacokinetic analysis': 'the absorption rate constant was
    # fixed (as 0.13 h-1 from previously published work')).
    lka <- fixed(log(0.13))
    label("First-order appearance rate constant of DeAQ (1/h)") # Table 3: Ka = 0.13 (fixed)
    lcl <- log(0.610)
    label("Apparent DeAQ clearance per kg (CL/F, L/kg/h)") # Table 3: CL/F = 0.610 (SE 0.038)
    lq <- log(0.680)
    label("Apparent intercompartmental clearance per kg (Q/F, L/kg/h)") # Table 3: Q = 0.680 (SE 0.306)
    lvc <- log(35.4)
    label("Apparent central volume per kg (Vc/F, L/kg)") # Table 3: Vcentral/F = 35.4 (SE 11.5)
    lvp <- log(87.9)
    label("Apparent peripheral volume per kg (Vp/F, L/kg)") # Table 3: Vperipheral/F = 87.9 (SE 17.1)

    # IIV is log-normal and the reported CV% is the square root of the variance
    # estimate (Methods), so the variance is 0.227^2. Source: Table 3, CL/F
    # inter-subject variation 22.7%.
    etalcl ~ 0.227^2

    expSd <- 0.457
    label("Residual SD on the log scale (additive on log concentration)") # Table 3: sigma = 0.457, 'residual additive error on the loge scale'
  })

  model({
    # Per-kilogram parameters scaled to absolute values for a dose in mg.
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * WT
    q <- exp(lq) * WT
    vc <- exp(lvc) * WT
    vp <- exp(lvp) * WT

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # central (mg) / vc (L) is mg/L; x1000 gives ng/mL.
    Cc <- 1000 * central / vc
    Cc ~ lnorm(expSd)
  })
}
