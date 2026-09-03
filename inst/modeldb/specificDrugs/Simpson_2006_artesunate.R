Simpson_2006_artesunate <- function() {
  description <- paste(
    "One-compartment population PK model of dihydroartemisinin (DHA), the principal",
    "active metabolite, following a single 10 mg/kg intra-rectal artesunate (ARS)",
    "suppository dose in 179 adults and children with moderately severe falciparum",
    "malaria pooled from two Phase II and three Phase III studies in Thailand, Ghana,",
    "Malawi and South Africa (Simpson 2006). DHA appears by a first-order process whose",
    "rate constant (0.2 /h) and lag time (0.14 h) were both fixed because the sparse,",
    "erratic data could not identify them; the appearance rate lumps rectal ARS",
    "absorption with ARS-to-DHA conversion and is slower than DHA elimination, so the",
    "profile is absorption rate limited (flip-flop). Apparent clearance and apparent",
    "volume are reported per kilogram: CL/F is a female value plus an additive male",
    "increment, and V/F is an affine function of body weight normalised to 70 kg.",
    "Artesunate itself could not be modelled - the ARS concentration-time data were too",
    "erratic to support a compartmental fit - so no parent PK is included here."
  )
  reference <- paste(
    "Simpson JA, Agbenyega T, Barnes KI, Di Perri G, Folb P, Gomes M, Krishna S,",
    "Krudsood S, Looareesuwan S, Mansor S, McIlleron H, Miller R, Molyneux M,",
    "Mwenechanya J, Navaratnam V, Nosten F, Olliaro P, Pang L, Ribeiro I, Tembo M,",
    "van Vugt M, Ward S, Weerasuriya K, Win K, White NJ. Population pharmacokinetics of",
    "artesunate and dihydroartemisinin following intra-rectal dosing of artesunate in",
    "malaria patients. PLoS Med. 2006;3(11):e444. doi:10.1371/journal.pmed.0030444"
  )
  vignette  <- "Simpson_2006_artesunate"
  units     <- list(
    time          = "h",
    dosing        = "mg",
    concentration = "mg/L"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. The depot holds the administered artesunate dose at the
  # rectal administration site; the paper assumed complete in-vivo conversion of
  # ARS to DHA and used the actual administered ARS dose (mg/kg) directly in the
  # DHA analysis, with no molar correction (Methods, "Population Pharmacokinetic
  # Modelling"). The mass leaving `depot` therefore enters `central` as DHA
  # one-for-one in mass units, and V/F and CL/F are expressed in ARS-dose terms.
  compartmentData <- list(
    depot   = list(analyte = "artesunate",         units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "dihydroartemisinin", units = "mg", specimen = "plasma",              verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Total body weight at enrolment",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Affine (not power) effect on the per-kilogram apparent volume of distribution,",
        "with the covariate entered as WT/70: V/F (L/kg) = 0.57 + 5.77 * (WT/70)",
        "(Simpson 2006 Table 5, 'Intercept' = 0.57 and 'Increase in V/F per unit (l/70 kg)'",
        "= 5.77). This reproduces every published typical value in Table 5 exactly",
        "(1.81 L/kg at 15 kg, 3.04 at 30 kg, 5.52 at 60 kg, 6.34 at 70 kg). The Methods",
        "state that weight was included 'using allometric scaling, i.e. (weight/70) for",
        "V/F and (weight/70)^0.75 for CL/F', but the final model in Table 5 retains only",
        "the linear (weight/70) term on V/F, added to an intercept; no weight effect on",
        "CL/F was retained ('No correlation was observed between eta_CL/F and body",
        "weight', Results). Because the parameter is per kilogram, absolute V/F (L) is",
        "WT * (0.57 + 5.77 * WT/70). Time-fixed at enrolment.",
        "Cohort weights ranged from 7.6 to 86 kg (Table 2)."
      ),
      source_name        = "weight"
    ),
    SEXF = list(
      description        = "Female sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 = female (the reference category in this encoding)",
      notes              = paste(
        "The paper's covariate is 'gender' and its estimated coefficient is the increase",
        "in CL/F for a MALE relative to a female (Simpson 2006 Table 5: 'Increase in CL/F",
        "(l/kg/h) for a male' = 1.14, SE 0.40). To store it under the canonical SEXF",
        "(1 = female), the male indicator is written (1 - SEXF), so",
        "CL/F (L/kg/h) = 2.03 + 1.14 * (1 - SEXF). Both endpoints are printed in Table 5",
        "(male 3.17, female 2.03) and 2.03 + 1.14 = 3.17 exactly, so neither the sign nor",
        "the reference category is ambiguous. Same encoding as Hennig_2013_tobra.R,",
        "Lee_2023_patritumab.R and Lane_2011_warfarin_s.R. No sex effect on V/F was",
        "retained ('There was no significant difference between males and females for the",
        "distribution of eta_V/F', Results). Time-fixed."
      ),
      source_name        = "gender"
    )
  )

  # Covariates screened by the source paper but not retained in the final model
  # (Simpson 2006 Results, "Relationship between DHA Pharmacokinetic Variables and
  # Covariates"). Documented for provenance only; not referenced in model().
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age at enrolment",
      units       = "years",
      type        = "continuous",
      notes       = paste(
        "Screened both continuously and dichotomised at 15 years. eta_V/F was higher in",
        "adults than children (p < 0.001) but 'the addition of baseline PCV and age did",
        "not significantly improve the objective function' once weight was in the model",
        "(Results, Table 5 paragraph); the univariate age association with V/F was",
        "confounded by weight. Not retained."
      )
    ),
    PCV = list(
      description = "Baseline packed cell volume (haematocrit)",
      units       = "percent",
      type        = "continuous",
      notes       = paste(
        "Positively correlated with eta_V/F (r = 0.22, p = 0.005) but did not",
        "significantly improve the objective function when added to the final model",
        "(Results). Not retained."
      )
    ),
    REGION_SEASIA = list(
      description = "Geographic group indicator (Southeast Asia versus Africa)",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "'There were no significant differences between patients from Southeast Asia and",
        "Africa, for both the distributions of eta_CL/F and eta_V/F' (Results).",
        "Not retained."
      )
    ),
    PARASITEMIA = list(
      description = "Baseline asexual Plasmodium falciparum parasitaemia",
      units       = "parasites/uL",
      type        = "continuous",
      notes       = "No correlation with either eta_CL/F or eta_V/F (Results). Not retained."
    ),
    LACTATE = list(
      description = "Baseline venous plasma lactate",
      units       = "mmol/L",
      type        = "continuous",
      notes       = paste(
        "No correlation with either eta_CL/F or eta_V/F (Results). Not measured at the",
        "Mae-Sot site, whose patients were excluded from the lactate covariate analysis",
        "(Methods, 'Covariates'). Not retained."
      )
    ),
    GLUC = list(
      description = "Baseline plasma glucose",
      units       = "mmol/L",
      type        = "continuous",
      notes       = paste(
        "No correlation with either eta_CL/F or eta_V/F (Results). Not measured at the",
        "Mae-Sot site (Methods, 'Covariates'). Not retained."
      )
    ),
    NSUPP = list(
      description = "Number of artesunate suppositories administered",
      units       = "count",
      type        = "count",
      notes       = paste(
        "'There were no significant differences between both the distributions of",
        "eta_CL/F and eta_V/F for the number of suppositories given (either one, two, or",
        "three or more)' (Results). Not retained."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 179L,
    n_studies      = 5L,
    n_subjects_pk  = "164 of 179 patients had posterior individual estimates of CL/F and V/F; 424 DHA concentrations from 164 patients entered the DHA model (three outliers > 6,000 ng/mL excluded).",
    age_range      = "11 months to 58 years (Table 1: Bangkok 16-50 y, Ghana 2-7 y, Mae-Sot 11 months-15 y, Malawi 16 months-10 y, South Africa 16-58 y)",
    weight_range   = "7.6 to 86 kg (Table 2 per-site ranges; site mean weights 14.2 to 61.9 kg)",
    sex_female_pct = 46.0,
    race_ethnicity = "Not reported by race; sites were Thai (Bangkok, Mae-Sot) and African (Ghana, Malawi, South Africa)",
    disease_state  = "Moderately severe Plasmodium falciparum malaria (asexual parasite density 0.1-27%) in patients unable to tolerate oral medication but without clinical or laboratory features of severe malaria; Mae-Sot enrolled only hyperparasitaemic patients (>= 4%) and Bangkok only patients with parasitaemia > 100,000/uL.",
    dose_range     = "Single intra-rectal artesunate dose targeting 10 mg/kg (50 mg and 200 mg suppositories, Mepha); per-site mean actual doses 9.2-11.0 mg/kg, individual range approximately 6.9-13.6 mg/kg (Table 2).",
    sampling       = "Phase II intensive (8 samples over 12 h in Ghana, 9 samples over 12 h in Bangkok); Phase III fixed 5 samples over 8 h in South Africa and randomised sparse sampling (3 samples in randomised time blocks over 24 h, later 12 h) in Malawi and Mae-Sot.",
    regions        = "Thailand (Bangkok, Mae-Sot), Ghana, Malawi, South Africa",
    notes          = paste(
      "Demographics from Simpson 2006 Tables 1 and 2; baseline laboratory data in Table 3.",
      "sex_female_pct is derived arithmetic, not a printed value: Table 2 reports percent",
      "male per site (Bangkok 42%, Ghana 55%, Mae-Sot 39%, Malawi 62%, South Africa 58%),",
      "which weights to 54.0% male across the 179 patients.",
      "Analysis was performed with the NLME procedure of S-PLUS 6.2, not NONMEM.",
      "Concentrations below the limit of quantification (50 ng/mL; 20 ng/mL at Mae-Sot) and",
      "above the assay range (3,200 ng/mL) were retained in the analysis rather than",
      "censored, and a residual-error covariate distinguishing those records did not",
      "improve the fit."
    )
  )

  ini({
    # ---- Fixed absorption/appearance parameters --------------------------
    # "A one-compartment model with first-order appearance and elimination
    # kinetics including lag time best fits the data (appearance rate was fixed
    # at 0.2/h and lag time at 0.14 h)" (Results, "DHA Pharmacokinetics");
    # "Due to the sparseness and erratic nature of the data, both the appearance
    # rate constant and lag time were fixed so that a satisfactory model fit
    # could be obtained" (Discussion).
    lka   <- fixed(log(0.2))
    label("First-order appearance rate constant of DHA (ka, 1/h)")  # Simpson 2006 Table 5: "k_a (/h) - fixed" = 0.2
    ltlag <- fixed(log(0.14))
    label("Appearance lag time (h)")  # Simpson 2006 Table 5: "Lag time (h) - fixed" = 0.14

    # ---- Apparent clearance (per kg) -------------------------------------
    # Table 5 prints both sex-specific typical values; 2.03 + 1.14 = 3.17
    # exactly, so lcl is the female reference and e_male_cl the male increment.
    lcl <- log(2.03)
    label("Apparent DHA clearance for a female (CL/F, L/kg/h)")  # Simpson 2006 Table 5: "CL/F (l/kg/h) for a: Female" = 2.03
    e_male_cl <- 1.14
    label("Additive increase in apparent DHA clearance for a male (CL/F, L/kg/h; applied as (1 - SEXF))")  # Simpson 2006 Table 5: "Increase in CL/F (l/kg/h) for a male" = 1.14 (SE 0.40); Abstract 95% CI 0.36-1.92

    # ---- Apparent volume of distribution (per kg) ------------------------
    # V/F (L/kg) = intercept + slope * (WT/70). Both terms come straight from
    # Table 5; together they reproduce the four printed typical values
    # (15 kg -> 1.81, 30 kg -> 3.04, 60 kg -> 5.52, 70 kg -> 6.34 L/kg).
    lvc <- log(0.57)
    label("Intercept of the apparent DHA volume of distribution (V/F, L/kg)")  # Simpson 2006 Table 5: "Intercept" = 0.57 (SE 1.68)
    e_wt_vc <- 5.77
    label("Increase in apparent DHA volume of distribution per unit of WT/70 (V/F, L/kg)")  # Simpson 2006 Table 5: "Increase in V/F per unit (l/70 kg)" = 5.77 (SE 0.55)

    # ---- Inter-individual variability ------------------------------------
    # "Inter-individual variability in the pharmacokinetic parameters was
    # modelled with log-normal error models" and "The magnitude of the
    # inter-individual variability was expressed as a coefficient of variation
    # (%CV; approximated by the square root of the variance estimate)"
    # (Methods). The reported %CV is therefore the log-scale SD directly, so the
    # nlmixr2 variance is (%CV/100)^2 -- NOT log(CV^2 + 1). "The data did not
    # support a full variance-covariance matrix for the random effects of CL/F
    # and V/F" (Results), so the two etas are uncorrelated.
    # No trailing comment may be attached to these two lines: rxode2 promotes a
    # trailing comment on an eta line into a label(), and the quotes a table
    # citation needs would break that parse (checkModelConventions test 1601).
    # Sources: Simpson 2006 Table 5 rows Inter-individual variability for CL/F
    # = 62 %CV (variance 0.3844) and for V/F = 75 %CV (variance 0.5625).
    etalcl ~ 0.62^2
    etalvc ~ 0.75^2

    # ---- Residual error ---------------------------------------------------
    # "The residual error model that gave the minimum value for the objective
    # function was a lognormal residual error model for all the data" (Results).
    # The reported 93% intra-individual variability is again a %CV approximated
    # by the square root of the variance estimate (Methods), i.e. the SD on the
    # log scale.
    expSd <- 0.93
    label("Log-normal residual SD for DHA plasma concentration (SD on the log scale)")  # Simpson 2006 Table 5: "Intra-individual variability" = 93 %CV
  })

  model({
    # 1. Typical (covariate-adjusted) per-kilogram disposition parameters.
    #    Both CL/F and V/F are reported per kilogram of body weight and both
    #    covariate effects are ADDITIVE on that per-kilogram scale, so the
    #    log-transformed ini() terms are the intercepts and the e_* terms are
    #    absolute (not fractional) offsets.
    #      CL/F (L/kg/h) = 2.03 + 1.14 * male
    #      V/F  (L/kg)   = 0.57 + 5.77 * (WT / 70)
    clkg <- (exp(lcl) + e_male_cl * (1 - SEXF)) * exp(etalcl)
    vckg <- (exp(lvc) + e_wt_vc * (WT / 70)) * exp(etalvc)

    # 2. Absolute individual parameters for a dose expressed in mg.
    cl <- clkg * WT
    vc <- vckg * WT

    ka   <- exp(lka)
    tlag <- exp(ltlag)

    # 3. Micro-constants. ka (0.2 /h) is smaller than kel for every patient in
    #    the studied weight range, so the terminal slope is the appearance rate:
    #    "The fixed value of the appearance rate constant was slower than the
    #    estimated elimination rate constant, suggesting that the ARS rectal
    #    capsules are absorption rate limited" (Discussion).
    kel <- cl / vc

    # 4. ODE system. Dose the administered artesunate amount (mg) into depot;
    #    complete conversion to DHA is assumed, with no molar correction, per
    #    Methods ("It was assumed that ARS was converted completely to DHA ...
    #    The exact dose (mg/kg) administered was used in the modelling").
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # 5. Appearance lag.
    alag(depot) <- tlag

    # 6. Observation. central (mg) / vc (L) gives Cc in mg/L, which is the unit
    #    declared in units$concentration above. Simpson 2006 plots DHA in ng/mL,
    #    so a reader comparing against Figure 2 must scale Cc by 1000; the
    #    vignette does that at the plotting step rather than here, so that the
    #    packaged model stays in the library's standard mg/L concentration unit.
    Cc <- central / vc
    Cc ~ lnorm(expSd)
  })
}
