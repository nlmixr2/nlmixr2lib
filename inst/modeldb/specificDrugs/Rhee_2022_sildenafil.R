Rhee_2022_sildenafil <- function() {
  description <- paste(
    "Joint parent + metabolite population PK model for oral sildenafil and",
    "its active metabolite N-desmethyl sildenafil (DMS) in 19 term and",
    "preterm infants with pulmonary arterial hypertension (Rhee 2022).",
    "Sildenafil is described by a one-compartment disposition with",
    "first-order absorption; its whole apparent clearance is assumed to form",
    "DMS (complete conversion, molar basis), which has its own",
    "one-compartment disposition. Current body weight enters both apparent",
    "clearances through power functions referenced to 3.14 kg. Correlated",
    "IIV on sildenafil V/F, sildenafil CL/F and DMS CL/F'; separate",
    "log-scale additive (exponential) residual errors per analyte."
  )
  reference <- paste(
    "Rhee SJ, Shin SH, Oh J, Jung YH, Choi CW, Kim HS, Yu KS.",
    "Population pharmacokinetic analysis of sildenafil in term and preterm",
    "infants with pulmonary arterial hypertension.",
    "Sci Rep. 2022;12:7393. doi:10.1038/s41598-022-11038-6."
  )
  vignette <- "Rhee_2022_sildenafil"
  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "ng/mL"
  )

  # The paper fitted the model on a molar scale (dose in umol, concentrations
  # in nmol/L; Methods 'Population pharmacokinetic analysis', molecular
  # weights 474.6 g/mol sildenafil and 460.6 g/mol DMS). Here dosing stays in
  # mg of sildenafil and each state holds the mass of its own analyte: the
  # sildenafil eliminated (mg) is converted to the mass of DMS formed by the
  # molecular-weight ratio, which is algebraically identical to the paper's
  # molar mass balance and returns both concentrations in ng/mL.
  compartmentData <- list(
    depot = list(analyte = "sildenafil", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "sildenafil", units = "mg", specimen = "plasma", verified = TRUE),
    central_ndmsil = list(analyte = "N-desmethyl sildenafil", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Current body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power effect on sildenafil CL/F (exponent 0.899) and on DMS CL/F'",
        "(exponent 1.34), both referenced to 3.14 kg as printed in the",
        "Table 2 covariate equations (the Table 1 cohort median is 3.18 kg).",
        "Observed range 0.79-4.09 kg. The paper does not say whether weight",
        "was time-varying over the up-to-6 samples per infant; supply the",
        "current weight at each record."
      ),
      source_name = "body weight"
    )
  )

  # Screened in the stepwise covariate search (paper Methods 'Population
  # pharmacokinetic analysis'; Results 'Final Population Pharmacokinetic
  # Model'; Discussion) but not retained in the final model.
  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Biological sex indicator, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Full model (Table 2): CL/F multiplied by theta_female = 1.78 for",
        "females. Statistically significant but excluded by the authors",
        "because there was no clear pharmacological explanation."
      )
    ),
    PNA = list(
      description = "Postnatal age",
      units = "months",
      type = "continuous",
      notes = paste(
        "Full model (Table 2): exponential-asymptotic maturation on",
        "sildenafil CL/F (theta_maturation = 4.28). Removed in backward",
        "elimination. The paper records postnatal age in days."
      )
    ),
    GA = list(
      description = "Gestational age at birth",
      units = "weeks",
      type = "continuous",
      notes = "Screened; not significant (Discussion)."
    ),
    PAGE = list(
      description = "Postmenstrual age",
      units = "months",
      type = "continuous",
      notes = paste(
        "Significant but strongly correlated with body weight (Pearson",
        "r = 0.803); body weight was retained instead (Discussion). The",
        "paper records postmenstrual age in weeks."
      )
    ),
    WT_BIRTH = list(
      description = "Body weight at birth",
      units = "kg",
      type = "continuous",
      notes = "Screened; not retained (Discussion). The paper records birth weight in g."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "IU/L",
      type = "continuous",
      notes = "Screened with a power function; not retained."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "IU/L",
      type = "continuous",
      notes = "Screened with a power function; not retained."
    ),
    BUN = list(
      description = "Blood urea nitrogen",
      units = "mg/dL",
      type = "continuous",
      notes = "Screened with a power function; not retained."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "mg/dL",
      type = "continuous",
      notes = "Screened with a power function; not retained."
    ),
    CONMED_CYP3A4_IND = list(
      description = "Concomitant CYP3A inducer (bosentan, phenobarbital)",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Four infants; no significant effect found (Discussion, limitations)."
      )
    ),
    CONMED_CYP3A4_INH = list(
      description = "Concomitant CYP inhibitor (fluconazole, ranitidine)",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Three infants; no significant effect found (Discussion, limitations)."
      )
    )
  )

  population <- list(
    species = "human",
    n_subjects = 19,
    n_studies = 1,
    n_observations = "99 PK samples (sildenafil and DMS each measured)",
    age_range = "postnatal age 5-98 days (median 11)",
    ga_range = "gestational age 24-41 weeks (median 36)",
    weight_range = "0.79-4.09 kg (median 3.18)",
    birth_weight_range = "590-4270 g (median 2340)",
    sex_female_pct = 47.4,
    race_ethnicity = "Korean (single country; not tabulated)",
    disease_state = paste(
      "Term and preterm neonates with pulmonary arterial hypertension in",
      "neonatal intensive care: PPHN 7, congenital heart disease 4,",
      "bronchopulmonary dysplasia 6, other 2."
    ),
    dose_range = paste(
      "Oral sildenafil 0.5 mg/kg four times daily for 2 weeks, increased to",
      "0.75 mg/kg four times daily when there was no echocardiographic",
      "improvement."
    ),
    regions = "Republic of Korea (Seoul National University Hospital and Bundang Hospital)",
    co_medication = paste(
      "Bosentan 3, fluconazole 2, phenobarbital 1, ranitidine 1 (Table 1)."
    ),
    notes = paste(
      "Demographics from paper Table 1. Open-label prospective study",
      "NCT02244528 (Feb 2015 - Jul 2016). Up to 6 opportunistic samples per",
      "infant, all after at least 4 doses. LC-MS/MS assay, LLOQ 1 ng/mL",
      "sildenafil and 0.5 ng/mL DMS. NONMEM 7.4, FOCE-I."
    )
  )

  ini({
    # Sildenafil -- paper Table 2, 'Final model' column.
    lka <- log(0.414)
    label("Sildenafil absorption rate constant (1/h)")
    # Table 2: KA = 0.414 1/h (RSE 10.7%).

    lvc <- log(19.8)
    label("Sildenafil apparent volume of distribution V/F (L)")
    # Table 2: V_Sil/F = 19.8 L (RSE 31.4%). Not weight-scaled.

    lcl <- log(10.1)
    label("Sildenafil apparent clearance CL/F at 3.14 kg (L/h)")
    # Table 2: theta_CL(Sil) = 10.1 L/h (RSE 24.7%).

    e_wt_cl <- 0.899
    label("Power exponent of body weight on sildenafil CL/F (unitless)")
    # Table 2: theta_weight (sildenafil) = 0.899 (RSE 7.0%).

    # N-desmethyl sildenafil -- paper Table 2, 'Final model' column.
    lvc_ndmsil <- log(1.78)
    label("DMS apparent volume of distribution V/F' (L)")
    # Table 2: V_DMS/F' = 1.78 L (RSE 11.2%). Not weight-scaled.

    lcl_ndmsil <- log(14.3)
    label("DMS apparent clearance CL/F' at 3.14 kg (L/h)")
    # Table 2: theta_CL(DMS) = 14.3 L/h (RSE 13.2%).

    e_wt_cl_ndmsil <- 1.34
    label("Power exponent of body weight on DMS CL/F' (unitless)")
    # Table 2: theta_weight (DMS) = 1.34 (RSE 26.1%).

    # IIV -- Table 2 reports %CV for exponential etas and the eta
    # correlations. Variances are omega^2 = log(1 + CV^2): 143.7% (V/F),
    # 63.5% (CL/F), 49.2% (DMS CL/F'). Covariances are r * omega_i * omega_j
    # with r = 0.295 (V, CL), 0.178 (CL, DMS CL), 0.0606 (V, DMS CL).
    # Lower-triangle order: var(cl); cov(vc,cl), var(vc); cov(cl_ndmsil,cl),
    # cov(cl_ndmsil,vc), var(cl_ndmsil).
    etalcl + etalvc + etalcl_ndmsil ~ c(
      0.338773,
      0.181716, 1.120040,
      0.048237, 0.029860, 0.216775
    )

    # Residual error -- additive on log-transformed concentrations
    # (Methods), i.e. exponential in linear space.
    expSd <- 0.553
    label("Log-scale residual SD for sildenafil (unitless)")
    # Table 2: residual variability (SD) sildenafil = 0.553 (RSE 7.7%).

    expSd_ndmsil <- 0.472
    label("Log-scale residual SD for DMS (unitless)")
    # Table 2: residual variability (SD) N-desmethyl sildenafil = 0.472 (RSE 15.9%).
  })

  model({
    # Table 2 covariate equations: CL = theta * (body weight / 3.14)^theta_weight
    wt_ref <- 3.14

    # Molecular weights (g/mol), paper Methods 'Population pharmacokinetic
    # analysis'; used only to convert the sildenafil mass eliminated into
    # the mass of DMS formed (the paper's molar-scale complete conversion).
    mw <- 474.6
    mw_ndmsil <- 460.6

    ka <- exp(lka)
    vc <- exp(lvc + etalvc)
    cl <- exp(lcl + etalcl) * (WT / wt_ref)^e_wt_cl

    vc_ndmsil <- exp(lvc_ndmsil)
    cl_ndmsil <- exp(lcl_ndmsil + etalcl_ndmsil) * (WT / wt_ref)^e_wt_cl_ndmsil

    kel <- cl / vc
    kel_ndmsil <- cl_ndmsil / vc_ndmsil

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    # Sildenafil is assumed to be completely metabolised to DMS (Methods).
    d/dt(central_ndmsil) <- kel * central * mw_ndmsil / mw - kel_ndmsil * central_ndmsil

    # mg/L -> ng/mL
    Cc <- 1000 * central / vc
    Cc_ndmsil <- 1000 * central_ndmsil / vc_ndmsil

    Cc ~ lnorm(expSd)
    Cc_ndmsil ~ lnorm(expSd_ndmsil)
  })
}
