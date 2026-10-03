NgougniPokem_2022_temocillin <- function() {
  description <- paste(
    "Three-compartment intravenous population PK model for UNBOUND temocillin in plasma and ascitic",
    "fluid of critically ill adults with septic shock associated with complicated intra-abdominal",
    "infection and ascitic fluid effusion, given a 2 g loading dose over 30 min followed by a",
    "continuous infusion of 6 g/24 h. Plasma disposition is two-compartment (central plus one",
    "peripheral compartment, parameterised with first-order distribution rate constants); a third",
    "compartment with its own volume represents the ascitic fluid, exchanging with the central",
    "compartment by first-order rate constants and drained by a non-renal clearance (the abdominal",
    "drain). Clearance from the central compartment is proportional to measured urinary creatinine",
    "clearance normalised to the cohort median of 39.9 mL/min. Fitted in Pmetrics with the",
    "non-parametric adaptive grid (NPAG) algorithm; the discrete joint density is approximated",
    "here by independent lognormal marginals centred on the published medians with variances",
    "from the published CV%, so the shape of the joint density is not recoverable from this",
    "encoding. Residual error is the Pmetrics lambda model on the published assay SD polynomial.",
    sep = " "
  )
  reference <- paste(
    "Ngougni Pokem P, Wittebole X, Collienne C, Rodriguez-Villalobos H, Tulkens PM, Elens L,",
    "Van Bambeke F, Laterre PF. Population Pharmacokinetics of Temocillin Administered by",
    "Continuous Infusion in Patients with Septic Shock Associated with Intra-Abdominal Infection",
    "and Ascitic Fluid Effusion. Antibiotics (Basel) 2022;11(7):898.",
    "doi:10.3390/antibiotics11070898. Parameter estimates from Table 4; model structure, covariate",
    "equation and error model from Supplementary Table S1 (Pmetrics model file).",
    sep = " "
  )
  vignette <- "NgougniPokem_2022_temocillin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # The ascitic-fluid compartment is a drained third compartment specific to
  # this paper's patient population; it has no canonical compartment name.
  paper_specific_compartments <- c("ascites")

  compartmentData <- list(
    central = list(analyte = "temocillin, unbound", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "temocillin, unbound", units = "mg", specimen = "plasma", verified = TRUE),
    ascites = list(analyte = "temocillin, unbound", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = paste(
        "Creatinine clearance MEASURED from a urine collection (the paper's CLCRurinary), raw",
        "mL/min, NOT BSA-normalized",
        sep = " "
      ),
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Table 1: median 39.90 mL/min (range 20.55-149.3) [local normal value > 78 mL/min].",
        "Enters the clearance from the central compartment as the through-origin ratio",
        "CL = CLi * (CRCL / 39.9) (Supplementary Table S1 #Sec block,",
        "'Ke = CLi*(CLCRurinary/39.9)/V'; Results section 2.4). There is no non-renal arm on",
        "the central clearance, so the model predicts zero central elimination at zero CRCL; the",
        "paper's simulations span 20-150 mL/min. Measured urinary rather than an estimating",
        "equation, so a Cockcroft-Gault or BSA-normalized value is not interchangeable.",
        sep = " "
      ),
      source_name = "CLCRurinary"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = paste(
        "Table 4 footnote b: individual Bayesian CL30 differed by sex in a post-hoc Mann-Whitney",
        "test (mean 3.04 L/h in females vs 1.38 L/h in males, p = 0.022), but sex was not added",
        "to the population model (Supplementary Table S1 #Cov block lists CLCRurinary only).",
        sep = " "
      )
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Table 1: median 67 kg (range 45-95). Screened (section 4.6.2) and reported with no",
        "statistical influence on the PK parameters (Table 4 footnote b; Results section 2.4).",
        sep = " "
      )
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Table 1: median 56 years (range 21-74). Screened (section 4.6.2) and reported with no",
        "statistical influence on the PK parameters (Table 4 footnote b).",
        sep = " "
      )
    )
  )

  population <- list(
    species = "human",
    n_subjects = 19L,
    n_studies = 1L,
    age_range = "21-74 years",
    age_median = "56 years",
    weight_range = "45-95 kg",
    weight_median = "67 kg",
    sex_female_pct = 68.4,
    race_ethnicity = "Not reported",
    disease_state = paste(
      "Critically ill adults in septic shock associated with complicated intra-abdominal",
      "infection and ascitic fluid effusion (11 spontaneous bacterial peritonitis, all cirrhotic;",
      "4 secondary peritonitis; 3 infected pancreatic necrosis; 1 liver abscess), all with a",
      "positive blood culture. SOFA 9 (4-14), APACHE II 18 (13-32).",
      sep = " "
    ),
    dose_range = paste(
      "Intravenous temocillin 2 g loading dose over 30 min, then continuous infusion of",
      "6 g/24 h; treatment duration 5 days (4-21)",
      sep = " "
    ),
    regions = "Belgium (single centre: Cliniques universitaires Saint-Luc, Brussels)",
    renal_function = paste(
      "Measured urinary creatinine clearance median 39.90 mL/min (range 20.55-149.3)",
      sep = " "
    ),
    albumin = "Plasma albumin median 22.30 g/L (13.70-30.80); ascitic fluid albumin 5.30 g/L (2.12-12.45)",
    n_concentrations = 114L,
    notes = paste(
      "Table 1. 114 unbound temocillin concentrations in plasma and ascitic fluid, sampled",
      "0.5-96 h after the start of treatment; ascitic fluid was collected via the drainage",
      "system. Unbound concentrations were measured by HPLC-MS/MS on ultrafiltrates (30 kDa",
      "cut-off), not computed from total concentrations. Unbound fraction 56.4% in plasma and",
      "57.4% in ascitic fluid; ascitic/plasma total AUC ratio 46.0% (Results section 2.3).",
      sep = " "
    )
  )

  ini({
    # ------------------------------------------------------------------------
    # Structural parameters: Table 4 (final covariate model, unbound
    # temocillin), MEDIAN column. Pmetrics NPAG estimates a discrete joint
    # density over support points; Table 4 summarises each marginal by mean,
    # SD, CV% and median (95% CI). The median is encoded as the median of a
    # lognormal marginal (log(x) is the mu of a lognormal whose median is x).
    # The paper's text quotes neither column, and re-simulating its Tables 5
    # and 6 does not discriminate between them (see vignette Assumptions).
    # The parameterisation follows the Supplementary Table S1 #Pri block:
    # V, CLi, K12, K21, K13, K31, V3, CL30.
    # ------------------------------------------------------------------------
    lvc <- log(13.90)
    label("Volume of the central compartment (L)")
    # Table 4, V row: mean 14.36, SD 4.18, CV 29.15%, median 13.90 (11.76-15.98).
    lcl <- log(2.56)
    label("Clearance from the central compartment at CRCL = 39.9 mL/min (L/h)")
    # Table 4, CLi row: mean 2.45, SD 0.91, CV 37.33%, median 2.56 (2.34-3.71).
    lk12 <- log(4.62)
    label("Rate constant, central to peripheral compartment (1/h)")
    # Table 4, K12 row: mean 4.95, SD 2.89, CV 58.47%, median 4.62 (2.93-6.53).
    lk21 <- log(5.85)
    label("Rate constant, peripheral to central compartment (1/h)")
    # Table 4, K21 row: mean 5.38, SD 3.62, CV 67.20%, median 5.85 (2.08-8.85).
    lk13 <- log(0.24)
    label("Rate constant, central to ascitic fluid compartment (1/h)")
    # Table 4, K13 row: mean 0.42, SD 0.47, CV 110.67%, median 0.24 (0.16-0.41).
    lk31 <- log(0.15)
    label("Rate constant, ascitic fluid to central compartment (1/h)")
    # Table 4, K31 row: mean 0.33, SD 0.46, CV 137.94%, median 0.15 (0.04-0.26).
    lv_ascites <- log(28.93)
    label("Volume of the ascitic fluid compartment (L)")
    # Table 4, V3 row: mean 30.00, SD 16.76, CV 55.87%, median 28.93 (15.83-42.77).
    lcl_ascites <- log(2.94)
    label("Non-renal (drain) clearance from the ascitic fluid compartment (L/h)")
    # Table 4, CL30 row: mean 2.91, SD 1.42, CV 48.75%, median 2.94 (2.11-3.60).

    # ------------------------------------------------------------------------
    # Between-subject variability. Table 4 CV% is the descriptive SD/mean of
    # the NPAG marginal (e.g. V: 4.18 / 14.36 = 29.1%). A lognormal with
    # omega^2 = log(CV^2 + 1) reproduces each CV% exactly. No correlations or
    # covariance matrix are published, so the etas are independent.
    # ------------------------------------------------------------------------
    # Table 4, V CV 29.15%: log(0.2915^2 + 1) = 0.08155441
    etalvc ~ 0.08155441
    # Table 4, CLi CV 37.33%: log(0.3733^2 + 1) = 0.1304605
    etalcl ~ 0.1304605
    # Table 4, K12 CV 58.47%: log(0.5847^2 + 1) = 0.2940672
    etalk12 ~ 0.2940672
    # Table 4, K21 CV 67.20%: log(0.6720^2 + 1) = 0.3726554
    etalk21 ~ 0.3726554
    # Table 4, K13 CV 110.67%: log(1.1067^2 + 1) = 0.7996602
    etalk13 ~ 0.7996602
    # Table 4, K31 CV 137.94%: log(1.3794^2 + 1) = 1.065657
    etalk31 ~ 1.065657
    # Table 4, V3 CV 55.87%: log(0.5587^2 + 1) = 0.2716637
    etalv_ascites ~ 0.2716637
    # Table 4, CL30 CV 48.75%: log(0.4875^2 + 1) = 0.2132195
    etalcl_ascites ~ 0.2132195

    # ------------------------------------------------------------------------
    # Residual error: Pmetrics lambda model. Supplementary Table S1 #Err block:
    #   L=2.26
    #   0.1,0.1,0,0   (plasma)
    #   0.1,0.1,0,0   (ascitic fluid)
    # i.e. assay SD = C0 + C1*Y with C0 = 0.1 mg/L and C1 = 0.1 for both
    # outputs, and lambda = 2.26 mg/L. Results section 2.4: 'The final Lambda
    # (L) error factor was set at 2.26 ... SD = 0.1 + 0.1Y ... for both plasma
    # and ascitic fluid.' Each observation is weighted by 1/Error^2 with
    # Error = (SD^2 + lambda^2)^0.5 (the Pmetrics additive lambda model; the
    # typeset formula lost the exponent on SD). All three values are stated
    # constants of the final model, hence fixed().
    # ------------------------------------------------------------------------
    addSd <- fixed(0.1)
    label("Assay-error polynomial intercept C0, plasma (mg/L)")
    # Supplementary Table S1 #Err, first polynomial row, C0 = 0.1.
    propSd <- fixed(0.1)
    label("Assay-error polynomial slope C1, plasma (fraction)")
    # Supplementary Table S1 #Err, first polynomial row, C1 = 0.1.
    addSd_Cascites <- fixed(0.1)
    label("Assay-error polynomial intercept C0, ascitic fluid (mg/L)")
    # Supplementary Table S1 #Err, second polynomial row, C0 = 0.1.
    propSd_Cascites <- fixed(0.1)
    label("Assay-error polynomial slope C1, ascitic fluid (fraction)")
    # Supplementary Table S1 #Err, second polynomial row, C1 = 0.1.
    lambdaSd <- fixed(2.26)
    label("Pmetrics lambda additive process noise, both outputs (mg/L)")
    # Supplementary Table S1 #Err, L = 2.26; Results section 2.4.
  })
  model({
    # 1. Individual parameters. Supplementary Table S1 #Sec block:
    #      Ke  = CLi*(CLCRurinary/39.9)/V
    #      K30 = CL30/V3
    #    so the central clearance is CLi scaled through the origin by measured
    #    urinary creatinine clearance relative to the cohort median 39.9 mL/min.
    vc <- exp(lvc + etalvc)
    cl <- exp(lcl + etalcl) * (CRCL / 39.9)
    k12 <- exp(lk12 + etalk12)
    k21 <- exp(lk21 + etalk21)
    k13 <- exp(lk13 + etalk13)
    k31 <- exp(lk31 + etalk31)
    v_ascites <- exp(lv_ascites + etalv_ascites)
    cl_ascites <- exp(lcl_ascites + etalcl_ascites)

    # 2. Micro-constants.
    kel <- cl / vc
    k30 <- cl_ascites / v_ascites

    # 3. ODE system. Supplementary Table S1 #Dif block:
    #      XP(1) = RATEIV(1) - (K12+K13+Ke)*X(1) + K21*X(2) + K31*X(3)
    #      XP(2) = K12*X(1) - K21*X(2)
    #      XP(3) = K13*X(1) - K31*X(3) - K30*X(3)
    #    X(1) = central, X(2) = peripheral1, X(3) = ascites. Doses are
    #    intravenous infusions into the central compartment.
    d/dt(central) <- -(k12 + k13 + kel) * central + k21 * peripheral1 + k31 * ascites
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(ascites) <- k13 * central - k31 * ascites - k30 * ascites

    # 4. Observation and error. Supplementary Table S1 #Out block:
    #      Y(1) = X(1)/V    (unbound plasma, mg/L)
    #      Y(2) = X(3)/V3   (unbound ascitic fluid, mg/L)
    #    Pmetrics lambda error: SD = sqrt((C0 + C1*Y)^2 + lambda^2).
    Cc <- central / vc
    Cascites <- ascites / v_ascites
    sdCc <- sqrt((addSd + propSd * Cc)^2 + lambdaSd^2)
    sdCascites <- sqrt((addSd_Cascites + propSd_Cascites * Cascites)^2 + lambdaSd^2)
    Cc ~ add(sdCc)
    Cascites ~ add(sdCascites)
  })
}
