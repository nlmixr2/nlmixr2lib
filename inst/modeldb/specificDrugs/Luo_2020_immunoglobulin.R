Luo_2020_immunoglobulin <- function() {
  description <- "Two-compartment population PK model with first-order subcutaneous absorption for polyclonal immunoglobulin G (IgPro10) in Japanese and non-Japanese patients with primary immunodeficiency (Luo 2020)"
  reference <- "Luo D, Baheti G, Tortorici MA, Hofmann J, Rojavin MA. Pharmacometric analysis of IgPro10 in Japanese and non-Japanese patients with primary immunodeficiency. Clin Ther. 2020;42(1):196-209.e195. doi:10.1016/j.clinthera.2019.11.012 -- parameter values transcribed from the secondary source: van der Zeeuw SL, van Tilburg SJ, Jacobs BC, Koch BCP, Dalm VASH, Crombag MBS, Preijers T. Population pharmacokinetics and pharmacodynamics of immunoglobulins: a systematic review. Clin Pharmacokinet. 2026;65(6):813-30. doi:10.1007/s40262-026-01641-5, Table 4 (reference 41)"
  vignette <- "vanderZeeuw_2026_immunoglobulin"
  units <- list(time = "day", dosing = "g", concentration = "g/L")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric power scaling on CL (estimated exponent 0.881) and on Vc (estimated exponent 0.501), reference weight 58.7 kg. Q and Vp carry no weight term in van der Zeeuw 2026 Table 4.",
      source_name = "BW"
    )
  )

  covariatesDataExcluded <- list(
    RACE_JAPANESE = list(
      description = "Japanese versus non-Japanese ethnicity",
      units = "(binary)",
      type = "binary",
      notes = "Ethnicity (Japanese / non-Japanese) was tested as a covariate -- it is the stated subject of the primary paper -- but was NOT retained in the final model: van der Zeeuw 2026 Table 3 lists 'BW, ethnicity (Japanese/non-Japanese), age, sex' under covariates tested and only 'BW on CL' under covariates in the final model.",
      source_name = "Ethnicity (Japanese / non-Japanese)"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "immunoglobulin G", units = "g", specimen = "administration site", verified = TRUE),
    central = list(analyte = "immunoglobulin G", units = "g", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "immunoglobulin G", units = "g", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 202L,
    n_studies = 7L,
    age_range = "3-81 years",
    age_median = "21 years",
    weight_range = "13-135 kg",
    weight_median = "58.7 kg",
    sex_female_pct = round(100 * 85 / 202, 1),
    race_ethnicity = "Japanese and non-Japanese (proportions not reported in the review)",
    disease_state = "Primary immunodeficiency (PID) on immunoglobulin replacement therapy",
    dose_range = "IVIg 13.3-913 mg/kg every 3 or 4 weeks; SCIg 26.7-379 mg/kg every 2 weeks",
    regions = "United States, Canada, Europe, Asia",
    notes = "Pooled analysis of NCT00419341, NCT00168025, NCT00322556, NCT00542997, NCT02711228 (2016-001631-12), NCT01199705 and NCT01458171 (van der Zeeuw 2026 Table 1). Baseline IgG not reported. The residual-error structure could not be transcribed: see the ini() note and the vignette Errata."
  )

  ini({
    # Structural parameters. van der Zeeuw 2026 Table 4, row 'Luo et al.
    # (2020) [41]'. Reference weight 58.7 kg as printed in the table.
    lcl     <- log(0.139); label("Clearance for a 58.7 kg patient (L/day)")             # van der Zeeuw 2026 Table 4: CL = 0.139 (BW/58.7)^0.881
    lvc     <- log(4.01);  label("Central volume of distribution for a 58.7 kg patient (L)")  # van der Zeeuw 2026 Table 4: Vc = 4.01 (BW/58.7)^0.501
    lq      <- log(0.300); label("Intercompartmental clearance (L/day)")                # van der Zeeuw 2026 Table 4: Q = 0.300
    lvp     <- log(3.51);  label("Peripheral volume of distribution (L)")               # van der Zeeuw 2026 Table 4: Vp = 3.51
    lka     <- log(0.506); label("First-order subcutaneous absorption rate constant (1/day)")  # van der Zeeuw 2026 Table 4: Ka = 0.506 (printed as '0.506 g/day'; see vignette Errata)
    lfdepot <- log(0.668); label("Subcutaneous bioavailability relative to intravenous (fraction)")  # van der Zeeuw 2026 Table 4: F1 = 66.8%

    # Allometric exponents, both estimated (van der Zeeuw 2026 sections
    # 3.2.1.3 and 3.2.1.4 list this model among those with estimated exponents).
    e_wt_cl <- 0.881; label("Allometric exponent on CL (unitless)")                     # van der Zeeuw 2026 Table 4: (BW/58.7)^0.881
    e_wt_vc <- 0.501; label("Allometric exponent on Vc (unitless)")                     # van der Zeeuw 2026 Table 4: (BW/58.7)^0.501

    # Endogenous IgG. van der Zeeuw 2026 section 3.2.1.6: 'In several studies,
    # the endogenous IgG concentration was fixed to 4 g/L [40, 41, 44, 54]' --
    # reference 41 is this model. Held constant, so wrapped in fixed().
    bl_igg  <- fixed(4); label("Endogenous (treatment-naive) IgG concentration (g/L)")  # van der Zeeuw 2026 section 3.2.1.6, held constant at 4 g/L

    # Inter-individual variability (apparent CV%, van der Zeeuw 2026 section
    # 2.3): omega^2 = log(1 + CV^2). This model reports IIV on six parameters,
    # more than any other in the review. The Vp value is extreme but is
    # transcribed as printed.
    etalcl     ~ 0.121991  # 36.02% CV; van der Zeeuw 2026 Table 4 IIV 'CL = 36.02'
    etalq      ~ 0.257975  # 54.25% CV; van der Zeeuw 2026 Table 4 IIV 'Q = 54.25'
    etalvc     ~ 0.614035  # 92.08% CV; van der Zeeuw 2026 Table 4 IIV 'Vc = 92.08'
    etalvp     ~ 2.450015  # 325.4% CV; van der Zeeuw 2026 Table 4 IIV 'Vp = 325.4'
    etalfdepot ~ 0.056002  # 24.0% CV; van der Zeeuw 2026 Table 4 IIV 'F1 = 24.0'
    etalka     ~ 0.995047  # 130.57% CV; van der Zeeuw 2026 Table 4 IIV 'Ka = 130.57'

    # Residual unexplained variability. van der Zeeuw 2026 Table 4 prints
    # 'Unclear' in BOTH residual columns for this row, with footnote b:
    # 'Multiple proportional errors are reported but unclear how they are
    # incorporated'. Set to zero rather than invented; see the vignette Errata.
    propSd <- fixed(0); label("Proportional residual error (fraction; ZERO - not transcribable from source)")  # van der Zeeuw 2026 Table 4: Prop = 'Unclear' (footnote b)
  })
  model({
    cl     <- exp(lcl     + etalcl)     * (WT / 58.7)^e_wt_cl
    vc     <- exp(lvc     + etalvc)     * (WT / 58.7)^e_wt_vc
    q      <- exp(lq      + etalq)
    vp     <- exp(lvp     + etalvp)
    ka     <- exp(lka     + etalka)
    fdepot <- exp(lfdepot + etalfdepot)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # `depot` receives subcutaneous doses; `central` receives intravenous
    # doses directly. States hold EXOGENOUS (therapeutic) IgG only.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    f(depot)          <-  fdepot

    # Observed total plasma IgG = exogenous concentration + endogenous baseline.
    Cc <- central / vc + bl_igg
    Cc ~ prop(propSd)
  })
}
