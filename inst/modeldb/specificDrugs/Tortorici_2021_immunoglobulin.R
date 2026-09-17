Tortorici_2021_immunoglobulin <- function() {
  description <- "Two-compartment population PK model with first-order subcutaneous absorption for polyclonal immunoglobulin G in chronic inflammatory demyelinating polyneuropathy (Tortorici 2021)"
  reference <- "Tortorici MA, Yuraszeck T, Cornblath D, Bril V, Hartung HP, Sobue G, et al. Pharmacometric analysis linking immunoglobulin exposure to clinical efficacy outcomes in chronic inflammatory demyelinating polyneuropathy. CPT Pharmacometrics Syst Pharmacol. 2021;10(8):839-50. doi:10.1002/psp4.12657 -- parameter values transcribed from the secondary source: van der Zeeuw SL, van Tilburg SJ, Jacobs BC, Koch BCP, Dalm VASH, Crombag MBS, Preijers T. Population pharmacokinetics and pharmacodynamics of immunoglobulins: a systematic review. Clin Pharmacokinet. 2026;65(6):813-30. doi:10.1007/s40262-026-01641-5, Table 4 (reference 43)"
  vignette <- "vanderZeeuw_2026_immunoglobulin"
  units <- list(time = "day", dosing = "g", concentration = "g/L")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric power scaling on all four disposition parameters, reference weight 82 kg: an estimated exponent of 0.615 on CL and Q, and 0.773 on Vc and Vp (van der Zeeuw 2026 Table 4 and sections 3.2.2.3-3.2.2.4).",
      source_name = "BW"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "immunoglobulin G", units = "g", specimen = "administration site", verified = TRUE),
    central = list(analyte = "immunoglobulin G", units = "g", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "immunoglobulin G", units = "g", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 235L,
    n_studies = 2L,
    age_range = "22-83 years",
    age_median = "58 years",
    weight_range = "42.3-133 kg",
    weight_median = "82 kg",
    sex_female_pct = round(100 * 86 / 235, 1),
    race_ethnicity = "Not reported",
    disease_state = "Chronic inflammatory demyelinating polyneuropathy (CIDP)",
    dose_range = "IVIg induction 2000 mg/kg over 2-5 days, followed by maintenance 1000 mg/kg; SCIg 200 or 400 mg/kg weekly",
    regions = "United States, Canada, Europe, Asia, Australia",
    notes = "Analysis of the PATH study (NCT01545076) and doi:10.1111/jns5.12017 (van der Zeeuw 2026 Table 1). Treatment-naive (endogenous) IgG median 12.5 g/L, range 5.6-33.0. Clearance in this CIDP cohort is the highest of any model in the review (0.435 L/day at 82 kg, approximately 0.39 L/day per 70 kg) -- roughly two to eight times the PID values -- which the review attributes to FcRn saturation at immunomodulatory doses (section 4). The companion exposure-response analysis of the INCAT disability score is described in section 3.3.2 but its parameter estimates (Emax, EC50, the probit thresholds) are NOT reported, so only the PK layer is extractable; see the vignette Errata."
  )

  ini({
    # Structural parameters. van der Zeeuw 2026 Table 4, row 'Tortorici et al.
    # (2021) [43]'. Reference weight 82 kg as printed in the table.
    lcl     <- log(0.435); label("Clearance for an 82 kg patient (L/day)")              # van der Zeeuw 2026 Table 4: CL = 0.435 (BW/82)^0.615
    lvc     <- log(4.69);  label("Central volume of distribution for an 82 kg patient (L)")  # van der Zeeuw 2026 Table 4: Vc = 4.69 (BW/82)^0.773
    lq      <- log(0.50);  label("Intercompartmental clearance for an 82 kg patient (L/day)")  # van der Zeeuw 2026 Table 4: Q = 0.50 (BW/82)^0.615
    lvp     <- log(1.87);  label("Peripheral volume of distribution for an 82 kg patient (L)")  # van der Zeeuw 2026 Table 4: Vp = 1.87 (BW/82)^0.773
    lfdepot <- log(0.824); label("Subcutaneous bioavailability relative to intravenous (fraction)")  # van der Zeeuw 2026 Table 4: F1 = 0.824

    # Absorption rate constant -- FIXED, not estimated. van der Zeeuw 2026
    # section 3.2.2.2: 'Tortorici et al. described the absorption of SCIg as a
    # first-order process, with a Ka of 0.439 day^-1 fixed from a previous
    # study because of limited data in the absorption phase [40, 43]'. The
    # prior study is Landersdorfer 2013, whose Table 4 Ka is indeed 0.439.
    lka     <- fixed(log(0.439)); label("First-order subcutaneous absorption rate constant (1/day)")  # van der Zeeuw 2026 section 3.2.2.2: fixed at 0.439 from Landersdorfer 2013

    # Allometric exponents, both estimated (van der Zeeuw 2026 section 3.2.2.4
    # reports the 0.615 clearance exponent in prose; section 3.2.2.3 reports
    # the 0.773 volume exponent).
    e_wt_cl_q  <- 0.615; label("Allometric exponent on CL and Q (unitless)")            # van der Zeeuw 2026 Table 4 and section 3.2.2.4: (BW/82)^0.615
    e_wt_vc_vp <- 0.773; label("Allometric exponent on Vc and Vp (unitless)")           # van der Zeeuw 2026 Table 4 and section 3.2.2.3: (BW/82)^0.773

    # Endogenous IgG. van der Zeeuw 2026 section 3.2.2.5: 'Tortorici et al.
    # modelled endogenous IgG based on observed treatment-naive IgG
    # concentrations, incorporating the residual error [43].' No point estimate
    # is given for a baseline PARAMETER, so the cohort's observed
    # treatment-naive median from Table 1 is used. Wrapped in fixed() because
    # it is an observed summary statistic, not an estimated parameter.
    bl_igg  <- fixed(12.5); label("Endogenous (treatment-naive) IgG concentration (g/L)")  # van der Zeeuw 2026 Table 1: treatment-naive/endogenous IgG median 12.5 g/L (5.6-33.0)

    # Inter-individual variability (apparent CV%, van der Zeeuw 2026 section
    # 2.3): omega^2 = log(1 + CV^2).
    etalcl ~ 0.075478  # 28.0% CV; van der Zeeuw 2026 Table 4 IIV 'CL = 28.0'
    etalvc ~ 0.051548  # 23.0% CV; van der Zeeuw 2026 Table 4 IIV 'Vc = 23.0'

    # Residual error. van der Zeeuw 2026 Table 4 prints '12.1%' in the
    # proportional column and '-' in the additive column.
    propSd <- 0.121; label("Proportional residual error (fraction)")                    # van der Zeeuw 2026 Table 4: Prop = 12.1%
  })
  model({
    cl     <- exp(lcl + etalcl) * (WT / 82)^e_wt_cl_q
    vc     <- exp(lvc + etalvc) * (WT / 82)^e_wt_vc_vp
    q      <- exp(lq)           * (WT / 82)^e_wt_cl_q
    vp     <- exp(lvp)          * (WT / 82)^e_wt_vc_vp
    ka     <- exp(lka)
    fdepot <- exp(lfdepot)

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
