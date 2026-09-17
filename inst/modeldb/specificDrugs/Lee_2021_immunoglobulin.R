Lee_2021_immunoglobulin <- function() {
  description <- "One-compartment population PK model for intravenous polyclonal immunoglobulin G in patients with predominantly antibody deficiencies, fitted to endogenous-subtracted (exogenous) IgG concentrations (Lee 2021)"
  reference <- "Lee JL, Mohd Saffian S, Makmor-Bakry M, Islahudin F, Alias H, Noh LM, et al. Population pharmacokinetic modelling of intravenous immunoglobulin in patients with predominantly antibody deficiencies. Br J Clin Pharmacol. 2021;87(7):2956-66. doi:10.1111/bcp.14710 -- parameter values transcribed from the secondary source: van der Zeeuw SL, van Tilburg SJ, Jacobs BC, Koch BCP, Dalm VASH, Crombag MBS, Preijers T. Population pharmacokinetics and pharmacodynamics of immunoglobulins: a systematic review. Clin Pharmacokinet. 2026;65(6):813-30. doi:10.1007/s40262-026-01641-5, Table 4 (reference 36)"
  vignette <- "vanderZeeuw_2026_immunoglobulin"
  units <- list(time = "day", dosing = "g", concentration = "g/L")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric power scaling on CL (estimated exponent 0.88) and on Vc (estimated exponent 0.66), reference weight 27 kg -- by far the lowest reference weight among the models in van der Zeeuw 2026, reflecting a small, predominantly paediatric Asian cohort (median 26.7 kg).",
      source_name = "BW"
    )
  )

  covariatesDataExcluded <- list(
    SNP_FCGRT_VNTR = list(
      description = "FcRn (FCGRT) variable-number-tandem-repeat genotype (VNTR2/3 heterozygous, VNTR3/3 homozygous)",
      units = "(categorical)",
      type = "categorical",
      notes = "Screened but NOT retained. van der Zeeuw 2026 section 3.2.1.6: the heterozygous polymorphism (VNTR2/3) was observed in four patients and the homozygous polymorphism (VNTR3/3) in six; 'Ultimately, FcRn polymorphisms were not identified as a significant covariate.' The review notes the small sample size may have driven this result. Documented here to preserve the provenance of the covariate screen; no point estimate is available.",
      source_name = "VNTR FcRn polymorphism"
    )
  )

  compartmentData <- list(
    central = list(analyte = "immunoglobulin G", units = "g", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 10L,
    n_studies = 1L,
    age_range = "3-64 years",
    age_median = "9.5 years",
    weight_range = "9.3-75 kg",
    weight_median = "26.7 kg",
    sex_female_pct = round(100 * 1 / 10, 1),
    race_ethnicity = "Asian (Malaysian cohort)",
    disease_state = "Predominantly antibody deficiency (a form of primary immunodeficiency); the cohort was mainly X-linked agammaglobulinaemia (XLA) patients",
    dose_range = "IVIg 360-600 mg/kg every 4 weeks",
    regions = "Asia",
    notes = "Single prospective observational study (van der Zeeuw 2026 Table 1). Treatment-naive (endogenous) IgG median 0.7 g/L, range 0.3-5.09 g/L -- the lowest endogenous concentrations of any cohort in the review, consistent with XLA. A one-compartment model was selected by Bayesian information criterion; van der Zeeuw 2026 section 3.2.1.3 attributes the lack of benefit from a two-compartment model to the small sample size and limited sampling spread."
  )

  ini({
    # Structural parameters. van der Zeeuw 2026 Table 4, row 'Lee et al.
    # (2021) [36]'. Reference weight 27 kg as printed in the table. IVIg only,
    # so no depot, Ka or bioavailability term.
    lcl     <- log(0.0624); label("Clearance for a 27 kg patient (L/day)")              # van der Zeeuw 2026 Table 4: CL = 0.0624 (BW/27 kg)^0.88
    lvc     <- log(2.77);   label("Central volume of distribution for a 27 kg patient (L)")  # van der Zeeuw 2026 Table 4: Vc = 2.77 (weight/27 kg)^0.66

    # Allometric exponents, both estimated (van der Zeeuw 2026 section 3.2.1.3
    # reports the Vc exponent of 0.66 in prose as well as in Table 4).
    e_wt_cl <- 0.88; label("Allometric exponent on CL (unitless)")                      # van der Zeeuw 2026 Table 4: (BW/27 kg)^0.88
    e_wt_vc <- 0.66; label("Allometric exponent on Vc (unitless)")                      # van der Zeeuw 2026 Table 4 and section 3.2.1.3: (weight/27 kg)^0.66

    # Endogenous IgG. Unlike every other model in the review, Lee 2021 removed
    # endogenous IgG from the DATA rather than modelling it: van der Zeeuw 2026
    # section 3.2.1.6, 'Lee et al. subtracted the measured endogenous IgG
    # concentrations prior to treatment from the total IgG concentration [36].
    # If pre-treatment data were unavailable, an endogenous IgG concentration
    # of 0.7 g/L was used.' The dependent variable is therefore EXOGENOUS IgG,
    # and the correct baseline offset inside the model is exactly zero. Users
    # wanting total IgG must add their own endogenous concentration.
    bl_igg  <- fixed(0); label("Endogenous IgG offset (g/L; ZERO - subtracted from the data before fitting)")  # van der Zeeuw 2026 section 3.2.1.6: endogenous IgG subtracted from the observations

    # Inter-individual variability (apparent CV%, van der Zeeuw 2026 section
    # 2.3): omega^2 = log(1 + CV^2). The Vc value is confirmed twice in the
    # source -- Table 4 and the prose of section 3.2.1.3.
    etalcl ~ 0.019328  # 13.97% CV; van der Zeeuw 2026 Table 4 IIV 'CL = 13.97'
    etalvc ~ 0.006380  # 8.00% CV; van der Zeeuw 2026 Table 4 and section 3.2.1.3 IIV 'Vc = 8.00'

    # Residual error. van der Zeeuw 2026 Table 4 prints '10.23%' in the
    # proportional column and '-' in the additive column.
    propSd <- 0.1023; label("Proportional residual error (fraction)")                   # van der Zeeuw 2026 Table 4: Prop = 10.23%
  })
  model({
    cl <- exp(lcl + etalcl) * (WT / 27)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 27)^e_wt_vc

    kel <- cl / vc

    # Intravenous administration only: doses go directly to `central`.
    d/dt(central) <- -kel * central

    # The observation is EXOGENOUS IgG (endogenous was subtracted from the
    # data); bl_igg is fixed at zero so this reduces to central / vc, while
    # keeping the same observation expression as its sibling models.
    Cc <- central / vc + bl_igg
    Cc ~ prop(propSd)
  })
}
