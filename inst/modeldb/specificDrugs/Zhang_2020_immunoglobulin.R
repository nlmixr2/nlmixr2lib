Zhang_2020_immunoglobulin <- function() {
  description <- "Two-compartment population PK model with first-order subcutaneous absorption for polyclonal immunoglobulin G (IgPro20, Hizentra) given weekly or biweekly in primary immunodeficiency (Zhang 2020)"
  reference <- "Zhang Y, Baheti G, Chapdelaine H, Hofmann J, Rojavin M, Tortorici M, et al. Population pharmacokinetic analysis of weekly and biweekly IgPro20 (Hizentra) dosing in patients with primary immunodeficiency. Int Immunopharmacol. 2020;81:106005. doi:10.1016/j.intimp.2019.106005 -- parameter values transcribed from the secondary source: van der Zeeuw SL, van Tilburg SJ, Jacobs BC, Koch BCP, Dalm VASH, Crombag MBS, Preijers T. Population pharmacokinetics and pharmacodynamics of immunoglobulins: a systematic review. Clin Pharmacokinet. 2026;65(6):813-30. doi:10.1007/s40262-026-01641-5, Table 4 (reference 44)"
  vignette <- "vanderZeeuw_2026_immunoglobulin"
  units <- list(time = "day", dosing = "g", concentration = "g/L")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric power scaling on CL (estimated exponent 0.768) and on Vc (estimated exponent 0.448), reference weight 66 kg. Q and Vp carry no weight term in van der Zeeuw 2026 Table 4.",
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
    n_subjects = 173L,
    n_studies = 5L,
    age_range = "3.0-81.0 years (per-study medians 18.0-32.0 years)",
    age_median = "NCT00419341: 32.0; NCT00168025: 25.0; NCT00322556: 23.0; NCT00542997: 18.0; NCT02711228: 19.0 years",
    weight_range = "13.0-135.0 kg",
    weight_median = "NCT00419341: 66.0; NCT00168025: 66.5; NCT00322556: 62.0; NCT00542997: 54.0; NCT02711228: 72.5 kg",
    sex_female_pct = round(100 * (27 + 34 + 39 + 16 + 9) / (22 + 27 + 46 + 34 + 26 + 39 + 35 + 16 + 8 + 9), 1),
    race_ethnicity = "Not reported",
    disease_state = "Primary immunodeficiency (PID) on immunoglobulin replacement therapy",
    dose_range = "IVIg 200-888 mg/kg every 3 or 4 weeks; SCIg 117.0-120.7 mg/kg weekly and 179.6-224.3 mg/kg every 2 weeks",
    regions = "United States, Canada, Europe",
    notes = "Pooled analysis of NCT00419341, NCT00168025, NCT00322556, NCT00542997 and NCT02711228 (van der Zeeuw 2026 Tables 1 and 2). Baseline IgG not reported."
  )

  ini({
    # Structural parameters. van der Zeeuw 2026 Table 4, row 'Zhang et al.
    # (2020) [44]'. Reference weight 66 kg as printed in the table.
    lcl     <- log(0.138); label("Clearance for a 66 kg patient (L/day)")               # van der Zeeuw 2026 Table 4: CL = 0.138 (BW/66)^0.768
    lvc     <- log(3.95);  label("Central volume of distribution for a 66 kg patient (L)")  # van der Zeeuw 2026 Table 4: Vc = 3.95 (BW/66)^0.448
    lq      <- log(0.260); label("Intercompartmental clearance (L/day)")                # van der Zeeuw 2026 Table 4: Q = 0.260
    lvp     <- log(4.44);  label("Peripheral volume of distribution (L)")               # van der Zeeuw 2026 Table 4: Vp = 4.44
    lka     <- log(0.444); label("First-order subcutaneous absorption rate constant (1/day)")  # van der Zeeuw 2026 Table 4: Ka = 0.444
    lfdepot <- log(0.676); label("Subcutaneous bioavailability relative to intravenous (fraction)")  # van der Zeeuw 2026 Table 4: F1 = 67.6%

    # Allometric exponents, both estimated (van der Zeeuw 2026 sections
    # 3.2.1.3 and 3.2.1.4 list this model among those with estimated exponents).
    e_wt_cl <- 0.768; label("Allometric exponent on CL (unitless)")                     # van der Zeeuw 2026 Table 4: (BW/66)^0.768
    e_wt_vc <- 0.448; label("Allometric exponent on Vc (unitless)")                     # van der Zeeuw 2026 Table 4: (BW/66)^0.448

    # Endogenous IgG. van der Zeeuw 2026 section 3.2.1.6: 'In several studies,
    # the endogenous IgG concentration was fixed to 4 g/L [40, 41, 44, 54]' --
    # reference 44 is this model. Held constant, so wrapped in fixed().
    bl_igg  <- fixed(4); label("Endogenous (treatment-naive) IgG concentration (g/L)")  # van der Zeeuw 2026 section 3.2.1.6, held constant at 4 g/L

    # Inter-individual variability (apparent CV%, van der Zeeuw 2026 section
    # 2.3): omega^2 = log(1 + CV^2).
    etalcl     ~ 0.121991  # 36.02% CV; van der Zeeuw 2026 Table 4 IIV 'CL = 36.02'
    etalq      ~ 0.320975  # 61.52% CV; van der Zeeuw 2026 Table 4 IIV 'Q = 61.52'
    etalvc     ~ 0.619019  # 92.58% CV; van der Zeeuw 2026 Table 4 IIV 'Vc = 92.58'
    etalvp     ~ 2.000029  # 252.77% CV; van der Zeeuw 2026 Table 4 IIV 'Vp = 252.77'
    etalfdepot ~ 0.099589  # 32.36% CV; van der Zeeuw 2026 Table 4 IIV 'F1 = 32.36'
    etalka     ~ 0.263020  # 54.85% CV; van der Zeeuw 2026 Table 4 IIV 'Ka = 54.85'

    # Combined residual error. van der Zeeuw 2026 Table 4 prints both terms as
    # VARIANCES ('sigma^2 = ...'), so each is entered here as its square root.
    # The review distinguishes sigma from sigma^2 elsewhere in the same table
    # (Fokkink 2022 is printed as 'sigma = 0.12'), so the squared notation is
    # taken at face value. See the vignette Errata.
    addSd  <- 0.079498; label("Additive residual error on total IgG (g/L)")             # van der Zeeuw 2026 Table 4: Add sigma^2 = 0.00632; sqrt = 0.079498
    propSd <- 0.100499; label("Proportional residual error (fraction)")                 # van der Zeeuw 2026 Table 4: Prop sigma^2 = 0.0101; sqrt = 0.100499
  })
  model({
    cl     <- exp(lcl     + etalcl)     * (WT / 66)^e_wt_cl
    vc     <- exp(lvc     + etalvc)     * (WT / 66)^e_wt_vc
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
    Cc ~ add(addSd) + prop(propSd)
  })
}
