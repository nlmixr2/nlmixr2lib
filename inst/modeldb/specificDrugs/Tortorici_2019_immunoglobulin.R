Tortorici_2019_immunoglobulin <- function() {
  description <- "Two-compartment population PK model for intravenous polyclonal immunoglobulin G (Privigen) in primary and secondary immunodeficiency, with a disease-type effect on central volume (Tortorici 2019)"
  reference   <- "Tortorici MA, Lawo JP, Weide R, Jochems J, Puli S, Hofmann J, et al. Privigen has similar pharmacokinetic properties in primary and secondary immune deficiency. Int Immunopharmacol. 2019;66:119-26. doi:10.1016/j.intimp.2018.11.013 -- parameter values transcribed from the secondary source: van der Zeeuw SL, van Tilburg SJ, Jacobs BC, Koch BCP, Dalm VASH, Crombag MBS, Preijers T. Population pharmacokinetics and pharmacodynamics of immunoglobulins: a systematic review. Clin Pharmacokinet. 2026;65(6):813-30. doi:10.1007/s40262-026-01641-5, Table 4 (reference 54)"
  vignette    <- "vanderZeeuw_2026_immunoglobulin"
  units       <- list(time = "day", dosing = "g", concentration = "g/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric power scaling on CL (estimated exponent 0.796) and on Vc (exponent 1.1), reference weight 72 kg. Q and Vp carry no weight term in van der Zeeuw 2026 Table 4.",
      source_name        = "BW"
    ),
    DIS_SAD = list(
      description        = "Secondary immunodeficiency indicator (1 = secondary immunodeficiency, 0 = primary immunodeficiency)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (primary immunodeficiency, PID)",
      notes              = "Time-fixed per subject. Tortorici 2019 is the only model in van der Zeeuw 2026 fitted to BOTH PID and SID data, and disease type was retained as a covariate on Vc only (van der Zeeuw 2026 section 3.2.1.6). The source term is 'secondary immunodeficiency (SID)'; the canonical column DIS_SAD (secondary antibody deficiency, reference category primary immunodeficiency) encodes the identical PID-versus-acquired-antibody-deficiency contrast and is reused here rather than minting a near-synonym. The SID cohort in this study was predominantly haematological-malignancy-associated. The authors concluded the Vc difference had no impact on overall IgG exposure (AUC0-28d).",
      source_name        = "Disease type (PID / SID)"
    )
  )

  compartmentData <- list(
    central     = list(analyte = "immunoglobulin G", units = "g", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "immunoglobulin G", units = "g", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 187L,
    n_studies      = 3L,
    age_range      = "PID mean (SD) 29.8 (20.3) years; SID mean (SD) 69.5 (10.4) years",
    age_median     = "Not reported (means given by disease group)",
    weight_range   = "PID mean (SD) 62.6 (26.4) kg; SID mean (SD) 76.8 (16.1) kg",
    weight_median  = "Not reported (means given by disease group)",
    sex_female_pct = round(100 * (41 + 39) / 187, 1),
    race_ethnicity = "Not reported",
    disease_state  = "Primary immunodeficiency (PID, n = 90) and secondary immunodeficiency (SID, n = 97)",
    dose_range     = "IVIg; PID 13.3-959.0 mg/kg, SID 90.9-678.0 mg/kg",
    regions        = "United States, Europe",
    notes          = "Pooled analysis of NCT00168025, NCT00322556 and non-interventional study NIS-Nr 182 (van der Zeeuw 2026 Table 1). Baseline IgG measured 28 +/- 2 days after an IgRT dose and therefore represents endogenous PLUS exogenous IgG: median 9.17 g/L (3.93-27.2) in PID and 6.15 g/L (2.05-17.1) in SID."
  )

  ini({
    # Structural parameters. van der Zeeuw 2026 Table 4, row 'Tortorici et al.
    # (2019) [54]'. Reference weight 72 kg as printed in the table. IVIg only,
    # so no depot, no Ka and no bioavailability term.
    lcl      <- log(0.152); label("Clearance for a 72 kg patient (L/day)")              # van der Zeeuw 2026 Table 4: CL = 0.152 (BW/72)^0.796
    lvc      <- log(3.08);  label("Central volume of distribution for a 72 kg PID patient (L)")  # van der Zeeuw 2026 Table 4: Vc PID = 3.08 (BW/72)^1.1
    lq       <- log(0.825); label("Intercompartmental clearance (L/day)")               # van der Zeeuw 2026 Table 4: Q = 0.825
    lvp      <- log(1.8);   label("Peripheral volume of distribution (L)")              # van der Zeeuw 2026 Table 4: Vp = 1.8

    # Allometric exponents, both estimated (van der Zeeuw 2026 section 3.2.1.4
    # lists 0.796 among the estimated CL exponents; section 3.2.1.3 lists
    # Tortorici 2019 among the five models with an estimated Vc exponent).
    e_wt_cl  <- 0.796; label("Allometric exponent on CL (unitless)")                    # van der Zeeuw 2026 Table 4: (BW/72)^0.796
    e_wt_vc  <- 1.1;   label("Allometric exponent on Vc (unitless)")                    # van der Zeeuw 2026 Table 4: (BW/72)^1.1

    # Disease-type effect on Vc, entered as the SID/PID ratio of the two
    # printed typical values so that both source numbers remain traceable.
    e_sid_vc <- 2.841; label("SID-vs-PID multiplicative ratio on Vc (Vc_SID / Vc_PID)")  # van der Zeeuw 2026 Table 4: Vc PID = 3.08, Vc SID = 8.75; 8.75/3.08 = 2.841

    # Endogenous IgG. van der Zeeuw 2026 section 3.2.1.6: 'In several studies,
    # the endogenous IgG concentration was fixed to 4 g/L [40, 41, 44, 54]' --
    # reference 54 is this model. Held constant, so wrapped in fixed().
    bl_igg   <- fixed(4); label("Endogenous (treatment-naive) IgG concentration (g/L)")  # van der Zeeuw 2026 section 3.2.1.6, held constant at 4 g/L

    # Inter-individual variability (apparent CV%, van der Zeeuw 2026 section
    # 2.3): omega^2 = log(1 + CV^2).
    etalcl ~ 0.248391  # 53.1% CV; van der Zeeuw 2026 Table 4 IIV 'CL = 53.1'
    etalvc ~ 0.450130  # 75.4% CV; van der Zeeuw 2026 Table 4 IIV 'Vc = 75.4'

    # Combined residual error. van der Zeeuw 2026 Table 4 prints the additive
    # term as a bare value in the 'Add (g/L)' column (an SD in g/L) and the
    # proportional term as a percentage.
    addSd  <- 0.930; label("Additive residual error on total IgG (g/L)")                # van der Zeeuw 2026 Table 4: Add = 0.930
    propSd <- 0.071; label("Proportional residual error (fraction)")                    # van der Zeeuw 2026 Table 4: Prop = 7.1%
  })
  model({
    cl <- exp(lcl + etalcl) * (WT / 72)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 72)^e_wt_vc * e_sid_vc^DIS_SAD
    q  <- exp(lq)
    vp <- exp(lvp)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Intravenous administration only: doses go directly to `central`.
    # States hold EXOGENOUS (therapeutic) IgG only.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Observed total plasma IgG = exogenous concentration + endogenous baseline.
    Cc <- central / vc + bl_igg
    Cc ~ add(addSd) + prop(propSd)
  })
}
