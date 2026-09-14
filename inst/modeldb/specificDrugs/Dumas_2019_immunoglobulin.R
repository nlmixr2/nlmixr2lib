Dumas_2019_immunoglobulin <- function() {
  description <- "One-compartment population PK model with first-order subcutaneous absorption for polyclonal immunoglobulin G (IVIg and SCIg 20%, Ig20Gly) in primary immunodeficiency (Dumas 2019)"
  reference   <- "Dumas T, Berry NS, Wolfsegger M, Jolles S, McCoy B, Yel L. Population pharmacokinetic modeling and simulation of immunoglobulin exposure with varying dosing intervals of subcutaneous immunoglobulin 20% (Ig20Gly) in patients with primary immunodeficiency diseases. Int Immunopharmacol. 2019;71:404-10. doi:10.1016/j.intimp.2019.03.043 -- parameter values transcribed from the secondary source: van der Zeeuw SL, van Tilburg SJ, Jacobs BC, Koch BCP, Dalm VASH, Crombag MBS, Preijers T. Population pharmacokinetics and pharmacodynamics of immunoglobulins: a systematic review. Clin Pharmacokinet. 2026;65(6):813-30. doi:10.1007/s40262-026-01641-5, Table 4 (reference 53)"
  vignette    <- "vanderZeeuw_2026_immunoglobulin"
  units       <- list(time = "day", dosing = "g", concentration = "g/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric power scaling on CL only, with reference weight 70 kg and an estimated exponent of 0.576. Dumas 2019 is listed among the 'BW on CL' models in van der Zeeuw 2026 section 3.2.1.4 but NOT among the 'BW on Vc' models in section 3.2.1.3, and its Table 4 Vc entry (4.01 L) carries no allometric term.",
      source_name        = "BW"
    )
  )

  compartmentData <- list(
    depot   = list(analyte = "immunoglobulin G", units = "g", specimen = "administration site", verified = TRUE),
    central = list(analyte = "immunoglobulin G", units = "g", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 102L,
    n_studies      = 2L,
    age_range      = "2.0-83 years",
    age_median     = "30.0 years",
    weight_range   = "13.2-161.8 kg",
    weight_median  = "63.7 kg",
    sex_female_pct = round(100 * 45 / 102, 1),
    race_ethnicity = "Not reported",
    disease_state  = "Primary immunodeficiency (PID) on immunoglobulin replacement therapy",
    dose_range     = "IVIg 300-1000 mg/kg every 3 or 4 weeks; SCIg (Ig20Gly) weekly-equivalent dose",
    regions        = "United States, Canada, Europe",
    notes          = "Pooled analysis of NCT01218438 and NCT01412385 (van der Zeeuw 2026 Table 1). Median baseline total IgG 9.9 g/L (range 3.4-19.9); patients were already on stable IVIg treatment, so this baseline represents endogenous PLUS exogenous IgG (van der Zeeuw 2026 section 3.2.1.6). A one-compartment model was used because samples were collected only at steady state (van der Zeeuw 2026 section 3.2.1.3)."
  )

  ini({
    # Structural parameters. van der Zeeuw 2026 Table 4, row 'Dumas et al.
    # (2019) [53]'. Reference weight 70 kg as printed in the table.
    lcl     <- log(0.09216); label("Clearance for a 70 kg patient (L/day)")             # van der Zeeuw 2026 Table 4: CL = 0.09216 (BW/70 kg)^0.576
    lvc     <- log(4.01);    label("Central volume of distribution (L)")                # van der Zeeuw 2026 Table 4: Vc = 4.01
    lka     <- log(0.096);   label("First-order subcutaneous absorption rate constant (1/day)")  # van der Zeeuw 2026 Table 4: Ka = 0.096
    lfdepot <- log(0.739);   label("Subcutaneous bioavailability relative to intravenous (fraction)")  # van der Zeeuw 2026 Table 4: F1 = 73.9%

    # Allometric exponent on CL -- estimated, not fixed (van der Zeeuw 2026
    # section 3.2.1.4 lists 0.576 among the six ESTIMATED exponents).
    e_wt_cl <- 0.576; label("Allometric exponent on CL (unitless)")                     # van der Zeeuw 2026 Table 4: (BW/70 kg)^0.576

    # Endogenous IgG. van der Zeeuw 2026 section 3.2.1.6: 'Two studies did not
    # report how they incorporated endogenous IgG concentrations [38, 53]' --
    # reference 53 is this model. The value below is NOT from the primary; it
    # is the review's own standardised simulation assumption for PID
    # (van der Zeeuw 2026 section 2.4: 'Endogenous IgG concentrations ... were
    # assumed to be ... 4.0 g/L for PID'), adopted so this model is comparable
    # with its siblings in the review's Figure 2. See the vignette Errata.
    bl_igg  <- fixed(4); label("Endogenous (treatment-naive) IgG concentration (g/L)")  # NOT reported for this model: van der Zeeuw 2026 section 2.4 simulation assumption for PID

    # Inter-individual variability. van der Zeeuw 2026 section 2.3 states that
    # all IIV values in the review were converted to apparent CV%, so
    # omega^2 = log(1 + CV^2).
    etalcl ~ 0.073414  # 27.6% CV; van der Zeeuw 2026 Table 4 IIV 'CL = 27.6'
    etalvc ~ 0.142943  # 39.2% CV; van der Zeeuw 2026 Table 4 IIV 'Vc = 39.2'

    # Residual error. van der Zeeuw 2026 Table 4 prints '5.3%' in the
    # proportional column and '-' in the additive column.
    propSd <- 0.053; label("Proportional residual error (fraction)")                    # van der Zeeuw 2026 Table 4: Prop = 5.3%
  })
  model({
    cl     <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    vc     <- exp(lvc + etalvc)
    ka     <- exp(lka)
    fdepot <- exp(lfdepot)

    kel <- cl / vc

    # `depot` receives subcutaneous doses; `central` receives intravenous
    # doses directly. States hold EXOGENOUS (therapeutic) IgG only.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central
    f(depot)      <-  fdepot

    # Observed total plasma IgG = exogenous concentration + endogenous baseline.
    Cc <- central / vc + bl_igg
    Cc ~ prop(propSd)
  })
}
