NavarroMora_2022_immunoglobulin <- function() {
  description <- "Two-compartment population PK model with first-order subcutaneous absorption for polyclonal immunoglobulin G given subcutaneously or intravenously in primary immunodeficiency (Navarro-Mora 2022)"
  reference   <- "Navarro-Mora G, Alberti JJ, Mondou E, Vilardell D, Vicente Torres J, Ayguasanosa J, et al. Pharmacokinetic modeling and simulation of subcutaneous and intravenous IgG dosing in patients with primary immunodeficiency diseases. Int Immunopharmacol. 2022;104:108472. doi:10.1016/j.intimp.2021.108472 -- parameter values transcribed from the secondary source: van der Zeeuw SL, van Tilburg SJ, Jacobs BC, Koch BCP, Dalm VASH, Crombag MBS, Preijers T. Population pharmacokinetics and pharmacodynamics of immunoglobulins: a systematic review. Clin Pharmacokinet. 2026;65(6):813-30. doi:10.1007/s40262-026-01641-5, Table 4 (reference 37)"
  vignette    <- "vanderZeeuw_2026_immunoglobulin"
  units       <- list(time = "day", dosing = "g", concentration = "g/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric power scaling on CL (estimated exponent 0.744), Vc (0.686) and Vp (1.04), reference weight 65.7 kg. This is the only model in van der Zeeuw 2026 with an ESTIMATED weight exponent on Vp: section 3.2.1.3, 'While Li et al. and Lee et al. fixed the allometric scaling component to 1, Navarro-Mora et al. estimated it to be 1.04'. Q carries no weight term in Table 4.",
      source_name        = "BW"
    )
  )

  compartmentData <- list(
    depot       = list(analyte = "immunoglobulin G", units = "g", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "immunoglobulin G", units = "g", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "immunoglobulin G", units = "g", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 95L,
    n_studies      = 3L,
    age_range      = "Not reported as a range; per-study means 10.8-42.5 years",
    age_median     = "NCT00389324: 42.5 (SD 15.8); NCT01465958: 10.8 (SD 3.7); NCT02604810: 36.8 (SD 21.36) years, all means",
    weight_range   = "16.7-153.0 kg",
    weight_median  = "65.7 kg",
    sex_female_pct = round(100 * 55 / 95, 1),
    race_ethnicity = "Not reported",
    disease_state  = "Primary immunodeficiency (PID) on immunoglobulin replacement therapy",
    dose_range     = "IVIg median 495 mg/kg (range 278-902) every 3 or 4 weeks; SCIg median 184.8 mg/kg (range 72.0-303.5) every 3 or 4 weeks",
    regions        = "United States, Canada",
    notes          = "Pooled analysis of NCT00389324, NCT01465958 and NCT02604810 (van der Zeeuw 2026 Tables 1 and 2). Baseline IgG not reported. The authors evaluated a fixed endogenous IgG of 4 g/L against a lower value of 1.5 g/L and retained 4 g/L because the change in PK parameters was minimal (section 3.2.1.6)."
  )

  ini({
    # Structural parameters. van der Zeeuw 2026 Table 4, row 'Navarro-Mora et
    # al. (2022) [37]'. Reference weight 65.7 kg as printed in the table.
    lcl     <- log(0.150); label("Clearance for a 65.7 kg patient (L/day)")             # van der Zeeuw 2026 Table 4: CL = 0.150 (BW/65.7)^0.744
    lvc     <- log(3.06);  label("Central volume of distribution for a 65.7 kg patient (L)")  # van der Zeeuw 2026 Table 4: Vc = 3.06 (BW/65.7)^0.686
    lq      <- log(0.474); label("Intercompartmental clearance (L/day)")                # van der Zeeuw 2026 Table 4: Q = 0.474
    lvp     <- log(1.93);  label("Peripheral volume of distribution for a 65.7 kg patient (L)")  # van der Zeeuw 2026 Table 4: Vp = 1.93 (BW/65.7)^1.04
    lka     <- log(0.246); label("First-order subcutaneous absorption rate constant (1/day)")  # van der Zeeuw 2026 Table 4: Ka = 0.246
    lfdepot <- log(0.705); label("Subcutaneous bioavailability relative to intravenous (fraction)")  # van der Zeeuw 2026 Table 4: F1 = 70.5%

    # Allometric exponents, all estimated (van der Zeeuw 2026 sections 3.2.1.3
    # and 3.2.1.4).
    e_wt_cl <- 0.744; label("Allometric exponent on CL (unitless)")                     # van der Zeeuw 2026 Table 4: (BW/65.7)^0.744
    e_wt_vc <- 0.686; label("Allometric exponent on Vc (unitless)")                     # van der Zeeuw 2026 Table 4: (BW/65.7)^0.686
    e_wt_vp <- 1.04;  label("Allometric exponent on Vp (unitless)")                     # van der Zeeuw 2026 Table 4 and section 3.2.1.3: (BW/65.7)^1.04

    # Endogenous IgG. van der Zeeuw 2026 section 3.2.1.6: 'Navarro-Mora et al.
    # evaluated a fixed value of 4 g/L and a lower value of 1.5 g/L. However,
    # because changes in PK parameters were minimal, endogenous IgG was fixed
    # at 4 g/L in the final model [37].'
    bl_igg  <- fixed(4); label("Endogenous (treatment-naive) IgG concentration (g/L)")  # van der Zeeuw 2026 section 3.2.1.6, held constant at 4 g/L in the final model

    # Inter-individual variability (apparent CV%, van der Zeeuw 2026 section
    # 2.3): omega^2 = log(1 + CV^2).
    etalcl ~ 0.080216  # 28.9% CV; van der Zeeuw 2026 Table 4 IIV 'CL = 28.9'
    etalvc ~ 0.036572  # 19.3% CV; van der Zeeuw 2026 Table 4 IIV 'Vc = 19.3'
    etalka ~ 0.517203  # 82.3% CV; van der Zeeuw 2026 Table 4 IIV 'Ka = 82.3'

    # Combined residual error. van der Zeeuw 2026 Table 4 prints both terms as
    # VARIANCES ('sigma^2 = ...'), so each is entered here as its square root.
    # The variance reading is the only one giving plausible magnitudes on this
    # row: read as standard deviations the proportional term would be 0.28%.
    addSd  <- 0.817313; label("Additive residual error on total IgG (g/L)")             # van der Zeeuw 2026 Table 4: Add sigma^2 = 0.668; sqrt = 0.817313
    propSd <- 0.053009; label("Proportional residual error (fraction)")                 # van der Zeeuw 2026 Table 4: Prop sigma^2 = 0.00281; sqrt = 0.053009
  })
  model({
    cl     <- exp(lcl + etalcl) * (WT / 65.7)^e_wt_cl
    vc     <- exp(lvc + etalvc) * (WT / 65.7)^e_wt_vc
    q      <- exp(lq)
    vp     <- exp(lvp)          * (WT / 65.7)^e_wt_vp
    ka     <- exp(lka + etalka)
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
    Cc ~ add(addSd) + prop(propSd)
  })
}
