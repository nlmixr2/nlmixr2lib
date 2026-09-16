Li_2022_immunoglobulin <- function() {
  description <- "Two-compartment population PK model with first-order subcutaneous absorption for polyclonal immunoglobulin G across intravenous, subcutaneous and hyaluronidase-facilitated subcutaneous products in primary immunodeficiency, scaled on lean body mass (Li 2022)"
  reference   <- "Li Z, Follman K, Freshwater E, Engler F, Yel L. Integrated population pharmacokinetics of immunoglobulin G following intravenous or subcutaneous administration of various immunoglobulin products in patients with primary immunodeficiencies. Int Immunopharmacol. 2022;113(Pt A):109331. doi:10.1016/j.intimp.2022.109331 -- parameter values transcribed from the secondary source: van der Zeeuw SL, van Tilburg SJ, Jacobs BC, Koch BCP, Dalm VASH, Crombag MBS, Preijers T. Population pharmacokinetics and pharmacodynamics of immunoglobulins: a systematic review. Clin Pharmacokinet. 2026;65(6):813-30. doi:10.1007/s40262-026-01641-5, Table 4 (reference 42)"
  vignette    <- "vanderZeeuw_2026_immunoglobulin"
  units       <- list(time = "day", dosing = "g", concentration = "g/L")

  covariateData <- list(
    LBM = list(
      description        = "Lean body mass",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reference 47 kg. This is the only model among the immunodeficiency models in van der Zeeuw 2026 that scales on LEAN body mass rather than total body weight (section 3.2.1.4: 'Notably, Li et al. derived the LBM and incorporated it as a covariate on CL instead of BW'). Exponents are FIXED: 0.75 on CL and Q, 1 on Vc and Vp (sections 3.2.1.3 and 3.2.1.4). The review derives LBM with the Boer formula (section 2.3); for a 70 kg, 170 cm male that gives approximately 55 kg.",
      source_name        = "LBM"
    ),
    FORM_IG_HYALURONIDASE = list(
      description        = "Hyaluronidase-facilitated subcutaneous immunoglobulin (fSCIg) product indicator (1 = fSCIg, 0 = unfacilitated SCIg)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (unfacilitated subcutaneous immunoglobulin, SCIg)",
      notes              = "Applies to subcutaneous doses only. Li 2022 is the only model in van der Zeeuw 2026 that includes fSCIg data (14.5% of included patients), and product type was the only covariate any of the reviewed models retained on bioavailability (section 3.2.1.2). Bioavailability is 70.5% without hyaluronidase and 79.4% with it.",
      source_name        = "IgG product / hyaluronidase yes-no"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Total body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Tested but not retained: van der Zeeuw 2026 Table 3 lists 'BW, BMI, LBM, sex, age, IgG product, hyaluronidase product yes/no' as covariates tested, and LBM was carried into the final model in place of total body weight.",
      source_name = "BW"
    )
  )

  compartmentData <- list(
    depot       = list(analyte = "immunoglobulin G", units = "g", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "immunoglobulin G", units = "g", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "immunoglobulin G", units = "g", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 340L,
    n_studies      = 8L,
    age_range      = "2.0-83.0 years",
    age_median     = "31.5 years",
    weight_range   = "11.9-162 kg",
    weight_median  = "66.0 kg",
    sex_female_pct = round(100 * 161 / 340, 1),
    race_ethnicity = "Not reported",
    disease_state  = "Primary immunodeficiency (PID) on immunoglobulin replacement therapy",
    dose_range     = "IVIg 300-1000 mg/kg every 3-4 weeks; SCIg weekly or Q4W equivalents (130-145% of the IVIg dose in some studies); fSCIg every 3-4 weeks",
    regions        = "United States, Canada, Europe",
    notes          = "Largest cohort in van der Zeeuw 2026: pooled analysis of NCT00814320, NCT01412385, NCT01218438, NCT00161993, NCT00157079, NCT00546871, NCT00782106 and NCT03277313 (Tables 1 and 2). Baseline IgG not reported in Table 1. Endogenous IgG was ESTIMATED rather than fixed in this model (section 3.2.1.6)."
  )

  ini({
    # Structural parameters. van der Zeeuw 2026 Table 4, row 'Li et al. (2022)
    # [42]'. Reference LEAN body mass 47 kg as printed in the table.
    lcl     <- log(0.183); label("Clearance at LBM 47 kg (L/day)")                      # van der Zeeuw 2026 Table 4: CL = 0.183 (LBM/47)^0.75
    lvc     <- log(3.01);  label("Central volume of distribution at LBM 47 kg (L)")     # van der Zeeuw 2026 Table 4: Vc = 3.01 (LBM/47)^1
    lq      <- log(0.353); label("Intercompartmental clearance at LBM 47 kg (L/day)")   # van der Zeeuw 2026 Table 4: Q = 0.353 (LBM/47)^0.75
    lvp     <- log(1.40);  label("Peripheral volume of distribution at LBM 47 kg (L)")  # van der Zeeuw 2026 Table 4: Vp = 1.40 (LBM/47)^1
    lka     <- log(0.395); label("First-order subcutaneous absorption rate constant (1/day)")  # van der Zeeuw 2026 Table 4: Ka = 0.395 day^-1
    lfdepot <- log(0.705); label("Unfacilitated subcutaneous bioavailability relative to intravenous (fraction)")  # van der Zeeuw 2026 Table 4: F1 without hyaluronidase = 70.5%

    # Allometric exponents on lean body mass -- both FIXED in the source
    # (van der Zeeuw 2026 sections 3.2.1.3 and 3.2.1.4: 'using a fixed
    # allometric scaling component of 0.75' for CL and Q; 'the exponents on Vc
    # and Vp were fixed to 1').
    e_lbm_cl_q  <- fixed(0.75); label("Allometric exponent on CL and Q (unitless)")     # van der Zeeuw 2026 sections 3.2.1.3-3.2.1.4: fixed at 0.75
    e_lbm_vc_vp <- fixed(1);    label("Allometric exponent on Vc and Vp (unitless)")    # van der Zeeuw 2026 section 3.2.1.3: exponents on Vc and Vp fixed to 1

    # Hyaluronidase effect on bioavailability, entered as the fSCIg/SCIg ratio
    # of the two printed bioavailability values so both source numbers remain
    # traceable.
    e_hyal_fdepot <- 1.1262; label("fSCIg-vs-SCIg multiplicative ratio on bioavailability (F_fSCIg / F_SCIg)")  # van der Zeeuw 2026 Table 4: F1 without hyaluronidase = 70.5%, with hyaluronidase = 79.4%; 79.4/70.5 = 1.1262

    # Endogenous IgG. van der Zeeuw 2026 section 3.2.1.6: 'Li et al. estimated
    # a typical value for endogenous IgG concentration [42]'. The value is not
    # in Table 4; the Discussion (section 4) supplies it: 'One popPK model did
    # estimate a parameter describing the endogenous IgG. However, the data
    # consisted mostly of patients who received IgG treatment prior to the
    # study, which may compromise the precision of the estimated value of
    # 6.15 g/L.' Entered as fixed() because no uncertainty is reported and the
    # value is recovered from prose rather than a parameter table. NOTE: 6.15
    # also coincides with the SID median baseline IgG reported for Tortorici
    # 2019 in Table 1; see the vignette Errata for why the Discussion
    # attribution to this model is nevertheless the better reading.
    bl_igg  <- fixed(6.15); label("Endogenous (treatment-naive) IgG concentration (g/L)")  # van der Zeeuw 2026 section 4 Discussion: estimated value of 6.15 g/L

    # Inter-individual variability (apparent CV%, van der Zeeuw 2026 section
    # 2.3): omega^2 = log(1 + CV^2).
    etalcl ~ 0.191183  # 45.9% CV; van der Zeeuw 2026 Table 4 IIV 'CL = 45.9'
    etalvc ~ 0.040773  # 20.4% CV; van der Zeeuw 2026 Table 4 IIV 'Vc = 20.4'

    # Combined residual error. van der Zeeuw 2026 Table 4 prints both terms as
    # VARIANCES ('sigma^2 = ...'), so each is entered here as its square root.
    # Table 4 footnote c attaches specifically to THIS row: 'Not clear whether
    # expressed as variance or standard deviation'. The squared notation is
    # taken at face value for consistency with the two sibling rows that use
    # it (Zhang 2020, Navarro-Mora 2022), where the variance reading is the
    # only one giving plausible magnitudes. Under the alternative
    # standard-deviation reading the values would be 0.962 g/L and 10.8%.
    # See the vignette Errata.
    addSd  <- 0.980816; label("Additive residual error on total IgG (g/L)")             # van der Zeeuw 2026 Table 4: Add sigma^2 = 0.962 (footnote c); sqrt = 0.980816
    propSd <- 0.328634; label("Proportional residual error (fraction)")                 # van der Zeeuw 2026 Table 4: Prop sigma^2 = 0.108 (footnote c); sqrt = 0.328634
  })
  model({
    cl     <- exp(lcl + etalcl) * (LBM / 47)^e_lbm_cl_q
    vc     <- exp(lvc + etalvc) * (LBM / 47)^e_lbm_vc_vp
    q      <- exp(lq)           * (LBM / 47)^e_lbm_cl_q
    vp     <- exp(lvp)          * (LBM / 47)^e_lbm_vc_vp
    ka     <- exp(lka)
    fdepot <- exp(lfdepot) * e_hyal_fdepot^FORM_IG_HYALURONIDASE

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # `depot` receives subcutaneous doses (SCIg or fSCIg, distinguished by
    # FORM_IG_HYALURONIDASE); `central` receives intravenous doses directly,
    # for which bioavailability is 1 by definition. States hold EXOGENOUS
    # (therapeutic) IgG only.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    f(depot)          <-  fdepot

    # Observed total plasma IgG = exogenous concentration + endogenous baseline.
    Cc <- central / vc + bl_igg
    Cc ~ add(addSd) + prop(propSd)
  })
}
