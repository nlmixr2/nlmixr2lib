Kyhl_2016_nalmefene <- function() {
  description <- "Population PK model for nalmefene in healthy volunteers (Kyhl 2016): two-compartment model with first-order absorption after oral dosing, separate absorption rates for tablet and solution formulations, and a link to mu-opioid receptor occupancy."
  reference <- "Kyhl LE, Li S, Faerch KU, Soegaard B, Larsen F, Areberg J. Population pharmacokinetics of nalmefene in healthy subjects and its relation to μ-opioid receptor occupancy. Br J Clin Pharmacol. 2016 Feb;81(2):290-300. doi: 10.1111/bcp.12805. Epub 2016 Jan 27. PMID: 26483076; PMCID: PMC4833148."
  vignette <- "Kyhl_2016_nalmefene"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Used (with sex and height) to derive LBM; not directly entered in the final PK model equations.",
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Linear (additive) effect on central volume, centered at 28 years.",
      source_name        = "AGE"
    ),
    LBM = list(
      description        = "Lean body mass",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on clearance with reference LBM 56.28 kg and exponent 0.626.",
      source_name        = "LBM"
    ),
    FED = list(
      description        = "Fed vs. fasted state at dosing",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (fasted)",
      notes              = "Multiplicative effect on oral bioavailability.",
      source_name        = "FED"
    ),
    TABLET = list(
      description        = "Tablet vs. solution formulation indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (solution)",
      notes              = "Switches between tablet and solution absorption rate constants.",
      source_name        = "TABLET"
    ),
    RIA_ASSAY = list(
      description        = "Assay type indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (LC-MS/MS)",
      notes              = "Switches the additive residual-error magnitude between radioimmunoassay (RIA) and LC-MS/MS measurements.",
      source_name        = "RIA_ASSAY"
    )
  )

  population <- list(
    n_subjects     = "TODO: from source paper",
    n_studies      = "TODO: from source paper",
    age_range      = "TODO: from source paper",
    age_median     = "TODO: from source paper (model centers AGE at 28 years)",
    weight_range   = "TODO: from source paper",
    weight_median  = "TODO: from source paper (LBM reference 56.28 kg)",
    sex_female_pct = "TODO: from source paper",
    race_ethnicity = "TODO: from source paper",
    disease_state  = "Healthy volunteers",
    dose_range     = "TODO: from source paper (oral tablet and oral solution formulations, fed and fasted)",
    regions        = "TODO: from source paper",
    notes          = "Pooled analysis of healthy-volunteer studies. Concentrations quantified by two assays (RIA and LC-MS/MS) with distinct additive residual-error magnitudes."
  )

  # Notes:
  # * Assuming that IIV% = sqrt(exp(variance) - 1)*100 (the typical, correct
  #   method, not one of the shortcut methods)
  # * Values were switched to log-scale from linear scale to standardize the
  #   model presentation with other models. This does not affect the model
  #   results.
  # * V2 and V3 were renamed to vc and vp to align with modeling used in
  #   `nlmixr2lib`
  # * All values are from Table 3
  ini({
    lka_tablet <- log(0.751); label("Absorption rate for oral tablet")
    etalka_tablet ~ sqrt(log((69.9/100)^2 + 1))
    lka_solution <- log(1.4)

    lcl <- log(60.4); label("Clearance (L/h)")
    etalcl ~ sqrt(log((18.7/100)^2 + 1))
    allo_lbm_cl <- 0.626; label("Lean body mass on CL (power)")
    lvc <- log(266); label("Volume of distribution, central compartment (L)")
    etalvc ~ sqrt(log((66.6/100)^2 + 1))
    e_age_vc <- -2.11; label("Effect of age on central volume of distribution (L/year)")
    lq <- log(109); label("Inter‐compartmental clearance (L/h)")
    etalq ~ sqrt(log((44.5/100)^2 + 1))
    lvp <- log(537); label("Volume of distribution, peripheral compartment (L)")
    etalvp ~ sqrt(log((45.8/100)^2 + 1))

    lf_oral <- log(0.406); label("Absolute oral bioavailability")
    etalf_oral ~ sqrt(log((20.1/100)^2 + 1))
    e_fed_f_oral <- 0.294

    propSd <- 0.094; label("Proportional residual error (fraction)")
    addSd_ria <- 0.0255; label("Additive residual error for radioimmunoassay (RIA)")
    addSd_lcms <- fixed(0); label("Additive residual error for LC-MS/MS")

    emax_mu_opioid <- 99.4; label("Maximum effect on mu-opioid receptor occupancy (percent)")
    lec50_mu_opioid <- log(0.338); label("Concentration producing 50% of maximum effect on mu-opioid receptor occupancy (ng/mL)")
  })
  model({
    ka_tablet <- exp(lka_tablet + etalka_tablet)
    ka_solution <- exp(lka_solution)
    ka <- ka_tablet*TABLET + ka_solution*(1 - TABLET)
    cl <- exp(lcl + etalcl) * (LBM/56.28)^allo_lbm_cl
    vc <- exp(lvc + etalvc) + e_age_vc*(AGE - 28)
    q <- exp(lq + etalq)
    vp <- exp(lvp + etalvp)

    f_oral <- exp(lf_oral + etalf_oral) * (1 + e_fed_f_oral*FED)

    kel <- cl/vc
    k12 <- q/vc
    k21 <- q/vp

    d/dt(depot) <- -ka*depot
    d/dt(central) <-  ka*depot - kel*central - k12*central + k21*peripheral1
    d/dt(peripheral1) <- k12*central - k21*peripheral1

    f(depot) <- f_oral
    # Unit conversion from mg/L to ng/mL
    Cc <- central / vc * 1000

    addSd <- addSd_ria*RIA_ASSAY + addSd_lcms*(1 - RIA_ASSAY)
    Cc ~ add(addSd) + prop(propSd)

    # These were not estimated within the primary model
    ec50_mu_opioid <- exp(lec50_mu_opioid)
    e_mu_opioid <- emax_mu_opioid*Cc/(ec50_mu_opioid + Cc)
  })
}
