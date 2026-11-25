Netterberg_2019_atezolizumab_itac <- function() {
  description <- paste(
    "PK / biomarker model for atezolizumab and plasma interferon-inducible",
    "T-cell alpha chemoattractant (ITAC, CXCL11) in relapsed or refractory",
    "non-small cell lung cancer (Netterberg 2019, phase I study PCD4989g). A",
    "two-compartment PK model with linear elimination drives, without an",
    "effect-compartment delay, a linear stimulation of ITAC release from a",
    "depletable precursor pool into plasma (pool indirect-response model).",
    "ITAC was evaluated as a predictor of tumor shrinkage but was not retained",
    "in the final tumor model. Time is in days. Companion model:",
    "modellib('Netterberg_2019_atezolizumab_il18_tgi').",
    sep = " "
  )
  reference <- paste(
    "Netterberg I, Li CC, Molinero L, Budha N, Sukumaran S, Stroh M,",
    "Jonsson EN, Friberg LE. A PK/PD Analysis of Circulating Biomarkers and",
    "Their Relationship to Tumor Response in Atezolizumab-Treated non-small",
    "Cell Lung Cancer Patients. Clin Pharmacol Ther. 2019;105(2):486-495.",
    "doi:10.1002/cpt.1198.",
    "PK parameters are in Supplementary Table S3; ITAC parameters in Table 1.",
    sep = " "
  )
  vignette <- "Netterberg_2019_atezolizumab"

  units <- list(
    time = "day",
    dosing = "mg",
    concentration = "ug/mL"
  )

  # The ITAC precursor pool and the ITAC plasma state are paper-specific
  # states.
  paper_specific_compartments <- c("itac_pool", "itac")

  covariateData <- list()

  compartmentData <- list(
    central = list(analyte = "atezolizumab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "atezolizumab", units = "mg", specimen = "not applicable", verified = TRUE),
    itac_pool = list(
      analyte = "interferon-inducible T-cell alpha chemoattractant (precursor pool)",
      units = "pg/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    itac = list(
      analyte = "interferon-inducible T-cell alpha chemoattractant",
      units = "pg/mL",
      specimen = "plasma",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 88L,
    n_studies = 1L,
    disease_state = "relapsed / refractory non-small cell lung cancer",
    dose_range = "atezolizumab 10, 15 or 20 mg/kg or 1,200 mg flat dose as a 30-60 min IV infusion every 3 weeks (one patient received 1 mg/kg)",
    regions = "not reported (phase Ia study PCD4989g, NCT01375842)",
    n_observations = "ITAC: 458 above LOQ (Supplementary Table S2)",
    notes = "The PK model was developed on 167 patients from PCD4989g who received 10, 15 or 20 mg/kg every 3 weeks; its individual estimates were used for the 88 NSCLC patients analysed for biomarkers (population PK parameters and data method). Age, weight, sex and race are not tabulated in Netterberg 2019."
  )

  ini({
    # Two-compartment PK (Supplementary Table S3)
    lcl <- log(0.233)
    label("Clearance CL (L/day)") # Supplementary Table S3: CL = 0.233 L/day (RSE 4.5%)
    lvc <- log(3.35)
    label("Central volume Vc (L)") # Supplementary Table S3: Vc = 3.35 L (RSE 2.3%)
    lq <- log(0.452)
    label("Intercompartmental clearance Q (L/day)") # Supplementary Table S3: Q = 0.452 L/day (RSE 8.6%)
    lvp <- log(3.36)
    label("Peripheral volume Vp (L)") # Supplementary Table S3: Vp = 3.36 L (RSE 7.8%)

    # ITAC pool indirect-response model (Table 1, ITAC model; Eq. 3)
    lkin_itac <- log(3.07)
    label("ITAC zero-order input rate Rin (pg/mL/day)") # Table 1: Rin = 3.07 pg/mL/day (RSE 24%)
    lkout_itac <- log(7.91e-3)
    label("ITAC first-order turnover rate constant kout (1/day)") # Table 1: kout = 7.91e-3 1/day (RSE 23%)
    lbl_itac_pool <- log(271)
    label("Baseline ITAC concentration in the precursor pool Pool0 (pg/mL)") # Table 1: Pool0 = 271 pg/mL (RSE 20%)
    lslope_itac <- log(0.0599)
    label("Linear effect of atezolizumab concentration on ITAC release (mL/ug)") # Table 1: Slope = 0.0599 mL/ug (RSE 35%)

    # IIV: CV% of log-normal etas; omega^2 = log(1 + CV^2)
    etalcl ~ 0.162430 # Supplementary Table S3: IIV CL 42% -> log(1 + 0.42^2)
    etalvc ~ 0.070363 # Supplementary Table S3: IIV Vc 27% -> log(1 + 0.27^2)
    etalkin_itac ~ 0.207339 # Table 1: IIV Rin (ITAC) 48% -> log(1 + 0.48^2)
    etalbl_itac_pool ~ 0.842643 # Table 1: IIV Pool0 (ITAC) 115% -> log(1 + 1.15^2)

    # Residual error (log-transformed data; additive and proportional error for
    # ITAC, proportional for the PK)
    propSd <- 0.187
    label("Proportional residual error, atezolizumab serum concentration (fraction)") # Supplementary Table S3: residual error 18.7% (RSE 6.6%)
    propSd_itac <- 0.283
    label("Proportional residual error, plasma ITAC (fraction)") # Table 1: ITAC proportional residual error 28.3% (RSE 3.9%)
    addSd_itac <- 51.8
    label("Additive residual error, plasma ITAC (pg/mL)") # Table 1: ITAC additive residual error 51.8 pg/mL (RSE 15%)
  })
  model({
    # Individual PK parameters
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q <- exp(lq)
    vp <- exp(lvp)
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Individual ITAC parameters. Baseline plasma ITAC is Rin / kout and the
    # pool-to-plasma release rate constant is krel = kout * BioM0 / Pool0
    # (Eq. 3).
    kin_itac <- exp(lkin_itac + etalkin_itac)
    kout_itac <- exp(lkout_itac)
    bl_itac_pool <- exp(lbl_itac_pool + etalbl_itac_pool)
    slope_itac <- exp(lslope_itac)
    bl_itac <- kin_itac / kout_itac
    krel_itac <- kout_itac * bl_itac / bl_itac_pool

    Cc <- central / vc
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Pool indirect-response model (Eq. 3) driven directly by the central
    # concentration (an effect compartment did not improve the ITAC fit)
    eff_itac <- slope_itac * Cc
    d/dt(itac_pool) <- kin_itac - krel_itac * itac_pool * (1 + eff_itac)
    d/dt(itac) <- krel_itac * itac_pool * (1 + eff_itac) - kout_itac * itac
    itac_pool(0) <- bl_itac_pool
    itac(0) <- bl_itac

    Cc ~ prop(propSd)
    itac ~ add(addSd_itac) + prop(propSd_itac)
  })
}
