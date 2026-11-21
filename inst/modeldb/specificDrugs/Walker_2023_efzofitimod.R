Walker_2023_efzofitimod <- function() {
  description <- paste(
    "Three-compartment population PK model for intravenous efzofitimod",
    "(an immunomodulatory neuropilin-2-directed Fc fusion protein) in",
    "healthy adult volunteers and adults with pulmonary sarcoidosis",
    "(Walker 2023), with first-order elimination from the central",
    "compartment, power effects of creatinine clearance and total",
    "bilirubin on clearance, power effects of body weight on the central",
    "and first peripheral volumes, and a separate log-scale residual",
    "error for the healthy-volunteer and patient studies."
  )
  reference <- paste(
    "Walker G, Adams R, Guy L, Chandrasekaran A, Kinnersley N, Ramesh P,",
    "Zhang L, Brown F, Niranjan V. Exposure-response analyses of",
    "efzofitimod in patients with pulmonary sarcoidosis. Front Pharmacol.",
    "2023;14:1258236. doi:10.3389/fphar.2023.1258236. PMCID PMC10580085.",
    "Parameter estimates are from Supplementary Table S1."
  )
  vignette <- "Walker_2023_efzofitimod"
  paper_specific_residual_sds <- c("expSdHv", "expSdPat")
  units <- list(time = "day", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Baseline body weight. Power effects on V1 (exponent 0.77) and V2 (exponent 1.52), referenced to 82 kg, the median weight of the pooled PK-evaluable population quoted in Walker 2023 Results section 3.2 ('relative to the median value of 82 kg'). Observed range 55-157 kg.",
      source_name = "BWT"
    ),
    CRCL = list(
      description = "Baseline creatinine clearance, raw (not BSA-normalised)",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on CL (exponent 0.71). Walker 2023 reports CrCL in mL/min with a 10th-90th percentile range of 103-187 mL/min but does not name the estimating equation or print the reference value of the power term; the model uses 139 mL/min, the geometric mean of the 10th and 90th percentiles (an estimate of the median under a log-symmetric distribution). Stored under the canonical CRCL column with raw mL/min units, following the raw Cockcroft-Gault precedents in the CRCL register entry.",
      source_name = "CRCL"
    ),
    TBILI = list(
      description = "Baseline total serum bilirubin",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on CL (exponent -0.17). Walker 2023 reports a 10th-90th percentile range of 3.2-19.6 umol/L but does not print the reference value of the power term; the model uses 7.9 umol/L, the geometric mean of the 10th and 90th percentiles (an estimate of the median under a log-symmetric distribution).",
      source_name = "BILI"
    ),
    DIS_HEALTHY = list(
      description = "Healthy-volunteer indicator: 1 = Phase 1 healthy-volunteer study (ATYR1923-C-001), 0 = Phase 1b/2a pulmonary-sarcoidosis study (ATYR1923-C-002)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (patient with pulmonary sarcoidosis)",
      notes = "Selects the study-specific log-additive residual variance (Supplementary Table S1 'epsilon (Phase 1)' and 'epsilon (Phase 2)'). Has no effect on any structural parameter. The two studies map one-to-one onto healthy volunteers and patients, so the study split is encoded on the healthy-participant indicator, as in Hughes_2022_brepocitinib.R.",
      source_name = "STUDY"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened on CL and the volumes; not significant (Walker 2023 Results 3.2: 'Race, sex, and age were not significant covariates')."
    ),
    SEXF = list(
      description = "Sex, 1 = female",
      units = "(binary)",
      type = "binary",
      notes = "Screened; not significant (Walker 2023 Results 3.2)."
    ),
    RACE_BLACK = list(
      description = "Black or African American race indicator",
      units = "(binary)",
      type = "binary",
      notes = "Race screened; not significant (Walker 2023 Results 3.2). Ethnicity was also screened and not retained."
    ),
    ALB = list(
      description = "Baseline serum albumin",
      units = "g/L",
      type = "continuous",
      notes = "Screened on CL and the volumes; not significant. ALP, ALT, AST, serum creatinine and eGFR were also screened on CL and not retained (Walker 2023 Methods 2.2 and Results 3.2)."
    )
  )

  compartmentData <- list(
    central = list(analyte = "efzofitimod", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "efzofitimod", units = "mg", specimen = "not applicable", verified = TRUE),
    peripheral2 = list(analyte = "efzofitimod", units = "mg", specimen = "not applicable", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 49L,
    n_studies = 2L,
    n_observations = 565L,
    age_range = "Phase 1 mean 23.8 (SD 3.94) years; Phase 1b/2a mean 52.3 (SD 10.3) years",
    weight_range = "55-157 kg",
    weight_median = "82 kg",
    sex_female_pct = 53.1,
    race_ethnicity = c(White = 61.2, Black = 22.4, Asian = 10.2, Other = 6.1),
    disease_state = "Healthy volunteers (Phase 1, n = 24) and patients with pulmonary sarcoidosis on oral corticosteroids (Phase 1b/2a, n = 25)",
    dose_range = "Phase 1: single 1-h IV infusion of 0.03, 0.1, 0.3, 1, 3 or 5 mg/kg. Phase 1b/2a: 1, 3 or 5 mg/kg IV every 4 weeks for six doses.",
    regions = "Australia (Phase 1), United States (Phase 1b/2a)",
    notes = "Demographics from Walker 2023 Table 3 (PK-evaluable population). Sex: 23 male (11 + 12), 26 female (13 + 13). Race: 30 White, 11 Black or African American, 5 Asian, 3 Other. Serum efzofitimod was measured by a validated immunoassay; rich sampling in Phase 1, pre-dose plus two post-dose samples on Day 1 and Week 20 in Phase 1b/2a."
  )

  ini({
    # Structural parameters at the reference covariates WT = 82 kg,
    # CRCL = 139 mL/min and TBILI = 7.9 umol/L.
    lcl <- log(1.68); label("Clearance CL (L/day)") # Supplementary Table S1 'CL' = 1.68 L/day (RSE 2.93%); Results 3.2
    lvc <- log(3.94); label("Central volume of distribution V1 (L)") # Supplementary Table S1 'V1' = 3.94 L (RSE 3.32%)
    lq <- log(2.4); label("Intercompartmental clearance Q2, central to peripheral1 (L/day)") # Supplementary Table S1 'Q2' = 2.4 L/day (RSE 6.96%)
    lvp <- log(3.3); label("First peripheral volume of distribution V2 (L)") # Supplementary Table S1 'V2' = 3.3 L (RSE 17.92%)
    lq2 <- log(1.2); label("Intercompartmental clearance Q3, central to peripheral2 (L/day)") # Supplementary Table S1 'Q3' = 1.2 L/day (RSE 13.76%)
    lvp2 <- log(7.43); label("Second peripheral volume of distribution V3 (L)") # Supplementary Table S1 'V3' = 7.43 L (RSE 6.19%)

    # Covariate effects (power form)
    e_tbili_cl <- -0.17; label("Power exponent of total bilirubin on CL (unitless)") # Supplementary Table S1 'BILI_CL' = -0.17 (RSE 23.02%)
    e_crcl_cl <- 0.71; label("Power exponent of creatinine clearance on CL (unitless)") # Supplementary Table S1 'CRCL_CL' = 0.71 (RSE 20.87%)
    e_wt_vc <- 0.77; label("Power exponent of body weight on V1 (unitless)") # Supplementary Table S1 'BWT_V1' = 0.77 (RSE 20.11%)
    e_wt_vp <- 1.52; label("Power exponent of body weight on V2 (unitless)") # Supplementary Table S1 'BWT_V2' = 1.52 (RSE 31.87%)

    # Between-subject variability: Supplementary Table S1 reports the
    # omega^2 variances directly, so they enter unchanged.
    etalcl ~ 0.03 # Supplementary Table S1 'omega2 (CL)' = 0.03 (RSE 45.98%)
    etalvc ~ 0.01 # Supplementary Table S1 'omega2 (V1)' = 0.01 (RSE 28.11%)
    etalq ~ 0.06 # Supplementary Table S1 'omega2 (Q2)' = 0.06 (RSE 44.05%)
    etalvp ~ 0.1 # Supplementary Table S1 'omega2 (V2)' = 0.1 (RSE 36.1%)
    etalq2 ~ fixed(0) # Supplementary Table S1 'omega2 (Q3)' = 0 FIX
    etalvp2 ~ 0.01 # Supplementary Table S1 'omega2 (V3)' = 0.01 (RSE 54.14%)

    # Residual error: log-additive, one variance per study. The table
    # reports the variances (0.01 and 1.63), so the log-scale SDs are their
    # square roots.
    expSdHv <- 0.1; label("Log-scale residual SD in the Phase 1 healthy-volunteer study") # Supplementary Table S1 'epsilon2 (Phase 1)' = 0.01 (RSE 17.94%) -> sqrt(0.01)
    expSdPat <- 1.2767; label("Log-scale residual SD in the Phase 1b/2a pulmonary-sarcoidosis study") # Supplementary Table S1 'epsilon2 (Phase 2)' = 1.63 (RSE 7.48%) -> sqrt(1.63) = 1.2767
  })

  model({
    # Individual parameters. Power covariate forms per Supplementary Table S1
    # ('Power of ... effect on ...'); reference values in the ini() header.
    cl <- exp(lcl + etalcl) * (CRCL / 139)^e_crcl_cl * (TBILI / 7.9)^e_tbili_cl
    vc <- exp(lvc + etalvc) * (WT / 82)^e_wt_vc
    q <- exp(lq + etalq)
    vp <- exp(lvp + etalvp) * (WT / 82)^e_wt_vp
    q2 <- exp(lq2 + etalq2)
    vp2 <- exp(lvp2 + etalvp2)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # Three-compartment disposition with first-order elimination from the
    # central compartment (Walker 2023 Figure 2). Doses are IV infusions
    # into central.
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1 - k13 * central + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2

    # mg / L = ug/mL; x 1000 gives ng/mL, the unit of the exposure metrics
    # (ng*day/mL) used throughout the paper.
    Cc <- 1000 * central / vc

    expSd <- expSdHv * DIS_HEALTHY + expSdPat * (1 - DIS_HEALTHY)
    Cc ~ lnorm(expSd)
  })
}
