Bardhi_2026_ampicillin_foal <- function() {
  description <- paste0(
    "Veterinary (horse, hospitalized neonatal foal). Two-compartment intravenous-bolus ",
    "population pharmacokinetic model with first-order elimination for ampicillin sodium in ",
    "12 critically ill hospitalized neonatal foals aged 9 to 194.5 h of life, each given ",
    "20 mg/kg intravenously every 6 h and sampled over the first 48 h of therapy ",
    "(205 plasma concentrations, Monolix 2024R1 SAEM). Postnatal age is the only retained ",
    "covariate and enters as a power function on clearance (exponent 0.40) and on the ",
    "peripheral volume (exponent 1.38), each normalised to the cohort weighted-mean age of ",
    "61.8 h. Body weight, serum albumin and serum creatinine were screened and not retained; ",
    "creatinine was the stronger clearance predictor but was rejected because a single ",
    "72-hour-old foal had a spurious 63.8 mg/L value (Bardhi 2026 Discussion), so they are ",
    "recorded in covariatesDataExcluded rather than in the model. The typical peripheral ",
    "volume (V2 = 22.12 L) and the between-subject SD of intercompartmental clearance ",
    "(omega_Q = 0.49) were held constant during estimation because their relative standard ",
    "errors exceeded the authors' 50% threshold. Residual variability uses the Monolix ",
    "'combined 1' form, in which the additive and proportional standard deviations add ",
    "linearly (sd = a + b * Cpred) rather than in quadrature; this maps onto rxode2's ",
    "combined1() error structure. Ampicillin is dosed directly into the `central` ",
    "compartment as an intravenous bolus. Clearances and volumes are absolute (L/h and L), ",
    "not weight-normalised: for the cohort median 48.5 kg foal the typical CL of 17.43 L/h ",
    "is 0.36 L/h/kg and Vss = Vc + Vp = 40.95 L is 0.84 L/kg. The paper's pharmacodynamic ",
    "layer is a probability-of-target-attainment analysis of fT > MIC computed by Monte ",
    "Carlo from this PK model using a fixed 15% plasma protein binding (unbound fraction ",
    "0.85); it adds no differential equations and is reproduced in the vignette rather than ",
    "in the model."
  )
  reference <- paste(
    "Bardhi A, Scala-Bertola J, Gehring R, Mariella J, Freccero F, Scarpellini R,",
    "Castagnetti C, Djerada Z, Barbarossa A. (2026).",
    "A population pharmacokinetic study of ampicillin therapy in hospitalized foals.",
    "Journal of Veterinary Internal Medicine 40(1):aalag021.",
    "doi:10.1093/jvimsj/aalag021.",
    sep = " "
  )
  vignette <- "Bardhi_2026_ampicillin_foal"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central     = list(analyte = "ampicillin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "ampicillin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    PNA = list(
      description        = "Postnatal age (chronological age since birth) at enrolment",
      units              = "months",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Time-fixed at enrolment. Bardhi 2026 reports foal age in HOURS ",
        "(cohort median 78.0 h, IQR 19.7-180.0 h, range 9-194.5 h) and centres the ",
        "covariate on the weighted-mean age of 61.8 h. The canonical PNA column carries ",
        "MONTHS (inst/references/covariate-columns.md), so model() converts back with ",
        "1 month = 30.4375 days = 730.5 h before forming the age ratio -- the same ",
        "reparameterisation used by Zhao_2018_omeprazole.R for a PNA reported in days. ",
        "The centring value 61.8 h therefore corresponds to 0.0846 months. Table 2 of ",
        "Bardhi 2026 prints '61.8 years'; that is a typographical slip for hours -- the ",
        "entire cohort is under 8 days old and both the Methods and Table 1 give age in h."
      ),
      source_name        = "Age (h)"
    )
  )

  # Screened by Bardhi 2026 (Materials and methods, 'Covariates analysis') but NOT
  # retained in the final model, so they must not appear in model(). Recorded here so
  # the paper's covariate screen is preserved without a 'declared but not referenced'
  # convention warning.
  covariatesDataExcluded <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Cohort median 48.5 kg (IQR 38.0-60.0), Bardhi 2026 Table 1. Tested one at a time ",
        "on every PK parameter and not retained; the published clearances and volumes are ",
        "absolute (L/h, L) rather than weight-normalised or allometrically scaled."
      ),
      source_name        = "TBW"
    ),
    ALB = list(
      description        = "Serum albumin measured before the first ampicillin dose",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Cohort median 28.9 g/L (IQR 26.5-30.7), Bardhi 2026 Table 1. Tested and not ",
        "retained. The Discussion notes that no sampled foal was hypoalbuminaemic, so a ",
        "protein-binding effect on distribution was considered unlikely."
      ),
      source_name        = "serum albumin"
    ),
    CREAT = list(
      description        = "Serum creatinine measured before the first ampicillin dose",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Cohort median 12.8 mg/L (IQR 9.5-14.7), Bardhi 2026 Table 1. Per the Discussion, ",
        "serum creatinine showed a STRONGER association with clearance than age ",
        "(Figure S1), but a single 63.8 mg/L value in the youngest foal -- probable ",
        "spurious hypercreatininaemia, which occurs in foals under 72 h old -- made the ",
        "association non-robust under bootstrapping, so the authors deliberately chose ",
        "age instead. Note the unit: mg/L, not the mg/dL common in human popPK."
      ),
      source_name        = "Scr"
    )
  )

  population <- list(
    species        = "horse (hospitalized neonatal foal)",
    n_subjects     = 12L,
    n_studies      = 1L,
    n_observations = 205L,
    age_range      = "9-194.5 h of life",
    age_median     = "78.0 h (IQR 19.7-180.0 h); mean 92.8 +/- 75.8 h; weighted mean 61.8 h",
    weight_range   = "IQR 38.0-60.0 kg",
    weight_median  = "48.5 kg",
    sex_female_pct = 50,
    breed          = "Standardbred (n = 6), Quarter horse (n = 3), Arabian horse (n = 2), mixed breed (n = 1)",
    disease_state  = paste0(
      "Critically ill hospitalized neonates. Diagnoses: omphalitis (n = 5), dysmaturity ",
      "(n = 1), neonatal isoerythrolysis (n = 1), umbilical cord haematoma (n = 1), total ",
      "failure of transfer of passive immunity with flexural deformity (n = 1), perinatal ",
      "asphyxia syndrome (n = 1), PAS with prematurity (n = 1), PAS with uroperitoneum ",
      "(n = 1). All foals survived. Most had some degree of dehydration or hypovolaemia and ",
      "received isotonic crystalloids."
    ),
    dose_range     = "20 mg/kg ampicillin sodium (Vetamplius, Fatro S.p.A.) intravenously every 6 h, dose expressed as ampicillin base",
    sampling       = paste0(
      "Blood drawn from a long-term jugular catheter by the push-pull technique at multiple ",
      "timepoints over the first 48 h of treatment; median 13 samples per foal (IQR 12-24). ",
      "Plasma assayed by UHPLC-MS/MS (Bardhi 2025), linear 0.3-100.0 mg/L, inter- and ",
      "intra-assay CV < 13%. One aberrantly high sample was excluded as probable ",
      "contamination and 17 below-LLOQ samples were treated as censored between 0 and LLOQ."
    ),
    regions        = "Italy (Veterinary Teaching Hospital, University of Bologna)",
    notes          = paste0(
      "Ethics: Animal Welfare Committee of the University of Bologna, Protocol No. 358467, ",
      "ID No. 4626. Cohort characteristics are Bardhi 2026 Table 1. Model selection favoured ",
      "2 compartments (BICc 1202.40 for the final model) over 1 compartment (1280.38) and ",
      "3 compartments (1232.47). Internal evaluation: 97% of observed values inside the ",
      "prediction intervals, all RSEs <= 44.8%, eta-shrinkage < 6.3%, 100% convergence over ",
      "1000 parametric bootstrap runs. Minimum inhibitory concentrations were measured in ",
      "the 7 pathogen-positive foals: Streptococcus equi subsp. zooepidemicus 0.06 mg/L ",
      "(n = 4), Clostridium perfringens 0.06 mg/L (n = 1), Bacillus licheniformis ",
      "0.25 mg/L (n = 1), Enterococcus faecalis 0.5 mg/L (n = 1)."
    )
  )

  ini({
    # =================================================================
    # Disposition -- Bardhi 2026 Table 2, 'Fixed effects' block.
    # Clearances are L/h and volumes are L for a whole foal; they are
    # NOT weight-normalised (body weight was screened and rejected as a
    # covariate).
    # =================================================================
    lcl <- log(17.43)
    label("Clearance (L/h)")  # Bardhi 2026 Table 2, Clpop = 17.43 L/h (rse-sa 17.74%, 95% CI 12.41-24.48)

    lvc <- log(18.83)
    label("Central compartment volume of distribution (L)")  # Bardhi 2026 Table 2, V1pop = 18.83 L (rse-sa 11.21%, 95% CI 15.15-23.41)

    lq <- log(16.79)
    label("Intercompartmental clearance (L/h)")  # Bardhi 2026 Table 2, Qpop = 16.79 L/h (rse-sa 32.31%, 95% CI 9.30-30.34)

    # V2 was HELD CONSTANT during the final estimation step: in the
    # initial 2-compartment model its RSE was 76.6%, above the authors'
    # 50% threshold, so once log(Age) had been added to V2 the typical
    # value was set to 22.12 L, which decreased BICc by 12.4 (Bardhi
    # 2026 Results, 'Population pharmacokinetic modeling'). Table 2
    # accordingly prints '22.12 (Fixed)' with no confidence interval and
    # no bootstrap row.
    lvp <- fixed(log(22.12))
    label("Peripheral compartment volume of distribution (L)")  # Bardhi 2026 Table 2, V2pop = 22.12 L (Fixed)

    # =================================================================
    # Postnatal-age covariate effects -- Bardhi 2026 Table 2 footnote:
    #   log(Cl) = log(Clpop) + beta_Cl_logAge * log(Age/61.8) + etaCl
    #   log(V2) = log(V2pop) + beta_V2_logAge * log(Age/61.8) + etaV2
    # Exponentiating makes each a power function of the age ratio, which
    # is how model() encodes it. Age is in HOURS in the source; see the
    # PNA covariateData note for the months conversion.
    # =================================================================
    e_pna_cl <- 0.40
    label("Power exponent of postnatal age on clearance (unitless)")  # Bardhi 2026 Table 2, beta_Cl_logAge = 0.40 (rse-sa 39.19%, 95% CI 0.092-0.7, Wald P = 2.56e-2)

    e_pna_vp <- 1.38
    label("Power exponent of postnatal age on peripheral volume (unitless)")  # Bardhi 2026 Table 2, beta_V2_logAge = 1.38 (rse-sa 44.8%, 95% CI 0.17-2.59, Wald P = 1.07e-2)

    # =================================================================
    # Between-subject variability -- Bardhi 2026 Table 2, 'SD of the
    # random effects' block. Monolix reports omega as a STANDARD
    # DEVIATION on the log scale ("omega was the standard deviation
    # estimated by Monolix", Materials and methods, 'Basic model
    # building'); nlmixr2's ini() takes VARIANCES, so each entry below
    # is the published omega squared. The variance-covariance matrix of
    # the random effects was diagonal, so there are no off-diagonals.
    # =================================================================
    etalcl ~ 0.3025   # Bardhi 2026 Table 2, omega_Cl = 0.55 (rse-sa 25.18%); variance = 0.55^2
    etalvc ~ 0.1024   # Bardhi 2026 Table 2, omega_V1 = 0.32 (rse-sa 32.50%); variance = 0.32^2

    # omega_Q had an RSE of 74.8% in the initial 2-compartment model,
    # above the authors' 50% threshold, so it was held at 0.49; that step
    # alone decreased BICc by 2.67 (Bardhi 2026 Results), and Table 2
    # accordingly prints it with no confidence interval and no bootstrap
    # row.
    etalq ~ fixed(0.2401)  # Bardhi 2026 Table 2, omega_Q = 0.49; variance = 0.49^2

    etalvp ~ 1.0816   # Bardhi 2026 Table 2, omega_V2 = 1.04 (rse-sa 39.53%); variance = 1.04^2

    # =================================================================
    # Residual error -- Bardhi 2026 Table 2, 'Error model parameters'.
    # The selected model is Monolix 'combined 1':
    #   Cobs = Cpred + (a + b * Cpred) * eps
    # i.e. the additive and proportional standard deviations add
    # LINEARLY. rxode2's default for add() + prop() is the combined2
    # (in-quadrature) form, so combined1() is declared explicitly in
    # model(). Combined 1 beat every alternative on BICc: +431.43 for
    # constant, +7.53 for proportional and +2.33 for combined 2
    # (Bardhi 2026 Results / Table S2).
    # =================================================================
    addSd <- 0.12
    label("Additive residual error SD (mg/L)")  # Bardhi 2026 Table 2, a = 0.12 (rse-sa 36.18%, 95% CI 0.063-0.23)

    propSd <- 0.40
    label("Proportional residual error (fraction)")  # Bardhi 2026 Table 2, b = 0.40 (rse-sa 8.30%, 95% CI 0.34-0.47)
  })

  model({
    # 1. Derived covariate term.
    # The canonical PNA covariate column is in MONTHS; Bardhi 2026
    # reports foal age in HOURS and centres on 61.8 h, so convert here.
    # 1 month = 30.4375 days = 730.5 h.
    pnaHours <- PNA * 730.5
    ageRatio <- pnaHours / 61.8

    # 2. Individual parameters. exp(beta * log(Age/61.8)) is the power
    # form (Age/61.8)^beta of the Table 2 footnote equations.
    cl <- exp(lcl + etalcl) * ageRatio^e_pna_cl
    vc <- exp(lvc + etalvc)
    q  <- exp(lq + etalq)
    vp <- exp(lvp + etalvp) * ageRatio^e_pna_vp

    # 3. Micro-constants.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 4. ODE system. Intravenous bolus with no absorption delay, dosed
    # directly into `central` (Bardhi 2026 Materials and methods:
    # "1-, 2-, and 3-compartment intravenous bolus models with no delay
    # in administration and first-order elimination").
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # 5. Observation and error.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd) + combined1()
  })
}
