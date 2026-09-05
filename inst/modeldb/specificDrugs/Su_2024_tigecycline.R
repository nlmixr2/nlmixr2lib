Su_2024_tigecycline <- function() {
  description <- "Two-compartment IV population PK model for tigecycline in critically ill adult patients, with creatinine clearance on clearance, body weight on both volumes, gamma-glutamyl transferase and total bilirubin on intercompartmental clearance, and albumin on peripheral volume (Su 2024)"
  reference <- "Su W, Song S, Liu J, Yu H, Feng B, Wu Y, Guo F, Yu Z. Population pharmacokinetics and individualized dosing of tigecycline for critically ill patients: a prospective study with intensive sampling. Front Pharmacol. 2024;15:1342947. doi:10.3389/fphar.2024.1342947"
  vignette <- "Su_2024_tigecycline"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Su 2024 Methods (plasma sampling,
  # LC-MS/MS assay) and the two-compartment structure of Table 2.
  compartmentData <- list(
    central     = list(analyte = "tigecycline", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "tigecycline", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description        = "Cockcroft-Gault creatinine clearance (raw, not BSA-normalized)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source column CCr, computed by the Cockcroft-Gault equation (Su 2024 Methods, 'Patient inclusion, drug administration and sample collection'); raw mL/min, NOT BSA-normalized to mL/min/1.73 m^2. Stored under the canonical CRCL column per inst/references/covariate-columns.md, which accepts raw mL/min when the source paper does not apply BSA normalization. Reference value 77.0 mL/min (population median, Table 1). Entered as an ADDITIVE LINEAR term with divisive normalization, matching the paper's Eq. 1 (theta_i = theta_pop + theta_cov * cov_i / cov_median): TVCL = exp(lcl) + e_crcl_cl * (CRCL / 77).",
      source_name        = "CCr"
    ),
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source column BW. Time-fixed baseline weight. Reference value 61 kg, taken from the divisor printed in the paper's final-model equations for V1 and V2 (the population median; Table 1 reports the MEAN as 62.3 +/- 12.4 kg, range 38-92.5). Enters as a power term on both V1 (exponent 1.95) and V2 (exponent 1.61) per the paper's Eq. 3.",
      source_name        = "BW"
    ),
    GGT = list(
      description        = "Serum gamma-glutamyl transferase activity",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source column GGT. Enters Q as a power term on the BASE-10 LOGARITHM of the covariate, normalized to the median of the log10 values: (log10(GGT) / 1.7)^0.956. The 1.7 divisor identifies the transform as base-10: log10(53) = 1.724 for the Table 1 median of 53 U/L, whereas ln(53) = 3.97. Domain restriction: the term requires GGT > 1 U/L so that log10(GGT) > 0 and the fractional power is real; the observed range was 8-1088 U/L.",
      source_name        = "GGT"
    ),
    TBILI = list(
      description        = "Total serum bilirubin",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source column TBIL, reported in SI units (umol/L) in Table 1. Enters Q as a power term on the BASE-10 LOGARITHM of the covariate, normalized to the median of the log10 values: (log10(TBIL) / 1.2)^-0.912. The 1.2 divisor identifies the transform as base-10: log10(15.7) = 1.196 for the Table 1 median of 15.7 umol/L, whereas ln(15.7) = 2.75. Domain restriction: requires TBIL > 1 umol/L; the observed range was 4-144 umol/L. The paper also tested TBIL on CL (Table S3 step 7) but eliminated it during backward elimination (step 11).",
      source_name        = "TBIL"
    ),
    ALB = list(
      description        = "Serum albumin",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source column ALB, reported in SI units (g/L) in Table 1. Enters V2 as a power term on the BASE-10 LOGARITHM of the covariate, normalized to the median of the log10 values: (log10(ALB) / 1.4)^4.52. The 1.4 divisor identifies the transform as base-10: log10(27.9) = 1.446 for the Table 1 median of 27.9 g/L, whereas ln(27.9) = 3.33. Domain restriction: requires ALB > 1 g/L; the observed range was 17.3-38.3 g/L. The large exponent makes V2 strongly sensitive to albumin across that range (roughly a 3-fold span), which the authors relate to the high plasma protein binding of tigecycline.",
      source_name        = "ALB"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 98L,
    n_studies      = 1L,
    n_observations = 751L,
    age_range      = "20-92 years",
    age_mean       = "63.4 years (SD 18.0)",
    weight_range   = "38-92.5 kg",
    weight_mean    = "62.3 kg (SD 12.4)",
    weight_median  = "61 kg (the normalization constant printed in the final-model V1 / V2 equations)",
    sex_female_pct = 33.7,
    race_ethnicity = "Not reported (single-center Chinese ICU cohort)",
    disease_state  = "Critically ill adults in a medical/surgical ICU receiving intermittent intravenous tigecycline. Infection site: intra-abdominal 37.8%, pulmonary 31.6%, other 22.4%, bloodstream 9.18%, trauma 7.14%, skin and soft tissue 6.12%. Predominant pathogens carbapenem-resistant A. baumannii (46.9%) and K. pneumoniae (24.5%).",
    dose_range     = "50 or 100 mg IV every 12 h as a 30-minute infusion (100 mg/day in 92.9% of patients, 200 mg/day in 7.14%)",
    regions        = "China (Sir Run Run Shaw Hospital, School of Medicine, Zhejiang University, Hangzhou)",
    renal_function = "Cockcroft-Gault creatinine clearance median 77.0 mL/min (IQR 51.1-142, range 6.40-338); patients on intermittent haemodialysis, peritoneal dialysis or CRRT were excluded (Supplementary Table S2)",
    hepatic_function = "Albumin median 27.9 g/L (IQR 25.8-30.3); GGT median 53 U/L (IQR 30-103); total bilirubin median 15.7 umol/L (IQR 10.0-30.4); ALT median 19 U/L; AST median 24 U/L",
    notes          = "Baseline demographics per Su 2024 Table 1. Single-center prospective study, July 2019 - July 2023. Intensive sampling: 8 samples per patient collected immediately before the seventh dose and at 0.5, 1, 2, 3, 4, 6 and 12 h post-dose, so every subject contributed a full steady-state dosing interval. Model built in NONMEM 7.5.0 with FOCEI; two-compartment structure selected over one-compartment on AIC (8242.914 vs 9325.681). Covariate screening (Supplementary Table S3) also tested age, BMI, sex, BUN, ALT, AST, ALP and the daily dose; DOSE on V1 and TBIL on CL entered forward inclusion but were removed during backward elimination, leaving the five covariates encoded here."
  )

  ini({
    # Structural parameters, Su 2024 Table 2 "Final model Estimate" column,
    # evaluated at the reference covariate values printed in the paper's
    # final-model equations (CCr 77 mL/min, BW 61 kg, log10(GGT) 1.7,
    # log10(TBIL) 1.2, log10(ALB) 1.4).
    #
    # NOTE on lcl: 3.09 L/h is the INTERCEPT of the paper's additive linear
    # CL ~ CCr model (Eq. 1), not the clearance of a typical patient. The
    # printed final-model equation is CL = 3.09 + (CCr/77) x 3.28, so the
    # typical value at the median CCr of 77 mL/min is 3.09 + 3.28 = 6.37 L/h.
    # The Results / Abstract text calls 3.09 L/h "the typical value of CL"
    # because it quotes the Table 2 row; the equation governs. See the
    # vignette "Assumptions and deviations" section for the corroboration.
    lcl <- log(3.09); label("Non-renal / baseline component of CL (intercept of the linear CCr model) (L/h)") # Su 2024 Table 2: CL = 3.09 L/h; final-model equation "CL (L/h) = 3.09+(CCr/77)x3.28"
    lvc <- log(32.1); label("Central volume of distribution V1 at 61 kg (L)")                                 # Su 2024 Table 2: V1 = 32.1 L; final-model equation "V1 (L) = 32.1x(BW/61)^1.95"
    lq  <- log(39.7); label("Intercompartmental clearance Q at median log10(GGT) and log10(TBIL) (L/h)")      # Su 2024 Table 2: Q = 39.7 L/h; final-model equation "Q (L/h) = 39.7x(logGGT/1.7)^0.956x(logTBIL/1.2)^-0.912"
    lvp <- log(113);  label("Peripheral volume of distribution V2 at 61 kg and median log10(ALB) (L)")        # Su 2024 Table 2: V2 = 113 L; final-model equation "V2 (L) = 113x(BW/61)^1.61x(logALB/1.4)^4.52"

    # Covariate effects, Su 2024 Table 2. CCr enters CL additively and
    # linearly (Eq. 1); the remaining four are power terms (Eq. 3), three of
    # them acting on the base-10 logarithm of the covariate.
    e_crcl_cl <- 3.28;   label("Renal CL slope per (CRCL / 77) (L/h)")                    # Su 2024 Table 2: theta_CLCR-CL = 3.28
    e_wt_vc   <- 1.95;   label("Power exponent on (WT / 61) for V1 (unitless)")           # Su 2024 Table 2: theta_BW-V1 = 1.95
    e_ggt_q   <- 0.956;  label("Power exponent on (log10(GGT) / 1.7) for Q (unitless)")   # Su 2024 Table 2: theta_GGT-Q = 0.956
    e_tbili_q <- -0.912; label("Power exponent on (log10(TBILI) / 1.2) for Q (unitless)") # Su 2024 Table 2: theta_TBIL-Q = -0.912
    e_wt_vp   <- 1.61;   label("Power exponent on (WT / 61) for V2 (unitless)")           # Su 2024 Table 2: theta_BW-V2 = 1.61
    e_alb_vp  <- 4.52;   label("Power exponent on (log10(ALB) / 1.4) for V2 (unitless)")  # Su 2024 Table 2: theta_ALB-V2 = 4.52

    # Inter-individual variability, Su 2024 Table 2 "Interindividual
    # variability" block, reported as omega (%). Exponential (log-normal)
    # IIV per Methods; converted to the internal variance scale via
    # omega^2 = log(CV^2 + 1). The paper states that IIV on V2 could not be
    # estimated ("it was not possible to estimate the interindividual
    # variability for V2"), so V2 carries no eta.
    etalcl ~ 0.070366 # log(0.270^2 + 1); Su 2024 Table 2 omega_CL = 27.0%
    etalvc ~ 0.422375 # log(0.725^2 + 1); Su 2024 Table 2 omega_V1 = 72.5%
    etalq  ~ 0.034735 # log(0.188^2 + 1); Su 2024 Table 2 omega_Q  = 18.8%

    # Residual error: proportional (Methods, "a proportional model [was] used
    # to describe ... residual variability"). Table 2 prints "sigma (%) 2.02".
    # Read literally as a 2.02% proportional CV that value is falsified by the
    # paper's own Figure 2 DV-vs-IPRED panel, whose scatter about the identity
    # line is an order of magnitude wider than +/-4%. It is encoded here as
    # the proportional-error VARIANCE expressed in percent (sigma^2 = 0.0202),
    # i.e. propSd = sqrt(0.0202) = 0.142, which matches that scatter. See the
    # vignette "Assumptions and deviations" section for the full adjudication.
    propSd <- 0.1421; label("Proportional residual error (fraction)") # Su 2024 Table 2: sigma = 2.02 (%), read as sigma^2 = 0.0202
  })
  model({
    # 1. Derived covariate terms. The paper's Q and V2 equations act on the
    # BASE-10 logarithm of the liver-chemistry covariates, each normalized by
    # the median of the log10 values rather than by the median covariate.
    # The base is identified by the printed divisors: log10(53) = 1.72 ~ 1.7,
    # log10(15.7) = 1.20 ~ 1.2, log10(27.9) = 1.45 ~ 1.4, against the Table 1
    # medians of 53 U/L, 15.7 umol/L and 27.9 g/L. Each requires its
    # covariate to exceed 1 in the units above so the ratio stays positive.
    lgGgt   <- log10(GGT)   / 1.7
    lgTbili <- log10(TBILI) / 1.2
    lgAlb   <- log10(ALB)   / 1.4

    # 2. Individual PK parameters. CL carries an additive linear covariate
    # term on CRCL with divisive normalization (TVCL = intercept + slope *
    # CRCL / 77), wrapped in exp(etalcl) so the IIV is log-normal on the
    # whole clearance.
    cl <- (exp(lcl) + e_crcl_cl * (CRCL / 77)) * exp(etalcl)
    vc <- exp(lvc + etalvc) * (WT / 61)^e_wt_vc
    q  <- exp(lq  + etalq)  * lgGgt^e_ggt_q * lgTbili^e_tbili_q
    vp <- exp(lvp)          * (WT / 61)^e_wt_vp * lgAlb^e_alb_vp

    # 3. Micro-constants
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # 4. ODE system. Tigecycline is given as a 30-minute intravenous infusion
    # into the central compartment.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # 5. Observation and error. Dose in mg, volumes in L -> central/vc has
    # units mg/L. Su 2024 plots concentrations in ug/L (= ng/mL), i.e.
    # 1000 x the value of Cc here.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
