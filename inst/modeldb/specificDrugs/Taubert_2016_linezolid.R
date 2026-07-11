Taubert_2016_linezolid <- function() {
  description <- paste(
    "Population PK model for linezolid in critically ill ICU patients",
    "(Taubert 2016). Two-compartment model with first-order absorption",
    "from a depot, first-order elimination, complete oral bioavailability",
    "(F = 1), and a combined proportional + additive residual error model.",
    "Covariate effects: body weight (power on Vc, exponent 1.31) and",
    "peritonitis (multiplier on Vc, factor 1.53); fibrinogen concentration",
    "(power on CL, exponent 0.04), serum lactate concentration (power on CL,",
    "exponent -0.21), and ARDS (multiplier on CL, factor 1.82). All",
    "continuous covariates are normalised to the patient-group-1 median",
    "(WT = 76 kg, FIB = 13.0 umol/L, LACT = 1.91 mmol/L)."
  )
  reference <- paste(
    "Taubert M, Zoller M, Maier B, Frechen S, Scharf C, Holdt L-M,",
    "Frey L, Vogeser M, Fuhr U, Zander J. (2016).",
    "Predictors of inadequate linezolid concentrations after standard",
    "dosing in critically ill patients.",
    "Antimicrob Agents Chemother 60(9):5254-5261.",
    "doi:10.1128/AAC.00356-16."
  )
  vignette <- "Taubert_2016_linezolid"
  units <- list(
    time          = "hour",
    dosing        = "mg",
    concentration = "mg/L"
  )

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power-of-ratio scaling on Vc (exponent 1.31) with reference 76 kg",
        "(patient-group-1 median, Taubert 2016 Table 1)."
      ),
      source_name        = "WT"
    ),
    FIB = list(
      description        = "Plasma fibrinogen concentration",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power-of-ratio scaling on CL (exponent 0.04) with reference",
        "13.0 umol/L (patient-group-1 median, Taubert 2016 Table 1).",
        "Time-varying covariate updated daily in the source data."
      ),
      source_name        = "Fibrinogen"
    ),
    LACT = list(
      description        = "Serum lactate concentration",
      units              = "mmol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power-of-ratio scaling on CL (exponent -0.21) with reference",
        "1.91 mmol/L (patient-group-1 median, Taubert 2016 Table 1).",
        "Time-varying covariate updated daily in the source data."
      ),
      source_name        = "Lactate"
    ),
    DIS_ARDS = list(
      description        = paste(
        "Acute respiratory distress syndrome status:",
        "1 = ARDS present, 0 = no ARDS."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "no ARDS = 0",
      notes              = paste(
        "Multiplicative factor on CL (1.82 when ARDS = 1, Taubert 2016",
        "page 5256). Time-varying covariate updated daily in the source",
        "data; ARDS was present in 15 of 52 patients (29 %) at baseline."
      ),
      source_name        = "ARDS"
    ),
    DIS_PERIT = list(
      description        = paste(
        "Peritonitis status: 1 = peritonitis present, 0 = no peritonitis."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "no peritonitis = 0",
      notes              = paste(
        "Multiplicative factor on Vc (1.53 when peritonitis = 1, Taubert",
        "2016 page 5256). Peritonitis was present in 9 of 52 patients",
        "(17 %) at baseline."
      ),
      source_name        = "Peritonitis"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 52L,
    n_studies      = 1L,
    age_range      = "28-84 years",
    age_median     = "58 years",
    weight_range   = "44-120 kg",
    weight_median  = "76 kg",
    sex_female_pct = 36.5,
    race_ethnicity = "Not reported (single-centre German ICU cohort)",
    disease_state  = paste(
      "Medical-surgical critically ill adults with clinically suspected or",
      "confirmed severe infections; principal infection sites pneumonia",
      "(67 %), peritonitis (17 %) and other; subgroups with ARDS (29 %),",
      "peritonitis (17 %), liver transplantation (13 %) and lung",
      "transplantation (29 %)."
    ),
    dose_range     = paste(
      "Linezolid 600 mg b.i.d., short-duration intravenous infusion",
      "(10-120 min) or oral; before day 1, patients had received 0 to 4",
      "linezolid doses."
    ),
    regions        = paste(
      "Single centre at the Department of Anesthesiology, Hospital of the",
      "Ludwig Maximilians University of Munich, Germany",
      "(ClinicalTrials.gov NCT01793012)."
    ),
    notes          = paste(
      "Patient group 1 (n = 52) is the linezolid PK fit cohort; patient",
      "group 2 (n = 134, treated with meropenem, piperacillin-tazobactam,",
      "cefepime or ciprofloxacin) supplied covariate values for the 67 000-",
      "subject simulation that drove the covariate-effect-size analysis.",
      "Demographics from Taubert 2016 Table 1; the modelled cohort is",
      "patient group 1.",
      "Median renal function (estimated CrCl, day 1) was 81 mL/min",
      "(range 5-293); 31 % were on continuous renal replacement therapy",
      "at baseline. APACHE II median 28 (range 9-38)."
    )
  )

  ini({
    # Structural PK parameters - reference values are the patient-group-1
    # medians (WT = 76 kg, FIB = 13.0 umol/L, LACT = 1.91 mmol/L). All
    # point estimates are bootstrap medians (n = 1000) per Taubert 2016
    # Supplemental Table S3.
    lka  <- log(1.72);  label("Absorption rate constant (1/h)")          # Suppl Table S3 (KA bootstrap median 1.72; 95 % CI 0.66-2.38)
    lcl  <- log(7.92);  label("Clearance at reference covariates (L/h)") # Suppl Table S3 (CL bootstrap median 7.92; 95 % CI 6.27-9.74)
    lvc  <- log(15);    label("Central volume at reference WT (L)")      # Suppl Table S3 (Vc bootstrap median 15.00; 95 % CI 8.58-20.67)
    lq   <- log(65.59); label("Inter-compartmental clearance (L/h)")     # Suppl Table S3 (Q  bootstrap median 65.59; 95 % CI 45.92-109.61)
    lvp  <- log(26.55); label("Peripheral volume (L)")                   # Suppl Table S3 (Vp bootstrap median 26.55; 95 % CI 21.24-32.63)

    # Covariate effects. The paper specifies the covariate equations
    # explicitly (page 5256, equation block):
    #   Vc = theta_Vc * exp(eta1) * (WT / 76)^1.31  * 1.53^DIS_PERIT
    #   CL = theta_CL * exp(eta2) * (FIB / 13.0)^0.04 *
    #        (LACT / 1.91)^-0.21 * 1.82^DIS_ARDS
    e_wt_vc      <- 1.31;  label("Exponent of WT/76 on Vc (unitless)")             # Suppl Table S3 (Weight on Vc median 1.31; 95 % CI 0.70-2.35)
    e_fib_cl     <- 0.04;  label("Exponent of FIB/13.0 on CL (unitless)")          # Suppl Table S3 (Fibrinogen on CL median 0.04; 95 % CI 0-0.08)
    e_lact_cl    <- -0.21; label("Exponent of LACT/1.91 on CL (unitless)")         # Suppl Table S3 (Lactate on CL median -0.21; 95 % CI -0.30 to -0.11)
    e_dis_ards_cl  <- log(1.82); label("Log-multiplier on CL when DIS_ARDS = 1")   # Suppl Table S3 (ARDS on CL median 1.82; 95 % CI 1.26-2.62)
    e_dis_perit_vc <- log(1.53); label("Log-multiplier on Vc when DIS_PERIT = 1")  # Suppl Table S3 (Peritonitis on Vc median 1.53; 95 % CI 1.16-2.11)

    # Inter-individual variability. Taubert 2016 page 5256 reports
    # "Interindividual variability terms were kept for Vc (37%) and CL
    # (58%) (coefficients of variation)." The variance entered into the
    # log-normal eta is omega^2 = log(1 + CV^2), so:
    #   omega^2_etalcl = log(1 + 0.58^2) = log(1.3364) = 0.2901
    #   omega^2_etalvc = log(1 + 0.37^2) = log(1.1369) = 0.1283
    # The supplement reports the same values (Suppl Table S3 omega^2_CL = 58%
    # [95 % CI 46-73] and omega^2_Vc = 37% [95 % CI 23-61]) and explicitly
    # labels the column as CV %.
    etalcl ~ 0.2901  # Suppl Table S3 (omega^2_CL reported as 58 % CV; log-normal variance = log(1 + 0.58^2))
    etalvc ~ 0.1283  # Suppl Table S3 (omega^2_Vc reported as 37 % CV; log-normal variance = log(1 + 0.37^2))

    # Residual error. Taubert 2016 Suppl Table S3 reports the combined
    # error model as variance estimates on the linear concentration scale:
    #   sigma^2_proportional = 0.115 (95 % CI 0.073-0.182)
    #   sigma^2_additive     = 0.005 (95 % CI 0.001-0.008) in (mg/L)^2
    # nlmixr2 ini() uses SDs, so propSd = sqrt(0.115) ~ 0.339 and
    # addSd = sqrt(0.005) ~ 0.0707 mg/L. The supplement was acquired
    # post-hoc via the PMC proof-of-work supplement route after the
    # initial extraction sidecar resolved with a missing-RUV fall-back
    # (fixed(0)); the published Table S3 values are used here because the
    # supplement is now on disk.
    propSd <- 0.339;  label("Proportional residual SD (fraction)")           # Suppl Table S3 (sigma^2_prop = 0.115; SD = sqrt(0.115))
    addSd  <- 0.0707; label("Additive residual SD (mg/L)")                   # Suppl Table S3 (sigma^2_add  = 0.005; SD = sqrt(0.005))
  })

  model({
    # Individual PK parameters - log-additive covariates per the
    # paper's published equations (page 5256) translated to the log
    # scale for log-normal IIV:
    #   log(Vc) = lvc + etalvc + e_wt_vc      * log(WT / 76) +
    #             e_dis_perit_vc * DIS_PERIT
    #   log(CL) = lcl + etalcl + e_fib_cl     * log(FIB / 13.0) +
    #             e_lact_cl * log(LACT / 1.91) + e_dis_ards_cl * DIS_ARDS
    vc <- exp(lvc + etalvc +
              e_wt_vc * log(WT / 76) +
              e_dis_perit_vc * DIS_PERIT)
    cl <- exp(lcl + etalcl +
              e_fib_cl * log(FIB / 13) +
              e_lact_cl * log(LACT / 1.91) +
              e_dis_ards_cl * DIS_ARDS)
    vp <- exp(lvp)
    q  <- exp(lq)
    ka <- exp(lka)

    # Micro-constants for the 2-compartment ODE system.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # 2-compartment PK with first-order absorption from a depot;
    # complete oral bioavailability (F = 1) per Taubert 2016 page 5256.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Observation: total linezolid plasma concentration (mg/L). The
    # paper's combined error model is additive + proportional in linear
    # (mg/L) space (Suppl Table S3).
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
