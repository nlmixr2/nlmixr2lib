Tien_2024_meropenem <- function() {
  description <- paste(
    "Meropenem in-bottle degradation coupled to a two-compartment",
    "population PK model in critically ill adults without CRRT (Tien",
    "2024). The infusion solution is an ODE state (infusion_container)",
    "dosed as a bolus when the bottle is hung; intact drug degrades by",
    "first-order kinetics and is pumped into central at constant volume",
    "rate over TINF, so the delivery rate is Dose * exp(-kdeg * t) / TINF.",
    "kdeg comes from a linear mixed-effects fit to in-vitro stability",
    "data (7 generic brands, 1 or 2 g/48 mL in saline, 25/30/37 degC,",
    "0-8 h) with categorical solution-temperature and solution-strength",
    "effects and correlated random intercept and slope per bottle series;",
    "the in-vitro fraction remaining is a second output with additive",
    "error on log(C/C0). The PK layer is the Ehmann 2019 ICU model as",
    "recoded by Tien: CL linear in Cockcroft-Gault creatinine clearance",
    "(centred 80.8 mL/min) up to an estimated hinge at 154 mL/min, power",
    "weight effect on Vc and linear albumin effect on Vp, with IIV on CL,",
    "Vc and Vp and one occasion of IOV on CL."
  )
  reference <- paste(
    "Tien NTN, Binh VN, Ha PTT, Lan DTN, Cho YS, Long NP, Shin JG, Anh NH,",
    "Quan TA, Tuan DN, Tiep NK, Thach PT, Anh NH, Hoa VD.",
    "Optimization of meropenem continuous",
    "infusion based on Monte Carlo simulation integrating with",
    "degradation study. PLoS One. 2024;19(12):e0313764.",
    "doi:10.1371/journal.pone.0313764. Stability model: Methods Eq. 1,",
    "Table 2 and the full-precision lmer output in S3 File; coupled",
    "PK/degradation equations: Methods 'Simulation procedure' and S4 File.",
    "PK parameters as recoded in S4 File from Ehmann L et al. Int J",
    "Antimicrob Agents. 2019;54(3):309-317,",
    "doi:10.1016/j.ijantimicag.2019.06.016."
  )
  vignette <- "Tien_2024_meropenem"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")
  dosing <- c("infusion_container", "central")

  covariateData <- list(
    CRCL = list(
      description = "Creatinine clearance (Cockcroft-Gault)",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Raw Cockcroft-Gault creatinine clearance in mL/min, NOT",
        "BSA-normalised (S4 File variable Clcr; Ehmann 2019 CLCR_CG).",
        "CL rises linearly with CRCL centred at 80.8 mL/min up to",
        "crcl_hinge (154 mL/min) and is constant above it. Tien simulated",
        "only CRCL <= 150 mL/min (Methods; Fig 6)."
      ),
      source_name = "Clcr"
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on Vc normalised to 70 kg (S4 File: Vc = tvVc * (Wt/70)^WT_Vc).",
      source_name = "Wt"
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "The source model uses g/dL centred at 2.8 g/dL (S4 File: Vp =",
        "tvVp * (1 + ALB_Vp * (Alb - 2.8))); the canonical g/L column is",
        "divided by 10 inside model(). Typical patient 28 g/L."
      ),
      source_name = "Alb"
    ),
    TINF = list(
      description = "Duration over which the current infusion bottle is delivered",
      units = "h",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-varying: set it on each record to the infusion length of the",
        "bottle most recently dosed into infusion_container (0.5 h for the",
        "500 mg loading dose, then 3, 4, 6 or 8 h in the paper). TINF must",
        "not exceed the bottle-renewal interval, since the input rate is",
        "computed from the most recent bottle only (podo/tad of",
        "infusion_container)."
      ),
      source_name = "duration"
    ),
    INFUSION_TEMP_30C = list(
      description = "Infusion solution held at 30 degC (1 = yes; 0 = 25 or 37 degC)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "Temperature entered the stability model as a categorical effect",
        "on kdeg with 25 degC as reference (Table 2). Mutually exclusive",
        "with INFUSION_TEMP_37C; both 0 = 25 degC. Only 25, 30 and 37 degC",
        "were studied, so no other temperature can be represented."
      ),
      source_name = "temp"
    ),
    INFUSION_TEMP_37C = list(
      description = "Infusion solution held at 37 degC (1 = yes; 0 = 25 or 30 degC)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "See INFUSION_TEMP_30C. Mutually exclusive with it.",
      source_name = "temp"
    ),
    FORM_MEROPENEM_SOLN2G48ML = list(
      description = "Infusion solution strength 2 g/48 mL (1) vs 1 g/48 mL (0)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "Meropenem in 0.9% sodium chloride; 1 g/48 mL is the reference",
        "(Table 2). Only these two strengths were studied."
      ),
      source_name = "conc"
    )
  )

  compartmentData <- list(
    # The prepared infusion solution before it enters the patient.
    infusion_container = list(analyte = "meropenem", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "meropenem", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "meropenem", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human (PK layer) + in vitro (infusion-solution stability)",
    n_subjects = 48,
    n_studies = 1,
    disease_state = paste(
      "Critically ill adults with severe infections, not on continuous",
      "renal replacement therapy (Ehmann 2019 PK cohort)."
    ),
    renal_function = paste(
      "Cockcroft-Gault creatinine clearance median 70.8 mL/min, 5th-95th",
      "percentiles 34.8-160 mL/min (Methods, citing Ehmann 2019)."
    ),
    n_observations = "1376 plasma samples (PK); 365 in-vitro stability observations",
    in_vitro = paste(
      "Seven meropenem generic brands reconstituted in 0.9% sodium chloride",
      "at 1 g/48 mL and 2 g/48 mL, held at 25, 30 and 37 degC and assayed",
      "by HPLC hourly for 8 h; 42 bottle series (7 brands x 3 temperatures",
      "x 2 strengths), 365 of 378 observations retained (Results; S2 and",
      "S3 Files). Brand was tested and not retained (Table 1, Model 4,",
      "p = 0.14)."
    ),
    regions = "Vietnam (stability study); the PK cohort's region is not stated in Tien 2024",
    notes = paste(
      "Typical patient of the deterministic simulations: CRCL 80.8 mL/min,",
      "WT 70 kg, albumin 2.8 g/dL. Monte Carlo simulations used 1000",
      "patients per regimen and condition, evaluated over treatment day 1."
    )
  )

  ini({
    # In-vitro stability (linear mixed-effects model fit3, S3 File lmer
    # summary; Table 2). lmer fits y = log(C/C0) with a NEGATIVE time
    # coefficient; kdeg here is its negation, so every effect on kdeg is
    # the negated lmer estimate.
    infusion_log_intercept <- -0.0002568
    label("Log fraction remaining at bottle preparation, log(C/C0) at t = 0 (unitless)") # S3 File fit3 '(Intercept) -0.0002568'; Table 2 mu -0.0003
    kdeg <- 0.0086083
    label("Degradation rate constant of meropenem in the infusion solution at 25 degC, 1 g/48 mL (1/h)") # S3 File fit3 'time -0.0086083'; Table 2 beta0 0.0086
    e_infusion_temp_30c_kdeg <- 0.0044789
    label("Additive increase in kdeg at 30 vs 25 degC (1/h)") # S3 File fit3 'time:temp30 degC -0.0044789'; Table 2 beta1 30 degC 0.0045
    e_infusion_temp_37c_kdeg <- 0.0185050
    label("Additive increase in kdeg at 37 vs 25 degC (1/h)") # S3 File fit3 'time:temp37 degC -0.0185050'; Table 2 beta1 37 degC 0.0185
    e_form_meropenem_soln2g48ml_kdeg <- 0.0041512
    label("Additive increase in kdeg for 2 g/48 mL vs 1 g/48 mL (1/h)") # S3 File fit3 'time:conc2g/48mL -0.0041512'; Table 2 beta2 0.0042

    # Ehmann 2019 two-compartment ICU model as recoded in S4 File
    # (Ehmann2019 RxODE block, @param line).
    lcl <- log(9.25)
    label("Clearance at CRCL 80.8 mL/min (L/h)") # S4 File tvCL = 9.25 (RSE 4.6%)
    lvc <- log(7.89)
    label("Central volume at 70 kg (L)") # S4 File tvVc = 7.89 (RSE 11.9%)
    lq <- log(28.4)
    label("Intercompartmental clearance (L/h)") # S4 File tvQ = 28.4 (RSE 16.1%)
    lvp <- log(16.1)
    label("Peripheral volume at albumin 2.8 g/dL (L)") # S4 File tvVp = 16.1 (RSE 7.4%)
    e_crcl_cl <- 0.00977
    label("Fractional change in CL per mL/min of CRCL from 80.8 mL/min (1/(mL/min))") # S4 File CLCR_CL = 0.00977 (RSE 9.2%)
    lcrcl_hinge <- log(154)
    label("CRCL above which CL no longer rises (mL/min)") # S4 File CLCR_INF = 154 (RSE 6.9%)
    e_wt_vc <- 0.945
    label("Power exponent of (WT/70) on Vc (unitless)") # S4 File WT_Vc = 0.945 (RSE 16.6%)
    e_alb_vp <- -0.202
    label("Fractional change in Vp per g/dL of albumin from 2.8 g/dL (1/(g/dL))") # S4 File ALB_Vp = -0.202 (RSE 36.6%)

    # In-vitro random intercept and slope per bottle series (S3 File fit3
    # random effects: variances 4.119e-05 and 2.149e-06, Corr 0.94). The
    # lmer correlation is between the intercept and the time coefficient;
    # kdeg is the negated time coefficient, so the covariance here is
    # -0.94 * 0.006418 * 0.001466 = -8.844e-06.
    etainfusion_log_intercept + etakdeg ~ c(4.119e-05, -8.844e-06, 2.149e-06)

    # IIV and IOV, log-normal; variance = log(CV^2 + 1) from the S4 File
    # Omega lotri (27.1, 31.5, 16.9 and 12.5 %CV).
    etalcl ~ 0.07086938
    etalvc ~ 0.09460539
    etalvp ~ 0.02816074
    etaiov_cl_1 ~ 0.01550419

    # Residual error
    propSd <- fixed(0)
    label("Proportional residual error on Cc (fraction); not reported in Tien 2024 or its supplements")
    expSd_infusion_fraction_remaining <- 0.010127
    label("Residual SD of log(C/C0) in the infusion solution (log scale)") # S3 File fit3 'Residual Std.Dev. 0.010127'; Table 2 sigma 0.0101
  })

  model({
    # Degradation rate constant of the hung bottle (Eq. 1; Table 2)
    kdeg_container <- kdeg +
      e_infusion_temp_30c_kdeg * INFUSION_TEMP_30C +
      e_infusion_temp_37c_kdeg * INFUSION_TEMP_37C +
      e_form_meropenem_soln2g48ml_kdeg * FORM_MEROPENEM_SOLN2G48ML +
      etakdeg
    log_intercept_container <- infusion_log_intercept + etainfusion_log_intercept

    # Time since the current bottle was hung. The bottle holds
    # A = Dose * exp(-kdeg * t) * (1 - t / TINF): intact drug decays while
    # the pump removes solution at a constant volume rate, so the drug
    # delivered is Dose * exp(-kdeg * t) / TINF while t < TINF (Methods,
    # 'Simulation procedure'; S4 File Depot/duration input).
    t_container <- tad(infusion_container)
    delivery <- 0
    if (t_container >= 0 && t_container < TINF) {
      delivery <- podo(infusion_container) * exp(-kdeg_container * t_container) / TINF
    }

    # Ehmann 2019 PK (S4 File @pre block). Above crcl_hinge CL is held at
    # its value at the hinge.
    crcl_hinge <- exp(lcrcl_hinge)
    cl_iiv <- exp(lcl + etalcl)
    cl <- cl_iiv * exp(etaiov_cl_1) * (1 + e_crcl_cl * (min(CRCL, crcl_hinge) - 80.8))
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    q <- exp(lq)
    vp <- exp(lvp + etalvp) * (1 + e_alb_vp * (ALB / 10 - 2.8))

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(infusion_container) <- -kdeg_container * infusion_container - delivery
    d/dt(central) <- delivery - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    Cc <- central / vc
    # In-vitro fraction of the initial solution concentration remaining
    # (Eq. 1 back-transformed); NA before the first bottle is hung.
    infusion_fraction_remaining <- exp(log_intercept_container - kdeg_container * t_container)

    Cc ~ prop(propSd)
    infusion_fraction_remaining ~ lnorm(expSd_infusion_fraction_remaining)
  })
}
