Stringer_2015_pioglitazone <- function() {
  description <- paste(
    "Population pharmacodynamic cascading indirect-response model for",
    "fasting plasma glucose (FPG) and glycosylated hemoglobin (HbA1c) in",
    "Japanese type 2 diabetes mellitus (T2DM) patients receiving",
    "pioglitazone or a non-thiazolidinedione oral glucose-lowering drug",
    "background over 2.5-4 years. FPG follows a zero-order production /",
    "first-order elimination turnover with a proportional disease-",
    "progression term on production; a time-driven Emax model stimulates",
    "FPG elimination to capture titration-related exposure ramp-up.",
    "HbA1c is driven by a power function of FPG plus a separate",
    "FPG-independent zero-order input with a linear-in-time disease-",
    "progression component. No PK data were collected; the drug effect is",
    "treatment-cohort- and time-driven, distinguished by different Emax",
    "and ET50 values in the pioglitazone vs control cohorts",
    "(Stringer 2015 Diabetes Technology and Therapeutics)."
  )
  reference <- paste(
    "Stringer F, DeJongh J, Enya K, Koumura E, Danhof M, Kaku K.",
    "Evaluation of the Long-Term Durability and Glycemic Control of",
    "Fasting Plasma Glucose and Glycosylated Hemoglobin for Pioglitazone",
    "in Japanese Patients with Type 2 Diabetes.",
    "Diabetes Technol Ther. 2015;17(3):215-222.",
    "doi:10.1089/dia.2014.0222"
  )
  vignette <- "Stringer_2015_pioglitazone"
  units <- list(
    time          = "day",
    dosing        = "none (drug effect is time-driven within treatment cohort; no PK data collected)",
    concentration = "FPG in mg/dL; HbA1c in % (NGSP)"
  )
  paper_specific_compartments <- c("fpg", "hba1c")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. analyte/specimen proposed by a local model from the
  # model description; units derived from the units block. verified = FALSE
  # means NOT checked against the source paper.
  compartmentData <- list(
    fpg   = list(analyte = "Fasting Plasma Glucose", units = NA_character_, specimen = "plasma", verified = FALSE),
    hba1c = list(analyte = "Glycosylated Hemoglobin", units = NA_character_, specimen = "blood cell", verified = FALSE)
  )

  covariateData <- list(
    SEXF = list(
      description        = "Sex (female indicator)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 (female)",
      notes              = paste(
        "Additive fractional multiplier on baseline FPG only.",
        "Male subjects (SEXF = 0) have baseline FPG that is",
        "(1 + sex_effect_fpg) times the female-typical value; the paper",
        "reports a +5% higher FPG baseline for males (Stringer 2015",
        "Table 2 'Gender on FPG BSL' and Results section 'Covariate",
        "analysis'). No effect on baseline HbA1c, KoutG, KoutH, Emax,",
        "or the disease-progression rates."
      ),
      source_name        = "SEX"
    ),
    TRT = list(
      description        = "Treatment-cohort indicator (0 = control non-TZD, 1 = pioglitazone with or without other non-TZD oral glucose-lowering drugs)",
      units              = "(categorical / integer-coded)",
      type               = "categorical",
      reference_category = "0 (control non-TZD)",
      notes              = paste(
        "Time-fixed per subject; carries the randomized-treatment",
        "assignment. Selects the drug-effect parameters of the",
        "time-driven Emax model (Stringer 2015 Eq. 5,",
        "DEF = Emax * TIME / (ET50 + TIME)): pioglitazone group has",
        "Emax = 17.3% and ET50 = 0 days FIX (effect essentially",
        "instantaneous by the 3-month first FPG sample); control",
        "group has Emax = 8.4% and ET50 = 49.2 days. The Emax IIV",
        "(omega^2 = 0.75) is shared between the two cohorts."
      ),
      source_name        = "TRT"
    )
  )

  population <- list(
    species                 = "human",
    n_subjects              = 587L,
    n_pio                   = 293L,
    n_control               = 294L,
    n_studies               = 1L,
    study_design            = "Multicenter, prospective, randomized, open-label, blinded-end-point study in Japanese T2DM patients (2.5-4 year treatment; median 3.14 years, maximum 3.9 years).",
    age_range               = "35.0-74.0 years (median 58.0; pooled cohorts)",
    weight_range            = "44.0-116.0 kg (median 68-69 kg)",
    bmi_range               = "18.5-42.6 kg/m^2 (median 26.2-26.5 kg/m^2)",
    sex_female_pct          = 37.8,
    disease_state           = "Type 2 diabetes mellitus (T2DM) with baseline HbA1c > 6.9% (NGSP); on non-TZD oral glucose-lowering background therapy; treatment adjusted to achieve HbA1c < 6.9% (NGSP)",
    baseline_fpg_median     = "153-157 mg/dL (pio: 153.0 [77.0-304.0]; control: 157.0 [81.0-371.0])",
    baseline_hba1c_median   = "7.6-7.9 % NGSP (pio: 7.9 [6.9-11.4]; control: 7.6 [6.9-11.8])",
    duration_of_diabetes    = "< 5 years in ~28%, > 5 years in ~72%",
    co_medication           = "Sulfonylureas (73.0-81.6%), alpha-glucosidase inhibitors (35.8-55.8%), biguanides (42.6-67.7%), rapid-acting insulin secretagogues (6.5-12.9%). No other TZD.",
    regions                 = "Japan",
    n_dropped_before_2p5yr  = 90L,
    notes                   = "Baseline demographics from Stringer 2015 Table 1. No PK data collected. The primary safety / tolerability results of the parent cardiovascular-outcome study are reported in the companion publications (Stringer 2015 refs 21, 22)."
  )

  ini({
    # ---- FPG structural parameters (Stringer 2015 Table 2) ----
    # Baseline FPG typical value is reported for females; males are +5% higher via sex_effect_fpg.
    lrbase_fpg <- log(156.0); label("Typical baseline FPG for females (mg/dL)")  # Table 2 BSL FPG (females) = 156.0
    lkout_fpg  <- log(0.0089); label("First-order FPG elimination rate constant KoutG (1/day)")  # Table 2 KoutG = 0.0089 days^-1

    # ---- HbA1c structural parameters (Stringer 2015 Table 2) ----
    # Baseline HbA1c is estimated with a Box-Cox (lambda ~ 3.28) transformation on the eta
    # for skewness; this implementation uses log-normal IIV as a practical approximation
    # (see the vignette 'Assumptions and deviations' section).
    lrbase_hba1c <- log(7.83); label("Typical baseline HbA1c (%)")  # Table 2 BSL HbA1c = 7.83
    lkout_hba1c  <- log(0.072); label("First-order HbA1c elimination rate constant KoutH (1/day)")  # Table 2 KoutH = 0.072 days^-1
    lkinz        <- log(0.28); label("Log zero-order FPG-independent input rate KinZ to HbA1c (% HbA1c / day)")  # Table 2 K in ZT = 0.28 days^-1
    c_exp        <- 1.91; label("Power exponent on FPG in HbA1c production (dimensionless)")  # Table 2 c = 1.91

    # ---- Drug effect model (Stringer 2015 Eq. 5; Table 2) ----
    # DEF = Emax * TIME / (ET50 + TIME) stimulates KoutG. Emax and ET50 are treatment-specific.
    lemax_pio  <- log(0.173); label("Maximum FPG-lowering fractional effect Emax for pioglitazone (dimensionless, i.e. 17.3%)")  # Table 2 Emax (Pioglitazone) = 17.3%
    lemax_ctrl <- log(0.084); label("Maximum FPG-lowering fractional effect Emax for control cohort (dimensionless, i.e. 8.4%)")  # Table 2 Emax (Control) = 8.4%
    let50_ctrl <- log(49.2);  label("Time-to-half-maximal drug effect ET50 for control cohort (days)")  # Table 2 ET50 (Control) = 49.2 days
    et50_pio   <- fixed(0);   label("Time-to-half-maximal drug effect ET50 for pioglitazone (days; 0 -- effect essentially instantaneous by the 3-month first FPG sample)")  # Table 2 ET50 (Pioglitazone) = 0 FIX

    # ---- Disease progression rates (Stringer 2015 Table 2; parameters expressed in years^-1) ----
    # These are converted to days^-1 inside model() by dividing by 365.25.
    fpgdp <- 0.017; label("FPG disease-progression proportional slope (1/year)")  # Table 2 FPGDP = 0.017 years^-1
    dpind <- 0.03;  label("HbA1c FPG-independent input linear-in-time disease-progression slope (1/year)")  # Table 2 DPind = 0.03 years^-1

    # ---- Sex covariate on FPG baseline (Stringer 2015 Table 2 / Covariate analysis) ----
    sex_effect_fpg <- 0.05; label("Fractional increase in baseline FPG for males vs females (dimensionless; +5% for males, SEXF = 0)")  # Table 2 Gender on FPG BSL = 0.05

    # ---- Inter-individual variability (Stringer 2015 Table 2; omega^2 values reported) ----
    # OMEGA BLOCK on the two baseline etas: BSL FPG omega^2 = 0.03; BSL HbA1c omega^2 = 0.01;
    # off-diagonal reported as "Correlation (omega^2 BSL HbA1c, omega^2 BSL FPG) = 0.01"
    # interpreted as the covariance element of the NONMEM OMEGA BLOCK (a correlation
    # coefficient of ~0.577 between the two baseline etas).
    etalrbase_fpg + etalrbase_hba1c ~ c(0.03,
                                        0.01, 0.01)  # omega^2_BSL_FPG=0.03; cov=0.01; omega^2_BSL_HbA1c=0.01
    etafpgdp   ~ 0.004  # omega^2_FPGDP = 0.004 (additive random effect on FPGDP)
    etalemax   ~ 0.75   # omega^2_Emax = 0.75 (log-normal; shared eta across both treatment cohorts)

    # ---- Residual error (Stringer 2015 Table 2; proportional model) ----
    propSd_fpg   <- 0.144; label("Proportional residual error for FPG (14.4% CV)")   # Table 2 residual error FPG = 14.4%
    propSd_hba1c <- 0.058; label("Proportional residual error for HbA1c (5.8% CV)")  # Table 2 residual error HbA1c = 5.8%
  })

  model({
    # ---- Time-unit conversion constants ----
    # KoutG, KoutH, ET50, KinZ are reported in days^-1 / days; FPGDP and DPind are reported
    # in years^-1. Model TIME is in days; the two annualized rates are scaled here.
    days_per_year <- 365.25

    # ---- Individual baselines ----
    # Baseline FPG: female typical value, +5% higher for males via SEXF.
    rbase_fpg_i <- exp(lrbase_fpg + etalrbase_fpg) *
                   (1 + sex_effect_fpg * (1 - SEXF))
    # Baseline HbA1c: log-normal IIV approximation of the paper's Box-Cox (lambda = 3.28)
    # eta transformation (see vignette Assumptions and deviations).
    rbase_hba1c_i <- exp(lrbase_hba1c + etalrbase_hba1c)

    # ---- Individual FPG disease-progression slope (additive random effect) ----
    fpgdp_i_year <- fpgdp + etafpgdp
    fpgdp_i_day  <- fpgdp_i_year / days_per_year
    dpind_day    <- dpind / days_per_year

    # ---- Individual Emax (per treatment cohort, shared log-normal eta) ----
    emax_pio_i  <- exp(lemax_pio  + etalemax)
    emax_ctrl_i <- exp(lemax_ctrl + etalemax)
    emax_i      <- emax_pio_i * TRT + emax_ctrl_i * (1 - TRT)

    # ---- Individual ET50 (pioglitazone fixed to 0; control on log scale) ----
    et50_i <- exp(let50_ctrl) * (1 - TRT) + et50_pio * TRT

    # ---- Rate constants ----
    kout_fpg   <- exp(lkout_fpg)
    kout_hba1c <- exp(lkout_hba1c)
    kinz       <- exp(lkinz)

    # ---- Derive KinH from HbA1c steady state at t = 0, no drug, no disease progression:
    # 0 = KinZ + KinH * BSL_FPG^c - KoutH * BSL_HbA1c
    #  => KinH = (KoutH * BSL_HbA1c - KinZ) / BSL_FPG^c
    # (Stringer 2015 Table 2 does not report KinH directly; it is back-calculated per subject
    # so each subject's HbA1c initial condition matches its cascading indirect-response steady
    # state given the correlated (BSL_FPG_i, BSL_HbA1c_i) draw.)
    kinh_i <- (kout_hba1c * rbase_hba1c_i - kinz) / (rbase_fpg_i^c_exp)

    # ---- FPG zero-order production with proportional disease progression (Eq. 1) ----
    kingdp <- rbase_fpg_i * kout_fpg * (1 + fpgdp_i_day * t)

    # ---- Drug effect (Eq. 5); 1e-10 avoids 0/0 for pio at t = 0 with ET50 = 0 ----
    def_effect <- emax_i * t / (et50_i + t + 1e-10)

    # ---- FPG-independent zero-order input to HbA1c with linear-in-time DP (Eq. 4) ----
    fpgind <- kinz * (1 + dpind_day * t)

    # ---- ODEs (Eqs. 2 and 3) ----
    d/dt(fpg)   <- kingdp - kout_fpg * (1 + def_effect) * fpg
    d/dt(hba1c) <- fpgind + kinh_i * (fpg^c_exp) - kout_hba1c * hba1c

    fpg(0)   <- rbase_fpg_i
    hba1c(0) <- rbase_hba1c_i

    # ---- Observations and residual error ----
    fpg   ~ prop(propSd_fpg)
    hba1c ~ prop(propSd_hba1c)
  })
}
