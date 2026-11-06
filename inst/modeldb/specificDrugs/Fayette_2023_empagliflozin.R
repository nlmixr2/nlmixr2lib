Fayette_2023_empagliflozin <- function() {
  description <- paste(
    "Two-compartment population PK with lagged first-order absorption for",
    "empagliflozin in adults with type 2 diabetes, coupled with turnover",
    "models for fasting plasma glucose (FPG) and glycated hemoglobin",
    "(HbA1c). Steady-state AUC (dose / CL/F) stimulates FPG removal through",
    "an Emax function; the FPG ratio to baseline raised to a power drives",
    "HbA1c production, and a decaying placebo term stimulates HbA1c removal.",
    "Refit of the Baron 2016 PK/PD model to five phase II trials (PK, n =",
    "1127) and two phase II trials (PD, n = 750), used to build Bayesian",
    "borrowing priors for the 24-week phase III target trial (Fayette 2023)."
  )
  reference <- paste(
    "Fayette L, Sailer MO, Perez-Pitarch A. Pharmacometrics enhanced",
    "Bayesian borrowing approach to improve clinical trial efficiency: Case",
    "of empagliflozin in type 2 diabetes. CPT Pharmacometrics Syst Pharmacol.",
    "2023;12(10):1386-1397. doi:10.1002/psp4.13035"
  )
  vignette <- "Fayette_2023_empagliflozin"
  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "Cc in nmol/L (converted from mg/L via MW 450.9 g/mol); glucose (FPG) in mmol/L; hba1c in % (NGSP)"
  )
  paper_specific_compartments <- c("hba1c")

  compartmentData <- list(
    depot = list(analyte = "empagliflozin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "empagliflozin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "empagliflozin", units = "mg", specimen = "plasma", verified = TRUE),
    glucose = list(analyte = "glucose", units = "mmol/L", specimen = "plasma", verified = TRUE),
    hba1c = list(analyte = "HbA1c", units = "%", specimen = "whole blood", verified = TRUE)
  )

  covariateData <- list(
    AGE = list(
      description = "Age at baseline",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Power-form effects centred at 50 years on CL/F, V2/F, V3/F and ka (Appendix S2 Table S4). The centring value is not printed in Fayette 2023; it is taken from the predecessor model (Baron 2016 Table S1).",
      source_name = "Age"
    ),
    BMI = list(
      description = "Body mass index at baseline",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Power-form effects on CL/F and V3/F centred at 25 kg/m^2 and on V2/F centred at 20 kg/m^2 (Table S4; centring values from Baron 2016 Table S1).",
      source_name = "BMI"
    ),
    SEXF = list(
      description = "Sex (female indicator)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Exponential effects exp(theta * SEXF) on CL/F, V2/F, V3/F and ka (Table S4). The paper's SEX column is read as a female indicator: every sign matches the female multipliers of Baron 2016 (e.g. CL/F exp(-0.127) = 0.88 vs 0.886).",
      source_name = "SEX"
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-Asian)",
      notes = "Exponential effects on CL/F, V2/F, V3/F, ka (Table S4) and on baseline FPG (Table S6 BFPGASIAN1).",
      source_name = "ASIAN race"
    ),
    SMOKE_CURRENT = list(
      description = "Current-smoker indicator (paired with SMOKE_NEVER)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = former smoker (when paired with SMOKE_NEVER = 0)",
      notes = "Table S4 reports CL/F effects for ex-smoker (0.0146) and current smoker (0.056) against a never-smoker reference. The model re-centres them on the canonical former-smoker reference: lcl absorbs +0.0146, e_smoke_never_cl = -0.0146 and e_smoke_current_cl = 0.056 - 0.0146, which is exact on the log scale.",
      source_name = "Current smoker"
    ),
    SMOKE_NEVER = list(
      description = "Never-smoker indicator (paired with SMOKE_CURRENT)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = former smoker (when paired with SMOKE_CURRENT = 0)",
      notes = "See SMOKE_CURRENT for the re-centring from the paper's never-smoker reference.",
      source_name = "Ex-smoker / Current smoker (never smoker = reference)"
    ),
    CRCL = list(
      description = "MDRD estimated glomerular filtration rate, BSA-normalised",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Power-form effect on CL/F centred at 100 mL/min/1.73 m^2 (Table S4; centring from Baron 2016). Values above 120 are capped at 120 inside the model, as the paper did for the covariate analysis (Table S2 footnote b).",
      source_name = "eGFR"
    ),
    TPRO = list(
      description = "Total serum protein",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Power-form effects centred at 70 g/L on CL/F, V2/F and V3/F (Table S4; centring from Baron 2016 Table S1, which prints the unit as g/dL although 70 is the g/L value).",
      source_name = "Total protein"
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Power-form effect on CL/F centred at 20 U/L (Table S4; centring from Baron 2016).",
      source_name = "ALT"
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Power-form effect on CL/F centred at 20 U/L (Table S4; centring from Baron 2016).",
      source_name = "AST"
    ),
    ALP = list(
      description = "Alkaline phosphatase",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Power-form effect on CL/F centred at 70 U/L (Table S4; centring from Baron 2016).",
      source_name = "AP"
    ),
    LDH = list(
      description = "Lactate dehydrogenase",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Power-form effect on CL/F centred at 160 U/L (Table S4; centring from Baron 2016).",
      source_name = "LDH"
    ),
    FPG = list(
      description = "Observed baseline fasting plasma glucose",
      units = "mmol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Uncentred exponential covariate on kHbA1c,in: exp(0.04 * FPG) (Table S6 KHINFPGB1; Appendix S1 Eq. for TVP). Distinct from the model-predicted baseline rbase_glucose, which sets the FPG state's initial condition.",
      source_name = "Observed baseline for FPG"
    ),
    DOSE_EMPA_MGD = list(
      description = "Once-daily empagliflozin dose",
      units = "mg/day",
      type = "continuous",
      reference_category = NULL,
      notes = "Drives the steady-state exposure AUCss = 2217.8 * DOSE / (CL/F) in nmol*h/L (Appendix S1). Set to 0 for placebo.",
      source_name = "DOSE"
    )
  )

  covariatesDataExcluded <- list(
    RACE_BLACK = list(
      description = "Black race indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-Black)",
      notes = "Table S4 lists BLACK effects on CL/F, V2/F, V3/F and ka, each reported as 0 with no standard error (the coefficients were held at zero), so the final model carries no Black-race effect.",
      source_name = "BLACK"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 1127L,
    n_subjects_pd = 750L,
    n_studies = 5L,
    n_pk_observations = 13693L,
    n_fpg_observations = 4036L,
    n_hba1c_observations = 4581L,
    age_mean_sd = "57.7 +/- 9.2 years (PK dataset); 58.0 +/- 9.3 years (PD dataset)",
    weight_mean_sd = "84.9 +/- 17.3 kg (PK); 85.8 +/- 17.1 kg (PD)",
    bmi_mean_sd = "29.8 +/- 4.8 kg/m^2 (PK); 30.3 +/- 4.7 kg/m^2 (PD)",
    sex_female_pct = 42.2,
    race_ethnicity = c(NonAsian = 78.6, Asian = 21.4),
    smoking = c(Never = 58.7, Former = 24.4, Current = 17.0),
    egfr_mean_sd = "90.7 +/- 18.6 mL/min/1.73 m^2 (MDRD; PK dataset)",
    disease_state = "Type 2 diabetes; 66.8% drug naive in the PK dataset, the rest on background oral antidiabetics.",
    dose_range = "Empagliflozin 1-100 mg once daily (PK); 1-50 mg once daily (PD), up to 12 weeks, plus placebo.",
    studies = "Five phase II trials (Heise 2013 x2, Kanada 2013, Ferrannini 2013, Rosenstock 2013); the PD model used the two with HbA1c and FPG (Table S1).",
    notes = "Demographics from Appendix S2 Table S2. The model was then used to simulate the phase III target trial (Roden 2013; drug-naive adults, placebo or 10 / 25 mg once daily for 24 weeks)."
  )

  ini({
    # ---- Structural PK (Table S4 point estimates; reported on the log scale) ----
    # Reference subject: 50-year-old non-Asian male former smoker, BMI 25
    # (20 for V2/F), eGFR 100, total protein 70 g/L, ALT 20, AST 20, AP 70,
    # LDH 160. The paper's CL/F reference is a never smoker; lcl adds the
    # ex-smoker coefficient so the canonical former-smoker reference holds.
    lcl <- 2.44 + 0.0146; label("Apparent clearance CL/F (log L/h), former-smoker reference") # Table S4 'log CL' = 2.44 plus 'Ex-smoker' = 0.0146
    lvc <- 1.19; label("Apparent central volume V2/F (log L)") # Table S4 'log V2' = 1.19
    lq <- 1.84; label("Apparent intercompartmental clearance Q/F (log L/h)") # Table S4 'log Q' = 1.84
    lvp <- 4.32; label("Apparent peripheral volume V3/F (log L)") # Table S4 'log V3' = 4.32
    lka <- -1.63; label("Absorption rate constant ka (log 1/h)") # Table S4 'log KA' = -1.63
    lalag <- fixed(log(0.5)); label("Absorption lag time (h)") # Table S4 ALAG1 = 0.5 h, footnote 'ALAG1 was fixed at 0.5 h'

    # ---- Covariate effects on CL/F (Table S4) ----
    e_age_cl <- -0.305; label("Power exponent of AGE/50 on CL/F") # Table S4 CL 'Age' = -0.305
    e_bmi_cl <- 0.249; label("Power exponent of BMI/25 on CL/F") # Table S4 CL 'BMI' = 0.249
    e_sexf_cl <- -0.127; label("Log-scale effect of female sex on CL/F") # Table S4 CL 'SEX' = -0.127
    e_smoke_never_cl <- -0.0146; label("Log-scale effect of never smoking (vs former) on CL/F") # Table S4 CL 'Ex-smoker' = 0.0146 vs never, re-centred
    e_smoke_current_cl <- 0.056 - 0.0146; label("Log-scale effect of current smoking (vs former) on CL/F") # Table S4 CL 'Current smoker' = 0.056 minus 'Ex-smoker' = 0.0146
    e_crcl_cl <- 0.27; label("Power exponent of eGFR/100 on CL/F") # Table S4 CL 'eGFR' = 0.27
    e_tpro_cl <- -0.242; label("Power exponent of TPRO/70 on CL/F") # Table S4 CL 'Total protein' = -0.242
    e_alt_cl <- -0.0132; label("Power exponent of ALT/20 on CL/F") # Table S4 CL 'ALT' = -0.0132
    e_ast_cl <- -0.0533; label("Power exponent of AST/20 on CL/F") # Table S4 CL 'AST' = -0.0533
    e_alp_cl <- -0.0325; label("Power exponent of ALP/70 on CL/F") # Table S4 CL 'AP' = -0.0325
    e_ldh_cl <- 0.0298; label("Power exponent of LDH/160 on CL/F") # Table S4 CL 'LDH' = 0.0298
    e_race_asian_cl <- -0.187; label("Log-scale effect of Asian race on CL/F") # Table S4 CL 'ASIAN race' = -0.187

    # ---- Covariate effects on V2/F (Table S4) ----
    e_age_vc <- 0.967; label("Power exponent of AGE/50 on V2/F") # Table S4 V2 'Age' = 0.967
    e_sexf_vc <- 0.26; label("Log-scale effect of female sex on V2/F") # Table S4 V2 'SEX' = 0.26
    e_tpro_vc <- 2.03; label("Power exponent of TPRO/70 on V2/F") # Table S4 V2 'Total protein' = 2.03
    e_bmi_vc <- 0.618; label("Power exponent of BMI/20 on V2/F") # Table S4 V2 'BMI' = 0.618
    e_race_asian_vc <- 0.0137; label("Log-scale effect of Asian race on V2/F") # Table S4 V2 'ASIAN race' = 0.0137

    # ---- Covariate effects on V3/F (Table S4) ----
    e_age_vp <- -0.065; label("Power exponent of AGE/50 on V3/F") # Table S4 V3 'Age' = -0.065
    e_sexf_vp <- -0.117; label("Log-scale effect of female sex on V3/F") # Table S4 V3 'SEX' = -0.117
    e_tpro_vp <- -0.196; label("Power exponent of TPRO/70 on V3/F") # Table S4 V3 'Total protein' = -0.196
    e_bmi_vp <- 0.36; label("Power exponent of BMI/25 on V3/F") # Table S4 V3 'BMI' = 0.36
    e_race_asian_vp <- -0.0609; label("Log-scale effect of Asian race on V3/F") # Table S4 V3 'ASIAN race' = -0.0609

    # ---- Covariate effects on ka (Table S4 cont.) ----
    e_age_ka <- 0.0845; label("Power exponent of AGE/50 on ka") # Table S4 KA 'Age' = 0.0845
    e_sexf_ka <- 0.164; label("Log-scale effect of female sex on ka") # Table S4 KA 'SEX' = 0.164
    e_race_asian_ka <- 0.225; label("Log-scale effect of Asian race on ka") # Table S4 KA 'ASIAN race' = 0.225

    # ---- PK IIV (Table S4; omega^2 with CV% beside it) ----
    etalcl + etalvp ~ c(
      0.106,
      0.0277, 0.0652
    ) # Table S4 omega2 CL = 0.106 (33.4% CV), omega CL,V3 = 0.0277, omega2 V3 = 0.0652 (26% CV)
    etalka ~ 0.0338 # Table S4 omega2 KA = 0.0338 (18.5% CV)

    # ---- PK residual error (Table S4) ----
    # Log-transformed concentrations; the reported CVs equal
    # sqrt(exp(sigma^2) - 1), so the error is log-normal. Studies 9-10 (the
    # sparsely sampled 12-week phase IIb trials) are encoded; Studies
    # 2-4-15 (intensive sampling) had sigma^2 = 0.0278 (16.8% CV).
    expSd <- sqrt(0.137); label("Log-normal residual SD for Cc, Studies 9-10") # Table S4 sigma2 Studies 9-10 = 0.137 (38.3% CV)

    # ---- FPG turnover (Table S6) ----
    lrbase_glucose <- log(9.61); label("Typical baseline FPG, non-Asian (mmol/L)") # Table S6 'Baseline FPG' = 9.61
    e_race_asian_rbase_glucose <- -0.0991; label("Log-scale effect of Asian race on baseline FPG") # Table S6 BFPGASIAN1 = -0.0991
    lkout_glucose <- fixed(log(0.0407)); label("First-order FPG removal rate constant kFPG,out (1/h)") # Table S6 'FPG elimination' = 0.0407, no RSE or CI (held fixed)
    lemax <- log(0.201); label("Maximal proportional stimulation of FPG removal GMAX (unitless)") # Table S6 GMAX = 0.201
    lauc50 <- log(718); label("Steady-state AUC giving half of GMAX (nmol*h/L)") # Table S6 GC50 = 718

    # ---- HbA1c turnover (Table S6) ----
    hba1c_limit <- 4.56; label("HbA1c boundary condition HbA1c,limit (%)") # Table S6 'HbA1c limit' = 4.56
    lkout_hba1c <- log(0.000641); label("First-order HbA1c removal rate constant kHbA1c,out (1/h)") # Table S6 'HbA1c elimination' = 0.000641
    lhba1c_ratio <- log(2.23); label("Ratio kHbA1c,in / kHbA1c,out at FPG = 0 (%)") # Table S6 'Ratio KHIN/KHOT' = 2.23
    e_fpg_hba1c_ratio <- 0.04; label("Exponential effect of observed baseline FPG on kHbA1c,in (L/mmol)") # Table S6 KHINFPGB1 = 0.04
    gamma_hba1c <- 1.21; label("Power of FPG / baseline FPG driving HbA1c production") # Table S6 POWER = 1.21
    pbo_hba1c <- 0.173; label("Typical amplitude of the placebo stimulation of HbA1c removal") # Table S6 PLCBM = 0.173
    lkpbo_hba1c <- log(0.000422); label("Decay rate constant of the placebo effect on HbA1c (1/h)") # Table S6 PLCBT = 0.000422

    # ---- PD IIV (Table S6) ----
    etalrbase_glucose + etalemax ~ c(
      0.0308,
      0.0519, 0.431
    ) # Table S6 ETA1 BFPG = 0.0308 (17.7% CV), Corr1 BFPG,GMAX = 0.0519 (covariance), ETA2 GMAX = 0.431 (73.4% CV)
    etapbo_hba1c ~ 0.431 # Table S6 ETA3 pcbo = 0.431; additive eta on the placebo amplitude
    etalhba1c_ratio ~ 0.0244 # Table S6 ETA4 kHbA1c,in = 0.0244 (15.7% CV)

    # ---- PD residual error (Table S6; exponential) ----
    expSd_glucose <- sqrt(0.0166); label("Log-normal residual SD for FPG") # Table S6 EPS1 = 0.0166 (12.9% CV)
    expSd_hba1c <- sqrt(0.00169); label("Log-normal residual SD for HbA1c") # Table S6 EPS2 = 0.00169 (4.1% CV)
  })

  model({
    # Empagliflozin MW 450.9 g/mol; 2217.8 nmol per mg (Appendix S1).
    nmol_per_mg <- 2217.8

    # Covariate centring values (not printed in Fayette 2023; taken from the
    # predecessor model, Baron 2016 Table S1).
    ref_age <- 50
    ref_bmi <- 25
    ref_bmi_vc <- 20
    ref_crcl <- 100
    ref_tpro <- 70
    ref_alt <- 20
    ref_ast <- 20
    ref_alp <- 70
    ref_ldh <- 160

    # eGFR above 120 mL/min/1.73 m^2 was set to 120 (Table S2 footnote b).
    crcl_cov <- min(CRCL, 120)

    # ---- Individual PK parameters ----
    cl <- exp(
      lcl +
        etalcl +
        e_sexf_cl * SEXF +
        e_race_asian_cl * RACE_ASIAN +
        e_smoke_never_cl * SMOKE_NEVER +
        e_smoke_current_cl * SMOKE_CURRENT
    ) *
      (AGE / ref_age)^e_age_cl *
      (BMI / ref_bmi)^e_bmi_cl *
      (crcl_cov / ref_crcl)^e_crcl_cl *
      (TPRO / ref_tpro)^e_tpro_cl *
      (ALT / ref_alt)^e_alt_cl *
      (AST / ref_ast)^e_ast_cl *
      (ALP / ref_alp)^e_alp_cl *
      (LDH / ref_ldh)^e_ldh_cl
    vc <- exp(lvc + e_sexf_vc * SEXF + e_race_asian_vc * RACE_ASIAN) *
      (AGE / ref_age)^e_age_vc *
      (TPRO / ref_tpro)^e_tpro_vc *
      (BMI / ref_bmi_vc)^e_bmi_vc
    vp <- exp(lvp + etalvp + e_sexf_vp * SEXF + e_race_asian_vp * RACE_ASIAN) *
      (AGE / ref_age)^e_age_vp *
      (TPRO / ref_tpro)^e_tpro_vp *
      (BMI / ref_bmi)^e_bmi_vp
    q <- exp(lq)
    ka <- exp(lka + etalka + e_sexf_ka * SEXF + e_race_asian_ka * RACE_ASIAN) *
      (AGE / ref_age)^e_age_ka

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    alag(depot) <- exp(lalag)

    Cc <- central / vc * nmol_per_mg

    # ---- Exposure driver (Appendix S1): AUCss = 2217.8 * DOSE / (CL/F) ----
    auc_ss <- nmol_per_mg * DOSE_EMPA_MGD / cl

    # ---- FPG turnover (Appendix S1 Eq. 1) ----
    # Removal is first order in FPG; kFPG,in = BFPG * kFPG,out keeps the
    # pre-treatment steady state.
    rbase_glucose <- exp(lrbase_glucose + etalrbase_glucose + e_race_asian_rbase_glucose * RACE_ASIAN)
    kout_glucose <- exp(lkout_glucose)
    kin_glucose <- rbase_glucose * kout_glucose
    emax <- exp(lemax + etalemax)
    auc50 <- exp(lauc50)
    stim <- emax * auc_ss / (auc50 + auc_ss)
    d/dt(glucose) <- kin_glucose - kout_glucose * glucose * (1 + stim)
    glucose(0) <- rbase_glucose

    # ---- HbA1c turnover (Appendix S1 Eq. 2) ----
    # kHbA1c,in = ratio * kHbA1c,out, with an uncentred exponential effect of
    # observed baseline FPG on the ratio; baseline HbA1c = ratio + limit.
    kout_hba1c <- exp(lkout_hba1c)
    hba1c_ratio <- exp(lhba1c_ratio + etalhba1c_ratio + e_fpg_hba1c_ratio * FPG)
    kin_hba1c <- hba1c_ratio * kout_hba1c
    kpbo_hba1c <- exp(lkpbo_hba1c)
    # The placebo amplitude switches on after the first dose (time 0).
    pbo_i <- (pbo_hba1c + etapbo_hba1c) * (time > 0)
    d/dt(hba1c) <- kin_hba1c * (glucose / rbase_glucose)^gamma_hba1c -
      kout_hba1c * (hba1c - hba1c_limit) * (1 + pbo_i * exp(-kpbo_hba1c * time))
    hba1c(0) <- hba1c_ratio + hba1c_limit

    Cc ~ lnorm(expSd)
    glucose ~ lnorm(expSd_glucose)
    hba1c ~ lnorm(expSd_hba1c)
  })
}
