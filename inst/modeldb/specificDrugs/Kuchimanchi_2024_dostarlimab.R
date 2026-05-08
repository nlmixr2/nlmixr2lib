Kuchimanchi_2024_dostarlimab <- function() {
  description <- "Two-compartment population PK model for dostarlimab (anti-PD-1 IgG4) with sigmoid I_max time-dependent clearance, fitted to GARNET (advanced solid tumours) plus RUBY Part 1 (primary advanced or recurrent endometrial cancer with carboplatin-paclitaxel) data (Kuchimanchi 2024)"
  reference <- "Kuchimanchi M, Jorgensen TL, Hanze E, et al. Population pharmacokinetics and exposure-response relationships of dostarlimab in primary advanced or recurrent endometrial cancer in part 1 of RUBY. Br J Clin Pharmacol. 2025;91(3):841-855. doi:10.1111/bcp.16325"
  vignette <- "Kuchimanchi_2024_dostarlimab"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight at baseline",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Baseline body weight. Allometric power scaling (WT / 70)^exponent on CL (exponent 0.523) and on Vc/Vp (shared exponent 0.48). Reference 70 kg (Kuchimanchi 2024 Methods, reference patient).",
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Age at baseline",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling (AGE / 64)^-0.238 on CL. Reference 64 years (median age of analysis set; Kuchimanchi 2024 Methods, reference patient).",
      source_name        = "AGE"
    ),
    ALB = list(
      description        = "Serum albumin (time-varying)",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying covariate. Power scaling (ALB / 39)^-0.922 on CL and (ALB / 39)^-0.132 on Vc. Reference 39 g/L (median of analysis set, Kuchimanchi 2024 Table 1; Methods reference patient).",
      source_name        = "ALB"
    ),
    ALT = list(
      description        = "Alanine aminotransferase (time-varying)",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying covariate. Power scaling (ALT / 18)^-0.0623 on CL. Reference 18 U/L per Methods reference-patient text and the equation denominator. Note the publication has an internal inconsistency: Table 1 reports median ALT of 17 U/L for the analysis set and the Figure 2 caption likewise lists ALT = 17 U/L; the Methods text and the structural equation use 18 U/L. The model uses 18 to match the equation; using 17 instead changes typical CL by less than 0.4%.",
      source_name        = "ALT"
    ),
    SEXF = list(
      description        = "Biological sex indicator (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 (female; the most common group, 82.0% of the analysis set)",
      notes              = "Female is the reference category (Kuchimanchi 2024: theta_CL_SEX and theta_Vc_SEX are 0 for females and estimated for males). Implemented as (1 + theta * (1 - SEXF)) so SEXF = 1 (female) gives factor 1 and SEXF = 0 (male) gives factor 1 + theta.",
      source_name        = "SEX (paper codes a male indicator; mapped to canonical SEXF via 1 - SEXM)"
    ),
    CONMED_CHEMO = list(
      description        = "Coadministration regimen: dostarlimab + platinum-doublet chemotherapy (carboplatin AUC 5 + paclitaxel 175 mg/m^2 Q3W)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (dostarlimab monotherapy; the reference regimen in this analysis)",
      notes              = "Indicator for whether dostarlimab is being given as part of combination chemotherapy (RUBY Part 1 cycles 1-6) versus monotherapy (GARNET, RUBY cycles 7+). Implemented multiplicatively as (1 + e_combo_cl * CONMED_CHEMO) with e_combo_cl = -0.0779, so combination therapy yields CL factor 0.9221 (~7.79% lower CL than monotherapy), matching Kuchimanchi 2024 Abstract and Table 2. The paper's Equation block writes the effect with an implicit indicator as '(1 - theta_CL_MONOTR)'; the most direct numerical match to the abstract phrasing 'CL was 7.79% lower in combination therapy' is achieved by treating the parameter as a fractional-change coefficient on a combination-therapy indicator.",
      source_name        = "MONOTR (the paper's monotherapy/combination indicator; mapped to canonical CONMED_CHEMO via the inverse value relation, see notes)"
    )
  )

  population <- list(
    n_subjects     = 868L,
    n_studies      = 2L,
    n_observations = 7957L,
    age_range      = "24-86 years",
    age_median     = "64 years",
    weight_range   = "34.0-182.0 kg",
    weight_median  = "73.0 kg",
    sex_female_pct = 82.0,
    race_ethnicity = c(White = 74.6, Black_AfricanAmerican = 5.4, Asian = 2.3,
                       AmericanIndianAlaskaNative = 0.6, NativeHawaiianPacificIslander = 0.1,
                       Other = 0.7, Unknown = 1.8, NotReported = 14.5),
    disease_state  = "Primary advanced or recurrent endometrial cancer (RUBY Part 1) plus advanced or recurrent solid tumours (GARNET: dMMR/MSI-H endometrial cancer, MMRp/MSS endometrial cancer, non-EC dMMR/MSI-H or POLE-mutant tumours, and NSCLC).",
    dose_range     = "GARNET: 1, 3 or 10 mg/kg IV Q2W (Part 1); 500 mg Q3W or 1000 mg Q6W (Part 2A); 500 mg Q3W x 4 cycles followed by 1000 mg Q6W (Part 2B). RUBY Part 1: 500 mg IV Q3W with carboplatin (AUC 5) + paclitaxel (175 mg/m^2) Q3W for 6 cycles, followed by 1000 mg IV Q6W for up to 3 years.",
    regions        = "Multinational; Europe 52.5%, North America 47.5%",
    coadministration = "Dostarlimab monotherapy in GARNET (n = 636) and in RUBY Part 1 cycles 7+; dostarlimab + carboplatin-paclitaxel in RUBY Part 1 cycles 1-6 (n = 233).",
    ada_positive_pct = 11.6,
    hepatic_impairment_pct = c(none = 88.8, mild = 10.6, moderate = 0.6, severe = 0),
    renal_impairment_pct   = c(none = 35.1, mild = 45.7, moderate = 18.9, severe = 0.3),
    notes          = "Baseline demographics from Kuchimanchi 2024 Table 1 (overall analysis set, N = 869; one patient excluded from PopPK so analysis set was 868 with 7,957 PK observations). Reference patient (Methods, 'PopPK model development'): female, 70 kg, 64 years, albumin 39 g/L, ALT 18 U/L (Methods text; Figure 2 caption reports ALT = 17 U/L)."
  )

  ini({
    # Structural PK parameters at the reference patient (female, 70 kg, 64 y,
    # ALB 39 g/L, ALT 18 U/L, monotherapy). Values are typical-value
    # population-PK parameter estimates from Kuchimanchi 2024 Table 2. Source
    # values reported in L/h (CL, Q) are converted to L/day so the time axis
    # matches the rest of the dostarlimab / mAb library and the time-dependent
    # CL equation, which is parameterised in days (T50 = 145 days).
    lcl  <- log(0.00732 * 24); label("Baseline clearance at t=0 for reference patient (CL_base, L/day)") # Kuchimanchi 2024 Table 2: CL = 0.00732 L/h
    lvc  <- log(3.09);         label("Central volume of distribution for reference patient (Vc_base, L)") # Kuchimanchi 2024 Table 2: Vc = 3.09 L
    lq   <- log(0.0191 * 24);  label("Intercompartmental clearance (Q, L/day)")                           # Kuchimanchi 2024 Table 2: Q = 0.0191 L/h
    lvp  <- log(2.48);         label("Peripheral volume of distribution for 70 kg reference (Vp, L)")     # Kuchimanchi 2024 Table 2: Vp = 2.48 L

    # Time-dependent clearance (sigmoid I_max function of time since first
    # dose; Kuchimanchi 2024 Results, equation block):
    #   CL(t) = CL_base * exp(I_max * t^Hill / (T50^Hill + t^Hill))
    # I_max is negative (a fractional decrease in log-CL at t >> T50). To
    # keep individual I_max strictly negative under log-normal IIV, the
    # magnitude is stored as log|I_max| and the sign is applied in model()
    # (same pattern as Melhem 2022 dostarlimab). The resulting maximum
    # reduction at t >> T50 is 1 - exp(I_max) = 1 - exp(-0.113) = 0.107
    # (10.7%), matching the paper's narrative.
    lImax <- log(0.113); label("log|I_max|; magnitude of the log-CL reduction at t >> T50 (unitless)") # Kuchimanchi 2024 Table 2: Imax = -0.113
    lt50  <- log(145);   label("log T50; time at half of I_max (days)")                                # Kuchimanchi 2024 Table 2: T50 = 145 days
    lhill <- log(7.05);  label("log Hill; sigmoid steepness coefficient (unitless)")                   # Kuchimanchi 2024 Table 2: Hill = 7.05

    # Allometric exponents on body weight (reference 70 kg).
    e_wt_cl    <- 0.523;  label("Allometric exponent of WT on CL (unitless)")              # Kuchimanchi 2024 Table 2: Effect of WT on CL = 0.523
    e_wt_vc_vp <- 0.48;   label("Shared allometric exponent of WT on Vc and Vp (unitless)") # Kuchimanchi 2024 Table 2: Effect of WT on Vc and Vp = 0.48

    # Continuous covariate effects on CL (power form: (cov / ref)^theta).
    e_age_cl <- -0.238;   label("AGE exponent on CL (power form (AGE/64)^theta)")        # Kuchimanchi 2024 Table 2: Effect of age on CL = -0.238
    e_alb_cl <- -0.922;   label("ALB exponent on CL (power form (ALB/39)^theta)")        # Kuchimanchi 2024 Table 2: Effect of ALB on CL = -0.922
    e_alt_cl <- -0.0623;  label("ALT exponent on CL (power form (ALT/18)^theta)")        # Kuchimanchi 2024 Table 2: Effect of ALT on CL = -0.0623
    e_alb_vc <- -0.132;   label("ALB exponent on Vc (power form (ALB/39)^theta)")        # Kuchimanchi 2024 Table 2: Effect of ALB on Vc = -0.132

    # Categorical sex effect on CL and Vc (paper codes a male indicator with
    # female as the reference; encoded here as (1 + theta * (1 - SEXF))).
    e_sex_cl <- 0.15;  label("Male fractional change on CL vs. female reference")        # Kuchimanchi 2024 Table 2: Effect of male on CL = 0.15
    e_sex_vc <- 0.137; label("Male fractional change on Vc vs. female reference")        # Kuchimanchi 2024 Table 2: Effect of male on Vc = 0.137

    # Categorical combination-therapy effect on CL. Kuchimanchi 2024 Table 2
    # reports "Effect of combination therapy on CL = -0.0779"; the abstract
    # describes this as combo CL being 7.79% lower than monotherapy. Encoded
    # as (1 + e_combo_cl * CONMED_CHEMO): combination therapy
    # (CONMED_CHEMO = 1) gives factor 1 - 0.0779 = 0.9221.
    e_combo_cl <- -0.0779; label("Combination-therapy fractional change on CL vs. monotherapy reference") # Kuchimanchi 2024 Table 2: Effect of combination therapy on CL = -0.0779

    # IIV - block on (CL, Vc) and independent eta on |I_max|. Kuchimanchi 2024
    # Table 2 reports: omega^2_CL = 0.0563, cov(CL, Vc) = 0.0193,
    # omega^2_Vc = 0.0278, omega^2_Imax = 0.903.
    etalcl + etalvc ~ c(0.0563,
                        0.0193, 0.0278)
    etalImax ~ 0.903

    # Residual error (combined additive + proportional). Kuchimanchi 2024
    # Table 2 reports two proportional residual errors (GARNET 0.16, RUBY
    # 0.246) and a single additive component (4.22 mg/L). The model file
    # carries the GARNET proportional error because GARNET is the larger,
    # primary dataset and contributes the recommended-therapeutic-dose PK
    # observations underpinning the model; users wanting the RUBY error
    # should set propSd = 0.246. The additive value is reported in mg/L
    # (interpreted as SD on the linear concentration scale).
    propSd <- 0.16; label("Proportional residual error SD (fraction; GARNET cohort)") # Kuchimanchi 2024 Table 2: Proportional error, GARNET = 0.16
    addSd  <- 4.22; label("Additive residual error SD (mg/L)")                         # Kuchimanchi 2024 Table 2: Additive error = 4.22 mg/L
  })

  model({
    # Allometric weight scaling (reference 70 kg).
    wt_cl    <- (WT / 70)^e_wt_cl
    wt_vc_vp <- (WT / 70)^e_wt_vc_vp

    # Continuous covariate effects (power form, Kuchimanchi 2024 final-model
    # equations).
    age_cl <- (AGE / 64)^e_age_cl
    alb_cl <- (ALB / 39)^e_alb_cl
    alt_cl <- (ALT / 18)^e_alt_cl
    alb_vc <- (ALB / 39)^e_alb_vc

    # Sex effect: female reference (SEXF = 1 -> multiplier 1); male
    # (SEXF = 0) -> multiplier 1 + theta.
    sex_cl <- 1 + e_sex_cl * (1 - SEXF)
    sex_vc <- 1 + e_sex_vc * (1 - SEXF)

    # Combination-therapy effect: monotherapy reference (CONMED_CHEMO = 0
    # -> multiplier 1); combo (CONMED_CHEMO = 1) -> multiplier 1 + theta.
    combo_cl <- 1 + e_combo_cl * CONMED_CHEMO

    # Time-dependent CL (Hill function of time since first dose; t in days).
    # I_max < 0; sign applied here to keep individual values strictly negative
    # under log-normal IIV on |I_max|.
    Imax_i <- -exp(lImax + etalImax)
    t50    <- exp(lt50)
    hill   <- exp(lhill)
    td_cl  <- exp(Imax_i * t^hill / (t50^hill + t^hill))

    # Individual PK parameters. CL_base is at t=0 for the reference patient;
    # td_cl folds in the time dependency.
    cl_base <- exp(lcl + etalcl) * wt_cl * age_cl * alb_cl * alt_cl * sex_cl * combo_cl
    cl      <- cl_base * td_cl
    vc      <- exp(lvc + etalvc) * wt_vc_vp * alb_vc * sex_vc
    vp      <- exp(lvp)          * wt_vc_vp
    q       <- exp(lq)

    # Two-compartment micro-constants.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # Concentration: dose in mg / volume in L = mg/L.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
