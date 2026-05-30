Melhem_2022_dostarlimab <- function() {
  description <- "Two-compartment population PK model for dostarlimab (anti-PD-1 IgG4) with time-dependent (sigmoid I_max) clearance in adults with advanced solid tumours (Melhem 2022)"
  reference <- "Melhem M, Hanze E, Lu Z, Venkatakrishnan K, Gupta N, Vugmeyster Y. Population pharmacokinetics and exposure-response of anti-programmed cell death protein-1 monoclonal antibody dostarlimab in advanced solid tumours. Br J Clin Pharmacol. 2022;88(9):4142-4154. doi:10.1111/bcp.15339"
  vignette <- "Melhem_2022_dostarlimab"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric power scaling (WT / 70)^exponent on CL (exponent 0.470) and on Vc/Vp (exponent 0.419). Reference 70 kg.",
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Age at baseline",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling (AGE / 64)^-0.227 on CL. Reference 64 years (median of analysis set; reference patient described in Melhem 2022 Methods, 'PopPK model development').",
      source_name        = "AGE"
    ),
    ALB = list(
      description        = "Serum albumin (time-varying)",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying covariate. Power scaling (ALB / 39)^-1.01 on CL and (ALB / 39)^-0.153 on Vc. Reference 39 g/L (median of analysis set, Table 2; the Methods text writes 'g dL-1' but this is a typo — Table 2 lists albumin in g/L with median 39).",
      source_name        = "ALB"
    ),
    ALT = list(
      description        = "Alanine aminotransferase (time-varying)",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying covariate. Power scaling (ALT / 18)^-0.0585 on CL. Reference 18 U/L (median of analysis set, Table 2).",
      source_name        = "ALT"
    ),
    SEXF = list(
      description        = "Biological sex indicator (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 (female; the most common group, 77.3% of the analysis set)",
      notes              = "Female is the reference category in Melhem 2022 (theta_CL_SEX and theta_Vc_SEX are 0 for females and estimated for males). Implemented here as multiplicative (1 + theta * (1 - SEXF)) so that the male effect (SEXF = 0) gives 1 + theta and female (SEXF = 1) gives 1.",
      source_name        = "SEX (paper codes male indicator; mapped to canonical SEXF via 1 - SEXM)"
    )
  )

  population <- list(
    n_subjects     = 546L,
    n_studies      = 1L,
    age_range      = "24-86 years",
    age_median     = "64 years",
    weight_range   = "34.0-182.0 kg",
    weight_median  = "71.4 kg",
    sex_female_pct = 77.3,
    race_ethnicity = c(White = 75.1, Black_AfricanAmerican = 3.5, Asian = 2.4, Other = 1.1, Unknown = 0.9, NotReported = 16.3),
    disease_state  = "Advanced or recurrent solid tumours (mismatch repair-proficient/deficient endometrial cancer, NSCLC, non-EC dMMR/MSI-H or POLE-mutant tumours, ovarian cancer)",
    dose_range     = "1, 3 or 10 mg/kg IV Q2W (Part 1); 500 mg Q3W or 1000 mg Q6W (Part 2A); 500 mg Q3W x 4 cycles followed by 1000 mg Q6W (Part 2B; recommended therapeutic dose)",
    regions        = "Multinational phase 1 (GARNET, NCT02715284)",
    ada_positive_pct = 18.5,
    hepatic_impairment_pct = c(none = 89.0, mild = 10.1, moderate = 0.9, severe = 0),
    renal_impairment_pct   = c(none = 38.3, mild = 43.0, moderate = 18.3, severe = 0.4),
    notes          = "Baseline demographics from Melhem 2022 Table 2 (analysis set, N = 546). 4783 PK observations used in the final model after excluding 21 records with |CWRES| > 5. Hepatic and renal impairment did not affect dostarlimab PK. Reference patient (Methods, last paragraph): female, 70 kg, 64 years, ALB 39 g/L, ALT 18 U/L."
  )

  ini({
    # Structural PK parameters at the reference patient (female, 70 kg, 64 y,
    # ALB 39 g/L, ALT 18 U/L). Values are the typical-value population PK
    # parameter estimates from Melhem 2022 Table 3.
    lcl  <- log(0.179);  label("Baseline clearance at t=0 for reference patient (CL_base, L/day)") # Melhem 2022 Table 3: CL = 0.179 L/d
    lvc  <- log(2.98);   label("Central volume of distribution for reference patient (Vc_base, L)") # Melhem 2022 Table 3: Vc = 2.98 L
    lq   <- log(0.547);  label("Intercompartmental clearance (Q, L/day)")                           # Melhem 2022 Table 3: Q = 0.547 L/d
    lvp  <- log(2.10);   label("Peripheral volume of distribution for 70 kg reference (Vp, L)")     # Melhem 2022 Table 3: Vp = 2.10 L

    # Time-dependent clearance (sigmoid I_max function of time since first
    # dose; Melhem 2022 Results, equation block):
    #   CL_timebase = CL_base * exp(I_max * t^Hill / (T50^Hill + t^Hill))
    # I_max is negative (a fractional decrease in log-CL at t >> T50). To
    # keep individual I_max strictly negative under log-normal IIV, the
    # magnitude is stored as log|I_max| and the sign is applied in model()
    # (same pattern as Masters_2022_avelumab). The resulting maximum
    # reduction at t >> T50 is 1 - exp(I_max) = 1 - exp(-0.161) = 0.149 (14.9 %),
    # matching the paper.
    limax <- log(0.161); label("log|I_max|; magnitude of the log-CL reduction at t >> T50 (unitless)") # Melhem 2022 Table 3: I_max = -0.161
    lt50  <- log(108);   label("log T50; time at half of I_max (days)")                                # Melhem 2022 Table 3: T50 = 108 d
    lhill <- log(5.29);  label("log Hill; sigmoid steepness coefficient (unitless)")                   # Melhem 2022 Table 3: Hill = 5.29

    # Allometric exponents on body weight (reference 70 kg).
    e_wt_cl    <- 0.470;  label("Allometric exponent of WT on CL (unitless)")              # Melhem 2022 Table 3: Effect of WT on CL = 0.470
    e_wt_vc_vp <- 0.419;  label("Shared allometric exponent of WT on Vc and Vp (unitless)") # Melhem 2022 Table 3: Effect of WT on Vc and Vp = 0.419

    # Continuous covariate effects on CL (power form: (cov / ref)^theta).
    e_age_cl <- -0.227;   label("AGE exponent on CL (power form (AGE/64)^theta)")        # Melhem 2022 Table 3: Effect of age on CL = -0.227
    e_alb_cl <- -1.01;    label("ALB exponent on CL (power form (ALB/39)^theta)")        # Melhem 2022 Table 3: Effect of ALB on CL = -1.01
    e_alt_cl <- -0.0585;  label("ALT exponent on CL (power form (ALT/18)^theta)")        # Melhem 2022 Table 3: Effect of ALT on CL = -0.0585
    e_alb_vc <- -0.153;   label("ALB exponent on Vc (power form (ALB/39)^theta)")        # Melhem 2022 Table 3: Effect of ALB on Vc = -0.153

    # Categorical sex effect on CL and Vc (paper codes a male indicator with
    # female as the reference; encoded here as (1 + theta * (1 - SEXF))).
    e_sex_cl <- 0.165; label("Male fractional change on CL vs. female reference")        # Melhem 2022 Table 3: Effect of male on CL = 0.165
    e_sex_vc <- 0.162; label("Male fractional change on Vc vs. female reference")        # Melhem 2022 Table 3: Effect of male on Vc = 0.162

    # IIV — block on (CL, Vc) and independent eta on |I_max|.
    # Melhem 2022 Table 3 reports the variance/covariance matrix elements:
    #   omega^2_CL          = 0.0551
    #   cov(CL, Vc)         = 0.0210  (correlation 0.557, consistent with the
    #                                  parenthetical 0.557 in the table)
    #   omega^2_Vc          = 0.0258
    #   omega^2_Imax        = 0.537
    etalcl + etalvc ~ c(0.0551,
                        0.0210, 0.0258)
    etalimax ~ 0.537

    # Residual error (combined additive + proportional). Melhem 2022 Table 3
    # reports proportional 0.133 (unitless) and additive 2.79 mg/L. The
    # additive value is reported with units mg/L (not mg^2/L^2), indicating
    # the table column is the standard deviation, not the variance; the
    # proportional value is treated identically as the SD (CV ~ 13.3 %),
    # which matches the residual-error magnitudes typical for mAb PK assays.
    propSd <- 0.133; label("Proportional residual error SD (fraction)")               # Melhem 2022 Table 3: Proportional residual variability = 0.133
    addSd  <- 2.79;  label("Additive residual error SD (mg/L)")                       # Melhem 2022 Table 3: Additive residual variability = 2.79 mg/L
  })

  model({
    # Allometric weight scaling (reference 70 kg).
    wt_cl    <- (WT / 70)^e_wt_cl
    wt_vc_vp <- (WT / 70)^e_wt_vc_vp

    # Continuous covariate effects (power form, Melhem 2022 final-model
    # equations).
    age_cl <- (AGE / 64)^e_age_cl
    alb_cl <- (ALB / 39)^e_alb_cl
    alt_cl <- (ALT / 18)^e_alt_cl
    alb_vc <- (ALB / 39)^e_alb_vc

    # Sex effect: female reference (SEXF = 1 -> multiplier 1); male
    # (SEXF = 0) -> multiplier 1 + theta.
    sex_cl <- 1 + e_sex_cl * (1 - SEXF)
    sex_vc <- 1 + e_sex_vc * (1 - SEXF)

    # Time-dependent CL (Hill function of time since first dose; t in days).
    # I_max < 0; sign applied here to keep individual values strictly negative
    # under log-normal IIV on |I_max|.
    imax_i <- -exp(limax + etalimax)
    t50    <- exp(lt50)
    hill   <- exp(lhill)
    td_cl  <- exp(imax_i * t^hill / (t50^hill + t^hill))

    # Individual PK parameters. CL_base is at t=0 for the reference patient;
    # td_cl folds in the time dependency.
    cl_base <- exp(lcl + etalcl) * wt_cl * age_cl * alb_cl * alt_cl * sex_cl
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
