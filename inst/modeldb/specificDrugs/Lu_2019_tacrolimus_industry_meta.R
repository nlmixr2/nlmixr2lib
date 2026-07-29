Lu_2019_tacrolimus_industry_meta <- function() {
  description <- paste0(
    "Industry meta-analysis. Two-compartment population PK model for ",
    "oral tacrolimus immediate-release (IR-T; Prograf, twice daily) and ",
    "prolonged-release (PR-T; Advagraf / Astagraf XL, once daily) ",
    "formulations in adult and paediatric liver, kidney, and heart ",
    "transplant recipients (Lu 2019). Pooled individual-patient data ",
    "from 8 Astellas Phase II studies (n = 408 patients, 23,176 whole-",
    "blood concentration records). Structural model: first-order ",
    "absorption with formulation-dependent Ka (PR-T ~50% slower than ",
    "IR-T), fixed absorption lag time, and two-compartment disposition ",
    "with first-order elimination. Significant covariates: Asian race ",
    "on CL/F (+59% vs Whites); log-AST on CL/F, Vc/F, Vp/F, and F1 ",
    "(power normalised at LAST = 3.15, i.e., AST ~= 23.3 IU/L); female ",
    "sex on Vc/F (-44.6% vs males); albumin on Vc/F and F1; and ",
    "Asian / Black race on F1 (Asians > Whites > Blacks). Type of ",
    "organ transplanted and adult-vs-paediatric population had no ",
    "significant effect on PK parameters."
  )
  reference <- paste0(
    "Lu Z, Bonate P, Keirns J. Population pharmacokinetics of ",
    "immediate- and prolonged-release tacrolimus formulations in ",
    "liver, kidney and heart transplant recipients. Br J Clin ",
    "Pharmacol. 2019;85(8):1692-1703. doi:10.1111/bcp.13952."
  )
  vignette <- "Lu_2019_tacrolimus_industry_meta"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    RACE_ASIAN = list(
      description        = "Asian race indicator (self-reported)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Asian; White or Black)",
      notes              = paste0(
        "Time-fixed per subject. Lu 2019 Methods 2.3.3: 'Race was ",
        "coded as a 3-category covariate in the analysis: White, ",
        "Black or Asian.' The Asian cohort in Lu 2019 (n = 44 of ",
        "408; Table 2) is heterogeneous - Japanese, Chinese and ",
        "other Far East groups pooled into a single indicator. ",
        "Effects: CL/F linear (1 + 0.59 * RACE_ASIAN), F1 linear ",
        "(1 + 0.25 * RACE_ASIAN + e_race_black_f * RACE_BLACK)."
      ),
      source_name        = "RACE (Asian level)"
    ),
    RACE_BLACK = list(
      description        = "Black / African American race indicator (self-reported)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Black; White or Asian)",
      notes              = paste0(
        "Time-fixed per subject. n = 24 of 408 patients (Table 2). ",
        "Effect retained only on F1 (-43.3% vs Whites); Black race ",
        "on CL/F was dropped during backward elimination due to ",
        "lack of precision (Lu 2019 Results 3.2)."
      ),
      source_name        = "RACE (Black level)"
    ),
    SEXF = list(
      description        = "Female sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = paste0(
        "Time-fixed per subject. 132 females / 408 (32.4%) per ",
        "Table 2. Linear effect on Vc/F: (1 - 0.446 * SEXF), i.e., ",
        "44.6% lower Vc/F in females (Table 3 'Sex on Vc/F' = ",
        "-0.446). The abstract narrative '55% lower' is a rounded ",
        "characterisation; the structural model parameter is ",
        "-0.446."
      ),
      source_name        = "SEX"
    ),
    ALB = list(
      description        = "Serum albumin concentration (time-varying)",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "SI units. Power scaling on Vc/F and F1 with reference ",
        "ALB = 39 g/L (Lu 2019 final-model equations; the same ",
        "value is used as the 'normal' ALB scenario in the Lu 2019 ",
        "simulation, Methods 2.3.5). Exponents: 1.03 on Vc/F, 1.04 ",
        "on F1."
      ),
      source_name        = "ALB"
    ),
    AST = list(
      description        = "Serum aspartate aminotransferase activity (time-varying)",
      units              = "IU/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Power scaling on CL/F, Vc/F, Vp/F, and F1. The Lu 2019 ",
        "final-model equations write the AST factor as ",
        "(LAST / 3.15)^theta with LAST = log(AST), 'normalized at ",
        "3.15 IU/L' (Lu 2019 Results 3.2). The exp(theta * (LAST - ",
        "3.15)) form gives the published behaviour - a 2.7-fold ",
        "(~e-fold) increase in AST shifts each parameter by ",
        "exp(theta), matching the paper's stated factors ",
        "(exp(-0.318) ~= 0.73 on CL/F, exp(1.73) ~= 5.6 on Vc/F, ",
        "exp(-0.945) ~= 0.39 on Vp/F, exp(0.74) ~= 2.1 on F1). ",
        "The reference AST is therefore AST_ref = exp(3.15) ~= ",
        "23.3 IU/L. AST simulation scenarios in Lu 2019 Methods ",
        "2.3.5: normal 25 IU/L, mild elevation 100 IU/L, moderate ",
        "elevation 400 IU/L."
      ),
      source_name        = "AST (entered as LAST = log(AST))"
    ),
    FORM_TAC_IR = list(
      description        = "Tacrolimus immediate-release vs prolonged-release formulation indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (prolonged-release, Advagraf / Astagraf XL)",
      notes              = paste0(
        "Per-subject (or per-occasion in cross-over studies). ",
        "Canonical convention: FORM_TAC_IR = 1 for the twice-daily ",
        "immediate-release tacrolimus formulation (Prograf), ",
        "FORM_TAC_IR = 0 for the once-daily prolonged-release ",
        "formulation (Advagraf / Astagraf XL). Lu 2019 uses the ",
        "opposite polarity internally (FORMULATION = 1 for PR-T); ",
        "model() derives form_pr = 1 - FORM_TAC_IR to match the ",
        "paper's published equation. Effect retained on Ka only: ",
        "Ka(PR-T) = 0.375 * 0.499 = 0.187 1/h (50% slower than ",
        "IR-T per Lu 2019 Results 3.2 and Table 3 'Prolonged-",
        "release tacrolimus on Ka' = 0.499). A formulation term on ",
        "F1 appears in the published F1 equation but is not ",
        "tabulated as a final-model fixed effect (Table 3 has no ",
        "row for Prolonged-release tacrolimus on F1); the abstract ",
        "characterises the PR-T:IR-T relative-bioavailability ",
        "geometric-mean ratio (95%, 90% CI 89-101%) as a posthoc ",
        "empirical-Bayes summary, not a structural parameter, so ",
        "the formulation effect on F1 is omitted here. See ",
        "vignette Assumptions and deviations."
      ),
      source_name        = "FORMULATION (1 = PR-T in Lu 2019)"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 408L,
    n_studies        = 8L,
    age_range        = "5 - 71 years (17 paediatric, 391 adult)",
    age_median       = "48 years",
    weight_range     = "18.5 - 148.5 kg",
    weight_median    = "74 kg",
    sex_female_pct   = 32.4,
    race_ethnicity   = c(White = 83.3, Asian = 10.8, Black = 5.9),
    disease_state    = paste0(
      "Adult and paediatric liver, kidney, and heart solid-organ ",
      "transplant recipients on immediate-release (Prograf) and/or ",
      "prolonged-release (Advagraf / Astagraf XL) oral tacrolimus ",
      "immunosuppression. Six conversion studies in stable ",
      "recipients converting from IR-T to PR-T (n = 265) and two ",
      "comparative studies in primary kidney or liver transplant ",
      "recipients (n = 143)."
    ),
    dose_range       = paste0(
      "Oral tacrolimus, dose individualised by trough monitoring ",
      "(target 5 - 20 ng/mL whole-blood trough); both formulations ",
      "given as 1 mg : 1 mg total-daily-dose conversion. Lu 2019 ",
      "simulation scenario used IR-T 5 mg BID and PR-T 10 mg QD."
    ),
    n_concentrations = 23176L,
    sampling_design  = paste0(
      "Rich PK sampling on PK-assessment days (predose and 0.5, 1, ",
      "2, 3, 4, 6, 8, 12, 12.5, 13, 14, 15, 16, 18, 20, 24 h post-",
      "dose), plus trough monitoring throughout the treatment ",
      "period. Whole-blood concentrations measured by LC-MS/MS in ",
      "seven studies and by immunoassay in one study (FJ-506E-",
      "KT01, Japanese cohort); lower limit of quantitation 0.1 ",
      "ng/mL."
    ),
    studies          = paste0(
      "8 Astellas Phase II studies pooled (Lu 2019 Table 1): ",
      "02-0-131 (adult kidney, n = 57); FG-506E-12-02 (adult ",
      "kidney, n = 60); FJ-506E-KT01 (Japanese adult kidney, ",
      "n = 35); 02-0-152 (adult liver, n = 51); FG-506-15-02 ",
      "(adult heart, n = 45); 03-0-160 (paediatric liver, n = 17, ",
      "mean age 9 years); FG-506E-12-01 (primary kidney, n = 66); ",
      "FG-506-11-01 (primary liver, n = 77)."
    ),
    regions          = "Multi-national (8 Astellas Phase II studies)",
    notes            = paste0(
      "Software: NONMEM v7.3 with FOCE-I; bootstrap n = 1000 ",
      "(509 successful runs, 491 minimization-terminated runs ",
      "skipped). Type of organ transplanted (kidney vs liver vs ",
      "heart) and adult-vs-paediatric population had no ",
      "significant effect on principal PK parameters. CYP3A4 / ",
      "CYP3A5 / P-gp genotype not included as covariates. Renal-",
      "function covariates not evaluated because urinary excretion ",
      "accounts for <2% of dose. Haematocrit and time-post-",
      "transplant noted as potentially relevant but not available ",
      "across all studies."
    )
  )

  ini({
    # --- Structural PK (Lu 2019 Table 3, final-model 'Value' column) ---
    # Reference subject for typical values: White, male, IR-T (Prograf),
    # AST = exp(3.15) ~= 23.3 IU/L (LAST = 3.15 normalisation point),
    # ALB = 39 g/L. CL/F, Q/F in L/h; Vc/F, Vp/F in L; Ka in 1/h;
    # ALAG1 in hours. Tacrolimus whole-blood concentrations in ng/mL.
    lcl     <- log(44.3)  ; label("Apparent oral clearance CL/F at reference (L/h)")                  # Lu 2019 Table 3 CL/F = 44.3 L/h (RSE 3.43%)
    lvc     <- log(110)   ; label("Apparent central volume Vc/F at reference (L)")                    # Lu 2019 Table 3 Vc/F = 110 L (RSE 10.55%)
    lq      <- log(131)   ; label("Apparent inter-compartmental clearance Q/F (L/h)")                 # Lu 2019 Table 3 Q/F = 131 L/h (RSE 5.42%)
    lvp     <- log(3180)  ; label("Apparent peripheral volume Vp/F at reference (L)")                 # Lu 2019 Table 3 Vp/F = 3180 L (RSE 7.39%)
    lka     <- log(0.375) ; label("Absorption rate constant Ka for IR-T (1/h)")                       # Lu 2019 Table 3 Ka = 0.375 1/h (RSE 4.48%) for immediate-release
    lfdepot <- log(1.51)  ; label("Relative bioavailability F1 at reference (unitless scaling)")     # Lu 2019 Table 3 F1 = 1.51 (RSE 2.96%); typical value at White / IR-T / ALB = 39 / AST = exp(3.15)
    ltlag   <- log(0.44)  ; label("Absorption lag time ALAG1 (h)")                                    # Lu 2019 Table 3 ALAG1 = 0.44 h (RSE 1.17%)

    # --- Covariate effects (Lu 2019 Table 3) ---
    # CL/F covariate effects.
    e_race_asian_cl <- 0.59    ; label("Linear effect of Asian race on CL/F: CL_typical * (1 + 0.59 * RACE_ASIAN)")            # Lu 2019 Table 3 'Asian race on CL/F' = 0.59 (+59% vs Whites). Black race on CL/F was dropped during backward elimination.
    e_ast_cl        <- -0.318  ; label("Power exponent of log-AST on CL/F: CL_typical * exp(-0.318 * (log(AST) - 3.15))")     # Lu 2019 Table 3 'AST on CL/F' = -0.318. Centring point LAST = 3.15 -> AST_ref = exp(3.15) ~= 23.3 IU/L.

    # Vc/F covariate effects.
    e_sex_vc <- -0.446 ; label("Linear effect of female sex on Vc/F: Vc_typical * (1 - 0.446 * SEXF)")                          # Lu 2019 Table 3 'Sex on Vc/F' = -0.446 (-44.6% vs males). Abstract '55% lower' is a rounded narrative; the structural parameter is -0.446.
    e_ast_vc <- 1.73   ; label("Power exponent of log-AST on Vc/F: Vc_typical * exp(1.73 * (log(AST) - 3.15))")                 # Lu 2019 Table 3 'AST on Vc/F' = 1.73.
    e_alb_vc <- 1.03   ; label("Power exponent of ALB on Vc/F: Vc_typical * (ALB / 39)^1.03")                                   # Lu 2019 Table 3 'ALB on Vc/F' = 1.03 with ALB_ref = 39 g/L.

    # Vp/F covariate effects.
    e_ast_vp <- -0.945 ; label("Power exponent of log-AST on Vp/F: Vp_typical * exp(-0.945 * (log(AST) - 3.15))")               # Lu 2019 Table 3 'AST on Vp/F' = -0.945.

    # Ka covariate effects.
    e_form_ka_pr <- 0.499 ; label("Multiplicative effect of prolonged-release formulation on Ka: Ka(PR-T) / Ka(IR-T) = 0.499")  # Lu 2019 Table 3 'Prolonged-release tacrolimus on Ka' = 0.499 (PR-T ~50% slower than IR-T). Applied to form_pr = 1 - FORM_TAC_IR.

    # F1 covariate effects.
    e_race_asian_f <- 0.25   ; label("Linear effect of Asian race on F1: F_typical * (1 + 0.25 * RACE_ASIAN + e_race_black_f * RACE_BLACK)") # Lu 2019 Table 3 'Asian race on F1' = 0.25.
    e_race_black_f <- -0.433 ; label("Linear effect of Black race on F1 (combined with Asian effect in race_f)")                              # Lu 2019 Table 3 'Black race on F1' = -0.433 (-43.3% vs Whites).
    e_ast_f        <- 0.74   ; label("Power exponent of log-AST on F1: F_typical * exp(0.74 * (log(AST) - 3.15))")                            # Lu 2019 Table 3 'AST on F1' = 0.74.
    e_alb_f        <- 1.04   ; label("Power exponent of ALB on F1: F_typical * (ALB / 39)^1.04")                                              # Lu 2019 Table 3 'ALB on F1' = 1.04 with ALB_ref = 39 g/L.

    # --- Inter-individual variability (Lu 2019 Table 3 IPV columns) ---
    # IPV reported as CV%; map to log-scale variance via omega^2 = log(CV^2 + 1)
    # (lognormal IIV on the exponentiated structural parameter, eta ~ N(0, omega^2)).
    # IIV on ALAG1 was not estimated in the final model (Lu 2019 Results 3.1).
    etalcl     ~ 0.0912   # Lu 2019 Table 3 IPV-CL = 30.9% CV -> log(1 + 0.309^2) = 0.0912
    etalvc     ~ 0.7531   # Lu 2019 Table 3 IPV-Vc = 106% CV  -> log(1 + 1.06^2)  = 0.7531
    etalq      ~ 0.1436   # Lu 2019 Table 3 IPV-Q  = 39.3% CV -> log(1 + 0.393^2) = 0.1436
    etalvp     ~ 0.6831   # Lu 2019 Table 3 IPV-Vp = 99% CV   -> log(1 + 0.99^2)  = 0.6831
    etalka     ~ 0.1187   # Lu 2019 Table 3 IPV-Ka = 35.5% CV -> log(1 + 0.355^2) = 0.1187
    etalfdepot ~ 0.0890   # Lu 2019 Table 3 IPV-F1 = 30.5% CV -> log(1 + 0.305^2) = 0.0890

    # --- Residual unexplained variability (Lu 2019 Table 3 RV rows) ---
    # The paper used log-transformed-both-sides (LTBS): NONMEM error model
    # Y = LOG(F) + EPS(1), so the residual coefficient maps directly to the
    # log-scale residual standard deviation used by Cc ~ lnorm(expSd).
    # Two assay-specific magnitudes: RV1 = 21.1% (LC-MS/MS, seven studies)
    # and RV2 = 15.8% (immunoassay, study FJ-506E-KT01 only). The package
    # ships the LC-MS/MS value as the default; the vignette documents the
    # immunoassay alternative for users simulating the Japanese cohort.
    expSd <- 0.211 ; label("Log-scale residual SD (LC-MS/MS assay)")  # Lu 2019 Table 3 RV1 = 21.1% (LC-MS/MS). RV2 = 15.8% for immunoassay (FJ-506E-KT01) is noted in the vignette.
  })

  model({
    # --- Reference values for power-model normalisation ---
    # AST_ref = exp(3.15) ~= 23.336 IU/L; the paper writes the AST factor
    # as (LAST/3.15)^theta with LAST = log(AST), 'normalized at 3.15 IU/L'.
    # The exp(theta * (LAST - 3.15)) form is mathematically equivalent and
    # reproduces the paper's stated factors (e.g. exp(-0.318) ~= 0.73 on
    # CL/F for an e-fold AST increase).
    # ALB_ref = 39 g/L (Lu 2019 final-model equations; also the 'normal'
    # ALB simulation scenario in Methods 2.3.5).

    # --- Polarity adjustment ---
    # Lu 2019 uses FORMULATION = 1 for prolonged-release (PR-T); the
    # canonical FORM_TAC_IR uses 1 for immediate-release (Prograf). Derive
    # Lu's PR-T indicator from the canonical column so the published Ka
    # equation reads naturally.
    form_pr <- 1 - FORM_TAC_IR

    # --- Derived covariate factors ---
    ast_cl <- exp(e_ast_cl * (log(AST) - 3.15))
    ast_vc <- exp(e_ast_vc * (log(AST) - 3.15))
    ast_vp <- exp(e_ast_vp * (log(AST) - 3.15))
    ast_f  <- exp(e_ast_f  * (log(AST) - 3.15))
    alb_vc <- (ALB / 39) ^ e_alb_vc
    alb_f  <- (ALB / 39) ^ e_alb_f

    race_cl <- 1 + e_race_asian_cl * RACE_ASIAN
    race_f  <- 1 + e_race_asian_f * RACE_ASIAN + e_race_black_f * RACE_BLACK
    sex_vc  <- 1 + e_sex_vc * SEXF
    ka_form <- 1 - (1 - e_form_ka_pr) * form_pr

    # --- Individual PK parameters ---
    cl   <- exp(lcl + etalcl) * race_cl * ast_cl
    vc   <- exp(lvc + etalvc) * sex_vc  * alb_vc * ast_vc
    q    <- exp(lq  + etalq)
    vp   <- exp(lvp + etalvp) * ast_vp
    ka   <- exp(lka + etalka) * ka_form
    tlag <- exp(ltlag)

    # --- Micro-constants ---
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # --- Two-compartment oral PK with first-order absorption + lag ---
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # --- Bioavailability (relative scaling factor) and lag time ---
    f(depot)    <- exp(lfdepot + etalfdepot) * race_f * alb_f * ast_f
    alag(depot) <- tlag

    # --- Observation: tacrolimus whole-blood concentration (ng/mL) ---
    # Dose in mg, central amount in mg, vc in L -> mg/L; multiply by 1000
    # to give ng/mL (the units used throughout Lu 2019).
    Cc <- central / vc * 1000
    Cc ~ lnorm(expSd)
  })
}
