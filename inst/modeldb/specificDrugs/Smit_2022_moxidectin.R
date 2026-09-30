Smit_2022_moxidectin <- function() {
  description <- paste0(
    "Two-compartment population PK model for oral moxidectin in 96 ",
    "Strongyloides stercoralis-infected adults in Laos (single 2-12 mg doses; ",
    "capillary whole-blood volumetric microsamples). First-order absorption ",
    "from a depot compartment after an absorption lag time, linear ",
    "elimination, allometric body-weight scaling (exponents fixed at 0.75 on ",
    "CL/F and Q/F and 1 on V1/F and V2/F, reference 70 kg) and a power effect ",
    "of age on the central volume (reference 44.25 years). Log-normal IIV on ",
    "all six structural parameters (lag-time IIV fixed) and proportional ",
    "residual error."
  )
  reference <- paste(
    "Smit C, Hofmann D, Sayasone S, Keiser J, Pfister M (2022).",
    "Characterization of the Population Pharmacokinetics of Moxidectin in",
    "Adults Infected with Strongyloides Stercoralis: Support for a Fixed-Dose",
    "Treatment Regimen. Clin Pharmacokinet 61:123-132.",
    "doi:10.1007/s40262-021-01048-4",
    sep = " "
  )
  vignette <- "Smit_2022_moxidectin"
  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "ng/mL"
  )

  compartmentData <- list(
    depot = list(analyte = "moxidectin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "moxidectin", units = "mg", specimen = "whole blood", verified = TRUE),
    peripheral1 = list(analyte = "moxidectin", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Time-fixed (single-dose study). Allometric scaling on all four ",
        "disposition parameters with exponents FIXED by the authors: 0.75 on ",
        "CL/F and Q/F, 1 on V1/F and V2/F (Smit 2022 Table 2 'WT on ... FIX' ",
        "rows and footnote a; supplementary MLXTRAN code beta_*_logtWT ",
        "method=FIXED). Centred on 70 kg (MLXTRAN logtWT = log(WT/70)). Cohort ",
        "median 56.2 kg, range 36.2-82.6 kg (Table 1)."
      ),
      source_name = "WT"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Power effect on V1/F only, centred on 44.2517 years (supplementary ",
        "MLXTRAN code logtAGE = log(AGE/44.2517); Table 2 footnote rounds the ",
        "reference to 44.3 years). Cohort median 45.0 years, range 22-65 ",
        "(Table 1). Adults only; the Discussion cautions against ",
        "extrapolating to children."
      ),
      source_name = "AGE"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Biological sex indicator, 1 = female, 0 = male.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = paste0(
        "Screened (source column SEX_F in the MLXTRAN data block). Sex on ",
        "V2/F was statistically significant (dOFV -9.8) but was NOT retained ",
        "because it gave no clear improvement in goodness of fit (Smit 2022 ",
        "Results 3.2). 37 of 96 (39%) female (Table 1)."
      ),
      source_name = "SEX_F"
    ),
    FFM = list(
      description = "Lean body weight (Janmahasatian equation), i.e. fat-free mass.",
      units = "kg",
      type = "continuous",
      notes = paste0(
        "Screened but not retained (Smit 2022 Methods 2.3: 'lean body weight ",
        "(LBW, based on the Janmahasatian equation)'). The Janmahasatian LBW ",
        "is the register's FFM quantity. Cohort values not reported."
      ),
      source_name = "LBW"
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Screened but not retained. Median 22.2 kg/m^2, range 17.0-32.3 (Table 1).",
      source_name = "BMI"
    ),
    HT = list(
      description = "Body height",
      units = "cm",
      type = "continuous",
      notes = "Recorded and screened but not retained. Median 158 cm, range 137-172 (Table 1).",
      source_name = "HT"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 96L,
    n_studies = 1L,
    age_range = "22-65 years",
    age_median = "45.0 years",
    weight_range = "36.2-82.6 kg",
    weight_median = "56.2 kg",
    sex_female_pct = 39,
    race_ethnicity = "not reported; participants resident in NamBak District, northern Laos",
    disease_state = paste0(
      "Strongyloides stercoralis-infected adults (baseline infection ",
      "intensity light 16.7%, moderate 39.6%, heavy 43.6%)"
    ),
    dose_range = paste0(
      "single oral dose of 2, 4, 6, 8, 10 or 12 mg moxidectin (2 mg tablets) ",
      "after a local lunch (n = 15, 16, 14, 21, 15, 15 per arm)"
    ),
    regions = "Laos",
    n_observations = paste0(
      "762 capillary whole-blood samples (Mitra 30 uL volumetric ",
      "microsamples) at 2, 4, 6, 7 h and 1, 3, 7, 28 days post-dose; 158 ",
      "(20.7%) below the 1.5 ng/mL LLOQ, handled with the Monolix SAEM ",
      "censoring extension (M3-like)"
    ),
    notes = paste0(
      "PK sub-study embedded in a phase IIa randomized, placebo-controlled, ",
      "dose-escalation trial (NCT04056325), November 2019 - March 2020. ",
      "Fitted with SAEM in Monolix 2019R2. Baseline characteristics are Smit ",
      "2022 Table 1; final estimates Table 2; MLXTRAN code in the Electronic ",
      "Supplementary Material."
    )
  )

  ini({
    # Structural parameters: Smit 2022 Table 2 'Population parameters'
    # (typical values for a 70 kg adult aged 44.3 years).
    ltlag <- log(1.64)
    label("Absorption lag time Tlag (h)")
    # Table 2: Tlag = 1.64 h (RSE 1.6%); bootstrap 1.63 (1.56-1.70).

    lka <- log(3.38)
    label("First-order absorption rate constant Ka (1/h)")
    # Table 2: Ka = 3.38 /h (RSE 11%); bootstrap 3.35 (2.53-4.44).

    lcl <- log(4.47)
    label("Apparent clearance CL/F (L/h) at 70 kg")
    # Table 2: CL/F pop = 4.47 L/h (RSE 7.6%); bootstrap 4.47 (3.63-5.39).

    lvc <- log(136)
    label("Apparent central volume V1/F (L) at 70 kg and 44.25 years")
    # Table 2: V1/F pop = 136 L (RSE 3.4%); bootstrap 136 (126-146).

    lq <- log(10.0)
    label("Apparent intercompartmental clearance Q/F (L/h) at 70 kg")
    # Table 2: Q/F pop = 10.0 L/h (RSE 4.5%); bootstrap 10.0 (9.10-11.0).

    lvp <- log(1172)
    label("Apparent peripheral volume V2/F (L) at 70 kg")
    # Table 2: V2/F pop = 1172 L (RSE 11%); bootstrap 1169 (875-1523).

    # Allometric exponents, all fixed (Table 2 'WT on ...' rows = FIX;
    # MLXTRAN beta_*_logtWT method=FIXED).
    e_wt_cl <- fixed(0.75)
    label("Allometric exponent of body weight on CL/F (unitless)")
    # Table 2: WT on CL/F = 0.75 FIX.

    e_wt_vc <- fixed(1)
    label("Allometric exponent of body weight on V1/F (unitless)")
    # Table 2: WT on V1/F = 1 FIX.

    e_wt_q <- fixed(0.75)
    label("Allometric exponent of body weight on Q/F (unitless)")
    # Table 2: WT on Q/F = 0.75 FIX.

    e_wt_vp <- fixed(1)
    label("Allometric exponent of body weight on V2/F (unitless)")
    # Table 2: WT on V2/F = 1 FIX.

    e_age_vc <- -0.422
    label("Power exponent of age on V1/F (unitless)")
    # Results 3.2: 'estimated exponent of -0.422 (95% CI -0.674 to -0.197)'.
    # Table 2 prints the same estimate rounded to -0.42 (RSE 30%). The
    # three-decimal value reproduces the Results' V1 of 199 L at 18 years
    # and 116 L at 65 years for a 70 kg adult.

    # IIV: Monolix log-normal random effects. Table 2 legend: 'omega
    # standard deviation of the interindividual variability parameter',
    # so the variances below are omega^2.
    etaltlag ~ fixed(0.01)
    # Table 2: omega Tlag = 0.1 FIX (SD) -> variance 0.01. MLXTRAN
    # omega_Tlag method=FIXED. Results 3.2 explains the fixing.

    etalka ~ 0.451584
    # Table 2: omega Ka = 0.672 (SD, RSE 14%) -> 0.672^2.

    etalcl ~ 0.390625
    # Table 2: omega CL/F = 0.625 (SD, RSE 8.6%) -> 0.625^2.

    etalvc ~ 0.088209
    # Table 2: omega V1/F = 0.297 (SD, RSE 9.2%) -> 0.297^2.

    etalq ~ 0.117649
    # Table 2: omega Q/F = 0.343 (SD, RSE 11%) -> 0.343^2.

    etalvp ~ 0.606841
    # Table 2: omega V2/F = 0.779 (SD, RSE 11%) -> 0.779^2.

    # Residual error: proportional only (MLXTRAN errorModel=proportional(b)).
    propSd <- 0.165
    label("Proportional residual error (fraction)")
    # Table 2: proportional error = 0.165 (RSE 4.9%), footnote b 'shown as
    # the standard deviation'.
  })

  model({
    # Covariate transforms (MLXTRAN [COVARIATE] block):
    # logtWT = log(WT/70); logtAGE = log(AGE/44.2517).
    wtref <- 70
    ageref <- 44.2517

    # Individual parameters (MLXTRAN [INDIVIDUAL] block; all logNormal).
    tlag <- exp(ltlag + etaltlag)
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) * (WT / wtref)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / wtref)^e_wt_vc * (AGE / ageref)^e_age_vc
    q <- exp(lq + etalq) * (WT / wtref)^e_wt_q
    vp <- exp(lvp + etalvp) * (WT / wtref)^e_wt_vp

    # Micro-constants (MLXTRAN model text: k10 = CL/V1, k12 = Q/V1, k21 = Q/V2).
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ODEs (MLXTRAN ddt_Ad, ddt_Ac, ddt_Ap).
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Lag time. The MLXTRAN code sets KA = 0 while t < Tlag, which for the
    # single dose of the study is identical to delaying the dose by Tlag.
    alag(depot) <- tlag

    # Dose in mg, volume in L: mg/L * 1000 = ng/mL.
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
