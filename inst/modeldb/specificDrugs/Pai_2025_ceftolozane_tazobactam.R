Pai_2025_ceftolozane_tazobactam <- function() {
  description <- "Joint one-compartment population PK model for ceftolozane and tazobactam given together as ceftolozane/tazobactam continuous infusion in hospitalized adults with Gram-negative infections (retrospective therapeutic drug monitoring cohort, Bologna, Italy). The combined ceftolozane + tazobactam dose enters a single dosing compartment that empties at 1/h and splits a fixed fraction 1/3 to tazobactam and 2/3 to ceftolozane (the 2:1 product ratio); each drug then has its own volume (fixed from the registration model) and clearance. Both clearances are power functions of the 2021 race-free CKD-EPI eGFR de-indexed to the individual Du Bois body surface area (mL/min), with strongly correlated random effects. Ceftolozane is the unsuffixed parent; tazobactam carries the _taz suffix."
  reference <- "Pai MP, Cojutti PG, Gatti M, Rinaldi M, Tonetti T, Siniscalchi A, Viale P, Pea F. Keeping the horse with the cart: are we underdosing tazobactam even when using continuous-infusion ceftolozane/tazobactam for effectively preventing resistance development by ESBL-producing Enterobacterales? Antimicrob Agents Chemother. 2025;69(12):e01215-25. doi:10.1128/aac.01215-25"
  vignette <- "Pai_2025_ceftolozane_tazobactam"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Doses are COMBINED ceftolozane + tazobactam amounts in mg, given into
  # `depot` (Supplement Figure S1: the infusion R0 enters A1, which is split
  # F to tazobactam and 1 - F to ceftolozane). A '3 g/1.5 g daily' regimen is
  # therefore 4500 mg per 24 h into `depot`; the 2 g/1 g loading dose is
  # 3000 mg over 1 h (Table 3 footnote: 'a 3 g loading dose (1 h infusion)').
  dosing <- "depot"

  compartmentData <- list(
    depot = list(analyte = "ceftolozane + tazobactam", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "ceftolozane", units = "mg", specimen = "plasma", verified = TRUE),
    central_taz = list(analyte = "tazobactam", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = "2021 race-free CKD-EPI creatinine eGFR de-indexed to the individual Du Bois body surface area (eGFR x BSA_DuBois / 1.73); raw mL/min, NOT BSA-normalized",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Source column ckdepi2021mlmin_bsa_dd (Supplement Table S1: 2021 CKD-EPI eGFR 'normalized",
        "to individual body surface area - Dubois Dubois equation', unit mL/min). The supplement",
        "STATA code computes it as CKDEPI2021 * bsa_dd / 1.73 with bsa_dd = ht^0.725 * 0.007184 *",
        "wt^0.425 (ht in cm, wt in kg) and serum creatinine in mg/dL. Selected as model 2 of",
        "Supplement Table S5 (AIC 4262.93; models 3-5, 8, 10 and 11 within 2 AIC points) and",
        "named in Table S6 as 'logtCkdepi2021mlmin_bsa_dd'. Reference value 60 mL/min (Table S5",
        "footnote a '(eGFR/60)'). Table 1 reports eGFR in mL/min median (IQR) 85.7 (42.4-107.5);",
        "Supplement Table S7 gives the five eGFR groups used for simulation. Patients on renal",
        "replacement therapy were excluded."
      ),
      source_name = "ckdepi2021mlmin_bsa_dd"
    )
  )

  # Supplement Tables S3-S4 univariate screens and Table S5 candidate models:
  # covariates tested on CLTOL and CLTAZ but not retained in the final model.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened (Supplement Table S3 CLTOL P = 0.862; Table S4 CLTAZ P = 0.603) but not retained. Cohort median (IQR) 66 (55-75) years (Table 1)."
    ),
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened together with ideal, adjusted, lean and dosing body weight, height, BMI and BSA (Supplement Tables S3-S4) but not retained; the Cockcroft-Gault and EKFC kidney-function alternatives in Table S5 were not selected. Cohort median (IQR) 70 (62-85) kg (Table 1)."
    ),
    SEXF = list(
      description = "Biological sex indicator, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      notes = "Screened (Supplement Table S3 CLTOL P = 0.313; Table S4 CLTAZ P = 0.509) but not retained. 61/139 female (Table 1)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 139L,
    n_studies = 1L,
    n_observations = "278 plasma samples, each assayed for both ceftolozane and tazobactam (Results)",
    age_range = "Adults; median (IQR) 66 (55-75) years, 26.6% aged 75 years or older (Table 1, Results)",
    weight_range = "Median (IQR) 70 (62-85) kg (Table 1)",
    bmi_range = "Median (IQR) 25.7 (22.5-29.0) kg/m^2; 18.7% with BMI >= 30 kg/m^2 (Table 1, Results)",
    sex_female_pct = 43.9,
    race_ethnicity = "Not reported (single Italian tertiary university hospital)",
    disease_state = "Hospitalized adults receiving continuous-infusion ceftolozane/tazobactam for Gram-negative infections with real-time TDM of both agents; patients on renal replacement therapy excluded (Methods).",
    renal_function = paste(
      "eGFR (2021 CKD-EPI) median (IQR) 84.2 (39.4-105.8) mL/min/1.73 m^2 and 85.7 (42.4-107.5)",
      "mL/min; Cockcroft-Gault eCLcr 76.1 (42.4-117.6) mL/min (Table 1). eGFR groups <15 / 15-29 /",
      "30-50 / 51-120 / >120 mL/min: n = 7 / 18 / 17 / 79 / 18 (Supplement Table S7)."
    ),
    dose_range = paste(
      "2 g/1 g ceftolozane/tazobactam loading dose over 1 h, then continuous infusion by eGFR",
      "(6 g/3 g daily if >= 50, 3 g/1.5 g if 30-49, 1.5 g/0.75 g if 15-29, 0.9 g/0.45 g if < 15",
      "mL/min/1.73 m^2), subsequently TDM-adjusted. Median (range) daily dose 3 g/1.5 g",
      "(1.5 g/0.75 g - 6 g/3 g) (Methods, Table 1)."
    ),
    regions = "Italy (IRCCS Azienda Ospedaliero-Universitaria di Bologna)",
    notes = paste(
      "Retrospective TDM cohort, November 2023 to June 2025. Steady-state total concentrations",
      "measured 48-72 h after start; free concentrations taken as 0.81 x total (ceftolozane) and",
      "0.70 x total (tazobactam). Median (IQR) fCss 33.3 (20.8-54.0) mg/L ceftolozane and 5.46",
      "(2.96-9.85) mg/L tazobactam (Table 1). LC-MS/MS assay, LLOQ 0.2 mg/L ceftolozane and",
      "0.1 mg/L tazobactam. Fitted with SAEM in Monolix 2024R1."
    )
  )

  ini({
    # Dosing-compartment outflow. Supplement Figure S1 writes
    # ddt_A1 = -F*A1 - (1-F)*A1, i.e. A1 empties with an implicit first-order
    # rate constant of 1/h, a fraction F to tazobactam and 1 - F to ceftolozane.
    lka <- fixed(log(1)); label("Dosing-compartment outflow rate constant (1/h)") # Supplement Figure S1 ODE 'ddt_A1 = -F*A1 - (1-F)*A1' (unit rate constant, not estimated)

    # Fraction of the combined dose that is tazobactam. Supplement Figure S1:
    # 'F represented in this model was set based on the known formulation
    # fraction of ceftolozane : tazobactam of 2:1' (initial-estimates table
    # F = 0.333); Table S6 prints the same value rounded as F_pop = 0.33.
    frac <- fixed(0.333); label("Fraction of the combined dose routed to tazobactam (fraction)") # Supplement Figure S1 table F = 0.333; Table S6 F_pop = 0.33, no SE (fixed)

    # Volumes, fixed from the Merck registration popPK model (Methods;
    # Supplement Figure S1 initial-estimates table).
    lvc <- fixed(log(32.1)); label("Ceftolozane volume V1 (L)") # Table S6 V1_pop = 32.1, no SE (fixed per Methods)
    lvc_taz <- fixed(log(43.3)); label("Tazobactam volume V2 (L)") # Table S6 V2_pop = 43.3, no SE (fixed per Methods)

    # Typical clearances at the reference eGFR of 60 mL/min.
    lcl <- log(2.78); label("Ceftolozane clearance CLTOL at CRCL = 60 mL/min (L/h)") # Table S6 Cltol_pop = 2.78 (RSE 5.74%)
    lcl_taz <- log(6.73); label("Tazobactam clearance CLTAZ at CRCL = 60 mL/min (L/h)") # Table S6 Cltaz_pop = 6.73 (RSE 6.68%)

    # Power exponents on (CRCL / 60); Monolix 'logt' covariate transform.
    e_crcl_cl <- 0.74; label("Power exponent of (CRCL/60) on ceftolozane CL (unitless)") # Table S6 beta_Cltol_logtCkdepi2021mlmin_bsa_dd = 0.74 (RSE 10.9%)
    e_crcl_cl_taz <- 0.84; label("Power exponent of (CRCL/60) on tazobactam CL (unitless)") # Table S6 beta_Cltaz_logtCkdepi2021mlmin_bsa_dd = 0.84 (RSE 11.1%)

    # IIV. Monolix reports omega as the SD of the log-normal random effect
    # with the CV in the next column. Variances use omega^2 = log(1 + CV^2)
    # from the four-significant-figure CVs (0.3686, 0.4751, 0.1764, 0.1444;
    # the rounded SDs squared agree to two significant figures). Covariance
    # CL-CL_taz = r * sqrt(var_CL * var_CL_taz) with r = 0.88.
    etalcl + etalcl_taz ~ c(
      log(1 + 0.6676^2),
      0.88 * sqrt(log(1 + 0.6676^2) * log(1 + 0.7798^2)),
      log(1 + 0.7798^2)
    ) # Table S6 omega_Cltol 0.61 (CV 66.76%), omega_Cltaz 0.69 (CV 77.98%), corr_Cltol_Cltaz 0.88
    etalvc ~ fixed(log(1 + 0.4392^2)) # Table S6 omega_V1 0.42 (CV 43.92%), no SE in Table S6
    etalvc_taz ~ fixed(log(1 + 0.3941^2)) # Table S6 omega_V2 0.38 (CV 39.41%), no SE in Table S6

    # Residual error. Table S5: Monolix 'Combined2' error, sd = sqrt(a^2 +
    # (b * f)^2), which is nlmixr2's add() + prop() with combined2(). Output 1
    # is ceftolozane, output 2 tazobactam.
    addSd <- 8.07; label("Ceftolozane additive residual error a1 (mg/L)") # Table S6 a1 = 8.07 (RSE 11.7%)
    propSd <- 0.23; label("Ceftolozane proportional residual error b1 (fraction)") # Table S6 b1 = 0.23 (RSE 10.7%)
    addSd_taz <- 1.3; label("Tazobactam additive residual error a2 (mg/L)") # Table S6 a2 = 1.3 (RSE 15.4%)
    propSd_taz <- 0.41; label("Tazobactam proportional residual error b2 (fraction)") # Table S6 b2 = 0.41 (RSE 6.89%)
  })

  model({
    # Individual parameters; clearances are power functions of the
    # BSA-de-indexed CKD-EPI 2021 eGFR centred at 60 mL/min.
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (CRCL / 60)^e_crcl_cl
    cl_taz <- exp(lcl_taz + etalcl_taz) * (CRCL / 60)^e_crcl_cl_taz
    vc <- exp(lvc + etalvc)
    vc_taz <- exp(lvc_taz + etalvc_taz)

    kel <- cl / vc
    kel_taz <- cl_taz / vc_taz

    # Supplement Figure S1: zero-order input of the combined dose into A1
    # (rate or dur on the dosing records), split F to tazobactam (A3) and
    # 1 - F to ceftolozane (A2).
    d/dt(depot) <- -ka * depot
    d/dt(central) <- (1 - frac) * ka * depot - kel * central
    d/dt(central_taz) <- frac * ka * depot - kel_taz * central_taz

    # Total plasma concentrations (mg/L). The paper's PK/PD targets use free
    # concentrations, 0.81 x Cc for ceftolozane and 0.70 x Cc_taz for
    # tazobactam (Methods); that conversion is applied in the vignette.
    Cc <- central / vc
    Cc_taz <- central_taz / vc_taz
    Cc ~ add(addSd) + prop(propSd) + combined2()
    Cc_taz ~ add(addSd_taz) + prop(propSd_taz) + combined2()
  })
}
