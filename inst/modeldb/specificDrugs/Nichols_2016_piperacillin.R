Nichols_2016_piperacillin <- function() {
  description <- "One-compartment population PK model for piperacillin in critically ill children (1-9 years) receiving extended-infusion piperacillin-tazobactam (Nichols 2016); IV zero-order input, first-order elimination, and a linear-additive effect of body weight on CL centered at the cohort median 18 kg."
  reference <- paste(
    "Nichols K, Chung EK, Knoderer CA, Buenger LE, Healy DP, Dees J,",
    "Crumby AS, Kays MB.",
    "Population pharmacokinetics and pharmacodynamics of extended-infusion",
    "piperacillin and tazobactam in critically ill children.",
    "Antimicrob Agents Chemother. 2016;60(1):522-531.",
    "doi:10.1128/AAC.02089-15.",
    sep = " "
  )
  vignette <- "Nichols_2016_piperacillin_tazobactam"
  units    <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Linear-additive effect on CL: TVCL = 3.51 + 0.0814 * (WT - 18) L/h,",
        "centered at the cohort median 18 kg (Nichols 2016 Results;",
        "individual weights 9.5-30.1 kg in Table 1). The paper does not",
        "model a weight effect on V."
      ),
      source_name        = "WT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 12,
    n_studies      = 1,
    age_range      = "12 months to 9 years (1-9 y in text; one 8-yr/9-yr cohort plus toddlers)",
    age_median     = "5 years (IQR 1.75-6.5)",
    weight_range   = "9.5-30.1 kg",
    weight_median  = "17.8 kg (IQR 11.4-20); model equation centers WT at 18 kg",
    sex_female_pct = 50,
    race_ethnicity = "Not reported",
    disease_state  = paste(
      "Critically ill children admitted to a pediatric intensive care unit",
      "with suspected or proven bacterial infection (most commonly pneumonia,",
      "VAP, sepsis, neutropenic fever); estimated GFR >= 60 mL/min/1.73 m^2",
      "(modified Schwartz)."
    ),
    dose_range     = paste(
      "100 mg/kg piperacillin (with 12.5 mg/kg tazobactam, 8:1 ratio) every 8",
      "hours by IV extended infusion over 4 hours (per institutional",
      "protocol). Maximum 3,000 mg piperacillin per dose. Observed doses",
      "1,070-3,375 mg of piperacillin per dose in this cohort (Table 1)."
    ),
    regions        = "United States (Riley Hospital for Children, Indianapolis, IN)",
    egfr_range     = "86-189 mL/min/1.73 m^2 (cohort median 103, IQR 96-111)",
    notes          = paste(
      "Twelve children sampled at steady state (6 samples per patient: pre-",
      "dose and at 2, 4 [end of infusion], 5, 6, and 8 hours after the start",
      "of the study dose). Patients had received a median of 5 prior doses",
      "(range 2-11) before the study dose. Patients receiving renal",
      "replacement therapy or with eGFR < 60 mL/min/1.73 m^2 were excluded.",
      "Demographics in Nichols 2016 Table 1."
    )
  )

  ini({
    # ===== Structural PK (Nichols 2016 Table 2 final piperacillin model) =====
    # Reference subject: WT = 18 kg (median centered).
    # TVCL = theta1 + theta3 * (WT - 18); TVV = theta2.
    lcl     <- log(3.51); label("Typical CL at WT = 18 kg (L/h)")  # Nichols 2016 Table 2: theta1 = 3.51 L/h (%SE 6.5)
    lvc     <- log(6.58); label("Typical V (L)")                    # Nichols 2016 Table 2: theta2 = 6.58 L (%SE 10.6)

    # Linear-additive WT effect on CL (Nichols 2016 Table 2 + Results final-model paragraph)
    e_wt_cl <- 0.0814;    label("Linear-additive effect of WT on CL (L/h/kg)")  # Nichols 2016 Table 2: theta3 = 0.0814 (%SE 45.1)

    # Reference covariate values (Nichols 2016 Results final-model paragraph)
    bw_ref  <- 18;        label("Reference body weight (kg, cohort centering value)")

    # ===== IIV (Nichols 2016 Table 2) =====
    # Exponential (log-normal) IIV; CV%-to-variance: omega^2 = log(1 + CV^2)
    # 17.3% -> log(1 + 0.173^2) = 0.029497
    # 25.2% -> log(1 + 0.252^2) = 0.061539
    # No correlation block: paper states "model did not support the correlation
    # between CL and V" (delta-OFV = -0.783).
    etalcl ~ 0.029497  # Nichols 2016 Table 2: omega_CL 17.3% CV (%SE 59.0, shrinkage 10.3%)
    etalvc ~ 0.061539  # Nichols 2016 Table 2: omega_V  25.2% CV (%SE 59.1, shrinkage 18.0%)

    # ===== Residual error (Nichols 2016 Table 2: proportional only) =====
    propSd <- 0.253; label("Proportional residual error (fraction)")  # Nichols 2016 Table 2: sigma_proportional 25.3% CV (%SE 28.7)
  })

  model({
    # ----- Individual PK parameters (Nichols 2016 Table 2 final piperacillin model) -----
    # TVCL = (theta1 + theta3 * (WT - 18)); CL_i = TVCL * exp(eta_CL)
    # IIV is multiplicative on the full typical-value CL, per the NONMEM
    # exponential-eta convention used throughout the paper.
    cl <- (exp(lcl) + e_wt_cl * (WT - bw_ref)) * exp(etalcl)
    vc <- exp(lvc + etalvc)

    # ----- Micro-constants -----
    kel <- cl / vc

    # ----- ODE system -----
    # IV piperacillin dosed into the central compartment (zero-order input
    # via the data-level RATE column; no depot compartment).
    d/dt(central) <- -kel * central

    # ----- Output -----
    # Piperacillin plasma concentration: dose in mg, vc in L -> mg/L.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
