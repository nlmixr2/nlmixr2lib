Nichols_2016_tazobactam <- function() {
  description <- "One-compartment population PK model for tazobactam in critically ill children (1-9 years) receiving extended-infusion piperacillin-tazobactam (Nichols 2016); IV zero-order input, first-order elimination, a multiplicative female-sex effect on CL, and a linear-additive WT effect on CL centered at the cohort median 18 kg."
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
        "Linear-additive effect on CL: TVCL = 3.43 * (1 - 0.285 * SEXF) +",
        "0.0676 * (WT - 18) L/h, centered at the cohort median 18 kg",
        "(Nichols 2016 Results; individual weights 9.5-30.1 kg in Table 1).",
        "The paper does not model a weight effect on V."
      ),
      source_name        = "WT"
    ),
    SEXF = list(
      description        = "Biological sex indicator (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = paste(
        "Multiplicative effect on CL: females have ~28.5% lower tazobactam",
        "CL than males (Nichols 2016 Table 2: theta3 = -0.285; entering the",
        "model as TVCL = theta1 * (1 + theta3 * SEXF) + theta4 * (WT - 18)).",
        "The source paper codes the indicator as 'sex' with 1 = female,",
        "0 = otherwise; this matches the canonical SEXF orientation with",
        "no value transformation. Reference category is 0 = male."
      ),
      source_name        = "sex"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 12,
    n_studies      = 1,
    age_range      = "12 months to 9 years (1-9 y in text)",
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
      "12.5 mg/kg tazobactam (with 100 mg/kg piperacillin, 8:1 ratio) every 8",
      "hours by IV extended infusion over 4 hours (per institutional",
      "protocol). Maximum 375 mg tazobactam per dose. Observed tazobactam",
      "doses 119-375 mg per dose in this cohort (calculated from the",
      "piperacillin doses in Table 1 / 8)."
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
    # ===== Structural PK (Nichols 2016 Table 2 final tazobactam model) =====
    # Reference subject: WT = 18 kg, male (SEXF = 0).
    # TVCL = theta1 * (1 + theta3 * SEXF) + theta4 * (WT - 18); TVV = theta2.
    lcl       <- log(3.43); label("Typical CL at WT = 18 kg, male (L/h)")  # Nichols 2016 Table 2: theta1 = 3.43 L/h (%SE 5.9)
    lvc       <- log(5.54); label("Typical V (L)")                          # Nichols 2016 Table 2: theta2 = 5.54 L (%SE 8.9)

    # Multiplicative female-sex effect on CL (Nichols 2016 Table 2)
    e_sexf_cl <- -0.285;    label("Multiplicative effect of female sex on CL (fraction)")  # Nichols 2016 Table 2: theta3 = -0.285 (%SE 20.9); Results: females have slower CL

    # Linear-additive WT effect on CL (Nichols 2016 Table 2)
    e_wt_cl   <- 0.0676;    label("Linear-additive effect of WT on CL (L/h/kg)")  # Nichols 2016 Table 2: theta4 = 0.0676 (%SE 38.6)

    # Reference covariate values (Nichols 2016 Results final-model paragraph)
    bw_ref    <- 18;        label("Reference body weight (kg, cohort centering value)")

    # ===== IIV (Nichols 2016 Table 2; IIV on V was not estimated) =====
    # Exponential (log-normal) IIV; CV%-to-variance: omega^2 = log(1 + CV^2)
    # 13.1% -> log(1 + 0.131^2) = 0.017028
    etalcl ~ 0.017028  # Nichols 2016 Table 2: omega_CL 13.1% CV (%SE 52.1, shrinkage 11.2%)

    # ===== Residual error (Nichols 2016 Table 2: combinational form) =====
    propSd <- 0.272; label("Proportional residual error (fraction)")  # Nichols 2016 Table 2: sigma_proportional 27.2% CV (%SE 35.2)
    addSd  <- 0.76;  label("Additive residual error (mg/L)")           # Nichols 2016 Table 2: sigma_additive 0.76 mg/L (%SE 47.8)
  })

  model({
    # ----- Individual PK parameters (Nichols 2016 Table 2 final tazobactam model) -----
    # TVCL = theta1 * (1 + theta3 * SEXF) + theta4 * (WT - 18);
    # with theta3 = -0.285 this reduces to TVCL = 3.43 * (1 - 0.285 * SEXF) +
    # 0.0676 * (WT - 18). IIV is multiplicative on the full typical-value CL,
    # per the NONMEM exponential-eta convention used throughout the paper.
    # No IIV on V (paper: omega_V dropped after addition did not significantly
    # decrease OFV; estimate had %SE > 100% and shrinkage > 30%).
    cl <- (exp(lcl) * (1 + e_sexf_cl * SEXF) + e_wt_cl * (WT - bw_ref)) * exp(etalcl)
    vc <- exp(lvc)

    # ----- Micro-constants -----
    kel <- cl / vc

    # ----- ODE system -----
    # IV tazobactam dosed into the central compartment (zero-order input
    # via the data-level RATE column; no depot compartment).
    d/dt(central) <- -kel * central

    # ----- Output -----
    # Tazobactam plasma concentration: dose in mg, vc in L -> mg/L.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
