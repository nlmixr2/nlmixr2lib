Zhang_2021_dupilumab <- function() {
  description <- "Two-compartment population PK model for dupilumab in adult and adolescent patients with asthma (Zhang 2021), with first-order SC absorption and parallel linear plus Michaelis-Menten elimination from the central compartment."
  reference <- "Zhang L, Gao Y, Li M, et al. Population pharmacokinetic analysis of dupilumab in adult and adolescent patients with asthma. CPT Pharmacometrics Syst Pharmacol. 2021;10(9):941-952. doi:10.1002/psp4.12667"
  vignette <- "Zhang_2021_dupilumab"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight (time-varying)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on Ke, V2 (central volume), and Vmax. Reference 78 kg, the median of the final asthma popPK dataset (Zhang 2021 Results, p. 945).",
      source_name        = "WT"
    ),
    ALB = list(
      description        = "Serum albumin (baseline)",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on V2. Reference 44 g/L, the median of the final asthma popPK dataset (Zhang 2021 Results, p. 945).",
      source_name        = "ALB"
    ),
    CRCL = list(
      description        = "Creatinine clearance normalized to body surface area (Cockcroft-Gault, BSA-normalized as 1.73 * CrCl / BSA)",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on Ke. Reference 111 mL/min/1.73 m^2, the median of the final asthma popPK dataset (Zhang 2021 Results, p. 945). Source column is CLCRN per Zhang 2021 Methods.",
      source_name        = "CLCRN"
    ),
    ADA_POS = list(
      description        = "Stationary anti-drug antibody (ADA) positivity indicator (positive at any time on study)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (ADA-negative; typical patient)",
      notes              = "Time-fixed per subject: 1 = positive at any time during the study, 0 = always negative. Multiplicative effect on Ke as Ke * (1 + 0.191 * ADA_POS). Per Zhang 2021 Results, the typical patient is ADA-negative; ADA-positive subjects make up 14.5% of the asthma development population (Table 2).",
      source_name        = "ADA"
    )
  )

  population <- list(
    n_subjects     = 2114L,
    n_observations = 14584L,
    n_studies      = 9L,
    age_range      = "12-83 years",
    age_median     = "48 years (asthma development cohort, N = 2114; healthy median 32, asthma median 49)",
    weight_range   = "32-186 kg (asthma); 52-95 kg (healthy)",
    weight_median  = "78 kg (final dataset median; healthy 77.3, asthma 80.0)",
    sex_female_pct = 59.6,
    race_ethnicity = "Race was a tested covariate; detailed race breakdown not reported in the main text Table 2.",
    disease_state  = "Pooled cohort of 202 healthy adults and 1912 patients with moderate-to-severe asthma (1844 adults + 68 adolescents aged 12 to <18 years).",
    dose_range     = "1-12 mg/kg IV and 75-600 mg SC single dose (Phase I); 200-300 mg SC q2w or q4w with 400-600 mg SC loading dose (Phase II/III maintenance).",
    regions        = "Multi-regional Phase I-III programme; 9 pooled studies (NCT01015027, NCT01484600, NCT01537653, NCT01537640, PKM14161, PKM14271, NCT01312961, NCT01854047, NCT02414854).",
    notes          = "Baseline demographics from Zhang 2021 Table 2 (final dataset N = 2114 with 14,584 dupilumab concentrations). Adolescents (12 to <18 years) N = 68 (3.2%). ADA-positive 14.5% of asthma cohort. Albumin median 44 g/L, CrCl-normalized median 111 mL/min/1.73 m^2. Eosinophil, FeNO, FEV1 percent-of-predicted-normal were tested but had no significant effect on dupilumab PK. The independent evaluation cohort (NCT02528214, N = 103, severe OCS-dependent asthma) was used for external validation only and is not included in n_subjects."
  )

  ini({
    # Structural PK parameters - Zhang 2021 Table 3 final-model estimates (reference covariate
    # values: 78 kg body weight, 44 g/L albumin, 111 mL/min/1.73 m^2 creatinine clearance,
    # ADA-negative).
    lka     <- log(0.263);  label("Absorption rate Ka (1/day)")                                    # Zhang 2021 Table 3, Ka row
    lkel    <- log(0.0418); label("Linear elimination rate Ke (1/day)")                            # Zhang 2021 Table 3, Ke row
    lvc     <- log(2.76);   label("Central compartment volume V2 (L)")                             # Zhang 2021 Table 3, V2 row
    lk12    <- log(0.0952); label("Central-to-peripheral rate K23 (1/day)")                        # Zhang 2021 Table 3, K23 row
    lk21    <- log(0.163);  label("Peripheral-to-central rate K32 (1/day)")                        # Zhang 2021 Table 3, K32 row
    lvmax   <- log(1.39);   label("Maximum target-mediated elimination rate Vmax (mg/L/day)")      # Zhang 2021 Table 3, Vmax row
    lkm     <- log(2.08);   label("Michaelis-Menten constant Km (mg/L)")                           # Zhang 2021 Table 3, Km row
    lfdepot <- log(0.609);  label("Subcutaneous bioavailability Fsc (fraction)")                   # Zhang 2021 Table 3, Fsc row

    # Covariate effects - Zhang 2021 Results (p. 945), final-model equations for V2, Vmax, Ke.
    e_wt_kel    <-  0.222;  label("Power exponent of WT/78 on Ke (unitless)")                      # Zhang 2021 Table 3: weight effect on Ke
    e_wt_vc     <-  0.667;  label("Power exponent of WT/78 on V2 (unitless)")                      # Zhang 2021 Table 3: weight effect on V2
    e_wt_vmax   <-  0.224;  label("Power exponent of WT/78 on Vmax (unitless)")                    # Zhang 2021 Table 3: weight effect on Vmax
    e_alb_vc    <- -0.484;  label("Power exponent of ALB/44 on V2 (unitless)")                     # Zhang 2021 Table 3: albumin effect on V2
    e_crcl_kel  <-  0.217;  label("Power exponent of CRCL/111 on Ke (unitless)")                   # Zhang 2021 Table 3: CLCRN effect on Ke
    e_ada_kel   <-  0.191;  label("Proportional multiplicative effect of ADA-positive on Ke (unitless)")  # Zhang 2021 Table 3: ADA effect on Ke; coded as Ke*(1 + 0.191*ADA)

    # Inter-individual variability: Zhang 2021 Table 3 reports omega^2 (variance on the log
    # scale) directly under the "Estimate" column; the percent in parentheses is the
    # small-variance approximation CV(%) ~ sqrt(omega^2) used by the authors. We use the
    # omega^2 values verbatim here.
    etalkel    ~ 0.0385    # Zhang 2021 Table 3: Ke IIV omega^2 = 0.0385 (CV ~19.6%, shrinkage 47.3%)
    etalvc     ~ 0.00834   # Zhang 2021 Table 3: V2 IIV omega^2 = 0.00834 (CV ~9.13%, shrinkage 47.7%)
    etalvmax   ~ 0.0589    # Zhang 2021 Table 3: Vmax IIV omega^2 = 0.0589 (CV ~24.3%, shrinkage 57.7%)
    etalka     ~ 0.243     # Zhang 2021 Table 3: Ka IIV omega^2 = 0.243 (CV ~49.2%, shrinkage 57.6%)
    etalfdepot ~ 0.132     # Zhang 2021 Table 3: Fsc IIV omega^2 = 0.132 (CV ~36.3%, shrinkage 36.3%)

    # Residual error: Zhang 2021 Table 3 reports sigma^2 as the "Estimate" with the
    # corresponding SD shown in parentheses under "(CV%)" / "SD". Combined proportional +
    # additive error model.
    propSd <- 0.197;  label("Proportional residual error (fraction)")                              # Zhang 2021 Table 3: proportional sigma^2 = 0.0388, SD = sqrt(0.0388) = 0.197 (~19.7% CV)
    addSd  <- 1.73;   label("Additive residual error (mg/L)")                                      # Zhang 2021 Table 3: additive sigma^2 = 2.98 mg^2/L^2, SD = sqrt(2.98) = 1.73 mg/L
  })
  model({
    # Individual PK parameters with Zhang 2021 covariate models (typical reference patient:
    # 78 kg, 44 g/L albumin, 111 mL/min/1.73 m^2 creatinine clearance, ADA-negative).
    kel <- exp(lkel + etalkel) *
           (1 + e_ada_kel * ADA_POS) *
           (WT / 78)^e_wt_kel *
           (CRCL / 111)^e_crcl_kel
    vc  <- exp(lvc + etalvc) *
           (WT / 78)^e_wt_vc *
           (ALB / 44)^e_alb_vc
    vmax <- exp(lvmax + etalvmax) *
            (WT / 78)^e_wt_vmax
    km <- exp(lkm)
    ka <- exp(lka + etalka)
    k12 <- exp(lk12)
    k21 <- exp(lk21)
    fdepot <- exp(lfdepot + etalfdepot)

    # Two-compartment PK with first-order SC absorption, parallel linear elimination (kel)
    # and Michaelis-Menten elimination (Vmax / Km) from the central compartment. Vmax is
    # reported in mg/L/day, so the MM rate in amount units (mg/day) is
    # vmax * vc * Cc / (km + Cc), equivalent to vmax * central / (km + Cc).
    Cc <- central / vc

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot -
                          kel * central -
                          vmax * central / (km + Cc) -
                          k12 * central +
                          k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    f(depot) <- fdepot

    Cc ~ prop(propSd) + add(addSd)
  })
}
