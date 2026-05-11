Okada_2025_rocatinlimab <- function() {
  description <- "Two-compartment population PK model with parallel linear and time-dependent saturable (Michaelis-Menten) clearance and first-order subcutaneous absorption for rocatinlimab (anti-OX40 mAb) in adults; covariates body weight, albumin, plaque-psoriasis disease state, and healthy-volunteer cohort indicator (Okada 2025)"
  reference <- "Okada H, Liao S, Khouri L, Liao L, Hruska MW, Nagata Y, Hasegawa M, Gewitz A, Marsteller D. Continuous-Time Markov Population PK/PD Modeling of Longitudinal EASI Categorical Score in Atopic Dermatitis Treated With Rocatinlimab, an Anti-OX40 Monoclonal Antibody. CPT Pharmacometrics Syst Pharmacol. 2025;14(10):1587-1597. doi:10.1002/psp4.70069"
  vignette <- "Okada_2025_rocatinlimab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-allometric scaling (WT/70)^exponent on CL (0.923), V1 (0.828), and Vmax (0.494). Reference 70 kg adult. Treated as baseline (time-fixed) in the Okada 2025 analysis.",
      source_name        = "WT"
    ),
    ALB = list(
      description        = "Serum albumin (baseline)",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling (ALB/44)^-1.30 on linear CL. Reference 44 g/L (population median). Source column ALBU; the supplement reports albumin in g/L (mean 43.8, median 44, range 31-53).",
      source_name        = "ALBU"
    ),
    DIS_PSORIASIS = list(
      description        = "Plaque-psoriasis disease-state indicator (1 = psoriasis patient, 0 = non-psoriasis)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-psoriasis: atopic dermatitis, ulcerative colitis, or healthy volunteer)",
      notes              = "Multiplicative shift on linear CL: CL is multiplied by (1 + e_psoriasis_cl) = (1 - 0.372) for psoriasis patients. Source column DIS in the NONMEM control stream encodes 0=healthy, 1=psoriasis, 2=UC, 3=AD; ingestion sets DIS_PSORIASIS = as.integer(DIS == 1).",
      source_name        = "DIS == 1"
    ),
    DIS_HEALTHY = list(
      description        = "Healthy-volunteer cohort indicator (1 = healthy volunteer, 0 = patient)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (any patient: AD, UC, or psoriasis)",
      notes              = "Multiplicative shift on Vmax: Vmax is multiplied by (1 + e_healthy_vmax) = (1 - 0.532) for healthy volunteers. Source column DIS in the NONMEM control stream uses DIS=0 for healthy; ingestion sets DIS_HEALTHY = as.integer(DIS == 0).",
      source_name        = "DIS == 0"
    )
  )

  population <- list(
    n_subjects     = 413L,
    n_studies      = 5L,
    age_range      = "18-89 years (combined; pooled across 5 studies)",
    age_median     = 35,
    weight_range   = "38.0-166.0 kg",
    weight_median  = 70.3,
    sex_female_pct = 32.7,
    race_ethnicity = c(
      White_not_Hispanic_or_Latino   = 43.9,
      Asian                          = 52.2,
      Black_or_African               = 2.2,
      White_Hispanic_or_Latino       = 1.2,
      American_Indian_Alaskan_Native = 0.5
    ),
    disease_state = "Pooled cohort: 64.4% atopic dermatitis (Studies 4 and 5), 13.8% ulcerative colitis (Studies 2 and 3 partial), 13.1% plaque psoriasis (Study 1), 8.7% healthy volunteers (Study 3). Target indication for the model is moderate-to-severe atopic dermatitis.",
    dose_range = "0.003-10 mg/kg IV; 1 and 3 mg/kg SC; 150, 300, and 600 mg flat SC (Q2W or Q4W)",
    regions = "Multinational (multi-site clinical-trial program; race composition 52% Asian, 44% White, 2% Black or African, 1% White Hispanic/Latino, 0.5% American Indian/Alaskan Native)",
    notes = "Adult atopic-dermatitis target population (Studies 4 and 5, n = 266 AD patients): weight median 67-68 kg (mean 68-71 kg, range 38-166), age median 31-34 years (range 18-89), 18-41% female, 65-100% Asian. Albumin median 44 g/L (range 31-53), creatinine clearance median 117-126 mL/min. The full PPK dataset (n = 413) pooled five studies in psoriasis, UC, AD, and healthy volunteers (Tables S1, S2, and Table 1 of Okada 2025)."
  )

  ini({
    # Structural parameters (final estimates from supplement Table S3)
    lcl     <- log(0.230);   label("Linear (nonspecific) clearance for a 70 kg AD/UC patient with albumin 44 g/L, non-psoriasis (CL, L/day)")
    lvc     <- log(3.30);    label("Central volume of distribution for a 70 kg subject (V1, L)")
    lvp     <- log(2.82);    label("Peripheral volume of distribution (V2, L)")
    lq      <- log(0.775);   label("Intercompartmental clearance (Q, L/day)")
    lvmax   <- log(0.968);   label("Maximum velocity of saturable Michaelis-Menten clearance for a 70 kg patient (Vmax, mg/day)")
    lkm     <- log(0.289);   label("Michaelis-Menten constant of rocatinlimab serum concentration (Km, ug/mL)")
    lkdes   <- log(0.00439); label("First-order rate of decline of Vmax over time (Kdes, 1/day)")
    lka     <- log(0.312);   label("First-order absorption rate after subcutaneous administration (ka, 1/day)")
    lfdepot <- log(0.855);   label("Subcutaneous bioavailability (F3, fraction)")

    # Power exponents and disease shifts (Table S3)
    e_wt_cl        <-  0.923; label("Body-weight power exponent on CL (unitless; reference 70 kg)")
    e_wt_vc        <-  0.828; label("Body-weight power exponent on V1 (unitless; reference 70 kg)")
    e_wt_vmax      <-  0.494; label("Body-weight power exponent on Vmax (unitless; reference 70 kg)")
    e_alb_cl       <- -1.30;  label("Serum-albumin power exponent on CL (unitless; reference 44 g/L)")
    e_psoriasis_cl <- -0.372; label("Multiplicative shift on CL for psoriasis vs non-psoriasis (fractional)")
    e_healthy_vmax      <- -0.532; label("Multiplicative shift on Vmax for healthy volunteers vs patients (fractional)")

    # Inter-individual variability (BSV); omega^2 = log(CV^2 + 1)
    etalcl   ~ 0.04600  # 21.7% CV (Table S3)
    etalvc   ~ 0.02949  # 17.3% CV
    etalvp   ~ 0.08075  # 29.0% CV
    etalvmax ~ 0.50151  # 80.7% CV
    etalkdes ~ 0.13356  # 37.8% CV
    etalka   ~ 0.16319  # 42.1% CV

    # Residual error: NONMEM additive on log scale -> proportional in nlmixr2 (linear scale).
    # Table S3 reports two additive-on-log-scale SDs (SIGMA fixed at 1; SD switched by patient/healthy):
    #   ADD ERR1 = 9.3 %CV (healthy), ADD ERR2 = 16.4 %CV (patients).
    # Population PK target is the AD patient cohort, so the patient SD is used as the single propSd.
    propSd <- 0.164; label("Proportional residual error for patients (fraction; ADD ERR2 from Table S3 = 16.4% CV)")
  })

  model({
    # Covariate effects (Table S3 footnote, supplement section 2.2)
    cl_cov   <- (WT / 70)^e_wt_cl * (ALB / 44)^e_alb_cl * (1 + e_psoriasis_cl * DIS_PSORIASIS)
    v1_cov   <- (WT / 70)^e_wt_vc
    vmax_cov <- (WT / 70)^e_wt_vmax * (1 + e_healthy_vmax * DIS_HEALTHY)

    # Individual parameters
    cl   <- exp(lcl + etalcl)     * cl_cov
    vc   <- exp(lvc + etalvc)     * v1_cov
    vp   <- exp(lvp + etalvp)
    q    <- exp(lq)
    vmax <- exp(lvmax + etalvmax) * vmax_cov
    km   <- exp(lkm)
    kdes <- exp(lkdes + etalkdes)
    ka   <- exp(lka + etalka)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Time-dependent saturable clearance: TDVM = Vmax * exp(-Kdes * t).
    # `t` is simulation time (days from t = 0); align dosing so the first
    # rocatinlimab dose is at t = 0 to reproduce the supplement's profile.
    tdvm <- vmax * exp(-kdes * t)

    Cc <- central / vc

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1 - tdvm * Cc / (km + Cc)
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # SC bioavailability applies to the depot compartment; IV doses
    # entered into `central` bypass it.
    f(depot) <- exp(lfdepot)

    Cc ~ prop(propSd)
  })
}
