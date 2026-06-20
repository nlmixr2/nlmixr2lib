Hwang_2023_monalizumab <- function() {
  description <- "Two-compartment population PK model for monalizumab (anti-CD94/NKG2A IgG4) in patients with advanced solid tumors or squamous cell carcinoma of the head and neck (Hwang 2023)"
  reference <- "Hwang M, Fan C, Yue MS, Zhou D, Paturel C, Andre P, Cheng L-Y, Mitchell P, Kourtesis P, Ruscica D, Das M, Morsli N, Ren S, Gibbs M, Phipps A, Song X. Population Pharmacokinetics of Monalizumab in Patients With Advanced Solid Tumors. J Clin Pharmacol. 2023;63(7):818-829. doi:10.1002/jcph.2220"
  vignette <- "Hwang_2023_monalizumab"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form effect on CL (exponent 0.626) and V1 (exponent 0.495) with reference 70.6 kg, the median of the pooled analysis population (Hwang 2023 Table 1, p. 821; covariate equation in Methods 'Covariate Considerations', p. 820: P_ki = theta_k * (X_ij / M(X_j))^theta_j). Source column BLWT (baseline body weight); baseline-only.",
      source_name        = "BLWT"
    ),
    ALB = list(
      description        = "Baseline serum albumin",
      units = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form effect on CL (exponent -1.20) and V1 (exponent -0.299) with reference 3.80 g/dL, the median of the pooled analysis population (Hwang 2023 Table 1, p. 821). Source column BALB (baseline albumin); Hwang 2023 reports albumin in g/dL throughout (Table 1, range 1.80-4.90 g/dL).",
      source_name        = "BALB"
    ),
    SEXF = list(
      description        = "Biological sex indicator (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male) — matches the paper's reference category. The Hwang 2023 Table 2 coefficient +0.111 means women have V1 ~11.1% higher than men.",
      notes              = "Proportional-shift effect on V1: V1 multiplied by (1 + 0.111)^SEXF (Hwang 2023 Table 2, row 'Coefficient in the proportional shift model of women on V1', p. 821). The categorical covariate equation is in Methods 'Covariate Considerations', p. 820: P_ki = theta_k * (1 + theta_j)^X_ij with X = 1 for the non-reference category. Storage convention SEXF (1 = female) matches the paper's encoding (women = 1), so no value transformation is needed.",
      source_name        = "SEX"
    ),
    SMOKE_CURRENT = list(
      description        = "1 if current smoker at baseline, 0 otherwise (former or never)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (former smoker, paired with SMOKE_NEVER = 0). Smoking-status reference in Hwang 2023 is former smoker, the most-common category (n = 319/507 = 62.9%, Table 1, p. 822).",
      notes              = "Proportional-shift effect on V1: V1 multiplied by (1 + 0.0484)^SMOKE_CURRENT (Hwang 2023 Table 2, row 'Coefficient in the proportional shift model of current smoker on V1', p. 821). Paired with SMOKE_NEVER to encode the 3-level smoking-status categorical (never / former / current); both indicators = 0 yields the former-smoker reference. See `inst/references/covariate-columns.md` for the canonical SMOKE_CURRENT entry.",
      source_name        = "Smoking status (Current vs. Former reference)"
    ),
    SMOKE_NEVER = list(
      description        = "1 if never smoker at baseline, 0 otherwise (former or current)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (former smoker, paired with SMOKE_CURRENT = 0).",
      notes              = "Proportional-shift effect on V1: V1 multiplied by (1 - 0.141)^SMOKE_NEVER (Hwang 2023 Table 2, row 'Coefficient in the proportional shift model of never smoker on V1', p. 821). Paired with SMOKE_CURRENT to encode the 3-level smoking-status categorical (never / former / current); both indicators = 0 yields the former-smoker reference.",
      source_name        = "Smoking status (Never vs. Former reference)"
    )
  )

  population <- list(
    n_subjects     = 507L,
    n_studies      = 2L,
    n_observations = 2842L,
    age_range      = "23.0-91.0 years",
    age_median     = "60.0 years",
    weight_range   = "37.4-154 kg",
    weight_median  = "70.6 kg",
    sex_female_pct = 36.9,
    race_ethnicity = c(
      White                                 = 66.9,
      `Black or African American`           = 5.1,
      Asian                                 = 11.0,
      `American Indian or Alaskan Native`   = 0.2,
      Other                                 = 3.2,
      Missing                               = 13.6
    ),
    disease_state  = "Advanced solid tumors (microsatellite-stable colorectal cancer 48.1%, squamous cell carcinoma of the head and neck 27.2%, ovarian cancer 8.5%, microsatellite-stable endometrial cancer 8.5%, non-small-cell lung cancer 3.9%, cervical cancer 3.2%, castration-resistant prostate cancer 0.6%)",
    dose_range     = "22.5-750 mg flat or 0.4-10 mg/kg weight-based, IV every 2 weeks (Q2W); alternative regimens 750 mg or 1500 mg IV every 4 weeks (Q4W)",
    regions        = "Multinational; the cohort includes a 40-patient East-Asia site, US, and EU sites",
    smoking_status = "Never 6.5%; Former 62.9%; Current 28.4%; Missing 2.2%",
    ada_status     = "Negative 48.5%; Positive 13.2%; Missing 38.3% (study IPH2201-203 still ongoing at data cutoff; ADA was tested but not retained as a covariate)",
    ecog_distribution = "ECOG 0 43.2%; ECOG 1 56.4%; Missing 0.4%",
    notes          = "Pooled from studies D419NC00001 (NCT02671435; N = 369) and IPH2201-203 (NCT02643550; N = 138); data cutoffs 2021-05-06 and 2021-05-13 respectively (Hwang 2023 Methods 'Study and Data Set', p. 819). Baseline demographics per Hwang 2023 Table 1, pp. 821-822. The final-model dataset excluded BLQ samples after the first dose using the Beal M1 method (224/3066 = 7.3% of all observations; Results 'Pharmacokinetic Data Set', p. 823). The 'Other' race category aggregates the four small categories per Table 1."
  )

  ini({
    # Structural parameters - typical values for a subject at the population
    # median covariates (WT = 70.6 kg, ALB = 3.80 g/dL, male, former smoker).
    # CL and Q reported in L/day; V1 and V2 in L. The model keeps time in days.
    lcl <- log(0.255); label("Clearance for the typical adult at population-median covariates (CL, L/day)") # Hwang 2023 Table 2, p. 821: CL = 0.255 L/day
    lvc <- log(3.58);  label("Central volume of distribution at population-median covariates (V1, L)")     # Hwang 2023 Table 2, p. 821: V1 = 3.58 L
    lq  <- log(0.567); label("Intercompartmental clearance (Q, L/day)")                                    # Hwang 2023 Table 2, p. 821: Q  = 0.567 L/day
    lvp <- log(2.78);  label("Peripheral volume of distribution (V2, L)")                                  # Hwang 2023 Table 2, p. 821: V2 = 2.78 L

    # Covariate effects on CL (Hwang 2023 Table 2; continuous power form per
    # Methods 'Covariate Considerations', p. 820: P = theta * (X / M(X))^theta_j).
    e_alb_cl <- -1.20;  label("Power exponent of ALB on CL (unitless; reference ALB 3.80 g/dL)")  # Hwang 2023 Table 2, p. 821: BALB on CL = -1.20
    e_wt_cl  <-  0.626; label("Power exponent of WT on CL (unitless; reference WT 70.6 kg)")      # Hwang 2023 Table 2, p. 821: BLWT on CL =  0.626

    # Covariate effects on V1 (Hwang 2023 Table 2). Continuous covariates use the
    # power form; categorical covariates use the proportional-shift form
    # P = theta * (1 + theta_j)^X (Methods 'Covariate Considerations', p. 820).
    e_alb_vc           <- -0.299;  label("Power exponent of ALB on V1 (unitless; reference ALB 3.80 g/dL)")           # Hwang 2023 Table 2, p. 821: BALB on V1 = -0.299
    e_wt_vc            <-  0.495;  label("Power exponent of WT on V1 (unitless; reference WT 70.6 kg)")               # Hwang 2023 Table 2, p. 821: BLWT on V1 =  0.495
    e_sexf_vc          <-  0.111;  label("Proportional-shift coefficient of women (SEXF = 1) on V1 (unitless)")       # Hwang 2023 Table 2, p. 821: women on V1 =  0.111
    e_smoke_current_vc <-  0.0484; label("Proportional-shift coefficient of current smoker on V1 (unitless)")         # Hwang 2023 Table 2, p. 821: current smoker on V1 =  0.0484
    e_smoke_never_vc   <- -0.141;  label("Proportional-shift coefficient of never smoker on V1 (unitless)")           # Hwang 2023 Table 2, p. 821: never smoker on V1   = -0.141

    # Inter-individual variability. CL, V1, V2 form a 3x3 log-normal block; the
    # paper does not estimate IIV on Q (Hwang 2023 Results 'Base Structural Model',
    # p. 823: 'Inclusion of IIV on intercompartmental clearance was not supported
    # by the data, yielding increased objective function value and an estimate
    # that was close to 0 with high relative standard error'). Variance /
    # covariance values (lower-triangular, row-major):
    #   row 1:  omega^2_CL                                = 0.136
    #   row 2:  cov(CL,V1),   omega^2_V1                  = 0.055,  0.0618
    #   row 3:  cov(CL,V2),   cov(V1,V2),   omega^2_V2    = 0.074,  0.0904, 0.318
    # Source: Hwang 2023 Table 2, p. 821 (rows 'Variance of IIV on CL/V1/V2' and
    # 'Covariance of CL and V1', 'Covariance of CL and V2', 'Covariance of V1 and V2').
    etalcl + etalvc + etalvp ~ c(
      0.136,
      0.055, 0.0618,
      0.074, 0.0904, 0.318
    )

    # Combined additive + proportional residual error model
    # Y_ij = C_ij * (1 + eps1_ij) + eps2_ij  (Hwang 2023 Methods 'Base Model
    # Development', p. 820). Hwang 2023 Table 2 reports VARIANCES of the residual
    # error terms (omega^2 form), so the SDs used by nlmixr2's add()/prop()
    # syntax are sqrt(variance).
    propSd <- sqrt(0.0772);  label("Proportional residual error SD (fraction)")  # Hwang 2023 Table 2, p. 821: Variance of proportional error = 0.0772
    addSd  <- sqrt(0.00766); label("Additive residual error SD (mg/L)")          # Hwang 2023 Table 2, p. 821: Variance of additive error    = 0.00766
  })
  model({
    # SI -> US-convention unit conversion (canonical ALB is in SI g/L per the
    # 2026-06-19 register standardization audit; the original calibration
    # used the g/dL reference value, so convert inline here).
    alb_gdL <- ALB * 0.1  # SI g/L -> US-convention g/dL (factor 0.1)

    # Covariate factors. Reference values are the population medians from
    # Hwang 2023 Table 1 (p. 821): alb_gdL = 3.80 g/dL, WT = 70.6 kg.
    alb_cl <- (alb_gdL / 3.80)^e_alb_cl
    wt_cl  <- (WT  / 70.6)^e_wt_cl

    alb_vc           <- (alb_gdL / 3.80)^e_alb_vc
    wt_vc            <- (WT  / 70.6)^e_wt_vc
    sexf_vc          <- (1 + e_sexf_vc)^SEXF
    smoke_current_vc <- (1 + e_smoke_current_vc)^SMOKE_CURRENT
    smoke_never_vc   <- (1 + e_smoke_never_vc)^SMOKE_NEVER

    # Individual PK parameters (Hwang 2023 Table 2 + covariate equations)
    cl <- exp(lcl + etalcl) * alb_cl * wt_cl
    vc <- exp(lvc + etalvc) * alb_vc * wt_vc * sexf_vc * smoke_current_vc * smoke_never_vc
    vp <- exp(lvp + etalvp)
    q  <- exp(lq)

    # Two-compartment IV micro-constants (Hwang 2023 Results 'Base Structural
    # Model', p. 823: 2-compartment with first-order elimination)
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # Dose in mg, vc in L => central / vc has units mg/L (= ug/mL).
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
