Yoneyama_2017_emicizumab <- function() {
  description <- "One-compartment population PK model with first-order subcutaneous absorption and elimination for emicizumab (ACE910), a bispecific anti-FIXa/FX humanized monoclonal antibody mimicking the cofactor function of activated factor VIII, in healthy male adult volunteers (Japanese and Caucasian) and Japanese male adult/adolescent patients with severe hemophilia A with or without factor VIII inhibitors (Yoneyama 2017). Body-weight allometric exponents are fixed (0.75 on CL/F, 1 on Vd/F) per Yoneyama 2017 Methods. Anti-emicizumab neutralizing antibody (ADA_POS) increases CL/F by a factor of exp(2.01) and the effect onsets 33.4 days post the first SC dose (NONMEM MTIME parameterisation). The companion repeated time-to-event (RTTE) bleeding-hazard model from Yoneyama 2017 Section 2.4 is not included here; nlmixr2lib does not currently support TTE models."
  reference   <- "Yoneyama K, Schmitt C, Kotani N, Levy GG, Kasai R, Iida S, Shima M, Kawanishi T. A Pharmacometric Approach to Substitute for a Conventional Dose-Finding Study in Rare Diseases: Example of Phase III Dose Selection for Emicizumab in Hemophilia A. Clin Pharmacokinet. 2018 May;57(5):613-619. doi:10.1007/s40262-017-0616-3 (online 6 December 2017; PMID 29209893)."
  vignette    <- "Yoneyama_2017_emicizumab"
  units       <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric power scaling with theoretical exponents fixed at 0.75 on CL/F and 1 on Vd/F per Yoneyama 2017 Methods. Reference weight 70 kg per Table 2 footnote b ('Standardized for a 70 kg, NAb-negative healthy volunteer'). Pooled-cohort range across healthy volunteers and patients: 41-87 kg (Table 1).",
      source_name        = "BW"
    ),
    DIS_HEALTHY = list(
      description        = "Healthy-participant cohort indicator: 1 = healthy male adult volunteer (single SC dose), 0 = Japanese male adult/adolescent patient with severe hemophilia A (multiple SC doses).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (severe hemophilia A patient)",
      notes              = "Time-fixed per subject. Yoneyama 2017 Eq. 2 (paper) encodes the patient-vs-volunteer effect with a PATIENT flag (1 = patient) and reports the typical CL/F = 0.222 L/day and Vd/F = 10.2 L standardized for a healthy volunteer (Table 2 footnote b). The canonical DIS_HEALTHY convention uses 0 = patient as reference; the model file shifts the structural typical values (lcl, lvc) to the patient state and negates the covariate coefficients so the model is mathematically identical: at DIS_HEALTHY = 0, CL/F = 0.280 L/day = 0.222 * exp(0.232) and Vd/F = 12.15 L = 10.2 * exp(0.175); at DIS_HEALTHY = 1, the exp(-0.232) and exp(-0.175) factors restore the paper's HV-typical 0.222 L/day and 10.2 L.",
      source_name        = "PATIENT"
    ),
    ADA_POS = list(
      description        = "Anti-emicizumab neutralizing antibody (NAb) positivity indicator: 1 = subject has developed neutralizing anti-emicizumab antibodies, 0 = NAb-negative.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (NAb-negative)",
      notes              = "Time-fixed per subject in this implementation. The NAb effect on CL/F is gated on time-since-first-dose >= 33.4 days via the structural parameter tNab (Yoneyama 2017 Table 2: 'Onset time of the effect of NAb positivity on CL/F following a single SC administration of emicizumab'; parameterised as MTIME in NONMEM). Of the 60 subjects in the analysis dataset only 1 Caucasian healthy volunteer developed NAbs (Table 1); the NAb effect estimate exp(2.01) = ~7.5x increase on CL/F therefore comes from one positive subject and should be interpreted with the wide 95% CI (1.26-2.74) reported in Table 2.",
      source_name        = "NAb"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 60L,
    n_studies        = 3L,
    age_range        = "12-58 years (healthy volunteers 21-44 years; patients 12-58 years)",
    age_median       = "Japanese HV 31; Caucasian HV 29; Cohort 1 patients 32; Cohort 2 patients 30; Cohort 3 patients 33 years",
    weight_range     = "41-87 kg (pooled across healthy volunteers and patients)",
    weight_median    = "Japanese HV 60; Caucasian HV 70; Cohort 1 60; Cohort 2 56; Cohort 3 66 kg",
    sex_female_pct   = 0,
    race_ethnicity   = "Healthy volunteers: 24 Japanese + 18 Caucasian male adults. Patients: 18 Japanese male adult/adolescent.",
    disease_state    = "Severe hemophilia A (factor VIII activity < 1 IU/dL) with or without FVIII inhibitors in the patient arm. Cohort 1 (n=6): 1 mg/kg SC loading followed by 0.3 mg/kg QW SC (4 FVIII-inhibitor-positive, 2 negative). Cohort 2 (n=6): 3 mg/kg SC loading followed by 1 mg/kg QW SC (4 FVIII-inhibitor-positive, 2 negative). Cohort 3 (n=6): 3 mg/kg QW SC (3 FVIII-inhibitor-positive, 3 negative). Healthy volunteers received single SC doses of 0.01, 0.1, 0.3, or 1 mg/kg (Japanese) or 0.1, 0.3, or 1 mg/kg (Caucasian); the 0.001 mg/kg Japanese cohort (n=6) was excluded because all plasma concentrations were below the lower limit of quantification (< 0.05 ug/mL).",
    dose_range       = "0.01-3 mg/kg SC (healthy volunteers, single doses 0.01, 0.1, 0.3, 1 mg/kg; patients, multiple doses 0.3-3 mg/kg QW with or without a loading dose). Phase III dosing regimens proposed from the simulations: 3 mg/kg QW loading for 4 weeks followed by 1.5 mg/kg QW, 3 mg/kg Q2W, or 6 mg/kg Q4W.",
    regions          = "Japan and Switzerland/EU (Japanese single-ascending-dose and patient studies; Caucasian single-ascending-dose study).",
    nab_pos_pct      = "1/60 (one Caucasian healthy volunteer developed NAbs)",
    fviii_inh_pos_pct = "11/18 patients (61%) had FVIII inhibitors at baseline.",
    notes            = "Demographic and cohort details from Yoneyama 2017 Table 1. The PopPK model pooled 42 healthy volunteers (single SC doses) with 18 severe hemophilia A patients (multiple SC doses). Phase III dose selection used this PopPK model jointly with the RTTE bleeding-hazard model (Section 2.4 of Yoneyama 2017, not implemented here). The plasma emicizumab concentration target identified for phase III was >= 45 ug/mL at steady-state trough."
  )

  ini({
    # ----------------------------------------------------------------------
    # Structural parameters - Yoneyama 2017 Table 2 final-model estimates.
    # Paper standardizes typical values for a 70 kg, NAb-negative HEALTHY
    # VOLUNTEER (Table 2 footnote b). The canonical DIS_HEALTHY register
    # uses 0 = patient as the reference; the structural typicals below are
    # shifted to that reference state so that:
    #   * at DIS_HEALTHY = 0, WT = 70, ADA_POS = 0 (patient typical):
    #       CL/F = exp(lcl) = 0.222 * exp(0.232) = 0.280 L/day
    #       Vd/F = exp(lvc) = 10.2  * exp(0.175) = 12.15 L
    #   * the e_dis_healthy_* coefficients in the model() block subtract
    #     0.232 / 0.175 (i.e., apply exp(-0.232) / exp(-0.175)) to recover
    #     the paper's HV typicals at DIS_HEALTHY = 1.
    # Half-life of absorption t_(1/2),abs = 1.56 days -> ka = ln(2) / 1.56
    #   = 0.4443 1/day; lka = log(0.4443).
    # ----------------------------------------------------------------------
    lka <- log(log(2) / 1.56)         ; label("First-order SC absorption rate ka (1/day); derived from t_(1/2),abs = 1.56 days")  # Yoneyama 2017 Table 2: t_(1/2),abs = 1.56 days
    lcl <- log(0.222) + 0.232         ; label("Apparent clearance CL/F (L/day); patient typical at WT = 70 kg, ADA-negative")     # Yoneyama 2017 Table 2: CL/F = 0.222 L/day (HV-typical) + log of paper's PATIENT effect 0.232
    lvc <- log(10.2) + 0.175          ; label("Apparent central volume Vd/F (L); patient typical at WT = 70 kg")                   # Yoneyama 2017 Table 2: Vd/F = 10.2 L (HV-typical) + log of paper's PATIENT effect 0.175

    # Allometric body-weight scaling - exponents fixed at theoretical values.
    e_wt_cl <- fixed(0.75)            ; label("Allometric exponent on CL/F (unitless)")                                            # Yoneyama 2017 Table 2: Effect of BW on CL/F = 0.75 (fixed)
    e_wt_vc <- fixed(1)               ; label("Allometric exponent on Vd/F (unitless)")                                            # Yoneyama 2017 Table 2: Effect of BW on Vd/F = 1 (fixed)

    # Categorical covariate effects (paper Eq. 2 form: P_i = P_TV * exp(Cov * theta)).
    # Coefficients are negated relative to the paper because DIS_HEALTHY = 1 - paper's PATIENT.
    e_dis_healthy_cl <- -0.232        ; label("Effect of DIS_HEALTHY = 1 on log-CL/F (unitless); -1 * paper PATIENT effect")        # Yoneyama 2017 Table 2: Effect of patient on CL/F = 0.232 (PATIENT flag = 1 - DIS_HEALTHY)
    e_dis_healthy_vc <- -0.175        ; label("Effect of DIS_HEALTHY = 1 on log-Vd/F (unitless); -1 * paper PATIENT effect")        # Yoneyama 2017 Table 2: Effect of patient on Vd/F = 0.175 (PATIENT flag = 1 - DIS_HEALTHY)

    # Anti-emicizumab NAb effect on CL/F - exponential form per paper Eq. 2,
    # with onset time tNab post the first SC dose (MTIME parameterisation in
    # NONMEM; see Yoneyama 2017 Methods 'PopPK Modeling').
    e_ada_pos_cl <- 2.01              ; label("Effect of ADA_POS = 1 (NAb-positive) on log-CL/F once the NAb effect is active (unitless)")  # Yoneyama 2017 Table 2: Effect of NAb positivity on CL/F = 2.01
    tNab         <- 33.4              ; label("Onset time (days) of the NAb effect on CL/F after the first SC dose of emicizumab")          # Yoneyama 2017 Table 2: Onset time of NAb effect = 33.4 days

    # ----------------------------------------------------------------------
    # Inter-individual variability - Yoneyama 2017 Table 2 (variances on a
    # log-normal scale per paper Methods 'Inter-individual variability
    # following a log-normal distribution was assumed for each structural
    # model parameter'). CL/F and Vd/F share a covariance term. IIV on the
    # absorption half-life t_(1/2),abs maps to IIV on lka with identical
    # variance magnitude because ka = ln(2) / t_(1/2),abs implies
    # eta_lka = -eta_l_t_half_abs, var(eta_lka) = var(eta_l_t_half_abs).
    # ----------------------------------------------------------------------
    etalcl + etalvc ~ c(0.0737,
                        0.0278, 0.0455)  # Yoneyama 2017 Table 2: Var(CL/F) = 0.0737; Cov(CL/F, Vd/F) = 0.0278; Var(Vd/F) = 0.0455
    etalka          ~ 0.502              # Yoneyama 2017 Table 2: Var(t_(1/2),abs) = 0.502; equal magnitude variance on lka (sign of eta is irrelevant for a zero-mean random variable)

    # ----------------------------------------------------------------------
    # Residual unexplained variability - Yoneyama 2017 Table 2 (combined
    # additive-plus-proportional error model; additive parameterised as SD
    # in ug/mL, proportional as CV).
    # ----------------------------------------------------------------------
    propSd <- 0.128                   ; label("Proportional residual error (fraction)")                                            # Yoneyama 2017 Table 2: Proportional error = 12.8%
    addSd  <- 0.0149                  ; label("Additive residual error (ug/mL)")                                                   # Yoneyama 2017 Table 2: Additive error = 0.0149 ug/mL
  })

  model({
    # NAb effect activation flag - 1 if subject is ADA_POS AND time since
    # simulation start (= time since first dose for a t = 0 first dose)
    # exceeds the estimated onset time tNab. This reproduces the NONMEM
    # MTIME/MPAST step that switches on the NAb effect 33.4 days after the
    # first SC administration (Yoneyama 2017 Methods 'PopPK Modeling').
    nab_active <- ADA_POS * (t > tNab)

    # Individual PK parameters - paper Eq. 1 (continuous covariate, power
    # form) and Eq. 2 (categorical covariate, exponential form) combined
    # multiplicatively.
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl + e_dis_healthy_cl * DIS_HEALTHY + e_ada_pos_cl * nab_active) *
          (WT / 70) ^ e_wt_cl
    vc <- exp(lvc + etalvc + e_dis_healthy_vc * DIS_HEALTHY) *
          (WT / 70) ^ e_wt_vc

    # Elimination rate constant for the one-compartment system.
    kel <- cl / vc

    # One-compartment SC PK with first-order absorption.
    d / dt(depot)   <- -ka * depot
    d / dt(central) <-  ka * depot - kel * central

    # Plasma concentration (ug/mL = mg/L). Doses in mg, central in mg,
    # vc in L, so central / vc has units mg/L = ug/mL with no conversion.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
