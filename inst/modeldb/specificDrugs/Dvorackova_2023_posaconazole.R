Dvorackova_2023_posaconazole <- function() {
  description <- "Population PK model for oral posaconazole tablets in adult lung-transplant recipients (Dvorackova 2023). One-compartment disposition with first-order absorption and first-order elimination, parameterised on the apparent (oral) scale as CL/F and Vd/F because bioavailability was not identifiable from the oral-only therapeutic-drug-monitoring data. The absorption rate constant Ka was fixed to 0.8 1/h (back-calculated from the tmax and half-life reported in the posaconazole SmPC) because all concentrations were sampled in the elimination phase; inter-individual variability on Ka was nevertheless estimated and is very large, reflecting that absorption is essentially unidentifiable from these data. Age is the only covariate retained in the final model, entering log-linearly (exponentially) on apparent clearance so that CL/F declines by about 0.9 percent per year of age. Residual variability is proportional. The model was used for Monte Carlo dose optimisation against the EUCAST trough targets of 0.7 mg/L for prophylaxis and 1.25 mg/L for therapy."
  reference   <- "Dvorackova E, Sima M, Zajacova A, Vyskocilova K, Kotowski T, Dunovska K, Klapkova E, Havlin J, Lischke R, Slanar O. Dosing Optimization of Posaconazole in Lung-Transplant Recipients Based on Population Pharmacokinetic Model. Antibiotics (Basel). 2023;12(9):1399. doi:10.3390/antibiotics12091399."
  vignette    <- "Dvorackova_2023_posaconazole"
  units       <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    AGE = list(
      description        = "Subject age at the start of posaconazole therapy",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject (recorded at the beginning of therapy; Dvorackova 2023 Table 1).",
        "The only covariate retained in the final model (Dvorackova 2023 Table 3: p = 0.012 for age on CL/F;",
        "all other tested covariates had p > 0.05). Enters log-linearly and UNCENTERED on apparent clearance:",
        "CL/F = CL/F_pop * exp(e_age_cl * AGE), per the printed model equation",
        "log(CL/F) = log(CL/F_pop) + beta_CL/F_age * age + eta_CL/F (Dvorackova 2023 Section 2.2).",
        "Because age is uncentered, the typical value lcl = log(8.8) is the mathematical intercept at AGE = 0,",
        "not a value attained by any study subject; at the cohort median age of 56 years the typical CL/F is",
        "8.8 * exp(-0.009 * 56) = 5.32 L/h. Observed age range 22-71 years; do not extrapolate outside it.",
        "See the vignette Assumptions and deviations section for the arithmetic confirming the exponential",
        "(rather than additive) covariate form against the paper's own simulated AUC values.",
        sep = " "
      ),
      source_name        = "Age"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 32L,
    n_studies      = 1L,
    age_range      = "22-71 years (median 56, IQR 48-61; Dvorackova 2023 Table 1)",
    age_median     = "56 years",
    weight_range   = "38-100 kg (median 69, IQR 60-83; Dvorackova 2023 Table 1)",
    weight_median  = "69 kg",
    sex_female_pct = 31.3,
    race_ethnicity = "Not reported (single-centre Czech cohort)",
    disease_state  = paste(
      "Adult lung-transplant recipients receiving oral posaconazole delayed-release tablets (Noxafil)",
      "for antifungal prophylaxis or for treatment of invasive fungal disease during the postoperative period.",
      "Indications for lung transplantation were interstitial lung disease (13), idiopathic pulmonary fibrosis",
      "or other fibrotic changes (7), chronic obstructive pulmonary disease (7), cystic fibrosis (3),",
      "bronchial asthma (1), and chronic aspergillosis (1).",
      sep = " "
    ),
    dose_range     = "100-400 mg once daily oral posaconazole tablets; all patients started at 300 mg once daily and were subsequently dose-adjusted by therapeutic drug monitoring.",
    regions        = "Czech Republic (Motol University Hospital, Prague), October 2020 to March 2023",
    n_observations = "80 serum posaconazole concentrations (1-12 per patient, mean 2.5, mode 2), all drawn in the elimination phase from day 4 of treatment onward (Dvorackova 2023 Sections 2.1 and 4.2).",
    renal_function = "Serum creatinine median 110 umol/L (range 53-301); CKD-EPI eGFR median 0.98 mL/s/1.73 m2 (range 0.32-2.16) (Dvorackova 2023 Table 1).",
    hepatic_function = "ALT median 0.72 ukat/L (range 0.24-4.97), AST median 0.35 ukat/L (range 0.15-1.68), GGT median 0.63 ukat/L (range 0.20-3.30) (Dvorackova 2023 Table 1).",
    co_medication  = "All patients received induction with basiliximab or antithymocyte globulin plus maintenance immunosuppression: tacrolimus (31) or cyclosporin A (1), prednisone (32), and mycophenolate mofetil (31). Gastric-pH-raising drugs (proton pump inhibitors or famotidine) were taken by 29 of 32 patients.",
    notes          = paste(
      "Prospective, observational, single-centre therapeutic-drug-monitoring study",
      "(ANZCTR ACTRN12622000997752; Ethics Committee EK-873/22).",
      "Posaconazole was quantified by a validated LC-MS/MS assay with a lower limit of quantification of 0.1 mg/L",
      "and a validated range of 0.1-20 mg/L (Dvorackova 2023 Section 4.3).",
      "Model estimated in Monolix Suite 2021R2 by SAEM; Monte Carlo dose simulations in Simulx 2021R2.",
      "Baseline demographics are in Dvorackova 2023 Table 1.",
      sep = " "
    )
  )

  ini({
    # ====================================================================
    # Structural PK parameters (Dvorackova 2023 Table 2, 'Fixed effects').
    # One-compartment model on the apparent (oral) scale: because
    # posaconazole was given only orally, bioavailability F could not be
    # estimated and the reported disposition parameters are CL/F and
    # Vd/F (Dvorackova 2023 Section 4.4 step 1). No separate
    # bioavailability term is therefore encoded: nominal dose enters the
    # depot compartment and F is absorbed into lcl and lvc.
    # ====================================================================

    # Ka was NOT estimated. Dvorackova 2023 Section 4.4 step 1: 'Ka was
    # fixed to a value of 0.8 h-1, which was calculated from t1/2 and
    # time to reach maximum plasma concentration (tmax) values reported
    # in SmPC'. Table 2 accordingly reports R.S.E. = NA for Ka_pop.
    lka <- fixed(log(0.8)); label("First-order absorption rate constant Ka (1/h, from the posaconazole SmPC)")  # Dvorackova 2023 Table 2: Ka_pop = 0.8 1/h, R.S.E. NA (fixed)

    lvc <- log(386.35);     label("Apparent volume of distribution Vd/F (L)")  # Dvorackova 2023 Table 2: Vd/F_pop = 386.35 L, R.S.E. 48.1 percent

    # lcl is the intercept of the log-linear age model at AGE = 0, i.e.
    # the typical CL/F extrapolated to age zero, NOT a value attained by
    # any subject (observed age range 22-71 y). At the cohort median age
    # of 56 y the typical CL/F is 8.8 * exp(-0.009 * 56) = 5.32 L/h.
    lcl <- log(8.8);        label("Apparent clearance CL/F at AGE = 0 (L/h; log-linear age-model intercept)")  # Dvorackova 2023 Table 2: CL/F_pop = 8.8 L/h, R.S.E. 32.0 percent

    # ====================================================================
    # Covariate effect: age on apparent clearance.
    #
    # FORM. The printed final model equation (Dvorackova 2023 Section 2.2)
    # is
    #     log(CL/F) = log(CL/F_pop) + beta_CL/F_age * age + eta_CL/F
    # i.e. an EXPONENTIAL (log-linear) effect of UNCENTERED age in years:
    #     CL/F = 8.8 * exp(-0.009 * AGE)
    # The abstract and Table 2 row label instead describe the effect
    # additively ('CL/F of 8.8 L/h decreases by 0.009 L/h with each year
    # of age'; Table 2 units column 'L/h per year'). The two readings are
    # not equivalent and the printed equation is the one implemented
    # here, per the standing convention that a printed equation
    # supersedes prose. The equation is independently corroborated by the
    # paper's own Monte Carlo results (Dvorackova 2023 Section 2.3): the
    # reported median steady-state AUC values imply CL/F = dose / AUC =
    # 300/51 = 400/68 = 5.88 L/h in the under-60 subgroup and
    # 200/40 = 300/60 = 5.00 L/h in the over-60 subgroup. Solving
    # 8.8 * exp(-0.009 * age) for those clearances gives ages of 44.8 and
    # 62.8 years respectively -- plausible subgroup medians for a cohort
    # with overall median age 56 y (IQR 48-61). The additive reading
    # would instead require ages of 324 and 422 years, and would produce
    # only a 0.4 percent change in CL/F across the entire 22-71 y range,
    # far too small to motivate the paper's 200/300/400 mg age-banded
    # dosing proposal. See the vignette for this check in full.
    # ====================================================================
    e_age_cl <- -0.009; label("Log-linear effect of age on apparent clearance (1/year; CL/F ~ exp(e_age_cl * AGE))")  # Dvorackova 2023 Table 2: beta_CL/F_age = -0.009, R.S.E. 63.0 percent

    # ====================================================================
    # Inter-individual variability (Dvorackova 2023 Table 2, 'Standard
    # deviation of the random effects'). Monolix reports omega directly
    # as the standard deviation of the log-scale random effect for
    # log-normally distributed parameters (Dvorackova 2023 Section 4.4:
    # 'The model parameters were assumed to be log-normally distributed'),
    # so the nlmixr2 variances below are simply omega^2 -- no CV-percent
    # to variance conversion is needed.
    #
    # NOTE on etalka. Ka's TYPICAL VALUE was fixed, but its random effect
    # was estimated and appears both in Table 2 (omega_Ka = 3.43, R.S.E.
    # 27.5 percent) and in the printed equation
    # log(Ka) = log(Ka_pop) + eta_Ka. The magnitude is extreme -- an
    # omega of 3.43 on the log scale spans roughly four orders of
    # magnitude either side of 0.8 1/h -- and is best read as the model
    # reporting that absorption is unidentifiable from these data, all of
    # which were sampled in the elimination phase (Dvorackova 2023
    # Section 4.2 and the study limitations in Section 3). It is encoded
    # faithfully here rather than trimmed. Because this is a
    # one-compartment model with linear elimination, ka does not affect
    # steady-state AUC over a dosing interval (AUC_tau = dose / (CL/F)
    # regardless of ka), so the practical impact is confined to the shape
    # of the absorption phase and a modest widening of the trough
    # distribution. Users refitting this model to data that do capture
    # the absorption phase should expect to re-estimate both lka and
    # etalka.
    # ====================================================================
    etalka ~ 3.43^2  # Dvorackova 2023 Table 2: omega_Ka   = 3.43, R.S.E. 27.5 percent; variance = 11.7649
    etalvc ~ 0.45^2  # Dvorackova 2023 Table 2: omega_Vd/F = 0.45, R.S.E. 47.7 percent; variance = 0.2025
    etalcl ~ 0.36^2  # Dvorackova 2023 Table 2: omega_CL/F = 0.36, R.S.E. 20.3 percent; variance = 0.1296

    # ====================================================================
    # Residual error. Dvorackova 2023 Section 2.2: a proportional error
    # model was selected over additive and combined alternatives.
    # ====================================================================
    propSd <- 0.29; label("Proportional residual error (fraction)")  # Dvorackova 2023 Table 2: Proportional error model parameter = 0.29, R.S.E. 11.3 percent
  })

  model({
    # ====================================================================
    # 1. Individual PK parameters
    #
    # Age enters CL/F log-linearly and uncentered, matching the printed
    # equation log(CL/F) = log(CL/F_pop) + beta_CL/F_age * age + eta_CL/F.
    # ====================================================================
    ka <- exp(lka + etalka)
    vc <- exp(lvc + etalvc)
    cl <- exp(lcl + e_age_cl * AGE + etalcl)

    # ====================================================================
    # 2. Micro-constants
    # ====================================================================
    kel <- cl / vc

    # ====================================================================
    # 3. ODE system -- one compartment with first-order absorption.
    # Oral dose records target cmt = depot with the nominal dose in mg;
    # bioavailability is absorbed into the apparent parameters CL/F and
    # Vd/F, so no f(depot) term is applied.
    # ====================================================================
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # ====================================================================
    # 4. Observation and residual error.
    # Cc (mg/L) = central (mg) / vc (L).
    # ====================================================================
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
