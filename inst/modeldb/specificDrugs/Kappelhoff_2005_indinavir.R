Kappelhoff_2005_indinavir <- function() {
  description <- paste(
    "One-compartment first-order-absorption popPK model for oral indinavir in",
    "HIV-1-infected adults, with multiplicative covariate effects of",
    "concomitant ritonavir (CL/F x 0.354) and concomitant NNRTI",
    "(efavirenz/nevirapine; CL/F x 1.41) on apparent clearance and of female",
    "sex on apparent bioavailability (F x 1.48). A 0.485 h absorption lag-time",
    "is applied only when ritonavir is co-administered (Kappelhoff 2005)."
  )
  reference <- paste(
    "Kappelhoff BS, Huitema ADR, Sankatsing SUC, Meenhorst PL, Van Gorp ECM,",
    "Mulder JW, Prins JM, Beijnen JH. Population pharmacokinetics of indinavir",
    "alone and in combination with ritonavir in HIV-1-infected patients.",
    "Br J Clin Pharmacol. 2005;60(3):276-86.",
    "doi:10.1111/j.1365-2125.2005.02436.x."
  )
  vignette <- "Kappelhoff_2005_indinavir"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    CONMED_RTV = list(
      description        = "Concomitant ritonavir co-administration indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant ritonavir)",
      notes              = paste(
        "Time-varying per occasion. Kappelhoff 2005 pooled ritonavir doses",
        "(100, 200, 400 mg) because the data supported maximal inhibition of",
        "indinavir CL at the lowest clinical dose; AUC50 was estimated to be",
        "very small in the AUC-based exploratory model. 288 of 443 occasions",
        "included ritonavir."
      ),
      source_name        = "RTV"
    ),
    CONMED_NNRTI = list(
      description        = "Concomitant NNRTI (efavirenz or nevirapine) co-administration indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant NNRTI)",
      notes              = paste(
        "Kappelhoff 2005 pooled efavirenz and nevirapine because including",
        "separate effects for each drug did not improve goodness-of-fit.",
        "35 of 443 occasions (7.9%) included an NNRTI."
      ),
      source_name        = "NNRTI"
    ),
    SEXF = list(
      description        = "Sex (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = paste(
        "Time-fixed per subject. Source paper uses SEX = 0 male, 1 female",
        "(Kappelhoff 2005 Table 2 footnote); same orientation as the canonical",
        "SEXF, so the column maps without inversion. 9 of 147 patients (6.1%)",
        "were female."
      ),
      source_name        = "SEX"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 147L,
    n_studies      = 1L,
    n_observations = 853L,
    age_median     = "40.3 years (IQR 34.9-47.1)",
    weight_median  = "73.0 kg (IQR 65.0-80.0)",
    sex_female_pct = 6.1,
    race_ethnicity = c(Caucasian = 82.3, Black = 9.5, Asian = 4.8, Latino = 3.4),
    disease_state  = paste(
      "HIV-1 infection on indinavir-containing antiretroviral therapy.",
      "Median CD4 cell count 380 / mm3 (IQR 220-575); median plasma",
      "log10 HIV-1 RNA 2.30 copies / mL (IQR 2.30-3.67). 5 of 138 (3.4%)",
      "with chronic hepatitis B and 8 of 138 (5.4%) with chronic hepatitis C."
    ),
    dose_range     = paste(
      "Indinavir 200-1400 mg per dose, given two or three times daily, with",
      "or without ritonavir 100-400 mg per dose twice daily. Most common",
      "regimens: 800 mg TID indinavir alone (112 occasions, 25.3%);",
      "800 mg BID indinavir + 100 mg BID ritonavir (201 occasions, 45.4%);",
      "400 mg BID indinavir + 400 mg BID ritonavir (37 occasions, 8.4%)."
    ),
    regions        = "Netherlands (Slotervaart Hospital and Academic Medical Centre, Amsterdam)",
    notes          = paste(
      "443 occasions across 147 patients. 45 patients contributed full PK",
      "profiles (8-12 timepoints per profile); the remainder contributed",
      "randomly timed therapeutic-drug-monitoring samples (2-3 samples /",
      "patient, range 1-18, follow-up 0-64 months). All samples at steady",
      "state, at least 2 weeks after initiation of the indinavir regimen.",
      "Plasma indinavir concentrations measured by LC-MS/MS over",
      "0.01-10 mg/L. Baseline demographics: Kappelhoff 2005 Table 1."
    )
  )

  ini({
    # Structural parameters - Kappelhoff 2005 Table 2 "Final model" column.
    lcl     <- log(46.8) ; label("Apparent clearance CL/F in the reference category (male, no RTV, no NNRTI) (L/h)")  # Table 2 Final CL/F = 46.8 L/h, RSE 5.75%
    lvc     <- log(82.3) ; label("Apparent volume of distribution V/F (L)")                                              # Table 2 Final V/F = 82.3 L, RSE 4.70%
    lka     <- log(2.62) ; label("First-order absorption rate constant ka (1/h)")                                        # Table 2 Final ka = 2.62 1/h, RSE 16.0%
    ltlag   <- log(0.485); label("Absorption lag-time when ritonavir is co-administered (h)")                            # Table 2 Final Lag-time = 0.485 h, RSE 1.79%; only estimated when indinavir and ritonavir were combined
    lfdepot <- fixed(log(1)); label("Apparent bioavailability F in the male reference category (anchor)")                # Table 2 equation block "F = 1 * 1.48^SEX" with male reference fixed at 1

    # Covariate effects - Kappelhoff 2005 Table 2 footnote equation:
    #   CL/F = 46.8 * 0.354^RTV * 1.41^NNRTI
    #   F    = 1 * 1.48^SEX
    e_conmed_rtv_cl    <- log(0.354); label("Effect of concomitant ritonavir on log(CL/F) (multiplicative 0.354)")     # Table 2 Final theta_ritonavir = 0.354, RSE 6.07%; -64.6% reduction in CL/F
    e_conmed_nnrti_cl  <- log(1.41) ; label("Effect of concomitant NNRTI on log(CL/F) (multiplicative 1.41)")          # Table 2 Final theta_concomitant_NNRTI = 1.41, RSE 4.78%; +41% increase in CL/F
    e_sexf_fdepot      <- log(1.48) ; label("Effect of female sex on log(F) (multiplicative 1.48)")                    # Table 2 Final theta_female = 1.48, RSE 16.7%; +48% increase in bioavailability

    # Inter-individual variability (omega^2 = log(CV^2 + 1) for the
    # exponential / log-normal error model used by Kappelhoff 2005).
    # Off-diagonal covariance = correlation * sqrt(var_CL * var_VC).
    # Final-model values: IIV CL/F = 24.2%, IIV V/F = 24.6%, correlation = 0.629.
    etalcl + etalvc ~ c(
      log(1 + 0.242^2),
      0.629 * sqrt(log(1 + 0.242^2) * log(1 + 0.246^2)),
      log(1 + 0.246^2)
    )  # Table 2 Final IIV CL/F 24.2%, IIV V/F 24.6%, correlation 0.629

    # Residual error (combined additive + proportional).
    addSd  <- 0.0491; label("Additive residual SD (mg/L)")           # Table 2 Final Additive error = 0.0491 mg/L, RSE 16.7%
    propSd <- 0.353 ; label("Proportional residual SD (fraction)")    # Table 2 Final Proportional error = 35.3%, RSE 6.18%
  })

  model({
    # Individual PK parameters. The exp(eta) factors multiply the typical
    # value; covariate effects enter as multiplicative powers per Table 2.
    cl     <- exp(lcl + e_conmed_rtv_cl * CONMED_RTV + e_conmed_nnrti_cl * CONMED_NNRTI + etalcl)
    vc     <- exp(lvc + etalvc)
    ka     <- exp(lka)
    fdepot <- exp(lfdepot + e_sexf_fdepot * SEXF)

    # Lag-time applies only when ritonavir is co-administered; setting
    # tlag = 0 in the RTV-free occasions reduces alag(depot) to 0.
    tlag   <- exp(ltlag) * CONMED_RTV

    # One-compartment first-order-absorption model.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - (cl / vc) * central

    f(depot)    <- fdepot
    alag(depot) <- tlag

    # Concentration: dose in mg, volume in L -> mg/L.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
