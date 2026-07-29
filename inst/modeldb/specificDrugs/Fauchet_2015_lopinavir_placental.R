Fauchet_2015_lopinavir_placental <- function() {
  description <- "One-compartment first-order-absorption population PK model for total lopinavir in HIV-infected pregnant and nonpregnant women with a maternal-to-fetal effect-compartment placental-transfer chain and a downstream fetal-to-amniotic-fluid distribution-and-elimination chain; a 39% pregnancy effect is applied multiplicatively to apparent maternal CL (Fauchet 2015 MFLA submodel)."
  reference <- "Fauchet F, Treluyer JM, Illamola SM, Pressiat C, Lui G, Valade E, Mandelbrot L, Lechedanec J, Delmas S, Blanche S, Warszawski J, Urien S, Tubiana R, Hirt D, for the ANRS 135 PRIMEVA Study Group. Population approach to analyze the pharmacokinetics of free and total lopinavir in HIV-infected pregnant women and consequences for dose adjustment. Antimicrob Agents Chemother. 2015;59(9):5727-5735. doi:10.1128/AAC.00863-15"
  vignette <- "Fauchet_2015_lopinavir"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  paper_specific_compartments <- c("fetal", "amniotic")

  covariateData <- list(
    PREG = list(
      description        = "Pregnancy indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = "1 = pregnant, 0 = nonpregnant. Applied multiplicatively to apparent maternal CL via beta_CL_ENC^PREG so that the published 1.39 multiplier (39% higher CL in pregnant women) is the verbatim source value (Fauchet 2015 Table 4 row 'beta_CL/ENC'; Results: 'Pregnant women were found to have a clearance that was 39% higher than that of nonpregnant women.').",
      source_name        = "PREG"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Tested in the covariate screen ('body weight, age, height, body mass index, gestational age, pregnancy, delivery, therapy group') but did not significantly improve the MFLA model after the pregnancy effect on CL was included (Fauchet 2015 Results, 'MFLA model' subsection: 'After inclusion of this factor in the model, no other covariate had a significant effect on the LPV pharmacokinetics').",
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Maternal age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Tested in the covariate screen but not retained after the pregnancy effect on CL was included (Fauchet 2015 Results, 'MFLA model' subsection).",
      source_name        = "AGE"
    ),
    GA = list(
      description        = "Gestational age at sampling",
      units              = "weeks",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Tested as a continuous covariate on CL via the paper's equation 3 (continuous form, applied only to pregnant subjects) but a binary PREG indicator gave a better fit (Fauchet 2015 Results, 'MFLA model' subsection). The canonical GA records gestational age at birth; the same column with time-of-sampling semantics is reused here because the orientation and units match.",
      source_name        = "GA"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 208L,
    n_studies      = 2L,
    n_observations = 527L,
    age_range      = "18.3-44.3 years",
    age_median     = "31.2 years",
    weight_range   = "45-122 kg",
    weight_median  = "76 kg",
    sex_female_pct = 100,
    race_ethnicity = "not reported in the source paper",
    disease_state  = "HIV-1 infection; pregnant, nonpregnant, and in-labour women treated for HIV infection or for prevention of mother-to-child transmission",
    ga_range       = "9-41 weeks (sampled gestational age across pregnant cohorts)",
    dose_range     = "LPV/r 400/100 mg twice daily oral (one subject 600 mg LPV BID); intravenous zidovudine during delivery for PRIMEVA-study women per French guidelines at the time",
    regions        = "France",
    notes          = "Pooled cohort across the ANRS 135 PRIMEVA randomised trial (n=103: 69 LPV/r monotherapy, 34 LPV/r + zidovudine/lamivudine triple therapy) and the Hospital Cochin (Paris) therapeutic drug monitoring cohort (n=105: 81 pregnant, 24 nonpregnant). MFLA submodel jointly fit 400 maternal + 79 cord blood + 48 amniotic fluid samples (paper Results); cohort demographics from Table 1."
  )

  ini({
    # Structural parameters from Fauchet 2015 Table 4 ('Population pharmacokinetic
    # parameters for lopinavir based on the MFLA model'). Apparent values absorb F;
    # CL and V are the maternal values for a nonpregnant reference subject.
    lka  <- fixed(log(0.255)); label("Absorption rate constant (1/h, fixed)")                  # Table 4 row 'Ka' = 0.255 /h (footnote c 'Fixed value'); Results 'MFLA model' subsection: 'Our data did not allow estimation of a Ka value. The stability of the model was improved for a Ka value fixed to 0.255'
    lcl  <- log(4.12);         label("Apparent nonpregnant maternal clearance (CL/F, L/h)")    # Table 4 row 'CL' = 4.12 L/h (nonpregnant reference; pregnancy applied via beta_CL_ENC)
    lvc  <- log(43.1);         label("Apparent maternal volume of distribution (V/F, L)")      # Table 4 row 'V' = 43.1 L

    # Placental-transfer rate constants (paper-mechanistic effect-compartment chain).
    # Fetal and amniotic states carry concentration (mg/L) so K_1F operates on Cc
    # rather than on the maternal amount, preserving the source's 'effect compartment
    # ... did not modify the compartmental model for the mother' formulation (Fauchet
    # 2015 Methods 'MFLA model' paragraph; Appendix differential system A(1)-A(4)).
    lk1f <- log(0.035);        label("Maternal-to-fetal transfer rate constant (1/h)")        # Table 4 row 'K_1F' = 0.035 /h
    lkfla<- log(0.289);        label("Fetal-to-amniotic-liquid transfer rate constant (1/h)") # Table 4 row 'K_FLA' = 0.289 /h
    lkla <- log(0.453);        label("Amniotic-liquid elimination rate constant (1/h)")       # Table 4 row 'K_LA' = 0.453 /h

    # Pregnancy effect on apparent maternal CL. Encoded as a power-of-factor multiplier
    # (CL_pregnant / CL_nonpregnant = e_preg_cl^PREG = 1.39 at PREG=1) so the published
    # 1.39 enters verbatim. Results: 'Pregnant women were found to have a clearance
    # that was 39% higher than that of nonpregnant women' (Table 4 row 'beta_CL/ENC').
    e_preg_cl <- 1.39; label("Multiplicative pregnancy effect on CL (CL_pregnant / CL_nonpregnant)")  # Table 4 row 'beta_CL/ENC' = 1.39

    # Inter-individual variability. Fauchet 2015 Table 4 footnote d notes 'omega, the
    # square root of the between-subject variance', so the variances entered here are
    # omega^2 of the reported SD values. IIV was estimated on CL, V, and K_1F per the
    # Results paragraph 'the between-subject variability was estimated for CL, V, and K_1F'.
    etalcl  ~ 0.02190  # Table 4 row 'omega_Cl'  = 0.148 -> variance = 0.148^2 = 0.02190
    etalvc  ~ 0.19625  # Table 4 row 'omega_V'   = 0.443 -> variance = 0.443^2 = 0.19625
    etalk1f ~ 0.78500  # Table 4 row 'omega_K1F' = 0.886 -> variance = 0.886^2 = 0.78500

    # Residual error -- separate proportional models for maternal, fetal, and amniotic
    # observations (Table 4 footnote f / Results 'MFLA model' subsection: 'A separate
    # proportional error was used for maternal, fetal, and amniotic fluid LPV
    # concentrations').
    propSd           <- 0.434; label("Proportional residual error on maternal LPV (fraction)")          # Table 4 row 'sigma_maternal'       = 0.434
    propSd_Cfetal    <- 0.455; label("Proportional residual error on fetal cord-blood LPV (fraction)")  # Table 4 row 'sigma_fetal'          = 0.455
    propSd_Camniotic <- 0.497; label("Proportional residual error on amniotic-fluid LPV (fraction)")    # Table 4 row 'sigma_amniotic liquid' = 0.497
  })

  model({
    # Individual maternal PK parameters.
    ka  <- exp(lka)
    cl  <- exp(lcl + etalcl) * e_preg_cl^PREG
    vc  <- exp(lvc + etalvc)
    kel <- cl / vc

    # Paper-mechanistic transfer rate constants.
    k1f <- exp(lk1f + etalk1f)
    kfla<- exp(lkfla)
    kla <- exp(lkla)

    # ODE system. Maternal central is in mg; fetal and amniotic states carry
    # concentration (mg/L) so K_1F * Cc has units (1/h) * (mg/L) = mg/L/h, consistent
    # with d(fetal)/dt in mg/L/h. The placental-transfer effect-compartment does NOT
    # subtract drug from the maternal compartment (per Fauchet 2015 Methods).
    d/dt(depot)    <- -ka * depot
    d/dt(central)  <-  ka * depot - kel * central
    Cc             <-  central / vc
    d/dt(fetal)    <-  k1f * Cc    - kfla * fetal
    d/dt(amniotic) <-  kfla * fetal - kla  * amniotic

    Cfetal    <- fetal
    Camniotic <- amniotic

    # Three-output residual error.
    Cc        ~ prop(propSd)
    Cfetal    ~ prop(propSd_Cfetal)
    Camniotic ~ prop(propSd_Camniotic)
  })
}
