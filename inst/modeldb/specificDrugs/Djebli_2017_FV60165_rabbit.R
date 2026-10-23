Djebli_2017_FV60165_rabbit <- function() {
  description <- "Preclinical (rabbit). Four-compartment population PK model of the src tyrosine kinase inhibitor FV-60165 in the anterior segment of the rabbit eye after topical instillation of its lipophilic prodrug FV-80228 as a 0.6% solution or nanosuspension. States: tear film (dose site, relative bioavailability F1), naso-lacrimal duct, cornea and aqueous humour, each with its own volume and observed concentration. Drug leaves the tears to the naso-lacrimal duct (formulation-dependent first-order rate) and to the cornea (first-order rate); cornea and aqueous humour exchange by an inter-compartmental clearance, and drug is cleared from the naso-lacrimal duct, the aqueous humour, and the cornea (cornea clearance fixed at 10% of aqueous-humour clearance). Each eye is one individual. Log-normal IIV on all ten estimated parameters (the paper's proportional IIV translated at equal CV); the paper's single proportional residual error (shared by the four matrices) is carried as four equal per-matrix parameters."
  reference <- "Djebli N, Khier S, Griguer F, Coutant AL, Tavernier A, Fabre G, Leriche C, Fabre D. Ocular Drug Distribution After Topical Administration: Population Pharmacokinetic Model in Rabbits. Eur J Drug Metab Pharmacokinet. 2017;42(1):59-68. doi:10.1007/s13318-016-0319-4"
  vignette <- "Djebli_2017_FV60165_rabbit"
  units <- list(time = "h", dosing = "ug", concentration = "ng/uL")

  covariateData <- list(
    FORM_FV60165_NANOSUSP = list(
      description = "Formulation indicator: 1 = prodrug instilled as the 0.6% nanosuspension, 0 = instilled as the 0.6% solution",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (0.6% solution in 7% 2-hydroxypropyl-beta-cyclodextrin, pH 5)",
      notes = "Djebli 2017 Results: k12 = theta_i * FORMULATION + theta_j * (1 - FORMULATION); Table 1 gives k12 = 0.259 1/h for the nanosuspension and 0.299 1/h for the solution. The paper does not state which level is coded 1; the two typical values are carried as separate parameters so the coding only selects between them.",
      source_name = "FORMULATION"
    )
  )

  compartmentData <- list(
    tear = list(analyte = "FV-60165", units = "ug", specimen = "tear fluid", verified = TRUE),
    nasolacrimal_duct = list(analyte = "FV-60165", units = "ug", specimen = "tissue", verified = TRUE),
    cornea = list(analyte = "FV-60165", units = "ug", specimen = "tissue", verified = TRUE),
    aqueous_humour = list(analyte = "FV-60165", units = "ug", specimen = "aqueous humour", verified = TRUE)
  )

  population <- list(
    species = "rabbit (HY79b pigmented)",
    n_subjects = 84L,
    n_studies = 1L,
    weight_range = "approximately 2-2.5 kg",
    sex_female_pct = 0,
    disease_state = "healthy",
    dose_range = "Two 30 uL drops of 0.6% prodrug FV-80228 (180 ug per drop) into each eye, 10 min apart (360 ug per eye); one third of each drop (60 ug) enters the model as the dose into the tears",
    regions = "France",
    notes = "42 male HY79b pigmented rabbits (Hypharm-Roussay); 21 received the solution and 21 the nanosuspension. Each eye was treated as an individual (84 eyes). Destructive sampling: n = 3 rabbits per time per formulation at 0.083, 0.25, 0.5, 1, 2, 4 and 8 h, one time point per eye; 336 drug concentrations in tear fluid, naso-lacrimal duct, cornea and aqueous humour (Djebli 2017 Sects. 2.2-2.3 and 3). Prodrug concentrations were mostly below the LLOQ and were not modelled."
  )

  ini({
    lfdepot <- log(0.300); label("Relative bioavailability into the tears F1 (fraction of the dosed amount)") # Table 1 row 'F1' = 0.300
    lk_tear_nasolacrimal_duct_soln <- log(0.299); label("Tear-to-naso-lacrimal-duct rate constant k12, solution (1/h)") # Table 1 row 'k12 (h-1) solution' = 0.299
    lk_tear_nasolacrimal_duct_nanosusp <- log(0.259); label("Tear-to-naso-lacrimal-duct rate constant k12, nanosuspension (1/h)") # Table 1 row 'k12 (h-1) nanosuspensions' = 0.259
    lk_tear_cornea <- log(0.340); label("Tear-to-cornea rate constant k13 (1/h)") # Table 1 row 'k13 (h-1)' = 0.340
    lq_cornea_aqueous_humour <- log(0.136); label("Cornea/aqueous-humour inter-compartmental clearance Q (uL/h)") # Table 1 row 'Q (uL/h)' = 0.136
    lcl_nasolacrimal_duct <- log(275); label("Clearance from the naso-lacrimal duct Clnl (uL/h)") # Table 1 row 'Clnl (uL/h)' = 275
    lcl_aqueous_humour <- log(231); label("Clearance from the aqueous humour Clah (uL/h)") # Table 1 row 'Clah (uL/h)' = 231
    f_cl_cornea <- fixed(0.1); label("Cornea clearance as a fraction of aqueous-humour clearance, Clc = 0.1 * Clah (unitless)") # Results and Fig. 1 legend: Clc fixed to 10% of Clah after a 25-run sensitivity analysis
    lv_tear <- log(155); label("Tear volume of distribution V1 (uL)") # Table 1 row 'V1 (uL)' = 155
    lv_nasolacrimal_duct <- log(1.24); label("Naso-lacrimal-duct volume of distribution V2 (uL)") # Table 1 row 'V2 (uL)' = 1.24
    lv_cornea <- log(42.8); label("Cornea volume of distribution V3 (uL)") # Table 1 row 'V3 (uL)' = 42.8
    lv_aqueous_humour <- log(218); label("Aqueous-humour volume of distribution V4 (uL)") # Table 1 row 'V4 (uL)' = 218

    # IIV: Table 1 prints CV% = 100 * omega of the paper's proportional form
    # P = TV * (1 + eta); the printed 95% CIs reproduce as
    # CV * sqrt(1 +/- 1.96 * RSE) with the RSE on the variance, which fixes
    # that scale. Translated to log-normal at equal CV: omega^2 = log(CV^2 + 1).
    etalfdepot ~ 0.31189 # Table 1 IIV 'F1' CV 60.5%; log(1 + 0.605^2)
    etalk_tear_nasolacrimal_duct ~ 0.49470 # Table 1 IIV 'k12' CV 80.0%, shared by both formulations; log(1 + 0.800^2)
    etalk_tear_cornea ~ 0.18669 # Table 1 IIV 'k13' CV 45.3%; log(1 + 0.453^2)
    etalq_cornea_aqueous_humour ~ 1.37597 # Table 1 IIV 'Q' CV 172%; log(1 + 1.72^2)
    etalcl_nasolacrimal_duct ~ 0.59832 # Table 1 IIV 'Clnl' CV 90.5%; log(1 + 0.905^2)
    etalcl_aqueous_humour ~ 0.61316 # Table 1 IIV 'Clah' CV 92.0%; log(1 + 0.920^2)
    etalv_tear ~ 1.20599 # Table 1 IIV 'V1' CV 153%; log(1 + 1.53^2)
    etalv_nasolacrimal_duct ~ 3.00093 # Table 1 IIV 'V2' CV 437%; log(1 + 4.37^2)
    etalv_cornea ~ 0.23201 # Table 1 IIV 'V3' CV 51.1%; log(1 + 0.511^2)
    etalv_aqueous_humour ~ 2.08164 # Table 1 IIV 'V4' CV 265%; log(1 + 2.65^2)

    # The paper estimated ONE proportional sigma for all four matrices; nlmixr2
    # needs one parameter per endpoint, so the same value is carried four times.
    propSd_Ctear <- 0.00491; label("Proportional residual error, tear fluid (fraction)") # Table 1 'sigma' = 0.491% (Results: 'the estimated proportional residual error was only about 0.491 %')
    propSd_Cnasolacrimal_duct <- 0.00491; label("Proportional residual error, naso-lacrimal duct (fraction)") # Table 1 'sigma' = 0.491%, shared across matrices
    propSd_Ccornea <- 0.00491; label("Proportional residual error, cornea (fraction)") # Table 1 'sigma' = 0.491%, shared across matrices
    propSd_Caqueous_humour <- 0.00491; label("Proportional residual error, aqueous humour (fraction)") # Table 1 'sigma' = 0.491%, shared across matrices
  })

  model({
    # Formulation selects between the two estimated k12 typical values
    # (Djebli 2017 Results: k12 = theta_i * FORMULATION + theta_j * (1 - FORMULATION)).
    lk12_form <- FORM_FV60165_NANOSUSP *
      lk_tear_nasolacrimal_duct_nanosusp +
      (1 - FORM_FV60165_NANOSUSP) * lk_tear_nasolacrimal_duct_soln

    fdepot <- exp(lfdepot + etalfdepot)
    k_tear_nasolacrimal_duct <- exp(lk12_form + etalk_tear_nasolacrimal_duct)
    k_tear_cornea <- exp(lk_tear_cornea + etalk_tear_cornea)
    q_cornea_aqueous_humour <- exp(lq_cornea_aqueous_humour + etalq_cornea_aqueous_humour)
    cl_nasolacrimal_duct <- exp(lcl_nasolacrimal_duct + etalcl_nasolacrimal_duct)
    cl_aqueous_humour <- exp(lcl_aqueous_humour + etalcl_aqueous_humour)
    cl_cornea <- f_cl_cornea * cl_aqueous_humour
    v_tear <- exp(lv_tear + etalv_tear)
    v_nasolacrimal_duct <- exp(lv_nasolacrimal_duct + etalv_nasolacrimal_duct)
    v_cornea <- exp(lv_cornea + etalv_cornea)
    v_aqueous_humour <- exp(lv_aqueous_humour + etalv_aqueous_humour)

    # Fig. 1 structure. Footnotes to Table 1: Q = k34 * V3 = k43 * V4,
    # Clnl = k20 * V2, Clah = k40 * V4.
    d/dt(tear) <- -(k_tear_nasolacrimal_duct + k_tear_cornea) * tear
    d/dt(nasolacrimal_duct) <- k_tear_nasolacrimal_duct * tear - cl_nasolacrimal_duct / v_nasolacrimal_duct * nasolacrimal_duct
    d/dt(cornea) <- k_tear_cornea * tear - (cl_cornea + q_cornea_aqueous_humour) / v_cornea * cornea + q_cornea_aqueous_humour / v_aqueous_humour * aqueous_humour
    d/dt(aqueous_humour) <- q_cornea_aqueous_humour / v_cornea * cornea - (q_cornea_aqueous_humour + cl_aqueous_humour) / v_aqueous_humour * aqueous_humour

    f(tear) <- fdepot

    # Amounts in ug and volumes in uL give ug/uL; x 1000 converts to the
    # paper's ng/uL.
    Ctear <- 1000 * tear / v_tear
    Cnasolacrimal_duct <- 1000 * nasolacrimal_duct / v_nasolacrimal_duct
    Ccornea <- 1000 * cornea / v_cornea
    Caqueous_humour <- 1000 * aqueous_humour / v_aqueous_humour

    Ctear ~ prop(propSd_Ctear)
    Cnasolacrimal_duct ~ prop(propSd_Cnasolacrimal_duct)
    Ccornea ~ prop(propSd_Ccornea)
    Caqueous_humour ~ prop(propSd_Caqueous_humour)
  })
}
