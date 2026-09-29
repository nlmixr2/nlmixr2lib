Brekkan_2018_omalizumab <- function() {
  description <- "Reduced quasi-equilibrium target-mediated drug disposition (TMDD) model for omalizumab and IgE, used as the reference model for optimal-design evaluation of reduced trial designs (Brekkan 2018). Three serum entities (free omalizumab, free IgE, and the omalizumab-IgE complex) are coupled through instantaneous-equilibrium binding with a concentration-ratio-dependent dissociation constant. Subcutaneous absorption is first-order. IgE is produced at a constant zero-order synthesis rate, so the pretreatment total-IgE state starts at its own steady state (ksyn * v_ige / cl_ige) rather than from a baseline-IgE covariate. This is the covariate-free reduction of the Hayashi 2007 model: body-weight and baseline-IgE covariate relationships and the correlation between random effects were removed so that the design optimization would not require integration over covariate distributions. Three observed quantities: total omalizumab (ug/mL), total IgE (ng/mL), and free IgE (ng/mL)."
  reference <- "Brekkan A, Jonsson S, Karlsson MO, Hooker AC. Reduced and optimized trial designs for drugs described by a target mediated drug disposition model. J Pharmacokinet Pharmacodyn. 2018;45(4):637-647. doi:10.1007/s10928-018-9594-9 (PMCID PMC6061097). Model structure and parameter values are given in Supplementary material Appendix 1 (Electronic Supplementary Material 10928_2018_9594_MOESM1_ESM.docx), Table 1A. The parent model, from which the covariate relationships and random-effect correlations were removed, is Hayashi N, Tsukamoto Y, Sallas WM, Lowe PJ. Br J Clin Pharmacol. 2007;63(5):548-561; see modellib('Hayashi_2007_omalizumab')."
  vignette <- "Brekkan_2018_omalizumab"
  units <- list(
    time = "day",
    dosing = "mg",
    concentration = "ug/mL (total omalizumab); ng/mL (free and total IgE)"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Supplement Appendix 1 names SC_D as 'the amount of drug
  # in the subcutaneous (SC) dosing compartment', OMA_T as 'the nanomolar
  # amount of total OMA' and IGE_T as 'the nanomolar amount of total IgE',
  # so analyte and units are taken from the source, not inferred.
  compartmentData <- list(
    depot = list(analyte = "omalizumab", units = "nmol", specimen = "administration site", verified = TRUE),
    central = list(analyte = "omalizumab", units = "nmol", specimen = "serum", verified = TRUE),
    total_target = list(analyte = "IgE", units = "nmol", specimen = "serum", verified = TRUE)
  )

  # This reduction has no covariates by construction. The two covariates of the
  # parent Hayashi 2007 model were deliberately removed and are recorded here
  # (documentation only) so the provenance of the reduction is not lost.
  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Power covariate on omalizumab CL and on the omalizumab / IgE volume in the parent Hayashi 2007 model. Removed in Brekkan 2018 ('The model was simplified by removing covariate relationships (body weight and baseline IgE levels) and correlation between parameters', Methods 'Population model'). Retained in modellib('Hayashi_2007_omalizumab').",
      source_name = "body weight"
    ),
    IGE = list(
      description = "Baseline serum total IgE concentration (pretreatment)",
      units = "ng/mL",
      type = "continuous",
      notes = "Power covariate on IgE CL and IgE production rate, and the source of the total-IgE initial condition, in the parent Hayashi 2007 model. Removed in Brekkan 2018; the pretreatment total-IgE state is instead the steady state of the zero-order synthesis / first-order loss balance, ksyn * v_ige / cl_ige. Retained in modellib('Hayashi_2007_omalizumab').",
      source_name = "baseline IgE level"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 48L,
    n_studies = 1L,
    n_observations = 1872L,
    disease_state = "Atopic patients receiving subcutaneous omalizumab; the reference design and the parameter values are those of the Hayashi 2007 Japanese atopic-asthma analysis, which pooled a single-dose study in healthy atopic Japanese volunteers and a multiple-dose study in Japanese seasonal allergic rhinitis patients",
    dose_range = "Single subcutaneous doses of 75, 150, 300 and 375 mg (four dose groups; Brekkan 2018 Table 1 footnote a)",
    regions = "Japan (underlying Hayashi 2007 analysis)",
    race_ethnicity = c(White = 0, Black = 0, Asian = 100, Other = 0),
    sex_female_pct = NA_real_,
    age_range = "Adults; Brekkan 2018 reports no demographic table, the design is specified only by dose group and sampling schedule",
    design_sampling_times = "0, 0.5, 1, 2, 4, 7, 10, 14, 28, 42, 56, 70, 84 days (Brekkan 2018 Table 1, reference design 1)",
    design_analytes = "Total omalizumab, total IgE and free IgE measured at each of the 13 sampling times, giving 39 observations per individual",
    notes = "n_subjects is not printed directly: Brekkan 2018 Table 1 reports 1872 total samples and 39 observations per individual for the reference design, giving 1872 / 39 = 48 individuals across the four dose groups. The 936 total samples of design 5 (two dose groups removed) confirm 12 individuals per dose group. This model was used to simulate and optimize trial designs, not refitted to new data; the parameter values are those of the parent Hayashi 2007 analysis (202 subjects, 2 studies) rounded as reported in Supplementary Appendix 1 Table 1A."
  )

  ini({
    # -----------------------------------------------------------------------
    # Structural parameters - Supplementary material Appendix 1, Table 1A.
    # The source table reports rates in 1/h and mL/h and volumes in mL; both
    # are converted to the nlmixr2lib day / L convention here
    # (x_per_day = x_per_h * 24; V_L = V_mL / 1000). The reduction carries no
    # covariates, so these are the typical values for every individual.
    # -----------------------------------------------------------------------
    lka <- log(0.02 * 24); label("Subcutaneous absorption rate constant for omalizumab (1/d; source: 0.02 1/h)") # Supplement Table 1A, ka
    lcl <- log(0.00732 * 24); label("Clearance of free omalizumab (L/d; source: 7.32 mL/h)") # Supplement Table 1A, CLOMA
    lvc <- log(5.9); label("Central volume of omalizumab (L; source: 5900 mL)") # Supplement Table 1A, VOMA
    lcl_ige <- log(0.071 * 24); label("Clearance of free IgE (L/d; source: 71 mL/h)") # Supplement Table 1A, CLIGE
    lv_ige <- log(5.9); label("Central volume of IgE (L; source: 5900 mL)") # Supplement Table 1A, VIGE
    lcl_complex <- log(0.00586 * 24); label("Clearance of the omalizumab-IgE complex (L/d; source: 5.86 mL/h)") # Supplement Table 1A, CLComp
    lvc_complex <- log(3.63); label("Central volume of the omalizumab-IgE complex (L; source: 3630 mL)") # Supplement Table 1A, VCOMP
    lksyn <- log(0.158 * 24); label("Zero-order IgE synthesis rate (nmol/d; source: 0.158 nmol/h)") # Supplement Table 1A, ksyn
    lkd0 <- log(1.07); label("Equilibrium dissociation constant when total omalizumab equals total IgE (nM; source: 0.00107 nmol/mL)") # Supplement Table 1A, Kd0
    alpha <- 0.157; label("Exponent scaling Kd by the ratio of total omalizumab to total IgE (unitless)") # Supplement Table 1A, alpha

    # -----------------------------------------------------------------------
    # Inter-individual variability. Supplement Appendix 1 footnote b:
    # 'Interindividual variability reported as coefficient of variation', and
    # the narrative states IIV 'is described according to log-normal
    # distributions'. For a log-normal parameter omega^2 = log(CV^2 + 1):
    #   CLOMA  CV 0.20 -> log(1 + 0.20^2) = 0.0392207
    #   VOMA   CV 0.13 -> log(1 + 0.13^2) = 0.0167588
    #   ka     CV 0.40 -> log(1 + 0.40^2) = 0.1484200
    #   CLIGE  CV 0.25 -> log(1 + 0.25^2) = 0.0606246
    #   ksyn   CV 0.23 -> log(1 + 0.23^2) = 0.0515483
    #   CLcomp CV 0.35 -> log(1 + 0.35^2) = 0.1155583
    #   Vcomp  CV 0.25 -> log(1 + 0.25^2) = 0.0606246
    # Table 1A also lists 'omega Kd0 = 0*', i.e. IIV on Kd0 fixed to zero, and
    # lists no IIV row at all for VIGE. Both are therefore carried without an
    # eta. They are omitted rather than written as `~ fixed(0)` because a
    # zero-variance diagonal makes OMEGA singular and breaks the Cholesky
    # sampler used by rxSolve. No off-diagonal terms: the correlation between
    # random effects present in the parent model was removed in this reduction
    # (Methods, 'Population model').
    # -----------------------------------------------------------------------
    etalka ~ 0.1484200 # Supplement Table 1A, 'omega Ka' CV 0.40
    etalcl ~ 0.0392207 # Supplement Table 1A, 'omega CLOMA' CV 0.20
    etalvc ~ 0.0167588 # Supplement Table 1A, 'omega VOMA' CV 0.13
    etalcl_ige ~ 0.0606246 # Supplement Table 1A, 'omega CLIGE' CV 0.25
    etalksyn ~ 0.0515483 # Supplement Table 1A, 'omega ksyn' CV 0.23
    etalcl_complex ~ 0.1155583 # Supplement Table 1A, 'omega CLcomp' CV 0.35
    etalvc_complex ~ 0.0606246 # Supplement Table 1A, 'omega Vcomp' CV 0.25

    # -----------------------------------------------------------------------
    # Residual error. Supplement Appendix 1 narrative: residual variability is
    # described 'with an additive error component in the log domain', and
    # footnote c reports the values as coefficients of variation. An additive
    # error on the log scale is proportional error in the linear scale, so
    # each value is encoded directly as the proportional-error SD via prop().
    # -----------------------------------------------------------------------
    propSd <- 0.17; label("Proportional residual error on total omalizumab (fraction)") # Supplement Table 1A, sigma OMAT
    propSd_totalIgE <- 0.21; label("Proportional residual error on total IgE (fraction)") # Supplement Table 1A, sigma IGET
    propSd_freeIgE <- 0.22; label("Proportional residual error on free IgE (fraction)") # Supplement Table 1A, sigma IGEF
  })

  model({
    # ------------------------------------------------------------------
    # Molecular weights. The supplement writes the ODE system in nmol and
    # reports the go/no-go thresholds in ng/mL, but does not tabulate the
    # molecular weights needed to move between the two. Both are taken from
    # the parent Hayashi 2007 analysis (Methods, page 552), and the IgE value
    # is confirmed by Brekkan 2018's own printed baseline: the steady-state
    # free-IgE concentration ksyn / cl_ige = 0.158 / 0.071 = 2.2254 nM
    # multiplied by 190 gives 422.82 ng/mL, exactly the baseline quoted in
    # Methods, 'Go/no-go decision'.
    # 1 nM * MW_kDa = MW_kDa ng/mL, and 1 mg = 1000 / MW_kDa nmol.
    # ------------------------------------------------------------------
    MWX <- 150 # omalizumab molecular weight (kDa = ng/nmol)
    MWE <- 190 # IgE molecular weight (kDa = ng/nmol)

    # ------------------------------------------------------------------
    # 1. Individual parameters. No covariate terms in this reduction.
    # ------------------------------------------------------------------
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    cl_ige <- exp(lcl_ige + etalcl_ige)
    v_ige <- exp(lv_ige)
    cl_complex <- exp(lcl_complex + etalcl_complex)
    vc_complex <- exp(lvc_complex + etalvc_complex)
    ksyn <- exp(lksyn + etalksyn)
    kd0 <- exp(lkd0)

    # ------------------------------------------------------------------
    # 2. Concentration-ratio-dependent dissociation constant
    #    (Supplement Appendix 1, equation 4):
    #      Kd = Kd0 * (OMA_T / IGE_T)^alpha
    #    alpha 'accounts for different complexes being formed at different
    #    concentrations of IgE and OMA'; when the ratio is 1, Kd = Kd0.
    #    At t = 0 there is no drug, so the ratio and hence Kd are 0 and no
    #    complex forms.
    # ------------------------------------------------------------------
    kd <- kd0 * (central / total_target)^alpha

    # ------------------------------------------------------------------
    # 3. Equilibrium-binding solution for the complex amount
    #    (Supplement Appendix 1, equation 5):
    #      COMP = 0.5 * (S - sqrt(S^2 - 4 * OMA_T * IGE_T))
    #    with S = Kd * V_OMA * V_IGE / V_COMP + OMA_T + IGE_T. Dimensions:
    #    Kd [nmol/L] * L * L / L = nmol, matching the two amounts.
    # ------------------------------------------------------------------
    S <- kd * vc * v_ige / vc_complex + central + total_target
    COMP <- 0.5 * (S - sqrt(S * S - 4 * central * total_target))

    # Free and complex concentrations in nM (Supplement equations 6-8).
    C_COMP <- COMP / vc_complex
    C_OMA_F <- (central - COMP) / vc
    C_IGE_F <- (total_target - COMP) / v_ige

    # ------------------------------------------------------------------
    # 4. ODE system (Supplement Appendix 1, equations 1-3).
    #    State units: depot, central and total_target in nmol.
    #    The depot-to-central transfer converts the mg dose to nmol via
    #    1000 / MWX; the supplement writes its equations with the dose
    #    already in nmol.
    # ------------------------------------------------------------------
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot * (1000 / MWX) - cl * C_OMA_F - cl_complex * C_COMP
    d/dt(total_target) <- ksyn - cl_ige * C_IGE_F - cl_complex * C_COMP

    # Pretreatment total IgE is the steady state of the zero-order synthesis
    # and first-order loss balance, since no drug is present at t = 0:
    #   ksyn = cl_ige * total_target(0) / v_ige.
    total_target(0) <- ksyn * v_ige / cl_ige

    # ------------------------------------------------------------------
    # 5. Observation outputs in assay units (Supplement equations 9-10).
    #    Total omalizumab in ug/mL = (C_OMA_F + C_COMP) [nM] * MWX / 1000.
    #    Total IgE in ng/mL        = (C_IGE_F + C_COMP) [nM] * MWE.
    #    Free IgE in ng/mL         = C_IGE_F [nM] * MWE.
    # ------------------------------------------------------------------
    Cc <- (C_OMA_F + C_COMP) * MWX / 1000
    totalIgE <- (C_IGE_F + C_COMP) * MWE
    freeIgE <- C_IGE_F * MWE

    Cc ~ prop(propSd)
    totalIgE ~ prop(propSd_totalIgE)
    freeIgE ~ prop(propSd_freeIgE)
  })
}
