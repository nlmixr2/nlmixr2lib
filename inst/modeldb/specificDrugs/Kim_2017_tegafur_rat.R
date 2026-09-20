Kim_2017_tegafur_rat <- function() {
  description <- paste(
    "Preclinical (rat). Population PK model for oral tegafur (given as the",
    "fluoropyrimidine combination S-1) and its active metabolite 5-FU in",
    "male Sprague-Dawley rats, with a herb-drug interaction arm for the",
    "traditional Korean polyherbal medicine Sipjeondaebo-tang (SDT, 1200",
    "mg/kg/day orally for seven consecutive days before the S-1 dose; Kim",
    "2017). Tegafur is dosed into an absorption-site compartment that",
    "drains by two competing first-order routes: intact tegafur into a",
    "two-compartment tegafur disposition model (Ka), and pre-systemic",
    "first-pass metabolism in gut and liver into an amount-only pool of",
    "the 5-FU precursor 5'-hydroxytegafur (Ka,Met). The precursor pool is",
    "additionally fed systemically by the fraction FMet of tegafur",
    "clearance and converts to 5-FU at KConv; 5-FU then follows its own",
    "two-compartment disposition. SDT pretreatment is carried as the",
    "binary covariate CONMED_SDT, which selects between separately",
    "estimated control and SDT-pretreated values of Ka, Ka,Met and the",
    "5-FU clearance; all other parameters are shared across the two arms.",
    "Repeated SDT dosing slowed tegafur absorption and raised 5-FU",
    "clearance 1.68-fold, reducing 5-FU exposure. Residual error is not",
    "reported in the source and is carried as fixed(0) -- see vignette",
    "Errata."
  )
  reference <- "Kim TH, Shin S, Shin JC, Bulitta JB, Weon KY, Yoo SD, Park GY, Jeong SW, Kwon DR, Min BS, Woo MH, Shin BS. Effect of Sipjeondaebo-Tang on the Pharmacokinetics of S-1, an Anticancer Agent, in Rats Evaluated by Population Pharmacokinetic Modeling. Molecules. 2017;22(9):1488. doi:10.3390/molecules22091488"
  vignette <- "Kim_2017_tegafur_rat"
  units <- list(
    time = "h",
    dosing = "ug/kg",
    concentration = "ng/mL"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Kim 2017 Section 4.5 equations 1-6
  # and Figure 2 (structural model diagram).
  compartmentData <- list(
    depot = list(analyte = "tegafur", units = "ug/kg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "tegafur", units = "ug/kg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "tegafur", units = "ug/kg", specimen = "tissue", verified = TRUE),
    precursor_5fu = list(
      analyte = "5'-hydroxytegafur (5-FU precursor)",
      units = "ug/kg",
      specimen = "not applicable",
      verified = TRUE
    ),
    central_5fu = list(analyte = "5-FU", units = "ug/kg", specimen = "plasma", verified = TRUE),
    peripheral1_5fu = list(analyte = "5-FU", units = "ug/kg", specimen = "tissue", verified = TRUE)
  )

  # XPre,5FU is an amount-only kinetic intermediate with no volume and no
  # measured concentration (Kim 2017 equation 4 acts on the amount). It is
  # whitelisted here pending operator ratification of `precursor_5fu` as a
  # canonical compartment -- see vignette Errata.
  paper_specific_compartments <- c("precursor_5fu")

  covariateData <- list(
    CONMED_SDT = list(
      description = "Sipjeondaebo-tang (SDT) pretreatment status: 1 = rat received SDT 1200 mg/kg orally once daily for seven consecutive days immediately before the S-1 dose; 0 = rat received the 1% CMC-Na vehicle on the same schedule.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (1% CMC-Na vehicle control)",
      notes = "Time-fixed arm indicator in the source parallel-group design. Selects between the separately estimated control and SDT-pretreated values of Ka, Ka,Met and CL5FU/F (Kim 2017 Table 2); all other structural parameters are shared. The modelled SDT arm is the seven-day REPEATED-dose pretreatment group -- a single SDT dose produced no significant pharmacokinetic change (Kim 2017 Section 2.1) and was not carried as a separate stratum.",
      source_name = "Con / Pretre (Kim 2017 Table 2 parameter subscripts)"
    )
  )

  population <- list(
    species = "rat (Sprague-Dawley)",
    n_subjects = 10L,
    n_studies = 2L,
    sex = "male",
    age_range = "8 weeks",
    weight_range = "280-300 g",
    disease_state = "Healthy male Sprague-Dawley rats; no tumour model.",
    dose_range = "Single oral gavage of S-1 providing tegafur 5 mg/kg, gimeracil 1.45 mg/kg and oteracil potassium 4.9 mg/kg in 7.5% DMSO, given 1 min after the final SDT (1200 mg/kg) or 1% CMC-Na dose. The co-fitted historical dataset used 5-FU 10 mg/kg intravenously.",
    regions = "Korea (Catholic University of Daegu; IACUC-2013-025).",
    n_observations = "Plasma sampled predose and at 0.25, 0.5, 1, 1.5, 2, 3, 4, 8, 12 and 24 h postdose from the jugular vein; tegafur, 5-FU and gimeracil assayed by LC/MS/MS (LLOQ 50, 10 and 50 ng/mL respectively).",
    notes = "The population fitted is the seven-day repeated-dose experiment: control (1% CMC-Na, n = 5) and SDT-pretreated (n = 5) rats, plus the plasma 5-FU profile after intravenous 5-FU 10 mg/kg from a previous study by the same group (Kim 2017 Section 4.5). Kim 2017 does not state whether the single-dose pretreatment groups (a further n = 5 per arm) also entered the fit; the covariate stratification and the reported SDT effect both track the repeated-dose arm, so only that arm is represented here -- see vignette Errata. Gimeracil was assayed and analysed noncompartmentally but was NOT part of the population model and is not represented in this file. Modelling software: importance-sampling MC-PEM in parallelised S-ADAPT 1.57 with SADAPT-TRAN."
  )

  ini({
    # ---- Tegafur absorption and competing pre-systemic first-pass ----
    # Kim 2017 equation 1: dXgut/dt = -(Ka + Ka,Met) * Xgut. The absorption
    # site drains by two parallel first-order routes, so Ka is the fraction
    # reaching tegafur central and Ka,Met the fraction converted to the
    # 5-FU precursor during gut and hepatic first pass. Both are estimated
    # separately in the control and SDT-pretreated arms (Kim 2017 Section
    # 2.3 and Table 2); the stratum suffixes follow the symmetric
    # stratum-suffix convention, so neither arm keeps the bare canonical.
    lka_ctl <- log(0.296)
    label("Tegafur absorption rate constant, vehicle control arm (Ka, 1/h)") # Kim 2017 Table 2, row 'Absorption rate constant for tegafur in control group': Ka,Con = 0.296 1/h
    lka_sdt <- log(0.197)
    label("Tegafur absorption rate constant, SDT-pretreated arm (Ka, 1/h)") # Kim 2017 Table 2, row 'Absorption rate constant for tegafur in SDT pretreatment group': Ka,Pretre = 0.197 1/h

    lk_precursor_5fu_form_ctl <- log(0.122)
    label("Pre-systemic formation rate constant of the 5-FU precursor from the absorption site, vehicle control arm (Ka,Met, 1/h)") # Kim 2017 Table 2, row 'Formation rate constant of 5-FU precursor from gut compartment in control group': Ka,Met,Con = 0.122 1/h
    lk_precursor_5fu_form_sdt <- log(0.0595)
    label("Pre-systemic formation rate constant of the 5-FU precursor from the absorption site, SDT-pretreated arm (Ka,Met, 1/h)") # Kim 2017 Table 2, row 'Formation rate constant of 5-FU precursor from gut compartment in SDT pretreatment group': Ka,Met,Pretre = 0.0595 1/h

    # ---- Precursor to 5-FU conversion ----
    # Kim 2017 equation 5: KConv * XPre,5FU. Shared across both arms.
    lk_5fu_form <- log(2.88)
    label("Conversion rate constant of the 5-FU precursor to 5-FU (KConv, 1/h)") # Kim 2017 Table 2, row 'Formation rate constant of 5-FU from 5-FU precursor': KConv = 2.88 1/h

    # ---- Tegafur disposition (apparent, /F) ----
    lcl <- log(0.0813)
    label("Tegafur apparent clearance (CLTeg/F, L/h/kg)") # Kim 2017 Table 2, row 'Clearance for tegafur': CLTeg/F = 0.0813 L/h/kg
    lq <- log(0.184)
    label("Tegafur apparent distribution clearance (CLdTeg/F, L/h/kg)") # Kim 2017 Table 2, row 'Distribution clearance for tegafur': CLdTeg/F = 0.184 L/h/kg
    lvc <- log(0.0464)
    label("Tegafur apparent central volume (V1,Teg/F, L/kg)") # Kim 2017 Table 2, row 'Central volume of distribution for tegafur': V1,Teg/F = 0.0464 L/kg
    lvp <- log(0.137)
    label("Tegafur apparent peripheral volume (V2,Teg/F, L/kg)") # Kim 2017 Table 2, row 'Peripheral volume of distribution for tegafur': V2,Teg/F = 0.137 L/kg

    # ---- Systemic routing of tegafur clearance into the precursor ----
    # Kim 2017 equation 4 adds CLTeg * FMet * C1,Teg to the precursor pool.
    # Note that equation 2 removes the FULL CLTeg from tegafur, so FMet
    # scales a shadow formation flux rather than splitting tegafur
    # elimination -- see vignette Errata.
    fm <- 0.342
    label("Fraction of tegafur clearance routed to 5-FU precursor formation (FMet, unitless)") # Kim 2017 Table 2, row 'Fraction of 5-FU clearance for 5-FU precursor formation': FMet = 0.342

    # ---- 5-FU disposition (apparent, /F) ----
    # CL5FU/F is the only disposition parameter estimated separately by arm.
    lcl_5fu_ctl <- log(3.52)
    label("5-FU apparent clearance, vehicle control arm (CL5FU/F, L/h/kg)") # Kim 2017 Table 2, row 'Clearance for 5-FU in control group': CL5FU,Con/F = 3.52 L/h/kg
    lcl_5fu_sdt <- log(5.93)
    label("5-FU apparent clearance, SDT-pretreated arm (CL5FU/F, L/h/kg)") # Kim 2017 Table 2, row 'Clearance for 5-FU in SDT pretreatment group': CL5FU,Pretre/F = 5.93 L/h/kg
    lq_5fu <- log(1.87)
    label("5-FU apparent distribution clearance (CLd5FU/F, L/h/kg)") # Kim 2017 Table 2, row 'Distribution clearance for 5-FU': CLd5FU/F = 1.87 L/h/kg
    lvc_5fu <- log(0.623)
    label("5-FU apparent central volume (V1,5FU/F, L/kg)") # Kim 2017 Table 2, row 'Central volume of distribution for 5-FU': V1,5FU/F = 0.623 L/kg
    lvp_5fu <- log(0.294)
    label("5-FU apparent peripheral volume (V2,5FU/F, L/kg)") # Kim 2017 Table 2, row 'Peripheral volume of distribution for 5-FU': V2,5FU/F = 0.294 L/kg

    # ---- Between-subject variability ----
    # Kim 2017 Table 2 reports one BSV per parameter in the parenthesised
    # 'Population Mean (BSV)' column, with Section 4.5 stating a log-normal
    # BSV distribution. The column carries bare decimals with no percent
    # sign, so the values are read here as log-scale variances (omega^2)
    # rather than coefficients of variation -- see vignette Errata for the
    # evidence and for what changes under the alternative reading.
    etalka_ctl ~ 0.0101 # Kim 2017 Table 2, Ka,Con BSV = 0.0101
    etalka_sdt ~ 0.0125 # Kim 2017 Table 2, Ka,Pretre BSV = 0.0125
    etalk_precursor_5fu_form_ctl ~ 0.308 # Kim 2017 Table 2, Ka,Met,Con BSV = 0.308
    etalk_precursor_5fu_form_sdt ~ 0.655 # Kim 2017 Table 2, Ka,Met,Pretre BSV = 0.655
    etalk_5fu_form ~ 0.747 # Kim 2017 Table 2, KConv BSV = 0.747
    etalcl ~ 0.00221 # Kim 2017 Table 2, CLTeg/F BSV = 0.00221
    etafm ~ 0.0162 # Kim 2017 Table 2, FMet BSV = 0.0162
    etalcl_5fu_ctl ~ 0.163 # Kim 2017 Table 2, CL5FU,Con/F BSV = 0.163
    etalcl_5fu_sdt ~ 0.025 # Kim 2017 Table 2, CL5FU,Pretre/F BSV = 0.025
    etalq ~ 0.105 # Kim 2017 Table 2, CLdTeg/F BSV = 0.105
    etalq_5fu ~ 0.168 # Kim 2017 Table 2, CLd5FU/F BSV = 0.168
    etalvc ~ 0.729 # Kim 2017 Table 2, V1,Teg/F BSV = 0.729
    etalvc_5fu ~ 0.291 # Kim 2017 Table 2, V1,5FU/F BSV = 0.291
    etalvp ~ 0.0141 # Kim 2017 Table 2, V2,Teg/F BSV = 0.0141
    etalvp_5fu ~ 0.184 # Kim 2017 Table 2, V2,5FU/F BSV = 0.184

    # ---- Residual unexplained variability ----
    # Kim 2017 Section 4.5 states 'Residual model with additive and
    # proportional error was used for tegafur and 5-FU concentrations' but
    # Table 2 reports no residual-error estimates and no other source on
    # disk supplies them. Carried as fixed(0) so the declared error
    # structure is preserved without inventing magnitudes.
    addSd <- fixed(0)
    label("Tegafur additive residual SD (ng/mL; not published)") # Kim 2017 Section 4.5 declares an additive + proportional residual; no value is reported
    propSd <- fixed(0)
    label("Tegafur proportional residual SD (unitless; not published)") # Kim 2017 Section 4.5 declares an additive + proportional residual; no value is reported
    addSd_5fu <- fixed(0)
    label("5-FU additive residual SD (ng/mL; not published)") # Kim 2017 Section 4.5 declares an additive + proportional residual; no value is reported
    propSd_5fu <- fixed(0)
    label("5-FU proportional residual SD (unitless; not published)") # Kim 2017 Section 4.5 declares an additive + proportional residual; no value is reported
  })

  model({
    # 1. Arm-specific parameters. CONMED_SDT switches between the control
    # and SDT-pretreated estimates of Kim 2017 Table 2.
    # Each individual log-parameter is formed on its own simple line so
    # rxode2 recognises the mu-referencing; the arm switch is applied
    # afterwards on the log scale.
    ilka_ctl <- lka_ctl + etalka_ctl
    ilka_sdt <- lka_sdt + etalka_sdt
    ilk_precursor_5fu_form_ctl <- lk_precursor_5fu_form_ctl + etalk_precursor_5fu_form_ctl
    ilk_precursor_5fu_form_sdt <- lk_precursor_5fu_form_sdt + etalk_precursor_5fu_form_sdt
    ilcl_5fu_ctl <- lcl_5fu_ctl + etalcl_5fu_ctl
    ilcl_5fu_sdt <- lcl_5fu_sdt + etalcl_5fu_sdt

    ka <- exp(ilka_ctl * (1 - CONMED_SDT) + ilka_sdt * CONMED_SDT)
    k_precursor_5fu_form <- exp(
      ilk_precursor_5fu_form_ctl * (1 - CONMED_SDT) +
        ilk_precursor_5fu_form_sdt * CONMED_SDT
    )
    cl_5fu <- exp(ilcl_5fu_ctl * (1 - CONMED_SDT) + ilcl_5fu_sdt * CONMED_SDT)

    # 2. Shared parameters.
    k_5fu_form <- exp(lk_5fu_form + etalk_5fu_form)
    cl <- exp(lcl + etalcl)
    q <- exp(lq + etalq)
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp + etalvp)
    fmet <- fm * exp(etafm)
    q_5fu <- exp(lq_5fu + etalq_5fu)
    vc_5fu <- exp(lvc_5fu + etalvc_5fu)
    vp_5fu <- exp(lvp_5fu + etalvp_5fu)

    # 3. Concentrations entering the clearance-parameterised ODEs
    # (Kim 2017 equations 2, 3, 5, 6 are written in terms of C1 and C2).
    Cc <- central / vc
    Cp <- peripheral1 / vp
    Cc_5fu <- central_5fu / vc_5fu
    Cp_5fu <- peripheral1_5fu / vp_5fu

    # 4. ODE system, Kim 2017 Section 4.5 equations 1-6 verbatim.
    # Eq 1: absorption site drains by both routes.
    d/dt(depot) <- -(ka + k_precursor_5fu_form) * depot
    # Eq 2: tegafur central.
    d/dt(central) <- ka * depot - (cl + q) * Cc + q * Cp
    # Eq 3: tegafur peripheral.
    d/dt(peripheral1) <- q * Cc - q * Cp
    # Eq 4: 5-FU precursor (5'-hydroxytegafur), fed pre-systemically from
    # the absorption site and systemically by the FMet fraction of tegafur
    # clearance, and drained by conversion to 5-FU.
    d/dt(precursor_5fu) <- k_precursor_5fu_form * depot + cl * fmet * Cc -
      k_5fu_form * precursor_5fu
    # Eq 5: 5-FU central.
    d/dt(central_5fu) <- k_5fu_form * precursor_5fu - (cl_5fu + q_5fu) * Cc_5fu +
      q_5fu * Cp_5fu
    # Eq 6: 5-FU peripheral.
    d/dt(peripheral1_5fu) <- q_5fu * Cc_5fu - q_5fu * Cp_5fu

    # 5. Observations. Kim 2017 Section 4.5: additive plus proportional
    # residual on both analytes (magnitudes not published; see ini()).
    Cc ~ add(addSd) + prop(propSd)
    Cc_5fu ~ add(addSd_5fu) + prop(propSd_5fu)
  })
}
