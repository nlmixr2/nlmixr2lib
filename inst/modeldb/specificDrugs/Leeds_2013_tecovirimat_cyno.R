Leeds_2013_tecovirimat_cyno <- function() {
  description <- "Preclinical (cynomolgus monkey). Two-compartment population PK model for oral tecovirimat (ST-246) in uninfected and monkeypox-virus-infected cynomolgus monkeys (Leeds 2013), with first-order absorption preceded by a lag time, allometric body-weight scaling (0.75 on CL/F and Q/F, 1 on Vc/F and Vp/F) about a 3.105 kg median weight, a weight-adjusted dose-level power effect on Ka, CL/F and Vc/F, an infection-status effect on Ka and CL/F, and two-occasion inter-occasion variability on Ka and CL/F. Companion human model: modellib('Leeds_2013_tecovirimat_human')."
  reference <- paste(
    "Leeds JM, Fenneteau F, Gosselin NH, Mouksassi MS, Kassir N, Marier JF,",
    "Chen Y, Grosenbach D, Frimm AE, Honeychurch KM, Chinsangaram J,",
    "Tyavanagimatt SR, Hruby DE, Jordan R. (2013). Pharmacokinetic and",
    "pharmacodynamic modeling to determine the dose of ST-246 to protect",
    "against smallpox in humans. Antimicrob Agents Chemother 57(3):1136-1143.",
    "doi:10.1128/AAC.00959-12.",
    sep = " "
  )
  vignette <- "Leeds_2013_tecovirimat"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot       = list(analyte = "tecovirimat", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "tecovirimat", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "tecovirimat", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric scaling with theory-based fixed exponents about the cohort median weight of 3.105 kg, which is the normalisation constant printed inside every covariate formula of Leeds 2013 Table 1 (NHP columns). Results 'NHP POP PK model development': 'Body weight was used as an allometric factor [wt/(median wt)^theta_eff] on CL/F, Vc/F, Q/F, and Vp/F, with theta_eff equal to 0.75 for clearance-related parameters and theta_eff equal to 1 for volume-related parameters'.",
      source_name        = "wt"
    ),
    DIS_INFECT_ACTIVE = list(
      description        = "Active monkeypox-virus (MPXV) infection episode indicator (1 = record falls after intravenous MPXV inoculation, 0 = uninfected).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (uninfected cynomolgus monkey)",
      notes              = "Time-varying per record. The paper's clinical criterion is experimental rather than symptom-based: animals were inoculated intravenously with approximately 5e7 PFU of MPXV per animal (Materials and Methods, 'In vivo study summaries (i) NHP studies'), so a record is infected if it follows inoculation. Study 2 is the design that identifies the effect within animal: 'Monkeys were administered a single dose of ST-246 10 days prior to infection to allow intra-animal comparison of the pharmacokinetics according to infection status', with ST-246 dosing then restarted 4 days after infection. Enters as a log-additive shift on Ka and on CL/F: Leeds 2013 Table 1 prints separate Infected and Uninfected rows for exactly those two parameters, and the Table 1 footnote states 'Infection was a covariate for both the Ka and CL/F, but the overall change in exposure was small'.",
      source_name        = "Status (Infected / Uninfected)"
    ),
    DOSE_TECOVIRIMAT_MGKG = list(
      description        = "Administered tecovirimat (ST-246) dose level, expressed per kg body weight.",
      units              = "mg/kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed within a dosing regimen. Leeds 2013 Table 1 footnote a defines the symbol: 'Dose, dose level in mg/kg'. Enters as a power function normalised to 10 mg/kg on Ka, CL/F and Vc/F, i.e. (dose/10)^theta, matching the printed Table 1 formulae. Study dose levels spanned 0.3 to 300 mg/kg once daily; the population PK dataset covers 0.3-30 mg/kg (study 1, uninfected) and 3-20 mg/kg (studies 2-6, infected). This is the dose LEVEL as a covariate label on the regimen, not the per-record amount: the amount actually administered is the event-table `amt` in mg.",
      source_name        = "dose"
    ),
    OCC = list(
      description        = "Integer-valued occasion / period indicator for inter-occasion-variability multiplexing.",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = "Two occasions are encoded: OCC = 1 for the pre-infection (uninfected) dosing occasion and OCC = 2 for the post-infection (infected) dosing occasion. Leeds 2013 does not state an occasion count; the two-occasion reading follows the paper's own attribution of the IOV to infection status -- the Table 1 title reports 'intraoccasion variability for those parameters altered by the infected state', IOV is printed on exactly the two rows (Ka, CL/F) that carry Infected / Uninfected sub-rows, and Results 'NHP POP PK model development' states 'the proposed two-compartment model demonstrated intraoccasion variability (IOV) for infection status'. The design that supplies both occasions within animal is study 2 (single dose 10 days before infection, then dosing from day 4 after infection). Decomposed inside `model()` into binary indicators `oc1` / `oc2` multiplexing the per-occasion IOV etas on log-Ka and log-CL.",
      source_name        = "OCC"
    )
  )

  population <- list(
    species        = "cynomolgus monkey (Macaca fascicularis)",
    n_subjects     = NA,
    n_studies      = 6L,
    weight_median  = "3.105 kg",
    sex_female_pct = NA,
    disease_state  = "Uninfected cynomolgus monkeys and cynomolgus monkeys with a nearly uniformly lethal systemic monkeypox-virus (MPXV) infection established by intravenous inoculation with approximately 5e7 PFU per animal. Studies 3-6 initiated treatment 3, 4 or 5 days after infection (at lesion onset in study 4).",
    dose_range     = "Oral ST-246 once daily. Study 1 (uninfected, GLP PK): 0.3, 3, 10, 20 and 30 mg/kg, n = 6 per dose (3 male, 3 female). Study 2 (infection-effect GLP study): 3, 10 or 20 mg/kg or vehicle, n = 6 per group, plus a single dose 10 days before infection. Studies 3-6 (efficacy, sparse PK): 0.3 to 300 mg/kg for 14 consecutive days.",
    notes          = "Population PK parameters from Leeds 2013 Table 1, NHP columns. The paper does not report the number of animals contributing to the PK analysis; it reports the concentration counts: 21 of 1,579 non-BQL preclinical plasma concentrations were excluded as outliers, and 1,558 plasma concentrations entered the NHP population PK analysis, none with |CWRES| > 4. Sex was not retained as a covariate in the NHP model (Results 'Human POP PK model development': gender affected Vc/F in humans but not in NHPs). Model fit in Phoenix NLME 6.2.0.416 with FOCE-ELS. The separate 96-animal survival / ROC pharmacodynamic analysis in the same paper is a nonparametric Kaplan-Meier and recursive-ROC analysis with no parametric hazard model, so it is not represented here; its exposure cutoffs are reproduced as validation targets in the vignette."
  )

  ini({
    # ---------------------------------------------------------------------
    # Structural PK parameters -- Leeds 2013 Table 1, 'NHP PK parameter /
    # Population estimate' column. Every NHP formula in Table 1 is printed
    # normalised to wt = 3.105 kg and, where a dose term appears, to
    # dose = 10 mg/kg, so the typical values below are the values at
    # WT = 3.105 kg and DOSE_TECOVIRIMAT_MGKG = 10 mg/kg. The Uninfected
    # row is taken as the reference; the infected shift is a separate
    # covariate-effect parameter below.
    # ---------------------------------------------------------------------
    lka   <- log(0.586)  ; label("Absorption rate constant Ka in uninfected animals at 10 mg/kg (1/h)")                    # Table 1 row 'Ka (h-1)' Uninfected: 0.586 * (dose/10)^0.160
    ltlag <- log(0.302)  ; label("Absorption lag time Tlag (h)")                                                           # Table 1 row 'Tlag (h)' NHP Population estimate = 0.302 (IIV reported as NA)
    lcl   <- log(2.827)  ; label("Apparent oral clearance CL/F in uninfected animals at WT = 3.105 kg and 10 mg/kg (L/h)") # Table 1 row 'CL/F (liters/h)' Uninfected: 2.827 * (wt/3.105)^0.75 * (dose/10)^0.093
    lvc   <- log(20.054) ; label("Apparent central volume of distribution Vc/F at WT = 3.105 kg and 10 mg/kg (L)")         # Table 1 row 'Vc/F (liters)': 20.054 * (wt/3.105)^1 * (dose/10)^0.623
    lq    <- log(3.244)  ; label("Apparent inter-compartmental clearance Q/F at WT = 3.105 kg (L/h)")                      # Table 1 row 'Q/F (liters/h)': 3.244 * (wt/3.105)^0.75
    lvp   <- log(13.34)  ; label("Apparent peripheral volume of distribution Vp/F at WT = 3.105 kg (L)")                   # Table 1 row 'Vp/F (liters)': 13.34 * (wt/3.105)

    # ---------------------------------------------------------------------
    # Allometric exponents. Theory-based and shared across the clearance
    # terms and across the volume terms, so one parameter serves each pair
    # (Results 'NHP POP PK model development': 'Body weight was used as an
    # allometric factor [wt/(median wt)^theta_eff] on CL/F, Vc/F, Q/F, and
    # Vp/F, with theta_eff equal to 0.75 for clearance-related parameters
    # and theta_eff equal to 1 for volume-related parameters (24)'). Cited
    # to an external allometry reference and printed without uncertainty in
    # Table 1, so encoded as structural rather than estimated.
    # ---------------------------------------------------------------------
    e_wt_cl_q  <- fixed(0.75) ; label("Allometric exponent on CL/F and Q/F with body weight (unitless)")     # Results 'NHP POP PK model development'; Table 1 NHP CL/F and Q/F formulae
    e_wt_vc_vp <- fixed(1.0)  ; label("Allometric exponent on Vc/F and Vp/F with body weight (unitless)")    # Results 'NHP POP PK model development'; Table 1 NHP Vc/F and Vp/F formulae

    # ---------------------------------------------------------------------
    # Infection-status effects. Table 1 prints the Infected and Uninfected
    # typical values side by side rather than a shift, so the shift is
    # recovered as the log ratio of the two printed values. The dose and
    # weight terms are identical in the two rows, so the ratio is exactly
    # the infection effect.
    # ---------------------------------------------------------------------
    e_dis_infect_active_ka <- log(0.868 / 0.586) ; label("Log-additive effect of active MPXV infection on Ka (unitless)")   # Table 1 row 'Ka (h-1)': Infected 0.868 vs Uninfected 0.586 (48.1% higher Ka when infected)
    e_dis_infect_active_cl <- log(2.809 / 2.827) ; label("Log-additive effect of active MPXV infection on CL/F (unitless)") # Table 1 row 'CL/F (liters/h)': Infected 2.809 vs Uninfected 2.827 (0.6% lower CL/F when infected)

    # ---------------------------------------------------------------------
    # Weight-adjusted dose-level effects, as power functions normalised to
    # 10 mg/kg (Results 'NHP POP PK model development': '... with covariates
    # of wt on Q/F and Vp/F and weight-adjusted dose on Ka, CL/F, and Vc/F').
    # Reported as estimates with no FIX signal, so left estimated.
    # ---------------------------------------------------------------------
    e_dose_tecovirimat_mgkg_ka <- 0.160 ; label("Power exponent on (dose / 10 mg/kg) for Ka (unitless)")    # Table 1 row 'Ka (h-1)': (dose/10)^0.160, both infection statuses
    e_dose_tecovirimat_mgkg_cl <- 0.093 ; label("Power exponent on (dose / 10 mg/kg) for CL/F (unitless)")  # Table 1 row 'CL/F (liters/h)': (dose/10)^0.093, both infection statuses
    e_dose_tecovirimat_mgkg_vc <- 0.623 ; label("Power exponent on (dose / 10 mg/kg) for Vc/F (unitless)")  # Table 1 row 'Vc/F (liters)': (dose/10)^0.623

    # ---------------------------------------------------------------------
    # Inter-individual variability. Table 1 column 'IIV (%)' gives one
    # percentage per parameter. Read as the standard deviation of the
    # log-scale random effect expressed as a percentage, i.e. the internal
    # variance is omega^2 = (IIV/100)^2. This is the same reading used for
    # the '% CV' IIV column of Jonsson_2011_ethambutol.R (also Antimicrob
    # Agents Chemother); the alternative exact log-normal back-transform
    # omega^2 = log(1 + (IIV/100)^2) is discussed in the vignette Errata and
    # changes omega by 2% on CL/F.
    #
    # No IIV is carried on Tlag: Table 1 reports 'NA' in the NHP IIV column
    # for that row.
    #
    # Leeds 2013 Results also reports 'correlation between CL/F and Vc/F'
    # for the NHP model, but Table 1 tabulates no covariance or correlation
    # coefficient, so the etas below are left uncorrelated. See vignette
    # Errata.
    # ---------------------------------------------------------------------
    etalka ~ 0.0121 # Table 1 'Ka (h-1)' NHP IIV = 11%;  omega^2 = 0.11^2
    etalcl ~ 0.0961 # Table 1 'CL/F' NHP IIV      = 31%;  omega^2 = 0.31^2
    etalvc ~ 0.2209 # Table 1 'Vc/F' NHP IIV      = 47%;  omega^2 = 0.47^2
    etalq  ~ 0.5625 # Table 1 'Q/F' NHP IIV       = 75%;  omega^2 = 0.75^2
    etalvp ~ 0.3025 # Table 1 'Vp/F' NHP IIV      = 55%;  omega^2 = 0.55^2

    # ---------------------------------------------------------------------
    # Inter-occasion variability on log-Ka and log-CL across the two
    # infection-status occasions (Table 1 column 'IOV (%)', populated for
    # exactly the Ka and CL/F rows). nlmixr2 has no equivalent of NONMEM's
    # `$OMEGA BLOCK(1) SAME`, so the second occasion's variance is fixed to
    # the first occasion's estimate, matching Chen_2023_nemonoxacin.R (the
    # other registered two-occasion IOV model) and Jonsson_2011_ethambutol.R.
    # ---------------------------------------------------------------------
    etaiov_ka_1 ~ 0.1024        # Table 1 'Ka (h-1)' NHP IOV = 32%; pi^2 = 0.32^2 (estimated)
    etaiov_ka_2 ~ fixed(0.1024) # SAME-equivalent: equal to the occasion-1 IOV variance
    etaiov_cl_1 ~ 0.0196        # Table 1 'CL/F' NHP IOV = 14%; pi^2 = 0.14^2 (estimated)
    etaiov_cl_2 ~ fixed(0.0196) # SAME-equivalent: equal to the occasion-1 IOV variance

    # ---------------------------------------------------------------------
    # Combined residual error. Table 1 'Error model' block, NHP column.
    # The additive term is printed in ug/liter and is converted to the
    # model's mg/L concentration scale by dividing by 1000.
    # ---------------------------------------------------------------------
    propSd <- 0.30         ; label("Proportional residual error (fraction)") # Table 1 'Error model / Proportional (%)' NHP = 30
    addSd  <- 0.133 / 1000 ; label("Additive residual error (mg/L)")         # Table 1 'Error model / Additive (ug/liter)' NHP = 0.133 ug/L
  })

  model({
    # Decompose the occasion column into binary indicators for IOV
    # multiplexing. OCC = 1 is the pre-infection (uninfected) occasion,
    # OCC = 2 the post-infection (infected) occasion; a single-occasion
    # simulation should pass OCC = 1 so the first IOV eta applies.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)

    iov_ka <- oc1 * etaiov_ka_1 + oc2 * etaiov_ka_2
    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2

    # Individual PK parameters. Allometric weight terms are normalised to
    # the 3.105 kg median weight and the dose terms to 10 mg/kg, exactly as
    # Table 1 prints them.
    ka <- exp(lka + e_dis_infect_active_ka * DIS_INFECT_ACTIVE + etalka + iov_ka) *
      (DOSE_TECOVIRIMAT_MGKG / 10)^e_dose_tecovirimat_mgkg_ka
    cl <- exp(lcl + e_dis_infect_active_cl * DIS_INFECT_ACTIVE + etalcl + iov_cl) *
      (WT / 3.105)^e_wt_cl_q * (DOSE_TECOVIRIMAT_MGKG / 10)^e_dose_tecovirimat_mgkg_cl
    vc <- exp(lvc + etalvc) * (WT / 3.105)^e_wt_vc_vp *
      (DOSE_TECOVIRIMAT_MGKG / 10)^e_dose_tecovirimat_mgkg_vc
    q  <- exp(lq + etalq) * (WT / 3.105)^e_wt_cl_q
    vp <- exp(lvp + etalvp) * (WT / 3.105)^e_wt_vc_vp

    tlag <- exp(ltlag)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Two-compartment model with first-order oral absorption from a depot.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Absorption lag applied to the dosing compartment.
    alag(depot) <- tlag

    # Dose units mg, Vc/F units L, so Cc units mg/L (= ug/mL).
    Cc <- central / vc

    Cc ~ add(addSd) + prop(propSd)
  })
}
