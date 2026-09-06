Krishna_2011_anacetrapib <- function() {
  description <- paste0(
    "Two-compartment population PK model with first-order absorption and an ",
    "absorption lag for oral anacetrapib, a cholesteryl ester transfer ",
    "protein (CETP) inhibitor, in healthy volunteers and patients with ",
    "dyslipidemia (Krishna 2011). Parameterised in absolute CL / Vc / Q / Vp ",
    "plus an explicitly modelled bioavailability F1, so this is a rare oral ",
    "model in which F1 is identified rather than absorbed into CL/F. F1 is ",
    "the product of three multiplicative terms (Eqs. 4-7): a saturable ",
    "reduction with increasing dose whose half-maximal dose is raised ",
    "six-fold by a high-fat meal; a four-level prandial factor (fasted, ",
    "low-fat, patient-selected AHA TLC diet, high-fat); and, for the ",
    "liquid-filled-capsule formulation only, an Emax increase with the ",
    "number of capsules per dose, which acts as a surrogate for the amount ",
    "of co-administered liquid surfactant. Inter-individual variability is a ",
    "correlated 2x2 block on CL and F1; residual error is proportional plus ",
    "additive. Two companion files carry the paper's exposure-response ",
    "layers: Krishna_2011_anacetrapib_hdlc.R and ",
    "Krishna_2011_anacetrapib_ldlc.R."
  )
  reference <- paste0(
    "Krishna R, Bergman AJ, Green M, Dockendorf MF, Wagner JA, Dykstra K. ",
    "Model-based development of anacetrapib, a novel cholesteryl ester ",
    "transfer protein inhibitor. AAPS J. 2011;13(2):179-190. ",
    "doi:10.1208/s12248-011-9254-0. Structural model: Table I and Eqs. 3-7."
  )
  vignette <- "Krishna_2011_anacetrapib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(
      analyte = "anacetrapib", units = "mg",
      specimen = "administration site", verified = TRUE
    ),
    central = list(
      analyte = "anacetrapib", units = "mg",
      specimen = "plasma", verified = TRUE
    ),
    peripheral1 = list(
      analyte = "anacetrapib", units = "mg",
      specimen = "plasma", verified = TRUE
    )
  )

  covariateData <- list(
    DOSE = list(
      description        = "Administered anacetrapib dose amount for the dose record",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Enters Eq. 6 as the saturable dose effect on bioavailability, ",
        "DG1 = 1 - Dmax * DOSE / (DOSE + dose_50 + e_fed_highfat_dose_50 * ",
        "FED_HIGHFAT). Krishna 2011 describes it as 'dose_i is the ",
        "anacetrapib dose for the ith subject', i.e. the nominal dose level ",
        "of the regimen the subject received, which for the once-daily ",
        "regimens studied equals the per-administration amount. Supply it as ",
        "a data column equal to the dose-record `amt` in mg rather than ",
        "relying on rxode2's podo(), which is silent about the difference ",
        "between the nominal and the F-adjusted amount. Higher doses have ",
        "LOWER bioavailability: the half-maximal dose is 55 mg fasted or on ",
        "a low-fat / patient-selected diet, so F1 has already fallen to 35% ",
        "of its zero-dose asymptote at the 100 mg phase III dose."
      ),
      source_name        = "dose"
    ),
    FED_HIGHFAT = list(
      description        = "High-fat-meal dosing indicator, 1 = dose taken after the standard high-fat breakfast, 0 = any other prandial state",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not a high-fat meal). The prandial reference level of this model is the FASTED state, which is encoded as all three fed indicators = 0, not by a separate column.",
      notes              = paste0(
        "Acts twice. In Eq. 5 it selects the prandial multiplier ",
        "e_fed_highfat_fdepot = 2.7 on F1, and in Eq. 6 it adds ",
        "e_fed_highfat_dose_50 = 274 mg to the half-maximal dose, raising it ",
        "from 55 mg to 329 mg (95% CI 123 to 535) and so mitigating the ",
        "dose-dependent fall in bioavailability. Krishna 2011 Clinical ",
        "Trial Simulations: the standard high-fat meal was 827 kcal with 57% ",
        "of calories from fat (two fried or scrambled eggs, two strips of ",
        "bacon, two slices of toast with two pats of butter, 113 g hash ",
        "browns, 240 mL whole milk), matching the FDA high-fat definition. ",
        "Mutually exclusive with FED_LOWFAT and FED_PATIENTSELECTED."
      ),
      source_name        = "I_HF"
    ),
    FED_LOWFAT = list(
      description        = "Low-fat-meal dosing indicator, 1 = dose taken after the standard low-fat breakfast, 0 = any other prandial state",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not a low-fat meal); the fasted reference is all three fed indicators = 0.",
      notes              = paste0(
        "Selects the prandial multiplier e_fed_lowfat_fdepot = 2.4 on F1 in ",
        "Eq. 5. Krishna 2011 Clinical Trial Simulations: the standard ",
        "low-fat meal was 373 kcal with 20% of calories from fat (two slices ",
        "of toasted white bread, one teaspoon low-fat margarine, one ",
        "tablespoon jelly, 5 oz skim milk, 5 oz orange juice). Mutually ",
        "exclusive with FED_HIGHFAT and FED_PATIENTSELECTED."
      ),
      source_name        = "I_LF"
    ),
    FED_PATIENTSELECTED = list(
      description        = "Patient-selected-meal dosing indicator, 1 = dose taken with a meal the participant chose under protocol dietary instruction, 0 = any other prandial state",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not a patient-selected meal); the fasted reference is all three fed indicators = 0.",
      notes              = paste0(
        "Selects the prandial multiplier e_fed_patientselected_fdepot = 3.3 ",
        "on F1 in Eq. 5 -- the LARGEST of the three prandial multipliers, ",
        "larger than the standard high-fat meal's 2.7. Krishna 2011 Clinical ",
        "Trial Simulations: the phase Ib and phase IIb trials instructed ",
        "participants to select their own meal conforming to the American ",
        "Heart Association's TLC diet, 'similar in fat and caloric content ",
        "to a low-fat meal'. The estimate is nonetheless kept separate from ",
        "FED_LOWFAT because the model resolves the two strata as distinct ",
        "coefficients (3.3 vs 2.4, a 37% difference in relative ",
        "bioavailability); collapsing them would misreproduce every ",
        "patient-selected prediction, including the 100 mg phase III ",
        "dose-selection case. Mutually exclusive with FED_HIGHFAT and ",
        "FED_LOWFAT."
      ),
      source_name        = "I_PB"
    ),
    NDOSEUNITS = list(
      description        = "Number of solid oral dosage units administered in the dose",
      units              = "(count)",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Number of liquid-filled capsules making up the dose. Krishna 2011 ",
        "Population PK Model: 'The number of capsules was considered as a ",
        "surrogate for the amount of liquid surfactant contained in the ",
        "capsule unit making up a given dose.' It is NOT recoverable from ",
        "DOSE, because capsule strength differs across the eight pooled ",
        "studies -- which is why the authors carried it as its own ",
        "covariate. It enters twice, both times gated on ",
        "FORM_ANACETRAPIB_LFC: as an Emax term in Eq. 7 (",
        "DG2 = 1 - ndoseunits_50 / (ndoseunits_50 + NDOSEUNITS), rising with ",
        "capsule count for every prandial state) and, in fasted subjects ",
        "only, as the exponential Eq. 5 term ",
        "exp(e_ndoseunits_fdepot * (NDOSEUNITS - 5.5)), centred on the ",
        "median 5.5 capsules among fasted subjects. Set NDOSEUNITS to the ",
        "physical unit count (1 for the single hot-melt-extruded tablet); ",
        "for tablet records the value is multiplied out by ",
        "FORM_ANACETRAPIB_LFC = 0 and does not affect the prediction."
      ),
      source_name        = "N_CAP"
    ),
    FORM_ANACETRAPIB_LFC = list(
      description        = "Anacetrapib Imwitor/Tween liquid-filled-capsule formulation indicator, 1 = liquid-filled capsule (LFC), 0 = hot-melt-extruded (HME) tablet",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (hot-melt-extruded tablet, the intended commercial and phase III formulation)",
      notes              = paste0(
        "Gates both NDOSEUNITS terms off for tablet records. Table I ",
        "contains no formulation theta of its own: Krishna 2011 states that ",
        "the phase IIb model update added 'formulation effects, capsule ",
        "effects', and the formulation effect IS the capsule effect, because ",
        "the capsule count is only defined for the LFC. The gating is not ",
        "printed in Eqs. 4-7 but is fixed by arithmetic and by the Methods ",
        "statement that the model assessed 'number of capsules per dose ",
        "(for LFC only)': the paper reports F1 = 0.35 for a fasted 100 mg ",
        "TABLET, and only DG2 = 1 with a unit fasted Eq. 5 term reproduces ",
        "it (1 - 100/155 = 0.3548). Carrying DG2 for the tablet instead ",
        "gives 0.21 at NDOSEUNITS = 1 and 0.32 at NDOSEUNITS = 5.5. The ",
        "same gating then reproduces the paper's three printed relative ",
        "bioavailabilities exactly -- see the vignette source-trace table. ",
        "Six of the eight pooled studies and 474 subjects used the LFC; two ",
        "trials (a 78-subject bridging study and a 24-subject food-effect ",
        "study) used the HME tablet."
      ),
      source_name        = "formulation type (HME tablets vs. LFC)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 576L,
    n_studies      = 8L,
    patients_pct   = 60,
    disease_state  = "Pooled healthy volunteers and patients with dyslipidemia; approximately 60% of the 576 subjects were patients (Krishna 2011 Results and Interpretation).",
    dose_range     = "Oral anacetrapib single and multiple doses spanning at least 10-400 mg once daily across the eight pooled phase I / phase IIb studies; the phase IIb atorvastatin co-administration arm used atorvastatin 20 mg once daily, and 100 mg once daily as the hot-melt-extruded tablet was the dose and formulation selected for phase III.",
    formulations   = "Imwitor/Tween liquid-filled capsule (LFC; six studies, 474 subjects) and hot-melt-extruded tablet (HME; a 78-subject bridging study and a 24-subject food-effect study).",
    prandial_states = "Four prandial strata: overnight fasted; standard low-fat breakfast (373 kcal, 20% fat); standard high-fat breakfast (827 kcal, 57% fat); and a patient-selected meal conforming to the American Heart Association TLC diet, used in the phase Ib and phase IIb trials.",
    notes          = paste0(
      "Per-study designs, doses, formulations and sampling schedules are in ",
      "Table IA of the Electronic Supplementary Material, which is not on ",
      "disk (see the vignette Errata). Every parameter of all three models ",
      "is in the main-text tables, so the missing supplement costs only ",
      "demographic detail: Krishna 2011 publishes no age, weight, sex or ",
      "race distribution for the popPK analysis set in the main text, and ",
      "no demographic covariate is in the model."
    )
  )

  ini({
    # =========================================================================
    # Krishna 2011 Table I (Population PK Model Parameter Estimates).
    #
    # NONMEM v5.0, FOCE INTER. The paper reports both variance-scale and
    # SD-scale-looking numbers in one table; the IIV rows are demonstrably
    # VARIANCES, because Table I also reports Cov(omega_CL, omega_F1) = 0.04
    # and a covariance is only defined against variances. Reading 0.098 and
    # 0.21 as SDs would imply a correlation of
    # 0.04 / sqrt(0.0096 * 0.0441) = 1.94, which is impossible. The residual
    # rows are read on the same variance scale by parallel construction, so
    # both are converted to the SD scale nlmixr2 expects -- see the vignette
    # Assumptions and deviations section, which also records the alternative
    # reading.
    # =========================================================================
    lka <- fixed(log(0.48)); label("Absorption rate constant (1/h)")                            # Table I: Ka_TV = 0.48 1/h, "Fixed"; footnote a: "Parameter was fixed at previously estimated value to achieve convergence (rounding errors were preventing convergence)"
    lcl <- log(7.6);         label("Systemic clearance from the central compartment (L/h)")      # Table I: CL_TV = 7.6 L/h (SE 1.0; %CV 13.4)
    lvc <- log(55);          label("Central volume of distribution (L)")                         # Table I: V2_TV = 55 L (SE 7; %CV 13.4)
    lq  <- log(5.3);         label("Intercompartmental clearance (L/h)")                         # Table I: Q_TV = 5.3 L/h (SE 0.7; %CV 13.9)
    lvp <- log(244);         label("Peripheral volume of distribution (L)")                      # Table I: V3_TV = 244 L (SE 33; %CV 13.6)
    ltlag <- log(0.918);     label("Absorption lag time (h)")                                    # Table I: Tlag = 0.918 h (SE 0.006; %CV 0.67)

    # Bioavailability is fully determined by the covariate terms of Eqs. 4-7:
    # Table I reports no F1 scale parameter, so the F1 = 1 anchor is explicit
    # and fixed here rather than left implicit.
    lfdepot <- fixed(log(1)); label("Bioavailability scale anchor for the depot (fraction)")     # Table I reports no F1 typical value; Eq. 4 gives F1_i = DG1 * DG2 * Feff_i * exp(eta_F1_i)

    # ---- Saturable dose effect on bioavailability (Eq. 6) ----
    # DG1 = 1 - dose_max * DOSE / (DOSE + dose_50 + e_fed_highfat_dose_50 * FED_HIGHFAT)
    # Same functional form as the prednisolone-on-F term of
    # Storset_2014_tacrolimus.R, whose <cov>_max / <cov>_50 naming is followed.
    dose_max <- fixed(1);          label("Maximum fractional reduction in bioavailability due to dose (unitless)")            # Table I: D_max = 1, "Fixed"
    dose_50  <- 55;                label("Anacetrapib dose giving half-maximal reduction in bioavailability (mg)")            # Table I: D50 = 55 mg (SE 5; %CV 10.0); Final Model Description gives 95% CI 44 to 66 mg
    e_fed_highfat_dose_50 <- 274;  label("Additive shift of a high-fat meal on the half-maximal dose (mg)")                   # Table I: theta_D = 274 (SE 99.5; %CV 36.3); Final Model Description: high-fat raises the half-maximal dose from 55 to 329 mg (95% CI 123 to 535)

    # ---- Prandial multipliers on bioavailability (Eq. 5) ----
    # Exactly one of the four prandial states applies per dose record; the
    # fasted state is the complement of the three fed indicators and its
    # multiplier is the exponential dosage-unit term rather than a theta.
    e_fed_highfat_fdepot         <- 2.7;  label("Bioavailability multiplier for a high-fat meal (unitless)")             # Table I: theta_HF = 2.7 (SE 0.3; %CV 11.7)
    e_fed_lowfat_fdepot          <- 2.4;  label("Bioavailability multiplier for a low-fat meal (unitless)")              # Table I: theta_LF = 2.4 (SE 0.4; %CV 14.3)
    e_fed_patientselected_fdepot <- 3.3;  label("Bioavailability multiplier for a patient-selected meal (unitless)")     # Table I: theta_PB = 3.3 (SE 0.4; %CV 13.1)

    # ---- Dosage-unit-count effects on bioavailability (Eqs. 5 and 7) ----
    e_ndoseunits_fdepot <- 0.07;  label("Exponential effect of dosage-unit count on bioavailability in fasted subjects (1/unit)")   # Table I: theta_Fast = 0.07 (SE 0.02; %CV 21.1)
    ndoseunits_50       <- 0.67;  label("Dosage-unit count giving half-maximal effect on bioavailability (count)")                  # Table I: Cap50 = 0.67 (SE 0.14; %CV 20.4)

    # ---- IIV: correlated 2x2 block on CL and F1 (Table I) ----
    # omega^2_CL = 0.098 -> 32.2 %CV; omega^2_F1 = 0.21 -> 48.3 %CV;
    # correlation = 0.04 / sqrt(0.098 * 0.21) = 0.279.
    # Table I: omega_CL = 0.098 (SE 9.4e-3), Cov(omega_CL, omega_F1) = 0.04
    # (SE 8.2e-3), omega_F1 = 0.21 (SE 0.02). The source trace sits above the
    # block rather than beside it: comments inside an omega c(...) are fragile
    # for both rxode2's parser and nlmixr2lib's source reconstructor.
    etalcl + etalfdepot ~ c(0.098,
                            0.04, 0.21)

    # ---- Residual error (Table I), converted from the variance scale ----
    propSd <- 0.428952; label("Proportional residual error (fraction)")   # Table I: sigma_proportional = 0.184 (SE 7.4e-3) read as a variance; sqrt(0.184) = 0.428952
    addSd  <- 5.692100; label("Additive residual error (ng/mL)")          # Table I: sigma_additive = 32.4 (SE 9.16) read as a variance; sqrt(32.4) = 5.692100
  })

  model({
    # ---- Individual disposition parameters (Eq. 3) ----
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc)
    q  <- exp(lq)
    vp <- exp(lvp)
    ka <- exp(lka)
    tlag <- exp(ltlag)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ---- Bioavailability (Eqs. 4-7) ----
    # Eq. 6: saturable reduction with dose; a high-fat meal adds 274 mg to the
    # half-maximal dose.
    dose50 <- dose_50 + e_fed_highfat_dose_50 * FED_HIGHFAT
    dg1 <- 1 - dose_max * DOSE / (DOSE + dose50)

    # Eq. 7: Emax increase with the number of dosage units. Defined only for
    # the liquid-filled capsule, so FORM_ANACETRAPIB_LFC = 0 collapses it to 1.
    dg2 <- 1 - FORM_ANACETRAPIB_LFC * ndoseunits_50 / (ndoseunits_50 + NDOSEUNITS)

    # Eq. 5: exactly one prandial state applies. The fasted state is the
    # complement of the three fed indicators and carries the exponential
    # dosage-unit term centred on the median 5.5 capsules; that term is also
    # capsule-only, so it collapses to 1 for tablet records.
    fasted <- 1 - FED_HIGHFAT - FED_LOWFAT - FED_PATIENTSELECTED
    feff <- fasted * exp(e_ndoseunits_fdepot * FORM_ANACETRAPIB_LFC * (NDOSEUNITS - 5.5)) +
      e_fed_highfat_fdepot * FED_HIGHFAT +
      e_fed_lowfat_fdepot * FED_LOWFAT +
      e_fed_patientselected_fdepot * FED_PATIENTSELECTED

    # Eq. 4
    fdepot <- exp(lfdepot + etalfdepot) * dg1 * dg2 * feff

    # ---- ODE system ----
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    f(depot)    <- fdepot
    alag(depot) <- tlag

    # Dose in mg, volumes in L, so central / vc is mg/L; x1000 gives ng/mL,
    # the units of Table I's additive residual and of both PD models' driver.
    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
