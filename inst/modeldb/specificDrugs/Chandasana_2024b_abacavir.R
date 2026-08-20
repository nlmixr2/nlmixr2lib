Chandasana_2024b_abacavir <- function() {
  description <- "Two-compartment oral population PK model with first-order absorption and elimination for abacavir in children with HIV-1 weighing 6 to 40 kg, applied without re-estimation to the ABC/DTG/3TC fixed-dose combination (dispersible tablet and tablet) in IMPAACT 2019 (Chandasana 2024)"
  reference <- paste(
    "Chandasana H, van Dijkman SC, Mehta R, Bush M, Rabie H, Flynn P,",
    "Cressey TR, Acosta EP, Brooks KM; for the IMPAACT 2019 Study Team.",
    "Population pharmacokinetic modeling of abacavir/dolutegravir/lamivudine",
    "to support a fixed-dose combination in children with HIV-1.",
    "Infect Dis Ther. 2024;13(8):1877-1891. doi:10.1007/s40121-024-01008-y.",
    "The abacavir model itself was developed in the US FDA clinical pharmacology",
    "review for Ziagen (abacavir sulfate) and Epivir (lamivudine)",
    "(Chandasana 2024 reference 16) and is reproduced in Chandasana 2024 Table 1;",
    "Chandasana 2024 applied it to IMPAACT 2019 without re-estimation",
    "(NONMEM MAXEVAL = 0 external validation).",
    sep = " "
  )
  vignette <- "Chandasana_2024_abacavir_dolutegravir_lamivudine_pediatric"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power (allometric) effect on CL/F and V2/F with reference weight 15.6 kg (Chandasana 2024 Table 1). Both exponents were FIXED in the source model. Weight range 4.6-61.3 kg in the pooled six-study model-development population; 8.15-39.30 kg in the IMPAACT 2019 external-validation cohort. Q/F and V3/F are NOT weight-scaled in the source model - Table 1 lists weight exponents only for CL/F and V2/F.",
      source_name        = "WT"
    ),
    STUDY_ARROW_PART2 = list(
      description        = "ARROW PK Substudy Part 2 study indicator (1 = record from ARROW PK Substudy Part 2, 0 = any other pooled study, including IMPAACT 2019)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (all pooled studies other than ARROW PK Substudy Part 2; relative bioavailability fixed at 1)",
      notes              = "A study-specific relative bioavailability term (F1) was included to describe the substantially higher observed abacavir exposure in the ARROW PK Substudy Part 2 despite administration of similar doses and formulations (Chandasana 2024 'ABC Pediatric PopPK Model'). The effect is formulation-specific within that substudy: 1.62 for the tablet and 1.75 for the oral solution (Chandasana 2024 Table 1), so it is applied jointly with FORM_SOLUTION. IMPAACT 2019 records carry STUDY_ARROW_PART2 = 0 and therefore relative bioavailability 1.",
      source_name        = "F1 study term"
    ),
    FORM_SOLUTION = list(
      description        = "Oral-solution formulation indicator (1 = abacavir oral solution, 0 = solid oral dosage form: tablet or dispersible tablet)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (solid oral dosage form)",
      notes              = "Used only in combination with STUDY_ARROW_PART2 to select between the two ARROW PK Substudy Part 2 relative-bioavailability estimates (tablet 1.62, solution 1.75; Chandasana 2024 Table 1). Outside that substudy the source model applies no formulation effect on abacavir bioavailability, so this covariate has no effect when STUDY_ARROW_PART2 = 0.",
      source_name        = "F1 formulation term"
    ),
    OCC = list(
      description        = "Integer-valued occasion / period indicator for inter-occasion variability on CL/F",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = "Chandasana 2024 Table 1 reports a single inter-occasion variability magnitude on CL/F (29.2% CV) shared across occasions, i.e. the NONMEM $OMEGA BLOCK(1) SAME idiom. Two occasions are encoded here (the minimum that makes IOV operational); set OCC = 1 for every record to reproduce the single steady-state occasion simulated in Chandasana 2024 Table 4. Users needing more occasions extend the oc<k> / etaiov_lcl_<k> pattern with additional ~ fix(0.08182) slots.",
      source_name        = "OCC"
    )
  )

  compartmentData <- list(
    depot       = list(analyte = "abacavir", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "abacavir", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "abacavir", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 169L,
    n_studies      = 6L,
    age_range      = "5 months to 13 years",
    weight_range   = "4.6-61.3 kg",
    disease_state  = "Children living with HIV-1 receiving oral abacavir",
    dose_range     = "Oral abacavir; in the IMPAACT 2019 confirmatory simulations 180 mg (>=6 to <10 kg), 240 mg (>=10 to <14 kg), 300 mg (>=14 to <20 kg) and 360 mg (>=20 to <25 kg) once daily as the ABC/DTG/3TC dispersible tablet, and 600 mg once daily as the ABC/DTG/3TC tablet (>=25 to <40 kg)",
    notes          = "Model-development population: a pooled analysis of six clinical studies reported in the US FDA clinical pharmacology review for Ziagen and Epivir (Chandasana 2024 reference 16), summarised in Chandasana 2024 'ABC Pediatric PopPK Model' and Table 1. External-validation cohort: IMPAACT 2019 (NCT03760458), an international phase I/II open-label study in children <12 years living with HIV-1 enrolled into five weight bands (>=6 to <10, >=10 to <14, >=14 to <20, >=20 to <25 and >=25 to <40 kg); 55 participants contributed 590 abacavir intensive and sparse PK samples, median (min-max) baseline age 6.0 (1.00-11.0) years and weight 17.00 (8.15-39.30) kg, 45.5% female, 67% Black and 31% Asian (Chandasana 2024 Results). The existing model was applied to the IMPAACT 2019 data with NONMEM MAXEVAL = 0, i.e. no parameter was re-estimated."
  )

  ini({
    # Structural parameters at the reference body weight of 15.6 kg, outside the
    # ARROW PK Substudy Part 2 (relative bioavailability 1).
    # All values from Chandasana 2024 Table 1 (previously reported ABC pediatric
    # PopPK parameter estimates; primary source is the FDA Ziagen/Epivir
    # clinical pharmacology review, Chandasana 2024 reference 16).
    lka <- log(0.85); label("Absorption rate constant Ka (1/h)")                                       # Chandasana 2024 Table 1: KA = 0.85 1/h (%RSE 2.31)
    lcl <- log(16.3); label("Apparent oral clearance CL/F at WT = 15.6 kg (L/h)")                      # Chandasana 2024 Table 1: CL/F = 16.3 L/h (%RSE 3.62)
    lvc <- log(10.1); label("Apparent central volume of distribution V2/F at WT = 15.6 kg (L)")        # Chandasana 2024 Table 1: V2/F = 10.1 L (%RSE 7.56)
    lq  <- log(1.69); label("Apparent intercompartmental clearance Q/F (L/h)")                         # Chandasana 2024 Table 1: Q/F = 1.69 L/h (%RSE 7.87)
    lvp <- log(23.0); label("Apparent peripheral volume of distribution V3/F (L)")                     # Chandasana 2024 Table 1: V3/F = 23.0 L (%RSE 17.4)

    # Allometric weight exponents; both reported as FIX in Chandasana 2024
    # Table 1, so they are encoded with fixed(). Table 1 lists weight exponents
    # only for CL/F and V2/F - Q/F and V3/F carry no weight effect.
    e_wt_cl <- fixed(0.794); label("Power exponent of body weight on CL/F, reference 15.6 kg (unitless)")   # Chandasana 2024 Table 1: CL/F~(WT/15.6) = 0.794 FIX
    e_wt_vc <- fixed(0.698); label("Power exponent of body weight on V2/F, reference 15.6 kg (unitless)")   # Chandasana 2024 Table 1: V/2F~(WT/15.6) = 0.698 FIX

    # Study-specific relative bioavailability for the ARROW PK Substudy Part 2,
    # encoded on the log scale so model() uses exp(effect * indicator), matching
    # the canonical (theta^indicator) parameterisation. Reference (every other
    # pooled study, including IMPAACT 2019) is relative bioavailability 1.
    e_arrow_tab_fdepot <- log(1.62); label("log relative bioavailability of the ARROW PK Substudy Part 2 tablet (unitless; exp = 1.62)")   # Chandasana 2024 Table 1: F, tablet ARROW PK Substudy Part 2 = 1.62 (%RSE 8.02)
    e_arrow_sol_fdepot <- log(1.75); label("log relative bioavailability of the ARROW PK Substudy Part 2 oral solution (unitless; exp = 1.75)") # Chandasana 2024 Table 1: F, solution ARROW PK Substudy Part 2 = 1.75 (%RSE 8.23)

    # IIV - Chandasana 2024 Table 1 reports interindividual variability as CV%,
    # converted to log-scale variance with omega^2 = log(1 + CV^2).
    # The source covariance matrix (eta correlations) is not reproduced in
    # Chandasana 2024 ("Further details about covariance matrix and full model
    # can be found in [16]"), so the etas are encoded as uncorrelated.
    etalcl ~ 0.12378  # Chandasana 2024 Table 1: IIV CL/F = 36.3% CV; log(1 + 0.363^2) = 0.12378
    etalvc ~ 0.23851  # Chandasana 2024 Table 1: IIV V2/F = 51.9% CV; log(1 + 0.519^2) = 0.23851
    etalq  ~ 0.37915  # Chandasana 2024 Table 1: IIV Q/F  = 67.9% CV; log(1 + 0.679^2) = 0.37915
    etalvp ~ 0.61224  # Chandasana 2024 Table 1: IIV V3/F = 91.9% CV; log(1 + 0.919^2) = 0.61224

    # IOV on CL/F - a single shared magnitude across occasions (NONMEM
    # $OMEGA BLOCK(1) SAME); the second and later occasions are fixed to the
    # first occasion's variance.
    etaiov_lcl_1 ~ 0.08182        # Chandasana 2024 Table 1: IOV-CL/F = 29.2% CV; log(1 + 0.292^2) = 0.08182
    etaiov_lcl_2 ~ fix(0.08182)   # Chandasana 2024 Table 1: same IOV-CL/F magnitude on every occasion

    # Residual error - proportional only.
    propSd <- 0.375; label("Proportional residual error (fraction)")  # Chandasana 2024 Table 1: proportional error = 37.5%
  })

  model({
    # Inter-occasion variability on CL/F, multiplexed by the occasion indicator.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    iov_cl <- oc1 * etaiov_lcl_1 + oc2 * etaiov_lcl_2

    # Study- and formulation-specific relative bioavailability. Both indicators
    # are 0 for IMPAACT 2019 and for every pooled study other than the ARROW PK
    # Substudy Part 2, giving relative bioavailability 1.
    frel <- exp(e_arrow_tab_fdepot * STUDY_ARROW_PART2 * (1 - FORM_SOLUTION) +
                e_arrow_sol_fdepot * STUDY_ARROW_PART2 * FORM_SOLUTION)

    # Individual PK parameters. Reference body weight 15.6 kg.
    ka <- exp(lka)
    cl <- exp(lcl + etalcl + iov_cl) * (WT / 15.6)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 15.6)^e_wt_vc
    q  <- exp(lq  + etalq)
    vp <- exp(lvp + etalvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    f(depot) <- frel

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
