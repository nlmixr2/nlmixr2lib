Chandasana_2024b_lamivudine <- function() {
  description <- "One-compartment oral population PK model with first-order absorption, an absorption lag time and formulation-dependent absolute bioavailability for lamivudine in children with HIV-1 weighing 6 to 40 kg, applied without re-estimation to the ABC/DTG/3TC fixed-dose combination (dispersible tablet and tablet) in IMPAACT 2019 (Chandasana 2024)"
  reference <- paste(
    "Chandasana H, van Dijkman SC, Mehta R, Bush M, Rabie H, Flynn P,",
    "Cressey TR, Acosta EP, Brooks KM; for the IMPAACT 2019 Study Team.",
    "Population pharmacokinetic modeling of abacavir/dolutegravir/lamivudine",
    "to support a fixed-dose combination in children with HIV-1.",
    "Infect Dis Ther. 2024;13(8):1877-1891. doi:10.1007/s40121-024-01008-y.",
    "The lamivudine model itself was developed in the US FDA clinical pharmacology",
    "review for Ziagen (abacavir sulfate) and Epivir (lamivudine)",
    "(Chandasana 2024 reference 16) and is reproduced in Chandasana 2024 Table 3;",
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
      notes              = "Power (allometric) effect on CL (exponent 0.758) and V (exponent 0.677) with reference weight 18.5 kg (Chandasana 2024 Table 3). Both exponents were estimated. Weight range 3.1-66.4 kg in the pooled six-study model-development population (reported as 5.1-66.4 kg in the Chandasana 2024 Discussion); 8.15-39.30 kg in the IMPAACT 2019 external-validation cohort.",
      source_name        = "WT"
    ),
    FORM_SOLUTION = list(
      description        = "Oral-solution formulation indicator (1 = lamivudine oral solution, 0 = solid oral dosage form: tablet, capsule or dispersible tablet)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (solid oral dosage form; absolute bioavailability 0.609)",
      notes              = "Chandasana 2024 Table 3 reports absolute bioavailability F1 = 0.496 for the oral solution and F1 = 0.609 for the tablet; the source model identified a higher absolute bioavailability for solid dosage forms (tablet and capsule) than for the oral solution, consistent with a lamivudine relative-bioavailability study in children (Chandasana 2024 '3TC Pediatric PopPK Model', citing Kasirye 2012). The ABC/DTG/3TC dispersible tablet used in IMPAACT 2019 is treated as a solid dosage form (F1 = 0.609): reproducing Chandasana 2024 Table 4 with F1 = 0.609 matches the published geometric-mean AUC0-24 to within 2.4% in every weight band, whereas F1 = 0.496 under-predicts it by about 26%.",
      source_name        = "formulation"
    ),
    OCC = list(
      description        = "Integer-valued occasion / period indicator for inter-occasion variability on CL, V and Ka",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = "Chandasana 2024 Table 3 reports a single inter-occasion variability magnitude per parameter (CL/F 24.9% CV, V/F 60.0% CV, Ka 19.7% CV) shared across occasions, i.e. the NONMEM $OMEGA BLOCK(1) SAME idiom. Two occasions are encoded here (the minimum that makes IOV operational); set OCC = 1 for every record to reproduce the single steady-state occasion simulated in Chandasana 2024 Table 4. Users needing more occasions extend the oc<k> / etaiov_*_<k> pattern with additional ~ fix(...) slots.",
      source_name        = "OCC"
    )
  )

  compartmentData <- list(
    depot   = list(analyte = "lamivudine", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "lamivudine", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 209L,
    n_studies      = 6L,
    age_range      = "4 months to 19 years",
    weight_range   = "3.1-66.4 kg",
    disease_state  = "Children living with HIV-1 receiving oral lamivudine",
    dose_range     = "Oral lamivudine; in the IMPAACT 2019 confirmatory simulations 90 mg (>=6 to <10 kg), 120 mg (>=10 to <14 kg), 150 mg (>=14 to <20 kg) and 180 mg (>=20 to <25 kg) once daily as the ABC/DTG/3TC dispersible tablet, and 300 mg once daily as the ABC/DTG/3TC tablet (>=25 to <40 kg)",
    notes          = "Model-development population: a pooled analysis of six clinical studies reported in the US FDA clinical pharmacology review for Ziagen and Epivir (Chandasana 2024 reference 16), summarised in Chandasana 2024 '3TC Pediatric PopPK Model' and Table 3. External-validation cohort: IMPAACT 2019 (NCT03760458), an international phase I/II open-label study in children <12 years living with HIV-1 enrolled into five weight bands (>=6 to <10, >=10 to <14, >=14 to <20, >=20 to <25 and >=25 to <40 kg); 55 participants contributed 597 lamivudine intensive and sparse PK samples, median (min-max) baseline age 6.0 (1.00-11.0) years and weight 17.00 (8.15-39.30) kg, 45.5% female, 67% Black and 31% Asian (Chandasana 2024 Results). The existing model was applied to the IMPAACT 2019 data with NONMEM MAXEVAL = 0, i.e. no parameter was re-estimated. The predefined exposure target for dose confirmation was a geometric-mean AUC0-24 of 6.3-26.5 ug*h/mL (Chandasana 2024 Methods)."
  )

  ini({
    # Structural parameters at the reference body weight of 18.5 kg.
    # All values from Chandasana 2024 Table 3 (previously reported 3TC pediatric
    # PopPK parameter estimates; primary source is the FDA Ziagen/Epivir
    # clinical pharmacology review, Chandasana 2024 reference 16).
    #
    # NOTE on the clearance and volume labels. Chandasana 2024 Table 3 labels
    # these "Apparent clearance, CL/F" and "Apparent central volume of
    # distribution, V/F", but the same table also reports an ABSOLUTE
    # bioavailability F1, which the source model applies separately. The two
    # therefore compose as CL/F = 9.16 / F1: reproducing Chandasana 2024 Table 4
    # as AUC0-24 = Dose * F1 / (9.16 * (WT/18.5)^0.758) matches the published
    # geometric means to within 2.4% across all five weight bands, whereas
    # treating 9.16 as already-apparent (F1 omitted) over-predicts them by ~64%.
    # The values are encoded here as true CL and V with F1 applied via f(depot).
    lka  <- log(2.08);  label("Absorption rate constant Ka (1/h)")                        # Chandasana 2024 Table 3: KA = 2.08 1/h (%RSE 9.76)
    ltlag <- log(0.297); label("Absorption lag time (h)")                                 # Chandasana 2024 Table 3: lag time ALAG1 = 0.297 h (%RSE 12.1)
    lcl  <- log(9.16);  label("Clearance CL at WT = 18.5 kg (L/h)")                       # Chandasana 2024 Table 3: CL/F = 9.16 L/h (%RSE 4.49)
    lvc  <- log(23.1);  label("Central volume of distribution V at WT = 18.5 kg (L)")     # Chandasana 2024 Table 3: V/F = 23.1 L (%RSE 4.68)

    # Absolute bioavailability, formulation-specific. Reference is a solid oral
    # dosage form (tablet / capsule / dispersible tablet).
    lfdepot     <- log(0.609); label("Absolute bioavailability F1, solid oral dosage form (unitless)")  # Chandasana 2024 Table 3: absolute bioavailability (F1) tablet PO = 0.609 (%RSE 5.35)
    lfdepot_sol <- log(0.496); label("Absolute bioavailability F1, oral solution (unitless)")           # Chandasana 2024 Table 3: absolute bioavailability (F1) solution PO = 0.496 (%RSE 5.36)

    # Allometric weight exponents; both estimated.
    e_wt_cl <- 0.758; label("Power exponent of body weight on CL, reference 18.5 kg (unitless)")  # Chandasana 2024 Table 3: CL/F~(WT/18.5) = 0.758 (%RSE 7.07)
    e_wt_vc <- 0.677; label("Power exponent of body weight on V, reference 18.5 kg (unitless)")   # Chandasana 2024 Table 3: VL/F~(WT/18.5) = 0.677 (%RSE 8.98)

    # IIV - Chandasana 2024 Table 3 reports interindividual variability as CV%,
    # converted to log-scale variance with omega^2 = log(1 + CV^2). The source
    # covariance matrix is not reproduced in Chandasana 2024 ("Further details
    # about covariance matrix and full model can be found in [16]"), so the etas
    # are encoded as uncorrelated.
    etalcl ~ 0.07862  # Chandasana 2024 Table 3: IIV CL/F = 28.6% CV; log(1 + 0.286^2) = 0.07862
    etalvc ~ 0.10159  # Chandasana 2024 Table 3: IIV V/F  = 32.7% CV; log(1 + 0.327^2) = 0.10159
    etalka ~ 0.46074  # Chandasana 2024 Table 3: IIV KA   = 76.5% CV; log(1 + 0.765^2) = 0.46074

    # IOV on CL, V and Ka - a single shared magnitude across occasions (NONMEM
    # $OMEGA BLOCK(1) SAME); occasions after the first are fixed to the first
    # occasion's variance.
    etaiov_lcl_1 ~ 0.06015        # Chandasana 2024 Table 3: IOV-CL/F = 24.9% CV; log(1 + 0.249^2) = 0.06015
    etaiov_lcl_2 ~ fix(0.06015)   # Chandasana 2024 Table 3: same IOV-CL/F magnitude on every occasion
    etaiov_lvc_1 ~ 0.30748        # Chandasana 2024 Table 3: IOV-V/F   = 60.0% CV; log(1 + 0.600^2) = 0.30748
    etaiov_lvc_2 ~ fix(0.30748)   # Chandasana 2024 Table 3: same IOV-V/F magnitude on every occasion
    etaiov_lka_1 ~ 0.03808        # Chandasana 2024 Table 3: IOV-KA    = 19.7% CV; log(1 + 0.197^2) = 0.03808
    etaiov_lka_2 ~ fix(0.03808)   # Chandasana 2024 Table 3: same IOV-KA magnitude on every occasion

    # Residual error. Chandasana 2024 Table 3 reports a single additive residual
    # error of 0.003 mg/L (= 0.003 ug/mL) together with a "Weighing factor for
    # residual error" of 4.72 (%RSE 7.22). The paper does not give the NONMEM
    # $ERROR construct in which that weighting factor is applied, and it cannot
    # be recovered from any other number in the paper, so only the additive term
    # is encoded here. See the vignette Assumptions and deviations section: this
    # is the one parameter of the three models that this secondary source does
    # not fully determine, and it affects only the residual error, not the
    # typical-value or between-subject predictions validated against Table 4.
    addSd <- 0.003; label("Additive residual error (ug/mL)")  # Chandasana 2024 Table 3: additive error = 0.003 mg/L
  })

  model({
    # Inter-occasion variability, multiplexed by the occasion indicator.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    iov_cl <- oc1 * etaiov_lcl_1 + oc2 * etaiov_lcl_2
    iov_vc <- oc1 * etaiov_lvc_1 + oc2 * etaiov_lvc_2
    iov_ka <- oc1 * etaiov_lka_1 + oc2 * etaiov_lka_2

    # Individual PK parameters. Reference body weight 18.5 kg.
    ka   <- exp(lka + etalka + iov_ka)
    tlag <- exp(ltlag)
    cl   <- exp(lcl + etalcl + iov_cl) * (WT / 18.5)^e_wt_cl
    vc   <- exp(lvc + etalvc + iov_vc) * (WT / 18.5)^e_wt_vc

    kel <- cl / vc

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Formulation-specific absolute bioavailability; the blend selects the
    # solid-dosage-form or oral-solution estimate.
    f(depot)    <- exp(lfdepot * (1 - FORM_SOLUTION) + lfdepot_sol * FORM_SOLUTION)
    alag(depot) <- tlag

    Cc <- central / vc
    Cc ~ add(addSd)
  })
}
