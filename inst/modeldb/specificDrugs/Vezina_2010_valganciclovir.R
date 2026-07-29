Vezina_2010_valganciclovir <- function() {
  description <- "One-compartment population PK model for ganciclovir following oral valganciclovir prophylaxis in pediatric solid organ transplant recipients at risk for Epstein-Barr virus disease (Vezina 2010). First-order absorption with no covariates retained in the final model; doses are mg of valganciclovir uncorrected for molecular weight, and the apparent CL/F and V/F absorb both oral bioavailability and the molar conversion from valganciclovir to ganciclovir."
  reference   <- "Vezina HE, Brundage RC, Nevins TE, Balfour HH Jr. The pharmacokinetics of valganciclovir prophylaxis in pediatric solid organ transplant patients at risk for Epstein-Barr virus disease. Clin Pharmacol Adv Appl. 2010;2:1-7. doi:10.2147/CPAA.S8341"
  vignette    <- "Vezina_2010_valganciclovir"
  units       <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list()

  population <- list(
    n_subjects     = 8L,
    n_studies      = 1L,
    n_observations = 43L,
    age_range      = "1.3 - 6.2 years (median 2.1 years)",
    weight_range   = "9.4 - 19.8 kg (median 14.1 kg)",
    height_range   = "73 - 107 cm (median 89.5 cm)",
    bsa_range      = "0.4 - 0.7 m^2 (median 0.6 m^2; computed as sqrt(weight_kg * height_cm / 3600))",
    crcl_range     = "61.9 - 127 mL/min/1.73 m^2 (median 106; Schwartz method)",
    sex_female_pct = 25,
    race_ethnicity = "Not reported in the publication.",
    disease_state  = "Pediatric solid organ transplant recipients (5 kidney, 3 liver) receiving valganciclovir prophylaxis for prevention of Epstein-Barr virus (EBV)-associated post-transplant lymphoproliferative disorder (PTLD). Subjects were at least six weeks post-transplant with two stable serum creatinine measurements (within +/- 0.2 mg/dL on consecutive occasions at least three days apart). Exclusion: ANC < 500/mm^3, platelets < 20,000/mm^3, hemoglobin < 6.5 g/dL.",
    dose_range     = "Oral valganciclovir suspension (90 mg/mL, compounded at the University of Minnesota Medical Center). Median (range) 11.1 (10.1-12.1) mg/kg every 12 hours, or 7.4 (5.3-11.3) mg/kg every 24 hours, dosed by the transplant center using weight-based valganciclovir adjusted for renal function.",
    regions        = "Single center, University of Minnesota Medical Center / General Clinical Research Center, Minneapolis, MN, USA.",
    co_medications = "Concomitant immunosuppression not separately tabulated; subjects were managed per standard post-transplant care.",
    notes          = "Prospective, open-label PK study. Sampling: a 12-hour intensive visit (0, 1, 2, 4, 6, 8, 12 h post-dose; subjects fasted overnight, dosed after a standardized 641 kcal breakfast) and/or two sparse visits during routine transplant clinic appointments. Two of eight subjects participated in both intensive and sparse visits; six participated only in sparse. 43 plasma samples used in the population PK analysis. Demographics from Table 1; pharmacokinetic results from Table 2. The Discussion notes a non-significant trend toward lower bioavailability (F < 65%) in subjects younger than three years; this signal was NOT incorporated into the final model and is not encoded here."
  )

  ini({
    # Structural PK parameters -- Vezina 2010 Table 2 final-model estimates.
    # One-compartment model with first-order oral absorption (NONMEM ADVAN2 /
    # TRANS2 in the Methods 'Pharmacokinetic analysis' section). Doses are
    # entered as milligrams of valganciclovir, NOT corrected for the molar
    # difference between valganciclovir and ganciclovir; the F factor inside
    # CL/F and V/F absorbs both the oral bioavailability and that molar
    # conversion (Methods: "Doses for valganciclovir were recorded as
    # milligrams of valganciclovir and not corrected for molecular weight
    # differences ... The model estimated parameters, apparent oral
    # clearance (CL/F), and apparent volume of distribution (V/F) implicitly
    # incorporate this difference into the bioavailability (F) term.").
    lcl <- log(7.33) ; label("Apparent oral clearance CL/F (L/h)")            # Vezina 2010 Table 2: CL/F = 7.33 L/h
    lvc <- log(35.1) ; label("Apparent volume of distribution V/F (L)")       # Vezina 2010 Table 2: V/F = 35.1 L
    lka <- log(0.85) ; label("First-order oral absorption rate constant Ka (1/h)")  # Vezina 2010 Table 2: Ka = 0.85

    # Inter-individual variability -- Vezina 2010 Table 2 reports IIV as a
    # CV%. NONMEM uses an exponential-on-log-scale parametrization (P_i =
    # P_pop * exp(eta)), so the log-scale variance is omega^2 = log(CV^2 + 1):
    #   CL/F  CV 36.3% -> log(0.363^2 + 1) = 0.1238
    #   V/F   CV 41.4% -> log(0.414^2 + 1) = 0.1582
    #   Ka    CV 74.3% -> log(0.743^2 + 1) = 0.4397
    # The Methods section ("Between-subject variabilities on CL/F, V/F, and
    # Ka were modeled using a proportional error model. This imposed a log
    # normal distribution on the parameters with results expressed as a
    # CV%.") confirms the lognormal IIV form.
    etalcl ~ 0.1238  # Vezina 2010 Table 2 IIV CL/F = 36.3%
    etalvc ~ 0.1582  # Vezina 2010 Table 2 IIV V/F = 41.4%
    etalka ~ 0.4397  # Vezina 2010 Table 2 IIV Ka  = 74.3%

    # Residual unexplained variability -- Vezina 2010 Methods ("Residual
    # unexplained variability was assumed to have a log normal distribution
    # and was also expressed as a CV%.") and Table 2 (RUV CV% = 33.5%). The
    # nlmixr2 prop() form Y = Cc * (1 + propSd * eps) approximates the
    # NONMEM lognormal residual at this CV magnitude (within ~0.6%); the CV
    # is taken as the fractional propSd.
    propSd <- 0.335 ; label("Proportional residual error (fraction)")  # Vezina 2010 Table 2: RUV CV = 33.5%
  })

  model({
    # Individual PK parameters with lognormal IIV. No covariates retained in
    # the final model (Vezina 2010 Results: "none of the patient-specific
    # covariates tested contributed statistically significant information
    # about ganciclovir pharmacokinetic parameters"; covariates evaluated
    # were age, weight, height, body surface area, Schwartz creatinine
    # clearance, transplant type, and sex).
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    ka <- exp(lka + etalka)

    kel <- cl / vc

    # One-compartment depot model. Dose lands in `depot` as mg of
    # valganciclovir; absorption to `central` is first-order (instantaneous
    # conversion to ganciclovir per Methods). F = 1 implicit because the
    # apparent CL/F and V/F already absorb both the true oral
    # bioavailability and the valganciclovir->ganciclovir molar conversion.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Internal mass / volume is mg / L = ug/mL; the assay reports plasma
    # ganciclovir in ng/mL (Methods 'Quantitation of ganciclovir in plasma',
    # Figure 1). Multiply by 1000 to match the assay scale and the units
    # field declared above.
    Cc <- central / vc * 1000
    Cc ~ prop(propSd)
  })
}
