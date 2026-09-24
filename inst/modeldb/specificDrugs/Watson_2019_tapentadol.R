Watson_2019_tapentadol <- function() {
  description <- "One-compartment population PK model for tapentadol oral solution in 92 children and adolescents aged 2 to <18 years with moderate-to-severe acute postsurgical pain (Watson 2019). First-order absorption from a depot compartment into the central compartment preceded by an absorption lag time, and first-order elimination; clearance and volume are apparent (CL/F, V/F) because only the oral solution was studied. Body weight is the only retained covariate, entering CL/F and V/F as a power function normalized to the 45 kg population median with exponents estimated (0.638 and 0.847) rather than fixed to the allometric 0.75 and 1. Inter-individual variability is exponential on CL/F, V/F and Ka with a full 3x3 correlation block, and the residual error is combined proportional plus additive. Age, sex, creatinine clearance, AST, ALT, ALP and bilirubin were screened by stepwise covariate modeling and none survived backward elimination, so the final model equals the weight-only base model."
  reference <- paste(
    "Watson E, Khandelwal A, Freijer J, van den Anker J, Lefeber C, Eerdekens M.",
    "Population pharmacokinetic modeling to facilitate dose selection of tapentadol",
    "in the pediatric population.",
    "J Pain Res. 2019;12:2835-2850.",
    "doi:10.2147/JPR.S208454",
    sep = " "
  )
  vignette <- "Watson_2019_tapentadol"

  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "ng/mL"
  )

  compartmentData <- list(
    depot = list(analyte = "tapentadol", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "tapentadol", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "The only covariate retained in the final model. Watson 2019 built body weight into",
        "the base model as the size descriptor on both CL/F and V/F using the continuous-covariate",
        "power function of the Methods, PTV = theta * (x / xref)^n, with xref set to the median of",
        "the covariate in the study population (45 kg; cohort mean (SD) 43 (19.7) kg, range",
        "12.7-80 kg per Table 1). Both exponents were estimated rather than fixed to the standard",
        "allometric 0.75 (CL) and 1 (V): the paper reports that estimating them gave a significant",
        "OFV reduction and better fit (Results, 'The final population model...'), yielding 0.638 for",
        "CL/F and 0.847 for V/F (Table 2). Because weight was already in the base model it was not",
        "part of the stepwise covariate screen. Treated as time-fixed here; the trials were",
        "single-dose so no within-subject weight change is modelled.",
        sep = " "
      ),
      source_name = "WT (Watson 2019 Table 2, rows 'Exponent CL-WT' and 'Exponent V-WT'); 'Bodyweight (kg)' in Table 1"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Subject age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Screened on both CL/F and V/F in the stepwise covariate model and not retained",
        "(Results, 'During the first iteration of the forward inclusion...'; Discussion,",
        "'Population PK model'). Watson 2019 interprets the absence of an age effect beyond",
        "weight as evidence that glucuronidation (the major tapentadol elimination route) is",
        "essentially mature by 2 years, so no maturation function is required above that age",
        "(Discussion, 'Maturation'). Cohort mean (SD) 11 (4.7) years, range 2 to <18 years",
        "(Table 1). The model is therefore not informative below 2 years of age.",
        sep = " "
      ),
      source_name = "Age (Watson 2019 Table 1)"
    ),
    SEXF = list(
      description = "Female sex indicator; 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = paste(
        "Screened on both CL/F and V/F using the categorical-covariate relationship of the",
        "Methods, PTV = theta_ref * (1 + sum(theta_j * I)), with male as the stated reference",
        "category; not retained (Results; Discussion, 'Population PK model'). The cohort was",
        "53.3% female / 46.7% male (Results, first paragraph; per-group percentages in Table 1).",
        sep = " "
      ),
      source_name = "% females (Watson 2019 Table 1); 'sex' in the covariate-analysis list"
    ),
    CRCL = list(
      description = "Creatinine clearance",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Screened on CL/F only (Methods, 'The effect of creatinine clearance (CRCL) on CL/F was",
        "also investigated'). It produced a significant OFV drop in the first forward-inclusion",
        "iteration but was not the most significant relationship and did not survive backward",
        "elimination, so it is absent from the final model (Results; Discussion, 'Population PK",
        "model'). No point estimate or baseline distribution is reported, so no value is carried",
        "here. Units not stated by the paper; recorded as the register default.",
        sep = " "
      ),
      source_name = "CRCL (Watson 2019 Methods, 'Covariate analysis')"
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Screened on both CL/F and V/F. AST on V/F was carried into the final forward-inclusion",
        "model but was then removed during backward deletion for failing the 1% significance",
        "criterion, so it does not appear in the final model (Results). No point estimate is",
        "reported for the discarded effect. Units not stated by the paper; recorded as the",
        "register canonical.",
        sep = " "
      ),
      source_name = "AST (Watson 2019 Methods, 'Covariate analysis')"
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Screened on both CL/F and V/F. ALT on CL/F gave the single most significant",
        "forward-inclusion step (delta OFV = 6.16) and was carried forward, but it was removed",
        "during backward deletion and is absent from the final model (Results). No point estimate",
        "is reported for the discarded effect. Units not stated by the paper; recorded as the",
        "register canonical.",
        sep = " "
      ),
      source_name = "ALT (Watson 2019 Methods, 'Covariate analysis'; Results, 'delta OFV =6.16')"
    ),
    ALP = list(
      description = "Alkaline phosphatase",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Screened on both CL/F and V/F; significant on CL/F and on V/F in the first",
        "forward-inclusion iteration but not retained in the final model (Results). No point",
        "estimate is reported. Units not stated by the paper; recorded as the register canonical.",
        sep = " "
      ),
      source_name = "ALP (Watson 2019 Methods, 'Covariate analysis')"
    ),
    TBILI = list(
      description = "Total bilirubin",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Screened on both CL/F and V/F and not retained; it is listed among the covariates",
        "with 'no relevant impact on the PK of tapentadol' (Methods, 'Covariate analysis';",
        "Discussion, 'Population PK model'). No point estimate or baseline distribution is",
        "reported. Units not stated by the paper; recorded as the register canonical.",
        sep = " "
      ),
      source_name = "bilirubin (Watson 2019 Methods, 'Covariate analysis')"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 92L,
    n_studies = 2L,
    n_observations = "424 tapentadol serum concentrations used for model building and evaluation, from 462 samples drawn in 109 treated patients (Watson 2019 Methods, 'Data set')",
    age_range = "2 to <18 years; group medians (range) 15 (12-17), 9 (6-11) and 3 (2-5) years; cohort mean (SD) 11 (4.7) years",
    age_median = "11 years (mean); stratified as 12 to <18 years (n=44), 6 to <12 years (n=33) and 2 to <6 years (n=15)",
    weight_range = "12.7-80 kg; group medians (range) 60 (41-80), 29.5 (20.2-58) and 16.3 (12.7-19.5) kg; cohort mean (SD) 43 (19.7) kg",
    weight_median = "45 kg (the covariate reference weight used for CL/F and V/F in Table 2)",
    sex_female_pct = 53.3,
    disease_state = "Moderate-to-severe acute postsurgical pain",
    dose_range = "Single tapentadol oral solution dose of 1.0 mg/kg body weight; the 4 mg/mL strength for patients <20 kg and the 20 mg/mL strength for patients >=20 kg, with total dose capped at 75 mg",
    regions = "Canada, Spain and the USA (NCT01134536); single-country sites for NCT01729728 (Quorum Review IRB, Seattle, USA)",
    notes = paste(
      "Pooled from two open-label single-dose phase 2 PK trials, NCT01729728 (n=56) and",
      "NCT01134536 (n=36); baseline demographics by trial and age stratum are in Table 1.",
      "Patients were recruited in a staggered fashion by descending age group. Of 109 treated",
      "patients, 17 (38 samples, 8%) were excluded for vomiting within 3 hrs of intake or for an",
      "incomplete dose. About 3% of samples were below the 0.2 ng/mL LOQ and were excluded rather",
      "than imputed or replaced. Body weight and age were highly correlated (r=0.92), which is why",
      "weight alone suffices as the size descriptor. Sampling intensity differed by age group:",
      "8 timed samples for 12 to <18 years, 4 windowed samples for 6 to <12 years, 2 windowed",
      "samples for 3-5 years and 4 timed samples for 2-year-olds. The tapentadol-O-glucuronide",
      "metabolite was assayed but is not modelled here; the paper models only the active moiety",
      "because the glucuronide is not analgesically active. Estimation used NONMEM 7.2 FOCE with",
      "interaction.",
      sep = " "
    )
  )

  ini({
    # ========================================================================
    # Structural parameters: one-compartment model with first-order absorption
    # from a depot compartment delayed by an absorption lag time, and
    # first-order elimination (Watson 2019 Results, 'The final population model
    # was best described as a 1-compartment model with linear oral absorption
    # from a dose compartment into the central compartment, with a lag time and
    # linear elimination'). CL/F and V/F are apparent: only the oral solution
    # was studied and no absolute bioavailability was estimated, so no separate
    # bioavailability term appears. Estimates are for the 45 kg reference
    # subject (Table 2 caption).
    # ========================================================================
    lka <- log(2.03); label("Absorption rate constant Ka (1/h)") # Watson 2019 Table 2: Ka = 2.03 1/h (RSE 16.5%; NONMEM 95% CI 1.373-2.687; bootstrap 1.599-3.263)
    ltlag <- log(0.247); label("Absorption lag time TLAG (h)") # Watson 2019 Table 2: TLAG = 0.247 h (RSE 0.7%; NONMEM 95% CI 0.243-0.251; bootstrap 0.245-0.273)
    lcl <- log(170); label("Apparent clearance CL/F at 45 kg (L/h)") # Watson 2019 Table 2: CL/F = 170 L/h (RSE 3.3%; NONMEM 95% CI 159.06-180.94; bootstrap 162.08-182.94)
    lvc <- log(685); label("Apparent volume of distribution V/F at 45 kg (L)") # Watson 2019 Table 2: V/F = 685 L (RSE 4.5%; NONMEM 95% CI 624.83-745.17; bootstrap 653.55-777.96)

    # ---- Body-weight power (allometric) exponents ---------------------------
    # Continuous-covariate power function from the Methods,
    #     PTV = theta * (x / xref)^n,
    # with xref the population median (45 kg). Both exponents were ESTIMATED,
    # not fixed to the canonical 0.75 / 1: Watson 2019 Results reports that
    # 'Estimated allometric scaling on CL/F and V/F was also found to give a
    # significant reduction in OFV and improvement in fit compared to fixing
    # the values to the standard 0.75 and 1', so neither is wrapped in fixed().
    e_wt_cl <- 0.638; label("Body-weight power exponent on CL/F (unitless)") # Watson 2019 Table 2: Exponent CL-WT = 0.638 (RSE 11.1%; NONMEM 95% CI 0.499-0.777; bootstrap 0.515-0.766)
    e_wt_vc <- 0.847; label("Body-weight power exponent on V/F (unitless)") # Watson 2019 Table 2: Exponent V-WT = 0.847 (RSE 10.2%; NONMEM 95% CI 0.678-1.016; bootstrap 0.718-1.029)

    # ========================================================================
    # Inter-individual variability
    # Watson 2019 Methods gives an exponential IIV model, Pi = Ptv * exp(eta_i),
    # with eta_i ~ N(0, omega^2). Table 2 reports the NONMEM $OMEGA entries
    # directly on the variance / covariance scale ('IIV CL/F (omega^2)',
    # 'Cov CL/F-V/F'), so the values below are used verbatim with no
    # log(1 + CV^2) conversion.
    #
    # Reported as a full correlation matrix between CL/F, V/F and Ka (Results,
    # 'Interindividual variability was best described using a full correlation
    # matrix between the parameters'), hence a single 3x3 block. The block is
    # positive definite (eigenvalues 1.993, 0.0678, 0.00157) and reproduces
    # every summary statistic the paper quotes in prose:
    #   sqrt(diag)  = 21.9% / 15.5% / 141.1%  vs reported 21.8% / 15.5% / 141.1%
    #   cov2cor     = 0.884 / 0.029 / -0.330  vs reported 0.88 / 0.03 / -0.33
    # ========================================================================
    etalcl + etalvc + etalka ~ c(
      0.048, # Watson 2019 Table 2, row 'IIV CL/F (omega^2)' = 0.048 (RSE 32.1%; bootstrap 0.025-0.088); sqrt = 21.9% IIV
      0.03, 0.024, # Watson 2019 Table 2, rows 'Cov CL/F-V/F' = 0.03 (RSE 46.1%) and 'IIV V/F (omega^2)' = 0.024 (RSE 61.5%); sqrt = 15.5% IIV; correlation 0.88
      0.009, -0.072, 1.99 # Watson 2019 Table 2, rows 'Cov CL/F-Ka' = 0.009 (RSE 614.5%), 'Cov V/F-Ka' = -0.072 (RSE 93.6%) and 'IIV Ka (omega^2)' = 1.99 (RSE 32.2%); sqrt = 141.1% IIV; correlations 0.03 and -0.33
    )

    # ========================================================================
    # Residual error
    # Watson 2019 Methods equation: Co,ij = Cp,ij * (1 + eps1,ij) + eps2,ij,
    # i.e. combined proportional (subscript 1) plus additive (subscript 2) --
    # 'the residual error was best described with a proportional and additive
    # error model' (Results).
    #
    # Both Table 2 entries are standard deviations, not variances: the table's
    # abbreviation list defines 'sigma, standard deviation' (and reserves
    # 'omega^2, variance' for the IIV rows), the additive term is given in
    # ng/mL rather than (ng/mL)^2, and the Results text reads the proportional
    # term straight off as a percentage -- 'The residual additive and
    # proportional error variabilities were 0.181 ng/mL and 32.9%'. A variance
    # reading would make the proportional error sqrt(0.329) = 57.4%, which
    # contradicts that sentence. So both are used as-is with no sqrt().
    # ========================================================================
    addSd <- 0.181; label("Additive residual error SD (ng/mL)") # Watson 2019 Table 2: Additive error = 0.181 ng/mL (RSE 39.1%; NONMEM 95% CI 0.042-0.32; bootstrap 0.036-0.415)
    propSd <- 0.329; label("Proportional residual error SD (fraction)") # Watson 2019 Table 2: Proportional error (sigma) = 0.329 (RSE 8.7%; NONMEM 95% CI 0.273-0.385; bootstrap 0.269-0.365); Results quote it as 32.9%
  })

  model({
    # ---- Individual parameters ---------------------------------------------
    # Body weight enters as the Methods power function normalized to the 45 kg
    # population median (Watson 2019 Table 2 caption: 'Estimates of CL/F and
    # V/F relate to a reference weight of 45 kg').
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) * (WT / 45)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 45)^e_wt_vc

    kel <- cl / vc

    # ---- ODE system --------------------------------------------------------
    # Dose records target the depot compartment. CL/F and V/F already absorb
    # the unknown oral bioavailability, so f(depot) is deliberately left at its
    # default of 1 rather than carrying a separate lfdepot term.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    alag(depot) <- exp(ltlag)

    # Amounts are in mg and vc is in L, so central / vc is mg/L; multiply by
    # 1000 to report ng/mL, the unit in which the concentrations were assayed
    # (LOQ 0.2 ng/mL) and in which the additive residual error was estimated.
    Cc <- (central / vc) * 1000

    Cc ~ add(addSd) + prop(propSd)
  })
}
