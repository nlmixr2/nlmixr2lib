Csajka_2004_indinavir <- function() {
  description <- "One-compartment first-order-absorption population PK model for oral indinavir 800 mg three-times-daily (alone) or 800 mg twice-daily with low-dose ritonavir in HIV-infected adults; concomitant ritonavir, sex, and body weight enter apparent oral clearance as multiplicative covariate effects (Csajka 2004)."
  reference <- "Csajka C, Marzolini C, Fattinger K, Decosterd LA, Telenti A, Biollaz J, Buclin T. Population pharmacokinetics of indinavir in patients infected with human immunodeficiency virus. Antimicrob Agents Chemother. 2004;48(9):3226-3232. doi:10.1128/aac.48.9.3226-3232.2004"
  vignette <- "Csajka_2004_indinavir"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Enters apparent oral clearance via a linear deviation from the 70 kg reference weight: (1 + e_wt_cl * (WT - 70) / 70). Cohort median 66.8 kg (range 41-116 kg) per Csajka 2004 Table 1.",
      source_name        = "BW"
    ),
    SEXF = list(
      description        = "Biological sex indicator (1 = female, 0 = male).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 (female); Csajka 2004's typical-value CL/F (32.4 L/h) is the female reference and the +30% male effect is applied to males.",
      notes              = "Csajka 2004 Table 3 footnote a encodes sex as a male-indicator (sex = 1 if male) and reports the male effect as theta_male = 0.30. The canonical SEXF (1 = female, 0 = male) inverts the values, so the effect is applied in model() as (1 + e_sex_cl * (1 - SEXF)), preserving Csajka 2004's female-reference CL/F = 32.4 L/h (32.4 * 1.30 = 42.1 L/h for males, matching the paper's reported 42.0 L/h within rounding).",
      source_name        = "sex"
    ),
    CONMED_RTV = list(
      description        = "Concomitant ritonavir (low-dose CYP3A4-inhibitor / PK-booster) coadministration indicator.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no ritonavir; indinavir 800 mg three-times-daily monotherapy).",
      notes              = "1 = subject receives ritonavir 100 mg twice-daily as a kinetic booster with indinavir 800 mg twice-daily; 0 = no ritonavir (indinavir 800 mg three-times-daily alone). Csajka 2004 Methods Study population: 177 of 239 patients received ritonavir. The covariate effect is applied multiplicatively to apparent oral clearance via (1 + e_rtv_cl * CONMED_RTV); e_rtv_cl = -0.63 (Table 3), so ritonavir reduces indinavir CL/F by 63% relative to the no-RTV reference.",
      source_name        = "PI"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 239L,
    n_studies       = 1L,
    n_observations  = 569L,
    n_visits        = 490L,
    age_range       = "16.3-73.4 years",
    age_median      = "40.1 years",
    weight_range    = "41-116 kg",
    weight_median   = "66.8 kg",
    height_range    = "150-194 cm",
    height_median   = "172 cm",
    sex_female_pct  = 29.3,
    race_ethnicity  = c(Caucasian = 92, Black = 5, Hispanic = 3, Asian = 0.8),
    disease_state   = "Adults infected with HIV-1 on combination antiretroviral therapy. CD4 12-1491 cells/mm^3 (median 433), HIV viral load 1-750000 copies/mm^3 (median 400) (Csajka 2004 Table 1).",
    dose_range      = "Oral indinavir 800 mg three-times-daily (62 of 239 patients, no ritonavir) or 800 mg twice-daily with ritonavir 100 mg twice-daily (177 of 239 patients). 21 of 239 had been already dose-reduced to 400 or 600 mg twice-daily for side effects prior to enrollment.",
    regions         = "Switzerland (University Hospital, Lausanne; University Hospital, Zurich).",
    notes           = "Single-centre observational cohort followed over a 40-month period. Sparse PK sampling (median 2 visits per patient, range 1-8) plus 7 patients with full 0-8 h time-course profiles. Indinavir quantified by reverse-phase HPLC (LLOQ 250 ug/L, linear to 10000 ug/L). NONMEM V / NM-TRAN II, first-order conditional estimation. Final population PK parameter estimates from Csajka 2004 Table 3; covariate-screening tests from Table 2."
  )

  ini({
    # Structural parameters from Csajka 2004 Table 3 (final-model column).
    # Reference subject: 70 kg female on indinavir 800 mg three-times-daily
    # without concomitant ritonavir.
    lcl <- log(32.4)
    label("Apparent oral clearance, female reference, no ritonavir (CL/F, L/h)")  # Csajka 2004 Table 3: CL_female = 32.4 L/h (95% CI 27.8-37.2)
    lvc <- log(65.7)
    label("Apparent volume of distribution (V/F, L)")                             # Csajka 2004 Table 3: V/F = 65.7 L (95% CI 55.9-76.0)
    lka <- log(1.0)
    label("First-order absorption rate constant (ka, 1/h)")                       # Csajka 2004 Table 3: Ka = 1.0 1/h (95% CI 0.80-1.41)

    # Bioavailability fixed at 1: indinavir was only given orally, so CL
    # and V are apparent (CL/F, V/F). Csajka 2004 Table 3 footnote h: "F
    # set to 1 because intravenous drug administration not available,
    # hence no 95% CI or % SE evaluable."
    lfdepot <- fixed(log(1.0))
    label("Oral bioavailability (FIXED to 1; no IV reference data)")              # Csajka 2004 Table 3 footnote h

    # Covariate effects on CL/F. Final-model multiplicative form (Csajka
    # 2004 Table 3 footnote a):
    #   CL/F = CL_female * (1 + theta_male) * (1 + theta_ritonavir)
    #                    * [1 + theta_BW * (BW - 70) / 70]
    #
    # Sex effect: Csajka codes sex(male) = 1 with female as reference;
    # the canonical SEXF (1 = female) inverts the values, so the effect
    # is applied via (1 - SEXF) in model() to preserve Csajka's female-
    # reference CL/F.
    e_sex_cl <- 0.30
    label("Fractional change in CL/F for male sex (multiplicative; applied via (1 - SEXF); unitless)")  # Csajka 2004 Table 3: theta_male = 0.30 (95% CI 0.14-0.49)

    # Concomitant ritonavir (CYP3A4 inhibitor / PK booster) multiplicative
    # decrease in CL/F when CONMED_RTV = 1.
    e_rtv_cl <- -0.63
    label("Fractional change in CL/F for concomitant ritonavir (multiplicative; unitless)")             # Csajka 2004 Table 3: theta_ritonavir = -0.63 (95% CI -0.60 to -0.64)

    # Body weight: linear deviation from the 70 kg reference weight,
    # applied multiplicatively to CL/F.
    e_wt_cl <- 0.16
    label("Proportionality coefficient relating CL/F to relative deviation of BW from 70 kg (unitless)") # Csajka 2004 Table 3: theta_BW = 0.16 (95% CI 0.12-0.35)

    # IIV. Csajka 2004 reports CV% on the natural-parameter scale
    # (exponential / log-normal IIV; Methods, "Model-based pharmacokinetic
    # analysis" paragraph). Convert to log-normal variance via
    # omega^2 = log(1 + CV^2), with CV as a fraction. The paper assigns
    # interpatient variability only to CL/F and ka (Results: "no
    # interpatient variability was an asset to V or F values"; Table 3
    # final-model interindividual-variability columns).
    # CL/F: CV = 39.0% -> omega^2 = log(1 + 0.39^2) = 0.14160
    # ka:   CV = 67.0% -> omega^2 = log(1 + 0.67^2) = 0.37067
    etalcl ~ 0.14160  # Csajka 2004 Table 3 final-model IIV CL/F = 39%
    etalka ~ 0.37067  # Csajka 2004 Table 3 final-model IIV ka  = 67%

    # Residual error: combined proportional + additive (Csajka 2004
    # Methods: "combined exponential and additive model was assigned to
    # the intrapatient (residual) variability"; Table 3 final-model rows
    # sigma(CV%) and sigma(SD; ug/L)). Concentration units in this model
    # are mg/L, so the paper's 670 ug/L additive SD becomes 0.670 mg/L.
    propSd <- 0.41
    label("Proportional residual error (fraction)")                               # Csajka 2004 Table 3: sigma(CV%) = 41%
    addSd <- 0.670
    label("Additive residual error (mg/L)")                                       # Csajka 2004 Table 3: sigma(SD) = 670 ug/L = 0.670 mg/L
  })

  model({
    # Sex effect: Csajka 2004 codes sex(male) = 1 with female as the
    # reference, so (1 - SEXF) reproduces the paper's male = 1 column
    # while SEXF (1 = female) remains the canonical storage convention.
    sex_male <- 1 - SEXF

    # Apparent oral clearance with covariate effects (Csajka 2004 final-
    # model equation, Table 3 footnote a):
    #   CL/F = exp(lcl + etalcl)
    #          * (1 + theta_male       * sex_male)
    #          * (1 + theta_ritonavir  * CONMED_RTV)
    #          * (1 + theta_BW         * (WT - 70) / 70)
    cl <- exp(lcl + etalcl) *
      (1 + e_sex_cl * sex_male) *
      (1 + e_rtv_cl * CONMED_RTV) *
      (1 + e_wt_cl * (WT - 70) / 70)

    # Apparent volume of distribution (no covariates and no IIV in the
    # final model).
    vc <- exp(lvc)

    # First-order absorption rate constant with IIV.
    ka <- exp(lka + etalka)

    # Bioavailability anchor (fixed to 1; oral-only data).
    fdepot <- exp(lfdepot)

    # Elimination micro-constant.
    kel <- cl / vc

    # ODE system: one-compartment open model with first-order absorption
    # from depot.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Bioavailability on the depot compartment (anchored to 1).
    f(depot) <- fdepot

    # Observation: oral CL/F and V/F absorb F into the apparent terms.
    # Dose in mg and V in L give Cc in mg/L (the paper reports
    # concentrations in ug/L; 1 mg/L = 1000 ug/L).
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
