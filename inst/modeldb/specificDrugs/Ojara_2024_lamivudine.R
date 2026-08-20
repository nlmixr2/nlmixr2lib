Ojara_2024_lamivudine <- function() {
  description <- "Lactation population PK model for oral lamivudine in breastfeeding Ugandan women living with HIV, fitted simultaneously to paired maternal plasma and breast-milk concentrations across a full 24-hour dosing interval. First-order absorption without lag time into a one-compartment plasma disposition model; plasma-to-breast-milk transfer is a partitioned effect compartment whose concentration state equilibrates at first-order rate ke0 towards ppc times the plasma concentration, where ppc is the estimated milk-to-plasma concentration ratio of 1.77. That structure reproduces the observed lag between the plasma and breast-milk concentration peaks (2 h vs 6 h). Correlated exponential interindividual variability on clearance and central volume, independent interindividual variability on the milk-to-plasma ratio, no interindividual variability on absorption, and separate proportional residual errors for plasma and for breast milk. Clearance and volume are apparent oral values. No covariates were retained in the final model."
  reference <- paste(
    "Ojara FW, Kawuma AN, Nakalema S, Kyohairwe I, Nakijoba R, Lamorde M,",
    "Pertinez H, Khoo S, Waitt C.",
    "Population pharmacokinetic modeling of paired plasma-breast milk",
    "lamivudine data for estimation of infant exposure in breastfeeding",
    "mother-infant pairs.",
    "CPT Pharmacometrics Syst Pharmacol. 2024;13(11):1978-1989.",
    "doi:10.1002/psp4.13274.",
    "Structural equations from Methods Equations 1-3; parameter estimates",
    "from Table 2.",
    sep = " "
  )
  vignette <- "Ojara_2024_lamivudine"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  # Issue #482. Doses are in mg and volumes in L, so central holds mg and the
  # derived Cc is mg/L = ug/mL; the paper reports observed concentrations in
  # ng/mL (multiply Cc and Cmilk by 1000). The `milk` state is an exception to
  # the "states hold amounts" rule: Ojara 2024 Equation 3 defines it as
  # "A3: concentration of lamivudine in the breast milk compartment", so no
  # milk volume enters the model at all (Discussion: the authors deliberately
  # replaced an earlier parameterisation that "fixed the breast milk volume to
  # a physiological value" with the volume-free effect-compartment concept).
  # verified = TRUE: analyte and specimen read off the Equation 1-3 variable
  # definitions and the Figure 1 model schematic.
  compartmentData <- list(
    depot   = list(analyte = "lamivudine", units = "mg",    specimen = "administration site", verified = TRUE),
    central = list(analyte = "lamivudine", units = "mg",    specimen = "plasma",              verified = TRUE),
    milk    = list(analyte = "lamivudine", units = "ug/mL", specimen = "milk",                verified = TRUE)
  )

  # The final model carries no covariates. Results, "Population
  # pharmacokinetic analysis": creatinine clearance on CL and weight on Vc
  # were statistically significant in the stepwise covariate analysis
  # (dOFV = 35.6, p < 0.01) but were imprecisely estimated (%RSE > 30%, Table
  # S1), "as such, the final plasma disposition model has no covariates
  # included". The screened-but-not-retained covariates are documented in
  # covariatesDataExcluded below.
  covariateData <- list()

  covariatesDataExcluded <- list(
    WT = list(
      description        = "Maternal body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened on both clearance and the central volume of distribution using linear, exponential and power functional forms with the covariate centred on the cohort median (Methods, 'Population pharmacokinetic analysis'). The weight effect on Vc reached the forward-inclusion and backward-elimination significance criteria but was imprecisely estimated (%RSE > 30%, Table S1) and was therefore not retained. Cohort median (range) 64.0 kg (50.0, 89.0) per Table 1.",
      source_name        = "Maternal weight"
    ),
    AGE = list(
      description        = "Maternal age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened on clearance only (Methods, 'Population pharmacokinetic analysis'); not significant and not retained. Cohort median (range) 30 years (19, 40) per Table 1.",
      source_name        = "Maternal age"
    ),
    CRCL = list(
      description        = "Creatinine clearance predicted from serum creatinine by the Cockcroft-Gault equation",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened on clearance (Methods, 'Population pharmacokinetic analysis'). Reached the stepwise significance criteria but was imprecisely estimated (%RSE > 30%, Table S1) and was therefore not retained; the Discussion attributes the imprecision to the small sample size and notes that earlier lamivudine publications did find a significant creatinine-clearance effect on clearance. Reported as mL/min (uncorrected for body surface area) because the Cockcroft-Gault equation returns an absolute clearance. Cohort median (range) 134.6 mL/min (89.2, 184.7) per Table 1.",
      source_name        = "CRCL"
    ),
    BMI = list(
      description        = "Maternal body mass index",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened on both clearance and the central volume of distribution (Methods, 'Population pharmacokinetic analysis'); not retained. Cohort median (range) 24.8 kg/m^2 (20.0, 30.5) per Table 1.",
      source_name        = "BMI"
    ),
    TPP = list(
      description        = "Time postpartum at the pharmacokinetic sampling visit",
      units              = "weeks",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened on the milk-to-plasma concentration ratio ppc: \"The postpartum day of PK sampling (visit day) was evaluated as a covariate on the milk-to-plasma ratio\" (Methods) and \"The postpartum day of PK sampling did not significantly affect the milk-to-plasma ratio\" (Results). The paper records this covariate in DAYS postpartum, not the canonical weeks: visit 1 median (range) 13.0 (7.00, 43.0) days and visit 2 73.5 (36.0, 84.0) days per Table 1 (the Results text gives the visit-2 median as 74.0 days). Divide by 7 to obtain the canonical unit. Not retained, so the choice of unit is documentation only.",
      source_name        = "Visit day"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 35,
    n_studies      = 1,
    age_range      = "19-40 years",
    age_median     = "30 years",
    weight_range   = "50.0-89.0 kg",
    weight_median  = "64.0 kg",
    sex_female_pct = 100,
    disease_state  = "Breastfeeding women living with HIV receiving first-line antiretroviral therapy, recruited antenatally and sampled postpartum. Lamivudine was given with nevirapine plus zidovudine (150 mg twice-daily group) or with efavirenz plus tenofovir disoproxil fumarate (300 mg once-daily group).",
    dose_range     = "Lamivudine 150 mg orally twice daily (12-hourly, n = 10) or 300 mg orally once daily (24-hourly, n = 25), at steady state.",
    regions        = "Uganda (Infectious Diseases Institute, Makerere University, and affiliated clinics in Kampala); enrolled 2016-2017.",
    renal_function = "Elevated creatinine clearance typical of the postpartum period: median (range) 134.6 mL/min (89.2, 184.7) by Cockcroft-Gault (Table 1).",
    infant_partner = "Each mother had one breastfed infant, freely breastfed throughout. Infant weight median (range) 3.60 kg (2.40, 5.58) at visit 1 and 5.17 kg (3.90, 7.13) at visit 2 (Table 1). Infant plasma concentrations were measured but were NOT part of the fitted model; they were used only to check the predicted infant steady-state exposure (Results, 'Estimation of infant breast milk exposure').",
    notes          = "Baseline demographics from Table 1. Mothers attended two of three possible visits at 1-2, 4-6 and 10-12 weeks postpartum; visit 1 median 13.0 days and visit 2 median 73.5 days postpartum. Sampling depended on the mother's habitual dosing time: 14 morning-dosing mothers gave plasma at 0, 1, 2, 4 and 8 h and breast milk at 0, 2, 4 and 8 h post directly-observed dose, while 21 evening-dosing mothers gave paired plasma and breast milk at 12, 16 and 20 h post self-reported dose. Breast milk was obtained by manual expression of 1 mL. The dataset contributing to this model comprises 248 maternal plasma and 256 breast-milk concentrations (plus 151 infant plasma concentrations that were not fitted). Concentrations were measured by a validated LC-MS/MS assay on dried blood and dried breast-milk spots with an LLOQ of 5 ng/mL for plasma and 16.6 ng/mL for breast milk. Steady state was assumed but not explicitly verified (Discussion, limitations). Estimation used FOCE with interaction in NONMEM 7.4.3; parameter precision came from a 1000-sample nonparametric bootstrap."
  )

  ini({
    # Structural parameters. Table 2 column "Estimates (%RSE)". CL and Vc are
    # APPARENT oral values: the study has no intravenous arm and Equations 1-2
    # deliver the whole dose to the depot with no bioavailability term, so F is
    # implicitly 1 and is not identifiable. The paper labels them "CL" and
    # "VC" without the /F qualifier.
    lcl <- log(19.4); label("Apparent oral clearance CL/F (L/h)")                            # Table 2: CL 19.4 L/h (4.0% RSE); bootstrap median 19.5 (17.5-21.5)
    lvc <- log(184); label("Apparent central volume of distribution Vc/F (L)")               # Table 2: VC 184 L (7.00% RSE); bootstrap median 180.1 (122.2-279.9)
    lka <- log(1.87); label("First-order absorption rate constant ka (1/h)")                 # Table 2: Ka 1.87 1/h (19.0% RSE); bootstrap median 1.85 (1.25-2.53)

    # Breast-milk partitioned effect compartment (Equation 3). ke0 sets how
    # fast the milk concentration equilibrates; ppc sets the level it
    # equilibrates to, i.e. the milk-to-plasma concentration ratio.
    lke0 <- log(0.245); label("Plasma-to-breast-milk equilibration rate constant Kcb (1/h)") # Table 2: Kcb 0.245 1/h (12.0% RSE); bootstrap median 0.251 (0.185-0.311)
    lppc <- log(1.77); label("Milk-to-plasma concentration ratio Rcb (unitless)")            # Table 2: Rcb 1.77 (6.00% RSE); bootstrap median 1.76 (1.58-1.95)

    # Interindividual variability. The Table 2 footnote defines the reported
    # %CV as "100 * sqrt(omega^2)", i.e. omega itself expressed as a
    # percentage, so variance = (%CV / 100)^2 -- NOT the log-normal
    # sqrt(exp(omega^2) - 1) transformation.
    #   IIV CL   20.9% -> 0.209^2 = 0.043681
    #   IIV Vc   76.4% -> 0.764^2 = 0.583696
    # The "IIV CL-VC, %CV 42.0%" row is a CORRELATION, not a covariance and
    # not a covariance-on-the-CV-scale. Reading 42.0% as 100 * sqrt(omega_12)
    # would give omega_12 = 0.1764 and a correlation of
    # 0.1764 / (0.209 * 0.764) = 1.105 > 1, which is not a positive
    # semi-definite covariance matrix; the same reading of the bootstrap row
    # (40.1 with CVs 20.2 and 71.5) also exceeds 1. Read as a correlation the
    # matrix is admissible, so covariance = 0.42 * 0.209 * 0.764 = 0.06706392.
    # Consistent with that row being the only one in Table 2 printed without a
    # %RSE, as a derived quantity would be.
    etalcl + etalvc ~ c(0.043681,
                        0.06706392, 0.583696)                                                # Table 2: IIV CL 20.9% CV (16.0% RSE), IIV Vc 76.4% CV (35.0% RSE), corr(CL,Vc) = 0.420
    etalppc ~ 0.025281                                                                       # Table 2: IIV Rcb 15.9% CV (19.0% RSE) -> 0.159^2

    # No IIV on ka. Results: "The available data did not support estimation of
    # interindividual variability on the absorption rate, resulting into a high
    # parameter shrinkage (i.e., 75%), as such the value was fixed to zero"
    # (21 of 35 women had no absorption-phase samples). Omitted rather than
    # written as `etalka ~ fixed(0)` because a zero-variance diagonal makes
    # OMEGA singular and breaks the Cholesky sampler used by rxSolve -- same
    # convention as Wattanakul_2024_primaquine_motherinfant.R.

    # Residual error. Separate proportional models for the two matrices; the
    # additive and combined models were also evaluated and rejected (Methods
    # and Results, "A proportional error model best characterized the residual
    # unexplained variability on plasma concentrations" / "in breast milk
    # concentrations").
    propSd <- 0.383; label("Proportional residual error, maternal plasma (fraction)")        # Table 2: RUV, PROP, PLASMA 38.3% CV (7.7% RSE); bootstrap median 38.7 (33.2-43.5)
    propSd_Cmilk <- 0.305; label("Proportional residual error, breast milk (fraction)")      # Table 2: RUV, PROP, BM 30.5% CV (8.6% RSE); bootstrap median 31.1 (25.3-36.7)
  })

  model({
    # 1. Individual parameters.
    ka <- exp(lka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    ke0 <- exp(lke0)
    ppc <- exp(lppc + etalppc)

    # 2. Micro-constant.
    kel <- cl / vc

    # 3. ODE system. Equations 1-2: the depot empties by first-order
    # absorption with no lag time (an absorption lag time gave
    # dOFV < 3.84 and transit-compartment absorption failed to estimate).
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    Cc <- central / vc

    # 4. Breast-milk partitioned effect compartment, Equation 3:
    #   dA3/dt = Kcb * (Rcb * A2 / VC - A3)
    # A3 is a CONCENTRATION, so no milk volume appears anywhere in the model.
    # At steady state over a complete dosing interval the integral of the
    # right-hand side is zero, so the average milk concentration equals ppc
    # times the average plasma concentration -- which is exactly the paper's
    # Equation 4 (ConcMilk = Rcb * ConcAve) used to derive the infant dose.
    d/dt(milk) <- ke0 * (ppc * Cc - milk)

    # 5. Observations and residual error.
    Cmilk <- milk
    Cc ~ prop(propSd)
    Cmilk ~ prop(propSd_Cmilk)
  })
}
