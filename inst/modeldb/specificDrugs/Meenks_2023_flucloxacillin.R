Meenks_2023_flucloxacillin <- function() {
  description <- "One-compartment integrated total/unbound population PK model for intravenous flucloxacillin in critically ill adults (Meenks 2023), developed by an iterative two-stage Bayesian procedure in EDSIM++/KINPOP++. The disposition is carried on the UNBOUND drug: unbound clearance is proportional to CKD-EPI creatinine clearance (CL = fr * CLcr) and the unbound volume scales linearly with body weight. Protein binding is linear; the measured TOTAL concentration is Cu / fu, with the unbound fraction a power function of serum albumin. Total and unbound concentrations are both observed endpoints, each with proportional residual error. Inter-individual variability on fr, V, fu and the albumin exponent."
  reference <- "Meenks SD, Punt N, le Noble JLML, Foudraine NA, Neef K, Janssen PKC. Target attainment and population pharmacokinetics of flucloxacillin in critically ill patients: a multicenter study. Crit Care. 2023;27(1):82. doi:10.1186/s13054-023-04353-5. PMCID PMC9982780. Covariate equations from the Results section 'Population PK model development'; parameter estimates from main-text Table 2 (final-model column); free-PK-mode structure, fr definition and allometric conventions from Additional file 1 (Supplementary e-Appendix, 'Population pharmacokinetic model development')."
  vignette <- "Meenks_2023_flucloxacillin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Linear (exponent 1) scaling of the unbound volume to 70 kg. Table 2 reports V in 'L/70 kg' and Additional file 1 states 'Regular allometric scaling was applied for the volumes of distribution, calculated in L/70 kg body weight'. Cohort median 76 kg (IQR 64-85) at ICU admission (Table 1).",
      source_name = "body weight"
    ),
    ALB = list(
      description = "Serum albumin concentration",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Reference 21.2 g/L, the literal denominator of the printed fu equation fu = 0.217 * (Alb / 21.2)^-ex (Results, 'Population PK model development'; 'The population estimate of serum albumin concentration was 21.2 g/L'). Measured by bromocresol purple colorimetric assay (Methods, Data collection). Mean 23.3 g/L (SD 8.2) at flucloxacillin sampling, 25.8 g/L (SD 8.5) at ICU admission (Table 1; Table S1). Enters on the unbound fraction only, so albumin changes TOTAL but not unbound concentrations.",
      source_name = "Alb"
    ),
    CRCL = list(
      description = "Estimated glomerular filtration rate by the CKD-EPI creatinine equation",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters linearly as CL (unbound) = fr * CLcr with no reference value (Results equation CLtot = fr * CLcr; Table 2 row 'fr, CL/CLcr'). The paper states the creatinine clearance was 'calculated using the CKD-EPI equation' and reports it in mL/min/1.73 m^2 (Table 1, Table S1); no de-normalisation by body surface area is described, so the BSA-normalised value is used as printed. model() converts mL/min to L/h by * 60 / 1000 so the dimensionless ratio fr returns L/h. Mean 66 (SD 42) at sampling, 68 (SD 42) at ICU admission.",
      source_name = "CLcr"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Used as a patient characteristic in the KINPOP++ fit (Additional file 1) but not a covariate on any parameter. Cohort median 69 years (IQR 54-76)."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = male",
      notes = "Not evaluated: 'No categorical covariate data were evaluated' (Additional file 1, Covariate models). 58% of the cohort was male."
    ),
    HT = list(
      description = "Body height",
      units = "cm",
      type = "continuous",
      reference_category = NULL,
      notes = "Recorded as a patient characteristic (Additional file 1) but not a covariate in the final model. Mean 172 cm (SD 10)."
    )
  )

  compartmentData <- list(
    central = list(
      analyte = "flucloxacillin (unbound)",
      units = "mg",
      specimen = "serum",
      verified = TRUE,
      notes = "Holds the UNBOUND flucloxacillin amount; dividing by V gives the unbound serum concentration Cu. Additional file 1: the model ran in EDSIM++ 'Free PK mode', in which 'total concentrations were calculated by ... C tot = C u / fu', and 'The PK parameters for unbound flucloxacillin were estimated from the unbound and total flucloxacillin serum concentrations'; fr is defined as 'the unbound renal clearance of flucloxacillin divided by the creatinine clearance' (Additional file 1; Table 2 footnote)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 31,
    n_studies = 1,
    age_range = "median 69 years (IQR 54-76)",
    age_median = "69 years",
    weight_range = "median 76 kg (IQR 64-85) at ICU admission; BMI median 25.3 kg/m^2 (IQR 23.4-28.9)",
    weight_median = "76 kg",
    sex_female_pct = 42,
    race_ethnicity = "Not reported (two Dutch ICUs)",
    disease_state = "Critically ill adults receiving intravenous flucloxacillin for (suspected) methicillin-sensitive Staphylococcus aureus infection; mean APACHE II score 21 (SD 10); 55% mechanically ventilated; 48% sepsis. Patients on renal replacement therapy or with liver cirrhosis were excluded.",
    renal_function = "CKD-EPI creatinine clearance mean 66 (SD 42) mL/min/1.73 m^2 at flucloxacillin sampling; 13% below 15 and 3% at or above 150 (Table S1).",
    albumin = "Serum albumin mean 23.3 g/L (SD 8.2) at sampling; 29% 10-19 g/L, 44% 20-34 g/L, 23% 35-55 g/L (Table S1).",
    dose_range = "4-12 g/day: continuous infusion of 6, 9 or 12 g per 24 h (74% of patients) or 15-30 min intermittent infusions of 1 g q6h, 1 g q4h or 2 g q4h (26%)",
    regions = "The Netherlands (VieCuri Medical Center Noord-Limburg, Venlo, and Maastricht University Medical Center+)",
    notes = "Prospective observational multicenter study, May 2017 - October 2019. 163 residual arterial blood samples (median 4 per patient, IQR 2-10) each analysed for total and unbound flucloxacillin by UPLC-MS/MS with ultrafiltration; total 1.3-668 mg/L, unbound 0.6-137 mg/L, unbound fraction median 22% (IQR 20-27, range 6-73%). Demographics in Table 1 and Additional file 1 Table S1. Population parameters by an iterative two-stage Bayesian procedure (EDSIM++ 2.04, KINPOP++ module); 1000-sample bootstrap. Final-model bias -27.1% and precision 53.3%."
  )

  ini({
    # ------------------------------------------------------------------
    # Table 2, 'Final model' column (Estimate, mean; RSE %). The EDSIM++
    # KINPOP++ iterative two-stage output reports each parameter as a
    # population mean with an inter-individual %CV, and Additional file 1
    # states 'Log-normal distribution was assumed for inter-individual PK
    # variance'. The mean is carried as the typical value and the CV is
    # converted to a log-scale variance by omega^2 = log(1 + CV^2), the
    # same convention as the MwPharm-built Daskapan_2019_darunavir and
    # Geers_2022_clozapine models. The vignette discusses the alternative
    # reading in which the printed mean is the arithmetic mean of the
    # log-normal.
    # ------------------------------------------------------------------

    # Unbound clearance ratio. The final model has no non-renal arm:
    # CL (unbound) = fr * CLcr (Results equation 'CLtot = fr * CLcr').
    # fr is dimensionless ('the unbound renal clearance of flucloxacillin
    # divided by the creatinine clearance', Table 2 footnote), so model()
    # converts CLcr from mL/min to L/h before multiplying.
    e_crcl_cl <- 19
    label("Ratio fr of unbound flucloxacillin clearance to CKD-EPI creatinine clearance (unitless)") # Table 2 final model fr (CL/CLcr) = 19 [RSE 13%]; bootstrap median 19 (95% CI 12-27)

    # Unbound volume of distribution, per 70 kg total body weight.
    lvc <- log(330)
    label("Unbound volume of distribution V for a 70 kg patient (L)") # Table 2 final model V = 330 L/70 kg [RSE 15%]; bootstrap median 321 (95% CI 204-476)

    # Unbound fraction at the reference albumin 21.2 g/L. The printed
    # covariate equation is fu = 0.217 * (Alb / 21.2)^-ex and the Results
    # restate 'The population estimate of the unbound fraction was 0.217'.
    # Table 2 lists fu = 25% (RSE 4.7%) for the final model; the equation
    # value is used here and the vignette discusses the difference.
    lfu <- log(0.217)
    label("Unbound fraction fu at serum albumin 21.2 g/L (fraction)") # Results 'Population PK model development' equation fu = 0.217 * (Alb/21.2)^-ex and text 'population estimate of the unbound fraction was 0.217'; Table 2 final model fu = 25% [RSE 4.7%], bootstrap 25 (21-30)

    # Albumin exponent on fu. Printed as (Alb / 21.2)^-ex with ex = 0.67,
    # so the signed exponent stored here is -0.67 (fu falls as albumin
    # rises). ex carries its own log-normal IIV in Table 2, applied
    # multiplicatively to the signed exponent in model().
    e_alb_fu <- -0.67
    label("Power exponent of (ALB / 21.2) on fu (unitless; printed as -ex)") # Table 2 final model ex = 0.67 [RSE 16%]; bootstrap median 0.63 (95% CI 0.29-1.14); sign from the printed equation (Alb/21.2)^-ex

    # ------------------------------------------------------------------
    # IIV (Table 2 final model, %CV; no covariance reported -- Methods,
    # Monte Carlo dosing simulations: 'Covariance between model
    # parameters was assumed to be absent'). omega^2 = log(1 + CV^2):
    #   fr 71% -> log(1 + 0.71^2) = 0.408195
    #   V  84% -> log(1 + 0.84^2) = 0.533917
    #   fu 26% -> log(1 + 0.26^2) = 0.065413
    #   ex 88% -> log(1 + 0.88^2) = 0.573462
    # ------------------------------------------------------------------
    etae_crcl_cl ~ 0.408195 # Table 2 final model 'IIV fr' = 71% (bootstrap 70%)
    etalvc ~ 0.533917 # Table 2 final model 'IIV V' = 84% (bootstrap 80%)
    etalfu ~ 0.065413 # Table 2 final model 'IIV fu' = 26% (bootstrap 23%)
    etae_alb_fu ~ 0.573462 # Table 2 final model 'IIV ex' = 88% (bootstrap 78%)

    # ------------------------------------------------------------------
    # Residual error: proportional on each endpoint (Additional file 1:
    # 'Residual variability ... was modelled for the proportional error
    # for the unbound and total flucloxacillin concentrations').
    # Cc (total) takes the suffix-free name.
    # ------------------------------------------------------------------
    propSd <- 0.42
    label("Proportional residual SD on total concentration Cc (fraction)") # Table 2 final model 'Proportional error, total flucloxacillin' = 42% [RSE 5.2%]; bootstrap 42 (26-52)
    propSd_Cu <- 0.35
    label("Proportional residual SD on unbound concentration Cu (fraction)") # Table 2 final model 'Proportional error, unbound flucloxacillin' = 35% [RSE 5.2%]; bootstrap 35 (25-42)
  })

  model({
    # 1. Individual unbound clearance: CL = fr * CLcr (Results equation),
    #    with CLcr converted from mL/min/1.73 m^2 to L/h (* 60 / 1000).
    fr <- e_crcl_cl * exp(etae_crcl_cl)
    cl <- fr * CRCL * 60 / 1000

    # 2. Unbound volume, linear in total body weight per 70 kg.
    vc <- exp(lvc + etalvc) * (WT / 70)

    # 3. Unbound fraction: fu = 0.217 * (Alb / 21.2)^-ex, with log-normal
    #    IIV on both fu and the exponent ex.
    ex_alb <- e_alb_fu * exp(etae_alb_fu)
    fu <- exp(lfu + etalfu) * (ALB / 21.2)^ex_alb

    # 4. One-compartment ODE on the UNBOUND amount; intravenous doses
    #    (intermittent or continuous infusion) enter central directly.
    kel <- cl / vc
    d/dt(central) <- -kel * central

    # 5. Observations. Free PK mode: Cu = amount / V and the total serum
    #    concentration Ctot = Cu / fu (Additional file 1, linear binding).
    Cu <- central / vc
    Cc <- Cu / fu

    Cc ~ prop(propSd)
    Cu ~ prop(propSd_Cu)
  })
}
