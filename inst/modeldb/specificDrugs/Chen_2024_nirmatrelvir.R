Chen_2024_nirmatrelvir <- function() {
  description <- "One-compartment first-order-absorption population PK model for oral nirmatrelvir given as the nirmatrelvir/ritonavir fixed combination in critically ill adults with COVID-19; apparent clearance depends on CKD-EPI creatinine clearance and on the co-administered ritonavir 12 h dosing-interval AUC through a doubly-centred power function, with correlated interindividual variability on CL/F and V/F (Chen 2024)."
  reference <- "Chen N, Yu X, Li L, Yang P, Dong R, Huang Y, Ling X, Shentu Q, Yu W, Jiang S. Target Attainment and Population Pharmacokinetics of Nirmatrelvir/Ritonavir in Critically Ill Adult Patients. Infect Drug Resist. 2024;17:4055-4065. doi:10.2147/IDR.S471918"
  vignette <- "Chen_2024_nirmatrelvir_ritonavir"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    depot   = list(analyte = "nirmatrelvir", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "nirmatrelvir", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description        = "Creatinine clearance estimated with the CKD-EPI equation",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Enters nirmatrelvir CL/F via the centred power form (CRCL / 80)^e_crcl_cl of Chen 2024 Equation 1, centred at 80 mL/min/1.73 m^2 (the rounded cohort mean of 78.5; Table 1). The paper labels the units inconsistently: the Table 3 footnote gives mL/min/1.73 m^2 while Table 1 and the Table 4 footnote give mL/min for the same CKD-EPI-derived quantity. The CKD-EPI equation returns a body-surface-area-normalised value, so mL/min/1.73 m^2 is used here, matching the Table 3 footnote that defines the covariate coefficient itself. Cohort value 78.5 +/- 34.4; the dosing-recommendation table spans 15 to >60. In combination with ritonavir the primary elimination pathway of nirmatrelvir shifts from liver to kidney, so CL/F falls as CRCL falls (Discussion).",
      source_name        = "CrCL"
    ),
    CONMED_RTV_AUC_12H = list(
      description        = "Ritonavir AUC over the 12 h (q12h) dosing interval for the co-administered 100 mg twice-daily ritonavir booster",
      units              = "mg*h/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Enters nirmatrelvir CL/F via the centred power form (CONMED_RTV_AUC_12H / 12.2)^e_rtv_auc_12h_cl of Chen 2024 Equation 1, centred at 12.2 mg*h/L. Computed per subject by Chen 2024 Equation 2 as AUC = DOSE / (CL/F) with DOSE = 100 mg (Figure 1 caption: 'AUC, area under curve of ritonavir base on 100mg') and CL/F the individual ritonavir apparent clearance from the companion ritonavir model in the same paper; see modellib('Chen_2024_ritonavir'). Simulated per-subject values are therefore obtained as 100 / cl_ritonavir. The Monte Carlo dosing simulations swept the 10th-90th percentiles of this covariate, and the Table 4 dose-recommendation grid spans 3.2 to 23.3 mg*h/L. Higher ritonavir exposure gives stronger CYP3A4/5 inhibition and hence lower nirmatrelvir CL/F, consistent with the negative exponent.",
      source_name        = "AUCRIT"
    )
  )

  # Screened during covariate model building but not retained in the final
  # nirmatrelvir model; the paper reports no point estimate for any of them.
  # Source: Chen 2024 Methods (Data Collection) for the screened list and
  # Results (Population PK Model Development) for the retained subset.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened but not retained; the Discussion states the number of young patients was insufficient to assess an age effect."
    ),
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened during stepwise covariate selection; not retained on CL/F or V/F."
    ),
    HT = list(
      description = "Height",
      units       = "cm",
      type        = "continuous",
      notes       = "Screened during stepwise covariate selection; not retained."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "categorical",
      notes       = "Screened during stepwise covariate selection; not retained. Cohort was 65% male (Table 1)."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Screened during stepwise covariate selection; not retained. Renal function entered the final model through CRCL instead."
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened during stepwise covariate selection; not retained."
    ),
    TPRO = list(
      description = "Total protein",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened during stepwise covariate selection; not retained."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened but not retained; the Discussion states the number of hepatic-impairment patients was insufficient to assess a liver-function effect."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened but not retained; see the ALT note."
    ),
    ALP = list(
      description = "Alkaline phosphatase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened during stepwise covariate selection; not retained."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Screened during stepwise covariate selection; not retained."
    ),
    APACHE_II = list(
      description = "Acute Physiology and Chronic Health Evaluation (APACHE) II score",
      units       = "(score)",
      type        = "continuous",
      notes       = "Screened during stepwise covariate selection; not retained. Cohort median 14 (IQR 11-21) per Table 1."
    ),
    SOFA = list(
      description = "Sequential Organ Failure Assessment score",
      units       = "(score)",
      type        = "continuous",
      notes       = "Screened during stepwise covariate selection; not retained. Cohort median 6 (IQR 3-9) per Table 1."
    ),
    CONMED_CYP3A4_INH = list(
      description = "Concomitant weak CYP3A4/5 inhibitor",
      units       = "(binary)",
      type        = "categorical",
      notes       = "Screened during stepwise covariate selection; not retained. 7 of 31 patients (22.6%) per Table 1. The Discussion notes that ritonavir, itself a strong CYP3A4/5 and P-gp inhibitor, may obscure the effect of other interacting co-medications."
    ),
    CONMED_CYP3A4_IND = list(
      description = "Concomitant weak CYP3A4/5 inducer",
      units       = "(binary)",
      type        = "categorical",
      notes       = "Screened during stepwise covariate selection; not retained. 10 of 31 patients (32.3%) per Table 1."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 31L,
    n_studies      = 1L,
    n_observations = 89L,
    age_median     = "69 years (IQR 63-77)",
    weight_median  = "60.6 kg (SD 12.8)",
    sex_female_pct = 35,
    disease_state  = "Critically ill adults (>=18 years) admitted to the intensive care unit with RT-PCR-confirmed SARS-CoV-2 infection and compatible pulmonary CT findings. Patients receiving continuous renal replacement therapy or extracorporeal membrane oxygenation were excluded. Median APACHE II score 14 (IQR 11-21); median SOFA score 6 (IQR 3-9). 81.2% of the cohort was older than 60 years (Discussion).",
    renal_function = "Creatinine clearance (CKD-EPI) 78.5 +/- 34.4 mL/min; serum creatinine 100.2 +/- 79.9 umol/L (Table 1).",
    dose_range     = "Oral nirmatrelvir/ritonavir 300 mg/100 mg twice daily in 27 of 31 patients (87%) and 150 mg/100 mg twice daily in 4 of 31 (13%). 29 of 31 (94%) received the tablets by nasal feeding and 2 (6%) orally, all for more than 5 days.",
    co_medication  = "Ritonavir 100 mg twice daily in every patient (the boosting component of the fixed combination). Concomitant CYP3A4/5 weak inhibitor in 7 of 31 (22.6%) and weak inducer in 10 of 31 (32.3%) per Table 1; named co-medications in the Discussion include voriconazole, dexamethasone, methylprednisolone and omeprazole.",
    regions        = "China (The First Affiliated Hospital, Zhejiang University School of Medicine, Hangzhou)",
    notes          = "Prospective observational study, January-June 2023. Two to three serial plasma samples per patient collected after the second dose; 89 plasma samples from 31 patients. Nirmatrelvir quantified by validated LC-MS/MS over 100.41-52991.51 ng/mL. Measured nirmatrelvir concentrations ranged 1214.31-21342.06 ng/mL. Model fit in Phoenix NLME 8.1 with FOCE; final OFV 95.0. Baseline demographics in Table 1; final parameter estimates and 1000-sample bootstrap in Table 3; the clearance covariate model is Equation 1. The efficacy target used in the paper's Monte Carlo simulations is the total in vitro EC90 of 292 ng/mL (0.292 mg/L), derived from the free EC90 of 90.5 ng/mL and approximately 70% protein binding; the safety ceiling is the human NOAEL of 79,700 ng/mL (79.7 mg/L)."
  )

  ini({
    # Structural parameters, Chen 2024 Table 3 (nirmatrelvir), "Final Model /
    # Estimate, Mean" column. The "CV, %" column of Table 3 is the relative
    # standard error of the estimate, not interindividual variability: the
    # Abstract reports these same estimates as "Mean (SD) ... 0.42 (0.10) h-1,
    # 36.5 (8.5) L, 3.6 (0.26) L/h", i.e. the estimate multiplied by the "CV, %"
    # value (0.42 * 24.9% = 0.105; 36.5 * 23.3% = 8.5; 3.6 * 7.1% = 0.256).
    lka <- log(0.42); label("Apparent first-order absorption rate constant (ka, 1/h)")  # Table 3: tvKa = 0.42 1/h (RSE 24.9%; bootstrap median 0.61, 95% CI 0.09-0.85)
    lvc <- log(36.5); label("Apparent central volume of distribution (V/F, L)")         # Table 3: tvV = 36.5 L (RSE 23.3%; bootstrap median 34.0, 95% CI 5.5-60.5)
    lcl <- log(3.6);  label("Apparent oral clearance at CRCL = 80 mL/min/1.73 m^2 and CONMED_RTV_AUC_12H = 12.2 mg*h/L (CL/F, L/h)")  # Table 3: tvCL = 3.6 L/h (RSE 7.1%; bootstrap median 3.2, 95% CI 2.6-3.7)

    # Covariate effects on CL/F. Chen 2024 Equation 1 (reproduced verbatim from
    # the publisher's equation image IDR-17-4055-e0001):
    #   CL = tvCL * (CrCL/80)^dCLdCrCL * (AUCRIT/12.2)^dCLdRIT * exp(nCL)
    # Both covariates enter as centred power functions. The paper writes the
    # ritonavir-AUC exponent as "dCLdRIT" in Equation 1 and as "dCLdAUCRIT" in
    # Table 3; these are the same parameter. The centring constants 80 and 12.2
    # appear only inside Equation 1 and are not repeated in any table.
    e_crcl_cl        <-  0.53; label("Power exponent of CKD-EPI creatinine clearance on nirmatrelvir CL/F (unitless)")             # Table 3: dCLdCrCL = 0.53 (RSE 21.7%; bootstrap median 0.58, 95% CI 0.23-1.15)
    e_rtv_auc_12h_cl <- -0.45; label("Power exponent of the ritonavir 12 h dosing-interval AUC on nirmatrelvir CL/F (unitless)")   # Table 3: dCLdAUCRIT = -0.45 (RSE magnitude 18.0%; bootstrap median -0.44, 95% CI -0.62 to -0.10)

    # Interindividual variability, exponential model on CL/F and V/F, with the
    # CL-V correlation carried in the off-diagonal element of the
    # variance-covariance matrix (Chen 2024 Results: "a correlation between CL
    # and V was observed and incorporated into the off-diagonal elements of the
    # variance-covariance matrix, resulting in a significantly improved OFV of
    # 95.0"). Table 3 reports variances directly (omega^2), so no CV%
    # back-transformation is needed. The off-diagonal covariance is
    # CorrV-CL * sqrt(omega2CL * omega2V) = 0.12 * sqrt(0.086 * 0.64) = 0.02815.
    # No interindividual variability was reported for ka.
    etalcl + etalvc ~ c(0.086,
                        0.02815, 0.64)  # Table 3: omega^2 CL = 0.086, omega^2 V = 0.64, CorrV-CL = 0.12 (bootstrap 0.086, 0.60, 0.12)

    # Residual variability. The paper selected "a one-compartment model with a
    # log-additive error option" (Results, Population PK Model Development); in
    # Phoenix NLME the log-additive error model is C = Cpred * exp(eps), which
    # is nlmixr2's lnorm() residual, so the reported stdev is the additive SD on
    # the log-transformed concentration scale.
    expSd <- 0.26; label("Additive residual SD on the log-transformed concentration scale (log-normal)")  # Table 3: stdev = 0.26 (RSE 17.4%; bootstrap median 0.26, 95% CI 0.16-0.37)
  })

  model({
    # Individual parameters (apparent, oral administration; F is not identifiable
    # separately and is subsumed into CL/F and V/F).
    # Chen 2024 Equation 1 for CL/F.
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) *
      (CRCL / 80)^e_crcl_cl *
      (CONMED_RTV_AUC_12H / 12.2)^e_rtv_auc_12h_cl
    vc <- exp(lvc + etalvc)

    # Micro-constant
    kel <- cl / vc

    # One-compartment model with first-order absorption
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    Cc <- central / vc

    Cc ~ lnorm(expSd)
  })
}
