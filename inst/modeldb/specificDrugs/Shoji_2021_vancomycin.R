Shoji_2021_vancomycin <- function() {
  description <- "One-compartment IV population PK model for vancomycin in pediatric liver transplant recipients aged 1-186 months (Shoji 2021), developed from 1,158 routine therapeutic-drug-monitoring concentrations in 161 children (270 treatment episodes) at a single Japanese pediatric transplant center. Clearance scales with body weight to the fixed 0.75 power and with power functions of serum creatinine (reference 0.16 mg/dL, exponent -0.70) and days from liver transplantation to the start of vancomycin (reference 17 days, exponent -0.09), so clearance is highest early after transplant and at low serum creatinine; volume of distribution is proportional to body weight. The final-model between-subject and residual variances are not reported (Table 3 prints them for the base model only) and are encoded as zero."
  reference <- "Shoji K, Saito J, Nakagawa H, Funaki T, Fukuda A, Sakamoto S, Kasahara M, Momper JD, Capparelli EV, Miyairi I. Population pharmacokinetics and dosing optimization of vancomycin in pediatric liver transplant recipients. Microbiol Spectr. 2021;9(2):e00460-21. doi:10.1128/Spectrum.00460-21"
  vignette <- "Shoji_2021_vancomycin"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Vancomycin was given as an intermittent intravenous
  # infusion (Shoji 2021 Discussion: "conventional intermittent intravenous
  # infusion regimens"), so the dose enters `central` directly. Shoji 2021
  # Methods measure "serum vancomycin concentrations" by immunoassay.
  compartmentData <- list(
    central = list(analyte = "vancomycin", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters clearance as the UNNORMALISED allometric term WT^0.75 and volume as WT^1 (Shoji 2021 Results equations and Table 3 footnotes b and c); both exponents were fixed before covariate evaluation (Methods). Table 1 median 9.1 kg (IQR 6.8-16.2, range 3.1-61.0) over 270 treatment episodes.",
      source_name = "wt"
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters clearance as the power term (CREAT/0.16)^-0.70; 0.16 mg/dL is the Table 1 cohort median (IQR 0.12-0.23, range 0.06-5.43). Univariate screening dOFV -435.09 (Table 2), the dominant clearance covariate. Patients on renal replacement therapy were excluded.",
      source_name = "sCr"
    ),
    POD = list(
      description = "Days from liver transplantation to the first day of the vancomycin treatment episode",
      units = "days",
      type = "continuous",
      reference_category = NULL,
      notes = "Shoji 2021 column DFLT, defined in Results as 'the number of days from the date of LT to the first day of vancomycin treatment'. It is therefore a per-episode CONSTANT fixed at the start of vancomycin, not a per-observation day count: the analysis used data within 14 days of vancomycin initiation, and the same DFLT applies to every concentration in the episode. Supply the value at the start of the episode on every row. Enters clearance as (POD/17)^-0.09; 17 days is the Table 1 cohort median (IQR 6-31, range 0-357). CAUTION: with a negative exponent the printed term is infinite at POD = 0, although Table 1 reports episodes starting on the day of transplant. Shoji 2021 does not say how those episodes were coded, so the model is encoded exactly as printed and needs POD > 0; the vignette uses POD >= 1.",
      source_name = "DFLT"
    )
  )

  # Screened in Shoji 2021's univariate and stepwise covariate analysis
  # (Methods 'Measurement of vancomycin concentrations and pharmacokinetics
  # analysis'; Table 2) and NOT retained in the final model.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Subject age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Table 1 median 13.3 months (IQR 7.6-53.5, range 1-186 months). Univariate dOFV -0.34 on CL and +2.64 on V (Table 2); not retained.",
      source_name = "Age"
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "male",
      notes = "Table 1: 72 male (44.7%) / 89 female (55.3%). Univariate dOFV +11.65 on CL and +11.95 on V (Table 2); not retained.",
      source_name = "Sex"
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Table 1 median 3.1 g/dL (IQR 2.8-3.3). Univariate dOFV +2.49 on CL and +1.36 on V (Table 2); not retained.",
      source_name = "Albumin"
    ),
    ALT = list(
      description = "Serum alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Table 1 median 54.7 U/L (IQR 27.4-134.9). Univariate dOFV +3.82 on CL and +1.90 on V (Table 2); not retained.",
      source_name = "ALT"
    ),
    DIS_BILIARY_ATRESIA = list(
      description = "Biliary atresia as the underlying disease requiring liver transplantation",
      units = "(binary)",
      type = "binary",
      reference_category = "other underlying disease",
      notes = "Table 2 footnote b: 'biliary atresia was set as a value of 1, and other disorders were set as a value of 0'. 84 of 161 patients (52.2%, Table 1). Univariate dOFV -18.86 on CL; at stepwise step 1 (after sCr on CL) adding it worsened the fit (dOFV +41.92, Table 2); not retained.",
      source_name = "Underlying diseases"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 161L,
    n_studies = 1L,
    n_sites = 1L,
    n_episodes = 270L,
    n_concentrations = 1158L,
    age_range = "1-186 months (Table 1; median 13.3 months, IQR 7.6-53.5); under 18 years by inclusion, with no neonates and few adolescents",
    weight_range = "3.1-61.0 kg (Table 1; median 9.1 kg, IQR 6.8-16.2)",
    sex_female_pct = 55.3,
    race_ethnicity = "Not reported (single Japanese center)",
    disease_state = "Pediatric liver transplant recipients with suspected or proven infection receiving intravenous vancomycin. Underlying diseases: biliary atresia 52.2%, metabolic disease 18.6%, fulminant hepatitis 16.1%, liver cirrhosis 4.3%, liver fibrosis 4.3%, vascular abnormalities 2.5%, liver tumor 1.9% (Table 1). Patients on renal replacement therapy were excluded.",
    dose_range = "Median 15.0 mg/kg/dose (IQR 14.0-15.0, range 4.3-37.0). Institutional starting regimen 15 mg/kg every 6 h for ages 1 month to 12 years and 15 mg/kg every 8 h for ages 13-17 years, adjusted to trough targets of 10-15 ug/mL (15-20 ug/mL in critically ill patients).",
    regions = "Japan (National Center for Child Health and Development, Tokyo)",
    renal_function = "Serum creatinine median 0.16 mg/dL (IQR 0.12-0.23, range 0.06-5.43).",
    co_medication = "Tacrolimus (trough level median 8.2 ug/mL, IQR 5.5-10.4).",
    notes = "Retrospective analysis of electronic-medical-record data, 2006-2014; data within 14 days of vancomycin initiation. Days from liver transplantation to vancomycin start median 17 (IQR 6-31, range 0-357). Serum concentrations measured by fluorescence polarization immunoassay (AxSYM, 2006-2010) and chemiluminescent immunoassay (Architect, 2010-2014). Fit in Phoenix NLME 8.2 by FOCE; evaluated by 1,000-sample bootstrap (100% convergence) and VPC."
  )

  ini({
    # Structural parameters, Shoji 2021 Table 3 'Final model estimate (RSE%)'
    # column, with the final-model equations in Results and Table 3 footnotes
    # b and c:
    #   CL (L/h) = theta_CL * wt^0.75 * (sCr/0.16)^theta_sCr * (DFLT/17)^theta_DFLT * exp(eta_CL)
    #   V  (L)   = theta_V * wt * exp(eta_V)
    # theta_CL is clearance per kg^0.75 at sCr = 0.16 mg/dL and DFLT = 17
    # days (Table 3 unit 'l/kg^0.75/h'); theta_V is volume per kg.
    lcl <- log(0.29); label("Clearance per kg^0.75 at CREAT = 0.16 mg/dL and POD = 17 days (L/h/kg^0.75)") # Shoji 2021 Table 3 theta_CL final 0.29 (RSE 4.13%); bootstrap median 0.29 (95% CI 0.26-0.32); base model 0.27
    lvc <- log(1.00); label("Volume of distribution per kg body weight (L/kg)") # Shoji 2021 Table 3 theta_V final 1.00 (RSE 5.97%); bootstrap median 1.00 (95% CI 0.89-1.14); base model 1.23

    # Allometric exponents fixed before covariate evaluation (Shoji 2021
    # Methods: 'The TVCL was scaled allometrically by subject weight
    # (weight^0.75), and typical value of the V (TVV) was also scaled by
    # subject weight (weight^1.0) before evaluation of other covariates').
    e_wt_cl <- fixed(0.75); label("Allometric exponent on WT for CL (unitless)") # Shoji 2021 Methods and Table 3 footnote b
    e_wt_vc <- fixed(1); label("Allometric exponent on WT for V (unitless)") # Shoji 2021 Methods and Table 3 footnote c

    # Covariate effects on clearance (estimated; Table 3 reports RSE and
    # bootstrap CI for both).
    e_creat_cl <- -0.70; label("Power exponent on (CREAT/0.16 mg/dL) for CL (unitless)") # Shoji 2021 Table 3 theta_sCr final -0.70 (RSE 8.23%); bootstrap median -0.69 (95% CI -0.79 to -0.57)
    e_pod_cl <- -0.09; label("Power exponent on (POD/17 days) for CL (unitless)") # Shoji 2021 Table 3 theta_DFL final -0.09 (RSE 34.8%); bootstrap median -0.09 (95% CI -0.15 to -0.03)

    # Between-subject variability. Shoji 2021 declares exponential IIV on CL
    # and V (Results equations exp(eta_CL), exp(eta_V)) and reports final-
    # model eta shrinkage of 22.0% (CL) and 27.0% (V), but Table 3 prints the
    # omega values ('% omega CL' 46.47, '% omega Vc' 48.14) ONLY in the
    # BASE-model column; the final-model and bootstrap columns are blank. The
    # base-model CL value predates the serum-creatinine effect (univariate
    # dOFV -435), so it overstates the final-model CL variability and would
    # give a chimeric model the authors never fitted. The final-model
    # variances are therefore encoded as zero rather than borrowed; see the
    # vignette Errata for the base-model values.
    etalcl ~ fixed(0) # Shoji 2021 Results (exp(eta_CL) declared); final-model magnitude not reported (Table 3 base model only: 46.47%)
    etalvc ~ fixed(0) # Shoji 2021 Results (exp(eta_V) declared); final-model magnitude not reported (Table 3 base model only: 48.14%)

    # Residual variability. Results: 'an additive/proportional error model
    # was selected'. Table 3 prints only a proportional term, 56.5%, and only
    # for the base model; no additive magnitude is printed anywhere and
    # epsilon shrinkage (4.55%) is the only final-model residual statistic.
    # Both magnitudes are encoded as zero.
    propSd <- fixed(0); label("Proportional residual SD (fraction; 0 -- final-model value not reported)") # Shoji 2021 Results combined error selected; Table 3 base model only: 56.5%
    addSd <- fixed(0); label("Additive residual SD (ug/mL; 0 -- not reported)") # Shoji 2021 Results combined error selected; no additive value printed
  })
  model({
    # Shoji 2021 Results final-model equations and Table 3 footnotes b-c.
    # WT enters unnormalised (theta_CL is per kg^0.75, theta_V per kg).
    # POD must be > 0: the negative exponent makes (POD/17)^e_pod_cl
    # infinite at POD = 0 (see covariateData notes).
    cl <- exp(lcl + etalcl) * WT^e_wt_cl * (CREAT / 0.16)^e_creat_cl * (POD / 17)^e_pod_cl
    vc <- exp(lvc + etalvc) * WT^e_wt_vc

    kel <- cl / vc

    d/dt(central) <- -kel * central

    # Dose in mg, volume in L, so central/vc is mg/L == ug/mL.
    Cc <- central / vc
    Cc ~ prop(propSd) + add(addSd)
  })
}
