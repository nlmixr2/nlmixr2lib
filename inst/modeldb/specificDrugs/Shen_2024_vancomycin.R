Shen_2024_vancomycin <- function() {
  description <- "One-compartment IV-infusion population PK model for vancomycin in 386 Southern Chinese children (Shen 2024). Clearance uses an age-cutoff structure at 2 years: a separate typical clearance AND a separate body-weight allometric exponent are estimated in each age stratum (2.59 L/h with exponent 0.38 for age > 2 years; 1.98 L/h with exponent 0.739 for age <= 2 years), both normalized to a 12 kg reference weight, with a shared Cockcroft-Gault creatinine-clearance power effect (exponent 0.517, reference 75 mL/min). Central volume scales linearly with body weight (22.4 L at 12 kg). Additive residual error is estimated separately for each of the two study centers. Both strata share the volume, the CLcr exponent and the clearance IIV in a single joint NONMEM fit."
  reference <- "Shen X, Li X, Lu J, Zhu J, He Y, Zhang Z, Chen Z, Zhang J, Fan X, Li W. Population pharmacokinetic analysis for dose regimen optimization of vancomycin in Southern Chinese children. CPT Pharmacometrics Syst Pharmacol. 2024;13(7):1201-1213. doi:10.1002/psp4.13151"
  vignette <- "Shen_2024_vancomycin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # The two additive residual SDs are center-specific (Shen 2024 Table 3
  # rows "RV1" / "RV2"); the canonical single `addSd` cannot express them.
  paper_specific_residual_sds <- c("addSdCenter1", "addSdCenter2")

  compartmentData <- list(
    central = list(analyte = "vancomycin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Shen 2024 Table 1 (Total model-building column): median 11.95 kg, range 0.88-49.5, mean 13.55 (SD 8.81). The reference weight 12 kg used in Eq. (7)-(9) is the cohort median rounded by the authors. Enters CL as (WT/12)^k with a stratum-specific exponent k (0.38 for age > 2 y, 0.739 for age <= 2 y; Table 3) and enters V as (WT/12)^1, i.e. simple linear weight proportionality (Eq. 9 and Eq. 1, the Holford size form).",
      source_name        = "BW"
    ),
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Shen 2024 Table 1 (Total): median 2.22 years, range 0.08-17.45, mean 3.42 (SD 3.48); 182 of 386 patients (47%) were younger than 2 years. Used ONLY as the 2-year cutoff that selects between the two clearance strata of Eq. (7) and Eq. (8) -- the paper's 'Model IV' age-cutoff structure. It does not enter the model as a continuous term: a continuous AGE effect (model 3) and a sigmoidal maturation function (Model III) were both tested and rejected in favour of the cutoff (Table 2, Supplementary Table 1).",
      source_name        = "Age"
    ),
    CRCL = list(
      description        = "Creatinine clearance by the Cockcroft-Gault formula (raw, NOT BSA-normalized)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Shen 2024 Methods 'Population pharmacokinetic model development': 'The CLcr was computed using the Cockroft-Gault formula.' Reported in raw mL/min, not normalized to 1.73 m^2. Table 1 (Total): median 74.16 mL/min, range 10.5-265.6, mean 79.95 (SD 43.33). The reference value 75 mL/min in Eq. (7)-(8) is the cohort median rounded by the authors. Enters CL as (CRCL/75)^0.517, shared across both age strata. Stored under canonical CRCL, which explicitly admits the raw Cockcroft-Gault variant (precedents Delattre_2010_amikacin.R, Chen_2023_nemonoxacin.R, Wada_2023_sparsentan.R); the assay form is documented here per the register's per-model requirement.",
      source_name        = "CLcr"
    ),
    STUDY_VANCO_CENTER2 = list(
      description        = "Shen 2024 vancomycin second-study-center indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (center N1, n = 210 patients / 272 concentrations)",
      notes              = "1 = subject enrolled at the second of the two study centers (Table 1 column N2, n = 176 patients / 249 concentrations); 0 = first center (column N1). Selects between the two additive residual-error magnitudes: RV1 = 4.64 mg/L and RV2 = 4.53 mg/L (Table 3). Shen 2024 Results: 'RV was best characterized by an additive error model for each of the two study centers.' The paper does NOT state which of Baoan Women's and Children's Hospital / Shenzhen Children's Hospital is N1 and which is N2, so the indicator is defined by the Table 1 column label rather than by hospital name. Note these are the two hospitals, not the two age strata: N1 has median age 0.95 y and N2 median 4.32 y, and both centers contribute patients on both sides of the 2-year cutoff.",
      source_name        = "Study center"
    )
  )

  # Covariates the paper screened but did not retain in the final model. They
  # are documentation only: no point estimate was published for any of them,
  # so none can be encoded.
  covariatesDataExcluded <- list(
    SEXF = list(
      description        = "Female sex indicator. Screened as a demographic covariate on CL (Methods 'Population pharmacokinetic model development') and rejected.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Shen 2024 Discussion: 'Other variables that produced no significant impact during the modeling process were Alb, sex, and concomitant therapies.' Table 1 (Total) reports 239 male / 147 female. No point estimate published."
    ),
    PNA = list(
      description        = "Postnatal age. Screened as a demographic covariate on CL and rejected.",
      units              = "weeks",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Supplementary Table 1 model 4: 'Add PNA on CL in model 1' gave OFV 2388.505 (dOFV -180.950), worse than the retained body-weight term (model 2, OFV 2320.382, dOFV -249.073), so PNA was not carried forward. Table 1 (Total): median 116.05 weeks, range 4.43-909.9. No point estimate published for the final model."
    ),
    ALB = list(
      description        = "Serum albumin. Screened as a hepatic-function covariate on CL and rejected.",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Shen 2024 Discussion lists Alb among the variables with no significant impact. Table 1 (Total): median 35.15 g/L, range 0.61-59.3. No point estimate published."
    ),
    BUN = list(
      description        = "Blood urea nitrogen. Screened as a renal-function covariate on CL; significant on its own but not selected.",
      units              = "mmol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Supplementary Table 1: adding Bun to model 2 (model 8) gave OFV 2269.018, marginally lower than adding CLcr (model 7, OFV 2272.810), but Results states 'the graphical diagnostics of the covariate-parameter relationship by the CLcr-adjusted model performed superior to that of the Bun', so CLcr was retained and Bun dropped. Table 1 (Total): median 3.18 mmol/L, range 0.7-25.5. No final-model point estimate published."
    ),
    CREAT = list(
      description        = "Serum creatinine. Screened as a renal-function covariate on CL and rejected in favour of the derived CLcr.",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Shen 2024 Discussion: 'Since Scr and Bun are known to be influenced by age, sex, muscle mass, and diet-limiting their utility as a marker of the GFR.' Table 1 (Total): median 26 umol/L, range 7.6-111. Serum creatinine is an input to the Cockcroft-Gault CRCL that IS retained, but does not enter the model separately. No point estimate published."
    ),
    CONMED_MEROPENEM = list(
      description        = "Concomitant meropenem indicator. Screened as a co-medication covariate on CL and rejected.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant meropenem)",
      notes              = "Screened because it was co-administered in more than 10% of patients (Methods). Table 1 (Total): 228 of 521 records (43.76%). Shen 2024 Discussion: concomitant therapies produced no significant impact. No point estimate published."
    ),
    CONMED_CEFOPERAZONE = list(
      description        = "Concomitant cefoperazone indicator. Screened as a co-medication covariate on CL and rejected.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant cefoperazone)",
      notes              = "Screened because it was co-administered in more than 10% of patients (Methods). Table 1 (Total): 135 (25.91%). No point estimate published."
    ),
    CONMED_CEFTRIAXONE = list(
      description        = "Concomitant ceftriaxone indicator. Screened as a co-medication covariate on CL and rejected.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant ceftriaxone)",
      notes              = "Screened because it was co-administered in more than 10% of patients (Methods). Table 1 (Total): 134 (25.72%). No point estimate published."
    ),
    CONMED_MANNITOL = list(
      description        = "Concomitant mannitol indicator. Screened as a co-medication covariate on CL and rejected.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant mannitol)",
      notes              = "Screened because it was co-administered in more than 10% of patients (Methods). Table 1 (Total): 66 (12.67%). No point estimate published."
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 386L,
    n_studies        = 2L,
    age_range        = "0.08-17.45 years (1 month to 17.45 years)",
    age_median       = "2.22 years (mean 3.42, SD 3.48); 182 of 386 patients (47%) were younger than 2 years",
    weight_range     = "0.88-49.5 kg",
    weight_median    = "11.95 kg (mean 13.55, SD 8.81)",
    sex_female_pct   = 38.1,
    race_ethnicity   = "Not reported. Southern Chinese pediatric cohort (Shenzhen, Guangdong province).",
    disease_state    = "Hospitalized children with Gram-positive infections who received intravenous vancomycin for more than 3 days. Exclusions: age under 1 month, hypersensitivity to vancomycin or its excipients, dialysis or blood purification during vancomycin treatment, Gram-positive bacterial colonization only. Principal pathogens were MRSA, methicillin-resistant coagulase-negative staphylococci, and Enterococcus species.",
    dose_range       = "Intravenous infusion over at least 60 min. Median vancomycin dose 43.24 mg/kg/day (range 12.3-79.4); the paper quotes 'a median daily dose of vancomycin was 43 mg/kg/day'.",
    regions          = "China (Shenzhen, Guangdong): Baoan Women's and Children's Hospital and Shenzhen Children's Hospital, 2016-2022.",
    renal_function   = "Cockcroft-Gault creatinine clearance median 74.16 mL/min (range 10.5-265.6, mean 79.95, SD 43.33); serum creatinine median 26 umol/L (range 7.6-111).",
    n_concentrations = 521L,
    notes            = "Baseline demographics from Shen 2024 Table 1 (Total model-building column). 386 patients contributing 521 vancomycin concentrations were used to build the model; a further 67 patients (76 concentrations) formed a randomly selected external validation set (15% of the 453 enrolled patients). Center N1 contributed 210 patients / 272 concentrations (median age 0.95 y), center N2 176 patients / 249 concentrations (median age 4.32 y). Routine therapeutic drug monitoring samples: troughs collected within 30 min before the 4th-5th dose, or peaks 0.5-1 h after the end of infusion; assayed by homogeneous enzyme immunoassay (Viva-E, Siemens) with within- and between-assay CV < 10%. Fit in NONMEM 7.5 with FOCE-I. Final model = Supplementary Table 1 model 7 (OFV 2272.810). Internal evaluation by GOF plots, VPC and a 1000-replicate bootstrap (996 successful; relative error < +/-2.6%); external evaluation gave MPE -14.08% and MAPE 24%. The paper also derives an efficacy target of AUC(0-24)/MIC >= 260 by logistic regression and ROC analysis on 180 culture-positive patients (29 treatment failures); that analysis is a statistical model, not a PK/PD structural model, and is not part of this extraction."
  )

  ini({
    # ---------------------------------------------------------------------
    # Structural parameters. Shen 2024 Table 3, final-model "Estimate"
    # column; the equations are Eq. (7)-(9) on p. 1205.
    #
    #   if age >  2 years  CL (L/h) = 2.59 * (BW/12)^0.38  * (CLcr/75)^0.517 * e^eta
    #   if age <= 2 years  CL (L/h) = 1.98 * (BW/12)^0.739 * (CLcr/75)^0.517 * e^eta
    #                       V (L)   = 22.4 * (BW/12)
    #
    # This is ONE jointly-fit model (single OFV 2272.810, Supplementary
    # Table 1 model 7): the volume, the CLcr exponent, the CL IIV and both
    # residual terms are shared across the two age strata. Only the typical
    # CL and the body-weight exponent are stratum-specific, so each of those
    # two quantities carries an explicit `_agele2` / `_agegt2` suffix
    # (symmetric stratum-suffix scheme; neither stratum is the reference).
    # ---------------------------------------------------------------------

    lcl_agegt2 <- log(2.59); label("Typical clearance for age > 2 years at WT = 12 kg, CRCL = 75 mL/min (L/h)")   # Shen 2024 Table 3: theta(CL age>2) = 2.59 L/h (RSE 6.7%, bootstrap 95% CI 2.30-2.88); Eq. (7)
    lcl_agele2 <- log(1.98); label("Typical clearance for age <= 2 years at WT = 12 kg, CRCL = 75 mL/min (L/h)")  # Shen 2024 Table 3: theta(CL age<=2) = 1.98 L/h (RSE 6.8%, bootstrap 95% CI 1.75-2.21); Eq. (8)
    lvc        <- log(22.4); label("Central volume of distribution at WT = 12 kg (L)")                            # Shen 2024 Table 3: theta(Vd) = 22.4 L (RSE 23.5%, bootstrap 95% CI 15.05-32.02); Eq. (9)

    # ---------------------------------------------------------------------
    # Covariate effects. The body-weight exponent on CL is estimated
    # separately in each age stratum; the CLcr exponent is shared.
    # ---------------------------------------------------------------------

    e_wt_cl_agegt2 <- 0.38;  label("Power exponent on (WT/12) for CL, age > 2 years (unitless)")   # Shen 2024 Table 3: theta(BW age>2)  = 0.38  (RSE 29.8%, bootstrap 95% CI 0.193-0.571); Eq. (7)
    e_wt_cl_agele2 <- 0.739; label("Power exponent on (WT/12) for CL, age <= 2 years (unitless)")  # Shen 2024 Table 3: theta(BW age<=2) = 0.739 (RSE 16.7%, bootstrap 95% CI 0.542-0.948); Eq. (8)
    e_crcl_cl      <- 0.517; label("Power exponent on (CRCL/75) for CL, both age strata (unitless)")  # Shen 2024 Table 3: theta(CLcr) = 0.517 (RSE 15.1%, bootstrap 95% CI 0.385-0.635); Eq. (7) and Eq. (8)

    # The weight exponent on V is structurally 1 (Eq. 9 writes V = 22.4 *
    # (BW/12), i.e. simple linear proportionality following the Holford size
    # model of Eq. 1). Table 3 reports no exponent for V, so it was not
    # estimated.
    e_wt_vc <- fixed(1.0); label("Power exponent on (WT/12) for V (unitless)")  # Shen 2024 Eq. (9): V = 22.4 * (BW/12), exponent 1 by construction; no theta in Table 3

    # ---------------------------------------------------------------------
    # Inter-individual variability. Exponential model on CL only (Methods:
    # 'an exponential error model was adopted to evaluate IIV of PK
    # parameters'), matching the explicit `* e^eta` in Eq. (7)-(8).
    #
    # Table 3 reports omega(CL) = 0.319 as an SD on the log scale, NOT a
    # variance: the row symbol is `omega (CL)` (not omega^2) and the same
    # table uses the parallel Greek symbol `sigma` for the two residual rows,
    # which are unambiguously SDs because they are quoted in mg/L. So
    # var(etalcl) = 0.319^2 = 0.101761, i.e. approximately 32.7% CV. See the
    # vignette Errata for the corroborating check against Figure 3.
    #
    # IIV on V was estimated at < 1e-3 and removed by the authors (Methods),
    # so there is deliberately no etalvc.
    # ---------------------------------------------------------------------

    etalcl ~ 0.101761  # 0.319^2  # Shen 2024 Table 3: omega(CL) = 0.319 (RSE 11.7%, bootstrap 95% CI 0.240-0.373); eta-shrinkage 43.1%

    # ---------------------------------------------------------------------
    # Residual error: additive, estimated separately for each study center
    # (Results: 'RV was best characterized by an additive error model for
    # each of the two study centers').
    # ---------------------------------------------------------------------

    addSdCenter1 <- 4.64; label("Additive residual error at study center 1 (mg/L)")  # Shen 2024 Table 3: sigma(additive1) = 4.64 mg/L (RSE 9.4%, bootstrap 95% CI 3.80-5.22); epsilon-shrinkage 13.3%
    addSdCenter2 <- 4.53; label("Additive residual error at study center 2 (mg/L)")  # Shen 2024 Table 3: sigma(additive2) = 4.53 mg/L (RSE 9.0%, bootstrap 95% CI 3.81-5.22); epsilon-shrinkage 12.4%
  })

  model({
    # ---- 1. Age-stratum indicator (Shen 2024 Eq. 7-8, cutoff at 2 years) ----
    is_agele2 <- (AGE <= 2.0)

    # ---- 2. Individual PK parameters ----
    # Typical CL and the body-weight exponent are stratum-specific; the CLcr
    # power effect and the exponential IIV are shared across strata.
    cl_wt <- exp(lcl_agele2) * (WT / 12.0)^e_wt_cl_agele2 * is_agele2 +
      exp(lcl_agegt2) * (WT / 12.0)^e_wt_cl_agegt2 * (1.0 - is_agele2)
    cl <- cl_wt * (CRCL / 75.0)^e_crcl_cl * exp(etalcl)

    vc <- exp(lvc) * (WT / 12.0)^e_wt_vc

    kel <- cl / vc

    # ---- 3. ODE system (one compartment, first-order elimination; NONMEM
    #         ADVAN1 TRANS2). Vancomycin is given as an IV infusion of at
    #         least 60 min, so dosing goes straight into `central`. ----
    d/dt(central) <- -kel * central

    # ---- 4. Observation and center-specific additive residual error ----
    # Dose in mg, vc in L -> central/vc has units mg/L.
    Cc <- central / vc

    addSd_i <- addSdCenter1 * (1.0 - STUDY_VANCO_CENTER2) +
      addSdCenter2 * STUDY_VANCO_CENTER2
    Cc ~ add(addSd_i)
  })
}
