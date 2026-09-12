Wang_2026_polymyxinB <- function() {
  description <- "Two-compartment intravenous population PK model for polymyxin B in critically ill adults, developed from a two-center Chinese ICU cohort sampled after at least the third dose (Wang 2026). CKD-EPI estimated glomerular filtration rate is the sole retained covariate, entering clearance as a power term normalized to the cohort median 42.88 mL/min/1.73 m^2 with exponent 0.43. Inter-individual variability on CL, V1 and Q; peripheral volume variability was fixed to zero. Combined proportional plus additive residual error."
  reference <- paste(
    "Wang Y, Wang X, Lei L, Sun W, Wu Z, Lan J, Chen J, Wang Y, Yao F,",
    "Hu L, Bai Y, Chen C.",
    "A multi-center study of population pharmacokinetics of polymyxin B",
    "in critically ill patients.",
    "Drug Des Devel Ther. 2026;20.",
    "doi:10.2147/DDDT.S521070. PMCID PMC13012298.",
    sep = " "
  )
  vignette <- "Wang_2026_polymyxinB"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Wang 2026: plasma polymyxin B was
  # quantified by HPLC-MS/MS (Methods, "Polymyxin B Administration and Sample
  # Collection") and the disposition model is the two-compartment structure of
  # the final-model equations on page 6 (V1 central, V2 peripheral).
  compartmentData <- list(
    central     = list(analyte = "polymyxinB", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "polymyxinB", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description        = "CKD-EPI estimated glomerular filtration rate (BSA-normalized)",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Wang 2026 Methods 'Study Design': 'CrCL was calculated using the Cockcroft-Gault formula, and eGFR was determined using the CKD-EPI formula.' The CKD-EPI eGFR is stored under the canonical CRCL column per inst/references/covariate-columns.md, whose CRCL entry explicitly accepts a CKD-EPI-estimated glomerular filtration rate in mL/min/1.73 m^2. Applied as a power covariate on CL, CL = 1.68 * (CRCL / 42.88)^0.43 * exp(eta_CL), which is the final-model equation printed on page 6 of the paper. IMPORTANT: the normalizing constant 42.88 mL/min/1.73 m^2 is the cohort MEDIAN eGFR and appears ONLY inside that typeset equation -- Table 1 reports the eGFR mean +/- SD (57.44 +/- 43.89) and never the median, and the equation is a vector graphic that plain text extraction drops entirely. See the vignette source-trace section. Both eGFR and Cockcroft-Gault CrCL were screened on CL in forward inclusion (eGFR dOFV -12.658; CrCL dOFV -5.705); eGFR gave the larger OFV drop and was the one retained, so the raw Cockcroft-Gault CrCL is deliberately NOT carried as a separate column here. Time-fixed per subject in this analysis.",
      source_name        = "eGFR"
    )
  )

  covariatesDataExcluded <- list(
    ECMO = list(
      description = "Extracorporeal membrane oxygenation support",
      units       = "(binary)",
      type        = "categorical",
      notes       = "Wang 2026 Methods: one of the three categorical covariates screened. Reached forward-inclusion significance on CL (dOFV -7.258, p < 0.05) and entered the full model, but was removed during backward elimination (p < 0.001 threshold). The Discussion attributes this to the very limited number of ECMO patients (n = 2, 3.57% of the cohort) and explicitly states the study 'could not robustly demonstrate an independent effect of ECMO on PB CL'."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Wang 2026 Methods: screened as a continuous covariate. Reached forward-inclusion significance on CL (dOFV -12.057, p < 0.05) and entered the full model alongside eGFR and ECMO, but was removed during backward elimination. Cohort value 40.21 +/- 78.90 U/L (Table 1)."
    ),
    BUN = list(
      description = "Blood urea nitrogen",
      units       = "mmol/L",
      type        = "continuous",
      notes       = "Wang 2026 Methods: screened as a continuous covariate. Reached forward-inclusion significance on CL (dOFV -9.301, p < 0.05) but was not retained in the final model. Cohort value 21.40 +/- 13.64 mmol/L (Table 1)."
    ),
    CRRT = list(
      description = "Continuous renal replacement therapy, delivered exclusively as continuous veno-venous hemofiltration (CVVH)",
      units       = "(binary)",
      type        = "categorical",
      notes       = "Wang 2026 Methods: recorded prospectively at each sampling time and time-aligned with each pharmacokinetic sample as a binary covariate (CRRT = 1 if ongoing at the exact sampling time, 0 otherwise), then tested with the categorical covariate structure of Eq.4. Not retained (CVVH dOFV -0.137, p > 0.05). 20 of 56 patients (35.71%) received CRRT, all via CVVH. A separate exploratory analysis found no significant correlation between effluent rate and drug clearance in the CVVH patients (absolute r = 0.158). The Discussion contrasts this null result with Hanafin et al, who reported significantly increased clearance under CVVHDF, and concludes the specific renal-replacement modality is a determinant that cannot be generalized."
    ),
    WT = list(
      description = "Total body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Wang 2026 Methods: screened as a continuous covariate. Not retained (dOFV -0.466, p > 0.05). Cohort value 60.73 +/- 10.78 kg (Table 1). Notable because the Discussion records that earlier polymyxin B analyses identified body weight as a covariate, and that Hanafin et al retained body weight on volume of distribution, whereas this final model did not."
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/L",
      type        = "continuous",
      notes       = "Wang 2026 Methods: screened as a continuous covariate. The Discussion reports it explicitly as the one named covariate tested on V1: 'During the covariate screening for factors influencing V1, no statistically significant covariate was identified, including ALB (dOFV -1.581, p > 0.05).' Cohort value 34.71 +/- 6.71 g/L (Table 1)."
    ),
    OTHER_SCREENED_COVARIATES = list(
      description = "Remaining candidate covariates screened and not retained",
      units       = "(various)",
      type        = "continuous",
      notes       = "Wang 2026 Methods lists the full screening set. Continuous: age, weight, BMI, ALT, AST, TP, ALB, TBIL, DBIL, BUN, eGFR, WBC, PLT, uric acid, Scr, and CrCL. Categorical: sex, CRRT status, and ECMO status. Only eGFR survived backward elimination. This entry groups age, BMI, AST, total protein, total bilirubin, direct bilirubin, white blood cell count, platelet count, uric acid, serum creatinine, Cockcroft-Gault creatinine clearance and sex, for which the paper reports no individual dOFV; the covariates that do carry a reported dOFV are given their own entries above."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 56L,
    n_studies      = 1L,
    age_mean       = "60.64 +/- 12.96 years",
    weight_mean    = "60.73 +/- 10.78 kg",
    sex_female_pct = 21.4,
    race_ethnicity = "Not reported (two-center Chinese ICU cohort)",
    disease_state  = "Critically ill adults (>= 18 years) in intensive care receiving intravenous polymyxin B sulfate. APACHE II 25.02 +/- 4.55. Pulmonary infection dominated (52 patients, 92.86%), with abdominal (4), sepsis (6), urinary tract (2), intracranial (1) and other (3) infections also present; patients could carry more than one infection type. Pathogens included Pseudomonas aeruginosa (39.29%), Klebsiella pneumoniae (26.79%) and Acinetobacter baumannii (25.00%). 20 patients (35.71%) received CRRT, exclusively CVVH; 2 patients (3.57%) received ECMO. Comorbidities: hypertension 33.93%, diabetes 16.07%, coronary disease 3.57%. Exclusions were pregnancy, allergy or intolerance to polymyxin B, and missing weight or renal function data.",
    dose_range     = "Intravenous polymyxin B sulfate, generally following the product label recommendation of 1.5-2.5 mg/kg per day (1 mg = 10,000 IU) divided into two doses. Regimens varied by treating clinician and were NOT adjusted for CRRT status or renal function; observed regimens included a 100 mg loading dose followed by 50 mg q12h, and maintenance doses of 50-100 mg q12h without a preceding load. Monte Carlo dosing simulations covered first and maintenance doses of 50-200 mg given q12h or q8h as 1-hour infusions.",
    regions        = "China (Guangdong Provincial People's Hospital and Maoming People's Hospital intensive care units), August 2020 to October 2022",
    renal_function = "CKD-EPI eGFR 57.44 +/- 43.89 mL/min/1.73 m^2 (cohort median 42.88, read from the final-model equation on page 6); Cockcroft-Gault creatinine clearance 150.79 +/- 145.79 mL/min; serum creatinine 190.18 +/- 154.35 umol/L; blood urea nitrogen 21.40 +/- 13.64 mmol/L. Renal function spanned the full range, and the dosing simulations were stratified into eGFR bands of <15, 15-30, 30-60, 60-90 and 90-130 mL/min/1.73 m^2.",
    n_observations = "350 polymyxin B plasma concentrations from 56 patients. Sampling began after at least the third dose, at seven time points per occasion: pre-dose (10 min before administration); 5 min, 1 h, 2 h, 4 h and 8 h after the end of infusion; and 10 min before the subsequent dose. Mean observed concentration 3.0774 +/- 2.1373 mg/L (Table 1).",
    notes          = "Prospective two-center study. NONMEM 7.3.0 with FOCE-I. Model evaluation used goodness-of-fit plots, prediction-corrected VPC (1000 replicates), NPDE and a 1000-replicate PsN bootstrap with a 93.5% success rate; all final estimates fell inside the bootstrap 95% CI (Table 2). Plasma quantified by HPLC-MS/MS using a previously published in-house method. The dosing recommendations of Table 3 target AUC24/MIC >= 50 with PTA >= 80% over MICs of 0.125-2 mg/L."
  )

  ini({
    # Structural parameters -- Wang 2026 Table 2 (final model) and the
    # final-model equations printed on page 6. The typical clearance is the
    # value at the cohort median eGFR of 42.88 mL/min/1.73 m^2.
    lcl <- log(1.68)  ; label("Clearance CL (L/h) at CRCL = 42.88 mL/min/1.73 m^2") # Table 2: CL = 1.68 L/h (RSE 11.4%; bootstrap median 1.66, 95% CI 1.23-2.03)
    lvc <- log(14.30) ; label("Central volume of distribution V1 (L)")               # Table 2: V1 = 14.30 L (RSE 7.9%; bootstrap median 14.15, 95% CI 12.42-16.17)
    lq  <- log(4.67)  ; label("Intercompartmental clearance Q (L/h)")                # Table 2: Q = 4.67 L/h (RSE 11.5%; bootstrap median 4.64, 95% CI 3.78-5.95)
    lvp <- log(48.84) ; label("Peripheral volume of distribution V2 (L)")            # Table 2: V2 = 48.84 L (RSE 23.5%; bootstrap median 50.64, 95% CI 27.85-81.98)

    # Covariate effect on CL -- Wang 2026 page 6 final model:
    #   CL_i = 1.68 * (eGFR/42.88)^0.43 * exp(eta_CL) L/h
    # This is the continuous power covariate form of Eq.3,
    # P_ij = P_tv,j * (COV/COV_median)^theta_j * exp(eta_j), with COV_median
    # the cohort median. Estimated (Table 2 gives RSE 16.6% and a bootstrap
    # 95% CI), so it is NOT wrapped in fixed().
    e_crcl_cl <- 0.43; label("Power exponent on (CRCL / 42.88 mL/min/1.73 m^2) for CL (unitless)") # Table 2: theta eGFR-CL = 0.43 (RSE 16.6%; bootstrap median 0.43, 95% CI 0.25-0.72)

    # Inter-individual variability. Wang 2026 Methods Eq.1 specifies an
    # exponential random-effects model, P_i = P_TV * exp(eta_i) with
    # eta ~ N(0, omega^2). Table 2 reports the IIV rows as 'omega_CL (%)',
    # 'omega_V1 (%)' and 'omega_Q (%)' -- the symbol named is omega itself,
    # not omega^2, so each percentage is read on the SD scale
    # (omega = pct/100) and squared to give the variance. Shrinkage is the
    # parenthesized figure in the Table 2 estimate column.
    # The same SD-scale reading of an 'omega (%)' row is used by the sibling
    # polymyxin B model Yang_2025_polymyxinB.R. See the vignette
    # 'Assumptions and deviations' section for the alternative exact-lognormal
    # reading and why it was not adopted.
    etalcl ~ 0.447561 # Table 2: omega CL = 66.9% (shrinkage 4%) -> omega = 0.669 -> variance 0.669^2
    etalvc ~ 0.236196 # Table 2: omega V1 = 48.6% (shrinkage 13%) -> omega = 0.486 -> variance 0.486^2
    etalq  ~ 0.352836 # Table 2: omega Q = 59.4% (shrinkage 29%) -> omega = 0.594 -> variance 0.594^2

    # Table 2 reports 'omega V2 (%) = 0 FIX', and the page 6 final-model
    # equation for V2 is written WITHOUT an exp(eta) term -- V2,i = 48.84 L --
    # unlike the CL, V1 and Q equations. Peripheral volume therefore carries no
    # inter-individual variability, encoded as an explicit zero-variance random
    # effect rather than omitted, so the fixed-to-zero status stays visible.
    etalvp ~ fixed(0) # Table 2: omega V2 = 0 FIX; page 6 equation V2,i = 48.84 L has no exp(eta)

    # Residual variability -- Wang 2026 Table 2, combined proportional plus
    # additive error (Eq.2: Y = F * (1 + eps1) + eps2). The Methods define
    # eps1 and eps2 as 'normally distributed with a mean of zero and variances
    # of sigma^2_prop and sigma^2_add', so the two Table 2 rows are NONMEM
    # $SIGMA VARIANCES and are converted to the SD scale here by taking square
    # roots. Reading them as SDs directly would imply a 1.36% proportional
    # error and a 0.0858 mg/L additive error against a mean observed
    # concentration of 3.0774 mg/L -- a total residual near 3%, far tighter
    # than an HPLC-MS/MS ICU popPK dataset supports and inconsistent with the
    # width of the Figure 2 pcVPC bands. See the vignette source-trace section.
    propSd <- 0.116619; label("Proportional residual error (fraction)") # Table 2: proportional error variance = 0.0136 (shrinkage 19%, RSE 38%) -> sqrt = 0.116619
    addSd  <- 0.292916; label("Additive residual error (mg/L)")         # Table 2: additive error variance = 0.0858 (shrinkage 19%, RSE 45%) -> sqrt = 0.292916
  })

  model({
    # Individual parameters -- Wang 2026 page 6 final-model equations. The
    # normalizing constant 42.88 mL/min/1.73 m^2 is the cohort median eGFR and
    # is printed only inside the typeset CL equation.
    cl <- exp(lcl + etalcl) * (CRCL / 42.88)^e_crcl_cl
    vc <- exp(lvc + etalvc)
    q  <- exp(lq + etalq)
    vp <- exp(lvp + etalvp)

    # Micro-constants
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Two-compartment disposition with first-order elimination from the
    # central compartment. Polymyxin B is given as an intravenous infusion, so
    # drug enters the central compartment directly (Methods: 'Enrolled patients
    # received intravenous infusions of PB sulfate'; the simulations of Table 3
    # use a 1-hour infusion).
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                  k12 * central - k21 * peripheral1

    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
