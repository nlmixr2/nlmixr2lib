Wang_2024_valproic_acid <- function() {
  description <- "One-compartment population PK model with first-order absorption for total plasma valproic acid in Chinese children with epilepsy (Wang 2024 final model), built to tailor the dose when switching between oral syrup and sustained-release tablets. Apparent clearance carries a power body-weight effect and a proportional female-sex effect; apparent volume carries a power body-weight effect. Formulation-specific absorption rate constants are FIXED from the literature (oral syrup 2.64 1/h reference, sustained-release tablet 0.46 1/h), because the therapeutic-drug-monitoring dataset is steady-state troughs only and contains no absorption-phase data."
  reference <- "Wang WJ, Li Y, Hu YH, Wang J, Zhang YY, Fan L, Dai HR, Guo HL, Ding XS, Chen F. Population pharmacokinetics of valproic acid in children with epilepsy: Implications for dose tailoring when switching from oral syrup to sustained-release tablets. CPT Pharmacometrics Syst Pharmacol. 2024;13(9):1555-1568. doi:10.1002/psp4.13191. PMCID PMC11533106. Final-model parameter estimates from Table 3; covariate equations from Equations 34 and 35."
  vignette <- "Wang_2024_valproic_acid"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    depot   = list(analyte = "valproic acid", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "valproic acid", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = "22 kg (the model-development cohort median; Wang 2024 Table 2)",
      notes              = "Power effect on both CL/F (exponent 0.717) and V/F (exponent 0.524), each normalised to the 22.0 kg cohort median. Cohort range 6.0-95.0 kg (Wang 2024 Table 2). The authors compared five body-weight / age scaling forms (fixed allometric, simple exponent, sigmoid maturation, body-weight-dependent exponent, age-dependent exponent; Methods Equations 6-13) and selected the SIMPLE EXPONENT model (Model II) with both exponents estimated: the two dependent-exponent models fitted better on OFV/AIC but their kmax and Hill terms had RSE > 60%, and TM50/kmax/Hill could not be estimated at all in the sigmoid maturation model. Note the estimated CL/F exponent 0.717 is close to, but was NOT fixed at, the allometric 0.75. Age is not in the final model because it correlates with body weight at r = 0.948 in this cohort (Wang 2024 Discussion).",
      source_name        = "BW"
    ),
    SEXF = list(
      description        = "Female sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "PROPORTIONAL (not exponential) effect on CL/F: the paper's Equation 34 is CL/F = 0.196 * (BW/22)^0.717 * (1 - 0.0436 * Sex) with 'Sex = 0 for male, Sex = 1 for female' stated immediately below the equation, so the paper's Sex indicator is already the canonical SEXF orientation and needs no recoding. Females therefore have 4.36% LOWER apparent clearance and hence 1/(1 - 0.0436) = 4.6% HIGHER exposure at an equivalent dose, which matches the Discussion's 'Compared with males, females exhibited 5-10% higher exposure of VPA at equivalent dose (Figure 2)'. Note the DIRECTION is opposite to the sibling Zhang_2024_valproic_acid.R, where women had 12.9% higher clearance; Zhang attributes its direction to the larger median body weight of the women in that adult-inclusive cohort. 194 of 498 model-development patients (39.0%) were female (Wang 2024 Table 2). Sex was retained on the strength of an OFV drop of 11.727 (p < 0.001) over the covariate-free Model II.",
      source_name        = "Sex"
    ),
    FORM_VPA_SR = list(
      description        = "Sustained-release valproic acid tablet formulation indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (oral syrup, the reference formulation in this cohort)",
      notes              = "Selects the FIXED sustained-release-tablet absorption rate constant Ka = 0.46 1/h; the oral-syrup reference is Ka = 2.64 1/h. Wang 2024 Methods, Base model: 'The fixed values for ka were 2.64 and 0.46 h-1 for syrup and SR tablets, respectively, due to a lack of absorption phase data', citing Mei 2018 (reference 19) - the same literature Ka pair that Zhang_2023_* and Zhang_2024_valproic_acid.R trace to Ding 2015. This cohort has only the two levels syrup and sustained-release tablet, so FORM_TABLET (the conventional immediate-release level used by Zhang_2023_valproic_acid_base.R) is not part of this model. Across the full 617-patient cohort 471 patients received syrup and 154 received SR tablets; 8 patients started on syrup and later switched, which is why the Table 2 per-arm percentages sum to slightly more than 100%.",
      source_name        = "Dosage forms"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Median 5.8 years, range 0.2-17.2 (Wang 2024 Table 2). Screened as a maturation covariate on CL/F in Models III (sigmoid maturation) and V (age-dependent exponent) but NOT retained: TM50, kmax and Hill 'were failed to be estimated in our study' (Discussion), and the strong body-weight/age correlation (r = 0.948) 'necessitates the exclusion of one of these variables, leading to the simplification of the model into a basic exponent model'."
    ),
    DOSE_VPA_MGKGD = list(
      description = "Total daily valproic acid dose per kg body weight",
      units       = "mg/kg/d",
      type        = "continuous",
      notes       = "Median 21.7 mg/kg/day, range 7.7-40.0 (Wang 2024 Table 2). Screened as a surrogate for concentration-dependent protein-binding saturation in Models VI (TDD simple exponent) and VII (dose-dependent Emax, with Emax and Hill fixed at 2.8 and 1.68 from Ding 2015). Model VII did lower the OFV to 7037.671, below the retained Model II, but was rejected because its typical CL/F of 0.052 L/h was 'unreasonable'. The Discussion further cautions that incorporating TDD in a therapeutic-drug-monitoring dataset is confounded by the TDM feedback effect. No nonlinear protein-binding term appears in the final model.",
      source_name = "TDD"
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened by forward inclusion as a continuous covariate on CL/F in linear, power and exponential forms (Methods Equations 23-25) and not retained. Albumin also enters Models VIII-X, the three explicit protein-binding models (one-binding-site with K = 15.5 1/mM and N = 1.98; Langmuir with Kd = 7.8 and Bm = 130 mg/L; linear non-saturable with Kd = 2.12, Bm = 67.3 mg/L and NS = 2.25), none of which improved on the linear Model II. Per-patient values are in Table S1, which is not reproduced in the main text."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened by forward inclusion as a hepatic-function covariate on CL/F and not retained (Wang 2024 Methods, Covariate model). Per-patient values are in Table S1, which is not reproduced in the main text."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened by forward inclusion as a hepatic-function covariate on CL/F and not retained (Wang 2024 Methods, Covariate model). Per-patient values are in Table S1, which is not reproduced in the main text."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units       = "umol/L",
      type        = "continuous",
      notes       = "A negative correlation with CL/F was seen during covariate evaluation but creatinine was NOT retained. The Discussion rejects it on physiological grounds: the patients' creatinine was almost entirely within the age-appropriate reference range, valproic acid is rarely excreted unchanged in urine, and creatinine 'is not a reliable surrogate for estimating the renal function in children' because it tracks muscle mass (Figure S4)."
    ),
    CYSC = list(
      description = "Serum cystatin C",
      units       = "mg/L",
      type        = "continuous",
      notes       = "Used in place of creatinine to derive the estimated glomerular filtration rate, because the study recorded no height and so could not apply a height-based paediatric equation. Not itself retained on CL/F.",
      source_name = "CysC"
    ),
    CRCL = list(
      description = "Estimated glomerular filtration rate",
      units       = "mL/min/1.73 m^2",
      type        = "continuous",
      notes       = "Derived three ways (Methods Equations 20-22: the KDIGO 2012 cystatin-C equation, Shull's equation, and a further cystatin-C-based equation) and screened on CL/F. 'Renal function did not impact VPA's CL/F' (Results, Covariate model); the Discussion adds that the KDIGO 2012 cystatin-C eGFR 'did not fulfill the inclusion/exclusion criteria'. The canonical column CRCL is the register's home for a renal-function covariate reported as an eGFR in mL/min/1.73 m^2; this model does not use it.",
      source_name = "eGFR"
    ),
    CONMED_OXC = list(
      description = "Concomitant oxcarbazepine indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "39 of 498 model-development patients (7.8%). Reached forward-inclusion significance on CL/F with an OFV drop of 7.981 (p < 0.05) but was removed in backward elimination (dOFV < 10.83). The Discussion attributes the non-retention to the small number of co-treated subjects.",
      source_name = "OXC"
    ),
    CONMED_LAMOTRIGINE = list(
      description = "Concomitant lamotrigine indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "55 of 498 model-development patients (11.0%). Forward-inclusion OFV drop of 8.329 (p < 0.05) on CL/F; removed in backward elimination (dOFV < 10.83).",
      source_name = "LTG"
    ),
    CONMED_CZP = list(
      description = "Concomitant clonazepam indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "39 of 498 model-development patients (7.8%). Forward-inclusion OFV drop of 7.049 (p < 0.05) on CL/F; removed in backward elimination (dOFV < 10.83).",
      source_name = "CZP"
    ),
    CONMED_PB = list(
      description = "Concomitant phenobarbital indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Only 1 of 498 model-development patients (0.2%), so phenobarbital fell below the paper's own '>5% of patients' screening threshold. The Discussion nonetheless reports a forward-inclusion OFV drop of 8.602 for concurrent phenobarbital therapy; it was removed in backward elimination.",
      source_name = "PB"
    ),
    CONMED_TPM = list(
      description = "Concomitant topiramate indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "51 of 498 model-development patients (10.2%), so topiramate passed the '>5% of patients' threshold and was screened as a covariate on CL/F, but the paper reports no OFV drop for it and it is not retained. 'In this study, little impact of concurrent ASMs on VPA's CL/F was observed.'",
      source_name = "TPM"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 498,
    n_studies      = 1,
    n_observations = 1138,
    age_range      = "0.2-17.2 years",
    age_median     = "5.8 years",
    weight_range   = "6.0-95.0 kg",
    weight_median  = "22.0 kg",
    sex_female_pct = 39.0,
    race_ethnicity = "Chinese (single-centre Nanjing cohort; sub-ethnicity not reported)",
    disease_state  = "Epilepsy, treated with valproic acid for at least 1 week before sampling",
    dose_range     = "7.7-40.0 mg/kg/day (median 21.7), given one to three times daily as oral syrup or sustained-release tablet",
    regions        = "China (Children's Hospital of Nanjing Medical University; January 2022 - March 2023)",
    co_medication  = "Concomitant antiseizure medications in the model-development set: levetiracetam 102 (20.5%), lamotrigine 55 (11.0%), topiramate 51 (10.2%), clonazepam 39 (7.8%), oxcarbazepine 39 (7.8%), perampanel 35 (7.0%), lacosamide 23 (4.6%), vigabatrin 5 (1.0%), zonisamide 5 (1.0%), phenobarbital 1 (0.2%). Only co-medications used by more than 5% of patients were screened as covariates, and none was retained in the final model. Patients taking any non-antiseizure co-medication were excluded at enrolment.",
    notes          = "1411 steady-state trough concentrations from 617 children were split 8:2 into a model-development set (498 patients, 1138 samples) and an external-evaluation set (119 patients, 273 samples); the numbers recorded in this metadata block are the MODEL-DEVELOPMENT set, which is what Table 3 was fitted to. External-evaluation set: 119 patients (78 male / 41 female), age median 6.2 years (0.5-13.8), body weight median 23.0 kg (6.0-75.0). Observed valproic acid concentrations 12.6-143.1 mg/L (median 68.6), all routine therapeutic-drug-monitoring TROUGHS - there is no absorption-phase sampling, which is why both Ka values are fixed from the literature and why interindividual variability on V/F 'was not informative enough to be calculated' and was dropped. Assay: enzyme multiplied immunoassay technique (EMIT; Viva-E, Siemens), calibration range 1.00-150 mg/L, within- and between-run coefficients of variation below 15%. Median 2 samples per patient (range 1-9). Model fitted in NONMEM 7.3.0 with PsN 5.2.6 using FOCE-I and the ADVAN2 TRANS2 subroutines. Baseline demographics: Wang 2024 Table 2."
  )

  ini({
    # ----------------------------------------------------------------
    # Absorption - FIXED to the literature Ka pair the authors adopted.
    # Wang 2024 Methods, Base model: "The fixed values for ka were 2.64
    # and 0.46 h-1 for syrup and SR tablets, respectively, due to a lack
    # of absorption phase data." Oral syrup is the reference
    # formulation, so lka is the syrup value and the SR indicator
    # carries the log-ratio shift (the same pattern as the sibling
    # Zhang_2024_valproic_acid.R, which fixes the identical pair).
    # ----------------------------------------------------------------
    lka <- fixed(log(2.64)); label("Absorption rate constant, oral syrup reference (1/h)")                                # Wang 2024 Methods Base model and Table 3 (Ka1 = 2.64, Fixed)
    e_form_vpa_sr_ka <- fixed(log(0.46 / 2.64)); label("Log-ratio shift on Ka for sustained-release tablet vs oral syrup") # Wang 2024 Methods Base model and Table 3 (Ka2 = 0.46, Fixed)

    # ----------------------------------------------------------------
    # Structural parameters, final model. Both are apparent (oral)
    # parameters describing TOTAL plasma valproic acid: "as
    # bioavailability (F) could not be determined, the CL and the Vd
    # were considered as the apparent CL (CL/F) and Vd (Vd/F)".
    #
    # Cross-check on the reference weight: at BW = 22 kg these give
    # CL/F = 0.196 / 22 = 0.0089 L/h/kg and V/F = 2.09 / 22 = 0.095
    # L/kg, which is exactly the paper's own summary immediately below
    # Equation 35 ("the typical value for the CL/F of VPA in children
    # at steady state obtained in our study was 0.009 L/h/kg and the
    # V/F was 0.1 L/kg"). This confirms 22 kg is the normalising
    # weight, which Equations 34-35 write only as the literal "22".
    # ----------------------------------------------------------------
    lcl <- log(0.196); label("Apparent clearance CL/F for a 22 kg male (L/h)")            # Wang 2024 Table 3 final model (CL/F = 0.196 L/h, RSE 3%, bootstrap 95% CI [0.183, 0.208]) and Equation 34
    lvc <- log(2.09);  label("Apparent volume of distribution V/F at 22 kg (L)")          # Wang 2024 Table 3 final model (V/F = 2.09 L, RSE 9%, bootstrap 95% CI [1.69, 2.49]) and Equation 35

    # ----------------------------------------------------------------
    # Covariate effects. Wang 2024 Equations 34 and 35:
    #   CL/F = 0.196 * (BW/22)^0.717 * (1 - 0.0436 * Sex)
    #   V/F  = 2.09  * (BW/22)^0.524
    # with "Sex = 0 for male, Sex = 1 for female".
    #
    # NOTE the sex term is the paper's PROPORTIONAL categorical form
    # (Methods Equation 27, Pi = TV(P) * (1 + theta * COV)), NOT the
    # exponential form (Equation 28). It is therefore applied as a
    # multiplicative (1 + e_sexf_cl * SEXF) factor in model() rather
    # than folded into the log-scale sum, so the encoded value is
    # exactly the tabulated -0.0436 with no back-transformation.
    # ----------------------------------------------------------------
    e_wt_cl   <-  0.717;   label("Power exponent on (WT/22) for CL/F (unitless)")         # Wang 2024 Table 3 "BW on CL/F" 0.717 (RSE 4%, bootstrap 95% CI [0.659, 0.775]); Equation 34
    e_wt_vc   <-  0.524;   label("Power exponent on (WT/22) for V/F (unitless)")          # Wang 2024 Table 3 "BW on V/F" 0.524 (RSE 15%, bootstrap 95% CI [0.361, 0.686]); Equation 35
    e_sexf_cl <- -0.0436;  label("Proportional shift on CL/F for female sex (unitless)")  # Wang 2024 Table 3 "Sex on CL/F" -0.0436 (RSE 29%, bootstrap 95% CI [-0.0685, -0.0187]); Equation 34

    # ----------------------------------------------------------------
    # IIV. Wang 2024 Table 3 reports a single "IIV on CL/F (%)" of
    # 11.9 (RSE 6%, bootstrap 95% CI [10.4, 13.2]) under the
    # exponential IIV model of Equation 1. The IIV on V/F "was not
    # informative enough to be calculated because only the Ctrough
    # samples of VPA were available, and then were excluded from the
    # model" (Results, Base model).
    #
    # The tabulated 11.9 is a CV%, so the internal variance is
    #   omega^2 = log(CV^2 + 1) = log(0.119^2 + 1) = 0.014062
    #
    # This is a load-bearing reading, because a variance-scale reading
    # (omega^2 = 0.119, i.e. CV = 35.5%) is also superficially
    # available. It is settled by the authors' own Figure 2, whose
    # boxes are "the median and IQR of Ctrough for each dosing
    # regimen" over 1000 simulated patients: the drawn boxes have
    # Q3/Q1 of roughly 1.2 (e.g. the 20 kg male / 40 mg/kg/day syrup
    # box spans about 92-112 mg/L about a median near 102), which for
    # a log-normal implies a log-scale SD near 0.15 - consistent with
    # a 11.9% CV plus this model's residual error, and flatly
    # inconsistent with a 35.5% CV, which would put the same panel's
    # Tukey whiskers below 10 and above 190 mg/L instead of the drawn
    # 75-140 mg/L.
    #
    # The zero-variance V/F term is OMITTED rather than written as
    # `etalvc ~ fixed(0)`: a zero diagonal makes OMEGA singular and
    # breaks the Cholesky sampler used by rxSolve (same treatment as
    # the sibling Zhang_2024_valproic_acid.R).
    # ----------------------------------------------------------------
    etalcl ~ 0.014062  # Wang 2024 Table 3 final model (IIV on CL/F = 11.9%, exponential per Equation 1)

    # ----------------------------------------------------------------
    # Residual variability - combined proportional and additive. Wang
    # 2024 Results, Base model: "The combined additive and proportional
    # model (Equation 5) showed the best fit for characterizing the
    # RV". Table 3 tabulates RVprop = 0.0151 (RSE 15%) and RVadd = 22.6
    # (RSE 35%), which are NONMEM $SIGMA VARIANCES, so the standard
    # deviations this file encodes are
    #   propSd = sqrt(0.0151) = 0.12288      (12.3%)
    #   addSd  = sqrt(22.6)   = 4.7539 mg/L
    #
    # The variance reading is confirmed by the assay: combining the two
    # at the cohort median trough of 68.6 mg/L gives
    #   sqrt((0.12288 * 68.6)^2 + 4.7539^2) / 68.6 = 14.1%
    # against an EMIT assay whose within- and between-run coefficients
    # of variation are "below 15%" (Methods, Bioassay of VPA). Reading
    # the tabulated numbers as standard deviations instead would put
    # the additive term alone at 22.6 mg/L, i.e. 33% of the median
    # trough, which no assay-limited TDM dataset supports.
    #
    # Equation 5 is TYPESET as Y = IPRED * (1 + eps1 + eps2), which
    # would make BOTH error terms proportional. That contradicts Table
    # 3's own row labels ("RVprop" and "RVadd") and is arithmetically
    # impossible for the tabulated magnitude - a proportional variance
    # of 22.6 is a 475% CV. It is a typesetting slip for the standard
    # NONMEM combined form Y = IPRED * (1 + eps1) + eps2, which is what
    # is encoded here.
    # ----------------------------------------------------------------
    propSd <- 0.12288; label("Proportional residual SD (fraction)")   # Wang 2024 Table 3 final model (RVprop = 0.0151 variance, RSE 15%, bootstrap 95% CI [0.0105, 0.0198])
    addSd  <- 4.7539;  label("Additive residual SD (mg/L)")           # Wang 2024 Table 3 final model (RVadd = 22.6 variance, RSE 35%, bootstrap 95% CI [6.58, 38.5])
  })

  model({
    # 1. Formulation-specific absorption rate constant. Oral syrup is
    #    the reference (FORM_VPA_SR = 0); both values are FIXED.
    ka <- exp(lka + e_form_vpa_sr_ka * FORM_VPA_SR)

    # 2. Apparent clearance, Equation 34. The body-weight power term is
    #    written on the log scale (exponent * log(ratio)), which is
    #    algebraically identical to the paper's (BW/22)^0.717 form; the
    #    sex term stays OUTSIDE the exponential because the paper uses
    #    the proportional categorical model (1 - 0.0436 * Sex), not the
    #    exponential one.
    cl <- exp(lcl + e_wt_cl * log(WT / 22) + etalcl) * (1 + e_sexf_cl * SEXF)

    # 3. Apparent volume of distribution, Equation 35. No
    #    interindividual variability (not estimable from trough-only
    #    data; see ini()).
    vc <- exp(lvc + e_wt_vc * log(WT / 22))

    # 4. Micro-constant
    kel <- cl / vc

    # 5. One-compartment ODE system with first-order oral absorption
    #    (NONMEM ADVAN2 TRANS2)
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # 6. Observation (total plasma valproic acid) and residual error
    Cc <- central / vc
    Cc ~ prop(propSd) + add(addSd)
  })
}
