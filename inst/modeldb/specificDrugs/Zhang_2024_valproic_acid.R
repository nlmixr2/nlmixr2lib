Zhang_2024_valproic_acid <- function() {
  description <- "One-compartment population PK model with first-order absorption for total plasma valproic acid in Chinese children and adults with epilepsy or after neurosurgery (Zhang 2024 final model). Apparent clearance carries body weight, serum creatinine, serum albumin, sex, concomitant carbapenem and concomitant enzyme-inducing antiepileptic drug effects; apparent volume carries body weight. Carbapenem coadministration multiplies CL/F by exp(1.50) = 4.48, the first population PK model to quantify this interaction. Formulation-specific absorption rate constants are FIXED from the literature (oral solution 2.64 1/h reference, sustained-release tablet 0.46 1/h)."
  reference <- "Zhang L, Wu R, Li X, Feng W, Zhao Z, Mei S. Combined carbapenem resulted in a 4.48-fold increase in valproic acid clearance: a population pharmacokinetic model in Chinese children and adults with epilepsy or after neurosurgery. Front Pharmacol. 2024;15:1423411. doi:10.3389/fphar.2024.1423411. PMCID PMC11581887. Final-model parameter estimates from Table 3 and the covariate equations in the Abstract and Results section 3.2."
  vignette <- "Zhang_2024_valproic_acid"
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
      reference_category = "60 kg (the cohort median; Zhang 2024 Results 3.2)",
      notes              = "Power effect on both CL/F (exponent 0.787) and V/F (exponent 0.751), each normalised to the 60.0 kg cohort median. Cohort range 5.50-120.00 kg (Zhang 2024 Table 1). The V/F exponent is reported only in the Abstract equation and in Results 3.2 ('0.75 for Vd and BW'); it is absent from Table 3, which tabulates only the CL/F exponent.",
      source_name        = "BW"
    ),
    CREAT = list(
      description        = "Serum creatinine",
      units              = "umol/L",
      type               = "continuous",
      reference_category = "50.3 umol/L (the cohort median; Zhang 2024 Results 3.2)",
      notes              = "Power effect on CL/F with a NEGATIVE exponent (-0.253): valproic acid clearance falls as creatinine rises. Cohort range 13.00-447.60 umol/L, with 33.86% of patients outside the 44.0-132 umol/L normal range (Zhang 2024 Discussion). Reported in SI units (umol/L), not mg/dL.",
      source_name        = "Cr"
    ),
    ALB = list(
      description        = "Serum albumin",
      units              = "g/L",
      type               = "continuous",
      reference_category = "39 g/L (the cohort median; Zhang 2024 Results 3.2)",
      notes              = "Power effect on CL/F with a NEGATIVE exponent (-0.873): apparent clearance of TOTAL valproic acid falls as albumin rises. Valproic acid is 90-95% albumin-bound, so a lower albumin leaves a larger unbound fraction and raises total-drug apparent clearance (Zhang 2024 Discussion). Cohort range 23.5-51.80 g/L, with 21.90% below the 28-54 g/L reference range. Reported in SI units (g/L).",
      source_name        = "ALB"
    ),
    SEXF = list(
      description        = "Female sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Exponential effect on CL/F: exp(0.121) = 1.129, i.e. women have 12.9% higher apparent clearance than men (Zhang 2024 Discussion states '12.9% higher'). 290 of 443 patients were female (Zhang 2024 Table 1). Note this is the OPPOSITE direction to several prior valproate studies, which the Discussion attributes to the larger median body weight of the women in this cohort (65 kg vs 57 kg).",
      source_name        = "gender"
    ),
    CONMED_CARBAPENEM = list(
      description        = "Concomitant carbapenem antibiotic indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant carbapenem)",
      notes              = "Exponential effect on CL/F: exp(1.50) = 4.482, a 448.2% increase, which is this paper's headline finding and the first time carbapenem coadministration has been retained as a covariate in a valproate population PK model. 22 of 443 patients (4.9%) were coadministered a carbapenem: 18 meropenem (4.1%) and 3 ertapenem (0.6%); the paper pools both agents under the single class flag CBP. Mechanism per the Discussion: inhibition of acylpeptide hydrolase reduces deglucuronidation of valproate glucuronide.",
      source_name        = "CBP"
    ),
    CONMED_EIAED = list(
      description        = "Concomitant enzyme-inducing antiepileptic drug indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant enzyme-inducing antiepileptic drug)",
      notes              = "This paper's 'enzyme inducer two (IND2)' definition: 1 = the patient takes at least one of oxcarbazepine, carbamazepine, phenobarbital or phenytoin. Note that OXCARBAZEPINE IS INCLUDED, which is wider than the classic carbamazepine / phenobarbital / phenytoin EIAED triad used by Rodrigues_2017_oxcarbazepine.R; broadening the definition to capture oxcarbazepine's modest CYP3A4 induction is an explicit contribution of this paper (Introduction). Exponential effect on CL/F: exp(0.15) = 1.162, i.e. clearance rises to 116% of the non-induced value. 91 of 443 patients (20.5%) were IND2-positive. The paper also recorded the narrow 'enzyme inducer one (IND1)' flag (carbamazepine / phenobarbital / phenytoin only, 35 patients / 7.9%; Zhang 2024 Table 1) but the final model retains IND2, not IND1; IND1 is not represented in this model file because it would collide with this same canonical column.",
      source_name        = "IND2"
    ),
    FORM_VPA_SR = list(
      description        = "Sustained-release valproic acid tablet formulation indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (oral solution, the reference formulation in this cohort)",
      notes              = "Selects the FIXED sustained-release-tablet absorption rate constant Ka = 0.46 1/h; the oral-solution reference is Ka = 2.64 1/h. Both values were fixed from Ding 2015 because the therapeutic-drug-monitoring dataset is almost entirely steady-state troughs and contains no absorption-phase data (Zhang 2024 Methods 2.4.1). This cohort has only the two levels solution and sustained-release tablet, so FORM_TABLET (the conventional immediate-release level used by Zhang_2023_valproic_acid_base.R) is not part of this model.",
      source_name        = "Dosage form"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Recorded for every patient (median 32.74 years, range 0.27-84.38; Zhang 2024 Table 1) and used to define the typical-patient profiles simulated in Table 4, but NOT retained on CL/F or V/F in the final model. Table 2 reports the stepwise procedure only for the six retained covariates, so the paper does not state the delta-OFV at which age was rejected."
    ),
    HT = list(
      description = "Body height",
      units       = "cm",
      type        = "continuous",
      notes       = "Recorded for every patient (median 162.00 cm, range 54.00-185.00; Zhang 2024 Table 1) but not retained in the final model and not reported in Table 2's stepwise procedure. The paper's column is spelled out as 'Height'; the canonical column is HT.",
      source_name = "Height"
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Recorded for every patient (median 16.20 U/L, range 0.00-436.20; Zhang 2024 Table 1) but not retained on CL/F. The Discussion reports that 36.4% of the carbapenem-coadministered patients had ALT above the reference range and reads this as a consequence of the interaction rather than as a clearance covariate."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Recorded for every patient (median 29.10 U/L, range 3.10-752.10; Zhang 2024 Table 1) but not retained on CL/F and not reported in Table 2's stepwise procedure."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 443,
    n_studies      = 1,
    n_observations = 615,
    age_range      = "0.27-84.38 years",
    age_median     = "32.74 years",
    weight_range   = "5.50-120.00 kg",
    weight_median  = "60.00 kg",
    sex_female_pct = 65.5,
    race_ethnicity = "Chinese (two-centre Beijing cohort; sub-ethnicity not reported)",
    disease_state  = "Epilepsy or post-neurosurgery seizure prophylaxis; 65.5% adults and 34.5% children",
    dose_range     = "80-3000 mg/day (median 800); 2.79-109.09 mg/kg/day (median 14.81); oral solution or sustained-release tablet given once or twice daily",
    regions        = "China (Beijing Tiantan Hospital and Beijing Children's Hospital, Capital Medical University; October 2016 - June 2022)",
    renal_function = "Serum creatinine 13.00-447.60 umol/L (median 50.30); 33.86% of patients outside the 44.0-132 umol/L normal range",
    co_medication  = "Only 50.11% received valproate monotherapy. Levetiracetam 137 (30.9%), oxcarbazepine 63 (14.2%), lamotrigine 28 (6.3%), clonazepam 26 (5.8%), topiramate 25 (5.6%), phenobarbital 23 (5.1%), carbamazepine 15 (3.3%), lacosamide 13 (2.9%), phenytoin 2 (0.4%), nitrazepam 2 (0.4%). Carbapenems 22 (4.9%): meropenem 18 (4.1%), ertapenem 3 (0.6%).",
    notes          = "615 total plasma valproic acid concentrations in 443 patients (290 female / 153 male), measured by fluorescence polarization immunoassay (Centaur XP, Siemens; quantitative range 1-150 mg/L). Observed concentrations 0.00-165.05 mg/L (median 63.01); 10 records (1.6%) were below the limit of quantification. Most samples are steady-state troughs from routine therapeutic drug monitoring, with no absorption-phase sampling - this is why Ka is fixed rather than estimated and why interindividual variability on V/F could not be estimated. All patients had received valproate for at least one month before sampling. Model fitted in Phoenix NLME 8.3 with FOCE-ELS. Baseline demographics: Zhang 2024 Table 1."
  )

  ini({
    # ----------------------------------------------------------------
    # Absorption - FIXED to the literature Ka pair the authors adopted.
    # Zhang 2024 Methods 2.4.1: "Due to the lack of available data
    # during the absorption phase, Ka was fixed at 0.46 and 2.64 h-1
    # for the sustained tablets and solutions, respectively, according
    # to previous studies (Ding et al., 2015)". Oral solution is the
    # reference formulation, so lka is the solution value and the SR
    # indicator carries the log-ratio shift (the same pattern as the
    # sibling Zhang_2023_valproic_acid_base.R).
    # ----------------------------------------------------------------
    lka <- fixed(log(2.64)); label("Absorption rate constant, oral solution reference (1/h)")                     # Zhang 2024 Methods 2.4.1 and Table 3 (Ka solutions = 2.64 1/h, FIXED)
    e_form_vpa_sr_ka <- fixed(log(0.46 / 2.64)); label("Log-ratio shift on Ka for sustained-release tablet vs oral solution") # Zhang 2024 Methods 2.4.1 and Table 3 (Ka sustained tablets = 0.46 1/h, FIXED)

    # ----------------------------------------------------------------
    # Structural parameters, final model. Both are apparent (oral)
    # parameters describing TOTAL plasma valproic acid.
    # ----------------------------------------------------------------
    lcl <- log(0.430); label("Apparent clearance CL/F at the reference covariate set (L/h)")   # Zhang 2024 Abstract CL equation and Results 3.2 ("0.430 (L/h) is the typical value of CL"); Table 3 final model tabulates 0.43 (3.22% RSE)
    lvc <- log(8.66);  label("Apparent volume of distribution V/F at 60 kg (L)")               # Zhang 2024 Abstract Vd equation and Table 3 final model (Vd = 8.66 L, 9.18% RSE)

    # ----------------------------------------------------------------
    # Covariate effects on CL/F. Zhang 2024 Abstract:
    #   CL (L/h) = 0.430 * (BW/60)^0.787 * (Cr/50.3)^-0.253
    #              * (ALB/39)^-0.873 * e^gender * e^CBP * e^IND2
    # with gender = 0.121 (female), CBP = 1.50, IND2 = 0.15, each 0 at
    # the reference level. The three continuous terms are power
    # functions of the median-normalised covariate; the three binary
    # terms are exponentials of the coefficient. Table 3 rounds the
    # same estimates to two significant figures (0.79, 0.12, -0.25,
    # -0.87, 1.50, 0.15); the Abstract equation carries the extra
    # digits and is used here.
    # ----------------------------------------------------------------
    e_wt_cl    <-  0.787; label("Power exponent on (WT/60) for CL/F (unitless)")               # Zhang 2024 Abstract CL equation (0.787); Table 3 "BW on CL" 0.79 (5.99% RSE, 95% CI 0.69-0.88)
    e_creat_cl <- -0.253; label("Power exponent on (CREAT/50.3) for CL/F (unitless)")          # Zhang 2024 Abstract CL equation (-0.253); Table 3 "Cr on CL" -0.25 (23.50% RSE, 95% CI -0.37 to -0.14)
    e_alb_cl   <- -0.873; label("Power exponent on (ALB/39) for CL/F (unitless)")              # Zhang 2024 Abstract CL equation (-0.873); Table 3 "ALB on CL" -0.87 (14.53% RSE, 95% CI -1.12 to -0.62)
    e_sexf_cl  <-  0.121; label("Log-scale shift on CL/F for female sex (unitless)")           # Zhang 2024 Abstract CL equation (gender = 0.121 when female); Table 3 "Gender on CL" 0.12 (30.25% RSE, 95% CI 0.049-0.19)
    e_conmed_carbapenem_cl <- 1.50; label("Log-scale shift on CL/F for concomitant carbapenem (unitless)")   # Zhang 2024 Abstract CL equation (CBP = 1.50 when combined with carbapenems); Table 3 "Mem on CL" 1.50 (9.79% RSE, 95% CI 1.21-1.79)
    e_conmed_eiaed_cl      <- 0.15; label("Log-scale shift on CL/F for concomitant enzyme-inducing antiepileptic (unitless)") # Zhang 2024 Abstract CL equation (IND2 = 0.15); Table 3 "Enzyme inducer 2 on CL" 0.15 (27.81% RSE, 95% CI 0.067-0.23)

    # ----------------------------------------------------------------
    # Covariate effect on V/F. Zhang 2024 Abstract:
    #   Vd (L) = 8.66 * (BW/60)^0.751
    # Results 3.2 quotes the same exponent as "0.75 for Vd and BW".
    # Table 3 does NOT tabulate this exponent (it lists only "BW on
    # CL"), so no RSE or confidence interval is available for it.
    # ----------------------------------------------------------------
    e_wt_vc <- 0.751; label("Power exponent on (WT/60) for V/F (unitless)")                    # Zhang 2024 Abstract Vd equation (0.751); Results 3.2 "0.75 for Vd and BW"; not tabulated in Table 3

    # ----------------------------------------------------------------
    # IIV. Zhang 2024 Table 3 reports a single exponential
    # interindividual variability of 23.33% for the final model, which
    # is the CL/F term - Results 3.2 states that interindividual
    # variability on Vd "needs to be fixed at zero" because sparse
    # trough-only data produced large shrinkage (Savic and Karlsson
    # 2009). Reported as a CV%, so the internal variance is
    #   omega^2 = log(CV^2 + 1) = log(0.2333^2 + 1) = 0.05300
    #
    # The zero-variance V/F term is OMITTED rather than written as
    # `etalvc ~ fixed(0)`: a zero diagonal makes OMEGA singular and
    # breaks the Cholesky sampler used by rxSolve (same treatment as
    # Wattanakul_2024_primaquine_motherinfant.R).
    # ----------------------------------------------------------------
    etalcl ~ 0.05300  # Zhang 2024 Table 3 final model (IIV exponential = 23.33%)

    # ----------------------------------------------------------------
    # Residual error - additive only. Zhang 2024 Results 3.2: "Random
    # residual variability was best described by the additive model".
    # Table 3 reports sigma (additive) = 17.68 with 6.63% RSE; the
    # Discussion confirms the scale and units by comparing "residual
    # variability was 17.68 mg/L" against prior studies' "3.11-17.3
    # mg/L", so the tabulated number is a standard deviation in mg/L,
    # not a variance.
    # ----------------------------------------------------------------
    addSd <- 17.68; label("Additive residual SD (mg/L)")                                       # Zhang 2024 Table 3 final model (sigma additive = 17.68, 6.63% RSE, 95% CI 15.38-19.98); units confirmed in the Discussion
  })

  model({
    # 1. Formulation-specific absorption rate constant. Oral solution
    #    is the reference (FORM_VPA_SR = 0); both values are FIXED.
    ka <- exp(lka + e_form_vpa_sr_ka * FORM_VPA_SR)

    # 2. Apparent clearance. Written on the log scale, which is
    #    algebraically identical to the paper's product form:
    #    the three power terms become exponent * log(ratio) and the
    #    three binary terms are already log-scale shifts.
    cl <- exp(lcl +
                e_wt_cl * log(WT / 60) +
                e_creat_cl * log(CREAT / 50.3) +
                e_alb_cl * log(ALB / 39) +
                e_sexf_cl * SEXF +
                e_conmed_carbapenem_cl * CONMED_CARBAPENEM +
                e_conmed_eiaed_cl * CONMED_EIAED +
                etalcl)

    # 3. Apparent volume of distribution. No interindividual
    #    variability (fixed at zero by the authors; see ini()).
    vc <- exp(lvc + e_wt_vc * log(WT / 60))

    # 4. Micro-constant
    kel <- cl / vc

    # 5. One-compartment ODE system with first-order oral absorption
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # 6. Observation (total plasma valproic acid) and residual error
    Cc <- central / vc
    Cc ~ add(addSd)
  })
}
