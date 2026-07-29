Chung_2013_vancomycin <- function() {
  description <- "One-compartment IV-infusion population PK model for vancomycin in Korean adults with normal serum creatinine (Chung 2013). CL and V are described by centered-linear additive deviations on age, total body weight, serum creatinine (CL only), and sex, plus a power-law effect of serum cystatin C on CL (reference 0.91 mg/L, exponent -0.78); cystatin C is the dominant CL covariate, accounting for ~62% of the inter-individual variability in CL even within the SCr <= 1.2 mg/dL inclusion window."
  reference <- "Chung JY, Jin SJ, Yoon JH, Song YG. Serum cystatin C is a major predictor of vancomycin clearance in a population pharmacokinetic analysis of patients with normal serum creatinine concentrations. J Korean Med Sci. 2013;28(1):48-54. doi:10.3346/jkms.2013.28.1.48"
  vignette <- "Chung_2013_vancomycin"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Chung 2013 Table 1: median 57 years (range 18-96). Centered at 57 in the structural model (Table 2 footnote: CL = CL_POP * (1 + theta_CLage * (AGE - 57)) * ...; V = V_POP * (1 + theta_Vage * (AGE - 57)) * ...). Linear-deviation form for both CL and V.",
      source_name        = "AGE"
    ),
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Chung 2013 Table 1: median 60.8 kg (range 27-140). Source paper symbol TBW (total body weight); stored under the canonical WT column. Centered at 60.8 in the structural model (Table 2 footnote: CL = ... * (1 + theta_CLTBW * (TBW - 60.8)) * ...; V = ... * (1 + theta_VTBW * (TBW - 60.8)) * ...). Linear-deviation form for both CL and V. Chung 2013 Discussion notes ideal/lean body weight did not improve fit and actual TBW is recommended for vancomycin TDM.",
      source_name        = "TBW"
    ),
    CREAT = list(
      description        = "Serum creatinine concentration (Jaffe kinetic method, Hitachi 7600)",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Chung 2013 Table 1: median 0.9 mg/dL (mean 0.80, range 0.39-1.2). Inclusion criterion SCr <= 1.2 mg/dL. Centered at 0.8 mg/dL in the structural model (Table 2 footnote: CL = ... * (1 + theta_CLSCr * (SCr - 0.8)) * ...). Linear-deviation form on CL only (no SCr effect on V). Source paper symbol SCr; stored under canonical CREAT (mg/dL).",
      source_name        = "SCr"
    ),
    CYSC = list(
      description        = "Serum cystatin C concentration (particle-enhanced immunoturbidimetric assay, Roche Cobas 6000)",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Chung 2013 Table 1: median 0.91 mg/L (mean 1.01, range 0.38-3.1). Reference ranges 0.57-0.97 mg/L (female) / 0.65-1.10 mg/L (male). Enters CL with a power-law form (Table 2 footnote: CL = ... * (Cystatin C / 0.91)^theta_CLcystatin) with reference 0.91 mg/L and exponent -0.780; described in Chung 2013 Results as the most influential covariate, reducing OFV by 428.3 and accounting for ~62% of the CL variability. No effect on V.",
      source_name        = "Cystatin C"
    ),
    SEXF = list(
      description        = "Biological sex indicator (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Chung 2013 Table 1: 400 male / 278 female (41% female). Enters both CL and V as a female-only additive multiplicative shift (Table 2 footnote: 'if female, apply 1 + theta_CLsex' on CL and '1 + theta_Vsex' on V), encoded as (1 + e * SEXF) so SEXF = 0 (male) recovers the typical-value reference. Chung 2013 Discussion: 'Females have 15% lower clearance, which is in line with the Cockcroft and Gault equation.'",
      source_name        = "SEX"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 678L,
    n_studies        = 1L,
    age_range        = "18-96 years",
    age_median       = "57 years",
    weight_range     = "27-140 kg",
    weight_median    = "60.8 kg",
    sex_female_pct   = 41.0,
    race_ethnicity   = "Korean (single-center cohort at Gangnam Severance Hospital, Seoul)",
    disease_state    = "Adult inpatients receiving intravenous vancomycin therapy with serum creatinine <= 1.2 mg/dL. 73% had ICU admission; 7% co-administered amikacin; 17% co-administered furosemide. Vancomycin concentrations measured for therapeutic drug monitoring after at least three doses (steady-state in most subjects).",
    dose_range       = "Intravenous vancomycin 250-4,500 mg/day (median 2,000 mg/day; 89% in 1,000-2,000 mg/day range); 1,000 mg dose dissolved in 100 mL saline and infused over 1 hour.",
    regions          = "Korea (Gangnam Severance Hospital, Seoul)",
    renal_function   = "Cystatin C median 0.91 mg/L (range 0.38-3.1); serum creatinine median 0.9 mg/dL (range 0.39-1.2 by inclusion criterion); BUN median 12.9 mg/dL (range 0-80).",
    n_concentrations = 1373L,
    notes            = "Patient characteristics from Chung 2013 Table 1 (June 2006 to May 2010). Sampling design: trough just before next infusion and peak 1 h after end of infusion. 2,000 mg/day was the median daily dose. Modeling done in NONMEM 7.1.0 with FOCE-INTER. Final model retained age, TBW, cystatin C, SCr, and sex on CL and age, TBW, and sex on V; the Cockcroft-Gault CLcr was rejected as redundant. Shrinkage: EBE_CL 5.8%, EBE_V 17.3%, IWRES 48.3%. Bootstrap (1000 replicates) confirmed parameter robustness."
  )

  ini({
    # Structural parameters (Chung 2013 Table 2 final-model column "Estimate (RSE)";
    # reference subject has AGE = 57, TBW = 60.8 kg, SCr = 0.8 mg/dL, CysC = 0.91 mg/L,
    # SEX = male).
    lcl <- log(4.90); label("Typical CL at reference covariates (L/h)") # Chung 2013 Table 2: CL_POP = 4.90 L/h (RSE 1.5%)
    lvc <- log(46.2); label("Typical V at reference covariates (L)")    # Chung 2013 Table 2: V_POP  = 46.2 L  (RSE 1.6%)

    # Covariate effects on CL (Chung 2013 Table 2 final-model values; functional
    # forms from the Table 2 footnote):
    #   CL = CL_POP * (1 + e_age_cl * (AGE - 57))
    #              * (1 + e_wt_cl  * (WT  - 60.8))
    #              * (1 + e_creat_cl * (CREAT - 0.8))
    #              * (CYSC / 0.91)^e_cysc_cl
    #              * (1 + e_sexf_cl * SEXF)
    e_age_cl   <- -0.00420; label("CL slope on (AGE - 57) (1/year)")              # Chung 2013 Table 2: theta_CLage     = -0.00420 (RSE 18.2%)
    e_wt_cl    <-  0.00997; label("CL slope on (WT - 60.8) (1/kg)")               # Chung 2013 Table 2: theta_CLTBW     =  0.00997 (RSE 8.8%)
    e_creat_cl <- -0.322;   label("CL slope on (CREAT - 0.8) (dL/mg)")            # Chung 2013 Table 2: theta_CLSCr     = -0.322   (RSE 17.1%)
    e_cysc_cl  <- -0.780;   label("CL power exponent on (CYSC / 0.91) (unitless)") # Chung 2013 Table 2: theta_CLcystatin = -0.780  (RSE 5.3%)
    e_sexf_cl  <- -0.150;   label("CL multiplicative shift in females (unitless)") # Chung 2013 Table 2: theta_CLsex     = -0.150   (RSE 13.5%)

    # Covariate effects on V (Chung 2013 Table 2; Table 2 footnote):
    #   V = V_POP * (1 + e_age_vc * (AGE - 57))
    #             * (1 + e_wt_vc  * (WT  - 60.8))
    #             * (1 + e_sexf_vc * SEXF)
    e_age_vc  <-  0.00580; label("V slope on (AGE - 57) (1/year)")     # Chung 2013 Table 2: theta_Vage  =  0.00580 (RSE 12.3%)
    e_wt_vc   <-  0.00661; label("V slope on (WT - 60.8) (1/kg)")      # Chung 2013 Table 2: theta_VTBW  =  0.00661 (RSE 15.9%)
    e_sexf_vc <- -0.119;   label("V multiplicative shift in females (unitless)") # Chung 2013 Table 2: theta_Vsex  = -0.119   (RSE 20.5%)

    # Inter-individual variability (Chung 2013 Table 2 final-model column reports
    # ISV as %CV; for log-normal etas omega^2 = log(CV^2 + 1)).
    etalcl ~ 0.05922 # log(0.247^2 + 1); 24.7% CV on CL (RSE 26.2%)
    etalvc ~ 0.06110 # log(0.251^2 + 1); 25.1% CV on V  (RSE 37.3%)

    # Combined additive + proportional residual error (Chung 2013 Table 2;
    # additive SD reported in mg/L, proportional reported as %CV).
    addSd  <- 1.40;   label("Additive residual error (mg/L)")        # Chung 2013 Table 2: additive    = 1.40 mg/L (RSE 12.9%)
    propSd <- 0.0639; label("Proportional residual error (fraction)") # Chung 2013 Table 2: proportional = 6.39% CV (RSE 22.8%)
  })
  model({
    # Individual PK parameters with centered-linear covariate effects on AGE, WT,
    # and CREAT (CL only), a power-law effect of CYSC on CL, and a female-only
    # multiplicative shift on both CL and V. Encoding follows Table 2 footnote
    # verbatim with the (1 + e * cov) form so the typical-value subject (AGE = 57,
    # WT = 60.8, CREAT = 0.8, CYSC = 0.91, SEXF = 0) recovers exp(lcl) and exp(lvc).
    cl <- exp(lcl + etalcl) *
          (1 + e_age_cl   * (AGE   - 57))   *
          (1 + e_wt_cl    * (WT    - 60.8)) *
          (1 + e_creat_cl * (CREAT - 0.8))  *
          (CYSC / 0.91)^e_cysc_cl           *
          (1 + e_sexf_cl  * SEXF)
    vc <- exp(lvc + etalvc) *
          (1 + e_age_vc * (AGE - 57))   *
          (1 + e_wt_vc  * (WT  - 60.8)) *
          (1 + e_sexf_vc * SEXF)

    kel <- cl / vc

    d/dt(central) <- -kel * central

    # Dose in mg, vc in L -> central/vc has units mg/L.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
