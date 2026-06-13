Tsuda_2010_tamsulosin <- function() {
  description <- "One-compartment population PK model for oral modified-release tamsulosin hydrochloride in paediatric patients (2-16 years) with neuropathic and non-neuropathic bladder (Tsuda 2010), with first-order absorption after a lag time, allometric (WT/70)^0.75 on apparent clearance and (WT/70)^1 on apparent central volume (allometric exponents fixed at theory values), a power-form alpha-1-acid glycoprotein (AAG/20 uM) effect on both CL/F and V/F, correlated inter-individual variability on CL/F and V/F, independent IIV on ka, and a combined additive + proportional residual error."
  reference <- "Tsuda Y, Tatami S, Yamamura N, Tadayasu Y, Sarashina A, Liesenfeld K-H, Staab A, Schaefer H-G, Ieiri I, Higuchi S. Population pharmacokinetics of tamsulosin hydrochloride in paediatric patients with neuropathic and non-neuropathic bladder. British Journal of Clinical Pharmacology. 2010;70(1):88-101. doi:10.1111/j.1365-2125.2010.03662.x"
  vignette <- "Tsuda_2010_tamsulosin"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric scaling on apparent clearance (exponent 0.75) and apparent volume (exponent 1) with reference 70 kg, per Tsuda 2010 Table 3 structural-model equations. Allometric exponents were estimated in the final model and reported at 0.768 (95% CI 0.598-0.938) for CL/F and 0.807 (95% CI 0.603-1.01) for V/F (Results 'Model building' Step 2); the 95% CIs included the theory values 0.75 and 1, so those were used in the final model. Cohort median 25.0 kg (range 12.1-92.2 kg, Table 2 step 2).",
      source_name        = "WT"
    ),
    AAG = list(
      description        = "Serum alpha-1-acid glycoprotein concentration",
      units              = "uM",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form covariate on both CL/F and V/F, applied as (AAG / 20.0)^theta with reference 20 uM (Tsuda 2010 Table 3 structural-model equation; reference value chosen by the authors close to the cohort median 19.95 uM in step 2 / 18.25 uM in step 1, Table 2). Cohort range 9.76-74.53 uM. The canonical-register AAG entry is in g/L (Kloft 2006 / Netterberg 2017 family); this model overrides per-model to uM so that the published exponents -0.844 (CL/F) and -0.663 (V/F) apply directly. To convert from g/L to uM use uM = (g/L) * 1000 / 41000 with the conventional 41 kDa AAG molecular weight (i.e. 20 uM ~ 0.82 g/L).",
      source_name        = "AAG"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Subject age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened in the stepwise covariate analysis (Tsuda 2010 Methods 'Step 2'); not retained in the final model. Authors interpret the null result as evidence that maturation of CYP3A4 / CYP2D6 and AAG concentration is complete before age 2 years (Discussion paragraph 5). Range 2-16 years."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened as 'gender' in the covariate analysis (Tsuda 2010 Methods 'Step 2'); not retained in the final model. Cohort 44% female (84 of 189 in step 2)."
    ),
    HT = list(
      description = "Height",
      units       = "cm",
      type        = "continuous",
      notes       = "Screened in the covariate analysis (Tsuda 2010 Methods 'Step 2'); not retained in the final model. Cohort median 125 cm (range 63-175)."
    ),
    BMI = list(
      description = "Body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Tested on V/F during forward inclusion (Tsuda 2010 Results 'Step 2': BMI on V/F entered the full model with OBJF drop -5.583) but was removed at backward elimination. Cohort median 17.41 kg/m^2 (range 10.65-36.32)."
    ),
    BSA = list(
      description = "Body surface area",
      units       = "m^2",
      type        = "continuous",
      notes       = "Screened in the covariate analysis (Tsuda 2010 Methods 'Step 2'); not retained in the final model. Cohort median 0.94 m^2 (range 0.43-2.03)."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units       = "mg/dL",
      type        = "continuous",
      notes       = "Screened in the covariate analysis (Tsuda 2010 Methods 'Step 2'); not retained in the final model. Cohort median 0.96 mg/dL (range 0.60-7.00)."
    ),
    CRCL = list(
      description = "Creatinine clearance (Cockcroft-Gault)",
      units       = "mL/min",
      type        = "continuous",
      notes       = "Screened in the covariate analysis (Tsuda 2010 Methods 'Step 2'); not retained in the final model. Cohort median 46.85 mL/min (range 5.76-197.5). The drug is largely metabolised (9-15% renal as parent per Discussion paragraph 1), so renal-function covariates were not expected to be primary drivers."
    ),
    GGT = list(
      description = "Gamma-glutamyl transferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Tested on CL/F during forward inclusion (Tsuda 2010 Results 'Step 2': GGT on CL/F entered the full model with OBJF drop -6.617) but was removed at backward elimination. Cohort median 29.25 U/L (range 10.71-130.0)."
    ),
    HGB = list(
      description = "Haemoglobin",
      units       = "g/dL",
      type        = "continuous",
      notes       = "Screened in the covariate analysis (Tsuda 2010 Methods 'Step 2'); not retained in the final model. Cohort median 14.29 g/dL (range 7.18-23.76)."
    ),
    ANTI_CHOLINERGIC = list(
      description = "Concomitant anti-cholinergic therapy indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened as a covariate of specific interest on CL/F and V/F (Tsuda 2010 Methods 'Step 2'); not retained in the final model. Prevalence 16.5% (percent of observation records) in step 2."
    ),
    NEUROGENIC = list(
      description = "Neurogenic-bladder patient indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened as a 'patient population' covariate of specific interest on CL/F and V/F (Tsuda 2010 Methods 'Step 2'); not retained in the final model. Cohort split: 45 non-neurogenic (trial 1) + 144 neurogenic (trials 2 and 3) = 189 in step 2."
    ),
    RACE = list(
      description = "Race / ethnic origin",
      units       = "(categorical)",
      type        = "categorical",
      notes       = "Screened as a covariate of specific interest on CL/F and V/F (Tsuda 2010 Methods 'Step 2'); not retained in the final model. Cohort distribution in step 2 (Table 2): White 87, Asian 83, Black or African American 6, American Indian or Alaska Native 13, Native Hawaiian or other Pacific Islander 0."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 189L,
    n_studies      = 3L,
    n_observations = 1082L,
    age_range      = "2-16 years (paediatric)",
    age_median     = "8 years",
    weight_range   = "12.1-92.2 kg",
    weight_median  = "25.0 kg",
    sex_female_pct = 44.4,
    race_ethnicity = "White 46%, Asian 44%, Black or African American 3%, American Indian or Alaska Native 7%, Native Hawaiian or other Pacific Islander 0% (step 2, n = 189; Tsuda 2010 Table 2).",
    disease_state  = "Paediatric patients with neuropathic (n = 144) and non-neuropathic (n = 45) bladder (Tsuda 2010 Table 2 step 2). Trial 1 non-neurogenic; Trials 2 and 3 neurogenic.",
    dose_range     = "Oral modified-release tamsulosin HCl capsule (sprinkled over applesauce or yogurt for paediatric administration in trials 2 and 3). Trial 1 single doses 0.1, 0.2, 0.4, or 0.8 mg by weight band. Trials 2 and 3 once-daily weight-based dosing at low (0.025-0.1 mg), medium (0.05-0.2 mg), or high (0.1-0.4 mg) per Table 1C.",
    regions        = "Multinational (Canada, USA, Germany, Philippines, Brazil, India, Italy, Korea, Mexico, Russia, South Africa, Spain, Ukraine; trial-investigator list in Appendix).",
    notes          = "Pooled analysis of three paediatric clinical trials (Table 1A): trial 1 (phase I single-dose PK / safety / tolerability in non-neurogenic bladder, 45 patients), trial 2 (phase II long-term safety / efficacy in neurogenic bladder with a PK sub-study, 29 patients), and trial 3 (phase IIb / III double-blind randomised dose-ranging placebo-controlled trial in neurogenic bladder, 115 patients). Step 1 base-model development used the frequent-sampling data from trials 1 and 2 (74 patients); step 2 covariate analysis added the sparse data from trial 3 (total 189 patients). Population baseline demographics summarised in Table 2 column 'Step 2 Trials 1, 2 and 3'."
  )

  ini({
    # =========================================================================
    # Structural disposition (Tsuda 2010 Table 3 'Parameter estimates from the
    # population pharmacokinetic final model'). The structural-model reference
    # is a 70 kg patient with AAG = 20 uM; the typical-value parameters
    # CL/F = 2.28 L/h and V/F = 37.5 L apply at that reference and the
    # allometric + AAG covariate factors below recreate individual-specific
    # values inside model(). Allometric exponents were estimated for the final
    # model (0.768 for CL/F and 0.807 for V/F) and the 95% CIs (0.598-0.938
    # and 0.603-1.01) included the theory values 0.75 and 1; the authors fixed
    # the exponents at theory values for the final model (Results 'Model
    # building' Step 2).
    # =========================================================================
    lcl   <- log(2.28)  ; label("Apparent oral clearance CL/F at WT 70 kg and AAG 20 uM (L/h)") # Tsuda 2010 Table 3 row 'CL/F (l h-1) = 2.28' (RSE 4.21%)
    lvc   <- log(37.5)  ; label("Apparent central volume of distribution V/F at WT 70 kg and AAG 20 uM (L)") # Tsuda 2010 Table 3 row 'V/F (l) = 37.5' (RSE 6.35%)
    lka   <- log(0.368) ; label("First-order absorption rate constant ka (1/h)")                # Tsuda 2010 Table 3 row 'Ka (h-1) = 0.368' (RSE 12.0%)
    ltlag <- log(0.957) ; label("Absorption lag time (h)")                                       # Tsuda 2010 Table 3 row 'ALAG1 (h) = 0.957' (RSE 1.14%)

    # =========================================================================
    # Covariate effects (Tsuda 2010 Table 3 'Structural model' equations).
    # Body-weight exponents are theory-fixed at the final-model step (see
    # comment above); AAG exponents are estimated power-form covariates.
    # =========================================================================
    e_wt_cl   <- fixed(0.75)   ; label("Allometric exponent of (WT/70) on CL/F (unitless; fixed at theory value)") # Tsuda 2010 Table 3 structural-model equation 'CL/F = theta_CL * (WT/70)^0.75 * ...'; Results 'Model building' Step 2: estimated 0.768 (95% CI 0.598-0.938), CI included 0.75 -> theory value used in final model
    e_wt_vc   <- fixed(1.0)    ; label("Allometric exponent of (WT/70) on V/F (unitless; fixed at theory value)")  # Tsuda 2010 Table 3 structural-model equation 'V/F = theta_V * (WT/70) * ...'; Results 'Model building' Step 2: estimated 0.807 (95% CI 0.603-1.01), CI included 1 -> theory value used in final model
    e_aag_cl  <- -0.844        ; label("Exponent of (AAG/20 uM) on CL/F (unitless; estimated)") # Tsuda 2010 Table 3 row 'theta_AAG_CL = -0.844' (RSE -15.1%)
    e_aag_vc  <- -0.663        ; label("Exponent of (AAG/20 uM) on V/F (unitless; estimated)")  # Tsuda 2010 Table 3 row 'theta_AAG_V = -0.663' (RSE -24.3%)

    # =========================================================================
    # Inter-individual variability (Tsuda 2010 Table 3 'IIV' rows). The paper
    # reports IIV as %CV with cohort covariance Cov_V/CL = 0.238 and the
    # interpreted correlation coefficient r = 0.715. The reported covariance
    # ties to the simpler approximation omega = %CV/100 rather than the exact
    # log-normal omega^2 = log(1 + CV^2): 0.715 * 0.544 * 0.612 = 0.238 (matches),
    # whereas 0.715 * sqrt(log(1+0.544^2)) * sqrt(log(1+0.612^2)) = 0.205 (no).
    # Therefore omega^2 values are taken as (CV/100)^2 directly:
    #   IIV CL/F  54.4% -> omega^2 = 0.544^2 = 0.2959
    #   IIV V/F   61.2% -> omega^2 = 0.612^2 = 0.3745
    #   IIV ka    117%  -> omega^2 = 1.17^2  = 1.3689
    # Correlated block on (CL/F, V/F): cov = 0.238 (Tsuda 2010 Table 3 row
    # 'Cov_V/CL = 0.238' with footnote 'translates to a coefficient of
    # correlation of 0.715'). IIV ka is independent (no Cov_Ka,CL or Cov_Ka,V
    # reported in Table 3).
    # =========================================================================
    etalcl + etalvc ~ c(0.2959,
                        0.238,    0.3745)  # Tsuda 2010 Table 3 IIV CL/F = 54.4%, IIV V/F = 61.2%, Cov_V/CL = 0.238 (correlation 0.715)
    etalka          ~ 1.3689                 # Tsuda 2010 Table 3 IIV ka = 117%

    # =========================================================================
    # Residual unexplained variability (Tsuda 2010 Table 3 'Residual random
    # effect model: Y = Yhat + Yhat * epsilon1 + epsilon2'). Combined
    # proportional + additive on the linear concentration scale.
    # =========================================================================
    propSd <- 0.284  ; label("Proportional residual error (fraction)")            # Tsuda 2010 Table 3 row 'Proportional residual variability (CV%) = 28.4' (RSE 6.13%)
    addSd  <- 0.178  ; label("Additive residual error (ng/mL)")                   # Tsuda 2010 Table 3 row 'Additive residual variability (SD ng ml-1) = 0.178' (RSE 27.0%)
  })

  model({
    # --- Individual PK parameters (Tsuda 2010 Table 3 structural-model
    # equations). Reference subject: WT 70 kg, AAG 20 uM.
    cl   <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl * (AAG / 20.0)^e_aag_cl
    vc   <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc * (AAG / 20.0)^e_aag_vc
    ka   <- exp(lka + etalka)
    tlag <- exp(ltlag)

    # --- One-compartment oral PK with first-order absorption + lag time.
    kel <- cl / vc

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    alag(depot) <- tlag

    # --- Observation: concentration in ng/mL. Internal units are dose in mg
    # and vc in L, so central/vc is in mg/L = ug/mL; multiply by 1000 to get
    # ng/mL to match the paper's reported concentration units.
    Cc <- central / vc * 1000
    Cc ~ add(addSd) + prop(propSd)
  })
}
