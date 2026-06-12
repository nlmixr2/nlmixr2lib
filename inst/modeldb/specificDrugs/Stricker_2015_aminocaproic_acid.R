Stricker_2015_aminocaproic_acid <- function() {
  description <- "Two-compartment IV population PK model for epsilon-aminocaproic acid (EACA) in infants undergoing craniofacial reconstruction and adolescents undergoing posterior spinal fusion surgery (Stricker 2015)"
  reference <- "Stricker PA, Gastonguay MR, Singh D, Fiadjoe JE, Sussman EM, Pruitt EY, Goebel TK, Zuppa AF. Population pharmacokinetics of epsilon-aminocaproic acid in adolescents undergoing posterior spinal fusion surgery. Br J Anaesth 2015;114(4):689-699. doi:10.1093/bja/aeu459"
  vignette <- "Stricker_2015_aminocaproic_acid"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Fixed allometric scaling with 70 kg reference; exponent 0.75 on CL and Q, 1.0 on V1 and V2 (Stricker 2015 Methods 'Full covariate model' and Table 5 caption).",
      source_name        = "WT"
    ),
    PNA = list(
      description        = "Postnatal age (chronological since birth)",
      units              = "months",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Drives the Emax-type CL maturation `PNA / (1.53 + PNA)` (Stricker 2015 Table 5 caption). Adolescent PSF cohort 158-211 months (~13-17.5 years); infant craniofacial cohort 6-25 months. Convert source ages from years/weeks to months before assignment. Effectively time-fixed within a single surgical observation window.",
      source_name        = "AGE"
    ),
    DIS_SCOL_IDIO = list(
      description        = "Idiopathic-scoliosis aetiology indicator (1 = idiopathic scoliosis or idiopathic kyphoscoliosis, 0 = otherwise)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not idiopathic-scoliosis; the reference complement here is the infant craniofacial-reconstruction cohort)",
      notes              = "Decomposed with DIS_SCOL_NONIDIO from a 3-level diagnosis / surgery-type categorical {craniofacial reference, idiopathic scoliosis, non-idiopathic scoliosis}; both indicators = 0 selects the craniofacial reference. Stricker 2015 Table 5 caption and Results p.694.",
      source_name        = "diagnosis"
    ),
    DIS_SCOL_NONIDIO = list(
      description        = "Non-idiopathic / syndromic scoliosis aetiology indicator (1 = non-idiopathic scoliosis: cerebral palsy, Marfan, Ehlers-Danlos, neurofibromatosis, spina bifida, congenital neuromuscular scoliosis, syringomyelia, 22q deletion, cortical dysgenesis; 0 = otherwise)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not non-idiopathic-scoliosis; the reference complement here is the infant craniofacial-reconstruction cohort)",
      notes              = "Decomposed with DIS_SCOL_IDIO from a 3-level diagnosis / surgery-type categorical {craniofacial reference, idiopathic scoliosis, non-idiopathic scoliosis}; both indicators = 0 selects the craniofacial reference. Stricker 2015 Table 5 caption and Table 1.",
      source_name        = "diagnosis"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 38,
    n_studies      = 2,
    age_range      = "6 months - 17.6 years (infant craniofacial cohort 27-107 weeks; adolescent PSF cohort 12-17 years)",
    weight_range   = "6.7-66.7 kg (infant cohort 6.7-11.8 kg; adolescent cohort 29.8-66.7 kg)",
    sex_female_pct = NA_real_,
    disease_state  = "Children undergoing major surgery with risk of significant blood loss: 18 infants undergoing craniofacial reconstruction for cranial synostosis, 10 adolescents with idiopathic scoliosis undergoing posterior spinal fusion, and 10 adolescents with non-idiopathic (syndromic) scoliosis undergoing posterior spinal fusion",
    dose_range     = "Craniofacial: 25 / 50 / 100 mg/kg IV loading bolus over 10 min followed by continuous IV infusion (CIVI) at 10 / 20 / 40 mg/kg/h. Adolescent PSF: 100 mg/kg IV loading bolus over 10 min followed by CIVI at 10 mg/kg/h until skin closure",
    regions        = "United States (Children's Hospital of Philadelphia)",
    notes          = "Pooled cohort: 20 adolescents from the PSF trial (NCT01408823, Stricker 2015) plus 18 infants from the earlier craniofacial reconstruction dose-escalation PK trial (Stricker 2013 BJA, Stricker 2015 reference 15). See Stricker 2015 Tables 1-3 for per-subject demographics. Race / ethnicity and sex distributions are not reported in the published tables."
  )

  ini({
    # Structural parameters at the 70 kg reference weight (Stricker 2015 Table 5 final-model column).
    # Units converted from the paper's mL/min/70kg to L/h/70kg by * 60 / 1000.
    lcl <- log(9.18);   label("Clearance for a 70 kg subject (CL, L/h)")                            # Stricker 2015 Table 5: 153 mL/min/70kg
    lvc <- log(8.78);   label("Central volume of distribution for a 70 kg subject (V1, L)")        # Stricker 2015 Table 5: 8.78 L/70kg
    lq  <- log(11.94);  label("Intercompartmental clearance for a 70 kg subject (Q, L/h)")         # Stricker 2015 Table 5: 199 mL/min/70kg
    lvp <- log(15.80);  label("Peripheral volume of distribution for a 70 kg subject (V2, L)")     # Stricker 2015 Table 5: 15.80 L/70kg

    # Allometric exponents, fixed at physiologic values per Stricker 2015 Methods 'Full covariate model'
    # ("an allometric power parameter ... fixed at 0.75 for clearances and at 1 for volumes").
    e_wt_cl_q  <- fixed(0.75); label("Shared allometric exponent on CL and Q (unitless)")  # Stricker 2015 Methods
    e_wt_vc_vp <- fixed(1.00); label("Shared allometric exponent on V1 and V2 (unitless)") # Stricker 2015 Methods

    # Emax-type postnatal-age maturation on CL: factor = PNA / (tm50_cl + PNA).
    # ~50% maturation at PNA = tm50_cl = 1.53 months and ~90% maturation by 15 months
    # (Stricker 2015 Results: "the model-estimated age at which 50% of full Cl was achieved was 1.53 months ...
    # 90% of full maturation Cl should occur at approximately 15 months").
    tm50_cl <- 1.53; label("Postnatal age at 50% CL maturation (months)")  # Stricker 2015 Table 5 'Age Cl 50%, months'

    # Diagnosis / surgery-type cohort effect on CL: multiplicative-power form with the craniofacial cohort
    # as the implicit reference (both indicators = 0 -> factor = 1.0). Authors note the effect is not
    # clinically relevant but retain it in the final model for magnitude / precision reporting
    # (Stricker 2015 Results p.694; Table 5 caption).
    e_dis_scol_idio_cl    <- 1.10; label("Multiplicative CL factor for idiopathic-scoliosis cohort (vs craniofacial reference)")    # Stricker 2015 Table 5 'Impact of idiopathic spines'
    e_dis_scol_nonidio_cl <- 0.97; label("Multiplicative CL factor for non-idiopathic-scoliosis cohort (vs craniofacial reference)") # Stricker 2015 Table 5 'Impact of non-idiopathic spines'

    # IIV: block-correlated on CL, V1, V2. No IIV on Q (Stricker 2015 Results:
    # "addition of an interindividual variability estimate on Q resulted in less precise parameter
    # estimates ... and therefore was not included in the final model").
    # Diagonal: omega^2 = (Table 5 reported value / 100)^2 since the caption reads
    # "Interindividual variability = (square root of variance)*100".
    #   CL : (0.2381)^2 = 0.0567
    #   V1 : (0.4990)^2 = 0.2490
    #   V2 : (0.2865)^2 = 0.0821
    # Off-diagonals are reported as variance-scale covariances directly in the Table 5 caption:
    #   cov(eta_CL, eta_V1) = 0.074
    #   cov(eta_CL, eta_V2) = 0.052
    #   cov(eta_V1, eta_V2) = 0.064
    etalcl + etalvc + etalvp ~ c(0.0567,
                                 0.074, 0.2490,
                                 0.052, 0.064, 0.0821)

    # Residual error (combined additive + proportional in linear-concentration space per Methods Eq.
    # "C_obs = C_pred * (1 + eps_P) + eps_A"). propSd / addSd are sqrt(reported variance) so that
    # the in-file values map directly to the published sigma^2 values.
    propSd <- sqrt(0.026); label("Proportional residual error (fraction)")  # Stricker 2015 Table 5: sigma^2 proportional = 0.026
    addSd  <- sqrt(0.673); label("Additive residual error (mg/L)")          # Stricker 2015 Table 5: sigma^2 additive = 0.673
  })

  model({
    # 1. Derived covariate terms
    # Emax-type postnatal-age maturation on CL: 0 at birth, 50% at PNA = 1.53 months,
    # ~90% by 15 months, ~1 in adolescents.
    age_factor_cl <- PNA / (tm50_cl + PNA)

    # Diagnosis / surgery-type cohort effect on CL (multiplicative power form;
    # both indicators 0 selects the craniofacial reference -> factor = 1.0).
    diag_cl <- (e_dis_scol_idio_cl ^ DIS_SCOL_IDIO) * (e_dis_scol_nonidio_cl ^ DIS_SCOL_NONIDIO)

    # 2. Individual PK parameters with fixed allometric weight scaling (reference 70 kg)
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl_q  * age_factor_cl * diag_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc_vp
    vp <- exp(lvp + etalvp) * (WT / 70)^e_wt_vc_vp
    q  <- exp(lq)           * (WT / 70)^e_wt_cl_q

    # 3. Micro-constants
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # 4. Two-compartment IV ODE system (no depot; dosing into central)
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # 5. Observation: plasma EACA concentration (dose mg, volume L -> mg/L = mcg/mL)
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
