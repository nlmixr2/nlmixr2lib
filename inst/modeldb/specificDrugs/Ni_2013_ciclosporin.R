Ni_2013_ciclosporin <- function() {
  description <- "One-compartment first-order-absorption population PK model for oral ciclosporin in Chinese children with aplastic anemia (Ni 2013)"
  reference <- "Ni SQ, Zhao W, Wang J, Zeng S, Chen SQ, Jacqz-Aigrain E, Zhao ZY. Population pharmacokinetics of ciclosporin in Chinese children with aplastic anemia: effects of weight, renal function and stanozolol administration. Acta Pharmacol Sin. 2013;34(7):969-975. doi:10.1038/aps.2013.9"
  vignette <- "Ni_2013_ciclosporin"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying allowed; used for allometric scaling with reference weight 70 kg (Ni 2013 Methods 'Covariate analysis').",
      source_name        = "WT"
    ),
    CREAT = list(
      description        = "Serum creatinine concentration",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Used in the renal-function factor RF = 1 - (CREAT/44) * e_creat_cl on CL/F. The denominator 44 umol/L is a cohort-derived scaling constant close to the cohort median CREAT of 44.2 umol/L (Ni 2013 Table 1).",
      source_name        = "CREA"
    ),
    CONMED_STANOZOLOL = list(
      description        = "Concomitant stanozolol coadministration indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no stanozolol coadministration)",
      notes              = "78 of 102 children were on stanozolol (Ni 2013 Table 1). Time-fixed per subject in the published analysis. Power-form multiplicative effect on CL/F: cl *= e_conmed_stanozolol_cl ^ CONMED_STANOZOLOL (Ni 2013 Table 2 'F comedication').",
      source_name        = "stanozolol indicator (no specific NONMEM column reported)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 102L,
    n_studies      = 1L,
    age_range      = "0.9-17.6 years",
    age_median     = "8.7 years (mean +/- SD 8.8 +/- 3.6)",
    weight_range   = "6.5-69.0 kg",
    weight_median  = "29.8 kg (mean +/- SD 31.3 +/- 12.7)",
    sex_female_pct = 52.9,
    race_ethnicity = "Chinese",
    disease_state  = "Chinese children with acquired or congenital aplastic anemia",
    dose_range     = "30-375 mg/day oral (median 150 mg/day; 5.3 +/- 1.7 mg/kg/day; administered BID, target trough C0 = 200 ng/mL)",
    regions        = "China (Children's Hospital of Zhejiang University, Hangzhou)",
    notes          = "Therapeutic drug monitoring (TDM) dataset, 2003-2009; 592 trough blood ciclosporin concentrations measured by FPIA (median sampling time 12.1 h post-dose, range 10.0-16.0 h). 78 of 102 patients were on stanozolol, 79 on prednisone; only stanozolol retained in the final covariate model. Baseline serum creatinine median 44.2 umol/L (range 5.8-120.0); see Ni 2013 Table 1."
  )

  ini({
    # Structural parameters; reference body weight 70 kg
    lka <- fixed(log(0.68))
    label("Absorption rate constant (Ka, 1/h); fixed from prior literature (Ni 2013 Results paragraph 1; Table 2)")  # Table 2 'Ka (Fixed)'
    lvc <- log(178)
    label("Apparent volume of distribution at WT = 70 kg (V/F, L); Ni 2013 Table 2 theta1")  # Table 2 theta1
    lcl <- log(31.5)
    label("Apparent oral clearance at WT = 70 kg, CREA = 0, no stanozolol (CL/F, L/h); Ni 2013 Table 2 theta2")  # Table 2 theta2

    # Allometric exponents fixed from West / Holford convention (Ni 2013 Methods; refs 11, 12)
    e_wt_cl <- fixed(0.75)
    label("Allometric (WT) exponent on CL/F (unitless); fixed (Ni 2013 Methods)")  # Methods; Table 2 equation
    e_wt_vc <- fixed(1.0)
    label("Allometric (WT) exponent on V/F (unitless); fixed (Ni 2013 Methods)")  # Methods; Table 2 equation

    # Covariate effects on CL/F
    e_creat_cl <- 0.0821
    label("Linear slope of (CREAT/44) inside the renal-function factor on CL/F (unitless); Ni 2013 Table 2 theta3")  # Table 2 theta3
    e_conmed_stanozolol_cl <- 0.83
    label("Multiplicative CL/F factor when stanozolol is coadministered (F_comedication; ratio); Ni 2013 Table 2")  # Table 2 'F comedication'

    # Inter-individual variability (omega^2 = log(CV^2 + 1))
    etalka ~ 0.29423
    # 58.5% CV on Ka (Ni 2013 Table 2); log(0.585^2 + 1) = 0.29423
    etalcl ~ 0.01625
    # 12.8% CV on CL/F (Ni 2013 Table 2); log(0.128^2 + 1) = 0.01625

    # Residual error (combined proportional + additive)
    propSd <- 0.224
    label("Proportional residual error (fraction); 22.4% per Ni 2013 Table 2")  # Table 2
    addSd <- 34.1
    label("Additive residual error (ng/mL); Ni 2013 Table 2")  # Table 2
  })

  model({
    # Renal-function multiplier on CL/F
    cl_rf <- 1 - (CREAT / 44) * e_creat_cl

    # Individual PK parameters
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl * cl_rf * e_conmed_stanozolol_cl^CONMED_STANOZOLOL
    vc <- exp(lvc) * (WT / 70)^e_wt_vc

    kel <- cl / vc

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Dose in mg, vc in L, so central/vc is in mg/L = ug/mL; *1000 -> ng/mL (paper units)
    Cc <- central / vc * 1000
    Cc ~ add(addSd) + prop(propSd)
  })
}
