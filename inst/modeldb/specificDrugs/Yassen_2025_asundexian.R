Yassen_2025_asundexian <- function() {
  description <- "Two-compartment population PK model with two transit absorption compartments for asundexian, an oral selective Factor XIa inhibitor, in healthy volunteers and adult patients at risk for thromboembolic / cardiovascular events (Yassen 2025)"
  reference <- paste(
    "Yassen A, Kanefendt F, Zisowsky J, Broeker A, Mundl H, Vis P,",
    "Garmann D, Berkhout J. Population Pharmacokinetics of Asundexian",
    "in People at Risk for Thromboembolic/Cardiovascular Events.",
    "CPT Pharmacometrics Syst Pharmacol. 2026;15:e70142.",
    "doi:10.1002/psp4.70142",
    sep = " "
  )
  vignette <- "Yassen_2025_asundexian"
  units <- list(time = "h", dosing = "mg", concentration = "ug/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on CL/F and Vc/F with reference 79 kg (typical participant in PACIFIC-STROKE).",
      source_name        = "BW"
    ),
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on CL/F and Vc/F with reference 68 years (population median).",
      source_name        = "AGE"
    ),
    SEXF = list(
      description        = "Biological sex indicator (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Proportional-shift effect on CL/F (-12.1%) and Vc/F (-14.0%) for females vs the male reference.",
      source_name        = "SEX"
    ),
    CRCL = list(
      description        = "CKD-EPI estimated glomerular filtration rate (eGFR)",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on CL/F with reference 77 mL/min/1.73 m^2 (typical participant in PACIFIC-STROKE). Computed from serum creatinine via the CKD-EPI formula per the source paper.",
      source_name        = "EGFR"
    ),
    CONMED_CYP3A4_INH = list(
      description        = "Concomitant CYP3A4 inhibitor coadministration indicator (1 = on inhibitor, 0 = none)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no CYP3A4 inhibitor coadministration)",
      notes              = "Pools weak and moderate CYP3A4 inhibitors into the CONMED_CYP3A4_INH = 1 category. Strong CYP3A4 inhibitors were a Phase II exclusion criterion; only 0.6% of participants received a strong inhibitor in Phase I, and the source paper combined all three strength tiers into a single covariate after testing weak / moderate alone showed no separable signal. Proportional-shift effect on CL/F (-5.3%) when CONMED_CYP3A4_INH = 1. Renamed from canonical CYP3A4_INH to CONMED_CYP3A4_INH on 2026-06-19 per the canonical-register standardization audit.",
      source_name        = "CYP3AI"
    )
  )

  population <- list(
    n_subjects     = 2914,
    n_studies      = 9,
    age_range      = "21-92 years",
    age_median     = "68 years",
    weight_range   = "34.6-178 kg",
    weight_median  = "80 kg",
    sex_female_pct = 30.8,
    race_ethnicity = c(White = 83.3, Asian = 15.0, Black = 0.769, Other = 0.879),
    disease_state  = "Pooled cohort: healthy volunteers (six Phase I studies including age/sex, multiple-dose, renal-impairment, hepatic-impairment, Japanese, and Chinese sub-studies) and adult participants at risk for thromboembolic / cardiovascular events from three Phase II studies (PACIFIC-STROKE, acute non-cardioembolic ischaemic stroke; PACIFIC-AMI, recent acute myocardial infarction; PACIFIC-AF, atrial fibrillation).",
    dose_range     = "10, 20, 25, 50, and 100 mg oral once daily across studies (Phase I 25/50/100 mg OD or BID; Phase II 10/20/50 mg OD).",
    regions        = "Multinational; Phase I included Japanese and Chinese sub-studies",
    notes          = "Baseline demographics from Tables 2 and 3 (n = 2914 across all nine pooled studies; 16,599 PK observations). eGFR is computed by CKD-EPI; CYP3A4 inhibitor coadministration pools weak + moderate strengths."
  )

  ini({
    # Reference participant (PACIFIC-STROKE typical): BW = 79 kg, age = 68 y,
    # eGFR = 77 mL/min/1.73 m^2, male, no CYP3A4 inhibitor.
    # Source: Yassen 2025 Table 4 (Estimate column) and section 3.2.

    lcl <- log(2.25);  label("Apparent clearance, reference participant (CL/F, L/h)")            # Table 4
    lvc <- log(35.3);  label("Apparent central volume of distribution, reference participant (Vc/F, L)")  # Table 4 -- model assumes Vp/F = Vc/F (overparameterisation constraint)
    lka <- log(1.87);  label("First-order absorption rate constant (Ka, 1/h)")                    # Table 4
    lq  <- log(27.4);  label("Apparent inter-compartmental clearance (Q/F, L/h)")                  # Table 4

    # Covariate effects on CL/F
    e_age_cl        <- -0.426;   label("Power exponent of age on CL/F (unitless)")                 # Table 4
    e_wt_cl         <-  0.396;   label("Power exponent of body weight on CL/F (unitless)")         # Table 4
    e_sexf_cl       <- -0.121;   label("Proportional shift in CL/F for females (unitless)")        # Table 4
    e_crcl_cl       <-  0.184;   label("Power exponent of eGFR (CKD-EPI) on CL/F (unitless)")      # Table 4
    e_cyp3a4_inh_cl <- -0.0531;  label("Proportional shift in CL/F with CYP3A4 inhibitor coadministration (unitless)")  # Table 4

    # Covariate effects on Vc/F
    e_age_vc  <- -0.176;  label("Power exponent of age on Vc/F (unitless)")               # Table 4
    e_wt_vc   <-  0.754;  label("Power exponent of body weight on Vc/F (unitless)")       # Table 4
    e_sexf_vc <- -0.14;   label("Proportional shift in Vc/F for females (unitless)")      # Table 4

    # IIV: 3x3 omega block (lower triangle order = var(CL), cov(CL,Vc),
    # var(Vc), cov(CL,Ka), cov(Vc,Ka), var(Ka)).
    # Variances on log scale derived from Table 4 CV%: omega^2 = log(1 + CV^2).
    #   CL/F CV = 30.1% -> 0.0867; Vc/F CV = 17.2% -> 0.0292; Ka CV = 49.2% -> 0.2168.
    # Covariances are reported directly on the log scale in Table 4.
    etalcl + etalvc + etalka ~ c(
      0.0867,
       0.0150,  0.0292,
      -0.0385, -0.00496, 0.2168
    )                                                                                 # Table 4

    # Residual error: NONMEM "additive on log scale" with sigma^2 = 0.0619 maps
    # to proportional error in nlmixr2 linear space; propSd = sqrt(0.0619) ~ 0.2488.
    propSd <- 0.2488; label("Proportional residual error (fraction)")               # Table 4
  })
  model({
    # Reference values for the typical participant (Table 4 footnote / section 3.2)
    ref_age  <- 68
    ref_wt   <- 79
    ref_crcl <- 77

    # Individual parameters
    cl <- exp(lcl + etalcl) *
          (AGE  / ref_age )^e_age_cl  *
          (WT   / ref_wt  )^e_wt_cl   *
          (1 + e_sexf_cl       * SEXF) *
          (CRCL / ref_crcl)^e_crcl_cl *
          (1 + e_cyp3a4_inh_cl * CONMED_CYP3A4_INH)
    vc <- exp(lvc + etalvc) *
          (AGE / ref_age)^e_age_vc *
          (WT  / ref_wt )^e_wt_vc  *
          (1 + e_sexf_vc * SEXF)
    ka <- exp(lka + etalka)
    q  <- exp(lq)
    vp <- vc  # paper constrains Vp/F = Vc/F to address model instability from imbalanced sampling

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Two-transit-compartment first-order absorption + 2-cmt linear disposition
    d/dt(depot)       <- -ka * depot
    d/dt(transit1)    <-  ka * depot     - ka  * transit1
    d/dt(transit2)    <-  ka * transit1  - ka  * transit2
    d/dt(central)     <-  ka * transit2  - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central  - k21 * peripheral1

    # Concentration: dose in mg, volumes in L -> central / vc gives mg/L;
    # multiply by 1000 to express asundexian concentration in ug/L per the source paper.
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
