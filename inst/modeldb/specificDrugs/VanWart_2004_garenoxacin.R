VanWart_2004_garenoxacin <- function() {
  description <- "One-compartment population pharmacokinetic model with first-order absorption and first-order elimination for oral garenoxacin (a des-F(6) quinolone) in adults with community-acquired respiratory tract infections (Van Wart 2004); CL/F covariates are creatinine clearance, ideal body weight, age, obesity (WT > 130% IBW), and concomitant pseudoephedrine; V/F covariates are body weight and male sex."
  reference <- "Van Wart S, Phillips L, Ludwig EA, Russo R, Gajjar DA, Bello A, Ambrose PG, Costanzo C, Grasela TH, Echols R, Grasela DM. Population pharmacokinetics and pharmacodynamics of garenoxacin in patients with community-acquired respiratory tract infections. Antimicrob Agents Chemother. 2004;48(12):4766-4777. doi:10.1128/AAC.48.12.4766-4777.2004"
  vignette <- "VanWart_2004_garenoxacin"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    CRCL = list(
      description        = "Cockcroft-Gault creatinine clearance (raw, not BSA-normalized)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Raw Cockcroft-Gault creatinine clearance in mL/min (NOT BSA-normalized to mL/min/1.73 m^2); stored under canonical CRCL per the precedent in Delattre_2010_amikacin.R. Used with power-form normalization (CRCL / 86.9)^0.436 on CL. Reference 86.9 mL/min is the development-cohort median (Van Wart 2004 Table 2). Cohort range 14.5-205 mL/min.",
      source_name        = "CRCL"
    ),
    IBW = list(
      description        = "Ideal body weight (Devine formula)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Ideal body weight by the Devine formula (Van Wart 2004 reference 30). Used in an additive linear deviation term 0.764 * (IBW - 64.2) on CL/F. Reference 64.2 kg is the development-cohort median (Table 2). Cohort range 20.9-96.9 kg.",
      source_name        = "IBW"
    ),
    AGE = list(
      description        = "Subject age at baseline",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Used in an additive linear deviation term 0.301 * (AGE - 49.5) on CL/F. Reference 49.5 years is the development-cohort median (Van Wart 2004 Table 2). Cohort range 18-88 years (outpatients aged 18 and over).",
      source_name        = "AGE"
    ),
    WT = list(
      description        = "Body weight at baseline",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Used in two places: (a) power-form allometric scaling on V/F: (WT / 79.3)^0.635 with reference 79.3 kg (development-cohort median, Van Wart 2004 Table 2); (b) derivation of the OBESE flag inside model() as WT > 1.3 * IBW (Van Wart 2004 Methods: 'obesity (defined as a WTKG greater than 130% of the IBW)'). Cohort range 34-178 kg.",
      source_name        = "WTKG"
    ),
    SEXF = list(
      description        = "Biological sex indicator (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male) -- the source paper uses male as the reference category for the additive V/F shift (V/F = 67.1 L for the typical 79.3 kg female).",
      notes              = "Van Wart 2004 encodes sex as SEXM (1 = male, 0 = female) with female as the V/F intercept; the additive shift +17.7 L applies to males. To store under the canonical SEXF (1 = female, 0 = male) while preserving the paper's female-reference V/F = 67.1 L at WT = 79.3 kg, the male shift is applied as 17.7 * (1 - SEXF) inside model(). Cohort 50% female (Table 2 development set).",
      source_name        = "SEXM"
    ),
    CONMED_PSEUDOEPHEDRINE = list(
      description        = "Concomitant pseudoephedrine coadministration (sample-level)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant pseudoephedrine at the sample)",
      notes              = "Van Wart 2004 PSEU_jk: 1 if the kth PK sample from the jth patient was collected during concomitant pseudoephedrine administration, 0 otherwise. Sample-level / time-varying within subject. Cohort prevalence 14% of patients (82/580; Table 1). Multiplicative effect on the structural CL/F covariate term: cl_structural *= (1 - 0.144 * CONMED_PSEUDOEPHEDRINE). Mechanism (Van Wart 2004 Discussion): pseudoephedrine is renally excreted by active tubular secretion in addition to glomerular filtration (renal CL 511-532 mL/min) and competes with garenoxacin's tubular secretion.",
      source_name        = "PSEU"
    )
  )

  population <- list(
    n_subjects     = 580L,
    n_studies      = 3L,
    n_observations = 1529L,
    age_range      = "18-88 years",
    age_mean       = "49.5 years",
    weight_range   = "34-178 kg",
    weight_mean    = "79.3 kg",
    ibw_range      = "20.9-96.9 kg",
    ibw_mean       = "64.2 kg",
    crcl_range     = "14.5-205 mL/min (Cockcroft-Gault, raw)",
    crcl_mean      = "86.9 mL/min",
    sex_female_pct = 50.3,
    race_ethnicity = "Caucasian 85%, Hispanic 10%, Black 6% (development set; Van Wart 2004 Table 2). Race was tested as a V/F covariate during forward selection but eliminated during backward stepwise (P > 0.00684).",
    disease_state  = "Adult outpatients (18+) with community-acquired respiratory tract infections: acute exacerbation of chronic bronchitis (ABECB), community-acquired pneumonia (CAP), or acute bacterial sinusitis (ABS).",
    dose_range     = "Oral garenoxacin 400 mg once daily for 5-10 days (5-day arm received placebo days 6-10 in the ABECB study; CAP and ABS were 10-day open-label).",
    regions        = "Three multinational phase II clinical trials sponsored by Bristol-Myers Squibb (T-3811ME / BMS-284756 program).",
    obesity_pct    = "34% obese (WTKG > 130% IBW) in development set",
    notes          = "Population PK analysis pooled across three phase II trials. PK sampling per protocol: predose (0 h), 2 h after first dose, and the day 3-5 visit. Plasma garenoxacin quantified by LC/MS/MS (LOQ 0.01 ug/mL, range 0.01-10 ug/mL). 80/20 random split into development (580 patients, 1529 samples) and validation (141 patients, 379 samples) data sets. Patients with serum creatinine > 2.0 mg/dL or hepatic enzymes >= 3x ULN were excluded; subjects with severe renal dysfunction were rare (n = 3 with CrCL 15-29 mL/min) and the model is not validated below CrCL 15 mL/min. Food effect (fed vs fasted) was tested and not significant (P = 0.8516; F_fed = 98.2%), so bioavailability is set to 1 across dosing states. Final population PK parameter estimates from Van Wart 2004 Table 4; covariate equations from Results 'Final model' section."
  )

  ini({
    # Structural parameters. Population means correspond to a reference
    # subject with CRCL = 86.9 mL/min, IBW = 64.2 kg, AGE = 49.5 yr,
    # WTKG = 79.3 kg, non-obese (WT <= 130% IBW), female (SEXF = 1), and
    # no concomitant pseudoephedrine -- i.e. the cohort-median female
    # (Van Wart 2004 Table 2 development set medians).

    # Absorption rate constant. No IIV (Van Wart 2004 Results: "the
    # interindividual variability of ka was removed from the model"
    # because the protocol sampling provided limited absorption-phase
    # information).
    lka <- log(2.41);  label("First-order absorption rate constant (Ka, 1/h)")               # Van Wart 2004 Table 4 final-model Ka

    # Apparent oral clearance intercept reported in mL/min (Van Wart 2004
    # Table 4). The model() block converts to L/h to match the time unit.
    lcl <- log(83.4);  label("Apparent oral clearance structural coefficient at reference covariates (CL/F, mL/min)") # Van Wart 2004 Table 4 final-model CL/F coefficient

    # Apparent volume of distribution structural coefficient for the
    # reference female; the male additive shift is applied separately.
    lvc <- log(67.1);  label("Apparent volume of distribution structural coefficient for a 79.3 kg female (V/F, L)") # Van Wart 2004 Table 4 final-model V/F coefficient

    # Covariate effects on CL/F (Van Wart 2004 Table 4 + Results 'Final
    # model' equation).
    e_crcl_cl  <-  0.436;  label("Power exponent of (CRCL / 86.9) on CL/F (unitless)")                # Van Wart 2004 Table 4: CRCL-CL power
    e_ibw_cl   <-  0.764;  label("Additive linear-deviation slope of (IBW - 64.2) on CL/F (mL/min per kg)") # Van Wart 2004 Table 4: IBW-CL slope
    e_obese_cl <- 10.9;    label("Additive shift in CL/F for obese subjects (WT > 130% IBW) (mL/min)")     # Van Wart 2004 Table 4: additive shift in CL for obesity
    e_age_cl   <-  0.301;  label("Additive linear-deviation slope of (AGE - 49.5) on CL/F (mL/min per year)") # Van Wart 2004 Table 4: Age-CL slope
    e_pseu_cl  <- -0.144;  label("Proportional shift in CL/F for concomitant pseudoephedrine (unitless)") # Van Wart 2004 Table 4: pseudoephedrine proportional shift

    # Covariate effects on V/F (Van Wart 2004 Table 4 + Results equation).
    e_wt_vc    <-  0.635;  label("Power exponent of (WT / 79.3) on V/F (unitless)")                    # Van Wart 2004 Table 4: WTKG-V power
    e_sexm_vc  <- 17.7;    label("Additive shift in V/F for males vs the female reference (L)")         # Van Wart 2004 Table 4: additive shift in V for males

    # Inter-individual variability (log-normal exponential model).
    # Reported as %CV in Van Wart 2004 Table 4; omega^2 = log(CV^2 + 1).
    # CL/F: 25.5% CV -> log(0.255^2 + 1) = 0.06300
    # V/F:  18.8% CV -> log(0.188^2 + 1) = 0.03473
    # Diagonal -- Van Wart 2004 Results: "the interindividual variability
    # terms for clearance and volume were not correlated".
    etalcl ~ 0.06300   # 25.5% CV (Van Wart 2004 Table 4)
    etalvc ~ 0.03473   # 18.8% CV (Van Wart 2004 Table 4)

    # Proportional residual error (Van Wart 2004 Table 4: 27.6% CV).
    propSd <- 0.276;   label("Proportional residual error on Cc (fraction)")                          # Van Wart 2004 Table 4 final-model residual variability
  })

  model({
    # ----- Derived covariate flags -----
    # OBESE = 1 if body weight exceeds 130% of ideal body weight (Van Wart
    # 2004 Methods: 'obesity (defined as a WTKG greater than 130% of the
    # IBW)'). Derived per subject from canonical WT and IBW; no separate
    # OBESE covariate column is required.
    obese_flag <- ifelse(WT > 1.3 * IBW, 1, 0)

    # SEXM = 1 if male, 0 if female (Van Wart 2004 paper notation). The
    # canonical SEXF stores 1 = female, so the male indicator is (1 - SEXF).
    sex_male <- 1 - SEXF

    # ----- Individual PK parameters -----
    # ka has no IIV.
    ka <- exp(lka)

    # CL/F covariate model. The eta is exponential on the whole structural
    # CL/F expression (Van Wart 2004 Methods: 'Interindividual variability
    # of all PK parameters was described with an exponential error model').
    # Structural CL/F in mL/min (paper units):
    #   CL_struct (mL/min)
    #     = [ 83.4 * (CRCL/86.9)^0.436
    #         + 0.764 * (IBW - 64.2)
    #         + 10.9  * OBESE
    #         + 0.301 * (AGE - 49.5)
    #       ] * (1 - 0.144 * PSEU)
    cl_mlmin <- (
      exp(lcl) * (CRCL / 86.9)^e_crcl_cl
      + e_ibw_cl   * (IBW - 64.2)
      + e_obese_cl * obese_flag
      + e_age_cl   * (AGE - 49.5)
    ) * (1 + e_pseu_cl * CONMED_PSEUDOEPHEDRINE) * exp(etalcl)

    # Convert apparent oral clearance from mL/min (paper unit) to L/h to
    # match the model's time = hour, volume = L convention.
    cl <- cl_mlmin * 60 / 1000

    # V/F covariate model. Structural V/F in L:
    #   V_struct (L) = 67.1 * (WTKG/79.3)^0.635 + 17.7 * SEXM
    vc <- (
      exp(lvc) * (WT / 79.3)^e_wt_vc
      + e_sexm_vc * sex_male
    ) * exp(etalvc)

    kel <- cl / vc

    # ----- ODEs: one-compartment open model with first-order absorption -----
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Observation. Dose in mg, vc in L -> central/vc has units mg/L
    # (equivalently ug/mL, the paper's assay unit).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
