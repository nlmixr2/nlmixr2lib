Huang_2021_cholecalciferol_pbpk <- function() {
  description <- paste(
    "PBPK (minimal whole-body, deSolve). Plasma vitamin D3 and",
    "25-hydroxyvitamin D3 in healthy adults given single or repeated daily",
    "oral vitamin D3, fitted to arm-mean data pooled from published trials",
    "(Huang 2021, final model Run009b). Vitamin D3 is absorbed first-order",
    "from a gastrointestinal depot into the liver and distributes between",
    "venous blood, arterial blood, liver and a lumped rest-of-body",
    "compartment; hepatic vitamin D3 clearance takes one value after a",
    "single dose and another under repeated daily dosing, and one third of",
    "it forms 25(OH)D3. 25(OH)D3 occupies its own four compartments and",
    "leaves the liver through a sigmoidal (Hill) clearance that rises with",
    "the 25(OH)D3 concentration. A constant endogenous vitamin D3 input",
    "ENDOG is back-calculated from the baseline plasma 25(OH)D3 so that the",
    "untreated system starts, and stays, at that baseline. Typical values",
    "only: the paper fitted arm means with a single parameter set, with no",
    "interindividual variability and no reported residual error."
  )
  reference <- paste(
    "Huang Z, You T. Personalise vitamin D3 using physiologically based",
    "pharmacokinetic modelling. CPT Pharmacometrics Syst Pharmacol.",
    "2021;10(7):723-734. doi:10.1002/psp4.12640.",
    sep = " "
  )
  vignette <- "Huang_2021_cholecalciferol_pbpk"

  # The deposited final-model script (supplement PSP4-10-723-s001.zip,
  # Final_model.R) works in hours, carries both analytes as nmol and converts
  # an oral dose in ug to nmol with DOSE/384.64*1000 before adding it to the
  # depot. Doses to this model are therefore in nmol (see the vignette).
  units <- list(
    time = "h",
    dosing = "nmol",
    concentration = "nmol/L"
  )

  compartmentData <- list(
    depot = list(
      analyte = "cholecalciferol",
      units = "nmol",
      specimen = "administration site",
      verified = TRUE
    ),
    venous = list(
      analyte = "cholecalciferol",
      units = "nmol",
      specimen = "plasma",
      verified = TRUE
    ),
    liver = list(
      analyte = "cholecalciferol",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    other = list(
      analyte = "cholecalciferol",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    arterial = list(
      analyte = "cholecalciferol",
      units = "nmol",
      specimen = "whole blood",
      verified = TRUE
    ),
    venous_25d3 = list(
      analyte = "25-hydroxyvitamin D3",
      units = "nmol",
      specimen = "plasma",
      verified = TRUE
    ),
    liver_25d3 = list(
      analyte = "25-hydroxyvitamin D3",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    other_25d3 = list(
      analyte = "25-hydroxyvitamin D3",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    arterial_25d3 = list(
      analyte = "25-hydroxyvitamin D3",
      units = "nmol",
      specimen = "whole blood",
      verified = TRUE
    )
  )

  covariateData <- list(
    D25OH_BL = list(
      description = "Baseline (pre-supplementation) plasma 25-hydroxyvitamin D concentration",
      units = "nmol/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-fixed. Huang 2021 Methods 'Model development': 'Steady state",
        "was assumed prior to dosing. Hence, ENDOG equals to the clearance",
        "rate of basal plasma 25(OH)D (D25BASE)', and every initial condition",
        "is expressed through D25BASE (supplement, Final Model section). Each",
        "published arm was simulated from that arm's mean baseline; the",
        "deposited script reads it as the Time = 0 concentration of the arm.",
        "Source symbol D25BASE."
      ),
      source_name = "D25BASE"
    ),
    REGI_QD = list(
      description = "Repeated once-daily vitamin D3 dosing (1) versus a single dose (0)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "Per-paper comparator: a SINGLE oral dose. Huang 2021 Methods: 'We",
        "assumed different values for vitamin D clearance under single dose",
        "(SCLH) and repeated daily doses (MCLH). This is similar to",
        "interoccasion variability'; supplement Final Model section: 'For",
        "single dose: CLH = SCLH. For repeated dose: CLH = MCLH.' The flag",
        "selects which hepatic vitamin D3 clearance applies for the whole",
        "simulation, including the pre-dose steady state that sets the",
        "vitamin D3 initial conditions. Only these two regimens were modelled",
        "(Methods: 'We excluded trials with other dosing frequencies')."
      ),
      source_name = "SCLH / MCLH"
    )
  )

  population <- list(
    species = "human",
    n_subjects = "307 (vitamin D3 PK) and 6484 (25(OH)D3 PK), arm means only",
    n_studies = "155 treatment arms from published trials, January 1970 to January 2019",
    age_range = "18 years and over (most vitamin D3 PK subjects 30-50 years)",
    sex_female_pct = NA,
    disease_state = "adults without disease or conditions known to alter vitamin D PK; normal renal function",
    dose_range = paste(
      "oral vitamin D3 single doses 70-50000 ug; repeated daily doses",
      "10-1250 ug/day; 25(OH)D3 7 and 20 ug/day"
    ),
    regions = "all continents except Antarctica; mostly USA, the Netherlands, UK and Canada",
    notes = paste(
      "Huang 2021 Results 'PK data collection': vitamin D3 PK from 307",
      "subjects in 13 arms of 6 trials (Table S2); 25(OH)D3 PK as 451 mean",
      "concentrations from 6484 subjects in 126 arms (Tables S3, S4).",
      "Training set: all vitamin D3 PK plus the 25(OH)D3 arms dosed 10 and",
      "100 ug/day (43 arms). Test set: 83 repeated-dose arms (12.5-1250",
      "ug/day), 16 high single-dose arms (1250-50000 ug) and two 25(OH)D3",
      "dosing arms. By mean BMI, 31 arms normal, 59 overweight, 18 obese, 47",
      "unknown. Plasma and serum treated as equivalent."
    )
  )

  ini({
    # ---------------------------------------------------------------
    # VITAMIN D3 (model Run004). Huang 2021 Table 1, posterior mean
    # (SD of the marginal MCMC posterior in parentheses). Supplement
    # Table S7 Run009b: 'Vitamin D PK: model Run004 with expected
    # parameter values'.
    # ---------------------------------------------------------------
    lka <- log(0.19)
    label("Vitamin D3 first-order absorption rate constant (log h^-1)") # Huang 2021 Table 1: Ka posterior 0.19 (0.05) h^-1; prior 0.20
    lkp_liver <- fixed(log(1))
    label("Vitamin D3 liver:plasma partition coefficient Kpl (unitless)") # Huang 2021 Table 1: Kpl = 1 (fixed)
    lkp_other <- log(0.09)
    label("Vitamin D3 rest-of-body:plasma partition coefficient Kprb (unitless)") # Huang 2021 Table 1: Kprb posterior 0.09 (0.02)
    lsclh <- log(0.32)
    label("Vitamin D3 hepatic clearance after a single dose SCLH (log L/h)") # Huang 2021 Table 1: SCLH posterior 0.32 (0.06); table unit 'h-1', but the ODE multiplies it by a concentration (see vignette)
    lmclh <- log(0.21)
    label("Vitamin D3 hepatic clearance under repeated daily dosing MCLH (log L/h)") # Huang 2021 Table 1: MCLH posterior 0.21 (0.004); table unit 'h-1', see vignette

    # ---------------------------------------------------------------
    # 25(OH)D3 (model Run009b). Huang 2021 Table 1 posterior means; they
    # equal the arithmetic means of exp() of the deposited MCMC chain
    # MCMC-Run009b.RDS (0.542, 86.33, 5.636, 0.03277).
    # ---------------------------------------------------------------
    fm_25d3 <- fixed(1 / 3)
    label("Fraction of hepatic vitamin D3 clearance forming 25(OH)D3 (unitless)") # Huang 2021 Table 1: Fm = 0.33 (fixed); the supplement ODEs and script use exactly 1/3 (formation x 1/3, ENDOG x 3)
    lkp_liver_25d3 <- fixed(log(1))
    label("25(OH)D3 liver:plasma partition coefficient Kp25l (unitless)") # Huang 2021 supplement Final Model ODEs carry Kp25l; value 1 from the deposited Final_model.R (Kp25l=1); not estimated (Table S7 Run009b fits only Kp25rb, C50, gamma, CLmax)
    lkp_other_25d3 <- log(0.54)
    label("25(OH)D3 rest-of-body:plasma partition coefficient Kp25rb (unitless)") # Huang 2021 Table 1: Kp25rb posterior 0.54 (0.15); prior 0.19
    lc50 <- log(86.3)
    label("25(OH)D3 concentration at half-maximal clearance C50 (log nmol/L)") # Huang 2021 Table 1: C50 posterior 86.3 (4.23) nmol/L
    lhill <- log(5.64)
    label("Hill exponent of the sigmoidal 25(OH)D3 clearance (unitless)") # Huang 2021 Table 1: gamma posterior 5.64 (1.28)
    lclmax <- log(0.033)
    label("Maximum 25(OH)D3 clearance CLmax (log L/h)") # Huang 2021 Table 1: CLmax posterior 0.033 (0.0023) L/h
  })

  model({
    # =================================================================
    # 1. PHYSIOLOGY - Huang 2021 supplement Table S6, 'PBPK models which
    #    do not consider adipose compartment (Run003, Run004, Run006,
    #    Run009a, Run009b, Run010)', lumped from Brown 1997 for a 70 kg
    #    adult (Table S5). Footnote: 'Vrb = 70L - Vl - Vven - Vart;
    #    Qrb = Qco - Ql.' Identical values in the deposited Final_model.R.
    # =================================================================
    q_co <- 312.000 # L/h  cardiac output
    q_liver <- 70.824 # L/h  liver blood flow (22.7% of cardiac output)
    q_other <- 241.176 # L/h  rest-of-body blood flow = Qco - Ql
    v_venous <- 4.20 # L
    v_liver <- 1.80 # L
    v_other <- 62.60 # L  = 70 - 1.80 - 4.20 - 1.40
    v_arterial <- 1.40 # L

    # =================================================================
    # 2. PARAMETERS
    # =================================================================
    ka <- exp(lka)
    kp_liver <- exp(lkp_liver)
    kp_other <- exp(lkp_other)
    sclh <- exp(lsclh)
    mclh <- exp(lmclh)
    kp_liver_25d3 <- exp(lkp_liver_25d3)
    kp_other_25d3 <- exp(lkp_other_25d3)
    c50 <- exp(lc50)
    hill <- exp(lhill)
    clmax <- exp(lclmax)

    # Supplement: 'For single dose: CLH = SCLH. For repeated dose:
    # CLH = MCLH.'
    clh <- REGI_QD * mclh + (1 - REGI_QD) * sclh

    # =================================================================
    # 3. ENDOGENOUS VITAMIN D3 INPUT (supplement Final Model, first
    #    equation; script term ... *D25BASE*3). The factor 3 is 1/Fm, so
    #    25(OH)D3 formation from ENDOG balances 25(OH)D3 clearance at the
    #    baseline concentration.
    # =================================================================
    endog <- clmax * D25OH_BL^hill / (c50^hill + D25OH_BL^hill) *
      D25OH_BL / fm_25d3

    # =================================================================
    # 4. VITAMIN D3 (supplement Final Model ODEs 1-5)
    # =================================================================
    cliver <- liver / (v_liver * kp_liver)
    cother <- other / (v_other * kp_other)
    cvenous <- venous / v_venous
    carterial <- arterial / v_arterial

    d/dt(depot) <- -ka * depot
    d/dt(venous) <- q_liver * cliver + q_other * cother - q_co * cvenous
    d/dt(liver) <- q_liver * carterial - q_liver * cliver - clh * cliver +
      endog + ka * depot
    d/dt(other) <- q_other * carterial - q_other * cother
    d/dt(arterial) <- q_co * cvenous - q_liver * carterial - q_other * carterial

    # =================================================================
    # 5. 25(OH)D3 (supplement Final Model ODEs 6-9). The supplement's
    #    clearance term divides the liver 25(OH)D3 amount by Kpl; the
    #    deposited script divides by Kp25l. Both are 1, so the choice is
    #    numerically immaterial; Kp25l is used as in the script.
    # =================================================================
    c25_liver <- liver_25d3 / (v_liver * kp_liver_25d3)
    c25_other <- other_25d3 / (v_other * kp_other_25d3)
    c25_venous <- venous_25d3 / v_venous
    c25_arterial <- arterial_25d3 / v_arterial

    cl_25d3 <- clmax * c25_liver^hill / (c50^hill + c25_liver^hill)

    d/dt(venous_25d3) <- q_liver * c25_liver + q_other * c25_other -
      q_co * c25_venous
    d/dt(liver_25d3) <- q_liver * c25_arterial - q_liver * c25_liver +
      fm_25d3 * clh * cliver - cl_25d3 * c25_liver
    d/dt(other_25d3) <- q_other * c25_arterial - q_other * c25_other
    d/dt(arterial_25d3) <- q_co * c25_venous - q_liver * c25_arterial -
      q_other * c25_arterial

    # =================================================================
    # 6. INITIAL CONDITIONS (supplement Final Model and Final_model.R):
    #    the pre-dose steady state sustained by ENDOG. The supplement
    #    prints A_l(0) with Vven; the script and the steady-state balance
    #    both require Vl (see the vignette Errata).
    # =================================================================
    depot(0) <- 0
    venous(0) <- endog / clh * v_venous
    liver(0) <- endog / clh * v_liver * kp_liver
    other(0) <- endog / clh * v_other * kp_other
    arterial(0) <- endog / clh * v_arterial
    venous_25d3(0) <- D25OH_BL * v_venous
    liver_25d3(0) <- D25OH_BL * v_liver * kp_liver_25d3
    other_25d3(0) <- D25OH_BL * v_other * kp_other_25d3
    arterial_25d3(0) <- D25OH_BL * v_arterial

    # =================================================================
    # 7. OBSERVATIONS - venous plasma concentrations (script output
    #    CP = A25_ven/Vven for 25(OH)D3; the vitamin D3 analogue for the
    #    Figure S4 vitamin D3 data). Deterministic: no residual error
    #    model was reported.
    # =================================================================
    Cc <- cvenous
    Cc_25d3 <- c25_venous
  })
}
