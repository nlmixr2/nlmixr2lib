Landersdorfer_2007_flucloxacillin <- function() {
  description <- paste(
    "Three-compartment population PK model for IV flucloxacillin in",
    "healthy adult volunteers (Landersdorfer 2007) with linear renal",
    "and non-renal elimination. The structural model splits total",
    "clearance into a renal arm (CL_R = 5.37 L/h) and a non-renal arm",
    "(CL_NR = 2.73 L/h); their sum reproduces the derived total",
    "clearance CL_T = 8.10 L/h reported in Table 2. The renal arm also",
    "drives a cumulative urinary excretion compartment that the paper",
    "fits jointly with plasma. Distribution uses a shallow peripheral",
    "(V_2 = 2.61 L, CLic_shallow = 15.3 L/h) and a deep peripheral",
    "(V_3 = 2.17 L, CLic_deep = 1.23 L/h); central volume V_1 = 4.79 L.",
    "Between-subject variability is reported as a full 5x5",
    "variance-covariance matrix (Table 3, natural-log scale) on CL_R,",
    "CL_NR, V_1, V_2, V_3; no BSV is included on the inter-compartmental",
    "clearances. Residual error is combined additive + proportional on",
    "both plasma concentrations (9.4% CV, 0.155 mg/L) and cumulative",
    "urinary amounts (20.9% CV, 1.04 mg). The 5-min infusion duration",
    "used in the study is supplied via dose records (DUR / RATE) rather",
    "than as a model parameter. No structural covariates were retained:",
    "the cohort was 10 healthy Caucasian adults (5 M / 5 F, weight",
    "52-83 kg, age 23-34 years) and demographics are not used inside",
    "the model. Monte Carlo dose-attainment simulations in the paper",
    "(continuous, 4-h, 0.5-h infusions) reuse these PK parameters",
    "together with 96% protein binding."
  )
  reference <- paste(
    "Landersdorfer CB, Kirkpatrick CMJ, Kinzig-Schippers M, Bulitta JB,",
    "Holzgrabe U, Drusano GL, Sorgel F.",
    "Population pharmacokinetics at two dose levels and pharmacodynamic",
    "profiling of flucloxacillin.",
    "Antimicrob Agents Chemother. 2007;51(9):3290-3297.",
    "doi:10.1128/AAC.01410-06"
  )
  vignette <- "Landersdorfer_2007_flucloxacillin"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = 10L,
    n_studies      = 1L,
    age_range      = "23-34 years (median 25)",
    weight_range   = "52-83 kg (median 71)",
    height_range   = "165-190 cm (median 178)",
    sex_female_pct = 50,
    race_ethnicity = "100% Caucasian (paper Methods, Study participants)",
    disease_state  = paste(
      "Healthy adult volunteers with normal renal and hepatic function;",
      "screened by physical exam, ECG, urinalysis, and drugs-of-abuse",
      "panel before entry (Landersdorfer 2007 Methods, Study participants)."
    ),
    dose_range     = paste(
      "Randomised two-way crossover: 500 mg and 1,000 mg single doses",
      "given as 5-min IV infusions, separated by a >=4-day washout",
      "(Landersdorfer 2007 Methods, Study design and drug administration).",
      "Monte Carlo simulation explored continuous-infusion, 4-h, and 0.5-h",
      "infusion regimens at daily doses of 6 g and 12 g (Methods, MCS)."
    ),
    regions        = "Germany (Institute for Biomedical and Pharmaceutical Research, Nuernberg-Heroldsberg)",
    notes          = paste(
      "Baseline demographics from Landersdorfer 2007 Results, Demographics",
      "(median weight 71 kg [52-83]; median height 178 cm [165-190];",
      "median age 25 years [23-34]; n=10, 5 M / 5 F). Plasma and urine",
      "sampling: blood drawn pre-infusion, end of infusion, then 5, 10,",
      "15, 20, 25, 45, 60, 75, 90 min and 2, 2.5, 3, 3.5, 4, 5, 6, 8, 24 h",
      "after end of infusion; urine collected pre-infusion, infusion start",
      "to 1 h after end of infusion, and over 1-2, 2-3, 3-4, 4-5, 5-6,",
      "6-8, 8-12, 12-24 h windows (Methods, Sampling schedule). NONMEM",
      "version V, release 1.1 with FOCE-I (ADVAN 11 / TRANS 4 analytical",
      "solution); the three-compartment model with linear elimination was",
      "preferred over saturable / parallel-saturable alternatives based on",
      "a 200-point OFV improvement and visual predictive checks."
    )
  )

  ini({
    # Structural clearance parameters -- Landersdorfer 2007 Table 2.
    # CL_T (8.10 L/h) is derived as CL_R + CL_NR, not estimated (footnote c).
    lcl_renal  <- log(5.37);  label("Renal clearance CL_R (L/h)")             # Landersdorfer 2007 Table 2 (CL_R)
    lcl_nonren <- log(2.73);  label("Non-renal clearance CL_NR (L/h)")         # Landersdorfer 2007 Table 2 (CL_NR)

    # Structural volume parameters -- Landersdorfer 2007 Table 2.
    # V_ss (9.57 L) is derived as V_1 + V_2 + V_3, not estimated (footnote c).
    lvc  <- log(4.79); label("Central volume V_1 (L)")                         # Landersdorfer 2007 Table 2 (V_1)
    lvp  <- log(2.61); label("Shallow peripheral volume V_2 (L)")              # Landersdorfer 2007 Table 2 (V_2)
    lvp2 <- log(2.17); label("Deep peripheral volume V_3 (L)")                 # Landersdorfer 2007 Table 2 (V_3)

    # Inter-compartmental clearances -- Landersdorfer 2007 Table 2.
    # The paper estimated point values but explicitly did NOT include BSV on
    # these terms (Table 2 footnote d; Methods, Individual PK model).
    lq  <- log(15.3); label("Inter-compartmental clearance central <-> shallow CLic_shallow (L/h)")  # Landersdorfer 2007 Table 2 (CLic_shallow)
    lq2 <- log(1.23); label("Inter-compartmental clearance central <-> deep CLic_deep (L/h)")        # Landersdorfer 2007 Table 2 (CLic_deep)

    # Inter-individual variability -- full 5x5 variance-covariance matrix on
    # natural-log scale (Landersdorfer 2007 Table 3). The paper reports BSV
    # CV% as sqrt(omega^2) directly (Methods, Individual PK model: "the
    # square root of the variance is an approximation of the apparent
    # coefficient of variation of a normal distribution on a log scale"), so
    # the Table 3 diagonal entries are the omega^2 variances to use here.
    # Order matches Table 3: CL_R, CL_NR, V_1, V_2, V_3 ->
    # etalcl_renal, etalcl_nonren, etalvc, etalvp, etalvp2.
    #   diag(omega^2) sanity-check vs Table 2 CV%:
    #     CL_R   sqrt(0.0343) = 0.185 -> 19%   (Table 2: 19%)
    #     CL_NR  sqrt(0.112)  = 0.335 -> 33%   (Table 2: 33%)
    #     V_1    sqrt(0.0282) = 0.168 -> 17%   (Table 2: 17%)
    #     V_2    sqrt(0.138)  = 0.371 -> 37%   (Table 2: 37%)
    #     V_3    sqrt(0.023)  = 0.152 -> 15%   (Table 2: 15%)
    etalcl_renal + etalcl_nonren + etalvc + etalvp + etalvp2 ~
      c(0.0343,
        0.0124,  0.112,
        0.00488, 0.00551, 0.0282,
        0.0415,  0.0641,  0.00168, 0.138,
        0.0172,  0.00804, 0.0175,  0.0351, 0.023)  # Landersdorfer 2007 Table 3

    # Residual error -- combined additive + proportional on the linear scale,
    # one pair per output (Landersdorfer 2007 Table 2, Methods, Observation
    # model).
    propSd          <- 0.094;  label("Proportional residual error on plasma (fraction)")              # Landersdorfer 2007 Table 2 (CV_C)
    addSd           <- 0.155;  label("Additive residual error on plasma (mg/L)")                       # Landersdorfer 2007 Table 2 (SD_C)
    propSd_urineAmt <- 0.209;  label("Proportional residual error on cumulative urinary amount (fraction)")  # Landersdorfer 2007 Table 2 (CV_AU)
    addSd_urineAmt  <- 1.04;   label("Additive residual error on cumulative urinary amount (mg)")            # Landersdorfer 2007 Table 2 (SD_AU)
  })

  model({
    # Individual PK parameters. Renal and non-renal clearance arms are
    # carried separately so the urine compartment can be driven by CL_R
    # alone (Landersdorfer 2007 Methods, Population PK analysis).
    cl_renal  <- exp(lcl_renal  + etalcl_renal)
    cl_nonren <- exp(lcl_nonren + etalcl_nonren)
    cl        <- cl_renal + cl_nonren                                          # Table 2 footnote c: CL_T derived
    vc        <- exp(lvc  + etalvc)
    vp        <- exp(lvp  + etalvp)
    vp2       <- exp(lvp2 + etalvp2)
    q         <- exp(lq)                                                       # Table 2 footnote d: no BSV
    q2        <- exp(lq2)                                                      # Table 2 footnote d: no BSV

    # Micro-constants
    kel  <- cl       / vc
    kelr <- cl_renal / vc
    k12  <- q  / vc
    k21  <- q  / vp
    k13  <- q2 / vc
    k31  <- q2 / vp2

    # ODE system. central, peripheral1 (shallow), peripheral2 (deep), and a
    # cumulative urinary excretion compartment driven by CL_R. Reproduces
    # the four-state system at Landersdorfer 2007 Methods, Population PK
    # analysis (equations for dX(1)/dt through dX(4)/dt).
    d/dt(central)     <- -(kel + k12 + k13) * central + k21 * peripheral1 + k31 * peripheral2
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(peripheral2) <-  k13 * central - k31 * peripheral2
    d/dt(urine)       <-  kelr * central

    # Observations. Dose in mg, vc in L -> Cc in mg/L.
    Cc       <- central / vc
    urineAmt <- urine

    Cc       ~ add(addSd) + prop(propSd)
    urineAmt ~ add(addSd_urineAmt) + prop(propSd_urineAmt)
  })
}
