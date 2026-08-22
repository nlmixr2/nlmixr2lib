Ohara_2014_warfarin_s <- function() {
  description <- "S-warfarin population PK/PD in Asian (Taiwanese Chinese) patients during warfarin induction therapy: a one-compartment first-order-absorption PK model for S-warfarin drives an indirect-response model in which S-warfarin inhibits synthesis of plasma normal prothrombin (NPT), and NPT depletion drives INR through a nonlinear power model (Ohara 2014). CYP2C9*3 genotype and body surface area are predictors of S-warfarin clearance; VKORC1 -1639G>A and CYP4F2*3 genotypes are predictors of IC50; baseline NPT is a predictor of the INR nonlinearity exponent lambda. The explicit covariate equations and the Chinese-cohort clinical-trial-simulation application are from Shi 2024. Dosing is expressed as the S-warfarin dose, i.e. half the racemic warfarin dose. R-warfarin is not described by this model."
  reference <- paste(
    "Ohara M, Takahashi H, Lee MTM, Wen M-S, Lee T-H, Chuang H-P, Luo C-H,",
    "Arima A, Onozuka A, Nagai R, Shiomi M, Mihara K, Morita T, Chen Y-T.",
    "Determinants of the Over-Anticoagulation Response during Warfarin Initiation",
    "Therapy in Asian Patients Based on Population Pharmacokinetic-Pharmacodynamic",
    "Analyses. PLoS One. 2014;9(8):e105891. doi:10.1371/journal.pone.0105891.",
    "PMID: 25148255. PMCID: PMC4141831.",
    "Structural equations from Methods Eqs 1-3; all parameter estimates, inter-individual",
    "and residual error from Table 2; NPT0 / INR0 baselines and cohort demographics from Table 1.",
    "The explicit covariate-model equations (CL(S), IC50 and lambda), the BSA reference of",
    "1.74 m^2 and the NPT0 centering value of 119 mg/L are taken from the application paper:",
    "Shi K, Deng J. Comparative performance of pharmacogenetics-based warfarin dosing",
    "algorithms in Chinese population: use of a pharmacokinetic/pharmacodynamic model to",
    "explore dosing regimen through clinical trial simulation.",
    "Pharmacogenet Genomics. 2024;34(8):275-284. doi:10.1097/FPC.0000000000000545.",
    "PMCID: PMC11424055 (Table 3 and its footnotes a-d)."
  )
  vignette <- "Ohara_2014_warfarin_s"

  # `npt` is the plasma normal-prothrombin turnover pool of Ohara 2014 Eq 2. It is
  # a concentration-scale biomarker state, not a drug compartment, so it has no
  # canonical equivalent in compartment-names.md.
  paper_specific_compartments <- c("npt")

  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Ohara 2014 reports both S-warfarin plasma concentration and NPT in ug/mL,
  # which is numerically identical to the mg/L used here (and used by Shi 2024
  # Table 3 for the same constants).
  compartmentData <- list(
    depot   = list(analyte = "S-warfarin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "S-warfarin", units = "mg", specimen = "plasma", verified = TRUE),
    npt     = list(analyte = "normal (fully carboxylated) prothrombin", units = "mg/L", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    BSA = list(
      description        = "Body surface area at baseline",
      units              = "m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed at baseline. Power model on CL(S) normalised to the cohort median BSA of 1.74 m^2 (Shi 2024 Table 3 footnote b; Ohara 2014 Table 1 reports BSA 1.74 +/- 0.18 m^2, n = 99). Ohara 2014 selected BSA over body weight and age, which were collinear with it.",
      source_name        = "BSA"
    ),
    CYP2C9_S3_COUNT = list(
      description        = "Count of CYP2C9*3 (rs1057910) reduced-function alleles per subject (0, 1 or 2)",
      units              = "(count, 0/1/2 alleles per subject)",
      type               = "continuous",
      reference_category = "*1/*1 (count 0) is the reference; CL(S) for *1/*1 is 0.240 L/h at BSA 1.74 m^2.",
      notes              = "Ohara 2014 and Shi 2024 both dichotomise this to a carriage flag rather than a per-allele dosage: Shi 2024 Table 3 footnote b states 'CYP2C9*3 = 0 in patients with CYP2C9*1/*1 and CYP2C9*3 = 1 in patients with CYP2C9*1/*3 or CYP2C9*3/*3'. The model therefore derives the carrier indicator as CYP2C9_S3_COUNT > 0, so a *3/*3 homozygote receives the same 0.543x CL(S) multiplier as a *1/*3 heterozygote. Ohara 2014 Table 1 cohort distribution (n = 99): 88 wild-type, 11 heterozygous, 0 homozygous (MAF 0.056).",
      source_name        = "CYP2C9*3"
    ),
    VKORC1_1639G_COUNT = list(
      description        = "Count of VKORC1 -1639G (rs9923231) alleles per subject (0, 1 or 2). The complementary -1639A count is 2 - VKORC1_1639G_COUNT.",
      units              = "(count, 0/1/2 alleles per subject)",
      type               = "continuous",
      reference_category = "-1639 A/A (G count 0) is the reference; typical IC50 for A/A is 0.0725 mg/L.",
      notes              = "The source papers label this indicator 'VKORC1*2', but per Shi 2024 Table 3 footnote c it flags carriage of the -1639G allele: 'VKORC1*2 = 0 in patients with VKORC1-1639 AA and VKORC1*2 = 1 in patients with VKORC1-1639 AG or VKORC1-1639 GG'. The model derives the carrier indicator as VKORC1_1639G_COUNT > 0. Note this is the opposite orientation to the Hamberg 2007 / Xia 2024 parameterisation, which treats G/G as the reference: here the warfarin-sensitive -1639A/A genotype is the reference and a G allele raises IC50 2.07-fold (i.e. reduces sensitivity), which is the same biology expressed from the other end. Ohara 2014 Table 1 cohort distribution (n = 99, reported as wild/hetero/homo for the -1639G>A substitution): 1 G/G, 17 A/G, 81 A/A (A allele frequency 0.904), a strongly A-predominant East-Asian cohort.",
      source_name        = "VKORC1*2 (rs9923231, -1639G>A)"
    ),
    SNP_CYP4F2_RS2108622_T_COUNT = list(
      description        = "Count of CYP4F2 c.1297C>T (rs2108622, p.V433M; the CYP4F2*3 allele) T alleles per subject (0, 1 or 2)",
      units              = "(count, 0/1/2 alleles per subject)",
      type               = "continuous",
      reference_category = "CYP4F2*1/*1 (T count 0) is the reference; typical IC50 for *1/*1 is 0.0725 mg/L.",
      notes              = "CYP4F2 is the vitamin-K1 oxidase; the *3 (433M) variant has reduced activity, raising hepatic vitamin K and so raising the warfarin concentration required for half-maximal inhibition of prothrombin synthesis. Ohara 2014 estimated a 1.30-fold IC50 increase in *3 carriers (Table 2). As with the other two genotypes the source papers dichotomise to carriage: Shi 2024 Table 3 footnote c gives 'CYP4F2*3 = 0' for *1/*1 and 1 otherwise, so the model derives the carrier indicator as SNP_CYP4F2_RS2108622_T_COUNT > 0. Ohara 2014 Table 1 cohort distribution (n = 99): 50 wild-type, 43 heterozygous, 6 homozygous (MAF 0.278). Shi 2024 assumed every simulated Chinese patient to be CYP4F2*1/*1, i.e. this column is 0 throughout their simulations (Table 3 footnote c).",
      source_name        = "CYP4F2*3"
    ),
    NPT_BASE = list(
      description        = "Baseline plasma normal (fully carboxylated) prothrombin concentration measured before the first warfarin dose",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject and load-bearing in three places: it is the initial condition of the npt state, it sets the zero-order synthesis rate kin = kout * NPT_BASE, and it is the covariate on the INR exponent lambda (lambda = 3.48 * exp(0.00588 * (NPT_BASE - 119)); Shi 2024 Table 3 footnote d, whose centering value of 119 mg/L is the cohort median). NPT_BASE is measured data rather than an estimated parameter: Ohara 2014 fitted no inter-individual variance term for it and excluded the three patients whose NPT0 was missing. Ohara 2014 Table 1 reports NPT0 = 118.2 +/- 22.1 ug/mL (= mg/L) in n = 99. Assayed by the carinactivase-1 method (Ohara 2014 Methods), which measures only the fully carboxylated, coagulation-competent fraction of prothrombin and therefore falls with warfarin exposure while total factor II is comparatively unchanged.",
      source_name        = "NPT0"
    ),
    INR_BASE = list(
      description        = "Pre-medication baseline INR measured before the first warfarin dose",
      units              = "(unitless ratio; INR has no units)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject. Enters Ohara 2014 Eq 3 as an additive constant so simulated INR returns to the subject's own baseline when NPT recovers to NPT_BASE. Ohara 2014 Table 1 reports INR0 = 1.05 +/- 0.10 in n = 99. Like NPT_BASE this is measured data, not an estimated parameter.",
      source_name        = "INR_Base (INR0)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 99L,
    n_studies      = 1L,
    age_range      = "not reported (mean 64.5 +/- 15.2 years)",
    age_median     = "64.5 years (mean)",
    weight_range   = "not reported (mean 68.4 +/- 12.4 kg)",
    weight_median  = "68.4 kg (mean)",
    sex_female_pct = 39.4,
    race_ethnicity = "Han Chinese (Taiwan)",
    disease_state  = "Adults starting warfarin induction therapy; indications atrial fibrillation 54.5%, stroke 29.3%, deep vein thrombosis 25.3%, pulmonary embolism 8.1%. Hypertension 66.7%, chronic kidney disease 16.2%, diabetes 19.2%, hepatic disease 10.1%.",
    dose_range     = "Racemic warfarin, mean starting dose 4.34 +/- 0.98 mg/day and mean maintenance dose 2.94 +/- 1.35 mg/day; target INR 2.0-3.0",
    regions        = "Taiwan (Chang Gung Memorial Hospital and Academia Sinica outpatient clinics)",
    cyp2c9_freq    = "CYP2C9*3 (rs1057910) wild/hetero/homo 88/11/0, MAF 0.056 (Ohara 2014 Table 1)",
    vkorc1_freq    = "VKORC1 -1639G>A (rs9923231) 1 G/G, 17 A/G, 81 A/A; A allele frequency 0.904 (Ohara 2014 Table 1)",
    cyp4f2_freq    = "CYP4F2*3 (rs2108622) wild/hetero/homo 50/43/6, MAF 0.278 (Ohara 2014 Table 1)",
    biomarkers     = "Baseline NPT 118.2 +/- 22.1 mg/L; baseline INR 1.05 +/- 0.10; peak INR during induction 2.25 +/- 0.88. 35 of 99 patients (35%) reached INR >= 4 during induction.",
    n_pd_records   = "NPT and INR analyses used n = 96 (three patients excluded for missing baseline NPT); the Cp(S) analysis used all n = 99",
    software       = "NONMEM via Wings for NONMEM 7.2.0; ADVAN6 for the NPT indirect-response model. Bootstrap: 1000 resamples for Cp(S) and INR, 100 for NPT.",
    notes          = paste(
      "Prospective randomized trial of genotype-guided (n = 77) versus standard (n = 22) warfarin",
      "initiation in Taiwan; recruitment September 2009 to December 2013, samples analysed July 2010",
      "to February 2012. See Ohara 2014 Table 1 for baseline demographics. Shi 2024 applied this model",
      "unchanged to simulated mainland-Chinese cohorts (n = 660 trial replication; Shi 2024 Table 1:",
      "age 67.4 +/- 10.1 y, weight 62.2 +/- 12.2 kg, height 161.9 +/- 8.0 cm, VKORC1 -1639 AA 80.2%,",
      "CYP2C9*1/*1 92.0%)."
    )
  )

  ini({
    # =====================================================================
    # PK: S-warfarin, 1 compartment with first-order absorption
    # (Ohara 2014 Eq 1; parameters Ohara 2014 Table 2 unless noted)
    # =====================================================================
    lfdepot <- fixed(log(1))    ; label("Oral bioavailability of S-warfarin (fraction)")        # Ohara 2014 Methods: "F is the bioavailability fixed at 1.0"
    lka     <- fixed(log(2))    ; label("Absorption rate constant (1/h)")                       # Ohara 2014 Methods: "Ka is the absorption rate constant fixed at 2 h-1"
    lvc     <- fixed(log(13.8)) ; label("Apparent volume of distribution of S-warfarin (L)")    # Ohara 2014 Methods: "Vd ... of S-warfarin fixed at 13.8 L" (from their ref 26)
    lcl     <- log(0.240)       ; label("Apparent oral clearance CL(S) at BSA 1.74 m^2, CYP2C9*1/*1 (L/h)")  # Ohara 2014 Table 2: 240 mL/h (95% CI 220, 260)

    # Covariate effects on CL(S) -- multiplicative, per Shi 2024 Table 3 footnote b:
    #   CL(S) = 0.24 * 0.543^CYP2C9*3 * (BSA_individual / 1.74)^2.14
    e_cyp2c9s3_cl <- 0.543 ; label("CL(S) ratio in CYP2C9*3 carriers vs *1/*1 (unitless)")  # Ohara 2014 Table 2 (95% CI 0.374, 0.712); Shi 2024 Table 3
    e_bsa_cl      <- 2.14  ; label("Power exponent for BSA on CL(S) (unitless)")            # Ohara 2014 Table 2 (95% CI 1.12, 3.16); Shi 2024 Table 3

    # =====================================================================
    # PD-1: normal prothrombin (NPT) indirect response, synthesis inhibition
    # (Ohara 2014 Eq 2; parameters Ohara 2014 Table 2)
    # =====================================================================
    lic50 <- log(0.0725) ; label("S-warfarin concentration giving 50% inhibition of NPT synthesis, VKORC1 -1639AA / CYP4F2*1/*1 (mg/L)")  # Ohara 2014 Table 2 (95% CI 0.0631, 0.0819)
    lkout <- log(0.0136) ; label("First-order NPT elimination rate constant (1/h)")                                                       # Ohara 2014 Table 2 (95% CI 0.0121, 0.0151)
    limax <- fixed(log(1)) ; label("Maximum fractional inhibition of NPT synthesis (unitless)")                                           # Ohara 2014 Methods: "IMax ... assumed to be 1.0 (complete inhibition of NPT synthesis)"

    # Covariate effects on IC50 -- multiplicative, per Shi 2024 Table 3 footnote c:
    #   IC50 = 0.0725 * 2.07^VKORC1*2 * 1.30^CYP4F2*3
    e_vkorc1_ic50 <- 2.07 ; label("IC50 ratio in VKORC1 -1639G carriers vs -1639AA (unitless)")  # Ohara 2014 Table 2 (95% CI 1.58, 2.56); Shi 2024 Table 3
    e_cyp4f2_ic50 <- 1.30 ; label("IC50 ratio in CYP4F2*3 carriers vs *1/*1 (unitless)")         # Ohara 2014 Table 2 (95% CI 1.07, 1.53); Shi 2024 Table 3

    # =====================================================================
    # PD-2: INR from fractional NPT depletion
    # (Ohara 2014 Eq 3; parameters Ohara 2014 Table 2)
    # =====================================================================
    linrmax <- fixed(log(5)) ; label("Maximum INR increase above baseline (unitless)")  # Ohara 2014 Table 2 row "INR Max (fixed)" = 5; Methods: set at 5 because 97.3% of observed INR were < 6
    llambda <- log(3.48)     ; label("Exponent of the NPT-depletion / INR relationship at NPT0 = 119 mg/L (unitless)")  # Ohara 2014 Table 2 (95% CI 3.30, 3.66)

    # Covariate effect on lambda -- exponential, per Shi 2024 Table 3 footnote d:
    #   lambda = 3.48 * exp(0.00588 * (NPT0_individual - 119))
    e_npt0_lambda <- 0.00588 ; label("Exponential coefficient for baseline NPT on lambda (L/mg)")  # Ohara 2014 Table 2 (95% CI 0.00304, 0.00872); Shi 2024 Table 3

    # =====================================================================
    # Inter-individual variability -- exponential on all four estimated
    # parameters. Ohara 2014 Table 2 reports these as "Inter-individual
    # error, omega X (%)", i.e. omega itself expressed as a percentage, so
    # the variance is (percent / 100)^2.
    # =====================================================================
    etalcl     ~ 0.15920 ; label("IIV on CL(S) (variance on log scale)")   # Ohara 2014 Table 2: omega CL(S) = 39.9% -> 0.399^2
    etalic50   ~ 0.14822 ; label("IIV on IC50 (variance on log scale)")    # Ohara 2014 Table 2: omega IC50  = 38.5% -> 0.385^2
    etalkout   ~ 0.20794 ; label("IIV on Kout (variance on log scale)")    # Ohara 2014 Table 2: omega Kout  = 45.6% -> 0.456^2
    etallambda ~ 0.05808 ; label("IIV on lambda (variance on log scale)")  # Ohara 2014 Table 2: omega lambda = 24.1% -> 0.241^2

    # =====================================================================
    # Residual unexplained variability -- one endpoint each
    # =====================================================================
    addSd   <- 0.0697 ; label("Additive residual error on S-warfarin plasma concentration (mg/L)")  # Ohara 2014 Table 2, PK residual error sigma (95% CI 0.0676, 0.0718)
    addSd_NPT  <- 12.2   ; label("Additive residual error on NPT (mg/L)")                              # Ohara 2014 Table 2, PD-1 residual error sigma (95% CI -16.6, 41.1)
    propSd_INR <- 0.247  ; label("Proportional residual error on INR (fraction)")                      # Ohara 2014 Table 2, PD-2 residual error sigma = 24.7% (95% CI 23.9, 25.5); "relative error model"
  })

  model({
    # -----------------------------------------------------------------
    # 1. Derived covariate terms
    #
    # All three genotype effects are carriage flags in the source papers,
    # not per-allele dosages: a homozygous variant carries the same
    # multiplier as a heterozygote (Shi 2024 Table 3 footnotes b and c).
    # -----------------------------------------------------------------
    cyp2c9s3Carrier <- (CYP2C9_S3_COUNT > 0)
    vkorc1gCarrier  <- (VKORC1_1639G_COUNT > 0)
    cyp4f2tCarrier  <- (SNP_CYP4F2_RS2108622_T_COUNT > 0)

    # -----------------------------------------------------------------
    # 2. Individual parameters
    # -----------------------------------------------------------------
    ka <- exp(lka)
    vc <- exp(lvc)
    cl <- exp(lcl + etalcl) * e_cyp2c9s3_cl^cyp2c9s3Carrier * (BSA / 1.74)^e_bsa_cl

    ic50 <- exp(lic50 + etalic50) * e_vkorc1_ic50^vkorc1gCarrier * e_cyp4f2_ic50^cyp4f2tCarrier
    kout <- exp(lkout + etalkout)
    imax <- exp(limax)

    inrmax <- exp(linrmax)
    lambda <- exp(llambda + etallambda) * exp(e_npt0_lambda * (NPT_BASE - 119))

    # -----------------------------------------------------------------
    # 3. Micro-constants
    #
    # Ohara 2014 Methods: "Kin is expressed as Kout * NPT0", so NPT is at
    # steady state at its measured baseline before the first dose.
    # -----------------------------------------------------------------
    kel <- cl / vc
    kin <- kout * NPT_BASE

    # -----------------------------------------------------------------
    # 4. ODE system
    #
    # Ohara 2014 Eq 1 is written as the analytic 1-compartment
    # first-order-absorption solution; it is encoded here as the
    # equivalent depot/central ODE pair.
    #
    # The S-warfarin concentration is written inline inside d/dt(npt)
    # rather than routed through the named `Cc` intermediate: a named
    # state expression used inside d/dt() can silently evaluate to zero
    # in the nlmixr2 model-function form.
    # -----------------------------------------------------------------
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central
    d/dt(npt)     <-  kin * (1 - imax * (central / vc) / (ic50 + central / vc)) - kout * npt

    npt(0) <- NPT_BASE

    # -----------------------------------------------------------------
    # 5. Bioavailability
    # -----------------------------------------------------------------
    f(depot) <- exp(lfdepot)

    # -----------------------------------------------------------------
    # 6. Observations and error models
    #
    # Ohara 2014 Eq 3. NPT can never exceed NPT_BASE under this model
    # (synthesis is only ever inhibited and npt starts at NPT_BASE), so
    # the depletion fraction is mathematically confined to [0, 1); the
    # max() guards only against a solver step returning a value a few
    # ulp above NPT_BASE, which would make a fractional power NaN.
    # -----------------------------------------------------------------
    Cc <- central / vc
    NPT <- npt
    nptDepletion <- max(1 - npt / NPT_BASE, 0)
    INR <- INR_BASE + inrmax * nptDepletion^lambda

    Cc  ~ add(addSd)
    NPT ~ add(addSd_NPT)
    INR ~ prop(propSd_INR)
  })
}
