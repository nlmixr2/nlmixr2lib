Moore_1999_lamivudine <- function() {
  description <- "One-compartment population PK model for oral lamivudine in HIV-1-infected adults pooled from the NUCA3001 and NUCA3002 phase III trials (Moore 1999); CL/F scales with a Cockcroft-Gault-style renal function index ((140 - AGE)/(CREAT * 100), * 0.85 if female) raised to an estimated power and with linear body weight, V/F and ka carry no covariates"
  reference <- "Moore KHP, Yuen GJ, Hussey EK, Pakes GE, Eron JJ Jr, Bartlett JA. Population pharmacokinetics of lamivudine in adult human immunodeficiency virus-infected patients enrolled in two phase III clinical trials. Antimicrob Agents Chemother. 1999;43(12):3025-3029. doi:10.1128/aac.43.12.3025"
  vignette <- "Moore_1999_lamivudine"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Baseline body weight. Linear scaling (WT/70) on CL/F outside the renal-function power term (Table 2 footnote d). Reference 70 kg.",
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Baseline age, entered into the Cockcroft-Gault-style renal-function index (140 - AGE)/(CREAT * 100). Cohort range 19-64 years, mean 35.7 years (Table 1).",
      source_name        = "AGE"
    ),
    CREAT = list(
      description        = "Serum creatinine",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Baseline serum creatinine in mg/dL (Table 1 mean 1.1 mg/dL, range 0.6-1.7). Enters the Cockcroft-Gault-style renal-function index (140 - AGE)/(CREAT * 100); the * 100 scaling centres the index near 1 for a typical adult with normal renal function (~0.95 at AGE=36, CREAT=1.1).",
      source_name        = "serum creatinine"
    ),
    SEXF = list(
      description        = "Biological sex indicator, 1 = female, 0 = male",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Multiplies the renal-function index by the Cockcroft-Gault female factor 0.85 (Table 2 footnote d female branch). Encoded inline as (1 - 0.15 * SEXF) so the male reference recovers a factor of 1.",
      source_name        = "gender"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 394,
    n_studies      = 2,
    age_range      = "19-64 years (mean 35.7)",
    age_median     = "mean 35.7 years",
    weight_range   = "37.2-138.5 kg (mean 74.2)",
    weight_median  = "mean 74.2 kg",
    sex_female_pct = 13,
    race_ethnicity = c(Caucasian = 62, Black = 15, Hispanic = 22, Other = 1),
    disease_state  = "HIV-1-infected adults; CDC class A 62%, B 32%, C 6%; mean CD4 306 cells/mm3; 81% HIV-1 RNA PCR < 50,000 copies/mL.",
    dose_range     = "Oral lamivudine 150 mg BID (two 75-mg tablets) or 300 mg BID (one 300-mg tablet), administered approximately 8 am and 4 pm; given alone or in combination with zidovudine 200 mg TID.",
    regions        = "United States (47 sites combined across NUCA3001 and NUCA3002).",
    n_observations = 1477,
    renal_function = "Mean calculated Cockcroft-Gault CrCl 97.4 mL/min; only 3 patients had CrCl < 60 mL/min.",
    notes          = "Pooled antiretroviral-therapy-naive (NUCA3001, n=245) and zidovudine-experienced (NUCA3002, n=149) cohorts with adequate dosing and sample-time documentation. Demographics in Table 1; the NUCA3001 / NUCA3002 trial designs are described in references 2 and 5 of the source paper."
  )

  ini({
    # Structural parameters at the reference patient (CRCL index ~1, WT = 70 kg, male)
    lka  <- log(4.65);  label("Absorption rate constant (ka, 1/h)")                     # Table 2 (theta_3 = 4.65 h^-1)
    lcl  <- log(25.1);  label("Apparent oral clearance reference (CL/F, L/h)")          # Table 2 (theta_1 = 25.1 L/h)
    lvc  <- log(128);   label("Apparent volume of distribution (V/F, L)")               # Table 2 (theta_2 = 128 L)

    # Covariate-effect exponent on the renal-function index ((140 - AGE)/(CREAT * 100))
    e_renal_cl <- 0.468; label("Power exponent of renal-function index on CL/F")        # Table 2 (theta_4 = 0.468)

    # Inter-individual variability; %CV = sqrt(omega^2) * 100 per Table 2 (35.2% and 31.1%)
    etalcl ~ 0.124   # Table 2 (eta on CL/F = 0.124, %CV = 35.2)
    etalvc ~ 0.0965  # Table 2 (eta on V/F = 0.0965, %CV = 31.1)

    # Residual error: paper Table 2 footnote f writes the error model as F * EXP(eps1)
    # with eps1 variance 0.163 (40.4% CV); encoded as lognormal residual on Cc.
    expSd <- 0.404; label("Log-normal residual error SD (exp scale)")                   # Table 2 (epsilon = 0.163, 40.4% CV)
  })
  model({
    # Cockcroft-Gault-style renal-function index, with female multiplicative factor 0.85
    # (Table 2 footnote d): male branch uses RFI = (140 - AGE)/(CREAT * 100); female
    # branch multiplies RFI by 0.85 inside the power term.
    rfi <- (140 - AGE) / (CREAT * 100) * (1 - 0.15 * SEXF)

    # PK parameters (no IIV or covariates on ka; CL/F scales with RFI^e_renal_cl * WT/70)
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * rfi^e_renal_cl * (WT / 70)
    vc <- exp(lvc + etalvc)

    kel <- cl / vc

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Concentration: dose in mg, volume in L -> mg/L (= ug/mL); paper assay reports ng/mL
    # so downstream comparison multiplies Cc by 1000.
    Cc <- central / vc
    Cc ~ lnorm(expSd)
  })
}
