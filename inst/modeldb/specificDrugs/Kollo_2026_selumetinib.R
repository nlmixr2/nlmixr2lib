Kollo_2026_selumetinib <- function() {
  description <- paste(
    "Nonparametric population PK model for oral selumetinib in children aged 5-16 years with inoperable",
    "neurofibromatosis type 1 or plexiform neurofibromas (Kollo 2026, model C).",
    "The authors' 'three-compartment' model is an absorption (depot) compartment feeding a central",
    "compartment with one peripheral compartment, parameterised in micro-constants: first-order absorption",
    "(Ka) after a lag time (Tlag) with relative bioavailability (F), first-order elimination (Ke) from the",
    "central compartment, and bidirectional distribution (KCP central-to-peripheral, KPC",
    "peripheral-to-central). Total body weight enters linearly on the central volume (V = V0 * WT), which",
    "outperformed allometric WT/70, BMI, BSA and age scaling.",
    "The model was fitted with the nonparametric adaptive grid (NPAG) algorithm of Pmetrics, so",
    "between-subject variability is a DISCRETE joint distribution of 22 posterior support points",
    "(Table S12) rather than a parametric omega matrix; the ini() values below are the population mean",
    "parameter vector (Table S14) and this file therefore carries NO inter-individual variability.",
    "The accompanying vignette reproduces the nonparametric distribution from the published support points."
  )
  reference <- paste(
    "Kollo Z, Kovacs J, Neely MN, Vasarhelyi B, Bruckner E, Szabo AJ, Garami M, Karvaly GB.",
    "A nonparametric population pharmacokinetic model of selumetinib in pediatric patients diagnosed with",
    "neurofibromatosis-I or plexiform neurofibromas.",
    "CPT Pharmacometrics Syst Pharmacol. 2026;15(1):e70156. doi:10.1002/psp4.70156"
  )
  vignette <- "Kollo_2026_selumetinib"
  units <- list(time = "h", dosing = "mg", concentration = "nmol/L")

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Baseline body weight measured at the sampling visit (Methods 2.2; individual values in Table S1).",
        "Enters LINEARLY and without a reference weight on the central volume: V = V0 * WT",
        "(Figure 4; Table S4 model C). Because the scaling is multiplicative rather than divisive there is",
        "no normalisation constant, so exp(lvc) carries units of L/kg. This linear form (model C) gave the",
        "best statistical performance of the 11 models compared; allometric WT/70 (model D), BMI (E),",
        "BSA (F) and age (B) scaling all performed worse (Tables S9 and S10)."
      ),
      source_name        = "WT"
    )
  )

  compartmentData <- list(
    depot       = list(analyte = "selumetinib", units = "nmol", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "selumetinib", units = "nmol", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "selumetinib", units = "nmol", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 22,
    n_studies      = 1,
    age_range      = "61-192 months (5.1-16.0 years)",
    age_median     = "140.5 months (11.7 years)",
    weight_range   = "17.0-70.0 kg",
    weight_median  = "39.35 kg",
    sex_female_pct = 45.5,
    disease_state  = "inoperable neurofibromatosis type I (NF-1) or plexiform neurofibromas",
    dose_range     = paste(
      "10-35 mg orally per administration, q12h fixed or morning/evening alternating doses;",
      "initial doses set from body surface area per the Koselugo label and later individualised on",
      "clinical response and tolerability (Table S1)"
    ),
    regions        = "Hungary (single centre, Semmelweis University, Budapest)",
    notes          = paste(
      "Twenty-eight children were recruited (July 2023-October 2024); 22 contributed 156 steady-state",
      "selumetinib concentrations (5-8 samples each) to model construction, and 6 further patients",
      "(S1-S5 and #18) contributed 10 samples to a limited external validation (Table S5).",
      "Patient #18 was excluded from model building because all of that patient's concentrations fell in",
      "the absorption phase. Baseline demographics are in Table 1 and Table S1; laboratory values in",
      "Table S6. All participants had normal liver and renal function and took no chronic co-medication.",
      "The cohort covers 60 percent of the Hungarian paediatric NF-1 population receiving selumetinib.",
      "NOTE: weight_median is 39.35 kg computed from the individual weights in Table S1; the Results text",
      "states 33.4 kg, which cannot be correct because it lies below both sex-specific medians reported in",
      "Table 1 (males 44.3 kg, females 36.1 kg). See the vignette Errata."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Population mean parameter vector of the final nonparametric model
    # (model C), Table S14 "Mean" column. Each value was independently
    # confirmed to equal the probability-weighted mean of the 22 posterior
    # support points in Table S12 (agreement to 3-4 significant figures;
    # the support-point probabilities sum to 1.000004).
    #
    # These are population MEANS of a discrete, markedly non-Gaussian and
    # partly multimodal distribution -- they are not the mode, and for the
    # skewed parameters they differ substantially from the published medians
    # (Table S14 "Median": Ka 18.30, Ke 0.49, V0 0.39, KCP 0.47, KPC 0.13,
    # Tlag 0.58, F 0.89). See the vignette for the full support-point table.
    # ------------------------------------------------------------------

    lka <- log(22.87)
    label("Absorption rate constant (Ka, 1/h)") # Table S14, Ka mean

    lkel <- log(0.58)
    label("Elimination rate constant from the central compartment (Ke, 1/h)") # Table S14, Ke mean

    lvc <- log(0.50)
    label("Central volume of distribution per kg of total body weight (V0, L/kg)") # Table S14, V0 mean

    lk12 <- log(1.32)
    label("Central-to-peripheral distribution rate constant (KCP, 1/h)") # Table S14, KCP mean

    lk21 <- log(0.78)
    label("Peripheral-to-central distribution rate constant (KPC, 1/h)") # Table S14, KPC mean

    ltlag <- log(0.62)
    label("Absorption lag time (Tlag1, h)") # Table S14, Tlag1 mean

    lfdepot <- log(0.81)
    label("Relative bioavailability (FA1, unitless)") # Table S14, FA1 mean

    # ------------------------------------------------------------------
    # Residual error.
    #
    # Pmetrics weights each observation by the assay-error polynomial
    # SD = 0.0372 * concentration + 0.879 nmol/L, established experimentally
    # by spiking 80 samples at 21 levels (Results 3.1 and Figure 5E; the
    # underlying data are Table S2). The final model used the MULTIPLICATIVE
    # error model, whose converged gamma was 2.474 (Results 3.1), so the total
    # residual SD is
    #     SD_total = gamma * (0.879 + 0.0372 * C)
    #              = 2.474 * 0.879 + 2.474 * 0.0372 * C
    #              = 2.1746 + 0.09203 * C   (nmol/L)
    # which is exactly an additive-plus-proportional error model in nlmixr2.
    # Both terms are therefore products of paper-reported numbers, not
    # separately reported estimates.
    # ------------------------------------------------------------------

    addSd <- 2.174646
    label("Additive residual error (nmol/L)") # gamma 2.474 x assay-error intercept 0.879 nmol/L

    propSd <- 0.0920328
    label("Proportional residual error (fraction)") # gamma 2.474 x assay-error slope 0.0372
  })

  model({
    # Selumetinib was assayed in nmol/L whereas doses are administered in mg,
    # so the amount entering the depot is converted from mg to nmol. The
    # molecular weight is NOT reported in the paper: 457.7 g/mol for
    # C17H15BrClFN4O3, from PubChem CID 10127622 (non-paper provenance; see
    # the vignette Errata). 1 mg selumetinib = 1e6 / 457.7 = 2184.8 nmol.
    mwSelumetinib <- 457.7
    mgToNmol <- 1e6 / mwSelumetinib

    # Individual parameters. No eta terms: the source model is nonparametric,
    # so between-subject variability is the discrete support-point
    # distribution of Table S12 and cannot be written as an omega matrix.
    # Imposing a log-normal omega here would contradict the paper's central
    # finding that parametric distributions fail to capture selumetinib
    # variability in this population.
    ka <- exp(lka)
    kel <- exp(lkel)
    k12 <- exp(lk12)
    k21 <- exp(lk21)
    tlag <- exp(ltlag)
    fdepot <- exp(lfdepot)

    # Total body weight enters linearly on the central volume: V = V0 * WT
    # (Figure 4; Table S4 model C). No reference weight is used because the
    # published scaling is multiplicative, not divisive.
    vc <- exp(lvc) * WT

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Relative bioavailability, carrying the mg -> nmol dose conversion.
    f(depot) <- fdepot * mgToNmol
    alag(depot) <- tlag

    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
