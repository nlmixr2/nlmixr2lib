Gupta_2006_peginterferon_alfa_2b <- function() {
  description <- paste(
    "One-compartment population PK model with first-order subcutaneous",
    "absorption for peginterferon alfa-2b (PEG-Intron) in adult patients",
    "with chronic myelogenous leukaemia (Gupta 2006). Apparent clearance",
    "declines over treatment time via an Emax-type function CL(t) = CL0 /",
    "(1 + (t / T50)^beta) with beta fixed to 1 in the final model, so",
    "CL(t) = CL0 / (1 + t / T50). Cockcroft-Gault creatinine clearance",
    "modifies baseline clearance via a power form. Exponential IIV on",
    "CL0, T50, and V; proportional residual error on plasma concentration."
  )
  reference <- paste(
    "Gupta S, Jen J, Kolz K, Cutler D.",
    "Dose selection and population pharmacokinetics of PEGIntron in",
    "patients with chronic myelogenous leukaemia.",
    "Br J Clin Pharmacol. 2006;63(3):292-299.",
    "doi:10.1111/j.1365-2125.2006.02757.x"
  )
  vignette <- "Gupta_2006_peginterferon_alfa_2b"
  units <- list(time = "day", dosing = "ug", concentration = "ng/mL")

  covariateData <- list(
    CRCL = list(
      description        = "Cockcroft-Gault creatinine clearance in raw mL/min (NOT BSA-normalized). Source column CLCR in Gupta 2006. Reference 120 mL/min matches the value at which Gupta 2006 Table 3 reports the typical CL0 = 44.1 L/day.",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed at baseline (the Methods section states 'The CLcr estimate was stable over the course of the study'). Used in a power form on baseline clearance: (CRCL / 120)^e_crcl_cl. Table 1's reported unit '(ml h-1)' is a typesetting error; the value 113 (53-223) is in mL/min, consistent with Table 3 footer's 'CLcr 120 ml min-1'. Source paper writes the column as CLcr; canonical column is CRCL with the raw-Cockcroft-Gault flavour (cf. Delattre 2010 amikacin and Bi 2017 peginterferon alfa-2a entries).",
      source_name        = "CLCR"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description        = "Body weight at baseline",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened in Gupta 2006 Table 2 (Model 2); deltaOFV = 3.686, p = 0.0549. Did not reach the alpha_ENTRY = 0.01 significance threshold; not retained in the final model."
    ),
    AGE = list(
      description        = "Age at baseline",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened in Gupta 2006 Table 2 (Model 4); deltaOFV = 3.326, p = 0.0682; not retained."
    ),
    SEXF = list(
      description        = "Female sex indicator (1 = female, 0 = male)",
      units              = "indicator",
      type               = "categorical",
      reference_category = "male",
      notes              = "Source paper uses a SEX 0/1 indicator; encoded canonically as SEXF (1 = female). Screened in Gupta 2006 Table 2 (Model 3); deltaOFV = 5.597, p = 0.018. Did not reach the alpha_ENTRY = 0.01 significance threshold; not retained in the final model."
    ),
    CREAT = list(
      description        = "Serum creatinine at baseline",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened in Gupta 2006 Table 2 (Model 5, SrCr); deltaOFV = 1.135, p = 0.287; not retained. Cohort serum creatinine was in the normal range (mean 0.8 mg/dL, range 0.4-1.1)."
    ),
    RACE_OTHER = list(
      description        = "Non-White race indicator (1 = Black/Other, 0 = White)",
      units              = "indicator",
      type               = "categorical",
      reference_category = "White",
      notes              = "Source paper encodes RACE as 0/1 (0 = White, 1 = other) - pooling Black + Other (3 + 19 of 137 patients). Screened in Gupta 2006 Table 2 (Model 7); deltaOFV = 1.858, p = 0.173; not retained."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 137L,
    n_studies      = 1L,
    n_observations = 624L,
    age_range      = "20-75 years (mean 51; Gupta 2006 Table 1)",
    weight_range   = "42-137 kg (mean 74.3; Gupta 2006 Table 1)",
    crcl_range     = "53-223 mL/min (mean 113; Gupta 2006 Table 1; published unit 'ml h-1' is a typesetting error - the values are mL/min)",
    srcr_range     = "0.4-1.1 mg/dL (mean 0.8; Gupta 2006 Table 1; all serum creatinine values within the normal range)",
    sex_female_pct = 43.06,
    race_ethnicity = "115/137 (83.9%) White, 3/137 (2.2%) Black, 19/137 (13.9%) Other (Gupta 2006 Table 1, N1/N2/N3 categories).",
    disease_state  = "Adults with newly diagnosed chronic-phase chronic myelogenous leukaemia (CML).",
    dose_range     = "Subcutaneous PEG-Intron once weekly; initial 6.0 ug/kg/week, allowed reductions to 4.5 or 3.0 ug/kg/week or discontinuation based on individual tolerability. Effective dose range across the cohort 3.0-6.0 ug/kg/week.",
    regions        = "United States (n = 26) and International (n = 111); randomised, multicentre, open-label, parallel-group Phase II/III trial.",
    notes          = "Sparse-sampling design with 624 observations from 137 patients. Pre-dose troughs at treatment weeks 4, 12, 24, 36, 48 plus single post-dose samples on weeks 12 (24 h), 24 (72 h), and 36 (120 h). Serum PEG-Intron assayed by a validated ECL immunoassay (LLOQ 50 pg/mL, linear range 50-2000 pg/mL, assay CV 12%; reference 9 of Gupta 2006). Structural one-compartment first-order absorption model and Ka = 1.9 /day were inherited from earlier intensive-sampling Phase I CML data and from prior CHC popPK analyses (references 8-10 of Gupta 2006)."
  )

  ini({
    # Structural parameters (Gupta 2006 Table 3 final model)
    lcl <- log(44.1)               # Gupta 2006 Table 3: theta_CL0 = 44.1 L/day at CLcr = 120 mL/min
    label("Typical baseline (maximal) apparent clearance CL0 at CLcr 120 mL/min (L/day)")
    lt50 <- log(23.8 * 28)         # Gupta 2006 Table 3: theta_T50 = 23.8 reported in units of 28 days -> 666.4 days
    label("Typical time at which apparent clearance has declined to half of CL0, T50 (days)")
    lvc <- log(149)                # Gupta 2006 Table 3: V = 149 L
    label("Typical apparent central volume V (L)")
    lka <- fixed(log(1.9))         # Gupta 2006 Table 3: Ka = 1.9 /day, FIXED at mean of Phase I CHC data (refs 9-10)
    label("First-order subcutaneous absorption rate constant Ka (1/day)")

    # Covariate effect on baseline clearance (Gupta 2006 Methods covariate equation; Table 3)
    e_crcl_cl <- 0.21              # Gupta 2006 Table 3: theta_CLcr = 0.21
    label("Power exponent of (CRCL / 120 mL/min) on baseline clearance CL0 (unitless)")

    # IIV - exponential / log-normal (Gupta 2006 Methods: eta_CL0 and eta_T50 ~ N(0, omega^2);
    # Table 3 omega(CV) column reports CV%. Internal variance is omega^2 = log(CV^2 + 1).
    # The Methods text explicitly enumerates only eta_CL0 and eta_T50, but Table 3 also
    # reports omega(CV) = 59% for V, so the final model carries IIV on V as well.
    etalcl  ~ 0.1156               # Gupta 2006 Table 3: omega(CV) CL0 = 35%  -> omega^2 = log(0.35^2 + 1) = 0.1156
    etalt50 ~ 0.7833               # Gupta 2006 Table 3: omega(CV) T50 = 109% -> omega^2 = log(1.09^2 + 1) = 0.7833
    etalvc  ~ 0.2986               # Gupta 2006 Table 3: omega(CV) V   = 59%  -> omega^2 = log(0.59^2 + 1) = 0.2986

    # Residual error: Gupta 2006 Table 3 reports sigma_epsilon = 41.4 +/- 84.2% with a
    # header unit '(l day-1)' that is a typesetting artifact (the residual is not a
    # clearance). The Methods text describes a combined 'multiplicative + additive'
    # error structure but Table 3 reports a single residual term. We encode the 41.4
    # value as a proportional residual CV (= 0.414) on plasma concentration - the
    # interpretation that is dimensionally consistent with sparsely-sampled popPK
    # practice and with the magnitude relative to the ECL immunoassay (LLOQ 50 pg/mL,
    # assay CV 12%, linear range 50-2000 pg/mL). The 84.2% is the relative standard
    # error on this estimate (high RSE reflects the sparse-sampling design).
    # See vignette Errata for the full rationale.
    propSd <- 0.414                # Gupta 2006 Table 3: sigma_epsilon = 41.4 (interpreted as proportional CV 41.4%)
    label("Proportional residual error on plasma concentration (fraction)")
  })

  model({
    # Individual baseline clearance CL0 with CLcr covariate effect (Gupta 2006 Methods
    # covariate equation):  ln(CL0_i) = ln(theta_CL0) + theta_CLcr * ln(CRCL_i / 120) + eta_CL0_i
    cl0 <- exp(lcl + etalcl) * (CRCL / 120)^e_crcl_cl

    # Individual half-decline time T50 (days)
    t50 <- exp(lt50 + etalt50)

    # Individual central volume V (L)
    vc <- exp(lvc + etalvc)

    # Absorption rate constant Ka (fixed; no IIV)
    ka <- exp(lka)

    # Time-varying apparent clearance: Gupta 2006 Methods Emax-type decline
    #   CL(t) = CL0 / (1 + (t / T50)^beta), with beta = 1 fixed in the final model.
    # 'time' is the model time variable (days from treatment start), matching the paper's t.
    cl <- cl0 / (1 + time / t50)

    kel <- cl / vc

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # PEG-Intron plasma concentration: dose entered in ug, vc in L => Cc in ug/L = ng/mL.
    # The published assay reports concentrations in pg/mL = 1000 x ng/mL (Gupta 2006
    # Figure 1 caption and Blood collection / assay section).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
