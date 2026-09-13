Lafaurie_2023_imipenem_recovery <- function() {
  description <- paste(
    "One-compartment IV population PK model for imipenem in 16 French",
    "neutropenic adults, AFTER neutropenia recovery (Lafaurie 2023). No",
    "covariates were retained. Inter-individual variability is exponential",
    "on clearance and volume, and residual error is combined proportional",
    "plus additive.",
    "The same study reported a second, separately-parameterised model for",
    "the same patients DURING the neutropenic phase, in which clearance and",
    "volume are both higher and the residual error is less than half as",
    "large; that phase is the companion model",
    "Lafaurie_2023_imipenem_neutropenia. The two are shipped as separate",
    "files because the review reports each with its own structural values,",
    "its own inter-individual variability AND its own residual error -- a",
    "per-phase residual error is the signature of two independent fits",
    "rather than one model carrying a phase covariate.",
    "Parameters transcribed from the Zhang 2025 imipenem population-PK",
    "systematic review (Tables 1-3 and Supplementary Table S1), not from the",
    "primary publication; re-verify against Lafaurie 2023 when the primary",
    "is obtained.",
    sep = " "
  )
  reference <- paste(
    "Lafaurie M, Burdet C, Hammas K, Goldwirt L, Bercot B, Sauvageon H,",
    "et al.",
    "Population pharmacokinetics and pharmacodynamics of imipenem in",
    "neutropenic adult patients.",
    "Infect Dis Now. 2023;53(1):104625. doi:10.1016/j.idnow.2022.09.020.",
    "Parameters transcribed from Zhang P, Zhao Y, Zhu J, Yang Y, Liang G,",
    "Wang X, Yu Z. Population pharmacokinetics of imipenem in different",
    "populations for individualized dosing: a systematic review.",
    "Front Pharmacol. 2025;16:1738055. doi:10.3389/fphar.2025.1738055",
    "(Tables 1-3, Supplementary Table S1).",
    sep = " "
  )
  vignette <- "Zhang_2025_imipenem_model_review"
  units    <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. verified = FALSE because the primary publication is
  # not on disk.
  compartmentData <- list(
    central = list(analyte = "imipenem", units = "mg", specimen = "plasma", verified = FALSE)
  )

  # No covariates. Zhang 2025 Table 3 records "NR" (no record) in every
  # column for this study -- covariate analysis method, covariates
  # screened, covariates incorporated, and formulation are all unreported
  # by the review. This is an absence of reporting in the secondary source
  # rather than positive evidence that the primary screened nothing.
  covariateData <- list()

  population <- list(
    species          = "human",
    n_subjects       = 16L,
    n_studies        = 1L,
    age_median       = "37 years (range 18.3-78.3)",
    weight_median    = "65.5 kg (range 48-101)",
    sex_female_pct   = 56.3,
    race_ethnicity   = NULL,
    disease_state    = paste(
      "Adult patients receiving imipenem-cilastatin, sampled after",
      "recovery from neutropenia. The same patients were sampled during",
      "the neutropenic phase; see the companion model",
      "Lafaurie_2023_imipenem_neutropenia."
    ),
    dose_range       = paste(
      "1000 mg imipenem intravenously every 8 h, or 500 mg every 6 h",
      "(Zhang 2025 Supplementary Table S1). The infusion duration is not",
      "reported by the review."
    ),
    regions          = "France",
    n_concentrations = NA_integer_,
    notes            = paste(
      "Prospective study (Zhang 2025 Table 1, study 16); 16 patients and",
      "118 samples IN TOTAL across both phases, sex split 7 male / 9",
      "female. The review does not report how the 118 samples split",
      "between the neutropenic and recovery phases, so n_concentrations is",
      "left NA for each phase model rather than guessing. Blood was",
      "sampled immediately before the dose and at 0.5, 1, 2, 3, 4, 6 and",
      "8 h afterwards, and assayed by HPLC-UV (Zhang 2025 Supplementary",
      "Table S1). Fitted in Monolix and evaluated by goodness-of-fit plots",
      "and NPDE (Zhang 2025 Table 2).",
      "Clearance after recovery (10.9 L/h) is 24% lower than during",
      "neutropenia (14.3 L/h) and the volume of distribution is 30%",
      "smaller (14.5 against 20.7 L), while the proportional residual",
      "error more than doubles (36.1% against 15.9%). Zhang 2025's",
      "Discussion does not comment on this study's phase contrast",
      "specifically, but the direction matches the augmented renal",
      "clearance and expanded extracellular fluid volume it describes for",
      "acutely unwell patients generally, with both resolving as the",
      "patient recovers.",
      "ALL PARAMETER VALUES ARE SECONDARY. They come from the Zhang 2025",
      "review's summary tables, not from Lafaurie 2023 itself."
    )
  )

  ini({
    # ===== Structural PK -- Zhang 2025 Table 2, Lafaurie et al. (2023)
    # row, 'After neutropenia recovery' block. =====
    lcl <- log(10.9); label("Clearance after neutropenia recovery (L/h)")            # Zhang 2025 Table 2 (Lafaurie 2023, after recovery): CL = 10.9 L/h
    lvc <- log(14.5); label("Volume of distribution after neutropenia recovery (L)") # Zhang 2025 Table 2 (Lafaurie 2023, after recovery): V = 14.5 L

    # ===== Inter-individual variability =====
    # SCALE CONVENTION. Zhang 2025 Table 2 prints IIV as a bare percentage
    # per parameter without stating the convention, and the column mixes at
    # least three conventions across the review's constituent studies (the
    # audit against the four already-extracted primaries is in the
    # vignette's 'Assumptions and deviations' section). Every model
    # transcribed from this review uses the same documented reading: the
    # printed percentage is an apparent CV of a log-normal random effect,
    # so omega^2 = log(1 + CV^2). Note that this study was fitted in
    # Monolix, whose parameter table conventionally reports omega as a
    # STANDARD DEVIATION rather than a CV -- the one other Monolix study
    # in this review whose primary is already extracted here, Couffignal
    # 2014, is exactly the case where the SD reading is correct. That makes
    # this model a likely site of the residual scale disagreement; at
    # 21.8% and 27.8% the two readings differ by 4-7% in variance, so the
    # practical consequence is small.
    etalcl ~ log(1 + 0.218^2)  # Zhang 2025 Table 2 (Lafaurie 2023, after recovery): IIV CL = 21.8%, read as an apparent CV
    etalvc ~ log(1 + 0.278^2)  # Zhang 2025 Table 2 (Lafaurie 2023, after recovery): IIV V = 27.8%, read as an apparent CV

    # ===== Residual error =====
    # Combined proportional plus additive. The additive term is printed in
    # mg/L, i.e. on the concentration scale, so it is taken as a standard
    # deviation, which is what nlmixr2's add() expects.
    propSd <- 0.361; label("Proportional residual error (fraction)")  # Zhang 2025 Table 2 (Lafaurie 2023, after recovery): Proportional = 36.1%
    addSd  <- 0.30;  label("Additive residual error (mg/L)")          # Zhang 2025 Table 2 (Lafaurie 2023, after recovery): Additive = 0.3 mg/L
  })

  model({
    # ----- Individual PK parameters -----
    # No covariates (Zhang 2025 Table 3 records 'NR' throughout for this
    # study), so each parameter is its typical value times an exponential
    # random effect.
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)

    # ----- Micro-constants -----
    kel <- cl / vc

    # ----- ODE system -----
    # Imipenem-cilastatin given as an IV infusion into the central
    # compartment; the infusion duration comes from the event table's
    # rate / dur column.
    d/dt(central) <- -kel * central

    # ----- Output -----
    Cc <- central / vc
    Cc ~ prop(propSd) + add(addSd)
  })
}
