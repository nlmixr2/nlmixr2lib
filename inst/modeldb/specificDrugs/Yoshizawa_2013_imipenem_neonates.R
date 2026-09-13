Yoshizawa_2013_imipenem_neonates <- function() {
  description <- paste(
    "One-compartment IV population PK model for imipenem in 60 Japanese",
    "neonates (Yoshizawa 2013). Total clearance is the sum of a renal and a",
    "non-renal arm, both expressed per kilogram of body weight, and the",
    "volume of distribution is likewise per kilogram; body weight therefore",
    "enters every disposition parameter as a linear scaler. Blood and urine",
    "were both assayed, which is what identifies the renal / non-renal",
    "split. Inter-individual variability is exponential on each clearance",
    "arm and residual error is combined proportional plus additive.",
    "The companion model for the same study's paediatric cohort is",
    "Yoshizawa_2013_imipenem_children.",
    "Parameters transcribed from the Zhang 2025 imipenem population-PK",
    "systematic review (Tables 1-3 and Supplementary Table S1), not from the",
    "primary publication; re-verify against Yoshizawa 2013 when the primary",
    "is obtained.",
    sep = " "
  )
  reference <- paste(
    "Yoshizawa K, Ikawa K, Ikeda K, Ohge H, Morikawa N.",
    "Population pharmacokinetic-pharmacodynamic target attainment analysis",
    "of imipenem plasma and urine data in neonates and children.",
    "Pediatr Infect Dis J. 2013;32(11):1208-1216.",
    "doi:10.1097/INF.0b013e31829b5880.",
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

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Body weight is not a screened-and-retained covariate in this model;",
        "it is a structural normalisation. Zhang 2025 Table 2 reports every",
        "disposition parameter for this cohort per kilogram -- CLr in",
        "L/h/kg, CLnr in L/h/kg and V in L/kg -- so each is multiplied by",
        "body weight to give the absolute parameter. That is a linear",
        "(exponent 1) weight scaling on clearance as well as on volume,",
        "which differs from the allometric 0.75 exponent used by the other",
        "paediatric imipenem models in this review (Dong 2019, Dao 2022);",
        "the per-kilogram form is what the review prints and is reproduced",
        "verbatim. Zhang 2025 Table 3 records the covariate analysis as 'NR'",
        "for this study and lists only age, gender and dose as screened, so",
        "no covariate effect beyond this normalisation can be reconstructed.",
        "Cohort mean weight 2.93 +/- 0.7 kg (Zhang 2025 Table 1)."
      ),
      source_name        = "BW"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 60L,
    n_studies        = 1L,
    age_mean         = "0.0288 +/- 0.0227 years (mean +/- SD), i.e. about 10.5 +/- 8.3 days",
    weight_mean      = "2.93 +/- 0.7 kg (mean +/- SD)",
    sex_female_pct   = 50.0,
    race_ethnicity   = NULL,
    disease_state    = "Neonates receiving imipenem-cilastatin",
    dose_range       = paste(
      "10.0-37.3 mg/kg imipenem intravenously (Zhang 2025 Supplementary",
      "Table S1). Neither the dosing interval nor the infusion duration is",
      "reported by the review."
    ),
    regions          = "Japan",
    n_concentrations = 443L,
    notes            = paste(
      "Retrospective study (Zhang 2025 Table 1, study 3). The neonatal",
      "cohort contributed 335 blood and 108 urinary samples (443 total);",
      "sex split 30 male / 30 female. Sampling times and assay method are",
      "recorded as 'NR' in Zhang 2025 Supplementary Table S1. Fitted in",
      "NONMEM and evaluated by goodness-of-fit plots, parameter sensitivity",
      "and leverage analyses, and VPC (Zhang 2025 Table 2).",
      "ERRATUM. The Zhang 2025 Table 2 footnote glosses the two clearance",
      "symbols as 'CLr, clearance of plasma drug; CLnr, clearance of urine",
      "drug'. That gloss is internally incoherent -- it makes the 'r' and",
      "'nr' subscripts meaningless and inverts the matrix each symbol would",
      "name. The subscripts, the fact that both blood and urine were",
      "assayed, and the universal convention in paediatric beta-lactam popPK",
      "all give the intended reading: CLr is RENAL clearance and CLnr is",
      "NON-RENAL clearance, summing to total clearance. The primary's own",
      "title, recovered from the review's reference list, settles it:",
      "'Population pharmacokinetic-pharmacodynamic target attainment",
      "analysis of imipenem PLASMA AND URINE DATA in neonates and",
      "children' -- a urine-informed renal clearance arm is exactly",
      "what that describes. That reading is used",
      "here. Because only the sum enters the plasma model, the split affects",
      "the interpretation of the two etas but not the predicted plasma",
      "concentrations.",
      "ALL PARAMETER VALUES ARE SECONDARY. They come from the Zhang 2025",
      "review's summary tables, not from Yoshizawa 2013 itself."
    )
  )

  ini({
    # ===== Structural PK -- Zhang 2025 Table 2, Yoshizawa et al. (2013)
    # row, 'Neonates' block. All three parameters are per kilogram. =====
    lcl_renal  <- log(0.0783); label("Renal clearance per kg (L/h/kg)")            # Zhang 2025 Table 2 (Yoshizawa 2013, neonates): CLr = 0.0783 L/h/kg
    lcl_nonren <- log(0.138);  label("Non-renal clearance per kg (L/h/kg)")        # Zhang 2025 Table 2 (Yoshizawa 2013, neonates): CLnr = 0.138 L/h/kg
    lvc        <- log(0.466);  label("Volume of distribution per kg (L/kg)")       # Zhang 2025 Table 2 (Yoshizawa 2013, neonates): V = 0.466 L/kg

    # Linear (exponent 1) body-weight scaling. Encoded as a fixed exponent
    # so the structural assumption is visible in the metadata rather than
    # buried in the per-kilogram units. See covariateData$WT.
    e_wt_cl <- fixed(1); label("Linear WT exponent on both clearance arms (unitless)")  # Zhang 2025 Table 2 (Yoshizawa 2013, neonates): CLr and CLnr reported in L/h/kg
    e_wt_vc <- fixed(1); label("Linear WT exponent on Vc (unitless)")                   # Zhang 2025 Table 2 (Yoshizawa 2013, neonates): V reported in L/kg

    # ===== Inter-individual variability =====
    # SCALE CONVENTION. Zhang 2025 Table 2 prints IIV as a bare percentage
    # per parameter without stating the convention, and the column mixes at
    # least three conventions across the review's constituent studies (the
    # audit against the four already-extracted primaries is in the
    # vignette's 'Assumptions and deviations' section). Every model
    # transcribed from this review uses the same documented reading: the
    # printed percentage is an apparent CV of a log-normal random effect,
    # so omega^2 = log(1 + CV^2).
    #
    # IIV is reported on the two clearance arms only; the review gives no
    # omega for V in the neonatal cohort, so none is declared. Writing
    # `etalvc ~ fixed(0)` instead would make OMEGA singular and break
    # simulation.
    etalcl_renal  ~ log(1 + 0.392^2)  # Zhang 2025 Table 2 (Yoshizawa 2013, neonates): IIV CLr = 39.2%, read as an apparent CV
    etalcl_nonren ~ log(1 + 0.334^2)  # Zhang 2025 Table 2 (Yoshizawa 2013, neonates): IIV CLnr = 33.4%, read as an apparent CV

    # ===== Residual error =====
    # Combined proportional plus additive. The additive term is printed in
    # mg/L, i.e. on the concentration scale, so it is taken as a standard
    # deviation, which is what nlmixr2's add() expects.
    propSd <- 0.252; label("Proportional residual error (fraction)")  # Zhang 2025 Table 2 (Yoshizawa 2013, neonates): Proportional = 25.2%
    addSd  <- 0.483; label("Additive residual error (mg/L)")          # Zhang 2025 Table 2 (Yoshizawa 2013, neonates): Additive = 0.483 mg/L
  })

  model({
    # ----- Individual PK parameters -----
    # Per-kilogram parameters multiplied by body weight. Total clearance is
    # the sum of the renal and non-renal arms; only the sum drives the
    # plasma model, since no urine compartment is carried.
    cl_renal  <- exp(lcl_renal  + etalcl_renal)  * WT^e_wt_cl
    cl_nonren <- exp(lcl_nonren + etalcl_nonren) * WT^e_wt_cl
    cl        <- cl_renal + cl_nonren
    vc        <- exp(lvc) * WT^e_wt_vc

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
