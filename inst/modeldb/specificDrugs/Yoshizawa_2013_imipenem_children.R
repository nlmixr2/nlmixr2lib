Yoshizawa_2013_imipenem_children <- function() {
  description <- paste(
    "Two-compartment IV population PK model for imipenem in 39 Japanese",
    "children (Yoshizawa 2013). Total clearance is the sum of a renal and a",
    "non-renal arm and every disposition parameter -- both clearance arms,",
    "both volumes and the intercompartmental clearance -- is expressed per",
    "kilogram of body weight, so body weight enters linearly throughout.",
    "Blood and urine were both assayed, which is what identifies the renal /",
    "non-renal split. Inter-individual variability is exponential on the two",
    "clearance arms and the central volume, and residual error is",
    "proportional.",
    "The companion model for the same study's neonatal cohort is",
    "Yoshizawa_2013_imipenem_neonates.",
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
    central     = list(analyte = "imipenem", units = "mg", specimen = "plasma", verified = FALSE),
    peripheral1 = list(analyte = "imipenem", units = "mg", specimen = "plasma", verified = FALSE)
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
        "disposition parameter for this cohort per kilogram -- CLr and CLnr",
        "in L/h/kg, V1 and V2 in L/kg, Q in L/h/kg -- so each is multiplied",
        "by body weight to give the absolute parameter. That is a linear",
        "(exponent 1) weight scaling on clearance as well as on volume,",
        "which differs from the allometric 0.75 exponent used by the other",
        "paediatric imipenem models in this review (Dong 2019, Dao 2022);",
        "the per-kilogram form is what the review prints and is reproduced",
        "verbatim. Zhang 2025 Table 3 records the covariate analysis as 'NR'",
        "for this study and lists only age, gender and dose as screened, so",
        "no covariate effect beyond this normalisation can be reconstructed.",
        "Cohort mean weight 29.5 +/- 10.9 kg (Zhang 2025 Table 1)."
      ),
      source_name        = "BW"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 39L,
    n_studies        = 1L,
    age_mean         = "9.61 +/- 3.16 years (mean +/- SD)",
    weight_mean      = "29.5 +/- 10.9 kg (mean +/- SD)",
    sex_female_pct   = 33.3,
    race_ethnicity   = NULL,
    disease_state    = "Children receiving imipenem-cilastatin",
    dose_range       = paste(
      "8.70-30.0 mg/kg imipenem intravenously (Zhang 2025 Supplementary",
      "Table S1). Neither the dosing interval nor the infusion duration is",
      "reported by the review."
    ),
    regions          = "Japan",
    n_concentrations = 385L,
    notes            = paste(
      "Retrospective study (Zhang 2025 Table 1, study 3). The paediatric",
      "cohort contributed 230 blood and 155 urinary samples (385 total);",
      "sex split 26 male / 13 female. Sampling times and assay method are",
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
      "concentrations. Note also that the renal / non-renal ordering",
      "reverses between the two cohorts of this study: in neonates the",
      "non-renal arm is the larger of the two (0.138 against 0.0783",
      "L/h/kg), while in children the renal arm dominates (0.187 against",
      "0.0711 L/h/kg). That is the direction renal maturation predicts and",
      "is a point in favour of the reading adopted above.",
      "ALL PARAMETER VALUES ARE SECONDARY. They come from the Zhang 2025",
      "review's summary tables, not from Yoshizawa 2013 itself."
    )
  )

  ini({
    # ===== Structural PK -- Zhang 2025 Table 2, Yoshizawa et al. (2013)
    # row, 'Children' block. All five parameters are per kilogram. =====
    lcl_renal  <- log(0.187);  label("Renal clearance per kg (L/h/kg)")                          # Zhang 2025 Table 2 (Yoshizawa 2013, children): CLr = 0.187 L/h/kg
    lcl_nonren <- log(0.0711); label("Non-renal clearance per kg (L/h/kg)")                      # Zhang 2025 Table 2 (Yoshizawa 2013, children): CLnr = 0.0711 L/h/kg
    lvc        <- log(0.203);  label("Central volume of distribution per kg (L/kg)")             # Zhang 2025 Table 2 (Yoshizawa 2013, children): V1 = 0.203 L/kg
    lvp        <- log(0.0569); label("Peripheral volume of distribution per kg (L/kg)")          # Zhang 2025 Table 2 (Yoshizawa 2013, children): V2 = 0.0569 L/kg
    lq         <- log(0.0621); label("Intercompartmental clearance per kg (L/h/kg)")             # Zhang 2025 Table 2 (Yoshizawa 2013, children): Q = 0.0621 L/h/kg

    # Linear (exponent 1) body-weight scaling on every disposition
    # parameter. Encoded as fixed exponents so the structural assumption is
    # visible in the metadata rather than buried in the per-kilogram units.
    e_wt_cl <- fixed(1); label("Linear WT exponent on both clearance arms and Q (unitless)")  # Zhang 2025 Table 2 (Yoshizawa 2013, children): CLr, CLnr and Q reported in L/h/kg
    e_wt_vc <- fixed(1); label("Linear WT exponent on Vc and Vp (unitless)")                  # Zhang 2025 Table 2 (Yoshizawa 2013, children): V1 and V2 reported in L/kg

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
    # IIV is reported on the two clearance arms and V1 only; the review
    # gives no omega for V2 or Q, so none is declared for them.
    etalcl_renal  ~ log(1 + 0.177^2)  # Zhang 2025 Table 2 (Yoshizawa 2013, children): IIV CLr = 17.7%, read as an apparent CV
    etalcl_nonren ~ log(1 + 0.395^2)  # Zhang 2025 Table 2 (Yoshizawa 2013, children): IIV CLnr = 39.5%, read as an apparent CV
    etalvc        ~ log(1 + 0.171^2)  # Zhang 2025 Table 2 (Yoshizawa 2013, children): IIV V1 = 17.1%, read as an apparent CV

    # ===== Residual error =====
    # Proportional only. Unlike the neonatal cohort of the same study, the
    # review reports no additive term for children.
    propSd <- 0.279; label("Proportional residual error (fraction)")  # Zhang 2025 Table 2 (Yoshizawa 2013, children): Proportional = 27.9%
  })

  model({
    # ----- Individual PK parameters -----
    # Per-kilogram parameters multiplied by body weight. Total clearance is
    # the sum of the renal and non-renal arms; only the sum drives the
    # plasma model, since no urine compartment is carried.
    cl_renal  <- exp(lcl_renal  + etalcl_renal)  * WT^e_wt_cl
    cl_nonren <- exp(lcl_nonren + etalcl_nonren) * WT^e_wt_cl
    cl        <- cl_renal + cl_nonren
    vc        <- exp(lvc + etalvc) * WT^e_wt_vc
    vp        <- exp(lvp)          * WT^e_wt_vc
    q         <- exp(lq)           * WT^e_wt_cl

    # ----- Micro-constants -----
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ----- ODE system -----
    # Imipenem-cilastatin given as an IV infusion into the central
    # compartment; the infusion duration comes from the event table's
    # rate / dur column.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # ----- Output -----
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
