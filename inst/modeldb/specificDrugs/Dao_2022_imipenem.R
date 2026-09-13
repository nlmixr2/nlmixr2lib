Dao_2022_imipenem <- function() {
  description <- paste(
    "One-compartment IV population PK model for imipenem in 82 Swiss",
    "neonates (Dao 2022). Clearance carries four covariates -- allometric",
    "body weight, and centred linear effects of postnatal age and",
    "gestational age, and an inverse power effect of serum creatinine --",
    "while the volume of distribution scales with body weight using the",
    "same 0.75 exponent as clearance. Inter-individual variability is",
    "exponential on clearance only, and residual error is combined",
    "proportional plus additive.",
    "Parameters transcribed from the Zhang 2025 imipenem population-PK",
    "systematic review (Tables 1-3 and Supplementary Table S1), not from the",
    "primary publication; re-verify against Dao 2022 when the primary is",
    "obtained.",
    sep = " "
  )
  reference <- paste(
    "Dao K, Fuchs A, Andre P, Giannoni E, Decosterd LA, Marchetti O, et al.",
    "Dosing strategies of imipenem in neonates based on pharmacometric",
    "modelling and simulation.",
    "J Antimicrob Chemother. 2022;77(2):457-465. doi:10.1093/jac/dkab394.",
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
        "Reference 1.16 kg, equal to the cohort median (range 0.5-4.1 kg;",
        "Zhang 2025 Table 1) -- a very-low-birth-weight neonatal cohort.",
        "Enters BOTH clearance and volume with the exponent 0.75 (Zhang",
        "2025 Table 3: CL = 0.21 * (BW/1.16)^0.75 * ...,",
        "V = 0.73 * (BW/1.16)^0.75).",
        "NOTE the volume exponent. The conventional allometric pairing is",
        "0.75 on clearance and 1 on volume; this model uses 0.75 on both.",
        "That is what the review prints, and it is reproduced verbatim",
        "rather than 'corrected' to 1. A 0.75 volume exponent is not",
        "unheard of in neonatal popPK but it is unusual enough that it",
        "should be re-verified against the primary, since a transcription",
        "slip in the secondary source is a live alternative. Recorded in",
        "the vignette Errata."
      ),
      source_name        = "BW"
    ),
    PNA = list(
      description        = "Postnatal age (chronological age since birth)",
      units              = "months",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "UNIT REPARAMETERISATION. The source expresses postnatal age in",
        "DAYS, centred on 21 days -- Zhang 2025 Table 1 reports the cohort",
        "PNA as 21 days (range 2.1-153) and Table 3's formula reads",
        "(1 + 0.22 * (PNA - 21) / 21). The canonical PNA column carries",
        "MONTHS per inst/references/covariate-columns.md, so model()",
        "recovers days with pna_days <- PNA * 30.4375 before forming the",
        "centred term, leaving the printed constants 0.22 and 21 exactly",
        "as the review gives them. This is the same reparameterisation",
        "that register entry records for Zhao 2018 (days) and Bardhi 2026",
        "(hours). The reference 21 days corresponds to 0.6899 months.",
        "The effect is a CENTRED LINEAR multiplier, not a power term, so",
        "it can in principle go non-positive: the factor reaches zero at",
        "PNA = 21 * (1 - 1/0.22) = -74.5 days, which is outside any",
        "achievable postnatal age, so the term is safe over its whole",
        "physical domain. Time-varying within subject."
      ),
      source_name        = "PNA"
    ),
    GA = list(
      description        = "Gestational age at birth",
      units              = "weeks",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Reference 26.9 weeks, equal to the cohort median (range 24.2-41.3",
        "weeks; Zhang 2025 Table 1) -- an extremely preterm cohort. Enters",
        "clearance as the centred linear multiplier",
        "(1 + 1.31 * (GA - 26.9) / 26.9) (Zhang 2025 Table 3).",
        "MAGNITUDE AND DOMAIN WARNING. The coefficient 1.31 is large: a",
        "term-born neonate at 40 weeks gets a factor of",
        "1 + 1.31 * (40 - 26.9)/26.9 = 1.64, i.e. 64% higher clearance than",
        "the 26.9-week reference, on top of the allometric weight effect.",
        "Because the form is linear rather than a power, the factor also",
        "reaches zero at GA = 26.9 * (1 - 1/1.31) = 6.4 weeks and would go",
        "NEGATIVE below that -- far outside any viable gestational age and",
        "outside the fitted range of 24.2-41.3 weeks, but a caller",
        "supplying nonsense would get a nonsense (or negative) clearance",
        "rather than an error. Do not extrapolate this term below the",
        "fitted range. Recorded in the vignette Errata. Time-fixed per",
        "subject."
      ),
      source_name        = "GA"
    ),
    CREAT = list(
      description        = "Serum creatinine",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Reference 46.6. The review does not state the unit, but 46.6",
        "umol/L (= 0.53 mg/dL) is a physiological neonatal serum",
        "creatinine while 46.6 mg/dL is not survivable, and the study is",
        "Swiss, where SI units are standard -- so umol/L is the only",
        "coherent reading. Recorded here explicitly because the canonical",
        "CREAT column accepts either umol/L or mg/dL and requires the",
        "per-model unit to be documented.",
        "Enters clearance as the INVERSE power term (46.6/SCr)^0.2 --",
        "note the reference is in the NUMERATOR, so higher creatinine",
        "lowers clearance, which is the expected direction for a renally",
        "eliminated drug. Zhang 2025 Results names this study as the only",
        "one of the 18 to retain serum creatinine itself (rather than a",
        "derived creatinine clearance) as a covariate. Time-varying."
      ),
      source_name        = "SCr"
    )
  )

  # Screened during covariate model building and not retained (Zhang 2025
  # Table 3, 'Covariates screened' column). The review prints no
  # coefficient for any of them, so none can be reconstructed.
  covariatesDataExcluded <- list(
    SEXF = list(description = "Female sex",           units = "(binary)", type = "binary",     notes = "Screened, not retained (Zhang 2025 Table 3). Cohort 44/82 female (Table 1)."),
    PMA  = list(description = "Postmenstrual age",    units = "weeks",    type = "continuous", notes = "Screened, not retained (Zhang 2025 Table 3); gestational age and postnatal age were retained separately instead. Cohort median 31 weeks, range 25.6-48.3 (Table 1). Not a registered canonical in inst/references/covariate-columns.md; recorded here as documentation only, since it is never referenced in model()."),
    SGA  = list(description = "Small for gestational age", units = "(binary)", type = "binary", notes = "Screened, not retained (Zhang 2025 Table 3). Not a registered canonical; recorded here as documentation only.")
  )

  population <- list(
    species          = "human",
    n_subjects       = 82L,
    n_studies        = 1L,
    age_median       = paste(
      "postnatal age 21 days (range 2.1-153); gestational age 26.9 weeks",
      "(range 24.2-41.3); postmenstrual age 31 weeks (range 25.6-48.3)"
    ),
    ga_range         = "24.2-41.3 weeks (median 26.9)",
    weight_median    = "1.16 kg (range 0.5-4.1)",
    sex_female_pct   = 53.7,
    race_ethnicity   = NULL,
    disease_state    = paste(
      "Neonates, predominantly extremely preterm, receiving",
      "imipenem-cilastatin. Together with the neonatal arm of Yoshizawa",
      "2013 this is one of only two neonatal cohorts among the 18 studies",
      "in the Zhang 2025 review, and it is the more preterm of the two."
    ),
    dose_range       = paste(
      "15-20 mg/kg imipenem intravenously every 8-12 h (Zhang 2025",
      "Supplementary Table S1). The infusion duration is not reported by",
      "the review."
    ),
    regions          = "Switzerland",
    n_concentrations = 173L,
    notes            = paste(
      "Retrospective study (Zhang 2025 Table 1, study 13); 82 patients,",
      "173 samples, sex split 38 male / 44 female. Samples were taken at",
      "Cmax (1-2 h after the start of infusion), at Cmin (at steady state,",
      "generally before the fourth dose), or both, and assayed by LC-MS/MS",
      "-- one of only three studies in the review not using HPLC-UV (Zhang",
      "2025 Supplementary Table S1). Covariates were selected by forward",
      "inclusion and backward elimination (thresholds not reported by the",
      "review). Concomitant treatments were also screened and not retained.",
      "Fitted in NONMEM and evaluated by bootstrap, pcVPC and",
      "goodness-of-fit plots (Zhang 2025 Table 2).",
      "The two age covariates enter as CENTRED LINEAR multipliers rather",
      "than power terms, which means both can in principle drive clearance",
      "to zero or below outside the fitted range; see the domain warnings",
      "in covariateData$GA and covariateData$PNA.",
      "ALL PARAMETER VALUES ARE SECONDARY. They come from the Zhang 2025",
      "review's summary tables, not from Dao 2022 itself."
    )
  )

  ini({
    # ===== Structural PK -- Zhang 2025 Table 2, Dao et al. (2022) row.
    # Typical values for the reference neonate: 1.16 kg, postnatal age 21
    # days, gestational age 26.9 weeks, serum creatinine 46.6 umol/L
    # (Zhang 2025 Table 3 formulations). =====
    lcl <- log(0.21); label("Clearance at the reference neonate (L/h)")               # Zhang 2025 Table 2 (Dao 2022): CL = 0.21 L/h; leading coefficient of the Table 3 CL formula
    lvc <- log(0.73); label("Volume of distribution at 1.16 kg (L)")                  # Zhang 2025 Table 2 (Dao 2022): V = 0.73 L; leading coefficient of the Table 3 V formula

    # ===== Covariate effects -- Zhang 2025 Table 3, Dao et al. (2022) =====
    #   CL = 0.21 * (BW/1.16)^0.75
    #             * (1 + 0.22 * (PNA - 21)/21)
    #             * (1 + 1.31 * (GA - 26.9)/26.9)
    #             * (46.6/SCr)^0.2
    #   V  = 0.73 * (BW/1.16)^0.75
    # The weight exponent 0.75 is the theory-based allometric value and may
    # well have been fixed rather than estimated, but the review reports no
    # standard errors and does not say, so it is left UNWRAPPED rather than
    # asserting a fixed() provenance the secondary source does not support.
    # Note that the SAME exponent is used on volume -- see covariateData$WT.
    e_wt_cl    <- 0.75; label("Allometric WT exponent on CL (unitless)")                     # Zhang 2025 Table 3 (Dao 2022): (BW/1.16)^0.75 in the CL formula
    e_wt_vc    <- 0.75; label("WT exponent on Vc (unitless)")                                # Zhang 2025 Table 3 (Dao 2022): (BW/1.16)^0.75 in the V formula
    e_pna_cl   <- 0.22; label("Centred linear coefficient of PNA on CL (unitless)")          # Zhang 2025 Table 3 (Dao 2022): (1 + 0.22 * (PNA - 21)/21), PNA in days
    e_ga_cl    <- 1.31; label("Centred linear coefficient of GA on CL (unitless)")           # Zhang 2025 Table 3 (Dao 2022): (1 + 1.31 * (GA - 26.9)/26.9), GA in weeks
    e_creat_cl <- 0.2;  label("Inverse power exponent of serum creatinine on CL (unitless)") # Zhang 2025 Table 3 (Dao 2022): (46.6/SCr)^0.2

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
    # IIV is reported on CL only; the review gives no omega for V.
    etalcl ~ log(1 + 0.20^2)  # Zhang 2025 Table 2 (Dao 2022): IIV CL = 20%, read as an apparent CV

    # ===== Residual error =====
    # Combined proportional plus additive. The additive term is printed in
    # mg/L, i.e. on the concentration scale, so it is taken as a standard
    # deviation, which is what nlmixr2's add() expects.
    propSd <- 0.37; label("Proportional residual error (fraction)")  # Zhang 2025 Table 2 (Dao 2022): Proportional = 37%
    addSd  <- 0.04; label("Additive residual error (mg/L)")          # Zhang 2025 Table 2 (Dao 2022): Additive = 0.04 mg/L
  })

  model({
    # ----- Unit reparameterisation -----
    # The canonical PNA column carries months; the source's covariate is
    # centred on 21 DAYS. Recover days so the printed constants can be used
    # verbatim. See covariateData$PNA.
    pna_days <- PNA * 30.4375  # canonical months -> days (1 month = 30.4375 days)

    # ----- Individual PK parameters -----
    cl <- exp(lcl + etalcl) *
          (WT / 1.16)^e_wt_cl *
          (1 + e_pna_cl * (pna_days - 21) / 21) *
          (1 + e_ga_cl  * (GA - 26.9) / 26.9) *
          (46.6 / CREAT)^e_creat_cl
    vc <- exp(lvc) * (WT / 1.16)^e_wt_vc

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
