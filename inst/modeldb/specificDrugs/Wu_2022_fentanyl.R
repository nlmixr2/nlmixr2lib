Wu_2022_fentanyl <- function() {
  description <- paste(
    "Two-compartment intravenous population PK model for fentanyl in preterm",
    "and term newborns (164 neonates, gestational age 23.9-42.3 weeks, birth",
    "weight 0.39-4.245 kg, postnatal age 0-68 days) pooled from a Helsinki",
    "continuous-infusion study and the Dutch multicentre DINO study.",
    "Clearance is a product of two power functions, one of birth weight",
    "(prenatal maturation, centred on 1.055 kg) and one of postnatal age in",
    "days plus 0.01 (postnatal maturation, centred on 0.99 days). Central",
    "volume scales with current body weight through a bodyweight-dependent",
    "exponent (BDE) that itself falls as a power function of current weight,",
    "both centred on 1.165 kg. Intercompartmental clearance and peripheral",
    "volume carry no covariates. Log-normal IIV on clearance and central",
    "volume; separate combined additive and proportional residual errors for",
    "the two pooled datasets."
  )
  reference <- paste(
    "Wu Y, Voller S, Flint RB, Simons SHP, Allegaert K, Fellman V, Knibbe CAJ",
    "(2022). Pre- and Postnatal Maturation are Important for Fentanyl Exposure",
    "in Preterm and Term Newborns: A Pooled Population Pharmacokinetic Study.",
    "Clin Pharmacokinet 61(3):401-412. doi:10.1007/s40262-021-01076-0.",
    "Electronic supplementary material 1 (NONMEM control stream) used for the",
    "structural form, the dataset-indicator orientation and the residual-error",
    "scale."
  )
  vignette <- "Wu_2022_fentanyl"
  units <- list(time = "h", dosing = "ug", concentration = "ug/L")

  covariateData <- list(
    WT_BIRTH = list(
      description = "Body weight at birth",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-fixed. Power effect on clearance centred on the 1.055 kg",
        "combined-cohort median (Wu 2022 Table 1: 1055 g, range 390-4245 g;",
        "Table 2 equation 'CL = TVCL x (BW/1055)^theta_BW x",
        "(PNA+0.01)^theta_PNA'). The source carries birth weight in grams; the",
        "canonical kg column enters as WT_BIRTH / 1.055, which is numerically",
        "identical to BW_g / 1055."
      ),
      source_name = "BW"
    ),
    PNA = list(
      description = "Postnatal age (chronological age since birth), time-varying",
      units = "months",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "The source analyses postnatal age in DAYS (Wu 2022 Table 1: median",
        "1.1 days at the start of treatment, range 0-68 days). The canonical",
        "PNA column is in months, so model() recovers days as",
        "PNA * 30.4375 before forming the power term (pna_days + 0.01)^0.505.",
        "The 0.01-day offset is the source's own (ESM control stream,",
        "'(PNA+0.01)**THETA(10)') and keeps the term finite on the day of",
        "birth; the Results state the typical clearance refers to a neonate of",
        "PNA 0.99 days, where (0.99 + 0.01)^theta = 1. Time-varying: postnatal",
        "age advances with the simulation clock, so supply it on every record",
        "of a dense time grid."
      ),
      source_name = "PNA"
    ),
    WT = list(
      description = "Current body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-varying. Enters central volume through the bodyweight-dependent",
        "exponent V1 = TVV1 x (CW/1165)^BDE, BDE = L1 x (CW/1165)^M, centred",
        "on the 1165 g combined-cohort median at the start of treatment (Wu",
        "2022 Table 1 and Eq. 3). The source carries current weight in grams;",
        "the canonical kg column enters as WT / 1.165. In the Helsinki dataset",
        "current weight was not recorded and was set equal to birth weight",
        "(all samples fell within the first week of life); in the DINO dataset",
        "missing values were linearly interpolated between measurements (Wu",
        "2022 Methods Section 2.4)."
      ),
      source_name = "WT"
    ),
    STUDY_DINO = list(
      description = "DINO study (dataset 2) record indicator; 1 = DINO, 0 = Helsinki dataset 1",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (Helsinki dataset 1)",
      notes = paste(
        "Selects the residual-error magnitudes only. Direct copy of the ESM",
        "control stream column ASY ('Dataset identification. 0 for dataset 1",
        "and 1 for dataset 2'), used as ADD = ASY*THETA(5) + (1-ASY)*THETA(6)",
        "and PRO = ASY*THETA(7) + (1-ASY)*THETA(8). Dataset 1 (Helsinki, 66",
        "neonates, arterial samples, radioimmunoassay with LLOQ 1 ug/L) versus",
        "dataset 2 (DINO study NCT02421068, 98 preterm neonates in four Dutch",
        "NICUs, arterial and capillary scavenge samples, LC-MS/MS with LLOQ",
        "0.3 ug/L). Set STUDY_DINO = 1 to simulate observations with the",
        "LC-MS/MS-era precision of the more recent dataset."
      ),
      source_name = "ASY"
    )
  )

  # Screened in the covariate analysis (Wu 2022 Methods Section 2.4) but not
  # retained in the final model. Documented here so the paper's covariate
  # screen is preserved without declaring covariates that model() never uses.
  covariatesDataExcluded <- list(
    GA = list(
      description = "Gestational age at birth",
      units = "weeks",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Table 1: median 28.95 weeks (range 23.90-42.30). Tested with PNA as",
        "a correlated pair on clearance; birth weight plus PNA was retained",
        "instead."
      ),
      source_name = "GA"
    ),
    PAGE = list(
      description = "Postmenstrual age",
      units = "weeks",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Computed in the control stream as PMA = PNA/7 + GA (weeks). Results",
        "Section 3.1: birth weight plus PNA on clearance was superior to PMA",
        "(dOFV = -91). Reported on the neonatal weeks scale rather than the",
        "register-default months."
      ),
      source_name = "PMA"
    ),
    SEXF = list(
      description = "Biological sex indicator, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = paste(
        "Control stream SEX: 1 = male, 2 = female (SEXF = SEX - 1). 61.6%",
        "male (Table 1). Tested as an additive shift; not retained."
      ),
      source_name = "SEX"
    ),
    SGA = list(
      description = "Small for gestational age indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (appropriate for gestational age)",
      notes = paste(
        "Control stream SGA: 1 = small for gestational age, 0 = appropriate.",
        "Tested as an additive shift; not retained."
      ),
      source_name = "SGA"
    )
  )

  compartmentData <- list(
    central = list(analyte = "fentanyl", units = "ug", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "fentanyl", units = "ug", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 164L,
    n_studies = 2L,
    age_range = "Postnatal age 0-68 days at the start of treatment; gestational age at birth 23.9-42.3 weeks",
    age_median = "Postnatal age 1.1 days at the start of treatment; gestational age 28.95 weeks",
    weight_range = "Birth weight 0.39-4.245 kg; current body weight at the start of treatment 0.39-4.245 kg",
    weight_median = "Birth weight 1.055 kg; current body weight 1.165 kg",
    sex_female_pct = 38.4,
    ga_range = "23.9-42.3 weeks (71 extremely preterm < 28 weeks, 81 very-to-late preterm 28-36 weeks, 12 term 37-42 weeks)",
    disease_state = paste(
      "Preterm and term newborns in neonatal intensive care receiving",
      "intravenous fentanyl for mechanical ventilation (Helsinki dataset),",
      "(re)intubation or analgosedation (DINO dataset)."
    ),
    dose_range = paste(
      "Helsinki dataset: 10.5 ug/kg over 1 h then a median 1.5 ug/kg/h",
      "continuous infusion for a median 58 h, with additional 1.5 ug/kg",
      "boluses in 11 patients. DINO dataset: physician-chosen boluses of",
      "0.5-3 ug/kg (median 2.1 ug/kg over a median 3 min) and/or infusions of",
      "0.5-3 ug/kg/h (median 1.0 ug/kg/h for a median 23 h)."
    ),
    regions = "Finland, Netherlands",
    notes = paste(
      "673 plasma concentrations (232 Helsinki, 441 DINO) from 164 neonates;",
      "Wu 2022 Table 1. BLOQ Helsinki samples (15.6%) were replaced by half the",
      "1 ug/L LLOQ; DINO BLOQ values were reported as measured."
    )
  )

  ini({
    # Structural parameters: Wu 2022 Table 2, 'Final model estimate (RSE %)'
    # column. The ESM control stream $THETA block holds the INITIAL values
    # (CL 0.315, V1 14, Q 0.358, V2 3.1, ...) and is used only for structure.
    lcl <- log(0.31)
    label("Clearance at birth weight 1.055 kg and postnatal age 0.99 days (L/h)") # Table 2 'TVCL' = 0.31 (10%)
    lvc <- log(10.6)
    label("Central volume at current body weight 1.165 kg (L)") # Table 2 'TVV1' = 10.6 (7%)
    lq <- log(0.573)
    label("Intercompartmental clearance (L/h)") # Table 2 'Q (L/h)' = 0.573 (35%)
    lvp <- log(3.37)
    label("Peripheral volume (L)") # Table 2 'V2 (L)' = 3.37 (23%)

    # Covariate effects. The printed Eq. 2 in the article body shows exponents
    # 1.57 and 0.502; Table 2 prints 1.47 and 0.505. Table 2 is used: its
    # values reproduce every fold-change the Results and Abstract report
    # (2^1.47 = 2.77 and 3^1.47 = 5.03 for BW 2000 and 3000 g vs 1000 g;
    # (7.01/1.01)^0.505 = 2.66, (14.01/1.01)^0.505 = 3.77,
    # (21.01/1.01)^0.505 = 4.63 for PNA 7, 14 and 21 vs 1 day), and the ESM
    # $THETA initial for BW on CL is also 1.47. See the vignette.
    e_wt_birth_cl <- 1.47
    label("Power exponent of birth weight (WT_BIRTH / 1.055 kg) on clearance (unitless)") # Table 2 'theta BW' = 1.47 (6%)
    e_pna_cl <- 0.505
    label("Power exponent of (postnatal age in days + 0.01) on clearance (unitless)") # Table 2 'theta PNA' = 0.505 (10%)

    # Bodyweight-dependent exponent on V1 (Eq. 3, Table 2):
    #   V1 = TVV1 * (CW/1165)^BDE,  BDE = L1 * (CW/1165)^M
    # At the reference weight BDE = L1, so L1 is the weight exponent on V1 at
    # 1.165 kg and M governs how that exponent changes with weight.
    e_wt_vc <- 1.56
    label("BDE intercept L1: exponent of (WT / 1.165 kg) on V1 at the reference weight (unitless)") # Table 2 'L1' = 1.56 (10%)
    bde_m_vc <- -0.417
    label("BDE slope M: power of (WT / 1.165 kg) on the V1 weight exponent (unitless)") # Table 2 'M' = -0.417 (32%)

    # IIV. Table 2 prints 'On CL (%) 44.4%' and 'On V1 (%) 45.6%'. Read as
    # 100*omega: the ESM $OMEGA initials 0.2 and 0.207 give 100*sqrt() =
    # 44.7% and 45.5%, within 1% of the printed finals, whereas the
    # sqrt(exp(omega^2)-1) reading of the same initials gives 47.1% and 48.0%
    # (5-6% away). The same group's Voller 2019 midazolam model was shown to
    # use the 100*omega convention. omega^2 = 0.444^2 and 0.456^2.
    etalcl ~ 0.197136 # Table 2 'On CL (%)' = 44.4% (RSE 9%, shrinkage 18%)
    etalvc ~ 0.207936 # Table 2 'On V1 (%)' = 45.6% (RSE 13%, shrinkage 30%)

    # Residual error: W = SQRT(ADD**2 + (PRO*IPRED)**2), Y = IPRED + EPS(1)*W
    # with $SIGMA 1 FIX (ESM control stream), so the thetas are SDs:
    # additive in ug/L and proportional as a fraction (Table 2 prints the
    # proportional rows as 0.23 and 0.361 under a '(%)' heading; they are
    # fractions, 23% and 36.1%). Separate pairs per dataset.
    addSd_helsinki <- 0.246
    label("Additive residual SD, Helsinki dataset 1 (ug/L)") # Table 2 'Additive (ug/L) on dataset 1' = 0.246 (27%)
    propSd_helsinki <- 0.23
    label("Proportional residual SD, Helsinki dataset 1 (fraction)") # Table 2 'Proportional (%) on dataset 1' = 0.23 (12%)
    addSd_dino <- 0.0297
    label("Additive residual SD, DINO dataset 2 (ug/L)") # Table 2 'Additive (ug/L) on dataset 2' = 0.0297 (26%)
    propSd_dino <- 0.361
    label("Proportional residual SD, DINO dataset 2 (fraction)") # Table 2 'Proportional (%) on dataset 2' = 0.361 (8%)
  })
  model({
    # Postnatal age: canonical months -> source days
    pna_days <- PNA * 30.4375

    # Eq. 2 / Table 2: CL = TVCL * (BW/1055)^thetaBW * (PNA+0.01)^thetaPNA
    cl <- exp(lcl + etalcl) * (WT_BIRTH / 1.055)^e_wt_birth_cl * (pna_days + 0.01)^e_pna_cl

    # Eq. 3 / Table 2: V1 = TVV1 * (CW/1165)^BDE, BDE = L1 * (CW/1165)^M
    bde_vc <- e_wt_vc * (WT / 1.165)^bde_m_vc
    vc <- exp(lvc + etalvc) * (WT / 1.165)^bde_vc

    q <- exp(lq)
    vp <- exp(lvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ADVAN3 TRANS4; all doses intravenous into central
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # ug / L = ug/L
    Cc <- central / vc

    addSd <- addSd_dino * STUDY_DINO + addSd_helsinki * (1 - STUDY_DINO)
    propSd <- propSd_dino * STUDY_DINO + propSd_helsinki * (1 - STUDY_DINO)
    Cc ~ add(addSd) + prop(propSd)
  })
}
