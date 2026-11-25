Samant_2018_ribociclib <- function() {
  description <- paste(
    "Two-compartment population PK model for oral ribociclib in adult patients",
    "with advanced cancer (Samant 2018; the Novartis registration popPK model",
    "of NDA 209092, pooled from studies X1101, X2101 and X2107). Absorption is",
    "zero-order directly into the central compartment over a fixed duration D1",
    "after an absorption lag time. Disposition is linear with apparent",
    "clearance and volumes (no bioavailability term estimated). Dose and body",
    "weight enter as power covariates on CL, Q and V2 (the peripheral volume),",
    "reproducing the more-than-dose-proportional exposure of ribociclib over",
    "50-1200 mg. IIV: V1 alone, plus a full 3x3 block on CL, Q and V2. Residual",
    "error is additive on the log scale. Proton-pump-inhibitor co-medication",
    "was tested as a covariate on relative bioavailability and not retained.",
    sep = " "
  )
  reference <- paste(
    "Samant TS, Dhuria S, Lu Y, Laisney M, Yang S, Grandeury A,",
    "Mueller-Zsigmondy M, Umehara K, Huth F, Miller M, Germa C, Elmeliegy M.",
    "Ribociclib Bioavailability Is Not Affected by Gastric pH Changes or Food",
    "Intake: In Silico and Clinical Evaluations. Clin Pharmacol Ther.",
    "2018;104(2):374-383. doi:10.1002/cpt.940.",
    "Final parameter estimates from Supplemental Table S3, identical to",
    "Table 55 of the FDA NDA 209092 (Kisqali) Multi-disciplinary Review and",
    "Evaluation (2017), which cites the applicant's population PK report Table 5-7.",
    sep = " "
  )
  vignette <- "Samant_2018_ribociclib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    central = list(analyte = "ribociclib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "ribociclib", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    DOSE = list(
      description = "Administered ribociclib dose level of the current dosing record",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power covariate on CL, Q and V2: P = TVP * (DOSE / 600)^e_dose_P.",
        "Samant 2018 Methods 'PopPK analysis' states the base model carried",
        "'dose as a structural covariate to account for nonlinearity of",
        "ribociclib PK'; Table S3 prints the exponents as 'Dose~CL', 'Dose~Q'",
        "and 'Dose~V2'. The reference dose is not printed in the paper, its",
        "supplement or the FDA review; 600 mg (the approved once-daily dose and",
        "the dose of the steady-state NCA in Table 3) is used, and the",
        "steady-state AUC check in the validation vignette supports it.",
        sep = " "
      ),
      source_name = "Dose"
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power covariate on CL, Q and V2: P = TVP * (WT / 70)^e_wt_P",
        "(Table S3 'BW~CL', 'BW~Q', 'BW~V2'). The reference weight is not",
        "printed; 70 kg is used. The FDA review reports a cohort median of",
        "72 kg (range 34-159) and describes the applicant's covariate",
        "simulation as 'a change of BW from 70 kg to 50 or 100 kg'.",
        sep = " "
      ),
      source_name = "BW"
    )
  )

  covariatesDataExcluded <- list(
    CONMED_PPI = list(
      description = "Concomitant proton-pump-inhibitor use (1 = PPI dose intensity >= 0.8)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "Tested on relative bioavailability in the full model and removed:",
        "binary classification F_rel = 0.95 (95% CI 0.82-1.09), tertiary",
        "classification 0.96 (0.80-1.11) (Samant 2018 Results 'PopPK",
        "analysis'). Not part of the final model.",
        sep = " "
      )
    )
  )

  population <- list(
    species = "human",
    n_subjects = 208,
    n_studies = 3,
    age_range = "median 60 years, range 23-84 (FDA NDA 209092 review)",
    weight_range = "median 72 kg, range 34-159 (FDA NDA 209092 review)",
    sex_female_pct = round(100 * 134 / 208, 1),
    race_ethnicity = c(Caucasian = 153, Asian = 23, Other = 32),
    disease_state = "Advanced solid tumours or lymphomas (X2101, X1101) and ER+ advanced breast cancer (X2107)",
    dose_range = "50-1200 mg orally once daily (X2101 dose escalation); 600 mg in X2107 and the approved regimen",
    regions = "Multinational; X1101 enrolled Japanese patients",
    notes = paste(
      "4,731 PK samples from 208 patients pooled from LEE011X1101, LEE011X2101",
      "and LEE011X2107 (Samant 2018 Methods 'PopPK analysis'; Table S1).",
      "52 patients had concomitant PPI use by the binary classification.",
      "External qualification on sparse data from A2301 (MONALEESA-2).",
      "Age, sex and race counts are from the FDA NDA 209092 review section",
      "13.4.3.1, which describes the same analysis.",
      sep = " "
    )
  )

  ini({
    # Apparent structural parameters at DOSE = 600 mg and WT = 70 kg.
    ltlag <- log(0.405)  ; label("Absorption lag time ALAG1 (h)")                                       # Table S3 'ALAG1 (h)' = 0.405
    ld1 <- fixed(log(3.4)) ; label("Duration of zero-order input into central D1 (h)")                 # Table S3 'D1 (h)' = 3.4 Fixed
    lcl <- log(26.433)   ; label("Apparent clearance CL/F (L/h)")                                       # Table S3 'CL (L/h)' = 26.433
    lq <- log(81.5)      ; label("Apparent intercompartmental clearance Q/F (L/h)")                     # Table S3 'Q (L/h)' = 81.5
    lvc <- log(239.492)  ; label("Apparent central volume V1/F (L)")                                    # Table S3 'V1 (L)' = 239.492
    lvp <- log(850.144)  ; label("Apparent peripheral volume V2/F (L)")                                 # Table S3 'V2 (L)' = 850.144

    # Dose power exponents (reference 600 mg)
    e_dose_cl <- -0.308  ; label("Power exponent on (DOSE / 600) for CL (unitless)")                    # Table S3 'Dose~CL' = -0.308
    e_dose_q <- -0.519   ; label("Power exponent on (DOSE / 600) for Q (unitless)")                     # Table S3 'Dose~Q' = -0.519
    e_dose_vp <- -0.344  ; label("Power exponent on (DOSE / 600) for V2 (unitless)")                    # Table S3 'Dose~V2' = -0.344

    # Body-weight power exponents (reference 70 kg)
    e_wt_cl <- 0.652     ; label("Power exponent on (WT / 70) for CL (unitless)")                       # Table S3 'BW~CL' = 0.652
    e_wt_q <- 1.216      ; label("Power exponent on (WT / 70) for Q (unitless)")                        # Table S3 'BW~Q' = 1.216
    e_wt_vp <- 0.773     ; label("Power exponent on (WT / 70) for V2 (unitless)")                       # Table S3 'BW~V2' = 0.773

    # IIV: log-normal variances. Table S3 footnote: all omega^2 are variances.
    etalvc ~ 0.925                                  # Table S3 omega2 V1 = 0.925
    etalcl + etalq + etalvp ~ c(
      0.263,                                        # Table S3 omega2 CL = 0.263
      0.124, 0.259,                                 # Table S3 omega2 CL-Q = 0.124; omega2 Q = 0.259
      0.192, 0.213, 0.312                           # Table S3 omega2 CL-V2 = 0.192; omega2 Q-V2 = 0.213; omega2 V2 = 0.312
    )

    # Residual error: additive on log-transformed concentrations, sigma^2 = 0.203
    expSd <- sqrt(0.203) ; label("Additive residual SD on the log scale (unitless)")                    # Table S3 sigma2 'additive in log domain' = 0.203
  })

  model({
    tlag <- exp(ltlag)
    d1 <- exp(ld1)

    cl <- exp(lcl + etalcl) * (DOSE / 600)^e_dose_cl * (WT / 70)^e_wt_cl
    q <- exp(lq + etalq) * (DOSE / 600)^e_dose_q * (WT / 70)^e_wt_q
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp + etalvp) * (DOSE / 600)^e_dose_vp * (WT / 70)^e_wt_vp

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- k21 * peripheral1 - k12 * central - kel * central
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Delayed zero-order input straight into central (NONMEM ALAG1 + D1).
    # Dose records must use rate = -2 so rxode2 takes the duration from dur().
    alag(central) <- tlag
    dur(central) <- d1

    # central in mg, vc in L -> mg/L; x1000 -> ng/mL
    Cc <- 1000 * central / vc
    Cc ~ lnorm(expSd)
  })
}
