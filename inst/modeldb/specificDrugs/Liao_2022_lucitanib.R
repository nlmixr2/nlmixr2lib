Liao_2022_lucitanib <- function() {
  description <- paste0(
    "Two-compartment population PK model for oral lucitanib (a VEGFR1-3 / ",
    "FGFR1-3 / PDGFRalpha/beta tyrosine kinase inhibitor) in 403 adults with ",
    "advanced cancers pooled from five phase 1/2 studies (Liao 2022). Absorption ",
    "is sequential: a zero-order release of duration D1 into the depot followed ",
    "by first-order absorption (ka) into the central compartment; elimination is ",
    "linear from central. D1 is 0.814 h for the hard gelatin capsule and 0.299-",
    "fold shorter (0.243 h) for the film-coated tablet. CL/F and Q/F scale with ",
    "(WT/70)^0.75 and Vc/F and Vp/F with (WT/70)^1 (exponents fixed). IIV on ",
    "CL/F and Vc/F (correlated), D1 and Vp/F; proportional residual error."
  )
  reference <- paste(
    "Liao M, Zhou J, Wride K, Lepley D, Cameron T, Sale M, Xiao J. (2022).",
    "Population Pharmacokinetic Modeling of Lucitanib in Patients with",
    "Advanced Cancer. European Journal of Drug Metabolism and",
    "Pharmacokinetics 47(5):711-723. doi:10.1007/s13318-022-00773-w.",
    sep = " "
  )
  vignette <- "Liao_2022_lucitanib"

  # Doses are in mg and volumes in L, so central / vc is mg/L; the 1000 factor
  # in model() reports Cc in ng/mL, the unit of the assay (LLOQ 2.00 ng/mL) and
  # of the control stream's scaling S2 = V2/1000.
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Fixed allometric scaling referenced to 70 kg: exponent 0.75 on CL/F and",
        "Q/F, 1 on Vc/F and Vp/F (Liao 2022 Section 2.2, Table 4 and the",
        "Supplemental File control stream). All covariates, weight included,",
        "were coded as time-varying in the analysis dataset (Section 2.4.2).",
        "Cohort median 67.5 kg, range 35.9-159 kg (Table 3)."
      ),
      source_name = "WT"
    ),
    FORM_TABLET = list(
      description = "Film-coated tablet formulation indicator (per dose record)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (hard gelatin capsule)",
      notes = paste(
        "NOTE the comparator is a CAPSULE, not the non-tablet oral liquid named as",
        "this canonical's default reference category (as for Lalovic 2020",
        "lemborexant, Zhu 2018 asunaprevir and Wada 2023 sparsentan). The",
        "control stream codes TVD1 = THETA(8) * THETA(10)**FORM with THETA(8) the",
        "capsule D1 (0.814 h) and THETA(10) the 'tablet to capsule ratio' 0.299",
        "(Table 4), so FORM = 1 is the tablet (0.243 h, as printed in the",
        "Abstract and Figure 1). A formulation effect on relative bioavailability",
        "F1 was tested and not retained (Supplemental Table 1 model 2); the",
        "control stream carries it as THETA(9) = 1 FIX. Patients in the",
        "CO-3810-025, FINESSE and E3810-II-02 studies received capsules and/or",
        "tablets; E-3810-I-01 and INES used capsules only (Table 1)."
      ),
      source_name = "FORM"
    )
  )

  covariatesDataExcluded <- list(
    CONMED_PGP_INH = list(
      description = "Concomitant P-glycoprotein inhibitor",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Significant on CL/F in forward addition (Supplemental Table 1 model 5)",
        "but removed at backward elimination (model 9, +8 < 10.828); carried as",
        "THETA(11) = 1 FIX in the final control stream (source column IPGP)."
      )
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/dL",
      type = "continuous",
      notes = paste(
        "Tested on Vc/F as a dichotomous ALB < 3.4 indicator (control stream",
        "ALBI); significant in forward addition (model 8) but removed at backward",
        "elimination (model 10, +9); carried as THETA(12) = 1 FIX. Table 3 prints",
        "the ALB unit as 'g/mL' (median 4.00), which is g/dL."
      )
    ),
    CONMED_PPI = list(
      description = "Concomitant proton pump inhibitor",
      units = "(binary)",
      type = "binary",
      notes = paste(
        "Effect on F not significant (model 4, -2); effect on ka statistically",
        "significant (model 13, -13.4) but judged not clinically meaningful from",
        "the VPCs (Supplemental Figure S5) and not included in the final model."
      )
    ),
    CONMED_CYP2C8_INH = list(
      description = "Concomitant CYP2C8 inhibitor",
      units = "(binary)",
      type = "binary",
      notes = "Effect on CL/F not significant (Supplemental Table 1 model 6, -0.6)."
    ),
    CONMED_CYP3A4_INH = list(
      description = "Concomitant CYP3A4 inhibitor",
      units = "(binary)",
      type = "binary",
      notes = "Visually inspected only, owing to limited data (Table 2; 11 patients)."
    ),
    RENAL_GROUP = list(
      description = "FDA renal function group (normal / mild / moderate)",
      units = "(categorical)",
      type = "categorical",
      notes = "Effect on CL/F not significant (Supplemental Table 1 model 7, -0.1)."
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Visually inspected; no significant trend (Table 2)."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      notes = "Visually inspected; no significant trend (Table 2)."
    )
  )

  compartmentData <- list(
    depot = list(
      analyte = "lucitanib",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "lucitanib",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "lucitanib",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 403,
    n_studies = 5,
    n_observations = 3540,
    age_range = "26-82 years",
    age_median = "55 years",
    weight_range = "35.9-159 kg",
    weight_median = "67.5 kg",
    sex_female_pct = 82.4,
    race_ethnicity = c(
      White = 84.1,
      Black = 4.0,
      Asian = 2.5,
      Other = 9.4
    ),
    disease_state = paste(
      "Advanced or metastatic solid tumours: mostly FGFR1/FGF-amplified or",
      "non-amplified ER+/HER2- metastatic breast cancer (CO-3810-025, FINESSE,",
      "INES), plus advanced solid tumours in the first-in-human study",
      "(E-3810-I-01) and advanced lung cancer (E3810-II-02)"
    ),
    dose_range = paste(
      "5-30 mg orally once daily continuously, or 15 mg once daily 5 days",
      "on / 2 days off or 21 days on / 7 days off, in 28-day cycles; INES",
      "combined 10 or 12.5 mg with fulvestrant 500 mg IM"
    ),
    formulation = "Film-coated tablet and/or hard gelatin capsule (immediate release)",
    renal_function = paste(
      "Creatinine clearance median 95.1 mL/min (range 34.5-268); renal function",
      "group normal:mild:moderate:missing 220:141:40:2"
    ),
    hepatic_function = "NCI-ODWG hepatic group normal:mild:moderate:missing 318:80:1:4",
    co_medication = "PPI 152 patients; P-gp inhibitor 21; CYP3A4 inhibitor 12; CYP2C8 inhibitor 14 (Table 3)",
    notes = paste(
      "Demographics from Liao 2022 Table 3 (all studies). Race percentages are",
      "computed from the Table 3 counts (White:Black:Asian:other",
      "339:16:10:38); sex male:female 71:332. 3540 PK records; two BLQ samples",
      "excluded (Section 3.1)."
    )
  )

  ini({
    # --- Structural parameters (Liao 2022 Table 4, final model) -----------------
    # Reference subject: 70 kg, capsule. All disposition parameters are
    # apparent (/F); the control stream fixes F1 = 1 for both formulations.
    lcl <- log(1.90); label("Apparent clearance CL/F at 70 kg (L/h)") # Table 4 'CL/F (L/h/70 kg)' = 1.90
    lvc <- log(63.8); label("Apparent central volume Vc/F at 70 kg (L)") # Table 4 'Vc/F (L/70 kg)' = 63.8
    lka <- log(4.86); label("First-order absorption rate constant ka (1/h)") # Table 4 'Ka (1/h)' = 4.86
    lq <- log(6.23); label("Apparent intercompartmental clearance Q/F at 70 kg (L/h)") # Table 4 'Q/F (L/h/70 kg)' = 6.23
    lvp <- log(69.1); label("Apparent peripheral volume Vp/F at 70 kg (L)") # Table 4 'Vp/F (L/70 kg)' = 69.1
    ld1 <- log(0.814); label("Duration of the zero-order release into the depot D1, capsule (h)") # Table 4 'Duration of constant release into depot for capsule, D1 (h)' = 0.814

    # --- Covariate effects -----------------------------------------------------
    e_form_tablet_d1 <- 0.299; label("Multiplicative factor on D1 for the tablet versus the capsule (unitless)") # Table 4 'Effect of formulation on D1 (tablet to capsule ratio)' = 0.299; 0.814 * 0.299 = 0.243 h, the tablet D1 of the Abstract and Figure 1
    e_wt_cl <- fixed(0.75); label("Allometric exponent of body weight on CL/F and Q/F (unitless)") # Table 4 'Fixed exponent of 0.75'; Section 2.2
    e_wt_vc <- fixed(1); label("Allometric exponent of body weight on Vc/F and Vp/F (unitless)") # Table 4 'Fixed exponent of 1'; Section 2.2

    # --- Inter-individual variability (Table 4) --------------------------------
    # Table 4 reports BSV as 'log-proportional, %' (a CV on the exponential
    # eta). Variances are omega^2 = log(1 + CV^2); the CL/F - Vc/F covariance
    # is 0.596 * sqrt(omega^2_CL * omega^2_Vc). No IIV on ka or Q/F (Section 3.2).
    etalcl + etalvc ~ c(
      0.205781,
      0.146895, 0.295201
    ) # Table 4 BSV CL/F 47.8% -> log(1 + 0.478^2); Vc/F 58.6% -> log(1 + 0.586^2); correlation 0.596
    etald1 ~ 0.308367 # Table 4 BSV D1 60.1% -> log(1 + 0.601^2)
    etalvp ~ 0.462657 # Table 4 BSV Vp/F 76.7% -> log(1 + 0.767^2)

    # --- Residual error (Table 4) ----------------------------------------------
    # The control stream codes W = SQRT(ADD**2 + PROP**2 * IPRED**2) with the
    # additive SD THETA(7) = 0 FIX, i.e. a proportional-only error model.
    propSd <- 0.336; label("Proportional residual SD (fraction)") # Table 4 'Residual error (proportional, %)' = 33.6%
  })

  model({
    # 1. Individual parameters (Supplemental File control stream $PK) -----------
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    ka <- exp(lka)
    q <- exp(lq) * (WT / 70)^e_wt_cl
    vp <- exp(lvp + etalvp) * (WT / 70)^e_wt_vc
    d1 <- exp(ld1 + etald1) * e_form_tablet_d1^FORM_TABLET

    # 2. Micro-constants ---------------------------------------------------------
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 3. ODE system (NONMEM ADVAN4 TRANS4: depot CMT 1, central CMT 2,
    #    peripheral CMT 3) ----------------------------------------------------------
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # 4. Absorption input ----------------------------------------------------------
    # Zero-order release of duration D1 into the depot, then first-order
    # absorption at ka. Dose records MUST carry rate = -2 for rxode2 to honour
    # dur(depot); without it the dose is a bolus and D1 is ignored.
    dur(depot) <- d1

    # 5. Observation and error -----------------------------------------------------
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
