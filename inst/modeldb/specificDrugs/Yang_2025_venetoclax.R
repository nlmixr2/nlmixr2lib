Yang_2025_venetoclax <- function() {
  description <- paste(
    "Two-compartment population PK model of venetoclax with first-order absorption,",
    "an absorption lag time and first-order elimination, in adults with hematologic",
    "malignancies receiving concomitant voriconazole (a strong CYP3A inhibitor).",
    "Every patient in the analysis dataset was co-administered voriconazole, so the",
    "CYP3A drug-drug interaction is baked into the typical values rather than carried",
    "as a covariate: apparent oral clearance is 1.31 L/h, roughly an order of magnitude",
    "below the 15.0-19.5 L/h reported for venetoclax alone and below the 2.2-3.6 L/h",
    "reported by PopPK models that treat strong CYP3A inhibition as a covariate.",
    "Serum albumin is the only retained covariate and acts on CL/F as a power function",
    "normalized to the cohort median of 38.4 g/L; higher albumin lowers apparent",
    "clearance of total (bound + unbound) venetoclax. With a terminal half-life near",
    "70 h, once-daily dosing does not reach steady state within two weeks. Venetoclax",
    "concentrations are in ug/mL.",
    sep = " "
  )
  reference <- paste(
    "Yang J, Wang H, Liu D, Cao W, Xing H, Wang P.",
    "A Population Pharmacokinetics Study of Venetoclax Concomitant with",
    "Voriconazole in Patients with Hematologic Malignancies.",
    "Drug Des Devel Ther. 2025;19:3681-3690. doi:10.2147/DDDT.S514173.",
    sep = " "
  )
  vignette <- "Yang_2025_venetoclax"
  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "ug/mL"
  )

  covariateData <- list(
    ALB = list(
      description        = "Serum albumin concentration.",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "The only covariate retained in the final model (Results 'PopPK Analysis';",
        "dOFV = 11.4, P < 0.001 on CL/F). Entered as a power function normalized to",
        "the development-cohort median albumin, CL/F = tvCL/F * (ALB / 38.4)^theta",
        "with theta = -1.49, so apparent clearance of total venetoclax FALLS as albumin",
        "rises. The paper attributes the direction to protein binding: venetoclax is",
        ">99% plasma-protein bound (fu < 0.01), so a lower albumin raises the free",
        "fraction available for clearance and hence the apparent clearance of the",
        "measured total concentration (Discussion). Reported in g/L in the source",
        "(Table 1 median 38.4, range 27.2-48.5), which is already the canonical SI unit,",
        "so no conversion is applied. Gender, age, weight, serum creatinine, creatinine",
        "clearance, serum proteins, alanine aminotransferase, aspartate aminotransferase",
        "and voriconazole concentration were all screened and none reached significance;",
        "they are recorded in covariatesDataExcluded.",
        sep = " "
      ),
      source_name        = "ALB"
    )
  )

  # Screened in the covariate analysis (Materials and Methods 'PopPK Modeling')
  # but NOT retained in the final model, and reported without a point estimate.
  # Documented here so the paper's covariate screen is preserved without
  # declaring covariates that model() never references.
  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Sex, coded female = 1 in the source (male = 0, female = 1).",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened on all PK parameters; no significant effect. Development cohort 43.3% female (Table 1)."
    ),
    AGE = list(
      description = "Age at baseline.",
      units       = "years",
      type        = "continuous",
      notes       = "Screened; no significant effect. Development cohort median 57.0 years (range 18.0-74.0, Table 1)."
    ),
    WT = list(
      description = "Body weight at baseline.",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened; no significant effect. Development cohort median 60.0 kg (range 40.0-100.0, Table 1)."
    ),
    SCR = list(
      description = "Serum creatinine at baseline.",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Screened; no significant effect. Development cohort median 58.5 umol/L (range 35.0-99.0, Table 1)."
    ),
    CRCL = list(
      description = "Creatinine clearance at baseline.",
      units       = "mL/min",
      type        = "continuous",
      notes       = "Screened; no significant effect. Development cohort median 104.5 mL/min (range 54.6-250.7, Table 1)."
    ),
    ALT = list(
      description = "Alanine aminotransferase at baseline.",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened; no significant effect. Development cohort median 18.0 U/L (range 5.0-141.0, Table 1)."
    ),
    AST = list(
      description = "Aspartate aminotransferase at baseline.",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened; no significant effect. Development cohort median 17.0 U/L (range 8.0-112.0, Table 1)."
    ),
    TPROT = list(
      description = "Total serum protein at baseline.",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened; no significant effect. Development cohort median 61.3 g/L (range 41.7-80.1, Table 1)."
    ),
    CONC_VORICONAZOLE = list(
      description = "Concomitant voriconazole plasma concentration.",
      units       = "ug/mL",
      type        = "continuous",
      notes       = paste(
        "Screened as a continuous covariate and found NOT significant on any venetoclax",
        "PK parameter (Results 'PopPK Analysis'). This is a load-bearing negative result:",
        "the CYP3A inhibition is already saturated at the voriconazole exposures observed,",
        "so the interaction is carried by the typical CL/F rather than being titrated by",
        "voriconazole concentration. The model therefore applies only to patients ON",
        "voriconazole and must not be used to predict venetoclax PK without it.",
        sep = " "
      )
    )
  )

  compartmentData <- list(
    depot       = list(analyte = "venetoclax", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "venetoclax", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "venetoclax", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 30L,
    n_studies      = 1L,
    n_observations = "261 venetoclax concentrations from the 30 development-dataset patients; a further 55 trough samples from 43 separate patients formed the external validation dataset (73 patients enrolled in total).",
    age_range      = "18.0-74.0 years",
    age_median     = "57.0 years (Table 1)",
    weight_range   = "40.0-100.0 kg",
    weight_median  = "60.0 kg (Table 1)",
    sex_female_pct = 43.3,
    race_ethnicity = "Not reported; single-centre Chinese cohort (Zhengzhou, Henan Province), so presumed predominantly Han Chinese.",
    disease_state  = paste(
      "Hematologic malignancies: acute myeloid leukemia 73.3%, chronic myelomonocytic",
      "leukemia 6.7%, acute lymphoblastic leukemia 6.7%, chronic lymphocytic leukemia 3.3%,",
      "mantle cell lymphoma 3.3%, myelodysplastic syndrome 3.3%, mixed phenotype acute",
      "leukemia 3.3% (Table 1). ECOG performance status 1 in 83.3%, 2 in 10.0%, 3 in 6.7%.",
      sep = " "
    ),
    dose_range     = paste(
      "Venetoclax 100 mg orally once daily without ramp-up in all but one patient, who",
      "reduced to 50 mg/day on day 3 for nausea and neutropenia. Simulations in the paper",
      "additionally explore 50 and 75 mg once daily.",
      sep = " "
    ),
    co_medication  = paste(
      "ALL patients received concomitant voriconazole, either continued at 200 mg twice",
      "daily or initiated at 400 mg twice daily on day 1 followed by 200 mg twice daily.",
      "Other concomitant drugs: azacitidine 60.0%, decitabine 30.0%, zanubrutinib 3.3%,",
      "none 6.7% (Table 1).",
      sep = " "
    ),
    renal_function = "Creatinine clearance median 104.5 mL/min (range 54.6-250.7); serum creatinine median 58.5 umol/L (range 35.0-99.0) (Table 1).",
    hepatic_function = "Alanine aminotransferase median 18.0 U/L (range 5.0-141.0); aspartate aminotransferase median 17.0 U/L (range 8.0-112.0); serum albumin median 38.4 g/L (range 27.2-48.5) (Table 1).",
    regions        = "China (single centre; First Affiliated Hospital of Zhengzhou University, Henan Province).",
    notes          = paste(
      "Prospective observational study conducted September 2022 to May 2023 (ethics",
      "approval KY-2022-0388). Development dataset used intensive sampling at 0-1 h",
      "pre-dose and 2, 4, 5, 6, 7, 8, 12 and 24 h post-dose on days 5-11; the external",
      "validation dataset used sparse pre-dose troughs on days 2-22. Model fitted with",
      "FOCE ELS in Phoenix NLME v8.3. Baseline demographics from Table 1.",
      sep = " "
    )
  )

  ini({
    # ======================================================================
    # Final model, Yang 2025 Table 2 ("Final Model" columns). All disposition
    # parameters are apparent (divided by the unknown oral bioavailability F),
    # so no separate bioavailability term is estimated. Estimation was by
    # FOCE ELS in Phoenix NLME v8.3 (Materials and Methods 'PopPK Modeling').
    # ======================================================================
    lka   <- log(0.11);   label("Apparent first-order absorption rate constant ka (1/h)")            # Table 2: tvka = 0.11 1/h (SE 0.04; RSE 32.48%; 95% CI 0.04 to 0.18; bootstrap median 0.11)
    lcl   <- log(1.31);   label("Apparent oral clearance CL/F (L/h) at the median albumin of 38.4 g/L")  # Table 2: tvCL/F = 1.31 L/h (SE 0.08; RSE 6.03%; 95% CI 1.15 to 1.46; bootstrap median 1.31)
    lvc   <- log(28.02);  label("Apparent central volume of distribution V/F (L)")                   # Table 2: tvV/F = 28.02 L (SE 8.39; RSE 29.96%; 95% CI 11.48 to 44.55; bootstrap median 27.13)
    lvp   <- log(87.26);  label("Apparent peripheral volume of distribution V2/F (L)")               # Table 2: tvV2/F = 87.26 L (SE 27.69; RSE 31.73%; 95% CI 32.72 to 141.80; bootstrap median 88.60)
    lq    <- log(5.29);   label("Apparent inter-compartmental clearance Q/F (L/h)")                  # Table 2: tvQ/F = 5.29 L/h (SE 1.46; RSE 27.57%; 95% CI 2.42 to 8.16; bootstrap median 5.21)
    ltlag <- log(3.36);   label("Absorption lag time Tlag (h)")                                      # Table 2: tvTlag = 3.36 h (SE 0.01; RSE 0.39%; 95% CI 3.33 to 3.38; bootstrap median 3.48)

    # ---- Retained covariate effect ---------------------------------------
    # Results 'PopPK Analysis', displayed equation:
    #     CL/F = tvCL/F * [ALB / 38.4]^(ALB on CL/F)
    # where 38.4 g/L is the median albumin of the development cohort (Table 1).
    e_alb_cl <- -1.49; label("Serum-albumin power exponent on CL/F (unitless)")                      # Table 2: ALB on CL/F = -1.49 (SE 0.33; RSE -22.19%; 95% CI -2.15 to -0.84; bootstrap median -1.49)

    # ======================================================================
    # Inter-individual variability
    # Materials and Methods 'PopPK Modeling': "An exponential error model was
    # used to measure between-subject variability", i.e. P_i = P_tv *
    # exp(eta_i). Table 2 reports the IIV terms as "omega^2", and the table
    # footnote confirms the entries are VARIANCES ("variance of
    # inter-individual variability"), so the printed values are used directly
    # on the log-normal variance scale with no CV conversion. Bracketed values
    # are the published eta shrinkages.
    #
    # No IIV on Q/F: "Random effect of Q/F was not taken into the model
    # because of shrinkage factor > 0.5" (Results 'PopPK Analysis').
    # ======================================================================
    etalka   ~ 0.09    # Table 2: omega^2 ka   = 0.09 (SE 0.04; RSE 44.44%; 95% CI 0.01 to 0.17) [shrinkage 36.55%]; bootstrap median 0.11
    etalvc   ~ 0.26    # Table 2: omega^2 V/F  = 0.26 (SE 0.12; RSE 46.15%; 95% CI 0.02 to 0.50) [shrinkage 36.18%]; bootstrap median 0.31
    etalvp   ~ 2.35    # Table 2: omega^2 V2/F = 2.35 (SE 0.72; RSE 30.64%; 95% CI 0.94 to 3.76) [shrinkage 14.34%]; bootstrap median 2.59
    etaltlag ~ 0.09    # Table 2: omega^2 Tlag = 0.09 (SE 0.04; RSE 44.44%; 95% CI 0.01 to 0.17) [shrinkage 24.03%]; bootstrap median 0.08

    # DEVIATION FROM THE PRINTED TABLE -- decimal-point correction.
    # Table 2 prints omega^2 CL/F as "0.82 +/- 0.40, RSE 48.78%, 95% CI 0.04
    # to 1.60" [shrinkage 28.68%]. That row is used here as 0.082, i.e. the
    # whole row divided by 10. Three independent lines of evidence:
    #
    #  1. The printed row is a UNIFORM 10x shift of a self-consistent row.
    #     0.40/0.82 = 48.78% and 0.82 +/- 1.96*0.40 = 0.04 to 1.60, but the
    #     identical relations hold for 0.082 +/- 0.040 (RSE 48.78%, 95% CI
    #     0.004 to 0.160) because RSE and the CI half-width scale with the
    #     estimate. Internal arithmetic therefore cannot distinguish the two
    #     readings -- it is consistent with either -- so it neither supports
    #     nor refutes the printed magnitude.
    #  2. The bootstrap column in the same row gives a median of 0.11 with a
    #     95% CI of 0.01 to 0.21. That interval excludes 0.82 outright and
    #     contains 0.082.
    #  3. The paper's own Monte Carlo simulation is reproduced only by 0.082.
    #     For 100 mg/day with voriconazole at the median albumin, the paper
    #     reports AUC24h on day 7 of 56.8 +/- 24.4 ug*h/mL (Discussion; 43%
    #     CV). Re-simulating this model over 1000 subjects gives 57.3 +/- 24.2
    #     (42% CV) with omega^2 CL/F = 0.082, versus 63.0 +/- 51.2 (81% CV)
    #     with 0.82 -- the printed value roughly doubles the published SD.
    #     The published observed AUC spread agrees with the corrected value
    #     too: 52.4 ug*h/mL median over a 19.2-121.9 range across 30 patients
    #     (Results) implies sd(log AUC) near 0.41-0.49, not the 0.91 that
    #     omega^2 = 0.82 would give.
    #
    # Correcting the decimal also removes an oddity in the printed table: at
    # 0.082 +/- 0.040 this row rounds to "0.08 +/- 0.04" at the table's two
    # decimal places, matching the ka and Tlag rows (0.09 +/- 0.04) and
    # explaining why the ka and CL/F bootstrap cells are identical.
    # Discussed in the vignette's Assumptions and deviations section.
    etalcl   ~ 0.082   # Table 2: omega^2 CL/F printed as 0.82 (SE 0.40; RSE 48.78%; 95% CI 0.04 to 1.60) [shrinkage 28.68%]; used as 0.082 -- see the note above; bootstrap median 0.11 (95% CI 0.01 to 0.21)

    # ======================================================================
    # Residual error
    # "The residual variability was well described by the proportional
    # residual error model" (Results 'PopPK Analysis'). Phoenix NLME reports
    # the proportional error as "stdev0", a standard deviation on the
    # fractional scale, so the estimate maps directly onto propSd with no
    # sqrt() conversion.
    # ======================================================================
    propSd <- 0.13; label("Proportional residual error on venetoclax Cc (fraction)")                 # Table 2: stdev0 = 0.13 (SE 0.01; RSE 5.99%; 95% CI 0.12 to 0.15) [shrinkage 14.55%]; bootstrap median 0.13
  })

  model({
    # ---- Individual parameters -------------------------------------------
    # Results 'PopPK Analysis': CL/F = tvCL/F * [ALB / 38.4]^(ALB on CL/F).
    ka   <- exp(lka + etalka)
    cl   <- exp(lcl + etalcl) * (ALB / 38.4)^e_alb_cl
    vc   <- exp(lvc + etalvc)
    vp   <- exp(lvp + etalvp)
    q    <- exp(lq)
    tlag <- exp(ltlag + etaltlag)

    # ---- Micro-constants -------------------------------------------------
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ---- Two-compartment ODE system with first-order absorption ----------
    # Materials and Methods 'PopPK Modeling' and Results 'PopPK Analysis':
    # two-compartment linear model with first-order absorption and
    # elimination and an absorption lag time.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                                k12 * central - k21 * peripheral1

    # ---- Absorption lag --------------------------------------------------
    alag(depot) <- tlag

    # ---- Observation and error model -------------------------------------
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
