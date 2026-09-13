Takechi_2025_nemolizumab <- function() {
  description <- paste0(
    "One-compartment population PK model with first-order subcutaneous ",
    "absorption and linear elimination for nemolizumab (humanized anti-",
    "interleukin-31 receptor A monoclonal antibody) in Japanese adult, ",
    "adolescent and paediatric patients with atopic dermatitis (AD) and in ",
    "patients with prurigo nodularis (PN). Takechi 2025 re-estimated an ",
    "existing adult-AD popPK model on a pooled dataset of 527 AD patients ",
    "(including 58 children aged 6-12 years) and then externally validated ",
    "the resulting model against 153 PN patients, concluding that AD and PN ",
    "share the same PK. Body weight scales apparent clearance allometrically ",
    "(exponent 0.75) and apparent volume proportionally (exponent 1), both ",
    "referenced to 70 kg; serum albumin scales CL/F through a power function ",
    "with an exponent of -1.52 referenced to 4.5 g/dL (45 g/L). ",
    "Interindividual variability is estimated on CL/F, V/F and ka, with a ",
    "CL/F-V/F covariance; the residual is log-normal with a 16.6% CV. ",
    "Companion models from the same paper: modellib('Takechi_2025_nemolizumab_ppnrs') ",
    "(population PD for weekly average Peak Pruritus NRS) and ",
    "modellib('Takechi_2025_nemolizumab_mbma_iga') (model-based meta-analysis ",
    "of IGA success rates for nemolizumab vs dupilumab)."
  )

  reference <- paste(
    "Takechi T, Shimizu J, Kabashima K, Ieiri I.",
    "Quantitative evaluation of nemolizumab pharmacokinetics and efficacy in",
    "prurigo nodularis: a population pharmacokinetics and model-based",
    "meta-analysis approach.",
    "Dermatol Ther (Heidelb). 2025;15(12):3615-3632.",
    "doi:10.1007/s13555-025-01554-4.",
    "Structural model: Results 'PopPK Analysis' paragraph 1 and the two",
    "displayed fixed-effect equations that follow Table 2.",
    "Parameter values: Table 2, 'PopPK' block, 'Original data / Estimate'",
    "column.",
    sep = " "
  )

  vignette <- "Takechi_2025_nemolizumab"

  units <- list(
    time          = "day",
    dosing        = "mg",
    concentration = "ng/mL"
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Allometric on both disposition parameters, referenced to 70 kg: ",
        "exponent 0.75 on CL/F and 1 (implicit, printed without an exponent) ",
        "on V/F. Neither exponent carries an %RSE in Table 2, so both are ",
        "encoded as fixed(). The pooled analysis population spans 16-151 kg ",
        "(Table 1, 'All' column), which is what makes the paediatric 6-12 ",
        "year AD cohort (median 25.3-32.9 kg) informative for the allometry."
      ),
      source_name        = "body weight (Takechi 2025 Results 'PopPK Analysis', CL/F and V/F equations)"
    ),
    ALB = list(
      description        = "Serum albumin concentration.",
      units              = "g/L (SI canonical); the source calibrated the effect on g/dL and this model converts inline",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Power effect on CL/F with exponent -1.52 referenced to 4.5 g/dL. ",
        "Per the covariate register, ALB is supplied in SI g/L, so model() ",
        "applies the documented inline conversion alb_gdL <- ALB * 0.1 before ",
        "the power term (identical to (ALB / 45)^-1.52). ",
        "UNIT-LABEL ERRATUM IN THE SOURCE: Takechi 2025 Table 1 heads the ",
        "albumin row 'Serum albumin (g/L)' but tabulates values of 4.3-4.5 ",
        "with a range of 2.0-5.2. Those are g/dL, not g/L (45 g/L albumin is ",
        "normal; 4.5 g/L is not survivable), and the model equation's own ",
        "reference of 4.5 confirms the g/dL reading. The Table 1 header is a ",
        "publication unit-label error; the equation is correct."
      ),
      source_name        = "ALB (Takechi 2025 Results 'PopPK Analysis', CL/F equation; Table 2 'Covariate effect of ALB')"
    )
  )

  compartmentData <- list(
    depot   = list(analyte = "nemolizumab", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "nemolizumab", units = "mg", specimen = "serum", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 680L,
    n_studies      = 7L,
    n_observations = 4389L,
    age_range      = "6-84 years",
    age_median     = "37 years",
    weight_range   = "16-151 kg",
    weight_median  = "62.8 kg",
    sex_female_pct = 43.2,
    race_ethnicity = c(
      White = 18.4,
      Black = 2.6,
      Asian = 78.7,
      Other = 0.3
    ),
    disease_state  = paste0(
      "Moderate-to-severe atopic dermatitis with moderate-to-severe pruritus ",
      "(six studies) and prurigo nodularis (one study). PN entry required a ",
      "diagnosis of more than 6 months, limb lesions, at least 20 bilateral ",
      "prurigo nodules, and inadequate response to high-potency topical ",
      "corticosteroids and oral antihistamines."
    ),
    dose_range     = paste0(
      "Subcutaneous. Phase 1 single doses 0.003-3 mg/kg; phase 2a 0.1, 0.5 or ",
      "2 mg/kg Q4W and 2 mg/kg Q8W; phase 3 flat doses of 60 mg Q4W, 30 mg ",
      "Q4W, and 30 mg Q4W with a 60 mg loading dose (Supplementary Table S1). ",
      "Only the initial treatment phase of each study contributed data."
    ),
    regions        = paste0(
      "Japan for six of the seven studies; the phase 2a study CIM003JG ",
      "(NCT01986933) was multinational and supplies all of the non-Asian ",
      "subjects in the pooled dataset."
    ),
    notes          = paste0(
      "Demographics are Takechi 2025 Table 1, 'All' column (the PopPK ",
      "analysis dataset), which reports medians and ranges. The race ",
      "percentages above are derived by summing the per-study White/Black/",
      "Asian/Other counts in Table 1, because Table 1 prints no pooled race ",
      "row; those counts total 125/18/535/2 = 680, matching the stated PopPK ",
      "subject count exactly. Supplementary Table S1 footnote a records that ",
      "placebo-arm subjects were excluded from the PopPK dataset, which is ",
      "why 680 PopPK subjects is smaller than the sum of the enrolled counts. ",
      "The 153 PN patients (905 observations) were held out as the external ",
      "validation set (Fig. 1, Fig. S5) rather than used for estimation."
    )
  )

  ini({
    # ========================================================================
    # All values are Takechi 2025 Table 2, 'PopPK' block, 'Original data /
    # Estimate' column, with the bootstrap median and 95% CI in the adjacent
    # columns used only as a transcription check.
    #
    # OMEGA SCALE. Table 2 prints the IIV rows as bare numbers with no CV%
    # label. They are VARIANCES, not SDs, and the 'Covariance for CL/F and
    # V/F' row proves it: reading 0.135 / 0.148 as variances gives a
    # correlation of 0.0896 / sqrt(0.135 * 0.148) = 0.634, which is a legal
    # correlation, whereas reading them as SDs gives 0.0896 / sqrt(0.135^2 *
    # 0.148^2) = 4.49, which is impossible. See the vignette source-trace.
    # ========================================================================

    lcl <- log(0.340)
    label("Log apparent clearance CL/F at 70 kg and 4.5 g/dL albumin (log L/day); back-transform CL/F = 0.340 L/day")
    # Table 2 CL/F = 0.340 L/day (%RSE 1.8; bootstrap median 0.339, 95% CI
    # 0.328-0.352). Reproduced in the displayed equation that follows Table 2:
    # CL/F (L/day) = 0.340 * (body weight / 70)^0.75 * (ALB / 4.5)^-1.52.

    lvc <- log(8.44)
    label("Log apparent central volume V/F at 70 kg (log L); back-transform V/F = 8.44 L")
    # Table 2 V/F = 8.44 L (%RSE 1.9; bootstrap median 8.42, 95% CI
    # 8.14-8.78). Reproduced in the displayed equation that follows Table 2:
    # V/F (L) = 8.44 * (body weight / 70).

    lka <- log(0.548)
    label("Log first-order subcutaneous absorption rate constant (log 1/day); back-transform ka = 0.548 1/day")
    # Table 2 ka = 0.548 1/day (%RSE 8.8; bootstrap median 0.547, 95% CI
    # 0.465-0.703).

    # ---- Allometric exponents -----------------------------------------------
    # Both are printed only inside the two displayed fixed-effect equations and
    # neither appears as an estimated row in Table 2, so neither has an %RSE or
    # a bootstrap CI. They are structural assumptions held constant, hence
    # fixed(). 0.75 is the theoretical allometric clearance exponent; the V/F
    # equation is printed WITHOUT an exponent, i.e. exactly proportional to
    # body weight, so the exponent is 1.
    e_wt_cl <- fixed(0.75)
    label("Allometric exponent on CL/F for body weight referenced to 70 kg (unitless)")
    # Takechi 2025 Results 'PopPK Analysis', displayed CL/F equation.

    e_wt_vc <- fixed(1)
    label("Allometric exponent on V/F for body weight referenced to 70 kg (unitless)")
    # Takechi 2025 Results 'PopPK Analysis', displayed V/F equation:
    # V/F (L) = 8.44 * (body weight / 70), i.e. an implicit exponent of 1.

    # ---- Covariate effect ---------------------------------------------------
    e_alb_cl <- -1.52
    label("Power exponent for serum albumin on CL/F referenced to 4.5 g/dL (unitless)")
    # Table 2 'Covariate effect of ALB' = -1.52 (%RSE 9.5; bootstrap median
    # -1.51, 95% CI -1.78 to -1.22). Negative: lower albumin gives higher
    # clearance, the usual direction for an IgG whose catabolism competes with
    # albumin for FcRn recycling.

    # ---- Interindividual variability ----------------------------------------
    # Table 2 IIV CL/F = 0.135 (%RSE 10.8), IIV V/F = 0.148 (%RSE 15.0),
    # Covariance for CL/F and V/F = 0.0896 (%RSE 18.8). Encoded as an
    # OMEGA block; see the omega-scale note at the top of ini().
    # Back-transformed CV%: CL/F sqrt(exp(0.135) - 1) = 38.1%,
    # V/F sqrt(exp(0.148) - 1) = 40.0%, correlation 0.634.
    etalcl + etalvc ~ c(
      0.135,
      0.0896, 0.148
    )

    etalka ~ 0.404
    # Table 2 IIV ka = 0.404 (%RSE 17.4; bootstrap median 0.406, 95% CI
    # 0.266-0.583). Variance; CV% = sqrt(exp(0.404) - 1) = 69.8%.

    # ---- Residual error -----------------------------------------------------
    expSd <- 0.164866
    label("Log-normal residual standard deviation on the log-concentration scale (unitless); equivalent to a 16.6% CV")
    # Table 2 'log normal error (CV%)' = 16.6 (%RSE 4.1; bootstrap median 16.5,
    # 95% CI 15.2-17.9). The paper reports the BACK-TRANSFORMED CV, so the
    # log-scale SD is sqrt(log(1 + 0.166^2)) = 0.164866. The 0.7% difference
    # from a naive expSd = 0.166 is immaterial but the exact inversion is used
    # so the model reproduces the printed CV% rather than approximating it.
  })

  model({
    # ---- 1. Derived covariate terms -----------------------------------------
    # The covariate register supplies ALB in SI g/L; Takechi 2025 calibrated
    # the power term against g/dL with a reference of 4.5. Convert inline per
    # the register's documented rule so the structural exponent stays aligned
    # with its original calibration.
    alb_gdL <- ALB * 0.1

    # ---- 2. Individual parameters -------------------------------------------
    # CL/F (L/day) = 0.340 * (WT / 70)^0.75 * (ALB / 4.5)^-1.52
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl * (alb_gdL / 4.5)^e_alb_cl
    # V/F (L) = 8.44 * (WT / 70)
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    ka <- exp(lka + etalka)

    # ---- 3. Micro-constants --------------------------------------------------
    kel <- cl / vc

    # ---- 4. ODE system -------------------------------------------------------
    # One compartment with first-order absorption from a subcutaneous depot and
    # linear elimination. F is not separately identifiable from a purely
    # subcutaneous dataset, so CL and V are apparent (CL/F, V/F) and no
    # bioavailability term is applied.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # ---- 5. Observation and error -------------------------------------------
    # Amounts are mg and V/F is L, so central / vc is mg/L; the source's assay
    # and every published concentration figure are in ng/mL, so scale by 1000.
    # Assay range 100-6400 ng/mL (Methods 'Pharmacokinetic and Pharmacodynamic
    # Assessments').
    Cc <- 1000 * central / vc
    Cc ~ lnorm(expSd)
  })
}
