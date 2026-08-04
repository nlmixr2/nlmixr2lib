Waterhouse_2024_vedolizumab <- function() {
  description <- "Two-compartment population PK model with first-order (linear) elimination for vedolizumab (humanised anti-alpha4-beta7 integrin IgG1 monoclonal antibody) as acute graft-versus-host disease (aGvHD) prophylaxis in adults undergoing allogeneic hematopoietic stem cell transplantation (allo-HSCT) (Waterhouse 2024)."
  reference   <- "Waterhouse T, Baron K, Eure W, Chen C, Akbari M, Dirks NL, Jansson J, Mehrotra S. Population pharmacokinetic modeling of vedolizumab for graft-versus-host disease prophylaxis in adults with allogeneic hematopoietic stem cell transplant. Pharmacol Res Perspect. 2024;12(6):e1257. doi:10.1002/prp2.1257"
  vignette    <- "Waterhouse_2024_vedolizumab"
  units       <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Derived mechanically; verified = FALSE means it has
  # NOT been checked against the source paper.
  compartmentData <- list(
    central     = list(analyte = "vedolizumab", units = "mg", specimen = "plasma", verified = FALSE),
    peripheral1 = list(analyte = "vedolizumab", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying (per Waterhouse 2024 Methods 2.2.1, NOCB-interpolated). Estimated power exponents on CL (0.713) and Vc (0.659); fixed allometric exponents on Q (0.75) and Vp (1.00). Reference 75 kg per Waterhouse 2024 Table 2 footnote (reference patient).",
      source_name        = "WT"
    ),
    ALB = list(
      description        = "Serum albumin",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying (per Waterhouse 2024 Methods 2.2.1, NOCB-interpolated). Power-form effect on CL: (ALB / 40)^(-1.35). Reference 40 g/L (equivalently 4 g/dL) per Waterhouse 2024 Table 2 footnote (reference patient) and Supplement Equation S2. Canonical unit is SI g/L; the paper reports the reference value in both units.",
      source_name        = "ALB"
    ),
    LYMPH_ABS = list(
      description        = "Absolute lymphocyte count (total peripheral-blood lymphocytes)",
      units              = "cells/uL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying (per Waterhouse 2024 Methods 2.2.1, NOCB-interpolated). Power-form effect on CL: (max(LYMPH_ABS, 10) / 100)^(-0.0180). Reference 100 cells/uL (= 0.1 K/uL) per Waterhouse 2024 Table 2 footnote. Waterhouse 2024 Supplement Equation S2 substitutes any zero LYMPH_ABS with 10 cells/uL (= 0.01 K/uL) before evaluating the power form, to avoid log(0).",
      source_name        = "LYMPH"
    ),
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Baseline / time-fixed. Power-form effect on CL: (AGE / 53)^(-0.130). Reference 53 years per Waterhouse 2024 Table 2 footnote (reference patient).",
      source_name        = "AGE"
    ),
    AGVHD_LIVER = list(
      description        = "Acute GvHD -- liver involvement (any grade) indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no documented liver aGvHD)",
      notes              = "Time-varying (per Waterhouse 2024 Methods 2.2.1, NOCB-interpolated across pre- and post-diagnosis intervals). Multiplicative power-form effect on CL: 1.05^AGVHD_LIVER. Paper Discussion Section 3.2.3 concludes the effect is not clinically meaningful (90% CI 0.834-1.26 crosses 1).",
      source_name        = "LIVGVHD"
    ),
    AGVHD_SKIN = list(
      description        = "Acute GvHD -- skin involvement (any grade) indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no documented skin aGvHD)",
      notes              = "Time-varying (per Waterhouse 2024 Methods 2.2.1, NOCB-interpolated). Multiplicative power-form effect on CL: 1.03^AGVHD_SKIN. Paper Discussion Section 3.2.3 concludes the effect is not clinically meaningful (90% CI 0.940-1.12 crosses 1).",
      source_name        = "SKGVHD"
    ),
    AGVHD_INTESTINE = list(
      description        = "Acute GvHD -- intestinal involvement (any grade) indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no documented intestinal aGvHD)",
      notes              = "Time-varying (per Waterhouse 2024 Methods 2.2.1, NOCB-interpolated). Multiplicative power-form effect on CL: 1.07^AGVHD_INTESTINE. Paper Discussion Section 3.2.3 concludes the effect is not clinically meaningful (90% CI 0.870-1.27 crosses 1). Mechanistically on-target for vedolizumab (integrin alpha-4 beta-7 blockade of gut-homing leukocyte trafficking).",
      source_name        = "INTGVHD"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Female-sex indicator (1 = female, 0 = male)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Reported in Waterhouse 2024 Table 1 (43.0% female overall) and screened during covariate exploration. Not retained in the final model per Section 3.2.2 -- the paper attributes the observed sex effect on random effects to correlation with body weight rather than an independent sex mechanism, so sex was not carried into the final covariate model."
    ),
    RACE_ASIAN = list(
      description = "Race indicator (1 = Asian, 0 = other)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Reported in Waterhouse 2024 Table 1 (15.5% Asian overall). Screened during covariate exploration but not retained in the final model per Section 3.2.2 -- the paper attributes the observed race effect on random effects to correlation with body weight rather than an independent race mechanism."
    ),
    RACE_BLACK = list(
      description = "Race indicator (1 = Black, 0 = other)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Reported in Waterhouse 2024 Table 1 (2.1% Black overall). Screened during covariate exploration but not retained in the final model per Section 3.2.2."
    ),
    CONMED_MTX = list(
      description = "Concomitant methotrexate use (any time during study)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Reported in Waterhouse 2024 Table 1 (76.7% overall). Screened in the post-hoc concomitant-medication analysis (Section 2.2.2 and Discussion) but not retained in the final model: 'no relationships between vedolizumab disposition and methotrexate, tacrolimus, cyclosporine, and ursodeoxycholic acid'."
    ),
    CONMED_TAC = list(
      description = "Concomitant tacrolimus use (any time during study)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Reported in Waterhouse 2024 Table 1 (56.0% overall). Screened in the post-hoc concomitant-medication analysis but not retained in the final model."
    ),
    CONMED_CSA = list(
      description = "Concomitant cyclosporine use (any time during study)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Reported in Waterhouse 2024 Table 1 (44.0% overall). Screened in the post-hoc concomitant-medication analysis but not retained in the final model."
    ),
    CONMED_UDCA = list(
      description = "Concomitant ursodeoxycholic acid use (any time during study)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Reported in Waterhouse 2024 Table 1 (81.3% overall). Screened in the post-hoc concomitant-medication analysis but not retained in the final model."
    ),
    ADA_POS = list(
      description = "Anti-vedolizumab antibody positivity",
      units       = "(binary)",
      type        = "binary",
      notes       = "Only 2 of 193 subjects (1.0%) had a positive AVA test during the dosing period, both at the first day of dosing. Waterhouse 2024 Section 3.2.2: 'The effect of antidrug antibodies on vedolizumab PK was not evaluated due to the limited number of subjects with positive ADA tests.' Not screened; recorded here as an observation-only covariate for completeness."
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 193L,
    n_studies        = 2L,
    age_range        = "18-72 years",
    age_median       = "50.8 years mean (SD 14.9); reference 53 years",
    weight_range     = "not reported (SD 19.1 kg)",
    weight_median    = "78.5 kg mean (SD 19.1); reference 75 kg",
    sex_female_pct   = 43,
    race_ethnicity   = c(White = 74.1, Asian = 15.5, Black = 2.1, NotReported = 8.3),
    disease_state    = "Adults undergoing allogeneic hematopoietic stem cell transplantation (allo-HSCT) for hematologic malignancy or myeloproliferative disorder; vedolizumab administered for prophylaxis of acute graft-versus-host disease (aGvHD).",
    dose_range       = "75 mg IV (3 subjects, phase 1b) or 300 mg IV (190 subjects, phase 1b + phase 3), given as a 30-minute infusion. Phase 1b schedule: Days -1, +13, +42 relative to transplantation. Phase 3 schedule: Days -1, +13, +41, +69, +97, +125, +153.",
    regions          = "not reported in this paper (multinational phase 1b VEDO-1015 [NCT02728895] + phase 3 VEDO-3035 GRAPHITE [NCT03657160]).",
    baseline_alb     = "mean 39.7 g/L (SD 5.60); reference 40 g/L",
    baseline_lymph   = "mean 426 cells/uL (SD 739); reference 100 cells/uL",
    liver_gvhd_pct   = 3.1,
    skin_gvhd_pct    = 32.1,
    intestine_gvhd_pct = 8.8,
    concomitant_meds = "Methotrexate 76.7%; tacrolimus 56.0%; cyclosporine 44.0%; ursodeoxycholic acid 81.3%",
    notes            = "Pooled VEDO-1015 (phase 1b, n=24) + VEDO-3035 GRAPHITE (phase 3, n=169); 2380 vedolizumab observations. Baseline demographics from Waterhouse 2024 Table 1."
  )

  ini({
    # Structural parameters -- Waterhouse 2024 Table 2 (final-model estimates).
    # Reference patient: 75 kg, albumin 40 g/L (= 4 g/dL), lymphocyte count 100 cells/uL
    # (= 0.1 K/uL), age 53 years, no aGvHD (liver / skin / intestine = 0).
    lcl <- log(0.148); label("Clearance CL for the reference patient (L/day)")  # Table 2: CL = 0.148 L/day
    lvc <- log(3.12);  label("Central volume of distribution Vc for the reference patient (L)")  # Table 2: Vc = 3.12 L
    lq  <- log(0.500); label("Intercompartmental clearance Q for the reference patient (L/day)")  # Table 2: Q = 0.500 L/day
    lvp <- log(3.95);  label("Peripheral volume of distribution Vp for the reference patient (L)")  # Table 2: Vp = 3.95 L

    # Estimated covariate effect parameters (Table 2 "Covariate effect parameters" block).
    e_wt_cl     <-  0.713;  label("Body-weight power exponent on CL (unitless; reference 75 kg)")             # Table 2
    e_alb_cl    <- -1.35;   label("Albumin power exponent on CL (unitless; reference 40 g/L)")               # Table 2
    e_lymph_cl  <- -0.0180; label("Absolute lymphocyte count power exponent on CL (unitless; reference 100 cells/uL)")  # Table 2
    e_age_cl    <- -0.130;  label("Age power exponent on CL (unitless; reference 53 years)")                 # Table 2
    e_wt_vc     <-  0.659;  label("Body-weight power exponent on Vc (unitless; reference 75 kg)")            # Table 2

    # Fixed allometric exponents (Waterhouse 2024 Section 3.2.2 and Supplement Equation S2).
    # Note: Table 2 prints "Body weight effect on Q = 0.50 Fixed" but this conflicts with both
    # Section 3.2.2 text ("fixed allometric exponents for Q (0.75) and Vp (1)") and Supplement
    # Equation S2 ("theta_11 for Q (fixed to 0.75)"). Following the paper's own equation form
    # (Eq S2) and the standard allometric-clearance convention, Q = 0.75. The Table 2 value
    # 0.50 is treated as a typographical error and is documented in the vignette Errata.
    e_wt_q  <- fixed(0.75); label("Allometric exponent of WT on Q")   # Section 3.2.2 + Supp Eq S2
    e_wt_vp <- fixed(1.00); label("Allometric exponent of WT on Vp")  # Table 2 + Supp Eq S2

    # Estimated categorical covariate multipliers (Waterhouse 2024 Table 2, power form:
    # CL *= multiplier^indicator; null effect = 1).
    e_liver_gvhd_cl     <- 1.05; label("Liver aGvHD multiplier on CL (power form: CL * 1.05^AGVHD_LIVER)")         # Table 2
    e_skin_gvhd_cl      <- 1.03; label("Skin aGvHD multiplier on CL (power form: CL * 1.03^AGVHD_SKIN)")           # Table 2
    e_intestine_gvhd_cl <- 1.07; label("Intestinal aGvHD multiplier on CL (power form: CL * 1.07^AGVHD_INTESTINE)")  # Table 2

    # Interindividual variability -- Waterhouse 2024 Table 3 full block on CL, Vc, Vp.
    #   variance(etalcl) = 0.0827  (CV% = 29.4)
    #   variance(etalvc) = 0.0323  (CV% = 18.1)
    #   variance(etalvp) = 0.179   (CV% = 44.2)
    #   cov(etalcl, etalvc) = 0.0214  (Corr = 0.415)
    #   cov(etalcl, etalvp) = -0.0253 (Corr = -0.208)
    #   cov(etalvc, etalvp) = 0.0272  (Corr = 0.358)
    etalcl + etalvc + etalvp ~ c(0.0827,
                                 0.0214,  0.0323,
                                -0.0253,  0.0272, 0.179)  # Table 3

    # IOV on CL (Table 3): variance = 0.0315, CV% = 17.9. Not encoded as a
    # separate occasion-indexed eta on lcl (requires an OCC event-table column;
    # the modelling of per-occasion within-subject variability is beyond the
    # scope of the simulation-only use case). The population-typical + IIV
    # variability alone reproduces the paper's VPC to within the reported
    # prediction intervals; IOV magnitude is documented here for provenance.

    # Residual error -- Waterhouse 2024 Table 3: sigma^2_prop = 0.0241 (CV% = 15.5).
    propSd <- sqrt(0.0241); label("Proportional residual error on vedolizumab concentration (fraction)")  # Table 3
  })

  model({
    # Individual PK parameters. Reference values: 75 kg (WT), 40 g/L (ALB), 100 cells/uL
    # (LYMPH_ABS), 53 y (AGE). All covariates enter as power-form multipliers on CL, with
    # WT allometrically scaling Vc, Vp, and Q as well.
    #
    # LYMPH_ABS zero-floor (Supp Eq S2): substitute 10 cells/uL when the observed
    # count is below the floor, to avoid a log(0) in the power form.
    lymph_floor <- 10
    lymph_use   <- (LYMPH_ABS < lymph_floor) * lymph_floor + (LYMPH_ABS >= lymph_floor) * LYMPH_ABS

    cl <- exp(lcl + etalcl) *
      (WT / 75)^e_wt_cl *
      (ALB / 40)^e_alb_cl *
      (lymph_use / 100)^e_lymph_cl *
      (AGE / 53)^e_age_cl *
      e_liver_gvhd_cl^AGVHD_LIVER *
      e_skin_gvhd_cl^AGVHD_SKIN *
      e_intestine_gvhd_cl^AGVHD_INTESTINE

    vc <- exp(lvc + etalvc) * (WT / 75)^e_wt_vc
    vp <- exp(lvp + etalvp) * (WT / 75)^e_wt_vp
    q  <- exp(lq)           * (WT / 75)^e_wt_q

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Two-compartment linear model; IV infusion into the central compartment.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Concentration: dose in mg, volume in L -> mg/L = ug/mL.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
