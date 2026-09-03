Kang_2020_adalimumab_phase1 <- function() {
  description <- paste(
    "Two-compartment population PK model with sequential zero- then",
    "first-order absorption and linear elimination for subcutaneous",
    "adalimumab in 324 healthy male subjects given a single 40 mg dose of",
    "the biosimilar adalimumab-adbm (Cyltezo), US-licensed Humira or",
    "EU-approved Humira in a 3-way bioequivalence Phase 1 study",
    "(NCT02045979; Kang 2020 Table 2). One apparent clearance was estimated",
    "across all three treatment arms. Clearance- and volume-related",
    "parameters are allometrically scaled by body weight with exponents",
    "fixed at 0.75 and 1. Apparent clearance carries two anti-drug-antibody",
    "terms: a factor for a negative ADA test and a power function of the ADA",
    "titre. IIV is a full 5x5 block on CL/F, Vc/F, Q/F, Vp/F and ka, and the",
    "residual error is combined additive plus proportional. This model",
    "supplied the informative priors used in the paper's two Phase 3 models",
    "(see modellib('Kang_2020_adalimumab_phase3_base') and",
    "modellib('Kang_2020_adalimumab_phase3_extension')).",
    sep = " "
  )
  reference <- paste(
    "Kang J, Eudy-Byrne RJ, Mondick J, Knebel W, Jayadeva G, Liesenfeld K-H.",
    "Population pharmacokinetics of adalimumab biosimilar adalimumab-adbm and",
    "reference product in healthy subjects and patients with rheumatoid",
    "arthritis to assess pharmacokinetic similarity.",
    "Br J Clin Pharmacol. 2020;86(11):2274-2285. doi:10.1111/bcp.14330.",
    "Parameter estimates from Table 2; structural and covariate model from",
    "Section 2.3, Section 3.2.1 and Equations 1-3.",
    sep = " "
  )
  vignette <- "Kang_2020_adalimumab"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot       = list(analyte = "adalimumab", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "adalimumab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "adalimumab", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Baseline body weight, normalised to a 70 kg reference. Allometric exponents were fixed at 0.75 on the clearance-related parameters (CL/F, Q/F) and 1 on the volume-related parameters (Vc/F, Vp/F) rather than estimated (Kang 2020 Section 3.2.1 and the Table 2 row labels).",
      source_name        = "WT"
    ),
    ADA_POS = list(
      description        = "Anti-drug-antibody positivity at the time of the PK sample",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 (ADA-positive). NOTE: this model's reference subject is ADA-POSITIVE at a titre of 16, not ADA-negative -- see notes.",
      notes              = "Time-varying. Kang 2020 Equation 5 writes the covariate as ADA-, an indicator taking the value 1 for a negative ADA test and 0 for a positive test; this file carries the canonical ADA_POS orientation and forms the paper's indicator inside model() as (1 - ADA_POS). The reference covariate set for the typical CL/F is an ADA-POSITIVE subject at a titre of 16 (Section 3.2.1), so exp(e_ada_neg_cl) = 0.421 is the multiplicative CL/F factor applied to ADA-negative subjects. Because 95% of subjects were ADA-negative at baseline, a typical ADA-negative 70 kg subject has CL/F = 0.0278 * 0.421 = 0.0117 L/h, which is the value comparable to published adalimumab clearances of about 0.3 L/day.",
      source_name        = "ADA-"
    ),
    ADA_TITER = list(
      description        = "Anti-drug-antibody titre, reciprocal-dilution convention",
      units              = "reciprocal dilution (powers of 2: 1, 2, 4, 8, 16, 32, ...)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying, matched in time to the PK sample. Determined by a 3-tier bridging electrochemiluminescence assay as the lowest 2-fold dilution still giving a positive response, so the values are generally powers of 2 (Kang 2020 Section 2.2). Zero-encoding convention: ADA-NEGATIVE records carry the REFERENCE titre of 16, not 0 and not 1, so that log(ADA_TITER / 16) vanishes and the entire ADA-negative effect is carried by the ADA_POS = 0 indicator. This is the only encoding consistent with the paper's own arithmetic: Table 2 reports (ADA-) * exp(theta7) = 0.421 and Section 3.2.1 states that ADA-negative subjects have 42.1% of the CL/F of subjects at a titre of 16, i.e. exactly the 0.421 factor with no additional titre contribution. model() applies a defensive guard so any ADA_TITER coding on ADA-negative records gives the same result. The titre of 16 was chosen as the reference because it was the commonly observed value and it aided convergence (Section 3.1).",
      source_name        = "ADA"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 324,
    n_studies      = 1,
    n_observations = 7255,
    age_range      = "18-55 years",
    age_median     = "mean 30.5 years (adalimumab-adbm arm) and 30.7 years (Humira arms)",
    weight_range   = "54.9-110 kg",
    weight_median  = "mean 79.4 kg (adalimumab-adbm arm) and 78.3 kg (Humira arms)",
    sex_female_pct = 0,
    race_ethnicity = c(White = 80, Black = 1, Asian = 7, `Native Hawaiian or other Pacific Islander` = 2, Other = 10),
    disease_state  = "healthy",
    dose_range     = "single 40 mg subcutaneous dose",
    regions        = "not reported",
    notes          = paste(
      "Phase 1 study NCT02045979, a randomised, double-blind, single-dose,",
      "parallel-arm, active-comparator 3-way bioequivalence study in healthy",
      "male subjects randomised 1:1:1 to adalimumab-adbm (n = 108),",
      "US-licensed Humira or EU-approved Humira (n = 216 combined).",
      "Baseline demographics are Kang 2020 Table 1, phase-1 columns; the",
      "concentration count is Section 3.1. All subjects were male. Serum",
      "adalimumab was measured by a validated ELISA with limits of",
      "quantification of 25-2000 ng/mL in neat human plasma (Section 2.2).",
      "Rich sampling on days 1-9, 14, 21, 28, 35, 44, 56 and 71.",
      sep = " "
    )
  )

  ini({
    # Structural parameters, at the reference covariate set of 70 kg body
    # weight and an ADA-POSITIVE subject with an ADA titre of 16
    # (Kang 2020 Section 3.2.1). All disposition parameters are apparent
    # (CL/F, Vc/F, Q/F, Vp/F) because no intravenous data were collected.
    lcl <- log(0.0278) ; label("Apparent clearance CL/F (L/h)")                             # Table 2: CL/F (WT/70)^0.75 = 0.0278 L/h (95% CI 0.0264, 0.0292)
    lvc <- log(2.49)   ; label("Apparent central volume of distribution Vc/F (L)")          # Table 2: Vc/F (WT/70)^1 = 2.49 L (95% CI 2.24, 2.78)
    lq  <- log(0.0708) ; label("Apparent inter-compartmental clearance Q/F (L/h)")          # Table 2: Q/F (WT/70)^0.75 = 0.0708 L/h (95% CI 0.0609, 0.0825)
    lvp <- log(3.74)   ; label("Apparent peripheral volume of distribution Vp/F (L)")       # Table 2: Vp/F (WT/70)^1 = 3.74 L (95% CI 3.53, 3.96)
    lka <- log(0.0108) ; label("First-order absorption rate constant ka (1/h)")             # Table 2: ka = 0.0108 /h (95% CI 0.00981, 0.0119)
    ld1 <- log(2.96)   ; label("Duration of the zero-order input into the depot D1 (h)")    # Table 2: D1 = 2.96 h (95% CI 2.75, 3.19)

    # Bioavailability is a structural anchor rather than a paper-estimated
    # value: Kang 2020 reports every disposition parameter in its apparent
    # form (CL/F, Vc/F, Q/F, Vp/F) and no intravenous arm was studied, so F
    # is not identifiable and is carried at 1.
    lfdepot <- fixed(log(1)) ; label("Subcutaneous bioavailability F (unitless)")           # Table 2 reports CL/F, Vc/F, Q/F and Vp/F throughout; F not estimated

    # Allometric exponents on body weight, held at the canonical 0.75 / 1
    # rather than estimated.
    e_wt_cl <- fixed(0.75) ; label("Allometric exponent on (WT/70) for CL/F (unitless)")    # Section 3.2.1: "coefficients fixed at 0.75 for clearance-related PK parameters and 1 for volume-related parameters"; Table 2 row label (WT/70)^0.75
    e_wt_q  <- fixed(0.75) ; label("Allometric exponent on (WT/70) for Q/F (unitless)")     # Section 3.2.1; Table 2 row label (WT/70)^0.75
    e_wt_vc <- fixed(1)    ; label("Allometric exponent on (WT/70) for Vc/F (unitless)")    # Section 3.2.1; Table 2 row label (WT/70)^1
    e_wt_vp <- fixed(1)    ; label("Allometric exponent on (WT/70) for Vp/F (unitless)")    # Section 3.2.1; Table 2 row label (WT/70)^1

    # Anti-drug-antibody effects on CL/F (Equation 5 terms theta7 and theta8).
    # The reference is an ADA-POSITIVE subject at a titre of 16.
    e_ada_neg_cl   <- log(0.421) ; label("Effect of a negative ADA test on log CL/F, versus an ADA titre of 16 (log-scale)")  # Table 2: (ADA-) * exp(theta7) = 0.421 (95% CI 0.399, 0.444); Section 3.2.1 "Subjects with negative ADA test had approximately 42.1% (39.9%, 44.4%) of CL/F in subjects with ADA titre value of 16"
    e_ada_titer_cl <- 0.242      ; label("Power exponent on (ADA_TITER/16) for CL/F (unitless)")                              # Table 2: (ADA/16)^theta8 = 0.242 (95% CI 0.229, 0.255)

    # IIV: full 5x5 block on log CL/F, log Vc/F, log Q/F, log Vp/F and log ka.
    # Table 2 reports variances (IIV rows) and covariances (COV rows)
    # directly, so no CV% back-transformation is needed. Ordering below is
    # lower-triangular row-major in the block order CL, Vc, Q, Vp, ka.
    # Row by row, all from Table 2:
    #   CL: IIV CL/F = 0.137 (CV 38.4%)
    #   Vc: COV Vc/F-CL/F = 0.0388 (corr 0.143); IIV Vc/F = 0.533 (CV 83.9%)
    #   Q:  COV Q/F-CL/F = 0.0306 (corr 0.0828); COV Q/F-Vc/F = -0.108
    #       (corr -0.149); IIV Q/F = 0.992 (CV 130%)
    #   Vp: COV Vp/F-CL/F = -0.0272 (corr -0.208); COV Vp/F-Vc/F = -0.0693
    #       (corr -0.270); COV Vp/F-Q/F = 0.217 (corr 0.618); IIV Vp/F = 0.124
    #       (CV 36.3%)
    #   ka: COV ka-CL/F = -0.0219 (corr -0.100); COV ka-Vc/F = 0.0861
    #       (corr 0.201); COV ka-Q/F = 0.404 (corr 0.690); COV ka-Vp/F = 0.101
    #       (corr 0.490); IIV ka = 0.345 (CV 64.2%)
    etalcl + etalvc + etalq + etalvp + etalka ~ c(
      0.137,
      0.0388,   0.533,
      0.0306,  -0.108,   0.992,
     -0.0272,  -0.0693,  0.217,  0.124,
     -0.0219,   0.0861,  0.404,  0.101,  0.345
    )                                              # Table 2, all IIV and COV rows; smallest eigenvalue 0.057, so the block is positive definite

    # Residual error: combined additive plus proportional on the linear
    # concentration scale (Equation 2). Table 2 reports VARIANCES; the SD
    # forms used here are their square roots and reproduce the CV% and SD
    # printed alongside each estimate.
    propSd <- 0.0817 ; label("Proportional residual error (fraction)")   # Table 2: proportional residual variance 0.00668 -> sqrt = 0.0817, matching the printed CV = 8.18%
    addSd  <- 0.141  ; label("Additive residual error (mg/L)")           # Table 2: additive residual variance 0.0198 -> sqrt = 0.141, matching the printed SD = 0.141 mg/L
  })

  model({
    # 1. Derived covariate terms.
    #
    # Kang 2020 Equation 5 carries the ADA effect as two additive terms on
    # log CL/F: theta7 * ADA- (an ADA-NEGATIVE indicator) and
    # theta8 * log(ADA / 16). The two terms partition the population, because
    # ADA-negative records hold the titre at its reference value of 16 in the
    # authors' data set so that the titre term vanishes there. The product
    # form below reproduces that exactly and is defensive: on an ADA-negative
    # record the ratio collapses to 1 regardless of how ADA_TITER is coded,
    # so a 0-coded titre column cannot produce log(0).
    adaTitreRatio <- ADA_POS * (ADA_TITER / 16) + (1 - ADA_POS)

    # 2. Individual parameters. Allometry is on the 70 kg reference subject.
    cl <- exp(lcl + etalcl +
                e_ada_neg_cl * (1 - ADA_POS) +
                e_ada_titer_cl * log(adaTitreRatio)) * (WT / 70)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    q  <- exp(lq  + etalq)  * (WT / 70)^e_wt_q
    vp <- exp(lvp + etalvp) * (WT / 70)^e_wt_vp
    ka <- exp(lka + etalka)
    d1 <- exp(ld1)

    # 3. Micro-constants.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 4. ODE system.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot + k21 * peripheral1 - k12 * central - kel * central
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # 5. Sequential zero- then first-order absorption (Section 3.2.1): the
    # subcutaneous dose enters the depot as a zero-order input over D1 hours
    # and then transfers to central first-order at ka. Dose records must set
    # rate = -2 so rxode2 takes the duration from dur(depot).
    f(depot)   <- exp(lfdepot)
    dur(depot) <- d1

    # 6. Observation. central is in mg and vc in L, so central / vc is mg/L,
    # the unit in which Kang 2020 reports the additive residual SD.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
