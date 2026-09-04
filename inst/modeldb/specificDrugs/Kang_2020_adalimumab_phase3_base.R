Kang_2020_adalimumab_phase3_base <- function() {
  description <- paste(
    "Two-compartment population PK model with sequential zero- then",
    "first-order absorption and linear elimination for subcutaneous",
    "adalimumab in 644 patients with active rheumatoid arthritis on",
    "background methotrexate, given 40 mg every 2 weeks of the biosimilar",
    "adalimumab-adbm (Cyltezo) or US-licensed Humira over 48 weeks in a",
    "Phase 3 equivalence study (NCT02137226; Kang 2020 Table 3). Fitted by",
    "full Bayesian MCMC in NONMEM with a noninformative prior on CL/F and",
    "informative priors on the remaining disposition parameters taken from",
    "the paper's Phase 1 model. A full covariate model places body weight",
    "(allometric, exponents fixed), anti-drug-antibody negativity, ADA titre,",
    "C-reactive protein, albumin and baseline rheumatoid factor on apparent",
    "clearance. IIV is a full 5x5 block on CL/F, Vc/F, Q/F, Vp/F and ka, and",
    "the residual error is combined additive plus proportional. Companion",
    "models: modellib('Kang_2020_adalimumab_phase1') and",
    "modellib('Kang_2020_adalimumab_phase3_extension').",
    sep = " "
  )
  reference <- paste(
    "Kang J, Eudy-Byrne RJ, Mondick J, Knebel W, Jayadeva G, Liesenfeld K-H.",
    "Population pharmacokinetics of adalimumab biosimilar adalimumab-adbm and",
    "reference product in healthy subjects and patients with rheumatoid",
    "arthritis to assess pharmacokinetic similarity.",
    "Br J Clin Pharmacol. 2020;86(11):2274-2285. doi:10.1111/bcp.14330.",
    "Parameter estimates from Table 3; structural and covariate model from",
    "Section 2.3, Section 3.3 and Equations 1-3 and 5.",
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
      notes              = "Baseline body weight, normalised to a 70 kg reference. Allometric exponents were fixed at 0.75 on the clearance-related parameters (CL/F, Q/F) and 1 on the volume-related parameters (Vc/F, Vp/F) rather than estimated (Kang 2020 Section 3.3 and the Table 3 row labels). Body weight was one of only two covariates whose 95% CI fell outside the +/-20% clinical-relevance band (Section 3.3 and Figure 3A).",
      source_name        = "WT"
    ),
    ADA_POS = list(
      description        = "Anti-drug-antibody positivity at the time of the PK sample",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 (ADA-positive). NOTE: this model's reference subject is ADA-POSITIVE at a titre of 16, not ADA-negative -- see notes.",
      notes              = "Time-varying. Kang 2020 Equation 5 writes the covariate as ADA-, an indicator taking the value 1 for a negative ADA test and 0 for a positive test; this file carries the canonical ADA_POS orientation and forms the paper's indicator inside model() as (1 - ADA_POS). The reference covariate set for the typical CL/F is an ADA-POSITIVE subject at a titre of 16 (Section 3.3), so exp(e_ada_neg_cl) = 0.560 is the multiplicative CL/F factor applied to ADA-negative subjects; Section 3.3 states this as clearance being approximately 44.0% lower in the absence of ADA. A typical ADA-negative 70 kg patient therefore has CL/F = 0.0244 * 0.560 = 0.0137 L/h (0.33 L/day). ADA appeared during treatment in 55% of subjects in this study (Section 3.1).",
      source_name        = "ADA-"
    ),
    ADA_TITER = list(
      description        = "Anti-drug-antibody titre, reciprocal-dilution convention",
      units              = "reciprocal dilution (powers of 2: 1, 2, 4, 8, 16, 32, ...)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying, matched in time to the PK sample. Determined by a 3-tier bridging electrochemiluminescence assay as the lowest 2-fold dilution still giving a positive response, so the values are generally powers of 2 (Kang 2020 Section 2.2). Zero-encoding convention: ADA-NEGATIVE records carry the REFERENCE titre of 16, not 0 and not 1, so that log(ADA_TITER / 16) vanishes and the entire ADA-negative effect is carried by the ADA_POS = 0 indicator. This is the only encoding consistent with the paper's own arithmetic in Section 3.3, which quotes the Table 3 factor 0.560 directly as the ADA-negative-versus-titre-16 clearance ratio. The titre effect is confirmed independently by the same paragraph: clearance is approximately 28% greater at a titre of 64 than at 16, and (64/16)^0.178 = 1.28. model() applies a defensive guard so any ADA_TITER coding on ADA-negative records gives the same result. ADA titre was the other covariate whose 95% CI fell outside the +/-20% clinical-relevance band (Figure 3A).",
      source_name        = "ADA"
    ),
    CRP = list(
      description        = "C-reactive protein, a marker of inflammatory disease activity",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying, normalised to a 3 mg/L reference (Kang 2020 Equation 5 and Section 3.3). Baseline values in this study ranged from 1 to 141 mg/L with arm means of 13.1-13.7 mg/L (Table 1). The effect is positive: apparent clearance rises with CRP. The authors judged it not clinically meaningful because the median effect size stayed within +/-20% of the typical CL/F (Section 3.3, Figure 3A).",
      source_name        = "CRP"
    ),
    ALB = list(
      description        = "Serum albumin",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying, normalised to a 43 g/L reference (Kang 2020 Equation 5 and Section 3.3). Baseline values in this study ranged from 32 to 52 g/L with arm means of 42.2-42.6 g/L (Table 1). The effect is negative: apparent clearance falls as albumin rises, which the authors attribute to higher albumin indicating more neonatal Fc receptor available to recycle IgG (Section 4.1). Judged not clinically meaningful (Section 3.3, Figure 3A). Reported in SI units in the source, matching the canonical g/L register unit; no conversion was applied.",
      source_name        = "ALB"
    ),
    RHEUMATOID_FACTOR = list(
      description        = "Baseline rheumatoid factor",
      units              = "IU/mL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-FIXED at its baseline value, unlike the time-varying CRP and ALB: Kang 2020 Equation 5 subscripts this covariate with i alone while CRP, ALB and the ADA terms carry an i,t subscript. Normalised to a 47 IU/mL reference. Baseline values in this study ranged from 5 to 6080 IU/mL with arm means of 129-182 IU/mL (Table 1); the reference of 47 is well below those means. The effect is positive: apparent clearance rises with rheumatoid factor. Judged not clinically meaningful (Section 3.3, Figure 3A). The source calls the column BRF for baseline rheumatoid factor.",
      source_name        = "BRF"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 644,
    n_studies      = 1,
    n_observations = 4342,
    age_range      = "21-80 years",
    age_median     = "arm means 52.9-53.3 years",
    weight_range   = "38.5-139 kg",
    weight_median  = "arm means 73.0-76.6 kg",
    sex_female_pct = 83,
    race_ethnicity = c(White = 95, Black = 2, Asian = 2, Other = 1),
    disease_state  = "active rheumatoid arthritis, all on background methotrexate",
    dose_range     = "40 mg subcutaneously every 2 weeks for 48 weeks",
    regions        = "not reported",
    notes          = paste(
      "Phase 3 base study NCT02137226, a randomised, double-blind,",
      "parallel-arm, multiple-dose, active-comparator equivalence trial.",
      "Subjects were randomised 1:1 to adalimumab-adbm or US-licensed Humira",
      "and re-randomised at week 24, giving three analysis arms: BI:BI",
      "(n = 323), Humira:BI (n = 146) and Humira:Humira (n = 175).",
      "Baseline demographics are Kang 2020 Table 1, phase-3-base columns; the",
      "concentration count is Section 3.1. Sampling was sparse, on days 1, 7,",
      "14, 28, 84, 168, 280, 336 and 406. Serum adalimumab was measured by a",
      "validated ELISA with limits of quantification of 25-2000 ng/mL",
      "(Section 2.2). ADA appeared during treatment in 55% of subjects.",
      sep = " "
    )
  )

  ini({
    # Structural parameters, at the reference covariate set of 70 kg body
    # weight, an ADA-POSITIVE subject with an ADA titre of 16, CRP 3 mg/L,
    # albumin 43 g/L and rheumatoid factor 47 IU/mL (Kang 2020 Section 3.3).
    # All disposition parameters are apparent (CL/F, Vc/F, Q/F, Vp/F) because
    # no intravenous data were collected.
    lcl <- log(0.0244) ; label("Apparent clearance CL/F (L/h)")                             # Table 3: CL/F (WT/70)^0.75 = 0.0244 L/h (95% CI 0.0233, 0.0257)
    lvc <- log(2.96)   ; label("Apparent central volume of distribution Vc/F (L)")          # Table 3: Vc/F (WT/70)^1 = 2.96 L (95% CI 2.64, 3.42)
    lq  <- log(0.0766) ; label("Apparent inter-compartmental clearance Q/F (L/h)")          # Table 3: Q/F (WT/70)^0.75 = 0.0766 L/h (95% CI 0.0672, 0.0872)
    lvp <- log(4.27)   ; label("Apparent peripheral volume of distribution Vp/F (L)")       # Table 3: Vp/F (WT/70)^1 = 4.27 L (95% CI 4.01, 4.53)
    lka <- log(0.0117) ; label("First-order absorption rate constant ka (1/h)")             # Table 3: ka = 0.0117 /h (95% CI 0.0108, 0.0129)
    ld1 <- log(2.82)   ; label("Duration of the zero-order input into the depot D1 (h)")    # Table 3: D1 = 2.82 h (95% CI 2.61, 3.03)

    # Bioavailability is a structural anchor rather than a paper-estimated
    # value: Kang 2020 reports every disposition parameter in its apparent
    # form and no intravenous arm was studied, so F is not identifiable and is
    # carried at 1. The paper did estimate inter-occasion variability around
    # relative F in this study (Table 3, IOV variance 0.0479, CV 22.1%); that
    # term is NOT encoded structurally here because the paper never defines
    # the occasions, and dosing ran every 2 weeks for 48 weeks. See the
    # vignette Errata.
    lfdepot <- fixed(log(1)) ; label("Subcutaneous bioavailability F (unitless)")           # Table 3 reports CL/F, Vc/F, Q/F and Vp/F throughout; F not estimated

    # Allometric exponents on body weight, held at the canonical 0.75 / 1
    # rather than estimated.
    e_wt_cl <- fixed(0.75) ; label("Allometric exponent on (WT/70) for CL/F (unitless)")    # Section 3.3 "allometrically scaled by WT with the fixed coefficients", per Section 3.2.1; Table 3 row label (WT/70)^0.75
    e_wt_q  <- fixed(0.75) ; label("Allometric exponent on (WT/70) for Q/F (unitless)")     # Section 3.3; Table 3 row label (WT/70)^0.75
    e_wt_vc <- fixed(1)    ; label("Allometric exponent on (WT/70) for Vc/F (unitless)")    # Section 3.3; Table 3 row label (WT/70)^1
    e_wt_vp <- fixed(1)    ; label("Allometric exponent on (WT/70) for Vp/F (unitless)")    # Section 3.3; Table 3 row label (WT/70)^1

    # Full covariate model on CL/F (Equation 5). The reference is an
    # ADA-POSITIVE subject at a titre of 16 with CRP 3 mg/L, ALB 43 g/L and
    # rheumatoid factor 47 IU/mL.
    e_ada_neg_cl   <- log(0.560) ; label("Effect of a negative ADA test on log CL/F, versus an ADA titre of 16 (log-scale)")  # Table 3: (ADA-) * exp(theta7) = 0.560 (95% CI 0.537, 0.580); Section 3.3 "clearances would be approximately 44.0% ... lower than that with the ADA titre value of 16"
    e_ada_titer_cl <- 0.178      ; label("Power exponent on (ADA_TITER/16) for CL/F (unitless)")                              # Table 3: (ADA/16)^theta8 = 0.178 (95% CI 0.163, 0.194); Section 3.3 cross-check (64/16)^0.178 = 1.28, the quoted 28% increase at a titre of 64
    e_crp_cl       <- 0.0867     ; label("Power exponent on (CRP/3) for CL/F (unitless)")                                     # Table 3: (CRP/3)^theta9 = 0.0867 (95% CI 0.0645, 0.108)
    e_alb_cl       <- -0.693     ; label("Power exponent on (ALB/43) for CL/F (unitless)")                                    # Table 3: (ALB/43)^theta10 = -0.693 (95% CI -1.03, -0.385)
    e_rf_cl        <- 0.0278     ; label("Power exponent on (RHEUMATOID_FACTOR/47) for CL/F (unitless)")                      # Table 3: (BRF/47)^theta11 = 0.0278 (95% CI 0.00488, 0.0432)

    # IIV: full 5x5 block on log CL/F, log Vc/F, log Q/F, log Vp/F and log ka.
    # Table 3 reports variances (IIV rows) and covariances (COV rows)
    # directly, so no CV% back-transformation is needed. Ordering below is
    # lower-triangular row-major in the block order CL, Vc, Q, Vp, ka.
    # Row by row, all from Table 3:
    #   CL: IIV CL/F = 0.122 (CV 36.0%)
    #   Vc: COV Vc/F-CL/F = -0.0139 (corr -0.0412); IIV Vc/F = 0.927 (CV 124%)
    #   Q:  COV Q/F-CL/F = 0.00861 (corr 0.0230); COV Q/F-Vc/F = -0.209
    #       (corr -0.202); IIV Q/F = 1.15 (CV 147%)
    #   Vp: COV Vp/F-CL/F = -0.0378 (corr -0.239); COV Vp/F-Vc/F = -0.233
    #       (corr -0.534); COV Vp/F-Q/F = 0.305 (corr 0.629); IIV Vp/F = 0.205
    #       (CV 47.7%)
    #   ka: COV ka-CL/F = -0.0468 (corr -0.179); COV ka-Vc/F = 0.142
    #       (corr 0.197); COV ka-Q/F = 0.560 (corr 0.698); COV ka-Vp/F = 0.140
    #       (corr 0.413); IIV ka = 0.561 (CV 86.7%)
    etalcl + etalvc + etalq + etalvp + etalka ~ c(
      0.122,
     -0.0139,   0.927,
      0.00861, -0.209,   1.15,
     -0.0378,  -0.233,   0.305,  0.205,
     -0.0468,   0.142,   0.560,  0.140,  0.561
    )                                              # Table 3, all IIV and COV rows; smallest eigenvalue 0.055, so the block is positive definite

    # Residual error: combined additive plus proportional on the linear
    # concentration scale (Equation 2). Table 3 reports VARIANCES; the SD
    # forms used here are their square roots and reproduce the CV% and SD
    # printed alongside each estimate.
    propSd <- 0.1058 ; label("Proportional residual error (fraction)")   # Table 3: proportional residual variance 0.0112 -> sqrt = 0.1058, matching the printed CV = 10.6%
    addSd  <- 0.727  ; label("Additive residual error (mg/L)")           # Table 3: additive residual variance 0.529 -> sqrt = 0.727, matching the printed SD = 0.727 mg/L
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

    # 2. Individual parameters. Equation 5 in full: the covariate terms are
    # additive on the log scale, i.e. multiplicative power functions of each
    # normalised covariate on CL/F. Allometry is on the 70 kg reference
    # subject; only CL/F carries the disease and immunogenicity covariates,
    # because the sparse Phase 3 data did not support covariates on the other
    # disposition parameters (Section 2.3).
    cl <- exp(lcl + etalcl +
                e_ada_neg_cl * (1 - ADA_POS) +
                e_ada_titer_cl * log(adaTitreRatio) +
                e_crp_cl * log(CRP / 3) +
                e_alb_cl * log(ALB / 43) +
                e_rf_cl  * log(RHEUMATOID_FACTOR / 47)) * (WT / 70)^e_wt_cl
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

    # 5. Sequential zero- then first-order absorption (Section 3.3, carried
    # over from the Phase 1 structural model): the subcutaneous dose enters
    # the depot as a zero-order input over D1 hours and then transfers to
    # central first-order at ka. Dose records must set rate = -2 so rxode2
    # takes the duration from dur(depot).
    f(depot)   <- exp(lfdepot)
    dur(depot) <- d1

    # 6. Observation. central is in mg and vc in L, so central / vc is mg/L,
    # the unit in which Kang 2020 reports the additive residual SD.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
