Isla_2024_fosfomycinCalcium <- function() {
  description <- paste(
    "Two-compartment population PK model with first-order absorption and an",
    "absorption lag time for oral fosfomycin calcium (Fosfocina) in 24 healthy",
    "adult women studied in a four-period randomized crossover bioavailability",
    "trial (500 mg capsule single dose, 1000 mg capsule single dose, 1000 mg",
    "oral suspension single dose, and 1000 mg capsules every 8 h for 3 days).",
    "All disposition parameters are apparent (CL/F, V1/F, Q/F, V2/F) because",
    "no intravenous arm was studied, so absolute bioavailability is not",
    "identifiable and F is absorbed into every volume and clearance term.",
    "Three covariate-parameter relationships are retained: raw Cockcroft-Gault",
    "creatinine clearance on CL/F as an exponential centered at the 108 mL/min",
    "cohort median, body weight on V1/F as a linear ratio to the 64 kg cohort",
    "median, and the oral-suspension formulation on both the absorption rate",
    "constant (1.17-fold faster than capsules) and the absorption lag time",
    "(0.84-fold, i.e. shorter, than capsules). Inter-individual variability is",
    "diagonal on CL/F, V1/F, ka and the lag time; inter-occasion variability",
    "across the four crossover periods is carried on CL/F and V1/F and is",
    "larger than the corresponding IIV for both. Residual error is combined",
    "additive plus proportional. Absorption is markedly rate-limiting",
    "(ka = 0.15 1/h against kel = CL/V1 = 0.97 1/h), so the disposition is",
    "flip-flop and the apparent terminal slope is governed by ka (Isla 2024).",
    sep = " "
  )
  reference <- paste(
    "Isla A, Alarcia-Lacalle A, Solinis MA, del Pozo-Rodriguez A, Abajo Z,",
    "Cabero M, Canut-Blasco A, Rodriguez-Gascon A.",
    "Population pharmacokinetics of oral fosfomycin calcium in healthy women.",
    "J Antimicrob Chemother. 2024;79(11):2891-2898.",
    "doi:10.1093/jac/dkae295.",
    sep = " "
  )
  vignette <- "Isla_2024_fosfomycinCalcium"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. verified = TRUE: checked against Isla 2024, whose only
  # measured matrix is plasma (Methods 'Data collection and drug assay':
  # blood samples centrifuged to plasma, fosfomycin quantified by HPLC-MS/MS).
  compartmentData <- list(
    depot       = list(analyte = "fosfomycin", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "fosfomycin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "fosfomycin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description        = "Creatinine clearance estimated by the Cockcroft-Gault equation, raw mL/min and NOT BSA-normalized",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Raw Cockcroft-Gault creatinine clearance in mL/min; Isla 2024 Table 1",
        "footnote (a) states the estimator explicitly and the table reports the",
        "value in mL/min with no BSA normalisation. Cohort mean 109.4, median",
        "108.0, SD 20.7, range 82.9-158.4 mL/min. Applied as an exponential",
        "effect centered on the 108 mL/min median:",
        "CL/F = theta_CL * exp(e_crcl_cl * (CRCL - 108)) with e_crcl_cl =",
        "0.0060 (Isla 2024 Table 2 row 'CL/F (L/h) = theta_CL * e^(theta_CLCR *",
        "(CLCR - 108))'). Note the centring constant is the cohort MEDIAN",
        "(108.0), not the mean (109.4); the Results text confirms it by quoting",
        "'the estimated CL/F value for a woman with CLCR of 108 mL/min is 23.7",
        "L/h', which is the untransformed theta_CL.",
        "Every participant had normal renal function by design (all values",
        "above 82.9 mL/min, women with renal failure excluded), so the model",
        "carries NO information about renal impairment even though CRCL is its",
        "sole clearance covariate; the paper's own Discussion says so and",
        "cautions against extrapolation. The exponential form also has no",
        "upper bound, so it extrapolates without limit above the fitted range.",
        sep = " "
      ),
      source_name        = "CLCR"
    ),
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Cohort mean 65.1, median 64.3, SD 9.6, range 51.7-94.8 kg (Isla 2024",
        "Table 1). Applied as a LINEAR ratio to a 64 kg reference on the",
        "apparent central volume only: V1/F = theta_V1 * (BW / 64) (Isla 2024",
        "Table 2 row 'V1/F (L) = theta_V1 * (BW/64)'). The exponent is",
        "structurally 1 -- the paper prints no exponent and estimated none --",
        "and is encoded here as e_wt_vc <- fixed(1) so the implicit exponent is",
        "explicit rather than hidden in the equation. The 64 is the rounded",
        "cohort median weight (64.3 kg).",
        "Body weight was NOT retained on clearance: the paper's stepwise",
        "covariate search kept CRCL on CL/F and weight on V1/F only, so there",
        "is no allometric term on CL/F or Q/F in this model. Body mass index",
        "was capped at 30 kg/m^2 by the inclusion criteria, so the fitted",
        "weight range excludes obesity and the paper explicitly cautions",
        "against extrapolating the volume term to obese women.",
        sep = " "
      ),
      source_name        = "BW"
    ),
    FORM_SYRUP = list(
      description        = "Oral suspension (Fosfocina 250 mg/5 mL) versus capsule (Fosfocina 500 mg) formulation indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (Fosfocina 500 mg hard capsule, the reference formulation carrying theta_KA1 and theta_TLAG1)",
      notes              = paste(
        "Per-dose-occasion indicator: every subject received both formulations",
        "in the randomized crossover, so the column varies within subject",
        "across study periods. 1 = the manufactured 250 mg/5 mL oral suspension",
        "(20 mL giving a 1000 mg dose); 0 = the manufactured 500 mg hard",
        "capsule. Both products are commercial Fosfocina from Laboratorios ERN",
        "S.A., so this is a two-distinct-drug-products contrast, which is what",
        "distinguishes FORM_SYRUP from FORM_SUSPENSION (the latter is reserved",
        "for the same solid drug product extemporaneously compounded into a",
        "liquid at bedside, e.g. Svensson 2018 bedaquiline).",
        "Applied as TWO multiplicative ratios, both encoded on the log scale",
        "because the paper parameterises them as ratios rather than as",
        "fractional changes (Isla 2024 Table 2):",
        "ka = exp(lka + eta) * exp(log(1.17) * FORM_SYRUP), so KA(suspension) =",
        "0.15 * 1.17 = 0.1755 ~ 0.18 1/h, matching the paper's own printed",
        "'KA (suspension) = 0.18' sub-row; and",
        "tlag = exp(ltlag + eta) * exp(log(0.84) * FORM_SYRUP), so",
        "TLAG(suspension) = 0.84 * 0.84 = 0.7056 ~ 0.70 h, matching the",
        "printed 'TLAG (suspension) = 0.70' sub-row. Note theta_TLAG1 and",
        "theta_TLAG2 coincidentally share the value 0.84 with different units",
        "(h and dimensionless respectively).",
        "Disposition (CL/F, V1/F, Q/F, V2/F) is NOT affected by formulation --",
        "the paper states this explicitly in the Discussion -- and no relative",
        "bioavailability term was estimated between the two products, so the",
        "suspension and the capsule share one apparent-clearance scale and the",
        "formulation effect is purely on absorption rate and lag.",
        sep = " "
      ),
      source_name        = "formulation type (capsule or suspension)"
    ),
    OCC = list(
      description        = "Study period / crossover occasion index (1-4) driving the inter-occasion variability",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = paste(
        "Integer 1-4 identifying which of the four randomized crossover",
        "treatment periods a record belongs to (Isla 2024 Methods 'Drug",
        "administration and dosing': 500 mg capsule single dose, 1000 mg",
        "capsule single dose, 1000 mg suspension single dose, 1000 mg capsules",
        "q8h for 3 days, each separated by a washout exceeding 1 week, with the",
        "sequence randomized per subject via Tables S1 and S2). The paper does",
        "not print the NONMEM occasion column name or its coding, so the 1-4",
        "mapping to the four treatments here is this package's convention; only",
        "the number of occasions and the fact that they are the study periods",
        "come from the source. Decomposed inside model() into oc1..oc4 binary",
        "indicators multiplexing the per-occasion IOV etas on CL/F and V1/F.",
        "Pass OCC = 1 for single-occasion data so the first IOV eta applies.",
        sep = " "
      ),
      source_name        = "study occasion"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 24L,
    n_studies      = 1L,
    n_observations = 1124L,
    age_range      = "19-49 years (mean 32, median 32, SD 9); protocol eligibility 18-55 years",
    weight_range   = "51.7-94.8 kg (mean 65.1, median 64.3, SD 9.6)",
    sex_female_pct = 100,
    disease_state  = paste(
      "Healthy adult women volunteers with no evidence of significant organic",
      "or psychiatric disease, normal clinical laboratory values, normal",
      "electrocardiogram and vital signs, and negative hepatitis B / hepatitis",
      "C / HIV serology. Non-smokers, not pregnant or breastfeeding, on no",
      "prescribed or over-the-counter medication for the preceding 14 days.",
      sep = " "
    ),
    renal_function = paste(
      "Uniformly normal to mildly augmented. Serum creatinine 0.64-0.92 mg/dL",
      "(mean 0.77); Cockcroft-Gault creatinine clearance 82.9-158.4 mL/min",
      "(mean 109.4, median 108.0). Women with renal failure were excluded, so",
      "the retained CRCL-on-CL/F covariate is fitted only over a normal-",
      "function range.",
      sep = " "
    ),
    hepatic_function = paste(
      "Normal. GOT 13-34 U/L, GPT 6-54 U/L, GGT 9-40 U/L, total plasma",
      "proteins 6.3-7.6 g/dL (Isla 2024 Table 1). Transaminases and total",
      "protein were screened as covariates and none was retained.",
      sep = " "
    ),
    dose_range     = paste(
      "Four treatments per subject in a randomized crossover with washout",
      "exceeding one week between periods: (i) 500 mg fosfomycin as one",
      "Fosfocina 500 mg capsule, single dose; (ii) 1000 mg as two Fosfocina",
      "500 mg capsules, single dose; (iii) 1000 mg as 20 mL of Fosfocina",
      "250 mg/5 mL oral suspension, single dose; (iv) 1000 mg as two Fosfocina",
      "500 mg capsules every 8 h for 3 days, with sampling after the last",
      "dose. All doses were given fasted with 200 mL of water, and water was",
      "withheld from 1 h before until 1 h after dosing.",
      sep = " "
    ),
    regions        = "Spain (single centre: Clinical Trial Unit, Araba University Hospital, Vitoria-Gasteiz).",
    notes          = paste(
      "Regulatory identifiers: AEMPS code PD7522.22, EudraCT 2020-001664-28.",
      "Thirteen plasma samples per subject per period (pre-dose and 1, 1.5, 2,",
      "2.5, 3, 3.5, 4, 4.5, 6, 8, 12 and 24 h), giving 1124 concentration-time",
      "records for model building. All 24 women completed the single-dose",
      "periods; four were excluded from the multiple-dose period (two dosing",
      "deviations, two sampling deviations), leaving n = 20 there, and one",
      "pre-dose multiple-dose sample was discarded as an analytical error.",
      "Assay: HPLC-MS/MS, linear 50 ng/mL (LLOQ) to 50000 ng/mL.",
      "Fitted in NONMEM 7.4 with FOCE-INTER; final model condition number",
      "8.08; evaluated by pcVPC (200 replicates) and a 1000-sample bootstrap",
      "of which 919 runs converged.",
      sep = " "
    )
  )

  ini({
    # ---- Structural fixed effects (Isla 2024 Table 2 'Estimate' column) ----
    # All disposition parameters are APPARENT (divided by the unknown oral
    # bioavailability F): the trial had no intravenous arm, so F is not
    # identifiable and no separate fdepot term is estimated. Reference subject
    # is a woman with CRCL = 108 mL/min and body weight 64 kg receiving the
    # capsule formulation.
    lcl   <- log(23.7)  ; label("Apparent oral clearance CL/F at CRCL 108 mL/min (L/h)")                             # Table 2 theta_CL = 23.7 L/h (RSE 5%; bootstrap median 23.1, 95% CI 20.5-25.6); Results text 'The estimated CL/F value for a woman with CLCR of 108 mL/min is 23.7 L/h'
    lvc   <- log(24.4)  ; label("Apparent central volume of distribution V1/F at 64 kg (L)")                         # Table 2 theta_V1 = 24.4 L (RSE 14%; bootstrap median 24.9, 95% CI 18.9-32.5); Results text 'the typical value for V1/F ... is 24.4 L'
    lq    <- log(4.04)  ; label("Apparent inter-compartmental clearance Q/F (L/h)")                                  # Table 2 theta_Q = 4.04 L/h (RSE 24%; bootstrap median 4.40, 95% CI 2.56-6.53); Results text 'inter-compartmental clearance (Q/F) ... 4.04 L/h'
    # The paper parameterises the peripheral volume as a RATIO to the central
    # volume rather than estimating it directly: Table 2 gives
    # 'Vss/F = (V1 * (1 + theta_V2))' with theta_V2 = 4.94, so the apparent
    # peripheral volume is V2/F = V1/F * theta_V2 = 24.4 * 4.94 = 120.536 L and
    # the steady-state volume is 24.4 * (1 + 4.94) = 144.936 L, reproducing the
    # Results text 'steady-state volume of distribution (Vss/F) ... 144.9 L'
    # and the Discussion's 'Vss/F (144.94 L)' to the digit. Because theta_V2
    # multiplies V1, the peripheral volume inherits V1's body-weight term (see
    # e_wt_vp below); it does NOT inherit V1's random effects, because every
    # equation in the Table 2 left-hand column is a typical-value relationship
    # with the etas tabulated separately underneath.
    lvp   <- log(24.4 * 4.94) ; label("Apparent peripheral volume of distribution V2/F at 64 kg (L)")                # Table 2 theta_V2 = 4.94 (RSE 38%; bootstrap median 5.13, 95% CI 3.14-16.21) as the V2/V1 ratio; V2/F = 24.4 * 4.94 = 120.536 L
    lka   <- log(0.15)  ; label("First-order absorption rate constant ka for the capsule (1/h)")                     # Table 2 theta_KA1 = 0.15 1/h (RSE 10%; bootstrap median 0.15, 95% CI 0.12-0.19)
    ltlag <- log(0.84)  ; label("Absorption lag time for the capsule (h)")                                           # Table 2 theta_TLAG1 = 0.84 h (RSE 2%; bootstrap median 0.84, 95% CI 0.81-0.86)

    # ---- Covariate effects (Isla 2024 Table 2 covariate rows) ----
    e_crcl_cl <- 0.0060 ; label("Exponential coefficient on centered Cockcroft-Gault creatinine clearance for CL/F (1/(mL/min))") # Table 2 theta_CLCR = 0.0060 (RSE 38%; bootstrap median 0.0062, 95% CI 0.0017-0.0107) in 'CL/F = theta_CL * e^(theta_CLCR * (CLCR - 108))'
    # The two weight exponents are structural, not estimated: Table 2 writes
    # 'V1/F (L) = theta_V1 * (BW/64)' with the ratio to the first power and no
    # exponent parameter anywhere in the table, so both are fixed at 1. The
    # peripheral exponent exists only because V2/F = V1/F * theta_V2 carries
    # the same (BW/64) factor through the ratio parameterisation above.
    e_wt_vc   <- fixed(1) ; label("Power exponent on (WT / 64 kg) for V1/F (unitless)")                              # Table 2 'V1/F (L) = theta_V1 * (BW/64)': ratio enters to the first power
    e_wt_vp   <- fixed(1) ; label("Power exponent on (WT / 64 kg) for V2/F (unitless)")                              # Table 2 'Vss/F = (V1 * (1 + theta_V2))': V2/F inherits V1/F's (BW/64) factor
    # Formulation effects are multiplicative RATIOS in the source, so they are
    # encoded on the log scale (the registered idiom -- see
    # Valle_2005_exemestane.R) rather than as (1 + fraction) terms.
    e_form_syrup_ka   <- log(1.17) ; label("Effect of the oral suspension on ka relative to the capsule (log scale)")    # Table 2 theta_KA2 = 1.17 (RSE 9%; bootstrap median 1.15, 95% CI 0.98-1.42) in 'KA (suspension) = theta_KA1 x theta_KA2'; 0.15 * 1.17 = 0.1755, matching the printed 'KA (suspension) = 0.18'
    e_form_syrup_tlag <- log(0.84) ; label("Effect of the oral suspension on the lag time relative to the capsule (log scale)") # Table 2 theta_TLAG2 = 0.84 (RSE 5%; bootstrap median 0.83, 95% CI 0.73-0.91) in 'TLAG (suspension) = theta_TLAG1 x theta_TLAG2'; 0.84 * 0.84 = 0.7056, matching the printed 'TLAG (suspension) = 0.70'

    # ---- Inter-individual variability (Isla 2024 Table 2 'IIV on ... (%)' rows) ----
    # Table 2 heads these rows with '(%)', and the Methods state IIV 'was
    # estimated assuming a log-normal distribution of the parameter values',
    # so the tabulated numbers are %CV on the log-normal scale. nlmixr2 takes
    # the log-scale VARIANCE, recovered as omega^2 = log(1 + CV^2):
    #   CL/F  15.7% -> log(1 + 0.157^2) = 0.0243501
    #   V1/F  29.4% -> log(1 + 0.294^2) = 0.0829026
    #   KA    32.7% -> log(1 + 0.327^2) = 0.1015895
    #   TLAG   4.8% -> log(1 + 0.048^2) = 0.0023013
    # This %CV reading is confirmed arithmetically by the paper's own Table 3:
    # the simulated cohort means there sit above the corresponding medians by
    # exactly the log-normal factor exp(sum(omega^2)/2) computed from these
    # variances plus the IOV variances below -- 1.0711 for CL/F and 1.2529 for
    # V1/F. Worked check for CL/F at CRCL 150 mL/min: median 23.7 *
    # exp(0.0060 * 42) = 30.49, times 1.0711 = 32.7 against Table 3's reported
    # mean of 32.8. For V1/F at 90 kg: median 24.4 * 90/64 = 34.31, times
    # 1.2529 = 43.0 against Table 3's reported mean of 42.7. Reading the
    # tabulated numbers as variances instead would inflate both by orders of
    # magnitude and reproduce neither column.
    etalcl   ~ 0.0243501  # Table 2 'IIV on CL/F (%)' = 15.7 (RSE 27%, eta-shrinkage 34%; bootstrap median 15.7, 95% CI 6.3-22.6) -> variance log(1 + 0.157^2)
    etalvc   ~ 0.0829026  # Table 2 'IIV on V1/F (%)' = 29.4 (RSE 47%, eta-shrinkage 43%; bootstrap median 30.8, 95% CI 10.5-51.6) -> variance log(1 + 0.294^2)
    etalka   ~ 0.1015895  # Table 2 'IIV on KA (%)'   = 32.7 (RSE 18%, eta-shrinkage 6.0%; bootstrap median 32.8, 95% CI 21.0-44.2) -> variance log(1 + 0.327^2)
    etaltlag ~ 0.0023013  # Table 2 'IIV on TLAG (%)' = 4.8 (RSE 34%, eta-shrinkage 44%; bootstrap median 5.1, 95% CI 2.1-7.5) -> variance log(1 + 0.048^2)

    # ---- Inter-occasion variability (Isla 2024 Table 2 'IOV on ... (%)' rows) ----
    # The Methods describe IOV as 'evaluated using exponential models ... to
    # assess differences in individual parameters across study occasions', and
    # the Results report that adding it reduced IIV on CL/F from 22.7% to
    # 20.4% and on V1/F from 42.5% to 32.1% and that it was retained. Same
    # %CV -> log-variance transform as the IIV block:
    #   CL/F  34.6% -> log(1 + 0.346^2) = 0.1130751
    #   V1/F  66.7% -> log(1 + 0.667^2) = 0.3680325
    # One variance per parameter is reported and it is shared by all four
    # crossover periods; nlmixr2 has no NONMEM `$OMEGA BLOCK(1) SAME`
    # shortcut, so occasion 1 carries the estimated variance and occasions
    # 2-4 fix it to the same value (the registered idiom -- see
    # Blackman_2026_methotrexate.R, Jonsson_2011_ethambutol.R).
    etaiov_cl_1 ~ 0.1130751       # Table 2 'IOV on CL/F (%)' = 34.6 (RSE 10%; bootstrap median 33.8, 95% CI 27.8-39.4) -> variance log(1 + 0.346^2) (estimated)
    etaiov_cl_2 ~ fix(0.1130751)  # SAME-equivalent: equal to the occasion-1 IOV variance
    etaiov_cl_3 ~ fix(0.1130751)  # SAME-equivalent: equal to the occasion-1 IOV variance
    etaiov_cl_4 ~ fix(0.1130751)  # SAME-equivalent: equal to the occasion-1 IOV variance
    etaiov_vc_1 ~ 0.3680325       # Table 2 'IOV on V1/F (%)' = 66.7 (RSE 14%; bootstrap median 63.9, 95% CI 43.2-83.3) -> variance log(1 + 0.667^2) (estimated)
    etaiov_vc_2 ~ fix(0.3680325)  # SAME-equivalent: equal to the occasion-1 IOV variance
    etaiov_vc_3 ~ fix(0.3680325)  # SAME-equivalent: equal to the occasion-1 IOV variance
    etaiov_vc_4 ~ fix(0.3680325)  # SAME-equivalent: equal to the occasion-1 IOV variance

    # ---- Residual error (Isla 2024 Table 2 'RE ...' rows) ----
    # Methods: 'Residual variability was explored according to additive,
    # proportional and combined (additive + proportional) error models'; the
    # final model reports one row of each, i.e. the combined model.
    # The proportional row is headed 'RE proportional (%)' but carries the
    # value 0.199, which is the NONMEM fraction (19.9%) rather than 0.199%;
    # the '(%)' in the header is a table artefact. A 0.199% proportional error
    # is irreconcilable with the 0.209 mg/L additive term on 1-12 mg/L
    # concentrations and with the epsilon-shrinkage of 8%.
    addSd  <- 0.209 ; label("Additive residual error standard deviation (mg/L)")                                     # Table 2 'RE additive (mg/L)' = 0.209 (RSE 17%, eps-shrinkage 8%; bootstrap median 0.201, 95% CI 0.128-0.309)
    propSd <- 0.199 ; label("Proportional residual error standard deviation (fraction)")                             # Table 2 'RE proportional (%)' = 0.199 (RSE 8%; bootstrap median 0.198, 95% CI 0.159-0.229) -- the value is the fraction 19.9%, see the note above
  })

  model({
    # 1. Decompose the integer crossover-period column into binary indicators
    #    to multiplex the four inter-occasion-variability etas on log-CL/F and
    #    log-V1/F. For single-occasion data pass OCC = 1.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)

    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2 + oc3 * etaiov_cl_3 + oc4 * etaiov_cl_4
    iov_vc <- oc1 * etaiov_vc_1 + oc2 * etaiov_vc_2 + oc3 * etaiov_vc_3 + oc4 * etaiov_vc_4

    # 2. Individual apparent PK parameters. Reference subject: CRCL 108 mL/min,
    #    body weight 64 kg, capsule formulation.
    #    - CL/F carries the exponential creatinine-clearance term centered at
    #      108 mL/min, its own eta, and the occasion-specific IOV eta.
    #    - V1/F carries the linear body-weight ratio, its own eta, and the
    #      occasion-specific IOV eta.
    #    - V2/F is the paper's V1 * theta_V2 ratio, so it carries the same
    #      body-weight ratio but no random effects of its own.
    #    - Q/F has neither a covariate nor a random effect in the final model.
    cl  <- exp(lcl + etalcl + iov_cl + e_crcl_cl * (CRCL - 108))
    vc  <- exp(lvc + etalvc + iov_vc) * (WT / 64)^e_wt_vc
    vp  <- exp(lvp) * (WT / 64)^e_wt_vp
    q   <- exp(lq)

    # 3. Absorption. The suspension absorbs 1.17-fold faster and starts
    #    0.84-fold sooner than the capsule; both effects are ratios, applied
    #    on the log scale.
    ka   <- exp(lka   + etalka   + e_form_syrup_ka   * FORM_SYRUP)
    tlag <- exp(ltlag + etaltlag + e_form_syrup_tlag * FORM_SYRUP)

    # 4. Micro-constants.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # 5. Two-compartment disposition with first-order absorption from an oral
    #    depot. Note ka (0.15 1/h for capsules) is far below kel
    #    (23.7 / 24.4 = 0.97 1/h), so absorption is rate-limiting and the
    #    profile is flip-flop: the apparent terminal slope reflects ka, not
    #    elimination.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central -
      k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # 6. Absorption lag on the depot input.
    alag(depot) <- tlag

    # 7. Observation. Dose in mg and vc in L give Cc in mg/L, the unit the
    #    paper reports plasma fosfomycin in throughout (Results: 'Cmax ranged
    #    from 1.1 to 5.2 mg/L with 500 mg capsules').
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
