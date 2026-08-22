Rosenborg_2025_salmeterol <- function() {
  description <- paste(
    "Empirical population PK model for inhaled salmeterol (SALM) in healthy",
    "adults given 3 x 50 ug (150 ug total) from a dry powder inhaler, fitted",
    "simultaneously across three two-way crossover bioequivalence studies",
    "(Rosenborg 2025 Model 4). Monophasic first-order absorption from a",
    "hypothetical pulmonary deposition site into a three-compartment apparent",
    "systemic disposition; the nominal inhaled dose is assumed to be deposited",
    "instantaneously and completely at the deposition site, so all clearances",
    "and volumes are apparent (CL/F, V1/F, ...) and carry an unknown absolute",
    "bioavailability < 1. The test product (Wixela Inhub) differs from the",
    "reference (Advair Diskus) in the absorption rate constant k41 and in",
    "relative bioavailability F4_rel; reference bioavailability is anchored at",
    "F = 1. Between-subject variability is carried on CL/F, Q2/F, Q3/F, V1/F,",
    "V2/F, V3/F, k41 and F4_rel, and a nested between-study random effect on",
    "F4_rel (indexed by SIDN) reproduces the paper's interstudy variability",
    "level. Companion models for the co-administered fluticasone propionate",
    "are Rosenborg_2025_fluticasone_300ug, Rosenborg_2025_fluticasone_750ug",
    "and Rosenborg_2025_fluticasone_1500ug.")
  reference <- paste(
    "Rosenborg J, Backman P, Bengtsson T, Haughie S.",
    "Relative Bioavailability of Inhaled Fluticasone Propionate and Salmeterol",
    "- is Population Pharmacokinetic Modelling a Relevant Alternative to a",
    "Non-Compartmental Approach?",
    "Drug Des Devel Ther. 2025;19:9653-9670. doi:10.2147/DDDT.S480189")
  vignette <- "Rosenborg_2025_fluticasone_salmeterol"
  units    <- list(time = "h", dosing = "ug", concentration = "ng/L")

  compartmentData <- list(
    depot       = list(analyte = "salmeterol", units = "ug", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "salmeterol", units = "ug", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "salmeterol", units = "ug", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "salmeterol", units = "ug", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    FORM_WIXELA_INHUB = list(
      description        = "Indicator that the inhalation was taken from the Wixela Inhub dry powder inhaler (test product), 1 = Wixela Inhub, 0 = Advair Diskus.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (Advair Diskus, the reference product; its relative extent of bioavailability is the anchor F = 1).",
      notes              = "Per-dose-record indicator. In this two-way crossover the same subject carries 1 on one period's dose row and 0 on the other's. The supplement's NONMEM code writes TREA = 1 for test and TREA = 2 for reference and derives TREA1 / TREA2 indicators from it; FORM_WIXELA_INHUB = TREA1. It gates two parameters (Rosenborg 2025 Figure 1): the monophasic absorption rate constant, K41 = EXP(TREA1*MU_7 + TREA2*MU_10 + ETA(7)), i.e. a separate typical value per product sharing one eta; and relative bioavailability, F4 = EXP(TREA1*(MU_11 + ETA(8) + ETA(9)) + TREA2*0), so reference doses get F = 1 exactly while test doses get F4_rel carrying both its subject-level and its study-level random effect.",
      source_name        = "TREA (supplement Sect. 1: 'Treatment alternative, test (TREA=1) and reference (TREA=2) formulation')"
    ),
    SIDN = list(
      description        = "Integer study index (1, 2 or 3) identifying which of the three crossover studies a subject's records belong to, used as the nesting level for the between-study random effect on relative bioavailability.",
      units              = "(count)",
      type               = "categorical",
      reference_category = "n/a -- SIDN is a nesting level, not a covariate with a reference category. It indexes ETA(9), declared in the supplement's NONMEM code as $LEVEL STUD=(9[1]) ; interstudy variability.",
      notes              = "The three studies are distinguished only by the fluticasone propionate strength inhaled alongside salmeterol (FP 3 x 100, 3 x 250 and 3 x 500 ug); the salmeterol dose was 3 x 50 ug in all three, which is why Model 4 pools them. SIDN maps directly onto the paper's STUD column (studies 1, 2 and 3). Because it is consumed by the ini() nesting syntax rather than by a model() expression, covariateData[['SIDN']] is deliberately not referenced anywhere in model(). Simulation caveat: an event table must carry at least two distinct SIDN values or rxSolve() fails, and omega must be passed explicitly because a nested omega is a list of matrices keyed by level.",
      source_name        = "STUD"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 198L,
    n_studies      = 3L,
    age_range      = "mean age 33.8 years (study 1), 37.7 years (study 2) and 35.7 years (study 3); Rosenborg 2025 Materials and Methods",
    sex_female_pct = 54.0,
    race_ethnicity = c(White = NA_real_, Black = NA_real_, Asian = NA_real_, Other = NA_real_),
    disease_state  = "Healthy adult volunteers.",
    dose_range     = "Three inhalations of salmeterol 50 ug (150 ug total) in every study, co-formulated with fluticasone propionate 100, 250 or 500 ug per inhalation; single dose per period, two periods with a 7-day washout.",
    regions        = "USA",
    n_observations = "7186 salmeterol plasma concentrations (Rosenborg 2025 Results).",
    notes          = "Three separate two-way crossover studies, each recruiting 66 healthy subjects (198 in total), of whom 61 to 65 per study completed both treatment periods. Panels were homogeneous, black and white, mean BMI approximately 26; study 1: 29 female / 37 male; study 2: 36 female / 30 male; study 3: 42 female / 24 male. Rosenborg 2025 states that 'neither demographic nor other covariates were considered in the evaluation of NCA-based results and therefore not in this alternative model analysis either', so product is the only fixed-effect covariate in the model. Sampling: pre-dose and 2, 5, 10, 15, 20, 30 and 45 min, and 1, 1.5, 2, 3, 4, 6, 8, 12, 24, 36 and 48 h post-dose; LLOQ 1 ng/L, with below-LLOQ values (55 samples from 24 h onwards across studies 1-3, Figure 5 caption) retained via the M3-style likelihood of Bauer 2019."
  )

  ini({
    # ==================================================================
    # Rosenborg 2025 Table 2, column 'SALM, Model 4 (RSE)' (studies 1-3
    # pooled, dose 150 ug); RSE quoted as the fraction reported in the
    # table.
    #
    # The supplement's NONMEM code (Fig. 1, 'NONMEM-code applied across
    # studies for SALM, Model 4') gives the parameterisation:
    #   $LEVEL STUD=(9[1])                           ; interstudy variability
    #   CL  = EXP(MU_1 + ETA(1))                     ; CL/F
    #   Q2  = EXP(MU_2 + ETA(2))                     ; Q2/F
    #   Q3  = EXP(MU_3 + ETA(3))                     ; Q3/F
    #   S1  = EXP(MU_4 + ETA(4))                     ; V1/F
    #   S2  = EXP(MU_5 + ETA(5))                     ; V2/F
    #   S3  = EXP(MU_6 + ETA(6))                     ; V3/F
    #   K41 = EXP(TREA1*MU_7 + TREA2*MU_10 + ETA(7))
    #   F4  = EXP(TREA1*(MU_11 + ETA(8) + ETA(9)) + TREA2*0)
    # so every parameter is estimated on the log scale. Unlike the FP
    # models, Q2/F and Q3/F carry random effects here (Table 2 reports
    # non-NA Q2/F_OMEGA_IIV and Q3/F_OMEGA_IIV only in the SALM column).
    #
    # Random-effect scale: Rosenborg 2025 Materials and Methods states
    # 'Random effects are presented as the standard deviation from an
    # associated typical parameter value, that is, the square root of the
    # estimated variance', so every *_OMEGA_* entry of Table 2 is an
    # omega (SD) and the variance below is its square.
    #
    # Which parameter carries the interstudy random effect: the
    # supplement's F4 line places ETA(9) on relative bioavailability, and
    # Figure 1 of the paper labels the SALM deposition compartment
    # 'F4_rel, ETA8/9, TREA' while Table 2 names the row
    # F4_rel_OMEGA_ISV. A sentence in Results / Population Analysis
    # instead describes 'the jointly estimated random effect of study on
    # apparent elimination and inter compartmental clearances of SALM';
    # the code, the schematic and the table agree against it, so ETA(9)
    # is encoded on F4_rel here (see the vignette's Errata).
    # ==================================================================

    # ---------------- Structural (typical-value) parameters ----------------
    lcl <- log(221.0)
    label("Log apparent elimination clearance CL/F (log(L/h))")                       # Rosenborg 2025 Table 2, Model 4: CL/F = 221.0 L/h (RSE 0.02500)
    lq <- log(1956)
    label("Log apparent distribution clearance Q2/F to peripheral1 (log(L/h))")       # Rosenborg 2025 Table 2, Model 4: Q2/F = 1956 L/h (RSE 0.03889)
    lq2 <- log(179.1)
    label("Log apparent distribution clearance Q3/F to peripheral2 (log(L/h))")       # Rosenborg 2025 Table 2, Model 4: Q3/F = 179.1 L/h (RSE 0.04003)
    lvc <- log(215.1)
    label("Log apparent central volume of distribution V1/F (log(L))")                # Rosenborg 2025 Table 2, Model 4: V1/F = 215.1 L (RSE 0.04616)
    lvp <- log(726.8)
    label("Log apparent first peripheral volume of distribution V2/F (log(L))")       # Rosenborg 2025 Table 2, Model 4: V2/F = 726.8 L (RSE 0.02669)
    lvp2 <- log(1494)
    label("Log apparent second peripheral volume of distribution V3/F (log(L))")      # Rosenborg 2025 Table 2, Model 4: V3/F = 1494 L (RSE 0.03437)

    # Monophasic absorption rate constant, one typical value per product
    # (supplement Fig. 1: K41 = EXP(TREA1*MU_7 + TREA2*MU_10 + ETA(7))).
    # Both strata carry an explicit suffix; neither keeps the bare lka.
    lka_test <- log(14.78)
    label("Log monophasic absorption rate constant k41 for the Wixela Inhub test product (log(1/h))")  # Rosenborg 2025 Table 2, Model 4: k41(test) = 14.78 1/h (RSE 0.04919)
    lka_ref <- log(20.42)
    label("Log monophasic absorption rate constant k41 for the Advair Diskus reference product (log(1/h))")  # Rosenborg 2025 Table 2, Model 4: k41(ref) = 20.42 1/h (RSE 0.01935)

    lfdepot <- log(1.028)
    label("Log relative extent of bioavailability F4_rel, test versus reference (unitless)")  # Rosenborg 2025 Table 2, Model 4: F4_rel = 1.028 (RSE 0.01442)

    # ---------------- Between-subject variability (Table 2 *_OMEGA_IIV rows) ----------------
    etalcl ~ 0.12229009                                                               # Rosenborg 2025 Table 2, Model 4: CL/F_OMEGA_IIV = 0.3497 (RSE 0.03261); 0.3497^2 = 0.12229009
    etalq ~ 0.26739241                                                                # Rosenborg 2025 Table 2, Model 4: Q2/F_OMEGA_IIV = 0.5171 (RSE 0.06330); 0.5171^2 = 0.26739241
    etalq2 ~ 0.28729600                                                               # Rosenborg 2025 Table 2, Model 4: Q3/F_OMEGA_IIV = 0.5360 (RSE 0.02264); 0.5360^2 = 0.287296
    etalvc ~ 0.19784704                                                               # Rosenborg 2025 Table 2, Model 4: V1/F_OMEGA_IIV = 0.4448 (RSE 0.05782); 0.4448^2 = 0.19784704
    etalvp ~ 0.11730625                                                               # Rosenborg 2025 Table 2, Model 4: V2/F_OMEGA_IIV = 0.3425 (RSE 0.03743); 0.3425^2 = 0.11730625
    etalvp2 ~ 0.22391824                                                              # Rosenborg 2025 Table 2, Model 4: V3/F_OMEGA_IIV = 0.4732 (RSE 0.01604); 0.4732^2 = 0.22391824
    etalka ~ 0.23040000                                                               # Rosenborg 2025 Table 2, Model 4: k41_OMEGA_IIV = 0.48000 (RSE 0.04152); 0.48^2 = 0.2304
    etalfdepot ~ 0.03297856                                                           # Rosenborg 2025 Table 2, Model 4: F4_rel_OMEGA_IIV = 0.1816 (RSE 0.06365); 0.1816^2 = 0.03297856

    # ---------------- Between-study (nested) variability (Table 2 *_OMEGA_ISV row) ----------------
    # ETA(9) of the supplement's F4 line, declared at the study level by
    # $LEVEL STUD=(9[1]). Precision is poor (RSE 0.5291) because only
    # three studies contribute.
    etalfdepot_study ~ 0.00085849 | SIDN                                              # Rosenborg 2025 Table 2, Model 4: F4_rel_OMEGA_ISV = 0.0293 (RSE 0.5291); 0.0293^2 = 0.00085849

    # ---------------- Residual variability (Table 2 CP_PROP / CP_ADD rows) ----------------
    # Supplement Fig. 1 $ERROR: W = SQRT(PLP**2*IPRED**2 + PLA**2),
    # Y = IPRED + SD*ERR(1) with PLP = CP_PROP and PLA = CP_ADD carried as
    # THETAs, i.e. nlmixr2's combined proportional-plus-additive form with
    # both magnitudes on the standard-deviation scale.
    propSd <- 0.1299
    label("Proportional residual error SD (fraction)")                                # Rosenborg 2025 Table 2, Model 4: CP_PROP = 0.1299 (RSE 0.01039)
    addSd <- 0.1444
    label("Additive residual error SD (ng/L)")                                        # Rosenborg 2025 Table 2, Model 4: CP_ADD = 0.1444 (RSE 0.1213)
  })

  model({
    # 1. Individual parameters. All eight subject-level etas of Model 4
    #    enter on the log scale.
    cl <- exp(lcl + etalcl)
    q <- exp(lq + etalq)
    q2 <- exp(lq2 + etalq2)
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp + etalvp)
    vp2 <- exp(lvp2 + etalvp2)

    # Product-specific absorption rate constant with a shared eta,
    # reproducing K41 = EXP(TREA1*MU_7 + TREA2*MU_10 + ETA(7)).
    ka <- exp(lka_test * FORM_WIXELA_INHUB + lka_ref * (1 - FORM_WIXELA_INHUB) + etalka)

    # 2. Micro-constants (supplement Fig. 1: K10 = CL/S1, K12 = Q2/S1,
    #    K21 = Q2/S2, K13 = Q3/S1, K31 = Q3/S3).
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # 3. ODE system. NONMEM ADVAN5 with
    #    $MODEL COMP=(DEFOBS1) COMP=(PERIPH1) COMP=(PERIPH2) COMP=(DEPOT1),
    #    i.e. cmt 1 = central, 2 = peripheral1, 3 = peripheral2,
    #    4 = the hypothetical pulmonary deposition site (depot), with the
    #    only depot transfer being k41 into central.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central -
      k12 * central + k21 * peripheral1 -
      k13 * central + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2

    # 4. Relative extent of bioavailability. Reference doses get F = 1
    #    exactly (F4 = EXP(TREA2*0)); test doses get F4_rel and carry both
    #    the subject-level (ETA8) and the study-level (ETA9) random
    #    effects on the test/reference ratio.
    f(depot) <- exp((lfdepot + etalfdepot + etalfdepot_study) * FORM_WIXELA_INHUB)

    # 5. Observation. The supplement's abbreviation list gives DOSE in ug
    #    and DV in ng/L while CL/F is in L/h and V1/F in L, so the amount
    #    in `central` is in ug and the factor 1000 converts ug/L to the
    #    ng/L scale on which the concentrations, the LLOQ of 1 ng/L and
    #    CP_ADD are reported.
    Cc <- 1000 * central / vc

    Cc ~ prop(propSd) + add(addSd)
  })
}
