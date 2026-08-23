Rosenborg_2025_fluticasone_750ug <- function() {
  description <- paste(
    "Empirical population PK model for inhaled fluticasone propionate (FP) in",
    "healthy adults given 3 x 250 ug (750 ug total) from a dry powder inhaler",
    "in a two-way crossover bioequivalence study (Rosenborg 2025 Model 2,",
    "study 2). Monophasic first-order absorption from a hypothetical",
    "pulmonary deposition site into a three-compartment apparent systemic",
    "disposition; the nominal inhaled dose is assumed to be deposited",
    "instantaneously and completely at the deposition site, so all clearances",
    "and volumes are apparent (CL/F, V1/F, ...) and carry an unknown absolute",
    "bioavailability < 1. The test product (Wixela Inhub) differs from the",
    "reference (Advair Diskus) in the absorption rate constant k41 and in",
    "relative bioavailability F4_rel; reference bioavailability is anchored at",
    "F = 1. Between-subject variability is carried on CL/F, V1/F, V2/F, V3/F,",
    "k41 and F4_rel; the distribution clearances Q2/F and Q3/F were fitted",
    "without random effects because of shrinkage. Fitted separately by",
    "study/dose - see Rosenborg_2025_fluticasone_300ug and",
    "Rosenborg_2025_fluticasone_1500ug for the other two FP doses and",
    "Rosenborg_2025_salmeterol for the co-administered salmeterol.")
  reference <- paste(
    "Rosenborg J, Backman P, Bengtsson T, Haughie S.",
    "Relative Bioavailability of Inhaled Fluticasone Propionate and Salmeterol",
    "- is Population Pharmacokinetic Modelling a Relevant Alternative to a",
    "Non-Compartmental Approach?",
    "Drug Des Devel Ther. 2025;19:9653-9670. doi:10.2147/DDDT.S480189")
  vignette <- "Rosenborg_2025_fluticasone_salmeterol"
  units    <- list(time = "h", dosing = "ug", concentration = "ng/L")

  compartmentData <- list(
    depot       = list(analyte = "fluticasone propionate", units = "ug", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "fluticasone propionate", units = "ug", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "fluticasone propionate", units = "ug", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "fluticasone propionate", units = "ug", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    FORM_WIXELA_INHUB = list(
      description        = "Indicator that the inhalation was taken from the Wixela Inhub dry powder inhaler (test product), 1 = Wixela Inhub, 0 = Advair Diskus.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (Advair Diskus, the reference product; its relative extent of bioavailability is the anchor F = 1).",
      notes              = "Per-dose-record indicator. In this two-way crossover the same subject carries 1 on one period's dose row and 0 on the other's. The supplement's NONMEM code writes TREA = 1 for test and TREA = 2 for reference and derives TREA1 / TREA2 indicators from it; FORM_WIXELA_INHUB = TREA1. It gates two parameters (Rosenborg 2025 Figure 1): the monophasic absorption rate constant, K41 = EXP(TREA1*MU_7 + TREA2*MU_10 + ETA(7)), i.e. a separate typical value per product sharing one eta; and relative bioavailability, F4 = EXP(TREA1*(MU_11 + ETA(8)) + TREA2*0), so reference doses get F = 1 exactly while test doses get F4_rel carrying its own inter-individual variability.",
      source_name        = "TREA (supplement Sect. 1: 'Treatment alternative, test (TREA=1) and reference (TREA=2) formulation')"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 66L,
    n_studies      = 1L,
    age_range      = "mean age 37.7 years (study 2; Rosenborg 2025 Materials and Methods)",
    sex_female_pct = 54.5,
    race_ethnicity = c(White = NA_real_, Black = NA_real_, Asian = NA_real_, Other = NA_real_),
    disease_state  = "Healthy adult volunteers.",
    dose_range     = "Three inhalations of fluticasone propionate 250 ug (750 ug total) plus salmeterol 50 ug (150 ug total), single dose per period, two periods with a 7-day washout.",
    regions        = "USA",
    n_observations = "7184 fluticasone propionate plasma concentrations across all three studies (Rosenborg 2025 Results); study 2 contributes the data for this model.",
    notes          = "Two-way crossover: each subject received both the test (Wixela Inhub) and the reference (Advair Diskus) product. 66 subjects were recruited, of whom 61 to 65 completed both treatment periods. Panels were homogeneous, black and white, mean BMI approximately 26. Rosenborg 2025 states that 'neither demographic nor other covariates were considered in the evaluation of NCA-based results and therefore not in this alternative model analysis either', so product is the only covariate in the model. Sampling: pre-dose and 2, 5, 10, 15, 20, 30 and 45 min, and 1, 1.5, 2, 3, 4, 6, 8, 12, 24, 36 and 48 h post-dose; LLOQ 1 ng/L, with below-LLOQ values retained via the M3-style likelihood of Bauer 2019 (one sample below LLOQ from 24 h onwards in this study, Figure 5 caption)."
  )

  ini({
    # ==================================================================
    # Rosenborg 2025 Table 2, column 'FP, Model 2 (RSE)' (study 2, dose
    # 750 ug); RSE quoted as the fraction reported in the table.
    #
    # The supplement's NONMEM code (Fig. 1, 'NONMEM-code applied by
    # study/dose for FP, Models 1-3') gives the parameterisation:
    #   CL  = EXP(MU_1 + ETA(1))                     ; CL/F
    #   Q2  = THETA(2) + ETA(2)                      ; Q2/F  (OMEGA fixed 0 here)
    #   Q3  = THETA(3) + ETA(3)                      ; Q3/F  (OMEGA fixed 0 here)
    #   S1  = EXP(MU_4 + ETA(4))                     ; V1/F
    #   S2  = EXP(MU_5 + ETA(5))                     ; V2/F
    #   S3  = EXP(MU_6 + ETA(6))                     ; V3/F
    #   K41 = EXP(TREA1*MU_7 + TREA2*MU_10 + ETA(7))
    #   F4  = EXP(TREA1*(MU_11 + ETA(8)) + TREA2*0)
    # so every parameter except Q2/F and Q3/F is estimated on the log
    # scale. Q2/F and Q3/F carry no random effect in Models 1-3
    # ('disregarding between subject variability of apparent distribution
    # clearances of FP owing to shrinkage', Results / Population
    # Analysis), which is why Table 2 reports NA for their IIV; with the
    # eta removed the linear THETA and the log-transformed canonical
    # lq / lq2 describe exactly the same typical value.
    #
    # Random-effect scale: Rosenborg 2025 Materials and Methods states
    # 'Random effects are presented as the standard deviation from an
    # associated typical parameter value, that is, the square root of the
    # estimated variance', so every *_OMEGA_IIV entry of Table 2 is an
    # omega (SD) and the variance below is its square.
    # ==================================================================

    # ---------------- Structural (typical-value) parameters ----------------
    lcl <- log(618.4)
    label("Log apparent elimination clearance CL/F (log(L/h))")                       # Rosenborg 2025 Table 2, Model 2: CL/F = 618.4 L/h (RSE 0.04010)
    lq <- log(1419)
    label("Log apparent distribution clearance Q2/F to peripheral1 (log(L/h))")       # Rosenborg 2025 Table 2, Model 2: Q2/F = 1419 L/h (RSE 0.1496)
    lq2 <- log(355.7)
    label("Log apparent distribution clearance Q3/F to peripheral2 (log(L/h))")       # Rosenborg 2025 Table 2, Model 2: Q3/F = 355.7 L/h (RSE 0.04612)
    lvc <- log(255.3)
    label("Log apparent central volume of distribution V1/F (log(L))")                # Rosenborg 2025 Table 2, Model 2: V1/F = 255.3 L (RSE 0.1079)
    lvp <- log(47.14)
    label("Log apparent first peripheral volume of distribution V2/F (log(L))")       # Rosenborg 2025 Table 2, Model 2: V2/F = 47.14 L (RSE 0.3307); the least precisely estimated parameter of the model, flagged in Results
    lvp2 <- log(4720)
    label("Log apparent second peripheral volume of distribution V3/F (log(L))")      # Rosenborg 2025 Table 2, Model 2: V3/F = 4720 L (RSE 0.05164)

    # Monophasic absorption rate constant, one typical value per product
    # (supplement Fig. 1: K41 = EXP(TREA1*MU_7 + TREA2*MU_10 + ETA(7))).
    # Both strata carry an explicit suffix; neither keeps the bare lka.
    lka_test <- log(0.2307)
    label("Log monophasic absorption rate constant k41 for the Wixela Inhub test product (log(1/h))")  # Rosenborg 2025 Table 2, Model 2: k41(test) = 0.2307 1/h (RSE 0.04234)
    lka_ref <- log(0.2525)
    label("Log monophasic absorption rate constant k41 for the Advair Diskus reference product (log(1/h))")  # Rosenborg 2025 Table 2, Model 2: k41(ref) = 0.2525 1/h (RSE 0.01656)

    lfdepot <- log(1.056)
    label("Log relative extent of bioavailability F4_rel, test versus reference (unitless)")  # Rosenborg 2025 Table 2, Model 2: F4_rel = 1.056 (RSE 0.04050)

    # ---------------- Between-subject variability (Table 2 *_OMEGA_IIV rows) ----------------
    etalcl ~ 0.09985600                                                               # Rosenborg 2025 Table 2, Model 2: CL/F_OMEGA_IIV = 0.3160 (RSE 0.1026); 0.3160^2 = 0.099856
    etalvc ~ 0.66292164                                                               # Rosenborg 2025 Table 2, Model 2: V1/F_OMEGA_IIV = 0.8142 (RSE 0.1015); 0.8142^2 = 0.66292164
    etalvp ~ 2.69288100                                                               # Rosenborg 2025 Table 2, Model 2: V2/F_OMEGA_IIV = 1.641 (RSE 0.1185); 1.641^2 = 2.692881. Results notes this eta distribution was 'remarkably skewed' after FP 750 and 1500 ug.
    etalvp2 ~ 0.09690769                                                              # Rosenborg 2025 Table 2, Model 2: V3/F_OMEGA_IIV = 0.3113 (RSE 0.2510); 0.3113^2 = 0.09690769
    etalka ~ 0.04791721                                                               # Rosenborg 2025 Table 2, Model 2: k41_OMEGA_IIV = 0.2189 (RSE 0.1274); 0.2189^2 = 0.04791721
    etalfdepot ~ 0.07529536                                                           # Rosenborg 2025 Table 2, Model 2: F4_rel_OMEGA_IIV = 0.2744 (RSE 0.1028); 0.2744^2 = 0.07529536
    # No etalq / etalq2: Table 2 reports Q2/F_OMEGA_IIV = Q3/F_OMEGA_IIV = NA for Models 1-3.

    # ---------------- Residual variability (Table 2 CP_PROP / CP_ADD rows) ----------------
    # Supplement Fig. 1 $ERROR: W = SQRT(PLP**2*IPRED**2 + PLA**2),
    # Y = IPRED + SD*ERR(1) with PLP = CP_PROP and PLA = CP_ADD carried as
    # THETAs, i.e. nlmixr2's combined proportional-plus-additive form with
    # both magnitudes on the standard-deviation scale.
    propSd <- 0.2261
    label("Proportional residual error SD (fraction)")                                # Rosenborg 2025 Table 2, Model 2: CP_PROP = 0.2261 (RSE 0.02080)
    addSd <- 0.5436
    label("Additive residual error SD (ng/L)")                                        # Rosenborg 2025 Table 2, Model 2: CP_ADD = 0.5436 (RSE 0.1655)
  })

  model({
    # 1. Individual parameters. Q2/F and Q3/F have no random effect in
    #    Models 1-3 (see the ini() note).
    cl <- exp(lcl + etalcl)
    q <- exp(lq)
    q2 <- exp(lq2)
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
    #    exactly (F4 = EXP(TREA2*0)); test doses get F4_rel and carry the
    #    inter-individual variability of the test/reference ratio.
    f(depot) <- exp((lfdepot + etalfdepot) * FORM_WIXELA_INHUB)

    # 5. Observation. The supplement's abbreviation list gives DOSE in ug
    #    and DV in ng/L while CL/F is in L/h and V1/F in L, so the amount
    #    in `central` is in ug and the factor 1000 converts ug/L to the
    #    ng/L scale on which the concentrations, the LLOQ of 1 ng/L and
    #    CP_ADD are reported.
    Cc <- 1000 * central / vc

    Cc ~ prop(propSd) + add(addSd)
  })
}
