Zou_2022_furmonertinib <- function() {
  description <- "Semi-mechanistic joint parent-metabolite population PK model for oral furmonertinib (AST2818, a third-generation irreversible EGFR TKI) and its active metabolite AST5902 (Zou 2022). The parent has two-transit-compartment absorption (rate constant ka shared for depot and both transits) feeding a two-compartment parent disposition (CL/F, Vc/F, Q/F, Vp/F). The parent is eliminated through the central compartment; a fraction Fm of the parent elimination becomes AST5902, which is described by a two-compartment metabolite (Clm/(F*Fm), Vcm/(F*Fm), and inter-compartment rate constants k67 and k76). Because no intravenous data were available, the absolute parent bioavailability F and the fraction Fm are non-identifiable and are absorbed into the apparent parameters. Autoinduction of furmonertinib metabolism (mediated by CYP3A4) is modelled as an indirect-response (IDR III) enzyme pool with unity baseline: d(A_ENZ)/dt = kENZ * (1 + S * Cc) - kENZ * A_ENZ, and the apparent parent clearance is CLbase/F * A_ENZ. Covariates: alkaline phosphatase (ALP, U/L) via power effect on CLbase/F and on Clm/(F*Fm) (median 77.2 U/L); body weight (WT, kg) via power effect on Clm/(F*Fm) (median 65 kg); and a categorical food-with-a-high-fat-meal effect on parent oral bioavailability (+22.4%) and on the fraction converted to AST5902 (-33.5%)."
  reference <- paste(
    "Zou HX, Zhang YF, Zhong DF, Jiang Y, Liu F, Zhao QY, Zuo Z,",
    "Zhang YF, Yan XY (2022).",
    "Effect of autoinduction and food on the pharmacokinetics of",
    "furmonertinib and its active metabolite characterized by a",
    "population pharmacokinetic model.",
    "Acta Pharmacologica Sinica 43(7):1865-1874.",
    "doi:10.1038/s41401-021-00798-y. PMID 34789919.",
    sep = " "
  )
  vignette <- "Zou_2022_furmonertinib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    ALP = list(
      description        = "Baseline serum alkaline phosphatase; liver-function marker used as a covariate on both parent (CLbase/F) and metabolite (Clm/(F*Fm)) apparent clearances via power scaling normalised to the cohort median 77.2 U/L (Zou 2022 Table 1 and Table 2).",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Cohort summary (Zou 2022 Table 1, N = 54): mean 88.1, SD 44.3, median 77.2, range 49.8-343 U/L. Two subjects had abnormally high ALP that appeared to drive the covariate significance in a sensitivity analysis (Discussion paragraph 3). Enters both CLbase/F and Clm/(F*Fm) with negative power exponents (larger ALP -> lower clearance). Reference 77.2 U/L for the power normalisation.",
      source_name        = "ALP"
    ),
    WT = list(
      description        = "Baseline body weight; covariate on the metabolite apparent clearance Clm/(F*Fm) via a power scaling normalised to the cohort median 65 kg (Zou 2022 Table 1 and Table 2).",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Cohort summary (Zou 2022 Table 1, N = 54): mean 66.5, SD 10.4, median 65, range 48-111 kg. Enters Clm/(F*Fm) only in the final model; the GAM-screened effects of body weight on other parameters were not retained by the stepwise SCM procedure (Results 'Covariate model' paragraph). Reference 65 kg for the power normalisation.",
      source_name        = "WT"
    ),
    FED_HIGHFAT = list(
      description        = "Categorical food-effect indicator: 1 = dose administered immediately after ingestion of a high-fat, high-calorie breakfast; 0 = fasted (overnight fast of at least 10 h). Time-fixed per dose record.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (fasted)",
      notes              = "Only Study 004 (n = 16 healthy males) tested the food effect via a crossover single-dose 80 mg design; Study 001 and Study 002 (NSCLC patients) dosed under fasted conditions (Zou 2022 Methods 'Study design' paragraphs 2-3). Modelled as two additive multiplicative effects (Zou 2022 covariate equations after Table 2): FTV = 1 + 0.224 * FOOD on parent oral bioavailability (F) and FmTV = 1 - 0.335 * FOOD on the fraction converted to AST5902 (Fm). The high-fat operational definition (FDA-style high-fat, high-calorie breakfast) motivates the FED_HIGHFAT canonical over the generic FED (only one meal-type was tested; the modeled effect quantifies a high-fat food-specific shift).",
      source_name        = "FOOD (Zou 2022 covariate equations)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 54L,
    n_studies      = 3L,
    age_range      = "21-68 years",
    age_median     = "51.5 years (mean 48, SD 13.9; Zou 2022 Table 1)",
    weight_range   = "48-111 kg",
    weight_median  = "65 kg (mean 66.5, SD 10.4; Zou 2022 Table 1)",
    sex_female_pct = 42.6,
    race_ethnicity = "Chinese (all studies conducted in China; Zou 2022 Methods paragraph 1)",
    disease_state  = "Pooled NSCLC patients with EGFR-sensitizing / T790M-resistance mutations (Study 001 dose escalation, n = 14; Study 002 dose expansion, n = 24) and healthy adult male volunteers (Study 004 food-effect crossover, n = 16). NSCLC patients had disease progression on prior first- or second-generation EGFR-TKI therapy.",
    dose_range     = "Furmonertinib 20, 40, 80, 160, or 240 mg orally once daily for 21 days per cycle in NSCLC patients (Study 001 and 002); single 80 mg oral dose (fed and fasted, two-period crossover) in healthy volunteers (Study 004).",
    regions        = "China (multicenter)",
    n_observations = "1,450 furmonertinib and 1,463 AST5902 plasma concentrations (LC/MS/MS; LLOQ 0.20 ng/mL parent, 0.050 ng/mL metabolite; Zou 2022 Methods 'Study design' paragraph 4 and Results 'Pharmacokinetic database description' paragraph 1).",
    notes          = "Baseline demographics from Zou 2022 Table 1. Hepatic dysfunction: 48 normal, 6 mild dysfunction. ClinicalTrials.gov identifiers NCT02973763 (Study 001), NCT03127449 (Study 002), NCT03926182 (Study 004). NONMEM 7.3 with FOCE-I; PsN 4.9.0 for SCM and bootstrap. Diagnostic and simulation post-processing in R 3.6.1. Furmonertinib is also known by its development codes AST2818 and alflutinib."
  )

  ini({
    # =========================================================================
    # Final-model structural parameters from Zou 2022 Table 2. Time unit is
    # HOURS. Furmonertinib doses are in milligrams; observed plasma
    # concentrations are in nanograms per millilitre (LLOQ 0.20 ng/mL parent,
    # 0.050 ng/mL AST5902 metabolite; Methods 'Study design' paragraph 4).
    # The observation equation in model() multiplies (central / vc) by 1000
    # to convert dose(mg)/Vc(L) = mg/L into ng/mL so the estimated S has
    # units of 1/(ng/mL) that match the paper's Cp scale.
    # =========================================================================

    # --- Parent (furmonertinib) apparent disposition ------------------------
    lcl_base <- log(70.5)   ; label("Apparent baseline clearance CLbase/F of furmonertinib (L/h)")                          # Zou 2022 Table 2: CLbase/F = 70.5 L/h, %RSE 6.41; bootstrap median 70.5 (61.3, 79.4)
    lvc      <- log(2897)   ; label("Apparent central volume of distribution Vc/F of furmonertinib (L)")                     # Zou 2022 Table 2: Vc/F = 2897 L, %RSE 6.35; bootstrap median 2837 (2344, 3216)
    lka      <- log(1.34)   ; label("Absorption and transit rate constant ka (1/h) shared across depot and two transit compartments") # Zou 2022 Table 2: ka = 1.34 1/h, %RSE 5.63; bootstrap median 1.32 (1.14, 1.48)
    lq       <- log(12.4)   ; label("Apparent inter-compartmental clearance Q/F of furmonertinib (L/h)")                     # Zou 2022 Table 2: Q/F = 12.4 L/h, %RSE 34.2; bootstrap median 12.6 (9.42, 53.0)
    lvp      <- log(1470)   ; label("Apparent peripheral volume of distribution Vp/F of furmonertinib (L)")                  # Zou 2022 Table 2: Vp/F = 1470 L, %RSE 26.7; bootstrap median 1499 (1113, 2561)

    # --- Autoinduction indirect-response enzyme pool (IDR III) --------------
    # The pool is unitless with A_ENZ(0) = 1. Baseline first-order production
    # rate = kENZ * A_ENZ(0) = kENZ so the pool sits at 1 with Cp = 0 (Zou
    # 2022 Methods 'Structural model building' paragraphs 4-5). Furmonertinib
    # plasma concentration Cp (ng/mL) stimulates the production rate; the
    # linear form (S * Cp) was retained after an initial Emax attempt
    # yielded extremely large Emax / SC50 (Zou 2022 Discussion paragraph 2).
    lkenz    <- log(0.00304); label("First-order rate constant for enzyme-pool degradation kENZ (1/h)")                       # Zou 2022 Table 2: kENZ = 0.00304 1/h, %RSE 6.56; bootstrap median 0.00311 (0.000952, 0.00560). log(2)/kENZ ~ 228 h, consistent with CYP3A4 half-life estimates in refs [27, 39-41].
    ls       <- log(0.0111) ; label("Slope of furmonertinib concentration on enzyme production rate S (1/(ng/mL))")           # Zou 2022 Table 2: S = 0.0111, %RSE 18.1; bootstrap median 0.0112 (0.00710, 0.0225). Units 1/(ng/mL) matched to the modelled Cc scale.

    # --- Metabolite (AST5902) apparent disposition --------------------------
    # Both Clm and Vcm are apparent; the (F * Fm) product is non-identifiable
    # from oral data alone and is absorbed into these apparent parameters
    # (Zou 2022 Methods 'Structural model building' paragraph 3). k67 and
    # k76 are distribution rate constants between the metabolite central and
    # peripheral compartments (Zou 2022 Fig. 1 legend and Methods paragraph 4).
    lcl_ast5902 <- log(119) ; label("Apparent clearance of AST5902 Clm/(F*Fm) (L/h)")                                          # Zou 2022 Table 2: Clm/(F*Fm) = 119 L/h, %RSE 3.56; bootstrap median 119 (111, 128)
    lvc_ast5902 <- log(291) ; label("Apparent central volume of distribution of AST5902 Vcm/(F*Fm) (L)")                       # Zou 2022 Table 2: Vcm/(F*Fm) = 291 L, %RSE 9.27; bootstrap median 286 (215, 341)
    lk67        <- log(0.952)  ; label("AST5902 central-to-peripheral distribution rate constant k67 (1/h)")                  # Zou 2022 Table 2: k67 = 0.952 1/h, %RSE 17.5; bootstrap median 0.972 (0.743, 1.67)
    lk76        <- log(0.0542) ; label("AST5902 peripheral-to-central distribution rate constant k76 (1/h)")                  # Zou 2022 Table 2: k76 = 0.0542 1/h, %RSE 9.20; bootstrap median 0.0540 (0.0471, 0.0711)

    # --- Covariate effects (power exponents; reference values are cohort medians) ---
    e_alp_cl_base    <- fixed(-0.505) ; label("Power exponent of ALP on CLbase/F (unitless); reference 77.2 U/L")             # Zou 2022 Table 2: theta_CLbase/F,ALP = -0.505, %RSE 17.5; bootstrap median -0.501 (-0.679, -0.235). Fixed here as a covariate coefficient without IIV, per convention.
    e_alp_cl_ast5902 <- fixed(-0.278) ; label("Power exponent of ALP on Clm/(F*Fm) (unitless); reference 77.2 U/L")           # Zou 2022 Table 2: theta_CLm/(F*Fm),ALP = -0.278, %RSE 32.7; bootstrap median -0.258 (-0.430, 0.0667)
    e_wt_cl_ast5902  <- fixed(0.622)  ; label("Power exponent of body weight on Clm/(F*Fm) (unitless); reference 65 kg")      # Zou 2022 Table 2: theta_CLm/(F*Fm),body weight = 0.622, %RSE 27.8; bootstrap median 0.629 (0.144, 0.956)
    e_fed_f          <- fixed(0.224)  ; label("Linear fed-vs-fasted coefficient on parent oral bioavailability F (unitless)") # Zou 2022 Table 2: theta_F,food = 0.224, %RSE 38.3; bootstrap median 0.225 (0.0500, 0.405). Applied as F = 1 + 0.224 * FED_HIGHFAT.
    e_fed_fm         <- fixed(-0.335) ; label("Linear fed-vs-fasted coefficient on fraction Fm of parent converted to AST5902 (unitless)") # Zou 2022 Table 2: theta_Fm,food = -0.335, %RSE 7.87; bootstrap median -0.332 (-0.382, -0.272). Applied as Fm_scaling = 1 - 0.335 * FED_HIGHFAT on the parent-to-metabolite formation flux.

    # =========================================================================
    # Inter-individual variability. Zou 2022 Table 2 reports the omega^2
    # variance directly for a log-normal (exponential) model:
    #     P_i = P_TV * exp(eta_i), eta_i ~ N(0, omega^2)
    # so the tabulated numbers are used verbatim as variances (no CV%
    # conversion). Shrinkage estimates are shown in parentheses in Table 2.
    # =========================================================================
    etalcl          ~ 0.0780   # Zou 2022 Table 2: omega^2 CL/F = 0.0780 (shrinkage 3.47%, %RSE 22.2)
    etalvc          ~ 0.144    # Zou 2022 Table 2: omega^2 Vc/F = 0.144 (shrinkage 3.39%, %RSE 27.1)
    etalka          ~ 0.161    # Zou 2022 Table 2: omega^2 ka = 0.161 (shrinkage 1.11%, %RSE 22.8)
    etalcl_ast5902  ~ 0.0485   # Zou 2022 Table 2: omega^2 Clm/(F*Fm) = 0.0485 (shrinkage 3.25%, %RSE 19.6)
    etalvc_ast5902  ~ 0.0970   # Zou 2022 Table 2: omega^2 Vcm/(F*Fm) = 0.0970 (shrinkage 11.4%, %RSE 29.6)

    # =========================================================================
    # Residual error. Zou 2022 Methods 'Stochastic model building' paragraph
    # 2: "an additive model of residual error was applied to log-transformed
    # data", i.e. Y_ij = ln(C_ij) + eps_ij, eps_ij ~ N(0, delta^2), which is
    # the additive-on-log <-> proportional-on-linear equivalence. Table 2
    # reports both delta values with "(CV%)" units, i.e. delta_parent = 0.336
    # (SD on log scale, ~ CV in linear space) and delta_metab = 0.275. Encoded
    # here as proportional residuals on the linear-scale plasma concentration.
    # =========================================================================
    propSd         <- 0.336   ; label("Proportional residual error on furmonertinib plasma concentration (fraction; additive-on-log-scale ~ CV in linear space)")   # Zou 2022 Table 2: delta_ADD_ERR = 33.6% (%RSE 8.05)
    propSd_ast5902 <- 0.275   ; label("Proportional residual error on AST5902 plasma concentration (fraction; additive-on-log-scale ~ CV in linear space)")         # Zou 2022 Table 2: delta_ADD_ERR_AST5902 = 27.5% (%RSE 8.54)
  })

  model({
    # ------------------------------------------------------------------------
    # 1. Individual structural parameters with covariate effects.
    # ------------------------------------------------------------------------
    # The apparent parent baseline clearance CLbase/F carries the ALP power
    # covariate (Zou 2022 Table 2 and covariate equation below Table 2:
    #   CLbase_TV = 70.5 * (ALP / 77.2) ^ (-0.505)).
    # The apparent metabolite clearance Clm/(F*Fm) carries both ALP and body
    # weight power covariates:
    #   CLm_TV = 119 * (ALP / 77.2) ^ (-0.278) * (WT / 65) ^ 0.622.
    # There is no reported IIV on Q/F, Vp/F, kENZ, S, k67, or k76.
    ka          <- exp(lka          + etalka)
    cl_base_iiv <- exp(lcl_base     + etalcl) * (ALP / 77.2) ^ e_alp_cl_base
    vc          <- exp(lvc          + etalvc)
    q           <- exp(lq)
    vp          <- exp(lvp)
    kenz        <- exp(lkenz)
    s_ind       <- exp(ls)
    cl_ast5902  <- exp(lcl_ast5902  + etalcl_ast5902) *
                     (ALP / 77.2) ^ e_alp_cl_ast5902 *
                     (WT  / 65)   ^ e_wt_cl_ast5902
    vc_ast5902  <- exp(lvc_ast5902  + etalvc_ast5902)
    k67         <- exp(lk67)
    k76         <- exp(lk76)

    # ------------------------------------------------------------------------
    # 2. Observations in ng/mL (dose in mg, volumes in L => central/vc in
    #    mg/L; multiply by 1000 to match the paper's ng/mL Cp scale so the
    #    autoinduction slope S = 0.0111 has units 1/(ng/mL) directly).
    # ------------------------------------------------------------------------
    Cc         <- 1000 * central         / vc
    Cc_ast5902 <- 1000 * central_ast5902 / vc_ast5902

    # ------------------------------------------------------------------------
    # 3. Autoinduction: time-varying apparent parent clearance CL/F(t) =
    #    CLbase/F * A_ENZ(t) with A_ENZ(0) = 1 (Zou 2022 Eq. 10).
    # ------------------------------------------------------------------------
    cl <- cl_base_iiv * enz_pool

    # ------------------------------------------------------------------------
    # 4. ODE system (Zou 2022 Fig. 1 and Methods 'Structural model building').
    #    Parent chain: depot -> transit1 -> transit2 -> central, all with
    #    the shared rate constant ka (Zou 2022 Methods paragraph 1: "The
    #    absorption and the transfer rate were quantified by ka"). Parent
    #    disposition is a two-compartment model (central + peripheral1).
    #    Metabolite disposition is a two-compartment model (central_ast5902
    #    + peripheral1_ast5902) with distribution rate constants k67 and
    #    k76 (Zou 2022 Methods paragraph 3).
    # ------------------------------------------------------------------------
    d/dt(depot)     <- -ka * depot
    d/dt(transit1)  <-  ka * depot    - ka * transit1
    d/dt(transit2)  <-  ka * transit1 - ka * transit2
    d/dt(central)   <-  ka * transit2 - (cl + q) / vc * central + q / vp * peripheral1
    d/dt(peripheral1) <- q / vc * central - q / vp * peripheral1

    # Metabolite formation flux enters central_ast5902 at rate (cl / vc) *
    # central (the total parent apparent elimination flux; in the (F * Fm)-
    # absorbed apparent-metabolite frame the fraction Fm is already carried
    # by the metabolite compartment's scaling). The food effect on Fm is
    # applied as a multiplicative modifier (1 + e_fed_fm * FED_HIGHFAT) on
    # this formation flux (Zou 2022 covariate equations below Table 2:
    # Fm_TV = 1 - 0.335 * FOOD).
    d/dt(central_ast5902)     <- (1 + e_fed_fm * FED_HIGHFAT) * cl / vc * central -
                                   cl_ast5902 / vc_ast5902 * central_ast5902 -
                                   k67 * central_ast5902 + k76 * peripheral1_ast5902
    d/dt(peripheral1_ast5902) <-  k67 * central_ast5902 - k76 * peripheral1_ast5902

    # Autoinduction enzyme pool (linear-slope IDR III form; Zou 2022 Eq. 11).
    d/dt(enz_pool) <- kenz * (1 + s_ind * Cc) - kenz * enz_pool
    enz_pool(0)    <- 1

    # ------------------------------------------------------------------------
    # 5. Parent oral bioavailability. F is anchored at 1 in the fasted state
    #    (all apparent parameters are reported at the fasted reference). Food
    #    shifts F by +22.4% (Zou 2022 covariate equations below Table 2:
    #    F_TV = 1 + 0.224 * FOOD).
    # ------------------------------------------------------------------------
    f(depot) <- 1 + e_fed_f * FED_HIGHFAT

    # ------------------------------------------------------------------------
    # 6. Observation model. Zou 2022 fit additive residual error on
    #    log-transformed data (Methods 'Stochastic model building' paragraph
    #    2); the CV% form reported in Table 2 corresponds to a proportional
    #    residual on linear-scale plasma concentration.
    # ------------------------------------------------------------------------
    Cc         ~ prop(propSd)
    Cc_ast5902 ~ prop(propSd_ast5902)
  })
}
