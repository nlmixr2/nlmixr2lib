Hand_2019_benzathine_benzylpenicillin_g <- function() {
  description <- "One-compartment population PK model for benzylpenicillin released from an intramuscular depot of benzathine benzylpenicillin G (Bicillin L-A) with two sequential first-order absorption stages (t1/2,abs-1 ~ 0.46 day for dissolution of the crystal suspension, t1/2,abs-2 ~ 8.9 day for the rate-limiting release into plasma) and a fixed, allometrically scaled elimination rate constant. Fat-free mass scales the apparent volume of distribution, and a body mass index at or above 25 kg/m^2 increases the slow absorption half-life by 86.5%, producing flip-flop kinetics in which absorption, not clearance, determines the observed terminal half-life. Developed from 256 dried-blood-spot benzylpenicillin concentrations collected over six monthly injection cycles in 18 children and adolescents receiving secondary prophylaxis for rheumatic heart disease (Hand 2019)."
  reference <- paste(
    "Hand RM, Salman S, Newall N, Vine J, Page-Sharp M, Bowen AC, Gray K,",
    "Baker A, Kado J, Joseph J, Marsh J, Ramsay J, Sika-Paotonu D, Batty KT,",
    "Manning L, Carapetis J. A population pharmacokinetic study of benzathine",
    "benzylpenicillin G administration in children and adolescents with",
    "rheumatic heart disease: new insights for improved secondary prophylaxis",
    "strategies. J Antimicrob Chemother. 2019;74(7):1984-1991.",
    "doi:10.1093/jac/dkz076.",
    "The fixed-elimination, absorption-limited structural approach is shared",
    "with the same group's later analyses; see",
    "modellib('Kado_2020_benzathine_benzylpenicillin_g') and",
    "modellib('Kado_2023_benzathine_benzylpenicillin_g').",
    sep = " "
  )
  vignette <- "Hand_2019_benzathine_benzylpenicillin_g"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  # Doses are expressed as milligrams of the benzathine benzylpenicillin G salt
  # (Methods, 'Clinical study procedures': "Bicillin L-A ... 2.3 mL containing
  # 900 mg (1.2 MIU) of BPG"), and the measured analyte is benzylpenicillin
  # (Methods, 'Measuring penicillin from DBS'). V is therefore an apparent
  # volume V/F that absorbs both the salt-to-penicillin mass conversion and the
  # unknown bioavailable fraction; no separate F is estimated. Concentrations
  # were assayed in dried blood spots but the paper analysed and reports them as
  # plasma concentrations (Methods, 'Pharmacokinetic modelling and
  # simulations': "Log_e plasma concentration-time datasets for
  # benzylpenicillin"), so `central` is annotated as plasma.
  compartmentData <- list(
    depot = list(
      analyte = "benzathine benzylpenicillin g",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    transit1 = list(
      analyte = "benzathine benzylpenicillin g",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "benzylpenicillin",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    )
  )

  covariateData <- list(
    FFM = list(
      description = "Fat-free mass, the allometric size descriptor for the apparent volume of distribution and for the fixed elimination rate constant.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Reference 70 kg: Table 3 reports kel as 'h^-1 . 70 kg^-1' and V as",
        "'L . 70 kg^-1'. Results: 'Fat-free mass was the best size parameter for",
        "allometric scaling on V.' The same size descriptor also carries the",
        "elimination rate constant. The paper does not restate that explicitly,",
        "but it is forced by its own equivalence claim: kel was fixed 'with",
        "allometric scaling with an exponential of -1/4 (equivalent to an",
        "exponential of 3/4 for CL and 1 for V)'. Because CL = kel * V, the",
        "exponents -1/4 and 1 sum to 3/4 only if kel and V are scaled by the SAME",
        "size descriptor; scaling V by fat-free mass and kel by total body weight",
        "would not produce a clean 3/4-power clearance and would contradict the",
        "paper's stated equivalence. Both are therefore scaled by FFM against a",
        "70 kg reference. Absorption parameters are NOT allometrically scaled.",
        "The paper derived FFM from body weight and body mass index using the",
        "published model of Anderson BJ, Holford NH, Drug Metab Pharmacokinet",
        "2009;24:25-36 (reference 20; Methods, 'Pharmacokinetic modelling and",
        "simulations': 'Fat-free mass was estimated from weight and BMI from a",
        "published model in children'). That reference reparameterises the",
        "Janmahasatian et al. equation already registered for this column;",
        "FFM = WHSmax * HT^2 * WT / (WHS50 * HT^2 + WT) is algebraically",
        "identical to WHSmax * WT / (WHS50 + BMI), with WHSmax = 42.92 and",
        "WHS50 = 30.93 for males and WHSmax = 37.99 and WHS50 = 35.98 for",
        "females. Hand 2019 itself prints none of those constants; a downstream",
        "user must either supply FFM directly or compute it from WT, HT and SEXF.",
        "Assumed time-fixed at baseline.",
        sep = " "
      ),
      source_name = "FFM"
    ),
    BMI = list(
      description = "Body mass index at baseline; dichotomised inside model() at the 25 kg/m^2 threshold the paper selected, and applied to the slow absorption half-life t1/2,abs-2.",
      units = "kg/m^2",
      type = "continuous",
      reference_category = "BMI < 25 kg/m^2 (the reference stratum, n = 10 of 18)",
      notes = paste(
        "The model consumes BMI as a CONTINUOUS column and forms the binary",
        "indicator internally, because the paper's covariate is explicitly a",
        "thresholded version of a continuous measurement: 'Although many body",
        "composition covariates were correlated with k_a-2, BMI as a categorical",
        "variable, with a threshold of >=25 kg/m^2, resulted in the best fit and",
        "was associated with an 86.5% increase in t1/2,abs-2' (Results). Applied",
        "as the linear multiplier t1/2,abs-2 = theta * (1 + 0.865 * [BMI >= 25]).",
        "BMI is also the second input to the fat-free-mass derivation described",
        "under FFM, so a downstream dataset that carries BMI and WT can produce",
        "both covariate columns. The threshold is the conventional",
        "overweight cut-off and is the same 25 kg/m^2 boundary the group used to",
        "stratify enrolment in its later phase 1 study",
        "(modellib('Kado_2023_benzathine_benzylpenicillin_g')). Assumed",
        "time-fixed at baseline. No other covariate relationship reached",
        "significance (Results: 'No other significant covariate relationships",
        "were identified').",
        sep = " "
      ),
      source_name = "BMI"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 18L,
    n_studies = 1L,
    age_range = "7.9-17.7 years",
    age_median = "14.1 years",
    weight_range = "29.9-149 kg",
    weight_median = "62.9 kg",
    height_range = "1.36-2.05 m",
    height_median = "1.61 m",
    bmi_range = "16.2-44.4 kg/m^2",
    bmi_median = "23.6 kg/m^2",
    sex_female_pct = 56,
    race_ethnicity = c(Aboriginal = 77.8, Maori = 11.1, Samoan = 11.1),
    disease_state = "children and adolescents with a history of acute rheumatic fever or established rheumatic heart disease, receiving monthly benzathine benzylpenicillin G as secondary prophylaxis; none had known established renal disease",
    dose_range = "900 mg (1.2 MIU) benzathine benzylpenicillin G by intramuscular injection into the upper outer gluteal quadrant, alternating sides, once every 28 days; all 18 analysed participants received the full 900 mg dose",
    regions = "Australia (metropolitan Perth, Western Australia)",
    notes = paste(
      "Longitudinal observational study conducted March-November 2017 through the",
      "Princess Margaret Hospital ambulatory care service. 22 participants were",
      "enrolled and 4 withdrew before contributing samples; 18 contributed 256",
      "benzylpenicillin concentrations, 16 (89%) with full datasets across all six",
      "monthly injection cycles. Samples were dried blood spots assayed by a",
      "validated LC-MS/MS method (LLOQ 0.0025 mg/L, LOD 0.001 mg/L); 25 (9.7%)",
      "concentrations were below the limit of quantification and were retained",
      "using the M3 likelihood method of Beal 2001. Intensive sampling on days 1,",
      "3, 6, 12 and 21 after the injection was performed in two of the six cycles,",
      "with an additional trough before each injection and unscheduled samples",
      "triggered by sore throat. Creatinine was not measured. Estimation was by",
      "NONMEM 7.2.0 LAPLACIAN with INTER on log_e-transformed concentrations.",
      "Baseline demographics: Table 1. Final estimates and bootstrap: Table 3.",
      "Individual post hoc absorption half-lives and exposure metrics by BMI",
      "stratum: Table 2. Model schematic: Figure 1. Goodness of fit: Figure 2.",
      "Prediction-corrected VPC stratified by BMI: Figure 3. Growth-chart-based",
      "dosing simulations: Figure 4.",
      sep = " "
    )
  )

  ini({
    # ==========================================================================
    # Disposition
    # ==========================================================================
    # kel was FIXED, not estimated: Table 3 gives the bootstrap column for kel as
    # "fixed", and Results explains why -- "Initial analysis using standard
    # compartmental modelling with various absorption models resulted in
    # estimates of elimination t1/2 that were much longer than previously
    # reported for benzylpenicillin. Therefore, the elimination rate constant was
    # fixed with allometric scaling with an exponential of -1/4 ... based on
    # previously published data in children receiving intravenous
    # benzylpenicillin." Table 3 reports it per hour; converted to 1/day because
    # every absorption parameter in this model is expressed in days. The same
    # 1.32 h^-1 . 70 kg^-1 value is carried forward by the group's later
    # analysis, modellib('Kado_2023_benzathine_benzylpenicillin_g').
    lkel <- fixed(log(1.32 * 24))
    label("Elimination rate constant at 70 kg fat-free mass (1/day)") # Table 3: kel = 1.32 h^-1 . 70 kg^-1, bootstrap 'fixed' -> 31.68 1/day

    lvc <- log(72.2)
    label("Apparent central volume of distribution V/F at 70 kg fat-free mass (L)") # Table 3: V = 72.2 L . 70 kg^-1 (bootstrap 72.0, 95% CI 64.0-84.2)

    # Allometric exponents were imposed a priori, not estimated (Results: the
    # exponent of -1/4 on kel is "equivalent to an exponential of 3/4 for CL and
    # 1 for V"), so both are fixed.
    e_ffm_vc <- fixed(1)
    label("Allometric exponent of fat-free mass on V/F (unitless), applied a priori") # Results: 'equivalent to an exponential of 3/4 for CL and 1 for V'
    e_ffm_kel <- fixed(-0.25)
    label("Allometric exponent of fat-free mass on kel (unitless), applied a priori") # Results: kel fixed 'with allometric scaling with an exponential of -1/4'

    # ==========================================================================
    # Absorption
    # ==========================================================================
    # Figure 1: Bolus --k_a-1--> Absorption --k_a-2--> V --k_el-->. Two
    # sequential FIRST-ORDER stages; Results states "First-order absorption for
    # both these stages performed better than models with zero-order process"
    # and "The addition of peripheral compartment(s) did not improve the model."
    # Both are parameterised by half-life, which model() converts to a rate
    # constant as k = log(2) / t1/2.
    lthalf_abs1 <- log(0.455)
    label("Fast absorption half-life from the injection depot, t1/2,abs-1 (day)") # Table 3: t1/2,abs-1 = 0.455 days (bootstrap 0.461, 95% CI 0.174-0.948)
    lthalf_abs2 <- log(8.88)
    label("Slow absorption half-life into the central compartment, t1/2,abs-2 (day)") # Table 3: t1/2,abs-2 = 8.88 days (bootstrap 8.79, 95% CI 5.71-12.5)

    # Linear fractional increase in t1/2,abs-2 for the higher-BMI stratum.
    # Estimated (the bootstrap column reports a median and a CI), not fixed.
    e_bmi_thalf_abs2 <- 0.865
    label("Fractional increase in t1/2,abs-2 when BMI is at or above 25 kg/m^2 (unitless)") # Table 3: 'increase in t1/2,abs-2 with BMI >=25 kg/m2 (%)' = 86.5 (bootstrap 86.8, 95% CI 33.4-198)

    # ==========================================================================
    # Inter-individual variability
    # ==========================================================================
    # Table 3 footnote: 'IIV, IOV and RV are presented as 100% x sqrt(variability
    # estimate)'. The tabulated percentage is therefore 100 x omega, so
    # omega^2 = (percentage / 100)^2 directly -- NOT log(1 + CV^2). Results:
    # 'A full covariance matrix model was used', so all three etas form one
    # block.
    #
    # Off-diagonals. The paper prints two of the three correlations:
    #   r(t1/2,abs-1, t1/2,abs-2) = -1     (Table 3, bootstrap column 'fixed')
    #   r(t1/2,abs-2, V)          = -0.746 (bootstrap -0.808, 95% CI -1.00 to -0.316)
    # The third is forced by the first: with r12 exactly -1 the two absorption
    # etas are perfectly collinear, eta1 = -(sd1/sd2) * eta2, so
    #   r(t1/2,abs-1, V) = -r(t1/2,abs-2, V) = +0.746.
    #
    # A correlation of exactly -1 makes the 3x3 block singular and rxode2's
    # Cholesky sampler cannot decompose it. Following the repository convention
    # for published perfect correlations, the OFF-DIAGONALS are scaled by 0.99
    # (so r12 becomes -0.99); the diagonal variances -- the published IIV values
    # -- are untouched. The resulting matrix is positive definite (smallest
    # eigenvalue 4.8e-03).
    #
    #   var(t1/2,abs-1) = 0.78^2 = 0.6084
    #   var(t1/2,abs-2) = 0.63^2 = 0.3969
    #   var(V)          = 0.26^2 = 0.0676
    #   cov(abs-1, abs-2) = 0.99 * (-1.000) * 0.78 * 0.63 = -0.486486
    #   cov(abs-1, V)     = 0.99 * (+0.746) * 0.78 * 0.26 = +0.149776
    #   cov(abs-2, V)     = 0.99 * (-0.746) * 0.63 * 0.26 = -0.120973
    etalthalf_abs1 + etalthalf_abs2 + etalvc ~ c(
      0.6084,
      -0.486486, 0.3969,
      0.149776, -0.120973, 0.0676
    )

    # ==========================================================================
    # Inter-occasion variability
    # ==========================================================================
    # Table 3: 'IOV in t1/2,abs-2' = 30% (shrinkage 46%; bootstrap 31, 95% CI
    # 20-48); omega^2 = 0.30^2 = 0.09. An occasion is one monthly injection
    # cycle, and the source data span six of them. rxode2 parses but cannot
    # simulate the native 'eta ~ var | occ' multi-level syntax, so this is
    # encoded as a single occasion-indexed eta (the repository's registered
    # etaiov_<param>_<occasion> form) that is drawn once per subject per solve.
    # A single-cycle simulation reproduces the source exactly; a multi-cycle
    # simulation reuses the occasion-1 draw across every cycle and therefore
    # under-represents within-subject cycle-to-cycle variation. See the
    # vignette's Assumptions and deviations section.
    etaiov_lthalf_abs2_1 ~ 0.09

    # ==========================================================================
    # Residual variability
    # ==========================================================================
    # Concentrations were analysed on the natural-log scale (Methods,
    # 'Pharmacokinetic modelling and simulations': 'Log_e plasma
    # concentration-time datasets ... were analysed'), so the additive
    # log-scale residual is a proportional error model on the linear scale.
    propSd <- 0.35
    label("Proportional residual error (fraction)") # Table 3: RV = 35% (shrinkage 13%; bootstrap 34, 95% CI 30-38)
  })

  model({
    # ---- Reference values ----------------------------------------------------
    ffmRef <- 70 # kg      -- Table 3 reports kel and V per 70 kg
    bmiCut <- 25 # kg/m^2  -- Results: 'BMI as a categorical variable, with a threshold of >=25 kg/m2'

    # ---- Covariate dichotomisation ------------------------------------------
    # The paper's covariate is the indicator [BMI >= 25 kg/m^2]; the model takes
    # continuous BMI and applies the published threshold here.
    bmiHigh <- 0
    if (BMI >= bmiCut) {
      bmiHigh <- 1
    }

    # ---- Disposition ---------------------------------------------------------
    # Volume carries the estimated typical value, the a priori exponent of 1 on
    # fat-free mass, and the only disposition IIV. kel is a fixed rate constant
    # carrying the a priori exponent of -1/4 on the same size descriptor and no
    # IIV, so clearance kel * V scales with the 3/4 power of fat-free mass, which
    # is exactly the equivalence the Results section asserts.
    vc <- exp(lvc + etalvc) * (FFM / ffmRef)^e_ffm_vc
    kel <- exp(lkel) * (FFM / ffmRef)^e_ffm_kel

    # ---- Absorption ----------------------------------------------------------
    # t1/2,abs-1 carries IIV only. t1/2,abs-2 carries IIV, IOV, and the linear
    # BMI effect applied to the typical value.
    # The IOV term is applied on a second line so that the mu-referenced line
    # carries exactly one subject-level random effect, as rxode2 requires.
    thalf_abs1 <- exp(lthalf_abs1 + etalthalf_abs1)
    thalf_abs2Base <- exp(lthalf_abs2 + etalthalf_abs2)
    thalf_abs2 <- thalf_abs2Base * exp(etaiov_lthalf_abs2_1) *
      (1 + e_bmi_thalf_abs2 * bmiHigh)

    ka1 <- log(2) / thalf_abs1
    ka2 <- log(2) / thalf_abs2

    # ---- Structure (Figure 1) ------------------------------------------------
    # depot ("Bolus") --ka1--> transit1 ("Absorption") --ka2--> central (V) --kel-->
    d/dt(depot) <- -ka1 * depot
    d/dt(transit1) <- ka1 * depot - ka2 * transit1
    d/dt(central) <- ka2 * transit1 - kel * central

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
