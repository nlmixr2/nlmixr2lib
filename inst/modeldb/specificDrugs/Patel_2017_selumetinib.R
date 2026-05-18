Patel_2017_selumetinib <- function() {
  description <- "Sequential two-compartment population PK model for oral selumetinib (AZD6244, ARRY-142886) and its active metabolite N-desmethyl-selumetinib in adults with advanced solid tumors pooled with children with recurrent low-grade glioma (Patel 2017). Selumetinib disposition uses sequential zero-order (release into the gut compartment over duration D1 with lag ALAG1) and first-order (rate Ka) absorption with bioavailability anchored at 1 under fasted conditions and reduced by an additive food-effect coefficient under fed conditions; D1 and ALAG1 carry additive food-effect coefficients. Body surface area (power on CL/F and Vc/F), age (power on Vc/F), and alanine aminotransferase (negative power on CL/F) modify selumetinib parameters; BSA (negative power) modifies the fraction metabolized to N-desmethyl-selumetinib. The metabolite is two-compartment with its central volume fixed equal to the parent central volume to resolve identifiability; metabolite clearance and intercompartmental clearance are apparent values."
  reference <- paste(
    "Patel YT, Daryani VM, Patel P, Zhou D, Fangusaro J, Carlile DJ, Martin PD, Aarons L, Stewart CF.",
    "Population pharmacokinetics of selumetinib and its metabolite N-desmethyl-selumetinib in adult patients with",
    "advanced solid tumors and children with low-grade gliomas.",
    "CPT Pharmacometrics Syst Pharmacol. 2017;6(5):305-314. doi:10.1002/psp4.12175"
  )
  vignette <- "Patel_2017_selumetinib"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    BSA = list(
      description        = "Body surface area (Gehan and George formula)",
      units              = "m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed at baseline. Power effects on selumetinib CL/F (exponent +0.923; positive correlation per text), Vc/F (exponent +1.24; positive correlation), and on fraction metabolized Fm (paper text 'negative correlation', encoded as (BSA/BSA_ref)^(-0.908) so the magnitude theta20 = +0.908 from Table 2 maps to the documented negative-correlation direction). BSA was computed via the Gehan and George formula per Patel 2017 Methods. Reference value 1.66 m^2 is the approximate pooled-cohort median across Studies 16+20+29 (n=105) derived from study-level medians in Table 1 (1.76, 1.97, 1.45); the paper does not state the exact normalization constant used in NONMEM (see vignette Assumptions).",
      source_name        = "BSA"
    ),
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed at baseline. Power effect on selumetinib Vc/F (exponent +0.327; positive correlation per text). Reference value 53 years is the approximate pooled-cohort median across Studies 16+20+29 (n=105) derived from study-level medians in Table 1 (60, 60, 13) with adults dominating the cohort; the paper does not state the exact normalization constant used in NONMEM (see vignette Assumptions).",
      source_name        = "AGE"
    ),
    ALT = list(
      description        = "Serum alanine aminotransferase activity (baseline)",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed at baseline. Power effect on selumetinib CL/F (paper text 'negative correlation' for ALT, encoded as (ALT/ALT_ref)^(-0.187) so the magnitude theta14 = +0.187 from Table 2 maps to the documented negative-correlation direction). Reference value 20 U/L is the approximate pooled-cohort median across Studies 16+20+29 (n=105) derived from study-level medians in Table 1 (17.5, 22, 22); the paper does not state the exact normalization constant used in NONMEM (see vignette Assumptions).",
      source_name        = "ALT"
    ),
    FED = list(
      description        = "Fed-state-at-dosing indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (fasted at dosing; bioavailability fixed at 1)",
      notes              = "Per-dose-record indicator. Additive effects on selumetinib absorption parameters under fed (high-fat-meal) condition: D1 increased by +4.09 h (theta9), ALAG1 increased by +0.348 h (theta10), and bioavailability decreased by an additive -0.117 (theta8) so F_fed = 1 - 0.117 = 0.883 (consistent with the paper's Figure 4 sensitivity analysis showing ~11% AUC0-Inf reduction under high-fat meal). Pediatric patients (Study 29) were modeled as fasted regardless of reported food intake per Patel 2017 Results because the reported pediatric food information was inadequate for modeling.",
      source_name        = "FED"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 105L,
    n_studies      = 3L,
    n_observations = 2418L,
    age_range      = "5.6-79 years",
    age_median     = "approximately 53 years (pooled across studies; adult-dominated)",
    weight_range   = "14-119 kg",
    sex_female_pct = 46.0,
    race_ethnicity = "Caucasian 92%, non-Caucasian 5%, missing 3% (Patel 2017 Table 1, n=105 model-development cohort)",
    disease_state  = "Pooled adults with advanced solid tumors (Study 16: non-small-cell lung cancer, n=42; Study 20: advanced solid malignancies, n=30) and children with recurrent or refractory low-grade glioma (Study 29, n=33). External validation in 44 additional pediatric LGG patients (Study 29B, NCT01089101).",
    dose_range     = "Adult: 75 mg b.i.d. oral selumetinib hydrogen-sulfate capsule (Study 16) or 75 mg single dose (Study 20); Pediatric: 25, 33, or 43 mg/m^2 b.i.d. (Study 29). External validation 25 mg/m^2 b.i.d. (Study 29B).",
    regions        = "Multinational (predominantly North America).",
    studies        = c("Study 16 (NSCLC adults)", "Study 20 (advanced solid malignancy adults)", "Study 29 (recurrent LGG children, Pediatric Brain Tumor Consortium)"),
    bsa_range_m2   = "0.64-2.46",
    notes          = "All patients received selumetinib as oral hydrogen-sulfate (Hyd-sulfate) capsule formulation. Per-study demographics are summarized in Patel 2017 Table 1. The pooled-cohort median values quoted in covariateData notes (BSA 1.66 m^2, age 53 years, ALT 20 U/L) are approximations from the study-level medians in Table 1; the paper does not state the exact NONMEM-internal normalization constants. The metabolite-to-parent AUC ratio on molar basis was approximately 0.15 in the upstream Banerji 2010 first-in-human study and Fm decreased by ~25% at steady state vs single dose in the present analysis. Inter-occasion variability (IOV) was estimated for D1, ALAG1, CL, V2, and Fm (Patel 2017 Table 2) but is not implemented in this model file because the simulation library encodes typical-value and IIV variability only; IOV is documented in vignette Assumptions and deviations."
  )

  ini({
    # Structural population parameters from Patel 2017 Table 2. Parent (selumetinib)
    # population parameters apply at the reference covariate set (BSA = 1.66 m^2,
    # AGE = 53 y, ALT = 20 U/L, FED = 0) and at single dose (the steady-state
    # reduction in Fm is documented in vignette Assumptions and deviations).
    # The Table 2 column header reads "h1 (nmol/hr)" and "h9 (nmol/hr)" but the
    # parameter description column says "Duration of zero-order drug input"; the
    # nmol/hr label is a Table 2 unit typo for these two duration entries. Both
    # values are interpreted as hours (only interpretation consistent with the
    # paper's Figure 1 absorption-model diagram and the Figure 4 ~11% food-effect
    # AUC reduction).
    ld1   <- log(0.622) ; label("Duration of zero-order drug input in the gut under fasted condition (D1, h)")  # Patel 2017 Table 2 theta1
    ltlag <- log(0.319) ; label("Absorption lag time under fasted condition (ALAG1, h)")                          # Patel 2017 Table 2 theta2
    lcl   <- log(13.5)  ; label("Apparent oral clearance of selumetinib at the reference covariate set (CL/F, L/h)") # Patel 2017 Table 2 theta3
    lvc   <- log(32.6)  ; label("Apparent selumetinib central volume of distribution at the reference covariate set (V2/F, L)") # Patel 2017 Table 2 theta4
    lvp   <- log(55.0)  ; label("Apparent selumetinib peripheral volume of distribution (V3/F, L)")               # Patel 2017 Table 2 theta5
    lq    <- log(8.2)   ; label("Apparent selumetinib inter-compartmental clearance (Q/F, L/h)")                  # Patel 2017 Table 2 theta6
    lka   <- log(3.7)   ; label("First-order absorption rate constant for selumetinib from the gut (Ka, 1/h)")    # Patel 2017 Table 2 theta7

    # Bioavailability anchor: F1 = 1 fixed under fasted condition per Patel 2017
    # Methods, "Bioavailability under the fasted condition was set to 1".
    lfdepot <- fixed(log(1.0)) ; label("Selumetinib bioavailability (F1) under fasted reference (fixed at 1)")    # Patel 2017 Methods; F1 fasted fixed at 1

    # Food-effect (FED = 1, high-fat meal) covariate coefficients applied as
    # linear additive shifts on selumetinib absorption parameters (Patel 2017
    # Methods: "food condition was evaluated as a structural model covariate on
    # the absorption parameters using a linear additive model"). The Table 2
    # estimates are magnitudes; the sign on the F effect is negative
    # (subtractive) because food decreases bioavailability (Figure 4 sensitivity
    # analysis reports ~11% AUC0-Inf reduction in adults, matching F_fed =
    # 1 - 0.117 = 0.883).
    e_fed_f       <- 0.117  ; label("Magnitude of food-effect bioavailability decrement (unitless; F_fed = 1 - e_fed_f * FED)") # Patel 2017 Table 2 theta8
    e_fed_d1      <- 4.09   ; label("Food-effect additive increase in zero-order drug input duration (h; D1_fed = D1 + e_fed_d1 * FED)") # Patel 2017 Table 2 theta9
    e_fed_tlag    <- 0.348  ; label("Food-effect additive increase in absorption lag time (h; ALAG1_fed = ALAG1 + e_fed_tlag * FED)") # Patel 2017 Table 2 theta10

    # Continuous-covariate power-form coefficients (magnitudes from Patel 2017
    # Table 2). Sign conventions per paper text:
    #   - BSA on CL/F: positive correlation -> (BSA/ref)^(+e_bsa_cl)
    #   - BSA on Vc/F: positive correlation -> (BSA/ref)^(+e_bsa_vc)
    #   - AGE on Vc/F: positive correlation -> (AGE/ref)^(+e_age_vc)
    #   - ALT on CL/F: negative correlation -> (ALT/ref)^(-e_alt_cl)
    #   - BSA on Fm  : negative correlation -> (BSA/ref)^(-e_bsa_fm)
    # Sign direction comes from Patel 2017 Results: "All significant covariates
    # showed a positive correlation with parameters except ALT" (selumetinib)
    # and "Covariate analysis identified BSA as a significant covariate on Fm
    # with a negative correlation" (metabolite).
    e_bsa_cl   <- 0.923 ; label("Magnitude of BSA power exponent on selumetinib CL/F (unitless)")                 # Patel 2017 Table 2 theta11
    e_bsa_vc   <- 1.24  ; label("Magnitude of BSA power exponent on selumetinib V2/F (unitless)")                 # Patel 2017 Table 2 theta12
    e_age_vc   <- 0.327 ; label("Magnitude of age power exponent on selumetinib V2/F (unitless)")                 # Patel 2017 Table 2 theta13
    e_alt_cl   <- 0.187 ; label("Magnitude of ALT power exponent on selumetinib CL/F (unitless; applied with negative sign)") # Patel 2017 Table 2 theta14

    # N-desmethyl-selumetinib (metabolite) population parameters from Patel 2017
    # Table 2. The metabolite central volume Vc_ndsel is fixed equal to
    # selumetinib Vc/F (V2/F) per Patel 2017 Methods, "volumes of the central
    # compartment for selumetinib and N-desmethyl-selumetinib were assumed equal
    # to avoid identifiability problem". Fm is reported as a single-dose
    # estimate that decreases by ~27.4% at steady-state (theta19), not modeled
    # here as a regime-switching effect -- see vignette Assumptions and
    # deviations. Fm > 1 is acknowledged in the paper Discussion to reflect a
    # likely V_ndsel < V_parent rather than a true >1 fractional conversion.
    lfm        <- log(1.37)  ; label("Selumetinib-to-N-desmethyl-selumetinib molar conversion coefficient at single dose (FM, unitless)") # Patel 2017 Table 2 theta15
    lcl_ndsel  <- log(240)   ; label("Apparent N-desmethyl-selumetinib clearance from central compartment at the reference covariate set (CL_Meta/F, L/h)") # Patel 2017 Table 2 theta16
    lq_ndsel   <- log(49.5)  ; label("Apparent N-desmethyl-selumetinib inter-compartmental clearance (Q_Meta/F, L/h)") # Patel 2017 Table 2 theta17
    lvp_ndsel  <- log(413)   ; label("Apparent N-desmethyl-selumetinib peripheral volume of distribution (V5/F, L)") # Patel 2017 Table 2 theta18

    # N-desmethyl-selumetinib covariate effects.
    e_bsa_fm  <- 0.908  ; label("Magnitude of BSA power exponent on Fm (unitless; applied with negative sign)") # Patel 2017 Table 2 theta20

    # Inter-individual variability. Patel 2017 Table 2 reports omega^2 (log-scale
    # variances) directly; the correlated block structure follows from the
    # reported Corr(eta_*, eta_*) entries (Patel 2017 Table 2 IIV section).
    # Selumetinib block {D1, ALAG1, V2, CL}:
    #   variances: omega^2_D1 = 0.171, omega^2_ALAG1 = 0.165, omega^2_V2 = 0.201,
    #              omega^2_CL = 0.070.
    #   covariances derived as r * sqrt(var_a * var_b):
    #     cov(D1, ALAG1)  = 0.664 * sqrt(0.171 * 0.165) = 0.1115
    #     cov(D1, V2)     = 0.820 * sqrt(0.171 * 0.201) = 0.1520
    #     cov(ALAG1, V2)  = 0.578 * sqrt(0.165 * 0.201) = 0.1053
    #     cov(CL, V2)     = 0.519 * sqrt(0.070 * 0.201) = 0.0616
    #     cov(D1, CL), cov(ALAG1, CL): not reported by Patel 2017 -> set to 0.
    etald1 + etaltlag + etalvc + etalcl ~ c(
      0.171,
      0.1115, 0.165,
      0.1520, 0.1053, 0.201,
      0.0000, 0.0000, 0.0616, 0.070
    )                                                                                                              # Patel 2017 Table 2 IIV block on {D1, ALAG1, V2, CL}
    # Selumetinib block {V3, Q}:
    #   variances: omega^2_V3 = 0.388, omega^2_Q = 0.295.
    #   covariance: cov(V3, Q) = 0.623 * sqrt(0.388 * 0.295) = 0.2108
    etalvp + etalq ~ c(0.388, 0.2108, 0.295)                                                                       # Patel 2017 Table 2 IIV block on {V3, Q}
    # Metabolite block {Fm, CL_Meta}:
    #   variances: omega^2_Fm = 0.162, omega^2_CL_Meta = 0.152.
    #   covariance: cov(Fm, CL_Meta) = 0.105 * sqrt(0.162 * 0.152) = 0.01648
    etalfm + etalcl_ndsel ~ c(0.162, 0.01648, 0.152)                                                              # Patel 2017 Table 2 IIV block on {Fm, CL_Meta}

    # Residual error. Patel 2017 Table 2 reports the combined additive +
    # proportional error variances on the linear-concentration scale (nmol/L).
    # Proportional SD = sqrt(sigma^2_prop). Additive SD = sqrt(0.63) = 0.794
    # nmol/L (FIX per paper). The model file converts the additive component to
    # ng/mL via the published selumetinib molecular weight 457.68 g/mol and the
    # N-desmethyl-selumetinib molecular weight 443.65 g/mol so the residuals
    # apply on the same ng/mL scale as the Cc / Cc_ndsel observations:
    #   addSd (ng/mL)        = 0.794 nM * 457.68 g/mol / 1000 = 0.3634 ng/mL
    #   addSd_ndsel (ng/mL)  = 0.794 nM * 443.65 g/mol / 1000 = 0.3523 ng/mL
    propSd       <- 0.3521 ; label("Proportional residual error on selumetinib (fraction; sqrt(0.124))")           # Patel 2017 Table 2 sigma^2_prop selumetinib = 0.124
    addSd        <- fixed(0.3634) ; label("Additive residual error on selumetinib (ng/mL; 0.794 nmol/L FIX * MW)") # Patel 2017 Table 2 sigma^2_add = 0.63 FIX
    propSd_ndsel <- 0.5367 ; label("Proportional residual error on N-desmethyl-selumetinib (fraction; sqrt(0.288))") # Patel 2017 Table 2 sigma^2_prop metabolite = 0.288
    addSd_ndsel  <- fixed(0.3523) ; label("Additive residual error on N-desmethyl-selumetinib (ng/mL; 0.794 nmol/L FIX * MW_ndsel)") # Patel 2017 Table 2 sigma^2_add = 0.63 FIX
  })

  model({
    # Reference covariate values used inside the power-covariate model. These
    # are pooled-cohort approximations derived from Patel 2017 Table 1 study-
    # level medians (105 patients across Studies 16+20+29) because the paper
    # does not state the exact NONMEM normalization constants used.
    ref_bsa <- 1.66
    ref_age <- 53
    ref_alt <- 20

    # Molecular weights used to convert mass-based simulated compartment
    # amounts (mg) into the molar-based 1:1 parent -> metabolite formation
    # flux. Selumetinib molecular formula C17H15BrClFN4O3, MW 457.68 g/mol
    # (PubChem CID 10127622). N-desmethyl-selumetinib molecular formula
    # C16H13BrClFN4O3, MW 443.65 g/mol (PubChem CID 11355684). The mass-
    # transfer factor (mw_ndsel / mw_parent) = 0.9694 converts the mass-based
    # parent-elimination flux into the mass-based metabolite-formation flux
    # so the molar 1:1 stoichiometry is preserved at the mass level.
    mw_parent <- 457.68
    mw_ndsel  <- 443.65

    # Individual selumetinib disposition parameters with covariate effects.
    # Power-form effects use the Patel 2017 sign convention documented in the
    # ini() block. BSA and AGE enter Vc/F with positive power exponents; BSA
    # enters CL/F with a positive exponent and ALT with a negative exponent
    # (per the paper text "All significant covariates showed a positive
    # correlation with parameters except ALT").
    cl <- exp(lcl + etalcl) *
          (BSA / ref_bsa)^e_bsa_cl *
          (ALT / ref_alt)^(-e_alt_cl)
    vc <- exp(lvc + etalvc) *
          (BSA / ref_bsa)^e_bsa_vc *
          (AGE / ref_age)^e_age_vc
    vp <- exp(lvp + etalvp)
    q  <- exp(lq  + etalq)
    ka <- exp(lka)

    # Sequential zero- and first-order absorption: under fasted (FED = 0) the
    # dose enters depot over duration D1 (h) starting at lag ALAG1 (h) and is
    # then absorbed first-order at rate Ka. Under fed (FED = 1) D1 and ALAG1
    # are increased additively and bioavailability is reduced by a subtractive
    # food coefficient (Patel 2017 Methods + Table 2). IIV for D1 and ALAG1 is
    # applied multiplicatively on the fasted value via exp(l<param> + eta);
    # the food increment is applied to the typical-value D1 and ALAG1 in
    # additive linear form because that is how the published table 2 theta_9
    # / theta_10 coefficients are reported.
    d1     <- exp(ld1   + etald1)   + e_fed_d1   * FED
    tlag   <- exp(ltlag + etaltlag) + e_fed_tlag * FED
    fdepot <- exp(lfdepot) * (1 - e_fed_f * FED)

    # Individual N-desmethyl-selumetinib disposition parameters. Fm > 1 is
    # accepted at face value as the published Table 2 theta15 single-dose
    # estimate; the steady-state reduction (Table 2 theta19, ~27.4%) is
    # documented in vignette Assumptions and deviations but is not applied as
    # a regime switch in the simulation because the paper does not specify
    # the time at which Fm transitions from the single-dose to the steady-
    # state value. BSA enters Fm with a negative power exponent per the paper
    # text "BSA as a significant covariate on Fm with a negative correlation".
    fm        <- exp(lfm + etalfm) *
                 (BSA / ref_bsa)^(-e_bsa_fm)
    cl_ndsel  <- exp(lcl_ndsel + etalcl_ndsel)
    vp_ndsel  <- exp(lvp_ndsel)
    q_ndsel   <- exp(lq_ndsel)
    vc_ndsel  <- vc

    # Micro-constants for the explicit ODE system.
    kel       <- cl / vc
    k23       <- q  / vc
    k32       <- q  / vp
    kel_ndsel <- cl_ndsel / vc_ndsel
    k45       <- q_ndsel  / vc_ndsel
    k54       <- q_ndsel  / vp_ndsel

    # ODE system (compartment amounts in mg). Selumetinib enters depot, is
    # absorbed first-order into central, distributes to peripheral1, and
    # eliminates from central at rate kel; the molar Fm fraction of the
    # parent elimination flux feeds the metabolite central compartment with
    # a mass-conserving (mw_ndsel / mw_parent) correction so that 1 mol
    # selumetinib eliminated produces Fm mol N-desmethyl-selumetinib at the
    # mass level. The metabolite then distributes between central_ndsel and
    # peripheral1_ndsel and eliminates from central_ndsel at rate kel_ndsel.
    d/dt(depot)             <- -ka * depot
    d/dt(central)           <-  ka * depot - kel * central - k23 * central + k32 * peripheral1
    d/dt(peripheral1)       <-  k23 * central - k32 * peripheral1
    d/dt(central_ndsel)     <-  fm * kel * central * (mw_ndsel / mw_parent) -
                                kel_ndsel * central_ndsel -
                                k45 * central_ndsel + k54 * peripheral1_ndsel
    d/dt(peripheral1_ndsel) <-  k45 * central_ndsel - k54 * peripheral1_ndsel

    # Sequential zero-order release + first-order absorption: rxode2 infuses
    # the dose into depot over duration `dur(depot) = d1` starting at lag
    # `alag(depot) = tlag`, after which depot drains first-order into central
    # at rate ka. Bioavailability F applies multiplicatively to the dose.
    dur(depot)  <- d1
    alag(depot) <- tlag
    f(depot)    <- fdepot

    # Plasma concentrations. Internal amount in mg, volume in L -> mg/L;
    # multiply by 1000 to report in ng/mL to align with clinical convention.
    # The Patel 2017 paper reports concentrations in nmol/L; the validation
    # vignette includes the conversion ng/mL = nmol/L * (MW / 1000).
    Cc       <- central       / vc       * 1000
    Cc_ndsel <- central_ndsel / vc_ndsel * 1000

    Cc       ~ prop(propSd)       + add(addSd)
    Cc_ndsel ~ prop(propSd_ndsel) + add(addSd_ndsel)
  })
}
