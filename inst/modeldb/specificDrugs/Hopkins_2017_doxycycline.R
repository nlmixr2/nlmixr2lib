Hopkins_2017_doxycycline <- function() {
  description <- "Two-compartment oral population PK model for doxycycline with two transit absorption compartments, fat-free-mass allometric scaling (CL exponent 0.75, V exponent 1.0, reference 70 kg FFM), and Doryx tablet (reference) / Doryx MPC delayed-release tablet / Doryx capsule formulation effects on relative bioavailability and absorption rate, plus a food (fed-status) effect on relative bioavailability and a formulation-dependent food effect on transit rate, plus a 14.4% increase in CL for female sex. Pooled from eight phase 1 healthy-volunteer trials (n = 178)."
  reference   <- paste(
    "Hopkins AM, Wojciechowski J, Abuhelwa AY, Mudge S, Upton RN, Foster DJR (2017).",
    "Population pharmacokinetic model of doxycycline plasma concentrations using pooled study data.",
    "Antimicrobial Agents and Chemotherapy 61(5):e02401-16.",
    "doi:10.1128/AAC.02401-16."
  )
  vignette <- "Hopkins_2017_doxycycline"
  units    <- list(time = "h", dosing = "mg", concentration = "ug/L")

  covariateData <- list(
    FFM = list(
      description        = "Fat-free mass at baseline; drives allometric scaling.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Reference 70 kg FFM. Janmahasatian formula (Hopkins 2017 Methods 'General modeling strategy' paragraph 5 citing Janmahasatian 2005). Exponent 0.75 on CL and CLP1; exponent 1.0 on V and VP1.",
      source_name        = "FFM"
    ),
    SEXF = list(
      description        = "Female-sex indicator (1 = female, 0 = male).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male).",
      notes              = "Multiplicative effect on CL: female sex increases CL by 14.4% (Hopkins 2017 Table 3 COVSEX = 0.144).",
      source_name        = "Sex"
    ),
    FED = list(
      description        = "Dose-record fed-vs-fasted indicator (1 = fed, 0 = fasted).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (fasted).",
      notes              = "Reduces relative bioavailability by 10.5% irrespective of formulation (Table 3 COVFEDF = 0.105 with the paper-reported -10.5% effect direction). Adds 0.203 h to the absorption lag (Table 3 FTLAG2). Reduces KTR by 20.9% for Doryx tablet / Doryx capsule (Table 3 COVFED = -0.209) and by 54.9% for Doryx MPC (Table 3 COVFED2 = -0.549).",
      source_name        = "Fed status"
    ),
    FORM_DOX_DORYX_MPC = list(
      description        = "Doryx MPC delayed-release tablet formulation indicator (1 = Doryx MPC, 0 = Doryx tablet / Doryx capsule).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (Doryx tablet; the reference formulation in Hopkins 2017).",
      notes              = "Relative bioavailability vs Doryx tablet = 0.863 (Table 3 F1MPC). Adds 0.115 h absorption lag (Table 3 ALAG1, shared with Doryx capsule). Strengthens the food effect on KTR (-54.9% in fed state vs -20.9% for tablet / capsule).",
      source_name        = "FMPC"
    ),
    FORM_CAPSULE = list(
      description        = "Doryx (conventional-release) capsule formulation indicator (1 = Doryx capsule, 0 = Doryx tablet / Doryx MPC).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (Doryx tablet).",
      notes              = "Relative bioavailability vs Doryx tablet = 0.978 (Table 3 F1CAP). Adds 0.115 h absorption lag (Table 3 ALAG1, shared with Doryx MPC). Capsule shares the food effect on KTR with the Doryx tablet (-20.9% in fed state).",
      source_name        = "FCAP"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 178L,
    n_studies      = 8L,
    age_range      = "18-73 years",
    age_median     = "28 years (mean)",
    weight_range   = "47.2-114.4 kg",
    weight_median  = "75.9 kg (mean)",
    sex_female_pct = 32.6,
    race_ethnicity = c(`Caucasian/White` = 55.6, Other = 44.4),
    disease_state  = "Healthy volunteers (no disease state).",
    dose_range     = "Single or multiple oral doses of 75 / 100 / 120 / 150 / 200 mg as Doryx tablet, Doryx MPC delayed-release tablet, or Doryx capsule.",
    regions        = "United States",
    notes          = "Pooled eight phase 1 trials conducted by Mayne Pharma International (study IDs 10, 20, 30, 40, 50, 60, 70, 80). 651 oral doses; 7,093 plasma observations; 12.7% BLOQ (handled via M1 method). Hopkins 2017 Table 1 (demographics) and Table 6 (study designs / LLOQs). Age mean 28 yr; FFM mean 56.3 kg (range 32.3-80.1)."
  )

  ini({
    # =========================================================================
    # Structural parameters from Hopkins 2017 Table 3 'Final model', referenced
    # to a 70 kg FFM adult. Allometric scaling is applied in model() with fixed
    # exponents 0.75 (clearance terms) and 1.0 (volume terms) per Methods
    # 'General modeling strategy' paragraph 5 (Janmahasatian FFM; Anderson &
    # Holford 2008 size principles). Doryx tablet is the reference formulation
    # (F = 1, no absorption lag, no MPC-specific food effect on KTR).
    # =========================================================================
    lcl   <- log(4.63)  ; label("Apparent clearance CL/F (L/h) at 70 kg FFM")             # Hopkins 2017 Table 3 CL (Theta 4) = 4.63 L/h/70 kg FFM (%SE 5.0)
    lvc   <- log(55.2)  ; label("Apparent central volume V/F (L) at 70 kg FFM")           # Hopkins 2017 Table 3 V (Theta 3) = 55.2 L/70 kg FFM (%SE 7.8)
    lvp   <- log(49.8)  ; label("Apparent peripheral volume VP1/F (L) at 70 kg FFM")      # Hopkins 2017 Table 3 VP1 (Theta 5) = 49.8 L/70 kg FFM (%SE 3.5)
    lq    <- log(11.3)  ; label("Apparent inter-compartmental clearance CLP1/F (L/h) at 70 kg FFM") # Hopkins 2017 Table 3 CLP1 (Theta 6) = 11.3 L/h/70 kg FFM (%SE 4.1)
    lktr  <- log(2.0)   ; label("Transit absorption rate constant KTR (1/h)")             # Hopkins 2017 Table 3 KTR (Theta 7) = 2.0 1/h (%SE 4.3)

    # =========================================================================
    # Formulation effects on relative bioavailability (Doryx tablet = reference
    # with F = 1). Encoded multiplicatively via log-shifts so the typical-value
    # f(depot) for Doryx tablet remains exactly 1.
    # =========================================================================
    lf1mpc <- log(0.863); label("Doryx MPC relative bioavailability vs Doryx tablet (unitless)")     # Hopkins 2017 Table 3 F1MPC (Theta 8) = 0.863 (%SE 4.9)
    lf1cap <- log(0.978); label("Doryx capsule relative bioavailability vs Doryx tablet (unitless)") # Hopkins 2017 Table 3 F1CAP (Theta 9) = 0.978 (%SE 3.8)

    # =========================================================================
    # Absorption lag terms. ALAG1 = 0.115 h applies whenever the subject is on
    # Doryx MPC or Doryx capsule (additive on the depot lag). FTLAG2 = 0.203 h
    # applies in the fed state regardless of formulation. The two lags are
    # additive on the depot lag time so a Doryx MPC fed subject carries
    # 0.115 + 0.203 = 0.318 h of total lag (Discussion paragraph 5).
    # =========================================================================
    e_alag1   <- 0.115  ; label("Absorption lag for Doryx MPC and Doryx capsule (h)")    # Hopkins 2017 Table 3 ALAG1 (Theta 12) = 0.115 h (%SE 24.4)
    e_ftlag2  <- 0.203  ; label("Additional absorption lag in the fed state (h)")        # Hopkins 2017 Table 3 FTLAG2 (Theta 13) = 0.203 h (%SE 14.2)

    # =========================================================================
    # Covariate effects (Hopkins 2017 Table 3 'Covariate effects' block).
    # =========================================================================
    e_fed_relf    <- 0.105  ; label("Magnitude of food-induced fractional decrease in relative bioavailability") # Hopkins 2017 Table 3 COVFEDF (Theta 19) = 0.105 (%SE 33.2); applied as RELF * (1 - e_fed_relf * FED) per the paper Discussion / Abstract: "fed status was observed to decrease F by 10.5%"
    e_sex_cl      <- 0.144  ; label("Fractional effect of female sex on CL")              # Hopkins 2017 Table 3 COVSEX (Theta 20) = 0.144 (%SE 20.4); applied as CL * (1 + e_sex_cl * SEXF), giving +14.4% CL in females
    e_fed_ktr_tab <- -0.209 ; label("Fractional effect of fed status on KTR for Doryx tablet and Doryx capsule") # Hopkins 2017 Table 3 COVFED (Theta 10) = -0.209 (%SE 26.3); applied as KTR * (1 + e_fed_ktr_tab * FED) when not on Doryx MPC
    e_fed_ktr_mpc <- -0.549 ; label("Fractional effect of fed status on KTR for Doryx MPC") # Hopkins 2017 Table 3 COVFED2 (Theta 11) = -0.549 (%SE 7.7); applied as KTR * (1 + e_fed_ktr_mpc * FED) on Doryx MPC

    # =========================================================================
    # Allometric exponents (fixed at the canonical Anderson-Holford values per
    # Hopkins 2017 Methods 'General modeling strategy' paragraph 5: "fixed
    # exponents of 0.75 for clearance parameters and 1 for volumes").
    # =========================================================================
    allo_cl <- fixed(0.75) ; label("Allometric exponent on CL and CLP1 (unitless)")       # Hopkins 2017 Methods paragraph 5; canonical Anderson-Holford 2008
    allo_v  <- fixed(1.0)  ; label("Allometric exponent on V and VP1 (unitless)")         # Hopkins 2017 Methods paragraph 5; canonical Anderson-Holford 2008

    # =========================================================================
    # Correlated between-subject / population variability on the four
    # disposition parameters (Hopkins 2017 Table 3 'Between-subject
    # variability' + 'Population variability' + 'Covariance' blocks). The
    # paper reports omega values as %CV (log-normal); on-diagonal variances
    # use omega^2 = log(1 + CV^2). The off-diagonal block-3 entries reported
    # in Table 3 under 'Covariance' have magnitudes incompatible with literal
    # covariances given the variance scale (the largest reported value 0.884
    # would imply a correlation >> 1) so they are interpreted as
    # correlations; covariances are reconstructed via cov(i,j) =
    # corr(i,j) * sqrt(var(i) * var(j)).
    # =========================================================================
    etalcl + etalktr + etalvp + etalvc ~ c(
      0.036592,
      0.022222,  0.076508,
      0.010422, -0.009717, 0.022545,
      0.050080,  0.088948, -0.004862, 0.132240
    ) # Hopkins 2017 Table 3: BSV_CL 19.3% CV; BSV_KTR 28.2% CV; PPV_VP1 15.1% CV; PPV_V 37.6% CV; correlations CL-KTR 0.420, CL-VP1 0.363, CL-V 0.720, KTR-VP1 -0.234, KTR-V 0.884, VP1-V -0.089

    # =========================================================================
    # Residual error (combined proportional + additive on linear scale per
    # Hopkins 2017 Methods 'General modeling strategy' paragraph 4, with the
    # NONMEM theta-on-prop / theta-on-add / epsilon-fixed-to-1 idiom).
    # =========================================================================
    propSd <- 0.196 ; label("Proportional residual error (fraction)")                     # Hopkins 2017 Table 3 RUV_PROP (Theta 1) = 19.6% CV (%SE 2.5)
    addSd  <- 19.8  ; label("Additive residual error (ug/L)")                             # Hopkins 2017 Table 3 RUV_ADD (Theta 2) = 19.8 ug/L (%SE 18.5)
  })

  model({
    # 1. Individual parameters with FFM allometric scaling and covariate
    #    effects.
    cl  <- exp(lcl + etalcl) * (FFM / 70) ^ allo_cl * (1 + e_sex_cl * SEXF)
    vc  <- exp(lvc + etalvc) * (FFM / 70) ^ allo_v
    vp  <- exp(lvp + etalvp) * (FFM / 70) ^ allo_v
    q   <- exp(lq)           * (FFM / 70) ^ allo_cl

    # 2. Transit-absorption rate with formulation-dependent food effect.
    #    KTR is reduced by 20.9% in fed state for Doryx tablet / Doryx
    #    capsule (e_fed_ktr_tab) and by 54.9% for Doryx MPC
    #    (e_fed_ktr_mpc). Using FORM_DOX_DORYX_MPC as the selector keeps the
    #    encoding unambiguous: when FORM_DOX_DORYX_MPC = 1 only the MPC effect
    #    fires, otherwise only the tablet/capsule effect fires.
    ktr_fed_effect <- e_fed_ktr_tab * (1 - FORM_DOX_DORYX_MPC) + e_fed_ktr_mpc * FORM_DOX_DORYX_MPC
    ktr <- exp(lktr + etalktr) * (1 + ktr_fed_effect * FED)

    # 3. Two-compartment disposition with two-transit absorption chain
    #    (Hopkins 2017 Fig. 1: dose -> T1 -> T2 -> central + peripheral1).
    #    Canonical naming: depot = T1 (receives the dose), transit1 = T2,
    #    central = V, peripheral1 = VP1. Each of depot, transit1, central
    #    drains at rate KTR per Methods 'General modeling strategy'
    #    paragraph 4 ('transit compartment' models with KTR as the
    #    first-order transit rate).
    d/dt(depot)       <- -ktr * depot
    d/dt(transit1)    <-  ktr * depot - ktr * transit1
    d/dt(central)     <-  ktr * transit1 - (cl / vc) * central - (q / vc) * central + (q / vp) * peripheral1
    d/dt(peripheral1) <-  (q / vc) * central - (q / vp) * peripheral1

    # 4. Bioavailability: Doryx tablet anchored at F = 1; Doryx MPC and
    #    Doryx capsule shift F via log-multiplicative terms; food applies
    #    a -10.5% shift on top, irrespective of formulation.
    f(depot) <- exp(lf1mpc * FORM_DOX_DORYX_MPC + lf1cap * FORM_CAPSULE) * (1 - e_fed_relf * FED)

    # 5. Absorption lag on the depot compartment. ALAG1 = 0.115 h adds
    #    whenever the formulation is Doryx MPC or Doryx capsule (the two
    #    indicators are mutually exclusive). FTLAG2 = 0.203 h adds in the
    #    fed state regardless of formulation. The two lags are additive.
    alag(depot) <- e_alag1 * (FORM_DOX_DORYX_MPC + FORM_CAPSULE) + e_ftlag2 * FED

    # 6. Observation in micrograms per liter (= ng/mL).
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd) + add(addSd)
  })
}
