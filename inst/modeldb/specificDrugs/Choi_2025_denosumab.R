Choi_2025_denosumab <- function() {
  description <- paste(
    "Two-compartment target-mediated drug disposition (TMDD) model with quasi-steady-state (QSS)",
    "approximation and first-order subcutaneous absorption for denosumab, coupled to an indirect-response",
    "(turnover) model in which free denosumab inhibits the first-order loss rate constant of lumbar-spine",
    "bone mineral density (BMD) through a sigmoid Imax function. Fitted by Choi 2025 to pooled individual",
    "data from a Phase I single-dose study in healthy male volunteers (SB16-1001) and a Phase III study in",
    "postmenopausal women with osteoporosis (SB16-3001), pooling the SB16 biosimilar with EU- and",
    "US-sourced reference denosumab. Study population (healthy volunteer vs patient) shifts absorption,",
    "baseline RANKL and inter-compartmental clearance; body weight enters Vc, Vp and CL as power terms;",
    "race shifts CL. The treatment-group (SB16 vs reference denosumab) effect on CL was retained by the",
    "authors for the comparative biosimilarity simulation despite not being statistically significant."
  )
  reference <- paste(
    "Choi S, Park S, Jung J, Baek S, Lim H-S.",
    "Population pharmacokinetics/pharmacodynamics analysis confirming biosimilarity of SB16 to reference",
    "denosumab. Front Pharmacol. 2025;16:1631034. doi:10.3389/fphar.2025.1631034"
  )
  vignette <- "Choi_2025_denosumab"
  units <- list(time = "h", dosing = "nmol", concentration = "nmol/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters CL/F, VC/F and VP/F as a normalised power term (WT/64)^theta (Choi 2025 Table 3 and",
        "Eq 14). The reference value of 64 kg is the pooled-cohort median printed inside the Table 3",
        "covariate-model expressions; the paper's Table 2 reports a pooled median of 66.3 kg, so 64 kg",
        "is the model's own centering constant and not a re-derived median."
      ),
      source_name        = "WT"
    ),
    DIS_HEALTHY = list(
      description        = "Healthy-participant cohort indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (postmenopausal women with osteoporosis, the Phase III SB16-3001 cohort)",
      notes              = paste(
        "1 = healthy male volunteer (Phase I SB16-1001), 0 = postmenopausal patient with osteoporosis",
        "(Phase III SB16-3001). Choi 2025 report a separate typical value per cohort for ka, R0 and Q/F",
        "(Table 3); this model takes the patient level as the reference so the covariate matches the",
        "DIS_HEALTHY canonical (reference = patient) and the clinically relevant target population."
      ),
      source_name        = "study population (HV vs PMO)"
    ),
    RACE_BLACK = list(
      description        = "Black / African American race indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (Caucasian; the CL/F reference level in Choi 2025 Table 3)",
      notes              = paste(
        "Choi 2025 pooled the two studies into Asian / Black / Caucasian only (Table 2), so RACE_BLACK",
        "and RACE_ASIAN together partition the cohort with Caucasian as the reference. All 46 Black",
        "subjects came from the Phase I healthy-volunteer study."
      ),
      source_name        = "Race"
    ),
    RACE_ASIAN = list(
      description        = "Asian race indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (Caucasian; the CL/F reference level in Choi 2025 Table 3)",
      notes              = paste(
        "See RACE_BLACK. RACE_BLACK and RACE_ASIAN are mutually exclusive in this cohort; setting both",
        "to 0 gives the Caucasian reference (84.78% of the pooled population)."
      ),
      source_name        = "Race"
    ),
    TRT_SB16 = list(
      description        = "SB16 biosimilar treatment-arm indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (reference denosumab, EU- and US-sourced pooled)",
      notes              = paste(
        "1 = SB16 (proposed denosumab biosimilar), 0 = reference denosumab. Choi 2025 found treatment",
        "group was NOT a statistically significant covariate and excluded it from the covariate-selected",
        "model, but deliberately re-introduced it on CL/F for the comparative biosimilarity simulation",
        "(Methods 2.3 and Results 3.3). The implemented CL/F ratio of SB16 to reference is 0.9982, i.e.",
        "a 0.18% difference; leaving this covariate at its reference level 0 reproduces the",
        "covariate-selected final model exactly."
      ),
      source_name        = "treatment group (SB16 vs DEN)"
    )
  )

  compartmentData <- list(
    depot = list(
      analyte = "denosumab", units = "nmol",
      specimen = "administration site", verified = TRUE
    ),
    central = list(
      analyte = "denosumab", units = "nmol",
      specimen = "serum", verified = TRUE
    ),
    peripheral1 = list(
      analyte = "denosumab", units = "nmol",
      specimen = "serum", verified = TRUE
    ),
    total_target = list(
      analyte = "RANKL", units = "nmol/L",
      specimen = "serum", verified = TRUE
    ),
    BMD_LS = list(
      analyte = "lumbar spine (L1-L4) areal bone mineral density", units = "g/cm^2",
      specimen = "not applicable", verified = TRUE
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 624L,
    n_studies      = 2L,
    age_range      = "28-81 years",
    age_median     = "63 years",
    weight_range   = "47.0-94.7 kg",
    weight_median  = "66.3 kg",
    sex_female_pct = 73.08,
    race_ethnicity = c(Caucasian = 84.78, Asian = 7.69, Black = 7.37),
    disease_state  = paste(
      "Pooled healthy male volunteers (Phase I SB16-1001, n = 168) and postmenopausal women with",
      "osteoporosis (Phase III SB16-3001, n = 456)"
    ),
    dose_range     = "60 mg subcutaneously; single dose (Phase I) or at months 0, 6 and 12 (Phase III)",
    regions        = "Republic of Korea, Czech Republic, Poland, Lithuania, Denmark",
    notes          = paste(
      "Demographics from Choi 2025 Table 2 (pooled N = 624 randomised). The PK dataset comprised 6,583",
      "serum denosumab concentrations from 615 subjects; the PD dataset comprised 1,716 lumbar-spine BMD",
      "measurements from the 456 Phase III patients. Assay LLOQ was 20 ng/mL; 1,129 of 4,262 post-dose",
      "samples (26.49%) were below the limit of quantification and treated as missing. Parameters were",
      "estimated by SAEM in Monolix Suite 2024R1 using a sequential population-PK-parameters-and-data",
      "(PPP&D) approach: the PK parameters were fixed from the PK model and the PD parameters estimated",
      "from PK and PD data simultaneously."
    )
  )

  ini({
    # ---- Absorption (Choi 2025 Table 3) ---------------------------------
    # Reference level is the postmenopausal-osteoporosis patient (DIS_HEALTHY = 0).
    lka <- log(0.0078)
    label("First-order subcutaneous absorption rate constant in patients (ka, 1/h)")  # Table 3, "ka_PMO" = 0.0078 1/h [4.96% RSE]
    e_healthy_ka <- 0.5849336
    label("Effect of healthy-volunteer status on ka (log-scale)")  # Table 3: ka(HV) 0.014 / ka_PMO 0.0078 = 1.795; Results 3.1 quotes the reciprocal ratio as 0.55

    # ---- Distribution (Choi 2025 Table 3) -------------------------------
    lvc <- log(1.58)
    label("Apparent central volume of distribution at 64 kg (VC/F, L)")  # Table 3, "VC/F" = 1.58 L [4.30% RSE]
    e_wt_vc <- 1.50
    label("Power exponent on (WT/64) for VC/F (unitless)")  # Table 3, "Covariate effect (theta) of body weight on VC/F" = 1.50 [13.78% RSE]
    lvp <- log(6.06)
    label("Apparent peripheral volume of distribution at 64 kg (VP/F, L)")  # Table 3, "VP/F" = 6.06 L [1.20% RSE]
    e_wt_vp <- 0.52
    label("Power exponent on (WT/64) for VP/F (unitless)")  # Table 3, "Covariate effect (theta) of body weight on VP/F" = 0.52 [11.25% RSE]
    lq <- log(0.20)
    label("Apparent inter-compartmental clearance in patients (Q/F, L/h)")  # Table 3, "Q/F_PMO" = 0.20 L/h [20.22% RSE]
    e_healthy_q <- 1.731656
    label("Effect of healthy-volunteer status on Q/F (log-scale)")  # Table 3: Q/F(HV) 1.13 / Q/F_PMO 0.20 = 5.65; Results 3.1 quotes the reciprocal ratio as 0.18

    # ---- Linear elimination (Choi 2025 Table 3) -------------------------
    lcl <- log(0.006)
    label("Apparent linear clearance in Caucasian subjects at 64 kg (CL/F, L/h)")  # Table 3, "CL/F in Caucasian" = 0.006 L/h [1.46% RSE]
    e_wt_cl <- 0.93
    label("Power exponent on (WT/64) for CL/F (unitless)")  # Table 3, "Covariate effect (theta) of body weight on CL/F" = 0.93 [7.98% RSE]
    e_black_cl <- 0.1397619
    label("Effect of Black race on CL/F (log-scale)")  # Table 3: CL/F Black 0.0069 / Caucasian 0.006 = 1.15; Results 3.1 quotes 1.14
    e_asian_cl <- 0.2097205
    label("Effect of Asian race on CL/F (log-scale)")  # Table 3: CL/F Asian 0.0074 / Caucasian 0.006 = 1.233; Results 3.1 quotes 1.22
    e_sb16_cl <- -0.001801622
    label("Effect of SB16 biosimilar treatment on CL/F (log-scale)")  # Results 3.3: "implemented CL/F ratio of SB16 to DEN was 0.9982"; log(0.9982)

    # ---- Target (RANKL) turnover and QSS binding (Choi 2025 Table 3) ----
    lrbase_target <- log(15.23)
    label("Baseline total target (RANKL) concentration in patients (R0, nmol/L)")  # Table 3, "R0_PMO" = 15.23 nmol/L [12.78% RSE]
    e_healthy_rbase_target <- -2.74347
    label("Effect of healthy-volunteer status on baseline RANKL (log-scale)")  # Table 3: R0(HV) 0.98 / R0_PMO 15.23 = 0.0644; Results 3.1 quotes the reciprocal ratio as 15.49
    lksyn <- log(0.01)
    label("Zero-order target (RANKL) synthesis rate constant (ksyn, nmol/L/h)")  # Table 3, "ksyn" = 0.01 1/h [2.86% RSE]; see model() for the ksyn/kdeg unit note
    lkint <- log(0.022)
    label("First-order denosumab-RANKL complex internalization rate constant (kint, 1/h)")  # Table 3, "kint" = 0.022 1/h [1.71% RSE]
    lkss <- log(1.56)
    label("Quasi-steady-state equilibrium constant (KSS = (koff + kint)/kon, nmol/L)")  # Table 3, "KSS" = 1.56 nmol/L [12.41% RSE]

    # ---- Lumbar-spine BMD turnover (Choi 2025 Table 4) ------------------
    lrbase_bmd <- log(0.76)
    label("Baseline lumbar spine BMD (BMD0, g/cm^2)")  # Table 4, "BMD0" = 0.76 g/cm^2 [0.46% RSE]
    lkout <- log(0.00018)
    label("First-order lumbar spine BMD loss rate constant (kout, 1/h)")  # Table 4, "kout" = 0.00018 1/h [9.5% RSE]
    logitimax <- -1.75
    label("Logit of the maximum fractional inhibition of kout (ImaxF, unitless)")  # Table 4, "ImaxF" = -1.75 [3.02% RSE]; Eq 10 Imax = exp(ImaxF)/(1 + exp(ImaxF)) = 0.148
    lic50 <- log(6.92)
    label("Free denosumab concentration giving 50% of Imax (IC50, nmol/L)")  # Table 4, "IC50" = 6.92 nmol/L [11.65% RSE]
    lhill <- log(0.17)
    label("Hill coefficient of the inhibitory sigmoid (unitless)")  # Table 4, "HILL" = 0.17 [4.9% RSE]

    # ---- Interindividual variability ------------------------------------
    # Choi 2025 Tables 3-4 report IIV as a coefficient of variation (%CV) for
    # exponentially-distributed (log-normal) parameters (Eq 11). Variances below
    # are omega^2 = log(CV^2 + 1).
    etalcl + etalvp ~ c(
      0.06732514,
      0.01724599, 0.02389254
    )  # Table 3: CL/F 26.39% CV, VP/F 15.55% CV, correlation CORR(Vp/F, CL/F) = 0.43 [12.35% RSE]
    etalka ~ 0.277644  # Table 3, ka 56.57% CV [5.05% RSE] (shared by both study populations)
    etalvc ~ 0.322493  # Table 3, VC/F 61.69% CV [5.68% RSE]
    etalq ~ 2.278396  # Table 3, Q/F 295.99% CV [6.97% RSE] (shared by both study populations)
    etalrbase_target ~ 1.254444  # Table 3, R0 158.3% CV [7.14% RSE] (shared by both study populations)
    etalksyn ~ 0.049214  # Table 3, ksyn 22.46% CV [11.92% RSE]
    etalkint ~ 0.006190  # Table 3, kint 7.88% CV [13.46% RSE]
    etalkss ~ 0.291022  # Table 3, KSS 58.12% CV [23.07% RSE]
    etalrbase_bmd ~ 0.306162  # Table 4, BMD0 59.85% CV [3.37% RSE]
    etalic50 ~ 0.009022  # Table 4, IC50 9.52% CV [17.45% RSE]

    # ---- Residual error --------------------------------------------------
    addSd <- 0.72
    label("Additive residual error on serum denosumab concentration (nmol/L)")  # Table 3, "eps_add" = 0.72 nmol/L [4.05% RSE]; Eq 12 combined error model
    propSd <- 0.07
    label("Proportional residual error on serum denosumab concentration (fraction)")  # Table 3, "eps_prop" = 0.07 [3.74% RSE]; Eq 12 combined error model
    addSd_BMD_LS <- 0.02
    label("Additive residual error on lumbar spine BMD (g/cm^2)")  # Table 4, "eps_add" = 0.02 g/cm^2 [2.06% RSE]; Eq 13 additive error model
  })

  model({
    # 1. Individual disposition parameters.
    #    Continuous covariates enter as normalised power terms centred on the
    #    64 kg reference printed in Choi 2025 Table 3 (Eq 14); categorical
    #    covariates enter exponentially (Eq 15).
    ka <- exp(lka + e_healthy_ka * DIS_HEALTHY + etalka)
    vc <- exp(lvc + etalvc) * (WT / 64)^e_wt_vc
    vp <- exp(lvp + etalvp) * (WT / 64)^e_wt_vp
    q <- exp(lq + e_healthy_q * DIS_HEALTHY + etalq)
    cl <- exp(
      lcl + e_black_cl * RACE_BLACK + e_asian_cl * RACE_ASIAN +
        e_sb16_cl * TRT_SB16 + etalcl
    ) * (WT / 64)^e_wt_cl

    # 2. Target-turnover and QSS binding parameters.
    #    Choi 2025 Eq 7 fixes the pre-dose target steady state as
    #    R0 = ksyn / kdeg, so kdeg is derived rather than estimated. Because R0
    #    carries the study-population effect, so does kdeg. Table 3 labels ksyn
    #    with units of 1/h, but Eq 4 adds ksyn to a concentration derivative and
    #    Eq 7 divides it by a first-order rate constant, so ksyn is a zero-order
    #    synthesis rate in nmol/L/h; the printed unit is a typographical slip.
    rbase_target <- exp(lrbase_target + e_healthy_rbase_target * DIS_HEALTHY + etalrbase_target)
    ksyn <- exp(lksyn + etalksyn)
    kdeg <- ksyn / rbase_target
    kint <- exp(lkint + etalkint)
    kss <- exp(lkss + etalkss)

    # 3. PD parameters. Eq 9 fixes kin from the pre-dose BMD steady state.
    rbase_bmd <- exp(lrbase_bmd + etalrbase_bmd)
    kout <- exp(lkout)
    kin <- kout * rbase_bmd
    imax <- expit(logitimax)
    ic50 <- exp(lic50 + etalic50)
    hill <- exp(lhill)

    # 4. Initial conditions: pre-dose target and BMD turnover are at steady
    #    state and no drug-bound target is present (Choi 2025, Eqs 7 and 9).
    total_target(0) <- rbase_target
    BMD_LS(0) <- rbase_bmd

    # 5. Quasi-steady-state binding in the central compartment.
    #    ctot is the total (free + bound) denosumab concentration - the quantity
    #    the ECLIA assay reports and the state Choi 2025 Eq 2 integrates. cfree
    #    is the positive root of the QSS binding quadratic (Eq 6) and complex is
    #    the bound target (Eq 5).
    ctot <- central / vc
    qssdisc <- ctot - total_target - kss
    cfree <- 0.5 * (qssdisc + sqrt(qssdisc * qssdisc + 4 * kss * ctot))
    complex <- total_target * cfree / (kss + cfree)

    # 6. ODE system.
    #    Eq 2 is printed with the (CL + Q) * C and Q * Ap/Vp terms lacking the
    #    1/VC that dimensional consistency requires: as printed those two terms
    #    carry units of nmol/h while dCtot/dt is nmol/L/h. The equations below
    #    are written on the amount scale (Eq 2 multiplied through by VC), which
    #    restores dimensional consistency and reproduces the two-compartment QSS
    #    TMDD form of Gibiansky et al. 2008 that Choi 2025 cite for this model.
    d/dt(depot) <- -ka * depot  # Eq 1
    d/dt(central) <- ka * depot - (cl + q) * cfree +
      q * (peripheral1 / vp) - kint * complex * vc  # Eq 2 (x VC)
    d/dt(peripheral1) <- q * (cfree - peripheral1 / vp)  # Eq 3
    d/dt(total_target) <- ksyn - kdeg * total_target -
      (kint - kdeg) * complex  # Eq 4

    # 7. Indirect response: free denosumab inhibits the BMD loss rate constant
    #    through a sigmoid Imax function (Eq 8).
    d/dt(BMD_LS) <- kin -
      kout * (1 - imax * cfree^hill / (ic50^hill + cfree^hill)) * BMD_LS

    # 8. Observations. Cc is the total serum denosumab concentration; BMD_LS is
    #    itself the ODE state, so event-table observation rows use
    #    cmt = "central" for PK and cmt = "BMD_LS" for BMD.
    Cc <- ctot
    Cc ~ add(addSd) + prop(propSd)  # Eq 12
    BMD_LS ~ add(addSd_BMD_LS)  # Eq 13
  })
}
