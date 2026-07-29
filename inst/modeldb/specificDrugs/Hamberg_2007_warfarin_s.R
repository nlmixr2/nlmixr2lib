Hamberg_2007_warfarin_s <- function() {
  description <- "S-warfarin population PK (2-compartment, first-order absorption) coupled to an inhibitory-Emax INR PD model with two parallel transit-compartment chains (6 + 1) driving the anticoagulant response (Hamberg 2007). CYP2C9 genotype and age are predictors for S-warfarin clearance; VKORC1 -1639G>A genotype is a predictor for INR sensitivity (EC50). R-warfarin is reported separately (modellib('Hamberg_2007_warfarin_r'))."
  reference <- paste(
    "Hamberg A-K, Dahl M-L, Barban M, Scordo MG, Wadelius M, Pengo V, Padrini R, Jonsson EN.",
    "A PK-PD Model for Predicting the Impact of Age, CYP2C9, and VKORC1 Genotype",
    "on Individualization of Warfarin Therapy.",
    "Clin Pharmacol Ther. 2007;81(4):529-538. doi:10.1038/sj.clpt.6100084.",
    "PMID: 17301738.",
    "PK parameters and CYP2C9/age effects from Table 2; PD parameters and VKORC1",
    "effects from Table 4; structural equations from the Appendix."
  )
  vignette <- "Hamberg_2007_warfarin_pkpd_pgx"
  paper_specific_compartments <- c(
    "coag_s1", "coag_s2", "coag_s3", "coag_s4", "coag_s5", "coag_s6",
    "coag_l1"
  )

  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    AGE = list(
      description        = "Subject age at baseline",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed at baseline. Linear scale on CL_S with reference age 71 years (cohort median); Hamberg 2007 Appendix Eq 4. CL_S decreases approximately 0.91% per year of age above 71.",
      source_name        = "AGE"
    ),
    CYP2C9_S1_COUNT = list(
      description        = "Count of CYP2C9*1 (wild-type) alleles per subject (0, 1, or 2)",
      units              = "(count, 0/1/2)",
      type               = "continuous",
      reference_category = "*1/*1 (S1_COUNT == 2) is the model's reference category; CL_S for *1/*1 is 0.314 L/h.",
      notes              = "Sum of CYP2C9_S1_COUNT + CYP2C9_S2_COUNT + CYP2C9_S3_COUNT = 2 per subject. Six diplotypes (*1/*1, *1/*2, *1/*3, *2/*2, *2/*3, *3/*3) map to the counts. Hamberg 2007 Table 2 reports per-diplotype CL_S reduction (vs *1/*1) which the model selects via diplotype indicator built from the three count covariates.",
      source_name        = "CYP2C9 genotype"
    ),
    CYP2C9_S2_COUNT = list(
      description        = "Count of CYP2C9*2 reduced-function alleles per subject (0, 1, or 2)",
      units              = "(count, 0/1/2)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "See CYP2C9_S1_COUNT. Hamberg 2007 cohort distribution (Table 1, n=150): *1/*2 17.3%, *2/*2 3.3%, *2/*3 2.7%.",
      source_name        = "CYP2C9 genotype"
    ),
    CYP2C9_S3_COUNT = list(
      description        = "Count of CYP2C9*3 reduced-function alleles per subject (0, 1, or 2)",
      units              = "(count, 0/1/2)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "See CYP2C9_S1_COUNT. Hamberg 2007 cohort distribution (Table 1, n=150): *1/*3 16.7%, *2/*3 2.7%, *3/*3 1.3%. *3 is the dominant reduced-function allele in this Italian cohort.",
      source_name        = "CYP2C9 genotype"
    ),
    VKORC1_1639G_COUNT = list(
      description        = "Count of VKORC1 -1639G alleles per subject (0, 1, or 2). The complementary -1639A count is 2 - VKORC1_1639G_COUNT.",
      units              = "(count, 0/1/2)",
      type               = "continuous",
      reference_category = "GG (G_COUNT == 2) is the EC50 reference category in the per-diplotype Hamberg parameterisation; typical EC50 = 4.61 mg/L.",
      notes              = "Hamberg 2007 cohort distribution (Table 1, Study I, n=56 with VKORC1 data): GG 39%, AG 41%, AA 20%. Per-diplotype EC50 (Table 4): GG 4.61, AG 3.02, AA 2.20 mg/L. The model uses diplotype indicators rather than a continuous per-allele dosage because the original paper reported three categorical EC50 values rather than two per-allele coefficients.",
      source_name        = "VKORC1 -1639G>A genotype"
    ),
    INR_BASE = list(
      description        = "Subject-specific baseline INR measured before warfarin administration",
      units              = "(unitless)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject (pre-medication INR). Added as an additive constant in the INR observation so the simulated INR returns to the subject's baseline when warfarin is withdrawn. Healthy untreated baseline is ~1.0; the Hamberg 2007 cohort had no centrally reported INR baseline distribution, so the default simulation value used here is 1.0.",
      source_name        = "BASE (Appendix Eq 11)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 150L,
    n_studies      = 2L,
    age_range      = "22-87 years",
    age_median     = "71 years",
    weight_range   = "45-120 kg",
    weight_median  = "80 kg",
    sex_female_pct = 34,
    race_ethnicity = "Italian (not stratified by race in source)",
    disease_state  = "Adults on long-term warfarin anticoagulant therapy for thromboembolic prophylaxis",
    dose_range     = "Single 10 mg racemic dose (Study I) and 6.25-78.75 mg/week maintenance (median 29.375 mg/week)",
    regions        = "Italy (Thrombosis Center, University of Padova)",
    cyp2c9_freq    = "*1/*1 58.7%, *1/*2 17.3%, *1/*3 16.7%, *2/*2 3.3%, *2/*3 2.7%, *3/*3 1.3% (n=150 total cohort; Hamberg 2007 Table 1)",
    vkorc1_freq    = "GG 39%, AG 41%, AA 20% (n=56 with VKORC1 typing; Study I only; Hamberg 2007 Table 1)",
    inr_steady     = "INR at steady state median 2.28 (range 1.36-4.37; Hamberg 2007 Table 1, total cohort)",
    n_pk_records   = "171 S-warfarin observations after single-dose + 150 after stable maintenance dosing (Study I + Study II)",
    n_pd_records   = "228 INR observations after single-dose + 56 at steady state (Study I only, n=57)",
    notes          = paste(
      "Pooled across Italian Studies I (n=57) and II (n=93). PK developed on the full 150-patient cohort;",
      "PD developed only on Study I (n=57) for which both single-dose plus steady-state INR and VKORC1 typing",
      "were available. See Hamberg 2007 Table 1."
    )
  )

  ini({
    # ============================================================
    # S-warfarin 2-cmt PK -- Hamberg 2007 Table 2 + Appendix
    # ============================================================
    lcls         <- log(0.314)        ; label("Apparent oral CL_S typical (L/h), CYP2C9 *1/*1, age 71 y")    # Table 2 (estimated, RSE 3.95%)
    lvc          <- log(13.8)         ; label("Apparent central volume V1 (L)")                              # Table 2 (estimated, RSE 3.64%)
    lq           <- log(0.131)        ; label("Apparent intercompartmental clearance Q (L/h)")               # Table 2 (estimated, RSE 16.3%)
    lvp          <- log(6.59)         ; label("Apparent peripheral volume V2 (L)")                           # Table 2 (estimated, RSE 24.1%)
    lka          <- fixed(log(2))     ; label("S-warfarin absorption rate Ka_S (1/h)")                       # Table 2 (fixed, sensitivity confirmed across 0.5-5; Methods)

    # Age effect: CL_S decreases ~0.91% per year of age above 71; Table 2 reports magnitude only,
    # prose ("CLS was reduced with increasing age, decreasing by approximately 9% per decade")
    # establishes the negative sign.
    e_age_cls    <- -0.00910          ; label("Age effect on CL_S (fractional change per year, ref 71 y)")   # Table 2 magnitude 0.91 (% change/year), RSE 28.9% (sign from prose: CL decreases with age)

    # CYP2C9 diplotype log-multipliers on CL_S, parameterised as fractional CL retention
    # (1 - %reduction) on the log scale so the diplotype indicators add to lcls.
    # Reference diplotype *1/*1 carries no parameter (multiplier = 1).
    e_cyp2c9_12_lcls <- log(1 - 0.315); label("CYP2C9 *1/*2 log multiplier on CL_S")                          # Table 2 (estimated, 31.5% reduction, RSE 13.8%)
    e_cyp2c9_13_lcls <- log(1 - 0.453); label("CYP2C9 *1/*3 log multiplier on CL_S")                          # Table 2 (estimated, 45.3% reduction, RSE 8.06%)
    e_cyp2c9_22_lcls <- log(1 - 0.722); label("CYP2C9 *2/*2 log multiplier on CL_S")                          # Table 2 (estimated, 72.2% reduction, RSE 8.70%)
    e_cyp2c9_23_lcls <- log(1 - 0.690); label("CYP2C9 *2/*3 log multiplier on CL_S")                          # Table 2 (estimated, 69.0% reduction, RSE 4.42%)
    e_cyp2c9_33_lcls <- log(1 - 0.852); label("CYP2C9 *3/*3 log multiplier on CL_S")                          # Table 2 (estimated, 85.2% reduction, RSE 0.91%)

    # ============================================================
    # INR PD -- Hamberg 2007 Table 4 + Appendix
    # ============================================================
    lemax        <- fixed(log(1))     ; label("Maximum anticoagulant effect Emax (fixed at 1 = 100% inhibition)")  # Table 4 (fixed; Methods)
    lhill        <- log(0.424)        ; label("Sigmoidicity factor gamma for the inhibitory Emax (unitless)")      # Table 4 (estimated, RSE 12.4%); paper symbol gamma
    lmtt1        <- log(11.6)         ; label("Mean transit time MTT1 (h) -- six-compartment 'long' chain (per-compartment)")  # Table 4 (estimated, RSE 4.65%)
    lmtt2        <- log(120)          ; label("Mean transit time MTT2 (h) -- single-compartment 'short' chain")              # Table 4 (estimated, RSE 23.0%)
    llambda      <- log(3.61)         ; label("INR-response exponent lambda (unitless, > 1 = supralinear)")                  # Table 4 (estimated, RSE 9.22%); paper symbol lambda
    linrmax      <- fixed(log(20))    ; label("Maximum INR increment above baseline (unitless)")                             # Appendix (fixed at 20; sensitivity confirmed across 15-25)

    # Per-VKORC1-diplotype typical EC50 values (mg/L), modelled categorically as Hamberg 2007 reported.
    lec50_gg     <- log(4.61)         ; label("Typical EC50 for VKORC1 GG (mg/L)")                            # Table 4 (estimated, RSE 41.4%)
    lec50_ga     <- log(3.02)         ; label("Typical EC50 for VKORC1 AG (mg/L)")                            # Table 4 (estimated, RSE 37.4%)
    lec50_aa     <- log(2.20)         ; label("Typical EC50 for VKORC1 AA (mg/L)")                            # Table 4 (estimated, RSE 38.6%)

    # ============================================================
    # IIV blocks
    #   PK: omega(CL_S) 31.0%, omega(V1) 26.2%, omega(V2) 99.1%, with covariance term between CL_S and V1.
    #     omega values from Hamberg 2007 Table 2 (CL_S RSE 18.6%, V1 RSE 32.9%, V2 RSE 54.5%).
    #     Variances on the log scale via log(1 + CV^2). Off-diagonal covariance is not numerically
    #     reported in Table 2; the starting value below preserves the structural form the authors used
    #     and will be re-estimated by the user.
    #   PD: omega(MTT1) 14.1%, omega(MTT2) 102%, omega(EC50) 40.9%, with covariance terms between all three.
    #     omega values from Hamberg 2007 Table 4 (MTT1 RSE 55.0%, MTT2 RSE 25.1%, EC50 RSE 19.1%).
    #     Off-diagonal covariances not numerically reported in Table 4; starting values preserve structure.
    # ============================================================
    etalcls + etalvc ~ c(0.0917, 0.025, 0.0664)
    etalvp   ~ 0.6845
    etalmtt1 + etalmtt2 + etalec50 ~ c(0.0197, 0.025, 0.7253, 0.020, 0.050, 0.1551)

    # ============================================================
    # Residual error
    # PK: Hamberg 2007 Table 2 reports two log-additive (~ constant-CV) residuals -- single-dose
    # (0.0908, RSE 16.4%) and steady-state (0.301, RSE 9.93%). The single registry parameter
    # below carries the LARGER (more conservative, steady-state) value because the rxode2 reserved
    # variable SS denotes steady-state DOSING (not observation phase) and cannot be used to switch
    # residual magnitudes per observation. See vignette Errata for the deviation note.
    # PD: single log-additive residual on INR.
    # ============================================================
    propSd     <- 0.301              ; label("S-warfarin proportional residual SD (fraction; Hamberg 2007 steady-state value 0.301; single-dose was 0.0908)")    # Table 2 (estimated, RSE 9.93%)
    propSd_INR <- 0.0325             ; label("INR proportional residual SD (fraction)")                              # Table 4 (estimated, RSE 11.2%)
  })

  model({
    # ---- Diplotype indicators built from CYP2C9 per-allele counts ----
    # Each indicator is 1 only for the named diplotype; *1/*1 carries no indicator (reference).
    is_cyp2c9_12 <- (CYP2C9_S1_COUNT == 1) * (CYP2C9_S2_COUNT == 1) * (CYP2C9_S3_COUNT == 0)
    is_cyp2c9_13 <- (CYP2C9_S1_COUNT == 1) * (CYP2C9_S2_COUNT == 0) * (CYP2C9_S3_COUNT == 1)
    is_cyp2c9_22 <- (CYP2C9_S2_COUNT == 2)
    is_cyp2c9_23 <- (CYP2C9_S1_COUNT == 0) * (CYP2C9_S2_COUNT == 1) * (CYP2C9_S3_COUNT == 1)
    is_cyp2c9_33 <- (CYP2C9_S3_COUNT == 2)

    cyp2c9_lfactor <- e_cyp2c9_12_lcls * is_cyp2c9_12 +
                      e_cyp2c9_13_lcls * is_cyp2c9_13 +
                      e_cyp2c9_22_lcls * is_cyp2c9_22 +
                      e_cyp2c9_23_lcls * is_cyp2c9_23 +
                      e_cyp2c9_33_lcls * is_cyp2c9_33

    # ---- VKORC1 diplotype EC50 selection ----
    is_vkorc1_gg <- (VKORC1_1639G_COUNT == 2)
    is_vkorc1_ga <- (VKORC1_1639G_COUNT == 1)
    is_vkorc1_aa <- (VKORC1_1639G_COUNT == 0)

    ec50_geno_log <- lec50_gg * is_vkorc1_gg +
                     lec50_ga * is_vkorc1_ga +
                     lec50_aa * is_vkorc1_aa

    # ---- Individual structural parameters ----
    cls   <- exp(lcls + cyp2c9_lfactor + etalcls) * (1 + e_age_cls * (AGE - 71))
    vc    <- exp(lvc  + etalvc)
    vp    <- exp(lvp  + etalvp)
    q     <- exp(lq)
    ka    <- exp(lka)

    emax  <- exp(lemax)
    hill  <- exp(lhill)
    mtt1  <- exp(lmtt1 + etalmtt1)
    mtt2  <- exp(lmtt2 + etalmtt2)
    lambda  <- exp(llambda)
    inrmax  <- exp(linrmax)
    ec50  <- exp(ec50_geno_log + etalec50)

    # ---- Micro-constants ----
    kel <- cls / vc
    k12 <- q   / vc
    k21 <- q   / vp

    # ---- S-warfarin 2-compartment PK ----
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                                k12 * central - k21 * peripheral1

    # ---- S-warfarin central concentration drives inhibitory Emax on coagulation-factor synthesis ----
    Cc          <- central / vc
    inhibition  <- emax * Cc^hill / (ec50^hill + Cc^hill)

    # ---- Two parallel transit-compartment chains (Appendix Eq 8) ----
    # "Long" chain: 6 compartments, MTT1 = 11.6 h per compartment (paper formula ktr1 = 1/MTT1);
    # represents factor VII-like fast coagulation factors. Naming follows the registry convention
    # of "_s" for the short per-compartment MTT chain.
    ktr1 <- 1 / mtt1
    d/dt(coag_s1) <- ktr1 * (1 - inhibition - coag_s1)
    d/dt(coag_s2) <- ktr1 * (coag_s1 - coag_s2)
    d/dt(coag_s3) <- ktr1 * (coag_s2 - coag_s3)
    d/dt(coag_s4) <- ktr1 * (coag_s3 - coag_s4)
    d/dt(coag_s5) <- ktr1 * (coag_s4 - coag_s5)
    d/dt(coag_s6) <- ktr1 * (coag_s5 - coag_s6)

    # "Short" chain: 1 compartment, MTT2 = 120 h (paper formula ktr2 = 1/MTT2);
    # represents factor II-like slow coagulation factors. Naming "_l" for the long-MTT chain.
    ktr2 <- 1 / mtt2
    d/dt(coag_l1) <- ktr2 * (1 - inhibition - coag_l1)

    # Baseline (no-drug) chain endpoints are 1 (unit synthesis at steady state); set all chain
    # initial conditions to 1 per Appendix ("The initial conditions for all partial derivatives ... were set to unity").
    coag_s1(0) <- 1
    coag_s2(0) <- 1
    coag_s3(0) <- 1
    coag_s4(0) <- 1
    coag_s5(0) <- 1
    coag_s6(0) <- 1
    coag_l1(0) <- 1

    # ---- INR observation (Appendix Eq 10) ----
    # The paper combines the two chain endpoints as a PRODUCT (not an arithmetic mean):
    #   log(INR_ij) = log(BASE_i + INRMAX * (1 - A6 * A7)^lambda) + eps_INR.
    # At baseline both chain endpoints are 1, giving 1 - A6*A7 = 0 and INR = BASE; under full
    # inhibition both end-compartments drift to 0 and INR -> BASE + INRMAX.
    anticoag_frac <- 1 - coag_s6 * coag_l1
    INR <- INR_BASE + inrmax * anticoag_frac^lambda

    # ---- Residual error ----
    Cc  ~ prop(propSd)
    INR ~ prop(propSd_INR)
  })
}
