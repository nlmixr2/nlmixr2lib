Hoefman_2021_asp8232 <- function() {
  description <- paste0(
    "Mechanism-based exposure-response (PD-only) model for the vascular ",
    "adhesion protein-1 (VAP-1) inhibitor ASP8232 in adults with diabetic ",
    "kidney disease. Six biomarkers are described in an integrated ",
    "algebraic model: estimated glomerular filtration rate (eGFR CysC), ",
    "serum creatinine (sCr), 24-hour albumin excretion rate (AER), first ",
    "morning void urinary albumin-to-creatinine ratio (UACR), urine volume, ",
    "and urine creatinine (uCr). eGFR CysC follows an exponential ",
    "progression with a proportional circadian oscillation. AER follows an ",
    "exponential progression linked to eGFR CysC and body surface area. ",
    "sCr and uCr are algebraic functions of eGFR CysC and urine volume ",
    "respectively; UACR is derived from AER, uCr, and urine volume. Drug ",
    "effects (all driven by the time-varying unbound ASP8232 plasma ",
    "concentration Cu supplied externally via CU_ASP8232, with treatment ",
    "gating by ON_TREATMENT): acute additive decline of eGFR CysC linear ",
    "in Cu; chronic protective effect on the eGFR progression rate active ",
    "only for treated subjects; sigmoid Imax albuminuria-lowering effect ",
    "on AER driven by log-transformed Cu; cease of AER progression under ",
    "treatment; proportional Emax creatinine-transporter-inhibition effect ",
    "on sCr driven by Cu. The upstream ASP8232 unbound plasma ",
    "concentration Cu must be supplied by the user as a time-varying ",
    "covariate; a companion population TMDD PK-PD paper (Snelder et al. ",
    "2021, doi:10.1007/s10928-020-09717-w) provides the Cu(t) source at ",
    "steady state Cu = 125.58 nM for 40 mg qd oral ASP8232 in a typical ",
    "DKD subject."
  )
  reference <- paste(
    "Hoefman S, Snelder N, van Noort M, Garcia-Hernandez A, Onkels H,",
    "Larsson TE, Bergmann KR.",
    "Mechanism-based modeling of the effect of a novel inhibitor of",
    "vascular adhesion protein-1 on albuminuria and renal function markers",
    "in patients with diabetic kidney disease.",
    "J Pharmacokinet Pharmacodyn. 2021;48(1):21-38.",
    "doi:10.1007/s10928-020-09716-x.",
    "Companion PK-VAP-1 model:",
    "Snelder et al. J Pharmacokinet Pharmacodyn. 2021;48(1):39-53,",
    "doi:10.1007/s10928-020-09717-w."
  )
  vignette <- "Hoefman_2021_asp8232"
  units <- list(
    time          = "hour (time after first dose)",
    dosing        = "n/a (PD-only; Cu supplied externally via CU_ASP8232)",
    concentration = paste0(
      "eGFR CysC (mL/min/1.73m^2); sCr (uM); AER (mg/24h); ",
      "UACR (mg/g); urine volume (L/24h); uCr (mM); Cu (nM)"
    )
  )

  covariateData <- list(
    BSA = list(
      description        = "Body surface area (m^2).",
      units              = "m^2",
      type               = "continuous",
      reference_category = "2.034 (median in the ALBUM DKD cohort)",
      notes              = paste0(
        "Enters the eGFR-to-AER filtration link (Eq. 4) and modifies typical ",
        "baseline urine volume and uCr (theta_29 and theta_30 in Table 1). ",
        "Reference value 2.034 m^2 is stated in Methods as the median value ",
        "for individuals in the dataset."
      ),
      source_name        = "BSA"
    ),
    SEXF = list(
      description        = "Sex indicator (1 = female, 0 = male).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = paste0(
        "theta_28 = 0.829 gives female-vs-male typical baseline eGFR CysC ",
        "(83% of male value); theta_25 = 0.770 gives female-vs-male typical ",
        "baseline sCr and uCr (77% of male value). Effect enters ",
        "multiplicatively: TVpar_i = TVpar * theta_sex^SEXF."
      ),
      source_name        = "SEX"
    ),
    AGE = list(
      description        = "Age (years).",
      units              = "years",
      type               = "continuous",
      reference_category = "70 (rounded typical value; the paper does not report a centring value; typical simulated subject is 69 years old)",
      notes              = paste0(
        "theta_26 = -0.595 (age effect on sCr). Encoded as ",
        "TVsCr_i = TVsCr * (AGE / 70)^theta_26 (power form; the paper does ",
        "not print an explicit reference age). The negative exponent gives ",
        "the paper-stated 'sCr decreases with age' direction."
      ),
      source_name        = "AGE"
    ),
    ALB = list(
      description        = "Baseline serum albumin (g/L).",
      units              = "g/L",
      type               = "continuous",
      reference_category = "42 (typical simulated subject baseline; Methods Simulations)",
      notes              = paste0(
        "theta_27 = -3.86 (baseline albumin effect on AER). Encoded as ",
        "TVAER_i = TVAER * (ALB / 42)^theta_27 (power form; the paper does ",
        "not print an explicit reference albumin). The negative exponent ",
        "gives the paper-stated 'AER decreases with increasing baseline ",
        "serum albumin' direction."
      ),
      source_name        = "ALB"
    ),
    ON_TREATMENT = list(
      description        = paste0(
        "1 = subject received ASP8232 40 mg qd for 12 weeks (active arm of ",
        "the ALBUM Phase 2 study); 0 = subject received placebo. Gates the ",
        "chronic eGFR slope effect and the AER-progression-cease effect ",
        "(both are treatment-group effects that are non-zero only for ",
        "ASP8232-treated subjects)."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (placebo arm)",
      notes              = paste0(
        "Canonical ON_TREATMENT covariate (per-subject binary treatment-arm ",
        "indicator). In this model the active arm is ASP8232 40 mg qd oral. ",
        "Per-subject time-fixed. The exposure-driven effects on eGFR CysC ",
        "(acute), AER (sigmoid Imax), and sCr (Emax) do NOT require ",
        "ON_TREATMENT because they collapse to zero when Cu = 0; ",
        "ON_TREATMENT is only used to switch on the two treatment-group ",
        "effects (chronic eGFR slope and AER progression cease)."
      ),
      source_name        = "ON_TREATMENT"
    ),
    CU_ASP8232 = list(
      description        = paste0(
        "Time-varying unbound ASP8232 plasma concentration (Cu; nM) as ",
        "computed by the companion population TMDD PK-PD model (Snelder et ",
        "al. 2021, doi:10.1007/s10928-020-09717-w). Supplied per row of the ",
        "event table by the user; the model file itself does not embed a PK ",
        "compartment for ASP8232."
      ),
      units              = "nM",
      type               = "continuous",
      reference_category = "0 (no drug); 125.58 nM steady-state value for 40 mg qd oral ASP8232 in a typical DKD subject (Methods Simulations)",
      notes              = paste0(
        "Time-varying per record. Enters three drug-effect terms: ",
        "(a) acute eGFR decline as h4 * Cu additive on eGFR CysC ",
        "(theta_4 = 0.00218; RSE 63%); ",
        "(b) sigmoid Imax albuminuria-lowering effect on AER driven by ",
        "log-transformed Cu, logCu = log(1 + 1000*Cu) with Cu in nM ",
        "(theta_9 = 95.9, theta_23 = 5.95, theta_24 = 10 [FIXED Hill]); ",
        "(c) proportional Emax creatinine-transporter-inhibition effect on ",
        "sCr, Cu / (theta_31 + Cu) with theta_31 = 52.9 nM [FIXED from a ",
        "separate multi-trial PopPK model, data on file per Discussion]. ",
        "For placebo subjects set CU_ASP8232 = 0 at every record."
      ),
      source_name        = "Cu"
    ),
    TCLOCK = list(
      description        = paste0(
        "Wall-clock time of day (hours; 0-24) at each observation, used to ",
        "compute the eGFR CysC circadian rhythm (Eq. 3). Time-varying per ",
        "record."
      ),
      units              = "hour of day",
      type               = "continuous",
      reference_category = "n/a (enters via a 24-hour periodic cosine function)",
      notes              = paste0(
        "Circadian factor: (1 + theta_32 * cos(2*pi * (TCLOCK + 24 - ",
        "theta_33) / 24)) with theta_32 = 0.0783 (amplitude) and ",
        "theta_33 = 10.5 h (time of maximum). For a subject sampled at a ",
        "consistent time of day (typical for the ALBUM design, morning ",
        "pre-dose), TCLOCK is approximately constant across observations. ",
        "If TCLOCK is unknown for a record, a sensible default is TCLOCK = ",
        "10.5 (matches theta_33; sets circadian factor to (1 + theta_32) = ",
        "1.0783 = peak of the circadian wave). Setting TCLOCK = 10.5 - 12 ",
        "= -1.5 (or +22.5) selects the trough."
      ),
      source_name        = "TCLOCK"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 120L,
    n_studies      = 1L,
    age_range      = "not fully tabulated in the paper; typical simulated subject 69 years old",
    age_median     = "not reported",
    weight_range   = "not reported",
    weight_median  = "not reported (BSA median 2.034 m^2 reported in Methods)",
    sex_female_pct = NA_real_,
    race_ethnicity = NULL,
    disease_state  = paste0(
      "Type 2 diabetes with chronic kidney disease and residual albuminuria ",
      "despite standard-of-care angiotensin-converting-enzyme inhibitor ",
      "(ACEi) or angiotensin-receptor blocker (ARB) therapy. Phase 2 ALBUM ",
      "trial (ClinicalTrials.gov NCT02358096)."
    ),
    dose_range     = "40 mg qd oral ASP8232 for 12 weeks (or matched placebo)",
    regions        = "not reported (multi-site parallel-group Phase 2)",
    notes          = paste0(
      "60 ASP8232-treated + 60 placebo DKD patients. Data pooled from the ",
      "12-week ALBUM Phase 2 trial and its 24-week follow-up. Typical ",
      "simulated subject (Methods Simulations): 69-year-old male, BSA 2.034 ",
      "m^2, serum albumin 42 g/L, baseline eGFR CysC 37.1 mL/min/1.73m^2, ",
      "baseline AER 983 mg/24h, baseline sCr 146 uM. IMPORTANT: this model ",
      "represents the PD/exposure-response layer only. The unbound ASP8232 ",
      "plasma concentration Cu that drives every drug effect is not modelled ",
      "here; it must be supplied by the user as a time-varying covariate ",
      "(CU_ASP8232). The companion multi-target TMDD PK-PD paper (Snelder et ",
      "al. 2021, doi:10.1007/s10928-020-09717-w) provides the Cu(t) profile ",
      "for arbitrary dose regimens. For steady-state Cu at 40 mg qd, use ",
      "125.58 nM (Methods Simulations)."
    )
  )

  ini({
    # --- Structural PD parameters (Table 1) ---
    # All population baseline values on log scale where log-normal IIV is used.

    lTVeGFR <- log(37.1)   ; label("Baseline eGFR CysC (mL/min/1.73m^2) for males, disregarding circadian variation")  # Table 1, theta_2 = 37.1, RSE 4.2%
    eGFRt   <- 0.226       ; label("eGFR CysC progression rate (unitless; multiplied by TAFD/10000)")                  # Table 1, theta_3 = 0.226, RSE 18%
    acute   <- 0.00218     ; label("Slope of acute eGFR CysC decline (per nM Cu; additive on eGFR CysC)")              # Table 1, theta_4 = 0.00218, RSE 63%
    chronic <- 0.807       ; label("Chronic eGFR slope treatment effect (unitless; multiplied by TAFD/10000)")         # Table 1, theta_5 = 0.807, RSE 49%
    lTVAER  <- log(983)    ; label("Baseline AER (mg/24h) for typical eGFR CysC and BSA")                              # Table 1, theta_7 = 983, RSE 9.6%
    aerFilt <- 0.0196      ; label("Slope of the filtration link between eGFR CysC and AER (per unit deviation of eGFRi*BSA/2.034 from TVeGFR)") # Table 1, theta_8 = 0.0196, RSE 7.7%
    lImax   <- log(95.9)   ; label("Imax of the albuminuria-lowering AER effect (unitless; sigmoid on log-transformed Cu)")  # Table 1, theta_9 = 95.9, RSE 39%
    lTVsCr  <- log(146)    ; label("Baseline sCr (uM) for males with typical baseline eGFR CysC")                      # Table 1, theta_11 = 146, RSE 2.9%
    sCrFilt <- 0.171       ; label("Slope of the filtration link between eGFR CysC and sCr (unitless; polynomial coefficient)")  # Table 1, theta_12 = 0.171, RSE 19%
    Emax    <- 0.0753      ; label("Emax of the creatinine-transporter-inhibition effect on sCr (unitless; proportional)")  # Table 1, theta_13 = 0.0753, RSE 17%
    lTVuvol <- log(2.01)   ; label("Baseline urine volume (L/24h) for typical BSA")                                    # Table 1, theta_15 = 2.01, RSE 2.7%
    lTVuCr  <- log(5.30)   ; label("Baseline uCr (mM) for males with typical baseline urine volume and BSA")           # Table 1, theta_17 = 5.30, RSE 3.2%
    volSlp  <- 0.304       ; label("Slope of the urine-volume-to-uCr link (per L deviation of uvol_i from TVuvol)")    # Table 1, theta_18 = 0.304, RSE 13%
    uacrSc  <- 0.858       ; label("UACR scale factor (unitless; ratio of FMV UACR to 24-h collection UACR)")          # Table 1, theta_20 = 0.858, RSE 2.5%
    AerProg <- 0.430       ; label("AER progression rate (unitless; multiplied by TAFD/10000; zero for ASP8232-treated subjects)")  # Table 1, theta_21 = 0.430, RSE 38%
    lIC50   <- log(5.95)   ; label("logIC50 of the albuminuria-lowering sigmoid Imax (unitless; on log-transformed Cu scale)")  # Table 1, theta_23 = 5.95, RSE 4.7%
    Hill    <- fixed(10)   ; label("Hill coefficient of the albuminuria-lowering sigmoid Imax (unitless; FIXED)")      # Table 1, theta_24 = 10 [FIXED]
    sexSCr  <- 0.770       ; label("Sex effect on baseline sCr and uCr (female = 0.770 * male)")                        # Table 1, theta_25 = 0.770, RSE 2.4%
    ageSCr  <- -0.595      ; label("Age power exponent on baseline sCr ((AGE/70)^ageSCr; negative gives paper-stated 'sCr decreases with age')")  # Table 1, theta_26 = -0.595, RSE 17%
    albAer  <- -3.86       ; label("Baseline-albumin power exponent on baseline AER ((ALB/42)^albAer; negative gives paper-stated 'AER decreases with increasing baseline serum albumin')")  # Table 1, theta_27 = -3.86, RSE 22%
    sexEGFR <- 0.829       ; label("Sex effect on baseline eGFR CysC (female = 0.829 * male)")                         # Table 1, theta_28 = 0.829, RSE 5.3%
    bsaUvol <- 0.732       ; label("BSA power exponent on baseline urine volume ((2.034/BSA)^bsaUvol; positive gives paper-stated 'urine volume decreases with BSA')")  # Table 1, theta_29 = 0.732, RSE 41%
    bsaUCr  <- 1.21        ; label("BSA power exponent on baseline uCr ((2.034/BSA)^bsaUCr; positive gives paper-stated 'uCr decreases with BSA')")  # Table 1, theta_30 = 1.21, RSE 19%
    EC50    <- fixed(52.9) ; label("EC50 of the creatinine-transporter-inhibition Emax on sCr (nM; FIXED from a separate multi-trial PopPK model, data on file)")  # Table 1, theta_31 = 52.9 [FIXED, data on file]
    Ampli   <- 0.0783      ; label("Amplitude of the eGFR CysC circadian rhythm (fraction of baseline)")               # Table 1, theta_32 = 0.0783, RSE 46%
    tmax    <- 10.5        ; label("Time of maximum of the eGFR CysC circadian wave (hours)")                          # Table 1, theta_33 = 10.5, RSE 15%

    # --- IIV (Table 2) ---
    # Correlated baseline eGFR / baseline AER: Table 2, variances 0.0816 and 0.523, covariance -0.112 (RSEs 14%, 17%, 26%)
    etalTVeGFR + etalTVAER ~ c(0.0816, -0.112, 0.523)
    etalTVsCr    ~ 0.0133  # Table 2, x^2 = 0.0133, RSE 13%
    etalTVuvol   ~ 0.0720  # Table 2, x^2 = 0.0720, RSE 23%
    etalTVuCr    ~ 0.0479  # Table 2, x^2 = 0.0479, RSE 17%
    etaeGFRt     ~ 0.0559  # Table 2, x^2 = 0.0559, RSE 22%
    etaAmpli     ~ 1.18    # Table 2, x^2 = 1.18, RSE 57%
    etaAerProg   ~ 1.01    # Table 2, x^2 = 1.01, RSE 23%
    etaEmax      ~ 0.467   # Table 2, x^2 = 0.467, RSE 43%

    # --- Residual error (Table 2) ---
    # Proportional residuals for eGFR CysC, sCr, urine volume, uCr; additive-on-log for AER and UACR.
    propSd_eGFRobs <- 0.0987 ; label("Proportional residual SD for eGFR CysC (fraction)")                               # Table 2, r_eGFR_CysC_proportional = 0.0987, RSE 3.2%
    propSd_sCrobs  <- 0.0857 ; label("Proportional residual SD for sCr (fraction)")                                     # Table 2, r_sCr_proportional = 0.0857, RSE 3.9%
    propSd_uvolobs <- 0.165  ; label("Proportional residual SD for urine volume (fraction)")                            # Table 2, r_urine_volume_proportional = 0.165, RSE 5.2%
    propSd_uCrobs  <- 0.304  ; label("Proportional residual SD for uCr (fraction)")                                     # Table 2, r_uCr_proportional = 0.304, RSE 4.6%
    addSd_AERobs   <- 0.424  ; label("Additive residual SD for log(AER) (log units); log-normal in linear space")        # Table 2, r_log(AER)_additive = 0.424, RSE 6.1%
    addSd_UACRobs  <- 0.342  ; label("Additive residual SD for log(UACR) (log units); log-normal in linear space")       # Table 2, r_log(UACR)_additive = 0.342, RSE 4.1%
  })

  model({
    # ------------------------------------------------------------------
    # 1. Individual typical baselines with covariate effects
    # ------------------------------------------------------------------
    TVeGFR_i <- exp(lTVeGFR) * (sexEGFR ^ SEXF)
    eGFR0    <- TVeGFR_i * exp(etalTVeGFR)

    TVAER_i  <- exp(lTVAER) * ((ALB / 42) ^ albAer)
    AER0     <- TVAER_i * exp(etalTVAER)

    TVsCr_i  <- exp(lTVsCr) * (sexSCr ^ SEXF) * ((AGE / 70) ^ ageSCr)
    sCr0     <- TVsCr_i * exp(etalTVsCr)

    TVuvol_i <- exp(lTVuvol) * ((2.034 / BSA) ^ bsaUvol)
    uvol_i   <- TVuvol_i * exp(etalTVuvol)

    # Per Eq. 7 (Hoefman 2021): uCr is proportional to uvol via volSlp; the
    # baseline typical uCr additionally carries sex and BSA effects
    # (theta_25 = sexSCr; theta_30 = bsaUCr). We combine both.
    TVuCr_i  <- exp(lTVuCr) * (sexSCr ^ SEXF) * ((2.034 / BSA) ^ bsaUCr)
    uCr_i    <- TVuCr_i * (1 + volSlp * (uvol_i - TVuvol_i)) * exp(etalTVuCr)

    # ------------------------------------------------------------------
    # 2. Time and drug-exposure inputs
    # ------------------------------------------------------------------
    # `time` is time after first dose (TAFD, hours). The model treats
    # `time` as TAFD directly -- construct the event table with time=0 at
    # the first observation / dose record.
    TAFD  <- time
    Cu    <- CU_ASP8232
    logCu <- log(1 + 1000 * Cu)     # per Table 1 caption: logCu = log(1 + 1000 * Cu)

    # ------------------------------------------------------------------
    # 3. Individual eGFR CysC progression + drug effects
    # ------------------------------------------------------------------
    # Progression rate is reduced by (1 - chronic * TAFD/10000) only for
    # ASP8232-treated subjects (ON_TREATMENT = 1); see vignette Errata for
    # the functional-form interpretation.
    eGFR_rate_i <- eGFRt * exp(etaeGFRt) * (1 - ON_TREATMENT * chronic * TAFD / 10000)

    # Circadian amplitude: individual random effect on amplitude
    Ampli_i <- Ampli * exp(etaAmpli)
    circ    <- 1 + Ampli_i * cos(2 * 3.14159265358979 * (TCLOCK + 24 - tmax) / 24)

    # Progression + circadian (Eq. 2 + Eq. 3)
    eGFR_prog <- eGFR0 * exp(-eGFR_rate_i * TAFD / 10000) * circ

    # Acute Cu-linear decline (Eq. 8; additive)
    eGFR_i <- eGFR_prog - acute * Cu

    # ------------------------------------------------------------------
    # 4. AER progression + drug effects (Eq. 4 + Eq. 9)
    # ------------------------------------------------------------------
    # AER progression rate: individual eta on log scale; zero for treated
    aerProg_i <- AerProg * exp(etaAerProg) * (1 - ON_TREATMENT)

    # Filtration link (individual-adjusted): (1 + aerFilt * (eGFR_i * BSA / 2.034 - TVeGFR_i))
    aer_filt  <- 1 + aerFilt * (eGFR_i * BSA / 2.034 - TVeGFR_i)

    # AER without drug-inhibition effect
    AER_prog  <- AER0 * aer_filt * exp(aerProg_i * TAFD / 10000)

    # Sigmoid Imax albuminuria-lowering effect on AER (Eq. 9)
    Imax_i    <- exp(lImax)                   # % (theta_9); enters as fraction/100 below
    IC50_i    <- exp(lIC50)
    imax_term <- (Imax_i / 100) * (logCu ^ Hill) / ((IC50_i ^ Hill) + (logCu ^ Hill))

    AER_i     <- AER_prog * (1 - imax_term)

    # ------------------------------------------------------------------
    # 5. sCr algebraic link + drug effect (Eq. 5 + Eq. 10)
    # ------------------------------------------------------------------
    egfr_ratio <- eGFR_i / TVeGFR_i
    sCr_link   <- 1 - sCrFilt * (egfr_ratio * egfr_ratio - 1) / egfr_ratio
    sCr_prog   <- sCr0 * sCr_link

    # Creatinine-transporter-inhibition Emax on sCr (Eq. 10)
    Emax_i     <- Emax * exp(etaEmax)
    emax_term  <- Emax_i * Cu / (EC50 + Cu)

    sCr_i      <- sCr_prog * (1 + emax_term)

    # ------------------------------------------------------------------
    # 6. UACR algebraic (Eq. 1)
    # ------------------------------------------------------------------
    UACR_i <- uacrSc * AER_i / (uvol_i * uCr_i)

    # ------------------------------------------------------------------
    # 7. Outputs and residual error
    # ------------------------------------------------------------------
    eGFRobs <- eGFR_i
    AERobs  <- log(AER_i)         # additive residual on log scale
    sCrobs  <- sCr_i
    uvolobs <- uvol_i
    uCrobs  <- uCr_i
    UACRobs <- log(UACR_i)        # additive residual on log scale

    eGFRobs ~ prop(propSd_eGFRobs)
    AERobs  ~ add(addSd_AERobs)
    sCrobs  ~ prop(propSd_sCrobs)
    uvolobs ~ prop(propSd_uvolobs)
    uCrobs  ~ prop(propSd_uCrobs)
    UACRobs ~ add(addSd_UACRobs)
  })
}
