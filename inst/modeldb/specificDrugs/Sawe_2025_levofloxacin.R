Sawe_2025_levofloxacin <- function() {
  description <- paste(
    "One-compartment population PK model for oral levofloxacin in South",
    "African adults treated for rifampicin-resistant tuberculosis (RR-TB),",
    "characterising the effect of third-trimester pregnancy (Sawe 2025;",
    "n = 47 pooled from two studies, 21 pregnant, 12 with matched",
    "antepartum / postpartum profiles). Savic transit-compartment",
    "absorption (analytical form, N fixed to 20, MTT = 1.07 h) feeds",
    "first-order absorption into a one-compartment disposition model.",
    "Clearance and volume are scaled allometrically by fat-free mass",
    "(Janmahasatian formula) with fixed exponents 0.75 and 1 and the",
    "cohort-median 39.4 kg as reference. Higher serum creatinine lowers",
    "clearance via a power function (exponent -0.367) centred on the",
    "cohort median 56.2 umol/L, and third-trimester pregnancy raises",
    "clearance by a further 38.1%. Between-subject variability is",
    "retained on CL only; between-occasion variability is carried on MTT,",
    "ka and F, with the BOV on F inflated 2.35-fold on the two occasions",
    "whose preceding dose was self-reported rather than observed.",
    "Bioavailability is fixed at 1.",
    sep = " "
  )
  reference <- paste(
    "Sawe S, Tsirizani L, Court R, Gausi K, Poswa A, Badat T, Wiesner L,",
    "Loveday M, Maartens G, Conradie F, Denti P.",
    "The effect of pregnancy on the population pharmacokinetics of",
    "levofloxacin in South Africans with rifampicin-resistant",
    "tuberculosis.",
    "Antimicrob Agents Chemother. 2025 May;69(5):e01626-24.",
    "doi:10.1128/aac.01626-24.",
    "Values taken from the corrected version posted 17 April 2025, which",
    "revised the Table 2 footnote defining the reported %CV.",
    sep = " "
  )
  vignette <- "Sawe_2025_levofloxacin"
  units    <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Sawe 2025 supplementary Fig. S1 and
  # the supplementary NONMEM control stream ($MODEL: COMP=(ABS DEFDOSE),
  # COMP=(CENTRAL DEFOBSERVATION); $ERROR IPRED = A(2)/V in mg/L).
  compartmentData <- list(
    depot   = list(analyte = "levofloxacin", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "levofloxacin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    FFM = list(
      description        = "Fat-free mass (Janmahasatian formula)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying: for the pregnant participants sampled at two",
        "visits, weight (and therefore FFM) was recorded at each visit and",
        "handled as a time-varying covariate (Sawe 2025 Methods,",
        "'Population pharmacokinetic modeling'). Derived in the",
        "supplementary material as",
        "FFM = WHSmax * HT^2 * WT / (WHS50 * HT^2 + WT) with",
        "WHSmax = 37.99 and WHS50 = 35.98 for females and",
        "WHSmax = 42.92 and WHS50 = 30.93 for males, HT in metres and WT",
        "in kg (the Janmahasatian et al. formula). FFM was the best size",
        "descriptor: including it improved the fit (dOFV 12.9) and beat",
        "total body weight (dOFV 3.46), which is why the sex difference in",
        "exposure is mediated through FFM rather than an explicit sex",
        "term. Cohort median FFM 39.4 kg (range 27.3-51.2, Table 1); 39.4",
        "kg is the normalisation reference per the Table 2 footnote 'b'.",
        "Allometric exponents fixed at 0.75 on CL and 1 on V.",
        sep = " "
      ),
      source_name        = "FFM"
    ),
    CREAT = list(
      description        = "Serum creatinine",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Used as a surrogate for glomerular filtration rate; the authors",
        "deliberately avoided creatinine-clearance estimating equations",
        "because these are reported to inconsistently underestimate renal",
        "function in pregnancy (Sawe 2025 Methods). Power effect on CL,",
        "CL_i = CL * (CREAT / 56.2)^e_creat_cl, where 56.2 umol/L is the",
        "cohort median (Table 1); the exponent is negative so higher serum",
        "creatinine gives lower clearance. Cohort range 25.3-110 umol/L",
        "(Table 1). Measured within +/- 2 weeks of the PK visit, not",
        "necessarily on the visit day; the authors caution against",
        "extrapolating outside the observed range. Median 45.3 umol/L in",
        "the third trimester versus 55.4 umol/L postpartum and 58.5",
        "umol/L in never-pregnant women (Discussion).",
        sep = " "
      ),
      source_name        = "CREATININE"
    ),
    PREG = list(
      description        = "Third-trimester pregnancy status indicator: 1 = pregnant, 0 = not pregnant",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not pregnant: pools the postpartum, never-pregnant-female and male records)",
      notes              = paste(
        "Time-varying within subject, unlike the time-fixed use of PREG in",
        "cohorts that enrol pregnant and non-pregnant women in parallel:",
        "12 of the 21 pregnant participants contributed matched antepartum",
        "and postpartum profiles, so the same subject carries PREG = 1 on",
        "the third-trimester visit and PREG = 0 on the postpartum visit.",
        "Pregnancy was tested as a three-level categorical covariate",
        "(pregnant / postpartum / never pregnant), but only the pregnant",
        "level was retained: the final model applies no separate postpartum",
        "effect, and the authors note this as an advantage of the model,",
        "since the FFM and serum-creatinine covariates alone reproduce the",
        "intermediate postpartum exposures (Discussion). Multiplicative",
        "fractional effect on CL, 1 + e_preg_cl * PREG, i.e. +38.1%",
        "(dOFV 58.4, 1 df, P < 0.001). Only third-trimester data were",
        "collected, so the model says nothing about earlier gestation.",
        "An alternative model that omits the serum-creatinine effect",
        "attributes a larger 53% increase to pregnancy (Table S3).",
        sep = " "
      ),
      source_name        = "PREGNANT"
    ),
    OCC = list(
      description        = "Sampling-occasion indicator: 1 = antepartum unobserved dose, 2 = antepartum observed dose, 3 = postpartum unobserved dose, 4 = postpartum observed dose",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = paste(
        "Sawe 2025 supplementary material, 'Implementation of between",
        "occasion variability and between subject variability': an",
        "occasion is a dosing event and its associated PK samples, while a",
        "visit is a PK-sampling event. Occasion 1 is the dose taken on the",
        "day before sampling and its trough sample; occasion 2 is the dose",
        "administered at the clinic on the sampling day and the subsequent",
        "profile. Non-pregnant participants had one visit, hence occasions",
        "1 and 2; pregnant participants had two visits, hence occasions 1",
        "and 2 antepartum and 3 and 4 postpartum. Decomposed inside",
        "model() into binary indicators oc1..oc4 that select the",
        "per-occasion BOV etas on log-MTT, log-ka and log-F. The",
        "odd-numbered occasions are exactly the control stream's UNOBS = 1",
        "records (the preceding dose was self-reported, with its timing",
        "imputed from the most recent reported dosing interval), and they",
        "carry the inflated BOV on F -- see the ini() comment on",
        "etaiov_fdepot_1. A user simulating a prospective regimen in which",
        "every dose is observed should pass an even occasion number",
        "(OCC = 2 or 4) so the un-inflated BOV on F applies.",
        sep = " "
      ),
      source_name        = "OCC"
    )
  )

  # Covariates that Sawe 2025 screened but did not retain in the final
  # model, plus the demographic inputs to the FFM derivation. Documentation
  # only: none of these is referenced in model().
  covariatesDataExcluded <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      notes              = paste(
        "Cohort median 58.0 kg (range 37.0-98.0, Table 1). Tested as the",
        "allometric size descriptor but beaten by fat-free mass",
        "(dOFV 3.46 for weight versus 12.9 for FFM), so it does not appear",
        "in the final model. Still required upstream as an input to the",
        "FFM formula documented under covariateData$FFM, and time-varying",
        "for the pregnant participants sampled twice.",
        sep = " "
      )
    ),
    HT = list(
      description        = "Height",
      units              = "m",
      type               = "continuous",
      notes              = paste(
        "Cohort median 1.60 m (range 1.46-1.88, Table 1). Not a covariate",
        "in the final model; required as an input to the FFM formula,",
        "which takes height in metres.",
        sep = " "
      )
    ),
    SEXF = list(
      description        = "Biological sex: 1 = female, 0 = male",
      units              = "(binary)",
      type               = "binary",
      notes              = paste(
        "33 of 47 participants (70%) female (Table 1). Tested as a",
        "covariate and not statistically significant, so absent from the",
        "final model; the authors attribute this to the sex difference in",
        "exposure already being captured by FFM (and serum creatinine),",
        "since women carry a larger proportion of essential body fat and",
        "hence less metabolically active mass (Discussion). Still required",
        "upstream to select the sex-specific WHSmax / WHS50 constants of",
        "the FFM formula.",
        sep = " "
      )
    ),
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      notes              = "Cohort median 32 years (range 19-51, Table 1). Reported as a population descriptor; not a covariate in the final model."
    ),
    ALB = list(
      description        = "Serum albumin",
      units              = "g/L",
      type               = "continuous",
      notes              = paste(
        "Cohort median 30.5 g/L (range 17.0-40.0, Table 1). Tested on",
        "clearance and not significant. The authors reason that albumin",
        "matters less here because levofloxacin is only moderately",
        "(24-38%) protein bound, and albumin concentrations were similar",
        "between pregnant and non-pregnant participants (Discussion).",
        sep = " "
      )
    ),
    HIV_POS = list(
      description        = "HIV-positive status indicator: 1 = HIV-positive, 0 = HIV-negative",
      units              = "(binary)",
      type               = "binary",
      notes              = paste(
        "31 of 47 participants (66%) living with HIV and on antiretroviral",
        "therapy, most commonly dolutegravir-based (n = 15, 48%; Table 1).",
        "Tested and not significant on levofloxacin pharmacokinetics.",
        sep = " "
      )
    ),
    EGA = list(
      description        = "Maternal estimated gestational age at the antepartum visit",
      units              = "weeks",
      type               = "continuous",
      notes              = paste(
        "Tested within the pregnant participants only, using linear,",
        "exponential and power functions centred on the data-set median,",
        "and found not significant (Sawe 2025 Methods and Results). All",
        "antepartum sampling was in the third trimester.",
        sep = " "
      )
    ),
    STUDY_BEAT = list(
      description        = "Source-study indicator: 1 = BEAT Tuberculosis (NCT04062201), 0 = King Dinuzulu Hospital observational cohort",
      units              = "(binary)",
      type               = "binary",
      notes              = paste(
        "Study and study site were both tested as covariates and neither",
        "was significant (Results), so the two pooled studies share one",
        "set of parameters.",
        sep = " "
      )
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 47L,
    n_studies       = 2L,
    n_observations  = 320L,
    n_profiles      = "19 antepartum, 14 postpartum, 12 non-pregnant female, 14 male",
    age_range       = "19-51 years (Table 1: median 32 years)",
    weight_range    = "37.0-98.0 kg (Table 1: median 58.0 kg)",
    height_range    = "1.46-1.88 m (Table 1: median 1.60 m)",
    ffm_range       = "27.3-51.2 kg (Table 1: median 39.4 kg)",
    creat_range     = "25.3-110 umol/L (Table 1: median 56.2 umol/L)",
    alb_range       = "17.0-40.0 g/L (Table 1: median 30.5 g/L)",
    sex_female_pct  = 70.2,
    race_ethnicity  = "38 (81%) black, 9 (19%) white (Table 1)",
    n_hiv_positive  = 31L,
    n_pregnant      = 21L,
    disease_state   = paste(
      "Rifampicin-resistant tuberculosis (RR-TB) on treatment. 21 of the",
      "33 female participants were pregnant, sampled in the third",
      "trimester; 12 of those contributed matched antepartum and",
      "postpartum profiles (postpartum sampling 4-8 weeks after delivery",
      "in BEAT Tuberculosis and about 6 weeks after delivery at King",
      "Dinuzulu Hospital)."
    ),
    dose_range      = paste(
      "Levofloxacin 750 or 1000 mg orally once daily by body-weight band",
      "per WHO guidance (34-50 kg and above 50 kg). A standard breakfast",
      "was given about 1 h before the observed dose in both studies."
    ),
    regions         = "South Africa (two sites for BEAT Tuberculosis; King Dinuzulu Hospital, Durban)",
    co_medication   = paste(
      "Antiretroviral therapy in 31 participants, commonly",
      "tenofovir/lamivudine/dolutegravir; RR-TB co-treatment most",
      "commonly linezolid (40, 85%), bedaquiline (44, 94%), clofazimine",
      "(47, 100%) and delamanid (13, 28%). Concomitant medication was",
      "tested as a covariate and was not significant."
    ),
    notes           = paste(
      "Pooled from BEAT Tuberculosis (ClinicalTrials.gov NCT04062201, a",
      "phase 3 RCT with a pregnancy PK sub-study; sampling pre-dose and at",
      "2, 4, 6, 8, 10 and 24 h, or pre-dose and 2, 4, 6, 8 and 24 h in the",
      "pregnancy sub-study) and a prospective observational cohort at King",
      "Dinuzulu Hospital, Durban (sampling pre-dose and 2, 4 and 6 h).",
      "Assayed by HPLC-MS/MS with an LLOQ of 0.0781 mg/L; 6 of 320",
      "samples (1.88%) were below the LLOQ, all pre-dose, and were imputed",
      "to LLOQ/2 per Beal's M6. Fitted in NONMEM 7.5.1 with FOCE-I;",
      "parameter uncertainty from sampling importance resampling. Observed",
      "median (IQR) AUC0-24 131 (108-170) mg*h/L and Cmax 11.3",
      "(9.68-13.7) mg/L."
    )
  )

  ini({
    # Structural PK parameters. Values are the final estimates of Sawe 2025
    # Table 2, quoted to the full precision of the $THETA block of the
    # supplementary NONMEM control stream (which agrees with Table 2 to
    # every printed digit). The typical values refer to a non-pregnant
    # participant with FFM 39.4 kg and serum creatinine 56.2 umol/L
    # (Table 2 footnote 'b' and the abstract).
    lcl  <- log(6.05859);  label("Typical apparent oral clearance CL/F at FFM 39.4 kg, serum creatinine 56.2 umol/L, not pregnant (L/h)")  # Sawe 2025 Table 2: CL = 6.06 L/h (95% CI 5.47-6.53); control stream $THETA 1 = 6.05859
    lvc  <- log(85.8803);  label("Typical apparent volume of distribution V/F at FFM 39.4 kg (L)")                                          # Sawe 2025 Table 2: V = 85.9 L (95% CI 80.6-91.7); control stream $THETA 2 = 85.8803
    lka  <- log(1.59498);  label("Absorption rate constant ka from the depot to central (1/h)")                                             # Sawe 2025 Table 2: Ka = 1.59 1/h (95% CI 1.11-2.40); control stream $THETA 3 = 1.59498
    lmtt <- log(1.07437);  label("Mean transit time of the Savic transit-absorption chain (h)")                                             # Sawe 2025 Table 2: MTT = 1.07 h (95% CI 0.771-1.32); control stream $THETA 7 = 1.07437

    # Fixed structural assumptions. The number of transit compartments was
    # fixed to 20 after a sensitivity analysis, to improve model stability;
    # the transit chain significantly beat a lag time (dOFV 5.87, 1 df).
    # Bioavailability is fixed at 1 -- levofloxacin bioavailability
    # approaches 100% (Introduction) -- so CL and V are apparent oral
    # values and the BOV on F below is identifiable only as relative
    # dose-to-dose variability.
    lnn     <- fixed(log(20)); label("Number of transit compartments N in the Savic chain (unitless)") # Sawe 2025 Table 2: NN = 20 (fixed); control stream $THETA 13 = 20 FIX
    lfdepot <- fixed(log(1));  label("Bioavailability F of the depot (unitless)")                      # Sawe 2025 Table 2: F = 1 (fixed); control stream $THETA 4 = 1 FIX

    # Allometric exponents on fat-free mass, both fixed to their
    # theory-based values rather than estimated (Sawe 2025 Methods:
    # "Allometric scaling of clearance (with a fixed exponent of 0.75) and
    # volume of distribution (with a fixed exponent of 1)"). Reference
    # FFM 39.4 kg, the cohort median.
    e_ffm_cl <- fixed(0.75); label("Allometric exponent on (FFM / 39.4) for CL (unitless)") # Sawe 2025 Methods; control stream ALLMCL_FFM = (FFM/TVFFM)**0.75, TVFFM = 39.4
    e_ffm_vc <- fixed(1);    label("Allometric exponent on (FFM / 39.4) for V (unitless)")  # Sawe 2025 Methods; control stream ALLMV_FFM  = (FFM/TVFFM)

    # Covariate effects on clearance.
    # Pregnancy is a FRACTIONAL multiplier, not a log-scale effect: the
    # control stream reads "IF (PREGNANT.EQ.1) preg_CL = 1 + THETA(12)",
    # so +38.1% corresponds to a coefficient of 0.380866 applied as
    # (1 + e_preg_cl * PREG), giving a CL multiplier of 1.381.
    e_preg_cl  <-  0.380866; label("Fractional change in CL during third-trimester pregnancy (PREG = 1) relative to not pregnant (unitless)") # Sawe 2025 Table 2: effect of pregnancy on CL = +38.1% (95% CI +23.4% to +57.1%); control stream $THETA 12 = 0.380866
    # Serum creatinine enters as a power exponent on the median-normalised
    # value, CL_i = CL * (sCr_i / sCr_median)^theta_sCr (Sawe 2025 Methods
    # equation, reproduced as supplementary inline equation m001).
    e_creat_cl <- -0.366504; label("Power exponent on (CREAT / 56.2) for CL (unitless)")                                                      # Sawe 2025 Table 2: effect of serum creatinine on CL = -0.367 (95% CI -0.493 to -0.104); control stream $THETA 15 = -0.366504

    # Between-subject variability. The final model retains BSV on CL only:
    # every other $OMEGA in the control stream (BSV on V, ka, F and on the
    # abandoned peripheral-compartment parameters) is "BLOCK(1) FIX 0", and
    # between-visit variability on CL was tested and not retained.
    #
    # NOTE ON THE VARIANCE SCALE. Table 2 footnote 'd' -- the footnote the
    # 17 April 2025 correction rewrote -- defines the reported percentages
    # as %CV = sqrt(omega^2) * 100, i.e. the reported number IS omega * 100
    # and NOT the log-normal CV. So the variances below are taken directly
    # from the control stream $OMEGA block and must NOT be converted via
    # log(CV^2 + 1). Each one reproduces its printed Table 2 percentage:
    # sqrt(0.0489868) = 22.1%, sqrt(0.211291) = 46.0% (Table 2 prints
    # 45.9%, a truncation of 45.97%), sqrt(0.732947) = 85.6%,
    # sqrt(0.0562092) = 23.7%.
    etalcl ~ 0.0489868 # Sawe 2025 Table 2: BSV on CL = 22.1% (95% CI 17.1-28.3); control stream $OMEGA 1 = 0.0489868

    # Between-occasion variability on the absorption parameters and on
    # bioavailability. BOV rather than BSV carries these because absorption
    # varies within a patient with gastric pH, gastric emptying, intestinal
    # motility and food (supplementary material). The control stream writes
    # each as "$OMEGA BLOCK(1)" followed by three "SAME" blocks, i.e. one
    # variance shared across the four occasions; that is encoded here as
    # the occasion-1 slot estimated and occasions 2-4 fixed equal to it.
    etaiov_mtt_1 ~ 0.211291        # Sawe 2025 Table 2: BOV on MTT = 45.9% (95% CI 30.5-70.3); control stream $OMEGA 21 = 0.211291 (occasion 1)
    etaiov_mtt_2 ~ fixed(0.211291) # occasion 2; $OMEGA BLOCK(1) SAME
    etaiov_mtt_3 ~ fixed(0.211291) # occasion 3; $OMEGA BLOCK(1) SAME
    etaiov_mtt_4 ~ fixed(0.211291) # occasion 4; $OMEGA BLOCK(1) SAME

    etaiov_ka_1 ~ 0.732947        # Sawe 2025 Table 2: BOV on Ka = 85.6% (95% CI 61.7-119); control stream $OMEGA 17 = 0.732947 (occasion 1)
    etaiov_ka_2 ~ fixed(0.732947) # occasion 2; $OMEGA BLOCK(1) SAME
    etaiov_ka_3 ~ fixed(0.732947) # occasion 3; $OMEGA BLOCK(1) SAME
    etaiov_ka_4 ~ fixed(0.732947) # occasion 4; $OMEGA BLOCK(1) SAME

    # BOV on F is the one place the four occasions do NOT share a variance.
    # Pre-dose concentrations were more variable than expected because the
    # dose preceding them was self-reported, with its timing imputed; the
    # model therefore inflates the BOV on F for those records. The control
    # stream does this by scaling the deviate,
    # "IF(UNOBS.EQ.1) BOVBIO = E_BOVF*BOVBIO" with E_BOVF = 2.35136, and
    # UNOBS = 1 is exactly the odd occasions (each visit's first occasion
    # is the previous day's unobserved dose and its trough sample --
    # supplementary material, 'Implementation of between occasion
    # variability...'). Scaling a normal deviate by a constant multiplies
    # its variance by the square of that constant, so the inflated slots
    # take 2.35136^2 * 0.0562092 = 0.310775 (equivalently 55.7% on the
    # Table 2 %CV scale). This is an exact re-expression of the published
    # parameterisation, not an approximation. Occasion 2 holds the
    # estimated variance; occasion 4 is fixed equal to it per SAME.
    etaiov_fdepot_1 ~ fixed(0.310775) # unobserved-dose occasion: 2.35136^2 * 0.0562092; Table 2 scaling factor 2.35 (95% CI 1.68-3.28), control stream $THETA 14 = 2.35136
    etaiov_fdepot_2 ~ 0.0562092       # Sawe 2025 Table 2: BOV on F = 23.7% (95% CI 19.3-28.1); control stream $OMEGA 13 = 0.0562092 (observed clinic dose)
    etaiov_fdepot_3 ~ fixed(0.310775) # unobserved-dose occasion, postpartum visit
    etaiov_fdepot_4 ~ fixed(0.0562092) # observed clinic dose, postpartum visit; $OMEGA BLOCK(1) SAME

    # Residual unexplained variability: combined additive + proportional on
    # the linear concentration scale (mg/L), matching the control stream's
    # W = SQRT(ADD**2 + PROP**2).
    #
    # The additive term is NOT simply the 0.244 mg/L of Table 2. The
    # control stream $ERROR sets "ADD = THETA(6) + (LLOQ*0.2)"
    # unconditionally on every record, with LLOQ = 0.0781 mg/L, so the
    # additive SD actually applied is 0.243823 + 0.2 * 0.0781 = 0.259443
    # mg/L. The same +20%-of-LLOQ floor appears in the sibling paediatric
    # model Denti_2018_levofloxacin.R. The two further additive inflations
    # in that $ERROR block -- +50% of LLOQ for the imputed BLQ records and
    # a huge fixed term that discards trailing BLQ records from the fit --
    # are M6 BLQ-handling devices conditional on CENS, apply to no
    # simulated observation, and are deliberately omitted. See the vignette
    # Errata.
    addSd  <- 0.259443;  label("Additive residual error (mg/L)")               # Sawe 2025 Table 2 additive error 0.244 mg/L ($THETA 6 = 0.243823) plus the control stream's unconditional 0.2 * LLOQ = 0.01562 mg/L
    propSd <- 0.0733236; label("Proportional residual error (fraction)")       # Sawe 2025 Table 2: proportional error = 7.33% (95% CI 6.49-8.10); control stream $THETA 5 = 0.0733236
  })

  model({
    # 1. Occasion indicators (binary decomposition of the OCC column).
    # Odd occasions are the visits' first occasion, whose preceding dose was
    # self-reported; even occasions are the observed clinic doses.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)

    # 2. Per-occasion BOV etas. Each absorption parameter carries exactly
    # one of its four etas, selected by the current occasion.
    iov_mtt    <- oc1 * etaiov_mtt_1    + oc2 * etaiov_mtt_2    + oc3 * etaiov_mtt_3    + oc4 * etaiov_mtt_4
    iov_ka     <- oc1 * etaiov_ka_1     + oc2 * etaiov_ka_2     + oc3 * etaiov_ka_3     + oc4 * etaiov_ka_4
    iov_fdepot <- oc1 * etaiov_fdepot_1 + oc2 * etaiov_fdepot_2 + oc3 * etaiov_fdepot_3 + oc4 * etaiov_fdepot_4

    # 3. Covariate effects on clearance.
    # Fractional pregnancy effect: multiplier 1.381 when PREG = 1.
    preg_cl  <- 1 + e_preg_cl * PREG
    # Power effect of serum creatinine, centred on the cohort median.
    creat_cl <- (CREAT / 56.2)^e_creat_cl

    # 4. Individual structural PK parameters. Allometric scaling on
    # fat-free mass with reference 39.4 kg; BSV on CL only.
    cl  <- exp(lcl + etalcl) * (FFM / 39.4)^e_ffm_cl * preg_cl * creat_cl
    vc  <- exp(lvc)          * (FFM / 39.4)^e_ffm_vc
    ka  <- exp(lka  + iov_ka)
    mtt <- exp(lmtt + iov_mtt)
    nn  <- exp(lnn)
    # Bioavailability is fixed at 1 as a typical value but still carries
    # the between-occasion variability described above.
    fbio <- exp(lfdepot + iov_fdepot)

    # 5. ODE system. Analytical Savic transit-compartment chain: rxode2's
    # transit(n, mtt, bio) emits the gamma-density input rate
    # bio * podo * ktr * (ktr * tad)^n * exp(-ktr * tad) / n!, with
    # ktr = (n + 1) / mtt -- exactly the control stream's
    # PIZZA = LOG(BIO*PD*KTR) - GAMLN(NN+1),
    # TRANSIT = EXP(PIZZA + NN*LOG(KTR*TEMPO) - KTR*TEMPO),
    # KTR = (NN+1)/MTT. The depot then drains to central via first-order
    # ka, matching DADT(1) = TRANSIT - KA*A(1) and
    # DADT(2) = KA*A(1) - K*A(2).
    #
    # Dose to `depot`. Do NOT add `f(depot) <- 0` here. That line is the
    # natural-looking analogue of the control stream's "F1 = 0" (and is
    # what several older transit models in this package carry), but under
    # rxUi with rxode2 5.1.7 it silently zeroes the transit input as well
    # and the model then absorbs nothing at all -- every concentration
    # comes back 0. Without it, transit() is already the only input
    # pathway: the bolus is not separately added, and mass balance is
    # exact (a single 1000 mg dose with F = 1 yields AUC(0-inf) * CL =
    # 999.98 mg, and the steady-state AUC identity in the validation
    # vignette holds to under 0.1%). The vignette's AUC identity check is
    # the standing regression gate for this.
    kel <- cl / vc
    d/dt(depot)   <- transit(nn, mtt, fbio) - ka * depot
    d/dt(central) <-                          ka * depot - kel * central

    # 6. Observation and error. Dose in mg and V in L give Cc in mg/L
    # (= ug/mL), the units the paper reports throughout; the control
    # stream's IPRED = A(2)/V agrees.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
