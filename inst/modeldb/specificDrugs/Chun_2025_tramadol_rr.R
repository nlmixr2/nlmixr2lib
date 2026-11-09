Chun_2025_tramadol_rr <- function() {
  description <- "Enantiomer-specific joint parent-metabolite population PK/PD model for (1R,2R)-tramadol, its active metabolite (1R,2R)-O-desmethyltramadol (M1) and (1R,2R)-N-desmethyltramadol (M2), and VAS pain score in adults with chronic neuropathic pain after oral racemic tramadol (Chun 2025). Two-compartment unbound tramadol with Weibull time-dependent absorption (rate 1 - exp(-(ra * (tad - tlag))^gam1), per hour) and a first-pass split of the dose: a fraction 1 - fp enters the tramadol depot and a fraction fp enters a first-pass depot that delivers M1 directly. Unbound tramadol is cleared by M1 formation (fm_m1 * CLF_M1), an 'other' route (0.3 * CL) and Michaelis-Menten formation of M2 driven by the central AMOUNT. M1 is one-compartment with enterohepatic circulation (biliary clearance into a gallbladder that empties between tmeal1 and tmeal1 + dge h after each dose). Total concentrations are unbound / fu with fu = 1 / (1 + kns); diabetes (type 1 or 2) scales the binding constants of tramadol and M1 and the IC50, BMI scales the unbound central volume, and the urinary metoprolol metabolic ratio (CYP2D6 activity) scales M1 formation and M2 clearance. Pain (VAS, 0-100 mm) follows an inhibitory Imax model driven by unbound M1 through an effect compartment. Amounts in nmol of the (1R,2R) enantiomer: dose BOTH depot and depot_m1 with the full enantiomer amount (half of the racemic dose, converted with 263.37 g/mol); f() splits it."
  reference <- "Chun D, Mehta P, Guzy S, Cicali B, Lauretti GR, Lanchote VL, Vozmediano V, De Moraes N. Enhanced Sensitivity to Tramadol in Diabetic Neuropathic Pain Compared to Nerve Compression Neuropathies: A Population PK/PD Model Analysis. CPT Pharmacometrics Syst Pharmacol. 2025;14(4):781-795. doi:10.1002/psp4.13315"
  vignette <- "Chun_2025_tramadol"
  units <- list(time = "h", dosing = "nmol", concentration = "nmol/L")

  covariateData <- list(
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on the apparent unbound central volume of tramadol, (BMI / 26.045)^e_bmi_vc; 26.045 kg/m^2 is the cohort median (Chun 2025 Methods, Covariate Analysis). Cohort mean (SD) 27.26 (6.51), range 17.2-38.9 kg/m^2 (Table S2 and Results).",
      source_name = "BMI"
    ),
    CYP2D6 = list(
      description = "CYP2D6 in vivo activity expressed as the metoprolol / alpha-hydroxymetoprolol metabolic ratio in 0-8 h urine after a 100 mg oral metoprolol probe dose",
      units = "(unitless ratio)",
      type = "continuous",
      reference_category = NULL,
      notes = "INVERSE orientation: a HIGHER ratio means LOWER CYP2D6 activity, so the negative power exponents on M1 formation clearance and M2 clearance give a positive relation with CYP2D6 activity. Enters as (CYP2D6 / 0.610)^e, where 0.610 is the cohort median (Chun 2025 Methods, Covariate Analysis; NONMEM column MRURINE in Tables S6/S8). Cohort mean (SD) 1.16 (1.66) (Table S2).",
      source_name = "MRURINE"
    ),
    DIS_DIAB = list(
      description = "Diabetes mellitus indicator (type 1 or type 2), 1 = diabetes, 0 = no diabetes",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no diabetes; neuropathic pain secondary to nerve compression)",
      notes = "Used jointly with DIS_DIAB_TYPE2 to rebuild the paper's three-level DMTYPE column (0 = no diabetes, 1 = type 1, 2 = type 2): type 1 = DIS_DIAB = 1 and DIS_DIAB_TYPE2 = 0; type 2 = DIS_DIAB = 1 and DIS_DIAB_TYPE2 = 1. Multiplies the tramadol and M1 binding constants and the IC50 by group-specific factors (Chun 2025 Equation 10, Table 1).",
      source_name = "DMTYPE"
    ),
    DIS_DIAB_TYPE2 = list(
      description = "Type 2 diabetes mellitus indicator, 1 = type 2 diabetes, 0 = type 1 diabetes or no diabetes",
      units = "(binary)",
      type = "binary",
      reference_category = "0",
      notes = "Separates type 2 from type 1 diabetes among DIS_DIAB = 1 subjects (DMTYPE = 2 in the source). The model reads type 2 as DIS_DIAB * DIS_DIAB_TYPE2, so DIS_DIAB_TYPE2 = 1 with DIS_DIAB = 0 has no effect.",
      source_name = "DMTYPE"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened (Chun 2025 Methods, Covariate Analysis); not retained.",
      source_name = "AGE"
    ),
    SEXF = list(
      description = "Sex (1 = female)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Screened; not retained.",
      source_name = "GENDER"
    ),
    FPG = list(
      description = "Fasting blood glucose",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened; not retained.",
      source_name = "FBG"
    ),
    HBA1C = list(
      description = "Glycated hemoglobin",
      units = "mmol/mol",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened; not retained (missing in 12 of 30 subjects per Table S2).",
      source_name = "HbA1C"
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened; not retained.",
      source_name = "CREATININE"
    ),
    BUN = list(
      description = "Serum urea",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened as serum urea; not retained.",
      source_name = "UREA"
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened; not retained. Table S2 prints the unit as g/L but the mean of 4.3 is in g/dL.",
      source_name = "ALBUMIN"
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened; not retained.",
      source_name = "AST"
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened; not retained.",
      source_name = "ALT"
    ),
    CYP3A4 = list(
      description = "CYP3A4 in vivo activity (midazolam clearance)",
      units = "mL/min/kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened; not retained.",
      source_name = "CLMDZ"
    ),
    CYP2B6 = list(
      description = "CYP2B6 genotype",
      units = "(categorical)",
      type = "categorical",
      reference_category = NULL,
      notes = "Screened as a genotype; not retained. The coding of the source column is not reported.",
      source_name = "CYP2B6GENO"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "(1R,2R)-tramadol", units = "nmol", specimen = "administration site", verified = TRUE),
    depot_m1 = list(
      analyte = "(1R,2R)-tramadol (first-pass fraction, delivered as (1R,2R)-O-desmethyltramadol)",
      units = "nmol",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(analyte = "(1R,2R)-tramadol unbound", units = "nmol", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "(1R,2R)-tramadol unbound", units = "nmol", specimen = "plasma", verified = TRUE),
    central_m1 = list(
      analyte = "(1R,2R)-O-desmethyltramadol (M1) unbound",
      units = "nmol",
      specimen = "plasma",
      verified = TRUE
    ),
    gallbladder = list(
      analyte = "(1R,2R)-O-desmethyltramadol (M1)",
      units = "nmol",
      specimen = "bile",
      verified = TRUE
    ),
    central_m2 = list(
      analyte = "(1R,2R)-N-desmethyltramadol (M2) unbound",
      units = "nmol",
      specimen = "plasma",
      verified = TRUE
    ),
    effect = list(
      analyte = "(1R,2R)-O-desmethyltramadol (M1) unbound, effect site",
      units = "nmol/L",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 30L,
    n_studies = 1L,
    n_observations = "2145 PK concentrations (total and unbound tramadol, M1 and M2 of both enantiomers); 402 pain scores from 27 subjects",
    age_range = "19-59 years",
    age_mean = "42.4 years (SD 11.5)",
    bmi_range = "17.2-38.9 kg/m^2",
    sex_female_pct = 53.3,
    race_ethnicity = "not reported in the source paper",
    disease_state = "Chronic neuropathic pain (more than 3 months, pain score above 40/100): diabetic neuropathy with type 1 diabetes (n = 9) or type 2 diabetes (n = 9, median HbA1c 10.7%), or neuropathic pain secondary to nerve compression (cervical or lumbar disc herniation, carpal tunnel syndrome; n = 12)",
    dose_range = "Single oral dose of 100 mg racemic tramadol (Tramal capsule) after an 8 h fast",
    regions = "Brazil (Hospital of the University of Sao Paulo)",
    notes = "Demographics from Chun 2025 Table S2 and Results 'Data'. Eight vomiting patients were excluded from all analyses and three non-responders from the PD analysis. CYP2D6 activity was phenotyped with a metoprolol probe (urinary metabolic ratio)."
  )

  ini({
    # Weibull time-dependent absorption of tramadol (Chun 2025 Equation 4; NONMEM
    # Table S6 (i): KA = 1 - EXP(-((TIME - GAMMA) / BETA)**ALPHA)). Table 1 swaps the
    # 'scale' and 'location' labels of gamma and beta; the code fixes the roles:
    # gamma is the location (lag) and beta the scale.
    lgam1 <- log(2.96); label("Weibull shape parameter alpha (unitless)") # Table 1 'alpha' = 2.96 (RSE 0.9%)
    lra <- log(1 / 0.762); label("Weibull rate-scaling parameter, reciprocal of the scale beta (1/h)") # Table 1 'beta' = 0.762 h (RSE 4.2%); ra = 1/beta
    ltlag <- log(0.147); label("Weibull location parameter gamma, absorption lag (h)") # Table 1 'gamma' = 0.147 h (RSE 4.5%)

    # Unbound tramadol disposition (apparent values)
    lvc <- log(63.4); label("Apparent unbound central volume of tramadol at BMI 26.045 kg/m^2, Vc,u/F (L)") # Table 1 'Vc parent,u/F' = 63.4 L (RSE 2.3%)
    e_bmi_vc <- 3.87; label("Power exponent of BMI on Vc,u/F (unitless)") # Table 1 'beta1' = 3.87 (RSE 8.5%); NONMEM applies it to V2 (central), not Vp as the Table 1 label says
    lvp <- log(439); label("Apparent unbound peripheral volume of tramadol, Vp,u/F (L)") # Table 1 'Vp peripheral,u/F' = 439 L (RSE 1.3%)
    lq <- log(299); label("Apparent unbound intercompartmental clearance of tramadol, Q,u/F (L/h)") # Table 1 'Q parent,u/F' = 299 L/h (RSE 0.7%)
    lcl <- log(36.4); label("Apparent unbound 'other' clearance of tramadol, CLother,u/F (L/h)") # Table 1 'CLother parent,u' = 36.4 L/h (RSE 0.7%)
    fm_other <- fixed(0.3); label("Fraction applied to CLother, the renally excreted unchanged share of the dose (unitless)") # Methods: 'CLother,u/F was fixed to 0.3'; Table S6 (iii) fOTHER = 0.3

    # Plasma protein binding of tramadol: Ctotal = Cunbound * (1 + Kb) (Equations 1-3)
    lkns <- log(0.718); label("Linear binding proportionality constant of tramadol, Kb = Cbound/Cunbound, no diabetes (unitless)") # Table 1 'Kb parent' = 0.718 (RSE 3.2%); fu = 0.582
    e_t1dm_kns <- 0.311; label("Multiplicative factor on tramadol Kb for type 1 diabetes (unitless)") # Table 1 'beta2' = 0.311 (RSE 24.6%)
    e_t2dm_kns <- 0.95; label("Multiplicative factor on tramadol Kb for type 2 diabetes (unitless)") # Table 1 'beta3' = 0.95 (RSE 1.0%)

    # M1 formation and first-pass input
    fm_m1 <- 0.538; label("Fraction of tramadol converted to M1, multiplier on CLF,M1 (unitless)") # Table 1 'fm1' = 0.538 (RSE 0.8%)
    lcl_form_m1 <- log(58.8); label("Apparent unbound formation clearance of M1 at CYP2D6 ratio 0.610, CLF,M1,u/F (L/h)") # Table 1 'CLF M1,u/F' = 58.8 L/h (RSE 1.6%)
    e_cyp2d6_cl_form_m1 <- -0.523; label("Power exponent of the urinary metoprolol metabolic ratio on CLF,M1,u/F (unitless)") # Table 1 'beta4' = -0.523 (RSE 25.2%)
    lfprat_m1 <- log(0.0661); label("Ratio of first-pass M1 input to bioavailable tramadol, FRAC; fp = FRAC/(1 + FRAC) (unitless)") # Table 1 'FRAC' = 0.0661 (RSE 2.7%); fp = 6.2%
    lka_m1 <- log(0.796); label("First-order rate of first-pass M1 input, Kmet1 (1/h)") # Table 1 'Kmet1' = 0.796 1/h (RSE 1.0%)

    # M1 disposition. Total M1 volume is scaled from rat with exponent 1 (Methods):
    # 0.1 L * (72.2 kg / 0.235 kg) = 30.7234 L; the unbound volume is that / fu_m1.
    lvc_m1 <- fixed(log(0.1 * 72.2 / 0.235)); label("Apparent total central volume of M1, allometrically scaled from rat (L)") # Methods: Vc R-M1,t,rat = 0.1 L, WT human 72.2 kg, WT rat 0.235 kg; Table S6 V6 = 30.7234/FUN2
    lcl_m1 <- log(42); label("Apparent unbound clearance of M1, CLM1,u/F (L/h)") # Table 1 'CLM1,u/F' = 42 L/h (RSE 0.6%)
    lkns_m1 <- log(3.07); label("Linear binding proportionality constant of M1, Kb,M1, no diabetes (unitless)") # Table 1 'Kb M1' = 3.07 (RSE 17.2%); fu = 0.246
    e_t1dm_kns_m1 <- 0.105; label("Multiplicative factor on M1 Kb for type 1 diabetes (unitless)") # Table 1 'beta5' = 0.105 (RSE 6.9%)
    e_t2dm_kns_m1 <- 0.821; label("Multiplicative factor on M1 Kb for type 2 diabetes (unitless)") # Table 1 'beta6' = 0.821 (RSE 7.5%)

    # Enterohepatic circulation of M1
    lcl_bile_m1 <- log(382); label("Apparent unbound biliary clearance of M1 into the gallbladder, CLbile,u/F (L/h)") # Table 1 'CLbile u/F' = 382 L/h (RSE 2.9%)
    lkehc_m1 <- log(0.0019); label("Gallbladder release rate constant of M1, Krel,u (1/h)") # Table 1 'Krel u' = 0.0019 1/h (RSE 4.9%)
    tmeal1 <- 2.85; label("Start of gallbladder release after each dose, STARTTIME (h)") # Table 1 'STARTTIME' = 2.85 h (RSE 1.3%)
    dge <- 6.15; label("Duration of gallbladder release, DURATION (h)") # Table 1 'DURATION' = 6.15 h (RSE 3.2%)

    # M2 formation (Michaelis-Menten on the central unbound AMOUNT) and disposition
    lvmax_m2 <- log(98.6); label("Maximum rate of M2 formation, Vmax (nmol/h)") # Table 1 'Vmax' = 98.6 nmol/h (RSE 0.3%)
    lkm_m2 <- log(16300); label("Michaelis-Menten constant of M2 formation, in nmol of central amount (nmol)") # Table 1 'Km' = 16,300 (RSE 1.2%); printed as nmol/L but NONMEM Table S6 (iii) divides by A(3), an amount
    lvc_m2 <- log(3.13); label("Apparent unbound volume of M2, VM2,u/F (L)") # Table 1 'VM2,u/F' = 3.13 L (RSE 0.7%)
    lcl_m2 <- log(7.75); label("Apparent unbound clearance of M2 at CYP2D6 ratio 0.610, CLM2,u/F (L/h)") # Table 1 'CLM2,u/F' = 7.75 L/h (RSE 1.9%)
    e_cyp2d6_cl_m2 <- -0.477; label("Power exponent of the urinary metoprolol metabolic ratio on CLM2,u/F (unitless)") # Table 1 'beta7' = -0.477 (RSE 10.1%)

    # Pain score (VAS fraction) inhibitory Imax model with an effect compartment
    lke0 <- log(0.0398); label("Effect-compartment equilibration rate constant, ke0 (1/h)") # Table 1 'ke0' = 0.0398 1/h (RSE 20.3%)
    logitimax <- qlogis(0.908); label("Logit of the maximum fractional pain reduction, Imax (unitless)") # Table 1 'Imax' = 0.908 (RSE 5.2%)
    logite0 <- qlogis(0.464); label("Logit of the baseline pain score as a fraction of the 100 mm VAS, PD0 (unitless)") # Table 1 'PD0' = 0.464 (RSE 3.6%)
    lic50 <- log(2.36); label("Unbound M1 concentration at half-maximal pain reduction, no diabetes, IC50 (nmol/L)") # Table 1 'IC50 Non-diabetic group' = 2.36 nmol/L (RSE 9.3%)
    e_t1dm_ic50 <- 0.229; label("Multiplicative factor on IC50 for type 1 diabetes (unitless)") # Table 1 'beta8' = 0.229 (RSE 12.3%)
    e_t2dm_ic50 <- 0.431; label("Multiplicative factor on IC50 for type 2 diabetes (unitless)") # Table 1 'beta9' = 0.431 (RSE 10.2%)

    # IIV (variances). Exponential on the PK parameters; additive on the logit
    # scale for Imax and PD0 (Equation 8). The eta on beta is carried on ra = 1/beta
    # (same variance, opposite sign of eta).
    etalgam1 ~ 0.1 # Table 1 'omega2 alpha' = 0.1 (RSE 29.7%)
    etalra ~ 0.441 # Table 1 'omega2 beta' = 0.441 (RSE 47.4%)
    etalvc ~ 1.29 # Table 1 'omega2 Vp parent,u/F' = 1.29; on V2 (central) in NONMEM Table S6 (i) and Table S3
    etalkns ~ 0.693 # Table 1 'omega2 Kb parent' = 0.693 (RSE 46.6%)
    etalcl_form_m1 ~ 0.726 # Table 1 'omega2 CLF M1,u/F' = 0.726 (RSE 73.8%)
    etalkns_m1 ~ 0.43 # Table 1 'omega2 Kb2' = 0.43 (RSE 86%)
    etalcl_bile_m1 ~ 0.487 # Table 1 'omega2 CLbile u/F' = 0.487 (RSE 83%)
    etalcl_m2 ~ 0.443 # Table 1 'omega2 CLM2,u/F' = 0.443 (RSE 60%)
    etalvmax_m2 ~ 0.415 # Table 1 'omega2 Vmax' = 0.415 (RSE 56.4%)
    etalkm_m2 ~ 1.99 # Table 1 'omega2 Km' = 1.99 (RSE 44.1%)
    etalogitimax ~ 8.97 # Table 1 'omega2 Imax' = 8.97 (RSE 63.3%), logit scale
    etalogite0 ~ 2.47 # Table 1 'omega2 PD0' = 2.47 (RSE 56.7%), logit scale

    # Residual error: proportional for concentrations (Table 1 sigma2 are
    # variances; SD = sqrt), additive on the logit scale for the pain score.
    propSd_Cunbound <- sqrt(0.044); label("Proportional residual SD, unbound (1R,2R)-tramadol (fraction)") # Table 1 'sigma2 prop,parent,u' = 0.044
    propSd <- sqrt(0.041); label("Proportional residual SD, total (1R,2R)-tramadol (fraction)") # Table 1 'sigma2 prop,parent,t' = 0.041
    propSd_Cunbound_m1 <- sqrt(0.181); label("Proportional residual SD, unbound (1R,2R)-M1 (fraction)") # Table 1 'sigma2 prop,M1,u' = 0.181
    propSd_m1 <- sqrt(0.217); label("Proportional residual SD, total (1R,2R)-M1 (fraction)") # Table 1 'sigma2 prop,M1,t' = 0.217 (row text says 'total tramadol')
    propSd_m2 <- sqrt(0.083); label("Proportional residual SD, unbound (1R,2R)-M2 (fraction)") # Table 1 'sigma2 prop,M2,u' = 0.083
    addSd_score <- sqrt(0.713); label("Additive residual SD of the pain score on the logit scale (unitless)") # Table 1 'sigma2 add' = 0.713 (RSE 8.0%)
  })

  model({
    # Diabetes groups from the paper's DMTYPE (0 none, 1 type 1, 2 type 2)
    t1dm <- DIS_DIAB * (1 - DIS_DIAB_TYPE2)
    t2dm <- DIS_DIAB * DIS_DIAB_TYPE2

    # Weibull absorption parameters
    gam1 <- exp(lgam1 + etalgam1)
    ra <- exp(lra + etalra)
    tlag <- exp(ltlag)

    # Unbound tramadol
    vc <- exp(lvc + etalvc) * (BMI / 26.045)^e_bmi_vc
    vp <- exp(lvp)
    q <- exp(lq)
    cl <- exp(lcl)
    kns <- exp(lkns + etalkns) * e_t1dm_kns^t1dm * e_t2dm_kns^t2dm
    fu <- 1 / (1 + kns)

    # M1
    cl_form_m1 <- exp(lcl_form_m1 + etalcl_form_m1) * (CYP2D6 / 0.610)^e_cyp2d6_cl_form_m1
    fprat_m1 <- exp(lfprat_m1)
    fp_m1 <- fprat_m1 / (1 + fprat_m1)
    ka_m1 <- exp(lka_m1)
    kns_m1 <- exp(lkns_m1 + etalkns_m1) * e_t1dm_kns_m1^t1dm * e_t2dm_kns_m1^t2dm
    fu_m1 <- 1 / (1 + kns_m1)
    vc_m1 <- exp(lvc_m1)
    vcu_m1 <- vc_m1 / fu_m1
    cl_m1 <- exp(lcl_m1)
    cl_bile_m1 <- exp(lcl_bile_m1 + etalcl_bile_m1)
    kehc_m1 <- exp(lkehc_m1)

    # M2
    vmax_m2 <- exp(lvmax_m2 + etalvmax_m2)
    km_m2 <- exp(lkm_m2 + etalkm_m2)
    vc_m2 <- exp(lvc_m2)
    cl_m2 <- exp(lcl_m2 + etalcl_m2) * (CYP2D6 / 0.610)^e_cyp2d6_cl_m2

    # Pain-score PD
    ke0 <- exp(lke0)
    logitimax_ind <- logitimax + etalogitimax
    imax <- expit(logitimax_ind)
    logite0_ind <- logite0 + etalogite0
    e0 <- expit(logite0_ind)
    ic50 <- exp(lic50) * e_t1dm_ic50^t1dm * e_t2dm_ic50^t2dm

    # Time-dependent absorption rate and gallbladder gate, both read against the
    # time after the most recent dose (the source fits single-dose data, where
    # this equals time since the dose).
    tdose <- tad()
    ka <- 0
    if (tdose > tlag) {
      ka <- 1 - exp(-(ra * (tdose - tlag))^gam1)
    }
    ehc <- 0
    if (tdose >= tmeal1 && tdose < tmeal1 + dge) {
      ehc <- 1
    }

    # Michaelis-Menten M2 formation on the unbound central amount (Table S6 (iii))
    form_m2 <- vmax_m2 * central / (km_m2 + central)

    d/dt(depot) <- -ka * depot
    d/dt(depot_m1) <- -ka_m1 * depot_m1
    d/dt(central) <- ka * depot - fm_m1 * cl_form_m1 / vc * central - q / vc * central + q / vp * peripheral1 -
      fm_other * cl / vc * central - form_m2
    d/dt(peripheral1) <- q / vc * central - q / vp * peripheral1
    d/dt(central_m1) <- ka_m1 * depot_m1 + fm_m1 * cl_form_m1 / vc * central - cl_m1 / vcu_m1 * central_m1 -
      cl_bile_m1 / vcu_m1 * central_m1 + ehc * kehc_m1 * gallbladder
    d/dt(gallbladder) <- cl_bile_m1 / vcu_m1 * central_m1 - ehc * kehc_m1 * gallbladder
    d/dt(central_m2) <- form_m2 - cl_m2 / vc_m2 * central_m2

    # First-pass split of the oral dose (Figure 1; Figure S19 mass balance)
    f(depot) <- 1 - fp_m1
    f(depot_m1) <- fp_m1

    # Concentrations (nmol/L). Total = unbound / fu, which equals the NONMEM
    # total-concentration states DADT(5) = DADT(3)/FUN and DADT(8) = DADT(6)/FUN2.
    Cunbound <- central / vc
    Cc <- Cunbound / fu
    Cunbound_m1 <- central_m1 / vcu_m1
    Cc_m1 <- Cunbound_m1 / fu_m1
    Cc_m2 <- central_m2 / vc_m2 # M2 was measured and modelled unbound only

    d/dt(effect) <- ke0 * (Cunbound_m1 - effect)

    # Pain score on the 0-100 mm VAS (the paper models the 0-1 fraction)
    score <- 100 * e0 * (1 - imax * effect / (ic50 + effect))

    Cunbound ~ prop(propSd_Cunbound)
    Cc ~ prop(propSd)
    Cunbound_m1 ~ prop(propSd_Cunbound_m1)
    Cc_m1 ~ prop(propSd_m1)
    Cc_m2 ~ prop(propSd_m2)
    score ~ logitNorm(addSd_score, 0, 100)
  })
}
