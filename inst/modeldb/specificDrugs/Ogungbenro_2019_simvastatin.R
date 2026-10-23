Ogungbenro_2019_simvastatin <- function() {
  description <- "Joint one-compartment population PK model for oral simvastatin (lactone parent, SV) and four plasma metabolites -- simvastatin acid (SVA), 6-hydroxymethyl simvastatin (HMSV), 6-hydroxymethyl simvastatin acid (HMSVA) and 3,5-dihydrodiol simvastatin (DHSV) -- in hyperlipidemic children and adolescents (Ogungbenro 2019). SV absorption is two parallel sequential zero- then first-order processes, the second lagged; each metabolite is a one-compartment species formed from its precursor's concentration and described by apparent (CL/V) rate constants because metabolite volumes were not identifiable. HMSV elimination is a saturable sigmoid (Hill) process and HMSVA is formed from both SVA and HMSV. SLCO1B1 c.521T>C (rs4149056, additive allele count) raises SVA formation and lowers its apparent volume; HMSVA and DHSV elimination decline with chronological age through a fractional Imax function. Doses in nmol and concentrations in nmol/L, matching the source. The Table 1 estimates are encoded as printed; they reproduce the paper's covariate relationships but give exposures 2- to 10-fold below its own VPC (Figure 3) -- see the vignette."
  reference <- paste(
    "Ogungbenro K, Wagner JB, Abdel-Rahman S, Leeder JS, Galetin A.",
    "A population pharmacokinetic model for simvastatin and its metabolites",
    "in children and adolescents.",
    "Eur J Clin Pharmacol. 2019;75(9):1227-1235.",
    "doi:10.1007/s00228-019-02697-y.",
    sep = " "
  )
  vignette <- "Ogungbenro_2019_simvastatin"
  units <- list(time = "h", dosing = "nmol", concentration = "nmol/L")

  # Log-scale residual-error parameters of the double-exponential error model
  # (Table 1 footnote c). They are thetas feeding a computed lnorm() SD, so the
  # canonical propSd / addSd matcher does not recognise them.
  paper_specific_residual_sds <- c(
    "expSd1Early",
    "expSd1Late",
    "expSd2",
    "mRes",
    "expSd1_acid",
    "expSd2_acid",
    "mRes_acid",
    "expSd1Early_hmsv",
    "expSd1Late_hmsv",
    "expSd2_hmsv",
    "mRes_hmsv",
    "expSd1_hmsva",
    "expSd2_hmsva",
    "mRes_hmsva",
    "expSd1Early_dhsv",
    "expSd1Late_dhsv",
    "expSd2_dhsv",
    "mRes_dhsv"
  )

  # Metabolite states hold CONCENTRATION (nmol/L): the paper estimated only
  # CL/V ratios because the metabolite volumes are structurally unidentifiable
  # (Methods, 'The structural and statistical model').
  compartmentData <- list(
    # Each depot receives the whole dose; the F1 / F2 split acts on its outflow.
    depot = list(analyte = "simvastatin", units = "nmol", specimen = "administration site", verified = TRUE),
    depot2 = list(analyte = "simvastatin", units = "nmol", specimen = "administration site", verified = TRUE),
    central = list(analyte = "simvastatin", units = "nmol", specimen = "plasma", verified = TRUE),
    central_acid = list(analyte = "simvastatin acid", units = "nmol/L", specimen = "plasma", verified = TRUE),
    central_hmsv = list(
      analyte = "6-hydroxymethyl simvastatin",
      units = "nmol/L",
      specimen = "plasma",
      verified = TRUE
    ),
    central_hmsva = list(
      analyte = "6-hydroxymethyl simvastatin acid",
      units = "nmol/L",
      specimen = "plasma",
      verified = TRUE
    ),
    central_dhsv = list(analyte = "3,5-dihydrodiol simvastatin", units = "nmol/L", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    AGE = list(
      description = "Chronological age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Drives the fractional Imax decline of HMSVA and DHSV apparent",
        "elimination (Eqs 3-4). The tabulated kel_hmsva and kel_dhsv are the",
        "AGE = 0 intercepts of those equations, not the value at a 14-year-old",
        "typical subject (see the vignette Errata on Table 1 footnote a).",
        "Cohort range 8.8-20.2 years (Supplemental Table S1)."
      ),
      source_name = "Age"
    ),
    SNP_SLCO1B1_RS4149056_HET = list(
      description = "SLCO1B1 c.521T>C heterozygous (TC) indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (wild-type TT, when SNP_SLCO1B1_RS4149056_HOM = 0)",
      notes = paste(
        "The paper uses an additive allele-count dummy rs56 = 0/1/2 for",
        "TT/TC/CC (Eqs 1-2); model() derives it as HET + 2 * HOM."
      ),
      source_name = "rs56 (rs4149056 c.521T>C)"
    ),
    SNP_SLCO1B1_RS4149056_HOM = list(
      description = "SLCO1B1 c.521T>C homozygous-variant (CC) indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (wild-type TT, when SNP_SLCO1B1_RS4149056_HET = 0)",
      notes = "Counts twice in the additive allele-count dummy rs56 (Eqs 1-2).",
      source_name = "rs56 (rs4149056 c.521T>C)"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Tested allometrically (fixed and estimated exponents) on formation/elimination and volume parameters; made the model worse and was not retained (Discussion)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 32L,
    n_studies = 1L,
    n_observations = "238 SV, 160 SVA, 230 HMSV, 207 HMSVA and 221 DHSV plasma samples; 7%, 38%, 10%, 19% and 14% below the LLOQ (Supplemental Table S1).",
    age_range = "8.8-20.2 years (mean 14.1)",
    weight_range = "30.2-149.2 kg (mean 79.4); most subjects overweight for age",
    height_range = "133.2-181.9 cm (mean 160.2)",
    bmi_range = "15.1-53.0 kg/m^2 (mean 30.2)",
    sex_female_pct = 53.1,
    race_ethnicity = "Majority Caucasian (30 of 32; Discussion)",
    disease_state = "Hyperlipidemic children and adolescents",
    dose_range = "Single oral dose of simvastatin: 10 mg for ages 8-17 years (n = 28), 20 mg for ages >= 18 years (n = 4). Sampling at 0, 0.5, 1, 1.5, 2, 3, 4, 6 and 8 h.",
    regions = "United States (Children's Mercy Hospital, Kansas City)",
    genotype = "SLCO1B1 c.521T>C (rs4149056) TT/TC/CC = 15/15/2; CYP3A5 rs776746 AA/AG/GG = 0/4/28 (Supplemental Table S1).",
    bloq_handling = "LLOQ 0.5 nmol/L for every analyte; BLQ samples retained as censored data (Beal M3).",
    notes = "Data from the genotype-stratified single-dose study of Wagner et al.; demographics from Supplemental Table S1. Fitted in NONMEM 7.4 (FOCE-I) on log-transformed concentrations."
  )

  ini({
    # Simvastatin (SV) -- Table 1. CL and V are apparent (oral) values.
    lcl <- log(1300); label("Apparent SV total elimination clearance CLSLe/F (L/h)") # Table 1 SV CLSLe/F = 1300 L/h (bootstrap 1210-1450)
    lvc <- log(110); label("Apparent SV volume of distribution VSL/F (L)") # Table 1 SV VSL/F = 110 L (bootstrap 98.7-114)

    # Two parallel absorption processes, each zero-order into its depot then
    # first-order out. Process 1 (depot) is the slow, unlagged one; process 2
    # (depot2) is the fast, lagged one.
    ld1 <- log(0.069); label("Zero-order duration of absorption process 1, D1 (h)") # Table 1 SV D1 = 0.069 (unit printed as h-1; Methods define D1 as a duration)
    ld2 <- log(0.39); label("Zero-order duration of absorption process 2, D2 (h)") # Table 1 SV D2 = 0.39 (unit printed as h-1; Methods define D2 as a duration)
    lka_slow <- log(0.030); label("First-order absorption rate constant of process 1, ka1 (1/h)") # Table 1 SV ka1 = 0.030 1/h
    lka_fast <- log(0.41); label("First-order absorption rate constant of process 2, ka2 (1/h)") # Table 1 SV ka2 = 0.41 1/h
    logitfrel <- -log(0.78); label("Logit of the dose fraction entering absorption process 1, F1 (logit units)") # Table 1 SV BA = 0.78 with footnote F1 = 1/(1 + BA), so logit(F1) = -log(BA) exactly
    ltlag <- log(0.18); label("Lag time of absorption process 2, ALAG (h)") # Table 1 SV ALAG = 0.18 h

    # Simvastatin acid (SVA) -- Table 1 and Eqs 1-2
    lkmet_acid <- log(0.043); label("SVA formation rate constant CLLA/VSVA for SLCO1B1 c.521TT (1/h)") # Table 1 SVA CLLA/VSVA = 0.043 1/h; Eq 1 intercept
    lkel_acid <- log(0.13); label("SVA elimination rate constant CLSVAe/VSVA for SLCO1B1 c.521TT (1/h)") # Table 1 SVA CLSVAe/VSVA = 0.13 1/h; Eq 2 intercept
    e_slco1b1_521_vc_acid <- -0.37; label("Fractional change in VSVA per c.521C allele (unitless)") # Table 1 theta561 = -0.37; Eqs 1-2
    e_slco1b1_521_kmet_acid <- 0.93; label("Fractional change in CLLA per c.521C allele (unitless)") # Table 1 theta562 = 0.93; Eq 1

    # 6-hydroxymethyl simvastatin (HMSV) -- Table 1
    lkmet_hmsv <- log(16.1); label("HMSV formation rate constant CLLH/VHSV (1/h)") # Table 1 HMSV CLLH/VHSV = 16.1 1/h
    lvmax_hmsv <- log(660); label("Maximum HMSV elimination rate CLEmax/VHSV (nmol/L/h)") # Table 1 HMSV CLEmax/VHSV = 660 nM/h
    lkm_hmsv <- log(20); label("HMSV concentration at half-maximal elimination rate CLEC50 (nmol/L)") # Table 1 HMSV CLEC50 = 20 nM
    hill_hmsv <- 0.86; label("Sigmoidicity of HMSV saturable elimination, gamma (unitless)") # Table 1 HMSV gamma = 0.86

    # 6-hydroxymethyl simvastatin acid (HMSVA) -- Table 1 and Eq 3
    lkmet_hmsva <- log(0.48); label("HMSVA formation rate constant from SVA, CLAH/VHSVA (1/h)") # Table 1 HMSVA CLAH/VHSVA = 0.48 1/h
    lkmet_hmsva_hmsv <- log(0.55); label("HMSVA formation rate constant from HMSV, CLHSVHSVA/VHSVA (1/h)") # Table 1 HMSVA CLHSVHSVA/VHSVA = 0.55 1/h
    lkel_hmsva <- log(4.91); label("HMSVA elimination rate constant CLHSVAe/VHSVA at AGE = 0 (1/h)") # Table 1 HMSVA CLHSVAe/VHSVA = 4.91 1/h; Eq 3 intercept
    e_age_kel_hmsva <- fixed(1); label("Maximal fractional decline of kel_hmsva with age, Fra (unitless)") # Table 1 thetaAGE1 'age (Fra) on CLHSVAe' = 1 (fixed); Eq 3
    lage50_kel_hmsva <- log(5.34); label("Age at half of the maximal decline of kel_hmsva (years)") # Table 1 thetaAGE2 'age (50%) on CLHSVAe' = 5.34 year; Eq 3

    # 3,5-dihydrodiol simvastatin (DHSV) -- Table 1 and Eq 4
    lkmet_dhsv <- log(4.11); label("DHSV formation rate constant CLLD/VDHSV (1/h)") # Table 1 DHSV CLLD/VDHSV = 4.11 1/h
    lkel_dhsv <- log(12.9); label("DHSV elimination rate constant CLDHSe/VDHSV at AGE = 0 (1/h)") # Table 1 DHSV CLDHSe/VDHSV = 12.9 1/h; Eq 4 intercept
    e_age_kel_dhsv <- 0.90; label("Maximal fractional decline of kel_dhsv with age, Fra (unitless)") # Table 1 thetaAGE3 'age (Fra) on CLDHSe' = 0.90; Eq 4
    lage50_kel_dhsv <- log(6.3); label("Age at half of the maximal decline of kel_dhsv (years)") # Table 1 thetaAGE4 'age (50%) on CLDHSe' = 6.3 year; Eq 4

    # Between-subject variability -- Table 1, reported as CV = sqrt(exp(omega^2) - 1)
    # (footnote b), so omega^2 = log(1 + CV^2). Diagonal: the paper estimated
    # variances only.
    etald2 ~ 1.664854 # Table 1 BSV D2 CV = 2.07 -> log(1 + 2.07^2)
    etalka_fast ~ 0.155378 # Table 1 BSV ka2 CV = 0.41 -> log(1 + 0.41^2)
    etaltlag ~ 1.358409 # Table 1 BSV ALAG CV = 1.70 -> log(1 + 1.70^2)
    etalcl ~ 0.334255 # Table 1 BSV CLSLe/F CV = 0.63 -> log(1 + 0.63^2)
    etalvc ~ 1.742026 # Table 1 BSV VSL/F CV = 2.17 -> log(1 + 2.17^2)
    etalkmet_acid ~ 0.653158 # Table 1 BSV CLLA/VSVA CV = 0.96 -> log(1 + 0.96^2)
    etalkel_acid ~ 0.484954 # Table 1 BSV CLSVAe/VSVA CV = 0.79 -> log(1 + 0.79^2)
    etalkmet_hmsv ~ 0.065413 # Table 1 BSV CLLH/VHSV CV = 0.26 -> log(1 + 0.26^2)
    etalkm_hmsv ~ 0.370805 # Table 1 BSV CLEC50 CV = 0.67 -> log(1 + 0.67^2)
    etalkmet_hmsva ~ 0.103369 # Table 1 BSV CLAH/VHSVA CV = 0.33 -> log(1 + 0.33^2)
    etalkmet_hmsva_hmsv ~ 0.231191 # Table 1 BSV CLHSVHSVA/VHSVA CV = 0.51 -> log(1 + 0.51^2)
    etalkel_hmsva ~ 0.031886 # Table 1 BSV CLHSVAe/VHSVA CV = 0.18 -> log(1 + 0.18^2)
    etalkmet_dhsv ~ 0.080750 # Table 1 BSV CLLD/VDHSV CV = 0.29 -> log(1 + 0.29^2)
    etalkel_dhsv ~ 0.075478 # Table 1 BSV CLDHSe/VDHSV CV = 0.28 -> log(1 + 0.28^2)

    # Residual error -- double-exponential model on log-transformed data
    # (Table 1 footnote c): ln(y) = ln(f + m) + f/(f + m) * eps1 + m/(f + m) * eps2.
    # Table 1 reports the VARIANCES of eps1 and eps2, so each SD below is
    # sqrt(variance). SV, HMSV and DHSV carry separate eps1 variances for
    # 0-2 h and 2-8 h after the dose.
    expSd1Early <- sqrt(0.38); label("SV log-scale SD of eps1, 0-2 h (unitless)") # Table 1 eps1SV = 0.38 (variance; 0-2 h)
    expSd1Late <- sqrt(0.32); label("SV log-scale SD of eps1, 2-8 h (unitless)") # Table 1 eps2SV = 0.32 (variance; 2-8 h)
    expSd2 <- sqrt(0.30); label("SV log-scale SD of eps2 (unitless)") # Table 1 eps3SV = 0.30 (variance)
    mRes <- 0.012; label("SV double-exponential error offset m (nmol/L)") # Table 1 mSV = 0.012
    expSd1_acid <- sqrt(0.26); label("SVA log-scale SD of eps1 (unitless)") # Table 1 eps1SVA = 0.26 (variance)
    expSd2_acid <- sqrt(0.16); label("SVA log-scale SD of eps2 (unitless)") # Table 1 eps2SVA = 0.16 (variance)
    mRes_acid <- 0.22; label("SVA double-exponential error offset m (nmol/L)") # Table 1 mSVA = 0.22
    expSd1Early_hmsv <- sqrt(0.29); label("HMSV log-scale SD of eps1, 0-2 h (unitless)") # Table 1 eps1HMSV = 0.29 (variance; 0-2 h)
    expSd1Late_hmsv <- sqrt(0.24); label("HMSV log-scale SD of eps1, 2-8 h (unitless)") # Table 1 eps2HMSV = 0.24 (variance; 2-8 h)
    expSd2_hmsv <- sqrt(0.12); label("HMSV log-scale SD of eps2 (unitless)") # Table 1 eps3HMSV = 0.12 (variance)
    mRes_hmsv <- fixed(0.0001); label("HMSV double-exponential error offset m (nmol/L)") # Table 1 mHMSV = 0.0001 (fixed)
    expSd1_hmsva <- sqrt(0.18); label("HMSVA log-scale SD of eps1 (unitless)") # Table 1 eps1HMSVA = 0.18 (variance)
    expSd2_hmsva <- sqrt(0.1); label("HMSVA log-scale SD of eps2 (unitless)") # Table 1 eps2HMSVA = 0.1 (variance)
    mRes_hmsva <- fixed(0.0001); label("HMSVA double-exponential error offset m (nmol/L)") # Table 1 mHMSVA = 0.0001 (fixed)
    expSd1Early_dhsv <- sqrt(0.24); label("DHSV log-scale SD of eps1, 0-2 h (unitless)") # Table 1 eps1DHSV = 0.24 (variance; 0-2 h)
    expSd1Late_dhsv <- sqrt(0.21); label("DHSV log-scale SD of eps1, 2-8 h (unitless)") # Table 1 eps2DHSV = 0.21 (variance; 2-8 h)
    expSd2_dhsv <- sqrt(0.086); label("DHSV log-scale SD of eps2 (unitless)") # Table 1 eps3DHSV = 0.086 (variance)
    mRes_dhsv <- fixed(0.0001); label("DHSV double-exponential error offset m (nmol/L)") # Table 1 mDHSV = 0.0001 (fixed)
  })

  model({
    # Additive allele-count dummy for SLCO1B1 c.521T>C: 0 = TT, 1 = TC, 2 = CC.
    rs56 <- SNP_SLCO1B1_RS4149056_HET + 2 * SNP_SLCO1B1_RS4149056_HOM

    # SV
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    kel <- cl / vc
    d1 <- exp(ld1)
    d2 <- exp(ld2 + etald2)
    ka_slow <- exp(lka_slow)
    ka_fast <- exp(lka_fast + etalka_fast)
    frel <- expit(logitfrel)
    tlag <- exp(ltlag + etaltlag)

    # SVA, Eqs 1-2: the c.521T>C effect on VSVA divides both apparent rates.
    f_vc_acid <- 1 + rs56 * e_slco1b1_521_vc_acid
    kmet_acid <- exp(lkmet_acid + etalkmet_acid) * (1 + rs56 * e_slco1b1_521_kmet_acid) / f_vc_acid
    kel_acid <- exp(lkel_acid + etalkel_acid) / f_vc_acid

    # HMSV
    kmet_hmsv <- exp(lkmet_hmsv + etalkmet_hmsv)
    vmax_hmsv <- exp(lvmax_hmsv)
    km_hmsv <- exp(lkm_hmsv + etalkm_hmsv)

    # HMSVA, Eq 3
    kmet_hmsva <- exp(lkmet_hmsva + etalkmet_hmsva)
    kmet_hmsva_hmsv <- exp(lkmet_hmsva_hmsv + etalkmet_hmsva_hmsv)
    age50_kel_hmsva <- exp(lage50_kel_hmsva)
    kel_hmsva <- exp(lkel_hmsva + etalkel_hmsva) * (1 - e_age_kel_hmsva * AGE / (age50_kel_hmsva + AGE))

    # DHSV, Eq 4
    kmet_dhsv <- exp(lkmet_dhsv + etalkmet_dhsv)
    age50_kel_dhsv <- exp(lage50_kel_dhsv)
    kel_dhsv <- exp(lkel_dhsv + etalkel_dhsv) * (1 - e_age_kel_dhsv * AGE / (age50_kel_dhsv + AGE))

    # SV amounts in nmol; metabolite states are concentrations in nmol/L,
    # formed at (CL_form / V_metab) * C_precursor. CLSLe and the metabolite
    # elimination terms are TOTAL elimination, so formation does not appear as
    # a separate loss from the precursor.
    #
    # Every dose is given twice, once to each depot, as a zero-order input
    # (rate = -2). The F1 / F2 split is applied to the first-order outflow
    # rather than through f(depot) / f(depot2): in rxode2 5.1.8, f() on two
    # depots combined with a lagged modelled duration is mis-solved (the lagged
    # depot goes negative) whenever the model branches on CMT, which every
    # multi-endpoint model does. The two forms give identical central and
    # metabolite profiles; the depots here hold the full, unsplit dose.
    d/dt(depot) <- -ka_slow * depot
    d/dt(depot2) <- -ka_fast * depot2
    d/dt(central) <- frel * ka_slow * depot + (1 - frel) * ka_fast * depot2 - kel * central
    Cc <- central / vc
    d/dt(central_acid) <- kmet_acid * Cc - kel_acid * central_acid
    d/dt(central_hmsv) <- kmet_hmsv * Cc -
      vmax_hmsv * central_hmsv^hill_hmsv / (km_hmsv^hill_hmsv + central_hmsv^hill_hmsv)
    d/dt(central_hmsva) <- kmet_hmsva * central_acid + kmet_hmsva_hmsv * central_hmsv -
      kel_hmsva * central_hmsva
    d/dt(central_dhsv) <- kmet_dhsv * Cc - kel_dhsv * central_dhsv

    dur(depot) <- d1
    dur(depot2) <- d2
    alag(depot2) <- tlag

    Cc_acid <- central_acid
    Cc_hmsv <- central_hmsv
    Cc_hmsva <- central_hmsva
    Cc_dhsv <- central_dhsv

    # Double-exponential residual error: log-scale SD
    # sqrt((f/(f+m))^2 * var(eps1) + (m/(f+m))^2 * var(eps2)).
    tpost <- tad(depot)
    early <- tpost <= 2
    w_sv <- Cc / (Cc + mRes)
    sd1_sv <- expSd1Early * early + expSd1Late * (1 - early)
    lnSd_sv <- sqrt((w_sv * sd1_sv)^2 + ((1 - w_sv) * expSd2)^2)
    w_acid <- Cc_acid / (Cc_acid + mRes_acid)
    lnSd_acid <- sqrt((w_acid * expSd1_acid)^2 + ((1 - w_acid) * expSd2_acid)^2)
    w_hmsv <- Cc_hmsv / (Cc_hmsv + mRes_hmsv)
    sd1_hmsv <- expSd1Early_hmsv * early + expSd1Late_hmsv * (1 - early)
    lnSd_hmsv <- sqrt((w_hmsv * sd1_hmsv)^2 + ((1 - w_hmsv) * expSd2_hmsv)^2)
    w_hmsva <- Cc_hmsva / (Cc_hmsva + mRes_hmsva)
    lnSd_hmsva <- sqrt((w_hmsva * expSd1_hmsva)^2 + ((1 - w_hmsva) * expSd2_hmsva)^2)
    w_dhsv <- Cc_dhsv / (Cc_dhsv + mRes_dhsv)
    sd1_dhsv <- expSd1Early_dhsv * early + expSd1Late_dhsv * (1 - early)
    lnSd_dhsv <- sqrt((w_dhsv * sd1_dhsv)^2 + ((1 - w_dhsv) * expSd2_dhsv)^2)

    Cc ~ lnorm(lnSd_sv)
    Cc_acid ~ lnorm(lnSd_acid)
    Cc_hmsv ~ lnorm(lnSd_hmsv)
    Cc_hmsva ~ lnorm(lnSd_hmsva)
    Cc_dhsv ~ lnorm(lnSd_dhsv)
  })
}
