Favie_2020_neonatalHypothermia <- function() {
  description <- "Integrated population PK model of seven drugs (morphine, midazolam, lidocaine, phenobarbital, amoxicillin, benzylpenicillin, gentamicin) and five metabolites (M3G, M6G, 1-hydroxymidazolam, hydroxymidazolam glucuronide, MEGX) in term encephalopathic neonates treated with therapeutic hypothermia (Favie 2020, PharmaCool). One- or two-compartment disposition per compound with parent-to-metabolite formation, fixed allometric birth-weight scaling, a fixed sigmoidal gestational-age maturation function, a linear postnatal-age (organ recovery) effect on clearance for high- and intermediate-clearance compounds, a linear body-temperature effect on clearance for intermediate-clearance compounds, and a common between-compound clearance random effect."
  reference <- paste(
    "Favie LMA, de Haan TR, Bijleveld YA, Rademaker CMA, Egberts TCG,",
    "Nuytemans DHGM, Mathot RAA, Groenendaal F, Huitema ADR.",
    "Prediction of Drug Exposure in Critically Ill Encephalopathic Neonates",
    "Treated With Therapeutic Hypothermia Based on a Pooled Population",
    "Pharmacokinetic Analysis of Seven Drugs and Five Metabolites.",
    "Clin Pharmacol Ther. 2020;108(5):1098-1106.",
    "doi:10.1002/cpt.1917.",
    sep = " "
  )
  vignette <- "Favie_2020_neonatalHypothermia"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT_BIRTH = list(
      description = "Birth weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Cohort mean 3.38 +/- 0.617 kg (Table 1). Allometric scaling",
        "relative to 3.5 kg with fixed exponents 0.75 on every clearance",
        "and intercompartmental clearance and 1 on every volume (Methods",
        "'Body size'; Table S1 footnote). The source abbreviates it 'BW'",
        "and describes it as body weight in the Methods, but Table 1,",
        "Figure 1 and Figure S2 all identify the covariate as birth weight."
      ),
      source_name = "BW"
    ),
    GA = list(
      description = "Gestational age at birth",
      units = "weeks",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Cohort mean 39.7 +/- 1.66 weeks, range 36-42 (Table 1, Methods).",
        "Enters every clearance through a fixed sigmoidal Hill maturation",
        "function (TM50 54.2 weeks, Hill 3.92) normalised to its value at",
        "40 weeks. The paper uses GA in place of postmenstrual age because",
        "data were collected only in the first 5 days of life."
      ),
      source_name = "GA"
    ),
    PNA = list(
      description = "Postnatal age (chronological since birth)",
      units = "months",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-varying. The source models PNA in HOURS of life (0-120 h in",
        "the data). Canonical PNA carries months, so model() recovers",
        "hours as PNA * 24 * 30.4375. Enters clearance linearly,",
        "(1 + slope * PNA_hours); the slope is 1.23 %/h for",
        "high-clearance compounds and 0.54 %/h for intermediate-clearance",
        "compounds, not applied to phenobarbital. Supply PNA on every",
        "record (or solve with linear covariate interpolation), because",
        "the slope is steep enough that a coarse step function biases",
        "clearance."
      ),
      source_name = "PNA"
    ),
    BODYTEMP = list(
      description = "Body temperature",
      units = "degC",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-varying; reference 36.5 degC. The source did not use",
        "measured temperature: it reconstructed each neonate's profile",
        "from the recorded start and end of therapeutic hypothermia as",
        "33.5 degC during hypothermia, rewarming at 0.4 degC/h (7.5 h)",
        "to 36.5 degC, then 36.5 degC (Methods 'Body temperature').",
        "Only the intermediate-clearance compounds (morphine, midazolam,",
        "1-hydroxymidazolam) carry the effect."
      ),
      source_name = "TEMP"
    )
  )

  compartmentData <- list(
    central = list(analyte = "morphine", units = "mg", specimen = "plasma", verified = TRUE),
    central_m3g = list(analyte = "morphine-3-glucuronide", units = "mg", specimen = "plasma", verified = TRUE),
    central_m6g = list(analyte = "morphine-6-glucuronide", units = "mg", specimen = "plasma", verified = TRUE),
    central_midazolam = list(analyte = "midazolam", units = "mg", specimen = "plasma", verified = TRUE),
    central_1ohm = list(analyte = "1-hydroxymidazolam", units = "mg", specimen = "plasma", verified = TRUE),
    central_hmg = list(analyte = "hydroxymidazolam glucuronide", units = "mg", specimen = "plasma", verified = TRUE),
    central_lidocaine = list(analyte = "lidocaine", units = "mg", specimen = "plasma", verified = TRUE),
    central_megx = list(analyte = "monoethylglycinexylidide", units = "mg", specimen = "plasma", verified = TRUE),
    central_phenobarbital = list(analyte = "phenobarbital", units = "mg", specimen = "plasma", verified = TRUE),
    central_amoxicillin = list(analyte = "amoxicillin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_amoxicillin = list(analyte = "amoxicillin", units = "mg", specimen = "plasma", verified = TRUE),
    central_benzylpenicillin = list(analyte = "benzylpenicillin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_benzylpenicillin = list(
      analyte = "benzylpenicillin",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    central_gentamicin = list(analyte = "gentamicin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_gentamicin = list(analyte = "gentamicin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 192L,
    n_studies = 1L,
    ga_range = "36-42 weeks (mean 39.7, SD 1.66)",
    weight_range = "birth weight mean 3.38 kg (SD 0.617)",
    pna_range = "0-120 h (samples on days 2-5 of life)",
    sex_female_pct = 38.5,
    disease_state = paste(
      "Term and near-term neonates with moderate or severe neonatal",
      "encephalopathy after perinatal asphyxia, treated with therapeutic",
      "hypothermia (33.5 degC for 72 h started within 6 h of birth,",
      "then rewarming at 0.4 degC/h)."
    ),
    dose_range = paste(
      "Clinical-care dosing of every drug; no study-protocol dosing",
      "(doses and regimens are described in the individual PharmaCool",
      "publications)."
    ),
    regions = "Netherlands and Belgium (12 level III NICUs)",
    per_drug_counts = paste(
      "Patients / samples (Table 1): morphine 180 / 534, amoxicillin",
      "125 / 1,280, midazolam 118 / 376, phenobarbital 113 / 378,",
      "gentamicin 47 / 471, benzylpenicillin 43 / 416, lidocaine 28 / 77.",
      "Metabolite samples are counted with their parent drug."
    ),
    notes = paste(
      "Prospective multicentre observational PharmaCool cohort",
      "(NTR2529). Every neonate appears once in the pooled dataset but",
      "may contribute data for several drugs."
    )
  )

  ini({
    # ---- Typical clearances (Table 3 and Table S1) ------------------------
    # All values are for a neonate of birth weight 3.5 kg, GA 40 weeks
    # (280 days), PNA 0 h and body temperature 36.5 degC (Table S1
    # footnote #). Metabolite clearances and volumes are apparent values
    # relative to the unestimated formation fraction (Table S1 footnote).
    lcl <- log(0.811); label("Morphine clearance CL (L/h)") # Table 3 / Table S1 Morphine Cl = 0.811 L/h
    lcl_m3g <- log(0.241); label("M3G apparent clearance CL/F (L/h)") # Table 3 / Table S1 M3G Cl = 0.241 L/h
    lcl_m6g <- log(0.765); label("M6G apparent clearance CL/F (L/h)") # Table 3 / Table S1 M6G Cl = 0.765 L/h
    lcl_midazolam <- log(0.511); label("Midazolam clearance CL (L/h)") # Table 3 / Table S1 Midazolam Cl = 0.511 L/h
    lcl_1ohm <- log(1.72); label("1-hydroxymidazolam apparent clearance CL/F (L/h)") # Table 3 / Table S1 OHM Cl = 1.72 L/h
    lcl_hmg <- log(0.111); label("Hydroxymidazolam glucuronide apparent clearance CL/F (L/h)") # Table 3 / Table S1 HMG Cl = 0.111 L/h
    lcl_lidocaine <- log(0.937); label("Lidocaine clearance CL (L/h)") # Table 3 / Table S1 Lidocaine Cl = 0.937 L/h
    lcl_megx <- log(1.51); label("MEGX apparent clearance CL/F (L/h)") # Table 3 / Table S1 MEGX Cl = 1.51 L/h
    lcl_phenobarbital <- log(0.00930); label("Phenobarbital clearance CL (L/h)") # Table 3 / Table S1 Phenobarbital Cl = 0.00930 L/h
    lcl_amoxicillin <- log(0.178); label("Amoxicillin clearance CL (L/h)") # Table 3 / Table S1 Amoxicillin Cl = 0.178 L/h
    lcl_benzylpenicillin <- log(0.359); label("Benzylpenicillin clearance CL (L/h)") # Table 3 / Table S1 Benzylpenicillin Cl = 0.359 L/h
    lcl_gentamicin <- log(0.108); label("Gentamicin clearance CL (L/h)") # Table 3 / Table S1 Gentamicin Cl = 0.108 L/h

    # ---- Intercompartmental clearances (Table S1) --------------------------
    lq_amoxicillin <- log(0.686); label("Amoxicillin intercompartmental clearance Q (L/h)") # Table S1 Amoxicillin Q = 0.686 L/h
    lq_benzylpenicillin <- log(0.178); label("Benzylpenicillin intercompartmental clearance Q (L/h)") # Table S1 Benzylpenicillin Q = 0.178 L/h
    lq_gentamicin <- log(0.158); label("Gentamicin intercompartmental clearance Q (L/h)") # Table S1 Gentamicin Q = 0.158 L/h

    # ---- Volumes, all fixed to the individual-drug models (Table S1) -------
    lvc <- fixed(log(8.88)); label("Morphine volume of distribution V (L)") # Table S1 Morphine V = 8.88 FIX
    lvc_m3g <- fixed(log(0.264)); label("M3G apparent volume V/F (L)") # Table S1 M3G V = 0.264 FIX
    lvc_m6g <- fixed(log(4.53)); label("M6G apparent volume V/F (L)") # Table S1 M6G V = 4.53 FIX
    lvc_midazolam <- fixed(log(5.42)); label("Midazolam volume of distribution V (L)") # Table S1 Midazolam V = 5.42 FIX
    lvc_1ohm <- fixed(log(4.18)); label("1-hydroxymidazolam apparent volume V/F (L)") # Table S1 OHM V = 4.18 FIX
    lvc_hmg <- fixed(log(1.06)); label("Hydroxymidazolam glucuronide apparent volume V/F (L)") # Table S1 HMG V = 1.06 FIX
    lvc_lidocaine <- fixed(log(10.9)); label("Lidocaine volume of distribution V (L)") # Table S1 Lidocaine V = 10.9 FIX
    lvc_megx <- fixed(log(4.59)); label("MEGX apparent volume V/F (L)") # Table S1 MEGX V = 4.59 FIX
    lvc_phenobarbital <- fixed(log(3.60)); label("Phenobarbital volume of distribution V (L)") # Table S1 Phenobarbital V = 3.60 FIX
    lvc_amoxicillin <- fixed(log(1.21)); label("Amoxicillin central volume Vc (L)") # Table S1 Amoxicillin Vc = 1.21 FIX
    lvp_amoxicillin <- fixed(log(1.21)); label("Amoxicillin peripheral volume Vp (L)") # Table S1 Amoxicillin Vp = 1.21 FIX
    lvc_benzylpenicillin <- fixed(log(2.08)); label("Benzylpenicillin central volume Vc (L)") # Table S1 Benzylpenicillin Vc = 2.08 FIX
    lvp_benzylpenicillin <- fixed(log(3.55)); label("Benzylpenicillin peripheral volume Vp (L)") # Table S1 Benzylpenicillin Vp = 3.55 FIX
    lvc_gentamicin <- fixed(log(1.63)); label("Gentamicin central volume Vc (L)") # Table S1 Gentamicin Vc = 1.63 FIX
    lvp_gentamicin <- fixed(log(1.52)); label("Gentamicin peripheral volume Vp (L)") # Table S1 Gentamicin Vp = 1.52 FIX

    # ---- Body size and maturation (Methods 'Body size', 'Maturation') -----
    e_wt_cl_q <- fixed(0.75); label("Allometric exponent of birth weight on every CL and Q (unitless)") # Methods 'Body size': exponent 0.75 for clearance
    e_wt_vc_vp <- fixed(1); label("Allometric exponent of birth weight on every volume (unitless)") # Methods 'Body size': exponent 1 for volume of distribution
    ga_tm50 <- fixed(54.2); label("Gestational age at 50% of mature clearance, TM50 (weeks)") # Methods 'Maturation': TM50 fixed to 54.2 weeks (Knosgaard et al.)
    ga_hill <- fixed(3.92); label("Hill coefficient of the gestational-age maturation function (unitless)") # Methods 'Maturation': Hill coefficient fixed to 3.92

    # ---- Postnatal-age (organ recovery) and temperature effects ------------
    e_pna_cl_high <- 0.0123; label("Fractional increase in clearance per hour of life, high-clearance compounds (1/h)") # Table 3 / Table S1 'PNA high clearance' = 1.23 %/h
    e_pna_cl_int <- 0.0054; label("Fractional increase in clearance per hour of life, intermediate-clearance compounds (1/h)") # Table 3 / Table S1 'PNA intermediate clearance' = 0.540 %/h
    e_bodytemp_cl <- 0.0683; label("Fractional change in clearance per degC of body temperature, intermediate-clearance compounds (1/degC)") # Table 3 / Table S1 'TEMP intermediate clearance' = 6.83 %/degC

    # ---- Common clearance random effect scaling (Table S1 'OMEGA structure')
    sd_ratio_cl <- fixed(1); label("Scale of the common clearance eta on morphine CL (unitless)") # Table S1 Morphine 'Common THETA on IIV Cl' = 1 FIX
    sd_ratio_cl_midazolam <- 1.46; label("Scale of the common clearance eta on midazolam CL (unitless)") # Table S1 Midazolam 'Common THETA on IIV Cl' = 1.46
    sd_ratio_cl_1ohm <- 1.38; label("Scale of the common clearance eta on 1-hydroxymidazolam CL (unitless)") # Table S1 OHM 'Common THETA on IIV Cl' = 1.38
    sd_ratio_cl_m3g <- 0.532; label("Scale of the common clearance eta on M3G CL (unitless)") # Table S1 M3G 'Common THETA on IIV Cl' = 0.532
    sd_ratio_cl_m6g <- 0.504; label("Scale of the common clearance eta on M6G CL (unitless)") # Table S1 M6G 'Common THETA on IIV Cl' = 0.504
    sd_ratio_cl_hmg <- 0.870; label("Scale of the common clearance eta on hydroxymidazolam glucuronide CL (unitless)") # Table S1 HMG 'Common THETA on IIV Cl' = 0.870
    sd_ratio_cl_amoxicillin <- 0.541; label("Scale of the common clearance eta on amoxicillin CL (unitless)") # Table S1 Amoxicillin 'Common THETA on IIV Clc' = 0.541
    sd_ratio_cl_benzylpenicillin <- 0.847; label("Scale of the common clearance eta on benzylpenicillin CL (unitless)") # Table S1 Benzylpenicillin 'Common THETA on IIV Clc' = 0.847
    sd_ratio_cl_gentamicin <- 0.327; label("Scale of the common clearance eta on gentamicin CL (unitless)") # Table S1 Gentamicin 'Common THETA on IIV Clc' = 0.327
    sd_ratio_e_pna_cl_int <- 0.765; label("Scale of the postnatal-age-effect eta for intermediate-clearance compounds (unitless)") # Table S1 'Scaling factor for intermediate clearance' = 0.765; correlation of PNA effects fixed to 100%

    # ---- Interindividual variability ---------------------------------------
    # Table S1 reports every IIV as a relative standard deviation (rsd),
    # which is sqrt(OMEGA): the predecessor morphine paper prints
    # 'variance 0.224 (rsd 47.3%)'. The variances below are rsd^2.
    etalcl_common ~ 0.190969 # Table S1 'Common OMEGA on Cl, rsd' 43.7% -> 0.437^2
    etalcl ~ 0.0841 # Table S1 Morphine 'IIV Cl, rsd' 29.0% -> 0.290^2
    # M3G and M6G: variances from Table S1 (49.6% and 52.1% rsd). The
    # covariance is not printed in Table S1, but Table 2 reports a 94.8%
    # M3G-M6G clearance correlation that the common eta alone (16.5%)
    # cannot produce; every other Table 2 cell is reproduced by the common
    # eta to within 0.4 points. The covariance is back-solved from Table 2
    # as 0.948 * sqrt(0.496^2 + (0.532 * 0.437)^2) *
    # sqrt(0.521^2 + (0.504 * 0.437)^2) - 0.532 * 0.504 * 0.437^2 = 0.242532
    # (compound-specific eta correlation 0.939). The Results text prints
    # 96.2% for the same correlation; the table value is used.
    etalcl_m3g + etalcl_m6g ~ c(0.246016, 0.242532, 0.271441) # Table S1 M3G/M6G 'IIV Cl, rsd' 49.6%/52.1%; covariance back-solved from Table 2 M3G-M6G 94.8%
    etalcl_midazolam ~ 0.3136 # Table S1 Midazolam 'IIV Cl, rsd' 56.0% -> 0.560^2
    etalcl_1ohm ~ 0.268324 # Table S1 OHM 'IIV Cl, rsd' 51.8% -> 0.518^2
    etalcl_hmg ~ 0.112225 # Table S1 HMG 'IIV Cl, rsd' 33.5% -> 0.335^2
    etalcl_lidocaine ~ 0.044944 # Table S1 Lidocaine 'IIV Cl, rsd' 21.2% -> 0.212^2
    etalcl_megx ~ 0.641601 # Table S1 MEGX 'IIV Cl, rsd' 80.1% -> 0.801^2
    etalcl_phenobarbital ~ 0.331776 # Table S1 Phenobarbital 'IIV Cl, rsd' 57.6% -> 0.576^2
    etalcl_amoxicillin ~ 0.163216 # Table S1 Amoxicillin 'IIV Cl, rsd' 40.4% -> 0.404^2
    etalcl_benzylpenicillin ~ 0.139876 # Table S1 Benzylpenicillin 'IIV Cl, rsd' 37.4% -> 0.374^2
    etalcl_gentamicin ~ 0.046656 # Table S1 Gentamicin 'IIV Cl, rsd' 21.6% -> 0.216^2
    etae_pna_cl_high ~ 0.512656 # Table S1 'IIV -- high clearance, rsd' 71.6% -> 0.716^2
    etalvc ~ fixed(0.463761) # Table S1 Morphine 'IIV V, rsd' 68.1 FIX -> 0.681^2
    etalvc_midazolam ~ fixed(0.933156) # Table S1 Midazolam 'IIV V, rsd' 96.6 FIX -> 0.966^2
    etalvc_hmg ~ fixed(0.837225) # Table S1 HMG 'IIV V, rsd' 91.5 FIX -> 0.915^2
    etalvc_lidocaine ~ fixed(0.525625) # Table S1 Lidocaine 'IIV V, rsd' 72.5 FIX -> 0.725^2
    etalvc_megx ~ fixed(0.868624) # Table S1 MEGX 'IIV V, rsd' 93.2 FIX -> 0.932^2
    etalvc_phenobarbital ~ fixed(0.0441) # Table S1 Phenobarbital 'IIV V, rsd' 21.0 FIX -> 0.210^2
    etalvc_amoxicillin ~ fixed(1.0609) # Table S1 Amoxicillin 'IIV Vc, rsd' 103 FIX -> 1.03^2
    etalvc_benzylpenicillin ~ fixed(0.670761) # Table S1 Benzylpenicillin 'IIV Vc, rsd' 81.9 FIX -> 0.819^2
    etalvp_benzylpenicillin ~ fixed(0.469225) # Table S1 Benzylpenicillin 'IIV Vp, rsd' 68.5 FIX -> 0.685^2
    etalvc_gentamicin ~ fixed(0.638401) # Table S1 Gentamicin 'IIV Vc, rsd' 79.9 FIX -> 0.799^2
    etalvp_gentamicin ~ fixed(0.729316) # Table S1 Gentamicin 'IIV Vp, rsd' 85.4 FIX -> 0.854^2

    # ---- Residual error (Table S1 'SIGMA structure') ----------------------
    # Proportional rsd = sqrt(SIGMA). The additive components were fixed
    # at LLOQ/2 for compounds with below-LLOQ data imputed at LLOQ/2
    # (Methods 'Population pharmacokinetic analysis').
    propSd <- 0.234; label("Morphine proportional residual error (fraction)") # Table S1 Morphine proportional rsd 23.4%
    propSd_m3g <- 0.192; label("M3G proportional residual error (fraction)") # Table S1 M3G proportional rsd 19.2%
    propSd_m6g <- 0.180; label("M6G proportional residual error (fraction)") # Table S1 M6G proportional rsd 18.0%
    propSd_midazolam <- 0.367; label("Midazolam proportional residual error (fraction)") # Table S1 Midazolam proportional rsd 36.7%
    addSd_midazolam <- fixed(0.01); label("Midazolam additive residual error (mg/L)") # Table S1 Midazolam additive 0.01 FIX mg/l
    propSd_1ohm <- 0.298; label("1-hydroxymidazolam proportional residual error (fraction)") # Table S1 OHM proportional rsd 29.8%
    addSd_1ohm <- fixed(0.01); label("1-hydroxymidazolam additive residual error (mg/L)") # Table S1 OHM additive 0.01 FIX mg/l
    propSd_hmg <- 0.272; label("Hydroxymidazolam glucuronide proportional residual error (fraction)") # Table S1 HMG proportional rsd 27.2%
    addSd_hmg <- fixed(0.01); label("Hydroxymidazolam glucuronide additive residual error (mg/L)") # Table S1 HMG additive 0.01 FIX mg/l
    propSd_lidocaine <- 0.223; label("Lidocaine proportional residual error (fraction)") # Table S1 Lidocaine proportional rsd 22.3%
    addSd_lidocaine <- fixed(0.1); label("Lidocaine additive residual error (mg/L)") # Table S1 Lidocaine additive 0.1 FIX mg/l
    propSd_megx <- 0.252; label("MEGX proportional residual error (fraction)") # Table S1 MEGX proportional rsd 25.2%
    addSd_megx <- fixed(0.1); label("MEGX additive residual error (mg/L)") # Table S1 MEGX additive 0.1 FIX mg/l
    propSd_phenobarbital <- 0.0925; label("Phenobarbital proportional residual error (fraction)") # Table S1 Phenobarbital proportional rsd 9.25%
    propSd_amoxicillin <- 0.231; label("Amoxicillin proportional residual error (fraction)") # Table S1 Amoxicillin proportional rsd 23.1%
    propSd_benzylpenicillin <- 0.362; label("Benzylpenicillin proportional residual error (fraction)") # Table S1 Benzylpenicillin proportional rsd 36.2%
    propSd_gentamicin <- 0.248; label("Gentamicin proportional residual error (fraction)") # Table S1 Gentamicin proportional rsd 24.8%
  })
  model({
    # Molecular weights (g/mol). The PharmaCool analyses converted doses
    # and concentrations of every parent-metabolite pair to umol and
    # umol/L (morphine: Favie 2019 PLoS One Methods; lidocaine: Favie
    # 2020 Br J Clin Pharmacol Methods; midazolam: the Figure S3 DV axes
    # reach 10 and 16 for midazolam and HMG, the observed maxima of
    # 3.25 and 8.34 mg/L expressed in umol/L). Parent elimination is
    # therefore converted mole-for-mole into metabolite mass here, so
    # every state is in mg and every concentration in mg/L. Morphine and
    # glucuronide weights are those printed in Favie 2019; the others are
    # standard molecular weights of the free bases.
    mw_morphine <- 285.3
    mw_m3g <- 461.5
    mw_m6g <- 461.5
    mw_midazolam <- 325.8
    mw_1ohm <- 341.8
    mw_hmg <- 517.9
    mw_lidocaine <- 234.3
    mw_megx <- 206.3

    # Covariate terms (Table S1 'Final model')
    pna_h <- PNA * 24 * 30.4375
    fsize_cl <- (WT_BIRTH / 3.5)^e_wt_cl_q
    fsize_v <- (WT_BIRTH / 3.5)^e_wt_vc_vp
    fmat <- (GA^ga_hill / (GA^ga_hill + ga_tm50^ga_hill)) /
      (40^ga_hill / (40^ga_hill + ga_tm50^ga_hill))
    ftemp <- 1 + e_bodytemp_cl * (BODYTEMP - 36.5)
    # Log-normal IIV on the postnatal-age slope, one eta shared by both
    # groups with the intermediate group's SD scaled by 0.765 (correlation
    # fixed to 100%).
    pna_slope_high <- e_pna_cl_high * exp(etae_pna_cl_high)
    pna_slope_int <- e_pna_cl_int * exp(sd_ratio_e_pna_cl_int * etae_pna_cl_high)
    fpna_high <- 1 + pna_slope_high * pna_h
    fpna_int <- 1 + pna_slope_int * pna_h

    # Clearances: intermediate-clearance group (PNA 0.54 %/h and
    # temperature): morphine, midazolam, 1-hydroxymidazolam.
    cl <- exp(lcl + etalcl + sd_ratio_cl * etalcl_common) * fsize_cl * fmat * fpna_int * ftemp
    cl_midazolam <- exp(lcl_midazolam + etalcl_midazolam + sd_ratio_cl_midazolam * etalcl_common) * fsize_cl * fmat * fpna_int * ftemp
    cl_1ohm <- exp(lcl_1ohm + etalcl_1ohm + sd_ratio_cl_1ohm * etalcl_common) * fsize_cl * fmat * fpna_int * ftemp
    # High-clearance group (PNA 1.23 %/h, no temperature effect).
    cl_m3g <- exp(lcl_m3g + etalcl_m3g + sd_ratio_cl_m3g * etalcl_common) * fsize_cl * fmat * fpna_high
    cl_m6g <- exp(lcl_m6g + etalcl_m6g + sd_ratio_cl_m6g * etalcl_common) * fsize_cl * fmat * fpna_high
    cl_hmg <- exp(lcl_hmg + etalcl_hmg + sd_ratio_cl_hmg * etalcl_common) * fsize_cl * fmat * fpna_high
    cl_amoxicillin <- exp(lcl_amoxicillin + etalcl_amoxicillin + sd_ratio_cl_amoxicillin * etalcl_common) * fsize_cl * fmat * fpna_high
    cl_benzylpenicillin <- exp(lcl_benzylpenicillin + etalcl_benzylpenicillin + sd_ratio_cl_benzylpenicillin * etalcl_common) * fsize_cl * fmat * fpna_high
    cl_gentamicin <- exp(lcl_gentamicin + etalcl_gentamicin + sd_ratio_cl_gentamicin * etalcl_common) * fsize_cl * fmat * fpna_high
    # Lidocaine and MEGX: high-clearance group, no common eta (the
    # correlation was not estimable from 28 neonates).
    cl_lidocaine <- exp(lcl_lidocaine + etalcl_lidocaine) * fsize_cl * fmat * fpna_high
    cl_megx <- exp(lcl_megx + etalcl_megx) * fsize_cl * fmat * fpna_high
    # Phenobarbital: no postnatal-age or temperature effect, no common eta.
    cl_phenobarbital <- exp(lcl_phenobarbital + etalcl_phenobarbital) * fsize_cl * fmat

    q_amoxicillin <- exp(lq_amoxicillin) * fsize_cl
    q_benzylpenicillin <- exp(lq_benzylpenicillin) * fsize_cl
    q_gentamicin <- exp(lq_gentamicin) * fsize_cl

    vc <- exp(lvc + etalvc) * fsize_v
    vc_m3g <- exp(lvc_m3g) * fsize_v
    vc_m6g <- exp(lvc_m6g) * fsize_v
    vc_midazolam <- exp(lvc_midazolam + etalvc_midazolam) * fsize_v
    vc_1ohm <- exp(lvc_1ohm) * fsize_v
    vc_hmg <- exp(lvc_hmg + etalvc_hmg) * fsize_v
    vc_lidocaine <- exp(lvc_lidocaine + etalvc_lidocaine) * fsize_v
    vc_megx <- exp(lvc_megx + etalvc_megx) * fsize_v
    vc_phenobarbital <- exp(lvc_phenobarbital + etalvc_phenobarbital) * fsize_v
    vc_amoxicillin <- exp(lvc_amoxicillin + etalvc_amoxicillin) * fsize_v
    vp_amoxicillin <- exp(lvp_amoxicillin) * fsize_v
    vc_benzylpenicillin <- exp(lvc_benzylpenicillin + etalvc_benzylpenicillin) * fsize_v
    vp_benzylpenicillin <- exp(lvp_benzylpenicillin + etalvp_benzylpenicillin) * fsize_v
    vc_gentamicin <- exp(lvc_gentamicin + etalvc_gentamicin) * fsize_v
    vp_gentamicin <- exp(lvp_gentamicin + etalvp_gentamicin) * fsize_v

    # Micro-constants
    kel <- cl / vc
    kel_m3g <- cl_m3g / vc_m3g
    kel_m6g <- cl_m6g / vc_m6g
    kel_midazolam <- cl_midazolam / vc_midazolam
    kel_1ohm <- cl_1ohm / vc_1ohm
    kel_hmg <- cl_hmg / vc_hmg
    kel_lidocaine <- cl_lidocaine / vc_lidocaine
    kel_megx <- cl_megx / vc_megx
    kel_phenobarbital <- cl_phenobarbital / vc_phenobarbital
    kel_amoxicillin <- cl_amoxicillin / vc_amoxicillin
    k12_amoxicillin <- q_amoxicillin / vc_amoxicillin
    k21_amoxicillin <- q_amoxicillin / vp_amoxicillin
    kel_benzylpenicillin <- cl_benzylpenicillin / vc_benzylpenicillin
    k12_benzylpenicillin <- q_benzylpenicillin / vc_benzylpenicillin
    k21_benzylpenicillin <- q_benzylpenicillin / vp_benzylpenicillin
    kel_gentamicin <- cl_gentamicin / vc_gentamicin
    k12_gentamicin <- q_gentamicin / vc_gentamicin
    k21_gentamicin <- q_gentamicin / vp_gentamicin

    # Morphine -> M3G and M6G. Each glucuronide receives the whole morphine
    # elimination flux; the unknown formation fractions are absorbed into
    # the apparent CL/F and V/F (Table S1 footnote).
    d/dt(central) <- -kel * central
    d/dt(central_m3g) <- kel * central * mw_m3g / mw_morphine - kel_m3g * central_m3g
    d/dt(central_m6g) <- kel * central * mw_m6g / mw_morphine - kel_m6g * central_m6g
    # Midazolam -> 1-hydroxymidazolam -> hydroxymidazolam glucuronide
    d/dt(central_midazolam) <- -kel_midazolam * central_midazolam
    d/dt(central_1ohm) <- kel_midazolam * central_midazolam * mw_1ohm / mw_midazolam - kel_1ohm * central_1ohm
    d/dt(central_hmg) <- kel_1ohm * central_1ohm * mw_hmg / mw_1ohm - kel_hmg * central_hmg
    # Lidocaine -> MEGX
    d/dt(central_lidocaine) <- -kel_lidocaine * central_lidocaine
    d/dt(central_megx) <- kel_lidocaine * central_lidocaine * mw_megx / mw_lidocaine - kel_megx * central_megx
    # Phenobarbital (one compartment)
    d/dt(central_phenobarbital) <- -kel_phenobarbital * central_phenobarbital
    # Two-compartment antibiotics
    d/dt(central_amoxicillin) <- -kel_amoxicillin * central_amoxicillin - k12_amoxicillin * central_amoxicillin + k21_amoxicillin * peripheral1_amoxicillin
    d/dt(peripheral1_amoxicillin) <- k12_amoxicillin * central_amoxicillin - k21_amoxicillin * peripheral1_amoxicillin
    d/dt(central_benzylpenicillin) <- -kel_benzylpenicillin * central_benzylpenicillin - k12_benzylpenicillin * central_benzylpenicillin + k21_benzylpenicillin * peripheral1_benzylpenicillin
    d/dt(peripheral1_benzylpenicillin) <- k12_benzylpenicillin * central_benzylpenicillin - k21_benzylpenicillin * peripheral1_benzylpenicillin
    d/dt(central_gentamicin) <- -kel_gentamicin * central_gentamicin - k12_gentamicin * central_gentamicin + k21_gentamicin * peripheral1_gentamicin
    d/dt(peripheral1_gentamicin) <- k12_gentamicin * central_gentamicin - k21_gentamicin * peripheral1_gentamicin

    Cc <- central / vc
    Cc_m3g <- central_m3g / vc_m3g
    Cc_m6g <- central_m6g / vc_m6g
    Cc_midazolam <- central_midazolam / vc_midazolam
    Cc_1ohm <- central_1ohm / vc_1ohm
    Cc_hmg <- central_hmg / vc_hmg
    Cc_lidocaine <- central_lidocaine / vc_lidocaine
    Cc_megx <- central_megx / vc_megx
    Cc_phenobarbital <- central_phenobarbital / vc_phenobarbital
    Cc_amoxicillin <- central_amoxicillin / vc_amoxicillin
    Cc_benzylpenicillin <- central_benzylpenicillin / vc_benzylpenicillin
    Cc_gentamicin <- central_gentamicin / vc_gentamicin

    Cc ~ prop(propSd)
    Cc_m3g ~ prop(propSd_m3g)
    Cc_m6g ~ prop(propSd_m6g)
    Cc_midazolam ~ add(addSd_midazolam) + prop(propSd_midazolam)
    Cc_1ohm ~ add(addSd_1ohm) + prop(propSd_1ohm)
    Cc_hmg ~ add(addSd_hmg) + prop(propSd_hmg)
    Cc_lidocaine ~ add(addSd_lidocaine) + prop(propSd_lidocaine)
    Cc_megx ~ add(addSd_megx) + prop(propSd_megx)
    Cc_phenobarbital ~ prop(propSd_phenobarbital)
    Cc_amoxicillin ~ prop(propSd_amoxicillin)
    Cc_benzylpenicillin ~ prop(propSd_benzylpenicillin)
    Cc_gentamicin ~ prop(propSd_gentamicin)
  })
}
