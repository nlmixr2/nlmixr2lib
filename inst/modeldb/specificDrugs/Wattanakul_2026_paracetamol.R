Wattanakul_2026_paracetamol <- function() {
  description <- "Two-compartment population PK model with first-order absorption for oral paracetamol (acetaminophen) 1 g every 6 h for 72 h in Malaysian patients aged 5 years and older with predominantly non-severe Plasmodium knowlesi malaria (PACKNOW trial; Wattanakul 2026). Relative bioavailability declines exponentially from 1 at the start of treatment to a minimum of 0.220 with a 10.2 h half-life, and carries inter-occasion variability across four treatment days. Allometric body-weight scaling (exponents 0.75 on CL/F and Q/F, 1 on Vc/F and Vp/F, reference 58 kg) and a linear age effect on ka. The structural parameters were estimated with an informative NONMEM $PRIOR on THETA from the Wattanakul 2016 falciparum-malaria model (modellib('Wattanakul_2016_paracetamol')); the values here are the posterior estimates reported for the knowlesi cohort."
  reference <- paste(
    "Wattanakul T, Cooper DJ, Plewes K, Grigg MJ, Rajahram GS, William T,",
    "Dondorp AM, Edstein MD, Birrell GW, Anstey NM, Hoglund RM, Barber BE,",
    "Tarning J. Population Pharmacokinetics and Pharmacodynamics of",
    "Paracetamol in Malaysian Patients With Plasmodium knowlesi Malaria.",
    "CPT Pharmacometrics Syst Pharmacol. 2026;15:e70283.",
    "doi:10.1002/psp4.70283.",
    "Prior PK model: Wattanakul T et al. Malaria Journal 2016;15:244,",
    "doi:10.1186/s12936-016-1283-9; see modellib('Wattanakul_2016_paracetamol').",
    sep = " "
  )
  vignette <- "Wattanakul_2026_paracetamol"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed enrollment value. Allometric scaling (WT/58)^0.75 on CL/F and Q/F and (WT/58)^1 on Vc/F and Vp/F. The reference 58 kg is the value printed in the supplement PK control stream ('TVCL = THETA(1)*(WT/58)**0.75'); the main text describes it as the study-population median.",
      source_name = "WT"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed enrollment value. Linear effect on ka centred at 35 years: ka = 0.917 * (1 - 0.0187 * (AGE - 35)) (supplement PK control stream 'KAAGE = (1 + THETA(9)*(AGE - 35))'). The study enrolled patients aged 5 years and older; the paracetamol-arm median (IQR) was 35 (26-49) years. The linear factor reaches zero at AGE = 88.5 years.",
      source_name = "AGE"
    ),
    OCC = list(
      description = "Treatment-day occasion index for inter-occasion variability on relative bioavailability.",
      units = "(integer)",
      type = "categorical",
      reference_category = NULL,
      notes = "Values 1, 2, 3 and 4 are treatment days 0, 1, 2 and 3 (Wattanakul 2026 Methods 2.2, 'between days random variability, that is, day 0, day 1, day 2, and day 3'). Decomposed inside model() into indicators oc1 .. oc4 that multiplex the four IOV etas, exactly as the supplement control stream builds IOVF1 from ETA(7)-ETA(10). Time-varying: set OCC on each dose record to the day that dose falls in (OCC = floor(time / 24) + 1 for a regimen starting at time 0). Source column OCCN.",
      source_name = "OCCN"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "paracetamol", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "paracetamol", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "paracetamol", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 183L,
    n_studies = 1L,
    n_observations = "2435 plasma paracetamol concentrations (453 intensive, 1982 sparse; 81.4% troughs); 2.0% below the LLOQ of 0.050 mg/L were omitted.",
    age_range = "5 years and older (inclusion criterion); paracetamol-arm IQR 26-49 years",
    age_median = "35 years",
    weight_range = "Paracetamol-arm IQR 53.5-68.4 kg",
    weight_median = "61.0 kg",
    sex_female_pct = 12.6,
    race_ethnicity = "Malaysian (Sabah); not further broken down.",
    disease_state = "Hospitalised patients with PCR-confirmed Plasmodium knowlesi monoinfection and recent fever, enrolled within 18 h of starting artesunate and/or artemether-lumefantrine; 7.1% met WHO research criteria for severe knowlesi malaria and 26.8% had acute kidney injury (KDIGO stage 1-3) at enrollment.",
    dose_range = "Paracetamol 1 g orally every 6 h for 72 h (directly observed).",
    regions = "Sabah, Malaysia (one tertiary and three district hospitals).",
    notes = "PACKNOW trial (NCT03056391), paracetamol arm only. Demographics from Wattanakul 2026 Table 1 (paracetamol column). Patients with measurable paracetamol at enrollment had the central compartment initialised from the observed enrollment concentration ('A_0(2) = INIT_DV*S2'), which this model does not reproduce; supply any pre-study amount as an initial condition when needed."
  )

  ini({
    # Wattanakul 2026 Table 2 (final PK model; posterior estimates under an
    # informative NONMEM $PRIOR on THETA from Wattanakul 2016). The same
    # values are the $THETA block of the supplement PK control stream.
    lka <- log(0.917)
    label("First-order absorption rate constant at AGE = 35 years (1/h)") # Table 2: ka = 0.917 1/h (RSE 18.3%)
    lcl <- log(9.80)
    label("Apparent elimination clearance CL/F at WT = 58 kg (L/h)") # Table 2: CL/F = 9.80 L/h (RSE 11.2%)
    lvc <- log(53.7)
    label("Apparent central volume Vc/F at WT = 58 kg (L)") # Table 2: Vc/F = 53.7 L (RSE 4.1%)
    lq <- log(10.7)
    label("Apparent intercompartmental clearance Q/F at WT = 58 kg (L/h)") # Table 2: Q/F = 10.7 L/h (RSE 15.2%)
    lvp <- log(33.2)
    label("Apparent peripheral volume Vp/F at WT = 58 kg (L)") # Table 2: Vp/F = 33.2 L (RSE 7.1%)
    lfdepot <- fixed(log(1))
    label("Relative oral bioavailability at time zero (fraction)") # Table 2: F = 1 Fixed; stream '(1) FIX ; 6. TVF1'

    # Time-dependent decline of relative bioavailability (Methods Eq. 3):
    # F(t) = F * (exp(-k * t) * (1 - FMIN) + FMIN), k = log(2) / F50
    lfdepot_thalf <- log(10.2)
    label("Half-life of the exponential decline in relative bioavailability, F50 (h)") # Table 2: F50 = 10.2 h (RSE 23.5%); stream 'KF1 = LOG(2)/HLF1'
    fdepot_min <- 0.220
    label("Minimum relative bioavailability reached by the time-dependent decline (fraction of F at time zero)") # Table 2: FMIN = 0.220 (RSE 16.6%)

    # Covariate effects
    e_age_ka <- -0.0187
    label("Linear effect of age on ka (fractional change per year from 35 years)") # Table 2: theta_AGE on ka = -0.0187 (RSE 8.1%); stream 'KAAGE = (1 + THETA(9)*(AGE - 35))'
    e_wt_cl_q <- fixed(0.75)
    label("Allometric exponent on CL/F and Q/F (unitless)") # Methods Eq. 5: n = 0.75 for all clearance parameters
    e_wt_vc_vp <- fixed(1)
    label("Allometric exponent on Vc/F and Vp/F (unitless)") # Methods Eq. 5: n = 1 for volume parameters

    # IIV (exponential). Variances are the supplement control-stream $OMEGA
    # values; each reproduces the Table 2 %CV through
    # CV = 100 * sqrt(exp(omega^2) - 1). IIV on CL and Q is '0 FIX' in the
    # stream and absent from Table 2, so no eta is carried for them.
    etalka ~ 2.13 # stream $OMEGA 5 IIV_KA = 2.13; sqrt(exp(2.13)-1) = 272% CV = Table 2 'IIV ka 272'
    etalvc ~ 0.338 # stream $OMEGA 2 IIV_V = 0.338; sqrt(exp(0.338)-1) = 63.4% CV = Table 2 'IIV VC 63.4'
    etalvp ~ 2.77 # stream $OMEGA 4 IIV_VP = 2.77; sqrt(exp(2.77)-1) = 387% CV = Table 2 'IIV VP 387'
    etalfdepot ~ 0.0662 # stream $OMEGA 6 IIV_F1 = 0.0662; sqrt(exp(0.0662)-1) = 26.2% CV = Table 2 'IIV F 26.2'

    # Inter-occasion variability on F across treatment days 0-3. The stream
    # writes $OMEGA BLOCK(1) 0.251 followed by three BLOCK(1) SAME, i.e. one
    # shared variance; each occasion gets its own eta with the variance fixed
    # to the shared value after occasion 1.
    etaiov_fdepot_1 ~ 0.251 # stream $OMEGA 7 IOV_F1_OCC1 = 0.251; sqrt(exp(0.251)-1) = 53.4% CV = Table 2 'IOV F 53.4'
    etaiov_fdepot_2 ~ fixed(0.251) # stream $OMEGA 8 'BLOCK (1) SAME'
    etaiov_fdepot_3 ~ fixed(0.251) # stream $OMEGA 9 'BLOCK (1) SAME'
    etaiov_fdepot_4 ~ fixed(0.251) # stream $OMEGA 10 'BLOCK (1) SAME'

    # Residual error: additive on the log scale. Table 2 reports
    # sigma = 0.290 as a VARIANCE (Table 2 abbreviations: 'variance of the
    # residual variability incorporated as an additive error on the
    # logarithmic scale'); stream '$SIGMA 0.29' and 'W = SQRT(SIGMA(1,1))'.
    expSd <- sqrt(0.290)
    label("Residual SD of additive error on log-concentration (log mg/L)") # Table 2: sigma = 0.290 (variance; RSE 8.3%)
  })

  model({
    # Occasion indicators and IOV on F (stream: IOVF1 = OCC1*ETA(7) + ...)
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    iov_fdepot <- oc1 * etaiov_fdepot_1 + oc2 * etaiov_fdepot_2 +
      oc3 * etaiov_fdepot_3 + oc4 * etaiov_fdepot_4

    # Individual parameters
    ka <- exp(lka + etalka) * (1 + e_age_ka * (AGE - 35))
    cl <- exp(lcl) * (WT / 58)^e_wt_cl_q
    vc <- exp(lvc + etalvc) * (WT / 58)^e_wt_vc_vp
    q <- exp(lq) * (WT / 58)^e_wt_cl_q
    vp <- exp(lvp + etalvp) * (WT / 58)^e_wt_vc_vp

    # Exponential decline of relative bioavailability with time since the
    # start of treatment (stream 'TEFFF1 = EXP(-KF1*TIME)*(1-MINF1)+MINF1').
    # Time zero is the first dose of the 72 h regimen.
    fdepot_thalf <- exp(lfdepot_thalf)
    kfdepot <- log(2) / fdepot_thalf
    fdepot_time <- exp(-kfdepot * t) * (1 - fdepot_min) + fdepot_min
    fdepot <- exp(lfdepot + etalfdepot + iov_fdepot) * fdepot_time

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    f(depot) <- fdepot

    Cc <- central / vc
    Cc ~ lnorm(expSd)
  })
}
