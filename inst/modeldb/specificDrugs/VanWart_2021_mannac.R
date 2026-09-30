VanWart_2021_mannac <- function() {
  description <- paste(
    "Semi-mechanistic joint population pharmacokinetic model of oral",
    "N-acetylmannosamine (ManNAc) and its metabolite N-acetylneuraminic acid",
    "(Neu5Ac, sialic acid) in adults with GNE myopathy (Van Wart 2021). ManNAc:",
    "one-compartment disposition with first-order absorption after a lag time,",
    "relative bioavailability falling with dose as a power function",
    "(F = 1 at 6 g), and a constant endogenous ManNAc production that holds",
    "the pre-dose baseline M0. Neu5Ac: an indirect-response production through",
    "a precursor compartment, both states draining at the Neu5Ac elimination",
    "rate constant kout, with production stimulated linearly by plasma ManNAc.",
    "The stimulation slope rises exponentially with time from SLP0 to SLPSS",
    "(first-order rate kinc), describing the increase in ManNAc-to-Neu5Ac",
    "conversion over the first week of dosing. Inter-occasion variability on",
    "ManNAc clearance over four occasions. No covariates were retained."
  )
  reference <- paste(
    "Van Wart S, Mager DE, Bednasz CJ, Huizing M, Carrillo N. Population",
    "Pharmacokinetic Model of N-acetylmannosamine (ManNAc) and",
    "N-acetylneuraminic acid (Neu5Ac) in Subjects with GNE Myopathy.",
    "Drugs R D. 2021;21(2):189-202. doi:10.1007/s40268-021-00343-6"
  )
  vignette <- "VanWart_2021_mannac"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # ManNAc is carried as an AMOUNT (mg) in depot and central, exactly as in the
  # deposited control stream (Online Resource 6, $MODEL DEPOT / MANNAC, dose in
  # mg). The two Neu5Ac states are CONCENTRATIONS (ng/mL): the stream carries
  # them as mg/L with IPRED = A(3) / (1/1000), and Van Wart 2021 Eqs. 3-4 write
  # them as concentrations with initial condition N0. No Neu5Ac volume is
  # defined anywhere, so none is asserted here.
  compartmentData <- list(
    depot = list(analyte = "ManNAc", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "ManNAc", units = "mg", specimen = "plasma", verified = TRUE),
    precursor1 = list(analyte = "Neu5Ac precursor (Pre-Neu5Ac)", units = "ng/mL", specimen = "plasma", verified = TRUE),
    central_neu5ac = list(analyte = "Neu5Ac", units = "ng/mL", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    DOSE_MANNAC_MG = list(
      description = "Administered ManNAc dose per administration",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Amount of a single administration, not the daily total: a 4 g TID subject carries 4000.",
        "Enters relative bioavailability as the power term (DOSE_MANNAC_MG / 6000)^e_dose_fdepot,",
        "reproducing the control stream line 'IF(DOSEMG.GT.0)F1=1*(DOSEMG/6000)**THETA(11)'",
        "(Online Resource 6). Rows with DOSE_MANNAC_MG = 0 (observation-only rows) leave F at 1,",
        "as in the stream; F only acts on dose rows, so set the column to the dose amount on every",
        "dose row. Studied doses 3, 4, 6 and 10 g."
      ),
      source_name = "DOSEMG"
    ),
    OCC = list(
      description = "Pharmacokinetic occasion for inter-occasion variability on ManNAc clearance",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "Four occasions, from the control stream $PK: OCC = 1 before Day 4 of treatment;",
        "OCC = 2 from Day 4 to before Day 30; OCC = 3 from Day 30 on the twice-daily regimen;",
        "OCC = 4 the three-times-daily (4 g Q8H) extension period at about Day 912. Any other",
        "value (e.g. 0) switches the IOV term off, so a typical-subject simulation can use OCC = 1."
      ),
      source_name = "OCC1-OCC4 (derived in $PK from DAY and TID)"
    )
  )

  # Screened in the formal covariate analysis and NOT retained: none reached
  # the forward-selection criterion (p < 0.01; Results 3.3.3 and Online
  # Resource 5, where the strongest candidate, weight on SLPSS, dropped the
  # objective function by only 4.212). Documentation only.
  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened (power on SLPSS, CL/F, V/F, kout), not retained. Mean 83.5 kg, range 49.3-115 (Table 2)."
    ),
    SEXF = list(
      description = "Female sex",
      units = "(binary)",
      type = "binary",
      notes = "Screened (proportional shift on V/F, kout), not retained. 52.9% female (Table 2)."
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened, not retained. Mean 41.3 years, range 25-65 (Table 2)."
    ),
    HT = list(
      description = "Height",
      units = "cm",
      type = "continuous",
      notes = "Screened, not retained. Mean 173 cm (Table 2)."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Screened, not retained. Mean 27.6 kg/m^2 (Table 2)."
    ),
    BSA = list(
      description = "Body surface area",
      units = "m^2",
      type = "continuous",
      notes = "Screened, not retained. Mean 2.01 m^2 (Table 2)."
    ),
    RACE_ASIAN = list(
      description = "Asian race",
      units = "(binary)",
      type = "binary",
      notes = "Race screened, not retained. 26.5% Asian, 70.6% Caucasian (Table 2)."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/dL",
      type = "continuous",
      notes = "Screened as a hepatic-injury marker, not retained. Mean 3.83 g/dL (Table 2)."
    ),
    CPK = list(
      description = "Creatine kinase",
      units = "U/L",
      type = "continuous",
      notes = "Screened as a hepatic-injury marker, not retained. Mean 236 U/L (Table 2)."
    ),
    CRCL = list(
      description = "Estimated glomerular filtration rate (CKD-EPI cystatin C equation)",
      units = "mL/min",
      type = "continuous",
      notes = "Screened, not retained. Cystatin C rather than creatinine because creatinine is unreliable at the low muscle mass of GNE myopathy (Methods 2.1). Mean 123 mL/min (Table 2)."
    ),
    GNE_MUTATION_TYPE = list(
      description = "GNE domain-mutation type (epimerase/epimerase, epimerase/kinase, kinase/kinase, deletion/kinase)",
      units = "(category)",
      type = "categorical",
      notes = "Screened as a disease-related index, not retained. 67.6% epimerase/kinase, 23.5% kinase/kinase (Table 2)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 34L,
    n_studies = 2L,
    age_range = "25-65 years",
    age_median = "39.5 years",
    weight_range = "49.3-115 kg",
    weight_median = "84.6 kg",
    sex_female_pct = 52.9,
    race_ethnicity = c(Caucasian = 70.6, Asian = 26.5, Other = 2.94),
    disease_state = "GNE myopathy (genetically confirmed biallelic pathogenic GNE variants); normal renal function (mean cystatin C eGFR 123 mL/min).",
    dose_range = paste(
      "Oral ManNAc, fasting. NIH 12-HG-0207 (NCT01634750): single doses of 3 g (n = 6),",
      "6 g (n = 8) or 10 g (n = 8) against placebo. NIH 15-HG-0068 (NCT02346461, n = 12):",
      "3 g or 6 g BID for Days 1-7, then 6 g BID to 30 months, with a 72 h washout after",
      "Day 90 and, after a 5-7 day interruption, 4 g TID (Q8H) on Days 912-914 (n = 8)."
    ),
    regions = "United States (NIH Clinical Center).",
    notes = paste(
      "Demographics from Van Wart 2021 Table 2 and Online Resource 1. 773 ManNAc and 777",
      "Neu5Ac concentrations in the original analysis, plus 72 of each added at the final",
      "stage; 4 ManNAc outliers excluded. Assay range 10-5000 ng/mL (ManNAc) and",
      "25-10,000 ng/mL (Neu5Ac). NONMEM 7.3 FOCE-I."
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # Final estimates: Van Wart 2021 Table 4 ('Final population
    # pharmacokinetic model parameter estimates'). The deposited control
    # stream (Online Resource 6) carries the final STRUCTURE but its $THETA /
    # $OMEGA records are the initial estimates handed to the final run
    # (e.g. VM 510 vs 506, SLP 0.000602 vs 0.000619), so Table 4 governs the
    # values and the stream governs the equations and units.
    # ---------------------------------------------------------------------
    lka <- log(0.256); label("First-order absorption rate constant of ManNAc, ka (1/h)") # Table 4 ka = 0.256 1/h (%SEM 15.2)
    lcl <- log(631); label("Apparent clearance of ManNAc, CLM/F (L/h)") # Table 4 CLM/F = 631 L/h (%SEM 14.8)
    lvc <- log(506); label("Apparent volume of distribution of ManNAc, VM/F (L)") # Table 4 VM/F = 506 L (%SEM 29.4)
    ltlag <- log(0.254); label("Absorption lag time of ManNAc, tlag (h)") # Table 4 tlag = 0.254 h (%SEM 26.4); stream bounds ALAG1 to (0, 0.5)
    lrbase <- log(61.1); label("Endogenous baseline plasma ManNAc concentration, M0 (ng/mL)") # Table 4 M0 = 61.1 ng/mL (%SEM 12.0)
    lrbase_neu5ac <- log(150); label("Endogenous baseline plasma Neu5Ac concentration, N0 (ng/mL)") # Table 4 N0 = 150 ng/mL (%SEM 5.71)
    lkout_neu5ac <- log(0.283); label("First-order Neu5Ac elimination rate constant, also the precursor transfer rate, kout (1/h)") # Table 4 kout = 0.283 1/h (%SEM 5.65)
    lslope0 <- log(0.000619); label("Initial linear stimulation slope of ManNAc-to-Neu5Ac conversion, SLP0 (mL/ng)") # Table 4 SLP0 = 0.000619 (ng/mL)^-1 (%SEM 29.1)
    lslope_ss <- log(0.00334); label("Steady-state linear stimulation slope of ManNAc-to-Neu5Ac conversion, SLPSS (mL/ng)") # Table 4 SLPSS = 0.00334 (ng/mL)^-1 (%SEM 35.0)
    lkinc <- log(0.0287); label("First-order rate constant of the rise from SLP0 to SLPSS, kinc (1/h)") # Table 4 kinc = 0.0287 1/h (%SEM 45.3)

    # Relative bioavailability: F = (dose / 6 g)^slope, F = 1 at 6 g (Table 4
    # 'F for 6 g dose = 1 Fixed'). Check against the Results: (3/6)^-0.405 =
    # 1.32 and (10/6)^-0.405 = 0.81, the printed range of F.
    lfdepot <- fixed(log(1)); label("Relative bioavailability of ManNAc at the 6 g reference dose (unitless)") # Table 4 'F for 6 g dose = 1, Fixed'
    e_dose_fdepot <- -0.405; label("Power exponent on (dose / 6000 mg) for relative bioavailability (unitless)") # Table 4 'F-Dose slope' = -0.405 (%SEM 39.0); stream F1 = (DOSEMG/6000)**THETA(11)

    # IIV: Table 4 prints omega^2 with a '%CV' equal to sqrt(omega^2)
    # (sqrt(0.0697) = 0.264 -> 26.4% CV), so the values are log-scale
    # variances. The kout IIV is FIXED to 0 in the final stream ('IIV on kout
    # ... no longer retained', Results 3.3.5) and is therefore omitted.
    etalka ~ 0.0697 # Table 4 omega2 for ka = 0.0697 (26.4% CV; %SEM 91.4)
    etalcl ~ 0.0636 # Table 4 omega2 for CLM/F = 0.0636 (25.2% CV; %SEM 93.2)
    etalvc ~ 0.120 # Table 4 omega2 for VM/F = 0.120 (34.6% CV; %SEM 170)
    etalrbase ~ 0.0966 # Table 4 omega2 for M0 = 0.0966 (31.1% CV; %SEM 43.5)
    etalrbase_neu5ac ~ 0.0439 # Table 4 omega2 for N0 = 0.0439 (21.0% CV; %SEM 55.4)
    etalslope_ss ~ 0.383 # Table 4 omega2 for SLPSS = 0.383 (61.9% CV; %SEM 130)

    # IOV on CLM/F: one variance shared by four occasions (stream $OMEGA
    # BLOCK(1) IOV_CLM1 followed by three BLOCK(1) SAME records).
    etaiov_cl_1 ~ 0.0580 # Table 4 IOV on CLM/F = 0.0580 (24.1% CV; %SEM 63.6)
    etaiov_cl_2 ~ fixed(0.0580) # SAME as occasion 1 (stream $OMEGA BLOCK(1) SAME)
    etaiov_cl_3 ~ fixed(0.0580) # SAME as occasion 1 (stream $OMEGA BLOCK(1) SAME)
    etaiov_cl_4 ~ fixed(0.0580) # SAME as occasion 1 (stream $OMEGA BLOCK(1) SAME)

    # Residual error: proportional (CCV) on each analyte. The stream's
    # $SIGMA is RV_CCV_M, RV_ADD_M = 0 FIXED, RV_CCV_N, RV_ADD_N = 0 FIXED,
    # so both additive components are fixed to zero. Table 4 prints the two
    # estimated variances as 'sigma2 CCV component 0.102' and 'Additive
    # component 0.0370'; the stream shows the second is the Neu5Ac CCV
    # (initial estimate RV_CCV_N = 0.0365), not an additive term.
    propSd <- 0.3194; label("Proportional residual error, plasma ManNAc (fraction)") # Table 4 sigma2 CCV component = 0.102 (31.9% CV); sqrt(0.102) = 0.3194; stream RV_CCV_M
    propSd_neu5ac <- 0.1924; label("Proportional residual error, plasma Neu5Ac (fraction)") # Table 4 row 'Additive component' = 0.0370 (19.2% CV); sqrt(0.0370) = 0.1924; stream RV_CCV_N
  })

  model({
    # 1. Inter-occasion variability on ManNAc clearance (stream $PK OCC1-OCC4)
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2 + oc3 * etaiov_cl_3 + oc4 * etaiov_cl_4

    # 2. Individual parameters
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl + iov_cl)
    vc <- exp(lvc + etalvc)
    tlag <- exp(ltlag)
    rbase <- exp(lrbase + etalrbase) # M0, ng/mL
    rbase_neu5ac <- exp(lrbase_neu5ac + etalrbase_neu5ac) # N0, ng/mL
    kout_neu5ac <- exp(lkout_neu5ac)
    slope0 <- exp(lslope0)
    slope_ss <- exp(lslope_ss + etalslope_ss)
    kinc <- exp(lkinc)

    # 3. Endogenous ManNAc production (paper k_syn = M0 * CLM/F). M0 is in
    #    ng/mL = ug/L, so M0 / 1000 is mg/L and ksyn is in mg/h (Table 4
    #    prints the typical value as 61.1 * 631 = 38,554 ug/h).
    ksyn <- rbase / 1000 * cl
    kel <- cl / vc

    # 4. Zero-order production into the Neu5Ac precursor, AS RUN. Van Wart
    #    2021 Eq. 5 writes k_pro = kout * N0 / (1 + SLP0 * M0), with M0 in
    #    ng/mL, and Table 4's k_pro = 40.9 ng/mL/h is that expression
    #    (0.283 * 150 / (1 + 0.000619 * 61.1) = 40.90). The fitted control
    #    stream instead evaluates KPRO = (KOUT*N0)/(1+SLP*M0) with
    #    M0 = THETA(4)/1000, i.e. M0 in mg/L while SLP is in mL/ng, so the
    #    stimulation term in k_pro is 1000-fold smaller than in the ODE's
    #    STIM = 1 + SLP * MC (MC in ng/mL). The published estimates were
    #    obtained with the as-run form, which is kept here; the /1000 below
    #    reproduces it. Consequence: k_pro is 42.45 rather than 40.9 ng/mL/h
    #    for the typical subject, and the Neu5Ac baseline drifts about 4% up
    #    from N0 before any dose. The vignette quantifies both forms.
    kpro <- kout_neu5ac * rbase_neu5ac / (1 + slope0 * rbase / 1000)

    # 5. Time-dependent stimulation slope, Eq. 6: SLP = SLPSS - (SLPSS -
    #    SLP0) * exp(-kinc * t). t is time since the start of the subject's
    #    record (NONMEM T in $DES), i.e. time since first dose when dosing
    #    starts at t = 0; it is NOT reset by a washout.
    slope <- slope_ss - (slope_ss - slope0) * exp(-kinc * t)

    # 6. ODEs (Eqs. 1-4; stream $DES). Cc is plasma ManNAc in ng/mL.
    Cc <- 1000 * central / vc
    stim <- 1 + slope * Cc
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ksyn + ka * depot - kel * central
    d/dt(precursor1) <- kpro * stim - kout_neu5ac * precursor1
    d/dt(central_neu5ac) <- kout_neu5ac * precursor1 - kout_neu5ac * central_neu5ac

    # 7. Initial conditions: endogenous steady state before dosing
    #    (stream A_INITIAL(2) = M0*VM with M0 in mg/L; A_INITIAL(3,4) = N0)
    central(0) <- rbase * vc / 1000
    precursor1(0) <- rbase_neu5ac
    central_neu5ac(0) <- rbase_neu5ac

    # 8. Absorption lag and dose-dependent relative bioavailability
    #    (stream: F1 = 1; IF(DOSEMG.GT.0) F1 = (DOSEMG/6000)**THETA(11))
    alag(depot) <- tlag
    fdose <- 1
    if (DOSE_MANNAC_MG > 0) fdose <- (DOSE_MANNAC_MG / 6000)^e_dose_fdepot
    f(depot) <- exp(lfdepot) * fdose

    # 9. Observations: plasma ManNAc (Cc) and plasma Neu5Ac (Cc_neu5ac), both
    #    ng/mL, each with a proportional residual error.
    Cc_neu5ac <- central_neu5ac
    Cc ~ prop(propSd)
    Cc_neu5ac ~ prop(propSd_neu5ac)
  })
}
