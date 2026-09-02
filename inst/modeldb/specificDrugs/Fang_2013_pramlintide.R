Fang_2013_pramlintide <- function() {
  description <- paste(
    "Mechanism-based population PK/PD model of pramlintide and postprandial",
    "plasma glucose in adults with type 1 diabetes mellitus (Fang 2013,",
    "a reanalysis of the Colburn 1996 AC137 crossover study). Pramlintide",
    "disposition is two-compartment with zero-order intravenous input and",
    "first-order elimination, a modelled 4 min duration for the nominal",
    "2 min bolus regimen and an absorption lag time. Postprandial glucose is",
    "a two-compartment indirect-response (turnover) system: a 50 g meal",
    "glucose dose enters an intestine compartment by zero-order gastric",
    "emptying of duration Tin and is absorbed into plasma at first-order",
    "rate ka, endogenous hepatic glucose production kin is balanced against",
    "first-order utilisation kout, and glucose distributes to a peripheral",
    "pool. Pramlintide acts by two mechanisms: a dose-independent maximal",
    "prolongation of the gastric-emptying duration (S = 1 + Smax whenever",
    "pramlintide is given) and a concentration-driven inhibitory Imax",
    "function on kin (I = 1 - Imax * Cp / (IC50 + Cp)), the latter active",
    "only from 60 min after dose initiation as coded in the published",
    "NONMEM PD control stream. Glucose distribution volumes and",
    "intercompartmental clearance are fixed to the Silber 2007 literature",
    "values. Two outputs: pramlintide plasma concentration Cc (pmol/L) and",
    "plasma glucose Gc (mg/dL)."
  )
  reference <- paste(
    "Fang J, Landersdorfer CB, Cirincione B, Jusko WJ (2013).",
    "Study Reanalysis Using a Mechanism-Based Pharmacokinetic/Pharmacodynamic",
    "Model of Pramlintide in Subjects with Type 1 Diabetes.",
    "The AAPS Journal 15(1):15-29.",
    "doi:10.1208/s12248-012-9409-7.",
    "PMCID PMC3535104.",
    "Original clinical study: Colburn WA, Gottlieb AB, Koda J, Kolterman OG",
    "(1996) J Clin Pharmacol 36(1):13-24.",
    sep = " "
  )
  vignette <- "Fang_2013_pramlintide"

  # Two dosed compartments: pramlintide into central and the 50 g meal glucose
  # into intestine. The registry only auto-detects depot / central, so both are
  # declared here.
  dosing <- c("central", "intestine")

  units <- list(
    time = "min",
    dosing = paste(
      "pmol of pramlintide into the central compartment;",
      "mg of glucose into the intestine compartment",
      "(the paper meal dose D = 50 g = 50000 mg)"
    ),
    concentration = paste(
      "pmol/L for pramlintide (Cc);",
      "mg/dL for plasma glucose (Gc)"
    )
  )

  # Peripheral glucose pool. Paper symbol Gp (Eq. 7); no canonical
  # compartment name exists for the peripheral limb of a glucose turnover
  # model. Precedent: Bosch_2024_cotadutide_qsp.R uses the same name.
  paper_specific_compartments <- c("glucose_per")

  covariateData <- list(
    ON_TREATMENT = list(
      description = paste(
        "Active pramlintide treatment indicator. 1 = the subject received",
        "pramlintide on this occasion, 0 = matching placebo. Drives the",
        "gastric-emptying prolongation S = 1 + Smax * ON_TREATMENT.",
        "Fang 2013 modelled the effect of pramlintide on the duration of",
        "glucose input as dose independent and maximal at every studied",
        "dose (30, 100 and 300 ug), so the switch is categorical rather",
        "than exposure driven. The inhibition of endogenous glucose",
        "production is exposure driven and needs no indicator."
      ),
      units = "(binary)",
      type = "binary",
      reference_category = "0 (placebo)",
      notes = paste(
        "Paper source column SIT in the NONMEM PD control stream",
        "(Appendix, p. 26): IF (SIT.EQ.0) THEN TVD3 = THETA(2) for placebo,",
        "IF (SIT.EQ.1.AND.GRP.EQ.k) THEN TVD3 = THETA(2) * (1 + SMAX) for",
        "each of the three dose groups k = 1, 2, 3. Time-varying across",
        "occasions in this crossover design: each subject contributed both",
        "placebo and active occasions."
      ),
      source_name = "SIT"
    ),
    OCC = list(
      description = paste(
        "Study occasion index, 1 to 4, following the Fang 2013 PD control",
        "stream coding: 1 = placebo infusion (or placebo infusion day 2),",
        "2 = drug infusion (or placebo infusion day 3), 3 = placebo bolus",
        "(or placebo bolus day 2), 4 = drug bolus (or placebo bolus day 3).",
        "Multiplexes the interoccasion-variability etas on baseline plasma",
        "glucose (four occasions) and on pramlintide clearance (two",
        "occasions, indexed by regimen)."
      ),
      units = "(integer index)",
      type = "count",
      reference_category = NULL,
      notes = paste(
        "The regimen is recoverable from OCC: OCC 1 and 2 are the 120 min",
        "infusion occasions, OCC 3 and 4 are the short (nominally 2 min,",
        "modelled 4 min) bolus occasions. The model uses that mapping for",
        "two purposes. First, the PK control stream (Appendix, p. 25)",
        "indexes interoccasion variability on CL by regimen, IF (OCC.EQ.0)",
        "BOVCL = ETA(5) for the bolus regimen and IF (OCC.EQ.1)",
        "BOVCL = ETA(6) for the infusion regimen with OMEGA BLOCK(1) SAME.",
        "Second, the same control stream assigns the lag time and the",
        "modelled infusion duration only inside IF (REGI.EQ.0), that is",
        "only for the bolus regimen, so the lag time here is multiplied by",
        "the bolus indicator (OCC >= 3). For a single-regimen simulation",
        "hold OCC at 1 or 2 (infusion) or at 3 or 4 (bolus). OCC does not",
        "by itself identify active treatment because subjects randomised to",
        "placebo throughout received placebo on every occasion; use",
        "ON_TREATMENT for that."
      ),
      source_name = "OCC"
    )
  )

  compartmentData <- list(
    central = list(
      analyte = "pramlintide", units = "pmol",
      specimen = "plasma", verified = TRUE
    ),
    peripheral1 = list(
      analyte = "pramlintide", units = "pmol",
      specimen = "plasma", verified = TRUE
    ),
    intestine = list(
      analyte = "glucose", units = "mg",
      specimen = "administration site", verified = TRUE
    ),
    glucose = list(
      analyte = "glucose", units = "mg",
      specimen = "plasma", verified = TRUE
    ),
    glucose_per = list(
      analyte = "glucose", units = "mg",
      specimen = "tissue", verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 25L,
    n_studies = 1L,
    age_range = "mean 29.6 years, SD 6.8 years",
    weight_range = "mean 170.4 lb, SD 17.0 lb (mean 77.3 kg, SD 7.7 kg)",
    sex_female_pct = 0,
    race_ethnicity = c(White = 96, Hispanic = 4),
    disease_state = paste(
      "Type 1 diabetes mellitus diagnosed 2 to 20 years before study entry,",
      "glycosylated haemoglobin 6.1 to 13.0 percent, basal C-peptide below",
      "0.6 ng/mL, stable glycaemic control on insulin for at least 2 weeks."
    ),
    dose_range = paste(
      "Pramlintide 30 ug (group 1), 100 ug (group 2) or 300 ug (group 3),",
      "each given intravenously both as a nominal 2 min bolus and as a",
      "120 min infusion in a two-period crossover, with matching placebo."
    ),
    regions = "United States (two study sites)",
    notes = paste(
      "Demographics from Fang 2013 Methods, Subjects. All subjects were",
      "male. 342 pramlintide concentration observations in 18 individuals",
      "were used for the PK model; individuals receiving placebo only were",
      "excluded from the PK fit because endogenous amylin was negligible or",
      "below the limit of detection. 1028 plasma glucose observations in all",
      "25 individuals, placebo and active, were used for the PD model.",
      "Two subjects per dose group received placebo throughout. Subjects",
      "took their usual pre-breakfast insulin about 30 min before dose",
      "initiation and breakfast was served 30 min after dose initiation;",
      "insulin was not measured and is not in the model. Pramlintide PK",
      "parameters were fixed at the individual empirical Bayes estimates",
      "during the PD estimation step, so the PK and PD random effects below",
      "come from two sequential NONMEM runs rather than one joint fit."
    )
  )

  ini({
    # -----------------------------------------------------------------
    # Pramlintide pharmacokinetics. Fang 2013 Table II, Estimate column.
    # -----------------------------------------------------------------
    lcl <- log(0.955); label("Pramlintide clearance CL (L/min)")                                    # Table II: CL = 0.955 (RSE 12.4 percent)
    lvc <- log(19.1); label("Pramlintide central volume V1 (L)")                                    # Table II: V1 = 19.1 (RSE 13.3 percent)
    lvp <- log(11.0); label("Pramlintide peripheral volume V2 (L)")                                 # Table II: V2 = 11.0 (RSE 11.3 percent)
    lq <- log(0.283); label("Pramlintide distribution clearance CLD (L/min)")                       # Table II: CL_D = 0.283 (RSE 12.7 percent)
    ld1 <- fixed(log(4)); label("Duration of the short intravenous input D1 (min)")                 # Table II: D1 = 4, no RSE; NONMEM PK control stream $THETA (4, FIX); 5 D1
    ltlag <- log(0.430); label("Pramlintide input lag time (min)")                                  # Table II: T_lag = 0.430 (RSE 10.9 percent)

    # -----------------------------------------------------------------
    # Glucose turnover and pramlintide effects.
    # Fang 2013 Table III, Estimate (bootstrap 95 percent CI) column.
    # -----------------------------------------------------------------
    ltin <- log(85.9); label("Duration of zero-order glucose input into the intestine Tin (min)")   # Table III: T_in = 85.9
    lkout <- log(0.0146); label("Glucose utilisation rate constant kout (1/min)")                   # Table III: k_out = 0.0146 (0.00987, 0.0271)
    lka <- log(0.0282); label("Glucose absorption rate constant from intestine ka (1/min)")         # Table III: k_a = 0.0282 (0.0167, 0.0371)
    lfglucose <- log(0.843); label("Oral bioavailability of the 50 g meal glucose dose F (fraction)")  # Table III: F = 0.843 (0.518, 1.24)
    lrbase <- log(161); label("Baseline plasma glucose concentration G0 (mg/dL)")                   # Table III: G_0 = 161 (135, 190)
    limax <- log(0.995); label("Maximum inhibition of endogenous glucose production Imax (fraction)")  # Table III: I_max = 0.995 (0.429, 1.0)
    lic50 <- log(23.8); label("Pramlintide concentration at half of Imax, IC50 (pmol/L)")           # Table III: IC_50 = 23.8 (1.53, 162)
    lsmax <- log(1.26); label("Maximum fractional prolongation of Tin by pramlintide Smax (unitless)")  # Table III: S_max = 1.26 (0.725, 2.98)

    # Glucose disposition, taken from Silber 2007 (Fang 2013 reference 20)
    # and held constant because the pramlintide data could not support
    # their estimation. Table III footnote a.
    lvgc <- fixed(log(9.33)); label("Central volume of distribution of glucose VGc (L; Silber 2007)")        # Table III footnote a: V_Gc = 9.33
    lvgp <- fixed(log(8.56)); label("Peripheral volume of distribution of glucose VGp (L; Silber 2007)")     # Table III footnote a: V_Gp = 8.56
    lqg <- fixed(log(0.442)); label("Intercompartmental clearance of glucose QG (L/min; Silber 2007)")       # Table III footnote a: Q_G = 0.442

    # -----------------------------------------------------------------
    # Random effects. Tables II and III report interindividual (IIV) and
    # interoccasion (IOV) variability as a standard deviation in percent,
    # so the variance entered here is (percent / 100)^2.
    # -----------------------------------------------------------------
    # CL and V1 were correlated; the covariance was estimated.
    # 0.195 * 0.233 * 0.487 = 0.0221268.
    etalcl + etalvc ~ c(
      0.038025,
      0.0221268, 0.054289
    )                                                                                              # Table II: IIV CL 19.5 percent, V1 23.3 percent, Corr CL_V1 0.487; $OMEGA BLOCK(2)
    etalvp ~ 0.381924                                                                              # Table II: IIV V2 61.8 percent
    etald1 ~ 0.138384                                                                              # Table II: IIV D1 37.2 percent

    # IOV on clearance, indexed by regimen. NONMEM PK control stream
    # $OMEGA BLOCK(1) 0.04 ; BOVCL followed by $OMEGA BLOCK(1) SAME, so
    # both occasions share one variance.
    etaiov_lcl_1 ~ 0.011664                                                                        # Table II: IOV CL 10.8 percent; bolus-regimen occasion
    etaiov_lcl_2 ~ fixed(0.011664)                                                                 # SAME as etaiov_lcl_1; infusion-regimen occasion

    etaltin ~ 0.157609                                                                             # Table III: IIV T_in 39.7 percent
    etalkout ~ 0.690561                                                                            # Table III: IIV k_out 83.1 percent (50.2, 122)
    etalfglucose ~ 0.235225                                                                        # Table III: IIV F 48.5 percent (0.478, 132)
    etalrbase ~ 0.147456                                                                           # Table III: IIV G_0 38.4 percent (21.5, 52.7)

    # IOV on baseline glucose over the four study occasions. NONMEM PD
    # control stream $OMEGA BLOCK(1) 0.04 ; BOVGLC1 followed by three
    # $OMEGA BLOCK(1) SAME lines.
    etaiov_lrbase_1 ~ 0.137641                                                                     # Table III: IOV G_0 37.1 percent (30.3, 42.8)
    etaiov_lrbase_2 ~ fixed(0.137641)                                                              # SAME as etaiov_lrbase_1
    etaiov_lrbase_3 ~ fixed(0.137641)                                                              # SAME as etaiov_lrbase_1
    etaiov_lrbase_4 ~ fixed(0.137641)                                                              # SAME as etaiov_lrbase_1

    # Fang 2013 fixed the IIV of CLD, ka, Imax, IC50 and Smax to zero
    # ($OMEGA 0.0 FIX in both control streams; Tables II and III show a
    # dash), so no eta is carried for those parameters.

    # -----------------------------------------------------------------
    # Residual error.
    # -----------------------------------------------------------------
    propSd <- 0.262; label("Proportional residual SD on pramlintide concentration (fraction)")      # Table II: Proportional residual variability 26.2 percent (RSE 3.85)
    expSd_Gc <- 0.157; label("Typical residual SD on log plasma glucose (log scale)")               # Table III: Residual error 15.7 (13.1, 18.1); PD control stream W = THETA(9) * EXP(ETA(13)), Y = LOG(IPRED) + W * EPS(1), $SIGMA 1 FIX
    etaexpSd_Gc ~ 0.156025                                                                         # Table III: IIV on the residual error 39.5 percent (11.3, 61.5)
  })

  model({
    # ---------------------------------------------------------------
    # 1. Occasion-derived indicators.
    #    OCC 1 and 2 are the 120 min infusion occasions, OCC 3 and 4 the
    #    short bolus occasions (see covariateData$OCC).
    # ---------------------------------------------------------------
    occ_bolus <- (OCC >= 3)
    occ_inf <- (OCC <= 2)
    iov_lcl <- occ_bolus * etaiov_lcl_1 + occ_inf * etaiov_lcl_2
    iov_lrbase <- (OCC == 1) * etaiov_lrbase_1 + (OCC == 2) * etaiov_lrbase_2 +
      (OCC == 3) * etaiov_lrbase_3 + (OCC == 4) * etaiov_lrbase_4

    # ---------------------------------------------------------------
    # 2. Individual parameters.
    # ---------------------------------------------------------------
    cl <- exp(lcl + etalcl + iov_lcl)
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp + etalvp)
    q <- exp(lq)
    d1 <- exp(ld1 + etald1)
    tlag <- exp(ltlag)

    tin <- exp(ltin + etaltin)
    kout <- exp(lkout + etalkout)
    ka <- exp(lka)
    fglucose <- exp(lfglucose + etalfglucose)
    rbase <- exp(lrbase + etalrbase + iov_lrbase)
    imax <- exp(limax)
    ic50 <- exp(lic50)
    smax <- exp(lsmax)
    vgc <- exp(lvgc)
    vgp <- exp(lvgp)
    qg <- exp(lqg)

    # ---------------------------------------------------------------
    # 3. Derived constants and drug effects.
    # ---------------------------------------------------------------
    # Glucose amounts are carried in mg and the glucose volumes are
    # reported in L, so the concentration conversion uses 10 dL/L
    # (the published control stream hard-codes 93.3 dL and 85.6 dL).
    vgc_dl <- vgc * 10
    vgp_dl <- vgp * 10
    kcp <- qg / vgc                                                # Eq. 6-7 QG/VGc; control stream 0.0473 1/min
    kpc <- qg / vgp                                                # Eq. 6-7 QG/VGp; control stream 0.0516 1/min
    gc0 <- rbase * vgc_dl                                          # Gc0, baseline central glucose amount (mg)
    kin <- kout * gc0                                              # Eq. 9: kin = kout * Gc0

    # Pramlintide plasma concentration drives the inhibition of kin.
    Cc <- central / vc

    # Eq. 8: I = 1 - Imax * Cp / (IC50 + Cp).
    inh <- imax * Cc / (ic50 + Cc)
    # The published PD control stream gates the inhibition on
    # IF (T.LT.60) THEN N = 0 ELSE N = 1, so kin is unaffected during the
    # first 60 min after dose initiation (30 min before and 30 min after
    # the meal). Eq. 6 as printed omits this gate; the control stream is
    # the implementation that produced the Table III estimates.
    if (t < 60) {
      nsw <- 0
    } else {
      nsw <- 1
    }
    ikin <- 1 - nsw * inh

    # Eq. 4: S = 1 + Smax, applied at every studied dose whenever
    # pramlintide was given and equal to 1 under placebo.
    s <- 1 + smax * ON_TREATMENT

    # ---------------------------------------------------------------
    # 4. ODE system. Compartment order matches the published PD control
    #    stream $MODEL: CENT, PERI, GUT, GLUC, GLUP.
    # ---------------------------------------------------------------
    d/dt(central) <- -(cl / vc) * central -
      (q / vc) * central + (q / vp) * peripheral1
    d/dt(peripheral1) <- (q / vc) * central - (q / vp) * peripheral1
    # Eq. 3: dInt/dt = k0 over [0, Tin * S] - ka * Int, Int(0) = 0.
    # The zero-order input rate k0 = D * F / (Tin * S) of Eq. 5 is
    # produced by the modelled duration and bioavailability below.
    d/dt(intestine) <- -ka * intestine
    # Eq. 6.
    d/dt(glucose) <- kin * ikin + ka * intestine - kout * glucose +
      kpc * glucose_per - kcp * glucose
    # Eq. 7.
    d/dt(glucose_per) <- kcp * glucose - kpc * glucose_per

    # ---------------------------------------------------------------
    # 5. Initial conditions, input durations, bioavailability and lag.
    # ---------------------------------------------------------------
    glucose(0) <- gc0
    # PD control stream A_0(5) = 0.917 * BGLC * 85.6. The 0.917 factor is
    # hard-coded there and is not tabulated; it puts the peripheral pool
    # slightly below the distribution equilibrium implied by Eq. 7.
    glucose_per(0) <- 0.917 * rbase * vgp_dl

    # Modelled 4 min duration for the nominal 2 min bolus regimen; supply
    # the dose with rate = -2 to use it. An explicit dur or rate on the
    # dosing record overrides this, which is how the 120 min infusion
    # regimen is given.
    dur(central) <- d1
    # The published PK control stream sets ALAG1 only inside
    # IF (REGI.EQ.0), that is only for the bolus regimen.
    alag(central) <- tlag * occ_bolus

    dur(intestine) <- tin * s
    f(intestine) <- fglucose

    # ---------------------------------------------------------------
    # 6. Observations.
    # ---------------------------------------------------------------
    Gc <- glucose / vgc_dl
    # PD control stream W = THETA(9) * EXP(ETA(13)) with SIGMA fixed to 1
    # and Y = LOG(IPRED) + W * EPS(1): additive residual on the log scale
    # with interindividual variability on its magnitude.
    w_Gc <- expSd_Gc * exp(etaexpSd_Gc)

    Cc ~ prop(propSd)
    Gc ~ lnorm(w_Gc)
  })
}
