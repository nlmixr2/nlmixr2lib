Hong_2015_moxifloxacin <- function() {
  description <- "Sequential population PK + PD (QT-interval) model for single-dose oral moxifloxacin (400 mg or 800 mg, Avelox tablets) in healthy adult Korean male volunteers (Hong 2015): a two-compartment first-order absorption PK model with a lag time and a dose-dependent absorption rate constant (different Ka for 400 mg vs 800 mg), followed by an individually corrected QT-interval PD model that adds two mixed-effect cosine circadian components (24 h and 6 h), a first-order-decaying placebo (water-intake) effect, and an Emax drug effect on QT prolongation."
  reference <- paste(
    "Hong T, Han S, Lee J, Jeon S, Park GJ, Park WS,",
    "Lim KS, Chung JY, Yu KS, Yim DS.",
    "Pharmacokinetic-pharmacodynamic analysis to evaluate the effect",
    "of moxifloxacin on QT interval prolongation in healthy Korean male",
    "subjects. Drug Des Devel Ther. 2015 Feb 26;9:1233-1245.",
    "doi:10.2147/DDDT.S79772. PMID 25750523.",
    sep = " "
  )
  vignette <- "Hong_2015_moxifloxacin"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    DOSE = list(
      description        = "Per-record administered moxifloxacin dose level used by the dose-dependent absorption term",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Hong 2015 estimated separate absorption rate constants for the 400 mg and 800 mg dose levels (Methods 'PK model' paragraph 2; Table 3 Ka1 = 16.7 1/h, Ka2 = 1.90 1/h). The packaged model uses a binary switch keyed on DOSE >= 600 mg: subjects with DOSE < 600 mg use the 400 mg ka and subjects with DOSE >= 600 mg use the 800 mg ka. For simulation set DOSE per subject (and per record if the dose changes within a subject) to the administered moxifloxacin dose in mg. Placebo simulations (where the model evaluates the circadian + placebo terms only) can use DOSE = 400 (the reference) since the value affects only the absorption rate.",
      source_name        = "DOSE"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 33L,
    n_studies      = 3L,
    age_range      = "20-40 years (mean 26.4, SD 4.8)",
    age_median     = "26 years (mean 26.4 years)",
    weight_range   = "approximately 60-78 kg (mean 68.3, SD 6.3)",
    weight_median  = "68.3 kg (mean)",
    sex_female_pct = 0,
    race_ethnicity = c(Korean = 100),
    disease_state  = "Healthy adult male Korean volunteers",
    dose_range     = "Single oral dose of moxifloxacin 400 mg or 800 mg as Avelox tablets with 240 mL water, in a three-way (placebo / 400 mg / 800 mg) William's-square crossover with 1-week washout between periods.",
    regions        = "Republic of Korea (three clinical-trial centres: Seoul St Mary's Hospital, Seoul National University Hospital, and Seoul National University Bundang Hospital).",
    notes          = "Demographics from Hong 2015 Table 1. Pooled dataset from three TQT studies (NCT01756521); 38 subjects enrolled, 33 completed and contributed PK data (660 plasma moxifloxacin concentrations). ECG / QT data from 27 subjects (6 lost to data archive errors): 810 baseline + 513 placebo + 540 drug-effect QT interval observations. No PK or PD covariate (age, height, weight, lean body mass, genotype) was retained in the final model."
  )

  ini({
    # =========================================================================
    # PK model (Hong 2015 Table 3) - two-compartment first-order absorption
    # with a lag time and dose-dependent absorption rate constant. Parameters
    # are apparent (CL/F, V2/F, Q/F, V3/F); the source assay is plasma
    # moxifloxacin (ng/mL) and the dose unit is mg, so the model returns Cc
    # in ng/mL via Cc = 1000 * central / vc with vc in L and dose in mg.
    # =========================================================================
    lka       <- log(16.7);  label("Absorption rate constant Ka1 at the reference 400 mg dose (1/h)")                                                   # Hong 2015 Table 3 Ka1 = 16.7 (%RSE 45.5)
    e_dose_ka <- log(1.90 / 16.7); label("Log-scale shift on ka for the 800 mg dose vs the 400 mg reference (unitless; ka_800 = exp(lka + e_dose_ka))") # Hong 2015 Table 3 Ka2 = 1.90 (%RSE 17.5); encoded relative to Ka1 so the source-paper Ka2 = exp(log(16.7) + log(1.90/16.7)) = 1.90 1/h
    lcl       <- log(11.8);  label("Apparent clearance CL/F (L/h)")                                                                                     # Hong 2015 Table 3 CL/F = 11.8 (%RSE 6.29)
    lvc       <- log(173);   label("Apparent central volume of distribution V2/F (L)")                                                                  # Hong 2015 Table 3 V2 = 173 (%RSE 3.17)
    lq        <- log(5.62);  label("Apparent inter-compartmental clearance Q/F (L/h)")                                                                  # Hong 2015 Table 3 Q = 5.62 (%RSE 19.2)
    lvp       <- log(47.1);  label("Apparent peripheral volume of distribution V3/F (L)")                                                               # Hong 2015 Table 3 V3 = 47.1 (%RSE 19.0)
    ltlag     <- log(0.46);  label("Absorption lag time alag (h)")                                                                                       # Hong 2015 Table 3 alag = 0.46 (%RSE 17.1)

    # IIV on CL/F and V2/F with reported correlation rho_cl-v2 = 0.92; IIV on
    # ka. IIV on Q, V3, alag was reported as 'not estimated' and is omitted.
    # omega^2 = log(1 + CV^2). Off-diagonal: cov = rho * sqrt(var_cl * var_vc).
    #   var_cl = log(1 + 0.270^2) = 0.07037
    #   var_vc = log(1 + 0.174^2) = 0.02983
    #   cov    = 0.92 * sqrt(0.07037 * 0.02983) = 0.04215
    etalcl + etalvc ~ c(0.07037, 0.04215, 0.02983)  # Hong 2015 Table 3 omega_CL = 27.0%, omega_V2 = 17.4%, rho_CL-V2 = 0.92
    etalka          ~ 0.6521                          # Hong 2015 Table 3 omega_Ka = 95.9% CV; log(1 + 0.959^2) = 0.6521

    propSd <- 0.140; label("Proportional residual error on plasma moxifloxacin concentration (fraction)")  # Hong 2015 Table 3 sigma_prop = 14.0% (%RSE 11.1)

    # =========================================================================
    # PD model - baseline + circadian QTc (Hong 2015 Table 4). The
    # individually-corrected QT interval (QTcI) is the modelled observation;
    # the heart-rate correction factor a is included for source-trace
    # completeness but is not required by the forward simulation (the
    # observation is QTcI, not QT). Circadian variation is the sum of two
    # cosines with periods 24 h and 6 h.
    # =========================================================================
    lqtcm   <- log(406);    label("Individual mesor (typical baseline) of corrected QT interval (ms)")              # Hong 2015 Table 4 QTcm = 406 (%RSE 0.69)
    la      <- log(0.26);   label("Individual heart-rate correction factor a in QTcI = QT * RR^a (unitless)")        # Hong 2015 Table 4 a = 0.26 (%RSE 5.45); used for QT to QTcI correction in the source dataset and retained for source-trace completeness
    lam24   <- log(0.0090); label("Amplitude of the 24-h circadian cosine, as a fraction of the mesor (unitless)")  # Hong 2015 Table 4 aM24 = 0.0090 (%RSE 19.0)
    lam6    <- log(0.0045); label("Amplitude of the 6-h circadian cosine, as a fraction of the mesor (unitless)")    # Hong 2015 Table 4 aM6 = 0.0045 (%RSE 16.7)
    ac24    <- 1.11;        label("Acrophase (phase shift) of the 24-h circadian cosine (h)")                       # Hong 2015 Table 4 AC24 = 1.11 (%RSE 53.4); time = 0 is 08:00 local time
    ac6     <- 4.65;        label("Acrophase (phase shift) of the 6-h circadian cosine (h)")                        # Hong 2015 Table 4 AC6 = 4.65 (%RSE 3.33)

    # =========================================================================
    # PD model - placebo effect (Hong 2015 Table 5 placebo block). Modeled
    # as a first-order-decaying effect that starts at PE at the time of the
    # 240 mL water intake and decays with rate constant k as time-after-dose
    # progresses. PE is reported as a negative number; the BSV is encoded as
    # additive on the linear scale because the source paper applies the etas
    # for the circadian / placebo model additively (Hong 2015 Results
    # 'Individualized baseline and circadian rhythm model' paragraph 1).
    # =========================================================================
    pe      <- -5.31;       label("Maximal value of the placebo (water-intake) effect on QTc (ms)")  # Hong 2015 Table 5 placebo block PE = -5.31 (%RSE 24.5)
    lkpl    <- log(0.06);   label("Elimination rate constant of the placebo effect (1/h)")            # Hong 2015 Table 5 placebo block k = 0.06 (%RSE 47.2)

    # =========================================================================
    # PD model - drug effect (Hong 2015 Table 5 drug block). The drug effect
    # is added on top of baseline + circadian + placebo:
    # E(C) = Emax * Cc / (EC50 + Cc) with Cc in ng/mL and Emax in ms.
    # =========================================================================
    lemax   <- log(34.7);   label("Maximum drug effect on QTc prolongation Emax (ms)")               # Hong 2015 Table 5 drug block Emax = 34.7 (%RSE 19.5)
    lec50   <- log(3920);   label("Drug concentration producing 50% of the maximum effect EC50 (ng/mL)") # Hong 2015 Table 5 drug block EC50 = 3920 (%RSE 29.3)

    # =========================================================================
    # PD IIV (Hong 2015 Tables 4 + 5).
    # - QTcm BSV 3.82% CV - log-normal multiplicative.
    # - a   BSV 19.6% CV  - log-normal multiplicative.
    # - aM24 BSV 93.0% CV - log-normal multiplicative on the fractional amplitude.
    # - aM6  IIV not estimated.
    # - AC24 reported as 206% (Hong 2015 Table 4); the source paper applies
    #   the circadian etas additively (linear-scale eta on ac24 in hours)
    #   so this implementation encodes the variance as (2.06 * 1.11)^2 =
    #   5.232 h^2 in additive form. See vignette Assumptions and deviations.
    # - AC6 IIV not estimated.
    # - PE BSV 96.9% CV - encoded as additive SD = 0.969 * |5.31| = 5.143 ms
    #   (variance 26.45 ms^2) because PE is negative and additive variability
    #   is the natural form (Hong 2015 places + eta_i on the circadian /
    #   placebo parameters).
    # - k  BSV 118% CV - log-normal multiplicative on the rate constant.
    # - Emax BSV 32.9% CV - log-normal multiplicative.
    # - EC50 IIV not estimated.
    # The QTcm interoccasional variability (0.88% CV; Table 4) is not encoded
    # because nlmixr2lib has no idiomatic encoding for IOV in distributed
    # models; see vignette Assumptions and deviations.
    # =========================================================================
    etalqtcm  ~ 0.001459       # 3.82% CV;  log(1 + 0.0382^2) = 0.001459
    etala     ~ 0.037808       # 19.6% CV;  log(1 + 0.196^2)  = 0.037808
    etalam24  ~ 0.622848       # 93.0% CV;  log(1 + 0.930^2)  = 0.622848
    etaac24   ~ 5.231836       # 206% (additive on hours); var = (2.06 * 1.11)^2 = 5.231836
    etape     ~ 26.453569      # 96.9% (additive on ms);  var = (0.969 * 5.31)^2 = 26.45357
    etalkpl   ~ 0.872446       # 118% CV;   log(1 + 1.18^2)   = 0.872446
    etalemax  ~ 0.102829       # 32.9% CV;  log(1 + 0.329^2)  = 0.102829

    propSd_QTc <- 0.0152; label("Proportional residual error on QTc (fraction)")  # Hong 2015 Table 5 drug-effect block sigma_prop = 1.52% (%RSE 39.4); the most-encompassing of the three reported residuals (baseline 1.35%, placebo 1.32%, drug 1.52%)
  })

  model({
    # -----------------------------------------------------------------------
    # 1. Dose-dependent absorption rate constant.
    # -----------------------------------------------------------------------
    # Hong 2015 Methods 'PK model' paragraph 2: 'Different absorption
    # characteristics per dose were assumed, so a different parameter was
    # estimated by dose in the absorption model.' Threshold 600 mg is the
    # midpoint of the two studied dose levels; subjects with DOSE < 600 mg
    # use Ka1 (= 16.7 1/h at 400 mg), subjects with DOSE >= 600 mg use Ka2
    # (= 1.90 1/h at 800 mg). The shift is encoded on the log scale via
    # e_dose_ka so that ka_400 = exp(lka) and ka_800 = exp(lka + e_dose_ka).
    is_dose_hi <- (DOSE >= 600)
    ka <- exp(lka + e_dose_ka * is_dose_hi + etalka)

    # -----------------------------------------------------------------------
    # 2. Individual PK parameters (no retained covariates in the final model).
    # -----------------------------------------------------------------------
    cl   <- exp(lcl + etalcl)
    vc   <- exp(lvc + etalvc)
    q    <- exp(lq)
    vp   <- exp(lvp)
    alag_dose <- exp(ltlag)

    # -----------------------------------------------------------------------
    # 3. Micro-constants and ODE system. Two-compartment disposition with
    #    first-order absorption from a depot compartment that is delayed by
    #    the absorption lag time alag.
    # -----------------------------------------------------------------------
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    alag(depot) <- alag_dose

    # -----------------------------------------------------------------------
    # 4. Plasma moxifloxacin concentration in ng/mL. Dose in mg and vc in L
    #    give central/vc in mg/L = ug/mL; multiplying by 1000 converts to
    #    ng/mL, the bioanalytical assay unit (LLOQ = 100 ng/mL).
    # -----------------------------------------------------------------------
    Cc <- 1000 * central / vc

    # -----------------------------------------------------------------------
    # 5. Individual PD parameters. The heart-rate correction factor a_hr is
    #    derived here for source-trace completeness; it is the individual
    #    exponent used by Hong 2015 to convert observed QT into QTcI on the
    #    raw data set, and is not required by the forward simulation of
    #    QTcI (the modelled observation).
    # -----------------------------------------------------------------------
    qtcm   <- exp(lqtcm + etalqtcm)
    a_hr   <- exp(la + etala)
    am24   <- exp(lam24 + etalam24)
    am6    <- exp(lam6)
    ac24_i <- ac24 + etaac24
    pe_i   <- pe + etape
    kpl    <- exp(lkpl + etalkpl)
    emax   <- exp(lemax + etalemax)
    ec50   <- exp(lec50)

    # -----------------------------------------------------------------------
    # 6. Baseline + circadian QTc (Hong 2015 Methods 'Circadian rhythm
    #    model'). The mesor QTcm is modulated by a 24-h cosine and a 6-h
    #    cosine; amplitudes are fractional (multiplicative of QTcm).
    #    pi_const = 3.141592653589793 (rxode2 has no built-in pi).
    # -----------------------------------------------------------------------
    pi_const <- 3.141592653589793

    qtb_c <- qtcm * (1
                     + am24 * cos(2 * pi_const * (t - ac24_i) / 24)
                     + am6  * cos(2 * pi_const * (t - ac6)    / 6))

    # -----------------------------------------------------------------------
    # 7. Placebo effect (Hong 2015 Methods 'Placebo effect model'). PE
    #    decays first-order from the time of placebo / drug intake. tad()
    #    is rxode2 time-since-most-recent-dose in the model's time unit
    #    (hours here). The placebo effect represents the response to the
    #    240 mL water intake that accompanies every dose event (active or
    #    placebo).
    # -----------------------------------------------------------------------
    pe_t <- pe_i * exp(-kpl * tad())

    # -----------------------------------------------------------------------
    # 8. Drug effect (Hong 2015 Methods 'Drug effect model'). Emax response
    #    on plasma concentration Cc (ng/mL).
    # -----------------------------------------------------------------------
    e_drug <- emax * Cc / (ec50 + Cc)

    # -----------------------------------------------------------------------
    # 9. Total individually-corrected QT interval observation (ms).
    # -----------------------------------------------------------------------
    QTc <- qtb_c + pe_t + e_drug

    Cc  ~ prop(propSd)
    QTc ~ prop(propSd_QTc)
  })
}
