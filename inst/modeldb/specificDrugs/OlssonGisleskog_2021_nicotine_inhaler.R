OlssonGisleskog_2021_nicotine_inhaler <- function() {
  description <- paste(
    "Population PK model for the nicotine inhaler (Nicorette Inhalator, 10",
    "and 15 mg cartridges, inhaled over 20-min sessions) in healthy adult",
    "smokers (Olsson Gisleskog 2021; 58 subjects, 3 studies). Disposition",
    "(three compartments with allometric weight scaling) and its IIV are",
    "fixed to the paper's intravenous model. The amount released per",
    "session (from the inhaler weight change) enters as a 20-min zero-order",
    "input split between the buccal cavity and the gut: a fraction Frsw",
    "(66.9%, logit-scale IOV) is swallowed and absorbed through a gut",
    "compartment and two transit compartments (ktr 5.26 1/h, IOV) with",
    "oral bioavailability fixed to 39.5%; the remainder is absorbed",
    "oromucosally (ka 0.753 1/h, F = 1 with IIV). Bioavailability of both",
    "routes is 66% higher in study 97NNIN024. A transient time-dependent",
    "clearance increase (23.2% maximum, IIV) describes the",
    "lower-than-expected accumulation on repeated dosing. Residual",
    "pre-study nicotine is a virtual 1-mg bolus into central at the start",
    "of the washout (bioavailability 7.31 mg).",
    sep = " "
  )
  reference <- paste(
    "Olsson Gisleskog PO, Perez Ruixo JJ, Westin A, Hansson AC, Soons PA.",
    "Nicotine Population Pharmacokinetics in Healthy Smokers After",
    "Intravenous, Oral, Buccal and Transdermal Administration.",
    "Clin Pharmacokinet. 2021;60(4):541-561.",
    "doi:10.1007/s40262-020-00960-5",
    sep = " "
  )
  vignette <- "OlssonGisleskog_2021_nicotine"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot_buccal = list(analyte = "nicotine", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "nicotine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "nicotine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "nicotine", units = "mg", specimen = "plasma", verified = TRUE),
    depot_oral = list(analyte = "nicotine", units = "mg", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "nicotine", units = "mg", specimen = "administration site", verified = TRUE),
    transit2 = list(analyte = "nicotine", units = "mg", specimen = "administration site", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric scaling to 70 kg with exponents 0.75 (CL, Q2, Q3) and 1 (volumes), fixed from the IV model.",
      source_name = "WT"
    ),
    STUDY_97NNIN024 = list(
      description = "Inhaler study 97NNIN024 indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (the other two inhaler studies)",
      notes = paste(
        "Multiplies both the buccal and the oral bioavailability by",
        "(1 + 0.660) (Table 5 'F increase, Study 97NNIN024'). Source flag:",
        "treatment codes TMT 28001 and 28002 (ESM FL_STUDY28)."
      ),
      source_name = "TMT (28001, 28002)"
    ),
    OCC = list(
      description = "Study period (occasion) index for inter-occasion variability",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "Integer 1-2, the study PERIOD column that indexes the IOV etas",
        "(ESM ETA(PERIOD_F2), ETA(PERIOD_FRSW), ETA(PERIOD_KTR)). Values",
        "outside 1-2 switch the IOV off."
      ),
      source_name = "PERIOD"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 58L,
    n_studies = 3L,
    age_range = "22-49 years",
    age_median = "29 years",
    weight_range = "43.0-101.0 kg",
    weight_median = "71.0 kg",
    sex_female_pct = 53.4,
    race_ethnicity = c(White = 98.3, Other = 1.7),
    disease_state = "healthy adult smokers (median 16 cigarettes/day)",
    dose_range = "Nicotine inhaler, about 2 mg released per 20-min inhalation session, single and repeated sessions",
    regions = "Sweden",
    notes = "Olsson Gisleskog 2021 Table 1 (inhaler row) and Table 2 (Inhaler column)."
  )

  ini({
    # Disposition fixed to the IV model (ESM inhaler control stream
    # $THETA ... FIX; Table 3).
    lcl <- fixed(log(67.4136))
    label("Clearance CL for a 70-kg subject (L/h)") # ESM $THETA 67.4136 FIX (Table 3 CL 67.4)
    lvc <- fixed(log(117.373))
    label("Central volume V1 for a 70-kg subject (L)") # ESM $THETA 117.373 FIX (Table 3 V1 117)
    lq <- fixed(log(38.615))
    label("Inter-compartmental flow Q2 to peripheral1 for a 70-kg subject (L/h)") # ESM $THETA 38.615 FIX (Table 3 Q2 38.6)
    lvp <- fixed(log(130.372))
    label("Peripheral volume V2 for a 70-kg subject (L)") # ESM $THETA 130.372 FIX (Table 3 V2 130)
    lq2 <- fixed(log(216.29))
    label("Inter-compartmental flow Q3 to peripheral2 for a 70-kg subject (L/h)") # ESM $THETA 216.29 FIX (Table 3 Q3 216)
    lvp2 <- fixed(log(53.4189))
    label("Peripheral volume V3 for a 70-kg subject (L)") # ESM $THETA 53.4189 FIX (Table 3 V3 53.4)
    e_wt_cl_q <- fixed(0.75)
    label("Allometric exponent on CL, Q2 and Q3 (unitless)") # Sect. 2.3.1
    e_wt_vc_vp <- fixed(1)
    label("Allometric exponent on V1, V2 and V3 (unitless)") # Sect. 2.3.1

    # Absorption (Table 5, Inhaler column).
    lka <- log(0.753)
    label("Oromucosal absorption rate constant ka (1/h)") # Table 5 'Ka (h-1)' 0.753
    logitfsw <- logit(0.669)
    label("Fraction of the dose swallowed, Frsw (logit)") # Table 5 'Frsw (%)' 66.9
    lktr <- log(5.26)
    label("Gut and transit rate constant Ktrg (1/h)") # Table 5 'Ktrg (h-1)' 5.26
    lfdepot_buccal <- fixed(log(1))
    label("Bioavailability of the oromucosally absorbed fraction (fraction)") # ESM $PK TVFBUCC = 1 with IIV; Sect. 2.3.4 'a F of 100%'
    lfdepot_oral <- fixed(log(0.395))
    label("Bioavailability of the swallowed fraction, from the oral model (fraction)") # ESM $PK ORAL_F = 0.395*EXP(ETA(6)); Table 4 F study 92NNBT005 39.5%
    e_study_97nnin024_f <- 0.660
    label("Fractional increase in buccal and oral bioavailability in study 97NNIN024 (fraction)") # Table 5 'F increase, Study 97NNIN024 (%)' 66.0

    # Time-dependent clearance change, Eq. 2 (Table 5).
    lcl_t50 <- log(2.06)
    label("Time after the first dose of 50% onset of the clearance increase, Start (h)") # Table 5 'Start (h)' 2.06
    lcl_time_dur <- log(4.47)
    label("Duration of the clearance increase, Duration (h)") # Table 5 'Duration (h)' 4.47
    lcl_time_max <- log(0.232)
    label("Maximal fractional clearance increase, Emax (fraction)") # Table 5 'Emax (%)' 23.2; ESM ICL_EFF = THETA(14)*EXP(ETA(10))
    lcl_time_hill <- fixed(log(80))
    label("Steepness of onset and offset of the clearance increase, pow (unitless)") # Table 5 'pow' '80 (fixed)'; ESM $THETA 80 FIX

    lfcentral <- log(7.31)
    label("Pre-washout nicotine dose: bioavailability of the virtual 1-mg central bolus (mg)") # Table 5 'Pre-washout nicotine dose (mg)' 7.31

    # IIV fixed from the IV model (ESM $OMEGA ... FIX).
    etalcl ~ fixed(0.0705245) # ESM $OMEGA 0.0705245 FIX (Table 3 IIV CL 27.0%)
    etalvc + etalvp2 + etalvp ~ fixed(c(
      0.381077,
      -0.230826, 0.450311,
      0, 0.5527, 1.83554
    )) # ESM $OMEGA BLOCK(3) FIX (Table 3)
    # Estimated IIV (Table 5; CV% = 100*sqrt(exp(omega^2)-1); variances from
    # the ESM output). The control stream's 'IIV_start' is 0 FIX and is
    # omitted here (a zero variance adds nothing and makes OMEGA singular).
    etalfdepot_buccal ~ 0.116 # Table 5 'IIV Fbuccal' 35.1% CV; ESM output OMEGA(5,5) 1.16E-01
    etalfdepot_oral ~ 0.0811 # Table 5 'IIV Foral' 29.1% CV; ESM output OMEGA(6,6) 8.11E-02
    etalfcentral ~ 0.140 # Table 5 'IIV pre-washout nicotine dose' 38.7% CV; ESM output OMEGA(7,7) 1.40E-01
    etalcl_time_dur ~ 0.0840 # Table 5 'IIV duration' 29.6% CV; ESM output OMEGA(9,9) 8.40E-02
    etalcl_time_max ~ 0.316 # Table 5 'IIV Emax' 60.9% CV; ESM output OMEGA(10,10) 3.16E-01

    # IOV over study periods 1-2; occasion 2 shares the occasion-1 variance
    # per $OMEGA BLOCK(1) SAME.
    etaiov_fcentral_1 ~ 0.0579 # Table 5 'IOV pre-washout nicotine dose' 24.4% CV; ESM output OMEGA(11,11) 5.79E-02
    etaiov_fcentral_2 ~ fixed(0.0579) # $OMEGA BLOCK(1) SAME
    etaiov_fsw_1 ~ 0.483 # Table 5 'IOV Frsw' 0.483 (logit scale); ESM output OMEGA(13,13) 4.83E-01
    etaiov_fsw_2 ~ fixed(0.483) # $OMEGA BLOCK(1) SAME
    etaiov_ktr_1 ~ 0.246 # Table 5 'IOV Ktr' 52.8% CV; ESM output OMEGA(15,15) 2.46E-01
    etaiov_ktr_2 ~ fixed(0.246) # $OMEGA BLOCK(1) SAME

    propSd <- 0.0761
    label("Proportional residual error (fraction)") # Table 5 'Proportional residual error (%)' 7.61
    addSd <- 0.179
    label("Additive residual error (ng/mL)") # Table 5 (continued) 'Additive residual error (ng/mL)' 0.179
  })
  model({
    # 1. Occasion indicators (ESM $ABBREVIATED REPLACE ETA(PERIOD_x)).
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    iov_fcentral <- oc1 * etaiov_fcentral_1 + oc2 * etaiov_fcentral_2
    iov_fsw <- oc1 * etaiov_fsw_1 + oc2 * etaiov_fsw_2
    iov_ktr <- oc1 * etaiov_ktr_1 + oc2 * etaiov_ktr_2

    # 2. Allometric scaling (ESM $PK).
    wt_cl <- (WT / 70)^e_wt_cl_q
    wt_v <- (WT / 70)^e_wt_vc_vp

    # 3. Time-dependent clearance change, Eq. 2, with T the time since the
    #    first dose of the period (ESM TSSP), taken as the time after the
    #    first gut dose. 1/(1 + (Start/T)^pow) is Eq. 2's
    #    T^pow/(Start^pow + T^pow) and avoids overflowing T^80; T is
    #    floored at 1e-6 h before the first dose, where the change is 0.
    tafd_inh <- tafd(depot_oral)
    if (is.na(tafd_inh)) tafd_inh <- 0
    if (tafd_inh < 1e-6) tafd_inh <- 1e-6
    cl_time_hill <- exp(lcl_time_hill)
    cl_t50 <- exp(lcl_t50)
    cl_time_dur <- exp(lcl_time_dur + etalcl_time_dur)
    cl_t_end <- cl_t50 + cl_time_dur
    cl_time_frac <- 1 / (1 + (cl_t50 / tafd_inh)^cl_time_hill) -
      1 / (1 + (cl_t_end / tafd_inh)^cl_time_hill)
    cl_time_max <- exp(lcl_time_max + etalcl_time_max)

    # 4. Individual parameters. cl is the total, time-varying clearance
    #    (ESM TVCL = THETA(1)*WTCL*(1+CL_EFF)).
    cl <- exp(lcl + etalcl) * wt_cl * (1 + cl_time_max * cl_time_frac)
    vc <- exp(lvc + etalvc) * wt_v
    q <- exp(lq) * wt_cl
    vp <- exp(lvp + etalvp) * wt_v
    q2 <- exp(lq2) * wt_cl
    vp2 <- exp(lvp2 + etalvp2) * wt_v
    ka <- exp(lka)
    ktr <- exp(lktr + iov_ktr)

    logit_fsw <- logitfsw + iov_fsw
    fsw <- expit(logit_fsw)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # 5. ODEs (ESM $DES: MOUTH -> CENTRAL at KA; GUT -> TRANS1 -> TRANS2 ->
    #    CENTRAL at KTR).
    d/dt(depot_buccal) <- -ka * depot_buccal
    d/dt(central) <- ka * depot_buccal + ktr * transit2 -
      kel * central - k12 * central + k21 * peripheral1 -
      k13 * central + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2
    d/dt(depot_oral) <- -ktr * depot_oral
    d/dt(transit1) <- ktr * depot_oral - ktr * transit1
    d/dt(transit2) <- ktr * transit1 - ktr * transit2

    # 6. Dosing. Each inhalation session is TWO zero-order records of amt =
    #    the released amount, one into depot_buccal and one into
    #    depot_oral, each at a fixed rate of amt / (20 min). As in NONMEM,
    #    bioavailability then scales the amount at that fixed rate, so the
    #    effective input duration is F * 20 min. ESM F1 = (1-FRSW)*BUCCAL_F
    #    and F5 = FRSW*ORAL_F. Only the virtual pre-washout smoking bolus
    #    enters central.
    f_study <- 1 + e_study_97nnin024_f * STUDY_97NNIN024
    fdepot_buccal <- exp(lfdepot_buccal + etalfdepot_buccal) * f_study
    fdepot_oral <- exp(lfdepot_oral + etalfdepot_oral) * f_study
    f(depot_buccal) <- fdepot_buccal * (1 - fsw)
    f(depot_oral) <- fdepot_oral * fsw
    fcentral <- exp(lfcentral + etalfcentral + iov_fcentral)
    f(central) <- fcentral

    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
