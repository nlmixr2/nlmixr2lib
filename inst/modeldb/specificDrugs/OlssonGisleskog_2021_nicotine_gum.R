OlssonGisleskog_2021_nicotine_gum <- function() {
  description <- paste(
    "Population PK model for nicotine chewing gum (Nicorette classic and",
    "Freshmint/Freshfruit coated gum, 2-6 mg, chewed for 30 min) in",
    "healthy adult smokers (Olsson Gisleskog 2021; 512 subjects, 14",
    "studies). Disposition (three compartments with allometric weight",
    "scaling) and its IIV are fixed to the paper's intravenous model.",
    "Nicotine is released from the gum by a first-order process whose rate",
    "is calculated from the individual amount released, and release stops",
    "at the end of the 0.5-h chewing period. A fraction Frsw (54.7%,",
    "logit-scale IOV) of the dose is swallowed, released at the same rate",
    "in the gut and absorbed through three transit compartments (ktr 5.53",
    "1/h) with oral bioavailability fixed to 39.5%; the remainder passes",
    "from the gum to the buccal cavity after a lag and is absorbed",
    "oromucosally (F = 1) with ka 26.5 1/h (Nicorette classic) or 9.29 1/h",
    "(Freshmint/Freshfruit) at 2 mg, scaled by (dose/2)^0.507. On repeated",
    "dosing a transient time-dependent decrease in oral bioavailability",
    "(75.2% maximum, logit-scale IIV) describes the lower-than-expected",
    "accumulation. Residual pre-study nicotine is a virtual 1-mg bolus into",
    "central at the start of the washout (bioavailability 6.93 mg).",
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
    depot = list(analyte = "nicotine", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "nicotine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "nicotine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "nicotine", units = "mg", specimen = "plasma", verified = TRUE),
    depot_oral = list(analyte = "nicotine", units = "mg", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "nicotine", units = "mg", specimen = "administration site", verified = TRUE),
    transit2 = list(analyte = "nicotine", units = "mg", specimen = "administration site", verified = TRUE),
    transit3 = list(analyte = "nicotine", units = "mg", specimen = "administration site", verified = TRUE),
    depot_buccal = list(analyte = "nicotine", units = "mg", specimen = "administration site", verified = TRUE)
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
    DOSE_NICOTINE_MG = list(
      description = "Nominal nicotine content of the chewing gum",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Labelled 2, 4 or 6 mg (ESM NDOSE). Enters the oromucosal ka as",
        "(DOSE_NICOTINE_MG/2)^0.507 and the release-rate calculation.",
        "Dose records carry the nominal content as amt (ESM ';AMT=NDOSE')."
      ),
      source_name = "NDOSE"
    ),
    DOSE_NICOTINE_RELEASED_MG = list(
      description = "Individual amount of nicotine released from the gum during chewing",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Nominal content minus residual nicotine assayed in the chewed gum",
        "(Sect. 2.2; ESM ADOSE). Sets the release rate constant",
        "Krel = -log(1 - ADOSE/NDOSE)/0.5 so that exactly this amount",
        "leaves the gum in the 0.5-h chewing period (Eq. 1); Krel = 12 1/h",
        "when ADOSE >= NDOSE. The paper reports an average Krel of 2.8 1/h",
        "(Table 5) and that on average 64-79% of the dose is released",
        "(Discussion)."
      ),
      source_name = "ADOSE"
    ),
    FORM_NICOTINE_FRESHMINT = list(
      description = "Freshmint/Freshfruit coated chewing gum indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (Nicorette classic chewing gum)",
      notes = paste(
        "Selects the oromucosal ka at 2 mg: 9.29 1/h (coated",
        "Freshmint/Freshfruit gum) vs 26.5 1/h (Nicorette classic). ESM",
        "FORMFL = 1 when FORM != 2, FORM 2 being Nicorette classic gum."
      ),
      source_name = "FORM (FORM.NE.2)"
    ),
    OCC = list(
      description = "Study period (occasion) index for inter-occasion variability",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "Integer 1-5, the study PERIOD column that indexes the IOV etas",
        "(ESM ETA(PERIOD_FRSW), ETA(PERIOD_F2)). Values outside 1-5",
        "switch the IOV off."
      ),
      source_name = "PERIOD"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 512L,
    n_studies = 14L,
    age_range = "18-50 years",
    age_median = "26 years",
    weight_range = "40.6-108.0 kg",
    weight_median = "71.5 kg",
    sex_female_pct = 48.0,
    race_ethnicity = c(White = 97.1, Asian = 0.6, Black = 0.6, Other = 0.6, Missing = 1.2),
    disease_state = "healthy adult smokers (median 20 cigarettes/day)",
    dose_range = "2, 4 and 6 mg nicotine chewing gum chewed for 30 min, single and repeated doses",
    regions = "Sweden",
    notes = "Olsson Gisleskog 2021 Table 1 (chewing gum row) and Table 2 (Gum column)."
  )

  ini({
    # Disposition fixed to the IV model (ESM gum control stream
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

    # Release and absorption (Table 5, Chewing gum column).
    ltlag <- log(0.0531)
    label("Lag time of release from the gum into the buccal cavity (h)") # Table 5 'Lag time (h)' 0.0531
    lka_classic <- log(26.5)
    label("Oromucosal absorption rate constant, Nicorette classic gum at 2 mg (1/h)") # Table 5 'Ka (h-1)' 26.5 (footnote d); ESM $PK TVKA reference NDOSE/2
    lka_freshmint <- log(9.29)
    label("Oromucosal absorption rate constant, Freshmint/Freshfruit gum at 2 mg (1/h)") # Table 5 'Ka (h-1)' 9.29 (footnote d)
    e_dose_nicotine_mg_ka <- 0.507
    label("Power exponent of (nominal dose/2 mg) on the oromucosal ka (unitless)") # Table 5 'Dose on Ka' 0.507; ESM TVKA = ...*(NDOSE/2)**POWKA
    logitfsw <- logit(0.547)
    label("Fraction of the dose swallowed, Frsw (logit)") # Table 5 'Frsw (%)' 54.7
    lktr <- log(5.53)
    label("Transit rate constant Ktrg of the swallowed fraction (1/h)") # Table 5 'Ktrg (h-1)' 5.53
    lfdepot_buccal <- fixed(log(1))
    label("Bioavailability of the oromucosally absorbed fraction (fraction)") # ESM $PK TVFBUCC = 1; Sect. 2.3.4 'a F of 100%'
    lfdepot_oral <- fixed(log(0.395))
    label("Bioavailability of the swallowed fraction, from the oral model (fraction)") # ESM $PK TVORAL_F = 0.395; Table 4 F study 92NNBT005 39.5%

    # Time-dependent decrease in oral bioavailability, Eq. 2 (Table 5).
    logitfdepot_oral_time_max <- logit(0.752)
    label("Maximal fractional decrease in oral bioavailability, Emax (logit)") # Table 5 'Emax (%)' -75.2 (a decrease); ESM output TH13 0.752 (F_EFF, logit IIV)
    lfdepot_oral_t50 <- log(0.636)
    label("Time after the first dose of 50% onset of the bioavailability decrease, Start (h)") # Table 5 'Start (h)' 0.636
    lfdepot_oral_time_dur <- log(13.9)
    label("Duration of the bioavailability decrease, Duration (h)") # Table 5 'Duration (h)' 13.9
    fdepot_oral_time_hill <- 6.75
    label("Steepness of onset and offset of the bioavailability decrease, pow (unitless)") # Table 5 'pow' 6.75

    lfcentral <- log(6.93)
    label("Pre-washout nicotine dose: bioavailability of the virtual 1-mg central bolus (mg)") # Table 5 'Pre-washout nicotine dose (mg)' 6.93

    # IIV fixed from the IV model (ESM $OMEGA ... FIX).
    etalcl ~ fixed(0.0705245) # ESM $OMEGA 0.0705245 FIX (Table 3 IIV CL 27.0%)
    etalvc + etalvp2 + etalvp ~ fixed(c(
      0.381077,
      -0.230826, 0.450311,
      0, 0.5527, 1.83554
    )) # ESM $OMEGA BLOCK(3) FIX (Table 3)
    # Estimated IIV (Table 5; CV% = 100*sqrt(exp(omega^2)-1); variances from
    # the ESM output).
    etaltlag ~ 0.132 # Table 5 'IIV lag time' 37.5% CV; ESM output OMEGA(5,5) 1.32E-01
    etalfdepot_oral ~ 0.411 # Table 5 'IIV Foral' 71.3% CV; ESM output OMEGA(6,6) 4.11E-01
    etalogitfdepot_oral_time_max ~ 2.67 # Table 5 'IIV Emax' 366 (logit-scale variance printed through the CV formula); ESM output OMEGA(7,7) 2.67E+00
    etalfdepot_oral_t50 ~ 0.275 # Table 5 'IIV start' 56.3% CV; ESM output OMEGA(8,8) 2.75E-01
    etalfdepot_oral_time_dur ~ 0.537 # Table 5 'IIV duration' 84.3% CV; ESM output OMEGA(9,9) 5.37E-01
    etalfcentral ~ 0.549 # Table 5 'IIV pre-washout nicotine dose' 85.5% CV; ESM output OMEGA(10,10) 5.49E-01

    # IOV over study periods 1-5; occasions 2-5 share the occasion-1
    # variance per $OMEGA BLOCK(1) SAME.
    etaiov_fsw_1 ~ 1.04 # Table 5 'IOV Frsw' 1.04 (logit scale); ESM output OMEGA(11,11) 1.04E+00
    etaiov_fsw_2 ~ fixed(1.04) # $OMEGA BLOCK(1) SAME
    etaiov_fsw_3 ~ fixed(1.04) # $OMEGA BLOCK(1) SAME
    etaiov_fsw_4 ~ fixed(1.04) # $OMEGA BLOCK(1) SAME
    etaiov_fsw_5 ~ fixed(1.04) # $OMEGA BLOCK(1) SAME
    etaiov_fcentral_1 ~ 0.164 # Table 5 'IOV pre-washout nicotine dose' 42.2% CV; ESM output OMEGA(16,16) 1.64E-01
    etaiov_fcentral_2 ~ fixed(0.164) # $OMEGA BLOCK(1) SAME
    etaiov_fcentral_3 ~ fixed(0.164) # $OMEGA BLOCK(1) SAME
    etaiov_fcentral_4 ~ fixed(0.164) # $OMEGA BLOCK(1) SAME
    etaiov_fcentral_5 ~ fixed(0.164) # $OMEGA BLOCK(1) SAME

    propSd <- 0.104
    label("Proportional residual error (fraction)") # Table 5 'Proportional residual error (%)' 10.4
    addSd <- 0.157
    label("Additive residual error (ng/mL)") # Table 5 (continued) 'Additive residual error (ng/mL)' 0.157
  })
  model({
    # 1. Occasion indicators (ESM $ABBREVIATED REPLACE ETA(PERIOD_x)).
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    oc5 <- (OCC == 5)
    iov_fsw <- oc1 * etaiov_fsw_1 + oc2 * etaiov_fsw_2 + oc3 * etaiov_fsw_3 +
      oc4 * etaiov_fsw_4 + oc5 * etaiov_fsw_5
    iov_fcentral <- oc1 * etaiov_fcentral_1 + oc2 * etaiov_fcentral_2 +
      oc3 * etaiov_fcentral_3 + oc4 * etaiov_fcentral_4 +
      oc5 * etaiov_fcentral_5

    # 2. Allometric scaling (ESM $PK).
    wt_cl <- (WT / 70)^e_wt_cl_q
    wt_v <- (WT / 70)^e_wt_vc_vp

    # 3. Individual parameters.
    cl <- exp(lcl + etalcl) * wt_cl
    vc <- exp(lvc + etalvc) * wt_v
    q <- exp(lq) * wt_cl
    vp <- exp(lvp + etalvp) * wt_v
    q2 <- exp(lq2) * wt_cl
    vp2 <- exp(lvp2 + etalvp2) * wt_v
    lka_form <- (1 - FORM_NICOTINE_FRESHMINT) * lka_classic + FORM_NICOTINE_FRESHMINT * lka_freshmint
    ka <- exp(lka_form) * (DOSE_NICOTINE_MG / 2)^e_dose_nicotine_mg_ka
    ktr <- exp(lktr)
    tlag <- exp(ltlag + etaltlag)

    # Release rate from the gum, Eq. 1: Krel = -log(1 - ADOSE/NDOSE)/0.5,
    # and 12 1/h when the released amount is not below the nominal dose
    # (ESM $PK KREL).
    frac_released <- DOSE_NICOTINE_RELEASED_MG / DOSE_NICOTINE_MG
    if (frac_released > 0.999999) frac_released <- 0.999999
    krel <- -log(1 - frac_released) / 0.5
    if (DOSE_NICOTINE_RELEASED_MG >= DOSE_NICOTINE_MG) krel <- 12

    # Release stops at the end of the 0.5-h chewing period that starts at
    # the dose (ESM ENDDOSE = DOSETIME + 0.5, KRELFL2 in $DES). The gut
    # dose record has no lag, so its tad() is the time since chewing began.
    tad_chew <- tad(depot_oral)
    if (is.na(tad_chew)) tad_chew <- 1e6
    chewing <- 0
    if (tad_chew <= 0.5) chewing <- 1

    # Fraction swallowed, logit-normal (ESM LTVFRSW, LFRSW).
    logit_fsw <- logitfsw + iov_fsw
    fsw <- expit(logit_fsw)

    # Time-dependent decrease in oral bioavailability, Eq. 2, with T the
    # time since the first dose of the period (ESM TSSP), taken as the time
    # after the first gut dose. 1/(1 + (Start/T)^pow) is Eq. 2's
    # T^pow/(Start^pow + T^pow); T is floored at 1e-6 h before the first
    # dose, where the change is 0 to machine precision.
    tafd_gum <- tafd(depot_oral)
    if (is.na(tafd_gum)) tafd_gum <- 0
    if (tafd_gum < 1e-6) tafd_gum <- 1e-6
    f_t50 <- exp(lfdepot_oral_t50 + etalfdepot_oral_t50)
    f_time_dur <- exp(lfdepot_oral_time_dur + etalfdepot_oral_time_dur)
    f_t_end <- f_t50 + f_time_dur
    f_time_frac <- 1 / (1 + (f_t50 / tafd_gum)^fdepot_oral_time_hill) -
      1 / (1 + (f_t_end / tafd_gum)^fdepot_oral_time_hill)
    logit_f_time_max <- logitfdepot_oral_time_max + etalogitfdepot_oral_time_max
    f_time_max <- expit(logit_f_time_max)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # 4. ODEs (ESM $DES: GUM -> MOUTH at KREL while chewing, MOUTH ->
    #    CENTRAL at KA; GUT -> TRANS1 at KREL while chewing, then TRANS1 ->
    #    TRANS2 -> TRANS3 -> CENTRAL at KTR).
    d/dt(depot) <- -krel * chewing * depot
    d/dt(central) <- ka * depot_buccal + ktr * transit3 -
      kel * central - k12 * central + k21 * peripheral1 -
      k13 * central + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2
    d/dt(depot_oral) <- -krel * chewing * depot_oral
    d/dt(transit1) <- krel * chewing * depot_oral - ktr * transit1
    d/dt(transit2) <- ktr * transit1 - ktr * transit2
    d/dt(transit3) <- ktr * transit2 - ktr * transit3
    d/dt(depot_buccal) <- krel * chewing * depot - ka * depot_buccal

    # 5. Dosing. Each gum is TWO bolus records of amt = the nominal content:
    #    one into depot (the gum in the mouth) and one into depot_oral (the
    #    swallowed share). ESM F1 = (1-FRSW)*BUCCAL_F and
    #    F5 = FRSW*ORAL_F with ORAL_F = 0.395*EXP(ETA(6))*(1-F_EFF). Only
    #    the virtual pre-washout smoking bolus enters central.
    fdepot_oral <- exp(lfdepot_oral + etalfdepot_oral) * (1 - f_time_frac * f_time_max)
    f(depot) <- exp(lfdepot_buccal) * (1 - fsw)
    f(depot_oral) <- fsw * fdepot_oral
    alag(depot) <- tlag
    fcentral <- exp(lfcentral + etalfcentral + iov_fcentral)
    f(central) <- fcentral

    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
