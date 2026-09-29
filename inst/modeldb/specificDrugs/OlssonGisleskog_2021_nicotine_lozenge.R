OlssonGisleskog_2021_nicotine_lozenge <- function() {
  description <- paste(
    "Population PK model for nicotine lozenges (Nicorette and NiQuitin, 2",
    "and 4 mg, dissolved in the mouth) in healthy adult smokers (Olsson",
    "Gisleskog 2021; 303 subjects, 6 studies). Disposition (three",
    "compartments with allometric weight scaling) and its IIV are fixed to",
    "the paper's intravenous model. Nicotine is released from the lozenge",
    "by a first-order process (Krel 6.41 1/h, IOV). A fraction Frsw of the",
    "dose is swallowed (68.8% for 4-mg NiQuitin, lower at 2 mg and higher",
    "for Nicorette lozenges; logit-scale IOV), released at the same rate in",
    "the gut and absorbed through three transit compartments (ktr 3.54",
    "1/h, IOV) with oral bioavailability fixed to 39.5%; the remainder",
    "passes after a lag into the buccal cavity and is absorbed",
    "oromucosally (ka 11.6 1/h, F = 1 with IIV). A transient",
    "time-dependent clearance increase (54.4% maximum) describes the",
    "lower-than-expected accumulation on repeated dosing. Residual",
    "pre-study nicotine is a virtual 1-mg bolus into central at the start",
    "of the washout (bioavailability 5.44 mg).",
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
      description = "Nominal nicotine content of the lozenge",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Labelled 2 or 4 mg (ESM NDOSE). A 2-mg lozenge shifts the logit",
        "of the fraction swallowed by -0.217 relative to 4 mg (ESM",
        "IF (NDOSE.EQ.2) EFF_FRSW2 = THETA(14)). Dose records carry the",
        "nominal content as amt."
      ),
      source_name = "NDOSE"
    ),
    FORM_NICOTINE_NICORETTE_LOZENGE = list(
      description = "Nicorette (vs NiQuitin) nicotine lozenge indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (NiQuitin lozenge, GlaxoSmithKline)",
      notes = paste(
        "1 = a Nicorette lozenge (McNeil AB). Adds 0.412 to the logit of",
        "the fraction swallowed (Table 5 'Nicorette on Frsw'). Source",
        "flag: treatment codes TMT 20001, 20003, 20005 and 20006 (ESM",
        "EFF_FRSW1)."
      ),
      source_name = "TMT (20001, 20003, 20005, 20006)"
    ),
    OCC = list(
      description = "Study period (occasion) index for inter-occasion variability",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "Integer 1-5, the study PERIOD column that indexes the IOV etas",
        "(ESM ETA(PERIOD_F2), ETA(PERIOD_FRSW), ETA(PERIOD_KREL),",
        "ETA(PERIOD_KTR)). Values outside 1-5 switch the IOV off."
      ),
      source_name = "PERIOD"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 303L,
    n_studies = 6L,
    age_range = "18-50 years",
    age_median = "28 years",
    weight_range = "52.2-105.0 kg",
    weight_median = "72.8 kg",
    sex_female_pct = 47.2,
    race_ethnicity = c(White = 65.0, Black = 0.3, Other = 0.3, Missing = 34.3),
    disease_state = "healthy adult smokers (median 20 cigarettes/day)",
    dose_range = "2 and 4 mg nicotine lozenge (Nicorette, NiQuitin), single and repeated doses",
    regions = "Sweden",
    notes = "Olsson Gisleskog 2021 Table 1 (lozenge row) and Table 2 (Lozenge column)."
  )

  ini({
    # Disposition fixed to the IV model (ESM lozenge control stream
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

    # Release and absorption (Table 5, Lozenge column).
    lkrel <- log(6.41)
    label("First-order release rate constant from the lozenge Krel (1/h)") # Table 5 'Krel (h-1)' 6.41
    ltlag <- log(0.0437)
    label("Lag time of release from the lozenge into the buccal cavity (h)") # Table 5 'Lag time (h)' 0.0437
    lka <- log(11.6)
    label("Oromucosal absorption rate constant ka (1/h)") # Table 5 'Ka (h-1)' 11.6
    logitfsw <- logit(0.688)
    label("Fraction of the dose swallowed for a 4-mg NiQuitin lozenge, Frsw (logit)") # Table 5 'Frsw (%)' 68.8 at 4 mg (footnote e); ESM output TH10 0.688
    e_form_nicotine_nicorette_lozenge_fsw <- 0.412
    label("Nicorette-lozenge shift in the logit of the fraction swallowed (logit units)") # Table 5 'Nicorette on Frsw' 0.412 (footnote g, additive on logit scale)
    e_dose2mg_fsw <- -0.217
    label("2-mg-lozenge shift in the logit of the fraction swallowed (logit units)") # ESM output TH14 -2.17E-01; reproduces Table 5 Frsw 64.0% at 2 mg
    lktr <- log(3.54)
    label("Transit rate constant Ktrg of the swallowed fraction (1/h)") # Table 5 'Ktrg (h-1)' 3.54
    lfdepot_buccal <- fixed(log(1))
    label("Bioavailability of the oromucosally absorbed fraction (fraction)") # ESM $PK TVFBUCC = 1 with IIV; Sect. 2.3.4 'a F of 100%'
    lfdepot_oral <- fixed(log(0.395))
    label("Bioavailability of the swallowed fraction, from the oral model (fraction)") # ESM $PK ORAL_F = 0.395*EXP(ETA(7)); Table 4 F study 92NNBT005 39.5%

    # Time-dependent clearance change, Eq. 2 (Table 5 and ESM output).
    lcl_t50 <- log(4.86)
    label("Time after the first dose of 50% onset of the clearance increase, Start (h)") # Table 5 'Start (h)' 4.86
    lcl_time_dur <- log(7.18)
    label("Duration of the clearance increase, Duration (h)") # ESM output TH17 7.18E+00 (SE 1.45, RSE 20.2%); Table 5 prints 20.2 in this cell, see vignette Errata
    lcl_time_max <- log(0.544)
    label("Maximal fractional clearance increase, Emax (fraction)") # Table 5 'Emax (%)' 54.4
    lcl_time_hill <- log(3.18)
    label("Steepness of onset and offset of the clearance increase, pow (unitless)") # ESM output TH18 3.18E+00 (SE 0.314, RSE 9.87%); Table 5 prints 9.88 in this cell, see vignette Errata

    lfcentral <- log(5.44)
    label("Pre-washout nicotine dose: bioavailability of the virtual 1-mg central bolus (mg)") # Table 5 'Pre-washout nicotine dose (mg)' 5.44

    # IIV fixed from the IV model (ESM $OMEGA ... FIX).
    etalcl ~ fixed(0.0705245) # ESM $OMEGA 0.0705245 FIX (Table 3 IIV CL 27.0%)
    etalvc + etalvp2 + etalvp ~ fixed(c(
      0.381077,
      -0.230826, 0.450311,
      0, 0.5527, 1.83554
    )) # ESM $OMEGA BLOCK(3) FIX (Table 3)
    # Estimated IIV (Table 5; CV% = 100*sqrt(exp(omega^2)-1); variances from
    # the ESM output).
    etaltlag ~ 0.0991 # Table 5 'IIV lag time' 32.3% CV; ESM output OMEGA(5,5) 9.91E-02
    etalka ~ 0.637 # Table 5 'IIV Ka' 94.4% CV; ESM output OMEGA(6,6) 6.37E-01
    etalfdepot_oral ~ 0.316 # Table 5 'IIV Foral' 61.0% CV; ESM output OMEGA(7,7) 3.16E-01
    etalfcentral ~ 0.320 # Table 5 'IIV pre-washout nicotine dose' 61.4% CV; ESM output OMEGA(8,8) 3.20E-01
    etalfdepot_buccal ~ 0.0533 # Table 5 'IIV Fbuccal' 23.4% CV; ESM output OMEGA(9,9) 5.33E-02
    etalcl_t50 ~ 0.171 # Table 5 'IIV start' 43.2% CV; ESM output OMEGA(10,10) 1.71E-01
    etalcl_time_dur ~ 0.428 # Table 5 'IIV duration' 73.1% CV; ESM output OMEGA(11,11) 4.28E-01

    # IOV over study periods 1-5; occasions 2-5 share the occasion-1
    # variance per $OMEGA BLOCK(1) SAME.
    etaiov_fcentral_1 ~ 0.171 # Table 5 'IOV pre-washout nicotine dose' 43.1% CV; ESM output OMEGA(12,12) 1.71E-01
    etaiov_fcentral_2 ~ fixed(0.171) # $OMEGA BLOCK(1) SAME
    etaiov_fcentral_3 ~ fixed(0.171) # $OMEGA BLOCK(1) SAME
    etaiov_fcentral_4 ~ fixed(0.171) # $OMEGA BLOCK(1) SAME
    etaiov_fcentral_5 ~ fixed(0.171) # $OMEGA BLOCK(1) SAME
    etaiov_fsw_1 ~ 0.245 # Table 5 'IOV Frsw' 0.245 (logit scale); ESM output OMEGA(17,17) 2.45E-01
    etaiov_fsw_2 ~ fixed(0.245) # $OMEGA BLOCK(1) SAME
    etaiov_fsw_3 ~ fixed(0.245) # $OMEGA BLOCK(1) SAME
    etaiov_fsw_4 ~ fixed(0.245) # $OMEGA BLOCK(1) SAME
    etaiov_fsw_5 ~ fixed(0.245) # $OMEGA BLOCK(1) SAME
    etaiov_krel_1 ~ 0.387 # Table 5 'IOV Krel' 68.7% CV; ESM output OMEGA(22,22) 3.87E-01
    etaiov_krel_2 ~ fixed(0.387) # $OMEGA BLOCK(1) SAME
    etaiov_krel_3 ~ fixed(0.387) # $OMEGA BLOCK(1) SAME
    etaiov_krel_4 ~ fixed(0.387) # $OMEGA BLOCK(1) SAME
    etaiov_krel_5 ~ fixed(0.387) # $OMEGA BLOCK(1) SAME
    etaiov_ktr_1 ~ 0.281 # Table 5 'IOV Ktr' 56.9% CV; ESM output OMEGA(27,27) 2.81E-01
    etaiov_ktr_2 ~ fixed(0.281) # $OMEGA BLOCK(1) SAME
    etaiov_ktr_3 ~ fixed(0.281) # $OMEGA BLOCK(1) SAME
    etaiov_ktr_4 ~ fixed(0.281) # $OMEGA BLOCK(1) SAME
    etaiov_ktr_5 ~ fixed(0.281) # $OMEGA BLOCK(1) SAME

    propSd <- 0.107
    label("Proportional residual error (fraction)") # Table 5 'Proportional residual error (%)' 10.7
    addSd <- 0.0851
    label("Additive residual error (ng/mL)") # Table 5 (continued) 'Additive residual error (ng/mL)' 0.0851
  })
  model({
    # 1. Occasion indicators (ESM $ABBREVIATED REPLACE ETA(PERIOD_x)).
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    oc5 <- (OCC == 5)
    iov_fcentral <- oc1 * etaiov_fcentral_1 + oc2 * etaiov_fcentral_2 +
      oc3 * etaiov_fcentral_3 + oc4 * etaiov_fcentral_4 +
      oc5 * etaiov_fcentral_5
    iov_fsw <- oc1 * etaiov_fsw_1 + oc2 * etaiov_fsw_2 + oc3 * etaiov_fsw_3 +
      oc4 * etaiov_fsw_4 + oc5 * etaiov_fsw_5
    iov_krel <- oc1 * etaiov_krel_1 + oc2 * etaiov_krel_2 +
      oc3 * etaiov_krel_3 + oc4 * etaiov_krel_4 + oc5 * etaiov_krel_5
    iov_ktr <- oc1 * etaiov_ktr_1 + oc2 * etaiov_ktr_2 + oc3 * etaiov_ktr_3 +
      oc4 * etaiov_ktr_4 + oc5 * etaiov_ktr_5

    # 2. Allometric scaling (ESM $PK).
    wt_cl <- (WT / 70)^e_wt_cl_q
    wt_v <- (WT / 70)^e_wt_vc_vp

    # 3. Time-dependent clearance change, Eq. 2, with T the time since the
    #    first dose of the period (ESM TSSP), taken as the time after the
    #    first gut dose (no lag on that compartment). 1/(1 + (Start/T)^pow)
    #    is Eq. 2's T^pow/(Start^pow + T^pow); T is floored at 1e-6 h before
    #    the first dose, where the change is 0 to machine precision.
    tafd_loz <- tafd(depot_oral)
    if (is.na(tafd_loz)) tafd_loz <- 0
    if (tafd_loz < 1e-6) tafd_loz <- 1e-6
    cl_time_max <- exp(lcl_time_max)
    cl_time_hill <- exp(lcl_time_hill)
    cl_t50 <- exp(lcl_t50 + etalcl_t50)
    cl_time_dur <- exp(lcl_time_dur + etalcl_time_dur)
    cl_t_end <- cl_t50 + cl_time_dur
    cl_time_frac <- 1 / (1 + (cl_t50 / tafd_loz)^cl_time_hill) -
      1 / (1 + (cl_t_end / tafd_loz)^cl_time_hill)

    # 4. Individual parameters. cl is the total, time-varying clearance.
    cl <- exp(lcl + etalcl) * wt_cl * (1 + cl_time_max * cl_time_frac)
    vc <- exp(lvc + etalvc) * wt_v
    q <- exp(lq) * wt_cl
    vp <- exp(lvp + etalvp) * wt_v
    q2 <- exp(lq2) * wt_cl
    vp2 <- exp(lvp2 + etalvp2) * wt_v
    krel <- exp(lkrel + iov_krel)
    ka <- exp(lka + etalka)
    ktr <- exp(lktr + iov_ktr)
    tlag <- exp(ltlag + etaltlag)

    # Fraction swallowed, logit-normal; reference NiQuitin 4 mg (ESM
    # EFF_FRSW1, EFF_FRSW2, LFRSW).
    dose2mg <- (DOSE_NICOTINE_MG == 2)
    logit_fsw <- logitfsw + e_form_nicotine_nicorette_lozenge_fsw * FORM_NICOTINE_NICORETTE_LOZENGE +
      e_dose2mg_fsw * dose2mg + iov_fsw
    fsw <- expit(logit_fsw)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # 5. ODEs (ESM $DES: LOZENGE -> MOUTH at KREL, MOUTH -> CENTRAL at KA;
    #    GUT -> TRANS1 at KREL, then TRANS1 -> TRANS2 -> TRANS3 -> CENTRAL
    #    at KTR).
    d/dt(depot) <- -krel * depot
    d/dt(central) <- ka * depot_buccal + ktr * transit3 -
      kel * central - k12 * central + k21 * peripheral1 -
      k13 * central + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2
    d/dt(depot_oral) <- -krel * depot_oral
    d/dt(transit1) <- krel * depot_oral - ktr * transit1
    d/dt(transit2) <- ktr * transit1 - ktr * transit2
    d/dt(transit3) <- ktr * transit2 - ktr * transit3
    d/dt(depot_buccal) <- krel * depot - ka * depot_buccal

    # 6. Dosing. Each lozenge is TWO bolus records of amt = the nominal
    #    content: one into depot (the lozenge in the mouth) and one into
    #    depot_oral (the swallowed share). ESM F1 = (1-FRSW)*BUCCAL_F and
    #    F5 = FRSW*ORAL_F. Only the virtual pre-washout smoking bolus
    #    enters central.
    fdepot_buccal <- exp(lfdepot_buccal + etalfdepot_buccal)
    fdepot_oral <- exp(lfdepot_oral + etalfdepot_oral)
    f(depot) <- fdepot_buccal * (1 - fsw)
    f(depot_oral) <- fdepot_oral * fsw
    alag(depot) <- tlag
    fcentral <- exp(lfcentral + etalfcentral + iov_fcentral)
    f(central) <- fcentral

    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
