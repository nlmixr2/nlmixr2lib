OlssonGisleskog_2021_nicotine_mouthspray <- function() {
  description <- paste(
    "Population PK model for nicotine oromucosal mouth spray (Nicorette",
    "QuickMist, 1 mg per spray) in healthy adult smokers (Olsson Gisleskog",
    "2021; 201 subjects, 6 studies, 1-4 mg). Disposition (three",
    "compartments with allometric weight scaling) and its IIV are fixed to",
    "the paper's intravenous model. The dose is delivered as a bolus split",
    "between the buccal cavity and the gut: a dose-dependent fraction Frsw",
    "(60.6% at 2 mg, rising with dose, logit-scale IOV) is swallowed and",
    "absorbed through a gut compartment and two transit compartments",
    "(ktr 3.70 1/h) with oral bioavailability fixed to 39.5%; the remainder",
    "is absorbed oromucosally (F = 1) after a 1.4-min lag with ka 15.9 1/h",
    "(buccal spraying) or 86.0 1/h (sublingual spraying). A transient",
    "time-dependent clearance increase (35.7% maximum) describes the",
    "lower-than-expected accumulation on repeated dosing. Residual",
    "pre-study nicotine is a virtual 1-mg bolus into central at the start",
    "of the washout (bioavailability 4.84 mg). IOV over up to six study",
    "periods on ka, Frsw, ktr and the pre-washout dose.",
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
    DOSE_NICOTINE_MG = list(
      description = "Nicotine dose of the mouth-spray administration",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Number of 1-mg sprays in the administration. Enters the fraction",
        "swallowed as (DOSE_NICOTINE_MG/2)^0.0928 (ESM TVFRSW =",
        "MFRSW*(ADOSE/2)**POWFRSW, reference 2 mg). The source column is",
        "the actual dose ADOSE, which for the spray equals the nominal",
        "dose. Must equal the amt of the administration's dose records."
      ),
      source_name = "ADOSE"
    ),
    ROUTE_SUBLINGUAL = list(
      description = "Sublingual (vs buccal) mouth-spray administration indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (spray directed into the buccal cavity, the standard administration)",
      notes = paste(
        "1 = the spray was directed under the tongue. Selects the",
        "oromucosal absorption rate constant (86.0 vs 15.9 1/h; Table 5",
        "footnote c). Source flag: treatment code TMT = 9003 (ESM",
        "FLAG_KA)."
      ),
      source_name = "TMT (TMT.EQ.9003)"
    ),
    OCC = list(
      description = "Study period (occasion) index for inter-occasion variability",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "Integer 1-6. The ESM control stream indexes the IOV etas by the",
        "study PERIOD column (ETA(PERIOD_KA) etc.). Each period starts",
        "with its own virtual pre-washout dose (EVID = 4 reset in the",
        "source data). Values outside 1-6 switch the IOV off."
      ),
      source_name = "PERIOD"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 201L,
    n_studies = 6L,
    age_range = "18-50 years",
    age_median = "27 years",
    weight_range = "49.4-105.6 kg",
    weight_median = "72.6 kg",
    sex_female_pct = 44.3,
    race_ethnicity = c(White = 99.0, Asian = 1.0),
    disease_state = "healthy adult smokers (median 20 cigarettes/day)",
    dose_range = "1, 2, 3 and 4 mg nicotine oromucosal spray, single and repeated doses",
    regions = "Sweden",
    notes = "Olsson Gisleskog 2021 Table 1 (mouth spray row) and Table 2 (Mouth spray column)."
  )

  ini({
    # Disposition fixed to the IV model (ESM mouth-spray control stream
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

    # Absorption (Table 5, Mouth spray column).
    lka_buccal <- log(15.9)
    label("Oromucosal absorption rate constant after buccal spraying (1/h)") # Table 5 'Ka (h-1)' 15.9 (footnote c, buccal)
    lka_sublingual <- log(86.0)
    label("Oromucosal absorption rate constant after sublingual spraying (1/h)") # Table 5 'Ka (h-1)' 86.0 (footnote c, sublingual)
    ltlag <- log(0.0230)
    label("Lag time of oromucosal absorption (h)") # Table 5 'Lag time (h)' 0.0230
    logitfsw <- logit(0.606)
    label("Fraction of the dose swallowed at a 2-mg dose, Frsw (logit)") # Table 5 'Frsw (%)' 60.6 at 2 mg (footnote e); ESM output TH10 0.606
    e_dose_nicotine_mg_fsw <- 0.0928
    label("Power exponent of (dose/2 mg) on the fraction swallowed (unitless)") # ESM output TH11 9.28E-02 (POWFRSW); reproduces Table 5 Frsw 64.6% at 4 mg
    lktr <- log(3.70)
    label("Gut and transit rate constant Ktrg (1/h)") # Table 5 'Ktrg (h-1)' 3.70
    lfdepot_buccal <- fixed(log(1))
    label("Bioavailability of the oromucosally absorbed fraction (fraction)") # ESM $THETA '1 FIX ; F total'; Sect. 2.3.4 'a F of 100%'
    lfdepot_oral <- fixed(log(0.395))
    label("Bioavailability of the swallowed fraction, from the oral model (fraction)") # ESM $PK ORAL_F = 0.395*EXP(ETA(6)); Table 4 F study 92NNBT005 39.5%; Sect. 3.4 'oral F fixed to 40%'

    # Time-dependent clearance change, Eq. 2 (Table 5).
    lcl_t50 <- log(3.37)
    label("Time after the first dose of 50% onset of the clearance increase, Start (h)") # Table 5 'Start (h)' 3.37
    lcl_time_dur <- log(8.78)
    label("Duration of the clearance increase, Duration (h)") # Table 5 'Duration (h)' 8.78
    lcl_time_max <- log(0.357)
    label("Maximal fractional clearance increase, Emax (fraction)") # Table 5 'Emax (%)' 35.7
    lcl_time_hill <- log(4.82)
    label("Steepness of onset and offset of the clearance increase, pow (unitless)") # Table 5 'pow' 4.82

    lfcentral <- log(4.84)
    label("Pre-washout nicotine dose: bioavailability of the virtual 1-mg central bolus (mg)") # Table 5 'Pre-washout nicotine dose (mg)' 4.84

    # IIV fixed from the IV model (ESM $OMEGA ... FIX).
    etalcl ~ fixed(0.0705245) # ESM $OMEGA 0.0705245 FIX (Table 3 IIV CL 27.0%)
    etalvc + etalvp2 + etalvp ~ fixed(c(
      0.381077,
      -0.230826, 0.450311,
      0, 0.5527, 1.83554
    )) # ESM $OMEGA BLOCK(3) FIX (Table 3)
    # Estimated IIV (Table 5; CV% = 100*sqrt(exp(omega^2)-1); variances from
    # the ESM output).
    etaltlag ~ 0.0472 # Table 5 'IIV lag time' 22.0% CV; ESM output OMEGA(5,5) 4.72E-02
    etalfdepot_oral ~ 0.307 # Table 5 'IIV Foral' 60.0% CV; ESM output OMEGA(6,6) 3.07E-01
    etalfcentral ~ 0.509 # Table 5 'IIV pre-washout nicotine dose' 81.5% CV; ESM output OMEGA(7,7) 5.09E-01
    etalcl_time_dur ~ 0.543 # Table 5 'IIV duration' 84.9% CV; ESM output OMEGA(8,8) 5.43E-01

    # IOV over study periods 1-6; occasions 2-6 share the occasion-1
    # variance per $OMEGA BLOCK(1) SAME.
    etaiov_fcentral_1 ~ 0.229 # Table 5 'IOV pre-washout nicotine dose' 50.7% CV; ESM output OMEGA(9,9) 2.29E-01
    etaiov_fcentral_2 ~ fixed(0.229) # $OMEGA BLOCK(1) SAME
    etaiov_fcentral_3 ~ fixed(0.229) # $OMEGA BLOCK(1) SAME
    etaiov_fcentral_4 ~ fixed(0.229) # $OMEGA BLOCK(1) SAME
    etaiov_fcentral_5 ~ fixed(0.229) # $OMEGA BLOCK(1) SAME
    etaiov_fcentral_6 ~ fixed(0.229) # $OMEGA BLOCK(1) SAME
    etaiov_ka_1 ~ 0.653 # Table 5 'IOV Ka' 96.0% CV; ESM output OMEGA(15,15) 6.53E-01
    etaiov_ka_2 ~ fixed(0.653) # $OMEGA BLOCK(1) SAME
    etaiov_ka_3 ~ fixed(0.653) # $OMEGA BLOCK(1) SAME
    etaiov_ka_4 ~ fixed(0.653) # $OMEGA BLOCK(1) SAME
    etaiov_ka_5 ~ fixed(0.653) # $OMEGA BLOCK(1) SAME
    etaiov_ka_6 ~ fixed(0.653) # $OMEGA BLOCK(1) SAME
    etaiov_fsw_1 ~ 0.294 # Table 5 'IOV Frsw' 0.294 (logit scale); ESM output OMEGA(21,21) 2.94E-01
    etaiov_fsw_2 ~ fixed(0.294) # $OMEGA BLOCK(1) SAME
    etaiov_fsw_3 ~ fixed(0.294) # $OMEGA BLOCK(1) SAME
    etaiov_fsw_4 ~ fixed(0.294) # $OMEGA BLOCK(1) SAME
    etaiov_fsw_5 ~ fixed(0.294) # $OMEGA BLOCK(1) SAME
    etaiov_fsw_6 ~ fixed(0.294) # $OMEGA BLOCK(1) SAME
    etaiov_ktr_1 ~ 0.353 # Table 5 'IOV Ktr' 65.1% CV; ESM output OMEGA(27,27) 3.53E-01
    etaiov_ktr_2 ~ fixed(0.353) # $OMEGA BLOCK(1) SAME
    etaiov_ktr_3 ~ fixed(0.353) # $OMEGA BLOCK(1) SAME
    etaiov_ktr_4 ~ fixed(0.353) # $OMEGA BLOCK(1) SAME
    etaiov_ktr_5 ~ fixed(0.353) # $OMEGA BLOCK(1) SAME
    etaiov_ktr_6 ~ fixed(0.353) # $OMEGA BLOCK(1) SAME

    propSd <- 0.101
    label("Proportional residual error (fraction)") # Table 5 'Proportional residual error (%)' 10.1
    addSd <- 0.136
    label("Additive residual error (ng/mL)") # Table 5 (continued) 'Additive residual error (ng/mL)' 0.136
  })
  model({
    # 1. Occasion indicators (ESM $ABBREVIATED REPLACE ETA(PERIOD_x)).
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    oc5 <- (OCC == 5)
    oc6 <- (OCC == 6)
    iov_fcentral <- oc1 * etaiov_fcentral_1 + oc2 * etaiov_fcentral_2 +
      oc3 * etaiov_fcentral_3 + oc4 * etaiov_fcentral_4 +
      oc5 * etaiov_fcentral_5 + oc6 * etaiov_fcentral_6
    iov_ka <- oc1 * etaiov_ka_1 + oc2 * etaiov_ka_2 + oc3 * etaiov_ka_3 +
      oc4 * etaiov_ka_4 + oc5 * etaiov_ka_5 + oc6 * etaiov_ka_6
    iov_fsw <- oc1 * etaiov_fsw_1 + oc2 * etaiov_fsw_2 + oc3 * etaiov_fsw_3 +
      oc4 * etaiov_fsw_4 + oc5 * etaiov_fsw_5 + oc6 * etaiov_fsw_6
    iov_ktr <- oc1 * etaiov_ktr_1 + oc2 * etaiov_ktr_2 + oc3 * etaiov_ktr_3 +
      oc4 * etaiov_ktr_4 + oc5 * etaiov_ktr_5 + oc6 * etaiov_ktr_6

    # 2. Allometric scaling (ESM $PK).
    wt_cl <- (WT / 70)^e_wt_cl_q
    wt_v <- (WT / 70)^e_wt_vc_vp

    # 3. Time-dependent clearance change, Eq. 2, with T the time since the
    #    first dose of the period (ESM TSSP), taken as the time after the
    #    first gut dose (no lag on that compartment). 1/(1 + (Start/T)^pow)
    #    is Eq. 2's T^pow/(Start^pow + T^pow); T is floored at 1e-6 h before
    #    the first dose, where the change is 0 to machine precision.
    tafd_spray <- tafd(depot_oral)
    if (is.na(tafd_spray)) tafd_spray <- 0
    if (tafd_spray < 1e-6) tafd_spray <- 1e-6
    cl_time_max <- exp(lcl_time_max)
    cl_time_hill <- exp(lcl_time_hill)
    cl_t50 <- exp(lcl_t50)
    cl_time_dur <- exp(lcl_time_dur + etalcl_time_dur)
    cl_t_end <- cl_t50 + cl_time_dur
    cl_time_frac <- 1 / (1 + (cl_t50 / tafd_spray)^cl_time_hill) -
      1 / (1 + (cl_t_end / tafd_spray)^cl_time_hill)

    # 4. Individual parameters. cl is the total, time-varying clearance.
    cl <- exp(lcl + etalcl) * wt_cl * (1 + cl_time_max * cl_time_frac)
    vc <- exp(lvc + etalvc) * wt_v
    q <- exp(lq) * wt_cl
    vp <- exp(lvp + etalvp) * wt_v
    q2 <- exp(lq2) * wt_cl
    vp2 <- exp(lvp2 + etalvp2) * wt_v
    lka_route <- (1 - ROUTE_SUBLINGUAL) * lka_buccal + ROUTE_SUBLINGUAL * lka_sublingual
    ka <- exp(lka_route + iov_ka)
    ktr <- exp(lktr + iov_ktr)
    tlag <- exp(ltlag + etaltlag)

    # Fraction swallowed: logit-normal with dose power on the typical value
    # (ESM TVFRSW, LTVFRSW, LFRSW).
    fsw_typ <- expit(logitfsw) * (DOSE_NICOTINE_MG / 2)^e_dose_nicotine_mg_fsw
    logit_fsw <- logit(fsw_typ) + iov_fsw
    fsw <- expit(logit_fsw)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # 5. ODEs (ESM $DES: ABS1 -> CENTRAL at KA; GUT -> TRANS1 -> TRANS2 ->
    #    CENTRAL at KTR2).
    d/dt(depot_buccal) <- -ka * depot_buccal
    d/dt(central) <- ka * depot_buccal + ktr * transit2 -
      kel * central - k12 * central + k21 * peripheral1 -
      k13 * central + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2
    d/dt(depot_oral) <- -ktr * depot_oral
    d/dt(transit1) <- ktr * depot_oral - ktr * transit1
    d/dt(transit2) <- ktr * transit1 - ktr * transit2

    # 6. Dosing. Each spray administration is TWO bolus records of the same
    #    amt: one into depot_buccal and one into depot_oral. The fraction
    #    swallowed and the route bioavailabilities are applied here (ESM
    #    F1 = BUCCAL_F*(1-FRSW), F5 = ORAL_F*FRSW). Only the virtual
    #    pre-washout smoking bolus enters central.
    fdepot_oral <- exp(lfdepot_oral + etalfdepot_oral)
    f(depot_buccal) <- exp(lfdepot_buccal) * (1 - fsw)
    f(depot_oral) <- fdepot_oral * fsw
    alag(depot_buccal) <- tlag
    fcentral <- exp(lfcentral + etalfcentral + iov_fcentral)
    f(central) <- fcentral

    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
