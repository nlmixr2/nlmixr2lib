OlssonGisleskog_2021_nicotine_transdermal <- function() {
  description <- paste(
    "Population PK model for transdermal nicotine patches (Nicorette patch",
    "5-15 mg/16 h and Nicorette Invisipatch 10-25 mg/16 h, worn 16 or 24",
    "h) in healthy adult smokers (Olsson Gisleskog 2021; 73 subjects, 3",
    "studies). Disposition (three compartments with allometric weight",
    "scaling) and its IIV are fixed to the paper's intravenous model.",
    "Nicotine leaves the patch by two parallel pathways: a fraction Fr1",
    "(40.0% Nicorette patch, 71.9% Invisipatch; logit-scale IIV) by a",
    "first-order release (Krel 0.146 1/h) that runs from application (after",
    "a 0.53-h lag for Invisipatch) for a fraction Frdur1 of 16 h (44.5% and",
    "96.2%), scaled so that the whole fraction is delivered in that window,",
    "and the remainder by a zero-order release from 4.06 h after",
    "application until patch removal. Released nicotine is absorbed through",
    "three transit compartments (ktr 3.62 1/h) with transdermal",
    "bioavailability 75.8% (IOV). Clearance is 11.6% higher from 25 h after",
    "the first application. Residual pre-study nicotine is a virtual 1-mg",
    "bolus into central at the start of the washout (bioavailability 3.82",
    "mg).",
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
    depot_td = list(analyte = "nicotine", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "nicotine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "nicotine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "nicotine", units = "mg", specimen = "plasma", verified = TRUE),
    transit1 = list(analyte = "nicotine", units = "mg", specimen = "administration site", verified = TRUE),
    transit2 = list(analyte = "nicotine", units = "mg", specimen = "administration site", verified = TRUE),
    transit3 = list(analyte = "nicotine", units = "mg", specimen = "administration site", verified = TRUE)
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
    FORM_NICOTINE_INVISIPATCH = list(
      description = "Nicorette Invisipatch (vs Nicorette patch) indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (Nicorette patch 5-15 mg/16 h)",
      notes = paste(
        "1 = Nicorette Invisipatch (NNTP). Selects the first-order release",
        "fraction Fr1 (71.9% vs 40.0%), its duration fraction Frdur1",
        "(96.2% vs 44.5% of 16 h) and a 0.53-h lag on the first-order",
        "pathway (0 for the Nicorette patch). ESM FORMFL = 1 when",
        "FORM = 8 (Nicorette patch)."
      ),
      source_name = "FORM (FORM.NE.8)"
    ),
    T_PATCH_WEAR = list(
      description = "Patch application (wear) duration",
      units = "h",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time from patch application to removal: 16 h normally, 16 or 24",
        "h in study A6431108 (Sect. 2.1). Sets the zero-order release",
        "duration to T_PATCH_WEAR minus the 4.06-h zero-order lag so that",
        "release ends at removal (ESM D5 = DUR - ALAG5). Supply it on the",
        "zero-order dose record."
      ),
      source_name = "DUR (data column DUR2 in the deposited dataset)"
    ),
    OCC = list(
      description = "Study period (occasion) index for inter-occasion variability",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "Integer 1-6, the study PERIOD column that indexes the IOV etas",
        "(ESM ETA(PERIOD_F2), ETA(PERIOD_FTOT)); all daily applications",
        "within a period share one occasion. Values outside 1-6 switch the",
        "IOV off."
      ),
      source_name = "PERIOD"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 73L,
    n_studies = 3L,
    age_range = "19-50 years",
    age_median = "24 years",
    weight_range = "43.2-112.8 kg",
    weight_median = "71.2 kg",
    sex_female_pct = 38.4,
    race_ethnicity = c(White = 100),
    disease_state = "healthy adult smokers (median 18 cigarettes/day)",
    dose_range = "5, 10, 15, 25 and 37 mg nicotine patches applied for 16 h (24 h in one arm), single and three daily applications",
    regions = "Sweden",
    notes = "Olsson Gisleskog 2021 Table 1 (transdermal patch row) and Table 2 (Transdermal column)."
  )

  ini({
    # Disposition fixed to the IV model (ESM transdermal control stream
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

    # Release and absorption (Table 6).
    lkrel <- log(0.146)
    label("First-order release rate constant from the patch Krel (1/h)") # Table 6 'Krel (h-1)' 0.146
    lktr <- log(3.62)
    label("Skin transit rate constant Ktrs (1/h)") # Table 6 'Ktrs (h-1)' 3.62
    logitffo_nicorette <- logit(0.400)
    label("Fraction released by the first-order pathway, Nicorette patch, Fr1 (logit)") # Table 6 'Fr1 Nicorette (%)' 40.0
    logitffo_invisipatch <- logit(0.719)
    label("Fraction released by the first-order pathway, Invisipatch, Fr1 (logit)") # Table 6 'Fr1 Invisipatch (%)' 71.9
    logitfrdur_nicorette <- logit(0.445)
    label("Duration of first-order release as a fraction of 16 h, Nicorette patch, Frdur1 (logit)") # Table 6 'Frdur1 Nicorette (%)' 44.5
    logitfrdur_invisipatch <- logit(0.962)
    label("Duration of first-order release as a fraction of 16 h, Invisipatch, Frdur1 (logit)") # Table 6 'Frdur1 Invisipatch (%)' 96.2
    lfdepot <- log(0.758)
    label("Absolute transdermal bioavailability F (fraction)") # Table 6 'F (%)' 75.8
    ltlag_invisipatch <- log(0.53)
    label("Lag time of the first-order release pathway, Invisipatch (h)") # Table 6 'Lag time1 (h)' 0.53 (ESM ALAG1 = THETA(14)*(1-FORMFL))
    ltlag2 <- log(4.06)
    label("Lag time of the zero-order release pathway, both patches (h)") # Table 6 'Lag time2 (h)' 4.06 (ESM ALAG5)
    lcl_time_max <- log(0.116)
    label("Fractional clearance increase from 25 h after the first application, CLch24 (fraction)") # Table 6 'CLch24 (%)' 11.6; ESM TVCL*(1+TFLAG*THETA(19)), TFLAG = TSSP > 25

    lfcentral <- log(3.82)
    label("Pre-washout nicotine dose: bioavailability of the virtual 1-mg central bolus (mg)") # Table 6 'Pre-washout nicotine dose (mg)' 3.82

    # IIV fixed from the IV model (ESM $OMEGA ... FIX).
    etalcl ~ fixed(0.0705245) # ESM $OMEGA 0.0705245 FIX (Table 3 IIV CL 27.0%)
    etalvc + etalvp2 + etalvp ~ fixed(c(
      0.381077,
      -0.230826, 0.450311,
      0, 0.5527, 1.83554
    )) # ESM $OMEGA BLOCK(3) FIX (Table 3)
    # Estimated IIV (Table 6; CV% = 100*sqrt(exp(omega^2)-1); variances from
    # the ESM output).
    etalkrel ~ 0.219 # Table 6 'IIV Krel' 49.5% CV; ESM output OMEGA(5,5) 2.19E-01
    etalktr ~ 0.835 # Table 6 'IIV KTR' 114% CV; ESM output OMEGA(6,6) 8.35E-01
    etalogitfrdur ~ 0.498 # Table 6 'IIV Frdur1' 0.498 (logit scale); ESM output OMEGA(7,7) 4.98E-01
    etalogitffo ~ 0.238 # Table 6 'IIV Fr1' 0.238 (logit scale); ESM output OMEGA(8,8) 2.38E-01
    etalfcentral ~ 0.906 # Table 6 'IIV pre-washout nicotine dose' 121% CV; ESM output OMEGA(9,9) 9.06E-01

    # IOV over study periods 1-6; occasions 2-6 share the occasion-1
    # variance per $OMEGA BLOCK(1) SAME.
    etaiov_fcentral_1 ~ 0.603 # Table 6 'IOV pre washout nicotine dose' 91% CV; ESM output OMEGA(10,10) 6.03E-01
    etaiov_fcentral_2 ~ fixed(0.603) # $OMEGA BLOCK(1) SAME
    etaiov_fcentral_3 ~ fixed(0.603) # $OMEGA BLOCK(1) SAME
    etaiov_fcentral_4 ~ fixed(0.603) # $OMEGA BLOCK(1) SAME
    etaiov_fcentral_5 ~ fixed(0.603) # $OMEGA BLOCK(1) SAME
    etaiov_fcentral_6 ~ fixed(0.603) # $OMEGA BLOCK(1) SAME
    etaiov_fdepot_1 ~ 0.0164 # Table 6 'IOV F' 12.9% CV; ESM output OMEGA(16,16) 1.64E-02
    etaiov_fdepot_2 ~ fixed(0.0164) # $OMEGA BLOCK(1) SAME
    etaiov_fdepot_3 ~ fixed(0.0164) # $OMEGA BLOCK(1) SAME
    etaiov_fdepot_4 ~ fixed(0.0164) # $OMEGA BLOCK(1) SAME
    etaiov_fdepot_5 ~ fixed(0.0164) # $OMEGA BLOCK(1) SAME
    etaiov_fdepot_6 ~ fixed(0.0164) # $OMEGA BLOCK(1) SAME

    propSd <- 0.19
    label("Proportional residual error (fraction)") # Table 6 'Proportional residual error' 0.19
    addSd <- 0.257
    label("Additive residual error (ng/mL)") # Table 6 'Additive residual error (ng/mL)' 0.257
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
    iov_fdepot <- oc1 * etaiov_fdepot_1 + oc2 * etaiov_fdepot_2 +
      oc3 * etaiov_fdepot_3 + oc4 * etaiov_fdepot_4 +
      oc5 * etaiov_fdepot_5 + oc6 * etaiov_fdepot_6

    # 2. Allometric scaling (ESM $PK).
    wt_cl <- (WT / 70)^e_wt_cl_q
    wt_v <- (WT / 70)^e_wt_vc_vp

    # 3. Product-specific release parameters (ESM TVFRD, TVFR1, ALAG1).
    tlag1 <- exp(ltlag_invisipatch) * FORM_NICOTINE_INVISIPATCH
    tlag2 <- exp(ltlag2)
    logit_frdur <- (1 - FORM_NICOTINE_INVISIPATCH) * logitfrdur_nicorette +
      FORM_NICOTINE_INVISIPATCH * logitfrdur_invisipatch + etalogitfrdur
    frdur <- expit(logit_frdur)
    dur_ffo <- frdur * 16
    logit_ffo <- (1 - FORM_NICOTINE_INVISIPATCH) * logitffo_nicorette +
      FORM_NICOTINE_INVISIPATCH * logitffo_invisipatch + etalogitffo
    ffo <- expit(logit_ffo)

    # 4. Time since application. rxode2's tad()/tafd() count from the
    #    lagged arrival of the depot_td dose, so tlag1 is added back to give
    #    the time since the patch was applied.
    tad_patch <- tad(depot_td)
    if (is.na(tad_patch)) tad_patch <- 1e6
    tad_patch <- tad_patch + tlag1
    tafd_patch <- tafd(depot_td)
    if (is.na(tafd_patch)) tafd_patch <- -1e6
    tafd_patch <- tafd_patch + tlag1

    # First-order release runs from application to application + Frdur1*16 h
    # (ESM MTIME windows and FL2); clearance is higher once more than 25 h
    # have passed since the first application of the period (ESM TFLAG).
    release_on <- 0
    if (tad_patch <= dur_ffo) release_on <- 1
    after25 <- 0
    if (tafd_patch > 25) after25 <- 1

    cl_time_max <- exp(lcl_time_max)

    # 5. Individual parameters. cl is the total, time-varying clearance.
    cl <- exp(lcl + etalcl) * wt_cl * (1 + cl_time_max * after25)
    vc <- exp(lvc + etalvc) * wt_v
    q <- exp(lq) * wt_cl
    vp <- exp(lvp + etalvp) * wt_v
    q2 <- exp(lq2) * wt_cl
    vp2 <- exp(lvp2 + etalvp2) * wt_v
    krel <- exp(lkrel + etalkrel)
    ktr <- exp(lktr + etalktr)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # 6. ODEs (ESM $DES: ABS1 -> TRANS1 at KA while the first-order window
    #    is open; TRANS1 -> TRANS2 -> TRANS3 -> CENTRAL at KTR).
    d/dt(depot_td) <- -krel * release_on * depot_td
    d/dt(central) <- ktr * transit3 -
      kel * central - k12 * central + k21 * peripheral1 -
      k13 * central + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2
    d/dt(transit1) <- krel * release_on * depot_td - ktr * transit1
    d/dt(transit2) <- ktr * transit1 - ktr * transit2
    d/dt(transit3) <- ktr * transit2 - ktr * transit3

    # 7. Dosing. Each patch is TWO records of amt = the nominal patch
    #    content: a bolus into depot_td (first-order pathway) and a
    #    zero-order record into transit1 with rate = -2 (modelled
    #    duration). ESM F1 = FTOT*FR1/(1-EXP(-KA*DD)) -- scaled so the whole
    #    first-order fraction is released by the end of its window -- and
    #    F5 = FTOT*(1-FR1), ALAG5 = 4.06 h, D5 = DUR - ALAG5. Only the
    #    virtual pre-washout smoking bolus enters central.
    fdepot <- exp(lfdepot + iov_fdepot)
    f(depot_td) <- fdepot * ffo / (1 - exp(-krel * dur_ffo))
    alag(depot_td) <- tlag1
    f(transit1) <- fdepot * (1 - ffo)
    alag(transit1) <- tlag2
    dur(transit1) <- T_PATCH_WEAR - tlag2
    fcentral <- exp(lfcentral + etalfcentral + iov_fcentral)
    f(central) <- fcentral

    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
