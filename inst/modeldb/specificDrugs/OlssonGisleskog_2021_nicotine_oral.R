OlssonGisleskog_2021_nicotine_oral <- function() {
  description <- paste(
    "Population PK model for orally ingested nicotine microtablets in",
    "healthy adult smokers (Olsson Gisleskog 2021; 26 subjects, 2 studies,",
    "2 and 6 mg). Disposition (three compartments with allometric weight",
    "scaling) and its IIV are fixed to the paper's intravenous model.",
    "Absorption differs by study: in 92NNBT005 (repeated dosing, tablets",
    "chewed) first-order absorption directly into the central compartment",
    "(ka 1.55 1/h, F 39.5%); in 93NNBT007 (single dose, tablets swallowed",
    "whole) first-order transfer at the same ka into a chain of three",
    "transit compartments (ktr 3.60 1/h, F 22.3%). A transient",
    "time-dependent increase in clearance during the repeated-dose day",
    "(77.3% maximum, 50% onset at 5.29 h after the first dose, lasting",
    "1.92 h) applies in 92NNBT005 only. Residual pre-study nicotine is a",
    "virtual 1-mg bolus into central at the start of the washout with",
    "bioavailability 4.92 mg.",
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
    STUDY_93NNBT007 = list(
      description = "Study 93NNBT007 indicator (single-dose microtablets swallowed whole)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (study 92NNBT005, repeated doses of microtablets chewed before swallowing)",
      notes = paste(
        "Selects the absorption structure and bioavailability: 1 = transit",
        "chain with F 22.3% and no time-dependent clearance change, 0 =",
        "direct first-order absorption with F 39.5% and the transient",
        "clearance increase. The paper attributes the difference to the",
        "microtablets being swallowed whole in 93NNBT007 but chewed in",
        "92NNBT005 (Sect. 2.1 and Discussion)."
      ),
      source_name = "STUDY (STUDY.EQ.6 for 93NNBT007, STUDY.EQ.5 for 92NNBT005)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 26L,
    n_studies = 2L,
    age_range = "24-47 years",
    age_median = "38 years",
    weight_range = "44.0-105.0 kg",
    weight_median = "62.5 kg",
    sex_female_pct = 53.8,
    race_ethnicity = c(Missing = 100),
    disease_state = "healthy adult smokers",
    dose_range = "2 and 6 mg nicotine microtablets, single (93NNBT007) and repeated (92NNBT005) oral doses",
    regions = "Sweden",
    notes = "Olsson Gisleskog 2021 Table 1 (oral microtablets row) and Table 2 (Oral column); race was not recorded."
  )

  ini({
    # Disposition fixed to the IV model (ESM oral control stream
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

    # Absorption (Table 4).
    lka <- log(1.55)
    label("First-order absorption rate constant ka (1/h)") # Table 4 'Ka (h-1)' 1.55
    lfdepot_92nnbt005 <- log(0.395)
    label("Oral bioavailability in study 92NNBT005 (fraction)") # Table 4 'F study 92NNBT005 (%)' 39.5
    lfdepot_93nnbt007 <- log(0.223)
    label("Oral bioavailability in study 93NNBT007 (fraction)") # Table 4 'F study 93NNBT007 (%)' 22.3
    lktr <- log(3.60)
    label("Transit rate constant ktr, study 93NNBT007 (1/h)") # Table 4 'Ktr study 93NNBT007 (h-1)' 3.60

    # Time-dependent clearance change, Eq. 2 (Table 4), study 92NNBT005.
    lcl_t50 <- log(5.29)
    label("Time after the first dose of 50% onset of the clearance increase, Start (h)") # Table 4 'Start (h)' 5.29
    lcl_time_dur <- log(1.92)
    label("Duration of the clearance increase, Duration (h)") # Table 4 'Duration (h)' 1.92
    lcl_time_max <- log(0.773)
    label("Maximal fractional clearance increase, Emax (fraction)") # Table 4 'Emax (%)' 77.3
    lcl_time_hill <- log(12.3)
    label("Steepness of onset and offset of the clearance increase, pow (unitless)") # Table 4 'pow' 12.3

    lfcentral <- log(4.92)
    label("Pre-washout nicotine dose: bioavailability of the virtual 1-mg central bolus (mg)") # Table 4 'Pre-washout nicotine dose (mg)' 4.92

    # IIV fixed from the IV model (ESM $OMEGA ... FIX).
    etalcl ~ fixed(0.0705245) # ESM $OMEGA 0.0705245 FIX (Table 3 IIV CL 27.0%)
    etalvc + etalvp2 + etalvp ~ fixed(c(
      0.381077,
      -0.230826, 0.450311,
      0, 0.5527, 1.83554
    )) # ESM $OMEGA BLOCK(3) FIX (Table 3)
    # Estimated IIV (Table 4 CV% = 100*sqrt(exp(omega^2)-1); variances from
    # the ESM output).
    etalfdepot ~ 0.0499 # Table 4 'IIV F' 22.6% CV; ESM output OMEGA(5,5) 4.99E-02
    etalktr ~ 0.186 # Table 4 'IIV Ktr' 45.2% CV; ESM output OMEGA(6,6) 1.86E-01
    etalfcentral ~ 0.610 # Table 4 'IIV pre-washout nicotine dose' 91.7% CV; ESM output OMEGA(7,7) 6.10E-01
    etalcl_time_dur ~ 0.340 # Table 4 'IIV duration' 63.7% CV; ESM output OMEGA(8,8) 3.40E-01

    propSd <- 0.0987
    label("Proportional residual error (fraction)") # Table 4 'Proportional residual error (%)' 9.87
    addSd <- 0.162
    label("Additive residual error (ng/mL)") # Table 4 'Additive residual error (ng/mL)' 0.162
  })
  model({
    # 1. Allometric scaling (ESM $PK).
    wt_cl <- (WT / 70)^e_wt_cl_q
    wt_v <- (WT / 70)^e_wt_vc_vp

    # 2. Time-dependent clearance change, Eq. 2:
    #    CL = CLbl * (1 + Emax * [T^pow/(Start^pow + T^pow) -
    #                             T^pow/((Start + Duration)^pow + T^pow)]),
    #    T = time since the first dose of the period (ESM TSSP), taken here
    #    as the time after the first oral dose. Written as
    #    1/(1 + (Start/T)^pow) for numerical stability; identical for T > 0.
    #    Before the first dose (T undefined or 0) T is floored at 1e-6 h,
    #    where both terms are 0 to machine precision (ESM TFL = 0 unless
    #    TSSP > 0). The change applies in 92NNBT005 only (ESM
    #    IF (STUDY.EQ.6) CL_EFF = 1).
    tafd_oral <- tafd(depot)
    if (is.na(tafd_oral)) tafd_oral <- 0
    if (tafd_oral < 1e-6) tafd_oral <- 1e-6
    cl_time_max <- exp(lcl_time_max)
    cl_time_hill <- exp(lcl_time_hill)
    cl_t50 <- exp(lcl_t50)
    cl_time_dur <- exp(lcl_time_dur + etalcl_time_dur)
    cl_t_end <- cl_t50 + cl_time_dur
    cl_time_frac <- 1 / (1 + (cl_t50 / tafd_oral)^cl_time_hill) -
      1 / (1 + (cl_t_end / tafd_oral)^cl_time_hill)
    cl_time_fac <- 1 + (1 - STUDY_93NNBT007) * cl_time_max * cl_time_frac

    # 3. Individual parameters. cl is the total, time-varying clearance.
    cl <- exp(lcl + etalcl) * wt_cl * cl_time_fac
    vc <- exp(lvc + etalvc) * wt_v
    q <- exp(lq) * wt_cl
    vp <- exp(lvp + etalvp) * wt_v
    q2 <- exp(lq2) * wt_cl
    vp2 <- exp(lvp2 + etalvp2) * wt_v
    ka <- exp(lka)
    ktr <- exp(lktr + etalktr)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # 4. Absorption (ESM $DES): 92NNBT005 absorbs depot -> central at ka;
    #    93NNBT007 transfers depot -> transit1 at ka, then three transits
    #    at ktr into central.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- (1 - STUDY_93NNBT007) * ka * depot +
      STUDY_93NNBT007 * ktr * transit3 -
      kel * central - k12 * central + k21 * peripheral1 -
      k13 * central + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2
    d/dt(transit1) <- STUDY_93NNBT007 * ka * depot - ktr * transit1
    d/dt(transit2) <- ktr * transit1 - ktr * transit2
    d/dt(transit3) <- ktr * transit2 - ktr * transit3

    # 5. Bioavailability. Only the virtual pre-washout smoking bolus is
    #    dosed into central, so f(central) is its bioavailability (the
    #    pre-washout dose in mg per 1-mg bolus).
    fdepot <- exp((1 - STUDY_93NNBT007) * lfdepot_92nnbt005 +
      STUDY_93NNBT007 * lfdepot_93nnbt007 + etalfdepot)
    f(depot) <- fdepot
    f(central) <- exp(lfcentral + etalfcentral)

    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
