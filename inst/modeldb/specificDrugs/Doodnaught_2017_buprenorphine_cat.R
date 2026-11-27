Doodnaught_2017_buprenorphine_cat <- function() {
  description <- "Preclinical/clinical veterinary (cat). Joint population PK-PD model for a high-concentration buprenorphine formulation (Simbadol, 1.8 mg/mL) and its metabolite norbuprenorphine in conscious cats after intravenous, buccal (oral transmucosal) and subcutaneous dosing, parameterised per kg body weight. Buprenorphine and norbuprenorphine are each two-compartment; norbuprenorphine is formed from central buprenorphine by a first-order conversion and, after buccal dosing, also enters through a first-pass gut-norbuprenorphine depot. Subcutaneous absorption is biphasic: an inverse-Gaussian input (fraction BIO of the bioavailable dose, peak at 5 min) plus a time-dependent first-order input whose rate constant rises logistically to its maximum. Thermal antinociception (thermal threshold, degrees C) is a sigmoid Emax function of an effect-compartment buprenorphine concentration (Doodnaught 2017)"
  reference <- paste(
    "Doodnaught GM, Monteiro BP, Benito J, Edge D, Beaudry F, Pelligand L,",
    "Steagall P. Pharmacokinetic and pharmacodynamic modelling after",
    "subcutaneous, intravenous and buccal administration of a high-concentration",
    "formulation of buprenorphine in conscious cats. PLoS One.",
    "2017;12(4):e0176443. doi:10.1371/journal.pone.0176443.",
    "Model structure from Figure 2 and S3 File (Appendix 3, equations 1-9).",
    sep = " "
  )
  vignette <- "Doodnaught_2017_buprenorphine_cat"
  # Every disposition parameter is published per kg body weight (CL in
  # L/kg/h, volumes in L/kg), so doses are given in mg/kg and compartment
  # amounts are carried in mg/kg. central / vc is then mg/L; the x1000 reports
  # ng/mL, the units of the assay (Methods 'Buprenorphine and norbuprenorphine
  # assay') and of the EC50.
  units <- list(time = "h", dosing = "mg/kg", concentration = "ng/mL")

  dosing <- c("central", "depot_buccal", "depot_norbup", "depot_sc")
  compartmentData <- list(
    depot_buccal = list(analyte = "buprenorphine", units = "mg/kg", specimen = "administration site", verified = TRUE),
    depot_sc = list(analyte = "buprenorphine", units = "mg/kg", specimen = "administration site", verified = TRUE),
    depot_norbup = list(analyte = "norbuprenorphine", units = "mg/kg", specimen = "tissue", verified = TRUE),
    central = list(analyte = "buprenorphine", units = "mg/kg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "buprenorphine", units = "mg/kg", specimen = "plasma", verified = TRUE),
    central_norbup = list(analyte = "norbuprenorphine", units = "mg/kg", specimen = "plasma", verified = TRUE),
    peripheral1_norbup = list(analyte = "norbuprenorphine", units = "mg/kg", specimen = "plasma", verified = TRUE),
    effect = list(analyte = "buprenorphine", units = "ng/mL", specimen = "not applicable", verified = TRUE)
  )

  # The thermal threshold (degrees C, Topcat WTT1 wireless device) is a novel
  # PD observation in the library; kept paper-specific until a second
  # thermal-nociception model arrives.
  paper_specific_compartments <- c("thermal_threshold")

  covariateData <- list()

  population <- list(
    species = "cat (domestic shorthair; Felis catus)",
    n_subjects = 6L,
    n_studies = 1L,
    n_administrations = 18L,
    weight_range = "4.9 +/- 0.7 kg (mean +/- SD)",
    sex_female_pct = 33.3,
    disease_state = "healthy adult laboratory cats (conscious)",
    dose_range = paste(
      "single dose of Simbadol (buprenorphine HCl 1.8 mg/mL): 0.12 mg/kg i.v.,",
      "0.24 mg/kg s.c. or 0.12 mg/kg buccal (oral transmucosal), randomised",
      "crossover with a 14-day minimum washout",
      sep = " "
    ),
    regions = "Canada (Universite de Montreal)",
    notes = paste(
      "Methods 'Animals': six domestic short-haired cats, four males and two females.",
      "Phase II crossover: every cat received all three routes, giving 18",
      "concentration-time profiles of each analyte. Plasma sampled from a jugular",
      "central venous catheter at 2 and 8 min and at 0.5, 1, 2, 4, 6, 8, 12, 24, 30, 36,",
      "48, 60 and 72 h; the buccal dose was given in the cheek contralateral to the",
      "catheter. Thermal threshold measured at baseline and at the same times from",
      "0.5 h, with a 55 degrees C safety cut-off. A saline (Phase I) arm contributed",
      "thermal thresholds only and is not part of the PK-PD fit.",
      sep = " "
    )
  )

  ini({
    # ---- Buprenorphine disposition (Results 'Pharmacokinetic modeling') ----
    # Shared by the three routes; IIV % follows each estimate in the text.
    lcl <- log(0.98); label("Buprenorphine clearance CL (L/kg/h)") # Results 'CL = 0.98 L kg-1 hour-1, 2.4%'
    lvc <- log(0.75); label("Buprenorphine central volume V1 (L/kg)") # Results 'V1 = 0.75 L kg-1, 11.3%'
    lq <- log(0.70); label("Buprenorphine intercompartmental clearance CL2 (L/kg/h)") # Results 'CL2 = 0.70 L kg-1 hour-1, 0%'
    lvp <- log(7.15); label("Buprenorphine peripheral volume V2 (L/kg)") # Results 'V2 = 7.15 L kg-1, 8.1%'

    # ---- Buccal (OTM) route ----
    # Bioavailabilities were estimated on the logit scale (S3 File eq 2).
    logitfdepot_buccal <- log(0.236 / (1 - 0.236)); label("Logit of buccal bioavailability FOTM, applied to both the buprenorphine and the first-pass norbuprenorphine input (unitless)") # Results 'FOTM was 23.6% (IIV 25%)'
    lka_buccal <- log(1.67); label("Buccal absorption rate constant kaOTM (1/h)") # Results 'kaOTM was 1.67 hour-1 (25.1%)'

    # ---- Subcutaneous route (S3 File eqs 5-7) ----
    logitfdepot_sc <- log(0.94 / (1 - 0.94)); label("Logit of total subcutaneous bioavailability FSC (unitless)") # Results 'FSC was 94% (IIV 23.4%)'
    logitfrel <- log(0.101 / (1 - 0.101)); label("Logit of the fraction BIO of the SC bioavailable dose absorbed by the inverse-Gaussian input (unitless)") # Results 'A proportion of 10.1% (IIV 113%) of FSC ... (BIO ...)'
    lmat_invgauss <- log(7.21); label("Mean input time MAT of the inverse-Gaussian SC input (h)") # Results 'MAT was 7.21 hours (3.5%)'
    lcv_invgauss <- log(5.46); label("Relative dispersion CV of the inverse-Gaussian SC input (unitless; squared in S3 File eq 5)") # Results 'variance of the input time (CV) was 5.46 (0.5%)'
    lka_sc <- log(0.062); label("Maximal time-dependent SC absorption rate constant kaSC (1/h)") # Results 'kaSC was 0.062 hour-1 (14.9%)'
    lka_t50 <- log(2.8); label("Time to 50% of the maximal SC absorption rate constant T50 (h)") # Results 'T50 was 2.8 hour (27%)'

    # ---- Norbuprenorphine (Results, norbuprenorphine paragraph) ----
    lkmet <- fixed(log(0.196)); label("Conversion rate constant of central buprenorphine to norbuprenorphine KaMet (1/h)") # Results 'KaMet = 0.196 hour-1, 19%'; estimated on the IV + SC data, 'then fixed' for the final 3-route fit
    lka_norbup <- log(0.626); label("First-pass norbuprenorphine absorption rate constant Kafirst-pass (1/h)") # Results 'Kafirst-pass 0.626 hour-1, 6.9%'
    lcl_norbup <- log(0.42); label("Norbuprenorphine clearance CLMet (L/kg/h)") # Results 'CLMet = 0.42 L kg-1 hour-1, IIV 4.3%'
    lvc_norbup <- log(0.323); label("Norbuprenorphine central volume V1Met (L/kg)") # Results 'V1Met = 0.323 L kg-1, 0.1%'
    lq_norbup <- log(7.662); label("Norbuprenorphine intercompartmental clearance CL2Met (L/kg/h)") # Results 'CL2 Met = 7.662 L kg-1 hour-1'; no random effect
    lvp_norbup <- log(4.69); label("Norbuprenorphine peripheral volume V2Met (L/kg)") # Results 'V2 Met = 4.69 L kg-1'; no random effect

    # ---- Thermal antinociception (Results 'PK-PD modeling', eq 1) ----
    lke0 <- log(0.52); label("Effect-compartment equilibration rate constant Ke0 (1/h)") # Results 'Ke0 was 0.52 hour-1 (IIV 0%)'
    le0 <- log(46.1); label("Baseline thermal threshold T0 (degrees C)") # Results 'T0 = 46.1 (IIV 0%)'
    lemax <- log(5.62); label("Maximal increase in thermal threshold Emax (degrees C)") # Results 'Emax was 5.62 degrees C (IIV 0%)'
    lec50 <- log(2.13); label("Effect-site buprenorphine concentration at half-maximal effect EC50 (ng/mL)") # Results 'EC50 ... was 2.13 ng mL-1 (IIV 447%)'
    lhill <- log(1.54); label("Hill coefficient n of the concentration-effect curve (unitless)") # Results 'slope (n) of 1.54 (91%)'

    # ---- IIV: CV% converted with S3 File eq 3, omega^2 = log(CV^2 + 1) ----
    # The logit-scale bioavailability etas use the same conversion. Parameters
    # reported with 'IIV 0%' (CL2, Ke0, T0, Emax) and those whose random
    # effects 'were not estimated' (CL2Met, V2Met) carry no eta.
    etalcl ~ 0.000575834 # Results 'CL ... 2.4%'
    etalvc ~ 0.0126882 # Results 'V1 ... 11.3%'
    etalvp ~ 0.00653957 # Results 'V2 ... 8.1%'
    etalogitfdepot_buccal ~ 0.0606246 # Results 'FOTM ... (IIV 25%)'
    etalka_buccal ~ 0.061096 # Results 'kaOTM ... (25.1%)'
    etalogitfdepot_sc ~ 0.0533095 # Results 'FSC ... (IIV 23.4%)'
    etalogitfrel ~ 0.822815 # Results 'BIO ... (IIV 113%)'
    etalmat_invgauss ~ 0.00122425 # Results 'MAT ... (3.5%)'
    etalcv_invgauss ~ 2.49997e-05 # Results 'CV ... (0.5%)'
    etalka_sc ~ 0.0219581 # Results 'kaSC ... (14.9%)'
    etalka_t50 ~ 0.0703653 # Results 'T50 ... (27%)'
    etalkmet ~ 0.0354637 # Results 'KaMet ... 19%'
    etalka_norbup ~ 0.0047497 # Results 'Kafirst-pass ... 6.9%'
    etalcl_norbup ~ 0.00184729 # Results 'CLMet ... IIV 4.3%'
    etalvc_norbup ~ 9.99999e-07 # Results 'V1Met ... 0.1%'
    etalec50 ~ 3.04361 # Results 'EC50 ... (IIV 447%)'
    etalhill ~ 0.603277 # Results 'slope (n) ... (91%)'

    # ---- Residual error ----
    # Figure 2 and Results: a proportional error for each analyte, common to
    # the three routes; S3 File Table 2: additive error for the PD model. No
    # magnitude is reported for any of them, so they are fixed to zero.
    propSd <- fixed(0); label("Proportional residual error, buprenorphine (fraction; not reported)") # Results 'common proportional residual error term'; value not reported
    propSd_norbup <- fixed(0); label("Proportional residual error, norbuprenorphine (fraction; not reported)") # Figure 2 'Residual: epsilon_Nor'; value not reported
    addSd_thermal_threshold <- fixed(0); label("Additive residual error, thermal threshold (degrees C; not reported)") # S3 File Table 2 'additive error model'; value not reported
  })

  model({
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q <- exp(lq)
    vp <- exp(lvp + etalvp)

    fdepot_buccal <- expit(logitfdepot_buccal + etalogitfdepot_buccal)
    ka_buccal <- exp(lka_buccal + etalka_buccal)

    fdepot_sc <- expit(logitfdepot_sc + etalogitfdepot_sc)
    frel <- expit(logitfrel + etalogitfrel)
    mat_invgauss <- exp(lmat_invgauss + etalmat_invgauss)
    cv_invgauss <- exp(lcv_invgauss + etalcv_invgauss)
    ka_sc <- exp(lka_sc + etalka_sc)
    ka_t50 <- exp(lka_t50 + etalka_t50)

    kmet <- exp(lkmet + etalkmet)
    ka_norbup <- exp(lka_norbup + etalka_norbup)
    cl_norbup <- exp(lcl_norbup + etalcl_norbup)
    vc_norbup <- exp(lvc_norbup + etalvc_norbup)
    q_norbup <- exp(lq_norbup)
    vp_norbup <- exp(lvp_norbup)

    ke0 <- exp(lke0)
    e0 <- exp(le0)
    emax <- exp(lemax)
    ec50 <- exp(lec50 + etalec50)
    hill <- exp(lhill + etalhill)

    # Time since, and amount of, the last SC dose. Both are NA before any SC
    # dose (IV and buccal arms), so they are zeroed there; the SC terms below
    # are then zero because depot_sc is empty and the dose amount is zero.
    tad_sc <- tad(depot_sc)
    if (is.na(tad_sc)) {
      tad_sc <- 0
    }
    dose_sc <- podo(depot_sc)
    if (is.na(dose_sc)) {
      dose_sc <- 0
    }

    # Inverse-Gaussian SC input rate (S3 File eq 5), carrying the fraction
    # FSC x BIO of the dose straight into central. The density is zero at the
    # dose time. With CV = 5.46 its mode is 0.08 h (S3 File eq 6; Results
    # 'initial peak ... at 0.08 hours'). Single-dose form: it reads only the
    # most recent SC dose.
    input_invgauss <- 0
    if (tad_sc > 0) {
      input_invgauss <- dose_sc * fdepot_sc * frel *
        sqrt(mat_invgauss / (2 * pi * cv_invgauss^2 * tad_sc^3)) *
        exp(-(tad_sc - mat_invgauss)^2 / (2 * cv_invgauss^2 * mat_invgauss * tad_sc))
    }

    # Time-dependent SC input (S3 File eq 7) for the remaining FSC x (1 - BIO):
    # a first-order rate constant rising logistically to kaSC with half-time
    # T50. The shape parameter gamma 'was approximated by the value of
    # clearance' (S3 File, after eq 7), i.e. the per-kg CL taken in hours.
    ka_sc_t <- ka_sc / (1 + exp(-(tad_sc - ka_t50) / cl))

    d/dt(depot_buccal) <- -ka_buccal * depot_buccal
    d/dt(depot_sc) <- -ka_sc_t * depot_sc
    d/dt(depot_norbup) <- -ka_norbup * depot_norbup
    # Buprenorphine leaves central only by CL; the conversion KaMet x central
    # to norbuprenorphine is a component of that clearance (Results 'total
    # body clearance'; the 12.3 h beta half-life follows from CL, V1, CL2, V2
    # alone).
    d/dt(central) <- ka_buccal * depot_buccal + input_invgauss + ka_sc_t * depot_sc -
      cl / vc * central - q / vc * central + q / vp * peripheral1
    d/dt(peripheral1) <- q / vc * central - q / vp * peripheral1
    d/dt(central_norbup) <- kmet * central + ka_norbup * depot_norbup -
      cl_norbup / vc_norbup * central_norbup -
      q_norbup / vc_norbup * central_norbup + q_norbup / vp_norbup * peripheral1_norbup
    d/dt(peripheral1_norbup) <- q_norbup / vc_norbup * central_norbup - q_norbup / vp_norbup * peripheral1_norbup

    # Buccal dose: FOTM applies to the buprenorphine depot and to the gut
    # norbuprenorphine (first-pass) depot alike (Results 'FOTM (parent and
    # metabolite OTM bioavailability)'), so a buccal dose is entered as two
    # records of the same amount, one into each depot.
    f(depot_buccal) <- fdepot_buccal
    f(depot_norbup) <- fdepot_buccal
    f(depot_sc) <- fdepot_sc * (1 - frel)

    Cc <- central / vc * 1000
    Cc_norbup <- central_norbup / vc_norbup * 1000

    # Effect compartment (S3 File eq 8) and sigmoid Emax (eq 1 / S3 File eq 9)
    d/dt(effect) <- ke0 * (Cc - effect)
    thermal_threshold <- e0 + emax * effect^hill / (ec50^hill + effect^hill)

    Cc ~ prop(propSd)
    Cc_norbup ~ prop(propSd_norbup)
    thermal_threshold ~ add(addSd_thermal_threshold)
  })
}
