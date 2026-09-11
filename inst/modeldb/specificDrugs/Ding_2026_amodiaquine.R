# Joint parent-metabolite population PK model for oral amodiaquine and its
# active metabolite desethylamodiaquine in patients with uncomplicated
# Plasmodium falciparum malaria treated with artemether-lumefantrine plus
# amodiaquine, pooled from the TRACII (NCT02453308) and TACT-CV
# (NCT03355664) trials (Ding 2026, Br J Clin Pharmacol 92(2):589-605;
# doi:10.1002/bcp.70301).

Ding_2026_amodiaquine <- function() {
  description <- paste(
    "Joint parent-metabolite population PK model for oral amodiaquine and",
    "its active CYP2C8-derived metabolite desethylamodiaquine in adults and",
    "children with acute uncomplicated Plasmodium falciparum malaria",
    "receiving artemether-lumefantrine plus amodiaquine as a triple",
    "artemisinin-based combination therapy (Ding 2026, pooled TRACII +",
    "TACT-CV, n = 302). First-order absorption feeds a two-compartment",
    "amodiaquine disposition model, with complete molar-corrected",
    "bioconversion to a three-compartment desethylamodiaquine disposition",
    "model. Allometric body-weight scaling on all apparent clearances",
    "(fixed exponent 0.75) and apparent volumes (fixed exponent 1.0) at a",
    "reference weight of 45 kg. Relative bioavailability is anchored at 1",
    "with inter-occasion variability on both bioavailability and the",
    "absorption rate constant. No covariate other than body weight was",
    "retained. Predictions are plasma amodiaquine and desethylamodiaquine",
    "concentrations in ng/mL.",
    sep = " "
  )
  reference <- paste(
    "Ding J, Hoglund RM, van der Pluijm RW, Callery JJ, Peto TJ, Tripura R,",
    "Das S, Nguyen HC, Promnarate C, Mukaka M, Dysoley L, Fanello C,",
    "Onyamboko MA, Anvikar AR, Mayxay M, Smithuis F, von Seidlein L,",
    "Dhorda M, Amaratunga C, Faiz MA, Ho DTN, White NJ, Day NPJ,",
    "Dondorp AM, Tarning J (2026). Population pharmacokinetics of",
    "artemether-lumefantrine plus amodiaquine in patients with",
    "uncomplicated Plasmodium falciparum malaria. British Journal of",
    "Clinical Pharmacology 92(2):589-605. doi:10.1002/bcp.70301.",
    sep = " "
  )
  vignette <- "Ding_2026_artemether_lumefantrine_amodiaquine"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Ding 2026 Methods ('PK sampling
  # scheme', 'Drug quantification': venous plasma amodiaquine and
  # desethylamodiaquine by LC-MS/MS) and Figure S5.
  compartmentData <- list(
    depot            = list(analyte = "amodiaquine",         units = "mg", specimen = "administration site", verified = TRUE),
    central          = list(analyte = "amodiaquine",         units = "mg", specimen = "plasma",              verified = TRUE),
    peripheral1      = list(analyte = "amodiaquine",         units = "mg", specimen = "plasma",              verified = TRUE),
    central_deaq     = list(analyte = "desethylamodiaquine", units = "mg", specimen = "plasma",              verified = TRUE),
    peripheral1_deaq = list(analyte = "desethylamodiaquine", units = "mg", specimen = "plasma",              verified = TRUE),
    peripheral2_deaq = list(analyte = "desethylamodiaquine", units = "mg", specimen = "plasma",              verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed at baseline. Ding 2026 Methods ('Covariates model'):",
        "'bodyweight was included on all clearance and volume parameters",
        "using a conventional allometric function with fixed exponents of",
        "0.75 and 1.0, respectively'. The Table 3 footnote fixes the",
        "reference: 'Population estimates are given for a typical adult",
        "patient weighing 45 kg with acute P. falciparum malaria', so 45 kg",
        "is the normalising constant encoded here. Unlike for artemether,",
        "the allometric term was strongly supported for amodiaquine",
        "(Results 3.1.2: delta-OFV = -189.216). Cohort median 41.5 kg",
        "(TRACII) and 52.2 kg (TACT-CV); range 9.0-98.8 kg (Table 1).",
        sep = " "
      ),
      source_name        = "BW"
    ),
    OCC = list(
      description        = "Dose occasion index, 1 to 6 across the six-dose amodiaquine regimen",
      units              = "(count)",
      type               = "categorical",
      reference_category = "1 (the first dose, at hour 0)",
      notes              = paste(
        "Integer occasion column taking value k on the interval starting at",
        "the k-th dose. Amodiaquine is given together with",
        "artemether-lumefantrine at 0, 8, 24, 36, 48 and 60 h (Methods,",
        "'Dosing regimen'; Table S2), so OCC = 1 for 0 <= t < 8 h, 2 for",
        "8 <= t < 24 h, 3 for 24 <= t < 36 h, 4 for 36 <= t < 48 h, 5 for",
        "48 <= t < 60 h and 6 for t >= 60 h. Used purely as the occasion",
        "grouping for the inter-occasion variability that Ding 2026",
        "Results 3.1.2 added on relative bioavailability",
        "(delta-OFV = -36.276) and on the absorption rate constant",
        "(delta-OFV = -26.552). Note that patients in the 5-14.9 kg weight",
        "band receive amodiaquine only at 0, 24 and 48 h (Table S2); the",
        "occasion index still advances at every artemether-lumefantrine",
        "dose time, and the occasions with no amodiaquine dose simply carry",
        "no absorption event.",
        sep = " "
      ),
      source_name        = "OCC"
    )
  )

  covariatesDataExcluded <- list(
    STUDY_TACTCV = list(
      description = "TACT-CV trial indicator (1 = TACT-CV, 0 = TRACII)",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "The only covariate the stepwise search identified for this model,",
        "and it was deliberately discarded: Ding 2026 Results 3.1.2, 'The",
        "only covariate identified in the covariate search was a difference",
        "in amodiaquine or desethylamodiaquine intercompartment clearance",
        "between trials (353% higher in the TACT-CV trial). This covariate",
        "was deemed implausible and was not retained in the final model.'",
        "The final model therefore has no study term. Retained here as",
        "documentation of the covariate screen; the same column IS retained",
        "in the sibling Ding_2026_lumefantrine.R.",
        sep = " "
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 302L,
    n_studies      = 2L,
    n_observations = paste(
      "352 amodiaquine concentrations from the 39 dense-PK patients (21",
      "TRACII, 18 TACT-CV) and 725 desethylamodiaquine concentrations from",
      "all 302 patients (150 TRACII, 152 TACT-CV) (Results 3.1.2)"
    ),
    age_range      = "17.0 years (TRACII) and 25.0 years (TACT-CV), medians (range 2.1-62.0 years) (Table 1)",
    weight_range   = "41.5 kg (TRACII) and 52.2 kg (TACT-CV), medians (range 9.0-98.8 kg) (Table 1)",
    sex_female_pct = 23.8,
    disease_state  = paste(
      "Acute uncomplicated Plasmodium falciparum malaria. Median admission",
      "asexual parasite count 52,500 parasites/uL (TRACII) and 21,500",
      "parasites/uL (TACT-CV); median admission body temperature 37.5 and",
      "37.7 degC respectively (Table 1)."
    ),
    dose_range     = paste(
      "Amodiaquine 150 mg base per tablet, target 10 mg base per kg per day",
      "given as a split dose twice daily together with",
      "artemether-lumefantrine at 0, 8, 24, 36, 48 and 60 h, directly",
      "observed. Tablets per dose by weight band: 0.5 at 0/24/48 h only",
      "(5-14.9 kg), 0.5 at all six times (15-24.9 kg), 1 at all six times",
      "(25-34.9 kg), 1.5 at all six times (>35 kg) (Table S2). Dense-PK",
      "cohort median amodiaquine dose 8.6-9.0 mg/kg/day (range 5.8-12.2)",
      "(Table S3)."
    ),
    regions        = paste(
      "TRACII (NCT02453308): Bangladesh, India, Myanmar, Democratic",
      "Republic of Congo and Lao PDR. TACT-CV (NCT03355664): western and",
      "eastern Cambodia and Vietnam. Dense PK sampling was feasible at one",
      "site per trial (Bangladesh and Vietnam)."
    ),
    notes          = paste(
      "Only patients randomised to the artemether-lumefantrine plus",
      "amodiaquine arm contribute. Amodiaquine parent data come from the",
      "dense-PK sub-cohort only; desethylamodiaquine data come from the",
      "full arm via sparse baseline / Day 7 / recurrence sampling.",
      "Concentrations below the lower limit of quantification were",
      "discarded (Beal M1): five amodiaquine samples within 96 h and 15",
      "desethylamodiaquine samples within 28 days. Amodiaquine samples",
      "collected on Day 7 and beyond were treated as missing because",
      "measurable Day 7 concentrations are implausible given the 5.2-13.7 h",
      "literature half-life and were judged spurious (Results 3.1.2)."
    )
  )

  ini({
    # ---- Absorption --------------------------------------------------
    # Ding 2026 Results 3.1.2: "The absorption of amodiaquine was
    # adequately described by a first-order absorption model, with no
    # further improvement in model fit with a transit compartment model."
    lka <- log(1.93)
    label("First-order absorption rate constant of amodiaquine from depot into central (1/h)")
    # Ding 2026 Table 3: Ka = 1.93 1/h (%RSE 20.2; SIR median 1.92,
    # 95% CI 1.4-2.91; eta shrinkage 68.4%)

    # ---- Amodiaquine disposition (two compartments) ------------------
    # Ding 2026 Table 3, "NONMEM estimates" column. Values are apparent
    # (relative to F = 1) and reported on the linear scale for a typical
    # adult patient weighing 45 kg; log() is applied here for the nlmixr2
    # internal log scale.
    lcl <- log(2250)
    label("Apparent amodiaquine elimination clearance CL/F at WT = 45 kg (L/h)")
    # Ding 2026 Table 3: CL/F AQ = 2250 L/h (%RSE 3.2; SIR median 2250,
    # 95% CI 2120-2390; eta shrinkage 74.0%). Discussion cross-checks this
    # against a previously reported pooled value of 2735 L/h.

    lvc <- log(12900)
    label("Apparent amodiaquine central volume of distribution Vc/F at WT = 45 kg (L)")
    # Ding 2026 Table 3: Vc/F AQ = 12,900 L (%RSE 8.6; SIR median 13,000,
    # 95% CI 10,700-15,100)

    lq <- log(3020)
    label("Apparent amodiaquine inter-compartmental clearance Q/F at WT = 45 kg (L/h)")
    # Ding 2026 Table 3: Q/F AQ = 3020 L/h (%RSE 7.6; SIR median 3010,
    # 95% CI 2640-3540)

    lvp <- log(27600)
    label("Apparent amodiaquine peripheral volume of distribution Vp/F at WT = 45 kg (L)")
    # Ding 2026 Table 3: Vp/F AQ = 27,600 L (%RSE 4.8; SIR median 27,700,
    # 95% CI 25,000-30,400)

    # ---- Desethylamodiaquine disposition (three compartments) --------
    lcl_deaq <- log(32.2)
    label("Apparent desethylamodiaquine elimination clearance CL/F at WT = 45 kg (L/h)")
    # Ding 2026 Table 3: CL/F DEAQ = 32.2 L/h (%RSE 3.5; SIR median 32.3,
    # 95% CI 30.1-34.4; eta shrinkage 37.8%). Discussion cross-checks this
    # against a previously reported pooled value of 30.1 L/h.

    lvc_deaq <- log(1260)
    label("Apparent desethylamodiaquine central volume of distribution Vc/F at WT = 45 kg (L)")
    # Ding 2026 Table 3: Vc/F DEAQ = 1260 L (%RSE 8.6; SIR median 1260,
    # 95% CI 1090-1500; eta shrinkage 68.8%)

    lq_deaq <- log(117)
    label("Apparent desethylamodiaquine inter-compartmental clearance Q1/F to the shallow peripheral compartment at WT = 45 kg (L/h)")
    # Ding 2026 Table 3: Q1/F DEAQ = 117 L/h (%RSE 17.5; SIR median 118,
    # 95% CI 82.7-161)

    lvp_deaq <- log(1640)
    label("Apparent desethylamodiaquine shallow peripheral volume of distribution Vp1/F at WT = 45 kg (L)")
    # Ding 2026 Table 3: Vp1/F DEAQ = 1640 L (%RSE 15.7; SIR median 1630,
    # 95% CI 1300-2250)

    lq2_deaq <- log(37.3)
    label("Apparent desethylamodiaquine inter-compartmental clearance Q2/F to the deep peripheral compartment at WT = 45 kg (L/h)")
    # Ding 2026 Table 3: Q2/F DEAQ = 37.3 L/h (%RSE 11.5; SIR median 37.5,
    # 95% CI 29.6-45.7)

    lvp2_deaq <- log(6440)
    label("Apparent desethylamodiaquine deep peripheral volume of distribution Vp2/F at WT = 45 kg (L)")
    # Ding 2026 Table 3: Vp2/F DEAQ = 6440 L (%RSE 7.9; SIR median 6420,
    # 95% CI 5570-7590)

    # ---- Relative bioavailability ------------------------------------
    # Ding 2026 Methods ('Population PK analysis'): "Relative
    # bioavailability (F) was fixed to unity in the population, allowing
    # for quantification of the IIV in the absorption process." Here the
    # absorption variability is carried entirely as inter-occasion
    # variability (Table 3 reports IOV, not IIV, on F).
    lfdepot <- fixed(log(1))
    label("Relative oral bioavailability of amodiaquine (unitless)")
    # Ding 2026 Table 3: F = 1 Fix

    # ---- Allometric exponents ----------------------------------------
    e_wt_cl <- fixed(0.75)
    label("Allometric WT exponent on all apparent clearance parameters (CL/F and Q/F of both analytes)")
    # Ding 2026 Methods, 'Covariates model': clearance exponent fixed 0.75

    e_wt_vc <- fixed(1.00)
    label("Allometric WT exponent on all apparent volume parameters (Vc/F and Vp/F of both analytes)")
    # Ding 2026 Methods, 'Covariates model': volume exponent fixed 1.0

    # ---- Inter-occasion variability ----------------------------------
    # Ding 2026 Table 3 footnote: "Coefficients of variation for
    # interindividual variability and interoccasion variability (IIV and
    # IOV) were calculated as 100 x (e^variance - 1)^(1/2)", so the
    # internal log-scale variance is recovered as omega^2 = log(CV^2 + 1).
    #
    #   Ka  IOV 254%  -> omega^2 = log(2.540^2 + 1) = 2.0084288
    #   F   IOV 27.5% -> omega^2 = log(0.275^2 + 1) = 0.0729019
    #
    # Six occasions, one per artemether-lumefantrine dose time. NONMEM fits
    # a single IOV magnitude shared across occasions via $OMEGA BLOCK(1)
    # SAME, which maps to one estimated slot followed by fixed() repeats.
    # The absorption-rate IOV is extreme (254% CV) because the dense-PK
    # design samples the amodiaquine absorption phase in only 39 patients;
    # the reported eta shrinkage on this slot is 68.4%.
    etaiov_ka_1 ~ 2.0084288
    # Ding 2026 Table 3: IOV on Ka = 254% CV (%RSE 29.8; SIR median 261,
    # 95% CI 157-560; eta shrinkage 68.4%)
    etaiov_ka_2 ~ fixed(2.0084288)
    etaiov_ka_3 ~ fixed(2.0084288)
    etaiov_ka_4 ~ fixed(2.0084288)
    etaiov_ka_5 ~ fixed(2.0084288)
    etaiov_ka_6 ~ fixed(2.0084288)

    etaiov_fdepot_1 ~ 0.0729019
    # Ding 2026 Table 3: IOV on F = 27.5% CV (%RSE 15.1; SIR median 27.5,
    # 95% CI 23.2-31.3; eta shrinkage 63.2%)
    etaiov_fdepot_2 ~ fixed(0.0729019)
    etaiov_fdepot_3 ~ fixed(0.0729019)
    etaiov_fdepot_4 ~ fixed(0.0729019)
    etaiov_fdepot_5 ~ fixed(0.0729019)
    etaiov_fdepot_6 ~ fixed(0.0729019)

    # ---- Inter-individual variability --------------------------------
    #   CL/F AQ    IIV 11.6% -> omega^2 = log(0.116^2 + 1) = 0.0133663
    #   CL/F DEAQ  IIV 30.6% -> omega^2 = log(0.306^2 + 1) = 0.0895079
    #   Vc/F DEAQ  IIV 42.5% -> omega^2 = log(0.425^2 + 1) = 0.1660440
    #
    # Table 3 reports no IIV on Ka, Vc/F AQ, Q/F AQ, Vp/F AQ, Q1/F DEAQ,
    # Vp1/F DEAQ, Q2/F DEAQ or Vp2/F DEAQ, so no eta slots are created for
    # those parameters.
    etalcl ~ 0.0133663
    # Ding 2026 Table 3: IIV on CL/F AQ = 11.6% CV (%RSE 35.8; SIR median
    # 11.8, 95% CI 7.0-15.4; eta shrinkage 74.0%)

    etalcl_deaq ~ 0.0895079
    # Ding 2026 Table 3: IIV on CL/F DEAQ = 30.6% CV (%RSE 21.3; SIR median
    # 30.7, 95% CI 24.0-37.5; eta shrinkage 37.8%)

    etalvc_deaq ~ 0.1660440
    # Ding 2026 Table 3: IIV on Vc/F DEAQ = 42.5% CV (%RSE 32.2; SIR median
    # 42.3, 95% CI 28.2-57.1; eta shrinkage 68.8%)

    # ---- Residual unexplained variability ----------------------------
    # Ding 2026 Methods ('Population PK analysis'): the residual was
    # "modelled as an additive error on log-transformed concentrations,
    # which is approximately equivalent to an exponential residual error on
    # an arithmetic scale", which maps to a proportional residual in linear
    # concentration space. The Table 3 footnote states "RUV is the residual
    # error variance", so the tabulated number is a variance and the SD is
    # its square root -- the same convention as the sibling
    # Ding_2024_amodiaquine.R.
    propSd <- sqrt(0.0676)
    label("Proportional residual SD for amodiaquine plasma concentration (SD on log scale)")
    # Ding 2026 Table 3: RUV AQ = 0.0676 (variance; %RSE 5.0; SIR median
    # 0.0676, 95% CI 0.0561-0.0827; epsilon shrinkage 15.5%)

    propSd_deaq <- sqrt(0.114)
    label("Proportional residual SD for desethylamodiaquine plasma concentration (SD on log scale)")
    # Ding 2026 Table 3: RUV DEAQ = 0.114 (variance; %RSE 3.7; SIR 95% CI
    # 0.0996-0.132; epsilon shrinkage 15.3%)
  })

  model({
    # Molecular weights of the free bases (g/mol). Ding 2026 Methods
    # ('Population PK analysis'): "Parent drugs were assumed to be
    # completely metabolized to their metabolites due to identifiability
    # issues with other model structures." The conversion factor is not
    # printed, so the mass flux leaving amodiaquine central is
    # molar-corrected before it enters desethylamodiaquine central,
    # matching the sibling Ding_2024_amodiaquine.R from the same group and
    # the WWARN amodiaquine model Ali_2018_amodiaquine.R, which states the
    # molar correction explicitly. The paper's own secondary parameters
    # support the correction: the published AUC ratio
    # AUC_DEAQ / AUC_AQ = 96.5 h*ug/mL / 1.53 h*ug/mL = 63.1 (a ratio that
    # cancels the dose and F entirely) is reproduced as
    # (CL_AQ / CL_DEAQ) * molarFactor = (2250 / 32.2) * 0.9212 = 64.4,
    # against 69.9 for an uncorrected mass-for-mass reading.
    mwAQ        <- 355.85
    mwDEAQ      <- 327.81
    molarFactor <- mwDEAQ / mwAQ

    # Occasion indicators for the inter-occasion variability slots
    # (Ding 2026 Results 3.1.2: IOV on relative bioavailability and on the
    # absorption rate constant).
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    oc5 <- (OCC == 5)
    oc6 <- (OCC == 6)

    iov_ka <- oc1 * etaiov_ka_1 + oc2 * etaiov_ka_2 + oc3 * etaiov_ka_3 +
              oc4 * etaiov_ka_4 + oc5 * etaiov_ka_5 + oc6 * etaiov_ka_6
    iov_fdepot <- oc1 * etaiov_fdepot_1 + oc2 * etaiov_fdepot_2 + oc3 * etaiov_fdepot_3 +
                  oc4 * etaiov_fdepot_4 + oc5 * etaiov_fdepot_5 + oc6 * etaiov_fdepot_6

    ka <- exp(lka + iov_ka)

    # Individual PK parameters. Allometric weight scaling on all apparent
    # clearances (exponent 0.75) and apparent volumes (exponent 1) centred
    # on the 45 kg reference of the Table 3 footnote.
    cl <- exp(lcl + etalcl) * (WT / 45)^e_wt_cl
    vc <- exp(lvc)          * (WT / 45)^e_wt_vc
    q  <- exp(lq)           * (WT / 45)^e_wt_cl
    vp <- exp(lvp)          * (WT / 45)^e_wt_vc

    cl_deaq  <- exp(lcl_deaq + etalcl_deaq) * (WT / 45)^e_wt_cl
    vc_deaq  <- exp(lvc_deaq + etalvc_deaq) * (WT / 45)^e_wt_vc
    q_deaq   <- exp(lq_deaq)                * (WT / 45)^e_wt_cl
    vp_deaq  <- exp(lvp_deaq)               * (WT / 45)^e_wt_vc
    q2_deaq  <- exp(lq2_deaq)               * (WT / 45)^e_wt_cl
    vp2_deaq <- exp(lvp2_deaq)              * (WT / 45)^e_wt_vc

    # Micro-rate constants (1/h).
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    kel_deaq <- cl_deaq  / vc_deaq
    k12_deaq <- q_deaq   / vc_deaq
    k21_deaq <- q_deaq   / vp_deaq
    k13_deaq <- q2_deaq  / vc_deaq
    k31_deaq <- q2_deaq  / vp2_deaq

    # ODE system (Ding 2026 Figure S5). Compartment amounts are in mg of
    # analyte base and volumes are in L, so amount/volume is mg/L and is
    # scaled to ng/mL below.
    d/dt(depot) <- -ka * depot

    # Amodiaquine central plus one peripheral compartment. The entire mass
    # flux leaving amodiaquine central by elimination (kel * central) is
    # routed to desethylamodiaquine central under the complete-conversion
    # assumption, after the molar correction.
    d/dt(central)     <- ka * depot - kel * central -
                         k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Desethylamodiaquine central plus two peripheral compartments.
    d/dt(central_deaq)     <- molarFactor * kel * central -
                              kel_deaq * central_deaq -
                              k12_deaq * central_deaq + k21_deaq * peripheral1_deaq -
                              k13_deaq * central_deaq + k31_deaq * peripheral2_deaq
    d/dt(peripheral1_deaq) <- k12_deaq * central_deaq - k21_deaq * peripheral1_deaq
    d/dt(peripheral2_deaq) <- k13_deaq * central_deaq - k31_deaq * peripheral2_deaq

    # Relative oral bioavailability on the amodiaquine dose. The population
    # anchor lfdepot = log(1) is fixed; all absorption variability is
    # carried as inter-occasion variability.
    f(depot) <- exp(lfdepot + iov_fdepot)

    # Plasma concentrations in ng/mL: amount (mg) / volume (L) is mg/L,
    # multiplied by 1000 to give ng/mL, the units used throughout Ding 2026
    # Table 3 and Figure 3.
    Cc      <- 1000 * central      / vc
    Cc_deaq <- 1000 * central_deaq / vc_deaq

    # Proportional residual error on the linear-concentration scale, the
    # linear-space equivalent of the paper's additive-on-log-scale error.
    Cc      ~ prop(propSd)
    Cc_deaq ~ prop(propSd_deaq)
  })
}
