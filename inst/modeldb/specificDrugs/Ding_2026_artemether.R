# Joint parent-metabolite population PK model for oral artemether and its
# active metabolite dihydroartemisinin in patients with uncomplicated
# Plasmodium falciparum malaria, pooled from the TRACII (NCT02453308) and
# TACT-CV (NCT03355664) trials (Ding 2026, Br J Clin Pharmacol
# 92(2):589-605; doi:10.1002/bcp.70301).

Ding_2026_artemether <- function() {
  description <- paste(
    "Joint parent-metabolite population PK model for oral artemether and",
    "its active metabolite dihydroartemisinin in adults and children with",
    "acute uncomplicated Plasmodium falciparum malaria, given the standard",
    "six-dose artemether-lumefantrine regimen alone or together with",
    "amodiaquine as a triple artemisinin-based combination therapy",
    "(Ding 2026, pooled TRACII + TACT-CV dense-PK cohorts, n = 79).",
    "Two-transit-compartment absorption feeds a two-compartment artemether",
    "disposition model whose apparent clearance rises linearly with dose",
    "occasion (an empirical description of CYP2B6 autoinduction), with",
    "complete molar-corrected bioconversion to a one-compartment",
    "dihydroartemisinin disposition model. Allometric body-weight scaling",
    "on all apparent clearances (fixed exponent 0.75) and apparent volumes",
    "(fixed exponent 1.0) at a reference weight of 45 kg. Relative",
    "bioavailability is anchored at 1 with inter-occasion variability on",
    "both bioavailability and mean transit time. Coadministered amodiaquine",
    "was not a significant covariate on any parameter. Predictions are",
    "plasma artemether and dihydroartemisinin concentrations in ng/mL.",
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
  # scheme', 'Drug quantification': venous plasma artemether and
  # dihydroartemisinin by LC-MS/MS) and Figure S1.
  compartmentData <- list(
    depot        = list(analyte = "artemether",         units = "mg", specimen = "administration site", verified = TRUE),
    transit1     = list(analyte = "artemether",         units = "mg", specimen = "administration site", verified = TRUE),
    transit2     = list(analyte = "artemether",         units = "mg", specimen = "administration site", verified = TRUE),
    central      = list(analyte = "artemether",         units = "mg", specimen = "plasma",              verified = TRUE),
    peripheral1  = list(analyte = "artemether",         units = "mg", specimen = "plasma",              verified = TRUE),
    central_dihydroart  = list(analyte = "dihydroartemisinin", units = "mg", specimen = "plasma",              verified = TRUE)
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
        "0.75 and 1.0, respectively'. The Table 2 footnote fixes the",
        "reference: 'Population estimates are given for a typical adult",
        "patient weighing 45 kg with acute P. falciparum malaria', so 45 kg",
        "is the normalising constant encoded here. Dense-PK cohort median",
        "50.0-52.3 kg (Table S3); full-trial range 9.0-98.8 kg (Table 1).",
        "Retained despite only a marginal improvement in fit",
        "(Results 3.1.1: delta-OFV = -0.826) 'due to the strong biological",
        "basis for this covariate'.",
        sep = " "
      ),
      source_name        = "BW"
    ),
    OCC = list(
      description        = "Dose occasion index, 1 to 6 across the six-dose artemether-lumefantrine regimen",
      units              = "(count)",
      type               = "categorical",
      reference_category = "1 (the first dose, at hour 0)",
      notes              = paste(
        "Integer occasion column taking value k on the interval starting at",
        "the k-th dose. Doses are given at 0, 8, 24, 36, 48 and 60 h",
        "(Ding 2026 Methods, 'Dosing regimen'; Table S1), so OCC = 1 for",
        "0 <= t < 8 h, 2 for 8 <= t < 24 h, 3 for 24 <= t < 36 h, 4 for",
        "36 <= t < 48 h, 5 for 48 <= t < 60 h and 6 for t >= 60 h. The",
        "column serves two distinct purposes in this model, which is why it",
        "is a covariate rather than a purely random-effect grouping: (i) the",
        "FIXED time-dependent clearance effect of the Table 2 footnote,",
        "'Time dependent CL is modelled as [(1 + theta*[OCC-1])*CL], where",
        "OCC is the dose occasion from 1 to 6', and (ii) the occasion",
        "grouping for the inter-occasion variability on mean transit time",
        "and relative bioavailability. Note that the time-dependent",
        "clearance term is explicitly EMPIRICAL and not extrapolable:",
        "Results 3.1.1 states 'this empirical time-dependent clearance model",
        "is not suitable for extrapolation to treatment durations beyond the",
        "standard 3 days'. Holding OCC at 6 beyond 60 h (as encoded in the",
        "companion vignette) freezes clearance at its last observed value",
        "rather than extrapolating the linear rise.",
        sep = " "
      ),
      source_name        = "OCC"
    )
  )

  covariatesDataExcluded <- list(
    CONMED_AMODIAQUINE = list(
      description = "Coadministration of amodiaquine (triple ACT versus artemether-lumefantrine alone)",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened as a binary drug-drug-interaction covariate on every PK",
        "parameter of artemether and dihydroartemisinin, and additionally",
        "assessed by a 500-bootstrap full covariate model, but NOT retained:",
        "Ding 2026 Results 3.1.1, 'Coadministration of amodiaquine did not",
        "affect the PK properties of artemether or dihydroartemisinin. This",
        "was further confirmed by the full covariate approach, which showed",
        "that the relative change associated with drug-drug interactions",
        "included zero for the main primary PK parameters (Figure S2)'. The",
        "final model therefore has no amodiaquine term, and this single",
        "model serves both the artemether-lumefantrine and the",
        "artemether-lumefantrine-amodiaquine arms.",
        sep = " "
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 79L,
    n_studies      = 2L,
    n_observations = paste(
      "626 plasma samples in total, including 42 artemether and 101",
      "dihydroartemisinin measurements below the lower limit of",
      "quantification, which were discarded (Beal M1; Results 3.1.1)"
    ),
    age_range      = "26.5-30.0 years, median by study/arm (range 8.0-58.4 years) (Table S3)",
    weight_range   = "50.0-52.3 kg, median by study/arm (range 18.0-78.0 kg) (Table S3)",
    sex_female_pct = 10.1,
    disease_state  = paste(
      "Acute uncomplicated Plasmodium falciparum malaria. Dense-PK cohort",
      "median admission asexual parasite count 5,597-45,000 parasites/uL and",
      "median admission body temperature 37.3-38.2 degC (Table S3)."
    ),
    dose_range     = paste(
      "Standard fixed-dose artemether-lumefantrine (20 mg artemether +",
      "120 mg lumefantrine per tablet) given orally as six doses over 3 days",
      "at 0, 8, 24, 36, 48 and 60 h, directly observed, with a fatty snack",
      "(TRACII) or 80 mL milk (TACT-CV). Tablets per dose by weight band:",
      "1 (5-14.9 kg), 2 (15-24.9 kg), 3 (25-34.9 kg), 4 (>35 kg)",
      "(Table S1). Dense-PK cohort median artemether dose 3.1-3.2",
      "mg/kg/day (range 2.1-4.4) (Table S3)."
    ),
    regions        = paste(
      "TRACII (NCT02453308) dense-PK sampling at one site in Bangladesh",
      "(n = 41); TACT-CV (NCT03355664) dense-PK sampling at one site in",
      "Vietnam (n = 38)."
    ),
    notes          = paste(
      "Only the dense-PK sub-cohorts contribute artemether and",
      "dihydroartemisinin data: samples at 1, 2, 4, 6, 8, 12, 24 and 64 h",
      "(plus 52 h in TRACII), because artemether and dihydroartemisinin were",
      "not quantified from Day 4 onwards given their short half-lives",
      "(Methods, 'PK sampling scheme'). Children below 20 kg were excluded",
      "from dense PK sampling. All patients except those at the Democratic",
      "Republic of Congo sites also received a single 0.25 mg/kg",
      "gametocytocidal dose of primaquine 24 h after the start of treatment."
    )
  )

  ini({
    # ---- Absorption --------------------------------------------------
    # Ding 2026 Results 3.1.1: "The absorption phase of artemether was best
    # described with two transit compartments, with an identical rate
    # constant between compartments." Table 2 reports a mean transit time
    # and a fixed transit-compartment count but NO separate absorption rate
    # constant, so the single rate governs all three transfers
    # (depot -> transit1 -> transit2 -> central) and the mean transit time
    # spans all three: ktr = (NN + 1) / MTT = 3 / MTT. This is the same
    # reading as the sibling MORU models Hoglund_2017_piperaquine.R and
    # Ali_2018_amodiaquine.R, and differs from Ding_2024_amodiaquine.R
    # only because that model DOES report a separate Ka. The reading is
    # confirmed numerically: ktr = 3/MTT reproduces the published typical
    # artemether Cmax to 6% (240 vs 256 ng/mL, Tmax 2.0 h), whereas
    # ktr = 2/MTT misses it by 25% (192 ng/mL). See the vignette source
    # trace.
    lmtt <- log(1.55)
    label("Mean transit time of the artemether absorption chain (h)")
    # Ding 2026 Table 2: Mean transit time = 1.55 h (%RSE 7.3; SIR median
    # 1.55, 95% CI 1.33-1.77)

    # ---- Artemether disposition (two compartments) -------------------
    # Ding 2026 Table 2, "NONMEM Estimates" column. Values are apparent
    # (relative to F = 1) and reported on the linear scale for a typical
    # adult patient weighing 45 kg; log() is applied here for the nlmixr2
    # internal log scale.
    lcl <- log(79.3)
    label("Apparent artemether elimination clearance CL/F on the first dose occasion at WT = 45 kg (L/h)")
    # Ding 2026 Table 2: CL/F ARM = 79.3 L/h (%RSE 4.8; SIR median 79.2,
    # 95% CI 71.4-86.6; eta shrinkage 14.1%)

    lvc <- log(141)
    label("Apparent artemether central volume of distribution Vc/F at WT = 45 kg (L)")
    # Ding 2026 Table 2: Vc/F ARM = 141 L (%RSE 7.1; SIR median 140,
    # 95% CI 122-161)

    lq <- log(21.7)
    label("Apparent artemether inter-compartmental clearance Q/F at WT = 45 kg (L/h)")
    # Ding 2026 Table 2: Q/F ARM = 21.7 L/h (%RSE 9.32; SIR median 21.5,
    # 95% CI 18.2-26.1)

    lvp <- log(283)
    label("Apparent artemether peripheral volume of distribution Vp/F at WT = 45 kg (L)")
    # Ding 2026 Table 2: Vp/F ARM = 283 L (%RSE 23; SIR median 279,
    # 95% CI 203-446)

    # ---- Dihydroartemisinin disposition (one compartment) ------------
    lcl_dihydroart <- log(255)
    label("Apparent dihydroartemisinin elimination clearance CL/F at WT = 45 kg (L/h)")
    # Ding 2026 Table 2: CL/F DHA = 255 L/h (%RSE 5.9; SIR median 254,
    # 95% CI 229-288; eta shrinkage 58.0%)

    lvc_dihydroart <- log(64.1)
    label("Apparent dihydroartemisinin central volume of distribution Vc/F at WT = 45 kg (L)")
    # Ding 2026 Table 2: Vc/F DHA = 64.1 L (%RSE 15.3; SIR median 64.3,
    # 95% CI 44.7-82.7; eta shrinkage 68.1%)

    # ---- Relative bioavailability ------------------------------------
    # Ding 2026 Methods ('Population PK analysis'): "Relative
    # bioavailability (F) was fixed to unity in the population, allowing
    # for quantification of the IIV in the absorption process." Here the
    # absorption variability is carried entirely as inter-occasion
    # variability (Table 2 reports IOV, not IIV, on F).
    lfdepot <- fixed(log(1))
    label("Relative oral bioavailability of artemether (unitless)")
    # Ding 2026 Table 2: F = 1 Fix

    # ---- Allometric exponents ----------------------------------------
    # Ding 2026 Methods ('Covariates model'): "bodyweight was included on
    # all clearance and volume parameters using a conventional allometric
    # function with fixed exponents of 0.75 and 1.0, respectively".
    e_wt_cl <- fixed(0.75)
    label("Allometric WT exponent on all apparent clearance parameters (CL/F ARM, Q/F ARM, CL/F DHA)")
    # Ding 2026 Methods, 'Covariates model': clearance exponent fixed 0.75

    e_wt_vc <- fixed(1.00)
    label("Allometric WT exponent on all apparent volume parameters (Vc/F ARM, Vp/F ARM, Vc/F DHA)")
    # Ding 2026 Methods, 'Covariates model': volume exponent fixed 1.0

    # ---- Time-dependent artemether clearance -------------------------
    # Ding 2026 Table 2 footnote gives the form verbatim: "Time dependent
    # CL is modelled as [(1 + theta*[OCC-1])*CL], where OCC is the dose
    # occasion from 1 to 6." Discussion: "we introduced a time-dependent
    # empirical parameter on clearance to describe the time-varying
    # clearance of artemether due to autoinduction of liver enzyme."
    e_occ_cl <- 0.551
    label("Linear fractional increase in apparent artemether clearance per dose occasion after the first (per occasion)")
    # Ding 2026 Table 2: Time dependency on CL = 0.551 (%RSE 16.7; SIR
    # median 0.546, 95% CI 0.387-0.736)

    # ---- Inter-occasion variability ----------------------------------
    # Ding 2026 Table 2 footnote: "Coefficients of inter-individual and
    # inter-occasion variability (IIV and IOV) were calculated as
    # 100 x (e^variance - 1)^(1/2)", so the internal log-scale variance is
    # recovered as omega^2 = log(CV^2 + 1).
    #
    #   MTT  IOV 54.4% -> omega^2 = log(0.544^2 + 1) = 0.2592332
    #   F    IOV 31.6% -> omega^2 = log(0.316^2 + 1) = 0.0951793
    #
    # Six occasions, one per dose (Table 2 footnote: 'OCC is the dose
    # occasion from 1 to 6'). NONMEM fits a single IOV magnitude shared
    # across occasions via $OMEGA BLOCK(1) SAME, which maps to one
    # estimated slot followed by fixed() repeats.
    etaiov_mtt_1 ~ 0.2592332
    # Ding 2026 Table 2: IOV on mean transit time = 54.4% CV (%RSE 16.3;
    # SIR median 54.2, 95% CI 46.5-65.7)
    etaiov_mtt_2 ~ fixed(0.2592332)
    etaiov_mtt_3 ~ fixed(0.2592332)
    etaiov_mtt_4 ~ fixed(0.2592332)
    etaiov_mtt_5 ~ fixed(0.2592332)
    etaiov_mtt_6 ~ fixed(0.2592332)

    etaiov_fdepot_1 ~ 0.0951793
    # Ding 2026 Table 2: IOV on F = 31.6% CV (%RSE 22.9; SIR median 31.7,
    # 95% CI 24.7-39.6; eta shrinkage 22.3%)
    etaiov_fdepot_2 ~ fixed(0.0951793)
    etaiov_fdepot_3 ~ fixed(0.0951793)
    etaiov_fdepot_4 ~ fixed(0.0951793)
    etaiov_fdepot_5 ~ fixed(0.0951793)
    etaiov_fdepot_6 ~ fixed(0.0951793)

    # ---- Inter-individual variability --------------------------------
    #   CL/F ARM  IIV 21.1% -> omega^2 = log(0.211^2 + 1) = 0.0435584
    #   CL/F DHA  IIV 41.7% -> omega^2 = log(0.417^2 + 1) = 0.1603222
    #
    # Table 2 reports no IIV on MTT, Vc/F ARM, Q/F ARM, Vp/F ARM or
    # Vc/F DHA, so no eta slots are created for those parameters.
    etalcl ~ 0.0435584
    # Ding 2026 Table 2: IIV on CL/F ARM = 21.1% CV (%RSE 29.1; SIR median
    # 21.5, 95% CI 16.1-27.7; eta shrinkage 14.1%)

    etalcl_dihydroart ~ 0.1603222
    # Ding 2026 Table 2: IIV on CL/F DHA = 41.7% CV (%RSE 22.2; SIR median
    # 42.2, 95% CI 33.4-52.0; eta shrinkage 58.0%)

    # ---- Residual unexplained variability ----------------------------
    # Ding 2026 Methods ('Population PK analysis'): "The residual
    # unexplained variability, assumed to be normally distributed with a
    # zero mean and variance sigma^2, was modelled as an additive error on
    # log-transformed concentrations, which is approximately equivalent to
    # an exponential residual error on an arithmetic scale." That
    # additive-on-log-scale residual maps to a proportional residual in
    # linear concentration space. The Table 2 footnote states "RUV is the
    # residual error variance", so the tabulated number is a variance and
    # the SD is its square root -- the same convention as the sibling
    # Ding_2024_amodiaquine.R and Hoglund_2017_piperaquine.R.
    propSd <- sqrt(0.229)
    label("Proportional residual SD for artemether plasma concentration (SD on log scale)")
    # Ding 2026 Table 2: RUV ARM = 0.229 (variance; %RSE 3.8; SIR median
    # 0.230, 95% CI 0.198-0.265; epsilon shrinkage 14.7%)

    propSd_dihydroart <- sqrt(0.262)
    label("Proportional residual SD for dihydroartemisinin plasma concentration (SD on log scale)")
    # Ding 2026 Table 2: RUV DHA = 0.262 (variance; %RSE 4.0; SIR median
    # 0.261, 95% CI 0.226-0.306; epsilon shrinkage 14.3%)
  })

  model({
    # Molecular weights (g/mol). Ding 2026 Methods ('Population PK
    # analysis') states that "Parent drugs were assumed to be completely
    # metabolized to their metabolites due to identifiability issues with
    # other model structures" but does not print the conversion factor. The
    # mass flux leaving artemether central is therefore molar-corrected
    # before it enters dihydroartemisinin central, matching the sibling
    # Ding_2024_amodiaquine.R from the same group. The paper's own
    # secondary parameters confirm the correction: with complete molar
    # conversion the published typical dihydroartemisinin AUC is
    # reproduced as Dose * (MW_DHA / MW_ARM) / CL_DHA = 480 * 0.95299 /
    # 255 = 1794 h*ng/mL against a published 1870, whereas an uncorrected
    # mass-for-mass reading gives 1882 h*ng/mL. Both are close because the
    # molar correction is only 4.7% here; the correction is applied for
    # consistency with the amodiaquine and lumefantrine siblings, where it
    # is decisive.
    mwARM       <- 298.38
    mwDHA       <- 284.35
    molarFactor <- mwDHA / mwARM

    # Occasion indicators (Ding 2026 Table 2 footnote: 'OCC is the dose
    # occasion from 1 to 6'). They drive both the inter-occasion
    # variability slots and the fixed time-dependent clearance term.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    oc5 <- (OCC == 5)
    oc6 <- (OCC == 6)

    iov_mtt <- oc1 * etaiov_mtt_1 + oc2 * etaiov_mtt_2 + oc3 * etaiov_mtt_3 +
               oc4 * etaiov_mtt_4 + oc5 * etaiov_mtt_5 + oc6 * etaiov_mtt_6
    iov_fdepot <- oc1 * etaiov_fdepot_1 + oc2 * etaiov_fdepot_2 + oc3 * etaiov_fdepot_3 +
                  oc4 * etaiov_fdepot_4 + oc5 * etaiov_fdepot_5 + oc6 * etaiov_fdepot_6

    # Absorption chain rate constant. NN = 2 transit compartments and the
    # single rate governs all NN + 1 = 3 transfers, so ktr = 3 / MTT (see
    # the lmtt annotation in ini()).
    mtt <- exp(lmtt + iov_mtt)
    ktr <- 3 / mtt

    # Individual PK parameters. Allometric weight scaling on all apparent
    # clearances (exponent 0.75) and apparent volumes (exponent 1) centred
    # on the 45 kg reference of the Table 2 footnote. Artemether clearance
    # additionally rises linearly with dose occasion.
    cl <- exp(lcl + etalcl) * (WT / 45)^e_wt_cl * (1 + e_occ_cl * (OCC - 1))
    vc <- exp(lvc)          * (WT / 45)^e_wt_vc
    q  <- exp(lq)           * (WT / 45)^e_wt_cl
    vp <- exp(lvp)          * (WT / 45)^e_wt_vc

    cl_dihydroart <- exp(lcl_dihydroart + etalcl_dihydroart) * (WT / 45)^e_wt_cl
    vc_dihydroart <- exp(lvc_dihydroart)              * (WT / 45)^e_wt_vc

    # Micro-rate constants (1/h).
    kel     <- cl / vc
    k12     <- q  / vc
    k21     <- q  / vp
    kel_dihydroart <- cl_dihydroart / vc_dihydroart

    # ODE system (Ding 2026 Figure S1). Compartment amounts are in mg of
    # analyte and volumes are in L, so amount/volume is mg/L and is scaled
    # to ng/mL below.
    #
    # Absorption: depot -> transit1 -> transit2 -> artemether central, all
    # three transfers at the same rate ktr.
    d/dt(depot)    <- -ktr * depot
    d/dt(transit1) <-  ktr * depot    - ktr * transit1
    d/dt(transit2) <-  ktr * transit1 - ktr * transit2

    # Artemether central plus one peripheral compartment. The entire mass
    # flux leaving artemether central by elimination (kel * central) is
    # routed to dihydroartemisinin central under the complete-conversion
    # assumption, after the molar correction.
    d/dt(central)     <- ktr * transit2 - kel * central -
                         k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Dihydroartemisinin central, one compartment. Its disposition is
    # formation-rate limited, which is why Ding 2026 Table 2 reports no
    # terminal half-life for dihydroartemisinin.
    d/dt(central_dihydroart) <- molarFactor * kel * central - kel_dihydroart * central_dihydroart

    # Relative oral bioavailability on the artemether dose. The population
    # anchor lfdepot = log(1) is fixed; all absorption variability is
    # carried as inter-occasion variability.
    f(depot) <- exp(lfdepot + iov_fdepot)

    # Plasma concentrations in ng/mL: amount (mg) / volume (L) is mg/L,
    # multiplied by 1000 to give ng/mL, the units used throughout Ding 2026
    # Table 2 and Figure 1.
    Cc     <- 1000 * central     / vc
    Cc_dihydroart <- 1000 * central_dihydroart / vc_dihydroart

    # Proportional residual error on the linear-concentration scale, the
    # linear-space equivalent of the paper's additive-on-log-scale error.
    Cc     ~ prop(propSd)
    Cc_dihydroart ~ prop(propSd_dihydroart)
  })
}
