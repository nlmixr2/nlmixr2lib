Zhang_2012_lopinavir_ritonavir_pedi <- function() {
  description <- paste(
    "Integrated one-compartment popPK model for oral lopinavir (LPV) and",
    "ritonavir (RTV) in 74 HIV-infected children (6 months to 4.5 years)",
    "treated with LPV/r oral solution with or without concomitant",
    "rifampicin-based antitubercular treatment (Zhang 2012). LPV uses a",
    "one-compartment model with first-order absorption; RTV uses a one-",
    "compartment model with a Savic-style 10-transit-compartment",
    "absorption chain followed by a separate first-order absorption step",
    "from the last transit to central. Apparent CL/F and V/F are",
    "allometrically scaled to the cohort median 10 kg with fixed exponents",
    "0.75 / 1. The dynamic LPV-RTV interaction is encoded as direct",
    "sigmoid-Emax inhibition of LPV apparent clearance by RTV plasma",
    "concentration (Emax = 0.9 fixed, EC50 = 0.0519 mg/L). Lopinavir",
    "bioavailability is modulated by concomitant rifampicin-based",
    "antitubercular treatment (-83.2% at the no-extra-ritonavir reference)",
    "and by the concomitant ritonavir dose in mg/kg (+2.1% per mg/kg",
    "above the 3 mg/kg reference). Ritonavir apparent clearance is",
    "+50% in subjects on rifampicin-based antitubercular treatment.",
    "Both drugs share random effects modelled as log-normal between-",
    "subject variability with selected inter-occasion variabilities",
    "folded in as BSV-equivalent (see vignette Assumptions and",
    "deviations). Residual error is proportional on the linear scale",
    "(implemented via NONMEM exponential error on log-transformed data)."
  )
  reference <- paste(
    "Zhang C, McIlleron H, Ren Y, van der Walt JS, Karlsson MO,",
    "Simonsson USH, Denti P.",
    "Population pharmacokinetics of lopinavir and ritonavir in combination",
    "with rifampicin-based antitubercular treatment in HIV-infected children.",
    "Antivir Ther. 2012;17(1):25-33.",
    "doi:10.3851/IMP1915.",
    sep = " "
  )
  vignette <- "Zhang_2012_lopinavir_ritonavir_pedi"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Used for allometric scaling at the cohort median 10 kg",
        "(Zhang 2012 Methods Equations 1 and 2 / Results 'Patients and",
        "data description' / Table 1: median body weight 10.2 kg, range",
        "5-17 kg). Exponent 0.75 on apparent CL/F (both LPV and RTV),",
        "exponent 1.0 on apparent V/F (both LPV and RTV). The Methods",
        "section 'Population pharmacokinetic analysis' paragraph 4",
        "states 'In order to account for size differences, allometric",
        "scaling based on the median body weight was tested and applied",
        "to apparent clearance (CL/F) and volume of distribution (V/F)'",
        "with the 10 kg reference appearing as the denominator of the",
        "weight ratio in Equations 1 and 2."
      ),
      source_name        = "WT"
    ),
    CONMED_RIF = list(
      description        = "Concomitant rifampicin-based antitubercular treatment indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no rifampicin-based antitubercular treatment)",
      notes              = paste(
        "1 = subject on rifampicin-based antitubercular treatment;",
        "0 = subject not on rifampicin. Set per occasion. Two effects",
        "in this model: (1) multiplicative -83.2% reduction in lopinavir",
        "relative bioavailability via the linear formula",
        "F_LPV = 1 - 0.832 * CONMED_RIF + 0.021 * (DOSE_RTV_MGKG - 3)",
        "(Zhang 2012 Equation 4 and Results 'Model description'",
        "paragraph 2); (2) +50% increase in ritonavir apparent clearance",
        "(12.8 L/h without rifampicin vs 19.1 L/h with rifampicin per",
        "Table 2 row 'Ritonavir CL/F'). The paper assumes the",
        "rifampicin enzyme-induction effect on CYP3A is at steady state",
        "after at least 2 weeks of co-administration so the indicator",
        "is treated as a steady-state binary switch with no induction",
        "lag (Methods 'Population pharmacokinetic analysis' paragraph 5:",
        "'the effect of rifampicin on enzyme induction could be assumed",
        "to be at steady state and within-day change could be",
        "neglected')."
      ),
      source_name        = "RIF"
    ),
    DOSE_RTV_MGKG = list(
      description        = "Concomitant ritonavir per-administration dose per kg body weight",
      units              = "mg/kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Per-dose ritonavir mg/kg given concomitantly with the lopinavir",
        "dose. Enters lopinavir bioavailability via the linear-shift",
        "regressor of Zhang 2012 Equation 4:",
        "F_LPV = 1 - 0.832 * CONMED_RIF + 0.021 * (DOSE_RTV_MGKG - 3)",
        "with DOSE_RTV_STD = 3 mg/kg (the median ritonavir dose in the",
        "no-rifampicin standard-LPV/r 4:1 reference arm, per the",
        "Methods 'Population pharmacokinetic analysis' paragraph 5).",
        "The standard 4:1 LPV/r at the median 11.6 mg/kg LPV dose",
        "(Table 1) gives DOSE_RTV_MGKG = 2.9 (~ 3 mg/kg reference,",
        "yielding F_LPV ~ 1). The super-boosted 1:1 LPV/r at the",
        "median 14 mg/kg LPV dose gives DOSE_RTV_MGKG = 14 and",
        "F_LPV = 0.40 (cohort estimate 40.5%, Table 2). The",
        "double-dose 4:1 LPV/r at the median 23 mg/kg LPV dose gives",
        "DOSE_RTV_MGKG = 5.75 and F_LPV = 0.23 (cohort estimate 22.6%,",
        "Table 2). The linear approximation is valid only inside the",
        "cohort-tested ritonavir-dose range of 2.9-14 mg/kg",
        "(Discussion paragraph 3: 'the relationship between ritonavir",
        "dose and bioavailability is probably quite complicated, in",
        "our model it was described using a linear proportionality.",
        "This choice was compelled by the limited range of ritonavir",
        "doses available in the dataset, and should not be used too",
        "far outside the tested range.')."
      ),
      source_name        = "DoseRTV"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 74L,
    n_studies      = 2L,
    age_range      = "6 months to 4.5 years (paediatric)",
    age_median     = "21 months",
    weight_range   = "5-17 kg",
    weight_median  = "10.2 kg (cohort) / 10 kg (allometric-scaling reference, Zhang 2012 Methods Equations 1 and 2)",
    sex_female_pct = 22.7,
    race_ethnicity = "South African paediatric HIV cohort (not further detailed in the paper).",
    disease_state  = "HIV-1 infection on lopinavir/ritonavir-based combination antiretroviral therapy (paediatric). Three sub-cohorts: 39 children without tuberculosis receiving standard LPV/r 4:1 oral solution every 12 h (median LPV dose 11.6 mg/kg); 15 children with HIV-associated tuberculosis receiving 'super-boosted' LPV (LPV/r 1:1 with extra ritonavir added to the standard LPV/r 4:1) plus rifampicin-based antitubercular treatment; 20 children with HIV-associated tuberculosis receiving doubled standard LPV/r 4:1 dose every 12 h plus rifampicin-based antitubercular treatment; 11 of the 'double-dose' children were re-sampled at least 4 weeks after completion of antitubercular treatment on standard LPV/r doses (counted within the 39-without-tuberculosis sample because no rifampicin was active at that occasion).",
    dose_range     = "LPV/r oral solution 230/57.5 mg/m^2 every 12 h in the standard reference cohort (median LPV dose 11.6 mg/kg, range 9.4-16.0; Table 1). Super-boosted cohort: extra ritonavir added to standard 4:1 LPV/r to a 1:1 ratio (median LPV 14.0 mg/kg, range 10.7-18.0). Double-dose cohort: double the standard LPV/r every 12 h (median LPV 23.0 mg/kg, range 13.8-29.5). Antituberculous regimen contained rifampicin 10 mg/kg/day. Simulation-based proposed dose recommendations during rifampicin co-treatment are 27 / 21 / 20 / 18 mg/kg every 8 h in the WHO weight bands 3-5.9 / 6-9.9 / 10-13.9 / 14-19.9 kg (Discussion paragraph 4 and Table 3).",
    regions        = "Cape Town, Stellenbosch, and Witwatersrand catchments (South Africa).",
    notes          = paste(
      "Pooled cohort from three antiretroviral clinics in South Africa.",
      "All samples were taken after at least 2 weeks of concurrent",
      "antitubercular and antiretroviral therapy to ensure pharmacokinetic",
      "steady state. About 5% of LPV/RTV samples were below the lower",
      "limit of quantification (LPV LLOQ 0.05 mg/L, RTV LLOQ 0.025 mg/L)",
      "and excluded; care was taken to ensure model predictions for the",
      "excluded samples were compatibly low. A total of 216 + 120 + 96 =",
      "432 LPV and RTV concentration measurements were available across",
      "the standard, super-boosted, and double-dose cohorts respectively.",
      "Baseline demographics from Zhang 2012 Table 1: median age 21 months",
      "(6 months to 4.5 years), median body weight 10.2 kg (5-17 kg),",
      "gender 34 males / 10 females (sex-female pct calculated as 10/44 =",
      "22.7%; Table 1 reports 34/10 across the included demographic table",
      "rows -- the 30 remaining children of the 74 total are not",
      "disaggregated by sex in the paper). Median height 79 cm",
      "(58-103 cm), median BSA 0.48 m^2 (0.28-0.69 m^2), median",
      "haemoglobin 10.7 g/L (5.7-29.7 g/L; the upper value is most likely",
      "a typo or unit ambiguity in the source PDF -- see vignette",
      "Assumptions and deviations), median albumin 38 g/L (29-47 g/L)."
    )
  )

  ini({
    # =====================================================================
    # Lopinavir (LPV, parent / substrate) structural parameters
    # Zhang 2012 Table 2 'Population pharmacokinetic parameter estimates'
    # column 'Final model estimates'. Values reported at the cohort median
    # 10 kg per Methods Equations 1 and 2.
    # =====================================================================
    lcl    <- log(4.18)
    label("LPV typical apparent clearance CL0/F at the 10 kg reference, without ritonavir (L/h)")  # Zhang 2012 Table 2 row 'Lopinavir CL/F = 4.18'; note Results 'Model description' paragraph 3: 'the typical clearance of LPV without ritonavir was 4.18 l/h, it should be kept in mind that this value is an extrapolation, since LPV was never given without ritonavir'
    lvc    <- log(11.6)
    label("LPV typical apparent volume of distribution V/F at the 10 kg reference (L)")            # Zhang 2012 Table 2 row 'Lopinavir V/F = 11.6'
    lka    <- log(0.74)
    label("LPV first-order absorption rate constant ka (1/h)")                                     # Zhang 2012 Table 2 row 'Lopinavir ka = 0.74'
    lfdepot <- fixed(log(1))
    label("LPV bioavailability anchor at the standard 4:1 LPV/r no-rifampicin reference (unitless, FIXED)")  # Standard 4:1 LPV/r at the median 3 mg/kg ritonavir without rifampicin co-administration is the F = 1 reference per Zhang 2012 Methods 'Population pharmacokinetic analysis' paragraph 5: 'The relative bioavailability in the control group (standard LPV/r dose, no rifampicin) was assumed as a reference (100%)'.

    # =====================================================================
    # Ritonavir (RTV, sibling-drug / perpetrator) structural parameters
    # Zhang 2012 Table 2 'Population pharmacokinetic parameter estimates'
    # column 'Final model estimates'. Values reported at the cohort median
    # 10 kg per Methods Equations 1 and 2. RTV CL/F differs by rifampicin
    # co-administration status (no-TB / after-TB vs with-TB).
    # =====================================================================
    lcl_rtv      <- log(12.8)
    label("RTV typical apparent clearance CL/F at the 10 kg reference, no rifampicin co-administration (L/h)")  # Zhang 2012 Table 2 row 'Ritonavir CL/F, No TB and after TB = 12.8'
    lcl_rtv_rif  <- log(19.1)
    label("RTV typical apparent clearance CL/F at the 10 kg reference, with rifampicin co-administration (L/h)")  # Zhang 2012 Table 2 row 'Ritonavir CL/F, With TB = 19.1'
    lvc_rtv      <- log(105)
    label("RTV typical apparent volume of distribution V/F at the 10 kg reference (L)")            # Zhang 2012 Table 2 row 'Ritonavir V/F = 105'
    lka_rtv      <- log(2.31)
    label("RTV first-order absorption rate constant ka from the final transit compartment to central (1/h)")    # Zhang 2012 Table 2 row 'Ritonavir ka = 2.31'
    lmtt_rtv     <- log(1.28)
    label("RTV mean transit time MTT through the 10-compartment Savic transit chain (h)")          # Zhang 2012 Table 2 row 'Ritonavir MTT = 1.28'
    nn_rtv       <- fixed(10)
    label("RTV number of Savic-style transit compartments (integer, unitless, FIXED)")             # Zhang 2012 Results 'Model description' paragraph 1: 'the absorption phase displayed more complex pharmacokinetics which was described best by a series of 10 transit compartments'; the paper does not state whether the chain length was estimated or fixed before the final fit, but Table 2 does not list NN as an estimated parameter; treated as fixed here

    # =====================================================================
    # Allometric exponents. The paper Methods 'Population pharmacokinetic
    # analysis' paragraph 4 / Equations 1 and 2 applies allometric scaling
    # at the 10 kg median body weight to apparent clearance and volume of
    # distribution of both drugs. The exponents are not separately listed
    # with bootstrap CIs in Table 2; they are the canonical 0.75 / 1
    # carried from Holford 1996 (Zhang 2012 references 10, 11), and treated
    # as fixed structural choices here.
    # =====================================================================
    e_wt_cl <- fixed(0.75)
    label("Allometric exponent on apparent CL/F for both LPV and RTV (unitless)")                  # Zhang 2012 Methods Equations 1 and 2, refs 10 and 11 (Holford / Anderson-Holford allometry)
    e_wt_vc <- fixed(1.0)
    label("Allometric exponent on apparent V/F for both LPV and RTV (unitless)")                   # Zhang 2012 Methods Equations 1 and 2, refs 10 and 11

    # =====================================================================
    # Lopinavir bioavailability covariate model (Zhang 2012 Equation 4):
    #   F_LPV = 1 - RIF * CONMED_RIF + SLP * (DOSE_RTV_MGKG - DOSE_RTV_STD)
    # with DOSE_RTV_STD = 3 mg/kg (median RTV dose in the standard-arm
    # reference). The reference is the standard 4:1 LPV/r at the median
    # RTV dose with no rifampicin: F_LPV = 1. Rifampicin alone (no extra
    # ritonavir) would give F_LPV = 1 - 0.832 = 0.168 (Discussion
    # paragraph 3: 'exposure to LPV would drop by 83.2% if antitubercular
    # treatment was concomitantly given without any further dose
    # adjustments').
    # =====================================================================
    e_rif_lpv <- 0.832
    label("Multiplicative -83.2% reduction in LPV relative bioavailability under rifampicin co-administration at the median 3 mg/kg ritonavir reference (unitless)")  # Zhang 2012 Table 2 row 'Lopinavir RIF on F = 0.832'; equation 4 coefficient 'RIF'
    e_dose_rtv_lpv <- 0.021
    label("Linear slope of LPV relative bioavailability on per-administration ritonavir dose in mg/kg (1/(mg/kg))")  # Zhang 2012 Table 2 row 'Lopinavir Slope = 0.021'; equation 4 coefficient 'SLP'

    # =====================================================================
    # LPV-RTV interaction: direct sigmoid-Emax inhibition of LPV apparent
    # clearance by ritonavir plasma concentration (Zhang 2012 Equation 3):
    #   CL_LPV(t) = CL0_LPV * (1 - Emax * C_RTV / (EC50 + C_RTV))
    # Emax was fixed at 0.9 because of numerical instability when all
    # interaction parameters were simultaneously estimated (Results
    # 'Model description' paragraph 4: 'the Emax was fixed to 0.9. This
    # value was estimated when ritonavir parameters were fixed and only
    # LPV parameters estimated').
    # =====================================================================
    emax_lpv <- fixed(0.9)
    label("Maximum fractional inhibition of LPV CL/F by ritonavir concentration (unitless, FIXED)")  # Zhang 2012 Table 2 row 'Lopinavir-ritonavir interaction Emax = 0.9 (fix)'
    ec50_lpv <- 0.0519
    label("Ritonavir plasma concentration producing 50% of Emax on LPV CL/F (mg/L)")                # Zhang 2012 Table 2 row 'Lopinavir-ritonavir interaction EC50 = 0.0519 mg/l'

    # =====================================================================
    # Inter-individual variability (Zhang 2012 Table 2 BSV/IOV rows).
    # IIV reported on the SD scale as approximate %CV; nlmixr2lib uses
    # variances on the log scale via omega^2 = log(1 + CV^2).
    #
    # LPV: only BSV on V (56.6%); IOV ka (76.2%) and IOV F (51.8%) are
    # folded in as BSV-equivalent following the Bienczak 2016 nevirapine
    # precedent (no separate BSV reported on the same parameters).
    # =====================================================================
    etalvc      ~ 0.27406   # Zhang 2012 Table 2 IIV V LPV = 56.6%; omega^2 = log(1 + 0.566^2) = 0.27406
    etalka      ~ 0.45078   # Zhang 2012 Table 2 IOV ka LPV = 76.2% folded in as BSV-equivalent; omega^2 = log(1 + 0.762^2) = 0.45078 (no separate BSV on LPV ka reported)
    etalfdepot  ~ 0.23475   # Zhang 2012 Table 2 IOV F LPV  = 51.8% folded in as BSV-equivalent; omega^2 = log(1 + 0.518^2) = 0.23475 (no separate BSV on LPV F reported)

    # =====================================================================
    # RTV: paper Results 'Model description' paragraph 1 reports 'A similar
    # solution was used for IIV in oral clearance and volume of
    # distribution of ritonavir' as for the ka_LPV / ka_RTV pair, i.e. the
    # two IIVs were estimated as proportional (high positive correlation).
    # Encoded here as a correlated IIV block with rho = 0.99 approximating
    # the perfect-proportionality assumption (see vignette Assumptions and
    # deviations). IOV CL RTV (41.6%) dropped because BSV is reported on
    # the same parameter (Bienczak 2016 precedent). IOV MTT and IOV ka of
    # RTV folded in as BSV-equivalent.
    # =====================================================================
    # Block entries below come from Zhang 2012 Table 2:
    #   var_cl_rtv = log(1 + 0.728^2) = 0.41888 (IIV CL RTV = 72.8%)
    #   var_vc_rtv = log(1 + 0.433^2) = 0.16758 (IIV V  RTV = 43.3%)
    #   cov_cl_vc  = rho * sqrt(var_cl_rtv * var_vc_rtv) with rho = 0.99
    #   approximating the paper's 'proportional' description of IIV CL_RTV
    #   and V_RTV (Results 'Model description' paragraph 1).
    etalcl_rtv + etalvc_rtv ~ c(
      0.41888,
      0.99 * sqrt(0.41888 * 0.16758),
      0.16758
    )
    etalmtt_rtv ~ 0.09128   # Zhang 2012 Table 2 IOV MTT RTV = 31.1% folded in as BSV-equivalent; omega^2 = log(1 + 0.311^2) = 0.09128 (no separate BSV on RTV MTT reported)
    etalka_rtv  ~ 0.65979   # Zhang 2012 Table 2 IOV ka  RTV = 98.1% folded in as BSV-equivalent; omega^2 = log(1 + 0.981^2) = 0.65979 (no separate BSV on RTV ka reported)

    # =====================================================================
    # Residual error. Zhang 2012 Methods 'Population pharmacokinetic
    # analysis' paragraph 3: 'The exponential error model was implemented
    # in NONMEM by log-transforming the data. Hence, using a first-order
    # approximation, the variability of the exponential model can be
    # considered as proportional to the observed value.' Encoded as
    # proportional on the linear scale per the paper's interpretation.
    # =====================================================================
    propSd     <- 0.304
    label("LPV proportional residual error on the linear scale (fraction)")                        # Zhang 2012 Table 2 row 'Lopinavir RUV = 0.304' (exponential error on log-data, equivalent to proportional on the linear scale)
    propSd_rtv <- 0.339
    label("RTV proportional residual error on the linear scale (fraction)")                        # Zhang 2012 Table 2 row 'Ritonavir RUV = 0.339' (exponential error on log-data, equivalent to proportional on the linear scale)
  })

  model({
    # ----------------------------------------------------------------------
    # 1. Individual ritonavir PK parameters (compute first so the RTV
    #    plasma concentration is available for the LPV CL inhibition term).
    # ----------------------------------------------------------------------
    # RTV typical CL/F depends on rifampicin status: the paper Table 2
    # reports two separate typical values (12.8 L/h no rifampicin, 19.1 L/h
    # with rifampicin). Encoded as a binary switch on the typical value
    # rather than a multiplicative coefficient on a baseline because the
    # paper does so (Results 'Model description' paragraph 2: 'different
    # typical values of clearance were estimated for the subjects with
    # and without antitubercular treatment'). Allometric scaling at 10 kg.
    cl_rtv_typ <- exp(lcl_rtv     + etalcl_rtv) * (1 - CONMED_RIF) +
                  exp(lcl_rtv_rif + etalcl_rtv) * CONMED_RIF
    cl_rtv     <- cl_rtv_typ * (WT / 10)^e_wt_cl
    vc_rtv     <- exp(lvc_rtv + etalvc_rtv) * (WT / 10)^e_wt_vc
    ka_rtv     <- exp(lka_rtv + etalka_rtv)
    mtt_rtv    <- exp(lmtt_rtv + etalmtt_rtv)
    ktr_rtv    <- nn_rtv / mtt_rtv      # transit-chain inter-compartment rate; Bienczak 2016 nevirapine precedent (depot -> transit_1 -> ... -> transit_NN with rate ktr; transit_NN -> central with rate ka)

    # ----------------------------------------------------------------------
    # 2. Sigmoid-Emax inhibition of LPV CL/F by ritonavir plasma
    #    concentration (Zhang 2012 Equation 3). When no ritonavir is
    #    present (Crtv = 0) the inhibition evaluates to zero and CL_LPV
    #    equals exp(lcl) = 4.18 L/h (the extrapolated no-RTV typical value).
    # ----------------------------------------------------------------------
    crtv  <- central_rtv / vc_rtv
    inhib <- emax_lpv * crtv / (ec50_lpv + crtv)

    # ----------------------------------------------------------------------
    # 3. Individual lopinavir PK parameters (CL/F inhibited by the
    #    ritonavir-driven sigmoid; allometric scaling at 10 kg).
    # ----------------------------------------------------------------------
    cl <- exp(lcl) * (WT / 10)^e_wt_cl * (1 - inhib)           # Zhang 2012 Table 2 reports no BSV on LPV CL/F (no etalcl term)
    vc <- exp(lvc + etalvc) * (WT / 10)^e_wt_vc
    ka <- exp(lka + etalka)

    # ----------------------------------------------------------------------
    # 4. Lopinavir bioavailability (Zhang 2012 Equation 4). Linear shift
    #    in DOSE_RTV_MGKG above the 3 mg/kg reference plus a -83.2%
    #    multiplicative reduction under rifampicin co-administration.
    #    Random effect on F (etalfdepot) carries the IOV-equivalent BSV
    #    folded in above; the random effect is on the log-multiplier of
    #    the typical bioavailability so it preserves the bounded-on-[0,1]
    #    behaviour reasonably in the simulated range.
    # ----------------------------------------------------------------------
    f_lpv_typ <- exp(lfdepot) *
                 (1 - e_rif_lpv * CONMED_RIF +
                      e_dose_rtv_lpv * (DOSE_RTV_MGKG - 3))
    f_lpv     <- f_lpv_typ * exp(etalfdepot)

    # ----------------------------------------------------------------------
    # 5. ODE system. Lopinavir: depot -> central with first-order ka and
    #    first-order systemic elimination CL/Vc. Ritonavir: depot_rtv ->
    #    transit1_rtv -> ... -> transit10_rtv -> central_rtv with shared
    #    transit-chain rate ktr_rtv on the depot-to-transit and inter-
    #    transit steps and a separate ka_rtv from the last transit into
    #    central_rtv (Bienczak 2016 nevirapine precedent for the
    #    Savic-style transit chain with a separate final ka). Apparent
    #    clearance CL/Vc for both drugs.
    # ----------------------------------------------------------------------
    d/dt(depot)         <- -ka * depot
    d/dt(central)       <-  ka * depot - (cl / vc) * central

    d/dt(depot_rtv)     <- -ktr_rtv * depot_rtv
    d/dt(transit1_rtv)  <-  ktr_rtv * depot_rtv     - ktr_rtv * transit1_rtv
    d/dt(transit2_rtv)  <-  ktr_rtv * transit1_rtv  - ktr_rtv * transit2_rtv
    d/dt(transit3_rtv)  <-  ktr_rtv * transit2_rtv  - ktr_rtv * transit3_rtv
    d/dt(transit4_rtv)  <-  ktr_rtv * transit3_rtv  - ktr_rtv * transit4_rtv
    d/dt(transit5_rtv)  <-  ktr_rtv * transit4_rtv  - ktr_rtv * transit5_rtv
    d/dt(transit6_rtv)  <-  ktr_rtv * transit5_rtv  - ktr_rtv * transit6_rtv
    d/dt(transit7_rtv)  <-  ktr_rtv * transit6_rtv  - ktr_rtv * transit7_rtv
    d/dt(transit8_rtv)  <-  ktr_rtv * transit7_rtv  - ktr_rtv * transit8_rtv
    d/dt(transit9_rtv)  <-  ktr_rtv * transit8_rtv  - ktr_rtv * transit9_rtv
    d/dt(transit10_rtv) <-  ktr_rtv * transit9_rtv  - ka_rtv  * transit10_rtv
    d/dt(central_rtv)   <-  ka_rtv  * transit10_rtv - (cl_rtv / vc_rtv) * central_rtv

    # ----------------------------------------------------------------------
    # 6. Bioavailability on the LPV depot. RTV bioavailability is the
    #    apparent F embedded in the CL/F and V/F estimates (the paper
    #    fits CL/F and V/F, not separate CL and V).
    # ----------------------------------------------------------------------
    f(depot) <- f_lpv

    # ----------------------------------------------------------------------
    # 7. Observation variables and proportional residual error.
    #    Cc     = LPV plasma concentration (mg/L)
    #    Cc_rtv = RTV plasma concentration (mg/L)
    # ----------------------------------------------------------------------
    Cc     <- central     / vc
    Cc_rtv <- central_rtv / vc_rtv

    Cc     ~ prop(propSd)
    Cc_rtv ~ prop(propSd_rtv)
  })
}
