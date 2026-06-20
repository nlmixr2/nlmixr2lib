Ramharter_2019_mefloquine <- function() {
  description <- paste(
    "Population PK model for the two mefloquine enantiomers and their",
    "carboxy-metabolite carboxymefloquine (CMQ) in pregnant African women",
    "receiving intermittent preventive treatment for malaria (Ramharter",
    "2019, MIPPAD trial, Gabon). Each parent enantiomer ((+)-mefloquine",
    "= '_r' = (11R, 2'S); (-)-mefloquine = '_s' = (11S, 2'R)) follows a",
    "two-compartment disposition with first-order oral absorption.",
    "Carboxymefloquine ('_cmq') is formed molar 1:1 from both parents via",
    "the apparent parent clearance and follows a two-compartment",
    "disposition; its first-order clearance is autoinduced by CMQ plasma",
    "concentration via a two-stage RNA + enzyme-pool turnover model",
    "(precursor1 = enzymatic-RNA precursor pool, precursor2 = metabolizing",
    "enzyme pool; both kdeg-driven, both at unit steady state in the",
    "absence of CMQ; CMQ clearance scales linearly with precursor2). A",
    "shared body-weight allometric exponent acts on the central volume of",
    "each parent enantiomer (reference 55 kg). The split-dose IPTp",
    "regimen (REGIMEN_SPLIT = 1: 7.5 mg/kg on two consecutive days) carries",
    "a small (+5%) bioavailability increment relative to the single-dose",
    "regimen (REGIMEN_SPLIT = 0: 15 mg/kg on day 1).",
    sep = " "
  )
  reference <- paste(
    "Ramharter M, Schwab M, Mombo-Ngoma G, Zoleko Manego R, Akerey-Diop D,",
    "Basra A, Mackanga J-R, Wurbel H, Wojtyniak J-G, Gonzalez R, Hofmann U,",
    "Geditz M, Matsiegui P-B, Kremsner PG, Menendez C, Kerb R, Lehr T. 2019.",
    "Population pharmacokinetics of mefloquine intermittent preventive",
    "treatment for malaria in pregnancy in Gabon.",
    "Antimicrob Agents Chemother 63:e01113-18.",
    "doi:10.1128/AAC.01113-18.",
    sep = " "
  )
  vignette <- "Ramharter_2019_mefloquine"
  units <- list(time = "hour", dosing = "nmol", concentration = "nmol/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight at baseline.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Drives allometric scaling of the central volume of distribution of",
        "BOTH parent enantiomers via a shared estimated exponent",
        "WTexponent = 1.33 centred on the cohort median 55 kg (Ramharter",
        "2019 Table 2 and Results 'Mefloquine PK model' paragraph). CMQ",
        "volumes and all clearances were not weight-scaled in the source",
        "covariate analysis. Baseline weight only; a more complex model",
        "incorporating gestational weight change did not improve fit (paper",
        "Results 'Mefloquine PK model')."
      ),
      source_name        = "WT"
    ),
    REGIMEN_SPLIT = list(
      description        = paste(
        "Binary indicator for the split-dose mefloquine IPTp regimen:",
        "1 = 7.5 mg/kg on two consecutive days (split-dose group);",
        "0 = 15 mg/kg on a single day (single-dose group)."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (single-dose 15 mg/kg on day 1)",
      notes              = paste(
        "Acts multiplicatively on the relative oral bioavailability of",
        "both parent enantiomer depots (paper Table 2: SPLIT = 1.05; 90% CI",
        "1.01-1.16). Encoded in the model as a +5% bioavailability",
        "increment when REGIMEN_SPLIT = 1 (i.e. e_regimen_split_fdepot =",
        "0.05). The paper attributes the effect to saturation of",
        "intestinal / hepatic uptake processes after the single high dose",
        "and concludes the magnitude 'may have only limited clinical",
        "significance' (Discussion paragraph 3). Time-fixed per subject."
      ),
      source_name        = "SPLIT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 263L,
    n_studies      = 1L,
    age_range      = "14-44 years (mean 24.0, SD 6.75, median 22.0)",
    weight_range   = "39.3-108 kg (mean 58.3, SD 10.8, median 55.9)",
    sex_female_pct = 100,
    race_ethnicity = "Sub-Saharan African (Gabonese; race / ethnicity composition not separately tabulated in source)",
    disease_state  = paste(
      "Asymptomatic, presumed-aparasitaemic pregnant women in malaria-endemic",
      "settings, attending the antenatal clinic for the first time between",
      "13 and 28 weeks of gestation (mean gestational age 17.8 weeks, SD 5.85,",
      "median 19.0, range 8-28). All HIV-negative (HIV infection was an",
      "exclusion criterion). Body mass index mean 23.2 kg/m^2 (SD 3.91, range",
      "15.6-41.7). Haemoglobin mean 10.2 g/dL (SD 1.38, range 5.7-14.7)."
    ),
    dose_range     = paste(
      "Mefloquine racemate 15 mg/kg total administered orally, split into one",
      "of two regimens per the MIPPAD trial design: single-dose (15 mg/kg on",
      "day 1; REGIMEN_SPLIT = 0) or split-dose (7.5 mg/kg on day 1 and 7.5",
      "mg/kg on day 2; REGIMEN_SPLIT = 1). Two IPTp administrations per",
      "subject at least 1 month apart during the second / third trimester.",
      "129 women allocated to the single-dose regimen, 134 to the split-dose",
      "regimen. 26 women attended only the first dosing interval. The user",
      "supplies two dose records per administration (one per enantiomer",
      "depot: cmt = depot_r and cmt = depot_s) each with amt set to half the",
      "racemate molar dose in nmol -- see vignette for the mg -> nmol",
      "conversion at the parent free-base molecular weight 378.31 g/mol."
    ),
    regions        = "Gabon (Lambarene and Fougamou study centres, MIPPAD trial NCT0081121)",
    notes          = paste(
      "Demographics from Ramharter 2019 Table 1 (n = 263; the abstract's n =",
      "264 is a typographical inconsistency in the source -- Table 1, Results",
      "paragraph 1, and the supplemental material all give n = 263, with the",
      "129 + 134 = 263 dose-arm split). Rich-PK sampling subgroup n = 37",
      "(5-14 samples per subject, median 11); sparse-PK subgroup n = 226",
      "(median 3.5 samples per subject). 924 (-)-MQ + 923 (+)-MQ + 924 CMQ",
      "plasma observations contributed to the fit. Cord-blood samples (n =",
      "126) were also fitted but with the cord-to-plasma factor (CORD)",
      "fixed at 1 and a separate set of residual-error magnitudes per",
      "matrix (paper Results 'Transplacental distribution'); the model file",
      "encodes only the plasma observation arm with its three plasma",
      "residual-error pairs (cord residual-error magnitudes are tabulated",
      "in the vignette Assumptions and deviations section, where the choice",
      "to drop cord observations is justified)."
    )
  )

  ini({
    # =========================================================================
    # (+)-mefloquine structural parameters (paper Table 2, "Fixed effects").
    # Suffix `_r` follows the registered R-enantiomer convention (Valitalo
    # 2017 ketorolac precedent); (+)-mefloquine has absolute configuration
    # (11R, 2'S) per Wong et al. (2017) Cell 169:73-85 -- the R designation
    # is keyed to the C-11 stereocentre. Apparent-F-relative volumes and
    # clearances (paper Table 2 column 'Description'): all rates per hour,
    # all volumes in L, concentrations in nmol/L.
    # =========================================================================
    lka_r  <- log(0.209)  ; label("(+)-MQ first-order absorption rate constant Ka (1/h)")            # Ramharter 2019 Table 2: Ka_(+)MQ = 0.209/h (90% CI 0.032-0.329)
    lvc_r  <- log(814)    ; label("(+)-MQ apparent central volume of distribution Vc/F (L)")         # Ramharter 2019 Table 2: Vcentral_(+)MQ/F = 814 L (90% CI 252-1230)
    lq_r   <- log(141)    ; label("(+)-MQ apparent inter-compartmental clearance Q/F (L/h)")          # Ramharter 2019 Table 2: Q_(+)MQ/F = 141 L/h (90% CI 22.2-182)
    lvp_r  <- log(1070)   ; label("(+)-MQ apparent peripheral volume of distribution Vp/F (L)")       # Ramharter 2019 Table 2: Vperipheral_(+)MQ/F = 1070 L (90% CI 808-1492)
    lcl_r  <- log(8.28)   ; label("(+)-MQ apparent metabolic clearance CL/F to CMQ (L/h)")            # Ramharter 2019 Table 2: CL_(+)MQ/F = 8.28 L/h (90% CI 7.55-9.28)

    # =========================================================================
    # (-)-mefloquine structural parameters (paper Table 2, "Fixed effects").
    # Suffix `_s` follows the registered S-enantiomer convention;
    # (-)-mefloquine has absolute configuration (11S, 2'R) per Wong et al.
    # (2017) Cell 169:73-85.
    # =========================================================================
    lka_s  <- log(0.157)  ; label("(-)-MQ first-order absorption rate constant Ka (1/h)")            # Ramharter 2019 Table 2: Ka_(-)MQ = 0.157/h (90% CI 0.023-0.303)
    lvc_s  <- log(476)    ; label("(-)-MQ apparent central volume of distribution Vc/F (L)")         # Ramharter 2019 Table 2: Vcentral_(-)MQ/F = 476 L (90% CI 132-801)
    lq_s   <- log(85.7)   ; label("(-)-MQ apparent inter-compartmental clearance Q/F (L/h)")          # Ramharter 2019 Table 2: Q_(-)MQ/F = 85.7 L/h (90% CI 12.8-120.5)
    lvp_s  <- log(860)    ; label("(-)-MQ apparent peripheral volume of distribution Vp/F (L)")       # Ramharter 2019 Table 2: Vperipheral_(-)MQ/F = 860 L (90% CI 657-1173)
    lcl_s  <- log(1.49)   ; label("(-)-MQ apparent metabolic clearance CL/F to CMQ (L/h)")            # Ramharter 2019 Table 2: CL_(-)MQ/F = 1.49 L/h (90% CI 1.38-1.60)

    # =========================================================================
    # Carboxymefloquine (CMQ) structural parameters (paper Table 2, "Fixed
    # effects"). CMQ is formed molar 1:1 from each parent enantiomer via the
    # parent CL/F and follows a two-compartment disposition with autoinduced
    # first-order clearance (see autoinduction parameters below).
    # =========================================================================
    lvc_cmq <- log(47.8)   ; label("CMQ apparent central volume of distribution Vc/F (L)")            # Ramharter 2019 Table 2: Vcentral_CMQ/F = 47.8 L (90% CI 7.6-62.0)
    lq_cmq  <- log(16.8)   ; label("CMQ apparent inter-compartmental clearance Q/F (L/h)")             # Ramharter 2019 Table 2: Q_CMQ/F = 16.8 L/h (90% CI 11.9-20.1)
    lvp_cmq <- log(968)    ; label("CMQ apparent peripheral volume of distribution Vp/F (L)")          # Ramharter 2019 Table 2: Vperipheral_CMQ/F = 968 L (90% CI 657-1063)
    lcl_cmq <- log(1.12)   ; label("CMQ baseline apparent clearance CL/F (L/h; scaled at run time by precursor2 = enzyme-pool amount)")  # Ramharter 2019 Table 2: CL_CMQ/F = 1.12 L/h (90% CI 0.78-2.58)

    # =========================================================================
    # CMQ autoinduction parameters. Two-stage turnover model (Ramharter 2019
    # Results "CMQ PK model"):
    #   precursor1 (enzymatic-RNA precursor pool) has zero-order synthesis
    #     induced by CMQ via a proportional Emax model; the maximum induction
    #     was estimated as a 3.38-fold INCREASE of the synthesis rate, so the
    #     standard (1 + Emax * Cc/(EC50 + Cc)) parameterisation gives a peak
    #     factor of 4.38 (paper text: "maximum induction was estimated as a
    #     3.38-fold increase of the synthesis rate").
    #   precursor2 (metabolizing enzyme pool) has zero-order synthesis
    #     proportional to the amount in precursor1; the half-life of the
    #     enzyme pool was estimated at 152 h (paper text), corresponding to
    #     kdeg = ln(2)/152 = 0.00456/h (paper Table 2 rounds to 0.00453/h).
    #   Both pools share kdeg (only one Kdeg value reported in Table 2);
    #   both pools are initialised at the steady-state value 1 (paper
    #   normalisation, so that CMQ clearance equals CL_CMQ * 1 at Cc_cmq = 0).
    #   EC50 is FIX (paper Table 2 column "90% CI" reads "FIX") because the
    #   data did not identify it.
    # =========================================================================
    lemax    <- log(3.38)              ; label("Maximum induction factor on the precursor1 synthesis rate (unitless multiplier; peak factor is 1 + Emax = 4.38)")  # Ramharter 2019 Table 2: Emax = 3.38 (90% CI 0.54-5.07)
    lec50    <- fixed(log(2.45))       ; label("CMQ plasma concentration achieving half-maximum induction effect EC50 (nmol/L); FIX in source")                     # Ramharter 2019 Table 2: EC50 = 2.45 nmol/L FIX
    lkdeg    <- log(0.00453)           ; label("First-order degradation rate constant of precursor1 and precursor2 enzyme-turnover pools (1/h); pool half-life 152 h")  # Ramharter 2019 Table 2: Kdeg = 0.00453/h (90% CI 0.00144-0.6433)

    # =========================================================================
    # Covariate effects (paper Table 2). The body-weight allometric exponent
    # is SHARED between Vc_(+)MQ and Vc_(-)MQ (paper text: "modelled
    # allometrically centered around the median body weight of 55 kg with an
    # estimated exponent of 1.33, where an increased body weight resulted in
    # increased volumes of distribution"); applied identically to both
    # central volumes via cov_wt_vc in model(). The split-dose effect is
    # multiplicative on the relative bioavailability of both parent depots.
    # =========================================================================
    e_wt_vc                  <- 1.33   ; label("Body-weight allometric exponent shared across the central volumes of (+)-MQ and (-)-MQ (unitless); reference 55 kg")  # Ramharter 2019 Table 2: WTexponent = 1.33 (90% CI 0.42-2.48)
    e_regimen_split_fdepot   <- 0.05   ; label("Proportional bioavailability increase under the split-dose IPTp regimen (unitless; +5% when REGIMEN_SPLIT = 1)")        # Ramharter 2019 Table 2: SPLIT = 1.05 (90% CI 1.01-1.16), encoded as (1 + 0.05 * REGIMEN_SPLIT)

    # =========================================================================
    # IIV. The paper estimated ONE eta on Vc shared across (+)-MQ and (-)-MQ
    # (paper Table 2: "IIV Vcentral,MQ"), ONE eta on CL shared across (+)-MQ
    # and (-)-MQ (paper Table 2: "IIV CLMQ"), and ONE eta on CL_CMQ. Reported
    # as %CV: omega^2 = log(1 + CV^2) for a log-normal random effect.
    #   Vc shared: %CV = 119  -> omega^2 = log(1 + 1.19^2) = log(2.4161) = 0.8821
    #   CL shared: %CV = 34.6 -> omega^2 = log(1 + 0.346^2) = log(1.1197) = 0.1132
    #   CL_CMQ:    %CV = 42.7 -> omega^2 = log(1 + 0.427^2) = log(1.1823) = 0.1675
    # The shared-across-enantiomers eta is implemented in nlmixr2 as a 2x2
    # block with var_r = var_s = cov = omega^2 (perfect correlation), so the
    # two registered etas (etalvc_r / etalvc_s, etalcl_r / etalcl_s) yield
    # identical realised values for each subject and the structure matches
    # the paper's single shared eta exactly.
    # =========================================================================
    etalvc_r + etalvc_s ~ c(0.8821, 0.8821, 0.8821)   # Ramharter 2019 Table 2: IIV Vcentral,MQ = 119 %CV (90% CI 80.3-189) -- shared across (+)-MQ and (-)-MQ
    etalcl_r + etalcl_s ~ c(0.1132, 0.1132, 0.1132)   # Ramharter 2019 Table 2: IIV CL,MQ = 34.6 %CV (90% CI 33.7-47.3) -- shared across (+)-MQ and (-)-MQ
    etalcl_cmq          ~ 0.1675                       # Ramharter 2019 Table 2: IIV CL,CMQ = 42.7 %CV (90% CI 39.4-55.6) -- independent CMQ-clearance eta

    # =========================================================================
    # Residual error (paper Table 2, "Random effects: residual variability").
    # Plasma-arm magnitudes only; cord-arm residual magnitudes are documented
    # in the vignette Assumptions and deviations section (CORD = 1 fixed, so
    # the cord prediction equals the plasma prediction with a different
    # residual error -- collapsing to the plasma arm omits the cord arm).
    # =========================================================================
    propSd_r   <- 0.385   ; label("(+)-MQ plasma proportional residual SD (fraction)")             # Ramharter 2019 Table 2: PRV (+)MQ plasma = 38.5% (90% CI 35.7-45.4)
    addSd_r    <- 12      ; label("(+)-MQ plasma additive residual SD (nmol/L)")                    # Ramharter 2019 Table 2: ARV (+)MQ plasma = +-12 nmol/L (90% CI 8.4-17.5)
    propSd_s   <- 0.299   ; label("(-)-MQ plasma proportional residual SD (fraction)")             # Ramharter 2019 Table 2: PRV (-)MQ plasma = 29.9% (90% CI 27.6-34.8)
    addSd_s    <- 53      ; label("(-)-MQ plasma additive residual SD (nmol/L)")                    # Ramharter 2019 Table 2: ARV (-)MQ plasma = +-53 nmol/L (90% CI 13.4-69.9)
    propSd_cmq <- 0.248   ; label("CMQ plasma proportional residual SD (fraction)")                # Ramharter 2019 Table 2: PRV CMQ plasma = 24.8% (90% CI 21.5-35.6)
    addSd_cmq  <- 63      ; label("CMQ plasma additive residual SD (nmol/L)")                       # Ramharter 2019 Table 2: ARV CMQ plasma = +-63 nmol/L (90% CI 8.6-94.8)
  })

  model({
    # -----------------------------------------------------------------------
    # 1. Reference body weight (paper Methods: cohort median 55 kg) and
    #    covariate factors. The body-weight allometric scaling acts on the
    #    central volume of BOTH parent enantiomers via the shared exponent
    #    e_wt_vc; the split-dose effect is a binary +5% bioavailability
    #    increment applied to BOTH parent depots when REGIMEN_SPLIT = 1.
    # -----------------------------------------------------------------------
    wt_ref       <- 55
    cov_wt_vc    <- (WT / wt_ref) ^ e_wt_vc
    cov_split_f  <- 1 + e_regimen_split_fdepot * REGIMEN_SPLIT

    # -----------------------------------------------------------------------
    # 2. Individual PK parameters. Shared etas across enantiomers are
    #    implemented via the perfect-correlation block in ini() -- etalvc_r
    #    and etalvc_s carry identical realisations, mirroring the paper's
    #    single shared eta on Vc_MQ; same for etalcl_r and etalcl_s on
    #    CL_MQ. etalcl_cmq is independent.
    # -----------------------------------------------------------------------
    ka_r   <- exp(lka_r)
    vc_r   <- exp(lvc_r + etalvc_r) * cov_wt_vc
    q_r    <- exp(lq_r)
    vp_r   <- exp(lvp_r)
    cl_r   <- exp(lcl_r + etalcl_r)

    ka_s   <- exp(lka_s)
    vc_s   <- exp(lvc_s + etalvc_s) * cov_wt_vc
    q_s    <- exp(lq_s)
    vp_s   <- exp(lvp_s)
    cl_s   <- exp(lcl_s + etalcl_s)

    vc_cmq <- exp(lvc_cmq)
    q_cmq  <- exp(lq_cmq)
    vp_cmq <- exp(lvp_cmq)
    cl_cmq <- exp(lcl_cmq + etalcl_cmq)

    emax   <- exp(lemax)
    ec50   <- exp(lec50)
    kdeg   <- exp(lkdeg)

    # -----------------------------------------------------------------------
    # 3. Plasma concentrations (observation variables; nmol/L).
    # -----------------------------------------------------------------------
    Cc_r   <- central_r   / vc_r
    Cc_s   <- central_s   / vc_s
    Cc_cmq <- central_cmq / vc_cmq

    # -----------------------------------------------------------------------
    # 4. (+)-mefloquine ODEs. The parent's elimination from central is the
    #    metabolic flux into CMQ central (cl_r * Cc_r); the same molar mass
    #    is added to d/dt(central_cmq) below (molar 1:1 stoichiometry; paper
    #    Methods: parent metabolised "by a first-order clearance process
    #    into CMQ"). Inter-compartmental transfer uses q_r/vc_r and q_r/vp_r.
    # -----------------------------------------------------------------------
    d/dt(depot_r)        <- -ka_r * depot_r
    d/dt(central_r)      <-  ka_r * depot_r -
                              cl_r * Cc_r -
                              (q_r / vc_r) * central_r +
                              (q_r / vp_r) * peripheral1_r
    d/dt(peripheral1_r)  <-  (q_r / vc_r) * central_r -
                              (q_r / vp_r) * peripheral1_r

    # -----------------------------------------------------------------------
    # 5. (-)-mefloquine ODEs. Same structural form as (+)-MQ above.
    # -----------------------------------------------------------------------
    d/dt(depot_s)        <- -ka_s * depot_s
    d/dt(central_s)      <-  ka_s * depot_s -
                              cl_s * Cc_s -
                              (q_s / vc_s) * central_s +
                              (q_s / vp_s) * peripheral1_s
    d/dt(peripheral1_s)  <-  (q_s / vc_s) * central_s -
                              (q_s / vp_s) * peripheral1_s

    # -----------------------------------------------------------------------
    # 6. Carboxymefloquine (CMQ) ODEs. Inflow is the molar 1:1 sum of the
    #    metabolic fluxes from each parent; outflow is the autoinduced
    #    first-order clearance cl_cmq * precursor2 multiplying Cc_cmq.
    # -----------------------------------------------------------------------
    d/dt(central_cmq)    <-  cl_r * Cc_r + cl_s * Cc_s -
                              cl_cmq * precursor2 * Cc_cmq -
                              (q_cmq / vc_cmq) * central_cmq +
                              (q_cmq / vp_cmq) * peripheral1_cmq
    d/dt(peripheral1_cmq) <- (q_cmq / vc_cmq) * central_cmq -
                              (q_cmq / vp_cmq) * peripheral1_cmq

    # -----------------------------------------------------------------------
    # 7. Two-stage enzyme-turnover model for CMQ autoinduction:
    #      precursor1 = enzymatic-RNA precursor pool (CMQ-induced synthesis)
    #      precursor2 = metabolizing enzyme pool (precursor1-driven synthesis;
    #                   linearly modulates cl_cmq above)
    #    Both pools share the same first-order degradation rate kdeg (paper
    #    reports a single Kdeg). With baseline Cc_cmq = 0 the induction
    #    factor (1 + Emax * Cc/(EC50 + Cc)) reduces to 1, and the steady-
    #    state value of each pool is 1; both states are initialised at 1.
    # -----------------------------------------------------------------------
    ind_factor      <- 1 + emax * Cc_cmq / (ec50 + Cc_cmq)
    d/dt(precursor1) <- kdeg * ind_factor - kdeg * precursor1
    d/dt(precursor2) <- kdeg * precursor1 - kdeg * precursor2

    precursor1(0) <- 1
    precursor2(0) <- 1

    # -----------------------------------------------------------------------
    # 8. Bioavailability / regimen effect. CL/F and Vc/F are apparent-F
    #    relative, so the baseline bioavailability is 1; the split-dose
    #    indicator carries the +5% increment uniformly across both parent
    #    depots. The user supplies two dose records per administration
    #    (cmt = depot_r and cmt = depot_s) each with amt = half the racemate
    #    molar dose in nmol.
    # -----------------------------------------------------------------------
    f(depot_r) <- cov_split_f
    f(depot_s) <- cov_split_f

    # -----------------------------------------------------------------------
    # 9. Plasma observations (cord-arm fits dropped from the model file;
    #    see vignette Assumptions and deviations for the CORD = 1 / per-
    #    matrix residual-error treatment in the source paper).
    # -----------------------------------------------------------------------
    Cc_r   ~ prop(propSd_r)   + add(addSd_r)
    Cc_s   ~ prop(propSd_s)   + add(addSd_s)
    Cc_cmq ~ prop(propSd_cmq) + add(addSd_cmq)
  })
}
