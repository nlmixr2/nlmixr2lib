Abuhelwa_2015_itraconazole <- function() {
  description <- "Population PK model for oral itraconazole and its active metabolite hydroxy-itraconazole in healthy adults (Abuhelwa 2015). Two-compartment parent with 4-transit-compartment Savic-style absorption and a one-compartment hydroxy-itraconazole metabolite eliminated by mixed linear and Michaelis-Menten kinetics. Encodes the SUBA-itraconazole vs Sporanox formulation effect on relative bioavailability (with formulation-dependent scaling of the F variability) and the fed-vs-fasted effect on both relative bioavailability and the transit-absorption rate constant; the metabolic conversion ratio fm is assumed = 1 so all parent clearance becomes metabolite, and the metabolite CL/V are apparent values scaled by the unknown fm."
  reference <- paste(
    "Abuhelwa AY, Foster DJR, Mudge S, Hayes D, Upton RN.",
    "Population pharmacokinetic modeling of itraconazole and",
    "hydroxyitraconazole for oral SUBA-itraconazole and Sporanox",
    "capsule formulations in healthy subjects in fed and fasted",
    "states. Antimicrob Agents Chemother. 2015;59(9):5681-5696.",
    "doi:10.1128/AAC.00973-15.",
    sep = " "
  )
  vignette <- "Abuhelwa_2015_itraconazole"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  # The FVAR random effect is a paper-mechanistic "common variability on
  # bioavailability" eta (Abuhelwa 2015 Methods 'Base model development of
  # single-dose data', stage 1; supplement Appendix S1 ETA(FVAR)). It is a
  # shared log-normal random effect on F1 with no corresponding fixed-effect
  # parameter (the typical value of FVAR is 1.0 by construction). Declare it
  # via paper_specific_etas so checkModelConventions() does not flag the
  # missing matching `fvar` typical-value parameter.
  paper_specific_etas <- c("etafvar")

  covariateData <- list(
    FED = list(
      description        = "Fed state at the time of itraconazole dose administration. 1 = fed (high-fat high-calorie breakfast eaten before dose); 0 = fasted (overnight fast of >= 10 h).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (fasted state).",
      notes              = "Per-dose-record indicator from the Abuhelwa 2015 crossover trials. Multiplicative effects on KTR (FEDKTR = -0.583, i.e. KTR_fed = 0.417 x KTR_fasted) and on relative bioavailability (FEDF = -0.269, i.e. F_fed = 0.731 x F_fasted), independent of formulation; see Abuhelwa 2015 Table 3 rows FEDKTR and FEDF.",
      source_name        = "FED"
    ),
    FORM_SUBA = list(
      description        = "SUBA-itraconazole vs Sporanox capsule formulation indicator. 1 = SUBA-itraconazole (solid dispersion in a pH-dependent polymeric matrix); 0 = Sporanox capsule (innovator product, the structural reference with relative bioavailability fixed to 1).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (Sporanox capsule; relative bioavailability F = 1, FVAR variability scaling ETASCALE = 1).",
      notes              = "Per-dose-record indicator from the Abuhelwa 2015 crossover trials (DRUG column in Appendix S1; 0 = Sporanox, 1 = SUBA-itraconazole). Multiplicative effect on relative bioavailability (FORMF = 0.729, i.e. F_SUBA = 1.729 x F_Sporanox = +73% relative bioavailability) and a separate multiplicative scaling of the FVAR random effect (ETASCALE = -0.213, i.e. F-variability for SUBA = 0.787 x F-variability for Sporanox) per Abuhelwa 2015 Table 3 rows FORMF and ETASCALE.",
      source_name        = "DRUG"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 238L,
    n_studies      = 7L,
    age_range      = "18-60 years (US studies median 38 y; UK studies median 28 y; combined Abuhelwa 2015 Table 2)",
    weight_range   = "49-108 kg (US studies 49-108 kg; UK studies 49.3-98.3 kg; combined Abuhelwa 2015 Table 2)",
    sex_female_pct = 60.9,
    race_ethnicity = "Combined US + UK cohort: 76% White (182 of 238), 11% Black (25), 0.4% Asian (1), 13% Other (30) per Abuhelwa 2015 Table 2.",
    disease_state  = "Healthy adult volunteers, non-tobacco-using, non-pregnant non-lactating, with no GI disease / malabsorption history and no CYP3A inhibitor or inducer exposure within 30 days of study entry.",
    dose_range     = "Oral single doses of 50, 65, 100, or 200 mg itraconazole (SUBA-itraconazole capsule or Sporanox capsule) in fed and fasted crossover periods; multidose phase used 100 or 200 mg SUBA-itraconazole / 200 or 400 mg Sporanox daily for 15 days under fed conditions (the multidose-phase dose- and time-dependent CL and F mechanisms are not implemented in this single-dose extraction; see vignette Assumptions and deviations).",
    regions        = "United States (5 trials, 154 subjects) and United Kingdom (2 trials, 84 subjects); Sporanox-source country (USA vs UK) is collinear with the study population indicator.",
    notes          = "Pooled phase I crossover bioequivalence data from seven Mayne Pharma International trials (MPG009, HGN007, HGN008, 10850702, 10850703, 10850705, 10850706). The data set contained 15,097 itraconazole plasma concentrations and 9,868 hydroxy-itraconazole plasma concentrations; BLQ samples (~12.5% of itraconazole, ~17.5% of hydroxy-itraconazole) were excluded via the M1 method. Sample collection: 0-72 h after dose for most studies, with HGN008 and MPG009 extending to 96-120 h."
  )

  ini({
    # ========================================================================
    # Structural parameters of the parent itraconazole 2-compartment model
    # (Abuhelwa 2015 Table 3 final estimates from the combined single- and
    # multidose itraconazole model). All clearance / volume parameters are
    # apparent (scaled by oral bioavailability F).
    # ========================================================================
    lcl  <- log(129);   label("Apparent oral clearance of itraconazole CL/F (L/h)")                                 # Abuhelwa 2015 Table 3 row 'CL/F (liters/h) = 129'
    lvc  <- log(861);   label("Apparent central volume of distribution of itraconazole V2/F (L)")                    # Abuhelwa 2015 Table 3 row 'V2/F (liters) = 861'
    lq   <- log(153);   label("Apparent inter-compartmental clearance of itraconazole Q/F (L/h)")                    # Abuhelwa 2015 Table 3 row 'Q/F (liters/h) = 153'
    lvp  <- log(2340);  label("Apparent peripheral volume of distribution of itraconazole V3/F (L)")                 # Abuhelwa 2015 Table 3 row 'V3/F (liters) = 2,340'
    lktr <- log(2.05);  label("Transit-compartment rate constant KTR for the 4-transit absorption chain (1/h)")      # Abuhelwa 2015 Table 3 row 'KTR (h-1) = 2.05'; mean transit time MTT = 5 / KTR = 2.44 h (population fasted, US-study reference)

    # ========================================================================
    # Formulation effect on bioavailability (Sporanox is the structural
    # reference with F = 1; SUBA-itraconazole has 1 + FORMF = 1.729 relative
    # bioavailability). The covariate FORM_SUBA = 1 selects the SUBA arm; for
    # Sporanox FORM_SUBA = 0 and the F factor degenerates to 1.
    # ========================================================================
    e_suba_f <- 0.729;  label("Multiplicative effect of SUBA-itraconazole vs Sporanox on relative bioavailability (unitless; F = 1 + e_suba_f x FORM_SUBA)")  # Abuhelwa 2015 Table 3 row 'FORMF = 0.729'; reported as 'SUBA-itraconazole had a relative bioavailability of 173% compared to that of Sporanox' (Results paragraph after Table 3)

    # ========================================================================
    # ETASCALE: SUBA-itraconazole has 1 + ETASCALE = 0.787 less F variability
    # than Sporanox (i.e. 21.3% less variable per the published Abstract). The
    # scaling enters the random effect on F as exp(etafvar * etascale_eff)
    # with etascale_eff = 1 + e_suba_etafvar x FORM_SUBA.
    # ========================================================================
    e_suba_etafvar <- -0.213;  label("Multiplicative scaling of the FVAR random effect for SUBA vs Sporanox (unitless; etascale = 1 + e_suba_etafvar x FORM_SUBA)")  # Abuhelwa 2015 Table 3 row 'ETASCALE = -0.213'; corresponds to 21.3% reduced F variability for SUBA per Abstract ('the bioavailability was 21% less variable between subjects')

    # ========================================================================
    # Fed-state effects (FED = 1 in fed periods, 0 in fasted periods). Both
    # effects are negative (food reduces F and slows the transit absorption
    # rate constant), enter as multiplicative (1 + e_fed_<param> x FED).
    # ========================================================================
    e_fed_f   <- -0.269;  label("Multiplicative effect of fed status on relative bioavailability (unitless; F factor = 1 + e_fed_f x FED)")    # Abuhelwa 2015 Table 3 row 'FEDF = -0.269'; reported as 'a decrease in the relative bioavailability to 73.1% of the fasted state' (Results, formulation section)
    e_fed_ktr <- -0.583;  label("Multiplicative effect of fed status on KTR (unitless; KTR factor = 1 + e_fed_ktr x FED)")                       # Abuhelwa 2015 Table 3 row 'FEDKTR = -0.583'; reported as 'a decrease in transit absorption rate constant (KTR) to 41.7% of the fasted state', corresponding to MTT increasing from 1.95 h fasted to 4.68 h fed

    # ========================================================================
    # Structural parameters of the hydroxy-itraconazole metabolite
    # (Abuhelwa 2015 Table 4 final estimates from the single-dose fit; the
    # paper's separate multidose TIMECLM factor is documented in the vignette
    # Assumptions and deviations section but not implemented here). All
    # metabolite CL / V values are apparent values scaled by the unknown
    # metabolic conversion ratio fm; the parent-to-metabolite stoichiometric
    # ratio is fixed to 1 in the model (every cleared parent molecule becomes
    # metabolite).
    # NOTE on V_max units: Abuhelwa 2015 Table 4 labels V_max units as
    # "liters/h", but the Appendix S1 NONMEM control stream states
    # ";UNITS ARE UG, L (NG/ML) AND H" and Appendix S2 implements the term
    # as ((VMAX * C8) / (KM + C8)) in DADT(8). With amount in micrograms and
    # concentration in ng/mL = ug/L, V_max must have mass/time units
    # (micrograms/h) for the V_max * C / (K_m + C) flux term to be in
    # micrograms/h. The natural reading of the published V_max value 403
    # is therefore 403 ug/h, not 403 L/h; the Table 4 unit label is a
    # presentation error and the value 403 is used here in micrograms/h.
    # ========================================================================
    lcl_ohi   <- log(45.6);  label("Apparent oral clearance of hydroxy-itraconazole CL_m/F_m (L/h)")                             # Abuhelwa 2015 Table 4 row 'CL_m/F_m (liters/h) = 45.6'
    lvc_ohi   <- log(5.28);  label("Apparent central volume of distribution of hydroxy-itraconazole V_1m/F_m (L)")                # Abuhelwa 2015 Table 4 row 'V_1m/F_m (liters) = 5.28'
    lvmax_ohi <- log(403);   label("Maximum elimination rate of the hydroxy-itraconazole Michaelis-Menten arm V_max (ug/h)")     # Abuhelwa 2015 Table 4 row 'V_max (liters/h) = 403' (units corrected to ug/h per Appendix S1 ';UNITS ARE UG, L (NG/ML) AND H' header; see vignette Assumptions and deviations)
    lkm_ohi   <- log(1.64);  label("Michaelis-Menten constant K_m (ng/mL)")                                                       # Abuhelwa 2015 Table 4 row 'K_m (ng/ml) = 1.64'

    # ========================================================================
    # IIV: parent has BSV on CL/F, V2/F, KTR and a "common variability on
    # bioavailability" parameter (FVAR) shared across all bioavailability-
    # affected parameters. Metabolite has BSV on CL_m only. The source
    # reports approximate %CV (footnote of Table 3 / Table 4); nlmixr2lib
    # uses variances on the log scale, converted via
    #     omega^2 = log(1 + CV^2)
    # ========================================================================
    etalcl     ~ 0.04681  # Abuhelwa 2015 Table 3 BSV CL/F = 22.1%; omega^2 = log(1 + 0.221^2) = 0.04681
    etalvc     ~ 0.09011  # Abuhelwa 2015 Table 3 BSV V2/F = 30.8%; omega^2 = log(1 + 0.308^2) = 0.09011
    etalktr    ~ 0.18120  # Abuhelwa 2015 Table 3 BSV KTR  = 44.8%; omega^2 = log(1 + 0.448^2) = 0.18120
    etafvar    ~ 0.27259  # Abuhelwa 2015 Table 3 FVAR     = 56.4% (common random effect on bioavailability); omega^2 = log(1 + 0.564^2) = 0.27259
    etalcl_ohi ~ 0.07497  # Abuhelwa 2015 Table 4 BSV CL_m/F_m = 27.9%; omega^2 = log(1 + 0.279^2) = 0.07497

    # ========================================================================
    # Residual error: the source uses a log-transform-both-sides (LTBS)
    # combined additive + proportional model (Methods Eq. 2 and Appendix S1),
    # which is equivalent in linear space to add() + prop() for the
    # concentration magnitudes used in this trial. The source reports
    # different residual-error magnitudes for single-dose fasted vs
    # single-dose fed vs multidose; this model uses the single-dose fasted
    # values for the parent (the primary bioequivalence-simulation scenario
    # in the published paper) and the single-dose values for the metabolite.
    # The single-dose fed and multidose values are documented as deviations
    # in the vignette Assumptions and deviations section.
    # ========================================================================
    propSd     <- 0.293;   label("Itraconazole proportional residual error, single-dose fasted (fraction)")   # Abuhelwa 2015 Table 3 row 'Proportional error-fasted (%CV) = 29.3'
    addSd      <- 0.168;   label("Itraconazole additive residual error, single-dose fasted (ng/mL)")            # Abuhelwa 2015 Table 3 row 'Additive error-fasted (ng/ml) = 0.168'
    propSd_ohi <- 0.387;   label("Hydroxy-itraconazole proportional residual error, single-dose (fraction)")    # Abuhelwa 2015 Table 4 row 'Single-dose studies / Proportional error (%CV) = 38.7'
    addSd_ohi  <- 0.329;   label("Hydroxy-itraconazole additive residual error, single-dose (ng/mL)")           # Abuhelwa 2015 Table 4 row 'Single-dose studies / Additive error (ng/ml) = 0.329'
  })

  model({
    # ----- 1. Derived covariate-effect terms ----------------------------------
    # Formulation effect on the typical-value relative bioavailability:
    #     F1_TV = 1 for Sporanox, = 1 + e_suba_f = 1.729 for SUBA-itraconazole.
    drugf_eff <- 1 + e_suba_f * FORM_SUBA

    # Fed-state effects on F and on KTR:
    #     fedf_eff   = 1 for fasted, = 1 + e_fed_f   = 0.731 for fed.
    #     fedktr_eff = 1 for fasted, = 1 + e_fed_ktr = 0.417 for fed.
    fedf_eff   <- 1 + e_fed_f   * FED
    fedktr_eff <- 1 + e_fed_ktr * FED

    # FVAR random-effect scaling: ETASCALE = 1 for Sporanox, = 1 + e_suba_etafvar = 0.787 for SUBA.
    etascale_eff <- 1 + e_suba_etafvar * FORM_SUBA

    # ----- 2. Individual PK parameters ---------------------------------------
    # Parent disposition.
    cl  <- exp(lcl  + etalcl)
    vc  <- exp(lvc  + etalvc)
    q   <- exp(lq)
    vp  <- exp(lvp)
    ktr <- exp(lktr + etalktr) * fedktr_eff

    # "Common variability on bioavailability" (Abuhelwa 2015 FVAR) -- shared
    # random effect on F across all formulations and fed states, scaled by
    # ETASCALE for SUBA-itraconazole.
    fvar <- exp(etafvar * etascale_eff)

    # Typical-value F times shared FVAR random effect (depot bioavailability).
    fbio <- drugf_eff * fedf_eff * fvar

    # Metabolite disposition.
    cl_ohi   <- exp(lcl_ohi + etalcl_ohi)
    vc_ohi   <- exp(lvc_ohi)
    vmax_ohi <- exp(lvmax_ohi)
    km_ohi   <- exp(lkm_ohi)

    # ----- 3. Micro-constants -------------------------------------------------
    kel <- cl / vc
    k23 <- q  / vc
    k32 <- q  / vp

    # ----- 4. ODE system ------------------------------------------------------
    # depot --(KTR)-> transit1 --(KTR)-> transit2 --(KTR)-> transit3 --(KTR)-> transit4 --(KTR)-> central
    # central <-> peripheral1 (k23 / k32), central --(kel)-> metabolite, metabolite --(linear CL + Michaelis-Menten)-> out
    d/dt(depot)       <- -ktr * depot
    d/dt(transit1)    <-  ktr * depot      - ktr * transit1
    d/dt(transit2)    <-  ktr * transit1   - ktr * transit2
    d/dt(transit3)    <-  ktr * transit2   - ktr * transit3
    d/dt(transit4)    <-  ktr * transit3   - ktr * transit4
    d/dt(central)     <-  ktr * transit4   - kel * central  - k23 * central + k32 * peripheral1
    d/dt(peripheral1) <-  k23 * central    - k32 * peripheral1

    # Metabolite: parent elimination flux enters as formation (fm = 1 assumed),
    # mixed linear (CL_m / V_m) + Michaelis-Menten (V_max C / (K_m + C)) elimination.
    cm <- central_ohi / vc_ohi
    d/dt(central_ohi) <- kel * central - (cl_ohi / vc_ohi) * central_ohi - vmax_ohi * cm / (km_ohi + cm)

    # ----- 5. Bioavailability and dose-unit conversion ------------------------
    # Convert the user-facing mg dose to ug-internal amounts so the
    # Michaelis-Menten V_max in ug/h (see lvmax_ohi label) is dimensionally
    # consistent with the linear elimination flux (cl_ohi/vc_ohi)*central_ohi
    # which is in ug/h when central_ohi is in ug. With internal amounts in ug
    # and volumes in L, the observation A/V is in ug/L = ng/mL natively,
    # matching the source paper's concentration units.
    f(depot) <- fbio * 1000

    # ----- 6. Observations and residual error --------------------------------
    Cc     <- central     / vc      # ug / L = ng/mL
    Cc_ohi <- central_ohi / vc_ohi  # ug / L = ng/mL

    Cc     ~ add(addSd)     + prop(propSd)
    Cc_ohi ~ add(addSd_ohi) + prop(propSd_ohi)
  })
}
