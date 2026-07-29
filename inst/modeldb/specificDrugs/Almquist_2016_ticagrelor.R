Almquist_2016_ticagrelor <- function() {
  description <- "Preclinical (mouse, C57Bl/6 male). Mechanistic interaction PK model for ticagrelor, its active metabolite (TAM, AR-C124910XX), and the ticagrelor-neutralising Fab antibody fragment MEDI2452 in mouse (Almquist 2016). Three-compartment disposition for ticagrelor and TAM (shared plasma V, tissue V1, V2; V1 in instantaneous equilibrium with V); MEDI2452 lives in plasma V only and reversibly binds the free fractions of ticagrelor and TAM with rate kon and dissociation constant Kd; both free MEDI2452 and the two MEDI2452-drug complexes are eliminated together at the Fab clearance Cl_f (no recycling). Naive-pooled fit (no IIV); multiplicative log-normal residual error on five plasma assays."
  reference   <- "Almquist J, Penney M, Pehrsson S, Sandinge AS, Janefeldt A, Maqbool S, Madalli S, Goodman J, Nylander S, Gennemark P. Unraveling the Pharmacokinetic Interaction of Ticagrelor and MEDI2452 (Ticagrelor Antidote) by Mathematical Modeling. CPT Pharmacometrics Syst Pharmacol. 2016;5(6):313-323. doi:10.1002/psp4.12089"
  vignette    <- "Almquist_2016_ticagrelor"
  units       <- list(time = "min", dosing = "nmol/kg", concentration = "nmol/L")

  covariateData <- list()

  population <- list(
    species        = "mouse (C57Bl/6 male)",
    n_subjects     = "Four validation/refinement studies in non-fasted male C57Bl/6 mice (Charles River, Sulzfeld, Germany; body weight 15-25 g); per-study group sizes summarised in Figure 1 and Supplementary Text S1.",
    n_studies      = 4,
    age_range      = "not reported",
    weight_range   = "15-25 g",
    sex_female_pct = 0,
    disease_state  = "Healthy mice (preclinical antidote pharmacology); also includes a Sprague-Dawley rat MEDI2452 pre-study (1000 mg/kg IV bolus) used to seed the allometrically scaled MEDI2452 parameters before re-estimation in mouse.",
    dose_range     = "Ticagrelor: IV infusion 240 ug/min/kg for 5 min then 30 ug/min/kg for 15 min (steady-state design used in studies 1-4); ticagrelor pre-study IV bolus 2000 ug/kg. MEDI2452: IV bolus 25-300 mg/kg (study 4 dose range).",
    regions        = "Preclinical (AstraZeneca R&D Molndal, Sweden; MedImmune, Cambridge, UK)",
    notes          = "Naive-pooled estimation in MATLAB (fminsearch) with maximum-likelihood multiplicative log-normal residual error; bootstrapping (N = 300) for parameter uncertainty. Final parameter values from Table 1 'Estimated value' column (refined model)."
  )

  ini({
    # Plasma and tissue volumes (body-weight-normalised, L/kg).
    # V (plasma) is fixed to the standard mouse plasma volume per Table 1.
    lvc  <- fixed(log(0.05));   label("Plasma volume of distribution V (L/kg, fixed: standard mouse plasma volume)")            # Almquist 2016 Table 1, V = 0.05 L/kg (not estimated)
    lvp  <- log(1.12);          label("Volume V1 of rapidly exchanging tissue (L/kg)")                                          # Almquist 2016 Table 1, V1 refined estimate
    lvp2 <- log(1.8);           label("Volume V2 of slowly exchanging peripheral tissue (L/kg)")                                # Almquist 2016 Table 1, V2 refined estimate

    # Distributional clearances for ticagrelor and TAM between V/V1/V2.
    # Cl_fast (V <-> V1) is fixed to a value much larger than the other
    # clearances in the system so that V and V1 act as instantaneously
    # equilibrating sub-compartments of the apparent central compartment
    # (Almquist 2016 assumption II).
    lq2 <- fixed(log(10));      label("Fast distributional clearance Cl_fast between V and V1 (L/min/kg, fixed: instantaneous equilibrium)")  # Almquist 2016 Table 1, Cl_fast = 10 L/min/kg (not estimated)
    lq  <- log(0.041);          label("Distributional clearance Cl_d between V and V2 (L/min/kg)")                              # Almquist 2016 Table 1, Cl_d refined estimate

    # Ticagrelor non-metabolic elimination and metabolic conversion to TAM.
    # TAM is assumed to follow the same Cl and the same distributional
    # clearances as ticagrelor (Almquist 2016 assumption IV); only one
    # 'lcl' parameter is needed for the two species.
    lcl     <- log(0.022);      label("Non-metabolic elimination clearance Cl, shared by ticagrelor and TAM (L/min/kg)")        # Almquist 2016 Table 1, Cl refined estimate
    lcl_met <- log(0.0080);     label("Metabolic clearance Cl_met from ticagrelor to TAM (L/min/kg)")                           # Almquist 2016 Table 1, Cl_met refined estimate

    # MEDI2452 Fab clearance: free MEDI2452 and both MEDI2452-drug
    # complexes are cleared at the same rate (assumption III - the
    # complex is eliminated together with MEDI2452, no recycling).
    lcl_target <- log(0.0025);  label("MEDI2452 Fab clearance Cl_f, applied to free MEDI2452 and to both Fab-drug complexes (L/min/kg)")  # Almquist 2016 Table 1, Cl_f refined estimate

    # MEDI2452 binding to free ticagrelor / free TAM.
    lkon <- log(0.11);          label("Second-order MEDI2452-drug association rate k_on (1/(nM*min))")                          # Almquist 2016 Table 1, k_on refined estimate
    lkd  <- fixed(log(0.02));   label("MEDI2452-drug dissociation constant K_d (nM, fixed at in vitro affinity)")               # Almquist 2016 Table 1, K_d = 0.02 nM (not estimated; Buchanan 2015 in vitro affinity)

    # Fraction of ticagrelor / TAM unbound to plasma protein. Held at the
    # internal AstraZeneca measurement (n = 38), not estimated in the PK fit.
    f_unbound <- fixed(0.0020); label("Fraction of ticagrelor / TAM unbound to plasma protein (unitless, fixed: internal data)") # Almquist 2016 Table 1, f = 0.0020 (not estimated)

    # Residual error: multiplicative log-normal per the Data Analysis
    # paragraph (y_obs = y_pred * exp(eps), eps ~ N(0, expSd^2)).
    # Table 1 reports r^2 = sigma^2 (log-scale variance); expSd = sqrt(r^2).
    expSd             <- sqrt(0.076); label("Log-scale residual SD for total ticagrelor in plasma (unitless)")                    # Almquist 2016 Table 1, r^2_tica = 0.076
    expSd_tam         <- sqrt(0.080); label("Log-scale residual SD for total TAM in plasma (unitless)")                           # Almquist 2016 Table 1, r^2_TAM = 0.080
    expSd_freeMEDI    <- sqrt(0.28);  label("Log-scale residual SD for free MEDI2452 in plasma (unitless)")                       # Almquist 2016 Table 1, r^2_MEDI = 0.28
    expSd_freeTicaObs <- sqrt(0.042); label("Log-scale residual SD for free ticagrelor measured ex vivo via the observation model (unitless)") # Almquist 2016 Table 1, r^2_freetica = 0.042
    expSd_freeTamObs  <- sqrt(0.060); label("Log-scale residual SD for free TAM measured ex vivo via the observation model (unitless)")        # Almquist 2016 Table 1, r^2_freeTAM = 0.060
  })

  model({
    vc  <- exp(lvc)
    vp  <- exp(lvp)
    vp2 <- exp(lvp2)
    q   <- exp(lq)
    q2  <- exp(lq2)
    cl     <- exp(lcl)
    cl_met <- exp(lcl_met)
    cl_target <- exp(lcl_target)
    kon <- exp(lkon)
    kd  <- exp(lkd)

    # Concentrations of each species in each physical compartment (nM).
    # State variables hold amounts in nmol/kg; dividing by the L/kg volume
    # yields nmol/L = nM. central / complex / complex_tam all share the
    # plasma volume vc; the ticagrelor / TAM peripheral compartments live
    # in their own vp / vp2.
    cTica    <- central         / vc
    cTica1   <- peripheral1     / vp
    cTica2   <- peripheral2     / vp2
    cTam     <- central_tam     / vc
    cTam1    <- peripheral1_tam / vp
    cTam2    <- peripheral2_tam / vp2
    cFab     <- target          / vc
    cFabTica <- complex         / vc
    cFabTam  <- complex_tam     / vc

    # Reversible MEDI2452 binding flux (amount/time). Almquist 2016 Eqs
    # 1, 4, 7-9: V * k_on * (f * X * F - K_d * X_complex), where X is the
    # total (free + protein-bound) ticagrelor or TAM in V.
    bind_tica <- vc * kon * (f_unbound * cTica * cFab - kd * cFabTica)
    bind_tam  <- vc * kon * (f_unbound * cTam  * cFab - kd * cFabTam)

    # Ticagrelor (parent). Almquist 2016 Eqs 1-3.
    d/dt(central)     <- -q2 * (cTica - cTica1) - q  * (cTica - cTica2) -
                          bind_tica - cl_met * cTica - cl * cTica
    d/dt(peripheral1) <-  q2 * (cTica - cTica1)
    d/dt(peripheral2) <-  q  * (cTica - cTica2)

    # TAM (ticagrelor active metabolite). Almquist 2016 Eqs 4-6.
    # Formation from ticagrelor enters the V compartment (the metabolic
    # step converts ticagrelor in V into TAM in V).
    d/dt(central_tam)     <- -q2 * (cTam - cTam1) - q  * (cTam - cTam2) -
                              bind_tam + cl_met * cTica - cl * cTam
    d/dt(peripheral1_tam) <-  q2 * (cTam - cTam1)
    d/dt(peripheral2_tam) <-  q  * (cTam - cTam2)

    # MEDI2452 free and bound species in V. Almquist 2016 Eqs 7-9.
    d/dt(target)      <- -bind_tica - bind_tam - cl_target * cFab
    d/dt(complex)     <-  bind_tica            - cl_target * cFabTica
    d/dt(complex_tam) <-  bind_tam             - cl_target * cFabTam

    # Observation model (Almquist 2016 "Observation model" paragraph and
    # Supplementary Text S2): when blood is sampled, all clearances are
    # interrupted and the binding equilibria are reached before bioanalysis
    # completes. With totals T_t (tica), T_m (TAM), T_f (Fab) conserved
    # and the equilibrium K_d = f * X_eq * Fab_eq / X_complex_eq, the free
    # MEDI2452 concentration Fab_eq at the moment of bioanalysis is the
    # positive root of the quadratic
    #   f * Fab_eq^2 + Fab_eq * (K_d + f*(T_t + T_m - T_f)) - T_f * K_d = 0,
    # and the measured free ticagrelor / TAM are then
    #   free_tica_obs = f * T_t / (1 + f * Fab_eq / K_d)
    #   free_tam_obs  = f * T_m / (1 + f * Fab_eq / K_d).
    # For free MEDI2452 the observation-model correction is marginal per
    # the paper (the in-vivo cFab is reported as the assay output).
    totalTica <- cTica + cFabTica
    totalTam  <- cTam  + cFabTam
    totalFab  <- cFab  + cFabTica + cFabTam

    a_eq <- f_unbound
    b_eq <- kd + f_unbound * (totalTica + totalTam - totalFab)
    c_eq <- -totalFab * kd
    discrim <- b_eq * b_eq - 4 * a_eq * c_eq
    fabEq <- (-b_eq + sqrt(discrim)) / (2 * a_eq)

    denomEq <- 1 + f_unbound * fabEq / kd
    ticaEq  <- totalTica / denomEq
    tamEq   <- totalTam  / denomEq

    # Plasma observations.
    Cc            <- totalTica            # total ticagrelor (nM)
    Cc_tam        <- totalTam             # total TAM (nM)
    freeMEDI      <- cFab                 # free MEDI2452 (nM); observation-model correction is marginal
    freeTicaObs   <- f_unbound * ticaEq   # free ticagrelor measured ex vivo (nM)
    freeTamObs    <- f_unbound * tamEq    # free TAM measured ex vivo (nM)

    Cc          ~ lnorm(expSd)
    Cc_tam      ~ lnorm(expSd_tam)
    freeMEDI    ~ lnorm(expSd_freeMEDI)
    freeTicaObs ~ lnorm(expSd_freeTicaObs)
    freeTamObs  ~ lnorm(expSd_freeTamObs)
  })
}
