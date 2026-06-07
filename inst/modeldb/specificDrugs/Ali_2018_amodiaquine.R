Ali_2018_amodiaquine <- function() {
  description <- "Joint parent-metabolite population PK model for oral amodiaquine and its CYP2C8-derived active metabolite desethylamodiaquine in adults and children with uncomplicated Plasmodium malaria, pooled across five WWARN cohorts (Burkina Faso, Ghana, Kenya, Uganda, Thailand). Two-transit absorption (NN = 2) into a 2-compartment amodiaquine disposition model with complete in-vivo conversion (with MW correction) to a 3-compartment desethylamodiaquine disposition model. Allometric body-weight scaling on CL/Q (exponent 0.75) and Vc/Vp (exponent 1.0) referenced at WT = 50 kg; sigmoidal postmenstrual-age maturation on both amodiaquine and desethylamodiaquine clearance; 22.4% lower bioavailability on the first daily dose relative to subsequent doses."
  reference <- paste(
    "Ali AM, Penny MA, Smith TA, Workman L, Sasi P, Adjei GO, Aweeka F,",
    "Kiechel JR, Jullien V, Rijken MJ, McGready R, Mwesigwa J,",
    "Kristensen K, Stepniewska K, Tarning J, Barnes KI, Denti P;",
    "WWARN Amodiaquine PK Study Group.",
    "Population pharmacokinetics of the antimalarial amodiaquine:",
    "a pooled analysis to optimize dosing.",
    "Antimicrob Agents Chemother. 2018;62(10):e02193-17.",
    "doi:10.1128/AAC.02193-17.",
    sep = " "
  )
  vignette <- "Ali_2018_amodiaquine"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed at baseline. Reference WT = 50 kg = pooled-cohort median (Ali 2018 Table 1). Allometric exponents 0.75 on CL and Q (for both AQ and DEAQ) and 1.0 on Vc and Vp (Ali 2018 Methods, 'Effect of body size and age').",
      source_name        = "WT"
    ),
    PAGE = list(
      description        = "Postmenstrual age (gestational age + postnatal age)",
      units              = "months",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Drives the sigmoidal maturation function on CL_AQ and CL_DEAQ. Per Ali 2018 Methods ('Effect of body size and age'), PMA is computed as PNA_months + 9 months assuming term gestation, because individual gestational ages were not available; for subjects with known GA the standard PAGE = GA_weeks/4.35 + PNA_months computation applies. Adults: PAGE = AGE_years*12 + 9.",
      source_name        = "PMA"
    ),
    CYCLE = list(
      description        = "Dose-number counter (1 = first daily dose of a 3-day course; 2 or 3 = second/third dose)",
      units              = "(count)",
      type               = "count",
      reference_category = NULL,
      notes              = "Used in a piecewise (CYCLE == 1) form to encode the first-day relative bioavailability reduction. Ali 2018 Table 3 reports an estimated -22.4% (95% CI -32.0% to -15.6%) reduction in F on day 1 relative to days 2 and 3 of the 3-day treatment course; the source paper attributes this to a transient malaria disease effect on absorption.",
      source_name        = "OCC"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 261,
    n_studies      = 6,
    n_pregnant     = 26,
    age_range      = "1-60 years",
    age_median     = "7.6 years",
    weight_range   = "6.5-93 kg",
    weight_median  = "21 kg",
    sex_female_pct = 53.3,
    disease_state  = "Uncomplicated Plasmodium falciparum or Plasmodium vivax malaria; includes 26 pregnant women in the second or third trimester.",
    dose_range     = "Oral amodiaquine (base, mg) targeting ~10 mg/kg once daily for 3 days, given alone or as artesunate + amodiaquine fixed-dose combination or as loose tablets; all administered doses are reported on the amodiaquine-base scale (Ali 2018 Materials and Methods).",
    regions        = "Burkina Faso, Ghana, Kenya, Uganda, Thailand (pooled via the WWARN Amodiaquine PK Study Group)",
    notes          = "Pooled patient-level data from five cohorts (six studies): Tarning 2008/2012 (pregnant women, Thailand, n = 26 + 7 resampled post-delivery), Jullien 2010 (adults, Uganda, n = 53), Mwesigwa 2010 (children 5-13 yr, Uganda, n = 20), Stepniewska 2009 (children 1-5 yr, Burkina Faso, n = 61), Adjei 2008 (children 1-14 yr, Ghana, n = 101). 36.4% of subjects are children under 5 years; 12.6% under 2 years. Baseline demographics summarised from Ali 2018 Table 1."
  )

  ini({
    # Amodiaquine structural parameters (apparent values; F is fixed at 1).
    # All clearances and volumes refer to a patient weighing WT = 50 kg,
    # the median weight in the pooled data set. Source: Ali 2018 Table 3.
    lka     <- log(0.589);  label("Absorption rate constant Ka (1/h); typical for WT = 50 kg")                                  # Table 3
    lmtt    <- log(0.236);  label("Mean transit time MTT (h)")                                                                  # Table 3
    lcl     <- log(2960);   label("Apparent amodiaquine clearance CL_AQ/F at WT = 50 kg (L/h)")                                 # Table 3
    lvc     <- log(13500);  label("Apparent amodiaquine central volume Vc_AQ/F at WT = 50 kg (L)")                              # Table 3
    lq      <- log(2310);   label("Apparent amodiaquine inter-compartmental clearance Q_AQ/F at WT = 50 kg (L/h)")              # Table 3
    lvp     <- log(22700);  label("Apparent amodiaquine peripheral volume Vp_AQ/F at WT = 50 kg (L)")                           # Table 3
    lfdepot <- fixed(log(1)); label("Reference relative bioavailability F (unitless; fixed at 1 per Ali 2018 Table 3)")          # Table 3, F = 1 Fixed

    # Desethylamodiaquine structural parameters (apparent values).
    lcl_deaq  <- log(32.6);  label("Apparent desethylamodiaquine clearance CL_DEAQ/F at WT = 50 kg (L/h)")                       # Table 3
    lvc_deaq  <- log(258);   label("Apparent desethylamodiaquine central volume Vc_DEAQ/F at WT = 50 kg (L)")                    # Table 3
    lq_deaq   <- log(154);   label("Apparent desethylamodiaquine inter-compartmental clearance Q1_DEAQ/F at WT = 50 kg (L/h); shallow peripheral exchange")  # Table 3
    lvp_deaq  <- log(2460);  label("Apparent desethylamodiaquine first peripheral volume Vp1_DEAQ/F at WT = 50 kg (L)")          # Table 3
    lq2_deaq  <- log(31.3);  label("Apparent desethylamodiaquine inter-compartmental clearance Q2_DEAQ/F at WT = 50 kg (L/h); deep peripheral exchange")  # Table 3
    lvp2_deaq <- log(5580);  label("Apparent desethylamodiaquine second peripheral volume Vp2_DEAQ/F at WT = 50 kg (L)")         # Table 3

    # Allometric exponents (fixed at canonical values 0.75 / 1.0).
    e_wt_cl <- fixed(0.75); label("Allometric exponent on CL_AQ, Q_AQ, CL_DEAQ, Q1_DEAQ, Q2_DEAQ (unitless; fixed)")  # Methods, 'Effect of body size and age'
    e_wt_vc <- fixed(1.0);  label("Allometric exponent on Vc_AQ, Vp_AQ, Vc_DEAQ, Vp1_DEAQ, Vp2_DEAQ (unitless; fixed)")  # Methods, 'Effect of body size and age'

    # Maturation parameters (estimated under weakly informative NONMEM priors;
    # Hill-type sigmoidal in postmenstrual age from conception).
    pma50_aq   <- 11.8; label("PMA at 50% AQ clearance maturation, time from conception (months)")    # Table 3
    hill_aq    <- 3.6;  label("Hill exponent for AQ clearance maturation (unitless)")                  # Table 3
    pma50_deaq <- 12.9; label("PMA at 50% DEAQ clearance maturation, time from conception (months)")   # Table 3
    hill_deaq  <- 3.22; label("Hill exponent for DEAQ clearance maturation (unitless)")                # Table 3

    # First-dose bioavailability reduction; applied via (CYCLE == 1).
    e_cycle1_fdepot <- -0.224; label("Fractional change in F on the first daily dose, applied multiplicatively as (1 + e_cycle1_fdepot * (CYCLE == 1)); -0.224 = 22.4% lower F on day 1")  # Table 3

    # Between-subject variability (BSV) on log-scale parameters. The paper
    # reports BSV as approximate CV%; convert via omega^2 = log(CV^2 + 1).
    etalcl      ~ 0.0987   # Ali 2018 Table 3, BSV_CL_AQ = 32.2% CV -> log(0.322^2 + 1) = 0.0987
    etalvc      ~ 0.2484   # Ali 2018 Table 3, BSV_Vc_AQ = 53.1% CV -> log(0.531^2 + 1) = 0.2484
    etalcl_deaq ~ 0.0392   # Ali 2018 Table 3, BSV_CL_DEAQ = 20.0% CV -> log(0.200^2 + 1) = 0.0392
    etalvc_deaq ~ 0.3727   # Ali 2018 Table 3, BSV_Vc_DEAQ = 67.2% CV -> log(0.672^2 + 1) = 0.3727

    # Residual error -- combined additive + proportional per analyte. Ali
    # 2018 fixed the additive error in each study to LLOQ/5 plus an
    # estimated parameter (significantly > 0 only for AQ). For a single
    # representative encoding, the most well-represented study LLOQ (Adjei
    # 2008 Ghana, n = 101 of 261 patients) of 10 ng/mL is used so that
    # LLOQ/5 = 2 ng/mL; the AQ extra estimated additive parameter (0.445
    # ng/mL) is added on top.
    addSd       <- 2.445;  label("Additive residual SD for amodiaquine (ng/mL; = LLOQ/5 with study-typical LLOQ = 10 ng/mL + 0.445 estimated)")  # Table 3, AQ Additive error: LLOQ/5 + 0.445
    propSd      <- 0.199;  label("Proportional residual SD for amodiaquine (CV fraction; = 19.9% per Ali 2018 Table 3)")  # Table 3, AQ Proportional error
    addSd_deaq  <- 2.0;    label("Additive residual SD for desethylamodiaquine (ng/mL; = LLOQ/5 with study-typical LLOQ = 10 ng/mL; estimated extra not different from zero per Ali 2018 Table 3 footnote d)")  # Table 3, DEAQ Additive error: LLOQ/5 Fixed
    propSd_deaq <- 0.242;  label("Proportional residual SD for desethylamodiaquine (CV fraction; = 24.2% per Ali 2018 Table 3)")  # Table 3, DEAQ Proportional error
  })

  model({
    # Molecular weight of the free base of each analyte (g/mol). The
    # source paper notes (Materials and Methods, 'Structural model') that
    # amodiaquine is assumed to be completely and irreversibly
    # metabolised to desethylamodiaquine and that a molar conversion
    # factor is applied to translate the mass flux of AQ leaving the AQ
    # central compartment into the mass flux of DEAQ entering the DEAQ
    # central compartment. The ratio mwDEAQ / mwAQ = 0.9212.
    mwAQ        <- 355.85
    mwDEAQ      <- 327.81
    molarFactor <- mwDEAQ / mwAQ

    # Sigmoidal Hill-type maturation of clearance on postmenstrual age
    # (months from conception); approaches 1 in adults and reaches 50%
    # at pma50_*.
    matAQ   <- PAGE^hill_aq   / (PAGE^hill_aq   + pma50_aq^hill_aq)
    matDEAQ <- PAGE^hill_deaq / (PAGE^hill_deaq + pma50_deaq^hill_deaq)

    # Individual PK parameters: allometric weight scaling (reference 50 kg)
    # combined with the maturation factor on each clearance.
    ka  <- exp(lka)
    mtt <- exp(lmtt)
    ktr <- 3 / mtt   # NN = 2 transit compartments fixed; KTR = (NN + 1) / MTT

    cl <- exp(lcl + etalcl) * (WT / 50)^e_wt_cl * matAQ
    vc <- exp(lvc + etalvc) * (WT / 50)^e_wt_vc
    q  <- exp(lq)           * (WT / 50)^e_wt_cl
    vp <- exp(lvp)          * (WT / 50)^e_wt_vc

    cl_deaq  <- exp(lcl_deaq + etalcl_deaq) * (WT / 50)^e_wt_cl * matDEAQ
    vc_deaq  <- exp(lvc_deaq + etalvc_deaq) * (WT / 50)^e_wt_vc
    q_deaq   <- exp(lq_deaq)   * (WT / 50)^e_wt_cl
    vp_deaq  <- exp(lvp_deaq)  * (WT / 50)^e_wt_vc
    q2_deaq  <- exp(lq2_deaq)  * (WT / 50)^e_wt_cl
    vp2_deaq <- exp(lvp2_deaq) * (WT / 50)^e_wt_vc

    # Micro-rate constants for the ODE system (1/h).
    kel_aq   <- cl       / vc
    k12_aq   <- q        / vc
    k21_aq   <- q        / vp
    kel_deaq <- cl_deaq  / vc_deaq
    k12_deaq <- q_deaq   / vc_deaq
    k21_deaq <- q_deaq   / vp_deaq
    k13_deaq <- q2_deaq  / vc_deaq
    k31_deaq <- q2_deaq  / vp2_deaq

    # ODE system. Compartment amounts are in mg of analyte base; volumes
    # are in L; concentrations are formed below as 1000 * amount / volume
    # to convert mg/L to ng/mL.
    #
    # Absorption chain: depot -> transit1 -> transit2 -> central (AQ);
    # NN = 2 transit compartments per Ali 2018 Table 3 (NN = 2.00,
    # 95% CI 1.09-6.31; the estimated value is treated as fixed at 2 for
    # implementation).
    d/dt(depot)    <- -ktr * depot
    d/dt(transit1) <-  ktr * depot    - ktr * transit1
    d/dt(transit2) <-  ktr * transit1 - ktr * transit2

    # Amodiaquine central + one peripheral. Mass flux leaving AQ central
    # equals (kel_aq * central) and is fully routed to the DEAQ central
    # compartment, after a mwDEAQ / mwAQ molar correction, under the
    # complete-conversion assumption.
    d/dt(central)     <-  ktr * transit2 - kel_aq * central -
                          k12_aq * central + k21_aq * peripheral1
    d/dt(peripheral1) <-  k12_aq * central - k21_aq * peripheral1

    # Desethylamodiaquine central + two peripherals.
    d/dt(central_deaq) <-  molarFactor * kel_aq * central -
                           kel_deaq * central_deaq -
                           k12_deaq * central_deaq + k21_deaq * peripheral1_deaq -
                           k13_deaq * central_deaq + k31_deaq * peripheral2_deaq
    d/dt(peripheral1_deaq) <- k12_deaq * central_deaq - k21_deaq * peripheral1_deaq
    d/dt(peripheral2_deaq) <- k13_deaq * central_deaq - k31_deaq * peripheral2_deaq

    # Bioavailability on the depot compartment. Reference F is fixed at 1
    # (per Ali 2018 Table 3 "F  1 Fixed"). The first-day effect drops F
    # to 1 + (-0.224) = 0.776 when CYCLE = 1; doses with CYCLE >= 2 keep
    # the reference value.
    f(depot) <- exp(lfdepot) * (1 + e_cycle1_fdepot * (CYCLE == 1))

    # Plasma observations in ng/mL: amount (mg) divided by volume (L)
    # gives mg/L = 10^3 ng/mL; multiply by 1000.
    Cc      <- 1000 * central      / vc
    Cc_deaq <- 1000 * central_deaq / vc_deaq

    Cc      ~ add(addSd)      + prop(propSd)
    Cc_deaq ~ add(addSd_deaq) + prop(propSd_deaq)
  })
}
