Wang_2010_romiplostim <- function() {
  description <- "Population PK/PD model for romiplostim in healthy subjects (Wang 2010 AAPS J). Pharmacodynamics-mediated drug disposition (PDMDD, a TMDD subtype) two-compartment quasi-equilibrium PK with first-order SC absorption, parallel linear (kel) and target-mediated (kint) elimination, coupled to a Krzyzanski-style cytokinetic precursor + platelet lifespan PD model with NP=10 megakaryocyte and NPLT=10 platelet age-compartments. Romiplostim free serum concentration stimulates platelet precursor production via a Hill function (Smax, SC50). The total c-Mpl receptor concentration is taken proportional to the circulating platelet count (Rtot = xi * PLT). Wang 2010 fit the model to MEAN PK and platelet-count data from 32 healthy subjects after single IV (0.3, 1, 10 ug/kg) or SC (0.1, 0.3, 1, 2 ug/kg) doses; no IIV was estimated (the population approach failed for this complex model)."
  reference <- "Wang YC, Krzyzanski W, Doshi S, Xiao JJ, Perez-Ruixo JJ, Chow AT. Pharmacodynamics-Mediated Drug Disposition (PDMDD) and Precursor Pool Lifespan Model for Single Dose of Romiplostim in Healthy Subjects. AAPS J. 2010 Dec;12(4):729-740. doi:10.1208/s12248-010-9234-9 (PMID 20963535)."
  vignette <- "Wang_2010_romiplostim"
  units <- list(time = "h", dosing = "ug", concentration = "ng/mL", platelet = "10^9/L")

  paper_specific_compartments <- c(
    "plt1", "plt2", "plt3", "plt4", "plt5",
    "plt6", "plt7", "plt8", "plt9", "plt10"
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Wang 2010 parameterizes PK in per-kg form (Vc reported in L/kg). Absolute central volume scales linearly: vc = exp(lvc) * WT. Dose amounts (amt column) are absolute (ug), so per-kg dosing in the source (e.g. 1 ug/kg in a 70 kg subject) translates to amt = 70 ug. No allometric exponent or further covariate effect is applied; the paper had no covariate analysis (mean-data fit, 32 healthy adults).",
      source_name        = "WT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 32L,
    n_studies      = 1L,
    age_range      = "18-50 years",
    sex_female_pct = NA,
    race_ethnicity = "Not reported in detail; described as demographically similar across cohorts (Wang 2010 Results / Subjects).",
    disease_state  = "Healthy non-obese volunteers (BMI < 30 kg/m^2). Eligible women surgically sterilized or postmenopausal.",
    dose_range     = "Single dose. IV: 0.3, 1, 10 ug/kg (4 subjects per group). SC: 0.1, 0.3, 1, 2 ug/kg (4 subjects per group except 2 ug/kg SC = 8). 16 placebo subjects across groups.",
    bmi_max        = "< 30 kg/m^2",
    samples        = "123 quantifiable serum concentrations in 17 participants (only IV cohorts and 2 ug/kg SC; lower SC doses were below the 0.018 ng/mL ELISA LOQ); 481 platelet counts from 32 participants.",
    regions        = "Not specified; phase I healthy-volunteer study, sponsored by Amgen.",
    notes          = "Wang 2010 fit MEAN concentration and MEAN platelet count by dose cohort (not individual data), so the published parameter table (Table I) is a typical-value fit without IIV. Baseline platelet counts were fixed at the predose mean for each dose level (Wang 2010 Statistical Model). lrbase in this model defaults to 250 x 10^9/L (a representative healthy baseline matching the Discussion's 'baseline platelet counts of 150-450 x 10^9 cells/L'); the vignette shows how to vary it. PK sampling at 0.083, 0.25, 0.5, 1, 2, 4, 8, 12 h and 1, 1.5, 2, 4, 7, 9, 10, 11, 12, 13, 14, 16, 19 days post-dose. Platelets at predose and days 2, 3, 5, 8, 10, 11, 12, 13, 14, 15, 17, 20, 28, 42."
  )

  ini({
    # ---------------------------------------------------------------------
    # PK structural parameters - Wang 2010 Table I.
    # Vc is reported per kg (L/kg); the model multiplies by WT in model().
    # Inter-compartmental rates kcp and kpc are reported directly (the
    # paper does not factor a Q / Vp); the model keeps that parameterisation
    # via lkcp / lkpc (precedent: Benkali_2010_tacrolimus.R).
    # ---------------------------------------------------------------------
    lka    <- log(0.0254);  label("First-order SC absorption rate ka (1/h)")                               # Wang 2010 Table I: ka = 0.0254 1/h (RSE 20%)
    lvc    <- log(0.0683);  label("Apparent central volume Vc per body weight (L/kg)")                     # Wang 2010 Table I: Vc = 0.0683 L/kg (RSE 21%)
    lkel   <- log(0.0382);  label("First-order linear elimination rate kel (1/h)")                         # Wang 2010 Table I: kel = 0.0382 1/h (RSE 48%)
    lkcp   <- log(0.0806);  label("Central-to-peripheral first-order rate kcp (1/h)")                      # Wang 2010 Table I: kcp = 0.0806 1/h (RSE 11%)
    lkpc   <- log(0.0148);  label("Peripheral-to-central first-order rate kpc (1/h)")                      # Wang 2010 Table I: kpc = 0.0148 1/h (RSE 27%)
    lfdepot <- log(0.499);  label("Absolute SC bioavailability F (unitless)")                              # Wang 2010 Table I: F = 0.499 (RSE 47%)

    # ---------------------------------------------------------------------
    # TMDD / QE-binding parameters - Wang 2010 Table I.
    # Rtot is computed as xi * PLT (with PLT in 10^9 cells/L and xi in
    # fg/platelet, the unit conversions give Rtot in ng/mL; see vignette
    # dimensional analysis).
    # ---------------------------------------------------------------------
    lkint  <- log(0.173);   label("Drug-receptor complex internalisation rate kint (1/h)")                 # Wang 2010 Table I: kint = 0.173 1/h (RSE 29%)
    lxi    <- log(0.0215);  label("c-Mpl receptor density per platelet, romiplostim weight equivalents (fg/platelet)") # Wang 2010 Table I: xi = 0.0215 fg/platelet (RSE 31%)
    lkd    <- log(0.131);   label("Equilibrium dissociation constant KD for romiplostim-c-Mpl binding (ng/mL)") # Wang 2010 Table I: KD = 0.131 ng/mL (RSE 129%)

    # ---------------------------------------------------------------------
    # PD - Hill stimulation of precursor production - Wang 2010 Table I.
    # Drives a (1 + Smax * Cc / (SC50 + Cc)) multiplier on kin (Eq 10).
    # ---------------------------------------------------------------------
    lec50  <- log(0.0520);  label("Romiplostim concentration for half-maximal precursor stimulation SC50 (ng/mL)") # Wang 2010 Table I: SC50 = 0.0520 ng/mL (RSE 16%)
    lemax  <- log(11.2);    label("Maximum stimulation Smax of precursor production rate (unitless)")      # Wang 2010 Table I: Smax = 11.2 (RSE 7.6%)

    # ---------------------------------------------------------------------
    # PD - lifespan parameters - Wang 2010 Table I.
    # NP = NPLT = 10 (Wang 2010 Methods, "Pharmacodynamic Model"); held
    # as numeric constants in model() rather than as estimable thetas.
    # ---------------------------------------------------------------------
    ltp    <- log(142);     label("Mean megakaryocyte precursor-cell lifespan TP (h)")                     # Wang 2010 Table I: TP = 142 h (~5.9 days; RSE 5.1%)
    ltplt  <- log(253);     label("Mean platelet lifespan TPLT (h)")                                       # Wang 2010 Table I: TPLT = 253 h (~10.5 days; RSE 1.6%)

    # ---------------------------------------------------------------------
    # Baseline platelet count (lrbase) - the source paper fixed PLT0 to
    # each cohort's predose mean (Wang 2010 Statistical Model). The default
    # here (250 x 10^9/L) is a representative healthy baseline within the
    # 150-450 x 10^9/L range cited in the Discussion (Wang 2010 Results).
    # Users simulating against a specific cohort should override lrbase to
    # that cohort's reported predose mean.
    # ---------------------------------------------------------------------
    lrbase <- log(250);     label("Baseline platelet count PLT0 (10^9 cells/L), drug-free steady state")   # Wang 2010 Discussion / Statistical Model; default mid-range healthy

    # ---------------------------------------------------------------------
    # Residual error - Wang 2010 Table I reports sigma^2 (variance) for a
    # proportional error model (Eq 18: Y_obs = Y_pred * (1 + eps)).
    # propSd = sqrt(sigma^2) is what nlmixr2 expects for prop() error.
    # ---------------------------------------------------------------------
    propSd     <- 0.3975;   label("Proportional residual error on Cc (fraction); sqrt(sigma^2_PK = 0.158)") # Wang 2010 Table I: sigma^2-PK = 0.158 (RSE 60%) -> sqrt = 0.3975
    propSd_PLT <- 0.0923;   label("Proportional residual error on PLT (fraction); sqrt(sigma^2_PD = 0.00851)") # Wang 2010 Table I: sigma^2-PD = 0.00851 (RSE 24%) -> sqrt = 0.0923
  })

  model({
    # ---------------------------------------------------------------------
    # Individual parameters (no IIV - Wang 2010 fit MEAN data; the
    # population approach failed for this complex PDMDD + lifespan model).
    # ---------------------------------------------------------------------
    ka    <- exp(lka)
    vc    <- exp(lvc) * WT          # L (per-kg parameter * body weight)
    kel   <- exp(lkel)
    kcp   <- exp(lkcp)
    kpc   <- exp(lkpc)
    fdepot <- exp(lfdepot)

    kint  <- exp(lkint)
    xi    <- exp(lxi)
    kd    <- exp(lkd)

    ec50  <- exp(lec50)
    emax  <- exp(lemax)
    tp    <- exp(ltp)
    tplt  <- exp(ltplt)
    rbase <- exp(lrbase)

    # ---------------------------------------------------------------------
    # Lifespan-chain constants - NP and NPLT compartments per the paper
    # (Wang 2010 Pharmacodynamic Model: "A cascade of NP = 10 age
    # compartments ... NPLT = 10 age compartments").
    # ---------------------------------------------------------------------
    np_per_tp     <- 10 / tp     # transit rate constant in the precursor chain (1/h)
    nplt_per_tplt <- 10 / tplt   # transit rate constant in the platelet chain (1/h)

    # ---------------------------------------------------------------------
    # Total circulating platelet count (sum of the 10 age compartments,
    # Wang 2010 Eq 14). Receptor concentration Rtot = xi * PLT (the unit
    # conversions cancel: fg/platelet * 10^9 platelets/L * 10^-6 ng/fg *
    # 10^-3 L/mL = ng/mL).
    # ---------------------------------------------------------------------
    PLT  <- plt1 + plt2 + plt3 + plt4 + plt5 + plt6 + plt7 + plt8 + plt9 + plt10
    rtot <- xi * PLT

    # ---------------------------------------------------------------------
    # Quasi-equilibrium TMDD relations (Wang 2010 Eq 5).
    # central holds the TOTAL romiplostim amount (free + bound DR) in
    # the central compartment; ctot is total drug concentration.
    # cfree solves cfree * (rtot - complex) = kd * complex with
    # ctot = cfree + complex, i.e. the standard Mager-Krzyzanski
    # quadratic; complex = ctot - cfree.
    # ---------------------------------------------------------------------
    ctot    <- central / vc
    qe_disc <- ctot - rtot - kd
    cfree   <- 0.5 * (qe_disc + sqrt(qe_disc * qe_disc + 4 * kd * ctot))
    complex <- ctot - cfree

    # ---------------------------------------------------------------------
    # PK ODEs (Wang 2010 Eqs 2-4; central tracks Atot = Ac + DR).
    # Distribution and linear elimination act on the FREE drug only
    # (kel * Ac, kcp * Ac, kpc * Ap); kint acts on the bound complex
    # (kint * DR). Ac = cfree * vc; DR = complex * vc; Ap = peripheral1.
    # ---------------------------------------------------------------------
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * cfree * vc -
                          kcp * cfree * vc + kpc * peripheral1 -
                          kint * complex * vc
    d/dt(peripheral1) <-  kcp * cfree * vc - kpc * peripheral1

    # SC bioavailability applied to the depot (IV doses enter `central`
    # directly with f = 1 by default).
    f(depot) <- fdepot

    # ---------------------------------------------------------------------
    # PD Hill stimulation of precursor production (Wang 2010 Eq 9).
    # stim drives kin -> kin * (1 + stim).
    # ---------------------------------------------------------------------
    stim <- emax * cfree / (ec50 + cfree)

    # ---------------------------------------------------------------------
    # Megakaryocyte precursor chain (Wang 2010 Eqs 10-11; NP = 10).
    # kin = rbase / tplt is the baseline precursor production rate
    # (Wang 2010 Eq 17). At steady state each precursor compartment
    # holds rbase * tp / (NP * tplt) cells/L; each platelet compartment
    # holds rbase / NPLT cells/L. The catenary chain transfer-rate
    # constants are NP/TP for precursors and NPLT/TPLT for platelets.
    # ---------------------------------------------------------------------
    kin <- rbase / tplt

    d/dt(precursor1)  <- kin * (1 + stim) - np_per_tp * precursor1
    d/dt(precursor2)  <- np_per_tp * precursor1 - np_per_tp * precursor2
    d/dt(precursor3)  <- np_per_tp * precursor2 - np_per_tp * precursor3
    d/dt(precursor4)  <- np_per_tp * precursor3 - np_per_tp * precursor4
    d/dt(precursor5)  <- np_per_tp * precursor4 - np_per_tp * precursor5
    d/dt(precursor6)  <- np_per_tp * precursor5 - np_per_tp * precursor6
    d/dt(precursor7)  <- np_per_tp * precursor6 - np_per_tp * precursor7
    d/dt(precursor8)  <- np_per_tp * precursor7 - np_per_tp * precursor8
    d/dt(precursor9)  <- np_per_tp * precursor8 - np_per_tp * precursor9
    d/dt(precursor10) <- np_per_tp * precursor9 - np_per_tp * precursor10

    # ---------------------------------------------------------------------
    # Platelet aging chain (Wang 2010 Eqs 12-13; NPLT = 10).
    # plt1 receives from precursor10 at rate np_per_tp; subsequent
    # platelet compartments transit at rate nplt_per_tplt.
    # ---------------------------------------------------------------------
    d/dt(plt1)  <- np_per_tp * precursor10 - nplt_per_tplt * plt1
    d/dt(plt2)  <- nplt_per_tplt * plt1  - nplt_per_tplt * plt2
    d/dt(plt3)  <- nplt_per_tplt * plt2  - nplt_per_tplt * plt3
    d/dt(plt4)  <- nplt_per_tplt * plt3  - nplt_per_tplt * plt4
    d/dt(plt5)  <- nplt_per_tplt * plt4  - nplt_per_tplt * plt5
    d/dt(plt6)  <- nplt_per_tplt * plt5  - nplt_per_tplt * plt6
    d/dt(plt7)  <- nplt_per_tplt * plt6  - nplt_per_tplt * plt7
    d/dt(plt8)  <- nplt_per_tplt * plt7  - nplt_per_tplt * plt8
    d/dt(plt9)  <- nplt_per_tplt * plt8  - nplt_per_tplt * plt9
    d/dt(plt10) <- nplt_per_tplt * plt9  - nplt_per_tplt * plt10

    # ---------------------------------------------------------------------
    # Initial conditions (drug-free steady state; Wang 2010 Eqs 7, 15, 16).
    # central / peripheral1 / depot default to 0 (no drug at t = 0).
    # ---------------------------------------------------------------------
    precursor1(0)  <- rbase * tp / (10 * tplt)
    precursor2(0)  <- rbase * tp / (10 * tplt)
    precursor3(0)  <- rbase * tp / (10 * tplt)
    precursor4(0)  <- rbase * tp / (10 * tplt)
    precursor5(0)  <- rbase * tp / (10 * tplt)
    precursor6(0)  <- rbase * tp / (10 * tplt)
    precursor7(0)  <- rbase * tp / (10 * tplt)
    precursor8(0)  <- rbase * tp / (10 * tplt)
    precursor9(0)  <- rbase * tp / (10 * tplt)
    precursor10(0) <- rbase * tp / (10 * tplt)

    plt1(0)  <- rbase / 10
    plt2(0)  <- rbase / 10
    plt3(0)  <- rbase / 10
    plt4(0)  <- rbase / 10
    plt5(0)  <- rbase / 10
    plt6(0)  <- rbase / 10
    plt7(0)  <- rbase / 10
    plt8(0)  <- rbase / 10
    plt9(0)  <- rbase / 10
    plt10(0) <- rbase / 10

    # ---------------------------------------------------------------------
    # Observations: Cc (free romiplostim serum concentration, ng/mL) and
    # PLT (total circulating platelet count, 10^9 cells/L). Both follow
    # proportional error (Wang 2010 Eq 18).
    # ---------------------------------------------------------------------
    Cc <- cfree
    Cc  ~ prop(propSd)
    PLT ~ prop(propSd_PLT)
  })
}
