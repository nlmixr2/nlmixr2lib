Yu_2022_ofatumumab <- function() {
  description <- "Population PK / B-cell-count model for subcutaneous ofatumumab in adults with relapsing multiple sclerosis (Yu 2022)"
  reference <- "Yu H, Graham G, David OJ, Kahn JM, Savelieva M, Pigeolet E, Das Gupta A, Pingili R, Willi R, Ramanathan K, Kieseier BC, Hach T, Aslanis V, Bagger Y, Ravenstijn P. Population Pharmacokinetic-B Cell Modeling for Ofatumumab in Patients with Relapsing Multiple Sclerosis. CNS Drugs. 2022;36(3):283-300. doi:10.1007/s40263-021-00895-w"
  vignette <- "Yu_2022_ofatumumab"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  # Final-model parameter estimates from Yu 2022 Table 3.
  # PK: quasi-steady-state TMDD with first-order SC absorption, two drug
  # compartments (central + peripheral), and a CD20 receptor compartment with
  # time-decaying synthesis rate ksyn(t) = ksyn_inf + (ksyn_0 - ksyn_inf) *
  # exp(-kdes * t / 365.25).
  # PD: indirect-response B-cell-count model with central + peripheral B cell
  # compartments and a sigmoid Emax stimulation of B-cell lysis driven by free
  # ofatumumab concentration.
  #
  # The QSS receptor parameters R0, ksyn0, ksyn_inf, KD are reported in nmol/L
  # in Yu 2022 Table 3. They are converted to drug-equivalent mg/L below using
  # the ofatumumab molecular weight (~149 kDa, IgG1) so that the QSS algebra
  # runs entirely on a single mg/L scale consistent with mg dosing and L
  # volumes (Lc = central / vc is natively mg/L). EC50 is already in mg/L.
  covariateData <- list(
    WT = list(
      description        = "Body weight (baseline; time-fixed)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on ka, Vc, ksyn0, CL, B0, kout per Yu 2022 covariate equations P_i = P_TV * (WT_i/70)^beta. Reference 70 kg (median weight in the pooled five-study cohort, Table 2).",
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Baseline age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on B0 per Yu 2022 covariate equation B0_i = B0 * (AGE_i/38)^beta. Reference 38 years (median age in the pooled five-study cohort, Table 2).",
      source_name        = "Age"
    ),
    BLBCELL = list(
      description        = "Baseline CD19+ B cell count (cells/uL) measured by FACS prior to first dose",
      units              = "cells/uL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on Emax per Yu 2022 covariate equation Emax_i = Emax * (Bcell0_i/200)^beta. Reference 200 cells/uL (median baseline CD19+ B cell count in the pooled cohort, Table 2).",
      source_name        = "Bcell0"
    ),
    ROUTE_IV = list(
      description        = "Intravenous administration indicator (1 = IV, 0 = SC)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (SC)",
      notes              = "Per-subject route covariate carrying paper's [Admin route = IV] effect on R0, CL, Q, ksyn_inf. Reference is SC. For simulation, set ROUTE_IV = 1 for IV cohorts (dose into central) and 0 for SC cohorts (dose into depot).",
      source_name        = "Admin route = IV"
    ),
    DEVICE_AI = list(
      description        = "SC autoinjector device indicator (1 = autoinjector, 0 = prefilled syringe)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (PFS)",
      notes              = "Per-subject device covariate carrying paper's [Formulation = AI] effect on k_e(P) and R0. Reference is PFS. Set to 0 for IV subjects since device is undefined for IV; IV effects are carried by ROUTE_IV.",
      source_name        = "Formulation = AI"
    ),
    STUDY_APLIOS = list(
      description        = "APLIOS (NCT03560739) study indicator (1 = APLIOS, 0 = other)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-APLIOS)",
      notes              = "Per-subject categorical effect on Emax per Yu 2022 covariate equation Emax_i = Emax * exp(beta * [Study_i = APLIOS]).",
      source_name        = "Study = APLIOS"
    ),
    STUDY_MIRROR = list(
      description        = "MIRROR (NCT01457924) study indicator (1 = MIRROR, 0 = other)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-MIRROR)",
      notes              = "Per-subject categorical effect on kout per Yu 2022 covariate equation kout_i = kout * exp(beta * [Study_i = MIRROR]).",
      source_name        = "Study = MIRROR"
    )
  )

  population <- list(
    n_subjects     = 1486L,
    n_studies      = 5L,
    age_range      = "18-56 years",
    age_median     = "38 years",
    weight_range   = "40.5-171.6 kg",
    weight_median  = "70.0 kg",
    sex_female_pct = 67.8,
    race_ethnicity = c(White = 91.1, Black = 2.4, Asian = 2.5, Other_AmInd_AlaskaNative = 3.2, Unknown = 0.7),
    disease_state  = "Adults with relapsing forms of multiple sclerosis (RMS), including relapsing-remitting MS (RRMS); EDSS 0-5.5 at screening; mean EDSS 2.8 (range 0-6).",
    dose_range     = "OMS115102: 100/300/700 mg IV at weeks 0+2 or 24+26 (700 mg arm excluded from PK analysis). MIRROR: 0/3/30/60 mg SC q12w or 60 mg SC q4w. APLIOS and ASCLEPIOS I/II: 20 mg SC q4w after three weekly 20-mg loading doses on days 1, 7, 14 (with the AI device used in 141 APLIOS patients; PFS otherwise).",
    regions        = "Multi-regional phase 2/3 programme (OMS115102, MIRROR, APLIOS, ASCLEPIOS I, ASCLEPIOS II).",
    studies        = "OMS115102 (N=25 IV RRMS), MIRROR (N=231 SC PFS RRMS), APLIOS (N=284 SC AI/PFS RMS), ASCLEPIOS I (N=465 SC PFS RMS), ASCLEPIOS II (N=481 SC PFS RMS); patients treated with ofatumumab and with B cell data who were included in the final model.",
    bcell_baseline_median = "200 cells/uL CD19+ (range 0-1520; pooled five-study cohort, Table 2)",
    edss_median    = "2.5 (range 0-6, baseline)",
    notes          = "Pooled phase 2/3 dataset assembled for the population PK-B cell analysis. PK dataset: 9168 plasma concentrations from 1440 patients (placebo and 700-mg cohorts excluded). PD dataset: 17,158 CD19+ B cell counts from 1486 patients. ASCLEPIOS B cell counts modelled as interval-censored (LLOQ 0/5/15/25 cells/uL bins). Anti-drug antibody incidence < 2 % with no detectable PK or B cell impact (Yu 2022 Results 'Data')."
  )

  ini({
    # ---- Structural PK (Yu 2022 Table 3 final-model estimates) ---------------
    lka      <- log(0.157);    label("Population first-order SC absorption rate ka (1/day)")                                  # Yu 2022 Table 3
    logitfdepot <- log(0.685 / (1 - 0.685)); label("Logit population SC bioavailability F (typical F = 0.685)")               # Yu 2022 Table 3
    lvc      <- log(2.62);     label("Population central volume of distribution Vc (L)")                                      # Yu 2022 Table 3
    lkep     <- log(1.31);     label("Population complex internalisation rate k_e(P) (1/day)")                                # Yu 2022 Table 3
    lcl      <- log(0.34);     label("Population free-drug clearance CL (L/day)")                                             # Yu 2022 Table 3
    lq       <- log(0.358);    label("Population intercompartmental clearance Q (L/day)")                                     # Yu 2022 Table 3
    lvp      <- fixed(log(2.8));        label("Peripheral volume of distribution Vp (L; FIXED per Ryman & Meibohm 2017)")     # Yu 2022 Table 3 (fixed)

    # ---- TMDD-QSS receptor and binding parameters ----------------------------
    # R0, ksyn0, ksyn_inf, and KD are reported in nmol/L (or nmol/L/day) in
    # Yu 2022 Table 3. Convert to drug-equivalent mg/L using ofatumumab MW
    # 149 kDa: c[mg/L] = c[nmol/L] * 149/1000 = c[nmol/L] * 0.149.
    #   R0       [mg/L]      = 32.5  * 0.149 = 4.8425
    #   ksyn0    [mg/L/day]  = 0.985 * 0.149 = 0.146765
    #   ksyn_inf [mg/L/day]  = 0.0554 * 0.149 = 0.0082546
    #   KD       [mg/L]      = 0.167 * 0.149 = 0.024883
    # k_off and kdes are 1/time and unit-invariant.
    lr0      <- log(4.8425);        label("Population baseline total CD20 receptor amount R0 in drug-equivalent mg/L (paper: 32.5 nmol/L)")  # Yu 2022 Table 3
    lksyn0   <- log(0.146765);      label("Population CD20 synthesis rate at t=0 ksyn0 in drug-equivalent mg/L/day (paper: 0.985 nmol/L/day)") # Yu 2022 Table 3
    lksyninf <- log(0.0082546);     label("Population CD20 synthesis rate at t=infinity ksyn_inf in drug-equivalent mg/L/day (paper: 0.0554 nmol/L/day)") # Yu 2022 Table 3
    lkd      <- fixed(log(0.024883));   label("Equilibrium dissociation constant KD in drug-equivalent mg/L (paper: 0.167 nmol/L; FIXED)")    # Yu 2022 Table 3 (fixed)
    lkoff    <- fixed(log(5.53));       label("Drug-target dissociation rate constant k_off (1/day; FIXED preclinical value)")                # Yu 2022 Table 3 (fixed)
    lkdes    <- log(2.58);          label("Time rate constant on receptor synthesis decay kdes (1/year)")                                     # Yu 2022 Table 3

    # ---- B cell PD parameters (Yu 2022 Table 3 final-model estimates) --------
    lb0      <- log(194);      label("Population baseline central-compartment B cell count B0 (cells/uL)")                  # Yu 2022 Table 3
    lemax    <- log(159);      label("Population maximum B-cell-lysis stimulatory effect Emax (unitless)")                  # Yu 2022 Table 3
    lec50    <- log(0.0057);   label("Population free-drug concentration producing 50% of Emax EC50 (mg/L)")                # Yu 2022 Table 3
    lhill   <- log(2.81);     label("Hill / sigmoidicity parameter for the lysis stimulatory function (unitless)")         # Yu 2022 Table 3
    lkout    <- log(0.0124);   label("Population B cell elimination rate kout (1/day)")                                     # Yu 2022 Table 3
    lqb      <- fixed(log(0.78));     label("Inter-B-cell-compartment flow QB (L/day; FIXED)")                              # Yu 2022 Table 3 (fixed; no IIV reported)
    lvb      <- log(3.7);      label("Peripheral B cell compartment volume Vb (L)")                                         # Yu 2022 Table 3

    # ---- Covariate effects (Yu 2022 Table 3 + covariate equations) -----------
    e_wt_ka       <- -0.457; label("Power exponent of (WT/70) on ka (unitless)")                                            # Yu 2022 Table 3
    e_wt_vc       <-  1.2;   label("Power exponent of (WT/70) on Vc (unitless)")                                            # Yu 2022 Table 3
    e_wt_ksyn0    <- -1.52;  label("Power exponent of (WT/70) on ksyn0 (unitless)")                                         # Yu 2022 Table 3
    e_wt_cl       <-  1.52;  label("Power exponent of (WT/70) on CL (unitless)")                                            # Yu 2022 Table 3
    e_wt_b0       <-  0.271; label("Power exponent of (WT/70) on B0 (unitless)")                                            # Yu 2022 Table 3
    e_wt_kout     <- -0.624; label("Power exponent of (WT/70) on kout (unitless)")                                          # Yu 2022 Table 3
    e_age_b0      <- -0.282; label("Power exponent of (AGE/38) on B0 (unitless)")                                           # Yu 2022 Table 3
    e_blbc_emax   <-  0.275; label("Power exponent of (BLBCELL/200) on Emax (unitless)")                                    # Yu 2022 Table 3
    e_iv_r0       <-  0.987; label("Exponential effect of IV administration route on R0 (unitless)")                        # Yu 2022 Table 3
    e_iv_cl       <- -1.07;  label("Exponential effect of IV administration route on CL (unitless)")                        # Yu 2022 Table 3
    e_iv_q        <- -2.31;  label("Exponential effect of IV administration route on Q (unitless)")                         # Yu 2022 Table 3
    e_iv_ksyninf  <-  2.49;  label("Exponential effect of IV administration route on ksyn_inf (unitless)")                  # Yu 2022 Table 3
    e_ai_kep      <-  0.713; label("Exponential effect of AI device on k_e(P) (unitless)")                                  # Yu 2022 Table 3
    e_ai_r0       <- -0.544; label("Exponential effect of AI device on R0 (unitless)")                                      # Yu 2022 Table 3
    e_aplios_emax <-  0.503; label("Exponential effect of APLIOS study on Emax (unitless)")                                 # Yu 2022 Table 3
    e_mirror_kout <- -0.554; label("Exponential effect of MIRROR study on kout (unitless)")                                 # Yu 2022 Table 3

    # ---- Inter-individual variability ---------------------------------------
    # Yu 2022 reports IIV as the SD of the log-normal random effect (or logit
    # for F). The nlmixr2 `etaXxx ~ value` syntax stores the variance, so SDs
    # are squared here. Three correlated blocks are reported in the paper text
    # (Methods "PKPD Modeling and Simulation", "Final Model Description"):
    #   Block 1 (PK): CL, ka, kdes
    #   Block 2 (PK): k_e(P), R0, ksyn_inf
    #   Block 3 (PD): Emax, kout, Vb
    # Off-diagonal covariances are computed from the reported correlations and
    # SDs as cov_ij = rho_ij * sd_i * sd_j.
    #
    # Block 3 typo handling: Yu 2022 Table 3 lists "Corr_kout_Vb" twice
    # (-0.336 and 0.423) and is missing a "Corr_kout_Emax" entry. Since the
    # block is over three parameters (Vb, kout, Emax), three pairwise
    # correlations are expected (Vb-Emax, kout-Vb, kout-Emax). The most
    # parsimonious interpretation is that the second 0.423 entry is a typo
    # for Corr_kout_Emax. The implied correlation matrix is positive
    # semidefinite (det ~= 0.855), supporting this reading. See the vignette's
    # "Assumptions and deviations" section.

    # Block 1 covariances (CL, ka, kdes):
    #   sd_ka = 0.652, sd_cl = 0.486, sd_kdes = 0.654
    #   cov(ka,  cl)   = -0.294 * 0.652 * 0.486 = -0.09316
    #   cov(ka,  kdes) =  0.433 * 0.652 * 0.654 =  0.18464
    #   cov(cl,  kdes) =  0.642 * 0.486 * 0.654 =  0.20407
    etalka + etalcl + etalkdes ~ c(0.652^2,
                                   -0.09316, 0.486^2,
                                    0.18464, 0.20407, 0.654^2)  # Yu 2022 Table 3 (Block 1: CL/ka/kdes)

    # Block 2 covariances (k_e(P), R0, ksyn_inf):
    #   sd_kep = 1.14, sd_r0 = 0.91, sd_ksyninf = 2.16
    #   cov(kep,     r0)        = -0.551 * 1.14 * 0.91 = -0.57164
    #   cov(kep,     ksyninf)   = -0.464 * 1.14 * 2.16 = -1.14266
    #   cov(r0,      ksyninf)   =  0.470 * 0.91 * 2.16 =  0.92395
    etalkep + etalr0 + etalksyninf ~ c(1.14^2,
                                       -0.57164, 0.91^2,
                                       -1.14266, 0.92395, 2.16^2)  # Yu 2022 Table 3 (Block 2: k_e(P)/R0/ksyn_inf)

    # Block 3 covariances (Vb, kout, Emax) -- typo-corrected:
    #   sd_vb = 1.37, sd_kout = 0.922, sd_emax = 0.587
    #   cov(vb,   emax) =  0.280 * 1.37  * 0.587 =  0.22513
    #   cov(vb,   kout) = -0.336 * 1.37  * 0.922 = -0.42446
    #   cov(kout, emax) =  0.423 * 0.922 * 0.587 =  0.22893
    etalvb + etalkout + etalemax ~ c(1.37^2,
                                     -0.42446, 0.922^2,
                                      0.22513, 0.22893, 0.587^2)  # Yu 2022 Table 3 (Block 3: Vb/kout/Emax; typo-corrected per vignette)

    etalogitfdepot ~ 0.531^2  # Yu 2022 Table 3 (logit-scale eta on F)
    etalvc         ~ 0.116^2  # Yu 2022 Table 3
    etalksyn0      ~ 0.0559^2 # Yu 2022 Table 3 (very small SD; high RSE 64.1%)
    etalq          ~ 0.705^2  # Yu 2022 Table 3
    etalb0         ~ 0.394^2  # Yu 2022 Table 3
    etalec50       ~ 0.927^2  # Yu 2022 Table 3
    etalhill      ~ 1.46^2   # Yu 2022 Table 3

    # ---- Residual error (Yu 2022 Table 3 final-model estimates) -------------
    propSd       <- 0.278;    label("Proportional residual error for ofatumumab plasma concentration (fraction)")  # Yu 2022 Table 3
    addSd        <- 0.0316;   label("Additive residual error for ofatumumab plasma concentration (mg/L)")          # Yu 2022 Table 3
    propSd_Bcell    <- 0.381;    label("Proportional residual error for B cell count (fraction)")                     # Yu 2022 Table 3
    addSd_Bcell     <- 0.153;    label("Additive residual error for B cell count (cells/uL)")                         # Yu 2022 Table 3
  })

  model({
    # ---- Individual structural parameters ----------------------------------
    # Covariate forms per Yu 2022 covariate equations (Methods, Final Model
    # Description). Continuous covariates enter as power scaling on the median
    # reference; categorical covariates enter as exp(beta * [indicator]).
    ka       <- exp(lka       + etalka)       * (WT / 70)^e_wt_ka
    logit_f  <- logitfdepot + etalogitfdepot
    fdepot   <- 1 / (1 + exp(-logit_f))
    vc       <- exp(lvc       + etalvc)       * (WT / 70)^e_wt_vc
    kep      <- exp(lkep      + etalkep)      * exp(e_ai_kep      * DEVICE_AI)
    cl       <- exp(lcl       + etalcl)       * (WT / 70)^e_wt_cl    * exp(e_iv_cl       * ROUTE_IV)
    q        <- exp(lq        + etalq)                               * exp(e_iv_q        * ROUTE_IV)
    vp       <- exp(lvp)
    r0       <- exp(lr0       + etalr0)       * exp(e_ai_r0        * DEVICE_AI) * exp(e_iv_r0       * ROUTE_IV)
    ksyn0    <- exp(lksyn0    + etalksyn0)    * (WT / 70)^e_wt_ksyn0
    ksyninf  <- exp(lksyninf  + etalksyninf)                         * exp(e_iv_ksyninf  * ROUTE_IV)
    kd       <- exp(lkd)
    koff     <- exp(lkoff)
    kdes     <- exp(lkdes     + etalkdes)
    b0_ind   <- exp(lb0       + etalb0)       * (WT / 70)^e_wt_b0    * (AGE / 38)^e_age_b0
    emax     <- exp(lemax     + etalemax)     * (BLBCELL / 200)^e_blbc_emax * exp(e_aplios_emax * STUDY_APLIOS)
    ec50     <- exp(lec50     + etalec50)
    hill    <- exp(lhill    + etalhill)
    kout     <- exp(lkout     + etalkout)     * (WT / 70)^e_wt_kout  * exp(e_mirror_kout * STUDY_MIRROR)
    qb       <- exp(lqb)
    vb       <- exp(lvb       + etalvb)

    # ---- TMDD-QSS algebra (Yu 2022 Eq. 'Lc = ...') --------------------------
    # All concentrations on a single drug-equivalent mg/L scale. Receptor
    # parameters lr0, lksyn0, lksyninf, lkd above are pre-converted from
    # nmol/L using MW 149 kDa.
    #   k_on  [(mg/L)^-1 day^-1] = k_off / KD
    #   Ks    [mg/L]             = (k_e(P) + k_off) / k_on
    #   k_deg [1/day]            = ksyn(t) / R0  (so dRtot/dt = 0 at t=0)
    kon      <- koff / kd
    ks       <- (kep + koff) / kon

    # Time-varying receptor synthesis (Yu 2022 'ksyn(t) = ksyn_inf + ...').
    # t is in days; kdes is in 1/year (paper uses /365.25).
    ksyn_t   <- ksyninf + (ksyn0 - ksyninf) * exp(-kdes * t / 365.25)
    kdeg_t   <- ksyn_t / r0

    # Total drug concentration in central from the 'central' state (drug
    # amount, mg) divided by Vc; peripheral drug concentration analogous.
    ltot     <- central / vc
    lp       <- peripheral1 / vp
    rtot     <- total_target  # state stored as concentration (drug-equivalent mg/L)

    # QSS quadratic for free drug Lc (Eq.: Lc = 0.5*(Ltot - Rtot - Ks +
    # sqrt((Ltot - Rtot - Ks)^2 + 4 Ltot Ks))).
    qsscore  <- ltot - rtot - ks
    lc       <- 0.5 * (qsscore + sqrt(qsscore * qsscore + 4 * ltot * ks))

    # Bound concentration (drug-equivalent mg/L), needed for the central ODE
    # internalisation term and for the receptor ODE.
    bound    <- ltot - lc

    # ---- ODE system --------------------------------------------------------
    # Standard mass-conserving 2-compartment QSS-TMDD form. Yu 2022's printed
    # equations parameterise kpc = Q/Vp in the central ODE (which would not
    # conserve mass when Vp != Vc); the standard form (Q/Vc for Lp's
    # contribution to the central concentration) is used here, matching what
    # Monolix's TMDD library evaluates internally. See vignette
    # "Assumptions and deviations".
    #
    # Drug central total amount (mg/day):
    #   ka * Adepot                          : SC absorption (after F applied via f(depot))
    #   - cl * lc                            : linear free-drug elimination
    #   - q  * (lc - lp)                     : 2-compartment free-drug exchange
    #   - kep * vc * bound                   : TMDD complex internalisation
    d/dt(depot)        <- -ka * depot
    d/dt(central)      <-  ka * depot - cl * lc - q * (lc - lp) - kep * vc * bound
    d/dt(peripheral1)  <-  q  * (lc - lp)
    # Total receptor (drug-equivalent mg/L) per Yu 2022 Eq. 'dRtot/dt = ...'
    d/dt(total_target) <-  ksyn_t - kdeg_t * rtot - (kep - kdeg_t) * bound

    # ---- B cell indirect-response model (Yu 2022 Eq. 'dB/dt = ...') --------
    # kin is fixed by the steady-state baseline B0:
    #   0 = kin - kout * B0  =>  kin = kout * B0
    # The 2-B-cell-compartment exchange terms cancel at steady state because
    # Bp(0) = B0 * Vb / Vc is the exact algebraic solution of dBp/dt = 0.
    kin      <- kout * b0_ind
    stim     <- emax * lc^hill / (ec50^hill + lc^hill)
    d/dt(bcell)         <- kin - kout * (1 + stim) * bcell - (qb / vc) * bcell + (qb / vb) * bcell_periph
    d/dt(bcell_periph)  <-                                    (qb / vc) * bcell - (qb / vb) * bcell_periph

    # ---- Initial conditions -----------------------------------------------
    # Drug compartments start empty; receptor at baseline R0; B cells at B0
    # and Bp(0) = B0 * Vb / Vc per Yu 2022 'B(0) = Bcell0, Bp(0) = ...'.
    total_target(0) <- r0
    bcell(0)        <- b0_ind
    bcell_periph(0) <- b0_ind * vb / vc

    # ---- Bioavailability ---------------------------------------------------
    # F applies to SC dosing into depot (logit-transformed in ini()). For IV
    # subjects, dose into 'central' directly; F = 1 by default for the central
    # compartment.
    f(depot) <- fdepot

    # ---- Observation variables and residual error -------------------------
    # Cc reported as total ofatumumab plasma concentration in mg/L (the assay
    # measures total drug per Methods 'B Cell-Depletion Target' and the QSS
    # model state Ltot). Bcell is the central-compartment CD19+ B cell count.
    Cc    <- ltot
    Bcell <- bcell
    Cc    ~ add(addSd)    + prop(propSd)
    Bcell ~ add(addSd_Bcell) + prop(propSd_Bcell)
  })
}
