Desmee_2015_PSA_mCRPC <- function() {
  description <- paste(
    "Mechanistic joint biomarker-survival model for prostate-specific antigen",
    "(PSA) kinetics under chemotherapy in metastatic castration-resistant",
    "prostate cancer (mCRPC). PSA is produced by a proliferating prostate-cell",
    "compartment and eliminated by a first-order process; chemotherapy blocks",
    "cell proliferation at a per-subject effectiveness eps until an individual",
    "escape time Tesc. The Weibull-baseline overall-survival hazard is",
    "log-linear in the current PSA value. This is a published-simulation-study",
    "model: parameter VALUES are pre-specified (inspired by one arm of the",
    "Tannock 2013 VENICE phase III trial), not estimated from real-data fits."
  )
  reference <- paste(
    "Desmee S, Mentre F, Veyrat-Follet C, Guedj J.",
    "Nonlinear mixed-effect models for prostate-specific antigen kinetics",
    "and link with survival in the context of metastatic prostate cancer:",
    "a comparison by simulation of two-stage and joint approaches.",
    "AAPS J. 2015 May;17(3):691-9. doi:10.1208/s12248-015-9745-5.",
    "Simulation parameter values inspired from one arm of the VENICE phase III",
    "trial (Tannock IF et al., Lancet Oncol. 2013;14:760-8;",
    "doi:10.1016/S1470-2045(13)70184-0).",
    sep = " "
  )
  vignette <- "Desmee_2015_PSA_mCRPC"
  paper_specific_compartments <- c("cells", "psa")

  units <- list(
    time          = "day",
    dosing        = "n/a (no explicit drug-dosing events; chemotherapy effect enters via the per-subject parameters eps and tesc)",
    concentration = "ng/mL (PSA) and 1/mL (cell density)"
  )

  covariateData <- list()
  # No covariates: in this simulation study the per-subject heterogeneity
  # is encoded entirely through eta on r, PSA0, eps (logit-scale), and Tesc.
  # The paper does not develop any covariate effects.

  population <- list(
    species        = "human (adult males with metastatic castration-resistant prostate cancer)",
    n_subjects     = 500L,
    n_studies      = 1L,
    age_range      = "not reported (the simulation parameters are inspired by the Tannock 2013 VENICE phase III arm; the simulation design itself does not specify a demographic distribution)",
    weight_range   = "not reported (not used; no allometric scaling in the model)",
    sex_female_pct = 0,
    disease_state  = "metastatic castration-resistant prostate cancer (mCRPC) under chemotherapy",
    dose_range     = "n/a (drug effect is encoded as a per-subject effectiveness eps held until an individual escape time tesc; there is no explicit dose or pharmacokinetic input)",
    regions        = "n/a (simulation study; no patient population was enrolled)",
    notes          = paste(
      "Simulation design (paper Methods, 'Simulation study' / 'Design'):",
      "M = 100 datasets of N = 500 patients each, with PSA measurements",
      "every 3 weeks for 2 years (last possible measurement at t = 735 day;",
      "max 36 observations per patient). Follow-up censored at t = 735 day;",
      "no dropout mechanism other than death. Four survival scenarios were",
      "simulated by varying beta and lambda jointly so that median end-of-",
      "study survival was approximately 25% under median PSA kinetics:",
      "No link (beta = 0, lambda = 580), Low link (beta = 0.005, lambda = 765),",
      "High link (beta = 0.02, lambda = 2150), and Short survival",
      "(beta = 0.02, lambda = 580). This model file uses the 'High link'",
      "scenario survival parameters as its canonical default; the vignette",
      "reproduces all four scenarios."
    )
  )

  ini({
    # ---------------------------------------------------------------
    # PSA-kinetics fixed effects (Table I; identical across all four
    # simulation scenarios because PSA dynamics do not depend on the
    # survival sub-model). All four are per-subject random parameters
    # with the listed transformation and omega.
    # ---------------------------------------------------------------
    lkg      <- log(0.05)             ; label("Prostate-cell proliferation rate r (1/day)")  # Table I: r = 0.05 day^-1, log-normal
    lpsa0    <- log(80)               ; label("Baseline PSA value PSA0 (ng/mL)")             # Table I: PSA0 = 80 ng/mL, log-normal
    logiteps <- log(0.3 / (1 - 0.3))  ; label("Logit of chemotherapy effectiveness epsilon (unitless on logit scale; back-transform in (0,1))")  # Table I: eps = 0.3, logit-normal
    ltesc    <- log(140)              ; label("Treatment-escape time Tesc (day)")            # Table I: Tesc = 140 day, log-normal

    # ---------------------------------------------------------------
    # Fixed mechanistic rate constants (paper Methods, "A Mechanistic
    # model for PSA kinetics"). Not estimated; held at literature
    # biological values, see refs [21] and [22].
    # ---------------------------------------------------------------
    ld_fix     <- fixed(log(0.046)) ; label("Cell death rate d (1/day) - fixed (Tu 1996; apoptotic-index 5% => half-life 15 days)")  # paper Methods, ref [21]
    ldelta_fix <- fixed(log(0.23))  ; label("PSA elimination rate delta (1/day) - fixed (Polascik 1999; PSA half-life ~3 days)")    # paper Methods, ref [22]

    # ---------------------------------------------------------------
    # Weibull baseline-hazard parameters and PSA-survival link (Table II).
    # This file uses the 'High link' scenario as the canonical default
    # (beta = 0.02, lambda = 2150, k = 1.5). The vignette swaps in the
    # other three scenarios for the multi-scenario figure replication.
    # ---------------------------------------------------------------
    llam_haz  <- log(2150) ; label("Weibull baseline-hazard scale parameter lambda (day) - 'High link' scenario")  # Table II: lambda = 2150 ('High link')
    lk_haz    <- log(1.5)  ; label("Weibull baseline-hazard shape parameter k (unitless)")                          # Table II: k = 1.5 (all scenarios)
    e_psa_haz <- 0.02      ; label("Hazard log-linear coefficient on current PSA, beta (1/(ng/mL)) - 'High link' scenario")  # Table II: beta = 0.02 ('High link')

    # ---------------------------------------------------------------
    # Inter-individual variability. The paper reports SDs (Table I,
    # right-most column); Omega is diagonal (paper Methods, "Statistical
    # model for PSA measurements"). Internal variance = SD^2.
    # ---------------------------------------------------------------
    etalkg      ~ 0.01    # Table I: omega = 0.1  on log(r)        => variance 0.01
    etalpsa0    ~ 0.36    # Table I: omega = 0.6  on log(PSA0)     => variance 0.36
    etalogiteps ~ 2.25    # Table I: omega = 1.5  on logit(eps)    => variance 2.25
    etaltesc    ~ 0.36    # Table I: omega = 0.6  on log(Tesc)     => variance 0.36

    # ---------------------------------------------------------------
    # Residual error on log(PSA + 1). Paper formula (4):
    #   y_ij = log(PSA(t_ij, psi_i) + 1) + sigma * e_ij,  e_ij ~ N(0,1)
    # This is additive error on the log(PSA + 1) transformed observation.
    # ---------------------------------------------------------------
    addSd <- 0.36 ; label("Additive residual SD on log(PSA + 1) scale (unitless)")  # Table I: sigma = 0.36
  })

  model({
    # 1. Individual PSA-kinetic parameters
    r        <- exp(lkg + etalkg)
    psa0     <- exp(lpsa0 + etalpsa0)
    le       <- logiteps + etalogiteps
    eps      <- exp(le) / (1 + exp(le))                # logit-normal back-transform => (0,1)
    tesc     <- exp(ltesc + etaltesc)

    # 2. Fixed mechanistic rate constants
    d_cell   <- exp(ld_fix)                            # 1/day
    delta    <- exp(ldelta_fix)                        # 1/day

    # 3. Time-varying chemotherapy effectiveness (paper Eq. 2)
    #    e(t) = eps  if t <= tesc
    #    e(t) = 0    if t >  tesc
    e_t      <- ifelse(t <= tesc, eps, 0)

    # 4. Mechanistic ODE system (paper Eq. 1)
    #    dC/dt   = (r * (1 - e(t)) - d) * C
    #    dPSA/dt = p * C - delta * PSA
    #
    #    The PSA production rate p (ng / (cell * day)) is NOT separately
    #    identifiable from PSA observations alone: it enters the PSA
    #    equation only through the product p * cells, and at quasi-steady
    #    state at treatment initiation p * cells(0) = delta * psa0.
    #    Setting p = 1 is a numerical convenience that leaves the PSA
    #    trajectory and the survival hazard unchanged; cells(0) is set
    #    by QSS to delta * psa0 / p = delta * psa0.
    p_psa    <- 1
    d/dt(cells) <- (r * (1 - e_t) - d_cell) * cells
    d/dt(psa)   <- p_psa * cells - delta * psa
    cells(0)    <- delta * psa0 / p_psa
    psa(0)      <- psa0

    # 5. Weibull baseline hazard with PSA-driven log-linear link
    #    (paper Eq. 5):
    #    h(t | PSA) = h0(t) * exp(beta * PSA(t)),
    #    h0(t)      = (k / lambda) * (t / lambda)^(k - 1).
    #
    #    A tiny offset del_t keeps the (t / lambda)^(k - 1) factor
    #    well-defined at t = 0 (here k > 1 so h0(0) = 0 in any case;
    #    the offset is a numerical safety net).
    lam_haz  <- exp(llam_haz)
    k_haz    <- exp(lk_haz)
    del_t    <- 1e-6
    hazard   <- (k_haz / lam_haz) *
                ((t + del_t) / lam_haz)^(k_haz - 1) *
                exp(e_psa_haz * psa)
    d/dt(cumhaz) <- hazard
    cumhaz(0)    <- 0
    sur          <- exp(-cumhaz)

    # 6. Observation and residual error (paper Eq. 4)
    logPSA1     <- log(psa + 1)
    logPSA1 ~ add(addSd)
  })
}
attr(Desmee_2015_PSA_mCRPC, "message") <-
  "This is a joint PSA-kinetics + Weibull-survival simulation model. The PSA-kinetics parameters (r, PSA0, eps, Tesc) are identical across all four published simulation scenarios; the survival parameters (lambda, k, beta) are set to the 'High link' scenario (beta=0.02, lambda=2150, k=1.5). The vignette reproduces all four scenarios by overriding lambda and beta. PSA observations are modelled on the log(PSA+1) scale; sur and hazard are derived outputs for forward simulation."
Desmee_2015_PSA_mCRPC
