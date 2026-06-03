Chang_2011_PF04455242_human <- function() {
  description <- "Two-compartment population PK and reduced direct-response PD model for PF-04455242 (kappa opioid receptor antagonist) in healthy adult volunteers (Chang 2011). PK is fit with zero-order oral absorption (duration D1) and lag time (ALAG1) into the central compartment; residual error uses the log-transform-both-sides (lognormal) form. PD is the reduced antagonism model (Eq. 11/12) that replaces the spiradoline PK with a deterministic Weibull-scaled placebo prolactin profile and predicts the time-matched prolactin response under PF-04455242 antagonism. Simulation time t = 0 must be aligned with the IM spiradoline challenge dose; PF-04455242 is dosed earlier (typically t = -1 h in the proof-of-mechanism study)."
  reference <- paste(
    "Chang C, Byon W, Lu Y, Jacobsen LK, Badura LL, Sawant-Basak A, Miller E,",
    "Liu J, Grimwood S, Wang EQ, Maurer TS.",
    "Quantitative PK-PD Model-Based Translational Pharmacology of a Novel",
    "Kappa Opioid Receptor Antagonist Between Rats and Humans.",
    "AAPS J. 2011;13(4):565-575. doi:10.1208/s12248-011-9296-3.",
    sep = " "
  )
  vignette <- "Chang_2011_PF04455242_prolactin"

  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = 69L,
    n_studies      = 3L,
    age_range      = "healthy adult",
    sex_female_pct = 0,
    disease_state  = "healthy male adult volunteers (proof-of-mechanism study); single- and multiple-dose PK studies pooled with the POM study for the population PK fit",
    dose_range     = paste(
      "PK pool: PF-04455242 0.5-30 mg PO single dose (n=18 subjects, 22 in reference (22))",
      "and PF-04455242 every-6-hour multiple dose followed by a single dose on day 7 (n=27 subjects, reference (23)).",
      "POM study: PF-04455242 18 or 30 mg PO at t=-1 h followed by spiradoline 3.2 ug/kg IM at t=0 (n=24 healthy male subjects, 8 per arm).",
      sep = " "
    ),
    regions        = "Pfizer Inc., USA (Clinical Pharmacology Primary Care Business Unit)",
    notes          = paste(
      "PK was fit on 943 PF-04455242 concentrations from the single- and multiple-dose phase 1 studies.",
      "The PD layer (this file) is the reduced antagonism model from Eq. 11 (paper Methods 'Clinical POM Prediction' and 'Clinical PK-PD Model') in which the spiradoline-driven STIMplacebo is described empirically by a Weibull function (Eq. 12) digitised from a placebo-arm prolactin time course.",
      "Sequential PK-PD: individual PF-04455242 concentrations were predicted using the population PK model and combined with the prolactin observations from the POM study.",
      sep = " "
    )
  )

  ini({
    # ------------------------------------------------------------------
    # PF-04455242 PK -- two-compartment, zero-order absorption into the
    # central compartment with duration D1 and absorption lag ALAG1.
    # Chang 2011 Table III (final population PK estimates).
    # ------------------------------------------------------------------
    lcl   <- log(54.6);  label("PF-04455242 apparent oral clearance CL/F (L/h)")             # Table III
    lvc   <- log(327);   label("PF-04455242 apparent central volume of distribution Vc/F (L)") # Table III
    ld1   <- log(1.24);  label("PF-04455242 zero-order absorption duration D1 (h)")          # Table III
    lvp   <- log(58.5);  label("PF-04455242 apparent peripheral volume Vp/F (L)")            # Table III
    lq    <- log(6.51);  label("PF-04455242 apparent intercompartmental clearance Q/F (L/h)") # Table III
    ltlag <- log(0.374); label("PF-04455242 absorption lag time ALAG1 (h)")                  # Table III

    # IIV on PK parameters. omega translated from reported CV% via
    # omega^2 = log(CV^2 + 1) (lognormal exponential random-effect model).
    # Block correlation between etalcl and etalvc: Chang 2011 Table III
    # reports cov(omega-CL, omega-Vc) = 48.5. The only physically-
    # consistent interpretation (correlation in [-1, 1]) is that 48.5 is
    # the correlation coefficient expressed as a percentage, i.e.,
    # correlation = 0.485. The log-scale covariance is then
    #   0.485 * sqrt(var_lcl * var_lvc)
    # = 0.485 * sqrt(log(0.537^2+1) * log(0.485^2+1))
    # = 0.485 * sqrt(0.2533 * 0.2113)
    # = 0.485 * 0.2314
    # = 0.1122
    # See vignette Errata for the alternative interpretations considered.
    etalcl + etalvc ~ c(0.2533, 0.1122, 0.2113)
    # Table III: omega-CL = 53.7%CV -> log(0.537^2+1) = 0.2533
    #            omega-Vc = 48.5%CV -> log(0.485^2+1) = 0.2113
    #            cov(omega-CL, omega-Vc) interpreted as r = 0.485 (see comment above)
    etald1   ~ 0.0896  # Table III, omega-D1   = 30.6%CV -> log(0.306^2+1) = 0.0896
    etalvp   ~ 0.0721  # Table III, omega-Vp   = 27.3%CV -> log(0.273^2+1) = 0.0721
    etalq    ~ 0.0721  # Table III, omega-Q    = 27.3%CV -> log(0.273^2+1) = 0.0721
    etaltlag ~ 0.0687  # Table III, omega-ALAG1 = 26.6%CV -> log(0.266^2+1) = 0.0687

    expSd <- 0.326; label("PF-04455242 log-scale residual SD (unitless, lnorm; LTBS)")  # Table III, 32.6% under the LTBS specification Eq. 14

    # ------------------------------------------------------------------
    # Prolactin PD -- reduced antagonism model (Chang 2011 Eq. 11):
    #   PRL = BL + STIMplacebo * (x^gamma - 1) /
    #              ((x - 1) * (1 + Cpf/Ki)^gamma + x^gamma - x)
    # with x = 2 (paper fixes the scaling factor x to 2 based on the
    # observed two-fold proportionality between spiradoline 1.6 and
    # 3.2 ug/kg IM prolactin elevation -- see Eq. 6 / Eq. 7 / Eq. 11
    # derivation in Methods 'Clinical POM Prediction').
    # STIMplacebo is the Weibull-scaled placebo prolactin elevation
    # above baseline (Eq. 12):
    #   STIMplacebo = BL * (kappa/lambda) * (TIME/lambda)^(kappa-1)
    #                 * exp(-(TIME/lambda)^kappa)
    # where TIME is the simulation time aligned with spiradoline
    # administration (t = 0 at the IM spiradoline dose).
    # ------------------------------------------------------------------
    lrbase  <- log(9.53);  label("Baseline plasma prolactin BL (ng/mL)")                          # Table IV
    lkappa  <- log(2.49);  label("Weibull placebo shape parameter kappa (unitless)")              # Table IV
    llambda <- log(0.943); label("Weibull placebo scale parameter lambda (h)")                    # Table IV
    lki     <- log(39.2);  label("PF-04455242 in vivo Ki at human KOR (ng/mL)")                   # Table IV
    lhill   <- log(1.5);   label("Hill coefficient (paper gamma) on the antagonism term (unitless)")  # Table IV (text mentions 1.66 once in Results, but Table IV and Discussion both report 1.5; see vignette Errata)
    x_scale <- fixed(2);   label("Scaling factor x relating spiradoline 3.2 vs 1.6 ug/kg prolactin elevation (unitless)")  # Methods 'Clinical POM Prediction'; fixed per the paper's derivation

    # IIV on baseline: Chang 2011 Table IV reports omega-BL = 3.1 (ng/mL)
    # -- the parenthesised unit indicates ADDITIVE IIV on the linear
    # scale, not a %CV. Encoded here as multiplicative log-scale IIV
    # using the equivalent CV = 3.1 / 9.53 = 32.5%, omega^2 =
    # log(0.325^2 + 1) = 0.1003. This approximates the paper's additive
    # IIV; the difference is < 5% for prolactin values within
    # +-2 SD of baseline and is documented in vignette Errata.
    etalrbase  ~ 0.1003   # Table IV, omega-BL = 3.1 ng/mL (additive); approximated as CV = 32.5% (see comment above)
    etallambda ~ 0.0388   # Table IV, omega-lambda = 19.9%CV -> log(0.199^2+1) = 0.0388

    propSd_PRL <- 0.175; label("Human prolactin proportional residual error (fraction)")  # Table IV, 17.5%
  })

  model({
    # ------------------------------------------------------------------
    # 1. Individual PK parameters.
    # ------------------------------------------------------------------
    cl   <- exp(lcl + etalcl)
    vc   <- exp(lvc + etalvc)
    d1   <- exp(ld1 + etald1)
    vp   <- exp(lvp + etalvp)
    q    <- exp(lq  + etalq)
    tlag <- exp(ltlag + etaltlag)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ------------------------------------------------------------------
    # 2. PD typical-value parameters.
    # ------------------------------------------------------------------
    rbase  <- exp(lrbase + etalrbase)
    kappa  <- exp(lkappa)
    lambda <- exp(llambda + etallambda)
    ki     <- exp(lki)
    hill   <- exp(lhill)

    # ------------------------------------------------------------------
    # 3. PK ODEs. Two-compartment with zero-order absorption into central
    # (duration d1, lag tlag). PF-04455242 doses enter cmt=central with
    # the duration / lag set per dose event via dur() and alag().
    # ------------------------------------------------------------------
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    dur(central)  <- d1
    alag(central) <- tlag

    # ------------------------------------------------------------------
    # 4. Plasma concentration and prolactin response. Internal central
    # is in mg (dose unit) and Vc/F is in L, so central/vc has units
    # mg/L = ug/mL. Multiply by 1000 to display in ng/mL, matching the
    # paper's reporting (Tables III/IV, Figures 6-7) and keeping Cc on
    # the same scale as the PD parameter Ki (ng/mL).
    # Cc      -- PF-04455242 plasma concentration (parent observation,
    #            ng/mL; lognormal residual under LTBS Eq. 14).
    # PRL     -- prolactin response (paper-named PD output, ng/mL). The
    #            reduced antagonism model (Eq. 11) and the Weibull
    #            placebo (Eq. 12) are evaluated at simulation time t
    #            aligned with the IM spiradoline dose (t = 0 = spiradoline
    #            administration). The PF-04455242 dose in the POM study
    #            was administered at t = -1 h.
    # ------------------------------------------------------------------
    Cc <- central / vc * 1000

    # Weibull placebo (Eq. 12). The (t > 0) gate keeps the formula at
    # zero for t <= 0 (where the Weibull PDF is undefined / pre-stimulus).
    stim_placebo <- (t > 0) * rbase * (kappa / lambda) *
                    (t / lambda)^(kappa - 1) *
                    exp(-(t / lambda)^kappa)

    # Antagonism factor (Eq. 11 dimensionless quotient). Reduces to 1
    # when Cpf = 0 (placebo arm) and shrinks toward zero as Cpf grows
    # large relative to Ki (full antagonism).
    antag <- (x_scale^hill - 1) /
             ((x_scale - 1) * (1 + Cc / ki)^hill + x_scale^hill - x_scale)

    PRL <- rbase + stim_placebo * antag

    Cc  ~ lnorm(expSd)
    PRL ~ prop(propSd_PRL)
  })
}
