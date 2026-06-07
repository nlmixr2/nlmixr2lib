Wang_2011_raltegravir <- function() {
  description <- paste(
    "Population PK model for plasma and intracellular (PBMC)",
    "raltegravir after a single 400 mg oral dose in six healthy male",
    "Singaporean volunteers (Wang 2011). Plasma PK is described by a",
    "one-compartment model with first-order elimination preceded by a",
    "chain of two transit compartments between depot and central",
    "(Kappelhoff 2005 transit-absorption parameterisation: MAT =",
    "(n + 1) / ktr with n = 2 transit compartments). Bioavailability",
    "F is implicit in the apparent CL/F and V/F. Intracellular (PBMC)",
    "raltegravir is described as an empirical partition of the",
    "predicted plasma concentration via the paper's accumulation ratio",
    "ACR (point estimate 11.2%) with its own inter-individual",
    "variability and exponential residual error (Wang 2011 Eq.:",
    "C_IC,obs = ACR * C_plasma,pred * exp(eps_IC)). The packaged model",
    "maps ACR to the canonical paper-named bare parameter `frac` and",
    "its log-transformed primary `lfrac`. No baseline covariates were",
    "retained in the final model. Note that the paper's term",
    "'accumulation ratio' is a misnomer in the conventional sense --",
    "ACR < 1 means raltegravir does NOT accumulate intracellularly,",
    "consistent with simple diffusion of unbound drug into PBMCs (the",
    "authors' Conclusions).",
    sep = " "
  )
  reference <- paste(
    "Wang L, Soon GH, Seng KY, Li J, Lee E, Yong EL, Goh BC, Flexner C,",
    "Lee L.",
    "Pharmacokinetic Modeling of Plasma and Intracellular Concentrations",
    "of Raltegravir in Healthy Volunteers.",
    "Antimicrobial Agents and Chemotherapy. 2011 Sep;55(9):4090-4095.",
    "doi:10.1128/AAC.00593-11.",
    sep = " "
  )
  vignette <- "Wang_2011_raltegravir"
  # The source paper reports raltegravir concentrations in nmol/L and doses
  # in mg; the packaged model uses mass-balanced mg / mg/L throughout
  # (matching Lee_2016_raltegravir). Conversions for cross-checking against
  # the paper's Table 1 NCA values use raltegravir molecular weight
  # 444.42 g/mol: 1 nmol/L = 444.42e-6 mg/L, so reported Cmax 2246 nmol/L =
  # 0.9982 mg/L and AUC0-inf 13119 nmol/L*h = 5.831 mg*h/L.
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = 6L,
    n_studies      = 1L,
    age_range      = "31-36 years",
    age_median     = "33.5 years",
    weight_range   = "56.5-87.1 kg",
    weight_median  = "67.4 kg",
    sex_female_pct = 0,
    disease_state  = paste(
      "Healthy male volunteers (HIV-negative). Inclusion required",
      "smoking <10 cigarettes/day, hemoglobin >10.9 g/dl, platelet count",
      ">=125,000/mm^3, creatinine clearance >=60 mL/min, and lipase or",
      "pancreatic amylase <1.1 x upper limit of normal. Subjects on any",
      "concomitant medication beyond aspirin / acetaminophen /",
      "chlorpheniramine / multivitamins were excluded.",
      sep = " "
    ),
    dose_range     = paste(
      "Single 400 mg raltegravir orally with water under fasted",
      "conditions on study day 1. Paired blood samples for plasma and",
      "PBMC were collected predose (0 h) and at 4, 8, 12, 24, and 48 h",
      "post-dose; additional plasma samples were collected at 0.5, 1,",
      "1.5, 2, 3, 5, 6, and 10 h post-dose. PBMC concentrations were",
      "computed assuming a cell volume of 0.4 pL.",
      sep = " "
    ),
    regions        = "Singapore (single centre, National University Health System)",
    age_eligibility = "21-65 years",
    notes          = paste(
      "Baseline demographics from Wang 2011 Results, 'Study population",
      "and safety' (page 4090) and 'Pharmacokinetic parameters' (page",
      "4091). Final population pharmacokinetic estimates are reported",
      "in Table 1 of the source. The visual predictive checks based on",
      "500 simulated profiles (Wang 2011 Fig. 4) covered 97% of the",
      "observed concentrations within the 95% prediction interval.",
      sep = " "
    )
  )

  ini({
    # ---- Structural fixed effects (Wang 2011 Table 1) ----
    # Number of transition compartments is reported as a structural
    # choice (n = 2), not as an estimated theta. It is encoded
    # implicitly in the ODE chain (depot -> transit1 -> transit2 ->
    # central) and in the relationship ktr = (n + 1) / mat = 3 / mat
    # inside model().
    lmat   <- log(0.887); label("Mean absorption time MAT (h)")                                # Wang 2011 Table 1: MAT = 0.887 h (RSE 27.6%)
    lcl    <- log(39.1);  label("Apparent oral clearance CL/F (L/h)")                          # Wang 2011 Table 1: CL/F = 39.1 L/h (RSE 20.8%)
    lvc    <- log(272);   label("Apparent oral central volume V/F (L)")                        # Wang 2011 Table 1: V/F = 272 L (RSE 21.8%)
    lfrac  <- log(0.112); label("Intracellular-to-plasma accumulation ratio ACR (unitless)")   # Wang 2011 Table 1: Accumulation ratio = 11.2% (RSE 35%)

    # ---- Inter-individual variability ----
    # The paper reports IIV as %CV on the log-normal scale; the
    # internal omega^2 on log scale is recovered as
    # omega^2 = log(CV^2 + 1).
    #   CV(MAT)   = 66.2%  -> omega^2 = log(1 + 0.662^2) = 0.36347
    #   CV(CL/F)  = 56.4%  -> omega^2 = log(1 + 0.564^2) = 0.27620
    #   CV(V/F)   = 56.7%  -> omega^2 = log(1 + 0.567^2) = 0.27877
    #   CV(ACR)   = 104%   -> omega^2 = log(1 + 1.04^2)  = 0.73315
    # Correlation between CL/F and V/F is reported as 0.559 (Table 1),
    # so the covariance entered into the lower-triangular block is
    #   cov(lcl, lvc) = 0.559 * sqrt(0.27620 * 0.27877) = 0.15506.
    etalcl + etalvc ~ c(0.27620,
                        0.15506, 0.27877)                                                       # Wang 2011 Table 1: IIV(CL/F)=56.4%, IIV(V/F)=56.7%, correlation = 0.559
    etalmat        ~ 0.36347                                                                    # Wang 2011 Table 1: IIV(MAT) = 66.2%
    etalfrac       ~ 0.73315                                                                    # Wang 2011 Table 1: IIV(accumulation ratio) = 104%

    # ---- Residual error -- exponential (log-additive) ----
    # The paper Methods state interindividual and residual unexplained
    # variability were estimated with an exponential error model. In
    # NONMEM this corresponds to Y = F * exp(eps), eps ~ N(0, sigma^2),
    # which maps to nlmixr2's lnorm(expSd) with expSd = sigma on the
    # log scale. Table 1 reports the percentages as SQRT(SIGMA) on log
    # scale (the standard NONMEM convention for $SIGMA in an
    # exponential error block): 86.7% -> expSd = 0.867; 64.9% ->
    # expSd_Cic = 0.649.
    expSd      <- 0.867; label("Exponential residual SD for plasma raltegravir (log scale)")          # Wang 2011 Table 1: Exponential error of plasma concns = 86.7% (RSE 15.4%)
    expSd_Cic  <- 0.649; label("Exponential residual SD for intracellular raltegravir (log scale)")   # Wang 2011 Table 1: Exponential error of intracellular concns = 64.9% (RSE 51.1%)
  })

  model({
    # ---- Individual parameters ----
    mat  <- exp(lmat + etalmat)
    cl   <- exp(lcl + etalcl)
    vc   <- exp(lvc + etalvc)
    frac <- exp(lfrac + etalfrac)

    # ---- Transit-chain rate constant (Kappelhoff 2005) ----
    # Wang 2011 Eq.: MAT = (n + 1) / ktr with n = 2 transit
    # compartments between depot and central. Inverting:
    # ktr = (n + 1) / MAT = 3 / MAT. The same first-order rate
    # constant drives each of the three transitions in the chain
    # (depot -> transit1, transit1 -> transit2, transit2 -> central).
    ktr <- 3 / mat

    # ---- Micro-constant ----
    kel <- cl / vc

    # ---- ODE system ----
    d/dt(depot)    <- -ktr * depot
    d/dt(transit1) <-  ktr * depot    - ktr * transit1
    d/dt(transit2) <-  ktr * transit1 - ktr * transit2
    d/dt(central)  <-  ktr * transit2 - kel * central

    # ---- Observations ----
    Cc  <- central / vc
    # Wang 2011 Eq.: C_IC,obs = ACR * C_plasma,pred * exp(eps_IC). The
    # intracellular concentration is an empirical multiplicative
    # transform of the plasma prediction, not a separate compartment.
    Cic <- frac * Cc

    Cc  ~ lnorm(expSd)
    Cic ~ lnorm(expSd_Cic)
  })
}
