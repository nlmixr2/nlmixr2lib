Othman_2007_carvedilol <- function() {
  description <- "Two-compartment population PK model for S(-)-carvedilol in healthy volunteers after oral administration of the immediate-release (IR) and the new controlled-release (CR) dosage forms of carvedilol (Othman 2007). Three parallel depot compartments encode the dosage-form-specific absorption: depot (CR, 3-stage time-varying KA), depot2 (IR morning, 2-stage), depot3 (IR evening, 2-stage). Diurnal variability in IR absorption is captured by separate morning and evening KAs and a lower IR-PM relative bioavailability."
  reference <- paste(
    "Othman AA, Tenero DM, Boyle DA, Eddington ND, Fossler MJ.",
    "Population pharmacokinetics of S(-)-carvedilol in healthy volunteers",
    "after administration of the immediate-release (IR) and the new",
    "controlled-release (CR) dosage forms of the racemate.",
    "AAPS J. 2007;9(2):E208-E218. doi:10.1208/aapsj0902023.",
    "Corrigendum: AAPS J. 2007;9(3) Article 37. doi:10.1208/aapsj0903037c",
    "(corrected alignment of Tables 1, 3, 4; numeric values are unchanged).",
    sep = " "
  )
  vignette <- "Othman_2007_carvedilol"

  paper_specific_etas <- c(
    "etalka_cr", "etalka_iram", "etalka_irpm",
    "etafrel",
    "etalka_cr_iov", "etafrel_cr_iov"
  )

  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = 96L,
    n_studies      = 3L,
    age_range      = "Healthy adults (specific demographic ranges not reported in the paper)",
    weight_range   = "Healthy adults (specific demographic ranges not reported in the paper)",
    sex_female_pct = NULL,
    race_ethnicity = "Not reported in the paper",
    disease_state  = "Healthy volunteers",
    dose_range     = paste(
      "IR: 25 mg of carvedilol (racemate, free base) every 12 hours for 2 doses;",
      "CR: single doses of 10, 20, 40, or 80 mg carvedilol phosphate",
      "(equivalent to 8.1, 16.2, 32.4, and 64.8 mg of carvedilol racemate free base)",
      "as a balanced 4-period crossover. All doses given under fed conditions."
    ),
    regions        = "Not reported",
    n_observations = 3328L,
    n_obs_ir       = 1270L,
    n_obs_cr       = 2058L,
    notes          = paste(
      "Three pooled studies (2 IR, 1 CR). Poor metabolizers of carvedilol were",
      "excluded by CYP2D6 genotyping. All doses administered under fed conditions",
      "(moderate-calorie, low-to-moderate-fat breakfast 30 minutes pre-dose).",
      "S(-)-carvedilol plasma concentrations measured by HPLC-MS/MS with LLOQ",
      "0.2 ng/mL. The analysis used S(-)-carvedilol (the beta-blocking enantiomer);",
      "doses entered NONMEM scaled to the administered S(-)-carvedilol free base",
      "amount (i.e., half the racemate dose). See Othman 2007 'Design of the Studies'",
      "and 'Population Analysis' sections (pages E209-E210)."
    )
  )

  ini({
    # Structural disposition parameters - corrigendum Table 1 (= main paper
    # Table 1; values are unchanged by the corrigendum, only the column
    # alignment was corrected). Final Model Point Estimate column of Table 3.
    lcl   <- log(149);    label("Apparent clearance CL/F (L/h)")                      # Table 1: 149 L/h (4.8% RSE)
    lvc   <- log(828);    label("Apparent central volume Vc/F (L)")                   # Table 1: 828 L (6.6% RSE)
    lvp   <- log(1150);   label("Apparent peripheral volume Vp/F (L)")                # Table 1: 1150 L (10.3% RSE)
    lq    <- log(94.7);   label("Apparent intercompartmental clearance Q/F (L/h)")    # Table 1: 94.7 L/h (5.3% RSE)

    # Time-varying first-order absorption rate constants - 3 stages for CR
    # (break points at 2 and 4 hours after CR dose) and 2 stages for each
    # of IR morning (IRAM) and IR evening (IRPM) (break point at 1 hour
    # after the respective IR dose). Source: Table 1 / Table 3.
    lka_cr_0to2   <- log(0.08);  label("Absorption rate constant CR, 0-2 h post-dose (1/h)")     # Table 1: 0.08 1/h (16.0% RSE)
    lka_cr_2to4   <- log(0.27);  label("Absorption rate constant CR, 2-4 h post-dose (1/h)")     # Table 1: 0.27 1/h (16.1% RSE)
    lka_cr_gt4    <- log(3.5);   label("Absorption rate constant CR, > 4 h post-dose (1/h)")     # Table 1: 3.5 1/h (17.7% RSE)
    lka_iram_0to1 <- log(0.92);  label("Absorption rate constant IR morning, 0-1 h post-dose (1/h)")  # Table 1: 0.92 1/h (21.5% RSE)
    lka_iram_gt1  <- log(8.79);  label("Absorption rate constant IR morning, > 1 h post-dose (1/h)")  # Table 1: 8.79 1/h (42.1% RSE)
    lka_irpm_0to1 <- log(0.42);  label("Absorption rate constant IR evening, 0-1 h post-dose (1/h)")  # Table 1: 0.42 1/h (26.7% RSE)
    lka_irpm_gt1  <- log(3.0);   label("Absorption rate constant IR evening, > 1 h post-dose (1/h)")  # Table 1: 3.0 1/h (31.9% RSE)

    # Relative bioavailability - IR morning (IRAM) is the reference (Frel = 1,
    # not estimated). Frel for CR and IR evening (IRPM) are estimated
    # relative to IRAM.
    lfrel_cr   <- log(0.76);  label("Relative bioavailability CR vs IR morning (fraction)")    # Table 1: 0.76 (7.4% RSE)
    lfrel_irpm <- log(0.80);  label("Relative bioavailability IR evening vs IR morning (fraction)")  # Table 1: 0.80 (3.2% RSE)

    # Absorption lag times - one common IR lag (estimated) and a separate
    # CR lag (fixed at 0.23 h based on sensitivity analysis; estimation of
    # the CR lag time gave rounding errors during NONMEM minimization).
    ltlag_ir <- log(0.20);          label("Absorption lag time for IR doses (h)")  # Table 1: 0.20 h (5.3% RSE)
    ltlag_cr <- fixed(log(0.23));   label("Absorption lag time for CR doses (h, fixed by sensitivity analysis)")  # Table 1: 0.23 h (Fixed)

    # Inter-individual variability (omega^2 on the log scale of each
    # exponentiated typical-value parameter; Eq. 1 in Othman 2007).
    # Variance values taken from the Final Model Point Estimate column of
    # Table 3. The %ISV in Table 1 is approximately omega * 100 (i.e.,
    # the standard deviation reported as a percent); for example
    # omega^2 = 0.14 -> %ISV ~ 37.4% matches Table 1's 37.1%.
    # Encoded here as separate (uncorrelated) IIV terms; Othman 2007
    # reports "ETA matrix plots did not show any trends of correlation
    # between the estimated intersubject variability parameters".
    etalvc       ~ 0.14   # Table 3: omega^2 Vc/F = 0.14
    etalka_cr    ~ 0.90   # Table 3: omega^2 KA,CR = 0.90 (one common variance for the 3 CR KA stages)
    etalka_iram  ~ 1.97   # Table 3: omega^2 KA,IR,AM = 1.97 (one common variance for the 2 IR-AM KA stages)
    etalka_irpm  ~ 3.73   # Table 3: omega^2 KA,IR,PM = 3.73 (one common variance for the 2 IR-PM KA stages)
    etafrel      ~ 0.11   # Table 3: omega^2 Frel = 0.11 (one common variance shared across Frel_CR, Frel_IRAM, Frel_IRPM)

    # Interoccasion variability on KA,CR and Frel,CR. The crossover CR
    # study has 4 occasions per subject (10, 20, 40, 80 mg). Othman 2007
    # uses Eq. 2 (P_ij = TVP * exp(eta_i + kappa_ij)). For the model
    # library these IOV terms are encoded as additional log-scale random
    # effects sampled per-subject by default (see Assumptions and
    # deviations in the vignette for the simplification); for per-occasion
    # simulation the user should re-sample these terms at each CR dose.
    etalka_cr_iov  ~ 1.29   # Table 3: pi^2 KA,CR = 1.29
    etafrel_cr_iov ~ 0.02   # Table 3: pi^2 Frel,CR = 0.02

    # Residual error - Othman 2007 used a log-normal exponential model
    # ln(Yobs) = ln(Ypred) + epsilon, epsilon ~ N(0, sigma^2 = 0.10).
    # This maps to the lnorm family in nlmixr2 with the log-scale SD
    # given by sqrt(sigma^2) = sqrt(0.10) ~ 0.3162. Othman 2007 also
    # added an IIV term on epsilon (Eq. 4, omega^2_RES = 0.05); that
    # eta-on-epsilon refinement is omitted here (see vignette
    # Assumptions and deviations for the rationale).
    expSd <- sqrt(0.10);  label("Log-scale residual SD for Cc (lognormal)")  # Table 3: sigma^2 = 0.10
  })

  model({
    # Per-subject structural PK parameters
    cl <- exp(lcl)
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp)
    q  <- exp(lq)

    # Time-varying KA per formulation. depot = CR, depot2 = IR morning,
    # depot3 = IR evening. tad(depot_n) returns time since the most
    # recent dose into that compartment, after the alag(...) shift.
    # The 3-stage CR KA switches at 2 h and 4 h (Table 1); each 2-stage
    # IR KA switches at 1 h. Below-break point assignment first; later
    # if-statements over-write when the threshold is reached so the most
    # recent applicable KA wins.
    # The IOV deviation (etalka_cr_iov) is factored out so each exp()
    # contains a single eta term (required for mu-referencing in
    # nlmixr2: rxode2 currently does not accept theta + eta1 + eta2).
    iov_ka_cr <- exp(etalka_cr_iov)
    ka_cr <- exp(lka_cr_0to2 + etalka_cr) * iov_ka_cr
    if (tad(depot) >= 2) ka_cr <- exp(lka_cr_2to4 + etalka_cr) * iov_ka_cr
    if (tad(depot) >= 4) ka_cr <- exp(lka_cr_gt4  + etalka_cr) * iov_ka_cr

    ka_iram <- exp(lka_iram_0to1 + etalka_iram)
    if (tad(depot2) >= 1) ka_iram <- exp(lka_iram_gt1 + etalka_iram)

    ka_irpm <- exp(lka_irpm_0to1 + etalka_irpm)
    if (tad(depot3) >= 1) ka_irpm <- exp(lka_irpm_gt1 + etalka_irpm)

    # Relative bioavailability factors. Frel_IRAM = 1 is the reference
    # (so f(depot2) carries only the eta_Frel deviation). Frel_CR and
    # Frel_IRPM are exponentiated typical values scaled by the same
    # eta_Frel; Frel_CR additionally carries an IOV deviation
    # (factored as iov_frel_cr to keep each exp() to a single eta).
    iov_frel_cr <- exp(etafrel_cr_iov)
    frel_cr   <- exp(lfrel_cr   + etafrel) * iov_frel_cr
    frel_iram <- exp(etafrel)
    frel_irpm <- exp(lfrel_irpm + etafrel)

    # Lag times
    tlag_cr <- exp(ltlag_cr)
    tlag_ir <- exp(ltlag_ir)

    # 2-compartment disposition micro-constants
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ODE system - three parallel oral absorption routes feeding a
    # single central compartment, with two-compartment linear
    # disposition.
    d/dt(depot)       <- -ka_cr   * depot
    d/dt(depot2)      <- -ka_iram * depot2
    d/dt(depot3)      <- -ka_irpm * depot3
    d/dt(central)     <-  ka_cr * depot + ka_iram * depot2 + ka_irpm * depot3 -
                          kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    f(depot)     <- frel_cr
    f(depot2)    <- frel_iram
    f(depot3)    <- frel_irpm

    alag(depot)  <- tlag_cr
    alag(depot2) <- tlag_ir
    alag(depot3) <- tlag_ir

    # Concentration: dose mg, vc in L gives mg/L; multiply by 1000 for
    # ng/mL (Othman 2007 reports Cc in ng/mL, e.g. Tables 4 and Cmax).
    Cc <- central / vc * 1000
    Cc ~ lnorm(expSd)
  })
}
