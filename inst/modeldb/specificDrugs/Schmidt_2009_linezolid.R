Schmidt_2009_linezolid <- function() {
  description <- paste(
    "In vitro (Staphylococcus aureus MRSA strain OC2878).",
    "Mechanism-based PD model of bacterial-killing time-kill curves for",
    "linezolid, the first-in-class FDA-approved oxazolidinone (Schmidt",
    "2009). Susceptibility-based two-subpopulation structure: an active",
    "self-replicating susceptible pool with logistic carrying-capacity",
    "limit and a dormant persister pool that is insusceptible to",
    "killing; first-order S->P conversion (P->S held fixed at 0),",
    "natural-death loss from both pools, exponential turn-on of growth",
    "and of drug-induced killing, and Emax killing of the susceptible",
    "subpopulation by the antibiotic. Linezolid was experimentally",
    "stable over 24 h so the published model carries kdeg = 0; for",
    "dynamic syringe-replacement experiments the user supplies the",
    "dilution-equivalent rate (~log(2)/5 1/h for a t1/2 ~ 5 h linezolid",
    "regimen) via rxSolve(..., params = c(kdeg = ...)). The same joint",
    "fit is shared with Schmidt_2009_rwj416457 (only EC50 and kdeg",
    "differ)."
  )
  reference <- paste(
    "Schmidt S, Sabarinath SN, Barbour A, Abbanat D, Manitpisitkul P,",
    "Sha S, Derendorf H. (2009). Pharmacokinetic-pharmacodynamic",
    "modeling of the in vitro activities of oxazolidinone antimicrobial",
    "agents against methicillin-resistant Staphylococcus aureus.",
    "Antimicrob Agents Chemother 53(12):5039-5045.",
    "doi:10.1128/AAC.00633-09."
  )
  vignette <- "Schmidt_2009_oxazolidinones"
  units <- list(
    time          = "hour",
    dosing        = "ug/mL (initial antibiotic concentration in MHB)",
    concentration = "log10 CFU/mL (bacterial Cc output); ug/mL (antibiotic central state)"
  )

  # bact_susceptible and bact_persister are paper-mechanistic states
  # (Schmidt 2009 Fig. 1: pools S and P). Declared so the convention
  # check does not flag them as non-canonical.
  paper_specific_compartment_pattern <- "^bact_"

  # No epidemiological covariates: this is an in vitro time-kill model.
  # The antibiotic concentration is a state (central), not a covariate.
  covariateData <- list()

  population <- list(
    species          = "in vitro (Staphylococcus aureus, MRSA strain OC2878)",
    n_subjects       = NA_integer_,
    n_studies        = 1L,
    organism         = "Methicillin-resistant Staphylococcus aureus, strain OC2878 (Johnson & Johnson Pharmaceutical R&D)",
    mic_values       = c(rwj416457 = "0.5 ug/mL", linezolid = "1.0 ug/mL"),
    system           = paste(
      "Two in vitro models simultaneously fit: (i) static time-kill curves",
      "in 50-mL cell-culture flasks at constant concentrations 0.25-16x MIC",
      "for 24 h; (ii) dynamic time-kill curves in the same flasks coupled to",
      "syringe systems that replaced 4.8 mL of broth every 2 h to simulate",
      "the linezolid human elimination half-life of ~5 h."
    ),
    medium           = "Mueller-Hinton broth (MHB), Difco; 20 mL per flask; 37 degC",
    inoculum         = "~5e5 CFU/mL after a 2-h pre-incubation with no antibiotic",
    dose_range       = "0.25-16x MIC, i.e. 0.25-16 ug/mL linezolid",
    drug_stability   = "Linezolid was completely stable over 24 h (Results 'Drug stability'); kdeg = 0",
    notes            = paste(
      "Triplicate static and dynamic experiments per drug per concentration.",
      "Parameters were estimated by a simultaneous fit of the same",
      "susceptibility-based model to RWJ-416457 and linezolid data with one",
      "shared set of structural and variance parameters (only EC50 differs",
      "between drugs); see Schmidt 2009 Table 1 footnote and Discussion",
      "paragraph 'allowed inferences about the potencies of RWJ-416457 and",
      "linezolid'."
    )
  )

  ini({
    # =================================================================
    # Bacterial growth, death, and subpopulation transfer (Schmidt 2009
    # Table 1; the model equations are stated in Materials and Methods
    # 'Mathematical modeling', Eqs 1-3).
    # =================================================================
    lks   <- log(1.19);   label("Log growth-rate constant ks (1/h)")                    # Table 1: ks = 1.19 +/- 0.061
    lkmax <- log(1.65);   label("Log maximum drug-induced kill-rate constant kmax (1/h)") # Table 1: kmax = 1.65 +/- 0.065

    # EC50 is the only structural parameter that differs between the two
    # oxazolidinones in the joint fit; this file carries the linezolid value.
    lec50 <- log(1.39);   label("Log linezolid EC50 (ug/mL)")                           # Table 1: EC50,LZD = 1.39 +/- 0.207

    # Carrying capacity (Nmax) and growth-onset delay (dgs) were fit on
    # the antibiotic-free growth-control data first and held constant in
    # the joint fit per Mathematical-modeling paragraph 4: "Nmax, dgs,
    # and kd were estimated from the growth control data and were
    # assumed to be constant." Their between-experiment IIVs are still
    # estimated, however -- see the variance block below.
    lnmax <- fixed(log(3.39e9));  label("Log maximum bacterial population Nmax (CFU/mL); FIXED from growth control") # Table 1: Nmax FIXED to 3.39e9
    ldgs  <- fixed(log(0.24));    label("Log growth-onset delay-rate dgs (1/h); FIXED from growth control")          # Table 1: dgs FIXED to 0.24

    # dks is estimated in the joint fit
    ldks  <- log(0.50);   label("Log kill-onset delay-rate dks (1/h)")                  # Table 1: dks = 0.50 +/- 0.049

    # Persister-pool transfer rates. Schmidt 2009 fixes kps to 0
    # (reference 22) because the experimental data did not identify a
    # back-transfer rate; kps therefore stays in linear scale so the
    # zero value can be encoded exactly.
    lksp  <- log(0.004);  label("Log susceptible-to-persister transfer rate ksp (1/h)") # Table 1: ksp = 0.004 +/- 0.002
    kps   <- fixed(0);    label("Persister-to-susceptible transfer rate kps (1/h); FIXED to 0 per Schmidt 2009 ref 22")

    # Natural-death rate of both pools (estimated from the growth-control
    # fit and then held fixed in the joint fit).
    lkd   <- fixed(log(0.015));   label("Log natural death-rate kd (1/h); FIXED from growth control")  # Table 1: kd FIXED to 0.015

    # Initial fraction of the inoculum in the susceptible pool F (the
    # rest 1-F seeds the persister pool). Logit-transformed to keep F in
    # (0,1) at simulation time.
    logitf <- log(0.83 / (1 - 0.83)); label("Logit-transformed initial susceptible fraction F (F=0.83)") # Table 1: F = 0.83 +/- 0.023

    # First-order antibiotic loss from the MHB medium. Linezolid was
    # experimentally stable over 24 h (Results 'Drug stability' /
    # Fig. 3B), so kdeg = 0 for static experiments. For dynamic
    # syringe-replacement experiments the user replaces this value with
    # the equivalent dilution rate (e.g.,
    # rxSolve(mod, ev, params = c(kdeg = log(2)/5))) so the central
    # state declines with the published linezolid t1/2 ~ 5 h. kdeg is
    # carried on the linear scale rather than via lkdeg because log(0)
    # is not finite; this matches the RWJ-416457 sibling file.
    kdeg  <- fixed(0); label("Drug degradation rate kdeg (1/h); FIXED to 0 (linezolid completely stable over 24 h per Results 'Drug stability')")

    # Initial inoculum -- the experimental setup uses ~5e5 CFU/mL at the
    # start of the antibiotic-exposure window (Materials and Methods
    # 'Organisms' / 'Static time-kill curves').
    lninit <- fixed(log(5e5)); label("Log initial total inoculum N0 (CFU/mL); FIXED from Methods")

    # =================================================================
    # Between-experiment IIV (Schmidt 2009 Table 1 'Variance model').
    # Reported as omega^2 (variance on the log scale) per NONMEM
    # convention; the bootstrap column in Table 1 prints the matching SD
    # (sqrt(omega^2)) which is how each value below was checked.
    # =================================================================
    etalks    ~ 0.013   # Table 1: eta(ks)   = 0.013 +/- 0.005 (bootstrap SD 0.11 ~ sqrt(0.013))
    etalnmax  ~ 0.219   # Table 1: eta(Nmax) = 0.219 +/- 0.082 (bootstrap SD 0.51 ~ sqrt(0.219))
    etaldgs   ~ 0.184   # Table 1: eta(dgs)  = 0.184 +/- 0.071 (bootstrap SD 0.43 ~ sqrt(0.184))
    etaldks   ~ 0.079   # Table 1: eta(dks)  = 0.079 +/- 0.029 (bootstrap SD 0.26 ~ sqrt(0.079))

    # =================================================================
    # Residual error. Schmidt 2009 fits a "log error model" (Data
    # analysis paragraph) on log-transformed bacterial counts; NONMEM
    # LOG() is natural log, so SIGMA = 0.29 (Table 1 'Residual
    # variability') is a variance on the natural-log scale (bootstrap SD
    # 0.53 ~ sqrt(0.29)). The output is reported in log10 CFU/mL for
    # microbiology readability, so addSd is rescaled from natural-log to
    # log10 by dividing the natural-log SD sqrt(0.29) = 0.539 by ln(10):
    #   addSd_log10 = sqrt(0.29) / log(10) = 0.234
    # See vignette Assumptions and deviations for the alternative
    # interpretation (the paper does not state the log base explicitly).
    # =================================================================
    addSd <- 0.234; label("Additive residual SD on log10 CFU/mL scale; converted from SIGMA = 0.29 (natural-log variance)")
  })

  model({
    # ---- Back-transform structural parameters from the log scale ----
    ks    <- exp(lks   + etalks)
    kmax  <- exp(lkmax)
    ec50  <- exp(lec50)
    nmax  <- exp(lnmax + etalnmax)
    dgs   <- exp(ldgs  + etaldgs)
    dks   <- exp(ldks  + etaldks)
    ksp   <- exp(lksp)
    kd    <- exp(lkd)
    fsusc <- exp(logitf) / (1 + exp(logitf))   # back-transform logit -> initial susceptible fraction
    ninit <- exp(lninit)

    # ---- Time-varying turn-on terms for growth and killing ----
    # Schmidt 2009 Eqs 2-3: growth and drug-induced killing each ramp
    # from 0 -> 1 via (1 - exp(-d*t)) to model the observed initial
    # lag phase before either process becomes effective.
    delay_growth <- 1 - exp(-dgs * t)
    delay_kill   <- 1 - exp(-dks * t)

    # ---- Logistic carrying-capacity factor (plat in (0,1)) ----
    ntot <- bact_susceptible + bact_persister
    plat <- 1 - ntot / nmax

    # ---- Emax drug-induced killing of the susceptible subpopulation ----
    kill <- kmax * central / (ec50 + central) * delay_kill

    # ---- Susceptible pool (Schmidt 2009 Eq 2) ----
    d/dt(bact_susceptible) <- ks  * plat * delay_growth * bact_susceptible -
                              kd  * bact_susceptible -
                              ksp * bact_susceptible +
                              kps * bact_persister -
                              kill * bact_susceptible

    # ---- Persister pool (Schmidt 2009 Eq 3) ----
    # Insusceptible to the antibiotic; only S<->P transfer and natural death apply.
    d/dt(bact_persister)   <- ksp * bact_susceptible -
                              kps * bact_persister -
                              kd  * bact_persister

    # ---- Antibiotic concentration in the MHB medium (first-order) ----
    # The central state numerically equals the broth concentration (no
    # explicit V because dose is delivered as a concentration into a
    # fixed-volume flask). For static linezolid experiments kdeg = 0
    # (the drug is stable); for dynamic experiments the user overrides
    # kdeg to log(2)/5 1/h to reproduce the syringe-dilution t1/2 ~ 5 h.
    d/dt(central) <- -kdeg * central

    # ---- Initial conditions (Schmidt 2009 Fig. 1 + F definition) ----
    bact_susceptible(0) <- ninit * fsusc
    bact_persister(0)   <- ninit * (1 - fsusc)

    # ---- Observation: log10 of total viable count ----
    Cc <- log10(ntot)
    Cc ~ add(addSd)
  })
}
