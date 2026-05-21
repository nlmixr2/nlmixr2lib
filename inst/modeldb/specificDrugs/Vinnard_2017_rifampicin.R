Vinnard_2017_rifampicin <- function() {
  description <- "One-compartment population PK model for oral rifampicin in HIV/TB patients in Botswana (Vinnard 2017), with a Savic 2007 analytical transit-compartment absorption chain feeding a virtual depot, oral bioavailability fixed at 1, between-subject variability on CL, F, MTT, and the (non-integer) number of transit compartments NN, and inter-occasion variability on F across two sampling visits (pre-ART vs after approximately 4 weeks of ART)."
  reference <- paste(
    "Vinnard C, Ravimohan S, Tamuhla N, Pasipanodya J, Srivastava S,",
    "Modongo C, Zetola NM, Weissman D, Gumbo T, Bisson GP. (2017).",
    "Markers of gut dysfunction do not explain low rifampicin bioavailability",
    "in HIV-associated TB. J Antimicrob Chemother 72(7):2020-2027.",
    "doi:10.1093/jac/dkx111.",
    sep = " "
  )
  vignette <- "Vinnard_2017_rifampicin"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    OCC = list(
      description        = "Integer-valued sampling-visit occasion for inter-occasion-variability multiplexing on bioavailability F.",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = "Two occasions in the Vinnard 2017 cohort: OCC = 1 is the first pharmacokinetic study visit (5-28 days after starting anti-TB therapy, prior to ART initiation) and OCC = 2 is the second visit (approximately 4 weeks after ART initiation). Decomposed inside model() into binary indicators oc1 and oc2 that multiplex the two IOV etas on log-F. Vinnard 2017 Methods 'Data collection' paragraph describes the two-visit design.",
      source_name        = "OCC"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 40L,
    n_studies      = 1L,
    age_range      = "21 years and older; median 32 years (IQR 27-43 years)",
    age_median     = "32 years",
    weight_range   = "Not directly tabulated; inferred from weight-based rifampicin dosing 7.9-12.5 mg/kg at doses of 300, 450, 600, or 750 mg (median 9.7 mg/kg)",
    sex_female_pct = NA_real_,
    race_ethnicity = "Citizens of Botswana (sub-Saharan African); detailed ancestry not reported.",
    disease_state  = "HIV-infected adults newly diagnosed with pulmonary TB and initiating standard first-line antitubercular therapy (isoniazid + rifampicin + ethambutol + pyrazinamide) under directly-observed therapy.",
    dose_range     = "Oral rifampicin once daily at 300, 450, 600, or 750 mg per WHO weight-based dosing bands (1, 19, 17, and 3 individuals respectively); intensive phase 2 months, intermittent phase additional 4 months on rifampicin + isoniazid. 7.9-12.5 mg/kg, median 9.7 mg/kg.",
    regions        = "Botswana (Gaborone; 22 public clinics and Princess Marina Hospital).",
    co_medication  = "First-line antitubercular fixed-dose combinations (isoniazid + ethambutol + pyrazinamide) at visit 1; tenofovir + emtricitabine + efavirenz-based ART at visit 2.",
    notes          = "40 HIV/TB patients completed visit 1 and 24 completed visit 2 (approximately 4 weeks after ART initiation). Median CD4 T-cell count 238 cells/uL (IQR 105-339). Exclusion criteria included pregnancy, renal insufficiency (CrCl < 50 mL/min), and hepatic dysfunction (ALT or AST > 3x ULN). Pharmacokinetic sampling at 0.3, 0.9, 2.2, 4.5, and 8 h post-dose. The paper evaluated I-FABP, sCD14, %CD38+DR+CD8+, IL-6, CD4 T-cell count, and HIV viral load as covariates on bioavailability and found none significant (Table 2)."
  )

  ini({
    # Structural PK parameters - final estimates from Vinnard 2017 Table 1
    # 'Final model parameter estimates for adult population pharmacokinetics
    # of rifampicin in HIV/TB patients'. The paper reports an oral
    # one-compartment model with first-order elimination and a Savic 2007
    # transit-compartment absorption chain (no separately estimated ka;
    # the chain itself delivers to the central compartment via the
    # KTR = (NN + 1) / MTT rate). All values are point estimates from the
    # 'Typical value (% RSE)' column.
    lcl  <- log(16.36); label("Apparent oral clearance CL/F (L/h)")                                       # Vinnard 2017 Table 1: oral apparent clearance (CL/F) = 16.36 L/h (RSE 6.27%)
    lvc  <- log(52.90); label("Apparent oral central volume of distribution V/F (L)")                     # Vinnard 2017 Table 1: oral apparent volume of distribution (V/F) = 52.90 L (RSE 5.10%)
    lmtt <- log(0.97);  label("Mean absorption transit time MTT (h)")                                     # Vinnard 2017 Table 1: mean transit time within compartments = 0.97 h (RSE 6.80%)
    lnn  <- log(4.36);  label("Number of absorption transit compartments NN (continuous, dimensionless)") # Vinnard 2017 Table 1: number of absorption transit compartments = 4.36 (RSE 9.01%)
    lfdepot <- fixed(log(1)); label("Oral bioavailability F (fixed at 1)")                                # Vinnard 2017 Table 1: oral bioavailability (F) = 1 (fixed)

    # Inter-individual variability (between-subject). Source reports BSV as CV%
    # in Table 1; converted to log-normal variance via omega^2 = log(1 + CV^2).
    etalcl     ~ 0.0998  # Vinnard 2017 Table 1: BSV CL/F = 32.4% CV; omega^2 = log(1 + 0.324^2) = 0.0998
    etalfdepot ~ 0.0436  # Vinnard 2017 Table 1: BSV F = 21.1% CV; omega^2 = log(1 + 0.211^2) = 0.0436
    etalnn     ~ 0.5517  # Vinnard 2017 Table 1: BSV NN = 85.8% CV; omega^2 = log(1 + 0.858^2) = 0.5517
    etalmtt    ~ 0.2702  # Vinnard 2017 Table 1: BSV MTT = 55.7% CV; omega^2 = log(1 + 0.557^2) = 0.2702

    # Inter-occasion variability on F across the two sampling visits
    # (pre-ART = OCC 1, post-ART approximately 4 weeks = OCC 2). The source
    # reports IOV F = 9.2% CV in Table 1. nlmixr2 has no NONMEM-style 'SAME'
    # shortcut so the second occasion's variance is fixed equal to the first
    # (matching the standard NONMEM $OMEGA BLOCK(1) + SAME pattern, e.g.
    # Wilkins_2008_rifampicin.R).
    etaiov_fdepot_1 ~ 0.00843        # Vinnard 2017 Table 1: IOV F = 9.2% CV; omega^2 = log(1 + 0.092^2) = 0.00843 (occasion 1, pre-ART)
    etaiov_fdepot_2 ~ fix(0.00843)   # occasion 2 (post-ART) - variance fixed equal to occasion 1 (single IOV variance shared across occasions in source)

    # Combined additive + proportional residual error (Vinnard 2017 Table 1).
    # Concentrations are in mg/L (= ug/mL).
    addSd  <- 0.16; label("Additive residual error (mg/L)")              # Vinnard 2017 Table 1: additive error (SD, mg/L) = 0.16 (RSE 0.05%)
    propSd <- 0.37; label("Proportional residual error (fraction, CV)")  # Vinnard 2017 Table 1: proportional error (% CV) = 0.37 (RSE 11.6%)
  })

  model({
    # 1. Decompose the integer-valued occasion column into binary indicators
    #    for IOV multiplexing on log-F. Records with OCC outside {1, 2} get
    #    zero IOV contribution (consistent with the Wilkins 2008 pattern).
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    iov_fdepot <- oc1 * etaiov_fdepot_1 + oc2 * etaiov_fdepot_2

    # 2. Individual PK parameters. No covariate effects on the structural
    #    model: Vinnard 2017 Table 2 explicitly evaluated I-FABP, sCD14,
    #    %CD38+DR+CD8+, IL-6, CD4 T-cell count, and HIV viral load on F and
    #    none significantly reduced the objective function value.
    cl     <- exp(lcl  + etalcl)
    vc     <- exp(lvc)                       # vc carries no IIV in the published model (Vinnard 2017 Table 1)
    mtt    <- exp(lmtt + etalmtt)
    nn     <- exp(lnn  + etalnn)
    fdepot <- exp(lfdepot + etalfdepot + iov_fdepot)

    kel <- cl / vc

    # 3. Absorption: rxode2 transit() implements the Savic 2007 analytical
    #    gamma-PDF input rate for non-integer NN. The source model has no
    #    separately estimated first-order absorption constant (Vinnard 2017
    #    Table 1 reports only NN and MTT for absorption); following the
    #    Wilkins 2008 / Tikiso 2021 / vanderWalt 2013 dapagliflozin pattern,
    #    we set ka >> (NN + 1)/MTT (here ka = 60 /h, so half-life in depot
    #    ~ 0.012 h, an order of magnitude faster than KTR ~ 5.5/h here) to
    #    collapse the depot's exponential tail and let the transit() output
    #    drive the central compartment input rate directly. f(depot) = 0
    #    suppresses the dose-event bolus into depot; the transit() function
    #    reads the raw dose amount from podo(depot) regardless of f(depot)
    #    so the full dose enters the absorption process through the
    #    analytical chain. Bioavailability fdepot enters via the bio
    #    argument of transit().
    ka <- 60

    d/dt(depot)   <- transit(nn, mtt, fdepot) - ka * depot
    d/dt(central) <- ka * depot - kel * central

    f(depot) <- 0

    # 4. Plasma rifampicin concentration. Dose mg / volume L = mg/L.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
