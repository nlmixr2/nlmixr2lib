Idkaidek_2011_ibuprofen_microgravity <- function() {
  description <- "One-compartment oral PK model with first-order absorption and linear elimination for a single 600 mg oral dose of ibuprofen in six healthy adult men studied in the SIMULATED-MICROGRAVITY arm (1-day antiorthostatic bed-rest position) of the Idkaidek 2011 crossover. Absorption is about three times faster than at normal gravity (Ka 2.24 vs 0.79 1/h, with inter-subject CV rising to 96%) and tmax is shorter, while exposure, elimination and bioavailability are unchanged - the paper's basis for concluding that no dose adjustment is needed in flight. Absorption rate constant from the paper's one-compartment fit; elimination rate constant from its noncompartmental analysis; apparent clearance and volume derived from the published AUC0-inf and Kel. Companion model at normal gravity: Idkaidek_2011_ibuprofen_normalGravity."
  reference <- "Idkaidek N, Arafat T. Effect of microgravity on the pharmacokinetics of ibuprofen in humans. J Clin Pharmacol. 2011;51(12):1685-1689. doi:10.1177/0091270010388652"
  vignette <- "Idkaidek_2011_ibuprofen"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    depot   = list(analyte = "ibuprofen", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "ibuprofen", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species        = "human",
    n_subjects     = 6,
    n_studies      = 1,
    age_range      = "18-45 years",
    sex_female_pct = 0,
    disease_state  = "Healthy adult male volunteers (Idkaidek 2011 'Human Participants'): body mass index 18.5-24.9 kg/m2, no clinically significant deviation from a normal medical condition on medical history, vital signs, physical examination, electrocardiogram or laboratory safety tests. All passed an antiorthostatic bed-rest tolerability test (able to eat and urinate while in the head-down position).",
    dose_range     = "Single 600 mg oral ibuprofen tablet with 240 mL water after a 10-hour overnight fast.",
    regions        = "Jordan (Jordan Center for Pharmaceutical Research, Al-Mowasah Hospital, Amman).",
    notes          = "Body position for this model file is the 1-day simulated-microgravity ANTIORTHOSTATIC BED-REST leg of a two-period sequential crossover with a 7-day washout; the same six men were also studied at normal gravity (see Idkaidek_2011_ibuprofen_normalGravity). Antiorthostatic (head-down-tilt) bed rest is the standard ground analogue for the cephalad fluid shift and altered gastric emptying / intestinal motility of spaceflight. Plasma was sampled at 0, 0.25, 0.5, 0.75, 1, 1.33, 1.66, 2, 2.5, 3, 4, 5, 6 and 8 h and assayed by a validated HPLC method with a linear range of 0.5-30 ug/mL. Saliva was also collected but ibuprofen was not detected in any saliva sample, so no salivary model exists."
  )

  ini({
    # Structural parameters. Idkaidek 2011 Table II, 'uG, Mean (CV%)' column.
    #
    # Ka and Kel are published directly. The apparent volume of distribution is
    # NOT published, so it is derived from two published quantities:
    #   CL/F = Dose / AUC0-inf = 600 mg / 120.77 ug*h/mL = 4.968 L/h
    #   V/F  = (CL/F) / Kel    = 4.968 / 0.36            = 13.80 L
    # (1 ug/mL = 1 mg/L, so 120.77 ug*h/mL = 120.77 mg*h/L and the quotient is
    # in L/h without further conversion.)
    lka <- log(2.24);  label("Absorption rate constant (1/h)")                            # Idkaidek 2011 Table II, Ka uG = 2.24 1/h (one-compartment data fit, Kinetica 2000; p = .03 vs 1G)
    lcl <- log(4.968); label("Apparent oral clearance CL/F (L/h)")                        # derived: 600 mg / AUC0-inf 120.77 ug*h/mL (Idkaidek 2011 Table II, uG)
    lvc <- log(13.80); label("Apparent central volume of distribution V/F (L)")           # derived: CL/F 4.968 L/h / Kel 0.36 1/h (Idkaidek 2011 Table II, uG)

    # IIV on Ka, from the published inter-subject CV%: omega^2 = log(CV^2 + 1)
    #   log(0.96^2 + 1) = 0.653158
    # The near-100% CV is a reported finding, not a transcription artefact:
    # Idkaidek 2011 Results notes that Ka "showed high intersubject variability,
    # reaching 96%" under the microgravity position.
    etalka ~ 0.653158                                                                  # Idkaidek 2011 Table II, Ka uG CV% = 96

    # IIV block on CL/F and V/F. Table II publishes a CV% for AUC0-inf (17%) and
    # for Kel (18%) but not for V/F, and CL/F and V/F are correlated, so the
    # block is computed from the individual concentration-time profiles in
    # Idkaidek 2011 Table I: per-subject AUC0-inf by linear trapezoid plus
    # Clast/lambda_z, lambda_z from log-linear regression on the 5, 6 and 8 h
    # points, then CL/F_i = 600/AUC0-inf_i and V/F_i = CL/F_i / lambda_z_i.
    # That reconstruction reproduces every uG cell of Table II exactly
    # (AUC0-inf 120.79 (17), Cmax 38.79 (21), tmax 1.14 (64), Kel 0.363 (18),
    # t1/2 1.96 (16)), so the derived block is on the same footing as the
    # published summaries. Resulting log-scale moments:
    #   var(log CL/F) = 0.034080  (geometric CV 19%, vs published AUC CV 17%)
    #   cov           = 0.004718  (correlation 0.377)
    #   var(log V/F)  = 0.004607
    # implying var(log Kel) = 0.029250, i.e. a Kel geometric CV of 17% against
    # the published 18%.
    etalcl + etalvc ~ c(0.034080,
                        0.004718, 0.004607)

    # The source fits each subject individually in Kinetica and reports no
    # residual-error estimate, so this is fixed at 0 per the standing policy for
    # unreported RUV (documented in the vignette Errata).
    propSd <- fixed(0); label("Proportional residual error (fraction; not reported in the source)")
  })

  model({
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)

    kel <- cl / vc

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
