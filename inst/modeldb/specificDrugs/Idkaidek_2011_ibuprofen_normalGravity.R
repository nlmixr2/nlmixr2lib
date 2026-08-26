Idkaidek_2011_ibuprofen_normalGravity <- function() {
  description <- "One-compartment oral PK model with first-order absorption and linear elimination for a single 600 mg oral dose of ibuprofen in six healthy adult men studied in the NORMAL-GRAVITY (1G, ambulatory) arm of the Idkaidek 2011 simulated-microgravity crossover. Absorption rate constant from the paper's one-compartment fit (Ka 0.79 1/h); elimination rate constant from the paper's noncompartmental analysis (Kel 0.39 1/h); apparent clearance and volume derived from the published AUC0-inf and Kel. Its companion model, Idkaidek_2011_ibuprofen_microgravity, holds the same six subjects under the antiorthostatic bed-rest (simulated microgravity) position, where absorption is roughly three times faster while exposure is unchanged."
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
    disease_state  = "Healthy adult male volunteers (Idkaidek 2011 'Human Participants'): body mass index 18.5-24.9 kg/m2, no clinically significant deviation from a normal medical condition on medical history, vital signs, physical examination, electrocardiogram or laboratory safety tests.",
    dose_range     = "Single 600 mg oral ibuprofen tablet with 240 mL water after a 10-hour overnight fast.",
    regions        = "Jordan (Jordan Center for Pharmaceutical Research, Al-Mowasah Hospital, Amman).",
    notes          = "Body position for this model file is the NORMAL-GRAVITY (1G) leg of a two-period sequential crossover with a 7-day washout; the same six men were also studied in the simulated-microgravity antiorthostatic bed-rest position (see Idkaidek_2011_ibuprofen_microgravity). Sample size of 6 was based on a previous ibuprofen intersubject PK variability of 25% giving 80% power at 95% confidence. Plasma was sampled at 0, 0.25, 0.5, 0.75, 1, 1.33, 1.66, 2, 2.5, 3, 4, 5, 6 and 8 h and assayed by a validated HPLC method with a linear range of 0.5-30 ug/mL. Saliva was also collected but ibuprofen was not detected in any saliva sample, so no salivary model exists."
  )

  ini({
    # Structural parameters. Idkaidek 2011 Table II, '1G, Mean (CV%)' column.
    #
    # Ka and Kel are published directly. The apparent volume of distribution is
    # NOT published, so it is derived from two published quantities:
    #   CL/F = Dose / AUC0-inf = 600 mg / 128.88 ug*h/mL = 4.655 L/h
    #   V/F  = (CL/F) / Kel    = 4.655 / 0.39            = 11.94 L
    # (1 ug/mL = 1 mg/L, so 128.88 ug*h/mL = 128.88 mg*h/L and the quotient is
    # in L/h without further conversion.)
    lka <- log(0.79);  label("Absorption rate constant (1/h)")                            # Idkaidek 2011 Table II, Ka 1G = 0.79 1/h (one-compartment data fit, Kinetica 2000)
    lcl <- log(4.655); label("Apparent oral clearance CL/F (L/h)")                        # derived: 600 mg / AUC0-inf 128.88 ug*h/mL (Idkaidek 2011 Table II, 1G)
    lvc <- log(11.94); label("Apparent central volume of distribution V/F (L)")           # derived: CL/F 4.655 L/h / Kel 0.39 1/h (Idkaidek 2011 Table II, 1G)

    # IIV on Ka, from the published inter-subject CV%: omega^2 = log(CV^2 + 1)
    #   log(0.35^2 + 1) = 0.115558
    etalka ~ 0.115558                                                                  # Idkaidek 2011 Table II, Ka 1G CV% = 35

    # IIV block on CL/F and V/F. Table II publishes a CV% for AUC0-inf (18%) and
    # for Kel (13%) but not for V/F, and CL/F and V/F are strongly correlated,
    # so the block is computed from the individual concentration-time profiles
    # in Idkaidek 2011 Table I: per-subject AUC0-inf by linear trapezoid plus
    # Clast/lambda_z, lambda_z from log-linear regression on the 5, 6 and 8 h
    # points, then CL/F_i = 600/AUC0-inf_i and V/F_i = CL/F_i / lambda_z_i.
    # That reconstruction reproduces every 1G cell of Table II exactly
    # (AUC0-inf 128.88 (18), Cmax 40.12 (20), tmax 1.37 (46), Kel 0.387 (13),
    # t1/2 1.82 (15)), so the derived block is on the same footing as the
    # published summaries. Resulting log-scale moments:
    #   var(log CL/F) = 0.037514  (geometric CV 20%, vs published AUC CV 18%)
    #   cov           = 0.015534  (correlation 0.700)
    #   var(log V/F)  = 0.013137
    # implying var(log Kel) = 0.019582, i.e. a Kel geometric CV of 14% against
    # the published 13%.
    etalcl + etalvc ~ c(0.037514,
                        0.015534, 0.013137)

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
