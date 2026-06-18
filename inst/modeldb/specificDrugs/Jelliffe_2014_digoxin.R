Jelliffe_2014_digoxin <- function() {
  description <- "Two-compartment population PK/PD model of digoxin in adults with first-order oral absorption, creatinine-clearance-dependent renal elimination, and a peripheral effect compartment normalized per body weight (Jelliffe et al. 2014, Ther Drug Monit; structural parameters carried from Reuning et al. 1973)."
  reference <- "Jelliffe R, Milman M, Schumitzky A, Bayard D, Van Guilder M. A Two-Compartment Population Pharmacokinetic-Pharmacodynamic Model of Digoxin in Adults, with Implications for Dosage. Ther Drug Monit. 2014 June;36(3):387-393. doi:10.1097/FTD.0000000000000023. PMCID: PMC4040260."
  vignette <- "Jelliffe_2014_digoxin"
  units <- list(time = "hour", dosing = "ug", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Per-kg model. The central volume Vc and the peripheral state amount-per-kg both scale linearly with WT (allometric exponent 1). Reuning 1973 used 70 kg as the assumed average adult body weight when computing Vc = 110 L (1.5714 L/kg).",
      source_name        = "WT"
    ),
    CRCL = list(
      description        = "Estimated creatinine clearance, BSA-normalized.",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Linear additive effect on the renal component of elimination: kel = knr + kr * CRCL. The paper estimated CRCL using the method of Jelliffe 2002 (Am J Nephrol 22:320-324). At a CRCL of 100 mL/min/1.73 m^2 the renal rate constant equals 0.0451 1/hr, equivalent to kr = 0.000451 1/hr per mL/min/1.73 m^2 (Methods, p3).",
      source_name        = "CCr"
    )
  )

  population <- list(
    n_subjects     = NA_integer_,
    n_studies      = 3,
    age_range      = "Adults",
    weight_range   = "Reference adult 70 kg (Reuning 1973 assumption)",
    sex_female_pct = NA_real_,
    race_ethnicity = NA_character_,
    disease_state  = "Adults requiring digoxin therapy (e.g., congestive heart failure with sinus rhythm; atrial fibrillation or flutter); model spans normal to anephric renal function via the CRCL covariate.",
    dose_range     = "Oral and intravenous; the paper illustrates loading doses of 0.875-1.149 mg given in 3 split doses 6 h apart followed by maintenance doses of 125-345 ug/day, but the structural model itself is dose-agnostic.",
    regions        = "Not specified",
    notes          = "The structural parameters (Vc, Knr, Kr, Kcp, Kpc) are carried from the three normal-renal-function studies pooled by Reuning, Sams, and Notari (J Clin Pharmacol 1973;13:127-141, Table 1, p 129; Vc = 110 L for a 70 kg adult, total elimination rate 0.0747 1/hr split 39% nonrenal / 61% renal). The absorption rate Ka and oral bioavailability F = 0.65 are stated by Jelliffe 2014 (Methods, p2). Per-parameter variability was assumed at 20% CV (Methods, p3). Jelliffe 2014 then converted the continuous parameter distributions into a 64-point discrete distribution via the maximum-entropy method of Milman, Jiang, and Jelliffe (Comput Biol Med 2001;31:197-214) for use in the BestDose / USC RightDose adaptive-control software."
  )

  ini({
    # Mean parameter values from Table 2 (corroborated by the Methods narrative on p2-3 and the
    # discrete support points of Table 1). Per-kg / per-CRCL forms are reproduced verbatim;
    # the per-WT scaling and Kel = Knr + Kr * CRCL composition are applied inside model().
    lka  <- log(0.6093);   label("First-order oral absorption rate (1/hr)")                                      # Table 2: Ka mean = 0.6093 (Table 1 support points 0.4874 and 0.7312)
    lvc  <- log(1.5714);   label("Central (serum) volume per kg body weight (L/kg)")                              # Table 2: Vc mean = 1.5714 L/kg (Methods p3: Vc = 110 L for a 70 kg adult)
    lknr <- log(0.0288);   label("Nonrenal first-order elimination rate constant (1/hr)")                         # Table 2: Knr mean = 0.0288 (Methods p3: 39% of overall 0.0747 1/hr)
    lkr  <- log(0.000451); label("Renal first-order elimination rate constant per CRCL unit (1/hr per mL/min/1.73 m^2)") # Table 2: Kr mean = 0.000451 (Methods p3: 0.0451 1/hr at CRCL = 100)
    lkcp <- log(0.56);     label("Central-to-peripheral first-order rate constant (1/hr)")                        # Table 2: Kcp mean = 0.56 (Methods p3, from Reuning 1973)
    lkpc <- log(0.15);     label("Peripheral-to-central first-order rate constant (1/hr)")                        # Table 2: Kpc mean = 0.15 (Methods p3, from Reuning 1973)
    lfdepot <- fixed(log(0.65)); label("Oral bioavailability (fraction)")                                          # Methods p2: oral bioavailability assumed 0.65

    # IIV: paper assumes 20% CV per parameter (Methods p3: "The variability in the parameter
    # distributions was assumed to be 20 percent"). Log-normal variance = log(1 + CV^2).
    etalka  ~ log(1 + 0.20^2)  # 20% CV on Ka  (Methods p3 / Table 2 SD 0.122   ~ 20% of 0.6093)
    etalvc  ~ log(1 + 0.20^2)  # 20% CV on Vc  (Methods p3 / Table 2 SD 0.314   = 20% of 1.5714)
    etalknr ~ log(1 + 0.20^2)  # 20% CV on Knr (Methods p3 / Table 2 SD 0.0058  ~ 20% of 0.0288)
    etalkr  ~ log(1 + 0.20^2)  # 20% CV on Kr  (Methods p3 / Table 2 SD 0.0001  ~ 22% of 0.000451)
    etalkcp ~ log(1 + 0.20^2)  # 20% CV on Kcp (Methods p3 / Table 2 SD 0.112   = 20% of 0.56)
    etalkpc ~ log(1 + 0.20^2)  # 20% CV on Kpc (Methods p3 / Table 2 SD 0.03    = 20% of 0.15)

    # Residual error not reported by Jelliffe 2014; a small fixed proportional term is supplied
    # so the observation Cc has a valid residual model in nlmixr2 syntax. Set to 1% to keep
    # typical-value simulations effectively deterministic.
    propSd <- fixed(0.01); label("Proportional residual error (placeholder; not reported in source)")  # not in source -- see vignette Assumptions and deviations
  })

  model({
    # Individual PK parameters
    ka  <- exp(lka  + etalka)
    vc  <- exp(lvc  + etalvc) * WT
    knr <- exp(lknr + etalknr)
    kr  <- exp(lkr  + etalkr)
    kcp <- exp(lkcp + etalkcp)
    kpc <- exp(lkpc + etalkpc)

    # Total central elimination rate constant: kel = knr + kr * CRCL
    kel <- knr + kr * CRCL

    # Mass-action ODEs (depot, central, peripheral1 carry amounts in ug)
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - kcp * central + kpc * peripheral1
    d/dt(peripheral1) <-                                kcp * central - kpc * peripheral1

    # Oral bioavailability applied to the depot compartment
    f(depot) <- exp(lfdepot)

    # Observation: serum concentration in ng/mL (dose ug / vc L = ug/L = ng/mL)
    Cc <- central / vc
    Cc ~ prop(propSd)

    # Peripheral compartment exposed in the paper's normalization (amount per kg, ug/kg)
    Cp_ugkg <- peripheral1 / WT
  })
}
