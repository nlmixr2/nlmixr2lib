Knibbe_2005_propofol_human <- function() {
  description <- paste(
    "Two-compartment intravenous population PK model for propofol",
    "in a 70 kg adult human, projected from male Wistar rat (0.25",
    "kg) parameters via the allometric power model with literature",
    "exponents 0.75 for clearances and 1 for volumes. Parameter",
    "values are taken from Knibbe 2005 Table 3 (column 'Scaled for",
    "humans (70 kg)'); inter- and intra-individual variability are",
    "inherited from the rat fit (Table 3, column 'Observed in the",
    "rat (250 g)') per the Methods text 'these human scaled",
    "pharmacokinetic parameters, together with ... intra- and",
    "interindividual variabilities estimated in the rat were used",
    "to simulate propofol concentrations'. The companion file",
    "Knibbe_2005_propofol_rat.R carries the rat-side parameters",
    "used as the scaling anchor. Knibbe 2005 demonstrated that",
    "concentrations simulated from this scaled-human model",
    "agreed (r^2 = 0.83, P < 0.0001) with concentrations observed",
    "in long-term-sedated critically ill patients (Figure 2).",
    sep = " "
  )
  reference <- paste(
    "Knibbe CAJ, Zuideveld KP, Aarts LPHJ, Kuks PFM, Danhof M.",
    "(2005). Allometric relationships between the pharmacokinetics",
    "of propofol in rats, children and adults.",
    "British Journal of Clinical Pharmacology 59(6):705-711.",
    "doi:10.1111/j.1365-2125.2005.02239.x.",
    sep = " "
  )
  vignette <- "Knibbe_2005_propofol"
  units    <- list(time = "minute", dosing = "mg", concentration = "mg/L")

  covariateData <- list()

  population <- list(
    species        = "human (allometric projection from male Wistar rat)",
    n_subjects     = NA_integer_,
    n_studies      = 0L,
    age_range      = "adult (70 kg reference subject)",
    weight_range   = "70 kg (reference subject only; scaled point estimate)",
    sex_female_pct = NA_real_,
    disease_state  = paste(
      "Reference 70 kg adult, intended in Knibbe 2005 as the",
      "starting point for first-in-human dosing simulations and",
      "compared in Figure 2 against observed concentrations in",
      "long-term-sedated critically ill patients (52-79 y, 70-96",
      "kg; Table 1 column 'Critically ill adult patients')."
    ),
    dose_range     = paste(
      "Validation dosing in Knibbe 2005 used the propofol regimens",
      "given to the long-term-sedated critically ill cohort (4",
      "mg/kg/h, with daily infusion-rate adjustments per Table 1).",
      "The model itself is dose-route agnostic for IV propofol",
      "delivery to the central compartment."
    ),
    regions        = "the Netherlands",
    notes          = paste(
      "Parameter values come from applying the allometric power",
      "model Y_human = Y_rat * (BW_human / BW_rat)^b to the",
      "rat point estimates of Table 3, using the canonical",
      "literature exponents b = 0.75 for CL and Q and b = 1 for",
      "V1 and V2 (Table 3 footnote). The table values can be",
      "recovered (within reporting precision) at BW_rat ~= 0.275",
      "kg and BW_human = 70 kg. Inter- and intra-individual",
      "variability are inherited unchanged from the rat fit because",
      "the simulation in Knibbe 2005 used the rat variability",
      "estimates (Methods, 'Scaling pharmacokinetic parameters ...",
      "and simulations in critically ill patients').",
      sep = " "
    )
  )

  ini({
    # Knibbe 2005 Table 3, column "Scaled for humans (70 kg)"
    # (page 708). Derived by Y_human = Y_rat * (BW_human / BW_rat)^b
    # with b = 0.75 for CL and Q and b = 1 for V1 and V2 (Table 3
    # footnote). Two-compartment IV PK with first-order elimination
    # from the central compartment.
    lcl <- log(1.63) ; label("Clearance (L/min)")                            # Knibbe 2005 Table 3 scaled-human: CL = 1.63 L/min
    lvc <- log(20.6) ; label("Central volume of distribution V1 (L)")        # Knibbe 2005 Table 3 scaled-human: V1 = 20.6 L
    lq  <- log(1.45) ; label("Inter-compartmental clearance Q (L/min)")      # Knibbe 2005 Table 3 scaled-human: Q  = 1.45 L/min
    lvp <- log(71.9) ; label("Peripheral volume of distribution V2 (L)")     # Knibbe 2005 Table 3 scaled-human: V2 = 71.9 L

    # Inter-individual variability inherited from the rat fit (Methods,
    # "intra- and interindividual variabilities estimated in the rat
    # were used"). Reported as %CV; internal log-normal variance is
    # omega^2 = log(CV^2 + 1).
    etalcl ~ 0.10937   # Knibbe 2005 Table 3 rat: IIV CL 34%, omega^2 = log(0.34^2 + 1)
    etalvc ~ 0.02226   # Knibbe 2005 Table 3 rat: IIV V1 15%, omega^2 = log(0.15^2 + 1)
    etalq  ~ 0.06541   # Knibbe 2005 Table 3 rat: IIV Q  26%, omega^2 = log(0.26^2 + 1)
    etalvp ~ 0.05154   # Knibbe 2005 Table 3 rat: IIV V2 23%, omega^2 = log(0.23^2 + 1)

    # Intra-individual residual error, also inherited from the rat
    # fit (Knibbe 2005 Eq. 4 constant-CV proportional model;
    # propSd = 0.199 reproduces the reported 19.9% rat residual CV).
    propSd <- 0.199 ; label("Proportional residual error (fraction)")        # Knibbe 2005 Table 3 rat: intra-individual variability 19.9%
  })

  model({
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q  <- exp(lq  + etalq)
    vp <- exp(lvp + etalvp)

    # Two-compartment IV PK with first-order elimination from the
    # central compartment. Doses land directly in `central` via the
    # cmt column of the user data set.
    d/dt(central)     <-  q / vp * peripheral1 - (cl + q) / vc * central
    d/dt(peripheral1) <-  q / vc * central     -  q       / vp * peripheral1

    # Whole-blood propofol concentration. Dose units mg, Vc units L
    # -> Cc units mg/L.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
