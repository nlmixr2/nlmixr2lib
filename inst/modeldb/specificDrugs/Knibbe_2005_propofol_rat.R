Knibbe_2005_propofol_rat <- function() {
  description <- paste(
    "Preclinical (rat). Two-compartment intravenous population PK",
    "model for propofol in male Wistar rats following a single 30",
    "mg/kg bolus delivered over 5 min, as reported in Table 3",
    "(column 'Observed in the rat (250 g)') of Knibbe 2005. The",
    "underlying NONMEM fit was performed by Knibbe et al. (reference",
    "11 of the paper) on 19 whole-blood samples from each of 22",
    "chronically instrumented rats; Knibbe 2005 reproduces those",
    "rat point estimates and uses them as the species anchor for an",
    "allometric scaling to humans (see the companion model file",
    "Knibbe_2005_propofol_human.R, which carries the human-projected",
    "parameters from Table 3 column 'Scaled for humans (70 kg)').",
    "Log-normal inter-individual variability on CL, V1, Q, V2 and a",
    "constant-CV proportional intra-individual residual error model.",
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
    species        = "rat (male Wistar)",
    n_subjects     = 22L,
    n_studies      = 1L,
    age_range      = "adult",
    weight_range   = "0.25-0.30 kg (reference weight 0.25 kg)",
    sex_female_pct = 0,
    disease_state  = paste(
      "Healthy chronically instrumented animals (jugular catheter for",
      "drug administration, arterial catheter for sampling);",
      "non-ventilated, experimental setting. Body temperature",
      "36.5-37.5 C; renal, hepatic, and cardiac function normal",
      "(Knibbe 2005 Table 1)."
    ),
    dose_range     = paste(
      "Single propofol bolus of 30 mg/kg administered over 5 minutes",
      "through the jugular catheter (Methods, Animals and patients).",
      "19 whole-blood samples per rat were drawn for HPLC-fluorescence",
      "assay over the concentration range 0.05-40 mg/L."
    ),
    regions        = "the Netherlands",
    notes          = paste(
      "The rat parameter values were originally fit in reference 11",
      "of Knibbe 2005 (the Knibbe et al. preclinical propofol popPK",
      "publication). Knibbe 2005 reproduces them in Table 3 as the",
      "rat-side anchor of the cross-species allometric analysis.",
      "Inter-individual variability is reported as %CV in Table 3",
      "with the footnote 'interindividual variability is expressed",
      "as %CV and equals the square root of the exponential variance",
      "of eta minus 1', i.e. CV = sqrt(exp(omega^2) - 1) and",
      "omega^2 = log(CV^2 + 1). Intra-individual variability is the",
      "proportional residual SD on the linear scale and is reported",
      "as %CV (= 19.9%).",
      sep = " "
    )
  )

  ini({
    # Knibbe 2005 Table 3, column "Observed in the rat (250 g)"
    # (page 708). Two-compartment IV PK with first-order elimination
    # from the central compartment, as fit in reference 11 of the
    # paper.
    lcl <- log(0.0261) ; label("Clearance (L/min)")                          # Knibbe 2005 Table 3 rat: CL = 0.0261 L/min (SE 0.00205)
    lvc <- log(0.0811) ; label("Central volume of distribution V1 (L)")      # Knibbe 2005 Table 3 rat: V1 = 0.0811 L (SE 0.00544)
    lq  <- log(0.0227) ; label("Inter-compartmental clearance Q (L/min)")    # Knibbe 2005 Table 3 rat: Q  = 0.0227 L/min (SE 0.00325)
    lvp <- log(0.291)  ; label("Peripheral volume of distribution V2 (L)")   # Knibbe 2005 Table 3 rat: V2 = 0.291 L (SE 0.0067)

    # Inter-individual variability. Reported as %CV in Table 3;
    # internal log-normal variance is omega^2 = log(CV^2 + 1).
    etalcl ~ 0.10937   # Knibbe 2005 Table 3 rat: IIV CL 34%, omega^2 = log(0.34^2 + 1)
    etalvc ~ 0.02226   # Knibbe 2005 Table 3 rat: IIV V1 15%, omega^2 = log(0.15^2 + 1)
    etalq  ~ 0.06541   # Knibbe 2005 Table 3 rat: IIV Q  26%, omega^2 = log(0.26^2 + 1)
    etalvp ~ 0.05154   # Knibbe 2005 Table 3 rat: IIV V2 23%, omega^2 = log(0.23^2 + 1)

    # Intra-individual residual error. Constant-CV proportional model
    # (Knibbe 2005 Eq. 4): c_ij = c_pred,ij * (1 + eps), with eps ~
    # N(0, sigma^2). The Table 3 footnote defines intra-individual
    # variability as CV = sqrt(sigma^2); propSd = 0.199 reproduces
    # the reported 19.9%.
    propSd <- 0.199 ; label("Proportional residual error (fraction)")        # Knibbe 2005 Table 3 rat: intra-individual variability 19.9%
  })

  model({
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q  <- exp(lq  + etalq)
    vp <- exp(lvp + etalvp)

    # Two-compartment IV PK with first-order elimination from the
    # central compartment. Doses land directly in `central` via the
    # cmt column of the user data set; rat doses in mg are obtained
    # as 30 mg/kg * body weight (kg).
    d/dt(central)     <-  q / vp * peripheral1 - (cl + q) / vc * central
    d/dt(peripheral1) <-  q / vc * central     -  q       / vp * peripheral1

    # Whole-blood propofol concentration. Dose units mg, Vc units L
    # -> Cc units mg/L, matching the HPLC-fluorescence assay range
    # 0.05-40 mg/L (Methods).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
