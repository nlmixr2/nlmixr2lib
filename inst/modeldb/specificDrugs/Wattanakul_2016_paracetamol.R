Wattanakul_2016_paracetamol <- function() {
  description <- paste(
    "Two-compartment population PK model for paracetamol (acetaminophen)",
    "administered as a single 600 mg dose by either intramuscular injection",
    "(zero-order absorption over DUR_IM) or oral syrup (first-order absorption",
    "with rate constant ka) in 21 adult Thai patients with uncomplicated",
    "Plasmodium falciparum malaria and fever > 38 C (Wattanakul 2016).",
    "Intramuscular bioavailability is fixed to F_IM = 1; the relative oral",
    "bioavailability is F_PO = 0.844 (95% CI 0.682-0.951). The depot",
    "compartment carries oral doses (f(depot) = F_PO) while intramuscular",
    "doses target central with rate = -2 to invoke the modeled",
    "dur(central) = DUR_IM. No covariates were retained: allometric scaling",
    "on body weight did not improve the fit and a stepwise covariate search",
    "(age, AST, ALT, bilirubin, BUN, creatinine, sex, hemoglobin, parasitaemia,",
    "systolic BP, temperature) found no significant effect at p < 0.05.",
    "Inter-individual variability for V_C and DUR_IM was estimated below 1%",
    "CV and fixed to zero in the source paper without changing the OFV;",
    "this model omits the corresponding etas accordingly."
  )
  reference <- paste(
    "Wattanakul T, Teerapong P, Plewes K, Newton PN, Chierakul W, Silamut K,",
    "Chotivanich K, Ruengweerayut R, White NJ, Dondorp AM, Tarning J (2016).",
    "Pharmacokinetic properties of intramuscular versus oral syrup paracetamol",
    "in Plasmodium falciparum malaria.",
    "Malaria Journal 15:244.",
    "doi:10.1186/s12936-016-1283-9.",
    sep = " "
  )
  vignette <- "Wattanakul_2016_paracetamol"
  units    <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    # Wattanakul 2016 Results, Pharmacokinetics paragraph 2: "Allometric scaling
    # of pharmacokinetic parameters did not improve model fit significantly.
    # Thus, body weight was not incorporated into the final model. The stepwise
    # covariate search showed no significant relationships in this population."
    # No covariate columns are required to simulate from this model.
  )

  population <- list(
    species        = "human",
    n_subjects     = 21L,
    n_studies      = 1L,
    age_range      = "15 to 54 years (median 25; IQR 22-37)",
    weight_range   = "47 to 70 kg (median 58; IQR 55-63)",
    sex_female_pct = 100 * (1 - 19 / 21),
    race_ethnicity = "Thai (Mae Sot Hospital, Tak Province, Thailand)",
    disease_state  = paste(
      "Adults with slide-confirmed uncomplicated Plasmodium falciparum",
      "malaria and aural temperature > 38 C (range 38.1-41.2 C);",
      "median parasitaemia ~47,500 parasites/uL.",
      "All patients also received intravenous artesunate plus oral",
      "doxycycline antimalarial therapy per the local protocol."
    ),
    dose_range     = paste(
      "Single 600 mg dose by either intramuscular injection",
      "(300 mg / 2 mL Partamol, two 2 mL injections in the anterior thigh)",
      "or oral syrup (Tylenol, given orally or via nasogastric tube),",
      "in a randomized open-label two-treatment crossover design with",
      "alternate-route dosing on day 1."
    ),
    sampling       = paste(
      "Pre-dose plus 0.5, 1.0, 1.5, 2, 3, 4, 6, 8, 10 and 12 h after each",
      "dose; 363 quantifiable paracetamol concentrations were included in",
      "the pharmacokinetic analysis."
    ),
    regions        = "Thailand (Mae Sot Hospital, Tak Province)",
    notes          = paste(
      "Baseline demographics from Wattanakul 2016 Table 1; study conducted",
      "May to June 2001. 19/21 (90%) were male. Patients on potentially",
      "interacting drugs or who had taken paracetamol within 12 h were",
      "excluded; severe-malaria criteria were also exclusionary at",
      "enrolment, though one patient in Group 2 was subsequently",
      "reclassified as severe malaria."
    )
  )

  ini({
    # Structural parameters: Wattanakul 2016 Table 2, "Population estimate" column.

    # Bioavailability. The depot compartment carries oral doses and the central
    # compartment carries IM doses. lfdepot encodes log(F_PO); F_IM is fixed at
    # 1 in the source paper and is left as the rxode2 default f(central) = 1.
    lfdepot <- log(0.844)
    label("Relative oral syrup bioavailability vs intramuscular (F_PO, fraction)")  # Table 2: F_PO = 0.844 (95% CI 0.682-0.951; 8.4% RSE)

    # First-order oral absorption rate.
    lka     <- log(4.15)
    label("Oral syrup absorption rate constant (ka_PO, 1/h)")  # Table 2: ka_PO = 4.15 1/h (95% CI 1.95-9.73; 44.5% RSE)

    # Zero-order IM input duration. Wattanakul 2016 Table 2 reports IIV (%CV)
    # < 1% for DUR_IM and so the source paper fixed the corresponding eta to
    # zero without changing the OFV; this model encodes that by omitting
    # etaldur_im. The DUR_IM population estimate itself is estimated (not
    # fixed) per Table 2.
    ldur_im <- log(0.689)
    label("IM zero-order input duration (DUR_IM, h)")  # Table 2: DUR_IM = 0.689 h (95% CI 0.621-0.784; 6.2% RSE)

    # Disposition parameters.
    lcl     <- log(10.7)
    label("Apparent elimination clearance (CL, L/h)")  # Table 2: CL = 10.7 L/h (95% CI 7.35-14.7; 16.9% RSE)
    lvc     <- log(45.5)
    label("Apparent central volume of distribution (V_C, L)")  # Table 2: V_C = 45.5 L (95% CI 36.7-51.5; 8.5% RSE)
    lq      <- log(10.3)
    label("Apparent intercompartmental clearance (Q, L/h)")  # Table 2: Q = 10.3 L/h (95% CI 4.80-20.1; 36.8% RSE)
    lvp     <- log(11.3)
    label("Apparent peripheral volume of distribution (V_P, L)")  # Table 2: V_P = 11.3 L (95% CI 5.01-29.0; 42.7% RSE)

    # Inter-individual variability. Wattanakul 2016 Table 2 reports %CV computed
    # as 100 * sqrt(exp(omega^2) - 1), so the internal-scale variance is
    # omega^2 = log(1 + (CV/100)^2) for each estimated eta. IIVs for V_C and
    # DUR_IM were estimated below 1% CV and fixed to zero in the source paper
    # without changing the OFV (Results, Pharmacokinetics paragraph 2); the
    # corresponding etas are omitted here.
    etalfdepot ~ log(1 + 2.87^2)  # Table 2: IIV F_PO = 287% CV (49.1% RSE)
    etalka     ~ log(1 + 2.32^2)  # Table 2: IIV ka_PO = 232% CV (49.7% RSE); note eta-shrinkage 35.8% in source analysis
    etalcl     ~ log(1 + 0.818^2) # Table 2: IIV CL = 81.8% CV (69.7% RSE)
    etalq      ~ log(1 + 0.774^2) # Table 2: IIV Q = 77.4% CV (44.0% RSE); note eta-shrinkage 54.9% in source analysis
    etalvp     ~ log(1 + 4.28^2)  # Table 2: IIV V_P = 428% CV (46.5% RSE)

    # Residual error. The source paper modelled natural-log plasma concentrations
    # with additive residual on the log scale (Methods, Population
    # pharmacokinetic and pharmacodynamic analysis: "The residual unexplained
    # variability was assumed to be additive on a logarithmic scale, essentially
    # equivalent to an exponential error on an arithmetic scale."). In nlmixr2
    # conventions this maps to proportional residual error in linear space;
    # propSd is the SD on the log scale, which equals the proportional CV to
    # first order.
    propSd <- sqrt(0.376)
    label("Proportional residual SD (SD on log scale)")  # Table 2: sigma = 0.376 (variance, 95% CI 0.316-0.436; 7.8% RSE); SD = sqrt(0.376) = 0.613
  })

  model({
    # Individual PK parameters.
    f_po   <- exp(lfdepot + etalfdepot)
    ka     <- exp(lka     + etalka)
    cl     <- exp(lcl     + etalcl)
    vc     <- exp(lvc)
    q      <- exp(lq      + etalq)
    vp     <- exp(lvp     + etalvp)
    dur_im <- exp(ldur_im)

    # Micro-constants.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ODE system. Oral syrup doses target depot for first-order absorption;
    # intramuscular doses target central with rate = -2 to invoke the modeled
    # dur(central) = DUR_IM (zero-order direct delivery to central). F_IM = 1
    # is the rxode2 default for f(central) so it is not specified explicitly.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-               k12 * central - k21 * peripheral1

    f(depot)     <- f_po
    dur(central) <- dur_im

    # Plasma concentration: dose in mg, volume in L -> mg/L.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
