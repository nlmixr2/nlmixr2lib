Nie_2025_sufentanil <- function() {
  description <- paste(
    "Two-compartment population PK model of epidural sufentanil in labouring women (Nie 2025;",
    "N = 41 primiparous women receiving sufentanil 0.3 ug/mL with 0.1% ropivacaine by",
    "patient-controlled epidural analgesia at six Chinese centres). Compartment 1 is the maternal",
    "central compartment and compartment 2 is the umbilical cord (placental circulation)",
    "compartment; the two exchange through a slow inter-compartmental clearance CL2. The epidural",
    "compartment could not be estimated from the sparse data and was assumed to have merged with",
    "the central compartment, so doses enter `central` directly and CL / V1 are apparent values",
    "(bioavailability assumed 1). Cord equilibration is far slower than maternal elimination",
    "(K21 = CL2/V2 gives a cord half-life of ~9.7 h against a ~2.0 h maternal elimination",
    "half-life), which reproduces the paper's central finding that umbilical cord concentrations",
    "decline only slowly once the epidural infusion is stopped. Inter-individual variability was",
    "estimated on CL and V1 only; IIV on CL2 and V2 did not improve the fit and was excluded.",
    "Residual error is proportional with a separate magnitude for maternal and umbilical cord",
    "observations.",
    sep = " "
  )
  reference <- paste(
    "Nie Y, Sun X, Cao R, Tang S, Zhou Q, Zhou M, Chen Z, Huang S.",
    "Population pharmacokinetic of epidural sufentanil in labouring women:",
    "a multicentric, prospective, observational study.",
    "Drug Des Devel Ther. 2025;19:971-980.",
    "doi:10.2147/DDDT.S500189.",
    sep = " "
  )
  vignette <- "Nie_2025_sufentanil"

  # `cord` is the umbilical cord / placental-circulation compartment named by
  # Nie 2025 Figure 1 ("V2, funicle umbilical cord compartment"). It is not a
  # generic `peripheral1`: it is an observed physiologic matrix with its own
  # residual-error magnitude, and its volume (0.187 L) is an umbilical cord
  # volume rather than a maternal tissue distribution volume. Declared
  # paper-specific following the placental-transfer precedents already in the
  # library -- Hirt_2007_nelfinavir.R (`cord_n` / `cord_m8`),
  # Fauchet_2015_lopinavir_placental.R (`fetal`) and
  # Ngamprasertwong_2016_propofol_sheep.R (`fetus`).
  paper_specific_compartments <- c("cord")

  # Nie 2025 reports CL and CL2 in L/h and V1 and V2 in L, and plots observed
  # concentrations in pg/mL (Figure 2 y-axis; LLOQ 1 pg/mL). Amount in ng
  # divided by volume in L gives ng/L, which is identical to pg/mL, so the
  # published parameter values are used unchanged with no conversion constant.
  # A 4.5 ug epidural loading dose is therefore entered as 4500 ng.
  units <- list(time = "h", dosing = "ng", concentration = "pg/mL")

  compartmentData <- list(
    central = list(analyte = "sufentanil", units = "ng", specimen = "plasma", verified = TRUE),
    cord    = list(analyte = "sufentanil", units = "ng", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 41L,
    n_studies      = 1L,
    age_range      = "20-40 years (inclusion criterion)",
    age_mean       = "27.2 +/- 2.9 years",
    weight_mean    = "70.2 +/- 7.5 kg",
    sex_female_pct = 100,
    disease_state  = paste(
      "Primiparous women with a singleton pregnancy requesting epidural labour analgesia,",
      "ASA physical status II-III, BMI 18.5-30 kg/m^2, cervical dilation 2 <= phi < 6 cm and a",
      "visual analogue scale score >= 50 mm before analgesia. Women with cardiopulmonary disease,",
      "insulin-treated diabetes, neuromuscular disease, alcohol or drug abuse, allergy or",
      "contraindications to epidural anaesthesia were excluded, as were women who did not require",
      "patient-controlled analgesia, who failed to respond to the loading dose, or whose epidural",
      "catheter became displaced.",
      sep = " "
    ),
    dose_range     = paste(
      "Epidural sufentanil 0.3 ug/mL mixed with 0.1% ropivacaine. Loading: 5 mL of the mixture",
      "with 1/400,000 epinephrine, then a 10 mL bolus 5 min later (15 mL, 4.5 ug total).",
      "Patient-controlled epidural analgesia from 30 min: continuous infusion 6 mL/h (1.8 ug/h),",
      "6 mL demand bolus, 15 min lockout, maximum 30 mL/h; maintained until delivery.",
      sep = " "
    ),
    regions        = paste(
      "China; six centres (Obstetrics and Gynecology Hospital of Fudan University, Shanghai;",
      "Chengdu Women's and Children's Central Hospital; First Affiliated Hospital of Guizhou",
      "University of Traditional Chinese Medicine, Guiyang; Jiangxi Maternal and Child Health",
      "Hospital, Nanchang; Fujian Provincial Maternity and Children's Hospital, Fuzhou;",
      "Affiliated Hospital of Qingdao Binhai University).",
      sep = " "
    ),
    n_observations = paste(
      "79 maternal central-compartment and 37 umbilical cord concentrations from 41 women",
      "(Results, 'Patient Characteristics'). Maternal venous samples were scheduled at 0.5 h and",
      "1 h after the epidural loading dose, at delivery and 2 h after delivery; a single umbilical",
      "venous sample was drawn immediately after cord clamping. Not all four maternal samples were",
      "mandatory, so the design is sparse and unbalanced. The Figure 4 caption quotes n = 40",
      "central and n = 35 cord observations at delivery, a subset of the totals above; five",
      "subjects had no umbilical cord measurement and were predicted from population mean",
      "parameters.",
      sep = " "
    ),
    notes          = paste(
      "Secondary PK analysis of a previously reported prospective trial comparing epidural",
      "nalbuphine-ropivacaine with sufentanil-ropivacaine (ChiCTR1800018810); enrolment November",
      "2018 to February 2019, concentration measurement February 2020. 90 parturients were",
      "recruited and 41 contributed to the final population PK analysis. Sufentanil was assayed by",
      "LC-MS/MS (Sciex Triple Quad 6500+) with sufentanil-d5 internal standard; LLOQ 1 pg/mL,",
      "inter- and intra-assay CV both < 15%. Estimation by FOCE with interaction in NONMEM 7.4.0;",
      "internal validation by 5000-run non-parametric bootstrap (PsN) and prediction-corrected VPC",
      "from 1000 simulated datasets. No covariate model was reported: Table 1 contains structural,",
      "inter-individual-variability and residual rows only, and the paper describes no covariate",
      "screening step. All women received concomitant epidural ropivacaine, which the authors note",
      "as a possible source of unmodelled pharmacokinetic interaction.",
      sep = " "
    )
  )

  ini({
    # Final-model population estimates, Nie 2025 Table 1 ("Population Parameter
    # Estimates for the Final Sufentanil Pharmacokinetic Model"). The
    # parenthesised percentages in that table are relative standard errors of
    # the estimates, not variability. NONMEM 7.4.0, FOCE+I.
    #
    # The structural model is defined in Methods, "Sufentanil Structural
    # Pharmacokinetic Model": A1 and A2 are the amounts in the central and
    # umbilical cord compartments, V1 and V2 their volumes, K12 = CL2/V1,
    # K21 = CL2/V2 and K10 = CL/V1.
    lcl    <- log(176)    ; label("Maternal central clearance (CL, L/h)")                                # Table 1: CL = 176 (RSE 9.25 percent)
    lvc    <- log(519)    ; label("Maternal central volume of distribution (V1, L)")                     # Table 1: V1 = 519 (RSE 6.45 percent)
    lqcord <- log(0.0134) ; label("Central-to-umbilical-cord inter-compartmental clearance (CL2, L/h)")  # Table 1: CL2 = 0.0134 (RSE 43.2 percent)
    lvcord <- log(0.187)  ; label("Umbilical cord volume of distribution (V2, L)")                       # Table 1: V2 = 0.187 (RSE 39.9 percent)

    # Inter-individual variability. Methods, "Population PK Model Development":
    # IIV was assumed log-normal and modelled with the exponential relationship
    # P_i = P_TV * exp(eta_i), where eta_i has mean zero and variance omega^2.
    # Table 1 reports the IIV magnitudes as %CV, so the internal variance is
    # omega^2 = log(CV^2 + 1).
    #
    # Results, "Population Pharmacokinetic Modeling": "The addition of IIV to
    # CL2 and V2 did not contribute to better model fitting and was therefore
    # not included in the final estimates" -- those etas are omitted entirely.
    etalcl ~ 0.255046  # Table 1: IIV of CL = 53.9 %CV (RSE 12 percent);   omega^2 = log(0.539^2 + 1) = 0.255046
    etalvc ~ 0.046429  # Table 1: IIV of V1 = 21.8 %CV (RSE 26.9 percent); omega^2 = log(0.218^2 + 1) = 0.046429

    # Residual error. Methods describe a combined proportional-plus-additive
    # model on linear-scale concentration, but state that "During the model
    # process, eps_add was removed due to its lesser significance". Table 1
    # therefore reports only the two proportional components, one per observed
    # matrix, expressed as percentages that equal the proportional-residual SD.
    propSd       <- 0.206 ; label("Proportional residual SD for maternal Cc (fraction)")               # Table 1: eps_prop1 = 20.6 percent (RSE 13.4 percent)
    propSd_Ccord <- 0.355 ; label("Proportional residual SD for umbilical cord Ccord (fraction)")      # Table 1: eps_prop2 = 35.5 percent (RSE 12 percent)
  })
  model({
    # Individual parameters. No covariate effects are applied: Nie 2025 reports
    # no covariate model. IIV is carried on CL and V1 only.
    cl    <- exp(lcl + etalcl)
    vc    <- exp(lvc + etalvc)
    qcord <- exp(lqcord)
    vcord <- exp(lvcord)

    # Micro-rate constants exactly as defined in Methods, "Sufentanil
    # Structural Pharmacokinetic Model".
    kel <- cl    / vc      # K10 = CL / V1
    k12 <- qcord / vc      # K12 = CL2 / V1
    k21 <- qcord / vcord   # K21 = CL2 / V2

    # Epidural doses are administered into `central`. The paper's Figure 1
    # draws a separate epidural compartment (Ve) upstream of the central
    # compartment, but Methods state that "due to limited data, the parameters
    # of the epidural compartment could not be estimated accurately.
    # Therefore, for the purpose of this analysis, it was assumed that the
    # epidural compartment merged with the central compartment." The fitted
    # model has no absorption step.
    d/dt(central) <- -kel * central - k12 * central + k21 * cord
    d/dt(cord)    <-  k12 * central - k21 * cord

    # Observations. Maternal venous plasma (Cc) is the canonical central
    # output; umbilical venous plasma (Ccord) is the paper-named output on the
    # umbilical cord compartment.
    Cc    <- central / vc
    Ccord <- cord    / vcord

    Cc    ~ prop(propSd)
    Ccord ~ prop(propSd_Ccord)
  })
}
