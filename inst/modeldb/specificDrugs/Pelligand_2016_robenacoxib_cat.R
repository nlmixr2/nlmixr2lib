Pelligand_2016_robenacoxib_cat <- function() {
  description <- "Preclinical/clinical veterinary (cat). Two-compartment population PK model for robenacoxib in cats, parameterised per kg body weight, pooling intravenous and subcutaneous dosing from eight studies. Subcutaneous absorption is a parallel mixed-order input: a fraction F0 of the bioavailable dose enters the central compartment by a zero-order process of duration Tk0 = 1.78 h and the remainder arrives first-order through a depot with ka = 0.68 1/h, giving flip-flop kinetics because ka is slower than the disposition terminal rate constant. General anaesthesia doubles the central volume of distribution; no other covariate was retained (Pelligand 2016)"
  reference <- paste(
    "Pelligand L, Soubret A, King JN, Elliott J, Mochel JP. Modeling of large",
    "pharmacokinetic data using nonlinear mixed-effects: a paradigm shift in",
    "veterinary pharmacology. A case study with robenacoxib in cats. CPT",
    "Pharmacometrics Syst Pharmacol. 2016;5(12):625-635. doi:10.1002/psp4.12141.",
    "Structural model transcribed from the MLXTRAN source listing deposited as",
    "Supplementary Data (PSP4-5-625-s008.txt); residual-error magnitudes and the",
    "anaesthesia effect on V1 digitised from the SAEM convergence traces of",
    "Supplementary Figure S4 (PSP4-5-625-s004.pdf).",
    sep = " "
  )
  vignette <- "Pelligand_2016_robenacoxib_cat"
  # Every structural parameter is published per kg body weight (CL in L/h/kg,
  # volumes in L/kg), so doses are given in mg/kg and compartment amounts are
  # carried in mg/kg. central/vc is then mg/L; the x1000 in Cc reports ng/mL to
  # match the bioanalytical method (Methods 'Analytical phase': HPLC-UV
  # 500-20,000 ng/mL, LC-MS 3-100 ng/mL) and the additive residual SD.
  units <- list(time = "h", dosing = "mg/kg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(
      analyte = "robenacoxib",
      units = "mg/kg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(analyte = "robenacoxib", units = "mg/kg", specimen = "whole blood", verified = TRUE),
    peripheral1 = list(analyte = "robenacoxib", units = "mg/kg", specimen = "whole blood", verified = TRUE)
  )

  covariateData <- list(
    ANESTH_GA = list(
      description = "Indicator for general anaesthesia during the pharmacokinetic observation window",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (conscious, unanaesthetised laboratory cat)",
      notes = paste(
        "Per-administration indicator. ANESTH_GA = 1 identifies the 36 clinical cats of the",
        "perioperative study, admitted for elective ovariohysterectomy and given robenacoxib",
        "s.c. under general anaesthesia (intramuscular buprenorphine + acepromazine",
        "premedication, propofol induction, isoflurane maintenance after endotracheal",
        "intubation; Methods 'Animal phase for sparse sampling dataset'). ANESTH_GA = 0",
        "identifies the 47 conscious laboratory cats of the seven CRA studies. The only",
        "covariate retained by the BIC backward elimination, acting on the central volume of",
        "distribution: Results 'Effects of demographics and anaesthesia on robenacoxib",
        "exposure' reports V1 = 0.33 L/kg in the clinical cats against 0.16 L/kg in the",
        "densely sampled laboratory cats, and Table 2 carries the two values as separate",
        "V1 rows labelled 'ANEST 5 0' and 'ANEST 5 1' (the '5' is a mis-rendered '=' in the",
        "converted text). The paper attributes the larger volume to volume expansion from",
        "intra-operative fluid therapy and to the vasodilatory anaesthetic agents",
        "(Discussion). The effect is a COHORT contrast, not a within-surgery time window:",
        "the ANESTH_GA = 1 blood samples were drawn at extubation and 2 h thereafter, i.e.",
        "post-operatively, so this column is NOT the canonical INTRAOP (whose reference",
        "category is explicitly 'pre- or post-operative'). Anaesthesia and sparse sampling",
        "are perfectly confounded in this design, which the paper flags as a limitation",
        "(Discussion, 'The outcome of the covariate analysis should be interpreted with",
        "caution').",
        sep = " "
      ),
      source_name = "ANEST"
    ),
    ROUTE_IV = list(
      description = "Indicator for intravenous administration",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (subcutaneous administration between the scapulae)",
      notes = paste(
        "Per-dose-record indicator (1 = intravenous bolus, 0 = subcutaneous). 23 cats were",
        "dosed i.v. and 74 s.c. over 97 administrations; 14 cats received both routes, so the",
        "column varies within subject and must be set per dose record. The MLXTRAN listing",
        "(Supplementary Data PSP4-5-625-s008.txt) declares both inputs on compartment 1:",
        "'iv(cmt=1,adm=2)' with no bioavailability term alongside the two s.c. inputs",
        "'oral(cmt=1,adm=1,Tk0=Tk0,p=Ftot*F0)' and 'oral(cmt=1,adm=1,ka=ka,p=Ftot*(1-F0))'.",
        "rxode2 applies one bioavailability per compartment rather than per administration",
        "type, so ROUTE_IV selects between them: f(central) collapses to 1 for the i.v. bolus",
        "and to Ftot*F0 for the zero-order half of the s.c. dose. An i.v. record must be a",
        "plain bolus into cmt = 'central' (no rate), which also makes rxode2 ignore the",
        "modelled dur(central). A s.c. administration needs TWO records at the same time:",
        "a bolus into cmt = 'depot' and a record into cmt = 'central' carrying rate = -2 so",
        "the modelled zero-order duration d1 is applied. ROUTE_IV = 0 on both s.c. records.",
        sep = " "
      ),
      source_name = "route of administration (IV)"
    )
  )

  # Weight, age and sex were screened on CL, V1 and Ftot and none was retained
  # by the BIC backward elimination (Results 'Effects of demographics and
  # anaesthesia on robenacoxib exposure'; Supplementary Figures S1, S5, S6). The
  # paper reports these screens only graphically (r^2 annotations on the
  # posterior-distribution scatterplots) and publishes no point estimate, so
  # they are documented rather than encoded.
  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Screened as t_WEIGHT0 (median-normalised and log-transformed, Methods 'Inclusion of",
        "covariate relationships') on CL, V1 and Ftot; not retained. Supplementary Figure S5",
        "reports r^2 = 0.0381 for log(CL) against t_WEIGHT0 and Supplementary Figure S1",
        "r^2 = 0.00911. Discussion: 'there is no need for dosing adjustment based on",
        "population demographics (age and bodyweight)', explicitly contrasted with the dog,",
        "where body weight did affect apparent robenacoxib CL and volume. Range 1.77-5.7 kg.",
        sep = " "
      ),
      source_name = "t_WEIGHT0"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Screened as t_AGE (median-normalised and log-transformed) on CL, V1 and Ftot; not",
        "retained. Supplementary Figure S5 reports r^2 = 0.000798 for log(CL) against t_AGE",
        "and Supplementary Figure S1 r^2 = 0.00356. Range 0.34-6.1 years.",
        sep = " "
      ),
      source_name = "t_AGE"
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = paste(
        "Screened as GENDER on CL, V1 and Ftot; not retained. Results: 'A modest but not",
        "significant difference was observed between the median CL of male and female cats'",
        "(Supplementary Figure S6). Not interpretable independently of anaesthesia in this",
        "design because all 36 perioperative cats were female (Discussion limitation).",
        sep = " "
      ),
      source_name = "GENDER"
    )
  )

  population <- list(
    species = "cat (domestic shorthair; Felis catus)",
    n_subjects = 83L,
    n_studies = 8L,
    n_administrations = 97L,
    n_observations = 652L,
    age_range = "0.34-6.1 years",
    weight_range = "1.77-5.7 kg",
    sex_female_pct = 69.9,
    disease_state = paste(
      "47 healthy conscious laboratory cats (seven CRA drug-development studies) plus 36",
      "clinical cats admitted for elective ovariohysterectomy and sampled perioperatively",
      "under general anaesthesia",
      sep = " "
    ),
    dose_range = "1.6-2.3 mg/kg single dose, intravenous or subcutaneous",
    regions = "France, Switzerland, United Kingdom",
    notes = paste(
      "Demographics from Table 1 (per-study listing) and Results 'Study demographics'.",
      "Laboratory cats: 22 female / 27 male, median body weight 3.9 kg (IQR 3.45-4.14,",
      "range 2.5-5.7), median age 1.32 years (IQR 1.0-1.55, range 0.97-6.1), 9-12 blood",
      "samples per cat. Clinical cats: 36 female, median body weight 2.73 kg (IQR",
      "2.41-2.99, range 1.77-4.0), median age 0.76 years (range 0.34-4.31), 1-2 samples per",
      "cat drawn 0.8-8.2 h after dosing. 23 cats received robenacoxib i.v. and 74 s.c., with",
      "14 receiving both routes at least one week apart. The female percentage is computed",
      "as (22 + 36) / 83; the paper's own laboratory-cat counts sum to 49 rather than 47.",
      "55 of 652 measured concentrations were below the limit of quantification and were",
      "handled by the Monolix equivalent of the NONMEM M3 method (Beal 2001).",
      sep = " "
    )
  )

  ini({
    # ---- Disposition (Table 2, 'IV model (initial fitting)') ----------------
    # The typical values of CL, V1, Q and V2 were estimated from the i.v. data
    # alone and then HELD FIXED for the sequential fit of the s.c. + i.v. data
    # (Table 2 footnote and Methods 'Data analysis and model evaluation'), which
    # is why Table 2 prints no relative standard error for them.
    lcl <- fixed(log(0.502))
    label("Clearance (L/h/kg)")
    # Results 'PKs': 'The systemic total body CL was estimated to be moderate
    # (0.502 L/kg/h)'; Table 2 prints the same value rounded to 0.50.

    lvc <- fixed(log(0.166))
    label("Central volume of distribution in a conscious cat (L/kg)")
    # Table 2 prints V1 (ANEST = 0) = 0.16 L/kg to two significant figures. Two
    # independent statements in the paper pin the unrounded value at 0.166:
    # (i) Results 'PKs' gives the population steady-state volume as 0.213 L/kg,
    # and Vss = V1 + V2, so V1 = 0.213 - 0.047 = 0.166; (ii) the beta_ANEST,V1
    # trace of Supplementary Figure S4 converges to 0.684, and V1 = 0.33 /
    # exp(0.684) = 0.1664. Using the displayed 0.16 instead would reproduce
    # neither the published Vss (it gives 0.207) nor the published V1 for an
    # anaesthetised cat.

    lvp <- fixed(log(0.047))
    label("Peripheral volume of distribution (L/kg)")
    # Table 2, row '-Peripheral compartment volume of distribution', V2 = 0.047
    # L/kg. Methods 'Data analysis and model evaluation' miscalls V2 'volume of
    # the central compartment'; Table 2 and the MLXTRAN 'peripheral(k12, k21,
    # amount=Ap)' block both make it the peripheral volume.

    lq <- fixed(log(0.065))
    label("Intercompartmental clearance (L/h/kg)")
    # Table 2, row '-Inter-compartmental clearance', Q = 0.065 L/h/kg.

    # ---- Subcutaneous absorption (Table 2, 'SC model (sequential fitting)') --
    lka <- log(0.68)
    label("First-order absorption rate constant (1/h)")
    # Table 2, row '-Absorption rate', ka = 0.68 1/h (RSE 5%).

    ld1 <- log(1.78)
    label("Duration of the zero-order absorption input (h)")
    # Table 2, row '-Absorption duration (0-order)', Tk0 = 1.78 h (RSE 5%).

    logitfdepot <- log(0.78 / (1 - 0.78))
    label("Logit of the total subcutaneous bioavailability Ftot (fraction)")
    # Table 2, row '-Bioavailability', Ftot = 0.78 (RSE 3%). Methods 'PK model
    # development': 'For the bioavailability model parameter, a logit transform
    # function was used to insure that the estimated value was bounded between
    # zero and one'. logit(0.78) = 1.2657.

    logitffo <- log(0.50 / (1 - 0.50))
    label("Logit of the fraction of the bioavailable dose absorbed first-order (fraction)")
    # Table 2 reports the complement, row '-Fraction absorbed through 0-order',
    # F0 = 0.50 (RSE 10%), so the first-order fraction is 1 - F0 = 0.50 and
    # logit(0.50) = 0. F0 is a PARALLEL SPLIT OF THE DOSE, not the sequential
    # zero-order time fraction 'fzo': the MLXTRAN listing applies it as the
    # bioavailability of one of two simultaneous inputs
    # ('oral(cmt=1,adm=1,Tk0=Tk0,p=Ftot*F0)' beside
    # 'oral(cmt=1,adm=1,ka=ka,p=Ftot*(1-F0))'), and the arithmetic test for the
    # sequential reading fails - fzo/(d1*(1-fzo)) = 0.50/(1.78*0.50) = 0.562,
    # which is not the published ka of 0.68.

    # ---- Covariate effect ---------------------------------------------------
    e_anesth_ga_vc <- log(0.33 / 0.166)
    label("Effect of general anaesthesia on the log central volume of distribution")
    # Table 2 carries two V1 rows, 0.16 L/kg at ANEST = 0 and 0.33 L/kg at
    # ANEST = 1 (RSE 14%); Results 'Effects of demographics and anaesthesia on
    # robenacoxib exposure' states the same pair. Monolix applies a categorical
    # covariate additively on the log scale of a log-normal parameter, so the
    # coefficient is log(0.33 / 0.166) = 0.6869. The beta_ANEST,V1 convergence
    # trace of Supplementary Figure S4 reads 0.684, agreeing to 0.003.

    # ---- Inter-individual variability --------------------------------------
    # Table 2's 'IIV (%)' column is omega x 100, the SD of the random effect on
    # the TRANSFORMED (log, or logit for Ftot and F0) scale - it is NOT a
    # back-transformed coefficient of variation. Settled against the supplement:
    # the eight omega convergence traces of Supplementary Figure S4 read 0.571,
    # 0.207, 0.355, 0.162, 0.415, 0.081, 0.0125 and 0.93 for Ftot, ka, Tk0, CL,
    # V1, Q, V2 and F0, each reproducing its Table 2 percentage to two decimal
    # places. Entries below are variances, omega^2.
    etalcl ~ 0.0256 # Table 2 IIV(CL) = 16%; 0.16^2
    etalvc ~ 0.1681 # Table 2 IIV(V1) = 41%; 0.41^2
    etalvp ~ 0.0001 # Table 2 IIV(V2) = 1%; 0.01^2
    etalq ~ 0.0064 # Table 2 IIV(Q) = 8%; 0.08^2
    etalogitfdepot ~ 0.3249 # Table 2 IIV(Ftot) = 57%; 0.57^2. Results says 'an interindividual variability of 3%', which is the RSE column, not the IIV column
    etalogitffo ~ 0.8649 # Table 2 IIV(F0) = 93%; 0.93^2. logit(1-x) = -logit(x), so the variance of the random effect is the same whether it is carried on F0 or on its complement

    # Results 'PKs': 'Adding a correlation between ka and Tk0 (corr 5 -0.24)
    # yielded the best objective function and BIC value and was, therefore,
    # included in the model structure' (the '5' is a mis-rendered '='). This is
    # the only correlation retained; Figure 5 shows the remaining random-effect
    # correlations to be flat when multiple posterior samples are used instead
    # of the EBEs. Covariance = -0.24 * 0.21 * 0.35 = -0.01764.
    etalka + etald1 ~ c(
      0.0441,
      -0.01764, 0.1225
    )

    # ---- Residual error -----------------------------------------------------
    # Methods 'PK model development' specifies a combined error on the linear
    # scale, y = F + a*eps_a + b*F*eps_b, with a and b as standard deviations.
    # Neither magnitude is printed in the text or in Table 2. Both are read off
    # the final-estimate markers of the SAEM convergence traces in
    # Supplementary Figure S4 (panels 'a' and 'b'). The digitisation is
    # validated against the eight parameters on the same page whose values ARE
    # tabulated - Ftot and the seven omegas above all reproduce Table 2 to two
    # decimal places - so the residual magnitudes carry the same accuracy.
    addSd <- 2.79
    label("Additive residual SD (ng/mL)")
    # Supplementary Figure S4 panel 'a' converges to 2.79 ng/mL. Consistent with
    # the 3 ng/mL lower limit of quantification of the LC-MS method (Methods
    # 'Analytical phase') and with the M3-equivalent BLQ likelihood term.

    propSd <- 0.21
    label("Proportional residual SD (fraction)")
    # Supplementary Figure S4 panel 'b' converges to 0.21.
  })

  model({
    # ---- Individual parameters ---------------------------------------------
    cl <- exp(lcl + etalcl)
    # Monolix enters a categorical covariate additively on the log scale, so the
    # anaesthetised cats carry exp(0.6869) = 1.99 times the conscious central
    # volume: 0.166 -> 0.330 L/kg.
    vc <- exp(lvc + e_anesth_ga_vc * ANESTH_GA + etalvc)
    vp <- exp(lvp + etalvp)
    q <- exp(lq + etalq)
    ka <- exp(lka + etalka)
    d1 <- exp(ld1 + etald1)

    # Total subcutaneous bioavailability and its split between the two parallel
    # inputs. ffo is the first-order share; 1 - ffo is the paper's F0, the share
    # entering the central compartment by the zero-order process.
    fdepot <- expit(logitfdepot + etalogitfdepot)
    ffo <- expit(logitffo + etalogitffo)

    # ---- Micro-constants (MLXTRAN 'PK:' block of PSP4-5-625-s008.txt) ------
    # k = CL/V1, k12 = Q/V1, k21 = Q/V2.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ---- ODE system ---------------------------------------------------------
    # Two-compartment mammillary disposition with first-order elimination from
    # central (Figure 1). `depot` carries only the first-order half of the
    # subcutaneous dose; the zero-order half is delivered straight into
    # `central` as a modelled-duration input, exactly as the MLXTRAN listing
    # declares both s.c. inputs on compartment 1.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot + k21 * peripheral1 - (kel + k12) * central
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # ---- Dose partitioning --------------------------------------------------
    # A subcutaneous administration is entered as two records at the same time:
    # a bolus into `depot` and a rate = -2 record into `central`. An
    # intravenous administration is a single plain bolus into `central`; with no
    # rate on that record rxode2 ignores dur(central), so the i.v. dose stays a
    # bolus and reads f(central) = 1.
    f(depot) <- fdepot * ffo
    f(central) <- ROUTE_IV + (1 - ROUTE_IV) * fdepot * (1 - ffo)
    dur(central) <- d1

    # ---- Observation --------------------------------------------------------
    # central in mg/kg and vc in L/kg give mg/L; x1000 reports ng/mL, the unit
    # of the bioanalytical method and of the additive residual SD.
    Cc <- central / vc * 1000
    Cc ~ add(addSd) + prop(propSd)
  })
}
