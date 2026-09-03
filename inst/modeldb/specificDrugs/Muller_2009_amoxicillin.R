Muller_2009_amoxicillin <- function() {
  description <- paste(
    "Five-compartment population PK model for intravenous amoxicillin in",
    "labouring women, the venous umbilical cord and the neonate, fitted",
    "simultaneously to maternal, arterial and venous umbilical cord and",
    "neonatal serum concentrations. Three maternal compartments (central",
    "V1 plus two peripheral compartments constrained to a common volume",
    "V2 = V3) exchange with a venous umbilical cord compartment V4 via",
    "the transplacental rate constants k14 / k41; a neonatal compartment",
    "V5 hangs off the cord compartment via k45 / k54 and eliminates drug",
    "from the system via k50. Central volume increases 4.2% per week of",
    "gestational age (Muller 2009)."
  )
  reference <- paste(
    "Muller AE, Oostvogel PM, DeJongh J, Mouton JW, Steegers EAP,",
    "Dorr PJ, Danhof M, Voskuyl RA. Pharmacokinetics of amoxicillin in",
    "maternal, umbilical cord, and neonatal sera.",
    "Antimicrob Agents Chemother. 2009;53(4):1574-1580.",
    "doi:10.1128/AAC.00119-08."
  )
  vignette <- "Muller_2009_amoxicillin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # The venous umbilical cord and neonatal states are paper-mechanistic
  # matrices of the maternal-fetal dyad; they are not members of any
  # canonical compartment family. Named after the matrices they hold
  # (Muller 2009 Fig. 2: "Umbilical Cord" V4, "Neonate" V5) rather than
  # after the ADVAN5 slot indices, following the precedent set by
  # Fauchet_2015_lopinavir_placental.R ("fetal", "amniotic").
  paper_specific_compartments <- c("cord_venous", "neonate")

  compartmentData <- list(
    central = list(
      analyte = "amoxicillin", units = "mg", specimen = "serum",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "amoxicillin", units = "mg", specimen = "serum",
      verified = TRUE
    ),
    peripheral2 = list(
      analyte = "amoxicillin", units = "mg", specimen = "serum",
      verified = TRUE
    ),
    # V4 in Muller 2009 Fig. 2. Fitted to the 28 venous umbilical cord
    # serum samples (Results paragraph 1; Fig. 3 middle panel).
    cord_venous = list(
      analyte = "amoxicillin", units = "mg", specimen = "serum",
      verified = TRUE
    ),
    # V5 in Muller 2009 Fig. 2. Fitted jointly to the 25 ARTERIAL
    # umbilical cord serum samples and the 14 neonatal heel-puncture
    # samples, which the authors deliberately pooled into one
    # compartment ("we adopted a simplified approach in which the
    # arterial umbilical cord and neonatal serum concentrations were
    # placed in the same compartment", Discussion paragraph 3;
    # Fig. 3 bottom panel legend).
    neonate = list(
      analyte = "amoxicillin", units = "mg", specimen = "serum",
      verified = TRUE
    )
  )

  covariateData <- list(
    GA = list(
      description        = paste(
        "Gestational age. The only covariate retained in the final",
        "model; it acts on the maternal central volume V1 through the",
        "linear-fractional relationship of Muller 2009 Results",
        "paragraph 4."
      ),
      units              = "weeks",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Centred on the study-population median of 36.8 weeks, stated",
        "verbatim in the covariate equation ('GA is the gestational age",
        "centered by its median value in the study population (36.8",
        "weeks)'). Consistent with the cohort mean amenorrhoea of 36 6/7",
        "weeks (= 36.857) in Table 1. Deliveries occurred at 30.0-42.4",
        "weeks (Results paragraph 1). Body mass index was also",
        "significant on V (dOFV = -8.5) but was dropped in favour of",
        "gestational age (dOFV = -10.5) because the two are correlated",
        "and adding both gave no further improvement (dOFV = -3.3)."
      ),
      source_name        = "GA"
    )
  )

  covariatesDataExcluded <- list(
    BMI = list(
      description = "Maternal body mass index.",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = paste(
        "Table 1: 29.4 +/- 5.3 kg/m^2 (n = 43). Screened on V and found",
        "significant on its own (dOFV = -8.5) but not retained: it is",
        "correlated with gestational age, which gave the larger drop",
        "(dOFV = -10.5), and the two together added nothing",
        "(dOFV = -3.3). Results paragraph 4. No coefficient is reported",
        "for the BMI effect, so it could not be implemented even as an",
        "alternative parameterisation."
      )
    ),
    WT = list(
      description = "Maternal body weight.",
      units       = "kg",
      type        = "continuous",
      notes       = paste(
        "Table 1: 79.4 +/- 13.9 kg (n = 43). Screened; maximum dOFV",
        "-0.5 (Results paragraph 4), not retained."
      )
    ),
    SBP = list(
      description = "Maternal blood pressure.",
      units       = "mmHg",
      type        = "continuous",
      notes       = paste(
        "Screened; maximum dOFV -0.3 (Results paragraph 4), not",
        "retained. The paper reports 'blood pressure' without stating",
        "systolic or diastolic; recorded in the standard workup",
        "(Methods 'Patients')."
      )
    ),
    HR = list(
      description = "Maternal pulse rate.",
      units       = "beats/min",
      type        = "continuous",
      notes       = "Screened; maximum dOFV -0.2 (Results paragraph 4), not retained."
    ),
    TEMP = list(
      description = "Maternal oral temperature.",
      units       = "degC",
      type        = "continuous",
      notes       = paste(
        "Screened but not evaluable: 'All models in which the",
        "temperature was incorporated resulted in running errors'",
        "(Results paragraph 4)."
      )
    ),
    OEDEMA = list(
      description = "Semi-quantitative oedema score, 0 (none) to 3 (above the knee).",
      units       = "(ordinal 0-3)",
      type        = "count",
      notes       = paste(
        "Table 1 reports the cohort distribution 29 / 12 / 3 across",
        "scores 0 / 1 / 2 (no patient scored 3). Screened on the",
        "parameters carrying interindividual variability (Methods 'PK",
        "analysis'); not retained and no coefficient reported."
      )
    ),
    WT_BIRTH = list(
      description = "Birth weight of the neonate.",
      units       = "g",
      type        = "continuous",
      notes       = paste(
        "Table 1: 2887.4 +/- 627.9 g (n = 46); range 1340-4470 g",
        "(Results paragraph 1). Screened on the maternal and neonatal",
        "parameters, including k45; not retained."
      )
    ),
    TWIN = list(
      description = "Twin versus singleton pregnancy indicator.",
      units       = "(binary)",
      type        = "binary",
      notes       = "Table 1: 3 twin pregnancies out of 44. Screened, not retained."
    ),
    DELIVERY_CESAREAN = list(
      description = "Mode of delivery: emergency caesarean section versus vaginal delivery.",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Table 1 reports 23 vaginal / 1 vacuum extraction / 3 emergency",
        "caesarean among the patients from whom umbilical cord serum was",
        "obtained. Screened, not retained."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 44L,
    n_studies      = 1L,
    n_observations = 971L,
    age_range      = "mean 30.0 +/- 6.85 years (maternal age, n = 44)",
    weight_range   = "mean 79.4 +/- 13.9 kg (n = 43)",
    sex_female_pct = 100,
    race_ethnicity = "not reported in the source paper",
    ga_range       = "30.0-42.4 weeks at delivery; cohort mean amenorrhoea 36 6/7 weeks (SD 2.7)",
    disease_state  = paste(
      "Women requiring intravenous amoxicillin or amoxicillin-clavulanic",
      "acid shortly before or during labour, for prevention of neonatal",
      "group B streptococcal disease (proven or unknown Streptococcus",
      "agalactiae carriage with recognised risk factors) or for suspected",
      "intra-amniotic infection. Part of the maternal data set overlaps",
      "the earlier PPROM cohort of Muller 2008 (416 shared maternal",
      "samples); no umbilical cord or neonatal sample was used in that",
      "earlier analysis."
    ),
    dose_range     = paste(
      "Amoxicillin 2 g IV over 30 min followed after 4 h by 1 g IV over",
      "15 min; or co-amoxiclav (1 g amoxicillin with 200 mg clavulanic",
      "acid) IV over 15 min every 8 h. All infusions at 50 mg/mL."
    ),
    regions        = "Netherlands (single centre, Medical Centre Haaglanden, The Hague)",
    notes          = paste(
      "904 maternal serum samples (3-41 per patient), 25 arterial and 28",
      "venous umbilical cord samples (53 total, from 44 women; both",
      "arterial and venous from 23 women; 4 cord samples from one twin",
      "pregnancy), and 14 neonatal heel-puncture samples from 13",
      "neonates taken 14.2-199.8 min after birth. Time from the last",
      "maternal dose to birth ranged from 24.4 to 316.8 min. Amoxicillin",
      "and co-amoxiclav patients were pooled because clavulanic acid does",
      "not alter intravenous amoxicillin PK (Discussion, final",
      "paragraph). Assay LLOQ 0.2 mg/L, between-run CV < 4%."
    )
  )

  ini({
    # Structural disposition parameters. Source: Muller 2009 Table 2,
    # 'Structural model parameters' block, Mean column. Every 95% CI in
    # that table reproduces as Mean +/- 1.96 * SE, which was used as an
    # internal transcription check.
    lcl  <- log(19.7);  label("Maternal clearance CL (L/h)")                                  # Table 2 row CL = 19.7 (SE 0.99, CI 17.8-21.6); also stated in the Abstract
    lvc  <- log(6.40);  label("Maternal central volume V1 at GA = 36.8 weeks (L)")            # Table 2 row V1 = 6.4 (SE 0.61, CI 5.2-7.6); Abstract '6.40 +/- 0.61 liters'
    lvp  <- log(5.88);  label("Maternal first peripheral volume V2 (L)")                      # Table 2 row V2 = 5.88 (SE 0.83, CI 4.2-7.5); Abstract '5.88 +/- 0.83 liters'
    lvp2 <- log(5.88);  label("Maternal second peripheral volume V3 (L), constrained V3 = V2") # Table 2 row V3 = 5.88, identical Mean/SE/CI to V2; Results paragraph 3: 'for the final analysis, V3 was assumed to be equal to V2'
    lq   <- log(56.6);  label("Maternal intercompartmental clearance Q1 between V1 and V2 (L/h)") # Table 2 row Q1 = 56.6 (SE 9.5, CI 38.0-75.2)
    lq2  <- log(10.7);  label("Maternal intercompartmental clearance Q2 between V1 and V3 (L/h)") # Table 2 row Q2 = 10.7 (SE 2.2, CI 6.3-15.1)

    # Transplacental and neonatal exchange rate constants. The paper
    # parameterises this half of the model with first-order rate
    # constants rather than clearances -- Methods 'PK analysis':
    # 'Because only venous cord blood enters the umbilical cord after
    # passage through the placenta, the antibiotic exchange between the
    # compartments might not be equal. Therefore, we used k values in our
    # model to describe the umbilical cord and neonatal data.' The
    # subscripts are the ADVAN5 slot indices of Fig. 2, kept verbatim as
    # the source's own parameter names: 1 = maternal central, 4 = venous
    # umbilical cord, 5 = neonate, 0 = eliminated from the system.
    lk14 <- log(0.76); label("Maternal central to venous umbilical cord transfer rate constant k14 (1/h)") # Table 2 row k14 = 0.76 (SE 0.28, CI 0.21-1.3)
    lk41 <- log(1.4);  label("Venous umbilical cord to maternal central transfer rate constant k41 (1/h)") # Table 2 row k41 = 1.4  (SE 0.31, CI 0.83-2.1)
    lk45 <- log(5.1);  label("Venous umbilical cord to neonate transfer rate constant k45 (1/h)")          # Table 2 row k45 = 5.1  (SE 2,    CI 1.1-9.1)
    lk54 <- log(1.4);  label("Neonate to venous umbilical cord transfer rate constant k54 (1/h)")          # Table 2 row k54 = 1.4  (SE 0.31, CI 0.83-2.1); identical Mean/SE/CI to k41, see the vignette Errata
    lk50 <- log(0.16); label("Neonatal elimination rate constant k50 (1/h)")                               # Table 2 row k50 = 0.16 (SE 0.033, CI 0.098-0.23)

    # Scaling volumes for the two fetal-side compartments. ADVAN5
    # estimates only rate constants for these slots, so the authors
    # back-calculated the volumes needed to turn amounts into
    # concentrations (Results paragraph 4: 'The values of V for the
    # venous umbilical cord and the neonate were calculated and were
    # found to be 3.4 liters for the venous umbilical cord and 11.9
    # liters for the neonate'). Both are reproduced by the flow-balance
    # identity implied by Fig. 2, which is the check that fixes the
    # direction of every k:
    #   V4 = k14 * V1 / k41 = 0.76 * 6.4 / 1.4 = 3.47 L  (paper: 3.4 L)
    #   V5 = k45 * V4 / k54 = 5.1 * 3.4  / 1.4 = 12.4 L  (paper: 11.9 L)
    # Both land inside the interval spanned by the two-significant-figure
    # rounding of the k values, so the reported 3.4 / 11.9 are carried
    # verbatim rather than recomputed. They are reported without a
    # standard error and are not estimated parameters, hence fixed().
    lv4 <- fixed(log(3.4));  label("Venous umbilical cord scaling volume V4 (L)")  # Results paragraph 4
    lv5 <- fixed(log(11.9)); label("Neonatal scaling volume V5 (L)")               # Results paragraph 4

    # Gestational-age effect on the maternal central volume, entered as
    # the additive-fractional form printed in Results paragraph 4:
    #   V1 = theta2 * [1 + theta12 * (GA - 36.8)] * exp(eta2)
    # 'An increase in V1 of 4.2% per week was found and was incorporated
    # into the model.' theta12 is not tabulated; 0.042 comes from that
    # sentence. Not reported with a standard error.
    e_ga_vc <- 0.042; label("Fractional change in V1 per week of gestational age above 36.8 weeks (1/week)")  # Results paragraph 4

    # Interindividual variability. Table 2's 'Variance model parameters'
    # block reports variances on the log scale: an exponential (log-normal)
    # model was used (Methods 'PK analysis': 'Individual estimates for PK
    # parameters were assumed to follow a log-normal distribution.
    # Therefore, an exponential distribution model was used to account for
    # interindividual variability'). Equivalent log-normal CV:
    #   sqrt(exp(0.076) - 1) = 28.1% on CL
    #   sqrt(exp(0.038) - 1) = 19.7% on V1
    # Results paragraph 3 states that 'A correlation between the random
    # parameters for interindividual variability was found and was
    # accounted for in the model', but the off-diagonal covariance is not
    # reported anywhere on disk, so a diagonal OMEGA is used here. See the
    # vignette 'Assumptions and deviations' section.
    etalcl ~ 0.076  # Table 2 row 'Interindividual variability in CL' (SE 0.026, CI 0.026-0.13)
    etalvc ~ 0.038  # Table 2 row 'Interindividual variability in V1' (SE 0.013, CI 0.014-0.063)

    # Residual error. A single proportional model was used for all three
    # matrices (Methods 'PK analysis': 'The residual variability ... was
    # described by using a proportional error model'; Results paragraph 3:
    # 'The residual error was found to be proportional to the blood
    # concentrations'). Table 2 reports the variance 0.03, so the
    # proportional SD is sqrt(0.03) = 0.1732 (17.3%).
    propSd          <- sqrt(0.03); label("Proportional residual SD on maternal Cc (fraction)")            # Table 2 row 'Residual variability' = 0.03 (SE 0.004, CI 0.022-0.037)
    propSd_Ccord    <- sqrt(0.03); label("Proportional residual SD on venous umbilical cord Ccord (fraction)")  # same Table 2 row; one residual-error model shared by all matrices
    propSd_Cneonate <- sqrt(0.03); label("Proportional residual SD on neonatal Cneonate (fraction)")      # same Table 2 row; one residual-error model shared by all matrices
  })

  model({
    # Maternal central volume with the gestational-age effect. Written
    # exactly as the source equation, with the eta multiplying the
    # covariate-adjusted typical value:
    #   V1 = theta2 * [1 + theta12 * (GA - 36.8)] * exp(eta2)
    # 36.8 weeks is the study-population median, printed in the equation
    # itself (Muller 2009 Results paragraph 4).
    vc  <- exp(lvc + etalvc) * (1 + e_ga_vc * (GA - 36.8))

    cl  <- exp(lcl + etalcl)
    vp  <- exp(lvp)
    vp2 <- exp(lvp2)
    q   <- exp(lq)
    q2  <- exp(lq2)
    v4  <- exp(lv4)
    v5  <- exp(lv5)

    # Maternal micro-constants.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # Transplacental and neonatal rate constants, estimated directly.
    k14 <- exp(lk14)
    k41 <- exp(lk41)
    k45 <- exp(lk45)
    k54 <- exp(lk54)
    k50 <- exp(lk50)

    # Five-compartment mass-transfer system exactly as drawn in Muller
    # 2009 Fig. 2: the maternal peripheral compartments V2 and V3 and the
    # venous umbilical cord V4 all hang off the maternal central V1; the
    # neonate V5 hangs off V4 and is the only fetal-side elimination
    # route ('Compartment 5 is attached to compartment 4, and the
    # antibiotics in compartment 5 are transferred back to compartment 4
    # and eliminated from the system', Results paragraph 3). Drug enters
    # V1 as an intravenous infusion set on the dose record.
    d/dt(central)     <- -(kel + k12 + k13 + k14) * central +
                          k21 * peripheral1 +
                          k31 * peripheral2 +
                          k41 * cord_venous
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(peripheral2) <-  k13 * central - k31 * peripheral2
    d/dt(cord_venous) <-  k14 * central - k41 * cord_venous -
                          k45 * cord_venous + k54 * neonate
    d/dt(neonate)     <-  k45 * cord_venous - k54 * neonate - k50 * neonate

    # Observations. Dose in mg divided by a volume in L gives mg/L, the
    # unit of every reported amoxicillin concentration in the paper.
    Cc       <- central     / vc
    Ccord    <- cord_venous / v4
    Cneonate <- neonate     / v5

    Cc       ~ prop(propSd)
    Ccord    ~ prop(propSd_Ccord)
    Cneonate ~ prop(propSd_Cneonate)
  })
}
