Olivo_2024_methotrexate <- function() {
  description <- paste(
    "Two-compartment population PK model for high-dose intravenous",
    "methotrexate (12 g/m^2 over a 4 h infusion) in Brazilian paediatric",
    "patients with osteosarcoma treated on the Brazilian Osteosarcoma",
    "Treatment Group (BOTG) protocol (Olivo 2024; n = 32 patients, 216",
    "cycles, 563 therapeutic-drug-monitoring concentrations). Linear",
    "first-order elimination from the central compartment. Serum creatinine",
    "scales clearance as a power function normalized to the cohort median",
    "0.58 mg/dL (exponent -0.192), and body surface area scales the central",
    "volume as a power function normalized to the cohort median 1.45 m^2",
    "(exponent 0.301). Correlated exponential between-subject variability",
    "on CL and Vc (correlation 94%), exponential between-occasion",
    "variability on CL shared across MTX cycles, and a proportional",
    "residual error. Built to anticipate leucovorin (folinic acid) rescue",
    "dose adjustments before the first monitored MTX concentration.",
    sep = " "
  )
  reference <- paste(
    "Olivo LB, de Oliveira Henz P, Wermann S, Dias BB, Porto GO,",
    "Pinhatti AV, Martins MD, Gregianin LJ, Costa TD, de Araujo BV",
    "(2024). Anticipating Leucovorin Rescue Therapy in Patients with",
    "Osteosarcoma through Methotrexate Population Pharmacokinetic Model.",
    "Pharmaceutics 16(9):1180. doi:10.3390/pharmaceutics16091180.",
    sep = " "
  )
  vignette <- "Olivo_2024_methotrexate"
  units    <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Olivo 2024 Results 3.3 ("a
  # two-compartment model with first-order elimination from the central
  # compartment parametrized in terms of clearance (CL), central compartment
  # volume (Vc), peripherical compartment volume (Vp), and
  # inter-compartmental clearance (Q)") and Methods 2.3 (MTX quantified in
  # serum by chemiluminescent microparticle immunoassay).
  compartmentData <- list(
    central     = list(analyte = "methotrexate", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "methotrexate", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    CREAT = list(
      description        = "Serum creatinine",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Source column SCr. Power-form effect on CL normalized to the cohort",
        "median, per Olivo 2024 Eq. 8:",
        "'CL_i = 14.8 x (SCr / 0.58)^-0.192 x e^(...)'. The 0.58 mg/dL",
        "reference is the cohort median serum creatinine reported in Table 2",
        "(median 0.58, range 0.17-3.2 mg/dL), consistent with Methods 2.3:",
        "'The analysis was performed using the linear, power, and exponential",
        "functions normalized by their respective medians.' Time-varying by",
        "cycle: Methods 2.1 records demographic and biochemical covariates",
        "'according to each cycle', missing continuous covariates were",
        "replaced by the nearest value by date for the patient, and",
        "time-varying continuous covariates were handled per Wahlby 2004.",
        "Adding SCr on CL lowered the OFV by 23.27 (p < 0.001) and reduced",
        "BSV and BOV on CL by 2.5% and 3.7% (Results 3.3). Creatinine",
        "clearance gave a similar OFV drop, so the primary biomarker SCr was",
        "retained instead; height (also used in the Schwartz CrCL equation)",
        "had no impact on BSV or BOV. Negative exponent: higher serum",
        "creatinine (worse renal function) lowers MTX clearance, which is",
        "mechanistically coherent because MTX elimination is mostly renal",
        "(Discussion)."
      ),
      source_name        = "SCr"
    ),
    BSA = list(
      description        = "Body surface area (Haycock equation)",
      units              = "m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power-form effect on Vc normalized to the cohort median, per Olivo",
        "2024 Eq. 9: 'Vc_i = 82.5 x (BSA / 1.45)^0.301 x e^(...)'. The",
        "1.45 m^2 reference is the cohort median BSA reported in Table 2",
        "(median 1.45, range 0.67-2.03 m^2). Table 2 footnote a: BSA",
        "'calculated with Haycock equation'. Adding BSA on Vc dropped the Vc",
        "variability by 11.3% (delta OFV 6.75, p < 0.01; Results 3.3 and",
        "Table S1). BSA is also the dosing descriptor for this protocol",
        "(12 g/m^2 per cycle), so a simulation must carry BSA for both the",
        "dose amount and the Vc covariate. Time-varying by cycle (Methods",
        "2.1). Note that the Conclusions section attributes the BSA effect to",
        "'Vp'; Table 3 (theta_BSA listed under Central Volume), Eq. 9, the",
        "Abstract and Results 3.3 all place it on Vc, which is what is",
        "encoded here."
      ),
      source_name        = "BSA"
    ),
    OCC = list(
      description        = "Integer-valued MTX cycle (occasion) indicator for between-occasion variability on CL",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = paste(
        "Olivo 2024 Methods 2.3 investigated between-occasion variability",
        "(BOV) and coded occasions per Ho Hui 2019 despite non-sequential",
        "cycles. Results 3.3: 'The inclusion of exponential BOV in CL (15.1%)",
        "decreased the OFV by 596.97 points ... Between-occasion variability",
        "was estimated as a single value for all different MTX occasions.'",
        "Because the BOV magnitude is a single shared variance, the occasion",
        "count only sets how many multiplexed slots exist. Twelve slots are",
        "provided here to cover a full BOTG course: Methods 2.1 states 'Each",
        "patient received 12 cycles of HDMTX by a 4-h infusion'. (The",
        "Discussion separately reports n = 21 cycles in the final dataset;",
        "extending past 12 occasions requires adding further etaiov_cl_<n>",
        "slots with the same fixed variance, following the",
        "Oosten_2016_fentanyl.R precedent.) Decomposed inside model() into",
        "binary indicators oc1..oc12; OCC = 0 or any value outside 1..12",
        "zeros every indicator and yields the typical-value CL with BSV only.",
        "For single-cycle simulation pass OCC = 1."
      ),
      source_name        = "OCC"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = paste(
        "Screened but not retained. Discussion: 'patient's age was not well",
        "distributed around all childhood lifetime stages, which probably",
        "explains why, in the current POPPK model, patient's age was not",
        "included as a covariate to explain PK parameter variability, as",
        "described in previous studies'. Cohort median 13.25 years (range",
        "5-18); 24 of 32 patients were adolescents (> 12 years)."
      )
    ),
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = paste(
        "Listed among the covariates tested in Methods 2.3 ('Covariates",
        "tested were weight, age, height, body surface area, body mass index,",
        "sex, race, serum creatinine, creatinine clearance, hepatic enzymes,",
        "hematocrit, and hemoglobin') but not retained in the final model;",
        "BSA carried the size effect on Vc. Table 2 median 47 kg (range",
        "13.80-85.50)."
      )
    ),
    HT = list(
      description = "Body height",
      units       = "cm",
      type        = "continuous",
      notes       = paste(
        "Screened on CL and not retained. Results 3.3: 'The investigation of",
        "height as a covariate for CL, which is also used in the Schwartz",
        "equation to calculate CrCL, had no impact on BSV or BOV.' Table 2",
        "median 159 cm (range 115-177)."
      )
    ),
    BMI = list(
      description = "Body mass index (Quetelet equation)",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = paste(
        "Tested per Methods 2.3 but not retained. Table 2 median 18.11",
        "kg/m^2 (range 10.43-28.56), footnote b 'calculated with Quetelet",
        "equation'."
      )
    ),
    SEXF = list(
      description = "Sex indicator (1 = female, 0 = male)",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Tested per Methods 2.3 but not retained. Results 3.1: 'Gender was",
        "well distributed in this population, comprising 18 males and 14",
        "females.'"
      )
    ),
    RACE_BLACK = list(
      description = "Black race indicator (1 = Black, 0 otherwise)",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Race was tested per Methods 2.3 but not retained. Results 3.1:",
        "'Around 78% of the patients were white, 18% black, and 3% were from",
        "another ethnicity.' The Discussion notes that cross-centre use of",
        "the model should be re-evaluated 'since Brazil is a country with",
        "significant ethnic diversity, which results in regional variation in",
        "patient's profiles'."
      )
    ),
    CRCL = list(
      description = "Creatinine clearance (Schwartz equation)",
      units       = "mL/min/1.73 m^2",
      type        = "continuous",
      notes       = paste(
        "Screened on CL and deliberately superseded by the primary biomarker.",
        "Results 3.3: 'The use of CrCL as a covariate on CL led to a similar",
        "decrease in the OFV as SCr, indicating that the primary biomarker",
        "was sufficient to explain the variability.' Table 2 median 190.28",
        "mL/min/1.73 m^2 (range 37.41-372.05), footnote c 'calculated with",
        "Schwartz equation'."
      )
    ),
    BUN = list(
      description = "Blood urea nitrogen",
      units       = "mg/dL",
      type        = "continuous",
      notes       = paste(
        "Retained as a CL covariate in the same group's earlier paediatric",
        "ALL model (Henz 2023) but NOT in this osteosarcoma model.",
        "Discussion: 'Differently from our model for MTX in ALL patients, in",
        "the present model, BUN did not explain MTX variability in CL,",
        "probably because 95% of the patients were well hydrated throughout",
        "the treatment.' Table 2 median 17 mg/dL (range 4-87)."
      )
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = paste(
        "Hepatic enzymes were tested but were not estimable here.",
        "Discussion: 'In the population used to build the present model,",
        "liver enzymes varied largely (from 6-8 to 1052 U/L), preventing the",
        "determination of their impact on PK parameters.' Table 2 median",
        "34 U/L (range 8-1052)."
      )
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = paste(
        "Tested with AST and not estimable for the same reason (Discussion).",
        "Table 2 median 47 U/L (range 6-1052)."
      )
    ),
    HCT = list(
      description = "Hematocrit",
      units       = "%",
      type        = "continuous",
      notes       = paste(
        "Tested per Methods 2.3 but not retained. Table 2 median 28.9%",
        "(range 10-52.20)."
      )
    ),
    HGB = list(
      description = "Hemoglobin",
      units       = "g/dL",
      type        = "continuous",
      notes       = paste(
        "Tested per Methods 2.3 but not retained. Table 2 median 9.6 g/dL",
        "(range 3.28-13.2)."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 32L,
    n_studies      = 1L,
    age_range      = "5 to 18 years",
    age_median     = "13.25 years (mean 13.25; 24 of 32 patients were adolescents > 12 years)",
    weight_range   = "13.80 to 85.50 kg",
    weight_median  = "47 kg",
    height_range   = "115 to 177 cm (median 159)",
    bsa_range      = "0.67 to 2.03 m^2 (median 1.45; Haycock equation)",
    sex_female_pct = 100 * 14 / 32,
    race_ethnicity = c(White = 78.1, Black = 18.8, Other = 3.1),
    disease_state  = paste(
      "Paediatric osteosarcoma treated on the Brazilian Osteosarcoma",
      "Treatment Group (BOTG) protocol: MAP chemotherapy (high-dose",
      "methotrexate + doxorubicin + cisplatin) for weeks 1-10, surgery in",
      "week 11 or 12, then a repeat of the 10-week cycle after surgery. Of",
      "the 216 analysed cycles, 96 were pre-tumour-resection and 120",
      "post-resection."
    ),
    renal_function = paste(
      "Serum creatinine median 0.58 mg/dL (range 0.17-3.2); Schwartz",
      "creatinine clearance median 190.28 mL/min/1.73 m^2 (range",
      "37.41-372.05). Patients received hydration fluids 3000 mL/m^2/day",
      "and sodium bicarbonate 45 mEq/m^2 for urine alkalization before",
      "(200 mL/m^2) and during (500 mL/m^2) each cycle; urinary pH was",
      "median 7.5 (range 7-9). Observations beyond 96 h were removed",
      "because in extensive TDM patients are selected for dialysis."
    ),
    dose_range     = paste(
      "High-dose methotrexate 12 g/m^2 per cycle as a 4 h intravenous",
      "infusion, 12 cycles per patient; realised dose median 11.9 g/m^2",
      "(range 5.9-12.9 g/m^2 depending on the patient's clinical",
      "condition)."
    ),
    regions        = "Brazil (Hospital de Clinicas de Porto Alegre, Rio Grande do Sul; southern Brazilian public reference hospital).",
    notes          = paste(
      "Retrospective therapeutic-drug-monitoring data collected between",
      "January 2015 and March 2023. Cycles were included if they had at",
      "least 2 observed MTX concentrations and the patient was monitored up",
      "to 96 h: 216 cycles and 563 serum MTX concentrations (median 3 per",
      "cycle, range 2-6); 147 cycles were monitored to 72 h after the end of",
      "infusion. Sampling time was missing for 65.3% of observations and was",
      "single-imputed using the mean 32.5 +/- 4.7 min difference between",
      "laboratory arrival time and reported sampling time (technique of Henz",
      "2023). Assay: chemiluminescent microparticle immunoassay (Alinity i",
      "Methotrexate Reagent Kit 09P48), LLOQ 0.013 mg/L, no observations",
      "below LLOQ. NONMEM 7.4 with FOCE-I plus PsN 4.9.0; evaluated by GOF",
      "plots, a 1000-simulation VPC, shrinkage (< 30%), %RSE and a 1000-run",
      "nonparametric bootstrap (853/1000 successful). Baseline demographics",
      "from Table 2; final parameter estimates from Table 3; the final",
      "individual-parameter covariate model from Eqs. 8-11. Two further",
      "screened covariates have no canonical column and are therefore not",
      "listed in covariatesDataExcluded: urinary pH (Table 2 median 7.5,",
      "range 7-9), excluded because 'urine pH values were mostly alkaline",
      "over the treatment period in all patients, impeding their use to",
      "explain CL variability among patients or cycles' (Discussion); and a",
      "dichotomous pre- versus post-tumour-resection period indicator",
      "(Methods 2.3), also not retained. Genetic polymorphisms were not",
      "investigated in this study. Because all patients followed the same",
      "chemotherapy protocol, interactions between the anticancer drugs",
      "could not be evaluated."
    )
  )

  ini({
    # Structural parameters -- Olivo 2024 Table 3 'Estimate' column, applied in
    # the individual-parameter forms printed as Eqs. 8-11:
    #
    #   Eq. 8   CL_i = 14.8  x (SCr / 0.58)^-0.192 x e^(0.0391 + 0.0228)   L/h
    #   Eq. 9   Vc_i = 82.5  x (BSA / 1.45)^0.301  x e^(0.0181)            L
    #   Eq. 10  Q_i  = 0.178                                              L/h
    #   Eq. 11  Vp_i = 5.72                                               L
    #
    # The exponential terms printed in Eqs. 8-9 are the omega variances of the
    # random effects (see the IIV / IOV block below), not covariate terms.
    lcl <- log(14.8);  label("Clearance at the reference serum creatinine of 0.58 mg/dL (L/h)")     # Table 3 TVCL = 14.8 L/h (RSE 19%; bootstrap median 14.4, 95% CI 10.7-19.3); Eq. 8
    lvc <- log(82.5);  label("Central volume of distribution at the reference BSA of 1.45 m^2 (L)") # Table 3 TVVC = 82.5 L (RSE 23%; bootstrap median 79.9, 95% CI 54.5-115.9); Eq. 9
    lq  <- log(0.178); label("Intercompartmental clearance (L/h)")                                  # Table 3 TVQ = 0.178 L/h (RSE 31%; bootstrap median 0.171, 95% CI 0.106-0.285); Eq. 10
    lvp <- log(5.72);  label("Peripheral volume of distribution (L)")                               # Table 3 TVVP = 5.72 L (RSE 35%; bootstrap median 5.42, 95% CI 3.1-10.5); Eq. 11

    # Covariate effects, both power-form and normalized by the cohort median of
    # the covariate (Methods 2.3: 'linear, power, and exponential functions
    # normalized by their respective medians').
    e_creat_cl <- -0.192; label("Power exponent of serum creatinine on CL (unitless)")       # Table 3 theta_SCr = -0.192 (RSE 31%; bootstrap median -0.201, 95% CI -0.32 to -0.08); Eq. 8
    e_bsa_vc   <-  0.301; label("Power exponent of body surface area on Vc (unitless)")      # Table 3 theta_BSA = 0.301 (RSE 36%; bootstrap median 0.298, 95% CI 0.025-0.481); Eq. 9

    # Between-subject variability. Table 3 reports omega as a percentage; those
    # percentages are the log-scale SDs rather than CV%, which the exponential
    # terms printed in Eqs. 8-9 confirm exactly: (0.198)^2 = 0.0391 and
    # (0.135)^2 = 0.0181 as printed, whereas log(1 + CV^2) would give 0.0385
    # and 0.0181. The variances below are therefore taken verbatim from
    # Eqs. 8-9. Results 3.3: 'BSV was exponentially added in CL (19.8%) and Vc
    # (13.5%), with a correlation of 94%.' The off-diagonal is
    # 0.94 x sqrt(0.0391 x 0.0181) = 0.0250.
    etalcl + etalvc ~ c(0.0391,
                        0.0250, 0.0181)   # Eq. 8 e^(0.0391 + ...) and Eq. 9 e^(0.0181); Table 3 omega BSV CL 19.8% (RSE 16%) and Vc 13.5% (RSE 29%), correlation CL-Vc 94%

    # Between-occasion variability on log-CL, a single shared variance across
    # MTX cycles (Results 3.3: 'Between-occasion variability was estimated as a
    # single value for all different MTX occasions'). nlmixr2 has no NONMEM
    # `$OMEGA BLOCK(1) SAME` shortcut, so occasions 2..12 each get their own eta
    # with the variance fixed equal to the occasion-1 estimate
    # (Jonsson_2011_ethambutol.R / Oosten_2016_fentanyl.R precedent).
    etaiov_cl_1  ~ 0.0228          # Eq. 8 e^(... + 0.0228); Table 3 omega BOV CL 15.1% (RSE 6%; bootstrap median 14.8, 95% CI 12.8-16.9) -> 0.151^2 = 0.0228
    etaiov_cl_2  ~ fixed(0.0228)   # SAME-equivalent: equal to the occasion-1 BOV variance
    etaiov_cl_3  ~ fixed(0.0228)
    etaiov_cl_4  ~ fixed(0.0228)
    etaiov_cl_5  ~ fixed(0.0228)
    etaiov_cl_6  ~ fixed(0.0228)
    etaiov_cl_7  ~ fixed(0.0228)
    etaiov_cl_8  ~ fixed(0.0228)
    etaiov_cl_9  ~ fixed(0.0228)
    etaiov_cl_10 ~ fixed(0.0228)
    etaiov_cl_11 ~ fixed(0.0228)
    etaiov_cl_12 ~ fixed(0.0228)

    # Residual error. Results 3.3: 'The residual variability was described by a
    # proportional error model.'
    propSd <- 0.309; label("Proportional residual error (fraction)")  # Table 3 'Residual Variability Proporcional Error' = 30.9% (RSE 4%; bootstrap median 30.3, 95% CI 27.5-33.1)
  })

  model({
    # Decompose the integer MTX-cycle column into binary occasion indicators
    # that multiplex the BOV etas onto log-CL. OCC = 1..12 selects the matching
    # per-occasion eta; OCC = 0 or any value outside 1..12 zeros every
    # indicator and leaves CL with between-subject variability only.
    oc1  <- (OCC == 1)
    oc2  <- (OCC == 2)
    oc3  <- (OCC == 3)
    oc4  <- (OCC == 4)
    oc5  <- (OCC == 5)
    oc6  <- (OCC == 6)
    oc7  <- (OCC == 7)
    oc8  <- (OCC == 8)
    oc9  <- (OCC == 9)
    oc10 <- (OCC == 10)
    oc11 <- (OCC == 11)
    oc12 <- (OCC == 12)
    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2 + oc3 * etaiov_cl_3 +
      oc4 * etaiov_cl_4 + oc5 * etaiov_cl_5 + oc6 * etaiov_cl_6 +
      oc7 * etaiov_cl_7 + oc8 * etaiov_cl_8 + oc9 * etaiov_cl_9 +
      oc10 * etaiov_cl_10 + oc11 * etaiov_cl_11 + oc12 * etaiov_cl_12

    # Individual PK parameters, Eqs. 8-11. Both covariate effects are power
    # functions of the covariate divided by its cohort median. Q and Vp carry
    # neither a covariate nor a random effect.
    cl <- exp(lcl + etalcl + iov_cl) * (CREAT / 0.58)^e_creat_cl
    vc <- exp(lvc + etalvc)          * (BSA / 1.45)^e_bsa_vc
    q  <- exp(lq)
    vp <- exp(lvp)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Two-compartment disposition; methotrexate is given as an intravenous
    # infusion directly into the central compartment (no absorption phase).
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Serum methotrexate concentration in mg/L (dose in mg, volumes in L).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
