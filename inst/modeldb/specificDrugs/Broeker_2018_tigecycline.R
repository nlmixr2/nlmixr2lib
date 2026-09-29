Broeker_2018_tigecycline <- function() {
  description <- paste(
    "Two-compartment population PK model for intravenous tigecycline in 11",
    "critically ill adults with acute kidney injury receiving continuous renal",
    "replacement therapy, 8 on continuous venovenous hemodialysis (CVVHD) and 3",
    "on continuous venovenous hemodiafiltration (CVVHDF) (Broeker 2018).",
    "Plasma and CRRT effluent concentrations were fitted SIMULTANEOUSLY, which",
    "is what separates the dialysis clearance from the physiological body",
    "clearance. Total elimination from the central compartment is the ADDITIVE",
    "sum of a body arm (CLbody, 18.3 L/h typical, encoded as lcl, carrying a",
    "median-normalised power effect of total bilirubin with a NEGATIVE exponent",
    "so that cholestatic patients clear tigecycline more slowly) and a dialysis",
    "arm whose typical value depends on the CRRT modality (CLdial 1.69 L/h for",
    "CVVHD, 2.71 L/h for CVVHDF). The dialysis arm is gated by",
    "RRT_CRRT_ACTIVE and selected by RRT_CVVHDF_STATUS, with CVVHD as the",
    "reference modality; only the CVVHD arm carries interindividual",
    "variability, an IIV on the CVVHDF arm having collapsed to zero during",
    "estimation. Besides the plasma concentration Cc the model returns the",
    "flow-normalised effluent concentration Ceffluent with its own proportional",
    "residual error, obtained by rearranging the paper's first-principle",
    "dialysis equations. Age, sex, serum creatinine and Cockcroft-Gault",
    "creatinine clearance were screened on the body clearance and not retained,",
    "and allometric scaling by total body weight did not improve the fit."
  )
  reference <- paste(
    "Broeker A, Wicha SG, Dorn C, Kratzer A, Schleibinger M, Kees F,",
    "Heininger A, Kees MG, Haeberle H.",
    "Tigecycline in critically ill patients on continuous renal replacement",
    "therapy: a population pharmacokinetic study.",
    "Crit Care. 2018;22(1):341.",
    "doi:10.1186/s13054-018-2278-4"
  )
  vignette <- "Broeker_2018_tigecycline"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "tigecycline", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "tigecycline", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    TBILI = list(
      description = "Total serum bilirubin",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "The ONLY covariate retained in the final model. Enters the body clearance as a",
        "median-normalised power term, Broeker 2018 Table 2 row heading:",
        "'Clearance (L/h) = theta_1 x (bilirubin/2.3)^theta_2' with theta_1 = 18.3 L/h and",
        "theta_2 = -0.29. UNITS: the paper reports bilirubin in the US convention mg/dL",
        "(Table 1 column header 'Bilirubin (mg/dL)') and the reference value 2.3 is",
        "therefore 2.3 mg/dL, stated in the Results as 'normalized by the population median",
        "of bilirubin, 2.3 mg/dL' and confirmed by the Table 1 median row. This column",
        "carries the register's canonical SI umol/L, so model() converts back inline with",
        "TBILI / 17.1 before forming the ratio; equivalently the reference is 39.3 umol/L.",
        "The exponent is NEGATIVE, so higher bilirubin means LOWER clearance -- Results:",
        "'lower bilirubin concentrations corresponded to higher clearances'. The effect is",
        "large because the cohort spans two orders of magnitude of bilirubin (0.7 to 43.3",
        "mg/dL, Table 1) in a population with four liver-failure or cirrhosis patients and",
        "two liver transplants: the Results quantify it as 'Individual clearance values",
        "varied from 9.3 L/h (10th percentile) to 19.1 L/h (90th percentile) depending on",
        "the bilirubin concentration (24 mg/dL to 1.8 mg/dL)', which the encoded term",
        "reproduces as 18.3 * (24/2.3)^-0.29 = 9.3 and 18.3 * (1.8/2.3)^-0.29 = 19.7.",
        "Including it dropped the objective function by 5.71 (p = 0.017) and cut the",
        "unexplained interindividual variability on clearance from 58.6% to 43.6%, so the",
        "43.6% in ini() is the POST-covariate value. The authors read bilirubin as a marker",
        "of hepatic function, tigecycline being predominantly biliary/faecally eliminated.",
        "Treated as a time-fixed baseline value (a single day-4 laboratory value per",
        "patient)."
      ),
      source_name = "bilirubin"
    ),
    RRT_CRRT_ACTIVE = list(
      description = "CRRT-active indicator (1 while the continuous renal replacement circuit is running, 0 otherwise)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (circuit not running)",
      notes = paste(
        "Gates the dialysis clearance arm and, with it, the effluent observable. EVERY",
        "subject in Broeker 2018 was on CRRT continuously for the whole sampling period --",
        "blood was drawn 'on day 4 of treatment with tigecycline after at least 24 h on",
        "CRRT' (Methods, Sampling and drug analysis) and the study has no off-CRRT",
        "subgroup -- so this column is identically 1 throughout the source analysis and the",
        "paper estimates no off-CRRT clearance. It is carried anyway, per the",
        "RRT_CRRT_EFFLUENT_FLOW register entry's instruction to pair the effluent flow with",
        "an on/off gate, so that a downstream user can simulate circuit interruption or",
        "discontinuation; setting it to 0 leaves the body clearance alone, which is the",
        "paper's CLbody and NOT a validated off-CRRT clearance for this population. The",
        "ACTIVE rather than the STATUS member of the RRT_<modality>_<kind> family is used",
        "because the quantity is physically time-varying (filter changes interrupt the",
        "circuit), matching ButraguenoLaiseca_2022_piperacillin."
      ),
      source_name = "(not a source column; implicitly 1 for all subjects)"
    ),
    RRT_CVVHDF_STATUS = list(
      description = "CRRT modality indicator (1 = continuous venovenous hemodiafiltration, 0 = continuous venovenous hemodialysis)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (CVVHD, the reference modality carrying 8 of the 11 subjects)",
      notes = paste(
        "Selects which of the two published dialysis clearances applies. Broeker 2018",
        "Table 2 estimates them as two separate parameters -- 'Dialysis clearance CVVHD",
        "(L/h) 1.69' and 'Dialysis clearance CVVHDF (L/h) 2.71' -- over a cohort that",
        "shares one body-PK model: 'Eleven patients ... receiving either continuous",
        "veno-venous hemodialysis (CVVHD, n = 8) or hemodiafiltration (CVVHDF, n = 3)'",
        "(Abstract, Methods). Per-subject modality assignment is recoverable from Table 1,",
        "whose footnotes read 'a Continuous veno-venous hemodialysis (CVVHD); b continuous",
        "veno-venous hemodiafiltration (CVVHDF)': patients 1, 2, 3, 6, 8, 9, 10 and 11 are",
        "CVVHD and patients 4, 5 and 7 are CVVHDF. MEANINGFUL ONLY WHEN RRT_CRRT_ACTIVE =",
        "1; the pair (RRT_CRRT_ACTIVE, RRT_CVVHDF_STATUS) is what distinguishes off-circuit",
        "from CVVHD, which a modality indicator alone cannot do. The two modalities differ",
        "physically in that CVVHDF adds convective transport through a post-dilution",
        "ultrafiltrate stream (QFil = 1 L/h, Methods) on top of CVVHD's purely diffusive",
        "dialysate stream, which is why its dialysis clearance is the larger of the two and",
        "why its mean saturation coefficient is higher (0.90 versus 0.79, Results). Only",
        "the CVVHD arm carries IIV: Table 2 prints a dash for the CVVHDF arm and the",
        "Results explain that 'an IIV for this method was not supported by the data (IIV",
        "tended to zero during estimation)', unsurprising at n = 3."
      ),
      source_name = "Table 1 footnote markers a / b"
    ),
    RRT_CRRT_EFFLUENT_FLOW = list(
      description = "Total effluent flow rate leaving the CRRT circuit (dialysate + ultrafiltrate)",
      units = "mL/h",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters ONLY the effluent observation equation, never the clearance model --",
        "Broeker 2018 estimates the dialysis clearance directly rather than as a sieving",
        "coefficient times a flow. Broeker 2018 Eq. 3 defines CLDial,CVVHD = QDial *",
        "(Ceff/CPla) and Eq. 4 defines CLDial,CVVHDF = (QDial + QFil) * (Ceff/CPla), where",
        "'QDial represents the dialysate flow rate, Ceff represents the concentration of",
        "tigecycline in the effluent, and CPla represents the pre-filter plasma",
        "concentration'. Rearranging for the observable gives Ceff = CLDial * CPla / Qeff",
        "with Qeff = QDial for CVVHD and QDial + QFil for CVVHDF -- exactly the form used",
        "by ButraguenoLaiseca_2022_piperacillin, so this column is the single Qeff",
        "denominator for both modalities. VALUES: the Methods prescribe the dialysate flow",
        "by weight band, 'Blood flow and dialysate flow were adjusted to body weight (< 90",
        "kg/> 90 kg; 100/120 mL/min and 2000/2500 mL/h, respectively)', and fix the",
        "CVVHDF ultrafiltration rate at 'QFil was 1 L/h'. The weight-band rule gives 2000",
        "mL/h for all eight CVVHD patients (all under 90 kg, Table 1) and 3000 mL/h for a",
        "CVVHDF patient under 90 kg. That reading is confirmed ARITHMETICALLY by the",
        "paper's own derived statistic: the published CVVHDF saturation coefficient is",
        "0.90 (Results) and 2.71 / 3.000 = 0.903. The CVVHD counterpart, 1.69 / 2.000 =",
        "0.845, sits inside the published mean of 0.79 plus or minus an SD of 0.36. The",
        "dialysate flow for the CVVHDF arm is not restated in the CVVHDF paragraph of the",
        "Methods and is assumed to follow the same weight band; see the vignette",
        "Assumptions section. Meaningful only when RRT_CRRT_ACTIVE = 1; converted to L/h",
        "inside model()."
      ),
      source_name = "QDial, QFil"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened on the body clearance and not retained (Broeker 2018 Methods, Pharmacometric analysis: 'Age, sex, serum creatinine, creatinine clearance (Cockcroft-Gault), and bilirubin were tested as covariates on the body clearance'; only bilirubin survived the likelihood-ratio criterion). Table 1: median 69 years, range 37 to 81. The protocol excluded patients over 85 or under 18."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      notes = "Screened on the body clearance and not retained (Methods, Pharmacometric analysis). The cohort is 10 male and 1 female (Table 1), so the covariate is effectively unidentifiable here regardless of any true effect."
    ),
    SCR = list(
      description = "Serum creatinine",
      units = "mg/dL",
      type = "continuous",
      notes = "Screened on the body clearance and not retained (Methods, Pharmacometric analysis). Table 1: median 1.2 mg/dL, range 0.5 to 2.4. Tigecycline is not renally eliminated to a meaningful extent, and in an anuric CRRT-dependent cohort serum creatinine reflects the dialysis prescription rather than native renal function."
    ),
    CRCL = list(
      description = "Cockcroft-Gault creatinine clearance",
      units = "mL/min",
      type = "continuous",
      notes = "Screened on the body clearance and not retained (Methods, Pharmacometric analysis, which names the Cockcroft-Gault estimator explicitly). Not tabulated per patient; by convention Cockcroft-Gault is not interpretable in a CRRT-dependent patient, which the null result is consistent with."
    ),
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      notes = "Tested as an allometric scalar on the structural parameters with both a fixed and a freely estimated exponent and NOT retained -- Methods: 'Allometric scaling models using total body weight with fixed and freely estimated scaling parameters were evaluated'; Results: 'Allometric scaling with a fixed exponent did not improve the model significantly and was not included.' Table 1: median 80 kg, range 68 to 104. Body weight nonetheless remains indirectly load-bearing through the CRRT prescription, because the Methods set the dialysate flow by a 90 kg weight band -- that pathway is carried by RRT_CRRT_EFFLUENT_FLOW, not by any covariate effect on a structural parameter."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 11L,
    n_studies = 1L,
    n_observations = 217L,
    age_range = "37 to 81 years",
    age_median = "69 years",
    weight_range = "68 to 104 kg",
    weight_median = "80 kg",
    sex_female_pct = 9.1,
    race_ethnicity = "Not reported.",
    disease_state = paste(
      "Critically ill adults in a 40-bed anaesthesiological ICU of a tertiary care",
      "hospital who required renal replacement therapy for acute kidney injury and were",
      "treated with tigecycline. Ten of the eleven were treated for complicated",
      "intra-abdominal infection and one for an Acinetobacter baumannii infection.",
      "Relevant co-conditions were liver failure or cirrhosis (four patients), liver",
      "transplantation (two), and extracorporeal membrane oxygenation (one). Two patients",
      "died during follow-up. APACHE II median 29, range 15 to 45 (Table 1). Exclusion",
      "criteria were age over 85 or under 18 years, severe liver insufficiency",
      "(Child-Pugh C), acute pancreatitis, concomitant anticoagulation therapy, and a",
      "history of tigecycline allergy. EudraCT 2012-005617-39."
    ),
    dose_range = paste(
      "Standard tigecycline dosing: a 100 mg intravenous loading dose followed by 50 mg",
      "twice daily (Methods, Setting and study population). Sampling took place on day 4",
      "of treatment, after at least 24 h on CRRT, and is therefore at steady state.",
      "Samples were drawn immediately before the start of the infusion (time 0), at 1 h",
      "(the end of infusion), and at 1.25, 1.5, 1.75, 2, 4, 6, 8 and 12 h, with effluent",
      "collected from the circuit effluent port at the same time points. The 1 h",
      "infusion duration is implied by the 'after 1 h (i.e., the end of infusion)'",
      "sampling description rather than stated as a prescription."
    ),
    regions = "Single centre, Germany (University Hospital Tuebingen / University Hospital Regensburg collaboration).",
    renal_function = paste(
      "All patients had acute kidney injury requiring continuous renal replacement",
      "therapy: 8 on CVVHD and 3 on CVVHDF, all using the Fresenius MultiFiltrate system",
      "with an Ultraflux AV 1000 S polysulfone membrane. CVVHD used Ci-Ca Dialysate K2",
      "with 4% sodium citrate regional anticoagulation at a median citrate flow of 176",
      "mL/h (under 3% of blood flow), targeting a post-filter ionised calcium of 0.25 to",
      "0.35 mmol/L. CVVHDF used multiBic fluid for both dialysis and post-filter",
      "(post-dilution) replacement with an ultrafiltration rate of 1 L/h and unfractionated",
      "heparin anticoagulation. Blood and dialysate flows were set by weight band: 100",
      "mL/min and 2000 mL/h under 90 kg, 120 mL/min and 2500 mL/h above. Serum creatinine",
      "median 1.2 mg/dL (0.5 to 2.4). No predilution correction of the CRRT clearance was",
      "applied, the citrate flow being very low relative to blood flow."
    ),
    hepatic_function = paste(
      "Deliberately broad, which is what identifies the bilirubin covariate. Total",
      "bilirubin median 2.3 mg/dL with a range of 0.7 to 43.3 mg/dL; albumin 2.1 to 3.1",
      "g/dL; total protein 3.5 to 6.4 g/dL (Table 1). Four patients had liver failure or",
      "cirrhosis and two were liver transplant recipients, while Child-Pugh C liver",
      "insufficiency was an exclusion criterion."
    ),
    notes = paste(
      "Baseline demographics from Broeker 2018 Table 1, which lists all eleven patients",
      "individually. A total of 109 blood and 108 effluent samples were used, after",
      "excluding two 12-h blood samples whose very high concentrations indicated the draw",
      "followed the start of the next infusion. Total and free tigecycline were measured",
      "by validated HPLC-UV, with a plasma limit of quantification of 0.05 mg/L and",
      "intra- and interassay imprecision under 6%; the corresponding effluent values were",
      "0.025 mg/L and under 9%. Free concentrations were determined by ultrafiltration at",
      "1, 2 and 12 h, giving a median unbound fraction of 61% (range 45 to 94%) -- the",
      "model is fitted to TOTAL concentrations and carries no protein-binding term. NONMEM",
      "7.4 with FOCEI executed via PsN 4.5.16, ADVAN1 and ADVAN3 routines; model selection",
      "by likelihood-ratio test (dOFV > 3.84), AIC for non-nested models, goodness-of-fit",
      "plots and visual predictive checks (n = 1000); parameter uncertainty from a",
      "nonparametric bootstrap (n = 1000). Shrinkage of the individual parameters was at",
      "most 26%."
    )
  )

  ini({
    # -----------------------------------------------------------------------
    # STRUCTURAL PARAMETERS -- Broeker 2018 Table 2, "Estimate" column, with
    # the nonparametric-bootstrap 95% CI quoted from the "95% CI" column.
    # A two-compartment model with first-order disposition was selected over
    # a one-compartment model by a dOFV of -113.77 (Results).
    #
    # Total elimination from the central compartment is the ADDITIVE sum of
    # the body arm below and a modality-specific dialysis arm. The paper's
    # own decomposition: 'the contribution of CRRT to tigecycline clearance
    # (CL) was only moderate (CLCVVHD: 1.69 L/h, CLCVVHDF: 2.71 L/h) in
    # comparison with CLbody (physiological part of the total clearance) of
    # 18.3 L/h' (Abstract).
    # -----------------------------------------------------------------------
    lcl <- log(18.3)
    label("Typical body clearance CLbody at the median bilirubin of 2.3 mg/dL (L/h)") # Table 2: theta_1 = 18.3 L/h (RSE 11.0%, bootstrap 95% CI 13.2, 22.7)
    lvc <- log(58.7)
    label("Typical central volume of distribution V1 (L)") # Table 2: V1 = 58.7 L (RSE 21.3%, 95% CI 29.3, 101.6)
    lvp <- log(154)
    label("Typical peripheral volume of distribution V2 (L)") # Table 2: V2 = 154 L (RSE 9.5%, 95% CI 124.3, 196.8)
    lq <- log(56.4)
    label("Typical distribution clearance Q (L/h)") # Table 2: Q = 56.4 L/h (RSE 15.3%, 95% CI 41.1, 76.6)

    # -----------------------------------------------------------------------
    # DIALYSIS CLEARANCE, one parameter per CRRT modality. Table 2 estimates
    # these as two independent parameters rather than as a reference value
    # plus a multiplier, so both printed numbers are transcribed verbatim and
    # RRT_CVVHDF_STATUS selects between them in model().
    #
    # These are the values the paper's first-principle dialysis model (Eq. 3
    # and Eq. 4) returns once the effluent measurements are added; the
    # descriptive saturation coefficients 0.79 (CVVHD) and 0.90 (CVVHDF)
    # quoted in the Results are Ceff/CPla ratios derived FROM these, not
    # separate estimated parameters.
    # -----------------------------------------------------------------------
    lcl_hemodialysis_cvvhd <- log(1.69)
    label("Typical dialysis clearance CLdial on CVVHD (L/h)") # Table 2: 'Dialysis clearance CVVHD (L/h)' = 1.69 (RSE 15.4%, 95% CI 1.26, 2.27)
    lcl_hemodialysis_cvvhdf <- log(2.71)
    label("Typical dialysis clearance CLdial on CVVHDF (L/h)") # Table 2: 'Dialysis clearance CVVHDF (L/h)' = 2.71 (RSE 8.9%, 95% CI 2.31, 3.16)

    # -----------------------------------------------------------------------
    # COVARIATE EFFECT -- Table 2 row heading:
    #   Clearance (L/h) = theta_1 x (bilirubin / 2.3)^theta_2
    # The reference 2.3 is the population MEDIAN bilirubin in mg/dL (Results;
    # Table 1 median row). The exponent is negative: higher bilirubin gives
    # lower clearance.
    # -----------------------------------------------------------------------
    e_tbili_cl <- -0.29
    label("Power exponent of median-normalised total bilirubin on body clearance (unitless)") # Table 2: theta_2 = -0.29 (RSE 33.1%, 95% CI -0.68, -0.10)

    # -----------------------------------------------------------------------
    # INTERINDIVIDUAL VARIABILITY -- Table 2, 'Interindividual variability
    # (%CV)' column. Eq. 1 gives the standard exponential model
    # P_k,i = theta_k * exp(eta_k,i), 'assuming log-normal distribution'.
    #
    # SCALE OF THE %CV COLUMN. The variances below are (CV/100)^2, i.e. the
    # printed %CV is read as 100 x the log-scale SD omega. Broeker 2018
    # prints no conversion footnote, but the Table 2 caption states that RSE
    # is 'reported on standard deviation scale for variability parameters',
    # and the same column reports the PROPORTIONAL residual errors (16.9%
    # and 40.6%) where a proportional sigma's %CV is unambiguously 100 x the
    # SD. Reading one column two ways would be inconsistent. The alternative
    # convention omega^2 = log(1 + CV^2) would give 0.174 / 0.802 / 0.161 /
    # 0.173 instead of the values below -- a difference of about 5% in omega
    # for the three moderate etas, though a larger one for V1. The vignette
    # Assumptions section records the choice.
    #
    # Only clearance, central volume and distribution clearance carry IIV:
    # 'The best model included IIV on clearance, central volume of
    # distribution and intercompartmental clearance' (Results). Table 2
    # prints a dash for V2 and for the CVVHDF dialysis arm. Shrinkage of the
    # individual parameters was at most 26%. The etas are independent: Table
    # 2 reports no off-diagonal element.
    #
    # The 43.6% on clearance is the POST-covariate value -- adding bilirubin
    # 'reduced the observed interindividual variability on CLbody from 58.6%
    # to 43.6%' (Abstract).
    # -----------------------------------------------------------------------
    etalcl ~ 0.436^2 # Table 2: IIV on clearance = 43.6 %CV
    etalvc ~ 1.109^2 # Table 2: IIV on V1 = 110.9 %CV
    etalq ~ 0.418^2 # Table 2: IIV on Q = 41.8 %CV
    etalcl_hemodialysis_cvvhd ~ 0.435^2 # Table 2: IIV on CVVHD dialysis clearance = 43.5 %CV

    # -----------------------------------------------------------------------
    # RESIDUAL ERROR -- Table 2, 'Residual variability' rows. Eq. 2 gives the
    # combined model Y_OBS = Y_PRED * (1 + eps_p) + eps_a, but the additive
    # component was dropped: 'A combined residual variability model
    # (proportional and additive) was not supported (additive error tended to
    # zero), so a proportional residual variability model was chosen'
    # (Results). Plasma and effluent carry separate magnitudes, the effluent
    # being the noisier matrix.
    # -----------------------------------------------------------------------
    propSd <- 0.169
    label("Proportional residual error for pre-filter plasma concentrations (fraction)") # Table 2: sigma proportional, pre-filter plasma = 16.9 %CV (RSE 16.1%, 95% CI 10.9, 21.7)
    propSd_Ceffluent <- 0.406
    label("Proportional residual error for effluent concentrations (fraction)") # Table 2: sigma proportional, effluent = 40.6 %CV (RSE 13.1%, 95% CI 30.4, 50.2)
  })

  model({
    # -------------------------------------------------------------------
    # 1. Body clearance and its bilirubin covariate.
    #
    # Table 2: Clearance (L/h) = theta_1 * (bilirubin / 2.3)^theta_2, with
    # the reference 2.3 in the paper's mg/dL. TBILI carries the register's
    # canonical SI umol/L, so it is converted back inline first.
    # -------------------------------------------------------------------
    tbili_mgdL <- TBILI / 17.1 # SI umol/L -> US-convention mg/dL (1 mg/dL = 17.1 umol/L)

    cl_body <- exp(lcl + etalcl) * (tbili_mgdL / 2.3)^e_tbili_cl

    # -------------------------------------------------------------------
    # 2. Dialysis clearance arm -- Table 2 gives one estimate per modality.
    #
    # RRT_CRRT_ACTIVE gates the whole arm off when the circuit is not
    # running; RRT_CVVHDF_STATUS then selects the modality, with CVVHD as
    # the reference. Only the CVVHD arm carries IIV -- for CVVHDF, 'an IIV
    # for this method was not supported by the data (IIV tended to zero
    # during estimation)' (Results).
    # -------------------------------------------------------------------
    cl_hemodialysis <-
      RRT_CRRT_ACTIVE * (
        (1 - RRT_CVVHDF_STATUS) *
          exp(lcl_hemodialysis_cvvhd + etalcl_hemodialysis_cvvhd) +
          RRT_CVVHDF_STATUS * exp(lcl_hemodialysis_cvvhdf)
      )

    # -------------------------------------------------------------------
    # 3. Total clearance, volumes and micro-constants.
    #
    # `cl` is deliberately the TOTAL elimination clearance (body + dialysis)
    # rather than the body arm alone, so that the cl / vc / q / vp quadruple
    # remains a faithful description of the system even if rxode2 elects to
    # solve it analytically.
    # -------------------------------------------------------------------
    cl <- cl_body + cl_hemodialysis
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp)
    q <- exp(lq + etalq)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 4. ODE system -- two-compartment, intravenous infusion into central
    #    (ADVAN3; Methods, Pharmacometric analysis).
    d/dt(central) <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # -------------------------------------------------------------------
    # 5. Observations -- plasma and effluent, fitted simultaneously.
    #
    # Cc is the pre-filter plasma concentration CPla of Eq. 3 and Eq. 4.
    #
    # EFFLUENT. The paper's first-principle dialysis equations are
    #   Eq. 3  CLdial,CVVHD  =  QDial          * (Ceff / CPla)
    #   Eq. 4  CLdial,CVVHDF = (QDial + QFil)  * (Ceff / CPla)
    # which rearrange to a single observable equation once the modality's
    # effluent flow is supplied as one column:
    #   Ceff = CLdial * CPla / Qeff
    # This is the same flow-normalised form used by
    # ButraguenoLaiseca_2022_piperacillin. RRT_CRRT_EFFLUENT_FLOW is stored
    # in the canonical mL/h and converted to the L/h in which clearances are
    # expressed; the floor keeps the denominator finite for a subject off
    # the circuit, where cl_hemodialysis is zero and the effluent
    # concentration is zero anyway.
    # -------------------------------------------------------------------
    effluent_flow <- max(RRT_CRRT_EFFLUENT_FLOW / 1000, 0.001)

    Cc <- central / vc
    Ceffluent <- cl_hemodialysis * Cc / effluent_flow

    Cc ~ prop(propSd)
    Ceffluent ~ prop(propSd_Ceffluent)
  })
}
