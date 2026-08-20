VarelaGonzalezAller_2025_fludarabine <- function() {
  description <- paste(
    "Three-compartment population PK model for intravenous fludarabine",
    "(measured as the circulating metabolite F-ara-A) given as",
    "lymphodepleting chemotherapy before CD19-directed CAR T-cell therapy in",
    "adults with relapsed/refractory large B-cell lymphoma",
    "(Varela-Gonzalez-Aller 2025; n = 56 patients, 348 serum samples from",
    "168 infusions). First-order elimination from the central compartment.",
    "Clearance is split into an additive non-renal plus renal pair of arms:",
    "the non-renal arm depends on which CAR T-cell product the patient was",
    "scheduled to receive (4.4 L/h axi-cel vs 3.9 L/h tisa-cel) and the",
    "renal arm scales linearly with creatinine clearance normalized to the",
    "96 mL/min cohort median (1.7 L/h at the median). Allometric body-weight",
    "scaling is fixed at 0.75 on all three clearance arms and 1 on all three",
    "volumes, with a 70 kg reference. Correlated inter-individual",
    "variability on CL and V1 (rho = 0.9); log-normal residual error. No",
    "exposure-response layer was reported, so only the population PK model",
    "is encoded here."
  )
  reference <- paste(
    "Varela-Gonzalez-Aller J, Sanchez-Salinas MA, Troconiz I, Iacoboni G,",
    "Alonso-Martinez C, Carreras-Soler MJ, Carpio C, Farriols-Danes A,",
    "Guerra-Gonzalez M, Rivas-Delgado A, Rivera Sanchez L, Feijoo S,",
    "Valdivia C, Barba P, Miarons M (2025). Towards Personalized",
    "Lymphodepletion: A Population Pharmacokinetic Fludarabine Model in",
    "Patients Receiving CAR T-Cell Therapy. Pharmaceutics 17(12):1592.",
    "doi:10.3390/pharmaceutics17121592.",
    sep = " "
  )
  vignette <- "VarelaGonzalezAller_2025_fludarabine"
  units    <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Dose amounts are F-ara-A equivalents in mg. Fludarabine is administered as
  # the monophosphate prodrug (F-ara-AMP) but the assay measures the
  # circulating dephosphorylated metabolite F-ara-A, so the paper converted
  # every dose before modelling -- Methods 2.1.1: 'As fludarabine is
  # administered as the monophosphate (f-ara-AMP), the corresponding f-ara-A
  # dose equivalents were obtained by multiplying the administered f-ara-AMP
  # doses by 0.78, according to the molecular weight ratio.' The 0.78 factor
  # therefore belongs in the event table, NOT in model(); passing an
  # unconverted F-ara-AMP amount would overstate exposure by 1/0.78 = 1.28x.
  compartmentData <- list(
    central     = list(analyte = "F-ara-A", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "F-ara-A", units = "mg", specimen = "serum", verified = TRUE),
    peripheral2 = list(analyte = "F-ara-A", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Source column WGT. Drives allometric scaling of all three clearance",
        "arms (exponent 3/4) and all three volumes (exponent 1) on a 70 kg",
        "reference, per the Table 2 footnote printed immediately above the",
        "table: 'f(WGT) = (WGT/70)^3/4 ; g(WGT) = WGT/70; WGT, Weight.'",
        "Table 2 reports no exponent estimate and no RSE for either",
        "exponent, and Methods 2.3.2 states that weight effects were",
        "'explored using power relationships with both estimated and fixed",
        "allometric [26] exponents' (reference 26 is Anderson & Holford",
        "2009, the source of the 0.75 / 1 theory-based pair), so both",
        "exponents are encoded as fixed(). Note that although Methods 2.3.2",
        "says continuous covariates were 'normalized by the median value of",
        "the study population', the printed f/g definitions normalize weight",
        "to 70 kg rather than to the 82.5 kg cohort median; the printed",
        "equation governs (and is confirmed by the Discussion back-check",
        "below). Cohort weight 79.6 +/- 12.9 kg, median 82.5, range 52-101",
        "(Table 1)."
      ),
      source_name        = "WGT"
    ),
    CRCL = list(
      description        = "Creatinine clearance, Cockcroft-Gault (Spanish Society of Nephrology calculator)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Source column CRCL. Methods 2.1.2: 'Creatinine clearance (CrCl) was",
        "estimated using the Cockcroft-Gault formula of the Spanish Society",
        "of Nephrology [19] at two time points (pre-lymphodepletion and prior",
        "to the third administration)', i.e. time-varying across the 3-day",
        "lymphodepletion course. The paper labels this quantity 'eGFR' in",
        "Table 1 with units mL/min/1.73 m^2, but Cockcroft-Gault returns a",
        "raw (non-BSA-normalized) mL/min value and Table 2 / the Discussion",
        "both quote it in plain mL/min ('a typical 70 kg patient with a",
        "median eGFR value (96 mL/min)'), so it is recorded here as the raw",
        "mL/min variant of the canonical CRCL column -- the same reading as",
        "the Chen_2023_nemonoxacin.R, Georges_2009_ceftazidime.R and",
        "Delattre_2010_amikacin.R precedents. Enters the renal clearance arm",
        "as a linear proportional term normalized to the 96 mL/min cohort",
        "median: Table 2 'CL_renal = theta_CRCL x CRCL/96' with theta_CRCL =",
        "1.7 L/h. Cohort eGFR mean 90, range 39.0-213.9 (Table 1); median 96",
        "(Discussion). Six patients received protocol dose reductions for",
        "impaired renal function (25% for eGFR 45-60, 40% for eGFR 30-45 per",
        "the internal protocol in Methods 2.1.1) -- those are dose changes in",
        "the event table, not covariate effects."
      ),
      source_name        = "CRCL"
    ),
    CONMED_TISACEL = list(
      description        = "Scheduled CAR T-cell product (1 = tisagenlecleucel, 0 = axicabtagene ciloleucel)",
      units              = "(binary)",
      type               = "binary",
      reference_category = paste(
        "0 (axicabtagene ciloleucel; the structural reference for the Table 2",
        "typical non-renal clearance of 4.4 L/h, corresponding to the paper's",
        "CART = 1 level)"
      ),
      notes              = paste(
        "Time-fixed per subject. The paper codes this as CAR_T with two",
        "levels (Table 2 footnote: 'CART = 1, Axi-cel CAR-T construct;",
        "CART = 2, Tisa-cel CAR-T construct') and selects the non-renal",
        "clearance intercept by level: theta_CART=1 = 4.4 L/h (RSE 1.2%),",
        "theta_CART=2 = 3.9 L/h (RSE 0.9%). Re-expressed here as an axi-cel",
        "referenced binary indicator with a multiplicative power-form factor",
        "e_tisacel_cl_nonren = 3.9/4.4 = 0.8864, which reproduces both",
        "printed thetas exactly. IMPORTANT semantic caveat: the CAR T-cell",
        "product is infused on day 0, whereas fludarabine is given on days -5",
        "to -3 and the entire PK sampling window closes 30 min before the",
        "cell infusion (Methods 2.1.1, 2.1.3). The product has therefore NOT",
        "been administered during any modelled observation, so this is a",
        "planned-treatment / cohort marker rather than a pharmacokinetic",
        "drug-drug interaction. The authors say so themselves (Discussion):",
        "'Although the inclusion of the CAR-T construct led to a",
        "statistically significant reduction in the OFV, the minimal impact",
        "on the typical parameter estimates and the behavior of the",
        "diagnostic plots indicates that this effect is unlikely to be",
        "clinically meaningful and may instead reflect unmeasured covariates",
        "or other patient-specific factors.' The product also determines the",
        "lymphodepletion dose (30 mg/m^2 axi-cel vs 25 mg/m^2 tisa-cel",
        "fludarabine phosphate); that is a dosing difference carried by the",
        "event table, not by this covariate. Cohort 38 axi-cel (68%) / 18",
        "tisa-cel (32%) (Table 1)."
      ),
      source_name        = "CAR_T"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 56L,
    n_studies      = 1L,
    age_range      = "23 to 82 years",
    age_median     = "59 years",
    weight_range   = "52 to 101 kg",
    weight_median  = "82.5 kg (mean 79.6 +/- 12.9)",
    sex_female_pct = 100 * 23 / 56,
    race_ethnicity = NULL,
    disease_state  = paste(
      "Relapsed/refractory large B-cell lymphoma (LBCL) after two or more",
      "lines of systemic therapy: diffuse large B-cell lymphoma 30 (54%),",
      "transformed from indolent lymphoma 20 (36%), high-grade B-cell",
      "lymphoma 3 (5%), other 3 (5%; T-cell/histiocyte-rich large B-cell",
      "lymphoma and intravascular large B-cell lymphoma). ECOG performance",
      "status 0 in 32 (57%), 1 in 20 (36%), 2 in 4 (7%). Two prior lines of",
      "treatment in 45 (80%), three or more in 11 (20%). Prior haematopoietic",
      "stem-cell transplantation in 21 (37%). Body mass index 25.2 +/- 3.5",
      "kg/m^2 (median 25.3, range 18.0-36.2), with obesity (BMI > 30)",
      "documented in two patients. Serum albumin median 4 g/dL (range",
      "2.7-4.8). Table 1."
    ),
    renal_function = paste(
      "Cockcroft-Gault renal function (labelled eGFR in Table 1) mean 90,",
      "range 39.0-213.9 mL/min; cohort median 96 mL/min (Discussion), which",
      "is the normalizing constant in the Table 2 renal clearance arm.",
      "Patients with impaired renal function received protocol-reduced",
      "fludarabine doses (25% reduction for eGFR 45-60 mL/min, 40% for",
      "30-45 mL/min per the hospital's internal protocol); six patients",
      "required a reduction (two axi-cel, four tisa-cel)."
    ),
    dose_range     = paste(
      "Fludarabine phosphate as a 30-minute intravenous infusion on days -5,",
      "-4 and -3 before CAR T-cell infusion, combined with cyclophosphamide:",
      "30 mg/m^2 fludarabine + 500 mg/m^2 cyclophosphamide for axi-cel",
      "recipients (cumulative 90 mg/m^2) and 25 mg/m^2 fludarabine + 250",
      "mg/m^2 cyclophosphamide for tisa-cel recipients (cumulative 75",
      "mg/m^2). Doses enter the model as F-ara-A equivalents, i.e. the",
      "administered F-ara-AMP dose multiplied by 0.78 (molecular-weight",
      "ratio, Methods 2.1.1)."
    ),
    regions        = "Spain (single centre: Hospital Universitari Vall d'Hebron, Barcelona).",
    notes          = paste(
      "Prospective, observational, single-centre study conducted January 2021",
      "to July 2022. Sixty patients enrolled, 56 analysed (four excluded for",
      "not completing the fludarabine regimen per protocol). 38 (68%)",
      "received axicabtagene ciloleucel and 18 (32%) tisagenlecleucel. 168",
      "fludarabine administrations yielded 400 blood samples, of which 348",
      "entered model development (samples below the LLOQ or above the upper",
      "LOQ without further dilution were excluded): 99 at 1.5 h, 97 at 2 h,",
      "94 at 24 h, 51 collected 30 min before CAR T-cell infusion, and only 7",
      "at 7 h post-infusion. Median 6 samples per patient (range 1-9).",
      "Sampling times follow the Monahan limited-sampling model (1, 2, 7 and",
      "24 h) as adapted to hospital practice. Assay: UPLC-MS/MS on serum",
      "after liquid-liquid extraction (Waters Acquity / Xevo TQ, m/z 286.10 >",
      "154.00, internal standard 2-chloroadenosine), calibration range 1-1000",
      "ng/mL with R^2 = 0.999, LLOQ 1 ng/mL. NONMEM 7.4.3 with FOCEI on",
      "log-transformed data; PsN 5.3.0, Pirana 11.1, R 4.2.3. Covariates",
      "selected by stepwise covariate modelling (forward 5%, backward 1%).",
      "Evaluated by prediction-corrected VPC over 1000 simulated datasets.",
      "Baseline demographics from Table 1; final parameter estimates and the",
      "covariate model from Table 2 and the equation block printed",
      "immediately above it. The 7 h sampling gap is an acknowledged",
      "limitation that reduces the precision of V3 and CLD3, which span the",
      "beta/gamma transition (Discussion). No pharmacogenetic data (notably",
      "rs2295890) were collected, and no exposure-response analysis was",
      "performed."
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = paste(
        "Tested but not retained. Discussion: 'Subsequent inclusion of eGFR",
        "in the model eliminated any residual body size-independent effects",
        "of age.' Cohort median 59 years, range 23-82 (Table 1)."
      )
    ),
    SEXF = list(
      description = "Sex indicator (1 = female, 0 = male)",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Collected (Methods 2.1.2) and available to the stepwise covariate",
        "search but not retained in the final model. 33 male (59%) / 23",
        "female (Table 1). The paper reports no covariate-screening table, so",
        "the objective-function change for sex is not recoverable."
      )
    ),
    BSA = list(
      description = "Body surface area",
      units       = "m^2",
      type        = "continuous",
      notes       = paste(
        "Collected as an alternative body-size descriptor (Methods 2.1.2) and",
        "relevant because dosing is per m^2, but body size entered the final",
        "model through allometric weight scaling instead: Discussion 'In our",
        "model, allometric scaling based on WGT best accounted for body size",
        "differences.' Not reported numerically in Table 1."
      )
    ),
    BMI = list(
      description = "Body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = paste(
        "Collected (Methods 2.1.2) and reported in Table 1 (25.2 +/- 3.5,",
        "median 25.3, range 18.0-36.2) but not retained; body size entered",
        "via allometric weight scaling."
      )
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/dL",
      type        = "continuous",
      notes       = paste(
        "Collected on days 1 and 3 of lymphodepletion and before CAR T-cell",
        "infusion (Methods 2.1.2) and reported in Table 1 (median 4, range",
        "2.7-4.8) but not retained in the final model."
      )
    ),
    CREAT = list(
      description = "Serum creatinine",
      units       = "mg/dL",
      type        = "continuous",
      notes       = paste(
        "Collected on days 1 and 3 of lymphodepletion and before CAR T-cell",
        "infusion (Methods 2.1.2). Enters the final model only indirectly, as",
        "the Cockcroft-Gault input to the retained CRCL covariate; it was not",
        "retained as a covariate in its own right."
      )
    ),
    LDH = list(
      description = "Serum lactate dehydrogenase",
      units       = "U/L",
      type        = "continuous",
      notes       = paste(
        "Collected on days 1 and 3 of lymphodepletion and before CAR T-cell",
        "infusion (Methods 2.1.2) but not retained. Not reported numerically",
        "in Table 1."
      )
    ),
    HGB = list(
      description = "Haemoglobin",
      units       = "g/dL",
      type        = "continuous",
      notes       = paste(
        "Collected on days 1 and 3 of lymphodepletion and before CAR T-cell",
        "infusion (Methods 2.1.2) but not retained. Not reported numerically",
        "in Table 1."
      )
    ),
    TX_HSCT = list(
      description = "Prior haematopoietic stem-cell transplantation indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Collected with its date (Methods 2.1.2) and reported in Table 1 (21",
        "patients, 37%) but not retained in the final model."
      )
    )
  )

  ini({
    # Final population estimates -- Varela-Gonzalez-Aller 2025 Table 2
    # ('Final population pharmacokinetic parameter estimates'), with the
    # covariate model as printed in the unnumbered equation block immediately
    # above Table 2 on page 9:
    #
    #   CL           = [CL_non-renal + CL_renal] x f(WGT)
    #   CL_non-renal = theta_CAR-T 1 , axi-cel
    #                  theta_CAR-T 2 , tisa-cel
    #   CL_renal     = theta_CRCL x CRCL/96
    #   f(WGT)       = (WGT/70)^3/4 ;  g(WGT) = WGT/70
    #
    # and the Table 2 'Description' column:
    #
    #   CL   (L/h) = [theta(CAR_T) + theta_CRCL x CRCL/96] x f(WGT)
    #   V1   (L)   = theta_V1   x g(WGT)
    #   V2   (L)   = theta_V2   x g(WGT)
    #   V3   (L)   = theta_V3   x g(WGT)
    #   CLD2 (L/h) = theta_CLD2 x f(WGT)
    #   CLD3 (L/h) = theta_CLD3 x f(WGT)
    #
    # Independent arithmetic check of the whole clearance block against the
    # Discussion, which states: 'For a typical 70 kg patient with a median
    # eGFR value (96 mL/min), the predicted CL was 6.1 L/h for axi-cel and
    # 5.6 L/h for tisa-cel.' At WT = 70 kg, f(WGT) = 1 and CRCL/96 = 1, so
    #   axi-cel : (4.4 + 1.7) x 1 = 6.1 L/h   (matches exactly)
    #   tisa-cel: (3.9 + 1.7) x 1 = 5.6 L/h   (matches exactly)
    # This confirms both the 70 kg allometric reference (a median-82.5 kg
    # reference would give 6.90 and 6.40 L/h before scaling, and 6.1 / 5.6
    # only at 82.5 kg, contradicting the quoted '70 kg patient') and the
    # additive -- rather than multiplicative -- combination of the two
    # clearance arms.

    # Non-renal clearance arm, referenced to the axi-cel cohort. Registered
    # multi-component CL canonical `lcl_nonren` (parameter-names.md
    # "Multi-component clearance"); the paper names this arm CL_non-renal.
    lcl_nonren <- log(4.4);  label("Non-renal clearance arm at WT = 70 kg, axi-cel (L/h)")   # Table 2 row 'CL': theta_CART=1 = 4.4 (RSE 1.2%)

    # Renal clearance arm at the cohort-median CRCL of 96 mL/min. Registered
    # multi-component CL canonical `lcl_renal`; the paper names this arm
    # CL_renal = theta_CRCL x CRCL/96, so theta_CRCL is by construction the
    # renal clearance of a patient at the median renal function.
    lcl_renal  <- log(1.7);  label("Renal clearance arm at CRCL = 96 mL/min and WT = 70 kg (L/h)") # Table 2 row 'CL': theta_CRCL = 1.7 (RSE 1.0%)

    lvc  <- log(41.2); label("Central volume V1 at WT = 70 kg (L)")            # Table 2 row 'V1':   theta_V1   = 41.2 (RSE 9.3%)
    lvp  <- log(14.5); label("Shallow peripheral volume V2 at WT = 70 kg (L)") # Table 2 row 'V2':   theta_V2   = 14.5 (RSE 14.0%)
    lvp2 <- log(10.8); label("Deep peripheral volume V3 at WT = 70 kg (L)")    # Table 2 row 'V3':   theta_V3   = 10.8 (RSE 3.0%)
    lq   <- log(4.8);  label("Intercompartmental clearance CLD2 (central-shallow) at WT = 70 kg (L/h)") # Table 2 row 'CLD2': theta_CLD2 = 4.8 (RSE 5.3%)
    lq2  <- log(3.6);  label("Intercompartmental clearance CLD3 (central-deep) at WT = 70 kg (L/h)")    # Table 2 row 'CLD3': theta_CLD3 = 3.6 (RSE 0.1%)

    # Allometric exponents, fixed a priori rather than estimated: Table 2
    # prints f(WGT) = (WGT/70)^3/4 and g(WGT) = WGT/70 as definitions with no
    # estimate row and no RSE, and Methods 2.3.2 explored 'power relationships
    # with both estimated and fixed allometric [26] exponents' (ref 26 =
    # Anderson & Holford 2009, the theory-based 0.75 / 1 pair).
    e_wt_cl_q  <- fixed(0.75); label("Allometric exponent of body weight on all clearance arms (unitless)") # Table 2 footnote 'f(WGT) = (WGT/70)^3/4'
    e_wt_vc_vp <- fixed(1.0);  label("Allometric exponent of body weight on all volumes (unitless)")        # Table 2 footnote 'g(WGT) = WGT/70'

    # CAR T-cell product effect on the non-renal clearance arm, re-expressed
    # from the paper's two-level categorical intercept onto an axi-cel
    # reference: 3.9 / 4.4 = 0.886364, so 4.4 x 0.886364 = 3.9 exactly.
    e_tisacel_cl_nonren <- 3.9 / 4.4; label("Tisagenlecleucel multiplicative factor on the non-renal clearance arm (power-form base)") # Table 2 row 'CL': theta_CART=2 = 3.9 (RSE 0.9%) relative to theta_CART=1 = 4.4

    # Inter-individual variability. Table 2 reports IIV as %CV and the note
    # printed above Table 2 gives the transform explicitly: 'Interindividual
    # variability (IIV) is expressed as coefficient of variation (CV, %)
    # calculated as sqrt(e^(omega^2) - 1) x 100, where omega^2 corresponds to
    # the variance of the random effects. The estimation of the correlation
    # coefficient between random effects was 0.9.' Inverting,
    # omega^2 = log(CV^2 + 1):
    #   CL: log(0.298^2 + 1) = 0.085080 -> 0.0851  (round-trips to 29.80% CV)
    #   V1: log(0.348^2 + 1) = 0.114314 -> 0.1143  (round-trips to 34.80% CV)
    # Off-diagonal from rho = 0.9: 0.9 x sqrt(0.0851 x 0.1143) = 0.088763
    # -> 0.0888 (round-trips to rho = 0.9004). Results 3.3 confirms the block
    # structure: 'Interindividual variability for CL and V1 was modeled using
    # a block covariance matrix ... The estimated covariance revealed a
    # moderate positive association between CL and V1 across the population
    # (correlation coefficient rho ~ 0.9).' The IIV multiplies total CL (i.e.
    # both arms combined), matching the Table 2 'Description' column, which
    # places the covariate expression inside a single CL row.
    #
    # No IIV on CLD2, CLD3, V2 or V3 -- Results 3.3: 'IIV terms for
    # intercompartmental clearances (CLD2 and CLD3) and peripheral volumes of
    # distribution (V2 and V3) were not included in the model.'
    etalcl + etalvc ~ c(0.0851,
                        0.0888, 0.1143)  # Table 2 'IIV on CL (%)' = 29.8 (RSE 25.1, shrinkage 8.9) and 'IIV on V1 (%)' = 34.8 (RSE 3.6, shrinkage 18.3); rho = 0.9 per the note above Table 2

    # Residual error. Methods 2.3: 'Data were logarithmically transformed for
    # analysis' and 2.3.1: 'The residual model structure was selected among
    # the additive, proportional, and combined error models in logarithmic
    # scale.' Table 2 reports the selected term as 'Error [log(ng/mL)] = 0.29
    # (RSE 18.6%, shrinkage 12.1%)' -- an additive SD on the natural-log scale,
    # i.e. O_ij = Y_ij x exp(eps_ij), which is nlmixr2's lnorm() error and NOT
    # add() on the linear scale. Encoded as expSd per parameter-names.md
    # (Billard_1995_methoxsalen_plasma.R precedent). Because the error is
    # additive on the log scale it is invariant to the mg/L-to-ng/mL scale
    # factor applied to Cc below.
    expSd <- 0.29; label("Log-normal residual error (log-scale SD)")  # Table 2 row 'Error [log(ng/mL)]' = 0.29 (RSE 18.6%)
  })

  model({
    # Clearance is the sum of a non-renal and a renal arm, each carrying the
    # same fixed allometric weight exponent, exactly as printed above Table 2:
    # CL = [CL_non-renal + CL_renal] x f(WGT). The correlated log-normal BSV
    # multiplies the combined total, matching the single CL row in Table 2.
    cl_nonren <- exp(lcl_nonren) * e_tisacel_cl_nonren^CONMED_TISACEL
    cl_renal  <- exp(lcl_renal)  * (CRCL / 96)
    cl        <- (cl_nonren + cl_renal) * (WT / 70)^e_wt_cl_q * exp(etalcl)

    # Distribution clearances and volumes. Volumes scale with g(WGT) = WGT/70;
    # distribution clearances with f(WGT) = (WGT/70)^(3/4).
    q   <- exp(lq)   * (WT / 70)^e_wt_cl_q
    q2  <- exp(lq2)  * (WT / 70)^e_wt_cl_q
    vc  <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc_vp
    vp  <- exp(lvp)  * (WT / 70)^e_wt_vc_vp
    vp2 <- exp(lvp2) * (WT / 70)^e_wt_vc_vp

    # Micro-rate constants for the three-compartment mammillary form.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # Three-compartment mammillary model with first-order elimination from the
    # central compartment (Results 3.3: 'A three-compartment model best
    # described the concentration data of fludarabine'; Discussion:
    # 'Fludarabine PK was best characterized using a three-compartment model
    # with first-order elimination'). All input is a 30-minute intravenous
    # infusion into the central compartment (Methods 2.1.1), so there is no
    # depot state and no bioavailability term.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1 - k13 * central + k31 * peripheral2
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1
    d/dt(peripheral2) <-                                                       k13 * central - k31 * peripheral2

    # central is in mg and vc in L, so central/vc is mg/L; multiply by 1000 to
    # obtain ng/mL (1 mg/L = 1000 ng/mL), the unit the paper reports throughout
    # (Table 2 'Error [log(ng/mL)]', the 1-1000 ng/mL calibration range in
    # Methods 2.2, and Figure 1).
    Cc <- 1000 * central / vc

    Cc ~ lnorm(expSd)
  })
}
