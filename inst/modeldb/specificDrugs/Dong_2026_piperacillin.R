Dong_2026_piperacillin <- function() {
  description <- paste(
    "One-compartment population PK model for intravenous piperacillin in",
    "critically ill adults (Dong 2026; n = 42 Chinese ICU patients, 117",
    "steady-state plasma concentrations spanning 1.25-376.34 mg/L). Linear",
    "first-order elimination from a single central compartment. Clearance",
    "carries two covariates -- cystatin-C-based CKD-EPI estimated glomerular",
    "filtration rate centred on 46.56 mL/min/1.73 m^2 (exponent 0.615) and",
    "total body weight centred on 70 kg (exponent 1.13) -- and central volume",
    "carries serum albumin centred on 34.8 g/L (exponent 1.21). Exponential",
    "inter-individual variability on CL and V; combined proportional-plus-",
    "additive residual error. Patients on continuous renal replacement therapy",
    "were excluded, so the model carries no information about extracorporeal",
    "clearance."
  )
  reference <- paste(
    "Dong Z, Shi H, Yang Y, Yi Q, Jiang Z, Li Y (2026).",
    "Population pharmacokinetics and dosing regimen optimization of",
    "piperacillin in critically ill patients.",
    "Drug Des Devel Ther 20.",
    "doi:10.2147/DDDT.S551307.",
    sep = " "
  )
  vignette <- "Dong_2026_piperacillin"
  units    <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    # Methods "Sample Collection": blood was centrifuged "to separate the
    # plasma" and "the plasma concentration of piperacillin was determined by
    # high-performance liquid chromatography". The assayed analyte is
    # piperacillin itself; tazobactam was not measured or modelled.
    central = list(analyte = "piperacillin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description        = paste(
        "Estimated glomerular filtration rate from the 2012 CKD-EPI",
        "cystatin-C equation, body-surface-area normalised to",
        "1.73 m^2."
      ),
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject. Power effect on CL, centred on 46.56",
        "mL/min/1.73 m^2, per the final-model equation printed on Dong 2026",
        "p. 6: CL (L/h) = 6.48 * e^eta * (eGFR/46.56)^0.615 *",
        "(total body weight/70)^1.13. The exponent 0.615 is the 'eGFR on CL'",
        "row of Table 3 (RSE 13.5%; bootstrap median 0.621, 5th-95th",
        "0.482-0.76), and the centring constant 46.56 equals the Table 1",
        "cohort median for '2012 CKD-EPI CYS-C' exactly (IQR 26.88-68.26).",
        "The choice of estimating equation is load-bearing rather than",
        "incidental: the authors screened Cockcroft-Gault, MDRD, CKD-EPIcr,",
        "CKD-EPIcys and CKD-EPIcr-cys as continuous covariates on CL and",
        "retained CKD-EPIcys because it produced the largest objective",
        "function drop, reasoning that in critically ill patients creatinine",
        "is confounded by muscle mass and systemic inflammation while",
        "cystatin C is not (Discussion). Note this cohort is renally",
        "IMPAIRED on the cystatin-C scale -- median 46.56, IQR 26.88-68.26 --",
        "which is far below the creatinine-based medians the same patients",
        "produce in Table 1 (2021 CKD-EPIcr 93.32; abbreviated MDRD 107.30),",
        "so this column is NOT interchangeable with a creatinine-based eGFR",
        "for the same subject. Patients on continuous renal replacement",
        "therapy were excluded. The paper's own target-attainment",
        "simulations extrapolate the term well above the fitted range, to",
        "eGFR levels at and beyond the 130 mL/min/1.73 m^2 augmented-renal-",
        "clearance threshold (Figure 6)."
      ),
      source_name        = "eGFR"
    ),
    WT = list(
      description        = "Total body weight.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject. Power effect on CL, centred on 70 kg, per",
        "the p. 6 final-model equation. The exponent 1.13 is the 'Total body",
        "weight on CL' row of Table 3 (RSE 13.1%; bootstrap median 1.12,",
        "5th-95th 0.73-1.4), and 70 kg equals the Table 1 cohort median",
        "exactly (IQR 65-75). Two points are worth carrying forward. First,",
        "this is an ESTIMATED exponent that came out near 1 rather than a",
        "fixed allometric 0.75; the bootstrap interval 0.73-1.4 covers both,",
        "but 1.13 is the value the typical CL of 6.48 L/h is conditioned on",
        "and is what is encoded. Second, weight enters CL only -- it does",
        "NOT scale V, which instead carries albumin. The authors note the",
        "weight-CL mechanism 'remains unclear' and speculate about increased",
        "renal blood flow in obesity and free-fatty-acid displacement of",
        "piperacillin from albumin (Discussion). The cohort IQR is narrow",
        "(65-75 kg) while the paper's simulations extrapolate to over 110 kg",
        "(Figure 7), so the term is doing substantial work outside the range",
        "it was fitted in."
      ),
      source_name        = "Total body weight"
    ),
    ALB = list(
      description        = "Serum albumin.",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject. Power effect on V, centred on 34.8 g/L, per",
        "the p. 6 final-model equation V (L) = 19 * e^eta * (albumin/34.8)^1.21.",
        "The exponent 1.21 is the 'Albumin on V' row of Table 3 and is the",
        "least precisely estimated parameter in the model (RSE 33.4%;",
        "bootstrap median 1.23, 5th-95th 0.53-2.06). The centring constant",
        "34.8 equals the Table 1 cohort median exactly (IQR 32.7-39.2).",
        "The SIGN is worth stating because it is the opposite of the usual",
        "expectation for a protein-bound drug: the exponent is POSITIVE, so",
        "LOWER albumin gives a SMALLER V, whereas the paper's own Discussion",
        "argues hypoalbuminaemia should raise V ('hypoalbuminemia reduces",
        "drug-albumin binding, leading to higher free drug concentrations",
        "that may distribute more readily into interstitial tissues').",
        "The fitted sign is what is encoded -- it is the value the model was",
        "estimated with and it is consistent across the paper's equation and",
        "table -- but a user should be aware the paper does not reconcile",
        "the discrepancy, and the 5th-95th bootstrap interval is wide enough",
        "that the effect is weakly identified. Note also that the assay",
        "measured TOTAL, not free, piperacillin (Limitations), so an albumin",
        "effect on the total-drug volume is not directly a protein-binding",
        "measurement."
      ),
      source_name        = "Albumin"
    )
  )

  # Screened during forward inclusion / backward elimination but not retained
  # in the final model (Methods, "Population Pharmacokinetic Model
  # Development"). The paper reports the screening list but no coefficient,
  # dOFV or point estimate for any of these, so nothing is encoded. The
  # non-retained eGFR equations (Cockcroft-Gault, abbreviated MDRD, modified
  # MDRD, MDRD CHN, 2021 CKD-EPIcr, 2021 CKD-EPIcr-cys) are not listed
  # separately here because they are alternative estimators of the SAME
  # canonical column that the model does retain (CRCL); the retained
  # estimating equation is documented in that entry's notes.
  covariatesDataExcluded <- list(
    SEXF = list(
      description        = "Female sex indicator.",
      units              = "unitless",
      type               = "categorical",
      reference_category = "male",
      notes              = "Screened as a covariate but not retained. Cohort was 30 male / 12 female (Table 1).",
      source_name        = "Sex"
    ),
    AGE = list(
      description        = "Age at enrolment.",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened but not retained. Cohort median 59 years (IQR 50.25-75.5), Table 1.",
      source_name        = "Age"
    ),
    BMI = list(
      description        = "Body mass index.",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened but not retained; total body weight was retained on CL",
        "instead. Cohort median 24.22 kg/m^2 (IQR 22.41-26.12), Table 1."
      ),
      source_name        = "BMI"
    ),
    APACHE_II = list(
      description        = "Acute Physiology and Chronic Health Evaluation II score.",
      units              = "unitless",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened but not retained. Cohort median 21.5 (IQR 19-26), Table 1.",
      source_name        = "APACHE II score"
    ),
    TBILI = list(
      description        = "Total bilirubin.",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened but not retained. Cohort median 10.1 umol/L (IQR 6.6-15.7), Table 1.",
      source_name        = "Total bilirubin"
    ),
    DBIL = list(
      description        = "Direct (conjugated) bilirubin.",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened but not retained. Cohort median 5.6 umol/L (IQR 4.1-10.6), Table 1.",
      source_name        = "Direct bilirubin"
    ),
    CREAT = list(
      description        = "Serum creatinine.",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened but not retained, both directly and by way of the",
        "creatinine-based eGFR equations. Cohort median 69 umol/L",
        "(IQR 48-130), Table 1."
      ),
      source_name        = "Serum creatinine"
    ),
    CYSC = list(
      description        = "Serum cystatin C.",
      units              = "mg/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Screened as a raw laboratory value but not retained; the model",
        "instead carries the cystatin-C-based CKD-EPI eGFR derived from it",
        "(see the CRCL entry). Cohort median 1.44 mg/L (IQR 1.02-2.19),",
        "Table 1."
      ),
      source_name        = "Cystatin C"
    ),
    BUN = list(
      description        = "Blood urea nitrogen.",
      units              = "mmol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened but not retained. Cohort median 8.9 mmol/L (IQR 6.2-14.6), Table 1.",
      source_name        = "Blood urea nitrogen"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 42L,
    n_studies      = 1L,
    n_observations = 117L,
    age_range      = "not reported; interquartile range 50.25-75.5 years",
    age_median     = "59 years",
    weight_range   = "not reported; interquartile range 65-75 kg",
    weight_median  = "70 kg",
    sex_female_pct = 100 * 12 / 42,
    race_ethnicity = c(Asian = 100),
    disease_state  = paste(
      "Adults admitted to the intensive care unit and treated with",
      "piperacillin. Pulmonary infection predominated (92.86%), followed by",
      "sepsis (14.29%), abdominal cavity infection (7.14%), skin and soft",
      "tissue infection (4.76%) and urinary, intestinal and multiple",
      "infections; 85.71% were on mechanical ventilation and 85.71% on a",
      "vasopressor (Table 2). Pseudomonas aeruginosa was the most common",
      "isolate (26.19%), then Klebsiella pneumoniae (14.29%), Escherichia",
      "coli and Stenotrophomonas maltophilia (11.90% each) and Acinetobacter",
      "baumannii (9.52%); no pathogen was identified in 16.67%. Median",
      "APACHE II score 21.5 (IQR 19-26). Patients on continuous renal",
      "replacement therapy, enrolled in another interventional trial, or",
      "pregnant were excluded."
    ),
    dose_range     = paste(
      "Intravenous piperacillin 4 g q8h by infusion in 39 of 42 patients",
      "(92.86%, 12 g/day), 4 g q8h by intravenous push in 1 (2.38%), and",
      "4 g q12h by infusion in 2 (4.76%, 8 g/day) -- Table 2. The clinical",
      "infusion duration is not stated; the paper's Monte Carlo simulations",
      "standardise intermittent infusions at 30 minutes and additionally",
      "explore continuous infusion and daily doses of 8, 12, 16, 20 and 24 g",
      "at q4h, q6h, q8h and q12h."
    ),
    renal_function = paste(
      "Markedly discordant between markers. Median 2012 CKD-EPI cystatin-C",
      "eGFR 46.56 mL/min/1.73 m^2 (IQR 26.88-68.26) -- the equation the",
      "final model uses -- against median 2021 CKD-EPI creatinine eGFR 93.32",
      "(IQR 54.64-108.68), abbreviated MDRD 107.30, modified MDRD 100.95,",
      "MDRD CHN 112.63, 2021 CKD-EPI creatinine-cystatin C 68.07 and",
      "Cockcroft-Gault creatinine clearance 85.46 mL/min (Table 1). Median",
      "serum creatinine 69 umol/L (IQR 48-130), cystatin C 1.44 mg/L (IQR",
      "1.02-2.19). Patients on continuous renal replacement therapy were",
      "excluded."
    ),
    sampling       = paste(
      "Sparse steady-state sampling after at least five piperacillin doses,",
      "2-3 samples per patient: a trough 30 minutes before the next",
      "infusion, a peak at the end of the infusion, and an optional",
      "mid-interval sample within the same dosing interval. Assayed by",
      "HPLC-UV at 218 nm over a linear range of 1.25-400 mg/L; observed",
      "concentrations spanned 1.25-376.34 mg/L. Samples below the limit of",
      "quantification (2.56% of the post-first-dose data) were retained and",
      "set to half the LLOQ."
    ),
    regions        = paste(
      "People's Republic of China (single centre; intensive care unit of the",
      "First Affiliated Hospital of Shandong First Medical University /",
      "Shandong Provincial Qianfoshan Hospital, Jinan)."
    ),
    notes          = paste(
      "Baseline demographics from Dong 2026 Tables 1 and 2. Prospective",
      "observational single-centre study running September 2021 to January",
      "2022 (ethics approval YXLL-KY-2021(050)). The final model was fit in",
      "NONMEM 7.4 via Pirana 2.9.7 with R 3.6.0, and evaluated by",
      "goodness-of-fit plots, a 1000-replicate bootstrap (all successful;",
      "every final estimate inside the 5th-95th bootstrap interval),",
      "normalised prediction distribution errors (mean 0.0161, variance",
      "0.954; Wilcoxon p = 0.93, Fisher variance p = 0.757, Shapiro-Wilk",
      "p = 0.592) and a prediction-corrected visual predictive check. The",
      "assay measured TOTAL plasma piperacillin; the unbound fraction of",
      "0.70 used in the paper's target-attainment simulations is a",
      "literature assumption, not a fitted parameter, and is therefore not",
      "encoded here. Tazobactam was not measured or modelled."
    )
  )

  ini({
    # Structural parameters. Typical values refer to the covariate reference
    # subject -- eGFR 46.56 mL/min/1.73 m^2, total body weight 70 kg and
    # serum albumin 34.8 g/L -- per the final-model equations printed on
    # Dong 2026 p. 6.

    lcl <- log(6.48); label("Clearance at eGFR 46.56 mL/min/1.73 m^2 and 70 kg (L/h)")  # Table 3, CL 6.48 (RSE 5.8%; bootstrap median 6.46, 5th-95th 5.82-7.1); also the leading coefficient of the p. 6 CL equation
    lvc <- log(19);   label("Central volume of distribution at albumin 34.8 g/L (L)")   # Table 3, V 19 (RSE 6.8%; bootstrap median 18.9, 5th-95th 17-21.4); also the leading coefficient of the p. 6 V equation

    # Covariate effects, all power functions of the covariate divided by its
    # cohort-median centring constant (p. 6 equations).
    e_crcl_cl <- 0.615; label("Power exponent of CKD-EPIcys eGFR on CL (unitless)")  # Table 3, 'eGFR on CL' 0.615 (RSE 13.5%; bootstrap median 0.621, 5th-95th 0.482-0.76)
    e_wt_cl   <- 1.13;  label("Power exponent of total body weight on CL (unitless)")  # Table 3, 'Total body weight on CL' 1.13 (RSE 13.1%; bootstrap median 1.12, 5th-95th 0.73-1.4)
    e_alb_vc  <- 1.21;  label("Power exponent of serum albumin on V (unitless)")  # Table 3, 'Albumin on V' 1.21 (RSE 33.4%; bootstrap median 1.23, 5th-95th 0.53-2.06)

    # Inter-individual variability, exponential on CL and V (Results: 'Inter-
    # individual variability of CL and V can be well fitted by an exponential
    # model').
    #
    # SCALE. Table 3 reports these rows as 'IIV-CL (%)' 33.8 and 'IIV-V (%)'
    # 26.3, a percentage column that on its own is ambiguous between a
    # log-scale SD and a lognormal CV. The paper's own printed equations
    # settle it: the p. 6 equations carry the literal factors e^0.114 on CL
    # and e^0.069 on V, which are the eta terms with the fitted variance
    # substituted for the random draw. Those constants match the SD reading
    # and only the SD reading --
    #   0.338^2 = 0.114244 -> 0.114   and   0.263^2 = 0.069169 -> 0.069,
    # whereas the CV reading would require omega^2 = log(1 + 0.338^2) = 0.108
    # and log(1 + 0.263^2) = 0.067. Read the other way round, an omega^2 of
    # 0.114 implies a CV of sqrt(exp(0.114) - 1) = 34.8%, not the tabulated
    # 33.8%. So the tabulated percentages are 100 * omega, and the variances
    # encoded below are the printed equation constants themselves.
    etalcl ~ 0.114  # p. 6 CL equation factor 'e^0.114'; = (0.338)^2 from Table 3 'IIV-CL (%)' 33.8 (RSE 14.6%; bootstrap median 32.7, 5th-95th 24.4-40.9)
    etalvc ~ 0.069  # p. 6 V equation factor 'e^0.069'; = (0.263)^2 from Table 3 'IIV-V (%)' 26.3 (RSE 24.8%; bootstrap median 24.2, 5th-95th 0.3-34.4)

    # Residual error: combined proportional plus additive (Results: 'residual
    # variability were best described by a combined model'). The Table 3
    # footnote names the two rows explicitly -- 'RSV_CV, proportional residual
    # variation; RSV_SD, additive type residual variation' -- so the 17.7%
    # row is a proportional SD and the 1.2 mg/L row an additive SD on the
    # mg/L concentration scale.
    propSd <- 0.177; label("Proportional residual error (fraction)")  # Table 3, RSV_CV 17.7% (RSE 15.9%; bootstrap median 17.7, 5th-95th 13-22.9)
    addSd  <- 1.2;   label("Additive residual error (mg/L)")          # Table 3, RSV_SD 1.2 mg/L (RSE 23.3%; bootstrap median 1.162, 5th-95th 0.528-1.723)
  })

  model({
    # Individual PK parameters. Both covariates on CL and the single
    # covariate on V enter as power functions centred on the cohort medians
    # (Dong 2026 p. 6 final-model equations).
    cl <- exp(lcl + etalcl) * (CRCL / 46.56)^e_crcl_cl * (WT / 70)^e_wt_cl
    vc <- exp(lvc + etalvc) * (ALB / 34.8)^e_alb_vc

    # Piperacillin is given intravenously; dose records target `central`
    # directly (as an infusion), so there is no depot compartment.
    d/dt(central) <- -(cl / vc) * central

    # Total plasma piperacillin concentration in mg/L.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
