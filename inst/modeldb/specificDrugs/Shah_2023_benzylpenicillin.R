Shah_2023_benzylpenicillin <- function() {
  description <- paste(
    "Two-compartment intravenous population PK model for benzylpenicillin",
    "in critically unwell adults, with a priori allometric body-weight",
    "scaling on all disposition parameters and a serum-creatinine power",
    "effect on clearance."
  )
  reference <- paste(
    "Shah RV, Kipper K, Baker EH, Barker CIS, Oldfield I, Philips BJ,",
    "Johnston A, Lipman J, Rhodes A, Basarab M, Sharland M, Almahdi S,",
    "Wake RM, Standing JF, Lonsdale DO (2023).",
    "Population Pharmacokinetic Study of Benzylpenicillin in Critically",
    "Unwell Adults. Antibiotics 12(4):643. doi:10.3390/antibiotics12040643.",
    sep = " "
  )
  vignette <- "Shah_2023_benzylpenicillin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Derived mechanically; verified = FALSE means it has
  # NOT been checked against the source paper.
  compartmentData <- list(
    central     = list(analyte = "benzylpenicillin", units = "mg", specimen = "plasma", verified = FALSE),
    peripheral1 = list(analyte = "benzylpenicillin", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Added a priori to all primary PK parameters by allometric scaling,",
        "normalised to a reference weight of 70 kg (Shah 2023 Section 4,",
        "Materials and Methods: 'Weight was added to primary pharmacokinetic",
        "parameters a priori using allometric scaling. Compartment volumes were",
        "scaled with a fixed exponent of 1, whereas clearance parameters were",
        "scaled to an allometric exponent of 0.75'). Table 2 reports every",
        "structural theta per 70 kg. Cohort median 70.0 kg (IQR 65.7-90.0,",
        "range 60.0-120.0; Table 1)."
      ),
      source_name        = "Weight (kg)"
    ),
    CREAT = list(
      description        = "Serum creatinine",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power effect on clearance normalised to 70 umol/L (Shah 2023",
        "Equation (1)). Table 1 prints the units as 'mmol/L', but the reported",
        "median of 70 (IQR 52-103.5, range 34-486) is only physiologically",
        "interpretable as umol/L -- 70 mmol/L is roughly a thousand-fold above",
        "any survivable serum creatinine. The 70 in Equation (1) is the cohort",
        "median from Table 1, so the normalisation constant and the covariate",
        "column must share the same unit; umol/L is used here. Time-varying in",
        "the source dataset (routine ICU biochemistry)."
      ),
      source_name        = "Serum creatinine"
    )
  )

  # Covariates screened in the stepwise covariate search (Shah 2023 Section 4
  # and Supplementary Table S1, runs 5-19) but NOT retained in the final model.
  # Documented here for provenance only; none is referenced in model(). These
  # names are documentation labels, not ratified canonical covariate columns.
  covariatesDataExcluded <- list(
    HT = list(
      description = "Height",
      units       = "cm",
      type        = "continuous",
      notes       = paste(
        "Tested on V1 (Table S1 run 6, dOFV 0), CL (run 7, dOFV -4.0) and V2",
        "(run 8, dOFV -0.2). The CL effect met the dOFV threshold but was",
        "rejected for a 5119% RSE on the covariate coefficient. Cohort median",
        "172.0 cm (IQR 163.5-178.0; Table 1)."
      )
    ),
    BMI = list(
      description = "Body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = paste(
        "Tested on V1 (Table S1 run 5, dOFV -0.8); rejected. Cohort median",
        "26.1 kg/m^2 (IQR 22.1-27.9; Table 1)."
      )
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "categorical",
      notes       = paste(
        "Tested on CL (Table S1 runs 17-18, dOFV -0.3 and 0) and V1 (run 19,",
        "dOFV 0); rejected. Cohort was 6 male : 6 female (Table 1). The source",
        "paper reports the covariate as 'SEX' without stating the reference",
        "level; because the effect was rejected, no coefficient or reference",
        "category needs to be resolved."
      )
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/L",
      type        = "continuous",
      notes       = paste(
        "Tested on V1 (Table S1 run 9, dOFV -0.6), Q (run 10, dOFV -1.6) and",
        "V2 (run 11, dOFV -2.3); all rejected. Cohort median 28 g/L (IQR",
        "24-34; Table 1)."
      )
    ),
    BODYTEMP = list(
      description = "Body temperature",
      units       = "degC",
      type        = "continuous",
      notes       = paste(
        "Tested on V1 (Table S1 run 12, dOFV -2.9), V2 (run 13, dOFV -10.2)",
        "and Q (run 14, dOFV -1.9). The V2 effect met the dOFV threshold but",
        "was rejected for a large increase in the uncertainty of the other",
        "parameters. Table 1 does not report the cohort temperature",
        "distribution."
      )
    ),
    APACHE_II = list(
      description = "Acute Physiology and Chronic Health Evaluation II score",
      units       = "points",
      type        = "continuous",
      notes       = paste(
        "Tested on V1 (Table S1 run 15, dOFV -0.8) and V2 (run 16, dOFV -0.2);",
        "rejected. Cohort median 14 points (IQR 12.5-18, range 5-23; Table 1).",
        "Documentation-only label -- APACHE II is not a ratified canonical",
        "covariate column in inst/references/covariate-columns.md, and it is",
        "deliberately not registered here because the final model does not use",
        "it. It is NOT the same instrument as the registered SAPS_II canonical."
      )
    ),
    CRCL = list(
      description = "Creatinine clearance",
      units       = "mL/min",
      type        = "continuous",
      notes       = paste(
        "Deliberately NOT tested by the authors. Shah 2023 Section 3:",
        "'unlike the two studies noted, we have chosen not to test CrCl as a",
        "covariate effect on clearance, as this was not directly measured",
        "during the study and estimates for CrCl are not validated in",
        "critically ill populations.' Recorded here so that a reader comparing",
        "this model with the CrCl-on-CL models of Bos et al. and",
        "Obrink-Hansen et al. sees why raw serum creatinine was used instead."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 12,
    n_studies      = 1,
    n_samples      = 77,
    age_range      = "25.7-71.7 years",
    age_median     = "57.7 years",
    weight_range   = "60.0-120.0 kg",
    weight_median  = "70.0 kg",
    height_range   = "150.0-188.0 cm",
    sex_female_pct = 50,
    race_ethnicity = c(
      `White British` = 50.0,
      `White Irish`   = 8.3,
      Asian           = 16.7,
      Caribbean       = 8.3,
      `Other/Not stated` = 16.7
    ),
    disease_state  = paste(
      "Critical illness requiring intensive care; treated infection sources",
      "were lower respiratory tract infection (7), skin/soft-tissue infection",
      "or abscess (5), infective endocarditis (1) and sepsis of unknown",
      "source (1)"
    ),
    dose_range     = paste(
      "1.2 g IV 4-hourly (8 patients), 1.2 g IV 6-hourly (2 patients),",
      "2.4 g IV 4-hourly (2 patients)"
    ),
    regions        = "United Kingdom (single centre: St George's Hospital, London)",
    renal_function = paste(
      "Serum creatinine median 70 umol/L (IQR 52-103.5, range 34-486);",
      "three participants had acute kidney injury. One patient received renal",
      "replacement therapy; only samples drawn after renal recovery and",
      "cessation of RRT were analysed"
    ),
    severity       = paste(
      "APACHE II median 14 points (IQR 12.5-18, range 5-23); 5 patients on",
      "vasopressors; 1 intubated and ventilated, 3 non-invasive ventilation,",
      "10 spontaneous ventilation (ventilation categories overlap over time)"
    ),
    notes          = paste(
      "Sub-study of the ABDose observational antibiotic PK/PD study",
      "(REC 14/LO/1999). Baseline demographics: Shah 2023 Table 1. Sampling",
      "was opportunistic within a dosing interval (Table 3). Twelve patients",
      "contributed 80 plasma samples, of which 77 entered the analysis.",
      "Model fit in NONMEM 7.5 using FOCE-i."
    )
  )

  ini({
    # Structural parameters: Shah 2023 Table 2, 'Fixed effects', reported as
    # typical values for a 70 kg adult with a serum creatinine of 70 umol/L.
    lcl <- log(23.1); label("Clearance (L/h/70 kg)")                      # Table 2, theta_CL = 23.1 L/h/70 kg (14% RSE)
    lvc <- log(15.1); label("Central volume of distribution (L/70 kg)")   # Table 2, theta_V1 = 15.1 L/70 kg (8% RSE)
    lq  <- log(11.1); label("Intercompartmental clearance (L/h/70 kg)")   # Table 2, theta_Q = 11.1 L/h/70 kg (50% RSE)
    lvp <- log(9.8);  label("Peripheral volume of distribution (L/70 kg)")# Table 2, theta_V2 = 9.8 L/70 kg (29% RSE)

    # Allometric exponents: fixed a priori, not estimated. Shah 2023 Section 4:
    # "Compartment volumes were scaled with a fixed exponent of 1, whereas
    # clearance parameters were scaled to an allometric exponent of 0.75".
    # A single exponent is shared across the two clearance parameters (CL, Q)
    # and a single exponent across the two volumes (V1, V2), hence the paired
    # e_<cov>_<param1>_<param2> canonical form.
    e_wt_cl_q  <- fixed(0.75); label("Allometric exponent on CL and Q (unitless)")   # Section 4, Materials and Methods
    e_wt_vc_vp <- fixed(1);    label("Allometric exponent on Vc and Vp (unitless)")  # Section 4, Materials and Methods

    # Covariate effect: power function of serum creatinine on clearance,
    # normalised to 70 umol/L. Shah 2023 Equation (1):
    #   CL = theta_CL * exp(eta_1) * (Creatinine / 70)^theta_creat
    e_creat_cl <- -0.916; label("Power exponent for serum creatinine on CL (unitless)")  # Table 2, theta_CREAT = -0.916 (18% RSE); Equation (1)

    # IIV: Shah 2023 Table 2, 'Random effects', reported as %CV. Converted to
    # the log-normal variance scale with omega^2 = log(CV^2 + 1).
    # No IIV on Q -- Section 2: "Parameters representing interindividual
    # variability were added to all parameters other than Q, for which this
    # parameter was found to be negligible."
    etalcl ~ log(0.420^2 + 1)  # Table 2, omega^2_1 CL = 42.0 %CV (43% RSE) -> 0.1624589
    etalvc ~ log(0.226^2 + 1)  # Table 2, omega^2_2 V1 = 22.6 %CV (42% RSE) -> 0.0498144
    etalvp ~ log(0.205^2 + 1)  # Table 2, omega^2_3 V2 = 20.5 %CV (60% RSE) -> 0.0411659

    # Residual error: Shah 2023 Table 2, 'Residual error', reported as sigma^2
    # variances; nlmixr2 parameterises the combined error by standard
    # deviations, so each is entered as sqrt(sigma^2).
    propSd <- sqrt(0.021); label("Proportional residual error (fraction)")  # Table 2, sigma^2_1 (proportional) = 0.021 (47% RSE) -> SD 0.1449
    addSd  <- sqrt(0.006); label("Additive residual error (mg/L)")          # Table 2, sigma^2_2 (additive)     = 0.006 (67% RSE) -> SD 0.0775
  })

  model({
    # Individual parameters. Reference weight 70 kg (Table 2 thetas are all
    # reported "/70 kg"); reference serum creatinine 70 umol/L (Equation (1)).
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl_q * (CREAT / 70)^e_creat_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc_vp
    q  <- exp(lq)           * (WT / 70)^e_wt_cl_q
    vp <- exp(lvp + etalvp) * (WT / 70)^e_wt_vc_vp

    # Micro-constants
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Two-compartment intravenous disposition (Shah 2023 Section 2: a
    # two-compartment structural model reduced the OFV by 104.4 relative to a
    # one-compartment model; Supplementary Table S1 runs 1-2).
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
