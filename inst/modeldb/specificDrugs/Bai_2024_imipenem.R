Bai_2024_imipenem <- function() {
  description <- paste(
    "Two-compartment intravenous population PK model for imipenem in",
    "critically ill adults with sepsis (Bai 2024; n = 51 Chinese ICU",
    "patients, 196 plasma samples). Zero-order infusion into a central",
    "compartment with first-order elimination and linear distribution to a",
    "single peripheral compartment. Cockcroft-Gault creatinine clearance",
    "enters clearance as a power function centred on 99.896 mL/min",
    "(exponent 0.473); it was the only covariate retained after forward",
    "inclusion and backward elimination. Inter-individual variability is",
    "exponential on all four disposition parameters and residual error is",
    "proportional. Imipenem was given as the fixed-ratio imipenem-cilastatin",
    "combination, but only imipenem was assayed and modelled."
  )
  reference <- paste(
    "Bai J, Wen A, Li Z, Li X, Duan M (2024).",
    "Population pharmacokinetics and dosing optimisation of imipenem in",
    "critically ill patients.",
    "Eur J Hosp Pharm 31(5):434-439.",
    "doi:10.1136/ejhpharm-2022-003403. PMCID PMC11347199.",
    sep = " "
  )
  vignette <- "Bai_2024_imipenem"
  units    <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Methods "Laboratory analysis": blood samples were
  # centrifuged and mixed 1:1 with a MOPS stabilising buffer, and "the plasma
  # imipenem concentration was determined by ... HPLC-UV". The assayed
  # analyte is imipenem itself; cilastatin was not measured.
  compartmentData <- list(
    central     = list(analyte = "imipenem", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "imipenem", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description        = paste(
        "Endogenous creatinine clearance computed with the Cockcroft-Gault",
        "equation in SI creatinine units, reported as raw mL/min and NOT",
        "normalised to 1.73 m^2 body surface area. Table 1 footnote gives the",
        "exact form used: CrCl = [(140 - age) * weight (kg)] /",
        "[0.818 * Cr (umol/L)], multiplied by 0.85 for women."
      ),
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject. Power effect on CL:",
        "CL (L/h) = 11.357 * (CrCl / 99.896)^0.473 * exp(eta), the paper's",
        "Eq. 8 (Results, 'Population pharmacokinetics'), with the exponent",
        "0.473 also tabulated as the 'fCrCl' row of Table 2 (RSE 5.181%; 95%",
        "CI 0.425-0.521; bootstrap median 0.468, 95% CI 0.375-0.590). Note",
        "that the centring constant 99.896 mL/min stated in the text",
        "immediately below Eq. 8 ('99.896 is the median value of CrCl') does",
        "NOT equal the Table 1 median of 99.34 mL/min; the equation constant",
        "is the one this model uses, because it is the value the fitted",
        "typical CL of 11.357 L/h is conditioned on. Fitted over an observed",
        "range of 17.80-256.22 mL/min (Table 1 and Results, 'Simulation');",
        "the paper states explicitly that 'the predication model can only be",
        "used within the range (17.80-256.22 mL/min)'. CrCl was the only",
        "covariate retained; introducing it dropped the objective function",
        "from 776.148 to 734.330 and the CL inter-individual variability",
        "from 38.562% to 35.748% CV."
      ),
      source_name        = "CrCl"
    )
  )

  # Screened as candidate covariates during forward inclusion but not retained
  # in the final model (Methods, 'Final model'; Results, 'Population
  # pharmacokinetics': "After forward inclusion-backward elimination of all
  # the candidate covariates, CL was markedly influenced by serum creatinine
  # clearance (CrCl)"). The Discussion notes that AGE, CREAT and WT also
  # tracked the individual parameter variation, but that CrCl produced by far
  # the largest OFV drop (>30) and that the others were correlated with it.
  # No point estimate is published for any of them, so nothing is encoded.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age.",
      units       = "years",
      type        = "continuous",
      notes       = paste(
        "Screened; Table 1 gives 56.45 +/- 18.76 years, median 56, range",
        "18-96. The Discussion states age 'was related to the individual",
        "variation in the parameters' (online supplemental figure 1) but it",
        "was not retained once CrCl entered, and age is itself an input to",
        "the Cockcroft-Gault CrCl."
      ),
      source_name = "AGE"
    ),
    WT = list(
      description = "Body weight.",
      units       = "kg",
      type        = "continuous",
      notes       = paste(
        "Screened; Table 1 gives 70.21 +/- 72.01 kg, median 69, range",
        "19.6-311.7. Not retained; weight is also an input to the",
        "Cockcroft-Gault CrCl. The reported standard deviation exceeds the",
        "mean and the 311.7 kg maximum is implausible for the cohort, so the",
        "Table 1 weight row should not be reused as a simulation",
        "distribution."
      ),
      source_name = "BW"
    ),
    CREAT = list(
      description = "Serum creatinine.",
      units       = "umol/L",
      type        = "continuous",
      notes       = paste(
        "Screened; Table 1 gives 94.35 +/- 187.94 umol/L, median 64.8, range",
        "32.8-883. Figure 2 shows the base-model CL ETA against both Cr and",
        "CrCl; the relationship disappeared in the final model after CrCl was",
        "included. Not retained as a separate covariate."
      ),
      source_name = "Cr"
    ),
    SEXF = list(
      description = "Female sex indicator.",
      units       = "unitless",
      type        = "categorical",
      notes       = paste(
        "Screened as an indicator variable (Methods, 'Final model'). Cohort",
        "33 male / 18 female (Table 1). Not retained; no coefficient is",
        "reported. Sex enters the Cockcroft-Gault CrCl via the 0.85 factor."
      ),
      source_name = "GNDR"
    ),
    AKI = list(
      description = "Acute kidney injury indicator.",
      units       = "unitless",
      type        = "categorical",
      notes       = paste(
        "Screened as an indicator variable. 20 of 51 patients (39.21%) met",
        "the KDIGO 2021 AKI criteria (Table 1). Not retained; no coefficient",
        "is reported."
      ),
      source_name = "AKI"
    ),
    APACHE_II = list(
      description = "Acute Physiology and Chronic Health Evaluation II score.",
      units       = "points",
      type        = "continuous",
      notes       = paste(
        "Screened; Table 1 gives 16.67 +/- 6.44, median 15, range 8-33. Not",
        "retained; no coefficient is reported."
      ),
      source_name = "APACHE"
    ),
    SOFA = list(
      description = "Sequential Organ Failure Assessment score.",
      units       = "points",
      type        = "continuous",
      notes       = paste(
        "Screened; Table 1 gives 6.78 +/- 5.06, median 5, range 2-19. Not",
        "retained; no coefficient is reported."
      ),
      source_name = "SOFA"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 51L,
    n_studies      = 1L,
    n_observations = 196L,
    age_range      = "18-96 years",
    age_median     = "56 years",
    weight_range   = "19.6-311.7 kg as tabulated (Table 1 reports 70.21 +/- 72.01 kg; the standard deviation exceeds the mean and the upper bound is implausible, so treat this row as unreliable)",
    weight_median  = "69 kg",
    sex_female_pct = 100 * 18 / 51,
    race_ethnicity = c(Asian = 100),
    disease_state  = paste(
      "Critically ill adults admitted to the intensive care unit and meeting",
      "the Sepsis 3.0 diagnostic criteria, treated with imipenem as part of",
      "their anti-infection therapy. Comorbid state at enrolment: acute",
      "kidney injury 20/51 (39.21%), hepatic insufficiency (Child-Pugh B or",
      "C) 26/51 (50.98%), septic shock 18/51 (35.29%), mechanical",
      "ventilation 26/51 (50.98%). Infection site: respiratory 27/51",
      "(52.94%), bloodstream 2/51, other 22/51 (43.14%). Pregnant patients",
      "and patients with an estimated survival time of 48 hours or less were",
      "excluded. Patients on continuous renal replacement therapy were not",
      "excluded by protocol, but the Discussion contrasts this cohort with",
      "the CRRT-containing cohort of Li et al. as a reason for the clearance",
      "difference between the two studies."
    ),
    dose_range     = paste(
      "Imipenem-cilastatin 500 mg/500 mg given intravenously; imipenem doses",
      "of 0.5 g q6h, 0.5 g q8h and 1 g q8h were the commonly prescribed",
      "regimens, with the regimen chosen by the attending physician. Each",
      "dose was diluted into 100 mL and infused by pump over 1 hour. Table 1",
      "gives a per-dose median of 500 mg (range 500-1000) and a daily-dose",
      "median of 2000 mg (range 1500-3000). The paper's Monte Carlo",
      "simulations span six regimens: 500 and 1000 mg q6h, q8h and q12h, all",
      "as 1-hour infusions."
    ),
    sampling       = paste(
      "Samples were drawn at or near steady state, at least 24 hours after",
      "the start of imipenem therapy. The nominal schedule was duplicate",
      "samples before the dose (0 h) and at 0.5, 1, 1.5, 2, 3, 4, 6 and 8",
      "hours after it; 7 patients provided samples at all time points and 44",
      "provided 3-4 time points, for 196 samples in total (3.84 per",
      "patient). Plasma was stabilised 1:1 with MOPS buffer (0.126 M, pH",
      "6.8) because imipenem hydrolyses rapidly in plasma."
    ),
    renal_function = paste(
      "Cockcroft-Gault creatinine clearance 104.59 +/- 60.95 mL/min, median",
      "99.34, range 17.80-256.22 mL/min (Table 1). Serum creatinine 94.35",
      "+/- 187.94 umol/L, median 64.8, range 32.8-883. The paper stratifies",
      "its dosing simulations into four CrCl bands: 17.80-30, 30-60, 60-90",
      "and 90-256.22 mL/min."
    ),
    regions        = "People's Republic of China (single centre; intensive care unit, Beijing Friendship Hospital, Capital Medical University, Beijing).",
    notes          = paste(
      "Baseline demographics from Bai 2024 Table 1; modelling details from",
      "Methods, 'Population pharmacokinetics/validation'. Prospective",
      "open-label study (ethics certificate 2018-P2-219-01). The model was",
      "fit in Phoenix NLME 8.0 with the first-order conditional",
      "estimation-least squares (FOCE-ELS) method, and evaluated with",
      "goodness-of-fit plots, a 1000-replicate visual predictive check and a",
      "1000-sample bootstrap. Imipenem was assayed by HPLC-UV at 298 nm with",
      "ceftazidime as the internal standard, LLOQ 0.3 ug/mL, linear over",
      "0.3-200.0 ug/mL. The unbound fraction used for the paper's fT>MIC",
      "target-attainment simulations is not stated, so no protein-binding",
      "term is encoded here."
    )
  )

  ini({
    # =========================================================================
    # Structural disposition parameters (Table 2, 'Model estimate' /
    # 'Estimate' column). Typical values refer to the covariate reference
    # subject, CrCl 99.896 mL/min (Eq. 8).
    # =========================================================================
    lcl <- log(11.357); label("Clearance (CL, L/h) at CrCl 99.896 mL/min")  # Table 2 row 'CL (L/h)' 11.357 (RSE 3.024%; 95% CI 10.679-12.035; bootstrap median 11.367, 95% CI 9.248-14.227); also the leading coefficient of Eq. 8
    lvc <- log(16.378); label("Central volume of distribution (Vc, L)")     # Table 2 row 'Vc (L)' 16.378 (RSE 2.283%; 95% CI 15.641-17.116; bootstrap median 16.460, 95% CI 5.110-32.423)
    lvp <- log(10.904); label("Peripheral volume of distribution (Vp, L)")  # Table 2 row 'Vp (L)' 10.904 (RSE 3.951%; 95% CI 10.054-11.754; bootstrap median 12.505, 95% CI 8.434-73.773)
    lq  <- log(7.645);  label("Intercompartmental clearance (Q, L/h)")      # Table 2 row 'Q (L/h)' 7.645 (RSE 4.086%; 95% CI 7.029-8.261; bootstrap median 7.539, 95% CI 2.919-17.858)

    # =========================================================================
    # Covariate effect. Eq. 7 gives the general continuous-covariate form
    # Pi = P * (COV / COV_median)^f * exp(eta); Eq. 8 instantiates it as
    # CL (L/h) = 11.357 * (CrCl / 99.896)^0.473 * exp(eta).
    # =========================================================================
    e_crcl_cl <- 0.473; label("Power exponent of CrCl on CL (unitless)")  # Table 2 row 'fCrCl' 0.473 (RSE 5.181%; 95% CI 0.425-0.521; bootstrap median 0.468, 95% CI 0.375-0.590); Eq. 8

    # =========================================================================
    # Inter-individual variability. Methods, 'Base model': "An exponential
    # variability model was selected to describe the inter-individual
    # variability (IIV): Pi = P * exp(eta_i) [Eq. 1] ... IIV is assumed to
    # follow a log-normal distribution, and the random variable eta_i is
    # normally distributed with a mean of zero and variance of omega^2."
    #
    # SCALE. Table 2 heads the column "Iiv (CV%)", so the tabulated numbers
    # are apparent coefficients of variation of the log-normal distribution,
    # not variances and not omega standard deviations expressed as percents.
    # Two internal cross-checks support that reading. (1) The companion
    # residual-error row is labelled "Residual error (proportional error,
    # CV%)" and carries 30.370 in the same column, which for a proportional
    # error model can only be a CV -- so the column is CV throughout.
    # (2) The Table 2 footnote and the Results text quote the CL entry as a
    # percentage twice: "The IIV decreased from 38.562% (base model) to
    # 35.748% (final model)". Converted with omega^2 = log(CV^2 + 1):
    #   Vc  CV 89.853% -> log(0.89853^2 + 1) = 0.5918651
    #   Vp  CV  8.319% -> log(0.08319^2 + 1) = 0.0068967
    #   Q   CV 24.453% -> log(0.24453^2 + 1) = 0.0580754
    #   CL  CV 35.748% -> log(0.35748^2 + 1) = 0.1202617
    # The paper reports a diagonal omega matrix; no correlations are given.
    # =========================================================================
    etalvc ~ 0.5918651  # CV 89.853% (Table 2 'Vc (L)' Iiv column; eta-shrinkage 0.466)
    etalvp ~ 0.0068967  # CV  8.319% (Table 2 'Vp (L)' Iiv column; eta-shrinkage 0.450)
    etalq  ~ 0.0580754  # CV 24.453% (Table 2 'Q (L/h)' Iiv column; eta-shrinkage 0.249)
    etalcl ~ 0.1202617  # CV 35.748% (Table 2 'CL (L/h)' Iiv column; eta-shrinkage 0.456)

    # =========================================================================
    # Residual error. Methods, 'Base model': "the proportional error model was
    # employed to calculate the residual error of the pharmacokinetic model:
    # Ci = C * (1 + eps) [Eq. 2] ... The proportional error was characterised
    # using eps, which is distributed with a mean of zero and variances of
    # sigma^2." Table 2 reports it as a CV%, i.e. a proportional SD of 30.370%.
    # =========================================================================
    propSd <- 0.30370; label("Proportional residual error (fraction)")  # Table 2 'Residual error (proportional error, CV%)' 30.370 (RSE 5.126%; 95% CI 27.298-33.442; bootstrap median 30.537, 95% CI 23.936-38.133)
  })

  model({
    # -----------------------------------------------------------------------
    # 1. Individual disposition parameters. Only CL carries a covariate
    #    (Eq. 8); Vc, Vp and Q are covariate-free with exponential IIV.
    # -----------------------------------------------------------------------
    cl <- exp(lcl + etalcl) * (CRCL / 99.896)^e_crcl_cl
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp + etalvp)
    q  <- exp(lq  + etalq)

    # -----------------------------------------------------------------------
    # 2. Micro-rate constants.
    # -----------------------------------------------------------------------
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # -----------------------------------------------------------------------
    # 3. Two-compartment intravenous disposition, Eqs. 3-6:
    #      dXc/dt = K0 - Q * (Cc - Cp) - CL * Cc,  Xc(0) = 0
    #      dXp/dt =      Q * (Cc - Cp),            Xp(0) = 0
    #      Cc = Xc / Vc,  Cp = Xp / Vp
    #    Substituting Cc = central/vc and Cp = peripheral1/vp turns the
    #    clearance form into the equivalent micro-constant form below. K0 is
    #    the zero-order infusion rate, supplied by the event table (imipenem
    #    was infused over 1 hour), so it does not appear as a model term.
    #    Doses in mg with volumes in L give Cc directly in mg/L (= ug/mL, the
    #    assay's reporting unit).
    # -----------------------------------------------------------------------
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # -----------------------------------------------------------------------
    # 4. Observation: total plasma imipenem concentration (Eq. 5).
    # -----------------------------------------------------------------------
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
