Fukumoto_2023_meropenem <- function() {
  description <- paste(
    "Two-compartment IV population PK model for meropenem in 31 adults with",
    "sepsis in a Japanese emergency center / ICU (Fukumoto 2023). Clearance",
    "scales as a power function of the BSA-normalized MEASURED creatinine",
    "clearance obtained from an 8-hour timed urine collection,",
    "CL = 13.6 * (CRCL / 87.6)^0.66 L/h, centered on the cohort-median",
    "87.6 mL/min/1.73 m^2. The measured CCr was the single retained covariate",
    "and outperformed Cockcroft-Gault CCr and eGFR, which underestimate renal",
    "function in septic patients with normal or augmented renal clearance.",
    "Interindividual variability was estimated on CL and V1 only; the random",
    "effects on Q and V2 were fixed at zero. Residual error is proportional."
  )
  reference <- paste(
    "Fukumoto S, Ohbayashi M, Okada A, Kohyama N, Tamatsukuri T, Inoue H,",
    "Kato A, Kotani T, Sagara H, Dohi K, Kogo M. (2023).",
    "Population pharmacokinetic model and dosing simulation of meropenem",
    "using measured creatinine clearance for patients with sepsis.",
    "Ther Drug Monit 45(3):392-399. doi:10.1097/FTD.0000000000001040"
  )
  vignette <- "Fukumoto_2023_meropenem"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Derived mechanically; verified = FALSE means it has
  # NOT been checked against the source paper.
  compartmentData <- list(
    central     = list(analyte = "meropenem", units = "mg", specimen = "plasma", verified = FALSE),
    peripheral1 = list(analyte = "meropenem", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    CRCL = list(
      description        = paste(
        "Measured creatinine clearance from an 8-hour timed urine",
        "collection, BSA-normalized to 1.73 m^2. NOT a creatinine-based",
        "estimating equation: urinary creatinine was assayed directly."
      ),
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Source column 'measured CCr'. Fukumoto 2023 Equation 1 (p. 393):",
        "measured CCr = (urine Cr [mg/dL] / SCr [mg/dL]) * urine output",
        "[mL/min] * (1.73 m^2 / BSA [m^2]). The trailing 1.73/BSA factor",
        "makes the value BSA-normalized despite the paper labeling the",
        "units simply 'mL/min'; it is therefore stored under the canonical",
        "CRCL column in its BSA-normalized sense per",
        "inst/references/covariate-columns.md (precedent:",
        "Xu_2019_sarilumab.R, MedellinGaribay_2015_gentamicin.R).",
        "Most patients had indwelling urinary catheters, so urine output",
        "was collected accurately; 8-hour rather than 24-hour collection",
        "was chosen to minimize circadian fluctuation of renal function",
        "(Discussion, p. 397).",
        "Cohort median 87.6 mL/min/1.73 m^2 (range 12.3-223), used as the",
        "centering reference in CL = theta2 * (CRCL/87.6)^theta5",
        "(Table 2, p. 395). 8/31 patients (26%) had CCr >= 130 (augmented",
        "renal clearance) and 3/31 (10%) had CCr <= 30 mL/min.",
        "The paper's central finding is that this measured CCr beat the",
        "Cockcroft-Gault estimate (CG-CCr, median 81.1 mL/min) and the",
        "Japanese eGFR equation (median 86.0 mL/min/1.73 m^2) on dOFV;",
        "CG-CCr significantly UNDERestimated renal function in the",
        "CCr >= 90 mL/min subgroup (Supplementary Table S1, not on disk).",
        "Treated as time-fixed per subject: a single 8-hour collection was",
        "paired with each patient's PK sampling occasion."
      ),
      source_name        = "measured CCr"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age at the time of blood collection",
      units       = "year",
      type        = "continuous",
      notes       = paste(
        "Screened on CL; reached significance in the forward-addition step",
        "(Fukumoto 2023 Results 'PPK Model', p. 395: 'Age, Alb, SCr, eGFR,",
        "measured CCr, CG-CCr, and SOFA score were selected as the",
        "significant covariates') but was not retained in the final model,",
        "which incorporates measured CCr on CL alone. No point estimate is",
        "published for the age effect (the screening dOFV table,",
        "Supplementary Table S2, is not on disk).",
        "Cohort median 72 years (range 18-94), Table 1."
      )
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/dL",
      type        = "continuous",
      notes       = paste(
        "Screened; significant on forward addition but not retained in the",
        "final model. No published point estimate. Cohort median 2.10 g/dL",
        "(range 1.50-3.10), Table 1 -- i.e. uniformly hypoalbuminemic, as",
        "expected in sepsis. Reported in g/dL by the source paper; the",
        "canonical ALB unit is g/L (multiply by 10)."
      )
    ),
    CREAT = list(
      description = "Serum creatinine",
      units       = "mg/dL",
      type        = "continuous",
      notes       = paste(
        "Screened; significant on forward addition but not retained. No",
        "published point estimate. Cohort median 0.69 mg/dL (range",
        "0.29-2.39), Table 1. The paper argues SCr is an unreliable renal",
        "marker in sepsis because it fluctuates with fluid resuscitation",
        "(Discussion, p. 397, citing Arulkumaran et al.), which motivates",
        "the measured-CCr covariate that was retained."
      )
    ),
    CRP = list(
      description = "C-reactive protein",
      units       = "mg/dL",
      type        = "continuous",
      notes       = paste(
        "Screened on CL and V1 as a marker of sepsis severity / vascular",
        "permeability; NOT significant and not retained. Cohort median",
        "11.0 mg/dL (range 0.72-41.0), Table 1. Note the Discussion",
        "(p. 397) quotes a median CRP of 9.25 mg/dL over the same range --",
        "an internal inconsistency in the source paper; Table 1 is taken as",
        "authoritative here."
      )
    ),
    SEXF = list(
      description = "Sex (1 = female, 0 = male)",
      units       = NULL,
      type        = "categorical",
      notes       = paste(
        "Screened; NOT significant and not retained. Cohort 25 men / 6",
        "women (19.4% female), Table 1."
      )
    ),
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = paste(
        "Screened on CL and V1; NOT significant and not retained -- so no",
        "allometric scaling appears in the final model and V1 / V2 / Q are",
        "absolute volumes and flows for the cohort. Cohort median 53.7 kg",
        "(range 35.3-91.7), Table 1. The Discussion (p. 397) notes body",
        "weight in critically ill patients is transiently inflated by fluid",
        "resuscitation, which is one reason weight-based Cockcroft-Gault",
        "performs poorly in this cohort."
      )
    ),
    SOFA = list(
      description = "Sequential Organ Failure Assessment score",
      units       = "points",
      type        = "continuous",
      notes       = paste(
        "Screened as a sepsis-severity covariate; significant on forward",
        "addition but not retained in the final model. No published point",
        "estimate. Cohort median 8 (range 2-15) on admission and 6 (range",
        "0-11) at the time of blood collection, Table 1. No canonical",
        "covariate column exists for SOFA in",
        "inst/references/covariate-columns.md; because the covariate is",
        "documentation-only here (never referenced in model()), no new",
        "canonical is proposed by this extraction."
      )
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 31L,
    n_studies        = 1L,
    n_concentrations = 100L,
    age_range        = "18-94 years",
    age_median       = "72 years",
    weight_range     = "35.3-91.7 kg",
    weight_median    = "53.7 kg",
    sex_female_pct   = 19.4,
    race_ethnicity   = "Not reported (single-center Japanese cohort)",
    disease_state    = paste(
      "Adults with sepsis admitted to the emergency center or intensive care",
      "unit, diagnosed by specialized physicians. Median SOFA score 8",
      "(range 2-15) on admission and 6 (range 0-11) at the time of blood",
      "collection. Median CRP 11.0 mg/dL (range 0.72-41.0) and WBC",
      "13.0 x 10^3/uL (range 2.70-30.6), indicating substantial systemic",
      "inflammation. Median serum albumin 2.10 g/dL (range 1.50-3.10).",
      "Vasopressors in 9/31 (29%); intravenous diuretics in 5/31 (16%);",
      "median daily intravenous fluid volume 1000 mL (range 200-3900).",
      "Exclusions: age < 18 years, hemodialysis, massive bleeding,",
      "pregnancy, and death during the urine collection."
    ),
    dose_range       = paste(
      "Meropenem IV 0.5 g every 8 or 12 hours, or 1 g every 8 or 12 hours,",
      "each administered as a 0.5-hour infusion. The regimen for each",
      "patient was chosen by the treating physician based on sepsis severity",
      "and renal function."
    ),
    regions          = paste(
      "Japan (single center: Showa University Hospital, Tokyo).",
      "Enrolled June 2016-August 2021."
    ),
    renal_function   = paste(
      "Broad range by design. Measured CCr (8-hour urine collection,",
      "BSA-normalized) median 87.6 mL/min/1.73 m^2 (range 12.3-223);",
      "Cockcroft-Gault CCr median 81.1 mL/min (range 14.0-246); eGFR median",
      "86.0 mL/min/1.73 m^2 (range 15.6-222); serum creatinine median",
      "0.69 mg/dL (range 0.29-2.39). 8/31 patients (26%) had augmented renal",
      "clearance (CCr >= 130 mL/min) and 3/31 (10%) had CCr <= 30 mL/min.",
      "Patients on hemodialysis were excluded."
    ),
    notes            = paste(
      "Demographics from Fukumoto 2023 Table 1 (p. 394). Prospective",
      "observational study; 32 patients enrolled, 1 excluded, 31 analyzed.",
      "Sparse sampling: 2-5 samples per patient (median 2) within a single",
      "dosing interval, at predose (trough) and 0.5, 1, 2, 4 and 8 hours",
      "post-dose, after at least 3 doses (i.e. at steady state).",
      "Total meropenem in plasma assayed by HPLC-UV at 300 nm with",
      "cefotaxime internal standard; calibration linear over 0.1-80 mcg/mL",
      "(R > 0.99), inter-measurement trueness and accuracy < 15%.",
      "Estimation: Phoenix NLME 8.1 (Certara), first-order conditional",
      "estimation-extended least squares (FOCE-ELS). Covariate selection by",
      "stepwise forward addition (P < 0.05) then backward elimination",
      "(P < 0.01). Validation: GOF plots, prediction-corrected VPC, and a",
      "1000-replicate bootstrap with a 100% success rate. Reported",
      "eta-shrinkage 41.4% for V1 and 7.75% for CL -- the V1 shrinkage is",
      "high and the authors flag the 127% RSE on IIV(V1) as a study",
      "limitation attributable to sparse sampling of the elimination phase."
    )
  )

  ini({
    # Structural parameters. Fukumoto 2023 Table 2 (p. 395), "Population
    # mean" block; the covariate model is printed verbatim in the row
    # headings of that table:
    #   V1 = theta1
    #   CL = theta2 * (measured CCr / 87.6)^theta5
    #   V2 = theta3
    #   Q  = theta4
    # Column order in Table 2 is: units | estimate (RSE %) | bootstrap
    # median | bootstrap 2.5% | bootstrap 97.5% | bias (%).
    lvc <- log(26.5); label("Central volume V1 (L)")                                # Table 2 theta1 = 26.5 L (RSE 13.6%); bootstrap median 26.6, 95% CI 21.4-36.0
    lcl <- log(13.6); label("Clearance CL at CRCL = 87.6 mL/min/1.73 m^2 (L/h)")    # Table 2 theta2 = 13.6 L/h (RSE 5.04%); bootstrap median 13.6, 95% CI 12.4-15.2
    lvp <- log(13.2); label("Peripheral volume V2 (L)")                             # Table 2 theta3 = 13.2 L (RSE 21.2%); bootstrap median 13.4, 95% CI 9.51-20.0
    lq  <- log(9.80); label("Intercompartmental clearance Q (L/h)")                 # Table 2 theta4 = 9.80 L/h (RSE 43.0%); bootstrap median 9.59, 95% CI 2.29-20.3

    # Covariate effect: power exponent on BSA-normalized measured creatinine
    # clearance, centered at the cohort-median 87.6 mL/min/1.73 m^2 (the
    # centering value is printed inside the Table 2 row heading and matches
    # the Table 1 median exactly).
    e_crcl_cl <- 0.66; label("Power exponent on (CRCL / 87.6) for CL (unitless)")   # Table 2 theta5 = 0.66 (RSE 11.9%); bootstrap median 0.67, 95% CI 0.52-0.84

    # Interindividual variability. Fukumoto 2023 Methods (p. 394): "The
    # interindividual variability was evaluated using an exponential error
    # model", i.e. theta_i = theta_TV * exp(eta_i). Table 2 reports the IIV
    # rows in units of "%", read here as %CV; for a log-normal eta the exact
    # variance is omega^2 = log(CV^2 + 1).
    etalvc ~ 0.051538  # log(0.230^2 + 1); Table 2 IIV V1 = 23.0% CV (RSE 127%); bootstrap median 21.5%, 95% CI 0.0003-56.7%. eta-shrinkage 41.4%.
    etalcl ~ 0.047686  # log(0.221^2 + 1); Table 2 IIV CL = 22.1% CV (RSE 32.4%); bootstrap median 21.1%, 95% CI 12.7-27.7%. eta-shrinkage 7.75%.

    # IIV on Q and V2 was estimated but did not improve the fit and was held
    # at zero in the final model. Fukumoto 2023 Results (p. 395): "The
    # effective interindividual variability on intercompartmental clearance
    # (Q) and peripheral Vd (V2) did not significantly improve the model fit.
    # Therefore, the random effects were fixed at zero."
    etalq  ~ fixed(0)  # Results p. 395: random effect on Q zero
    etalvp ~ fixed(0)  # Results p. 395: random effect on V2 zero

    # Residual error. Fukumoto 2023 Results (p. 395): "The final residual
    # error model was described using the proportional error model."
    propSd <- 0.194; label("Proportional residual error (fraction)")                # Table 2 "Residual variability / Proportional (%)" = 19.4% (RSE 13.1%); bootstrap median 18.4%, 95% CI 13.0-23.4%
  })

  model({
    # Individual PK parameters. CL carries a power-function effect of the
    # BSA-normalized measured creatinine clearance centered at the cohort
    # median, exactly as printed in the Table 2 row heading
    # CL = theta2 * (measured CCr / 87.6)^theta5. No weight / allometric
    # scaling appears in the final model -- body weight was screened and not
    # retained -- so V1, V2 and Q are absolute for this cohort.
    cl <- exp(lcl + etalcl) * (CRCL / 87.6)^e_crcl_cl
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp + etalvp)
    q  <- exp(lq  + etalq)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Two-compartment ODE system. Meropenem is given intravenously; the dose
    # enters `central` with an infusion duration supplied at simulation time
    # (0.5 h in the clinical study; 0.5 / 3 / 8 / 12 h in the paper's Monte
    # Carlo dosing simulations).
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                  k12 * central - k21 * peripheral1

    # Dose in mg and volumes in L give Cc in mg/L, which equals the mcg/mL
    # units the paper reports concentrations and MICs in. Meropenem plasma
    # protein binding is ~2%, so total and free concentrations are
    # operationally interchangeable for the % T > MIC target.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
