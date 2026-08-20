Morath_2025_apixaban <- function() {
  description <- paste0(
    "One-compartment population pharmacokinetic model with first-order ",
    "absorption for oral apixaban in adults with postoperative atrial ",
    "fibrillation after cardiac surgery (Morath 2025), quantifying the ",
    "amiodarone drug-drug interaction. Apparent oral clearance CL/F = 3.05 ",
    "L/h at a creatinine clearance of 100 mL/min without amiodarone, scaled ",
    "by a power function of Cockcroft-Gault creatinine clearance ",
    "(CRCL/100)^0.279 and reduced by a multiplicative factor of 0.679 (a ",
    "32% reduction) during concomitant amiodarone therapy. Apparent volume ",
    "of distribution Vd/F = 23.7 L and absorption rate constant ka = 0.652 ",
    "1/h, neither carrying covariates. Interindividual variability was ",
    "supported on CL/F only (29.4% CV) and not on Vd/F or ka; ",
    "interoccasion variability was not supported on any parameter. ",
    "Residual variability is proportional (31.4% CV). Apixaban dose was not ",
    "a significant covariate on bioavailability."
  )
  reference <- paste0(
    "Morath B, Foerster KI, Chiriac U, Zaradzki M, Hoppe-Tichy T, ",
    "Schrey D, Burhenne J, Czock D, Karck M, Haefeli WE, Wicha SG. ",
    "Effect of amiodarone on apixaban exposure in patients after cardiac ",
    "surgery - a population pharmacokinetic study. ",
    "Clin Pharmacokinet. 2025;64(8):1191-1201. ",
    "doi:10.1007/s40262-025-01534-z."
  )
  vignette <- "Morath_2025_apixaban"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Morath 2025 Methods sections 2.3
  # (apixaban plasma concentrations by LC-MS/MS) and 3.1 (one-compartment
  # model with first-order absorption of oral apixaban tablets).
  compartmentData <- list(
    depot   = list(analyte = "apixaban", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "apixaban", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description        = "Creatinine clearance (Cockcroft-Gault)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Raw Cockcroft-Gault creatinine clearance in mL/min, NOT ",
        "BSA-normalised (Morath 2025 Methods 2.4: 'creatinine clearance ",
        "(CrCL) according to Cockcroft and Gault's equation on CL/F'). ",
        "Cohort 76.2 +/- 33.3 mL/min overall (control cohort 77.3 +/- 38.8, ",
        "interaction cohort 74.8 +/- 25.3; Table 1). Reference value 100 ",
        "mL/min in the CL/F equation (CrCL/100)^0.279 (Morath 2025 ",
        "Discussion equation and Table 2 COV-CL-CrCL); the paper centres to ",
        "a standard value of 100 rather than to the cohort median, per ",
        "Methods 2.4 ('centered around the median covariate value of the ",
        "population or centered to standard values'). The exposure ",
        "simulations (Table S1) span CrCL 15-100 mL/min, which is wider ",
        "than the fitting cohort; extrapolation below about 30 mL/min rests ",
        "on the power model rather than on observed data. The raw ",
        "Cockcroft-Gault scale follows the register's existing raw-CrCl ",
        "precedent (Delattre 2010 amikacin, Chen 2023 nemonoxacin, Wada ",
        "2023 sparsentan, Shu 2024 posaconazole, Ueshima 2018 apixaban)."
      ),
      source_name        = "CrCL"
    ),
    CONMED_AMIO = list(
      description        = "Concomitant amiodarone therapy indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant amiodarone)",
      notes              = paste0(
        "1 = subject is receiving concomitant amiodarone for postoperative ",
        "atrial fibrillation, 0 = apixaban alone. Cohort prevalence 15/33 ",
        "(45%) (Morath 2025 Results). May be time-varying: amiodarone start ",
        "and stop were captured per patient and are shaded in ",
        "Supplementary Fig. S1. IMPLEMENTATION NOTE - the source paper did ",
        "not enter amiodarone as a plain indicator column during ",
        "estimation. Methods 2.4 states that an established amiodarone PK ",
        "model was encoded inside the same NONMEM control stream, typical ",
        "amiodarone profiles were simulated, and a threshold amiodarone ",
        "concentration for onset of inhibition was to be estimated. That ",
        "threshold was not estimable and 'the threshold was set to a small ",
        "value of 0.001 mg/L' (Results 3.1), which is roughly three orders ",
        "of magnitude below any therapeutic amiodarone concentration, so ",
        "the inhibition is switched on for the entire duration of ",
        "amiodarone therapy. The paper also reports 'no difference between ",
        "amiodarone loading or maintenance doses were observed'. The final ",
        "model is therefore written by the authors themselves in binary ",
        "form: 'Individual CL/F = 3.05 x (CrCL/100)^0.279 x 0.679 (for ",
        "amiodarone comedication)' (Morath 2025 Discussion). This model ",
        "file encodes that published equation directly, so the simulated ",
        "amiodarone PK layer and its 0.001 mg/L threshold are not ",
        "reproduced - they carry no information beyond the on/off switch. ",
        "Per the register convention the per-paper definition is 'any ",
        "concurrent amiodarone use', independent of dose or loading status."
      ),
      source_name        = "amiodarone comedication"
    )
  )

  # Covariates screened by Morath 2025 but NOT retained in the final model.
  # Documentation only: these are not referenced in model(). Recorded here so
  # the provenance of the paper's covariate screen is preserved.
  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = paste0(
        "Tested on both CL/F and Vd/F using linear and power models ",
        "(Morath 2025 Methods 2.4) and not retained. Discussion: ",
        "'Previously observed covariates such as age or body weight were ",
        "not confirmed in our patients'. Cohort 83.0 kg median, range ",
        "61.3-136.0 kg (Table 1). No point estimate is reported, so no ",
        "effect can be encoded."
      )
    ),
    DOSE = list(
      description = "Administered apixaban dose level (2.5 or 5 mg b.i.d.)",
      units       = "mg",
      type        = "categorical",
      notes       = paste0(
        "Tested as a categorical covariate on bioavailability F and not ",
        "retained (dOFV -3.386, P = 0.065; Morath 2025 Results 3.1). The ",
        "authors flag this as a power-limited negative result: 'the small ",
        "sample size and few samples in the absorption phase may have ",
        "prevented the estimation of a significant effect of apixaban dose ",
        "on bioavailability. Hence, the absence of such an effect cannot be ",
        "excluded' (Limitations). Contrasts with Cirincione 2018, which did ",
        "find a dose effect on F."
      )
    ),
    CONMED_METAMIZOLE = list(
      description = "Concomitant metamizole (dipyrone) therapy indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste0(
        "Tested on CL/F and not retained (Morath 2025 Results 3.1: 'There ",
        "was no interaction between metamizole, rosuvastatin, or ",
        "atorvastatin with CL/F'). Cohort prevalence 31/33 (93.9%), which ",
        "leaves almost no unexposed reference group and so very little ",
        "power to resolve an effect (Table 1)."
      )
    ),
    CONMED_ROSUVASTATIN = list(
      description = "Concomitant rosuvastatin therapy indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste0(
        "Tested on CL/F as a potential BCRP-competing comedication and not ",
        "retained (Morath 2025 Results 3.1). Cohort prevalence 5/33 ",
        "(15.1%) (Table 1)."
      )
    ),
    CONMED_ATORVASTATIN = list(
      description = "Concomitant atorvastatin therapy indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste0(
        "Tested on CL/F as a potential CYP3A4- / transporter-competing ",
        "comedication and not retained (Morath 2025 Results 3.1). Cohort ",
        "prevalence 19/33 (57.6%) (Table 1)."
      )
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 33L,
    n_studies        = 1L,
    n_observations   = 76L,
    age_range        = "32.0-87.0 years (median 72.0)",
    age_median       = "72.0 years",
    weight_range     = "61.3-136.0 kg (median 83.0)",
    weight_median    = "83.0 kg",
    sex_female_pct   = 12.1,
    disease_state    = paste0(
      "Adults (age >= 18 years) with postoperative atrial fibrillation ",
      "after cardiac surgery at a German tertiary-care university hospital, ",
      "receiving oral apixaban for stroke prophylaxis either alone ",
      "(control cohort, n = 18) or with concomitant amiodarone ",
      "(interaction cohort, n = 15). Most frequent surgery was coronary ",
      "artery bypass graft (15/33, 45.5%), followed by aortic conduits for ",
      "aortic aneurysm (10/33, 30.3%). Comorbidities: arterial ",
      "hypertension 84.8%, hyperlipidemia 72.7%, type 2 diabetes mellitus ",
      "42.4%, past acute coronary syndrome 27.7%, previous stroke 12.1%. ",
      "Median left ventricular ejection fraction 53% (range 20-70%). ",
      "Median body mass index 27.6 +/- 4.45 kg/m^2. No patient had severe ",
      "liver impairment. Patients with a mechanical valve prosthesis were ",
      "excluded. No CYP inhibitors or inducers per Flockhart's table were ",
      "prescribed."
    ),
    dose_range       = paste0(
      "Oral apixaban 5 mg or 2.5 mg twice daily (b.i.d.) per the drug ",
      "label. Concentrations were measured in residual blood at the time ",
      "of routine clinical sampling and are therefore not restricted to ",
      "trough samples, giving a spread of times after dose. Apixaban ",
      "plasma concentrations by LC-MS/MS, lower limit of quantification ",
      "0.5 ng/mL; the assay met ICH M10 bioanalytical validation ",
      "recommendations."
    ),
    regions          = "Germany (single-centre prospective study, Heidelberg University Hospital; ethics approval S-523/2021)",
    co_medication    = "Concomitant amiodarone in 15/33 (45%); metamizole 93.9%, atorvastatin 57.6%, rosuvastatin 15.1%",
    creatinine_clearance = "76.2 +/- 33.3 mL/min overall (Cockcroft-Gault; raw mL/min, not BSA-normalised)",
    notes            = paste0(
      "Baseline demographics in Table 1; final-model parameter estimates ",
      "in Table 2. Software: NONMEM 7.5 with FOCE-I, executed through ",
      "Perl-speaks-NONMEM 5.0.0; exposure simulations in Simulx 2023R1. ",
      "Parameter uncertainty was assessed by log-likelihood-profiling-based ",
      "sampling-importance-resampling (LLP-SIR) rather than bootstrap, ",
      "because bootstrap is unreliable in small-N datasets; Table 2 ",
      "therefore reports 95% confidence intervals and not RSEs. Model ",
      "evaluation used a prediction-corrected VPC with 500 simulations ",
      "(Fig. 1), PRED/IPRED-versus-DV plots (Fig. 2), and NPDE plots ",
      "(Fig. 3). A two-compartment model did not improve the fit (dOFV ",
      "-2.09). Covariates tested and NOT retained: body weight (on CL/F ",
      "and Vd/F), apixaban dose (on F), and comedications metamizole, ",
      "rosuvastatin and atorvastatin (on CL/F) - see ",
      "covariatesDataExcluded. Sample size was planned for 15 patients per ",
      "group (80% power, effect size 1.5, CV 0.37, alpha 0.05). Limitations ",
      "noted by the authors: real-world residual-blood sampling, sparse ",
      "absorption-phase data (hence the imprecise ka), small N, and the PK ",
      "of the active amiodarone metabolite N-desethylamiodarone was not ",
      "modelled, so the long-term-therapy setting may not be fully ",
      "represented."
    )
  )

  ini({
    # Structural parameters - Morath 2025 Table 2 final estimates, with
    # LLP-SIR 95% confidence intervals given in the trailing comments.

    # Absorption rate constant. The authors flag ka as the one imprecisely
    # estimated structural parameter (Results 3.1: 'Apart from ka and the
    # proportional residual variability, the model parameters were estimated
    # with adequate precision given the smallN dataset'; Limitations:
    # 'limited data for the absorption phase impacted precision of the
    # estimated parameter ka'). It was estimated, not fixed.
    lka <- log(0.652); label("Absorption rate constant ka (1/h)")  # Morath 2025 Table 2: ka = 0.652 /h (95% CI 0.338-1.345)

    # Apparent oral clearance at the reference covariate values, i.e. CrCL =
    # 100 mL/min and no amiodarone comedication. This is the 3.05 L/h
    # intercept of the published covariate equation 'Individual CL/F = 3.05 x
    # (CrCL/100)^0.279 x 0.679 (for amiodarone comedication)'.
    lcl <- log(3.05); label("Apparent oral clearance CL/F (L/h) at CRCL = 100 mL/min without amiodarone")  # Morath 2025 Table 2: CL/F = 3.05 L/h (95% CI 2.54-3.71)

    # Apparent volume of distribution. No covariates were retained on Vd/F.
    lvc <- log(23.7); label("Apparent volume of distribution Vd/F (L)")  # Morath 2025 Table 2: Vd/F = 23.7 L (95% CI 16.11-31.7)

    # Power exponent on (CRCL / 100) for CL/F. Retained after backward
    # elimination (dOFV -6.758, P = 0.0093).
    e_crcl_cl <- 0.279; label("Power exponent on (CRCL / 100) for CL/F (unitless)")  # Morath 2025 Table 2: COV-CL-CrCL = 0.279 (95% CI 0.073-0.506)

    # Multiplicative factor on CL/F during concomitant amiodarone therapy.
    # Retained after backward elimination (dOFV -14.327, P = 0.00015).
    # Encoded as the LOG of the published multiplicative factor 0.679 so that
    # model() can apply it as exp(e_amio_cl * CONMED_AMIO): the exponent
    # collapses to 0 (factor 1) without amiodarone and to log(0.679) (factor
    # 0.679) with amiodarone. Log encoding also keeps the parameter
    # unconstrained if the model is re-fitted, and matches the sibling
    # apixaban model Ueshima_2018_apixaban.R. The reduction is about 32%,
    # which the Discussion rounds to 'reduced apixaban CL/F by approximately
    # 30%'.
    e_amio_cl <- log(0.679); label("Log of the multiplicative factor on CL/F for concomitant amiodarone")  # Morath 2025 Table 2: COV-CL/F-amiodarone = 0.679 (95% CI 0.556-0.820)

    # Interindividual variability - Morath 2025 Table 2. IIV was supported on
    # CL/F only, not on Vd/F or ka (Results 3.1), and interoccasion
    # variability was not supported on any parameter; the model therefore
    # carries a single eta. Table 2 reports IIV CL/F = 29.4% CV with the
    # footnote '%CV for IIV calculated by sqrt(e^omega^2 - 1)', so the
    # log-scale variance is recovered by inverting that transform:
    # omega^2 = log(CV^2 + 1) = log(0.294^2 + 1) = 0.0829 (omega = 0.288).
    etalcl ~ 0.0829  # Morath 2025 Table 2: IIV CL/F = 29.4% CV (95% CI 20.4-43.7); omega^2 = log(0.294^2 + 1) = 0.0829

    # Residual unexplained variability - proportional error, reported as
    # 31.4% CV. The authors note this parameter, like ka, was estimated with
    # limited precision (Results 3.1). A proportional model was selected over
    # additive and combined alternatives (Methods 2.4).
    propSd <- 0.314; label("Proportional residual error (fraction)")  # Morath 2025 Table 2: Prop. RUV = 31.4% CV (95% CI 25.8-39.3)
  })

  model({
    # 1. Individual parameters.
    # CL/F carries the two retained covariate effects and the single eta.
    # The published equation is 'Individual CL/F = 3.05 x (CrCL/100)^0.279 x
    # 0.679 (for amiodarone comedication)' (Morath 2025 Discussion), with
    # log-normal IIV applied multiplicatively via exp(etalcl).
    cl <- exp(lcl + etalcl) * (CRCL / 100)^e_crcl_cl * exp(e_amio_cl * CONMED_AMIO)

    # No IIV and no covariates on Vd/F or ka (Morath 2025 Results 3.1).
    vc <- exp(lvc)
    ka <- exp(lka)

    # 2. Micro-constants.
    kel <- cl / vc

    # 3. One-compartment ODE system with first-order oral absorption.
    # Bioavailability F is not identifiable from oral-only data and is
    # absorbed into the apparent parameters CL/F and Vd/F, so no f(depot)
    # term is applied. Apixaban dose was tested on F and not retained.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # 4. Observation and error.
    # central is in mg and vc in L, so central / vc is mg/L; multiply by 1000
    # to report ng/mL (1 mg/L = 1 ug/mL = 1000 ng/mL), matching the paper's
    # reporting units for concentrations, the 0.5 ng/mL LLOQ, the 230 ng/mL
    # bleeding-risk threshold and the ng*h/mL AUCss of Table S1.
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
