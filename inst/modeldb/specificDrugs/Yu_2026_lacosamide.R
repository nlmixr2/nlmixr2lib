Yu_2026_lacosamide <- function() {
  description <- "One-compartment population PK model for oral lacosamide in adult patients with epilepsy, with concomitant carbamazepine, sex and creatinine clearance on apparent clearance"
  reference <- paste(
    "Yu L, Mao F, Chen S, Yu K, Hu Y, Chen J, Hu W, Yu Z, Dai H (2026).",
    "Development and validation of a population pharmacokinetic model for",
    "lacosamide in adult patients with epilepsy to inform precision dosing.",
    "BMC Pharmacology and Toxicology.",
    "doi:10.1186/s40360-026-01114-2. PMCID: PMC13067655.",
    sep = " "
  )
  vignette <- "Yu_2026_lacosamide"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Enters only the apparent volume of distribution, which Yu 2026 fixed at",
        "0.6 L/kg (Table 2), so V/F = 0.6 * WT and the weight exponent on V/F is",
        "structurally 1 rather than an estimated allometric exponent. Body weight",
        "was screened as a covariate on CL/F (Results 'PPK model development';",
        "Table S1 tests SEX, AGE, CRCL and carbamazepine) but was NOT retained on",
        "CL/F in the final model. Baseline value; median 65 kg, range 36.0-130 kg",
        "(Table 1).",
        sep = " "
      ),
      source_name        = "WT"
    ),
    SEXF = list(
      description        = "Female sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = paste(
        "Yu 2026 Results (equation gloss following the CL/F equation) states",
        "'SEX is 1 for females and 0 for males', so the paper's SEX column maps",
        "directly onto the canonical SEXF with no value inversion. Females carry",
        "0.875-fold the apparent clearance of males, i.e. about 12.5 percent lower",
        "CL/F and correspondingly higher trough concentrations, consistent with the",
        "Discussion citation of Markoula et al. Cohort was 101 of 180 female",
        "(Table 1).",
        sep = " "
      ),
      source_name        = "SEX"
    ),
    CRCL = list(
      description        = "Creatinine clearance (Cockcroft-Gault)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "RAW Cockcroft-Gault creatinine clearance in mL/min, NOT BSA-normalized to",
        "mL/min/1.73 m^2 (Methods 'Patient inclusion': 'Creatinine clearance (CRCL)",
        "was calculated according to the Cockcroft-Gault equation'). Supply raw",
        "mL/min values, as in Delattre_2010_amikacin.R and Chen_2023_nemonoxacin.R;",
        "supplying a BSA-normalized value would silently rescale the renal term.",
        "Power-model reference is 119 mL/min, taken from the printed CL/F equation",
        "(Results; Table 2 footnote 'CL=THETA(1)*THETA(4)**CBP*THETA(5)**SEX*",
        "(CRCL/119)**THETA(6)'). Note this is not identical to the baseline median",
        "of 115 mL/min in Table 1: Table 1 footnote states 'All data were presented",
        "as baseline value' whereas the model was centred on the median across the",
        "294 observation records, so CRCL is treated as time-varying at the",
        "observation level. Baseline range 25.2-308 mL/min (Table 1).",
        sep = " "
      ),
      source_name        = "CLCR"
    ),
    CONMED_CBZ = list(
      description        = "Concomitant carbamazepine indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant carbamazepine)",
      notes              = paste(
        "1 = subject is receiving concomitant carbamazepine, an enzyme inducer that",
        "raises lacosamide apparent clearance 1.48-fold (Table 2), 0 = otherwise.",
        "28 of 180 subjects (15.6 percent) received carbamazepine (Table 1).",
        "Applied as a power form on the binary indicator, 1.48^CONMED_CBZ, matching",
        "the NONMEM control-stream line in the Table 2 footnote and the",
        "Hashimoto_1994_zonisamide.R precedent; for a 0/1 indicator this is",
        "numerically identical to a 1 + 0.48 * CONMED_CBZ fractional form. Yu 2026",
        "screened the individual concomitant antiseizure medications separately",
        "(Table 1 lists ten) and retained only carbamazepine, so the drug-specific",
        "canonical is used rather than the class-level CONMED_EIAED. The",
        "carbamazepine DOSE was not collected (Discussion, study limitations), so",
        "the effect is an all-or-none indicator only. Time-varying if carbamazepine",
        "starts or stops within the observation window.",
        sep = " "
      ),
      source_name        = "CBZ"
    )
  )

  compartmentData <- list(
    depot   = list(analyte = "lacosamide", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "lacosamide", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species          = "human",
    n_subjects       = 180,
    n_observations   = 294,
    n_studies        = 1,
    age_range        = "18-82.8 years",
    age_median       = "32.9 years (IQR 24.9-51.2)",
    weight_range     = "36.0-130 kg",
    weight_median    = "65 kg (IQR 55.3-75.0)",
    sex_female_pct   = 56.1,
    disease_state    = "epilepsy",
    dose_range       = "50-225 mg/day oral (median 100 mg/day, IQR 100-150)",
    renal_function   = "creatinine clearance median 115 mL/min (IQR 95.5-137), range 25.2-308 mL/min (Cockcroft-Gault)",
    co_medication    = paste(
      "Concomitant antiseizure medications (Table 1, percent of 180 subjects):",
      "levetiracetam 44.4, oxcarbazepine 28.3, valproic acid 26.7,",
      "lamotrigine 23.3, carbamazepine 15.6, perampanel 9.4, clonazepam 5.0,",
      "topiramate 3.9, phenobarbital 1.7, brivaracetam 0.6. Only carbamazepine",
      "was retained as a covariate in the final model.",
      sep = " "
    ),
    regions          = "China (single centre: Second Affiliated Hospital, Zhejiang University School of Medicine, Hangzhou)",
    notes            = paste(
      "Retrospective therapeutic drug monitoring cohort of inpatients treated",
      "April 2022 to February 2024 (Methods 'Patient inclusion'); baseline",
      "demographics in Table 1. ALL 294 records were trough concentrations",
      "(Results 'Patient inclusion and characteristics'), which is why Ka and V",
      "were fixed from the literature rather than estimated and why no",
      "interindividual variability on V could be identified (Discussion, study",
      "limitations). Sex percentage computed as 101/180 = 56.1 percent from the",
      "Table 1 counts; Table 1 prints 56.6 percent, which is inconsistent with its",
      "own male count of 79 (43.9 percent) summing to 180.",
      sep = " "
    )
  )

  ini({
    # Structural parameters. Only CL/F was estimated; Ka and V/F were fixed from
    # the literature because the dataset contains trough concentrations only
    # (Methods 'Base model': "owing to the concentration distribution, Ka was
    # fixed at 6.47 h-1, and V was fixed at 0.6 L/kg").
    lka <- fixed(log(6.47)); label("Absorption rate constant, literature value (1/h)")  # Table 2: Ka = 6.47 h-1 FIXED
    lcl <- log(1.86); label("Apparent clearance CL/F in males, no concomitant carbamazepine, at CRCL 119 mL/min (L/h)")  # Table 2: CL/F = 1.86 L/h (RSE 6.13%)
    lvc <- fixed(log(0.6)); label("Apparent volume of distribution V/F per kg body weight, literature value (L/kg)")  # Table 2: V = 0.6 L/kg FIXED

    # Covariate effects on CL/F. Table 2 footnote gives the NONMEM control-stream
    # form: CL = THETA(1)*THETA(4)**CBP*THETA(5)**SEX*(CRCL/119)**THETA(6).
    e_conmed_cbz_cl <- 1.48;  label("Multiplicative factor on CL/F when CONMED_CBZ = 1 (unitless)")  # Table 2: THETA(4) = 1.48 (RSE 4.71%, 95% CI 1.34-1.62)
    e_sexf_cl       <- 0.875; label("Multiplicative factor on CL/F when SEXF = 1 (unitless)")        # Table 2: THETA(5) = 0.875 (RSE 4.31%, 95% CI 0.801-0.949)
    e_crcl_cl       <- 0.311; label("Power exponent of (CRCL/119) on CL/F (unitless)")               # Table 2: THETA(6) = 0.311 (RSE 24.4%, 95% CI 0.162-0.460)

    # Interindividual variability. Table 2 reports the raw NONMEM $OMEGA element
    # (a VARIANCE on the exponential-IIV log scale), not a standard deviation:
    # 0.041 with SE 0.0133 gives RSE 32.4%, matching the printed RSE, and NONMEM
    # reports $OMEGA on the variance scale. Corresponds to
    # sqrt(exp(0.041) - 1) = 20.5% CV. No IIV on V/F could be estimated because
    # only trough concentrations were available (Discussion, study limitations).
    etalcl ~ 0.041  # Table 2: omega_CL = 0.041 (RSE 32.4%, shrinkage 21.4%); bootstrap 0.0392 (Table 3)

    # Residual error. Table 2 labels the row "sigma (proportional)" and reports
    # the raw NONMEM $SIGMA element, i.e. a VARIANCE; nlmixr2's prop() takes a
    # standard deviation, so propSd = sqrt(0.066) = 0.257 (25.7% proportional
    # error). Same convention as Luu_2017_nusinersen.R.
    propSd <- sqrt(0.066); label("Proportional residual SD on Cc (fraction; sqrt(NONMEM $SIGMA = 0.066))")  # Table 2: sigma = 0.066 (RSE 26.8%, shrinkage 18%); bootstrap 0.0662 (Table 3)
  })

  model({
    # Individual parameters. CL/F carries all three retained covariates; the
    # power form on the binary indicators reproduces the Table 2 footnote
    # control-stream line exactly.
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) *
      e_conmed_cbz_cl^CONMED_CBZ *
      e_sexf_cl^SEXF *
      (CRCL / 119)^e_crcl_cl
    # V/F was fixed per kilogram of body weight, so the weight exponent is
    # structurally 1 and is not an estimated allometric term.
    vc <- exp(lvc) * WT

    kel <- cl / vc

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
