Staatz_2005_vancomycin <- function() {
  description <- "One-compartment IV population PK model for vancomycin in adult cardiothoracic-surgery patients with unstable renal function (Staatz 2005). Clearance scales linearly with raw Cockcroft-Gault creatinine clearance centred at the population median (57 mL/min); volume of distribution scales linearly with body weight (typical-value reported per kg). Vancomycin did not benefit from the Wahlby 2004 baseline-CrCl + change-from-baseline (BCOV+DCOV) decomposition in the paper -- the simpler covariate form was retained as the final vancomycin model -- so this implementation reproduces the paper's published vancomycin final model exactly."
  reference <- paste(
    "Staatz CE, Byrne C, Thomson AH.",
    "Population pharmacokinetic modelling of gentamicin and vancomycin in",
    "patients with unstable renal function following cardiothoracic surgery.",
    "Br J Clin Pharmacol. 2006;61(2):164-176",
    "(OnlineEarly 13 December 2005).",
    "doi:10.1111/j.1365-2125.2005.02547.x",
    sep = " "
  )
  vignette <- "Staatz_2005_gentamicin_vancomycin"
  units    <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Staatz 2005 Table 1: median 74 kg (range 44-110) in the vancomycin model-building cohort; median 69 kg (range 47-112) in the test cohort. Enters the structural model as a linear scalar on V (Table 2 footer for the vancomycin final model: V = theta x WT). No reference-weight normalisation; the typical-value V is reported directly in L/kg.",
      source_name        = "WT"
    ),
    CRCL = list(
      description        = "Cockcroft-Gault creatinine clearance (raw, not BSA-normalized)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source column CL_Cr. Computed by the Cockcroft-Gault equation in raw mL/min (NOT BSA-normalised to mL/min/1.73 m^2); Cr_Se measurements below 60 umol/L were set to 60 umol/L per the paper's Methods (better CL estimates in patients with low creatinine production). Stored under the canonical CRCL column per inst/references/covariate-columns.md (CRCL accepts raw mL/min when the source paper does not apply BSA normalization, with the per-model description recording the assay form). Reference value 57 mL/min corresponds to the CL_Cr_median used in the vancomycin final model (Staatz 2005 Table 2 bold/final row: 'CL = theta1 x (1 + theta2 x (CL_Cr - 57))'). CRCL ranged from 12 to 172 mL/min in the vancomycin model-building cohort (Table 1).",
      source_name        = "CL_Cr"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 139L,
    n_studies      = 1L,
    age_range      = "17-87 years (model-building); 36-92 years (test)",
    age_median     = "66 years (model-building); 68 years (test)",
    weight_range   = "44-110 kg (model-building); 47-112 kg (test)",
    weight_median  = "74 kg (model-building); 69 kg (test)",
    sex_female_pct = 29.5,
    race_ethnicity = "Not reported (UK single-centre cardiothoracic surgery cohort)",
    disease_state  = "Adults receiving intravenous vancomycin for postoperative sepsis or related infection following cardiothoracic surgery, with unstable renal function. 78% had cardiac surgery, 6% thoracic surgery, 2% wound infection; the remainder did not have admission reason recorded.",
    dose_range     = "Intravenous vancomycin infusion over 0.2-4.42 h (median 2 h). Median dose 1000 mg (range 120-2000) in the model-building set; 1000 mg (range 500-2000) in the test set. Dosage adjusted to maintain trough concentrations of 5-10 mg/L or (from 2002 onwards) 10-15 mg/L.",
    regions        = "United Kingdom (Western Infirmary, Glasgow; Cardiothoracic Surgery Unit). Data collected January 1998 - September 2003 (model-building) and September 2003 - August 2004 (test).",
    renal_function = "Cockcroft-Gault creatinine clearance median 60 mL/min (range 12-172 in model-building); median 49 mL/min (range 18-173 in test). Max within-subject CL_Cr change median 11.5 mL/min (range 0-93.3, model-building) and 10.6 mL/min (0-84.3, test). Cr_Se < 60 umol/L (the lower limit of the reference range) was set to 60 umol/L before computing CL_Cr.",
    n_concentrations = 559L,
    notes          = paste(
      "Baseline demographics from Staatz 2005 Table 1. 102 patients in the model-building set",
      "(408 vancomycin concentrations, 1-19 samples per patient, median 3), 37 patients in the",
      "test set (151 concentrations, 1-13 samples per patient, median 4). Combined cohort N = 139",
      "is reported as the sum (the paper notes 42 patients received both gentamicin and vancomycin",
      "but does not specify the overlap between vancomycin model-building and test sets, so the",
      "subject totals are read directly from Table 1).",
      "Female% computed as 41/139 from the combined male/female counts (71+27=98 / 31+10=41).",
      "Cr_Se concentrations were typically measured every 1-3 days; missing values were imputed",
      "by linear daily interpolation (method A) for the primary analysis. Data fit using NONMEM",
      "V 1.1 with FOCE-I. A mono-exponential elimination model with combined residual error and",
      "no IIV covariance between CL and V was retained as the final structural model; the",
      "extended Wahlby 2004 BCOV+DCOV decomposition was tested but produced no clinically",
      "meaningful improvement and was not retained."
    )
  )

  ini({
    # Structural parameters (Staatz 2005 Table 4 'All data' column for the
    # vancomycin final model). Volume of distribution is reported per kg (L/kg).
    lcl <- log(2.97);  label("Clearance at CRCL = 57 mL/min (L/h)")            # Table 4 'All data': theta1 = 2.97 L/h
    lvc <- log(1.24);  label("Volume of distribution per kg body weight (L/kg)") # Table 4 'All data': V = 1.24 L/kg

    # Covariate effect on CL: linear centred-deviation in raw CL_Cr (mL/min).
    # This is the paper's final vancomycin form exactly (Table 2 vancomycin
    # final row: 'CL = theta1 x (1 + theta2 x (CL_Cr - 57))'; Table 4 'All
    # data': theta2 = 0.0205 per mL/min). No BCOV+DCOV split is used for
    # vancomycin because the extended covariate model produced no improvement
    # in IIV in CL and was driven by a single outlier patient (paper Results
    # paragraph on the extended-model analysis of the vancomycin data).
    e_crcl_cl <- 0.0205; label("Linear effect of (CRCL - 57 mL/min) on CL")    # Table 4 'All data': theta2 = 0.0205 per mL/min

    # Inter-individual variability (Staatz 2005 Table 4 'All data' column).
    # The paper retains IIV on CL and V; no IIV covariance was modelled
    # ("no covariance between interindividual variabilities in CL and volume
    # of distribution (V)"). Log-normal IIV: omega^2 = log(CV^2 + 1).
    etalcl ~ 0.07037   # log(0.27^2 + 1); 27% CV on CL (Table 4 'All data': etaCL = 27%)
    etalvc ~ 0.12188   # log(0.36^2 + 1); 36% CV on V (Table 4 'All data': etaV = 36%)
    # Note: Table 4 'All data' reports eta_V = 36% (RSE 24%); the model-building
    # column reports 33% and the test column reports 33%. The combined-data
    # value 36% is used here to match the published all-data fit.

    # Residual error: combined additive + proportional (Staatz 2005 Table 4
    # 'All data' column).
    addSd  <- 1.6;  label("Additive residual error (mg/L)")                    # Table 4 'All data': Add = 1.6 mg/L
    propSd <- 0.15; label("Proportional residual error (fraction)")            # Table 4 'All data': Prop = 15%
  })

  model({
    # Individual PK parameters. CL has a linear centred-deviation effect on
    # CRCL (raw Cockcroft-Gault mL/min) centred at 57 (the CL_Cr_median used
    # in the vancomycin final model). V scales linearly with body weight
    # (typical-value reported per kg).
    cl <- exp(lcl + etalcl) * (1 + e_crcl_cl * (CRCL - 57))
    vc <- exp(lvc + etalvc) * WT

    kel <- cl / vc

    # One-compartment IV PK. Doses enter the central compartment directly
    # (IV infusion over 0.2-4.42 h per Staatz 2005 Methods).
    d/dt(central) <- -kel * central

    # Dose in mg, vc in L -> central / vc has units mg/L.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
