Staatz_2005_gentamicin <- function() {
  description <- "Two-compartment IV population PK model for gentamicin in adult cardiothoracic-surgery patients with unstable renal function (Staatz 2005). Clearance scales linearly with raw Cockcroft-Gault creatinine clearance centred at the population baseline median (63 mL/min); central and peripheral volumes scale linearly with body weight; intercompartmental clearance is a population constant. Operator-resolved sidecar (request-001) replaced the paper's Wahlby 2004 baseline-CrCl + change-from-baseline (BCOV+DCOV) decomposition with the simpler CrCl-only covariate form to avoid adding a new canonical baseline-CrCl column; for stable-CrCl subjects the published final-model parameters reproduce the paper's CL exactly because DCOV is zero (see vignette Errata)."
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
      notes              = "Staatz 2005 Table 1: median 72 kg (range 39-138) in the combined gentamicin cohort. Enters the structural model as a linear scalar on both V1 and V2 (Table 2 footer for the gentamicin final model: V1 = theta x WT, V2 = theta x WT, Q = theta). No reference-weight normalisation; the typical-value V1 and V2 are reported directly in L/kg.",
      source_name        = "WT"
    ),
    CRCL = list(
      description        = "Cockcroft-Gault creatinine clearance (raw, not BSA-normalized)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source column CL_Cr. Computed by the Cockcroft-Gault equation in raw mL/min (NOT BSA-normalised to mL/min/1.73 m^2); Cr_Se measurements below 60 umol/L were set to 60 umol/L per the paper's Methods (better CL estimates in patients with low creatinine production). Stored under the canonical CRCL column per inst/references/covariate-columns.md (CRCL accepts raw mL/min when the source paper does not apply BSA normalization, with the per-model description recording the assay form). Reference value 63 mL/min corresponds to the BCOV_median of the gentamicin model-building data set (Staatz 2005 Table 2 bold/final row: 'CL = theta1 x (1 + theta2 x (BCOV - 63)) + theta3 x DCOV'); the combined-data BCOV median is not explicitly restated in Table 4, so 63 is used consistent with the published final-model parameter values reported in the same Table 4 'All data' column. CRCL ranged from 9 to 169 mL/min in the gentamicin model-building cohort (Table 1).",
      source_name        = "CL_Cr"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 135L,
    n_studies      = 1L,
    age_range      = "29-83 years (model-building); 31-82 years (test)",
    age_median     = "62 years (model-building); 67 years (test)",
    weight_range   = "39-138 kg (model-building); 47-119 kg (test)",
    weight_median  = "72 kg (model-building and test)",
    sex_female_pct = 28.9,
    race_ethnicity = "Not reported (UK single-centre cardiothoracic surgery cohort)",
    disease_state  = "Adults receiving intravenous gentamicin for postoperative sepsis or related infection following cardiothoracic surgery, with unstable renal function. 79% had cardiac surgery, 4% thoracic surgery, 1% wound infection, 1% native valve endocarditis; the remainder did not have admission reason recorded.",
    dose_range     = "Intravenous gentamicin: short infusion over 10-30 min or slow bolus over 2-3 min. Median dose 120 mg (range 60-300) in the model-building set; 160 mg (range 40-300) in the test set. Therapy adjusted to achieve 1 h post-dose peak concentrations of 4-6 or 7-10 mg/L and troughs <2 mg/L.",
    regions        = "United Kingdom (Western Infirmary, Glasgow; Cardiothoracic Surgery Unit). Data collected January 1998 - September 2003 (model-building) and September 2003 - August 2004 (test).",
    renal_function = "Cockcroft-Gault creatinine clearance median 57 mL/min (range 9-169 in model-building); median 55 mL/min (range 21-171 in test). Max within-subject CL_Cr change median 15.1 mL/min (range 0-95.1, model-building) and 18.3 mL/min (1.0-85.6, test). Cr_Se < 60 umol/L (the lower limit of the reference range) was set to 60 umol/L before computing CL_Cr.",
    n_concentrations = 550L,
    notes          = paste(
      "Baseline demographics from Staatz 2005 Table 1. 96 patients in the model-building set",
      "(365 gentamicin concentrations, 1-20 samples per patient, median 2), 39 patients in the",
      "test set (185 concentrations, 1-19 samples per patient, median 4). Combined cohort N = 135",
      "is reported as the sum (the paper notes 42 patients received both gentamicin and vancomycin",
      "but does not specify the overlap between gentamicin model-building and test sets, so the",
      "subject totals are read directly from Table 1).",
      "Female% computed as 39/135 from the combined male/female counts (69+27=96 / 27+12=39).",
      "Cr_Se concentrations were typically measured every 1-3 days; missing values were imputed",
      "by linear daily interpolation (method A) for the primary analysis. Data fit using NONMEM",
      "V 1.1 with FOCE-I."
    )
  )

  ini({
    # Structural parameters (Staatz 2005 Table 4 'All data' column for the
    # gentamicin final model). Volumes are reported per kg (L/kg);
    # intercompartmental clearance Q is a population value (L/h, no weight scaling).
    lcl <- log(2.81);  label("Clearance at CRCL = 63 mL/min (L/h)")            # Table 4 'All data': theta1 = 2.81 L/h
    lvc <- log(0.251); label("Central volume per kg body weight (L/kg)")       # Table 4 'All data': V1 = 0.251 L/kg
    lvp <- log(0.156); label("Peripheral volume per kg body weight (L/kg)")    # Table 4 'All data': V2 = 0.156 L/kg
    lq  <- log(1.45);  label("Intercompartmental clearance (L/h)")             # Table 4 'All data': Q = 1.45 L/h

    # Covariate effect on CL: linear centred-deviation in raw CL_Cr (mL/min),
    # using the published final-model BCOV coefficient (theta2 = 0.0150)
    # collapsed to operate on the time-varying CRCL alone. The paper's full
    # final form is CL = theta1 x (1 + theta2 x (BCOV - 63) + theta3 x DCOV);
    # this library implements CL = theta1 x (1 + theta2 x (CRCL - 63)),
    # equivalent for stable-CrCl subjects (DCOV = 0) and recommended by the
    # operator-resolved sidecar (request-001). The theta3 = 0.0174 DCOV term
    # is documented in vignette Errata but not encoded here.
    e_crcl_cl <- 0.0150; label("Linear effect of (CRCL - 63 mL/min) on CL")    # Table 4 'All data': theta2 = 0.0150 per mL/min

    # Inter-individual variability (Staatz 2005 Table 4 'All data' column).
    # The paper retains IIV on CL and V2 only; IIV on V1 and Q was removed
    # because retaining it did not improve OFV and destabilised the fit
    # ("Removal of IIV in the central compartment volume (V1) and
    # intercompartmental clearance (Q) further stabilized the model with no
    # change in OFV"). Log-normal IIV: omega^2 = log(CV^2 + 1).
    etalcl ~ 0.07037   # log(0.27^2 + 1); 27% CV on CL (Table 4 'All data': etaCL = 27%)
    etalvp ~ 0.45506   # log(0.84^2 + 1); 84% CV on V2 (Table 4 'All data': etaV2 = 84%)

    # Residual error: combined additive + proportional (Staatz 2005 Table 4
    # 'All data' column).
    addSd  <- 0.13;  label("Additive residual error (mg/L)")                   # Table 4 'All data': Add = 0.13 mg/L
    propSd <- 0.24;  label("Proportional residual error (fraction)")           # Table 4 'All data': Prop = 24%
  })

  model({
    # Individual PK parameters. CL has a linear centred-deviation effect on
    # CRCL (raw Cockcroft-Gault mL/min) centred at 63 (the BCOV_median of the
    # gentamicin model-building data set, applied to combined-data parameters
    # consistent with the published final model when DCOV = 0). V1 and V2
    # scale linearly with body weight (typical-value reported per kg).
    cl <- exp(lcl + etalcl) * (1 + e_crcl_cl * (CRCL - 63))
    vc <- exp(lvc) * WT
    vp <- exp(lvp + etalvp) * WT
    q  <- exp(lq)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Two-compartment IV PK. Doses enter the central compartment directly
    # (short IV infusion or slow bolus per Staatz 2005 Methods).
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # Dose in mg, vc in L -> central / vc has units mg/L.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
