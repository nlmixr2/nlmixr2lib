Lanoiselee_2021_cefazolin <- function() {
  description <- "Two-compartment population PK model for intravenous bolus cefazolin given as antibiotic prophylaxis to adults undergoing primary total hip arthroplasty (Lanoiselee 2021, n = 100 patients, 29% obese, 484 total plasma concentrations). Elimination clearance scales with CKD-EPI estimated creatinine clearance through an estimated power function centred at 80 mL/min/1.73 m^2; no body-size descriptor (total body weight, BMI, lean body weight) was retained. Log-normal between-subject variability on CL, Vc, Q and Vp with a CL-Vc correlation; proportional residual error. Estimated in Monolix 4.3 by SAEM."
  reference <- "Lanoiselee J, Chaux R, Hodin S, Bourayou S, Gibert A, Philippot R, Molliex S, Zufferey PJ, Delavenne X, Ollier E. Population pharmacokinetic model of cefazolin in total hip arthroplasty. Sci Rep. 2021;11:19763. doi:10.1038/s41598-021-99162-7. PMCID PMC8492877. ClinicalTrials.gov NCT02252497 (PORTO study). Parameter estimates from Table 2 (Model 2, final model); covariate equation and error model from the Methods display equations."
  vignette <- "Lanoiselee_2021_cefazolin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Lanoiselee 2021 Methods ('Sample collection and drug assay'): TOTAL plasma
  # cefazolin concentrations by LC-MS/MS (Q-Exactive Plus), linear range
  # 5-500 mg/L, LLOQ 5 mg/L. No protein-binding sub-model, so both states hold
  # total cefazolin.
  compartmentData <- list(
    central = list(analyte = "cefazolin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "cefazolin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = "Creatinine clearance estimated by the CKD-EPI equation (BSA-normalized)",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Lanoiselee 2021 Methods: 'CrCLi is centred on 80 to provide an estimation of CL for the typical CrCL value in our population'; CL_i = CL_POP * (CrCL_i / 80)^theta_CrCL * exp(eta_i). Table 1 reports the covariate in mL/min/1.73 m^2 (mean 83, range 17-129; 81 patients > 60, 8 between 30 and 60, 1 < 30), which is the native CKD-EPI unit; the Table 2 and Table 3 footnotes label it 'mL/min' but the Table 3 column header and the Figure 4 caption repeat mL/min/1.73 m^2, so the BSA-normalized reading is used here. Cockcroft-Gault creatinine clearance was tested as an alternative renal-function estimator and CKD-EPI was the one retained (Results: 'interpatient variability for parameter CL was best explained by CrCL according to the CKD-EPI formula'). Adding the covariate reduced the between-subject SD on CL from 0.39 to 0.32 (an 18% decrease) and the BIC by 68.51 points.",
      source_name = "CrCL"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested in the stepwise BIC-based covariate search (Lanoiselee 2021 Methods) and not retained; cohort mean 67 years, range 24-91 (Table 1)."
    ),
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested as a body-size descriptor (Methods: 'total body weight (TBW)') and not retained; cohort mean 76 kg, range 48-123 (Table 1). The Discussion states 'weight did not affect cefazolin concentrations' despite 29% of the cohort being obese."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested as a body-size descriptor and not retained. Table 1: 71 patients < 30, 24 between 30 and 35, 5 > 35 kg/m^2."
    ),
    LBM = list(
      description = "Lean body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested as a body-size descriptor (Methods: 'lean body weight (LBW)') and not retained. The formula used to compute it is not stated."
    ),
    SEXF = list(
      description = "Sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Tested and not retained; 49 female, 51 male (Table 1)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 100L,
    n_studies = 1L,
    age_range = "24-91 years (mean 67; Table 1)",
    weight_range = "48-123 kg (mean 76; Table 1)",
    sex_female_pct = 49,
    race_ethnicity = "Not reported (single-centre French cohort)",
    disease_state = "Adults undergoing cementless primary unilateral total hip arthroplasty (lateral position, posterior approach) receiving cefazolin surgical prophylaxis. 29% obese (BMI > 30 kg/m^2); 9% renal impairment (8 with CKD-EPI CrCL 30-60, 1 with CrCL < 30 mL/min/1.73 m^2).",
    dose_range = "2000 mg intravenous bolus at anaesthesia induction (n = 96); 3000 mg (n = 1); 4000 mg (n = 3) when BMI > 35 kg/m^2 and total body weight > 100 kg. Protocol redosing of 1000 mg if surgery exceeded 4 h (Table 1, Methods).",
    regions = "France (University Hospital of Saint-Etienne, single centre)",
    renal_function = "CKD-EPI CrCL mean 83 mL/min/1.73 m^2, range 17-129 (Results, Table 1).",
    notes = "Residual blood from the PORTO randomized trial of tranexamic acid dosing in THA (NCT02252497; EudraCT 2013-000791-15), enrolled April 2014 - December 2015. Sampling at 3 min and 20 min after the bolus, at the end of surgery, and 3 h and 8 h after the bolus; 484 concentrations analysed. Mean time from injection to incision 0.85 h (range 0.18-2.23); mean surgery duration 1.16 h (range 0.49-2.48)."
  )

  ini({
    # Structural parameters -- Lanoiselee 2021 Table 2, Model 2 (final model).
    # Table 2 prints the clearance equation as CL (L/h) = theta1 * (CrCL/80)^theta2.
    lcl <- log(2.86); label("Clearance at CRCL = 80 mL/min/1.73 m^2 (L/h)") # Table 2 Model 2: theta1 = 2.86 (RSE 3.26%)
    lvc <- log(5.2); label("Central volume of distribution (L)") # Table 2 Model 2: Vc = 5.2 L (RSE 6.89%)
    lq <- log(10.9); label("Intercompartmental clearance (L/h)") # Table 2 Model 2: Q = 10.9 L/h (RSE 11.9%)
    lvp <- log(4.56); label("Peripheral volume of distribution (L)") # Table 2 Model 2: Vp = 4.56 L (RSE 3.07%)

    e_crcl_cl <- 0.79; label("Power exponent on (CRCL/80) for clearance (unitless)") # Table 2 Model 2: theta2 = 0.79 (RSE 10.6%)

    # Between-subject variability. Table 2 prints the Omega rows as 100 x the
    # Monolix log-scale SD (omega_CL = 0.32, etc.), despite the footnote calling
    # them variances. The SD reading is adjudicated by the paper's own Table 3:
    # with omega as SD the model reproduces PTA(C > 20 mg/L at 4 h) = 94.5% and
    # 99.8% at CrCL 120 and 90, whereas the variance reading gives ~80% and ~92%.
    # Variances below are the squared SDs; covariance = 0.83 * 0.32 * 0.57.
    etalcl + etalvc ~ c(0.1024, 0.151392, 0.3249) # Table 2 Model 2: Omega CL 32 (RSE 7.46%), Omega Vc 57 (RSE 9%), correlation 0.83 (RSE 5.81%)
    etalq ~ 0.4356 # Table 2 Model 2: Omega Q 66 (RSE 15%)
    etalvp ~ 0.01 # Table 2 Model 2: Omega Vp 10 (RSE 71.8%)

    # Residual error. Methods display equation: Obs = F + (a + b * F) * eps,
    # eps ~ N(0, 1); Results: 'Residual variability was best described by a
    # proportional error model' (a = 0). Monolix b is an SD on the fraction
    # scale, printed x 100 in Table 2.
    propSd <- 0.12; label("Proportional residual error (fraction)") # Table 2 Model 2: proportional residual 12 (RSE 4.9%)
  })

  model({
    # Clearance: power function of CKD-EPI CrCL centred at 80 (Methods equation).
    cl <- exp(lcl + etalcl) * (CRCL / 80)^e_crcl_cl
    vc <- exp(lvc + etalvc)
    q <- exp(lq + etalq)
    vp <- exp(lvp + etalvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
