deRouw_2024_pemetrexed_ckdepicr <- function() {
  description <- "Three-compartment population PK model for intravenous pemetrexed in adults with lung cancer (de Rouw 2024), with total clearance a linear function of the creatinine-based CKD-EPI eGFR de-indexed to the patient's body surface area (mL/min): CL = 1.66 + 0.0358 * CRCL (L/h). This is one of six renal-biomarker covariate models the paper compares (Table 2); the CKD-EPI creatinine-cystatin C model (deRouw_2024_pemetrexed_ckdepicrcys) fitted best. The paper reports only the clearance equation and the IIV on clearance; the volumes, intercompartmental clearances, their IIV and the residual error are carried from the de Rouw 2021 base model, which the 2024 analysis refitted to its own data without reporting the refit values."
  reference <- paste(
    "de Rouw N, Beunders R, Hartmann O, Schulte J, Boosman RJ, Derijks HJ, Burger DM, van den Heuvel MM, Hilbrands LB, Pickkers P, ter Heine R. A comparison of the renal function biomarkers serum creatinine, pro-enkephalin and cystatin C to predict clearance of pemetrexed. Cancer Chemother Pharmacol. 2024;94(6):799-806. doi:10.1007/s00280-024-04717-w.",
    "Volumes, intercompartmental clearances, their IIV and the residual error: de Rouw N, Boosman RJ, Huitema ADR, Hilbrands LB, Svensson EM, Derijks HJ, van den Heuvel MM, Burger DM, ter Heine R. Rethinking the application of pemetrexed for patients with renal impairment: a pharmacokinetic analysis. Clin Pharmacokinet. 2021;60(5):649-654. doi:10.1007/s40262-020-00972-1 (Supplementary Table S1, base model)."
  )
  vignette <- "deRouw_2024_pemetrexed"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    CRCL = list(
      description = "Creatinine-based CKD-EPI estimated glomerular filtration rate, de-indexed to the patient's body surface area",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = "RAW, NOT BSA-normalized (mL/min, not mL/min/1.73 m^2). de Rouw 2024 Methods: 'All outcomes of eGFR were expressed in ml/min and, thus, adjusted for individual body surface area (BSA)', i.e. CKD-EPI eGFR x BSA / 1.73 (the de Rouw 2021 supplement used the Du Bois BSA for the same de-indexing). Cohort median 97.0 mL/min (IQR 85.4-104.1), Table 1. Enters CL as an UNCENTRED linear term, so CL at CRCL = 0 is the 1.66 L/h intercept.",
      source_name = "CKD-EPICR"
    ),
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric scaling of the volumes (exponent 1) and intercompartmental clearances (exponent 0.75) to 70 kg, carried with those parameters from the de Rouw 2021 base model (supplement section 1.1 and its NONMEM $PK block). Clearance itself is not weight-scaled: the 2024 Table 2 equation has no weight term. Cohort median 77.1 kg (IQR 63.7-84.8), Table 1.",
      source_name = "WT"
    )
  )

  compartmentData <- list(
    central = list(analyte = "pemetrexed", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "pemetrexed", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "pemetrexed", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 66,
    n_studies = 2,
    age_median = "65 years (IQR 59-71)",
    weight_median = "77.1 kg (IQR 63.7-84.8)",
    sex_female_pct = 50,
    disease_state = "lung cancer treated with pemetrexed",
    dose_range = "intravenous pemetrexed; the dose is not restated in de Rouw 2024",
    renal_function = "CKD-EPI creatinine eGFR median 97.0 mL/min (IQR 85.4-104.1); five patients below 60 mL/min. CKD-EPI cystatin C 75.1 (59.8-92.2), CKD-EPI creatinine-cystatin C 89.0 (73.3-99.5), eGFR-PENK-creatinine 112.5 (98.9-122.6) mL/min; plasma pro-enkephalin 37.3 (30.8-46.1) pmol/L. All Table 1.",
    notes = "Patients from two multicentre clinical trials (NCT03655821, NCT03656549) with paired creatinine, cystatin C and pro-enkephalin at baseline; 378 pemetrexed plasma concentrations, at least 4 per patient. 19 patients (28.4%) were sarcopenic by fat-free mass index. Baseline characteristics are de Rouw 2024 Table 1."
  )

  ini({
    # Clearance. de Rouw 2024 Table 2, row 'CKD-EPICR (mL/min)':
    #   Clpred = 1.66 + (0.0358 . CKD-EPISCR)   (L/h; eGFR in mL/min)
    # The intercept is the clearance at zero eGFR, hence lcl_nonren; the
    # slope is the renal arm per mL/min of eGFR.
    lcl_nonren <- log(1.66); label("Clearance intercept at zero eGFR (L/h)") # Table 2 CKD-EPICR row: intercept 1.66
    e_crcl_cl_renal <- 0.0358; label("Slope of clearance on CKD-EPI creatinine eGFR (L/h per mL/min)") # Table 2 CKD-EPICR row: slope 0.0358

    # Volumes and intercompartmental clearances: de Rouw 2021 supplement
    # Table S1, 'Base' column, the model the 2024 analysis refitted
    # ('A previously developed three-compartment pharmacokinetic model
    # without covariates was fitted to the data [5]'). The 2024 refit values
    # are not reported. Values are for a 70 kg patient.
    lvc <- log(6.70); label("Central volume of distribution V1 for 70 kg (L)") # de Rouw 2021 Table S1 Base: V1 = 6.70 L
    lq <- log(6.56); label("Intercompartmental clearance Q1 for 70 kg (L/h)") # de Rouw 2021 Table S1 Base: Q1 = 6.56, printed as mL/min, read as L/h
    lvp <- log(8.01); label("First peripheral volume V2 for 70 kg (L)") # de Rouw 2021 Table S1 Base: V2 = 8.01 L
    lq2 <- log(0.04); label("Intercompartmental clearance Q2 for 70 kg (L/h)") # de Rouw 2021 Table S1 Base: Q2 = 0.04, printed as mL/min, read as L/h
    lvp2 <- log(1.25); label("Second peripheral volume V3 for 70 kg (L)") # de Rouw 2021 Table S1 Base: V3 = 1.25 L

    e_wt_vc <- fixed(1); label("Allometric exponent of WT/70 on V1, V2 and V3 (unitless)") # de Rouw 2021 supplement 1.1: exponent 1 for all volumes
    e_wt_q <- fixed(0.75); label("Allometric exponent of WT/70 on Q1 and Q2 (unitless)") # de Rouw 2021 supplement 1.1: exponent 0.75 for flows

    # IIV. de Rouw reports IIV as sqrt(omega^2) x 100 (the 2021 supplement
    # pairs 'IIV CL 21.0%' with $OMEGA 0.0442 = 0.210^2), so each variance is
    # (percent / 100)^2. IIV on CL: de Rouw 2024 Table 3 for this model.
    etalcl ~ 0.026244 # de Rouw 2024 Table 3 'CKD-EPICR' IIV 16.2%: 0.162^2
    # V1/V2/V3 block: de Rouw 2021 Table S1 Base column, IIV 33.6 / 31.2 /
    # 48.1 % and correlations V1-V2 32.5 %, V1-V3 63.1 %, V2-V3 61.9 %;
    # covariance = correlation x SD1 x SD2.
    etalvc + etalvp + etalvp2 ~ c(
      0.112896,
      0.034070, 0.097344,
      0.101979, 0.092894, 0.231361
    )

    propSd <- 0.235; label("Proportional residual error (fraction)") # de Rouw 2021 Table S1 Base: residual prop. error 23.5%
  })

  model({
    # Clearance as printed in de Rouw 2024 Table 2, with log-normal IIV on
    # the total (de Rouw 2021 NONMEM code: CL = TVCL * EXP(ETA(1))).
    cl <- (exp(lcl_nonren) + e_crcl_cl_renal * CRCL) * exp(etalcl)

    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    vp <- exp(lvp + etalvp) * (WT / 70)^e_wt_vc
    vp2 <- exp(lvp2 + etalvp2) * (WT / 70)^e_wt_vc
    q <- exp(lq) * (WT / 70)^e_wt_q
    q2 <- exp(lq2) * (WT / 70)^e_wt_q

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1 - k13 * central + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
