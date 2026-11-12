Koomen_2023_tacrolimus <- function() {
  description <- "Joint one-compartment population PK model for plasma (Cc) and whole-blood (Cblood) tacrolimus trough concentrations in adult kidney and lung transplant recipients (Koomen 2023). Plasma is described by first-order absorption (ka and F fixed) and first-order elimination; whole blood is the haematocrit-weighted sum of the plasma concentration and an erythrocyte-bound concentration that saturates in plasma (Bmax, Kd). Bioavailability depends on daily dose (per formulation: Prograf, Advagraf, Envarsus), food intake and time after transplantation (Prograf only); apparent CL on haematocrit, eGFR and transplant type; V on haptoglobin; Bmax on haematocrit, total cholesterol, haptoglobin, vitamin B12, transplant type and concomitant mycophenolic acid and vitamin D3."
  reference <- "Koomen JV, Knobbe TJ, Zijp TR, Kremer D, Gan CT, Verschuuren EAM, Bakker SJL, Touw DJ, Colin PJ; TransplantLines Investigators. A Joint Pharmacokinetic Model for the Simultaneous Description of Plasma and Whole Blood Tacrolimus Concentrations in Kidney and Lung Transplant Recipients. Clin Pharmacokinet. 2023;62(8):1117. doi:10.1007/s40262-023-01259-x"
  vignette <- "Koomen_2023_tacrolimus"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    DOSE_TAC_MGD = list(
      description = "Total daily tacrolimus dose",
      units = "mg/day",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying (TDM-titrated). Enters the bioavailability equation (Koomen 2023 Eq. 8) as a formulation-specific power of Dose / 2.5 mg. The paper does not state whether 'Dose' is the daily or the per-administration dose; the daily dose is the reading that reproduces the paper's own reference-subject simulation (3 mg Prograf -> plasma trough 0.09 ng/mL and whole blood 4.4 ng/mL, Results 'Simulations'), with 3 mg/day given as 1.5 mg twice daily. The per-administration reading gives 0.15-0.19 ng/mL. See the vignette.",
      source_name = "Dose"
    ),
    FORM_TAC_IR = list(
      description = "Immediate-release twice-daily tacrolimus (Prograf) indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (Advagraf or Envarsus once-daily products)",
      notes = "Encodes the three oral products jointly with FORM_TAC_ENVARSUS: Prograf = (FORM_TAC_IR = 1, FORM_TAC_ENVARSUS = 0); Advagraf = (0, 0); Envarsus = (0, 1). Selects the Prograf dose exponent, and gates the fasted-state and time-after-transplantation effects on F, which Koomen 2023 Eq. 8 applies to Prograf only. The paper writes 'Prograft' throughout.",
      source_name = "drug product (Prograft = 1, Other = 0)"
    ),
    FORM_TAC_ENVARSUS = list(
      description = "Envarsus (LCP-Tac) once-daily tacrolimus indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (Prograf or Advagraf)",
      notes = "Selects the Envarsus dose exponent on F (Koomen 2023 Eq. 8 / Table 2). Advagraf is the level with both FORM_TAC_IR and FORM_TAC_ENVARSUS equal to 0.",
      source_name = "formulation"
    ),
    FED = list(
      description = "Fed-state indicator at the dose (1 = fed / non-fasted, 0 = fasted)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (fasted)",
      notes = "Koomen 2023 codes the opposite indicator ('fasted', 1 = fasted); the model uses fasted = 1 - FED. The food effect (F lower by 11% when fasted) applies to Prograf only (Eq. 8; Table 2 'Prograft, fasted').",
      source_name = "fasted (inverted)"
    ),
    POD = list(
      description = "Time after transplantation",
      units = "days",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying. Koomen 2023 Eq. 8 writes the term as 1 + exp(-0.32 * Time / 1000) for Prograf. Table 2 labels the covariate 'Time after transplantation (days)', but only a time in HOURS reproduces the paper's own month-3 versus month-6 reference-subject simulations (ESM Figs 1 and 2: plasma ~0.115 vs ~0.093 ng/mL, a 1.2-fold ratio; days give a 1.5% change). The model therefore converts POD to hours (POD * 24) before applying the published coefficient. See the vignette.",
      source_name = "Time"
    ),
    HCT = list(
      description = "Haematocrit",
      units = "%",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying (last observation carried forward in the source). Power effects on CL (reference 39%) and on Bmax (reference 39%), and the weighting factor between plasma and erythrocyte-bound concentration in the whole-blood equation (Eq. 2), where it is used as a fraction (HCT / 100).",
      source_name = "Haematocrit"
    ),
    CRCL = list(
      description = "Estimated glomerular filtration rate (CKD-EPI), BSA-normalised",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Creatinine-based CKD-EPI eGFR (Koomen 2023 Methods). Power effect on CL with reference 51 mL/min/1.73 m^2.",
      source_name = "eGFR"
    ),
    TX_LUNG = list(
      description = "Lung transplant recipient indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (kidney transplant recipient)",
      notes = "The cohort contains kidney (reference) and lung recipients only. Proportional effects on CL (-22%) and Bmax (+10%).",
      source_name = "Transplant type (0 = kidney, 1 = lung)"
    ),
    HAPTOGLOBIN = list(
      description = "Serum haptoglobin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effects on V (reference 1.4 g/L, exponent 0.13) and on Bmax (reference 1.4 g/L, exponent -0.05).",
      source_name = "Haptoglobin"
    ),
    TCHOL = list(
      description = "Total serum cholesterol",
      units = "mmol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on Bmax with reference 4.6 mmol/L.",
      source_name = "Cholesterol"
    ),
    VITAMIN_B12 = list(
      description = "Serum vitamin B12 (cobalamin)",
      units = "pmol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on Bmax with reference 320 pmol/L.",
      source_name = "Vitamin B12"
    ),
    CONMED_MPA = list(
      description = "Concomitant mycophenolic acid (mycophenolate mofetil or mycophenolate sodium) indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no mycophenolic acid)",
      notes = "Proportional effect on Bmax (+14% when used). 83.7% of the cohort used mycophenolic acid (Table 1).",
      source_name = "Mycophenolic acid (0 = not used, 1 = used)"
    ),
    CONMED_COLECALCIFEROL = list(
      description = "Concomitant vitamin D3 (colecalciferol) indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no vitamin D3)",
      notes = "Proportional effect on Bmax (+5% when used; Table 2 'Vitamin D use'). Eq. 9 is typeset as '1 x 0.05 + VitaminD', which is read as 1 + 0.05 * VitaminD per the 'Proportional' structure in Table 2.",
      source_name = "Vitamin D (0 = not used, 1 = used)"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Tested on CL and V in the first covariate step (ESM Table 1) but not retained."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      notes = "Tested on V (ESM Table 1) but not retained; the Discussion notes albumin was not significant."
    ),
    CONMED_CCB = list(
      description = "Concomitant calcium-channel blocker",
      units = "(binary)",
      type = "binary",
      notes = "Tested on CL (ESM Table 1) but not retained."
    ),
    CONMED_PPI = list(
      description = "Concomitant proton-pump inhibitor",
      units = "(binary)",
      type = "binary",
      notes = "Tested on F (ESM Table 1) but not retained."
    )
  )

  compartmentData <- list(
    depot = list(analyte = "tacrolimus", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "tacrolimus", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 1060L,
    n_studies = 1L,
    n_observations = "1973 plasma and 1961 whole-blood tacrolimus trough concentrations (plasma 1670 kidney / 303 lung; whole blood 1666 kidney / 295 lung); 2.5 samples per patient on average (range 1-6).",
    age_range = "Adults (>= 18 years); mean 54.7 (SD 13.6) years",
    weight_range = "Mean 81.0 (SD 16.3) kg",
    sex_female_pct = 39.5,
    race_ethnicity = c(Caucasian = 91.0, Asian = 2.0, Black = 1.3, Other_or_missing = 5.7),
    disease_state = "Adult kidney (n = 845) and lung (n = 215) transplant recipients from the TransplantLines Biobank and Cohort Study, sampled at outpatient visits 2-409 days after transplantation at first sample (median 47 days); haematocrit 39.0 (SD 5.5)%, eGFR 50.2 (SD 19.2) mL/min/1.73 m^2.",
    dose_range = "Oral tacrolimus titrated by whole-blood TDM (kidney trough targets 6-10 ng/mL at month 3 and 4-6 ng/mL from month 6; lung 7-10 ng/mL from month 3). Doses 2.9 (SD 1.9) mg. Prograf twice daily (74.7% of plasma samples), Advagraf once daily (18.5%), Envarsus once daily (5.4%).",
    regions = "Netherlands (University Medical Center Groningen)",
    co_medication = "Prednisolone/prednisone 97.8%, proton-pump inhibitor 86.4%, mycophenolic acid 83.7%, colecalciferol 39.3%, calcium-channel blocker 37.6% (Table 1).",
    notes = "Koomen 2023 Table 1. Samples drawn at 08:00 after an overnight fast with tacrolimus taken at 10:00 (and 22:00 for twice-daily regimens), so sampling was assumed at 10 h (twice daily) or 22 h (once daily) after the last dose. Plasma processed at 21 degrees C; the plasma concentrations and the derived plasma target ranges depend on that centrifugation temperature (Discussion)."
  )

  ini({
    # Structural parameters (Koomen 2023 Table 2 and Eqs. 5-10). CL and V are
    # apparent plasma parameters relative to the fixed F of 0.25.
    lka <- fixed(log(3.5)); label("First-order absorption rate constant ka (1/h)") # Eq. 5 and Methods: ka fixed to 3.5 1/h
    lcl <- log(257.6); label("Apparent plasma clearance CL at the reference covariates (L/h)") # Table 2 'Apparent plasma CL (L/h)' 257.6 (RSE 3.6%)
    lvc <- log(4213.2); label("Apparent plasma volume of distribution V at haptoglobin 1.4 g/L (L)") # Table 2 'Apparent plasma Vd (L)' 4213.2 (RSE 8.3%)
    lfdepot <- fixed(log(0.25)); label("Oral bioavailability F at the reference covariates (fraction)") # Eq. 8 and Methods: F fixed to 25%
    lbmax <- log(35.7); label("Maximum erythrocyte-bound concentration Bmax at the reference covariates (ng/mL)") # Table 2 'Bmax (ng/mL)' 35.7 (RSE 6.7%)
    lkd <- log(0.24); label("Plasma concentration giving half-maximal erythrocyte binding, Kd (ng/mL)") # Table 2 'C50 (ng/mL)' 0.24 (RSE 10.4%); Eq. 10 Kd

    # Covariate effects on F (Eq. 8 and Table 2)
    e_dose_tac_fdepot_prograf <- -0.64; label("Power exponent of daily dose / 2.5 mg on F, Prograf (unitless)") # Table 2 F Dose 'Prograft' power -0.64 (RSE 4.7%)
    e_dose_tac_fdepot_envarsus <- -0.74; label("Power exponent of daily dose / 2.5 mg on F, Envarsus (unitless)") # Table 2 F Dose 'Envarsus' power -0.74 (RSE 9.7%)
    e_dose_tac_fdepot_advagraf <- -0.62; label("Power exponent of daily dose / 2.5 mg on F, Advagraf (unitless)") # Table 2 F Dose 'Advagraf' power -0.62 (RSE 9.4%)
    e_pod_fdepot <- -0.32; label("Exponential coefficient of time after transplantation (per 1000 h) in the Prograf F term 1 + exp(coef * time / 1000) (unitless)") # Table 2 'Time after transplantation' exponential -0.32 (RSE 8.9%); Eq. 8
    e_fasted_fdepot <- -0.11; label("Proportional change in F when fasted, Prograf only (fraction)") # Table 2 'Prograft, fasted' proportional -0.11 (RSE 22.4%)

    # Covariate effects on CL (Eq. 6 and Table 2)
    e_hct_cl <- 0.49; label("Power exponent of haematocrit / 39% on CL (unitless)") # Table 2 CL Haematocrit power 0.49 (RSE 16.4%)
    e_crcl_cl <- 0.07; label("Power exponent of eGFR / 51 mL/min/1.73 m^2 on CL (unitless)") # Table 2 CL eGFR power 0.07 (RSE 40.2%)
    e_tx_lung_cl <- -0.22; label("Proportional change in CL for lung transplant recipients (fraction)") # Table 2 CL 'Transplantation type (lung recipient)' -0.22 (RSE 10.9%)

    # Covariate effect on V (Eq. 7 and Table 2)
    e_haptoglobin_vc <- 0.13; label("Power exponent of haptoglobin / 1.4 g/L on V (unitless)") # Table 2 V Haptoglobin power 0.13 (RSE 31.2%)

    # Covariate effects on Bmax (Eq. 9 and Table 2)
    e_hct_bmax <- -0.43; label("Power exponent of haematocrit / 39% on Bmax (unitless)") # Table 2 Bmax Haematocrit power -0.43 (RSE 14.6%)
    e_tchol_bmax <- -0.18; label("Power exponent of total cholesterol / 4.6 mmol/L on Bmax (unitless)") # Table 2 Bmax 'Total cholesterol' power -0.18 (RSE 16.4%)
    e_haptoglobin_bmax <- -0.05; label("Power exponent of haptoglobin / 1.4 g/L on Bmax (unitless)") # Table 2 Bmax Haptoglobin power -0.05 (RSE 22.6%)
    e_vitamin_b12_bmax <- 0.04; label("Power exponent of vitamin B12 / 320 pmol/L on Bmax (unitless)") # Table 2 Bmax 'Vitamin B12' power 0.04 (RSE 32.2%)
    e_tx_lung_bmax <- 0.10; label("Proportional change in Bmax for lung transplant recipients (fraction)") # Table 2 Bmax 'Transplantation type (lung recipient)' 0.10 (RSE 21.9%)
    e_conmed_mpa_bmax <- 0.14; label("Proportional change in Bmax with concomitant mycophenolic acid (fraction)") # Table 2 Bmax 'Mycophenolic acid' 0.14 (RSE 16.3%)
    e_conmed_colecalciferol_bmax <- 0.05; label("Proportional change in Bmax with concomitant vitamin D3 (fraction)") # Table 2 Bmax 'Vitamin D use' 0.05 (RSE 28.3%)

    # Between-subject variability (Table 2; variances on the log scale, with
    # CV% = sqrt(exp(omega^2) - 1)). The covariance -0.02 gives r = -0.39.
    etalfdepot + etalbmax ~ c(0.13, -0.02, 0.02) # Table 2 BSV F 0.13 (37.3%), BSV Bmax 0.02 (14.4%), covariance -0.02 (r = -0.39)

    # Residual variability (Table 2; proportional on both matrices). The table
    # gives variances 0.14 and 0.06 and CV% = sqrt(variance); the printed CV%
    # carries more digits than the rounded variance and is used here. The
    # plasma / whole-blood residual covariance (0.06, r = 0.60; NONMEM L2) is
    # not encoded; see the vignette.
    propSd <- 0.371; label("Proportional residual error, plasma Cc (fraction)") # Table 2 'Residual variability for plasma' 0.14 (CV 37.1%)
    propSd_Cblood <- 0.253; label("Proportional residual error, whole blood Cblood (fraction)") # Table 2 'Residual variability for whole blood' 0.06 (CV 25.3%)
  })

  model({
    # Formulation indicators (Prograf / Advagraf / Envarsus) and fasted state
    form_prograf <- FORM_TAC_IR
    form_envarsus <- FORM_TAC_ENVARSUS
    form_advagraf <- (1 - FORM_TAC_IR) * (1 - FORM_TAC_ENVARSUS)
    fasted <- 1 - FED

    # Eq. 8: bioavailability. Dose is the daily dose; time after transplantation
    # enters in hours (POD * 24) divided by 1000. For Advagraf and Envarsus the
    # time term is 1 + exp(0) = 2, as printed.
    e_dose_tac_fdepot <- form_prograf * e_dose_tac_fdepot_prograf +
      form_envarsus * e_dose_tac_fdepot_envarsus +
      form_advagraf * e_dose_tac_fdepot_advagraf
    tx_time_h <- POD * 24
    fdepot <- exp(lfdepot + etalfdepot) *
      (DOSE_TAC_MGD / 2.5)^e_dose_tac_fdepot *
      (1 + e_fasted_fdepot * fasted * form_prograf) *
      (1 + exp(e_pod_fdepot * tx_time_h / 1000 * form_prograf))

    # Eq. 6 (CL), Eq. 7 (V), Eq. 5 (ka)
    ka <- exp(lka)
    cl <- exp(lcl) * (HCT / 39)^e_hct_cl * (CRCL / 51)^e_crcl_cl *
      (1 + e_tx_lung_cl * TX_LUNG)
    vc <- exp(lvc) * (HAPTOGLOBIN / 1.4)^e_haptoglobin_vc
    kel <- cl / vc

    # Eq. 9 (Bmax) and Eq. 10 (Kd)
    bmax <- exp(lbmax + etalbmax) * (HCT / 39)^e_hct_bmax *
      (TCHOL / 4.6)^e_tchol_bmax *
      (HAPTOGLOBIN / 1.4)^e_haptoglobin_bmax *
      (VITAMIN_B12 / 320)^e_vitamin_b12_bmax *
      (1 + e_conmed_mpa_bmax * CONMED_MPA) *
      (1 + e_tx_lung_bmax * TX_LUNG) *
      (1 + e_conmed_colecalciferol_bmax * CONMED_COLECALCIFEROL)
    kd <- exp(lkd)

    # Eq. 1: one-compartment plasma PK with first-order absorption
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    f(depot) <- fdepot

    # Plasma concentration: mg / L -> ng/mL
    Cc <- 1000 * central / vc

    # Eqs. 2 and 4: erythrocyte-bound concentration saturating in plasma, and
    # whole blood as the haematocrit-weighted sum of plasma and bound
    # concentrations, assumed in equilibrium at sampling
    Cbound <- Cc * bmax / (kd + Cc)
    Cblood <- Cc * (1 - HCT / 100) + Cbound * HCT / 100

    Cc ~ prop(propSd)
    Cblood ~ prop(propSd_Cblood)
  })
}
