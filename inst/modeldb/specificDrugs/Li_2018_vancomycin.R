Li_2018_vancomycin <- function() {
  description <- "One-compartment IV-infusion population PK model for vancomycin in critically ill Chinese ICU neonates (Li 2018). CL scales allometrically with body weight (reference 2.9 kg, exponent 1.55) and as an inverse power of serum creatinine (reference 23.3 umol/L, exponent 0.337 on the SCr_ref/SCr ratio). V scales allometrically with body weight (reference 2.9 kg, exponent 1.05). IIV is on CL only; residual error is proportional."
  reference <- "Li Z, Liu Y, Jiao Z, Qiu G, Huang J, Xiao Y, Wu S, Wang C, Hu W, Sun H. Population Pharmacokinetics of Vancomycin in Chinese ICU Neonates: Initial Dosage Recommendations. Front Pharmacol. 2018;9:603. doi:10.3389/fphar.2018.00603"
  vignette <- "Li_2018_vancomycin"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight (current; time-varying as recorded in NICU charts)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Li 2018 Table 1: median 2.74 kg, mean 2.87 kg (range 1.4-5.6). Used with allometric scaling reference 2.9 kg in both CL and V structural equations (Li 2018 Table 3 footnote: CL = theta1 * (WT/2.9)^theta2 * (23.3/Scr)^theta3; V = theta4 * (WT/2.9)^theta5).",
      source_name        = "WT"
    ),
    CREAT = list(
      description        = "Serum creatinine (enzymatic method, Hitachi 7180 Automatic Analyzer)",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Li 2018 Table 1: median 28.3 umol/L, mean 23.2 umol/L (range 5.85-61.6). Reference 23.3 umol/L in the CL structural equation. Enters CL as a power on the SCr_ref/SCr ratio (Li 2018 Table 3 footnote: CL = ... * (23.3/Scr)^theta3) so higher SCr reduces CL. No effect on V.",
      source_name        = "Scr"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 80L,
    n_studies        = 1L,
    age_range        = "PNA 4-126 days (median 24); GA 25.7-41.1 weeks at birth; PMA 29-47.1 weeks",
    age_median       = "PNA 24 days; PMA 40 weeks",
    weight_range     = "1.4-5.6 kg",
    weight_median    = "2.74 kg",
    sex_female_pct   = 32.5,
    race_ethnicity   = "Chinese (single-center cohort at Shanghai Children's Hospital, Shanghai Jiao Tong University)",
    disease_state    = "Critically ill neonates in the neonatal ICU; 59% preterm; 57% with respiratory tract infections. Patients with extracorporeal membrane oxygenation or continuous renal replacement therapy were excluded.",
    dose_range       = "Intravenous vancomycin 10-15 mg/kg every 8 h (q8h) or every 12 h (q12h) administered as a 2-h infusion per local protocols.",
    regions          = "China (Shanghai)",
    renal_function   = "Serum creatinine median 28.3 umol/L (range 5.85-61.6); BUN median 4.1 mmol/L (range 0.4-28.5); albumin median 32 g/L (range 21.6-46.8).",
    n_concentrations = 165L,
    notes            = "Patient characteristics from Li 2018 Table 1 (January 2013 to December 2016 enrolment). Sampling design: peak 1 h after end of 2-h infusion and trough half an hour before the next dose, both obtained after at least four repeated doses (75 trough and 90 peak observations). Modelling done in NONMEM 7.4 with FOCE-I (eta-eps interaction). Volume-of-distribution IIV had RSE > 50% with the sparse peak/trough design and was not estimated."
  )

  ini({
    # Structural parameters (Li 2018 Table 3 final-model "Estimate"; reference
    # subject has WT = 2.9 kg and Scr = 23.3 umol/L, per the Table 3 footnote
    # final-model equations).
    lcl <- log(0.309); label("Typical CL at reference covariates (L/h)") # Li 2018 Table 3: theta1 = 0.309 (RSE 5%)
    lvc <- log(2.63);  label("Typical V at reference covariates (L)")    # Li 2018 Table 3: theta4 = 2.63  (RSE 8%)

    # Covariate effects on CL (Li 2018 Table 3 footnote final-model equation:
    #   CL = theta1 * (WT/2.9)^theta2 * (23.3/Scr)^theta3 )
    e_wt_cl    <- 1.55;  label("Allometric exponent on WT for CL (unitless)")           # Li 2018 Table 3: theta2 = 1.55  (RSE 10%)
    e_creat_cl <- 0.337; label("Power exponent on (CREAT_ref / CREAT) for CL (unitless)") # Li 2018 Table 3: theta3 = 0.337 (RSE 40%)

    # Covariate effect on V (Li 2018 Table 3 footnote final-model equation:
    #   V = theta4 * (WT/2.9)^theta5 )
    e_wt_vc    <- 1.05;  label("Allometric exponent on WT for V (unitless)")  # Li 2018 Table 3: theta5 = 1.05  (RSE 27%)

    # Inter-individual variability (Li 2018 Table 3 final-model column reports
    # IIV as %CV; for log-normal etas omega^2 = log(CV^2 + 1)). V IIV was not
    # estimated because the sparse peak/trough sampling design produced RSE
    # > 50% (Li 2018 Results "Model Building").
    etalcl ~ 0.1343 # log(0.379^2 + 1); 37.9% CV on CL (RSE 26%)

    # Proportional residual error (Li 2018 Results "Model Building": "The
    # residual unexplained variability was described best by a proportional
    # model"; Table 3 reports 37.5% CV with RSE 19%).
    propSd <- 0.375; label("Proportional residual error (fraction)") # Li 2018 Table 3: 37.5% CV (RSE 19%)
  })
  model({
    # Individual PK parameters. Encoding follows Table 3 footnote verbatim:
    # the typical-value subject (WT = 2.9 kg, CREAT = 23.3 umol/L) recovers
    # exp(lcl) and exp(lvc). The Scr factor is the inverse ratio
    # (23.3 / CREAT)^e_creat_cl so that elevated SCr reduces CL.
    cl <- exp(lcl + etalcl) *
          (WT / 2.9)^e_wt_cl *
          (23.3 / CREAT)^e_creat_cl
    vc <- exp(lvc) *
          (WT / 2.9)^e_wt_vc

    kel <- cl / vc

    d/dt(central) <- -kel * central

    # Dose in mg, vc in L -> central/vc has units mg/L.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
