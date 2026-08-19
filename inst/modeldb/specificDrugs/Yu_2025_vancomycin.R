Yu_2025_vancomycin <- function() {
  description <- "One-compartment IV population PK model for vancomycin in Chinese pediatric patients (birth to 15 years) receiving intermittent intravenous infusions and monitored by steady-state trough concentrations (Yu 2025). Clearance scales allometrically with body weight (fixed exponent 0.75, reference 70 kg) and as a power function of Schwartz-estimated GFR (estimated exponent 0.812, reference 173.42 mL/min/1.73 m^2); volume of distribution scales linearly with body weight (fixed exponent 1, reference 70 kg). Interindividual variability was estimable only on clearance; the proportional residual error was fixed at 30%."
  reference <- "Yu B, Mei K, Zhan D, Tang Q, Cai H, Zhang R. Establishment of a Vancomycin Population Pharmacokinetic Model for Pediatric Patients Based on the Non-Linear Mixed-Effects Model. Drugs R D. 2025;25:309-320. doi:10.1007/s40268-025-00523-8"
  vignette <- "Yu_2025_vancomycin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Vancomycin was given as a 60-minute intravenous
  # intermittent infusion (Yu 2025 Sect. 2.2), so the dose enters `central`
  # directly and there is no depot state. The specimen is verified: Sect. 2.2
  # states that venous blood was collected into EDTA tubes, centrifuged, and
  # the drug concentration determined in the PLASMA supernatant by enzyme
  # amplification immunoassay.
  compartmentData <- list(
    central = list(analyte = "vancomycin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Actual body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Yu 2025 Table 1: mean 21.49 kg, median 24 kg (range 1.25-74). Enters both CL and V allometrically with a 70 kg reference and exponents fixed at 0.75 and 1.0 respectively (Yu 2025 Eqs. 8, 9, 12, 13; Sect. 2.4.1 states the exponents were fixed). Note that the 70 kg reference is far outside this cohort's weight range, so the typical values in Table 3 (CL 8.22 L/h, V 113 L) describe a hypothetical 70 kg subject rather than a typical study child; the paper makes this explicit in Sect. 4 by dividing 8.22 L/h by 70 kg to obtain 0.117 L/h/kg for comparison with the international literature.",
      source_name        = "WT"
    ),
    CRCL = list(
      description        = "Estimated glomerular filtration rate from the height-based Schwartz equation, BSA-normalized: eGFR (mL/min/1.73 m^2) = k * height (cm) / serum creatinine (mg/dL)",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Yu 2025 Eq. 2 with k = 0.33 for preterm infants younger than 1 year, 0.45 for infants younger than 1 year, 0.55 for children 1-12 years, 0.55 for girls older than 12 years, and 0.7 for boys older than 12 years; serum creatinine was converted from umol/L to mg/dL by Eq. 1 (1 mg/dL = 88.4 umol/L). Yu 2025 Table 1: mean 163.83, median 173.42 mL/min/1.73 m^2 (range 31.9-260.82). The reference value 173.42 in Eq. 12 is the cohort median, per Sect. 2.4.2 ('continuous covariates were standardized by the median'). Stored under canonical CRCL, which covers BSA-normalized creatinine-based GFR estimates; the assay form here is the height-based (original, not bedside) Schwartz estimate. Renal strata used in the Sect. 2.6 Monte Carlo simulation: normal >= 90, mild 60-89, moderate 30-59 mL/min/1.73 m^2; children with eGFR <= 30 were excluded because they usually received renal replacement therapy.",
      source_name        = "eGFR"
    )
  )

  # Screened in the Yu 2025 stepwise covariate search (Table 2) but NOT retained
  # in the final model, so they are documentation only and are not referenced in
  # model(). eGFR was the sole covariate surviving backward elimination: it gave
  # the largest forward-step drop (dOFV -52.62) and its removal raised the OFV by
  # 52.62 (Table 2 model 11, p < 0.01). Age, ALT and BUN each passed the forward
  # alpha = 0.05 screen on their own but added nothing once eGFR was in the model
  # (dOFV -1.15, -0.38 and -1.15 respectively); sex failed even the forward screen.
  covariatesDataExcluded <- list(
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Yu 2025 Table 1: mean 5.17 years, median 4 (range 0.0027-15). Table 2 model 2: an exponential effect on CL gave dOFV -21.28 (p < 0.05) in the univariate forward step, the second-largest single-covariate drop after eGFR, but model 7 (CL-eGFR-Age) added only -1.15 on top of eGFR (p > 0.05). Not retained.",
      source_name        = "Age"
    ),
    ALT = list(
      description        = "Alanine aminotransferase",
      units              = "U/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Yu 2025 Table 1: mean 52.13 U/L, median 21.7 (range 5-805). Table 2 model 3: dOFV -4.60 (p < 0.05) univariately; model 8 (CL-eGFR-ALT) added -0.38 (p > 0.05). Not retained.",
      source_name        = "ALT"
    ),
    BUN = list(
      description        = "Blood urea nitrogen",
      units              = "mmol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Yu 2025 Table 1: mean 4.66 mmol/L, median 4.2 (range 0.3-25.4). Table 2 model 4: dOFV -5.01 (p < 0.05) univariately; model 9 (CL-eGFR-BUN) added -1.15 (p > 0.05). Not retained.",
      source_name        = "BUN"
    ),
    SEXF = list(
      description        = "Female sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "male",
      notes              = "Yu 2025 Table 1: 41 of 100 patients (41%) female. The only categorical covariate screened; Sect. 2.4.2 Eq. 11 gives the linear form used for binary covariates. Table 2 model 5: dOFV -1.67 (p > 0.05), failing even the forward alpha = 0.05 screen; model 10 (CL-eGFR-Sex) added -1.35 (p > 0.05). Not retained. Sex nevertheless enters the model indirectly through the Schwartz k coefficient for children older than 12 years (0.55 for girls, 0.7 for boys).",
      source_name        = "Sex"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 100L,
    n_studies        = 1L,
    n_sites          = 1L,
    n_concentrations = 124L,
    age_range        = "0.0027-15 years (eligibility was younger than 18 years)",
    age_median       = "4 years (mean 5.17)",
    weight_range     = "1.25-74 kg",
    weight_median    = "24 kg (mean 21.49)",
    height_median    = "106.5 cm (mean 107.98; range 37-178)",
    sex_female_pct   = 41,
    race_ethnicity   = "Not reported; single-center Chinese cohort",
    disease_state    = "Pediatric inpatients with confirmed or suspected Gram-positive infection receiving intravenous vancomycin. Excluded: renal replacement therapy, extracorporeal membrane oxygenation, undetectable vancomycin concentrations, or incomplete clinical data.",
    dose_range       = "Clinician-chosen empirical regimens; vancomycin 0.5 g/vial diluted in 0.9% saline or 5% glucose and given as a 60-minute intravenous intermittent infusion. The Sect. 2.6 Monte Carlo simulation used the Chinese 2020 guideline dose of 60 mg/kg/day as 20 mg/kg every 8 h or 15 mg/kg every 6 h.",
    regions          = "China (Anhui Provincial Children's Hospital, Hefei, Anhui)",
    renal_function   = "Schwartz eGFR median 173.42 mL/min/1.73 m^2 (mean 163.83, range 31.9-260.82); serum creatinine median 30.1 umol/L (mean 51.24, range 12.6-997.2)",
    notes            = "Single-center retrospective therapeutic-drug-monitoring study, September 2021 to November 2023 (Yu 2025 Sect. 2.1). Sampling was sparse and dominated by steady-state troughs drawn 30 minutes before the fourth dose per the Chinese Expert Consensus on Therapeutic Drug Monitoring in Children, giving only 124 concentrations from 100 patients (median slightly above one sample per subject). Plasma vancomycin was measured by enzyme amplification immunoassay on a Siemens Viva-ProE analyzer; the lower limit of quantification was 2 mg/L and measurable samples below it were retained in the modeling. Model fit in NONMEM 7.4.1 with FOCE-I, PsN 5.1.2 for covariate screening and validation; stepwise covariate selection by forward inclusion at alpha = 0.05 (dOFV 3.84) then backward elimination at alpha = 0.01 (dOFV 6.63). Validation was internal only (goodness of fit, 1000-replicate bootstrap with 997 successful runs, and a 200-replicate VPC); no external dataset was used. The paper reports that the VPC 95th percentile was poorly reproduced and recommends against using the model to predict concentrations above 20 mg/L."
  )

  ini({
    # Structural typical values (Yu 2025 Table 3). Both are reported at the
    # 70 kg allometric reference of Eqs. 12 and 13, not at the cohort median
    # weight; Sect. 4 confirms this by computing 8.22 / 70 = 0.117 L/h/kg.
    lcl <- log(8.22); label("Clearance at WT=70 kg and eGFR=173.42 mL/min/1.73 m^2 (L/h)")  # Yu 2025 Table 3: theta_CL 8.22 L/h (RSE 9%, bootstrap median 8.271, 2.5-97.5th 6.92-9.88)
    lvc <- log(113);  label("Volume of distribution at WT=70 kg (L)")                        # Yu 2025 Table 3: theta_V 113 L (RSE 18%, bootstrap median 113.082, 2.5-97.5th 81.19-172.72)

    # Allometric exponents. Yu 2025 Sect. 2.4.1: "Clearance and V were scaled
    # using fixed exponents of 0.75 and 1.0, respectively", and neither appears
    # in the Table 3 estimate list, so both are fixed structural assumptions.
    e_wt_cl <- fixed(0.75); label("Allometric exponent on (WT/70) for CL (unitless)")  # Yu 2025 Eqs. 8 and 12: (WT/70)^0.75
    e_wt_vc <- fixed(1);    label("Allometric exponent on (WT/70) for V (unitless)")   # Yu 2025 Eqs. 9 and 13: (WT/70)^1

    # Renal-function effect on clearance. Estimated (Table 3 reports an RSE),
    # so it is not wrapped in fixed().
    e_crcl_cl <- 0.812; label("Power exponent on (CRCL/173.42) for CL (unitless)")  # Yu 2025 Table 3: theta_GFR 0.812 (RSE 15%, bootstrap median 0.8162, 2.5-97.5th 0.55-1.06)

    # Interindividual variability. Yu 2025 Eq. 3 gives the exponential model
    # P_i = theta_P * exp(eta_i) with eta_i ~ N(0, omega^2), and Table 3 reports
    # the row as "Inter-individual variability omega_CL (%)" = 56.5. The row
    # names omega itself expressed as a percentage, so omega = 0.565 and the
    # variance is 0.565^2 = 0.319225. See the vignette Errata for the discarded
    # alternative reading (omega_CL as %CV = sqrt(exp(omega^2) - 1) * 100, which
    # would give omega = 0.5265) and for the RSE argument that rules out the
    # value being a variance.
    etalcl ~ 0.319225  # 0.565^2; Yu 2025 Table 3: omega_CL 56.5% (RSE 10%, bootstrap median 55.5%, 2.5-97.5th 44.7-68.9%)

    # Residual variability. Yu 2025 Sect. 3.2.1 selected the proportional model
    # (Eq. 5, Obs = Pred * (1 + eps)) over additive, combined and exponential
    # alternatives. Sect. 4 states the residual magnitude could not be estimated
    # stably from 124 sparse trough samples, so after a sensitivity analysis over
    # the candidate values 0.2, 0.3 and 0.4 it was fixed at 0.3.
    propSd <- fixed(0.3); label("Proportional residual SD (fraction of prediction)")  # Yu 2025 Table 3: residual unexplained variability 30% (fixed)
  })
  model({
    # Yu 2025 Eqs. 12 and 13:
    #   CL = 8.22 * (WT/70)^0.75 * (GFR/173.42)^0.812
    #   V  = 113  * (WT/70)
    # Interindividual variability enters CL only: Sect. 3.2.1 states that the
    # estimated interindividual variability in V was too small and was dropped
    # to stabilize the fit, so there is no etalvc.
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl * (CRCL / 173.42)^e_crcl_cl
    vc <- exp(lvc)          * (WT / 70)^e_wt_vc

    kel <- cl / vc

    d/dt(central) <- -kel * central

    # Dose in mg, volume in L, so central/vc is mg/L, the unit of the Yu 2025
    # observed concentrations (Table 1: vancomycin concentration in mg/L).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
