Wi_2017_teicoplanin <- function() {
  description <- "Two-compartment IV bolus population PK model for teicoplanin in adult patients receiving venoarterial extracorporeal membrane oxygenation (VA-ECMO) for cardiogenic shock, with binary within-subject ECMO indicators on the central volume of distribution (V1) and inter-compartmental clearance (Q) and a binary CRRT indicator on the peripheral volume of distribution (V2) (Wi 2017)"
  reference   <- "Wi J, Noh H, Min KL, Yang S, Jin BH, Hahn J, Bae SK, Kim J, Park MS, Choi D, Chang MJ. Population pharmacokinetics and dose optimization of teicoplanin during venoarterial extracorporeal membrane oxygenation. Antimicrob Agents Chemother. 2017;61(9):e01015-17. doi:10.1128/AAC.01015-17"
  vignette    <- "Wi_2017_teicoplanin"
  units       <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    ECMO_STATUS = list(
      description        = "Within-subject binary indicator for active venoarterial extracorporeal membrane oxygenation (VA-ECMO) cardiopulmonary support during the modeled record",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = "Source column ECMO. 1 = record sampled while the subject was on VA-ECMO; 0 = record sampled after the subject was weaned off VA-ECMO. Naturally time-varying: 4 of 10 study subjects survived ECMO and provided paired post-weaning samples that serve as their own controls (Wi 2017 Results, 'Patient characteristics'). Stored under the canonical ECMO_STATUS column per inst/references/covariate-columns.md. The cohort is exclusively VA-ECMO; the column does not distinguish VA from VV.",
      source_name        = "ECMO"
    ),
    CRRT_STATUS = list(
      description        = "Subject-level binary indicator for concomitant continuous renal replacement therapy (continuous venovenous hemodiafiltration via Prismaflex) during the teicoplanin sampling period",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = "Source column CRRT. 1 = subject was concomitantly receiving continuous venovenous hemodiafiltration via Prismaflex; 0 = no CRRT (Wi 2017 Materials and Methods, 'Study procedures'). Treated as time-fixed at the subject level: 5 of 10 study patients received CRRT throughout the sampling window. Stored under the canonical CRRT_STATUS column per inst/references/covariate-columns.md (distinct from HEMODIAL which is intermittent-hemodialysis-only).",
      source_name        = "CRRT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 10L,
    n_studies      = 1L,
    age_range      = "19-77 years (Wi 2017 Results: 'median age was 62.5 years (range, 19 to 77 years)')",
    age_median     = "62.5 years",
    weight_range   = "41-87 kg (Wi 2017 Results: 'median body weight was 67.5 kg (range, 41 to 87 kg)')",
    weight_median  = "67.5 kg",
    sex_female_pct = 30,
    race_ethnicity = "Not reported (single-centre Korean tertiary referral cardiac ICU, presumed predominantly Korean)",
    disease_state  = "Critically ill adults with refractory cardiogenic shock requiring VA-ECMO support. Indications for VA-ECMO: acute myocardial infarction (n = 7), myocarditis (n = 2), valvular heart disease (n = 1). All received teicoplanin for infection prophylaxis or treatment. Five subjects received concomitant CRRT via Prismaflex hemodiafiltration.",
    dose_range     = "Teicoplanin IV bolus. Per-protocol regimen during the study: loading dose 400 mg q12h x 3 doses then maintenance 400 mg q24h (study regimen A in the dosing-optimization simulations). The published dosing simulations (Wi 2017 Table 3) compared eight regimens with loading doses 400-1,200 mg q12h x 3 and maintenance doses 400-1,000 mg q24h; the authors recommend higher doses (regimen B/D for mild-to-moderate infections, regimen F/H for severe infections) than the standard regimen A to achieve >50% PTA at 72 h.",
    regions        = "Korea (Severance Hospital, Yonsei University Health System, Seoul; cardiac ICU)",
    ecmo_circuit   = "VA-ECMO via centrifugal pump (Capiox SP-101, Terumo) and X-coated circuit (Capiox EBS, Terumo); peripheral femoral vein-femoral artery cannulation. Mean ECMO duration 6.82 days (range 1.88-12.1); mean blood flow rate 2,256 mL/min (range 1,650-2,520).",
    crrt_modality  = "Continuous venovenous hemodiafiltration via Prismaflex (Gambro Inc., Meyzieu, France) when applied.",
    renal_function = "Serum creatinine 0.6-3.5 mg/dL across the cohort (Wi 2017 Table 1). Cockcroft-Gault CrCL was tested as a covariate on CL and not retained.",
    screened_covariates = "Tested-but-not-retained covariates (Wi 2017 Methods, 'Population pharmacokinetic analysis'): age, sex, body weight, serum albumin, serum creatinine, blood urea nitrogen, Cockcroft-Gault creatinine clearance, ECMO blood flow rate. Only the within-subject ECMO indicator (on V1 and Q) and the binary CRRT indicator (on V2) were retained in the final model.",
    notes          = "Baseline demographics per Wi 2017 Table 1 (n = 10 adults). Study registered as ClinicalTrials.gov NCT02581280; IRB approval 4-2014-0919 (Yonsei University). Sampling: 8 timepoints per subject per crossover period (pre-dose, 5 min, 1, 2, 3, 6, 12, 24 h after teicoplanin administration on day 2 of ECMO and on day 2 after ECMO weaning when applicable); 99 total concentrations. Assay: HPLC-MS/MS (Shimadzu LCMS-8050); LLOQ 2.0 mg/L, linear range 2-150 mg/L, CV <15% at QC samples (2, 6, 12, 120 mg/L)."
  )

  ini({
    # Structural fixed-effect parameters from Wi 2017 Table 2 (Population estimate column;
    # final population PK model, n = 10 with within-subject ECMO crossover).
    lcl       <- log(0.95); label("Clearance (L/h)")                   # Wi 2017 Table 2: CL = 0.95 L/h (RSE 29.6%)
    lvc       <- log(15.7); label("Central volume V1 (L)")             # Wi 2017 Table 2: V1 = 15.7 L (RSE 17.3%)
    lq        <- log(5.57); label("Intercompartmental clearance Q (L/h)") # Wi 2017 Table 2: Q = 5.57 L/h (RSE 43.1%)
    lvp       <- log(71.7); label("Peripheral volume V2 (L)")          # Wi 2017 Table 2: V2 = 71.7 L (RSE 28.1%)

    # Covariate effects on V1, Q, V2 (Wi 2017 Table 2 and final-model equation).
    # Encoded so that the model() expressions reproduce the paper's printed forms
    # V1 = 15.7 * (1 - 0.34 * ECMO), Q = 5.57 * (1 - 0.50 * ECMO), V2 = 71.7 * 1.5^CRRT
    # verbatim. The leading "(1 -" form is preserved here for direct comparison with
    # the published equation; the sign convention is therefore: e_ecmo_vc and e_ecmo_q
    # are positive fractional decrements (larger value = larger decrement on V1 / Q).
    e_ecmo_vc <- 0.34;      label("Fractional decrement in V1 when on ECMO (unitless)") # Wi 2017 Table 2: theta_ECMO on V1 = -0.34 (RSE 16.3%); paper formula V1 = TVV1 * (1 - 0.34 * ECMO)
    e_ecmo_q  <- 0.50;      label("Fractional decrement in Q when on ECMO (unitless)")  # Wi 2017 Table 2: theta_ECMO on Q = -0.50 (RSE 74.6%); paper formula Q = TVQ * (1 - 0.50 * ECMO)
    e_crrt_vp <- 1.50;      label("Power-model CRRT multiplier on V2 (unitless)")       # Wi 2017 Table 2: theta_CRRT on V2 = 1.50 (RSE 25.3%); paper formula V2 = TVV2 * 1.5^CRRT (multiplier 1.5 raised to the CRRT indicator)

    # Between-subject variability from Wi 2017 Table 2 (random-effects block; reported
    # values are omega^2 of the exponential IIV model as printed by NONMEM $OMEGA).
    # No off-diagonal covariances were published despite the Methods statement that a
    # variance-covariance matrix was estimated; independent etas are used here.
    # No IIV is reported for Q or V2 -- only CL, V1 (Vc), and Q variability rows appear
    # in Table 2, but the Q row reports omega^2 = 0.15 with RSE 52.3%, so Q does carry
    # an estimated IIV. V2 does not.
    etalcl ~ 0.34  # Wi 2017 Table 2: omega^2 = 0.34 on CL (RSE 463%); corresponds to approx 65% CV on linear scale via sqrt(exp(0.34) - 1)
    etalvc ~ 0.13  # Wi 2017 Table 2: omega^2 = 0.13 on V1 (RSE 54.7%); corresponds to approx 37% CV on linear scale
    etalq  ~ 0.15  # Wi 2017 Table 2: omega^2 = 0.15 on Q  (RSE 52.3%); corresponds to approx 40% CV on linear scale

    # Combined additive + proportional residual error on the linear-concentration scale
    # (Wi 2017 Methods: "Residual variability was described using the following combined
    # additive and proportional model" and Table 2 random-effects 'Residual variability' rows).
    addSd  <- 3.32;  label("Additive residual error (mg/L)")            # Wi 2017 Table 2: additive = 3.32 ug/mL (RSE 84.3%); ug/mL is numerically equivalent to mg/L
    propSd <- 0.243; label("Proportional residual error (fraction)")    # Wi 2017 Table 2: proportional = 24.3% (RSE 4.12%); printed as a percentage, encoded here as a fraction
  })
  model({
    # Individual PK parameters (log-normal IIV; reference values are the typical
    # population values from Wi 2017 final model, expressed before any covariate effect).
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc) * (1 - e_ecmo_vc * ECMO_STATUS)
    q  <- exp(lq  + etalq)  * (1 - e_ecmo_q  * ECMO_STATUS)
    vp <- exp(lvp)          * e_crrt_vp^CRRT_STATUS

    # Micro-constants for the explicit two-compartment ODE.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # IV bolus into central; first-order distribution to peripheral1 and first-order
    # elimination from central. Dose in mg, volumes in L -> central / vc has units mg/L.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
