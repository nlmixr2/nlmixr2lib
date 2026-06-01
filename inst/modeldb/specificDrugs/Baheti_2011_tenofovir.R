Baheti_2011_tenofovir <- function() {
  description <- "Two-compartment first-order-absorption population PK model for plasma tenofovir (TFV) in HIV-1-infected adults on once-daily tenofovir disoproxil fumarate (TDF) coupled with a stimulatory indirect-response (Dayneka 1993) model for intracellular tenofovir diphosphate (TFV-DP) in peripheral blood mononuclear cells; plasma TFV drives TFV-DP formation through a sigmoidal Emax stimulation function. Creatinine clearance enters CL/F and Vc/F via a power covariate. Fitted sequentially (PK first, PD with PK individual post-hoc Bayes estimates fixed)."
  reference <- "Baheti G, Kiser JJ, Havens PL, Fletcher CV. Plasma and intracellular population pharmacokinetic analysis of tenofovir in HIV-1-infected patients. Antimicrob Agents Chemother. 2011;55(11):5294-5299. doi:10.1128/AAC.05317-11"
  vignette <- "Baheti_2011_tenofovir"
  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "ng/mL"
  )
  # Dose note: doses entering the depot are TFV-equivalent mg, not TDF-equivalent
  # mg. The paper fitted apparent CL/F = 42 L/h with dose corrected to TFV per
  # the Jullien 2005 convention the authors cite (Baheti 2011 Discussion: "Jullien
  # et al. reported a population estimate for TFV CL/F of 50.5 liters/h (after
  # correction of TFV dose from TFV disoproxil to TFV)"). To simulate a 300 mg
  # TDF oral dose, convert with the molecular-weight ratio MW(TFV)/MW(TDF) =
  # 287.213 / 635.51 = 0.452, so amt = 300 * 0.452 = 135.6 mg TFV-equivalent.
  # TFV-DP is measured in fmol per 10^6 PBMCs (paper-specific units, not a
  # volume-based concentration); see TFVDP output below.

  covariateData <- list(
    CRCL = list(
      description = "Estimated creatinine clearance by the Cockcroft-Gault method (raw mL/min, not BSA-normalized).",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = "Baheti 2011 Methods 'Population pharmacokinetic modeling': estimated CrCL by the Cockcroft-Gault method (Cockcroft DW, Gault MH. Nephron 1976;16:31-41) and entered as a power covariate on CL/F and Vc/F. The reference (centering) value used here is 108 mL/min, the overall median CrCL of the pooled cohort (Baheti 2011 Table 1 'CrCL (ml/min)' Overall median). The CRCL canonical column is BSA-normalized in the register; Baheti 2011 uses the raw Cockcroft-Gault value without BSA-normalization, matching the CLCR source-alias precedent in Delattre_2010_amikacin.R.",
      source_name = "CrCL"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Subject age in years at study entry.",
      units = "year",
      type = "continuous",
      notes = "Tested in univariate analysis (Baheti 2011 Methods 'Population pharmacokinetic modeling', Table S1 in supplemental). Significant on CL/F and Vc/F in univariate analysis but dropped during stepwise forward addition: 'age and CrCL were significant on CL/F and for final model development only CrCL on CL/F was utilized, as age is a variable in calculation of CrCL'."
    ),
    WT = list(
      description = "Body weight in kg at study entry.",
      units = "kg",
      type = "continuous",
      notes = "Screened in univariate analysis but not retained in the final model (Baheti 2011 Methods 'Population pharmacokinetic modeling', and Discussion: 'Sex, weight, TBIL and concomitant PI medications did not affect TFV CL/F or V/F')."
    ),
    SEXF = list(
      description = "Subject sex indicator (1 = female).",
      units = "(binary)",
      type = "binary",
      notes = "Screened (Baheti 2011 Table 1 19 of 55 female, 35%) but not retained in the final model. Discussion: 'Sex ... did not affect TFV CL/F or V/F'."
    ),
    TBIL = list(
      description = "Total bilirubin in mg/dL at study entry.",
      units = "mg/dL",
      type = "continuous",
      notes = "Baheti 2011 Table 1 (overall median 1.1, range 0.2-4.7). Screened but not retained: Discussion 'TBIL ... did not affect TFV CL/F or V/F'."
    ),
    PI = list(
      description = "Concomitant boosted protease inhibitor co-medication (atazanavir/ritonavir or lopinavir/ritonavir) indicator.",
      units = "(binary)",
      type = "binary",
      notes = "Baheti 2011 Table 1 (ATV/r 50% of ATN056 cohort, LPV/r 50% of 1427 cohort, no-PI 50% of 1427). Screened in univariate analysis and significant on CL/F, Vc/F, Vp/F. Not retained after backward elimination -- Discussion: 'concomitant PI medications did not affect TFV CL/F or V/F'."
    ),
    RACE_BLACK = list(
      description = "Self-reported Black or African American race indicator (Baheti 2011 Table 1 Race/ethnicity row 'Black/African American').",
      units = "(binary)",
      type = "binary",
      notes = "Baheti 2011 Table 1 23/55 (42%) Black/African American. Race was screened among demographic covariates but not retained in the final model."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 55,
    n_studies = 2,
    age_range = "18.6-60 years; pooled-cohort median 33 (mean 32); cohort 1427 median 41.5 (mean 39.5), range 25-60; cohort ATN056 median 22.8 (mean 22.9), range 18.6-25",
    weight_range = "38.7-131.6 kg; pooled-cohort median 79.1 (mean 74.7); cohort 1427 median 79.6, range 38.7-121.4; cohort ATN056 median 78.6, range 46.9-131.6",
    sex_female_pct = 34.5,
    race_ethnicity = c(Black = 42, White = 45, Asian = 2, Other = 9, Unknown = 2),
    hispanic_latino_pct = 13,
    disease_state = "HIV-1 infection on stable combination antiretroviral therapy; near-fully suppressed (all but 6 subjects had HIV-1 RNA below 400 copies/mL at PK sampling); pooled across two parent studies",
    dose_range = "300 mg TDF (tenofovir disoproxil fumarate) PO once daily following a meal; steady-state sampling over the 24-h dosing interval",
    regions = "United States (single-site protocol 1427 and multi-site Adolescent Trials Network protocol ATN056)",
    renal_function = "median Cockcroft-Gault CrCL 108 mL/min (range 43.2-227.1); cohort 1427 median 95.1 (43.2-154.9), cohort ATN056 median 127.1 (65.7-227.1)",
    notes = "Pooled from two HIV-1 cohorts: 1427 (n = 30 adults aged 25-60, LPV/r vs no-PI arms) and ATN056 (n = 25 adolescents and young adults aged 18.6-25, ATV/r). Steady-state PK sampling on 8 (ATN056) or 11 (1427) timepoints over 24 h; intracellular PBMC TFV-DP at 3 timepoints per subject (predose, 5, 24 h for 1427; 1, 4, 24 h for ATN056) in 51 of 55 subjects. Counts: 529 plasma TFV samples + 151 intracellular TFV-DP samples. Demographics from Baheti 2011 Table 1; covariate screen results from Methods 'Population pharmacokinetic modeling' / Table S1 in supplemental."
  )

  ini({
    # Structural plasma PK (Baheti 2011 Table 2 'Final model' column).
    lka <- log(1.05)
    label("First-order absorption rate constant ka (1/h)")  # Table 2 Final model 'ka (h-1)' = 1.05 (95% CI 0.618-1.481); no IIV (Results 'the variance component for ... ka was not estimated')
    lcl <- log(42.0)
    label("Apparent clearance CL/F at reference CRCL = 108 mL/min (L/h)")  # Table 2 Final model 'CL/F' = 42.0 (95% CI 38.2-45.8)
    lvc <- log(277)
    label("Apparent central volume Vc/F at reference CRCL = 108 mL/min (L)")  # Table 2 Final model 'Vc/F' = 277 (95% CI 183-370)
    lq <- log(182)
    label("Apparent inter-compartmental clearance Q/F (L/h)")  # Table 2 Final model 'Q/F' = 182 (95% CI 150-213); no IIV (Results 'the variance component for ... Q/F ... was not estimated')
    lvp <- log(436)
    label("Apparent peripheral volume Vp/F (L)")  # Table 2 Final model 'Vp/F' = 436 (95% CI 348-523)

    # Covariate effects on plasma PK (CrCL as power function on the natural
    # scale, centered around the pooled-cohort median 108 mL/min per
    # Baheti 2011 Methods 'centered-around-median or a power function'
    # approach and Methods 'covariate analysis' / Table 2).
    e_crcl_cl <- 0.489
    label("CrCL power exponent on CL/F (unitless)")  # Table 2 Final model 'CrCL ~ CL/F' = 0.489 (95% CI 0.273-0.704)
    e_crcl_vc <- 1.01
    label("CrCL power exponent on Vc/F (unitless)")  # Table 2 Final model 'CrCL ~ Vc/F' = 1.01 (95% CI 0.58-1.43)

    # Intracellular TFV-DP (Baheti 2011 Table 3 'Final model' column; sequential
    # fit with plasma PK individual Bayes estimates fixed per Methods 'A
    # sequential analysis approach was utilized').
    lkin <- log(0.276)
    label("Zero-order TFV-DP production rate constant kin (1/h)")  # Table 3 Final model 'k in (h-1)' = 0.276 (bootstrap 95% CI 0.0145-1.56)
    lkout <- log(0.00808)
    label("First-order TFV-DP elimination rate constant kout (1/h)")  # Table 3 Final model 'k out (h-1)' = 0.00808 (bootstrap 95% CI 0.0007-0.0372); implies typical half-life ln(2)/kout = 86 h (Discussion '87-h TFV-DP half-life')
    lec50 <- log(99.9)
    label("Plasma TFV concentration at half-maximum TFV-DP stimulation EC50 (ng/mL)")  # Table 3 Final model 'EC 50 (ng/ml)' = 99.9 (bootstrap 95% CI 1.000-403); abstract rounds to 100 ng/mL
    lemax <- log(300)
    label("Maximal intracellular TFV-DP concentration Emax at saturating plasma TFV (fmol/10^6 cells)")  # Table 3 Final model 'E max (fmol/10^6 cells)' = 300 (bootstrap 95% CI 4.000-484)

    # IIV. Plasma PK: 3x3 block on CL/F, Vc/F, Vp/F with off-diagonal
    # covariances estimated (Baheti 2011 Methods 'the off-diagonal elements of
    # the covariance matrix were estimated'); ka and Q have no IIV (Results
    # 'Due to data limitations, the variance component for apparent
    # intercompartmental clearance (Q/F) and absorption rate (ka) was not
    # estimated'). Reported IIV percentages are sqrt(omega^2) per the
    # exponential parameterization (Results 'An exponential error model
    # described interpatient random effects'); diagonal entries:
    #   IIV CL/F = 33.5%  -> var = 0.335^2 = 0.112225
    #   IIV Vc/F = 64.8%  -> var = 0.648^2 = 0.419904
    #   IIV Vp/F = 46.5%  -> var = 0.465^2 = 0.216225
    # Off-diagonals from Table 2 Final model:
    #   Cov CL/F-Vc/F =  0.118  (R = 0.553)
    #   Cov Vc/F-Vp/F = -0.0415 (R = -0.139)
    #   Cov CL/F-Vp/F =  0.113  (R = 0.731)
    # Cross-check: R(CL,Vc) = 0.118 / (0.335 * 0.648) = 0.544 (Table 2 reports 0.553);
    # R(Vc,Vp) = -0.0415 / (0.648 * 0.465) = -0.138 (Table 2 reports -0.139);
    # R(CL,Vp) = 0.113 / (0.335 * 0.465) = 0.725 (Table 2 reports 0.731). Within rounding.
    etalcl + etalvc + etalvp ~ c(
      0.112225,
      0.118, 0.419904,
      0.113, -0.0415, 0.216225
    )

    # TFV-DP IIV on EC50 and Emax only (Baheti 2011 Results 'Interindividual
    # variability was described only for EC 50 and E max'). Reported IIV
    # percentages (sqrt of variance):
    #   IIV EC50 = 106%  -> var = 1.06^2 = 1.1236  (Table 3 Final model 'IIV EC50' = 106%)
    #   IIV Emax = 82.3% -> var = 0.823^2 = 0.677329  (Table 3 Final model 'IIV E max' = 82.3%)
    etalec50 ~ 1.1236
    etalemax ~ 0.677329

    # Residual error -- proportional for both plasma TFV and intracellular
    # TFV-DP (Baheti 2011 Results 'Residual error was modeled as a
    # proportional error model' for plasma; Results 'The residual error model
    # was a proportional model' for TFV-DP).
    propSd <- 0.183
    label("Proportional residual SD on plasma TFV (fraction)")  # Table 2 Final model 'RUV (% CV)' = 18.3% (95% CI 16.0-20.3); abstract rounds to 18%
    propSd_TFVDP <- 0.567
    label("Proportional residual SD on intracellular TFV-DP (fraction)")  # Table 3 Final model 'RUV (% CV)' = 56.7% (bootstrap 95% CI 46-63)
  })

  model({
    # Reference CrCL for the power covariate (Baheti 2011 Table 1 overall
    # cohort median; Methods 'centered-around-median ... approach').
    crcl_ref <- 108

    # Typical-value and individual plasma PK parameters. ka and Q have no
    # IIV (paper Results).
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (CRCL / crcl_ref)^e_crcl_cl
    vc <- exp(lvc + etalvc) * (CRCL / crcl_ref)^e_crcl_vc
    q  <- exp(lq)
    vp <- exp(lvp + etalvp)

    # TFV-DP typical-value and individual parameters. kin and kout have no IIV.
    kin   <- exp(lkin)
    kout  <- exp(lkout)
    ec50  <- exp(lec50 + etalec50)
    emax  <- exp(lemax + etalemax)

    # Micro-constants for the two-compartment PK.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Two-compartment first-order absorption PK (Baheti 2011 Results
    # 'two-compartment model with first-order absorption ... ADVAN4 TRANS4'
    # and Fig. 1 schematic). Doses enter `depot` in TFV-equivalent mg
    # (see header dose note).
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Plasma TFV concentration. central is in mg, vc is in L, so
    # central / vc is in mg/L = 1000 ng/mL.
    Cc <- 1000 * central / vc

    # Intracellular TFV-DP indirect response (Dayneka stimulation of kin).
    # Baheti 2011 Methods writes:
    #     S(t) = 1 + Emax * Cp / (EC50 + Cp)
    #     dR/dt = kin * S(t) - kout * R
    # interpreted with Emax in fmol/10^6 cells representing "the maximal
    # intracellular concentration of TFV-DP" (Methods text after the
    # equation, and Fig. 1 caption). Reading S(t)'s Emax as a dimensionless
    # stimulation factor would give R_ss > 4000 fmol/10^6 cells at typical
    # plasma TFV concentrations (Cp ~ 72-135 ng/mL), inconsistent with
    # Baheti 2011 Discussion ('would predict TFV-DP concentrations of 128 to
    # 174 fmol/10^6 cells using a typical trough (72 ng/ml) or steady-state
    # average (135 ng/ml) plasma concentration of TFV') and with the observed
    # mean of 90.7 fmol/10^6 cells. The internally-consistent encoding
    # treats Emax as the asymptote of R at saturating Cp:
    #     dR/dt = kin + (kout * Emax - kin) * Cp / (EC50 + Cp) - kout * R
    # which gives R_ss(Cp -> 0) = kin/kout = 34 fmol/10^6 cells (baseline)
    # and R_ss(Cp -> inf) = Emax = 300 fmol/10^6 cells. At Cp = 72 the
    # predicted R_ss is 145 fmol/10^6 cells; at Cp = 135 it is 187
    # fmol/10^6 cells -- matching the Baheti 2011 Discussion range 128-174.
    # See vignette 'Assumptions and deviations' for the full derivation.
    stim_frac <- Cc / (ec50 + Cc)
    d/dt(effect) <- kin + (kout * emax - kin) * stim_frac - kout * effect

    # Initial intracellular TFV-DP at simulation start (Cp = 0). Pre-treatment
    # baseline equilibrium of the indirect-response model. Subjects in the
    # actual study were at steady state on TDF; simulating from t = 0 with a
    # TDF dose lets the long TFV-DP half-life (~ 87 h) accumulate to the
    # individual-subject steady state.
    effect(0) <- kin / kout

    # TFV-DP output for residual error. fmol/10^6 cells (paper-mechanistic
    # units, not a volume-based concentration; not transformable to ng/mL
    # without a cell-volume assumption -- Baheti 2011 Methods 'No conversion
    # of TFV-DP concentrations to the same units as plasma was made').
    TFVDP <- effect

    Cc    ~ prop(propSd)
    TFVDP ~ prop(propSd_TFVDP)
  })
}
