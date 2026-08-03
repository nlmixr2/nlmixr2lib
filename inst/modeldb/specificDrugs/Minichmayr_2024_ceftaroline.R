Minichmayr_2024_ceftaroline <- function() {
  description <- "Two-compartment population PK model with linear elimination for ceftaroline (the active metabolite of the prodrug ceftaroline fosamil) in 12 healthy male volunteers receiving either the standard dosing regimen (ceftaroline fosamil 600 mg every 12 h as a 1 h intravenous infusion) or the approved intensified regimen (600 mg every 8 h as a 2 h intravenous infusion). Doses are entered as milligrams of the administered prodrug ceftaroline fosamil; the model converts to ceftaroline via the fixed ceftaroline-to-ceftaroline-fosamil molar-mass ratio of 0.883 applied as the bioavailability of the central compartment, reflecting the paper's finding that prodrug-to-drug conversion is complete and effectively instantaneous (a depot / k_trans conversion compartment worsened the fit, and a sensitivity analysis favoured k_trans > 1000/h). Exponential interindividual variability on CL, Vc and Vp; combined additive and proportional residual error. No covariates were retained: creatinine clearance was screened on CL as linear and power functions and via a previously published CL-CLcr relationship, but was not significant in this small, homogeneous cohort with unimpaired renal function. The model underpins the paper's systematic comparison of the three components of dosing intensification (dosing interval, infusion duration and total daily dose) on the PK/PD index fT>MIC against Staphylococcus aureus."
  reference <- "Minichmayr IK, Wicha SG, Matzneller P, Kloft C, Zeitlinger M. Impact of key components of intensified ceftaroline dosing on pharmacokinetic/pharmacodynamic target attainment. Clin Pharmacokinet. 2024;63(1):121-131. doi:10.1007/s40262-023-01325-4. Final population PK parameter estimates are Table 1; the structural-model selection narrative (two-compartment, linear elimination, complete and immediate prodrug conversion, no CLcr covariate) is Results section 3.1; the ceftaroline-to-ceftaroline-fosamil molar-mass ratio of 0.883 and the 20% plasma protein binding used for fT>MIC are Methods sections 2.2 and 2.3; published median fT>MIC values across the eight systematically varied dosing regimens are Tables 2 and 3. The electronic supplementary material (ESM) adds baseline demographics (Table S1), the objective function values of the key candidate models explored during model development (Table S2), and the complete 72-cell grid of median fT>MIC across six MIC values, two dosing intervals, two total daily doses and three infusion durations including the 3 h infusions omitted from the main text (Table S3)."
  vignette <- "Minichmayr_2024_ceftaroline"
  units <- list(
    time          = "hour",
    dosing        = "mg (ceftaroline fosamil, the administered prodrug)",
    concentration = "mg/L (total plasma ceftaroline)"
  )

  # No covariates were retained in the final model. Creatinine clearance was
  # the only covariate formally screened (on CL); see covariatesDataExcluded.
  covariateData <- list()

  # Screened but not retained. Creatinine clearance (Cockcroft-Gault) was
  # investigated as a covariate on CL using linear and power functions as well
  # as a CL-CLcr relationship previously identified for healthy subjects
  # (Methods 2.2), but "Inclusion of a CLCR-CL covariate-parameter relationship
  # was not statistically significant (delta-OFV < 1), irrespective of its
  # functional form" (Results 3.1). The paper attributes this to the small,
  # homogeneous cohort with a narrow range of unimpaired renal function
  # (CLcr_Cockcroft-Gault 93.2-165 mL/min; Discussion). No point estimate is
  # reported for any of the screened functional forms, so the effect cannot be
  # encoded even as a fixed term.
  covariatesDataExcluded <- list(
    CRCL = list(
      description = "Creatinine clearance estimated by the Cockcroft-Gault equation, reported in raw mL/min (not BSA-normalized) for this cohort",
      units       = "mL/min",
      type        = "continuous",
      notes       = "Screened as a covariate on CL with linear and power functions and with a published healthy-subject CL-CLcr relationship; not statistically significant (delta-OFV < 1) and not retained in the final model. Median 145 mL/min (range 93.2-165). Source column name in the study data was CLCR.",
      source_name = "CLCR"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 12L,
    n_studies        = 1L,
    n_observations   = 274L,
    age_range        = "22-50 years (main text); 5th-95th percentile 23.1-45.6 years (ESM Table S1)",
    age_median       = "27.5 years",
    weight_range     = "63-106 kg (main text); 5th-95th percentile 64.4-96.1 kg (ESM Table S1)",
    weight_median    = "74.5 kg",
    height_median    = "183 cm (5th-95th percentile 173-194 cm; ESM Table S1)",
    bsa_median       = "1.92 m^2 (5th-95th percentile 1.79-2.28 m^2; ESM Table S1)",
    bmi_range        = "19.4-27.6 kg/m^2 (main text); 5th-95th percentile 19.4-25.9 kg/m^2 (ESM Table S1)",
    bmi_median       = "23.0 kg/m^2",
    sex_female_pct   = 0,
    disease_state    = "healthy volunteers (no renal impairment)",
    renal_function   = "CLcr (Cockcroft-Gault) median 145 mL/min, range 93.2-165 mL/min (main text; 5th-95th percentile 95.4-158 mL/min, ESM Table S1); serum creatinine median 0.88 mg/dL (5th-95th percentile 0.76-1.04 mg/dL, ESM Table S1); no renal impairment",
    dose_range       = "Ceftaroline fosamil 600 mg every 12 h as a 1 h intravenous infusion (standard dosing, n = 6) or 600 mg every 8 h as a 2 h intravenous infusion (intensified dosing, n = 6)",
    regions          = "Austria (Medical University of Vienna)",
    notes            = paste(
      "Prospective, open-label study; EudraCT 2012-005134-11 (Study Population and Methods 2.1; demographics in ESM Table S1, which reports medians with 5th-95th percentiles, and in the main-text Methods narrative, which reports medians with full ranges).",
      "All 12 participants were male, so sex was not a candidate covariate.",
      "Rich sampling of total plasma ceftaroline on two occasions: after the first dose and after three (q12h group) or four (q8h group) repeated doses; n = 274 concentrations total.",
      "Interoccasion variability between the two sampling days was not supported by the data (delta-OFV <= 5.23 for IOV on CL of 3%), and no correlations between interindividual random effects were supported (Results 3.1).",
      "Estimated with NONMEM 7.3.0 using FOCE with interaction, assisted by PsN 4.7.0.",
      "Plasma protein binding of ceftaroline was assumed to be 20% (unbound fraction fu = 0.8) when converting total concentrations to the unbound concentrations that drive fT>MIC (Methods 2.3). fu is NOT applied inside model() because the model was fit to TOTAL plasma concentrations; multiply Cc by 0.8 to obtain the unbound concentration used for PK/PD target attainment.",
      "Alternative structures explored and rejected: three-compartment (comparable fit, less precise and less plausible estimates); a depot compartment with a first-order prodrug conversion rate k_trans (increased OFV; sensitivity analysis favoured k_trans > 1000/h, i.e. virtually immediate conversion); Michaelis-Menten and mixed linear plus Michaelis-Menten elimination (no significant OFV decrease, RSE > 80%).",
      paste(
        "ESM Table S2 gives the objective function values behind those choices; the selected model (A: two-compartment, IIV on CL/V1/V2, combined additive plus proportional RUV) has OFV 39.498, against",
        "339.608 for the same model without IIV, 462.005 for one compartment without IIV and 336.511 for three compartments without IIV;",
        "38.918 for a linear CLcr-CL covariate relationship (delta-OFV 0.58, not significant against the 3.84 threshold);",
        "39.487 for Michaelis-Menten elimination; 34.268 / 39.506 / 37.248 for interoccasion variability on CL / V1 / V2 respectively (the CL case is the delta-OFV 5.23 quoted in Results 3.1);",
        "and 50.339 with k_trans fixed from a prodrug clearance of 228 L/h versus 39.693 with k_trans fixed to 1000/h, while estimating k_trans outright did not converge (k_trans exceeded 1e4/h)."
      )
    )
  )

  ini({
    # Structural disposition parameters. Table 1 reports these on the linear
    # scale with RSE and LLP-SIR 95% CIs; they are log-transformed here per
    # nlmixr2lib convention. All four were estimated (no FIX flags reported).
    lcl <- log(10.9);  label("Clearance (L/h)")                          # Minichmayr 2024 Table 1 (CL = 10.9 L/h; RSE 4.70%; 95% CI 10.1-12.1)
    lvc <- log(15.3);  label("Central volume of distribution (L)")       # Minichmayr 2024 Table 1 (V_C = 15.3 L; RSE 5.44%; 95% CI 13.8-16.7)
    lvp <- log(7.82);  label("Peripheral volume of distribution (L)")    # Minichmayr 2024 Table 1 (V_P = 7.82 L; RSE 8.90%; 95% CI 6.86-8.84)
    lq  <- log(4.82);  label("Intercompartmental clearance (L/h)")       # Minichmayr 2024 Table 1 (Q = 4.82 L/h; RSE 16.5%; 95% CI 3.91-5.97)

    # Interindividual variability. Table 1 reports omega on the %CV scale for
    # an exponential (lognormal) IIV model, so the variance on the log scale is
    # omega^2 = log(1 + CV^2):
    #   CL: log(1 + 0.156^2) = 0.0240446
    #   Vc: log(1 + 0.130^2) = 0.0167588
    #   Vp: log(1 + 0.166^2) = 0.0271832
    # The alternative reading of the %CV column as the approximate SD scale
    # (omega^2 = CV^2) gives 0.024336 / 0.016900 / 0.027556 -- within 1.4% of
    # the values used here, so the two readings are numerically equivalent at
    # this magnitude. Table 1 footnote a notes that the RSEs of these omegas
    # (not the omegas themselves) are on an approximate standard-deviation
    # scale. Correlations between the random effects were not supported by the
    # data (Results 3.1), so this is a diagonal omega.
    etalcl ~ 0.0240446   # Minichmayr 2024 Table 1 (omega_CL = 15.6 %CV; RSE 16.7%; 95% CI 11.1-22.7)
    etalvc ~ 0.0167588   # Minichmayr 2024 Table 1 (omega_VC = 13.0 %CV; RSE 19.3%; 95% CI 7.57-21.7)
    etalvp ~ 0.0271832   # Minichmayr 2024 Table 1 (omega_VP = 16.6 %CV; RSE 24.9%; 95% CI 10.1-28.9)

    # Combined additive plus proportional residual error, reported directly on
    # the SD scale (sigma_additive in mg/L; sigma_proportional as %CV, i.e. a
    # fraction of 0.135). A combined model was selected over additive-only and
    # proportional-only alternatives (Results 3.1).
    addSd  <- 0.0441; label("Additive residual error SD (mg/L)")             # Minichmayr 2024 Table 1 (sigma_additive = 0.0441 mg/L; RSE 44.2%; 95% CI 0.0200-0.0734)
    propSd <- 0.135;  label("Proportional residual error SD (fraction)")     # Minichmayr 2024 Table 1 (sigma_proportional = 13.5 %CV; RSE 12.7%; 95% CI 12.3-15.1)
  })

  model({
    # Ceftaroline-to-ceftaroline-fosamil molar-mass ratio (unitless). Doses are
    # entered as mg of the administered prodrug ceftaroline fosamil; this
    # stoichiometric factor converts the dosed amount to the ceftaroline
    # amount that the disposition parameters and the measured concentrations
    # refer to (Methods 2.2: "The difference between the administered and
    # measured entities (CPT-F and CPT) was considered using their ratio of
    # molar masses (CPT-to-CPT-F ratio: 0.883)"). Conversion was found to be
    # complete and immediate, so no depot / conversion compartment is used
    # (Results 3.1).
    cptRatio <- 0.883

    # Individual disposition parameters (exponential IIV).
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp + etalvp)
    q  <- exp(lq)

    # Micro-constants for the two-compartment system.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Two-compartment disposition with linear elimination. Ceftaroline fosamil
    # is given as an intravenous infusion, so doses go directly to central
    # (use rate = or dur = on the dose records).
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                  k12 * central - k21 * peripheral1

    # Prodrug-to-drug molar conversion, applied as the bioavailability of the
    # central compartment so that the user doses ceftaroline fosamil in mg.
    f(central) <- cptRatio

    # Total plasma ceftaroline concentration (mg/L). Unbound concentration for
    # fT>MIC is 0.8 * Cc (20% plasma protein binding, Methods 2.3); the model
    # itself was fit to total concentrations.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
