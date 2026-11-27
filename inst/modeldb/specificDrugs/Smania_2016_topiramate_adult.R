Smania_2016_topiramate_adult <- function() {
  description <- "Population PK-PD model of adjunctive topiramate (TPM) in adults with refractory partial-onset seizures, as implemented by Smania 2016 from the Girgis 2010 model to derive the adult prior for Bayesian pediatric trial designs. PK: two-compartment model with first-order absorption (fixed-effect transfer rate constants k12 and k21) whose clearance scales with body weight (power, centred at 69.9 kg) and age (exponential, centred at 31.4 years) and is 47.9% higher on adjunctive therapy than on topiramate monotherapy; central volume scales with body weight. IIV on CL, Vc and ka; combined proportional plus additive residual error. PD: the adult exposure-response is linear in the steady-state trough concentration, Y = b0 + b1 * Cmin + eps, where Y = log(100 * (S - B) / B + 110) is the log-transformed translated percent change in 28-day seizure frequency from baseline (S treatment, B baseline); placebo (Cmin = 0) gives Y = b0, a 24.6% seizure reduction. The concomitant-antiepileptic-specific clearance factors of the Girgis 2010 model were switched off in the Smania simulation code and are not included."
  reference <- paste(
    "Smania G, Baiardi P, Ceci A, Cella M, Magni P. Model-Based Assessment of Alternative Study Designs in Pediatric Trials. Part II: Bayesian Approaches. CPT Pharmacometrics Syst Pharmacol. 2016;5(8):402. doi:10.1002/psp4.12092. PMCID: PMC4999603.",
    "Smania G, Baiardi P, Ceci A, Magni P, Cella M. Model-Based Assessment of Alternative Study Designs in Pediatric Trials. Part I: Frequentist Approaches. CPT Pharmacometrics Syst Pharmacol. 2016;5(6):305-312. doi:10.1002/psp4.12083. PMCID: PMC5131885.",
    "Parameter values originate from: Girgis IG, et al. Pharmacokinetic-pharmacodynamic assessment of topiramate dosing regimens for children with epilepsy 2 to <10 years of age. Epilepsia. 2010;51:1954-1962. doi:10.1111/j.1528-1167.2010.02598.x."
  )
  vignette <- "Smania_2016_topiramate"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on CL ((WT/69.9)^0.453) and on Vc ((WT/69.9)^1.14). The 69.9 kg reference is the centring value of the pooled adult-plus-pediatric Girgis 2010 PK model (Smania 2016 Part I Supplementary Material 2, NONMEM $PK FCWT and FVWT).",
      source_name = "BW"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Exponential effect on CL, exp(-0.00306 * (AGE - 31.4)), centred at 31.4 years (Smania 2016 Part I Supplementary Material 2, NONMEM $PK FCAGE).",
      source_name = "AGE"
    ),
    CONMED_AED = list(
      description = "Adjunctive-therapy indicator: 1 = topiramate given on top of the patient's current antiepileptic treatment, 0 = topiramate monotherapy",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (topiramate monotherapy)",
      notes = "Multiplicative effect on CL, (1 + 0.479 * CONMED_AED). Table S1 of Smania 2016 Part I labels theta2 'CLSTA (effect of adjuvant)'; the qualifying concomitant antiepileptic drugs are not enumerated in any source on disk. The adult PD data were from adjunctive-therapy trials in refractory epilepsy (Smania 2016 Part I Supplementary Material 1); Smania 2016 did not simulate adult trials.",
      source_name = "ADJ"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "topiramate", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "topiramate", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "topiramate", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 663,
    n_studies = 8,
    age_range = "adults (range not reported in the Smania 2016 sources); the PK model was fitted to pooled pediatric and adult data",
    weight_range = "not reported in the Smania 2016 sources",
    disease_state = "Treatment-refractory epilepsy; topiramate given as adjunctive therapy",
    dose_range = "Average adult regimen used by Smania 2016 to derive the adult treatment effect: 150 mg oral topiramate twice daily",
    notes = "Smania 2016 Part I Supplementary Material 1: PD data came from eight clinical trials in refractory epileptic patients (pediatric and adult). Smania 2016 Part II Supplementary Table S1: pediatric PD fit n = 115, adult PD fit n = 663. The PK model is shared with the pediatric PK-PD model (Smania_2016_topiramate_pediatric); only the PD intercept and slope differ. Smania 2016 Part II used the adult slope with a population trough of 8 mg/L to set the adult treatment effect (Supplementary Material 5, BD_simulations.R). Sex was not a covariate."
  )

  ini({
    # PK -- Smania 2016 Part I Supplementary Table S1 ('adapted from Girgis et al.') and
    # Supplementary Material 2 (NONMEM $THETA); references 69.9 kg and 31.4 years.
    lcl <- log(1.21)
    label("Clearance on topiramate monotherapy at 69.9 kg and 31.4 years (L/h)") # Part I Table S1 'CLSTM (baseline clearance monotherapy) (theta1)' = 1.21 (1.2 %SE); NONMEM THETA(1)
    e_conmed_aed_cl <- 0.479
    label("Fractional increase in CL on adjunctive therapy (unitless)") # Part I Table S1 'CLSTA (effect of adjuvant) (theta2)' = 0.479 (25.3 %SE); CL = theta1 * (1 + ADJ * theta2)
    e_wt_cl <- 0.453
    label("Power exponent on WT/69.9 for CL (unitless)") # Part I Table S1 'FCWT (effect of weight) (theta3)' = 0.453 (9.0 %SE)
    e_age_cl <- -0.00306
    label("Exponential age coefficient on CL, exp(e * (AGE - 31.4)) (1/year)") # Part I Table S1 'FCAGE (effect of age) (theta4)' = -0.00306 (30.9 %SE)
    lvc <- log(4.61)
    label("Central volume of distribution at 69.9 kg (L)") # Part I Table S1 'VST (theta8)' = 4.61 (33.2 %SE)
    e_wt_vc <- 1.14
    label("Power exponent on WT/69.9 for Vc (unitless)") # Part I Table S1 'FVWT (effect of weight) (theta9)' = 1.14 (19.1 %SE)
    lka <- log(0.105)
    label("First-order absorption rate constant (1/h)") # Part I Table S1 'Ka (h-1) (theta10)' = 0.105 (27.0 %SE)
    lk12 <- log(0.577)
    label("Transfer rate constant central -> peripheral1 (1/h)") # Part I Table S1 'K23 (h-1) (theta11)' = 0.577 (16.7 %SE); NONMEM ADVAN4 K23
    lk21 <- log(0.0586)
    label("Transfer rate constant peripheral1 -> central (1/h)") # Part I Table S1 'K32 (h-1) (theta12)' = 0.0586 (23.6 %SE); NONMEM ADVAN4 K32

    # PD -- adult column of Smania 2016 Part II Supplementary Table S1 (n = 663)
    e0_seizY <- 4.4469
    label("Placebo response on the Y = log(percent change + 110) scale (unitless)") # Part II Table S1 'Placebo effect beta0', Adults = 4.4469 (3.13 %SE)
    slope_seizY <- -0.0627
    label("Change in Y per unit steady-state trough concentration (L/mg)") # Part II Table S1 'Concentration effect beta1', Adults = -0.0627 (0.97 %SE); BD_simulations.R beta1_a

    # IIV -- NONMEM $OMEGA of Part I Supplementary Material 2, entered as (%CV/100)^2
    etalcl ~ 0.07441984 # $OMEGA 'CL --> (27.28/100)^2'; Part I Table S1 BSV 27.28 %
    etalvc ~ 1.350244 # $OMEGA 'V --> (116.2/100)^2'; Part I Table S1 BSV 116.2 %
    etalka ~ 0.04990756 # $OMEGA 'KA --> (22.34/100)^2'; Part I Table S1 BSV 22.34 %

    # Residual error
    propSd <- 0.2546
    label("Proportional residual error on topiramate concentration (fraction)") # Part I Table S1 'CCV residual error (%CV)' = 25.46; $SIGMA (25.46/100)^2
    addSd <- 0.1797
    label("Additive residual error on topiramate concentration (mg/L)") # Part I Table S1 'Additive residual error (mg/L)' = 0.1797; $SIGMA 0.1797^2
    addSd_seizY <- 0.751664
    label("Additive residual error on Y (unitless)") # Part II Table S1 sigma, Adults = 0.7517 (same as pediatric); Part I Table S1 'sigma_eps' = 0.751664; BD_simulations.R sigma_epsilon and 'sA = sP'
  })

  model({
    # Individual PK parameters (NONMEM $PK of Part I Supplementary Material 2)
    cl <- exp(lcl + etalcl) * (1 + e_conmed_aed_cl * CONMED_AED) *
      (WT / 69.9)^e_wt_cl * exp(e_age_cl * (AGE - 31.4))
    vc <- exp(lvc + etalvc) * (WT / 69.9)^e_wt_vc
    ka <- exp(lka + etalka)
    k12 <- exp(lk12)
    k21 <- exp(lk21)
    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    Cc <- central / vc

    # Exposure-response: Y = b0 + b1 * Cmin. The source model uses the
    # steady-state trough concentration, so seizY reproduces the published
    # endpoint only on observation rows at a steady-state pre-dose trough.
    # Placebo subjects (no topiramate) have Cc = 0 and Y = b0.
    seizY <- e0_seizY + slope_seizY * Cc
    # Percent change in 28-day seizure frequency from baseline implied by Y
    # (negative = reduction): Y = log(100 * (S - B) / B + 110).
    seizPctChg <- exp(seizY) - 110

    Cc ~ prop(propSd) + add(addSd)
    seizY ~ add(addSd_seizY)
  })
}
