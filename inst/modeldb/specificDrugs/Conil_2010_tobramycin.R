Conil_2010_tobramycin <- function() {
  description <- "Two-compartment IV population PK model for tobramycin in adult ICU patients receiving once-daily aminoglycoside therapy for nosocomial Gram-negative infections (Conil 2010); additive linear covariate effects of Cockcroft-Gault creatinine clearance and height on CL, with no IIV on Q or V2."
  reference <- "Conil JM, Georges B, Ruiz S, Rival T, Seguin T, Cougot P, Fourcade O, Houin G, Saivin S. Tobramycin disposition in ICU patients receiving a once daily regimen: population approach and dosage simulations. Br J Clin Pharmacol. 2011;71(1):61-71. doi:10.1111/j.1365-2125.2010.03793.x"
  vignette <- "Conil_2010_tobramycin"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    CRCL = list(
      description        = "Cockcroft-Gault creatinine clearance (raw, not BSA-normalized)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source column COCK in the NONMEM control stream. Computed by the Cockcroft-Gault equation in raw mL/min (NOT BSA-normalized to mL/min/1.73 m^2). Stored under the canonical CRCL column per inst/references/covariate-columns.md (raw-Cockcroft precedent: Delattre_2010_amikacin.R). Reference value 94 mL/min (population median, Conil 2010 Methods - 'Simulation and dosage regimen propositions'). Used in an additive linear covariate term on CL: TVCL = exp(lcl) + e_crcl_cl * (CRCL - 94).",
      source_name        = "COCK"
    ),
    HT = list(
      description        = "Body height at baseline",
      units              = "cm",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source column HEIG in the NONMEM control stream. Reference value 172 cm (population median, Conil 2010 Methods - 'Simulation and dosage regimen propositions'). Used in an additive linear covariate term on CL: TVCL = exp(lcl) + e_ht_cl * (HT - 172). Per the Conil 2010 Discussion, height correlated more strongly with CL than total body weight or ideal body weight; in ICU patients body weight is biased by oedema, and the authors retained height as a more stable size descriptor of clearance.",
      source_name        = "HEIG"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 49L,
    n_studies      = 1L,
    n_observations = 277L,
    age_range      = "not reported (excluded < 15 years)",
    age_median     = "61.1 +/- 16.5 years (mean +/- SD)",
    weight_range   = "not reported",
    weight_median  = "78.7 +/- 18.9 kg (mean +/- SD)",
    height_median  = "172 +/- 9 cm (mean +/- SD)",
    sex_female_pct = 18,
    race_ethnicity = "Not reported (French university-hospital ICU population)",
    disease_state  = "Nosocomial Gram-negative infections in adult ICU patients (poly-trauma 41%, post-surgical 14%, medical/respiratory 45%); all patients were mechanically ventilated and haemodynamically stable; SAPS I 15 +/- 4 and SAPS II 56 +/- 14",
    dose_range     = "5 mg/kg tobramycin IV infusion over 30 minutes once daily for 3 to 5 days, given in combination with a beta-lactam (ceftazidime or imipenem)",
    regions        = "France (single centre: Toulouse-Rangueil University Hospital, ICU)",
    renal_function = "Cockcroft creatinine clearance 108 +/- 55 mL/min (raw, not BSA-normalized); Robert creatinine clearance 75 +/- 29 mL/min; serum creatinine 84.2 +/- 37.2 umol/L",
    exclusion      = "Pregnancy, age < 15 years, drug allergies / intolerance to aminoglycosides, oligo-anuric renal failure, cochlear problems",
    notes          = "Demographics from Conil 2010 Table 1 (Total n=49 column). 32 patients (182 concentrations) entered the model-building set; 17 patients (95 concentrations) formed the qualification set evaluated by npde. After qualification the two populations were pooled and final parameters were re-estimated on all 49 patients - this is the model encoded here (Conil 2010 Table 2, 'Final model (n=49)' column). Study period October 2005 - December 2007. Indication-driven combination antibiotic therapy is documented but the structural PK model is for tobramycin alone."
  )

  ini({
    # Structural parameters (Conil 2010 Table 2, 'Final model (n=49)' column).
    # The typical CL is parameterised as a linear-additive model centred at the
    # population reference covariates (CRCL = 94 mL/min, HT = 172 cm). The
    # intercept exp(lcl) = THETA(1) = 3.83 L/h is the typical CL at reference;
    # additive deviations from reference are added through e_crcl_cl and e_ht_cl
    # below. Q and V2 carry no IIV in the final model (the paper fixed Omega 3
    # and Omega 4 to zero).
    lcl <- log(3.83);  label("Typical CL at the reference patient (CRCL=94 mL/min, HT=172 cm) (L/h)") # Conil 2010 Table 2: THETA(1) = 3.83 L/h (Final model n=49)
    lvc <- log(25.5);  label("Central volume of distribution V1 (L)")                                  # Conil 2010 Table 2: THETA(4) = 25.50 L (Final model n=49)
    lq  <- log(4.74);  label("Inter-compartmental clearance Q (L/h)")                                  # Conil 2010 Table 2: THETA(5) = 4.74 L/h (Final model n=49)
    lvp <- log(30.6);  label("Peripheral volume of distribution V2 (L)")                               # Conil 2010 Table 2: THETA(6) = 30.60 L (Final model n=49)

    # Covariate effects on CL: linear-additive deviations from the reference
    # covariate values, applied to the typical (intercept) CL exp(lcl).
    # Form (Conil 2010 Table 2 equation): TVCL = THETA1 + THETA2*(COCK-94) + THETA3*(HEIG-172).
    e_crcl_cl <- 0.020; label("Slope of CL per (CRCL - 94 mL/min) (L/h per mL/min)") # Conil 2010 Table 2: THETA(2) = 0.020 L/h per mL/min (Final model n=49)
    e_ht_cl   <- 0.052; label("Slope of CL per (HT - 172 cm) (L/h per cm)")          # Conil 2010 Table 2: THETA(3) = 0.052 L/h per cm   (Final model n=49)

    # Inter-individual variability (Conil 2010 Table 2, 'Final model' column).
    # Omega values are variances on the log-scale because the paper applies
    # IIV multiplicatively as EXP(ETA): CL = TVCL * exp(eta1); V1 = TVV1 * exp(eta2).
    # Q and V2 IIV were fixed to zero in the final model and are not included.
    etalcl ~ 0.095 # Conil 2010 Table 2: Omega 1 = 0.095 on log-scale (Final model n=49); sqrt(exp(0.095)-1) ~= 31% CV, close to reported 29%
    etalvc ~ 0.045 # Conil 2010 Table 2: Omega 2 = 0.045 on log-scale (Final model n=49); sqrt(exp(0.045)-1) ~= 21% CV

    # Residual error (proportional only). Conil 2010 Table 2 reports Sigma 1
    # = 0.055 in the final model; in the NONMEM convention Y = F*(1 + EPS) this
    # is the variance of EPS, so the SD on the proportional scale is sqrt(0.055).
    propSd <- sqrt(0.055); label("Proportional residual error SD (fraction)") # Conil 2010 Table 2: Sigma 1 = 0.055 (Final model n=49); SD = sqrt(0.055) ~= 0.234 (CV 23%)
  })
  model({
    # Individual PK parameters. CL takes the additive linear-deviation form on
    # CRCL and HT; eta on CL is applied multiplicatively via exp(etalcl).
    # V1 has log-normal IIV; Q and V2 have no IIV.
    cl <- (exp(lcl) + e_crcl_cl * (CRCL - 94) + e_ht_cl * (HT - 172)) * exp(etalcl)
    vc <- exp(lvc + etalvc)
    q  <- exp(lq)
    vp <- exp(lvp)

    # Two-compartment IV micro-constants.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # Dose in mg, V in L -> central / vc has units mg/L.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
