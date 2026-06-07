Rower_2017_tacrolimus <- function() {
  description <- paste0(
    "One-compartment population pharmacokinetic model for oral / enteral ",
    "tacrolimus in paediatric heart transplant recipients (Rower 2017): ",
    "first-order absorption with fixed Ka = 3.43 1/h; AGE power effect on ",
    "apparent volume with exponent 0.775 and reference 5.7 years; ",
    "creatinine-clearance power effect on apparent elimination rate with ",
    "exponent 0.850 and reference 122.4 mL/min/1.73 m^2; concomitant ",
    "fluconazole reduces apparent elimination by 34%. Originally ",
    "parameterised in NONMEM ADVAN2 TRANS1 on (ke, V); converted here to ",
    "the canonical (CL/F, V/F) form via CL/F = ke * V, so the AGE effect ",
    "propagates to CL/F with the same exponent as on V/F."
  )
  reference <- paste0(
    "Rower JE, Stockmann C, Linakis MW, Kumar SS, Liu X, Korgenski EK, ",
    "Sherwin CMT, Molina KM. Predicting tacrolimus concentrations in ",
    "children receiving a heart transplant using a population ",
    "pharmacokinetic model. BMJ Paediatrics Open. 2017;1:e000147. ",
    "doi:10.1136/bmjpo-2017-000147."
  )
  vignette <- "Rower_2017_tacrolimus"
  units <- list(time = "hour", dosing = "mg", concentration = "ug/L")

  covariateData <- list(
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Power-form effect on apparent volume of distribution with ",
        "exponent 0.775 and reference 5.7 years (model-building cohort ",
        "median; Table 1). In the canonical (CL/F, V/F) parameterisation ",
        "the same AGE exponent appears on CL/F because the source paper ",
        "estimated (ke, V) and CL/F = ke * V; see the model description."
      ),
      source_name        = "AGE"
    ),
    CRCL = list(
      description        = "Creatinine clearance, BSA-normalised (bedside Schwartz)",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste0(
        "Power-form effect on apparent elimination rate (and therefore ",
        "on CL/F in the canonical form) with exponent 0.850 and reference ",
        "122.4 mL/min/1.73 m^2 (model-building cohort median; Table 1). ",
        "Calculated using the bedside Schwartz equation (Methods, PK ",
        "modelling section). Time-varying within subject; missing values ",
        "carried forward or back up to 48 hours, beyond which the ",
        "population median was imputed (Methods, PK modelling section)."
      ),
      source_name        = "CrCL"
    ),
    CONMED_AZOLE = list(
      description        = "Concomitant azole antifungal therapy indicator (fluconazole in the source cohort)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant azole antifungal)",
      notes              = paste0(
        "Fluconazole coadministration reduces tacrolimus apparent ",
        "elimination by 34% (Table 2 final model: ke = 0.0408 /h without ",
        "fluconazole and 0.0268 /h with fluconazole; fractional change ",
        "-0.343, encoded as the e_azole_cl coefficient). 15 of 30 subjects ",
        "in the model-building cohort received fluconazole (Table 1). The ",
        "source paper labels this indicator FLUC (fluconazole only); the ",
        "canonical register entry CONMED_AZOLE generalises across azole ",
        "antifungals -- interpret it as fluconazole-specific in this model. ",
        "Time-varying within subject."
      ),
      source_name        = "FLUC"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 30L,
    n_studies      = 1L,
    age_range      = "0.1-17.7 years (model building); 0.3-18.4 years (validation cohort)",
    age_median     = "5.7 years (model building); 2.0 years (validation)",
    weight_range   = "7.0-77.2 kg (model building); 4.9-63.0 kg (validation)",
    weight_median  = "28.9 kg (model building); 11.2 kg (validation)",
    sex_female_pct = 36.7,
    race_ethnicity = c(White_Caucasian = 93.3, Black_African_American = 3.3, Other = 3.3),
    disease_state  = paste0(
      "Paediatric heart transplant recipients, within the first 6 weeks ",
      "post-transplant; on tacrolimus + mycophenolate immunosuppression ",
      "with milrinone for post-transplant cardiac support."
    ),
    dose_range     = paste0(
      "Oral / enteral immediate-release tacrolimus, ",
      "0.02-0.49 mg/kg/day (model-building cohort; typically split twice ",
      "daily). Administered either orally or via nasogastric / nasojejunal ",
      "tube; the source paper did not detect an administration-route effect."
    ),
    regions        = "Single centre, Primary Children's Hospital, Salt Lake City, Utah, USA",
    n_concentrations_modelbuild = 395L,
    n_concentrations_validation = 330L,
    creatinine_clearance_baseline = "median 122.4 mL/min/1.73 m^2 (range 15.6-442.2; bedside Schwartz)",
    fluconazole_coadministered_pct_modelbuild = 50.0,
    transplant_indications = c(CongenitalHeartDisease = 14L, Cardiomyopathy = 16L, Arrhythmia = 0L),
    notes          = paste0(
      "Retrospective inpatient data, 2007-2015. Model-building cohort = 30 ",
      "(2007-2013); external-validation cohort = 18 (2014-2015). ",
      "Bioanalytical assay: whole-blood tacrolimus by LC-MS/MS at ARUP ",
      "Laboratories, linear range 1-40 ng/mL (= ug/L). Sample times measured ",
      "relative to the first dose of tacrolimus. About 40% of observed ",
      "troughs were within the target therapeutic range of 12-16 ug/L in ",
      "both cohorts (Table 1). CYP3A5 genotype was not available."
    )
  )

  ini({
    # ----- Structural PK (Rower 2017 Table 2 final-model column) -----
    # Source paper used NONMEM ADVAN2 TRANS1 (parameterised on ke and V).
    # Re-expressed here in the canonical (CL/F, V/F) form: CL/F = ke * V.
    # Reference covariates: AGE = 5.7 years, CRCL = 122.4 mL/min/1.73 m^2,
    # CONMED_AZOLE = 0. F/bioavailability was not estimated; defaults to 1.
    lka <- fixed(log(3.43))                                                 ; label("Absorption rate (Ka, 1/h)")  # Rower 2017 Table 2 final model: Ka fixed at 3.43/h (the estimate from the base model, Results paragraph 1 of "Population PK model")
    lcl <- log(0.0408 * 233)                                                ; label("Apparent clearance at reference covariates (CL/F, L/h)")  # CL/F = ke_typ * V_typ = 0.0408 * 233 = 9.5064 L/h; both factors from Rower 2017 Table 2 final-model column
    lvc <- log(233)                                                         ; label("Apparent central volume at reference age (V/F, L)")  # Rower 2017 Table 2 final-model V = 233 L (RSE 17%)

    # ----- Covariate exponents (Rower 2017 Table 2 final-model column) -----
    # Both exponents are estimated point estimates with reported RSE -- not fixed.
    e_age_cl_vc <- 0.775                                                    ; label("Power exponent of AGE on CL/F and V/F (unitless)")  # Rower 2017 Table 2 "Volume: age exponent" = 0.775 (RSE 13%); applies to V only in the source (ke, V) form; appears on CL/F as well in the canonical (CL/F, V/F) form because CL/F = ke * V
    e_crcl_cl   <- 0.850                                                    ; label("Power exponent of CRCL on CL/F (unitless)")  # Rower 2017 Table 2 "Elimination rate: creatinine clearance exponent" = 0.850 (RSE 24%); acts on ke (and hence on CL/F) only; V is not affected by CRCL
    e_azole_cl  <- -0.343                                                   ; label("Fractional change in CL/F from concomitant azole antifungal (unitless)")  # Computed as ke_with_fluc / ke_typ - 1 = 0.0268 / 0.0408 - 1 = -0.343; Rower 2017 Table 2 final-model "Fluconazole elimination rate" 0.0268/h (RSE 5%) vs typical ke 0.0408/h; equivalent to the 34% reduction in elimination reported in Results

    # ----- IIV (Rower 2017 Table 2 final-model column) -----
    # The source paper reports diagonal omegas on the (ke, V) parameterisation:
    # omega^2(log ke) = 0.262 (RSE 40%); omega^2(log V) = 0.329 (RSE 35%).
    # Converting (ke, V) to (CL/F, V/F): log(CL/F) = log(ke) + log(V), so under
    # the source's diagonal assumption (eta_ke independent of eta_V),
    #   etalcl = eta_ke + eta_V  =>  Var(etalcl) = 0.262 + 0.329 = 0.591;
    #   etalvc = eta_V           =>  Var(etalvc) = 0.329;
    #   Cov(etalcl, etalvc) = Var(eta_V) = 0.329.
    # This block correlation reproduces the source's per-subject ke and V
    # distributions exactly (the structural ke and V are preserved; only the
    # decomposition surface changes). Eta-shrinkage in the source was 16% on
    # ke and 14% on V (Discussion paragraph 1).
    etalcl + etalvc ~ c(0.591, 0.329, 0.329)                                # Derived from Rower 2017 Table 2 final-model omega^2(ke) = 0.262 and omega^2(V) = 0.329 (both diagonal); algebra in the block comment above

    # ----- Residual error (Rower 2017 Table 2 final-model column) -----
    # Additive error model selected on the basis of model stability (Results
    # paragraph 1 of "Population PK model").
    addSd <- 3.69                                                           ; label("Additive residual SD (ug/L)")  # Rower 2017 Table 2 final-model additive RUV SD = 3.69 ug/L (RSE 13%)
  })

  model({
    # 1. Derived covariate scalings. Reference values are the model-building
    # cohort medians from Table 1: AGE_ref = 5.7 years, CRCL_ref = 122.4
    # mL/min/1.73 m^2. Methods describes the continuous-covariate equation as
    # theta_pop_star = theta_m * (COV_i / COV_median)^gamma; the fluconazole
    # categorical follows theta_pop_star = theta_pop * (1 + Delta * gamma).
    age_eff   <- (AGE  / 5.7  )^e_age_cl_vc
    crcl_eff  <- (CRCL / 122.4)^e_crcl_cl
    azole_eff <- 1 + e_azole_cl * CONMED_AZOLE

    # 2. Individual PK parameters.
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * crcl_eff * azole_eff * age_eff
    vc <- exp(lvc + etalvc) * age_eff

    # 3. Micro-constants.
    kel <- cl / vc

    # 4. ODE system (one-compartment with first-order oral / enteral absorption).
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # 5. Bioavailability -- F not estimated; defaults to 1 (Methods, PK modelling).

    # 6. Observation and error. Dose in mg, central amount in mg, vc in L
    # gives mg/L; multiply by 1000 to report whole-blood tacrolimus in ug/L
    # (= ng/mL), the bioanalytical assay output (assay linear 1-40 ng/mL;
    # Methods, Data collection).
    Cc <- 1000 * central / vc
    Cc ~ add(addSd)
  })
}
