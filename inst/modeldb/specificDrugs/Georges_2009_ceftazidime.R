Georges_2009_ceftazidime <- function() {
  description <- "Two-compartment IV population PK model for ceftazidime in critically ill adults (ICU). Total clearance is an additive linear function of MDRD-estimated glomerular filtration rate; central volume V1 is selected by mechanical-ventilation status; peripheral volume V2 is selected by ICU admission etiology (polytrauma, postsurgical, or medical)."
  reference   <- "Georges B, Conil J-M, Seguin T, Ruiz S, Minville V, Cougot P, Decun J-F, Gonzalez H, Houin G, Fourcade O, Saivin S. Population pharmacokinetics of ceftazidime in intensive care unit patients: influence of glomerular filtration rate, mechanical ventilation, and reason for admission. Antimicrob Agents Chemother. 2009;53(10):4483-4489. doi:10.1128/AAC.00430-09"
  vignette    <- "Georges_2009_ceftazidime"
  units       <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    CRCL = list(
      description        = "MDRD-estimated glomerular filtration rate (raw, not BSA-normalized)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source column MDRD (Modification of Diet in Renal Disease formula). Georges 2009 Methods, Patients section states '... the modification of the diet in renal disease (MDRD) (20, 21) were also calculated.' Population mean MDRD = 121 +/- 55 mL/min (Table 1, total n=72; range across the simulation scenario in Fig. 3 is 30-180 mL/min). The paper does not BSA-normalize MDRD when entering the model and the structural equation 'TVCL = theta1 + theta2 x MDRD, with MDRD in ml/min' treats it as raw mL/min. Stored under the canonical CRCL column per inst/references/covariate-columns.md (CRCL accepts raw mL/min when the source paper does not apply BSA normalization, with the per-model description recording the assay form). The effect on CL is additive linear (TVCL = exp(lcl) + e_crcl_cl * CRCL), no centering or divisive normalization.",
      source_name        = "MDRD"
    ),
    MECH_VENT = list(
      description        = "Mechanical-ventilation indicator at study entry",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not mechanically ventilated)",
      notes              = "Treatment-status flag captured at ICU admission. Georges 2009 Table 1 reports 12/60 no/yes for the total cohort (n=72). Selector effect on V1 in the final model: TVV1 = theta3 when MECH_VENT = 0 (18.9 L) and TVV1 = theta4 when MECH_VENT = 1 (9.02 L). The paper's Discussion attributes the V1 decrease to positive-pressure-ventilation-induced increases in plasma renin activity, aldosterone, and antidiuretic hormone (references 2 and 40 of the source paper).",
      source_name        = "mechanical ventilation"
    ),
    ICU_ADM_POLYTRAUMA = list(
      description        = "ICU admission etiology indicator: polytrauma",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (other admission etiology -- postsurgical or medical)",
      notes              = "One of three mutually-exclusive ICU admission etiology indicators (ICU_ADM_POLYTRAUMA + ICU_ADM_POSTSURG + ICU_ADM_MEDICAL = 1 for every subject). Georges 2009 Table 1 reports 27/72 polytrauma admissions (16 in group 1 model-building + 11 in group 2 validation). Selector effect on V2: TVV2 = theta6 (57.1 L) when ICU_ADM_POLYTRAUMA = 1; the multiplicative form uses ICU_ADM_MEDICAL = 1 as the structural reference (smallest V2 = 13.6 L per Discussion 'patients with a medical reason for admission presented a volume of distribution in the same order of magnitude as those of healthy subjects'). The reference-category indicator ICU_ADM_MEDICAL is not declared in this model's covariateData (it is implicit: ICU_ADM_MEDICAL = 1 - ICU_ADM_POLYTRAUMA - ICU_ADM_POSTSURG) and is registered as a canonical column in inst/references/covariate-columns.md for completeness.",
      source_name        = "Admission for polytrauma"
    ),
    ICU_ADM_POSTSURG = list(
      description        = "ICU admission etiology indicator: postsurgical",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (other admission etiology -- polytrauma or medical)",
      notes              = "One of three mutually-exclusive ICU admission etiology indicators. Georges 2009 Table 1 reports 19/72 postsurgical admissions (16 + 3). Selector effect on V2: TVV2 = theta7 (25.7 L) when ICU_ADM_POSTSURG = 1. See ICU_ADM_POLYTRAUMA notes for the partition convention and the reason ICU_ADM_MEDICAL is not declared here.",
      source_name        = "Admission for postsurgical"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 72L,
    n_studies      = 1L,
    age_range      = "Adults (>=18 years); mean 58 +/- 17 years (Table 1, total n=72)",
    weight_range   = "Mean 76.8 +/- 15.8 kg (Table 1)",
    height_range   = "Mean 172 +/- 7 cm (Table 1)",
    sex_female_pct = 15.3,
    race_ethnicity = "Not reported (French single-centre ICU cohort, Toulouse)",
    disease_state  = "ICU adults with Pseudomonas aeruginosa nosocomial pneumonia or bacteremia presumed sensitive to ceftazidime; mechanically ventilated in 60/72 subjects; admission etiology polytrauma (27), postsurgical (19), or medical (26).",
    dose_range     = "Three regimens: 2 g ceftazidime IV over 30 min q8h (intermittent, n=22); 6 g/day continuous infusion via syringe pump (n=22); 2 g IV loading dose over 30 min followed by 6 g/day continuous infusion (n=28).",
    regions        = "France (Toulouse, single centre, Rangueil University Hospital ICU)",
    severity_scores = "Mean SAPS I 13.9 +/- 4.4, mean SAPS II 47.8 +/- 15.6 (Table 1)",
    renal_function  = "MDRD-eGFR mean 121 +/- 55 mL/min (Table 1; simulation range 30-180 mL/min, Fig. 3)",
    notes          = "Baseline demographics and disease characteristics per Georges 2009 Table 1. Group 1 (n=49, 300 concentrations) used for model-building; group 2 (n=23, 143 concentrations) used for predictive validation. The two groups were then pooled (total n=72, 443 concentrations) to fit the 'final' / 'total' model whose parameter estimates are reproduced in this model file. Per the abstract: 'The mean pharmacokinetic parameters were as follows: CL, 5.48 liters/h, 40%; V1, 10.48 liters, 34%; V2, 32.12 liters, 59%; total volume, 42.60 liters, 45%; and intercompartmental clearance, 16.19 liters/h, 42%' -- these are cohort empirical means with CV%, not the structural thetas (which are encoded below)."
  )

  ini({
    # Structural parameters from Georges 2009 Table 2 'Final model (n = 72)'
    # column. The published parameter equations (Results, 'Population model'
    # section) are:
    #   TVCL = theta1 + theta2 * MDRD              (MDRD in mL/min)
    #   TVV1 = theta3 if MECH_VENT = 0; theta4 if MECH_VENT = 1
    #   TVQ  = theta5
    #   TVV2 = theta6 (polytrauma); theta7 (postsurgical); theta8 (medical)
    # Final-model point estimates: theta1=2.24, theta2=0.024, theta3=18.90,
    # theta4=9.02, theta5=15.20, theta6=57.10, theta7=25.70, theta8=13.60.

    # CL intercept (theta1). Additive linear CL covariate model: the
    # exponentiated lcl is the CL component independent of renal function,
    # and e_crcl_cl is the slope per mL/min of MDRD. Pattern follows the
    # Delattre 2010 amikacin model file in the same registry.
    lcl <- log(2.24);  label("CL intercept theta1 (L/h)")                                  # Georges 2009 Table 2 Final model theta1 = 2.24 L/h

    # V1 reference value at MECH_VENT = 0 (theta3).
    lvc <- log(18.90); label("Central V1 at MECH_VENT = 0 (L)")                            # Georges 2009 Table 2 Final model theta3 = 18.90 L

    # V2 reference value at the ICU_ADM_MEDICAL stratum (theta8); the other
    # two strata enter as positive log-shifts on this reference.
    lvp <- log(13.60); label("Peripheral V2 at ICU_ADM_MEDICAL = 1 (L)")                   # Georges 2009 Table 2 Final model theta8 = 13.60 L

    # Intercompartmental clearance (theta5); no covariates.
    lq  <- log(15.20); label("Intercompartmental clearance Q (L/h)")                       # Georges 2009 Table 2 Final model theta5 = 15.20 L/h

    # Covariate effects (Georges 2009 Table 2 Final model column).
    e_crcl_cl              <- 0.024;             label("Slope of additive linear MDRD effect on CL (L/h per mL/min)")  # Georges 2009 Table 2 Final model theta2 = 0.024
    e_mech_vent_vc         <- log(9.02 / 18.90); label("Log effect of MECH_VENT = 1 on V1 (= log(theta4 / theta3))")   # Georges 2009 Table 2 Final model: theta4 = 9.02 L vs theta3 = 18.90 L
    e_icu_adm_polytrauma_vp <- log(57.10 / 13.60); label("Log effect of ICU_ADM_POLYTRAUMA = 1 on V2 (= log(theta6 / theta8))") # Georges 2009 Table 2 Final model: theta6 = 57.10 L vs theta8 = 13.60 L
    e_icu_adm_postsurg_vp   <- log(25.70 / 13.60); label("Log effect of ICU_ADM_POSTSURG = 1 on V2 (= log(theta7 / theta8))")   # Georges 2009 Table 2 Final model: theta7 = 25.70 L vs theta8 = 13.60 L

    # Inter-individual variability (Georges 2009 Table 2 Final model
    # omegas). Values are reported as exponential-IIV variances on the
    # multiplicative log-normal eta (NONMEM $OMEGA block); CV(%) values
    # in Table 2 parentheses are the RSE of the omega estimate, not a
    # transformed CV of the parameter.
    etalcl ~ 0.09  # Georges 2009 Table 2 Final model Omega CL = 0.09 -> CV(CL) = sqrt(exp(0.09) - 1) = 30.7%
    etalvc ~ 0.12  # Georges 2009 Table 2 Final model Omega V1 = 0.12 -> CV(V1) = sqrt(exp(0.12) - 1) = 35.7%
    etalvp ~ 0.11  # Georges 2009 Table 2 Final model Omega V2 = 0.11 -> CV(V2) = sqrt(exp(0.11) - 1) = 34.1%
    etalq  ~ 0.50  # Georges 2009 Table 2 Final model Omega Q  = 0.50 -> CV(Q)  = sqrt(exp(0.50) - 1) = 80.5%

    # Residual error: proportional only (Georges 2009 Results: 'A
    # proportional-error model was the most accurate for residual and
    # interpatient variability.'). Table 2 final-model Sigma = 0.05
    # is the NONMEM $SIGMA variance; the linear-scale proportional SD
    # is sqrt(0.05) ~= 0.2236 (i.e., ~ 22.4% proportional CV on Cc).
    propSd <- sqrt(0.05); label("Proportional residual SD on Cc (fraction)")  # Georges 2009 Table 2 Final model Sigma = 0.05 -> propSd = sqrt(0.05) = 0.2236
  })

  model({
    # Individual PK parameters. CL has an additive linear covariate term on
    # raw (non-BSA-normalized) MDRD-eGFR per the paper's structural form
    # 'TVCL = theta1 + theta2 * MDRD'; the multiplicative log-normal eta is
    # applied to the combined intercept + slope, matching the Delattre
    # 2010 amikacin precedent. V1 and V2 use the conventional baseline-
    # plus-log-effect parameterization with reference categories
    # MECH_VENT = 0 (V1) and ICU_ADM_MEDICAL = 1 (V2); on those references
    # exp(lvc) and exp(lvp) recover theta3 = 18.90 L and theta8 = 13.60 L
    # exactly, and the positive log-shifts recover the larger thetas for
    # ventilated patients (theta4 on V1) and polytrauma / postsurgical
    # patients (theta6 / theta7 on V2).
    cl <- (exp(lcl) + e_crcl_cl * CRCL) * exp(etalcl)
    vc <- exp(lvc + e_mech_vent_vc * MECH_VENT + etalvc)
    vp <- exp(lvp + e_icu_adm_polytrauma_vp * ICU_ADM_POLYTRAUMA +
                    e_icu_adm_postsurg_vp   * ICU_ADM_POSTSURG +
                    etalvp)
    q  <- exp(lq  + etalq)

    # Two-compartment IV disposition micro-constants.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ODE system. Doses are administered into the central compartment
    # via the data event table (30-minute IV infusions for intermittent
    # 2 g doses; continuous-infusion rate for the 6 g/day arms).
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # Ceftazidime serum concentration in mg/L (dose in mg, volumes in L).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
