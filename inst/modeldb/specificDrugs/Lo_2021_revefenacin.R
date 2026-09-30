Lo_2021_revefenacin <- function() {
  description <- paste(
    "Parent + metabolite population pharmacokinetic model for nebulized",
    "revefenacin (a lung-selective long-acting muscarinic antagonist) and its",
    "major hydrolysis metabolite THRX-195518 in 935 patients with chronic",
    "obstructive pulmonary disease from three phase II and two phase III",
    "studies (Lo 2021). Revefenacin is described by a two-compartment model",
    "with first-order absorption (ka fixed at 200 1/h) from a dosing depot",
    "representing the lung, with relative bioavailability depending on the",
    "dose level (power) and lower in phase II Study 0059. THRX-195518 is",
    "formed from a fixed 21% of the individual revefenacin clearance and is",
    "described by a two-compartment model whose central volume shares the",
    "metabolite-clearance random effect through a scale factor. Retained",
    "covariates: age on revefenacin CL/F, body weight on revefenacin Q/F, age",
    "on THRX-195518 CL/F and body weight on the formed fraction. The paper",
    "fitted the metabolite sequentially on the parent's post hoc estimates;",
    "this file couples the two analytes in one rxode2 model. Residual-error",
    "magnitudes are not reported and are held at zero."
  )
  reference <- paste(
    "Lo A, Borin MT, Bourdet DL. Population Pharmacokinetics of Revefenacin",
    "in Patients with Chronic Obstructive Pulmonary Disease.",
    "Clin Pharmacokinet. 2021;60(3):391-401. doi:10.1007/s40262-020-00938-3.",
    "Parameter estimates are from Table 2; covariate equations from Methods",
    "2.4; covariate medians (age 64 years, weight 81 kg) from Results 3.4.",
    "The omega variances, the population median weight (81.2 kg) and the",
    "shared-eta form of the CLmet/V3 'correlation' are confirmed against the",
    "FDA Clinical Pharmacology Review of NDA 210598 (Yupelri, 2018),",
    "Tables 4.1.2.2.2 and 4.1.2.4.1."
  )
  vignette <- "Lo_2021_revefenacin"
  units <- list(time = "h", dosing = "ug", concentration = "ng/mL")

  covariateData <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power effect normalized to the population median of 64 years (Lo",
        "2021 Results 3.4; FDA review Table 4.1.2.2.2 median 64.0) on both",
        "revefenacin CL/F and THRX-195518 CL/F. Analysis-population range",
        "41-88 years, mean 63.5 (SD 8.72) (Table 1)."
      ),
      source_name = "AGE"
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power effect normalized to the population median on revefenacin Q/F",
        "and on the fraction of revefenacin clearance forming THRX-195518.",
        "Lo 2021 quotes the median as 81 kg (Results 3.4); the FDA review",
        "(Table 4.1.2.2.2) prints it as 81.2 kg, which is used here.",
        "Analysis-population range 38.5-192 kg, mean 83.3 (SD 21.8)."
      ),
      source_name = "WT"
    ),
    DOSE_REVEFENACIN_UG = list(
      description = "Nominal nebulized revefenacin dose of the current administration",
      units = "ug",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power effect on relative bioavailability, (DOSE / 175)^0.0987. The",
        "paper normalizes continuous covariates to the population median but",
        "does not print the median dose; 175 ug is assumed (see vignette",
        "Assumptions and deviations). Must equal the amt of the dosing",
        "record. Studied levels 22, 44, 88, 175, 350 and 700 ug."
      ),
      source_name = "Dose"
    ),
    STUDY_0059 = list(
      description = "Phase II Study 0059 (Lo 2021 'Study 1', NCT03064113) indicator; 1 = Study 0059, 0 = otherwise",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (all other studies; use 0 to simulate the marketed product)",
      notes = paste(
        "Multiplies relative bioavailability by 0.553 (categorical form",
        "theta_eff^K_ind, Methods 2.4) to reflect the lower exposures",
        "observed in that single-dose crossover study; the FDA reviewer",
        "suggests a different nebulizer as a possible cause."
      ),
      source_name = "Study 1"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Biological sex, 1 = female",
      units = "(binary)",
      type = "binary",
      notes = "Tested in the stepwise covariate analysis (Methods 2.4) and not retained."
    ),
    CRCL = list(
      description = "Creatinine clearance",
      units = "mL/min",
      type = "continuous",
      notes = "Tested in the stepwise covariate analysis (Methods 2.4) and not retained; range 22-151 mL/min."
    ),
    SMOKE = list(
      description = "Current smoker",
      units = "(binary)",
      type = "binary",
      notes = "Tested in the stepwise covariate analysis and not retained; 46% current smokers."
    ),
    FEV1_BL = list(
      description = "Baseline forced expiratory volume in 1 s",
      units = "mL",
      type = "continuous",
      notes = "Tested on CL/F, V1/F, Q/F and V2/F and not retained."
    )
  )

  compartmentData <- list(
    depot = list(analyte = "revefenacin", units = "ug", specimen = "administration site", verified = TRUE),
    central = list(analyte = "revefenacin", units = "ug", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "revefenacin", units = "ug", specimen = "plasma", verified = TRUE),
    central_thrx195518 = list(
      analyte = "THRX-195518 (revefenacin metabolite)",
      units = "ug",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1_thrx195518 = list(
      analyte = "THRX-195518 (revefenacin metabolite)",
      units = "ug",
      specimen = "plasma",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 935L,
    n_studies = 5L,
    n_observations = "10043 revefenacin and 10717 THRX-195518 measurable plasma concentrations (Results 3)",
    age_range = "41-88 years",
    age_median = "64 years",
    weight_range = "38.5-192 kg",
    weight_median = "81.2 kg",
    sex_female_pct = 47.8,
    race_ethnicity = c(White = 90.3),
    disease_state = "Moderate to very severe chronic obstructive pulmonary disease",
    dose_range = "22-700 ug once daily by jet nebulizer for 1, 7 or 28 days (phase II); 88 or 175 ug once daily for 12 weeks (phase III)",
    renal_function = "Estimated creatinine clearance 22-151 mL/min, mean 71.7 (SD 20.7)",
    co_medication = "32.5% concomitant LABA/ICS therapy",
    regions = "Multinational (USA, New Zealand, South Africa, UK and others)",
    notes = paste(
      "Baseline demographics from Lo 2021 Table 1 and Results 3.1 (488 men,",
      "447 women). Studies: 0059 (Study 1, single dose 350/700 ug, n = 32),",
      "0091 (Study 2, 22-700 ug for 7 days, n = 61), 0117 (Study 3, 44-350",
      "ug for 28 days, n = 34 with PK), 0126 and 0127 (Studies 4 and 5, phase",
      "III, 88 or 175 ug for 12 weeks, n = 808 with PK). Patients with",
      "moderate to severe hepatic impairment were excluded."
    )
  )

  ini({
    # ==================================================================
    # REVEFENACIN (parent). Lo 2021 Table 2. All parameters are apparent
    # (per unit of the reference relative bioavailability) because every
    # study dosed the nebulized inhalation solution.
    # ==================================================================
    lcl <- log(668)
    label("Revefenacin apparent clearance CL/F (L/h)") # Table 2 CL/F = 668 L/h (RSE 3.17%)
    lvc <- log(867)
    label("Revefenacin apparent central volume V1/F (L)") # Table 2 V1/F = 867 L (RSE 3.77%)
    lq <- log(2607)
    label("Revefenacin apparent intercompartmental clearance Q/F (L/h)") # Table 2 Q/F = 2607 L/h (RSE 2.51%)
    lvp <- log(15495)
    label("Revefenacin apparent peripheral volume V2/F (L)") # Table 2 V2/F = 15,495 L (RSE 4.88%)
    lka <- fixed(log(200))
    label("First-order absorption rate constant from the lung depot ka (1/h)") # Table 2 Ka = 200 (fixed); units printed 'L/h', Results 3.2 and 3.5 give 200/h

    # Relative bioavailability is anchored at 1 for the reference dose (175
    # ug) outside Study 0059; only the covariate terms and the eta move it.
    lfdepot <- fixed(log(1))
    label("Relative bioavailability of the nebulized dose at the reference dose (unitless)") # Results 3.2: F1 carries the Study 1 and dose covariate terms only; reference value not estimated

    e_study0059_fdepot <- 0.553
    label("Multiplicative factor on relative bioavailability for Study 0059 (unitless)") # Table 2 'Study 1 effect on F1' = 0.553 (RSE 7.09%)
    e_dose_fdepot <- 0.0987
    label("Power exponent of dose / 175 ug on relative bioavailability (unitless)") # Table 2 'Dose effect on F1' = 0.0987 (RSE 4.80%)
    e_age_cl <- -0.559
    label("Power exponent of age / 64 years on revefenacin CL/F (unitless)") # Table 2 'Age effect on CL/F' = -0.559 (RSE 26.0%)
    e_wt_q <- 0.485
    label("Power exponent of weight / 81.2 kg on revefenacin Q/F (unitless)") # Table 2 'Weight effect on Q/F' = 0.485 (RSE 12.4%)

    # Table 2 prints IIV as omega SD x 100 (the NONMEM '%CV' convention):
    # the FDA review Table 4.1.2.4.1 lists the variances 0.316, 0.0722,
    # 0.0962, 0.272 and 0.114, which are exactly the squares of 0.562,
    # 0.269, 0.310, 0.522 and 0.337.
    etalcl ~ 0.316 # Table 2 IIV 56.2%; FDA review 'Revefenacin Apparent Clearance Variance' = 0.316
    etalvc ~ 0.0722 # Table 2 IIV 26.9%; FDA review 'Apparent Central Vd Variance' = 0.0722
    etalq ~ 0.0962 # Table 2 IIV 31.0%; FDA review 'Intercompartmental Clearance Variance' = 0.0962
    etalvp ~ 0.272 # Table 2 IIV 52.2%; FDA review 'Apparent Peripheral Vd Variance' = 0.272
    etalfdepot ~ 0.114 # Table 2 IIV 33.7% (row 'Study 1 effect on F1'; Results 3.2 places it on F1); FDA review 'Bioavailability Variance' = 0.114

    # ==================================================================
    # THRX-195518 (metabolite). Lo 2021 Table 2.
    # ==================================================================
    lcl_thrx195518 <- log(53.2)
    label("THRX-195518 apparent clearance CLmet/F (L/h)") # Table 2 CLmet/F = 53.2 L/h (RSE 1.84%)
    lvc_thrx195518 <- log(20.4)
    label("THRX-195518 apparent central volume V3/F (L)") # Table 2 V3/F = 20.4 L (RSE 3.62%)
    lq_thrx195518 <- log(36.3)
    label("THRX-195518 apparent intercompartmental clearance Qmet/F (L/h)") # Table 2 Qmet/F = 36.3 L/h (RSE 4.31%)
    lvp_thrx195518 <- log(35.8)
    label("THRX-195518 apparent peripheral volume V4/F (L)") # Table 2 V4/F = 35.8 L (RSE 5.44%)
    fm <- fixed(0.21)
    label("Fraction of revefenacin clearance forming THRX-195518 at the reference weight (unitless)") # Table 2 Fmet = 0.21 (fixed), from the human ADME study recovery (Results 3.3)

    # Table 2 'Correlation between CLmet and V3' = 1.45 is not a correlation
    # coefficient but the scale of a shared random effect: V3 carries
    # 1.45 x eta(CLmet). 1.45 x 36.0% = 52.2%, the IIV Table 2 prints for V3,
    # and the FDA review describes it as 'an additional THETA term to
    # represent this correlation'.
    vc_thrx195518_eta_scale <- 1.45
    label("Scale of the THRX-195518 clearance random effect carried by V3/F (unitless)") # Table 2 'Correlation between CLmet and V3' = 1.45 (RSE 7.46%)
    e_age_cl_thrx195518 <- -0.777
    label("Power exponent of age / 64 years on THRX-195518 CLmet/F (unitless)") # Table 2 'Age effect on CLmet/F' = -0.777 (RSE 13.1%)
    e_wt_fm <- -0.406
    label("Power exponent of weight / 81.2 kg on the fraction metabolized Fmet (unitless)") # Table 2 'Weight effect on Fmet' = -0.406 (RSE 18.1%)

    etalcl_thrx195518 ~ 0.13 # Table 2 IIV 36.0%; FDA review 'Metabolite Apparent Clearance Variance' = 0.13

    # ==================================================================
    # RESIDUAL ERROR. Results 3.2 / 3.3: combined additive + proportional
    # error for the phase II data and a separate proportional error for the
    # phase III data, for each analyte. Neither the paper, its supplement
    # nor the FDA review prints the magnitudes, so they are held at zero;
    # the model therefore simulates individual predictions only.
    # ==================================================================
    propSd <- fixed(0)
    label("Revefenacin proportional residual SD (fraction; not reported, held at zero)") # not reported in Lo 2021, its ESM or the FDA review
    addSd <- fixed(0)
    label("Revefenacin additive residual SD (ng/mL; not reported, held at zero)") # not reported in Lo 2021, its ESM or the FDA review
    propSd_thrx195518 <- fixed(0)
    label("THRX-195518 proportional residual SD (fraction; not reported, held at zero)") # not reported in Lo 2021, its ESM or the FDA review
    addSd_thrx195518 <- fixed(0)
    label("THRX-195518 additive residual SD (ng/mL; not reported, held at zero)") # not reported in Lo 2021, its ESM or the FDA review
  })
  model({
    # Covariate equation (Methods 2.4): theta_i = theta_typ x (Cov / median)^theta_eff
    # for continuous covariates and theta_i = theta_typ x theta_eff^K for binary.
    age_ref <- AGE / 64
    wt_ref <- WT / 81.2
    dose_ref <- DOSE_REVEFENACIN_UG / 175

    cl <- exp(lcl + etalcl) * age_ref^e_age_cl
    vc <- exp(lvc + etalvc)
    q <- exp(lq + etalq) * wt_ref^e_wt_q
    vp <- exp(lvp + etalvp)
    ka <- exp(lka)
    fdepot <- exp(lfdepot + etalfdepot) * dose_ref^e_dose_fdepot * e_study0059_fdepot^STUDY_0059

    cl_thrx195518 <- exp(lcl_thrx195518 + etalcl_thrx195518) * age_ref^e_age_cl_thrx195518
    vc_thrx195518 <- exp(lvc_thrx195518 + vc_thrx195518_eta_scale * etalcl_thrx195518)
    q_thrx195518 <- exp(lq_thrx195518)
    vp_thrx195518 <- exp(lvp_thrx195518)
    fmet <- fm * wt_ref^e_wt_fm

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    kel_thrx195518 <- cl_thrx195518 / vc_thrx195518
    k12_thrx195518 <- q_thrx195518 / vc_thrx195518
    k21_thrx195518 <- q_thrx195518 / vp_thrx195518

    # Fig. 1: Fmet x CL/F of the eliminated revefenacin forms THRX-195518 and
    # (1 - Fmet) x CL/F leaves by other routes. Formation is mass for mass
    # (the amide-to-acid hydrolysis changes molecular weight by < 0.2%).
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(central_thrx195518) <- fmet * kel * central - kel_thrx195518 * central_thrx195518 -
      k12_thrx195518 * central_thrx195518 + k21_thrx195518 * peripheral1_thrx195518
    d/dt(peripheral1_thrx195518) <- k12_thrx195518 * central_thrx195518 - k21_thrx195518 * peripheral1_thrx195518

    f(depot) <- fdepot

    # ug / L = ng / mL
    Cc <- central / vc
    Cc_thrx195518 <- central_thrx195518 / vc_thrx195518

    Cc ~ add(addSd) + prop(propSd)
    Cc_thrx195518 ~ add(addSd_thrx195518) + prop(propSd_thrx195518)
  })
}
