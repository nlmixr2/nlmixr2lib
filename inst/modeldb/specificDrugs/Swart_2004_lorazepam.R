Swart_2004_lorazepam <- function() {
  description <- "Two-compartment IV population PK model for lorazepam by continuous infusion in mechanically ventilated critically ill adult ICU patients. Clearance is selected by chronic alcohol-abuse status (a flat 0.74 L/h for alcohol-abuse subjects; a PEEP-adjusted 4.13 - (PEEP - 5) * 0.417 L/h otherwise). Steady-state volume of distribution decreases linearly with age above 58 years. Fitted by NONMEM V in the Swart 2004 learning cohort (n = 28)."
  reference   <- "Swart EL, Zuideveld KP, de Jongh J, Danhof M, Thijs LG, Strack van Schijndel RJM. Comparative population pharmacokinetics of lorazepam and midazolam during long-term continuous infusion in critically ill patients. Br J Clin Pharmacol. 2004;57(2):135-145. doi:10.1046/j.1365-2125.2003.01957.x"
  vignette    <- "Swart_2004_lorazepam_midazolam"
  units       <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    PEEP = list(
      description        = "Positive end expiratory pressure applied during positive-pressure mechanical ventilation",
      units              = "mmHg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying within an ICU stay as ventilator settings are adjusted; Swart 2004 Methods lists PEEP among the time-dependent covariates entered against the weighted-residual plots for lorazepam. Learning-group cohort mean 5.3 +/- 2.5 mmHg (range 0-17). Reference / centring value = 5 mmHg (the paper's expression is CL = 4.13 - (PEEP - 5) * 0.417 for the no-alcohol-abuse stratum).",
      source_name        = "PEEP"
    ),
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject. Learning-group cohort mean 58 +/- 17 years (range 21-84). Reference / centring value = 58 years (the paper's expression is Vss = 156 - (age - 58) * 2.07).",
      source_name        = "age"
    ),
    ALCOHOL_ABUSE = list(
      description        = "Chronic alcohol abuse indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no chronic alcohol abuse)",
      notes              = "Defined per Swart 2004 as chronic use of more than 6 units per day (Methods, 'alcohol abuse'). Time-fixed per subject, captured at study or ICU admission. Selector effect on CL: alcohol-abuse subjects take a flat CL = 0.74 L/h with no PEEP effect; non-alcohol-abuse subjects take CL = 4.13 - (PEEP - 5) * 0.417 L/h.",
      source_name        = "alcohol abuse"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 28L,
    n_studies      = 1L,
    age_range      = "18-84 years (learning group); mean 58 +/- 17 (range 21-84)",
    weight_range   = "Mean 80 +/- 25 kg (range 40-175); learning group",
    sex_female_pct = 39.3,
    race_ethnicity = "Not reported (single-centre Dutch ICU cohort, Amsterdam)",
    disease_state  = "Critically ill mechanically-ventilated adult ICU patients (medical ICU); admission diagnostic groups per Table 1: cardio-respiratory insufficiency (19), sepsis (5), postoperative (2), trauma (1), miscellaneous (1). APACHE II on admission mean 18 +/- 7 (range 6-36).",
    dose_range     = "Continuous IV infusion via volumetric pump at 0.16 mg/mL; typical starting rate 2 mL/h, adjusted 0.5-5.3 mL/h (~ 0.08-0.85 mg/h) to the desired Addenbrooke sedation level. Learning-group amount 18 +/- 13.9 mg/day (range 3.8-63.2); mean duration 149 +/- 157 h (24-572 h).",
    regions        = "The Netherlands (Vrije Universiteit Medical Center, Amsterdam)",
    notes          = "Two-treatment, open-label, randomized, parallel-group study (learning group); n = 28 evaluable of 66 initially enrolled with lorazepam or midazolam (17 excluded per Methods). Sampled 0, 15, 30, 45 min and 1, 2, 4, 8 h from infusion start, then every 24 h until infusion stop; post-stop at 0, 15, 30 min, 1, 2, 8, 24 h. Analytical: HPLC-UV, LLOQ 10 ng/mL, linear range 10-1000 ng/mL, inter- / intra-assay CV < 10%. 344 measured concentrations in the learning group. Evaluation cohort (n = 31, 120 concentrations, higher infusion concentration 0.33 mg/mL) is used for external validation of the model reported here; the fitted parameters in the ini() below are from Swart 2004 Table 4 'Model with covariates' column."
  )

  ini({
    # Structural parameters from Swart 2004 Table 4 'Model with covariates'
    # column for lorazepam. The published parameter equations are:
    #   CL (no alcohol abuse) = 4.13 - (PEEP - 5) * 0.417   L/h
    #   CL (alcohol abuse)    = 0.74                        L/h  (flat, no PEEP effect)
    #   V   (central)         = 0.743                       L
    #   Vss (steady-state)    = 156 - (age - 58) * 2.07     L
    #   Q                     = 36.3                        L/h
    # V2 = Vss - V is derived inside model().

    # Baseline CL for the non-alcohol-abuse stratum at the reference PEEP = 5 mmHg.
    # The PEEP effect is an additive linear deviation from this baseline, so lcl
    # here is the log of the intercept at PEEP = 5.
    lcl_noalc <- log(4.13);   label("CL intercept at PEEP = 5, no alcohol abuse (L/h)")   # Swart 2004 Table 4 Model-with-covariates: 4.13

    # Baseline CL for the alcohol-abuse stratum (flat; no PEEP effect).
    lcl_alc   <- log(0.74);   label("CL for alcohol-abuse stratum (L/h)")                 # Swart 2004 Table 4 Model-with-covariates: 0.74

    # Central compartment volume V1 (called V in Table 4).
    lvc       <- log(0.743);  label("Central V1 (L)")                                     # Swart 2004 Table 4 Model-with-covariates: 0.743

    # Peripheral compartment volume V2 = Vss - V1 at the reference AGE = 58 years.
    # From Table 4: Vss(58) = 156 L, V1 = 0.743 L, so V2(58) = 156 - 0.743 = 155.257 L.
    # The age slope on Vss (-2.07 L per year above 58) applies directly to V2 because
    # V1 has no age effect.
    lvp       <- log(155.257); label("Peripheral V2 at AGE = 58 (L)")                     # Swart 2004 Table 4: Vss(58) - V1 = 156 - 0.743 = 155.257

    # Intercompartmental clearance (no covariates in the final model).
    lq        <- log(36.3);   label("Intercompartmental clearance Q (L/h)")               # Swart 2004 Table 4 Model-with-covariates: 36.3

    # Covariate effects (linear-scale slopes; Swart 2004 Table 4 Model-with-covariates).
    e_peep_cl <- 0.417;   label("Additive linear slope of (PEEP - 5) on CL for no-alcohol stratum (L/h per mmHg, subtractive)")  # Swart 2004 Table 4: coefficient in CL = 4.13 - (PEEP - 5) * 0.417
    e_age_vp  <- 2.07;    label("Additive linear slope of (AGE - 58) on Vss / V2 (L per year, subtractive)")                     # Swart 2004 Table 4: coefficient in Vss = 156 - (age - 58) * 2.07

    # Inter-individual variability (Swart 2004 Table 4 Model-with-covariates CV %).
    # The paper assumes a proportional, log-normal distribution for IIV (Methods,
    # step 1). Reported CV (%) is converted to a NONMEM $OMEGA-style variance via
    # omega^2 = ln(1 + CV^2) with CV expressed as a fraction. Table 4 gives two
    # separate CV values for the two alcohol strata (63 % no-alcohol vs 389 %
    # alcohol), so two separate etas are fitted on CL, one per stratum.
    etalcl_noalc ~ 0.3341   # Swart 2004 Table 4 no-alcohol CL CV = 63 %: omega^2 = ln(1 + 0.63^2) = ln(1.3969) = 0.3341
    etalcl_alc   ~ 2.7811   # Swart 2004 Table 4 alcohol-abuse CL CV = 389 %: omega^2 = ln(1 + 3.89^2) = ln(16.1321) = 2.7811
    etalvc       ~ 0.9314   # Swart 2004 Table 4 V CV = 124 %: omega^2 = ln(1 + 1.24^2) = ln(2.5376) = 0.9314
    etalvp       ~ 0.8624   # Swart 2004 Table 4 Vss CV = 117 % (applied to V2 since V1 has its own eta): omega^2 = ln(1 + 1.17^2) = ln(2.3689) = 0.8624
    etalq        ~ 0.4463   # Swart 2004 Table 4 Q CV = 75 %: omega^2 = ln(1 + 0.75^2) = ln(1.5625) = 0.4463

    # Residual error. Swart 2004 Results (Model building for lorazepam): 'A
    # proportional model best described the intra-individual residual
    # variability. The intra-individual residual variability was ... 15 % in the
    # covariate-adjusted model.'
    propSd <- 0.154; label("Proportional residual SD on Cc (fraction)")  # Swart 2004 Table 4 Model-with-covariates: 15.4 %
  })

  model({
    # Individual PK parameters. Two separate etas on CL, one per alcohol-abuse
    # stratum, per Swart 2004 Table 4 (which reports two distinct CV % values,
    # 63 % vs 389 %); a subject enters only the eta of their stratum via the
    # selector. The PEEP effect is an additive linear deviation from the
    # no-alcohol reference at PEEP = 5, and applies only in the no-alcohol
    # stratum. V2 uses the AGE-dependent Vss with the V1 offset absorbed into
    # lvp so exp(lvp) recovers V2(58) = 155.257 L.
    cl_noalc <- (exp(lcl_noalc) - (PEEP - 5) * e_peep_cl) * exp(etalcl_noalc)
    cl_alc   <- exp(lcl_alc + etalcl_alc)
    cl       <- ifelse(ALCOHOL_ABUSE == 1, cl_alc, cl_noalc)
    vc       <- exp(lvc + etalvc)
    vp       <- (exp(lvp) - (AGE - 58) * e_age_vp) * exp(etalvp)
    q        <- exp(lq + etalq)

    # Two-compartment IV disposition micro-constants.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ODE system. Doses (continuous IV infusion) are administered into the
    # central compartment via the data event table (infusion rate and duration).
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # Lorazepam plasma concentration (mg/L; dose in mg, volumes in L).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
