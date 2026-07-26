Swart_2004_midazolam <- function() {
  description <- "Two-compartment IV population PK model for midazolam by continuous infusion in mechanically ventilated critically ill adult ICU patients. Clearance is selected by chronic alcohol-abuse status and decreases linearly with age above 57 years; intercompartmental clearance decreases linearly with APACHE II score above 26. Fitted by NONMEM V in the Swart 2004 learning cohort (n = 21)."
  reference   <- "Swart EL, Zuideveld KP, de Jongh J, Danhof M, Thijs LG, Strack van Schijndel RJM. Comparative population pharmacokinetics of lorazepam and midazolam during long-term continuous infusion in critically ill patients. Br J Clin Pharmacol. 2004;57(2):135-145. doi:10.1046/j.1365-2125.2003.01957.x"
  vignette    <- "Swart_2004_lorazepam_midazolam"
  units       <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject. Learning-group cohort mean 57 +/- 16 years (range 21-84). Reference / centring value = 57 years (the paper's CL expressions are 11.3 - (age - 57) * 0.145 (no alcohol) and 7.3 - (age - 57) * 0.145 (alcohol abuse)).",
      source_name        = "age"
    ),
    APACHE_II = list(
      description        = "Acute Physiology and Chronic Health Evaluation II score at ICU admission",
      units              = "points",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject; captured within the first 24 h of ICU admission per Knaus 1985. Learning-group cohort mean 26 +/- 9 points (range 6-34). Reference / centring value = 26 points (the paper's expression is Q = 40.8 - (APACHE - 26) * 2.75).",
      source_name        = "APACHE"
    ),
    ALCOHOL_ABUSE = list(
      description        = "Chronic alcohol abuse indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no chronic alcohol abuse)",
      notes              = "Defined per Swart 2004 as chronic use of more than 6 units per day (Methods, 'alcohol abuse'). Time-fixed per subject, captured at study or ICU admission. Selector effect on CL: alcohol-abuse subjects take CL = 7.3 - (AGE - 57) * 0.145 L/h; non-alcohol-abuse subjects take CL = 11.3 - (AGE - 57) * 0.145 L/h (same age slope, ~4 L/h lower baseline).",
      source_name        = "alcohol abuse"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 21L,
    n_studies      = 1L,
    age_range      = "18-84 years (learning group); mean 57 +/- 16 (range 21-84)",
    weight_range   = "Mean 71 +/- 13 kg (range 40-90); learning group",
    sex_female_pct = 38.1,
    race_ethnicity = "Not reported (single-centre Dutch ICU cohort, Amsterdam)",
    disease_state  = "Critically ill mechanically-ventilated adult ICU patients (medical ICU); admission diagnostic groups per Table 1: cardio-respiratory insufficiency (14), sepsis (6), postoperative (1), trauma (0), miscellaneous (0). APACHE II on admission mean 26 +/- 9 (range 6-34).",
    dose_range     = "Continuous IV infusion via volumetric pump at 5 mg/mL undiluted; typical starting rate 2 mL/h (10 mg/h), adjusted 0.5-16 mL/h (~ 2.5-80 mg/h) to the desired Addenbrooke sedation level. Learning-group amount 497 +/- 417 mg/day (range 119-1682); mean duration 134 +/- 172 h (24-715 h).",
    regions        = "The Netherlands (Vrije Universiteit Medical Center, Amsterdam)",
    notes          = "Two-treatment, open-label, randomized, parallel-group study (learning group); n = 21 evaluable of 66 initially enrolled with lorazepam or midazolam (17 excluded per Methods). Sampling schedule matches the lorazepam arm. Analytical: HPLC-UV, LLOQ 10 ng/mL, linear range 10-10000 ng/mL, inter- / intra-assay CV < 10%. 494 measured concentrations in the learning group. Evaluation cohort (n = 33, higher rate of interacting CYP3A4-inhibitor comedication than the learning group) is used for external validation of the model reported here; the fitted parameters in the ini() below are from Swart 2004 Table 5 'Model with covariates' column. 1-hydroxymidazolam concentrations were also measured but 'No model could be identified to describe the pharmacokinetics of 1-hydroxymidazolam, because only a few data points were available and the concentrations of this metabolite were low' (Results); this file therefore encodes only the parent midazolam PK."
  )

  ini({
    # Structural parameters from Swart 2004 Table 5 'Model with covariates'
    # column for midazolam. The published parameter equations are:
    #   CL (no alcohol abuse) = 11.3 - (age - 57) * 0.145        L/h
    #   CL (alcohol abuse)    = 7.3  - (age - 57) * 0.145        L/h
    #   V   (central)         = 7.15                             L
    #   Vss (steady-state)    = 431                              L  (no age effect)
    #   Q                     = 40.8 - (APACHE - 26) * 2.75      L/h
    # V2 = Vss - V is derived inside model().

    # Baseline CL for the non-alcohol-abuse stratum at the reference AGE = 57.
    lcl_noalc <- log(11.3);  label("CL intercept at AGE = 57, no alcohol abuse (L/h)")   # Swart 2004 Table 5 Model-with-covariates: 11.3

    # Baseline CL for the alcohol-abuse stratum at the reference AGE = 57 (same age slope).
    lcl_alc   <- log(7.3);   label("CL intercept at AGE = 57, alcohol abuse (L/h)")      # Swart 2004 Table 5 Model-with-covariates: 7.3

    # Central compartment volume V1 (called V in Table 5).
    lvc       <- log(7.15);  label("Central V1 (L)")                                     # Swart 2004 Table 5 Model-with-covariates: 7.15

    # Peripheral compartment volume V2 = Vss - V1. Table 5: Vss = 431 L, V1 = 7.15 L,
    # V2 = 431 - 7.15 = 423.85 L. Vss has no age dependence in the midazolam model.
    lvp       <- log(423.85); label("Peripheral V2 (L)")                                 # Swart 2004 Table 5: Vss - V1 = 431 - 7.15 = 423.85

    # Baseline Q at the reference APACHE_II = 26. The APACHE effect is an additive
    # linear deviation from this baseline.
    lq        <- log(40.8);  label("Intercompartmental clearance Q at APACHE_II = 26 (L/h)")  # Swart 2004 Table 5 Model-with-covariates: 40.8

    # Covariate effects (linear-scale slopes; Swart 2004 Table 5 Model-with-covariates).
    e_age_cl     <- 0.145; label("Additive linear slope of (AGE - 57) on CL (L/h per year, subtractive; same slope in both alcohol strata)")  # Swart 2004 Table 5: coefficient 0.145 in CL formulas for both alcohol strata
    e_apache_q   <- 2.75;  label("Additive linear slope of (APACHE_II - 26) on Q (L/h per point, subtractive)")                              # Swart 2004 Table 5: coefficient in Q = 40.8 - (APACHE - 26) * 2.75

    # Inter-individual variability (Swart 2004 Table 5 Model-with-covariates CV %).
    # The paper assumes a proportional, log-normal distribution for IIV (Methods,
    # step 1). Reported CV (%) is converted via omega^2 = ln(1 + CV^2) with CV
    # expressed as a fraction. Table 5 gives a single CV (%) = 77 for both alcohol
    # strata of CL (Table 5 CL rows both list '77'), so a single eta on CL is fitted
    # (in contrast to Swart 2004 lorazepam where two CL etas are needed).
    etalcl ~ 0.4654   # Swart 2004 Table 5 CL CV = 77 % (same in both alcohol strata): omega^2 = ln(1 + 0.77^2) = ln(1.5929) = 0.4654
    etalvc ~ 1.8023   # Swart 2004 Table 5 V CV = 225 %: omega^2 = ln(1 + 2.25^2) = ln(6.0625) = 1.8023
    etalvp ~ 2.6103   # Swart 2004 Table 5 Vss CV = 355 % (applied to V2 since V1 has its own eta): omega^2 = ln(1 + 3.55^2) = ln(13.6025) = 2.6103
    etalq  ~ 0.2561   # Swart 2004 Table 5 Q CV = 54 %: omega^2 = ln(1 + 0.54^2) = ln(1.2916) = 0.2561

    # Residual error. Swart 2004 Results (Model building for midazolam): 'An additive
    # plus proportional model best described the intra-individual residual
    # variability, values of which were ... 31 % plus 32 ng/ml in the covariate-
    # adjusted model.' 32 ng/mL = 0.032 mg/L (matching the file's mg/L concentration
    # unit).
    propSd <- 0.309; label("Proportional residual SD on Cc (fraction)")   # Swart 2004 Table 5 Model-with-covariates: 30.9 %
    addSd  <- 0.032; label("Additive residual SD on Cc (mg/L)")           # Swart 2004 Table 5 Model-with-covariates: 32 ng/mL = 0.032 mg/L
  })

  model({
    # Individual PK parameters. A single eta on CL applies to both alcohol strata
    # (Table 5 reports the same 77 % CV for both strata). The alcohol-abuse
    # selector picks the correct intercept, and the AGE effect uses the same
    # slope in both strata. The APACHE_II effect is an additive linear deviation
    # on Q from the baseline at APACHE_II = 26.
    cl_noalc <- (exp(lcl_noalc) - (AGE - 57) * e_age_cl)
    cl_alc   <- (exp(lcl_alc)   - (AGE - 57) * e_age_cl)
    cl       <- ifelse(ALCOHOL_ABUSE == 1, cl_alc, cl_noalc) * exp(etalcl)
    vc       <- exp(lvc + etalvc)
    vp       <- exp(lvp + etalvp)
    q        <- (exp(lq) - (APACHE_II - 26) * e_apache_q) * exp(etalq)

    # Two-compartment IV disposition micro-constants.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ODE system. Doses (continuous IV infusion) are administered into the
    # central compartment via the data event table (infusion rate and duration).
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # Midazolam plasma concentration (mg/L; dose in mg, volumes in L).
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
