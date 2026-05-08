Archary_2019_abacavir <- function() {
  description <- "Two-compartment population PK model for abacavir in severely malnourished HIV-infected children (Archary 2019); CL/F steps up between day 1 and day 14 of antiretroviral treatment and bioavailability is 31% higher in the early-ART arm"
  reference <- "Archary M, McIlleron H, Bobat R, LaRussa P, Sibaya T, Wiesner L, Hennig S. Population pharmacokinetics of abacavir and lamivudine in severely malnourished human immunodeficiency virus-infected children in relation to treatment outcomes. Br J Clin Pharmacol. 2019;85(8):1881-1890. doi:10.1111/bcp.13998"
  vignette <- "Archary_2019_abacavir"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying; used for allometric scaling on CL/F (exponent 0.75) and Vc/F (exponent 1) with reference weight 7 kg (population median).",
      source_name        = "WT"
    ),
    DAY14 = list(
      description        = "Day-14-post-ART landmark indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (day 1 / pre-rehabilitation)",
      notes              = "Within-subject step indicator: 0 on day 1 of antiretroviral treatment (acute malnutrition baseline), 1 on day 14 (post-nutritional-rehabilitation steady state). Gates the day-1 (3.33) vs day-14 (5.86) typical CL/F per 7 kg step reported in Table 3 of the source.",
      source_name        = "DAY14"
    ),
    EARLY_ART = list(
      description        = "Early-vs-delayed antiretroviral-treatment-initiation arm indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (delayed-ART arm; F = 1 typical-value reference)",
      notes              = "Time-fixed per subject (MATCH-trial randomization). 1 = early ART (initiated within 14 days of admission); 0 = delayed ART (initiated after nutritional recovery, > 14 days). Adds +31% to bioavailability F per Table 3.",
      source_name        = "EARLY"
    )
  )

  population <- list(
    n_subjects     = 75,
    n_studies      = 1,
    age_range      = "0.1-10.8 years (median 1.4)",
    age_median     = "1.4 years",
    weight_range   = "1.88-19.6 kg",
    weight_median  = "7 kg",
    sex_female_pct = 41,
    race_ethnicity = "South African paediatric cohort (race not stratified in source).",
    disease_state  = "Severely malnourished HIV-infected children (weight-for-length Z-scores < -3, mid-upper arm circumference < 115 mm, or peripheral oedema) initiating antiretroviral treatment.",
    dose_range     = "WHO weight-band paediatric oral abacavir + lamivudine + LPV/r dosing; liquid formulation for children < 14 kg, solid formulation > 14 kg (only 2 patients received solid formulations).",
    regions        = "South Africa (King Edward VIII Hospital, Durban).",
    n_observations = 623,
    notes          = "MATCH (Malnutrition and ART Timing in Children with HIV) trial, PACTR201609001751384; 75 patients with day-1 abacavir + lamivudine concentrations, 69 of whom had day-14 samples; 623 abacavir concentrations sampled 0.4-12.4 h post-dose. Demographics summarised in Table 1 of the source."
  )

  ini({
    # Structural parameters; reference weight 7 kg = population median
    lka  <- log(0.97);   label("Absorption rate constant (ka, 1/h)")                           # Table 3
    lcl  <- log(3.33);   label("Apparent clearance day 1 at 7 kg reference (CL/F, L/h)")        # Table 3 (day-1 typical value)
    lvc  <- log(4.63);   label("Apparent central volume at 7 kg reference (Vc/F, L)")           # Table 3
    lq   <- log(0.63);   label("Apparent intercompartmental clearance (Q/F, L/h)")              # Table 3
    lvp  <- log(1.65);   label("Apparent peripheral volume (Vp/F, L)")                          # Table 3
    lfdepot <- fix(log(1)); label("Log baseline bioavailability F for delayed-ART reference (fixed at log(1) = 0)") # Table 3 (delayed-arm typical F = 1)

    # Allometric exponents (paper-fixed per Methods Section 2.3)
    e_wt_cl_q  <- 0.75; label("Allometric exponent on CL/F and Q/F (unitless; fixed)")          # Methods 2.3
    e_wt_vc_vp <- 1;    label("Allometric exponent on Vc/F and Vp/F (unitless; fixed)")         # Methods 2.3

    # DAY14 effect on CL/F (multiplicative additive shift; encodes the day-1 -> day-14 step)
    # 5.86 / 3.33 - 1 = 0.760, i.e. CL/F rises by 76% on day 14 vs day 1
    e_day14_cl <- 0.760; label("Effect of DAY14 (post-ART rehabilitation) on CL/F (fraction)")  # Table 3 (day-14 5.86 vs day-1 3.33 L/h per 7 kg)

    # EARLY_ART effect on bioavailability (additive shift; +31% in early-ART arm)
    e_earlyart_f <- 0.31; label("Effect of EARLY_ART arm on bioavailability F (fraction)")      # Table 3

    # IIV (omega^2 = log(CV^2 + 1))
    # Source reports separate per-arm IIV CL/F (29.2% delayed, 62.5% early) and IOV CL/F = 52%;
    # the model file uses the delayed-arm IIV (29.2%) as a single typical-value omega and
    # documents the early-arm and IOV inflations in the vignette Errata.
    etalcl     ~ 0.0820  # log(1 + 0.292^2); CV 29.2% delayed-arm IIV per Table 3
    # Source reports IIV F = 21.4% only in the early-ART arm; the model file applies it to
    # all subjects (small bias for delayed-arm subjects whose F = 1) and documents this
    # simplification in the vignette Errata.
    etalfdepot ~ 0.0448  # log(1 + 0.214^2); CV 21.4% early-arm IIV F per Table 3

    # Residual error
    propSd <- 0.362; label("Proportional residual error (fraction)")                            # Table 3
  })
  model({
    # PK parameters (allometric weight scaling, day-14 CL step)
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (WT / 7)^e_wt_cl_q  * (1 + e_day14_cl * DAY14)
    vc <- exp(lvc)          * (WT / 7)^e_wt_vc_vp
    vp <- exp(lvp)
    q  <- exp(lq)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Bioavailability (additive shift for early-ART arm; multiplicative log-normal IIV)
    f(depot) <- exp(lfdepot + etalfdepot) * (1 + e_earlyart_f * EARLY_ART)

    # Concentration: dose in mg, volume in L -> mg/L (= ug/mL)
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
