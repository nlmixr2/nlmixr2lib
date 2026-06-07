Hussein_2005_levocetirizine <- function() {
  description <- "One-compartment population PK model with first-order absorption and first-order elimination for orally administered levocetirizine in atopic young children (12-48 months, 8-20 kg) receiving 0.125 mg/kg twice-daily levocetirizine (administered as 0.25 mg/kg twice-daily racemic cetirizine) for 18 months in the ETAC study (Hussein 2005). CL/F and V/F are linear functions of body weight (CL/F = 0.244 + 0.0442 * WT L/h; V/F = 0.639 * WT L). The absorption rate constant ka is parameterised as ka = theta_ka + CL/V to guard against flip-flop kinetics, with theta_ka = 1.140 1/h and CL/V contributing on average less than 5% to ka. Residual variability is additive with two concentration-dependent magnitudes: 53.5 ng/mL for Cc <= 400 ng/mL and 316 ng/mL for Cc > 400 ng/mL (the 400 ng/mL threshold was selected by sensitivity analysis and has no clinical or therapeutic implication). Bioavailability is anchored at F = 1 here; the paper additionally estimated F_noncomp = 0.281 applied to 12% of records flagged as suspected noncompliance and recorded in the vignette Assumptions and deviations."
  reference <- "Hussein Z, Pitsiu M, Majid O, Aarons L, de Longueville M, Stockis A, on behalf of the ETAC Study Group. Retrospective population pharmacokinetics of levocetirizine in atopic children receiving cetirizine: the ETAC study. Br J Clin Pharmacol. 2005;59(1):28-37. doi:10.1111/j.1365-2125.2005.02242.x"
  vignette <- "Hussein_2005_levocetirizine"
  paper_specific_residual_sds <- c("addSdLow", "addSdHigh")
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight at the visit (time-varying within subject across the 18-month follow-up)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Linear effects on CL/F and V/F in the final model (Hussein 2005 Table 3): CL/F = 0.244 + 0.0442 * WT (L/h), V/F = 0.639 * WT (L). The covariate is time-varying (Hussein 2005 Methods, weight recorded at each visit; visit means 11.9 / 13.9 / 15.0 kg at months 3 / 12 / 18). Range across all visits 8.2-20.5 kg.",
      source_name        = "WT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 343L,
    n_studies      = 1L,
    age_range      = "12-48 months (14.4-46.3 months across the three visits)",
    age_median     = "Visit-by-visit medians 19.1 / 27.8 / 34.6 months at months 3 / 12 / 18",
    weight_range   = "8.2-20.5 kg across visits",
    weight_median  = "Visit-by-visit medians 11.8 / 13.7 / 14.8 kg at months 3 / 12 / 18",
    sex_female_pct = NA_real_,
    race_ethnicity = "Not reported",
    disease_state  = "Atopic young children at high risk of developing asthma but not yet affected; enrolled in the ETAC (Early Treatment of the Atopic Child) randomized, double-blind, parallel-group, placebo-controlled trial. Concomitant exposures recorded at each visit included corticosteroids, penicillins, macrolides, hydroxyzine, and diarrhoea / gastro-enteritis; allergic sensitization (aeroallergen IgE) and severe-allergy markers (eosinophil count) were also captured. Severe allergy, aeroallergen sensitization, gender, age, BSA, and creatinine clearance were screened as covariates but did not reach significance in the final model (Hussein 2005 Tables 2-3 and Results).",
    dose_range     = "Cetirizine dihydrochloride (50% levocetirizine + 50% dextrocetirizine) 0.25 mg/kg twice daily PO for 18 months as oral drops (10 mg/mL formulation). Modelled dose is the levocetirizine enantiomer, 0.125 mg/kg twice daily.",
    regions        = "Multi-centre European trial (sites in the ETAC programme)",
    samples        = "753 plasma levocetirizine concentration records used in the population analysis (66 children with 1 record, 144 with 2, 133 with 3); ~75% of samples taken 1-5 h post-dose. 3.9% of measured concentrations were excluded for sampling time < 0 h or > 12 h relative to dose, or concentration <= 12 ng/mL (dilution-corrected lower limit of determination).",
    follow_up      = "18-month treatment period with PK sampling at months 3, 12, and 18 (visits 3, 6, 8); plasma drawn at each visit with the date and hour of the last dose recorded.",
    estimation     = "NONMEM V level 1.1; first-order (FO) estimation method (FOCE also evaluated but did not improve precision).",
    notes          = "Hussein 2005 Br J Clin Pharmacol 59(1):28-37, doi:10.1111/j.1365-2125.2005.02242.x. Baseline demographics in Table 1; final-model parameter estimates in Table 3. ETAC enrolled 343 of 399 cetirizine-treated children with usable plasma data; visit-3 / visit-6 / visit-8 included 248 / 279 / 226 children. Levocetirizine plasma assay: chiral HPLC with tandem MS detection; lower limit of determination 12 ng/mL after 4-fold plasma dilution."
  )

  ini({
    # Structural fixed effects from Hussein 2005 Table 3 (final-model
    # estimates). The reported relationships are linear in body weight:
    #   CL/F (L/h) = q_CL + q_CL,WT * WT  (intercept 0.244 + slope 0.0442 / kg)
    #   V/F  (L)   = q_V,WT * WT          (no intercept; slope 0.639 / kg)
    # CL/F is therefore a linear-additive form with the intercept exp(lcl)
    # at WT = 0 and the per-kg slope e_wt_cl; V/F is a multiplicative form
    # where exp(lvc) is V/F per kg of body weight and the per-subject volume
    # is exp(lvc) * WT.
    lcl <- log(0.244); label("Intercept of CL/F at WT = 0 (theta_CL, L/h)")            # Hussein 2005 Table 3 final-model q_CL
    lvc <- log(0.639); label("V/F per kg of body weight (theta_V/F per WT, L/kg)")     # Hussein 2005 Table 3 final-model q_V (no intercept)

    # Per-kg slope for the linear CL/F covariate effect (units: L/h/kg).
    e_wt_cl <- 0.0442; label("Slope of CL/F on body weight (theta_CL,WT, L/h/kg)")     # Hussein 2005 Table 3 final-model q_CL,WT

    # ka was modelled as ka = theta_ka + CL/V to guard against flip-flop;
    # theta_ka is the additive offset reported in Table 3 with %CV(IIV) =
    # 105.4. CL/V contributes on average less than 5% to ka (Hussein 2005
    # Table 3 footnote a).
    lka <- log(1.140); label("Additive offset on ka in ka = theta_ka + CL/V (theta_ka, 1/h)") # Hussein 2005 Table 3 final-model q_ka

    # Bioavailability fixed at 1.0 for the compliant population. Hussein 2005
    # additionally estimated F_noncomp = 0.281 (Table 3, q_F1) for 12% of
    # records flagged as suspected noncompliance; that effect is documented
    # in the vignette Assumptions and deviations but is not encoded here
    # because the noncompliance flag is a per-record retrospective
    # data-quality indicator rather than a prospective subject covariate.
    lfdepot <- fixed(log(1)); label("Bioavailability anchor for the compliant population (unitless)") # Hussein 2005 structural anchor (F1 = 1 for compliant)

    # Inter-patient variability (Hussein 2005 Table 3 footnote b): the
    # reported %CV approximates omega (square root of variance, x100) under a
    # first-order expansion of the exponential error model. Variances are
    # therefore (CV/100)^2 directly.
    etalcl ~ 0.05954    # Hussein 2005 Table 3 IIV(CL/F) = 24.4% CV -> 0.244^2
    etalvc ~ 0.02161    # Hussein 2005 Table 3 IIV(V/F)  = 14.7% CV -> 0.147^2
    etalka ~ 1.11092    # Hussein 2005 Table 3 IIV(ka)   = 105.4% CV -> 1.054^2

    # Residual variability is additive with two concentration-dependent
    # magnitudes (Hussein 2005 Table 3): 53.5 ng/mL for Cc <= 400 ng/mL and
    # 316 ng/mL for Cc > 400 ng/mL. The 400 ng/mL threshold was selected by
    # sensitivity analysis (Hussein 2005 Results, "BASE 1" derivation) and
    # has no clinical or therapeutic implication. The two SDs are switched
    # at simulation time in model() by comparing the predicted Cc against
    # the threshold; the per-bin parameters are declared via
    # paper_specific_residual_sds at the top of the function body.
    addSdLow  <- 53.5; label("Additive residual SD when predicted Cc <= 400 ng/mL (ng/mL)")  # Hussein 2005 Table 3 final-model additive RUV (Cc <= 400)
    addSdHigh <- 316;  label("Additive residual SD when predicted Cc > 400 ng/mL (ng/mL)")   # Hussein 2005 Table 3 final-model additive RUV (Cc > 400)
  })

  model({
    # Individual structural parameters. CL/F is linear-additive in body
    # weight with a non-zero intercept; V/F is linear-multiplicative without
    # an intercept (the per-kg slope is exp(lvc)). IIV is applied
    # multiplicatively on each per-subject parameter following the paper's
    # exponential error model.
    cl <- (exp(lcl) + e_wt_cl * WT) * exp(etalcl)
    vc <- exp(lvc + etalvc) * WT

    # ka = theta_ka + CL/V (Hussein 2005 Table 3 footnote a); theta_ka
    # carries the log-normal IIV with 105.4% CV. The CL/V additive
    # contribution uses the same per-subject cl and vc above so the
    # flip-flop guard tracks each subject's disposition.
    theta_ka <- exp(lka + etalka)
    ka       <- theta_ka + cl / vc

    # First-order elimination rate constant.
    kel <- cl / vc

    # One-compartment open model with first-order absorption from the depot.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Bioavailability for oral dosing (compliant population anchor F = 1).
    f(depot) <- exp(lfdepot)

    # Levocetirizine concentration in the central compartment. Dose in mg,
    # vc in L, so central / vc = mg/L = ug/mL; multiply by 1000 to express
    # in ng/mL (the source paper's concentration unit and the LLD scale).
    Cc <- central / vc * 1000

    # Concentration-stratified additive residual error. Predicted Cc <= 400
    # ng/mL uses addSdLow; predicted Cc > 400 ng/mL uses addSdHigh. The
    # arithmetic switch is equivalent to ifelse(Cc <= 400, addSdLow,
    # addSdHigh) but stays inside the rxode2 expression grammar.
    low_bin <- (Cc <= 400)
    addSd   <- addSdLow * low_bin + addSdHigh * (1 - low_bin)
    Cc ~ add(addSd)
  })
}
