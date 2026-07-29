Reijers_2016_trastuzumab <- function() {
  description <- "Three-compartment population PK model with parallel linear and Michaelis-Menten nonlinear elimination from the central compartment for intravenous trastuzumab in healthy male volunteers from a phase I biosimilarity trial of the FTMB biosimilar vs Herceptin reference product (Reijers 2016, combined model on all dose levels 0.49-6.44 mg/kg); covariates are lean body mass on central volume of distribution V1 and BMI on the linear elimination rate constant ke."
  reference <- "Reijers JAA, van Donge T, Schepers FML, Burggraaf J, Stevens J. Use of population approach non-linear mixed effects models in the evaluation of biosimilarity of monoclonal antibodies. Eur J Clin Pharmacol. 2016;72(11):1343-1352. doi:10.1007/s00228-016-2101-6"
  vignette <- "Reijers_2016_trastuzumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    LBM = list(
      description        = "Lean body mass at baseline",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline value. Linear (unit-slope) ratio effect on central volume V1 per Reijers 2016 Online Resource Equation 1: V1_i = V1_pop * (LBW_i / LBW_median) * exp(eta_i). Reference value 61 kg approximates the cohort-weighted mean of the combined-model dose groups in Table 1 (the paper does not state the analysis-dataset median directly). Source column LBW (lean body weight) is an alias of canonical LBM (lean body mass) and the values transfer without transformation.",
      source_name        = "LBW"
    ),
    BMI = list(
      description        = "Body mass index at baseline",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline value. Linear (unit-slope) ratio effect on the linear elimination rate constant ke per Reijers 2016 Online Resource Equation 2: ke_i = ke_pop * (BMI_i / BMI_median) * exp(eta_i). Reference value 23 kg/m^2 approximates the cohort-weighted mean of the combined-model dose groups in Table 1.",
      source_name        = "BMI"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 110L,
    n_studies       = 1L,
    n_observations  = 2159L,
    study           = "Phase I randomised, single-dose, parallel-group bioequivalence trial of FTMB (Synthon BV biosimilar candidate, test) vs Herceptin (EU-licenced reference) trastuzumab products, preceded by a placebo-controlled dose-escalation part (Wisman 2014).",
    age_range       = "18-45 years",
    age_summary     = "Mean (SD) 24-33 (2.3-9.1) years by dose group, per Reijers 2016 Table 1.",
    weight_summary  = "Mean (SD) total body weight 72.0-79.5 (7.5-12.6) kg by dose group, per Reijers 2016 Table 1.",
    weight_median   = "Combined-model cohort weighted mean total body weight is approximately 77 kg (computed from Table 1 dose-group means).",
    lbm_summary     = "Mean (SD) lean body mass 57.5-62.6 (3.8-8.4) kg by dose group, per Reijers 2016 Table 1.",
    lbm_median      = "Reference LBM = 61 kg (cohort-weighted mean across dose groups).",
    bmi_summary     = "Mean (SD) BMI 21.2-23.5 (2.1-3.3) kg/m^2 by dose group, per Reijers 2016 Table 1.",
    bmi_median      = "Reference BMI = 23 kg/m^2 (cohort-weighted mean across dose groups).",
    sex_female_pct  = 0,
    disease_state   = "Healthy male volunteers deemed healthy after full medical screening.",
    dose_range      = "Single 90-min IV infusion in 250 mL 0.9% NaCl. Dose-escalation part: 0.5 mg/kg (n=6), 1.5 mg/kg (n=6), 3 mg/kg (n=6) of test product (actual doses 0.49, 1.48, 2.96 mg/kg). Bioequivalence part: 6 mg/kg test (n=46, actual 5.96 mg/kg) and 6 mg/kg reference (n=46, actual 6.44 mg/kg).",
    regions         = "The Netherlands (Centre for Human Drug Research, Leiden).",
    reference_subject = "61 kg LBM, BMI 23 kg/m^2 (the reference subject implied by the cohort-weighted means used as covariate-equation references).",
    notes           = "Baseline demographics in Reijers 2016 Table 1 (n=110 across 5 treatment-arm strata). Trastuzumab serum concentrations quantified by ELISA with LLOQ 0.060 ug/mL; 1247 observations on the test product and 912 on the reference product (n=2159 total), of which 143 (test) and 51 (reference) were <LLOQ. All pre-dose <LLOQ values were set to zero before analysis; post-dose <LLOQ values were excluded. Drug product (test vs reference) was not a statistically significant covariate on any model parameter (max decrease in OFV 5.80, p>0.01), supporting biosimilarity. The combined model uses all 110 participants and all five dose strata; this file encodes that combined model. The paper additionally reports two separate models (test-only and reference-only at 6 mg/kg) for biosimilarity assessment; those models are not packaged here because the combined model is the published recommendation for general use and the separate models share the structural form."
  )

  ini({
    # Structural parameters (Reijers 2016 Table 2, combined-model final
    # estimates). The paper parameterises the model with:
    #   - linear elimination rate constant ke (1/h in the paper)
    #   - parallel Michaelis-Menten elimination with Vmax (ug/h) and
    #     Km (ug/L) on the central compartment
    #   - three-compartment distribution with central volume V1 (L) and
    #     two peripheral volumes V2, V3 (L), with intercompartmental
    #     clearances Q1, Q2 (L/h)
    # Conversion to nlmixr2lib canonical units (time = day, dose = mg,
    # concentration = ug/mL):
    #   Q1 (L/day)  = 2.91e-3 L/h * 24 h/day = 0.0698 L/day
    #   Q2 (L/day)  = 4.34e-2 L/h * 24 h/day = 1.0416 L/day
    #   Vmax (mg/day) = 178 ug/h * 24 h/day / 1000 = 4.272 mg/day
    #   Km (ug/mL)  = 937 ug/L / 1000 = 0.937 ug/mL (= 0.937 mg/L)
    #   ke (1/day)  = 2.20e-3 1/h * 24 h/day = 0.0528 1/day
    # Concentration in the central compartment is Cc = central / vc with
    # dose in mg and volumes in L, so Cc is in mg/L = ug/mL. Vmax is in
    # mg/day and Km is in mg/L so the MM elimination flux is in mg/day,
    # matching the linear flux kel * central.
    lvc   <- log(3.28);   label("Central volume V1 for the reference subject (L)")                    # Reijers 2016 Table 2, combined model
    lvp   <- log(1.89);   label("First peripheral volume V2 (L)")                                     # Reijers 2016 Table 2, combined model
    lvp2  <- log(1.96);   label("Second peripheral volume V3 (L)")                                    # Reijers 2016 Table 2, combined model
    lq    <- log(0.0698); label("Intercompartmental clearance Q1 to peripheral1 (L/day)")             # Reijers 2016 Table 2 (2.91e-3 L/h * 24)
    lq2   <- log(1.0416); label("Intercompartmental clearance Q2 to peripheral2 (L/day)")             # Reijers 2016 Table 2 (4.34e-2 L/h * 24)
    lvmax <- log(4.272);  label("Maximum nonlinear (Michaelis-Menten) elimination rate Vmax (mg/day)")# Reijers 2016 Table 2 (178 ug/h * 24/1000)
    lkm   <- log(0.937);  label("Michaelis-Menten constant Km (ug/mL = mg/L)")                        # Reijers 2016 Table 2 (937 ug/L / 1000)
    lkel  <- log(0.0528); label("Linear elimination rate constant ke for the reference subject (1/day)") # Reijers 2016 Table 2 (2.20e-3 1/h * 24)

    # Covariate effects (Reijers 2016 Online Resource Equation 1, Equation 2).
    # The supplement parameterises both effects as a unit-slope ratio:
    #   V1_i = V1_pop * (LBW_i / LBW_median) * exp(eta_i)
    #   ke_i = ke_pop * (BMI_i / BMI_median) * exp(eta_i)
    # i.e. the covariate exponent is structurally 1 (linear in the
    # covariate ratio) rather than an estimated power. Encoded with
    # fixed(1) to document the slope-1 form explicitly.
    e_lbm_vc  <- fixed(1); label("Power exponent of LBM on V1 (unitless, fixed at 1 by SI Eq. 1)")     # Reijers 2016 Online Resource Eq. 1
    e_bmi_kel <- fixed(1); label("Power exponent of BMI on ke (unitless, fixed at 1 by SI Eq. 2)")    # Reijers 2016 Online Resource Eq. 2

    # Inter-individual variability (Reijers 2016 Table 2 combined model).
    # Table 2 reports omega^2 directly; the CV% column is a derived
    # back-transformation (omega^2 = log(CV^2 + 1)). The reported omega^2
    # and CV% pairs are internally consistent:
    #   omega^2 V1 = 0.0217 (CV 14.8%)
    #   omega^2 KM = 0.121  (CV 35.9%)
    #   omega^2 ke = 0.0292 (CV 17.2%)
    # The paper notes "An omega block was required to correct for the
    # parameter correlation between KM and ke in the model" but neither
    # the main paper nor Online Resource reports the off-diagonal
    # covariance / correlation; the diagonal omega block here is the
    # closest faithful encoding given on-disk data. See the vignette's
    # Assumptions and deviations section.
    etalvc  ~ 0.0217   # V1 omega^2 -- Reijers 2016 Table 2 combined model
    etalkm  ~ 0.121    # KM omega^2 -- Reijers 2016 Table 2 combined model
    etalkel ~ 0.0292   # ke omega^2 -- Reijers 2016 Table 2 combined model

    # Residual error: combined proportional + additive (Reijers 2016
    # Results: "A combined residual error structure proved best fit for
    # purpose."). Table 2 reports the variances on the internal
    # ug/L scale:
    #   sigma^2 proportional = 0.0222  -> propSd = sqrt(0.0222) = 0.149 (unitless fraction; 14.9% CV consistent with Results)
    #   sigma^2 additive     = 1520    -> addSd  = sqrt(1520) ug/L = 38.99 ug/L = 0.039 ug/mL
    # The additive SD converted to canonical ug/mL units is 0.039 ug/mL
    # which is below the assay LLOQ (0.060 ug/mL) -- internally consistent.
    propSd <- 0.149;  label("Proportional residual error (fraction)")  # Reijers 2016 Table 2 (sigma^2 = 0.0222)
    addSd  <- 0.039;  label("Additive residual error (ug/mL)")         # Reijers 2016 Table 2 (sigma^2 = 1520 (ug/L)^2 = 1.52e-3 (ug/mL)^2)
  })

  model({
    # Individual parameters with covariate scaling per Reijers 2016
    # Online Resource Equations 1 and 2. Reijers 2016 reports a (KM, ke)
    # omega block but does not state the off-diagonal covariance in
    # either the main paper or the SI; etalkm and etalkel are therefore
    # encoded as independent log-normal random effects here (see the
    # vignette's Assumptions and deviations section).
    vc   <- exp(lvc + etalvc) * (LBM / 61)^e_lbm_vc
    vp   <- exp(lvp)
    vp2  <- exp(lvp2)
    q    <- exp(lq)
    q2   <- exp(lq2)
    vmax <- exp(lvmax)
    km   <- exp(lkm + etalkm)
    kel  <- exp(lkel + etalkel) * (BMI / 23)^e_bmi_kel

    # Micro-constants for the three-compartment model.
    k12 <- q  / vc
    k21 <- q  / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    Cc <- central / vc

    # Three-compartment IV model with parallel linear (rate constant kel)
    # and Michaelis-Menten elimination (Vmax, Km) from the central
    # compartment (Reijers 2016 Figure 1):
    #   linear elimination flux = kel * central (mg/day)
    #   MM elimination flux     = Vmax * Cc / (Km + Cc) (mg/day)
    # All trastuzumab doses enter the central compartment directly
    # (90-min IV infusion in the source data; encoded as bolus into
    # central in the validation vignette since the infusion is short
    # relative to the multi-week disposition half-life).
    d/dt(central)     <- -kel * central -
                          vmax * Cc / (km + Cc) -
                          k12 * central -
                          k13 * central +
                          k21 * peripheral1 +
                          k31 * peripheral2
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(peripheral2) <-  k13 * central - k31 * peripheral2

    Cc ~ add(addSd) + prop(propSd)
  })
}
