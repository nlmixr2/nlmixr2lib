Xu_2020_daratumumab <- function() {
  description <- "Two-compartment population PK model for intravenous daratumumab (anti-CD38 IgG1k) in adults with multiple myeloma, with parallel linear and Michaelis-Menten eliminations from the central compartment. The maximum velocity of the saturable (target-mediated) elimination decays mono-exponentially from its baseline value at first-order rate KDES, mimicking depletion of the CD38 target over weekly 16 mg/kg therapy (Xu 2020 MMY1001 D-Kd / D-KRd cohorts)."
  reference <- "Xu XS, Moreau P, Usmani SZ, et al. Split First Dose Administration of Intravenous Daratumumab for the Treatment of Multiple Myeloma (MM): Clinical and Population Pharmacokinetic Analyses. Adv Ther. 2020;37(4):1464-1478. doi:10.1007/s12325-020-01247-8"
  vignette <- "Xu_2020_daratumumab"
  units <- list(time = "hour", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight at baseline",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power covariate on CL (exponent 0.451) and on V1 (exponent 0.375). Reference 78.6 kg (Xu 2020 Online Resource 6 footnote).",
      source_name        = "WT"
    ),
    ALB = list(
      description        = "Baseline serum albumin",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power covariate on linear CL (exponent -1.149). Reference 37.0 g/L (Xu 2020 Online Resource 6 footnote). Units inferred from the reference value (a value of 37 corresponds to g/L; albumin reported in g/dL would be around 3.7).",
      source_name        = "ALB"
    ),
    MM_NIGG = list(
      description        = "Multiple-myeloma immunoglobulin type indicator: 1 = non-IgG MM, 0 = IgG MM",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 (non-IgG MM) -- in Xu 2020 the typical-value CL is anchored to non-IgG MM and an IgG-MM patient receives an additive shift of +0.806 (i.e., 80.6% higher CL); this is the opposite reference orientation from Fau 2020 isatuximab. The canonical MM_NIGG column semantics are preserved (1 = non-IgG, 0 = IgG); only the model's reference category differs.",
      notes              = "Additive (linear) shift on CL: TPMMCL = 1 + e_igg_cl * (1 - MM_NIGG), with e_igg_cl = 0.806 from Xu 2020 Online Resource 6. IgG MM patients have 80.6% higher linear clearance than non-IgG MM patients, hypothesised to reflect FcRn-mediated competition between endogenous IgG M-protein and the therapeutic IgG1k mAb.",
      source_name        = "Type of MM (IgG vs non-IgG)"
    ),
    SEXF = list(
      description        = "Biological sex indicator: 1 = female, 0 = male",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Additive (linear) shift on V1: SEXV1 = 1 + e_sexf_vc * SEXF, with e_sexf_vc = -0.205 from Xu 2020 Online Resource 6. Females have 20.5% lower central volume than males.",
      source_name        = "Sex"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 107L,
    n_studies      = 1L,
    age_range      = "34-85 years",
    age_median     = "65 years",
    weight_range   = "45.0-160.8 kg",
    weight_median  = "70.0 kg (D-Kd cohort); 79.9 kg (D-KRd cohort) per Xu 2020 Table 1",
    sex_female_pct = 45.8,
    race_ethnicity = c(White = 81.3, Black_or_African_American = 3.7, Asian = 2.8, American_Indian_or_Alaska_Native = 0.9, Not_reported = 11.2),
    disease_state  = "Multiple myeloma; D-Kd cohort (n=85) relapsed/refractory MM with 1-3 prior lines of therapy including bortezomib and an IMiD; D-KRd cohort (n=22) newly diagnosed MM. ECOG 0/1/2: ~42% / ~52% / ~8%.",
    dose_range     = "Daratumumab 16 mg/kg IV; first dose given as either a single 16-mg/kg infusion on Cycle 1 Day 1 (D-Kd n=10) or as two 8-mg/kg infusions on Cycle 1 Days 1 and 2 (D-Kd n=75, D-KRd n=22). Subsequent doses 16 mg/kg weekly for Cycles 1-2, every 2 weeks for Cycles 3-6, every 4 weeks thereafter.",
    regions        = "Multinational (MMY1001 phase 1b sites in France, Spain, and the United States; ClinicalTrials.gov NCT01998971).",
    notes          = "MMY1001 D-Kd and D-KRd cohorts; baseline demographics from Xu 2020 Table 1. Reference covariates for the typical-value equations (Online Resource 6 footnote): WT 78.6 kg, ALB 37.0 g/L, non-IgG MM, male. Model objective function value -1577.363; final-model condition number 27.80."
  )

  ini({
    # Structural typical-value PK parameters for the reference patient
    # (Xu 2020 Online Resource 6; reference: WT 78.6 kg, ALB 37.0 g/L,
    # non-IgG MM, male).
    lcl   <- log(0.00485); label("Linear clearance CL at reference covariates (L/hour)") # Xu 2020 Online Resource 6: CL 0.00485 L/hr, RSE 10%
    lvc   <- log(4.09);    label("Central volume of distribution V1 at reference (L)")   # Xu 2020 Online Resource 6: V1 4.09 L, RSE 3.5%
    lq    <- log(0.0642);  label("Intercompartmental clearance Q (L/hour)")              # Xu 2020 Online Resource 6: Q 0.0642 L/hr, RSE 8.0%
    lvp   <- log(3.06);    label("Peripheral volume of distribution V2 (L)")             # Xu 2020 Online Resource 6: V2 3.06 L, RSE 9.7%

    # Parallel Michaelis-Menten (target-mediated) elimination from the
    # central compartment. Vmax (mg/hour) decays mono-exponentially from
    # its baseline at first-order rate KDES (1/hour); KM is fixed.
    lvmax <- log(2.08);    label("Baseline Michaelis-Menten Vmax of nonlinear CL (mg/hour)") # Xu 2020 Online Resource 6: Vmax 2.08 mg/hr, RSE 16.7%
    lkdes <- log(0.0013);  label("First-order rate of decrease of Vmax (1/hour)")            # Xu 2020 Online Resource 6: KDES 0.0013 1/hr, RSE 17.8%
    lkm   <- fixed(log(0.93)); label("Michaelis-Menten Km of nonlinear CL (ug/mL)")          # Xu 2020 Online Resource 6: KM 0.93 ug/mL, fixed

    # Covariate effects. Power-form on continuous covariates; additive
    # (linear) shifts on binary covariates as parameterised in Xu 2020
    # (Online Resource 6 footnote): TVCL formula and TVV1 formula.
    e_wt_cl  <- 0.451;  label("Power exponent of WT/78.6 on linear CL (unitless)")              # Xu 2020 Online Resource 6: WT on CL, RSE 50.1%
    e_alb_cl <- -1.149; label("Power exponent of ALB/37.0 on linear CL (unitless)")             # Xu 2020 Online Resource 6: serum albumin on CL, RSE 27.2%
    e_igg_cl <- 0.806;  label("Additive shift on CL for IgG MM (vs non-IgG MM reference)")       # Xu 2020 Online Resource 6: Type of MM (IgG vs non-IgG) on CL, RSE 29.8%
    e_wt_vc  <- 0.375;  label("Power exponent of WT/78.6 on V1 (unitless)")                     # Xu 2020 Online Resource 6: WT on V1, RSE 38.1%
    e_sexf_vc <- -0.205; label("Additive shift on V1 for female sex (vs male reference)")       # Xu 2020 Online Resource 6: Sex on V1, RSE 22.2%

    # Inter-individual variability. Xu 2020 Online Resource 6 reports IIV
    # as % CV; convert to log-normal variance via omega^2 = log(CV^2 + 1).
    # CL    40.7%   -> log(0.407^2 + 1) = 0.1531
    # V1    21.8%   -> log(0.218^2 + 1) = 0.04643
    # Vmax  71.3%   -> log(0.713^2 + 1) = 0.4111
    # KDES  43.4%   -> log(0.434^2 + 1) = 0.1726
    # No IIV reported for Q, V2, or KM (KM is fixed).
    etalcl   ~ 0.1531  # Xu 2020 Online Resource 6: omega(CL)   40.7% CV, RSE 10.4%
    etalvc   ~ 0.04643 # Xu 2020 Online Resource 6: omega(V1)   21.8% CV, RSE 10.3%
    etalvmax ~ 0.4111  # Xu 2020 Online Resource 6: omega(Vmax) 71.3% CV, RSE 23.8%
    etalkdes ~ 0.1726  # Xu 2020 Online Resource 6: omega(KDES) 43.4% CV, RSE 45.6%

    # Residual error. Xu 2020 reports "Additive error term on the
    # log-scale (% CV) = 13.8%, RSE 1.3%". Additive error on the
    # log-transformed observation in NONMEM is equivalent to a
    # proportional error in linear space with SD ~ CV/100 for small
    # variances; encode as propSd = 0.138.
    propSd <- 0.138;   label("Proportional residual error (fraction)") # Xu 2020 Online Resource 6: additive error on log-scale 13.8% CV
  })
  model({
    # Individual PK parameters with covariate effects (Xu 2020 final
    # covariate model; Online Resource 6 footnote equations):
    #   TVCL = 0.00485 * (WT/78.6)^0.451 * (ALB/37.0)^-1.149 * TPMMCL
    #        TPMMCL = 1 for non-IgG MM (MM_NIGG = 1)
    #        TPMMCL = 1.806 for IgG MM  (MM_NIGG = 0)
    #   TVV1 = 4.09 * (WT/78.6)^0.375 * SEXV1
    #        SEXV1 = 1 for male   (SEXF = 0)
    #        SEXV1 = 0.795 for female (SEXF = 1)
    cl <- exp(lcl + etalcl) *
            (WT  / 78.6)^e_wt_cl *
            (ALB / 37.0)^e_alb_cl *
            (1 + e_igg_cl * (1 - MM_NIGG))
    vc <- exp(lvc + etalvc) *
            (WT / 78.6)^e_wt_vc *
            (1 + e_sexf_vc * SEXF)
    q    <- exp(lq)
    vp   <- exp(lvp)
    vmax <- exp(lvmax + etalvmax)
    kdes <- exp(lkdes + etalkdes)
    km   <- exp(lkm)

    # Time-varying Vmax: first-order decay from the baseline value at
    # rate KDES. KDES = 0.0013 /hr corresponds to a Vmax half-life of
    # ln(2) / 0.0013 ~ 533 hours (~22 days), capturing the gradual
    # depletion of CD38 target burden over weekly therapy.
    vmax_t <- vmax * exp(-kdes * t)

    # Two-compartment IV-input PK with parallel linear and (time-varying)
    # Michaelis-Menten eliminations from the central compartment. The
    # Michaelis-Menten term is in concentration form: Vmax is reported
    # in mg/h and is therefore a mass-rate that depends on Cc = central / vc.
    Cc <- central / vc

    d/dt(central)     <- -(cl / vc) * central -
                          vmax_t * Cc / (km + Cc) -
                          (q / vc) * central +
                          (q / vp) * peripheral1
    d/dt(peripheral1) <-  (q / vc) * central -
                          (q / vp) * peripheral1

    Cc ~ prop(propSd)
  })
}
