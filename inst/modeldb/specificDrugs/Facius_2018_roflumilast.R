Facius_2018_roflumilast <- function() {
  description <- paste(
    "Integrated population PK model for oral roflumilast and its primary",
    "active metabolite roflumilast N-oxide in adult patients with severe",
    "chronic obstructive pulmonary disease (COPD) (Facius 2018). The",
    "structural model is the joint parent-metabolite model previously",
    "developed by Lahu 2010 on 21 phase I + 2 phase II/III studies: a",
    "two-compartment parent disposition with first-order absorption and a",
    "shared lag time, and a two-compartment N-oxide disposition with",
    "first-order absorption from a separate pre-systemic dose compartment",
    "(relative bioavailability F5) plus complete first-order conversion",
    "from the parent central compartment. All structural disposition",
    "parameters and the F5 / KAm-to-KAp ratio are fixed to the Lahu 2010",
    "base-model estimates re-applied to OPTIMIZE via a Bayesian feedback",
    "MAXEVAL = 0 step; only the phase II-III dichotomous patient effects",
    "(on KA, parent CL, N-oxide CL, and N-oxide central volume), the",
    "covariate effects, the between-subject variability on parent and",
    "N-oxide clearance (with a Box-Cox-shape transformation), and the",
    "log-additive residual errors were estimated on the combined OPTIMIZE",
    "and REACT phase III dataset of 1238 + 461 patients. Covariates",
    "retained are body weight on all volume terms and on N-oxide CL,",
    "smoking (current vs not-current) on parent and N-oxide CL, age on",
    "parent and N-oxide CL, and sex on N-oxide CL. tPDE4i (total",
    "phosphodiesterase-4 inhibitory activity), the exposure metric used",
    "in the paper's downstream PK/adverse-event and PK/time-to-event",
    "models, is a per-dosing-interval summary derived from the predicted",
    "average plasma concentrations and is not embedded in the ODE",
    "system; see the vignette for the derivation.")
  reference <- paste(
    "Facius A, Marostica E, Gardiner P, Watz H, Lahu G.",
    "Pharmacokinetic and Pharmacodynamic Modelling to Characterize the",
    "Tolerability of Alternative Up-Titration Regimens of Roflumilast in",
    "Patients with Chronic Obstructive Pulmonary Disease.",
    "Clin Pharmacokinet. 2018;57(8):1029-1038.",
    "doi:10.1007/s40262-018-0671-4.",
    "Structural disposition parameters are fixed from Lahu 2010 (the",
    "Facius 2018 base model); see modellib('Lahu_2010_roflumilast') for an",
    "independent extraction of the upstream paper.")
  vignette <- "Facius_2018_roflumilast"
  units <- list(time = "h", dosing = "ug", concentration = "ug/L")

  covariateData <- list(
    WT = list(
      description        = "Total body weight at baseline (constant within an individual in the source dataset).",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric power scaling (WT/70)^1.22 on all volume terms (V2, V4, V3, V6) per Facius 2018 Online Resource Table S2 footnote e. Additional power scaling (WT/70)^0.273 on N-oxide CL. Reference weight is 70 kg (combined OPTIMIZE+REACT mean is 75.3 +/- 17.7 kg per Table 1).",
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Subject age at baseline (constant within an individual in the source dataset).",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling (AGE/60)^-0.611 on parent CL and (AGE/60)^-0.531 on N-oxide CL per Facius 2018 Online Resource Table S2 footnote g. Reference age is 60 years (combined OPTIMIZE+REACT mean is 64.4 +/- 8.2 years per Table 1; note this differs from the 40-year reference used by Lahu 2010 on the broader phase I-III dataset).",
      source_name        = "Age"
    ),
    SEXF = list(
      description        = "Biological sex indicator, 1 = female, 0 = male.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male).",
      notes              = "Facius 2018 Online Resource Table S2 footnote h codes the source covariate as SEX = 1 for male and 0 for female, with effect on N-oxide CL = -0.112 (i.e., males have 11.2% lower N-oxide CL than females). The canonical SEXF column inverts this; the model() block applies the source coefficient via the male indicator (1 - SEXF). Sex was not retained on parent CL or on any volume term.",
      source_name        = "Sex"
    ),
    SMOKE = list(
      description        = "Current-smoker indicator, 1 = current smoker at baseline, 0 = former/never smoker.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-current-smoker).",
      notes              = "Facius 2018 Online Resource Table S2 footnote f: SMK = 1 for current smokers and 0 for former/never smokers; linear additive effects +15.1% on parent CL and +15.1% on N-oxide CL. Mechanism is CYP1A2 induction in smokers (main text Discussion).",
      source_name        = "Smoking"
    ),
    DIS_COPD = list(
      description        = "Chronic obstructive pulmonary disease patient indicator, 1 = phase II/III COPD patient, 0 = phase I healthy volunteer.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (healthy volunteer; the implicit reference category from the Lahu 2010 base-model fit on phase I data).",
      notes              = "Dichotomous phase II-III patient effects re-estimated on the combined OPTIMIZE+REACT dataset (Facius 2018 Online Resource Table S2): -73.3% on KA (so KA in patients is 26.7% of the healthy-volunteer KA), -55.2% on parent CL, -24.4% on N-oxide CL, -20.7% on N-oxide V3. No phase II-III effect on parent volumes (V2, V4), N-oxide peripheral V6, parent or N-oxide Q, or on F5 / KAm-to-KAp ratio. The OPTIMIZE+REACT studies enrol only severe-COPD patients, so the model's most useful operating point is DIS_COPD = 1; DIS_COPD = 0 returns the phase I healthy-volunteer reference predictions.",
      source_name        = "Phase"
    )
  )

  covariatesDataExcluded <- list(
    RACE_ASIAN = list(
      description        = "Asian-vs-other race indicator, tested in the Facius 2018 formal covariate analysis but not retained.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Asian).",
      notes              = "Facius 2018 Online Resource Methods Section 2 page 3: 'Due to the low number of patients in the race categories different from Asian and White, the covariate race was tested as a binary relation (i.e. RASIA = 1 if race = Asian, RASIA = 0 if race different from Asian).' The covariate was screened on parent and N-oxide clearance (the two parameters with retained BSV) under the standard forward-inclusion / backward-deletion procedure but was not retained in the final model and is not present in Online Resource Table S2. Documented here for provenance only; the model does not use this column.",
      source_name        = "RASIA"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 1699L,
    n_studies      = 2L,
    age_range      = "40-92 years",
    age_median     = "64.0 years",
    weight_range   = "33.5-160 kg",
    weight_median  = "74.0 kg",
    sex_female_pct = 25.0,
    race_ethnicity = c(White = 92.9, Asian = 5.3, Black = 0.8, Other = 0.6, Hispanic = 0.4),
    disease_state  = "severe-to-very-severe chronic obstructive pulmonary disease (COPD) with chronic productive cough, current or former smokers, post-bronchodilator FEV1 <= 50% of predicted, FEV1/FVC < 70%. OPTIMIZE patients (n = 1238) received background standard-of-care LABA/ICS +/- LAMA +/- theophylline; REACT patients (n = 461) received background ICS/LABA +/- LAMA.",
    dose_range     = "Oral roflumilast 250 microgram (lower up-titration arm), 500 microgram every other day (alternate up-titration arm), or 500 microgram once daily (full-dose arm) for the initial 4-week up-titration phase in OPTIMIZE; 500 microgram once daily for 8 weeks of maintenance in OPTIMIZE and for 52 weeks in REACT. The OPTIMIZE down-titration sub-cohort received 250 microgram once daily after main-study discontinuation. The model file uses microgram (ug) as the dosing unit and ug/L as the concentration unit (numerically equivalent to ng/mL).",
    regions        = "Multinational; OPTIMIZE (NCT02165826) and REACT (NCT01329029) phase III studies.",
    n_observations = "OPTIMIZE: 18,983 quantifiable plasma samples from 1238 patients (5878 from 250 ug OD up-titration arm, 5835 from 500 ug EOD up-titration arm, 5556 from 500 ug OD arm, 1714 from 250 ug OD open-label down-titration). REACT: 3176 quantifiable plasma samples from 461 patients in the roflumilast 500 ug OD arm. Combined dataset: roflumilast 11,005 samples (9416 OPTIMIZE + 1589 REACT) and roflumilast N-oxide 11,140 samples (9553 OPTIMIZE + 1587 REACT) per Facius 2018 Table 1. Concentrations were quantified by validated HPLC with tandem mass-spectrometric detection.",
    notes          = "OPTIMIZE was a randomized, double-blind, three-arm, parallel-group phase III trial of three 4-week up-titration regimens (250 ug OD, 500 ug EOD, 500 ug OD) before escalation to 500 ug OD maintenance for 8 weeks. REACT was a randomized, double-blind phase III trial of roflumilast 500 ug OD vs placebo for 52 weeks. The integrated popPK model was built on the existing Lahu 2010 base model: structural disposition parameters (Tlag, KAp, CLp, V2, Qp, V4, CLm, V3, F5, KAm/KAp ratio, V6, Qm) were carried over unchanged via a Bayesian-feedback (MAXEVAL = 0) step, and only the phase II-III patient effects, BSV (with Box-Cox transformation on the parent and N-oxide clearance etas), residual error, and re-tested covariates (age, sex, weight, smoking, COPD status) were re-estimated on the combined OPTIMIZE+REACT dataset. The Box-Cox shape parameter (lambda = 0.704) is documented in the model file and discussed in the vignette Errata, but the model() block uses a plain log-normal eta to keep the file rxode2-compatible. Asian-vs-other race (RASIA) was tested in the formal covariate analysis but not retained; see covariatesDataExcluded."
  )

  ini({
    # ------------------------------------------------------------------
    # STRUCTURAL DISPOSITION PARAMETERS (fixed from Lahu 2010 base model;
    # re-applied to OPTIMIZE via Bayesian feedback MAXEVAL = 0 per
    # Facius 2018 Online Resource Methods Section 2 page 2). Values are
    # the "Fixed effects" block of Online Resource Table S2.
    # ------------------------------------------------------------------

    ltlag <- fixed(log(0.227))
    label("Absorption lag time, shared parent and N-oxide depot (h)")          # Facius 2018 Online Resource Table S2: Tlag = 0.227 h (Fixed)

    lka <- fixed(log(3.36))
    label("Parent first-order absorption rate constant at healthy reference (1/h)")  # Facius 2018 Online Resource Table S2: KAp = 3.36 1/h (Fixed)

    lcl <- fixed(log(12.6))
    label("Parent apparent clearance CLp/F at healthy reference (L/h)")        # Facius 2018 Online Resource Table S2: CLp = 12.6 L/h (Fixed)

    lvc <- fixed(log(63.9))
    label("Parent apparent central volume V2/F at 70 kg reference (L)")        # Facius 2018 Online Resource Table S2: V2 = 63.9 L (Fixed)

    lq <- fixed(log(15.6))
    label("Parent apparent inter-compartmental clearance Qp/F (L/h)")          # Facius 2018 Online Resource Table S2: Qp = 15.6 L/h (Fixed)

    lvp <- fixed(log(171))
    label("Parent apparent peripheral volume V4/F at 70 kg reference (L)")     # Facius 2018 Online Resource Table S2: V4 = 171 L (Fixed)

    lcl_noxide <- fixed(log(1.26))
    label("N-oxide apparent clearance CLm/F at healthy / 70 kg / 60 y / female reference (L/h)")  # Facius 2018 Online Resource Table S2: CLm = 1.26 L/h (Fixed)

    lvc_noxide <- fixed(log(13.9))
    label("N-oxide apparent central volume V3/F at healthy / 70 kg reference (L)")  # Facius 2018 Online Resource Table S2: V3 = 13.9 L (Fixed)

    lq_noxide <- fixed(log(4.83))
    label("N-oxide apparent inter-compartmental clearance Qm/F (L/h)")         # Facius 2018 Online Resource Table S2: Qm = 4.83 L/h (Fixed)

    lvp_noxide <- fixed(log(12.4))
    label("N-oxide apparent peripheral volume V6/F at 70 kg reference (L)")    # Facius 2018 Online Resource Table S2: V6 = 12.4 L (Fixed)

    lf5 <- fixed(log(0.0649))
    label("N-oxide pre-systemic depot relative bioavailability F5 (unitless; applied with MW correction MWm/MWp = 419.21/403.22 in model())")  # Facius 2018 Online Resource Table S2: F5 = 0.0649 (Fixed)

    lratio_ka <- fixed(log(0.635))
    label("Ratio of N-oxide depot KAm to parent depot KAp (unitless)")         # Facius 2018 Online Resource Table S2: RatioKAm/KAp = 0.635 (Fixed)

    # ------------------------------------------------------------------
    # PHASE II-III DICHOTOMOUS PATIENT EFFECTS (estimated on OPTIMIZE+REACT)
    # Applied multiplicatively as (1 + effect * DIS_COPD) so DIS_COPD = 0
    # returns the phase I healthy-volunteer reference and DIS_COPD = 1
    # returns the OPTIMIZE+REACT patient operating point.
    # ------------------------------------------------------------------

    e_pat_ka <- -0.733
    label("Phase II-III patient effect on KA (unitless; -73.3%)")              # Facius 2018 Online Resource Table S2 Phase II-III effects: Effect on KA = -73.3% (SE 1.73, 95% CI -76.7% to -69.9%)

    e_pat_cl <- -0.552
    label("Phase II-III patient effect on parent CL (unitless; -55.2%)")       # Facius 2018 Online Resource Table S2 Phase II-III effects: Effect on CLp = -55.2% (SE 1.08, 95% CI -57.3% to -53.1%)

    e_pat_cl_noxide <- -0.244
    label("Phase II-III patient effect on N-oxide CL (unitless; -24.4%)")      # Facius 2018 Online Resource Table S2 Phase II-III effects: Effect on CLm = -24.4% (SE 1.46, 95% CI -27.3% to -21.5%)

    e_pat_vc_noxide <- -0.207
    label("Phase II-III patient effect on N-oxide central V3 (unitless; -20.7%)")  # Facius 2018 Online Resource Table S2 Phase II-III effects: Effect on V3 = -20.7% (SE 8.47, 95% CI -37.3% to -4.1%)

    # Note: a Box-Cox shape parameter (lambda = 0.704; SE 0.0651, 95% CI
    # 0.576 to 0.832) was fitted on the parent and N-oxide CL etas in
    # the source (Facius 2018 Online Resource Table S2 Phase II-III
    # effects row "Box-Cox shape parameter"). The model file uses
    # unmodified log-normal etas (lambda = 0 equivalent) so the
    # parameter is not encoded here; see population$notes and the
    # vignette Errata for the deviation.

    # ------------------------------------------------------------------
    # COVARIATE EFFECTS (re-estimated on OPTIMIZE+REACT)
    # ------------------------------------------------------------------

    allov <- 1.22
    label("Allometric power exponent of (WT/70) on all volume terms (V2, V4, V3, V6) (unitless)")  # Facius 2018 Online Resource Table S2 Covariate effects: WT on all volume terms = 1.22 (SE 0.121, 95% CI 0.983 to 1.46)

    e_smk_cl <- 0.151
    label("Linear coefficient for current-smoker on parent CL (unitless; +15.1%)")  # Facius 2018 Online Resource Table S2 Covariate effects: Smok on CLp = 0.151 (SE 0.0372, 95% CI 0.0781 to 0.224)

    e_wt_cl_noxide <- 0.273
    label("Power exponent of (WT/70) on N-oxide CL (unitless)")                # Facius 2018 Online Resource Table S2 Covariate effects: WT on CLm = 0.273 (SE 0.0406, 95% CI 0.193 to 0.353)

    e_smk_cl_noxide <- 0.151
    label("Linear coefficient for current-smoker on N-oxide CL (unitless; +15.1%)")  # Facius 2018 Online Resource Table S2 Covariate effects: Smok on CLm = 0.151 (SE 0.0253, 95% CI 0.101 to 0.201)

    e_age_cl <- -0.611
    label("Power exponent of (AGE/60) on parent CL (unitless)")                # Facius 2018 Online Resource Table S2 Covariate effects: AGE on CLp = -0.611 (SE 0.138, 95% CI -0.881 to -0.341)

    e_age_cl_noxide <- -0.531
    label("Power exponent of (AGE/60) on N-oxide CL (unitless)")               # Facius 2018 Online Resource Table S2 Covariate effects: AGE on CLm = -0.531 (SE 0.0829, 95% CI -0.693 to -0.369)

    e_sex_m_cl_noxide <- -0.112
    label("Linear coefficient for male-vs-female on N-oxide CL (unitless; -11.2%; canonical SEXF inverts via (1 - SEXF))")  # Facius 2018 Online Resource Table S2 Covariate effects: Sex on CLm = -0.112 (SE 0.0186, 95% CI -0.148 to -0.0755)

    # ------------------------------------------------------------------
    # INTER-INDIVIDUAL VARIABILITY (re-estimated on OPTIMIZE+REACT)
    # ------------------------------------------------------------------
    # Facius 2018 Online Resource Table S2 reports omega^2 directly:
    # var(etalcl) = 0.227, var(etalcl_noxide) = 0.273, cov = 0.163
    # (translates to a 50.5% / 56.0% lognormal CV on the parent and
    # N-oxide CL respectively, with correlation rho =
    # 0.163 / sqrt(0.227 * 0.273) = 0.655). A Box-Cox shape parameter
    # (lambda = 0.704) was fitted on these etas in the source; the model
    # file uses unmodified log-normal etas (see vignette Errata).

    etalcl + etalcl_noxide ~ c(0.227,
                               0.163, 0.273)                                    # Facius 2018 Online Resource Table S2 Random effects: omega^2(CLp) = 0.227 (CV 50.5%), cov(CLp, CLm) = 0.163, omega^2(CLm) = 0.273 (CV 56.0%)

    # ------------------------------------------------------------------
    # RESIDUAL ERROR (re-estimated on OPTIMIZE+REACT; log-additive)
    # ------------------------------------------------------------------
    # The source residuals are reported as "SD in log scale" -- i.e.
    # epsilon ~ N(0, sigma^2) added to log(Cc), which maps to a
    # log-normal residual error in nlmixr2 (Cc ~ lnorm(expSd)).

    expSd <- 0.566
    label("Parent log-additive residual SD on log-concentration (unitless)")    # Facius 2018 Online Resource Table S2 Residual error: sigma_p,add (log-scale) = 0.566 (SE 0.0116, 95% CI 0.543 to 0.589)

    expSd_noxide <- 0.448
    label("N-oxide log-additive residual SD on log-concentration (unitless)")   # Facius 2018 Online Resource Table S2 Residual error: sigma_m,add (log-scale) = 0.448 (SE 0.0172, 95% CI 0.414 to 0.482)
  })

  model({
    # Reference values for the power / linear effect equations
    # (Facius 2018 Online Resource Table S2 footnotes e, g).
    ref_wt   <- 70
    ref_age  <- 60

    # Molecular weight ratio applied to the pre-systemic N-oxide depot
    # bioavailability so that a user-supplied roflumilast dose (in ug)
    # routed to depot_noxide enters as the equivalent N-oxide mass.
    # Facius 2018 Online Resource Methods Section 3 page 3:
    # "The dose was obtained by multiplying the roflumilast dose by the
    # ratio of the molecular weights (MW): MW metabolite / MW parent
    # (i.e. 419.21/403.22)."
    mw_ratio <- 419.21 / 403.22

    # ------------------------------------------------------------------
    # PARENT (roflumilast) individual parameters
    # ------------------------------------------------------------------

    tlag_dep <- exp(ltlag)
    ka       <- exp(lka) *
                (1 + e_pat_ka * DIS_COPD)
    cl       <- exp(lcl + etalcl) *
                (1 + e_pat_cl * DIS_COPD) *
                (1 + e_smk_cl * SMOKE) *
                (AGE / ref_age)^e_age_cl
    vc       <- exp(lvc) * (WT / ref_wt)^allov
    q        <- exp(lq)
    vp       <- exp(lvp) * (WT / ref_wt)^allov

    # ------------------------------------------------------------------
    # METABOLITE (roflumilast N-oxide) individual parameters
    # ------------------------------------------------------------------

    ka_noxide <- ka * exp(lratio_ka)
    cl_noxide <- exp(lcl_noxide + etalcl_noxide) *
                 (1 + e_pat_cl_noxide * DIS_COPD) *
                 (WT  / ref_wt)^e_wt_cl_noxide *
                 (1 + e_smk_cl_noxide * SMOKE) *
                 (AGE / ref_age)^e_age_cl_noxide *
                 (1 + e_sex_m_cl_noxide * (1 - SEXF))
    vc_noxide <- exp(lvc_noxide) *
                 (1 + e_pat_vc_noxide * DIS_COPD) *
                 (WT / ref_wt)^allov
    q_noxide  <- exp(lq_noxide)
    vp_noxide <- exp(lvp_noxide) * (WT / ref_wt)^allov

    # Micro-constants for the explicit ODE system.
    kel_p <- cl        / vc
    k12_p <- q         / vc
    k21_p <- q         / vp
    kel_m <- cl_noxide / vc_noxide
    k12_m <- q_noxide  / vc_noxide
    k21_m <- q_noxide  / vp_noxide

    # ODE system. Roflumilast follows a two-compartment disposition with
    # first-order absorption from the parent depot. Pre-systemic N-oxide
    # formation enters the N-oxide central compartment via a separate
    # first-order absorption from depot_noxide (rate constant KAm =
    # KAp * ratio); post-systemic N-oxide formation is the full mass
    # rate of parent elimination (kel_p * central) scaled by MWm/MWp,
    # since the source assumes complete conversion of roflumilast to
    # N-oxide (Online Resource Methods Section 1 page 2). N-oxide
    # disposition is two-compartment with first-order elimination.

    d/dt(depot)          <- -ka * depot
    d/dt(central)        <-  ka * depot -
                             kel_p * central -
                             k12_p * central + k21_p * peripheral1
    d/dt(peripheral1)    <-  k12_p * central - k21_p * peripheral1

    d/dt(depot_noxide)   <- -ka_noxide * depot_noxide
    d/dt(central_noxide) <-  ka_noxide * depot_noxide +
                             kel_p * central * mw_ratio -
                             kel_m * central_noxide -
                             k12_m * central_noxide + k21_m * peripheral2
    d/dt(peripheral2)    <-  k12_m * central_noxide - k21_m * peripheral2

    # Absorption modifiers. The parent depot has F1 = 1 (parent absolute
    # bioavailability fixed to 1, since no IV roflumilast data; Online
    # Resource Methods Section 1 page 2). The N-oxide pre-systemic depot
    # has F = F5 * MWm/MWp so that a user-supplied dose in roflumilast
    # micrograms is rescaled to the N-oxide-equivalent mass actually
    # absorbed pre-systemically.
    lag(depot)        <- tlag_dep
    lag(depot_noxide) <- tlag_dep
    f(depot)          <- 1.0
    f(depot_noxide)   <- exp(lf5) * mw_ratio

    # Observations. Dose is in ug, volumes in L, so central / vc gives
    # ug/L = ng/mL (the units used in Facius 2018 main-text Figure 2
    # axis labels).
    Cc        <- central        / vc
    Cc_noxide <- central_noxide / vc_noxide

    Cc        ~ lnorm(expSd)
    Cc_noxide ~ lnorm(expSd_noxide)
  })
}
