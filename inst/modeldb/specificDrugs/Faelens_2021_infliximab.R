Faelens_2021_infliximab <- function() {
  description <- "One-compartment IV population PK model of infliximab in adults with moderate-to-severe ulcerative colitis (Faelens 2021 adapted model; baseline-covariate-only re-fit of Dreesen 2019)"
  reference <- "Faelens R, Wang Z, Bouillon T, Declerck P, Ferrante M, Vermeire S, Dreesen E. Model-Informed Precision Dosing during Infliximab Induction Therapy Reduces Variability in Exposure and Endoscopic Improvement between Patients with Ulcerative Colitis. Pharmaceutics. 2021;13(10):1623. doi:10.3390/pharmaceutics13101623"
  vignette <- "Faelens_2021_infliximab"
  units <- list(time = "day", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    SCORE_MAYO_E = list(
      description        = "Mayo endoscopic subscore at baseline (1, 2, or 3)",
      units              = "(score, 1-3)",
      type               = "categorical",
      reference_category = "2 (moderate active UC)",
      notes              = "Categorical effect on the elimination rate constant KE: typical KE values are 0.0422 /day for Mayo 1, 0.0463 /day for Mayo 2 (reference), and 0.0570 /day for Mayo 3 per Faelens 2021 supplement Table S1 (Adapted Model column). The original NONMEM dataset uses source column MPRE and additionally codes a sentinel `MPRE = -99` for missing values whose typical KE was carried over from the original Dreesen 2019 model with an unconverged initial estimate; that level is out-of-domain for this library implementation, which only supports Mayo 1/2/3.",
      source_name        = "MPRE"
    ),
    CONMED_STEROID = list(
      description        = "Baseline corticosteroid use (binary indicator)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no baseline corticosteroid use)",
      notes              = "Multiplicative fold-change on V (and therefore on CL since CL = KE * V): V is multiplied by 1.30 when CONMED_STEROID = 1 per Faelens 2021 supplement Table S1 (Adapted Model column; THETA(5) in the NONMEM control stream). Source column name `CS`; renamed to canonical CONMED_STEROID per inst/references/covariate-columns.md.",
      source_name        = "CS"
    ),
    DISEXT_EP = list(
      description        = "Extensive colitis at baseline (binary indicator)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not extensive colitis)",
      notes              = "Multiplicative fold-change on V (and therefore on CL): V is multiplied by 1.25 when DISEXT_EP = 1 per Faelens 2021 supplement Table S1 (Adapted Model column; THETA(8) in the NONMEM control stream). Source column name `EXTCOL`; renamed to canonical DISEXT_EP per inst/references/covariate-columns.md. The Faelens 2021 dataset uses a binary EXTCOL indicator (presence of extensive colitis) without a separate `other` category, so the paired DISEXT_OTHER canonical does not apply.",
      source_name        = "EXTCOL"
    ),
    FFM = list(
      description        = "Fat-free mass (Janmahasatian formula from WT, HT, SEX)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power scaling on V with reference 52 kg per Faelens 2021 supplement Table S1 (Adapted Model column; THETA(7) = 0.517 in the NONMEM control stream). FFM is computed in the NONMEM $PK block via the Janmahasatian (2005) formula: FFM_male = 9.27e3 * WT / (6.68e3 + 216 * BMI), FFM_female = 9.27e3 * WT / (8.78e3 + 244 * BMI), with BMI = WT / HT^2 (HT in metres). Reference 52 kg corresponds to the cohort-typical FFM in the original Dreesen 2019 dataset.",
      source_name        = "FFM"
    )
  )

  population <- list(
    n_subjects     = 204L,
    n_studies      = 1L,
    age_range      = "Adults; full demographics in source paper Dreesen 2019 (BJCP 85:782-795).",
    weight_range   = "Adults; full demographics in source paper Dreesen 2019.",
    sex_female_pct = NA_real_,
    race_ethnicity = "Not reported in Faelens 2021; single-center cohort (KU Leuven IBD Biobank, Belgium).",
    disease_state  = "Adults with moderate-to-severe ulcerative colitis (baseline Mayo endoscopic subscore typically 2 or 3).",
    dose_range     = "5 mg/kg or 10 mg/kg IV infliximab induction at days 0, 14, and 42 (label dosing or escalated dosing scenarios in the simulation paper).",
    regions        = "Single-center cohort, Leuven, Belgium (IBD Biobank B322201213950/S53684).",
    n_pk_samples   = 583L,
    notes          = "Model fit to 583 PK samples from 204 patients with UC (Faelens 2021 Section 2.1; original cohort in Dreesen 2019). The Faelens 2021 'Adapted Model' is a re-fit of the Dreesen 2019 popPK model that drops time-varying covariates (CRP, serum albumin) so that simulations of higher-than-observed doses are not biased by dose-dependent feedback through acute-phase proteins; the baseline covariates (Mayo endoscopic subscore, corticosteroid use, extensive colitis, FFM) are retained. Excluded interoccasion variability on KE (6.70% CV per occasion in the source) for nlmixr2lib portability; see vignette Assumptions section."
  )

  ini({
    # Structural parameters -- typical values for the reference patient
    # (baseline Mayo endoscopic subscore = 2, no corticosteroid use, no extensive
    # colitis, FFM = 52 kg) per Faelens 2021 supplement Table S1 (Adapted Model
    # column). The source paper parametrizes the one-compartment model with KE
    # and V; nlmixr2lib's standard parameterization is CL and Vc, so:
    #   typical CL = KE_Mayo2 * V_typ = 0.0463 /d * 6.97 L = 0.32271 L/d
    #   typical Vc = V_typ              = 6.97 L
    lcl <- log(0.32271); label("Clearance for the reference patient (CL, L/day)")            # Faelens 2021 supplement Table S1 (Adapted Model); KE_Mayo2 = 0.0463 /d, V_typ = 6.97 L; CL = KE * V
    lvc <- log(6.97);    label("Central volume of distribution for the reference patient (Vc, L)")  # Faelens 2021 supplement Table S1 (Adapted Model); THETA(6) = TVV = 6.97 L

    # Baseline Mayo endoscopic subscore effect on KE (and therefore on CL).
    # Source models the effect as a categorical lookup, with separate typical
    # KE for Mayo 1, 2, 3; encoded here as the log fold-change on CL relative
    # to the reference Mayo 2 category. KE values from Faelens 2021 supplement
    # Table S1 (Adapted Model column): KE_Mayo1 = 0.0422, KE_Mayo2 = 0.0463,
    # KE_Mayo3 = 0.0570 (all in /day):
    #   e_mayo1_cl = log(0.0422 / 0.0463) = -0.09275
    #   e_mayo3_cl = log(0.0570 / 0.0463) =  0.20785
    e_mayo1_cl <- -0.09275; label("Effect of baseline Mayo endoscopic subscore = 1 on CL (log fold-change vs Mayo 2; KE_Mayo1/KE_Mayo2 = 0.0422/0.0463)")  # Faelens 2021 supplement Table S1
    e_mayo3_cl <-  0.20785; label("Effect of baseline Mayo endoscopic subscore = 3 on CL (log fold-change vs Mayo 2; KE_Mayo3/KE_Mayo2 = 0.0570/0.0463)")  # Faelens 2021 supplement Table S1

    # Covariate effects on V (also propagate to CL since CL = KE * V).
    # Source covariate forms per Faelens 2021 supplement NONMEM control stream:
    #   TVV = THETA(6) * THETA(5)^CONMED_STEROID * (FFM/52)^THETA(7) * THETA(8)^DISEXT_EP
    # with THETA(5) = 1.30, THETA(7) = 0.517, THETA(8) = 1.25.
    e_conmed_steroid_vc   <- 1.30;  label("Multiplicative fold-change on Vc for baseline corticosteroid use (Vc * e_conmed_steroid_vc^CONMED_STEROID)")           # Faelens 2021 supplement Table S1; THETA(5)
    e_ffm_vc       <- 0.517; label("Power exponent of fat-free mass on Vc with reference FFM 52 kg ((FFM/52)^e_ffm_vc)")                       # Faelens 2021 supplement Table S1; THETA(7)
    e_disext_ep_vc <- 1.25;  label("Multiplicative fold-change on Vc for extensive colitis at baseline (Vc * e_disext_ep_vc^DISEXT_EP)")       # Faelens 2021 supplement Table S1; THETA(8)

    # Inter-individual variability. Source has independent ETAs on KE and V
    # (no $OMEGA BLOCK between ETA(1) and ETA(2) in the NONMEM control stream).
    # Reparameterised to (CL, Vc) using log(CL) = log(KE) + log(V):
    #   var(etalcl)         = var(ETA_KE) + var(ETA_V)        = 0.10577 + 0.05417 = 0.15994
    #   var(etalvc)         = var(ETA_V)                       = 0.05417
    #   cov(etalcl, etalvc) = var(ETA_V)                       = 0.05417
    # Source CV% values from supplement Table S1 (Adapted Model column):
    #   KE: 33.4% CV -> omega^2 = log(0.334^2 + 1) = 0.10577
    #   V : 23.6% CV -> omega^2 = log(0.236^2 + 1) = 0.05417
    # The induced correlation between log(CL) and log(Vc) is
    #   rho = cov / sqrt(var_cl * var_vc) = 0.05417 / sqrt(0.15994 * 0.05417) = 0.582.
    etalcl + etalvc ~ c(0.15994,
                        0.05417, 0.05417)   # Faelens 2021 supplement Table S1 IIV (33.4% CV on KE, 23.6% CV on V; reparameterized to CL, Vc)

    # Residual error (combined additive + proportional). The NONMEM $ERROR
    # block is Y = IPRED * (1 + ERR(1)) + ERR(2), which in nlmixr2 maps to
    # add(addSd) + prop(propSd). Final estimates from supplement Table S1
    # (Adapted Model column): proportional 32.9% CV, additive 0.300 mg/L FIX.
    addSd  <- 0.300;  label("Additive residual error (mg/L)")         # Faelens 2021 supplement Table S1 (Adapted Model); FIX
    propSd <- 0.329;  label("Proportional residual error (fraction)") # Faelens 2021 supplement Table S1 (Adapted Model); 32.9% CV
  })
  model({
    # Categorical effect of baseline Mayo endoscopic subscore on CL.
    # Source uses a 4-way categorical lookup on KE for MPRE in {1, 2, 3, missing};
    # this implementation supports MPRE in {1, 2, 3}, with Mayo 2 as reference.
    mayo_cl <- exp(e_mayo1_cl * (SCORE_MAYO_E == 1) + e_mayo3_cl * (SCORE_MAYO_E == 3))

    # Covariate effects on Vc (and propagated to CL since CL = KE * V).
    cov_vc <- (e_conmed_steroid_vc^CONMED_STEROID) * ((FFM / 52)^e_ffm_vc) * (e_disext_ep_vc^DISEXT_EP)

    # Individual PK parameters.
    cl <- exp(lcl + etalcl) * mayo_cl * cov_vc
    vc <- exp(lvc + etalvc) * cov_vc

    kel <- cl / vc

    d/dt(central) <- -kel * central

    # Concentration: dose in mg, volume in L -> mg/L
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
