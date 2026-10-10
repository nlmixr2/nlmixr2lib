Li_2026_epcoritamab <- function() {
  description <- "Repeated time-to-event (RTTE) model for the hazard of Grade >= 2 cytokine release syndrome (CRS) after subcutaneous epcoritamab (CD3xCD20 bispecific antibody) in adults with relapsed or refractory aggressive or indolent B cell non-Hodgkin lymphoma (Li 2026). The hazard is the product of a stimulatory sigmoid Emax function of plasma epcoritamab concentration (SMAX, S50, SHILL) and a tolerance moderator: a turnover pool starting at 1 whose zero-order input is inhibited by concentration (Imax = 1, I50, Hill = 1) with kin = kout, so the hazard falls as exposure accumulates. Prior CAR T cell therapy lowers SMAX; Cycle 1 prophylaxis with intravenous fluids or dexamethasone, or with both, raises S50. No between-subject variability on the hazard. The concentration comes from the embedded two-compartment quasi-steady-state TMDD population PK model of Li 2025 (modellib('Li_2025_epcoritamab')), reproduced unchanged; body weight and age act only through that PK model. The model exposes the instantaneous hazard (per day), the cumulative hazard and the survival function sur = exp(-cumhaz); 1 - sur is the probability of at least one Grade >= 2 CRS event."
  reference <- paste(
    "Li T, Tredennick A, Polhamus D, Putnins M, Liu S, Sanghavi K, Thalhauser CJ,",
    "Parikh A, Noorani B, Mohamed MEF, Le Gallo C, Elliott B, Gupta M, Xu S.",
    "Epcoritamab Step-Up Dosing Regimen Selection and Optimization Using Repeated",
    "Time-to-Event Modeling for Cytokine Release Syndrome Risk Mitigation.",
    "Clin Pharmacol Ther. 2026;120(2):542-551. doi:10.1002/cpt.70362.",
    "Hazard equations in the Supplementary Methods; final estimates in Table 2.",
    "Embedded PK model: Li T, Gibiansky L, Parikh A, et al. Population",
    "Pharmacokinetics of Epcoritamab Following Subcutaneous Administration in",
    "Relapsed or Refractory B Cell Non-Hodgkin Lymphoma. Clin Pharmacokinet.",
    "2025;64(1):127-141. doi:10.1007/s40262-024-01464-2 (reference 15 of the",
    "CRS paper; modellib('Li_2025_epcoritamab')).",
    sep = " "
  )
  vignette <- "Li_2026_epcoritamab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  # etaRUV is the random effect on the magnitude of the PK residual error,
  # carried over from Li 2025 (supplementary Fig. S1, eta9); it has no
  # canonical IIV name.
  paper_specific_etas <- c("etaRUV")

  covariateData <- list(
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Acts on the embedded Li 2025 PK model only (power covariate on CL/F, Q/F, Vc/F and Vp/F, normalised to 75 kg), and so on the CRS hazard only through concentration. The CRS paper simulated virtual patients with body weight drawn from the observed pooled-population distribution (Methods, 'Model-based simulation'); Table S1 reports a median of 72.80 kg (range 45.99-109.41) in the 600 patients.",
      source_name = "WT"
    ),
    AGE = list(
      description = "Baseline age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Acts on the embedded Li 2025 PK model only (power covariate on ka normalised to 65 years). Table S1 median 66 years (range 33-83).",
      source_name = "AGE"
    ),
    PRIOR_CART = list(
      description = "Prior chimeric antigen receptor (CAR) T cell therapy, 1 = yes",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no prior CAR T cell therapy)",
      notes = "Time-fixed per subject. Multiplies SMAX by exp(-0.964) = 0.382 (a 61.8% lower maximum stimulation of the CRS hazard; Table 2 theta9). 125 of 600 patients (20.8%) had prior CAR T cell therapy (Table 1). The paper's Table 2 also lists effects of prior CAR T cell therapy on Kin (theta8) and on S50 (theta10), both fixed to 0 and therefore not part of the final model.",
      source_name = "prior CAR T cell therapy (yes/no)"
    ),
    CONMED_DEXAMETHASONE = list(
      description = "Dexamethasone given as the corticosteroid for CRS prophylaxis with the first four epcoritamab doses in Cycle 1, 1 = yes",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no dexamethasone prophylaxis; includes patients premedicated with prednisolone)",
      notes = "Time-fixed per subject: a patient-level category describing Cycle 1 prophylaxis, applied for the whole simulated period, as in the source. Patients who received prednisolone and no dexamethasone are 0 (Methods, 'Study design and input data'). Dexamethasone or IV fluids given to MANAGE a CRS event were not used for the categorisation. Combined with CONMED_IV_FLUIDS into the paper's three-level prophylaxis covariate: neither (reference), exactly one (S50 x exp(1.19) = 3.28, Table 2 theta12) or both (S50 x exp(1.63) = 5.10, theta13).",
      source_name = "CRS prophylaxis (dexamethasone/IV fluids) category, Table 1"
    ),
    CONMED_IV_FLUIDS = list(
      description = "Intravenous fluid hydration given as CRS prophylaxis with the first four epcoritamab doses in Cycle 1, 1 = yes",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no prophylactic IV fluids)",
      notes = "Time-fixed per subject. See CONMED_DEXAMETHASONE for how the two indicators combine into the paper's three-level Cycle 1 prophylaxis covariate on S50. IV fluids were evaluated only in the EPCORE NHL-1 dose-optimisation part (Table 1).",
      source_name = "CRS prophylaxis (dexamethasone/IV fluids) category, Table 1"
    )
  )

  covariatesDataExcluded <- list(
    TUMTP_INHL = list(
      description = "Indolent (versus aggressive) non-Hodgkin lymphoma indicator, 1 = iNHL",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (aggressive NHL)",
      notes = "Disease type (aNHL/iNHL) was considered as a potential modifier of the CRS hazard (Methods) but was not retained in the final model; Table 2 carries no disease-type parameter. The paper simulates aNHL and iNHL separately (Figures 3-5) with the same hazard parameters. 331 aNHL (55.2%), 269 iNHL (44.8%) (Table 1)."
    ),
    STUDY_NHL3 = list(
      description = "EPCORE NHL-3 (Japan) study indicator, 1 = NHL-3",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (EPCORE NHL-1)",
      notes = "Effect of study = NHL-3 on S50 (Table 2 theta11) was fixed to 0 and not part of the final model; the study effect 'did not demonstrate statistical significance' (Results)."
    )
  )

  compartmentData <- list(
    depot = list(analyte = "epcoritamab", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(
      analyte = "epcoritamab (total: free + target-bound)",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(analyte = "epcoritamab", units = "mg", specimen = "plasma", verified = TRUE),
    moderator1 = list(
      analyte = "tolerance moderator of the CRS hazard (EFF_Inh, unitless fraction)",
      units = NA_character_,
      specimen = "not applicable",
      verified = TRUE
    ),
    cumhaz = list(
      analyte = "cumulative hazard of Grade >= 2 CRS",
      units = NA_character_,
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 600L,
    n_studies = 2L,
    age_range = "33-83 years",
    age_median = "66 years",
    weight_range = "45.99-109.41 kg",
    weight_median = "72.80 kg",
    disease_state = "Relapsed or refractory B cell non-Hodgkin lymphoma: aggressive NHL including (D)LBCL 331 (55.2%), indolent NHL including FL 269 (44.8%); prior CAR T cell therapy 125 (20.8%)",
    dose_range = "Subcutaneous epcoritamab in 28-day cycles; priming, intermediate and full doses 0.0128-60 mg in dose escalation; 2-step-up 0.16/0.8/48 mg (expansion) and 3-step-up 0.16/0.8/3 or 6/48 mg (FL optimisation); full doses QW in cycles 1-3, Q2W in cycles 4-9 and Q4W from cycle 10",
    co_medication = "Cycle 1 CRS prophylaxis: no dexamethasone and no IV fluids 373 (62.2%), IV fluids only 30 (5.0%), dexamethasone only 91 (15.2%), dexamethasone and IV fluids 106 (17.7%); prednisolone 100 mg was the default corticosteroid outside the optimisation part",
    regions = "EPCORE NHL-1 (US/EU) and EPCORE NHL-3 (Japan)",
    notes = "Pooled EPCORE NHL-1 (NCT03625037; dose escalation/expansion n = 364, dose optimisation n = 172) and EPCORE NHL-3 (NCT04542824; n = 64) monotherapy patients (Table 1). Age and weight from supplementary Table S1 (median [range]). The embedded PK model was developed in 327 of these patients (Li 2025). Endpoint: Grade >= 2 CRS events (repeated)."
  )

  ini({
    # ---- Embedded PK model: Li 2025 (Clin Pharmacokinet 64:127), Table 2 ----
    # Reproduced unchanged from modellib('Li_2025_epcoritamab'). The CRS paper
    # used this PPK model with individual EBEs for calibration and typical PK
    # values with sampled covariates for simulation (Methods).
    lcl <- log(0.481); label("Apparent nonspecific (linear) clearance CL/F (L/day)") # Li 2025 Table 2: CL/F = 0.481 L/day
    lq <- log(0.488); label("Apparent intercompartmental clearance Q/F (L/day)") # Li 2025 Table 2: Q/F = 0.488 L/day
    lvc <- log(9.33); label("Apparent central volume of distribution Vc/F (L)") # Li 2025 Table 2: Vc/F = 9.33 L
    lvp <- log(14.1); label("Apparent peripheral volume of distribution Vp/F (L)") # Li 2025 Table 2: Vp/F = 14.1 L
    lka <- log(0.584); label("First-order subcutaneous absorption rate constant ka (1/day)") # Li 2025 Table 2: ka = 0.584 1/day

    lrbase_target <- log(2.03); label("Baseline (constant) total target concentration BASE (ug/mL)") # Li 2025 Table 2: BASE = 2.03 ug/mL
    lkss <- log(0.214); label("Quasi-steady-state constant KSS (ug/mL)") # Li 2025 Table 2: KSS = 0.214 ug/mL
    lkint <- log(0.0278); label("Internalization rate constant of the drug-target complex kint (1/day)") # Li 2025 Table 2: kint = 0.0278 1/day

    e_wt_cl <- 0.875; label("Power exponent of (WT/75) on CL/F (unitless)") # Li 2025 Table 2: weight on CL/F = 0.875
    e_wt_q <- fixed(0.75); label("Power exponent of (WT/75) on Q/F (unitless)") # Li 2025 Table 2: weight on Q/F = 0.75 (fixed)
    e_wt_vc <- 0.603; label("Power exponent of (WT/75) on Vc/F (unitless)") # Li 2025 Table 2: weight on Vc/F = 0.603
    e_wt_vp <- fixed(1); label("Power exponent of (WT/75) on Vp/F (unitless)") # Li 2025 Table 2: weight on Vp/F = 1 (fixed)
    e_age_ka <- -0.503; label("Power exponent of (AGE/65) on ka (unitless)") # Li 2025 Table 2: age on ka = -0.503

    # Li 2025 Table 2 IIV (CV%), omega^2 = log(1 + CV^2)
    etalcl ~ 0.063959 # Li 2025 Table 2: IIV CL/F CV 25.7%
    etalq ~ 0.56355 # Li 2025 Table 2: IIV Q/F CV 87.0%
    etalvc ~ 0.092893 # Li 2025 Table 2: IIV Vc/F CV 31.2%
    etalvp ~ 1.0615 # Li 2025 Table 2: IIV Vp/F CV 137.5%
    etalka ~ 0.26176 # Li 2025 Table 2: IIV ka CV 54.7%
    etalrbase_target ~ 0.32885 # Li 2025 Table 2: IIV BASE CV 62.4%
    etalkss ~ 0.54082 # Li 2025 Table 2: IIV KSS CV 84.7%
    etalkint ~ 0.45109 # Li 2025 Table 2: IIV kint CV 75.5%
    etaRUV ~ 0.043963 # Li 2025 Table 2: random effect on residual-error magnitude CV 21.2%

    propSd <- 0.189; label("Proportional residual error of the PK model (fraction)") # Li 2025 Table 2: proportional residual error = 0.189
    addSd <- 0.0133; label("Additive residual error of the PK model (ug/mL)") # Li 2025 Table 2: additive residual error = 0.0133 ug/mL

    # ---- Grade >= 2 CRS RTTE hazard: Li 2026 Table 2 ----
    # Table 2 'Estimate' column is on the log scale; 'Transformed estimate' is
    # exp(Estimate). Parameters printed as 0 (-) (theta8, theta10, theta11)
    # were not used in the final model and are omitted. All four IIV variances
    # (omega1,1-omega4,4) are 0, so the hazard has no random effects.
    lkin_moderator1 <- -0.832; label("Turnover rate constant of the tolerance moderator, Kin = Kout (1/day)") # Table 2 theta1 = -0.832 -> Kin = 0.435 1/day
    lemax <- 0.0539; label("Maximum stimulated Grade >= 2 CRS hazard SMAX (1/day)") # Table 2 theta2 = 0.0539 -> SMAX = 1.06
    lec50 <- -2.36; label("Epcoritamab concentration giving half-maximal hazard stimulation S50 (ug/mL)") # Table 2 theta3 = -2.36 -> S50 = 0.0943 mg/L
    lhill <- -0.124; label("Hill coefficient of the hazard stimulation SHILL (unitless)") # Table 2 theta7 = -0.124 -> SHILL = 0.884
    lic50 <- -7.72; label("Epcoritamab concentration giving half-maximal inhibition of moderator input I50 (ug/mL)") # Table 2 theta4 = -7.72 -> I50 = 0.000446 mg/L
    limax <- fixed(log(1)); label("Maximum inhibition of the moderator input IMAX (fraction)") # Table 2 theta5 = 1.0 (-), used at the fixed value 1
    lhill_inh <- fixed(log(1)); label("Hill coefficient of the moderator-input inhibition IHILL (unitless)") # Table 2 theta6 = 1.0 (-), used at the fixed value 1

    e_prior_cart_emax <- -0.964; label("Log-scale effect of prior CAR T cell therapy on SMAX (unitless)") # Table 2 theta9 = -0.964 -> 0.382-fold
    e_proph_either_ec50 <- 1.19; label("Log-scale effect of Cycle 1 IV fluids OR dexamethasone (exactly one) on S50 (unitless)") # Table 2 theta12 = 1.19 -> 3.28-fold
    e_proph_both_ec50 <- 1.63; label("Log-scale effect of Cycle 1 IV fluids AND dexamethasone on S50 (unitless)") # Table 2 theta13 = 1.63 -> 5.10-fold
  })
  model({
    # ---- Embedded Li 2025 PK model (supplementary Fig. S1 of Li 2025) ----
    cl <- exp(lcl + etalcl) * (WT / 75)^e_wt_cl
    q <- exp(lq + etalq) * (WT / 75)^e_wt_q
    vc <- exp(lvc + etalvc) * (WT / 75)^e_wt_vc
    vp <- exp(lvp + etalvp) * (WT / 75)^e_wt_vp
    ka <- exp(lka + etalka) * (AGE / 65)^e_age_ka
    rbase_target <- exp(lrbase_target + etalrbase_target)
    kss <- exp(lkss + etalkss)
    kint <- exp(lkint + etalkint)

    # QSS free concentration, rationalised so it does not cancel to 0 at low
    # total concentration (identical algebra to Li_2025_epcoritamab).
    ctot <- central / vc
    qss_a <- ctot - rbase_target - kss
    cfree <- 2 * kss * ctot / (sqrt(qss_a^2 + 4 * kss * ctot) - qss_a)

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - cl * cfree - q * cfree + q / vp * peripheral1 -
      kint * rbase_target * cfree * vc / (kss + cfree)
    d/dt(peripheral1) <- q * cfree - q / vp * peripheral1

    # ---- Grade >= 2 CRS hazard (main-text Model development; Supplementary
    # Methods 'STIM and EFF_Inh models'; Figure 1) ----
    # Cp is the PK model's predicted plasma (free) epcoritamab concentration,
    # mg/L = ug/mL, the quantity the PK model was fitted to.
    cp_crs <- max(cfree, 0)

    # Cycle 1 prophylaxis: the paper's three-level covariate (neither / IV
    # fluids or dexamethasone / both) rebuilt from two indicators.
    proph_both <- CONMED_DEXAMETHASONE * CONMED_IV_FLUIDS
    proph_either <- CONMED_DEXAMETHASONE + CONMED_IV_FLUIDS - 2 * proph_both

    # Covariates enter as an exponentiated linear combination (Methods).
    kin_moderator1 <- exp(lkin_moderator1)
    emax <- exp(lemax + e_prior_cart_emax * PRIOR_CART)
    ec50 <- exp(lec50 + e_proph_either_ec50 * proph_either + e_proph_both_ec50 * proph_both)
    hill <- exp(lhill)
    ic50 <- exp(lic50)
    imax <- exp(limax)
    hill_inh <- exp(lhill_inh)

    # STIM(t) = SMAX * Cp^SHILL / (S50^SHILL + Cp^SHILL)
    stim <- emax * cp_crs^hill / (ec50^hill + cp_crs^hill)

    # dEFF_Inh/dt = Kin * (1 - IMAX * Cp^IHILL / (I50^IHILL + Cp^IHILL)) -
    #               Kout * EFF_Inh, with Kin = Kout so EFF_Inh -> 1 without drug.
    inh <- imax * cp_crs^hill_inh / (ic50^hill_inh + cp_crs^hill_inh)
    d/dt(moderator1) <- kin_moderator1 * (1 - inh) - kin_moderator1 * moderator1
    moderator1(0) <- 1

    # h(t) = STIM(t) x EFF_Inh(t)  (events per day)
    hazard <- stim * moderator1
    d/dt(cumhaz) <- hazard
    cumhaz(0) <- 0
    # Probability of no Grade >= 2 CRS event since time 0; 1 - sur is the
    # probability of at least one event (the quantity in Figures 3-5).
    sur <- exp(-cumhaz)

    # ---- PK observation (Li 2025) ----
    Cc <- cfree
    propSdi <- propSd * exp(etaRUV)
    addSdi <- addSd * exp(etaRUV)
    Cc ~ add(addSdi) + prop(propSdi) + combined2()
  })
}
