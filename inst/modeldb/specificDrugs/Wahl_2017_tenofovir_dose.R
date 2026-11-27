Wahl_2017_tenofovir_dose <- function() {
  description <- paste0(
    "Preclinical (mouse). Static landmark dose-response model for the ",
    "fraction of BLT-humanized mice protected from vaginal HIV-1 ",
    "acquisition by systemic tenofovir disoproxil fumarate (TDF) ",
    "pre-exposure prophylaxis (Wahl 2017). This is the dose-response arm ",
    "of a five-model PK-PD comparison (dose, plasma TFV, female ",
    "reproductive tract [FRT] TFV, FRT TFV-diphosphate, cervicovaginal ",
    "lavage TFV); the plasma and dose models gave the most precise fits ",
    "and plasma TFV was judged the best predictor of efficacy. The ",
    "protected fraction is a fractional sigmoidal Hill (Emax) function of ",
    "the administered TDF dose in mg/kg, E0 + (Emax - E0) * ",
    "DOSE^H / (EC50^H + DOSE^H) with Emax fixed to 1 (Wahl 2017 Eqn 2, ",
    "Table 2 'Dose' column). There is no PK layer and no ODE: the TDF ",
    "dose is supplied as the DOSE_TDF_MGKG covariate column, and the ",
    "probability of protection is emitted as the landmark output ",
    "prob_protection. The models were fitted by maximum likelihood in ",
    "ADAPT5, with the PD (protection) data from BLT mice and the exposure ",
    "data from parallel BALB/c mice given the same TDF doses."
  )
  reference <- paste(
    "Wahl A, Ho PT, Denton PW, Garrett KL, Hudgens MG, Swartz G, O'Neill C,",
    "Veronese F, Kashuba AD, Garcia JV.",
    "Predicting HIV Pre-exposure Prophylaxis Efficacy for Women using a",
    "Preclinical Pharmacokinetic-Pharmacodynamic In Vivo Model.",
    "Sci Rep. 2017 Jan 31;7:41098. doi:10.1038/srep41098. PMCID PMC5286499.",
    sep = " "
  )
  vignette <- "Wahl_2017_tenofovir_prep"

  units <- list(
    time = "n/a (static landmark dose-response model; no time dimension)",
    dosing = "n/a (no dose events; the TDF dose enters as the DOSE_TDF_MGKG covariate column)",
    concentration = "DOSE_TDF_MGKG in mg/kg; prob_protection is a probability (0-1)"
  )

  covariateData <- list(
    DOSE_TDF_MGKG = list(
      description = "Daily systemic (intraperitoneal) tenofovir disoproxil fumarate dose that each mouse received, in mg of TDF per kg body weight.",
      units = "mg/kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "The administered TDF dose is the exposure driver of this ",
        "dose-response arm. Wahl 2017 used four dose levels (20, 50, 140, ",
        "300 mg/kg) given once daily; BLT mice were dosed for seven days ",
        "with vaginal HIV challenge 3 h after the third dose (Fig. 1a, ",
        "Methods 'Experimental design'). A value of 0 is an untreated ",
        "control mouse, giving prob_protection = E0 (the baseline ",
        "protected fraction). Member of the DOSE_<DRUG>_<UNITS> ",
        "auto-approve covariate family; 'TDF' is tenofovir disoproxil ",
        "fumarate (the prodrug actually dosed), distinct from the ",
        "registered DOSE_TAF_MG (tenofovir alafenamide)."
      ),
      source_name = "TDF Dose (mg/kg)"
    )
  )

  population <- list(
    species = "mouse (BLT humanized, protection endpoint; exposures from parallel BALB/c mice)",
    n_subjects = "78 BLT mice across the four TDF arms contributed the protection outcomes (20 mg/kg n=14, 50 mg/kg n=12, 140 mg/kg n=13, 300 mg/kg n=11, plus 28 untreated controls); exposures were measured in parallel BALB/c mice (n=8 per dose)",
    n_studies = 1L,
    n_observations = "Binary protection outcome per BLT mouse (protected = absence of HIV-RNA in plasma at every time point and of HIV-DNA in blood and tissues at necropsy)",
    age_range = "Not reported (young adult laboratory mice)",
    sex_female_pct = 100,
    race_ethnicity = "n/a (mouse)",
    disease_state = "HIV-1 pre-exposure prophylaxis challenge model. BLT (bone marrow/liver/thymus) humanized NSG mice reconstituted with a human immune system, susceptible to vaginal HIV-1(JR-CSF) infection; BALB/c mice (immune-deficient parental strain, not HIV-susceptible) used as the exposure surrogate.",
    dose_range = "TDF 20, 50, 140, or 300 mg/kg intraperitoneally once daily; TDF solubilized in saline at 10 mg/mL",
    regions = "USA (University of North Carolina at Chapel hill)",
    notes = paste0(
      "Observed protection by dose (Table 1, uninfected / total): 20 ",
      "mg/kg 7/14 (50%), 50 mg/kg 8/12 (67%), 140 mg/kg 11/13 (85%), 300 ",
      "mg/kg 11/11 (100%), untreated controls 7/28 (25%). The 50 mg/kg ",
      "dose gave a 70% risk reduction comparable to the ~70% seen in ",
      "highly adherent women in Partners PrEP, and produced plasma TFV ",
      "(~370 ng/mL) similar to human exposures. HIV-1(JR-CSF) challenge: ",
      "3.5e5 TCIU in 20 uL intravaginally."
    )
  )

  ini({
    # ================================================================
    # Wahl 2017 Eqn 2 (fractional sigmoidal Hill / Emax model) fitted to
    # the BLT-mouse protection data against administered TDF dose, final
    # ADAPT5 estimates from Table 2, 'Dose (mg/kg)' column. Values are
    # mean (relative standard error [%SE]).
    #
    # E0 is the baseline protected fraction at zero dose; it is a
    # probability in (0, 1) and is kept on the natural scale (as in
    # Chen_2023_tilmicosin.R, where E0 is likewise an unconstrained-sign
    # baseline effect). EC50 and the Hill coefficient are kept positive
    # via the canonical log transform. Emax was FIXED to 1 for every one
    # of the paper's models ('Emax was fixed to one for all models',
    # Methods 'PD model development').
    # ================================================================
    e0 <- 0.2502
    label("Baseline protected fraction at zero TDF dose E0 (probability)")  # Wahl 2017 Table 2, Dose column E0 = 0.2502 (%SE 3.67)
    lemax <- fixed(log(1))
    label("Log maximum protected fraction Emax (unitless); Wahl 2017 held Emax at 1 for all models")  # Wahl 2017 Methods 'PD model development': 'Emax was fixed to one for all models'
    lec50 <- log(38.7)
    label("Log TDF dose giving half-maximal protection EC50 (mg/kg)")  # Wahl 2017 Table 2, Dose column EC50 = 38.7 mg/kg (%SE 7.628)
    lhill <- log(1.138)
    label("Log Hill coefficient H of the dose-response curve (unitless)")  # Wahl 2017 Table 2, Dose column H = 1.138 (%SE 10.4); estimated for the dose and plasma models
    # ================================================================
    # Wahl 2017 fitted 'proportional residual error variance models ...
    # with CV% of the assays as initial estimates' (Methods 'PD model
    # development'), but Table 2 reports no final residual-error estimate,
    # and the packaged model emits a deterministic typical-value
    # probability. A tiny placeholder additive residual is exposed so the
    # nlmixr2 observation machinery accepts the model, following the
    # prob_<endpoint> landmark convention (Komatsu_2016_penicillinG_*,
    # Fukae_2024_valemetostat_* families). Not from the source.
    # ================================================================
    addSd_prob_protection <- fixed(0.001)
    label("Placeholder additive residual SD on the typical-value protected fraction; Wahl 2017 used a proportional error whose final estimate it did not report")  # not from source; see vignette Assumptions and deviations
  })

  model({
    emax <- exp(lemax)
    ec50 <- exp(lec50)
    hill <- exp(lhill)

    # Fractional sigmoidal Hill (Emax) model, Wahl 2017 Eqn 2:
    # E = E0 + (Emax - E0) * DV^H / (EC50^H + DV^H), Emax = 1.
    # Driver DV is the administered TDF dose (mg/kg); at DOSE = 0 the
    # ratio is 0 and prob_protection = E0.
    prob_protection <-
      e0 + (emax - e0) * DOSE_TDF_MGKG^hill / (ec50^hill + DOSE_TDF_MGKG^hill)

    # Deterministic probability of protection from vaginal HIV
    # acquisition. Downstream callers can sample binary protection with
    # rbinom(n, 1, prob_protection) on the rxSolve output.
    prob_protection ~ add(addSd_prob_protection)
  })
}
