Wahl_2017_tenofovir_cvl <- function() {
  description <- paste0(
    "Preclinical (mouse). Static landmark exposure-response model for the ",
    "fraction of BLT-humanized mice protected from vaginal HIV-1 ",
    "acquisition by systemic tenofovir disoproxil fumarate (TDF) ",
    "pre-exposure prophylaxis, driven by the cervicovaginal lavage (CVL) ",
    "tenofovir (TFV) concentration (Wahl 2017). This is the CVL arm of a ",
    "five-model PK-PD comparison (dose, plasma TFV, female reproductive ",
    "tract [FRT] TFV, FRT TFV-diphosphate, CVL TFV) and the least precise ",
    "of the five, reflecting the large (expected) variability in CVL ",
    "concentrations. The protected fraction is a fractional sigmoidal ",
    "Hill (Emax) function of the CVL TFV concentration, E0 + (Emax - E0) ",
    "* C^H / (EC50^H + C^H) with Emax fixed to 1 and the Hill coefficient ",
    "H FIXED to 1 (Wahl 2017 Eqn 2, Table 2 'CVL' column). There is no PK ",
    "layer and no ODE: the CVL TFV concentration enters as the ",
    "CONC_TFV_CVL_NGML covariate column, and the probability of ",
    "protection is emitted as the landmark output prob_protection. No ",
    "BALB/c-to-BLT conversion factor was generated for CVL because of the ",
    "large sample variability, so the BALB/c CVL concentration is used ",
    "directly. Fitted by maximum likelihood in ADAPT5."
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
    time = "n/a (static landmark exposure-response model; no time dimension)",
    dosing = "n/a (no dose events; the CVL TFV concentration enters as the CONC_TFV_CVL_NGML covariate column)",
    concentration = "CONC_TFV_CVL_NGML in ng/mL; prob_protection is a probability (0-1)"
  )

  covariateData <- list(
    CONC_TFV_CVL_NGML = list(
      description = "Cervicovaginal lavage (CVL) fluid tenofovir (TFV) concentration driving the protected-fraction model, in ng/mL.",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Landmark CVL TFV concentration at the time of vaginal HIV ",
        "challenge, measured in parallel BALB/c mice (Wahl 2017 Fig. 3, ",
        "Table S4). No BLT:BALB/c conversion factor was generated for CVL ",
        "because of the large (expected) variability in sample ",
        "concentrations (Results 'Comparative TFV and TFVdp ",
        "concentrations'), so the BALB/c value is used directly. BALB/c ",
        "median CVL TFV by TDF dose (Table S4): 20 mg/kg 323.5, 50 mg/kg ",
        "266.3, 140 mg/kg 43900, 300 mg/kg 5870 ng/mL (highly variable ",
        "and non-monotone in dose). A value of 0 gives prob_protection = ",
        "E0. Member of the applied-concentration CONC_<DRUG>_<UNITS> ",
        "family; the matrix (CVL) is carried in the name to distinguish ",
        "it from the plasma TFV driver of this same paper ",
        "(CONC_TFV_PLASMA_NGML), which has identical units."
      ),
      source_name = "CVL TFV (ng/ml)"
    )
  )

  population <- list(
    species = "mouse (BLT humanized, protection endpoint; CVL TFV exposures from parallel BALB/c mice)",
    n_subjects = "Protection from 78 BLT mice across four TDF arms; CVL TFV from BALB/c mice (n=8 per dose)",
    n_studies = 1L,
    n_observations = "Binary protection outcome per BLT mouse related to the cohort-median BALB/c CVL TFV at the matching TDF dose",
    age_range = "Not reported (young adult laboratory mice)",
    sex_female_pct = 100,
    race_ethnicity = "n/a (mouse)",
    disease_state = "HIV-1 pre-exposure prophylaxis challenge model. BLT humanized NSG mice (HIV-susceptible) for protection; BALB/c mice as the CVL-fluid exposure surrogate.",
    dose_range = "TDF 20, 50, 140, or 300 mg/kg intraperitoneally once daily",
    regions = "USA (University of North Carolina at Chapel Hill)",
    notes = paste0(
      "CVL TFV concentrations were highly variable and not dose ",
      "proportional; the CVL model was the least precise of the five ",
      "(largest %SE on EC50), consistent with the known 1-log ",
      "variability of human CVL data (Results and Discussion)."
    )
  )

  ini({
    # ================================================================
    # Wahl 2017 Eqn 2 (fractional sigmoidal Hill / Emax) fitted to the
    # BLT-mouse protection data against the CVL TFV concentration, final
    # ADAPT5 estimates from Table 2, 'CVL (ng/ml)' column. Values are
    # mean (relative standard error [%SE]). E0 natural scale; EC50
    # log-transformed. Emax fixed to 1 (all models). The Hill coefficient
    # H is FIXED to 1 for the tissue and CVL models (Table 2 reports H =
    # NA; estimating it gave false minima, Methods 'PD model
    # development').
    # ================================================================
    e0 <- 0.4016
    label("Baseline protected fraction at zero CVL TFV E0 (probability)")  # Wahl 2017 Table 2, CVL column E0 = 0.4016 (%SE 8.918)
    lemax <- fixed(log(1))
    label("Log maximum protected fraction Emax (unitless); Wahl 2017 held Emax at 1 for all models")  # Wahl 2017 Methods 'PD model development': 'Emax was fixed to one for all models'
    lec50 <- log(238.7)
    label("Log CVL TFV concentration giving half-maximal protection EC50 (ng/mL)")  # Wahl 2017 Table 2, CVL column EC50 = 238.7 ng/mL (%SE 51.53)
    lhill <- fixed(log(1))
    label("Log Hill coefficient H of the CVL TFV concentration-response curve (unitless)")  # Wahl 2017 Table 2 reports H = NA; Methods: H fixed to 1 for tissue and CVL models (estimating it gave false minima)
    # ================================================================
    # Placeholder additive residual for the deterministic landmark
    # probability; Wahl 2017 used a proportional error whose final
    # estimate it did not report. Not from source.
    # ================================================================
    addSd_prob_protection <- fixed(0.001)
    label("Placeholder additive residual SD on the typical-value protected fraction; Wahl 2017 used a proportional error whose final estimate it did not report")  # not from source; see vignette Assumptions and deviations
  })

  model({
    emax <- exp(lemax)
    ec50 <- exp(lec50)
    hill <- exp(lhill)

    # Fractional sigmoidal Hill (Emax) model, Wahl 2017 Eqn 2:
    # E = E0 + (Emax - E0) * DV^H / (EC50^H + DV^H), Emax = 1, H = 1.
    # Driver DV is the CVL TFV concentration (ng/mL).
    prob_protection <-
      e0 + (emax - e0) * CONC_TFV_CVL_NGML^hill /
        (ec50^hill + CONC_TFV_CVL_NGML^hill)

    prob_protection ~ add(addSd_prob_protection)
  })
}
