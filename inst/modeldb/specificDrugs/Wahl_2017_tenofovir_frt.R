Wahl_2017_tenofovir_frt <- function() {
  description <- paste0(
    "Preclinical (mouse). Static landmark exposure-response model for the ",
    "fraction of BLT-humanized mice protected from vaginal HIV-1 ",
    "acquisition by systemic tenofovir disoproxil fumarate (TDF) ",
    "pre-exposure prophylaxis, driven by the female reproductive tract ",
    "(FRT) tissue tenofovir (TFV) concentration (Wahl 2017). This is the ",
    "FRT-TFV arm of a five-model PK-PD comparison (dose, plasma TFV, FRT ",
    "TFV, FRT TFV-diphosphate, cervicovaginal lavage TFV). The protected ",
    "fraction is a fractional sigmoidal Hill (Emax) function of the FRT ",
    "TFV concentration, E0 + (Emax - E0) * C^H / (EC50^H + C^H) with Emax ",
    "fixed to 1 and the Hill coefficient H FIXED to 1 (Wahl 2017 Eqn 2, ",
    "Table 2 'FRT TFV' column; H was fixed to 1 for the tissue and CVL ",
    "models because estimating it gave false minima). There is no PK ",
    "layer and no ODE: the FRT TFV concentration enters as the ",
    "CONC_TFV_FRT_NGG covariate column (on the BLT scale), and the ",
    "probability of protection is emitted as the landmark output ",
    "prob_protection. The model was fitted to BLT-equivalent FRT TFV ",
    "obtained by scaling the parallel BALB/c FRT TFV by 0.78 (the ratio ",
    "of BLT to BALB/c median tissue concentrations). Fitted by maximum ",
    "likelihood in ADAPT5."
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
    dosing = "n/a (no dose events; the FRT TFV concentration enters as the CONC_TFV_FRT_NGG covariate column)",
    concentration = "CONC_TFV_FRT_NGG in ng/g of FRT tissue; prob_protection is a probability (0-1)"
  )

  covariateData <- list(
    CONC_TFV_FRT_NGG = list(
      description = "Female reproductive tract (FRT) tissue tenofovir (TFV) concentration driving the protected-fraction model, on the BLT-mouse scale, in ng/g.",
      units = "ng/g",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Landmark FRT (uterus + cervix + vagina) TFV concentration at the ",
        "time of vaginal HIV challenge. The model was FITTED on ",
        "BLT-equivalent concentrations, obtained by multiplying the ",
        "parallel BALB/c FRT TFV by the BLT:BALB/c median ratio 0.78 ",
        "(Wahl 2017 Results 'Comparative TFV and TFVdp concentrations', ",
        "and Fig. 4c which scales FRT TFV by 0.78). This covariate ",
        "therefore carries the BLT-scale concentration: to reproduce a ",
        "dose cohort from the BALB/c measurements in Table S4, multiply ",
        "the BALB/c median by 0.78. BALB/c median FRT TFV by TDF dose ",
        "(Table S4): 20 mg/kg 2804, 50 mg/kg 2780, 140 mg/kg 12034, 300 ",
        "mg/kg 22496.5 ng/g. A value of 0 gives prob_protection = E0. ",
        "Member of the applied-concentration CONC_<DRUG>_<UNITS> family ",
        "(founding example CONC_TFVDP_FMOLMC, Jayachandran 2021); the ",
        "matrix (FRT) is carried in the name to distinguish it from the ",
        "plasma and CVL TFV drivers of this same paper."
      ),
      source_name = "FRT TFV (ng/g)"
    )
  )

  population <- list(
    species = "mouse (BLT humanized, protection endpoint; FRT TFV exposures from parallel BALB/c mice, scaled to the BLT level)",
    n_subjects = "Protection from 78 BLT mice across four TDF arms; FRT TFV from BALB/c mice (n=8 per dose)",
    n_studies = 1L,
    n_observations = "Binary protection outcome per BLT mouse related to the BLT-scaled cohort-median FRT TFV at the matching TDF dose",
    age_range = "Not reported (young adult laboratory mice)",
    sex_female_pct = 100,
    race_ethnicity = "n/a (mouse)",
    disease_state = "HIV-1 pre-exposure prophylaxis challenge model. BLT humanized NSG mice (HIV-susceptible) for protection; BALB/c mice as the FRT-tissue exposure surrogate.",
    dose_range = "TDF 20, 50, 140, or 300 mg/kg intraperitoneally once daily",
    regions = "USA (University of North Carolina at Chapel Hill)",
    notes = paste0(
      "FRT TFV concentrations did not increase dose-proportionally and ",
      "were higher in tissue than in plasma; the tissue models show ",
      "considerable overlap in the concentrations providing 84-100% ",
      "protection (Results 'PK-PD analysis of tenofovir')."
    )
  )

  ini({
    # ================================================================
    # Wahl 2017 Eqn 2 (fractional sigmoidal Hill / Emax) fitted to the
    # BLT-mouse protection data against the BLT-scale FRT TFV
    # concentration, final ADAPT5 estimates from Table 2, 'FRT TFV
    # (ng/g)' column. Values are mean (relative standard error [%SE]).
    # E0 natural scale; EC50 log-transformed. Emax fixed to 1 (all
    # models). The Hill coefficient H is FIXED to 1 for the tissue and
    # CVL models: 'H was fixed to one for TFV in tissue and CVL because
    # the estimation routine resulted in false minima (H < 0.5)' (Methods
    # 'PD model development'); Table 2 reports H as NA for this column.
    # ================================================================
    e0 <- 0.3151
    label("Baseline protected fraction at zero FRT TFV E0 (probability)")  # Wahl 2017 Table 2, FRT TFV column E0 = 0.3151 (%SE 7.299)
    lemax <- fixed(log(1))
    label("Log maximum protected fraction Emax (unitless); Wahl 2017 held Emax at 1 for all models")  # Wahl 2017 Methods 'PD model development': 'Emax was fixed to one for all models'
    lec50 <- log(4101)
    label("Log FRT TFV concentration giving half-maximal protection EC50 (ng/g)")  # Wahl 2017 Table 2, FRT TFV column EC50 = 4101 ng/g (%SE 21.18)
    lhill <- fixed(log(1))
    label("Log Hill coefficient H of the FRT TFV concentration-response curve (unitless)")  # Wahl 2017 Table 2 reports H = NA; Methods: H fixed to 1 for tissue and CVL models (estimating it gave false minima)
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
    # Driver DV is the BLT-scale FRT TFV concentration (ng/g).
    prob_protection <-
      e0 + (emax - e0) * CONC_TFV_FRT_NGG^hill /
        (ec50^hill + CONC_TFV_FRT_NGG^hill)

    prob_protection ~ add(addSd_prob_protection)
  })
}
