Wahl_2017_tenofovir_plasma <- function() {
  description <- paste0(
    "Preclinical (mouse). Static landmark exposure-response model for the ",
    "fraction of BLT-humanized mice protected from vaginal HIV-1 ",
    "acquisition by systemic tenofovir disoproxil fumarate (TDF) ",
    "pre-exposure prophylaxis, driven by the PLASMA tenofovir (TFV) ",
    "concentration (Wahl 2017). This is the plasma arm of a five-model ",
    "PK-PD comparison (dose, plasma TFV, female reproductive tract [FRT] ",
    "TFV, FRT TFV-diphosphate, cervicovaginal lavage TFV); plasma TFV was ",
    "judged the best predictor of protective efficacy across the dose ",
    "range. The protected fraction is a fractional sigmoidal Hill (Emax) ",
    "function of the plasma TFV concentration, E0 + (Emax - E0) * ",
    "C^H / (EC50^H + C^H) with Emax fixed to 1 (Wahl 2017 Eqn 2, Table 2 ",
    "'Plasma' column). There is no PK layer and no ODE: the plasma TFV ",
    "concentration is supplied as the CONC_TFV_PLASMA_NGML covariate ",
    "column, and the probability of protection is emitted as the landmark ",
    "output prob_protection. BLT and BALB/c plasma TFV concentrations ",
    "differed by less than 9% (not significant), so no BALB/c-to-BLT ",
    "conversion factor was applied to this driver. Fitted by maximum ",
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
    dosing = "n/a (no dose events; the plasma TFV concentration enters as the CONC_TFV_PLASMA_NGML covariate column)",
    concentration = "CONC_TFV_PLASMA_NGML in ng/mL; prob_protection is a probability (0-1)"
  )

  covariateData <- list(
    CONC_TFV_PLASMA_NGML = list(
      description = "Plasma tenofovir (TFV) concentration driving the protected-fraction model, in ng/mL.",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Landmark plasma TFV concentration measured 3 h after the third ",
        "daily TDF dose (the time of vaginal HIV challenge), from BALB/c ",
        "mice given the same TDF doses (Wahl 2017 Fig. 3, Table S4). ",
        "BLT and BALB/c plasma TFV differed by <9% (p = 0.43), so the ",
        "BALB/c value is used directly with NO conversion factor, unlike ",
        "the FRT TFV and FRT TFVdp drivers. BALB/c median plasma TFV by ",
        "TDF dose (Table S4): 20 mg/kg 184.5, 50 mg/kg 370.5, 140 mg/kg ",
        "1500, 300 mg/kg 2410 ng/mL. A value of 0 gives prob_protection = ",
        "E0. Member of the applied-concentration CONC_<DRUG>_<UNITS> ",
        "family (founding example CONC_TFVDP_FMOLMC, Jayachandran 2021); ",
        "the matrix (plasma) is carried in the name because this paper ",
        "also uses a cervicovaginal-lavage TFV driver of the same units ",
        "(CONC_TFV_CVL_NGML), which the family's bare drug+units spelling ",
        "could not disambiguate."
      ),
      source_name = "Plasma TFV (ng/ml)"
    )
  )

  population <- list(
    species = "mouse (BLT humanized, protection endpoint; plasma TFV exposures from parallel BALB/c mice)",
    n_subjects = "Protection from 78 BLT mice across four TDF arms; plasma TFV from BALB/c mice (n=8 per dose)",
    n_studies = 1L,
    n_observations = "Binary protection outcome per BLT mouse related to the cohort-median BALB/c plasma TFV at the matching TDF dose",
    age_range = "Not reported (young adult laboratory mice)",
    sex_female_pct = 100,
    race_ethnicity = "n/a (mouse)",
    disease_state = "HIV-1 pre-exposure prophylaxis challenge model. BLT humanized NSG mice (HIV-susceptible) for protection; BALB/c mice as the plasma-exposure surrogate.",
    dose_range = "TDF 20, 50, 140, or 300 mg/kg intraperitoneally once daily",
    regions = "USA (University of North Carolina at Chapel Hill)",
    notes = paste0(
      "At 50 mg/kg TDF the plasma TFV (~370 ng/mL) mimicked human ",
      "exposures and gave a ~70% risk reduction, comparable to highly ",
      "adherent women in Partners PrEP. Predicted protection at the ",
      "cohort-median plasma TFV concentrations: 58, 67, 84, 88% for 20, ",
      "50, 140, 300 mg/kg (Results 'PK-PD analysis of tenofovir')."
    )
  )

  ini({
    # ================================================================
    # Wahl 2017 Eqn 2 (fractional sigmoidal Hill / Emax) fitted to the
    # BLT-mouse protection data against the plasma TFV concentration,
    # final ADAPT5 estimates from Table 2, 'Plasma (ng/ml)' column.
    # Values are mean (relative standard error [%SE]). E0 is a baseline
    # probability on the natural scale; EC50 and H are log-transformed to
    # stay positive. Emax fixed to 1 for all of the paper's models.
    # ================================================================
    e0 <- 0.2580
    label("Baseline protected fraction at zero plasma TFV E0 (probability)")  # Wahl 2017 Table 2, Plasma column E0 = 0.2580 (%SE 13.24)
    lemax <- fixed(log(1))
    label("Log maximum protected fraction Emax (unitless); Wahl 2017 held Emax at 1 for all models")  # Wahl 2017 Methods 'PD model development': 'Emax was fixed to one for all models'
    lec50 <- log(259)
    label("Log plasma TFV concentration giving half-maximal protection EC50 (ng/mL)")  # Wahl 2017 Table 2, Plasma column EC50 = 259 ng/mL (%SE 19.6)
    lhill <- log(0.7449)
    label("Log Hill coefficient H of the plasma concentration-response curve (unitless)")  # Wahl 2017 Table 2, Plasma column H = 0.7449 (%SE 14.98); estimated for the dose and plasma models
    # ================================================================
    # Placeholder additive residual for the deterministic landmark
    # probability; Wahl 2017 used a proportional error whose final
    # estimate it did not report. Not from source (prob_<endpoint>
    # landmark convention).
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
    # Driver DV is the plasma TFV concentration (ng/mL).
    prob_protection <-
      e0 + (emax - e0) * CONC_TFV_PLASMA_NGML^hill /
        (ec50^hill + CONC_TFV_PLASMA_NGML^hill)

    prob_protection ~ add(addSd_prob_protection)
  })
}
