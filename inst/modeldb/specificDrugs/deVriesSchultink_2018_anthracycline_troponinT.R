deVriesSchultink_2018_anthracycline_troponinT <- function() {
  description <- paste(
    "Kinetic-pharmacodynamic (K-PD) direct-effect model for serum",
    "high-sensitive cardiac troponin T (hs-TnT) in early-breast-cancer",
    "patients receiving an adjuvant anthracycline regimen (de Vries",
    "Schultink 2018). The dose enters a virtual K-PD body amount",
    "compartment with first-order elimination at rate kel; the linear",
    "direct effect TRP = TRP0 * (1 + SLOPE * Aant) raises serum troponin",
    "T above a population baseline TRP0 in proportion to the current",
    "K-PD amount. The proportional SLOPE is anthracycline-type dependent:",
    "epirubicin produces a 0.524-fold smaller effect than the doxorubicin",
    "reference, encoded jointly by the CONMED_DOXORUBICIN and",
    "CONMED_EPIRUBICIN indicators. No other covariates retained.",
    "Companion file `deVriesSchultink_2018_trastuzumab_LVEF.R` consumes",
    "the per-subject peak troponin T from this model as a covariate."
  )
  reference <- paste(
    "de Vries Schultink AHM, Boekhout AH, Gietema JA, Burylo AM, Dorlo",
    "TPC, van Hasselt JGC, Schellens JHM, Huitema ADR. Pharmacodynamic",
    "modeling of cardiac biomarkers in breast cancer patients treated",
    "with anthracycline and trastuzumab regimens. J Pharmacokinet",
    "Pharmacodyn. 2018;45(3):431-442. doi:10.1007/s10928-018-9579-8"
  )
  vignette <- "deVriesSchultink_2018_cardiotoxicity"
  units <- list(time = "day", dosing = "mg", concentration = "ng/L")

  covariateData <- list(
    CONMED_DOXORUBICIN = list(
      description        = "Indicator that the subject's anthracycline cycles used doxorubicin",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (subject not on doxorubicin; paired with CONMED_EPIRUBICIN = 1 in this cohort)",
      notes              = "Time-fixed per subject in this study (each patient received exactly one anthracycline class for the full 2-6 cycle course). The de Vries Schultink 2018 dataset records the anthracycline class as a categorical 'doxorubicin' / 'epirubicin'; the canonical covariate register splits this into the two binary indicators CONMED_DOXORUBICIN and CONMED_EPIRUBICIN per the CONMED_<drug> pattern. Doxorubicin is the reference category in the paper's model (paper's anthracycline-type-on-SLOPE coefficient is the epirubicin-vs-doxorubicin contrast). Per Table 1, doxorubicin n = 181 (87.9%), epirubicin n = 25 (12.1%).",
      source_name        = "ANTH_TYPE (= 'doxorubicin')"
    ),
    CONMED_EPIRUBICIN = list(
      description        = "Indicator that the subject's anthracycline cycles used epirubicin",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (subject not on epirubicin)",
      notes              = "See CONMED_DOXORUBICIN. The two indicators are mutually exclusive in this cohort (exactly one is 1 per subject). The epirubicin-on-SLOPE multiplicative effect is 0.524 (de Vries Schultink 2018 Table 2, 'Proportional anthracycline-type effect on SLOPE', RSE 17.5%); the doxorubicin-driven SLOPE is the population reference. The text-worked example (paper p. 438) confirms: 100 mg doxorubicin raises troponin T from 3 to 5.6 ng/L (an increment of 2.6 ng/L driven by SLOPE * dose at t = 0), and an equivalent epirubicin dose raises troponin T from 3 to 4.3 ng/L (an increment of 1.3 ng/L = 0.524 * 2.6). The encoding here is SLOPE_i = SLOPE_pop * (1 + e_anth_slope * CONMED_EPIRUBICIN - e_anth_slope * CONMED_DOXORUBICIN) = SLOPE_pop * (CONMED_DOXORUBICIN + 0.524 * CONMED_EPIRUBICIN) when the two indicators are mutually exclusive.",
      source_name        = "ANTH_TYPE (= 'epirubicin')"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 190L,
    n_studies      = 1L,
    n_observations = 1230L,
    age_range      = "25-69 years",
    age_median     = "50 years",
    weight_range   = NA_character_,
    weight_median  = NA_character_,
    sex_female_pct = 100,
    race_ethnicity = "Not explicitly tabulated in de Vries Schultink 2018; clinical-trial cohort recruited at Dutch and University-of-Groningen sites (predominantly European white).",
    disease_state  = "HER2-positive early (stage I-III) breast cancer; cardiac-marker-evaluable subset (n = 190 of 206 randomised) of the parent candesartan-vs-placebo cardioprotection trial (Boekhout 2016, JAMA Oncology 2:1030-1037, doi:10.1001/jamaoncol.2016.1726). Patients had favorable baseline cardiac history (max age 69 y; only 13% > 60 y) and no significant prior cardiac disease.",
    dose_range     = "Adjuvant anthracycline as 4 [2-6] cycles of either doxorubicin (median absolute dose 110 mg/cycle [75-150], n = 181) or epirubicin (median 170 mg/cycle [100-200], n = 25), administered IV at three-week intervals; cumulative doses below the lifetime cardiotoxicity threshold (< 550 mg/m^2 for doxorubicin and < 950 mg/m^2 for epirubicin).",
    regions        = "Netherlands (Antoni van Leeuwenhoek / NKI Amsterdam, University Medical Center Groningen).",
    notes          = "Cohort and patient characteristics from de Vries Schultink 2018 Table 1 (the 'troponin T measurements available' subset, n = 190 of 206). 4.6% of the 1230 troponin T observations were below the assay LLOQ (3 ng/L for the Roche Modular E hs-TnT assay); BLQ values were imputed as LLOQ/2 in the source dataset per Methods."
  )

  ini({
    # Structural parameters (de Vries Schultink 2018 Table 2, anthracycline-
    # troponin T model rows). The K-PD body amount compartment is named
    # depot_kpd (canonical) and the dose enters the compartment directly in
    # mg. Troponin T is computed as an algebraic output troponin_t = TRP0 *
    # (1 + SLOPE * depot_kpd) per the paper's Eq. block on p. 435:
    #     d/dt(Aant) = -ke * Aant
    #     TRP = TRP0 * (1 + SLOPE * Aant)
    # Doxorubicin is the reference anthracycline; epirubicin scales SLOPE by
    # the multiplicative factor of 0.524 (Table 2, "Proportional
    # anthracycline-type effect on SLOPE", RSE 17.5%).
    ltrp0  <- log(4.72);    label("Baseline serum troponin T concentration (ng/L)")             # Table 2: TRP0 = 4.72 ng/L, RSE 3.5%
    lkel   <- log(8.49e-3); label("K-PD body-amount elimination rate constant (1/day)")          # Table 2: ke = 8.49e-3 /day, RSE 4.0%
    lslope <- log(8.84e-3); label("Doxorubicin direct-effect SLOPE on troponin T (1 / mg of K-PD anthracycline amount)") # Table 2: SLOPE = 8.84e-3 (ng/L)^-1, but the equation reads TRP = TRP0 * (1 + SLOPE * Aant) so SLOPE has units 1/Aant = 1/mg

    # Anthracycline-type covariate effect on SLOPE: encoded as the
    # epirubicin-vs-doxorubicin multiplicative ratio. Since the paper's
    # Table 2 reports 0.524 (i.e., epirubicin SLOPE is 0.524-fold of
    # doxorubicin SLOPE), the covariate-effect parameter `e_epirubicin_slope`
    # carries the absolute multiplier 0.524 directly. SLOPE_i is then
    # selected via the binary indicators with doxorubicin as reference:
    #     slope_i = slope_pop * (CONMED_DOXORUBICIN + e_epirubicin_slope * CONMED_EPIRUBICIN)
    # When CONMED_DOXORUBICIN = 1 (and CONMED_EPIRUBICIN = 0) the multiplier
    # is 1; when CONMED_EPIRUBICIN = 1 (and CONMED_DOXORUBICIN = 0) the
    # multiplier is 0.524.
    e_epirubicin_slope <- 0.524; label("Epirubicin-vs-doxorubicin multiplicative ratio on SLOPE (unitless; doxorubicin is the reference)") # Table 2: anthracycline-type effect 0.524, RSE 17.5%

    # Inter-individual variability (Table 2; CV%). Paper reports
    # log-normal IIV on SLOPE (57.7% CV, shrinkage 31.0%) and TRP0
    # (39.2% CV, shrinkage 12.6%). No off-diagonal covariance reported.
    # omega^2 = log(CV^2 + 1) under the log-normal convention.
    etalslope ~ 0.298017                                                                          # Table 2: omega_SLOPE = 57.7% CV -> omega^2 = log(0.577^2 + 1)
    etaltrp0  ~ 0.142567                                                                          # Table 2: omega_TRP0  = 39.2% CV -> omega^2 = log(0.392^2 + 1)

    # Residual error: proportional 30.1% on troponin T (Table 2,
    # shrinkage 11.2%). Single-output model -> use the bare-suffix
    # `propSd` form per the canonical convention.
    propSd <- 0.301; label("Proportional residual SD on troponin T (fraction)")                  # Table 2: proportional residual error = 30.1%
  })

  model({
    # Individual parameters (log-normal IIV; population point estimates
    # exponentiated to the linear scale). Doxorubicin is the reference
    # category; epirubicin scales SLOPE by the 0.524 multiplier.
    trp0    <- exp(ltrp0 + etaltrp0)
    kel     <- exp(lkel)
    slope_i <- exp(lslope + etalslope) *
                 (CONMED_DOXORUBICIN + e_epirubicin_slope * CONMED_EPIRUBICIN)

    # K-PD virtual body amount compartment with first-order elimination
    # (de Vries Schultink 2018 p. 435, dAant/dt equation). The dose enters
    # depot_kpd directly in mg per cycle of anthracycline.
    d/dt(depot_kpd) <- -kel * depot_kpd

    # Direct-effect troponin T (de Vries Schultink 2018 p. 435,
    # TRP = TRP0 * (1 + SLOPE * Aant) equation). Algebraic output with no
    # ODE state -- troponin_t is the observable concentration in ng/L.
    troponin_t <- trp0 * (1 + slope_i * depot_kpd)

    troponin_t ~ prop(propSd)
  })
}
