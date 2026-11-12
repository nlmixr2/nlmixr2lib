Zecchin_2026_depemokimab_dose_time <- function() {
  description <- "Bayesian nonlinear mixed-effects dose-time response model of log blood eosinophil count after a single subcutaneous depemokimab dose (2-300 mg) or placebo in adults with mild-to-moderate asthma (Zecchin 2026, Table 1). Algebraic (no PK): the depemokimab arm follows a sigmoid in log dose whose half-maximal dose rises exponentially with time after the dose (log ED50 at Day 182 plus KE x (day - 182)), the placebo arm is a separate log-linear baseline regression, both arms use log baseline BEC as a covariate, the 2 mg dose has its own relative potency F, and the subject-level intercept, Emax and log ED50 random effects are correlated. Valid for a single dose only (time = days after the dose)."
  reference <- "Zecchin C, Schalkwijk S, Pouliquen IJ, Berges A, Bird N, Follows R, Austin D. Straight to Phase III: Model-Informed Approach Speeds Depemokimab Clinical Development in Interleukin-5-Driven Diseases. Clin Pharmacol Ther. 2026;119(4):1005. doi:10.1002/cpt.70183"
  vignette <- "Zecchin_2026_depemokimab"
  units <- list(time = "day", dosing = "mg", concentration = "10^9 cells/L")

  covariateData <- list(
    EOS = list(
      description = "Baseline blood eosinophil count (time-fixed per subject)",
      units = "cells/uL",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters both arms as LNBASE = log(baseline BEC in GI/L) = log(EOS / 1000) (Zecchin 2026 Table 1 footnote: 'LNBASE, log(Baseline)'; the GI/L scale of LNEOS is read from the Figure 1 y-axis 'Absolute BEC, GI/L'). The FTiH cohort was enrolled at BEC >= 200 cells/uL at screening.",
      source_name = "LNBASE (log baseline BEC, GI/L)"
    ),
    DOSE_DEPEMOKIMAB_MG = list(
      description = "Single subcutaneous depemokimab dose of the subject's arm, in mg; 0 = placebo",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "Enters as LNDOSE = log(dose in mg). A value of 0 selects the placebo equation (LNEOS = b1 + E1 + E3 * LNBASE). The 2 mg dose additionally carries the relative potency F (Table 1 'F (2 mg dose)'), applied as -log(F) inside the sigmoid. Phase I doses: 2, 10, 30, 100, 300 mg.",
      source_name = "DOSE (LNDOSE = log(DOSE))"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 48L,
    n_studies = 1L,
    age_range = "19-65 years (mean 44.0, SD 11.2; Thorsted 2026 Table 3, same 48 patients)",
    weight_range = "57.8-107.0 kg (mean 82.9, SD 11.6; Thorsted 2026 Table 3)",
    sex_female_pct = NA_real_,
    race_ethnicity = NA_character_,
    disease_state = "Mild-to-moderate asthma on low-to-moderate inhaled corticosteroid (+/- LABA) with blood eosinophil count >= 200 cells/uL at screening (baseline mean 366, SD 147 cells/uL; Thorsted 2026 Table 3)",
    dose_range = "Single subcutaneous depemokimab 2, 10, 30, 100 or 300 mg, or placebo",
    regions = "Europe (first-time-in-human Study 205722)",
    notes = "622 BEC measurements from 48 patients in the depemokimab first-time-in-human Phase I trial; the model was fitted to the 293 observations from Day 57 onward (return-to-baseline phase) and externally checked on the 47 Day 29 observations (Zecchin 2026 Results). Estimation: SAS Bayesian procedure (MCMC, 10,000 retained samples); the values here are the posterior means of Table 1."
  )

  ini({
    # Posterior means, Zecchin 2026 Table 1 (N = 10,000 MCMC samples each).
    # LNEOS = log(BEC in GI/L); LNBASE = log(baseline BEC in GI/L); LNDOSE = log(dose in mg).
    e0 <- -0.2601; label("Active-arm intercept of log BEC away from the drug effect, E0 (log GI/L)") # Table 1 'E0' mean -0.2601 (SD 0.3054)
    e0_placebo <- -0.6583; label("Placebo-arm intercept of log BEC, E1 (log GI/L)") # Table 1 'E1' mean -0.6583 (SD 0.4209)
    e_eos_e0 <- 0.5128; label("Slope of log baseline BEC in the active arm, E2 (unitless)") # Table 1 'E2' mean 0.5128 (SD 0.1995)
    e_eos_e0_placebo <- 0.5084; label("Slope of log baseline BEC in the placebo arm, E3 (unitless)") # Table 1 'E3' mean 0.5084 (SD 0.3865); median 0.3425
    emax <- -2.6047; label("Intercept of log BEC at full drug effect, Emax (log GI/L)") # Table 1 'E max' mean -2.6047 (SD 0.2451)
    led50 <- 2.9513; label("Log half-maximal dose at Day 182 after the dose, LNED50 (log mg)") # Table 1 'LNED50 (Day 182)' mean 2.9513 (SD 0.2093); ED50 = exp(2.9513) = 19.13 mg
    ked50 <- 0.0207; label("Rate of increase of log ED50 with time after the dose, KE (1/day)") # Table 1 'KE (/day)' mean 0.0207 (SD 0.00263); T1/2 = log(2)/KE = 33.5 days
    lhill <- log(1.2358); label("Hill slope of the sigmoid in log dose, SLP (unitless)") # Table 1 'Hill slope' mean 1.2358 (SD 0.2222)
    lfrel_dose2mg <- log(2.3730); label("Relative potency (bioavailability) of the 2 mg dose, F (unitless)") # Table 1 'F (2 mg dose)' mean 2.3730 (SD 0.9708)

    # Between-patient random effects b1 (intercept, both arms), b2 (Emax) and
    # b3 (LNED50); Table 1 S2B rows are the variances / covariances.
    etae0 + etaemax + etaled50 ~ c(
      0.2096,
      -0.2058, 0.5194,
      0.0550, -0.2120, 0.3850
    ) # Table 1: S2B1 (E0) 0.2096; S2B1B2 = S2B2B1 -0.2058; S2B2 (E max) 0.5194; S2B1B3 = S2B3B1 0.0550; S2B2B3 = S2B3B2 -0.2120; S2B3 (LNED50) 0.3850

    # EPS is additive on log BEC with variance 0.0804 (Table 1 'S2Y (EPS)'),
    # reported as BOV = sqrt(exp(0.0804) - 1) = 28.9%; SD = sqrt(0.0804).
    expSd <- 0.28355; label("Residual SD of log BEC (log GI/L)") # Table 1 'S2Y (EPS)' mean 0.0804 is a variance -> sqrt(0.0804) = 0.28355
  })

  model({
    lnbase <- log(EOS / 1000)
    hill <- exp(lhill)
    frel_dose2mg <- exp(lfrel_dose2mg)

    # Relative potency F applies to the 2 mg dose only (Table 1 footnote
    # 'F, relative bioavailability of 2 mg dose'); every other dose has F = 1.
    frel <- 1
    if (DOSE_DEPEMOKIMAB_MG == 2) {
      frel <- frel_dose2mg
    }

    # Table 1 'Model' row: time after the dose in days (TDAYS = t), centred on
    # Day 182 for LNED50.
    # b1 (etae0) shifts both arms; b2 (etaemax) and b3 (etaled50) act on
    # the active arm only.
    e0i <- e0 + etae0
    b1 <- e0i - e0
    emaxi <- emax + etaemax
    led50i <- led50 + etaled50
    if (DOSE_DEPEMOKIMAB_MG > 0) {
      lneos <- e0i + e_eos_e0 * lnbase +
        (emaxi - e0) /
          (1 + exp(hill * (led50i + ked50 * (t - 182) - log(frel) - log(DOSE_DEPEMOKIMAB_MG))))
    } else {
      lneos <- b1 + e0_placebo + e_eos_e0_placebo * lnbase
    }

    blood_eosinophils <- exp(lneos)
    blood_eosinophils ~ lnorm(expSd)
  })
}
