Homkham_2019_efavirenz <- function() {
  description <- paste(
    "Population PK-PD HIV-1 viral-dynamics model of efavirenz in Thai",
    "HIV-infected children starting a first-line efavirenz-based",
    "antiretroviral regimen (PHPT Observational Cohort, Homkham 2019).",
    "A three-state target-cell model (Ribeiro and Perelson) of uninfected",
    "CD4+ target cells `target` (T), productively infected cells `infected`",
    "(T*) and free virions `virus` (V), fitted in Monolix to HIV-1 RNA load",
    "and CD4 cell count over the first three years of treatment on a",
    "month time base. Before treatment the system sits at its chronic",
    "infection steady state; efavirenz inhibits the infection rate beta",
    "through an Imax function I = 1 - Imax * Cav / (IC50 + Cav) of the",
    "individual average efavirenz plasma concentration Cav = (dose / tau) /",
    "CL, supplied as the CAV covariate (mg/L) from the upstream Homkham",
    "2015 population PK model; Imax is fixed to 1 and the infected-cell",
    "death rate delta is fixed to 15.2 per month. Lognormal between-subject",
    "variability on lambda, d (correlated), beta, p, c and IC50; additive",
    "residual error on log10 HIV-1 RNA and on CD4 count. No covariate was",
    "retained on any PD parameter.",
    sep = " "
  )
  reference <- paste(
    "Homkham N, Cressey TR, Bouazza N, Ingsrisawang L, Techakunakorn P,",
    "Mekmullica J, Borkird T, Puangsombat A, Na-Rajsima S, Treluyer JM,",
    "Urien S, Jourdain G. Role of efavirenz plasma concentrations on",
    "long-term HIV suppression and immune restoration in HIV-infected",
    "children. PLoS One. 2019;14(5):e0216868.",
    "doi:10.1371/journal.pone.0216868.",
    "Structural equations are Eqs 1-4 and the pre-treatment steady state",
    "of Methods 'Population pharmacokinetic analysis'; parameter estimates",
    "are Table 2.",
    sep = " "
  )
  vignette <- "Homkham_2019_efavirenz"
  units <- list(
    time = "month",
    dosing = "mg",
    concentration = paste(
      "No drug-concentration state: efavirenz exposure enters as the CAV",
      "covariate in mg/L. log10_viral_load is log10 HIV-1 RNA copies/mL;",
      "target (the CD4 cell count output) is cells/mm3."
    )
  )

  covariateData <- list(
    CAV = list(
      description = "Individual average steady-state efavirenz plasma concentration over the once-daily dosing interval",
      units = "mg/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Cav = (dose / dosing interval) / CL_i, with CL_i the individual",
        "efavirenz clearance estimated from the previously published",
        "population PK model of the same cohort (Homkham 2015, Pediatr",
        "Infect Dis J 34:1214; one-compartment model with two transit",
        "absorption compartments, allometric weight scaling and a",
        "fast/slow-metaboliser mixture on CL/F). One constant value per",
        "child from the start of treatment (time 0); set CAV = 0 before",
        "treatment, which leaves the system at its untreated chronic",
        "infection steady state. Cohort median 1.7 mg/L (IQR 1.3 to 2.2;",
        "Results 'Efavirenz plasma concentrations and viral replication')."
      ),
      source_name = "Cav"
    )
  )

  compartmentData <- list(
    target = list(
      analyte = "uninfected CD4+ T lymphocytes (target cells)",
      units = "cells/mm3",
      specimen = "whole blood",
      verified = TRUE
    ),
    infected = list(
      analyte = "productively HIV-1-infected CD4+ T lymphocytes",
      units = "cells/mm3",
      specimen = "whole blood",
      verified = TRUE
    ),
    virus = list(
      analyte = "free HIV-1 virions (HIV-1 RNA)",
      units = "copies/mL",
      specimen = "plasma",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 87L,
    n_studies = 1L,
    age_range = "children; IQR 5 to 10 years",
    age_median = "8 years",
    weight_range = "IQR 14 to 23 kg",
    weight_median = "17 kg",
    sex_female_pct = 49,
    race_ethnicity = "Thai",
    disease_state = paste(
      "HIV-1-infected, antiretroviral-naive children starting first-line",
      "efavirenz plus two NRTIs (zidovudine n = 57 or stavudine n = 30,",
      "plus lamivudine n = 84 or didanosine n = 3); baseline HIV-1 RNA",
      "median 5.1 log10 copies/mL (IQR 4.6 to 5.4), CD4 count median 71",
      "cells/mm3 (IQR 23 to 262)"
    ),
    dose_range = paste(
      "Efavirenz once daily in the evening by US FDA weight bands;",
      "median 13.9 mg/kg (IQR 13.0 to 16.2)"
    ),
    regions = "Thailand",
    notes = paste(
      "Program for HIV Prevention and Treatment (PHPT) Observational Cohort",
      "(NCT00433030), first-line efavirenz started between May 2002 and",
      "March 2010; children on interacting drugs such as rifampicin were",
      "excluded. Median follow-up 24 months (IQR 12 to 34); 443 HIV-1 RNA",
      "loads and 447 CD4 counts. Baseline characteristics are Table 1",
      "(whose printed age IQR '5 to 6' contradicts the abstract's '5 to",
      "10'; the abstract value is used here)."
    )
  )

  ini({
    # Fixed effects (Table 2). Rates are per month, the Monolix time base.
    lkin <- log(89.2); label("Production rate of uninfected target cells lambda (cells/mm3/month)") # Table 2 'lambda' = 89.2 cells/mm3/month (RSE 8%)
    lkout <- log(0.1); label("Death rate constant of uninfected target cells d (1/month)") # Table 2 'd' = 0.1 (RSE 11%); Table 2 prints the unit as cells/mm3/month, but Fig 1 and Eq 1 make d a first-order rate constant
    lbeta <- log(6.7e-6); label("Infection rate constant of target cells beta (mL/copy/month)") # Table 2 'beta x10^-6' = 6.7 (RSE 17%)
    ldelta <- fixed(log(15.2)); label("Death rate constant of infected cells delta (1/month)") # Table 2 'delta (fixed)' = 15.2 per month; Results: fixed from Perelson 1996
    lp <- log(6.5e4); label("Virion production rate per infected cell p ((copies/mL)/(cells/mm3)/month)") # Table 2 'p x10^4' = 6.5 (RSE 10%)
    lc <- log(3.3); label("Clearance rate constant of free virions c (1/month)") # Table 2 'c' = 3.3 per month (RSE 10%)
    imax <- fixed(1); label("Maximum fractional inhibition of beta by efavirenz Imax (fraction)") # Table 2 'Imax (fixed)' = 1.0; Results: fixed to 100% (not identifiable)
    lic50 <- log(0.1); label("Efavirenz Cav at half-maximal inhibition IC50 (mg/L)") # Table 2 'IC50' = 0.1 mg/L (RSE 47%)

    # Between-subject variability. Table 2 omegas are Monolix standard
    # deviations of the lognormal random effects, so variance = omega^2;
    # 'omega lambda,d' = 0.6 is the correlation, so the covariance is
    # 0.6 * 0.6 * 0.6 = 0.216 (read as a covariance it would imply a
    # correlation of 1.67, which is impossible).
    etalkin + etalkout ~ c(
      0.36,
      0.216, 0.36
    ) # Table 2 'omega lambda' = 0.6, 'omega d' = 0.6, 'omega lambda,d' = 0.6 (correlation)
    etalbeta ~ 1.0 # Table 2 'omega beta' = 1.0 -> variance 1.0
    etalp ~ 0.16 # Table 2 'omega p' = 0.4 -> variance 0.16
    etalc ~ 0.25 # Table 2 'omega c' = 0.5 -> variance 0.25
    etalic50 ~ 4.41 # Table 2 'omega IC50' = 2.1 -> variance 4.41

    # Residual error (Table 2): additive on both outputs.
    addSd_log10_viral_load <- 0.7; label("Additive residual error on log10 HIV-1 RNA (log10 copies/mL)") # Table 2 'sigma VL' = 0.7 log10 copies/mL (RSE 4%)
    addSd_target <- 176; label("Additive residual error on CD4 cell count (cells/mm3)") # Table 2 'sigma CD4' = 176.0 cells/mm3 (RSE 4%)
  })

  model({
    # Individual parameters
    kin <- exp(lkin + etalkin)
    kout <- exp(lkout + etalkout)
    beta <- exp(lbeta + etalbeta)
    delta <- exp(ldelta)
    p <- exp(lp + etalp)
    cvir <- exp(lc + etalc)
    ic50 <- exp(lic50 + etalic50)

    # Eq 4: efavirenz inhibition of the infection rate; CAV = 0 before
    # treatment gives inh = 1 (no drug effect).
    inh <- 1 - imax * CAV / (ic50 + CAV)

    # Untreated chronic-infection steady state (Methods, after Fig 1):
    # T0 = c * delta / (beta * p), T*0 = c * VL0 / p,
    # VL0 = lambda * p / (c * delta) - d / beta. VL0 is positive only when
    # the basic reproductive ratio lambda * beta * p / (d * c * delta) > 1.
    target_bl <- cvir * delta / (beta * p)
    virus_bl <- kin * p / (cvir * delta) - kout / beta
    infected_bl <- cvir * virus_bl / p
    target(0) <- target_bl
    infected(0) <- infected_bl
    virus(0) <- virus_bl

    # Eqs 1-3
    d/dt(target) <- kin - kout * target - inh * beta * target * virus
    d/dt(infected) <- inh * beta * target * virus - delta * infected
    d/dt(virus) <- p * infected - cvir * virus

    # Inhibitory score ISEFV (Eq 5), percent
    isefv <- 100 * imax * CAV / (ic50 + CAV)

    # Observations. The paper does not say whether the CD4 cell count is T
    # or T + T*; T is used (T* is about 5 cells/mm3 at baseline and near
    # zero on suppressive treatment, against a 176 cells/mm3 residual SD). The
    # viral load is floored at 1e-4 copies/mL before the log10 transform:
    # under suppression V decays toward 0, and the floor (5.7 log10 below
    # the 50 copies/mL quantification limit) keeps the output finite
    # without changing any prediction above that limit.
    log10_viral_load <- log10(max(virus, 1e-4))
    log10_viral_load ~ add(addSd_log10_viral_load)
    target ~ add(addSd_target)
  })
}
