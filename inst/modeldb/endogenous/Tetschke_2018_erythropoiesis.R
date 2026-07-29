Tetschke_2018_erythropoiesis <- function() {
  description <- "Three-compartment population mixed-effects model for human erythropoiesis (red blood cell regeneration after a phlebotomy / blood donation) in healthy adults"
  reference   <- "Tetschke M, Lilienthal P, Pottgiesser T, Fischer T, Schalk E, Sager S. Mathematical Modeling of Red Blood Cell Count Dynamics after Blood Loss. Processes. 2018;6(9):157. doi:10.3390/pr6090157"
  vignette    <- "Tetschke_2018_erythropoiesis"

  covariateData <- list(
    THB_MASS = list(
      description        = "Subject baseline (steady-state) total hemoglobin mass in grams. Plasma-volume-independent measurement obtained by the optimised CO-rebreathing method (Schmidt 2005, Pottgiesser 2008); paper symbol Base.",
      units              = "g",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Subject-level steady-state baseline. Tetschke 2018 estimates per-subject Base from the arithmetic mean of pre-donation tHb measurements (Section 3.3); the population in Pottgiesser 2008 (29 healthy adult males, 30 +/- 10 years, 76.6 +/- 11.2 kg) had Base ~870 g (Tetschke 2018 Table A2 mean across 28 estimable subjects). Reference value 885.42 g (Tetschke 2018 Table 1 example subject). The model uses THB_MASS to (i) derive the constant inflow term X0 = alpha * THB_MASS, (ii) define the steady-state initial conditions of all three states (Eq. 3), and (iii) set the negative-feedback set point in Eq. 2.",
      source_name        = "Base"
    )
  )

  population <- list(
    n_subjects     = 29,
    n_studies      = 1,
    age_range      = "30 +/- 10 years (mean +/- SD), all adult males",
    weight_range   = "76.6 +/- 11.2 kg (mean +/- SD)",
    sex_female_pct = 0,
    disease_state  = "Healthy adult male volunteers undergoing a single 1-unit standard erythrocyte-concentrate blood donation",
    dose_range     = "Single 1-unit blood donation; per-subject hemoglobin loss inferred from the difference between pre-donation Base and the minimum tHb value in the eight days following donation (Tetschke 2018 Section 3.3; Appendix A.1)",
    regions        = "Germany (Pottgiesser et al. 2008 Transfusion 48:1390 dataset)",
    notes          = "Population NLME estimation in NONMEM 7.4 FOCE-I across 276 tHb observations from 29 subjects (Tetschke 2018 Section 4.2). Diagonal OMEGA, exponential (log-normal) IIV on beta and gamma, additive residual error model; the additive residual SD is not reported in the manuscript so the packaged model omits residual error and is intended for typical-value + IIV simulation rather than refitting. Height 181 +/- 7 cm. Total hemoglobin (tHb) was measured by the optimised CO-rebreathing method (Schmidt 2005)."
  )

  ini({
    lbeta  <- log(1.02); label("Cell maturation amplification factor (beta, dimensionless)")  # Tetschke 2018 Section 4.2 (population fixed effect, SE 0.151)
    lgamma <- log(0.46); label("EPO feedback amplification factor (gamma, dimensionless)")    # Tetschke 2018 Section 4.2 (population fixed effect, SE 0.0651)

    etalbeta  ~ 0.294   # Tetschke 2018 Section 4.2 (diagonal OMEGA, log-normal IIV variance, SE 0.125)
    etalgamma ~ 0.346   # Tetschke 2018 Section 4.2 (diagonal OMEGA, log-normal IIV variance, SE 0.148)
  })

  model({
    # Mechanism constants (Tetschke 2018 Table 1; Assumption 12, Section 2.2):
    # k1 and k2 are reciprocals of the 8-day EPO-proliferating and 6-day non-EPO-proliferating
    # phase durations; alpha is the reciprocal of the 120-day mature-erythrocyte lifespan.
    k1    <- 1 / 8     # 1/day, transit out of EPO-proliferating compartment
    k2    <- 1 / 6     # 1/day, transit out of non-EPO-proliferating compartment
    alpha <- 1 / 120   # 1/day, mature-erythrocyte apoptosis rate

    # Individual amplification factors with log-normal IIV
    beta  <- exp(lbeta + etalbeta)
    gamma <- exp(lgamma + etalgamma)

    # Constant inflow X0 derived from the steady-state condition X0 = k1 * x1_bar = alpha * Base
    # (Tetschke 2018 Eq. 3 with x1_bar = (alpha/k1) * Base).
    X0 <- alpha * THB_MASS

    # Negative feedback function on the EPO-proliferating compartment (Tetschke 2018 Eq. 2)
    Fb <- gamma * (THB_MASS - thb) / THB_MASS

    # ODE system (Tetschke 2018 Eq. 1):
    # precursor1 -- EPO-proliferating progenitors (paper x1, dimensionless cell-count units)
    # precursor2 -- non-EPO-proliferating progenitors (paper x2, dimensionless)
    # thb        -- mature-erythrocyte total hemoglobin mass (paper x3, grams)
    d/dt(precursor1) <- beta * (X0 - k1 * precursor1) + Fb * precursor1
    d/dt(precursor2) <- beta * (k1 * precursor1 - k2 * precursor2)
    d/dt(thb)        <- beta * (k2 * precursor2 - alpha * thb)

    # Steady-state initial conditions (Tetschke 2018 Eq. 3)
    precursor1(0) <- alpha / k1 * THB_MASS
    precursor2(0) <- alpha / k2 * THB_MASS
    thb(0)        <- THB_MASS
  })
}
