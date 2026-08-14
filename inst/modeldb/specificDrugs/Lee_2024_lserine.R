Lee_2024_lserine <- function() {
  description <- "Two-compartment population PK with zero-order absorption, linked through an effect compartment to a linear drug-effect / linear disease-progression model of the Korean-Vineland Adaptive Behavior Scale-II Adaptive Behavior Composite (K-VABS-II-ABC) score, for AST-001 (L-serine) in pediatric patients with autism spectrum disorder (Lee 2024)"
  reference <- "Lee S, Hwang S-K, Cho J-S, Ryu HC, Chung J-Y. Population pharmacokinetic and pharmacodynamic model guided weight-tiered dose of AST-001 in pediatric patients with autism spectrum disorder. Front Pharmacol. 2024;15:1452526. doi:10.3389/fphar.2024.1452526 -- PK layer fixed from the upstream healthy-adult model: Lee S, Hwang SK, Nam HS, Cho JS, Chung JY. Population pharmacokinetic model of AST-001, L-isomer of serine, combining endogenous production and exogenous administration in healthy subjects. Front Pharmacol. 2022;13:891227. doi:10.3389/fphar.2022.891227"
  vignette <- "Lee_2024_lserine"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    central = list(
      analyte = "L-serine (AST-001, exogenous only)", units = "mg",
      specimen = "plasma", verified = TRUE
    ),
    peripheral1 = list(
      analyte = "L-serine (AST-001, exogenous only)", units = "mg",
      specimen = "plasma", verified = TRUE
    ),
    effect = list(
      analyte = "L-serine effect-site concentration (Ce)", units = "ug/mL",
      specimen = "not applicable", verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Empirical allometric scaling of the fixed upstream adult PK parameters to",
        "pediatric size, normalized to a 70 kg adult. Exponent 0.75 on CL/F and Q/F",
        "and 1 on V1/F and V2/F, both fixed rather than estimated (Lee 2024 'Pharmacokinetic",
        "model'; Lee 2022 Equation 1)."
      ),
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Only covariate retained in the final model. Power function on the baseline",
        "K-VABS-II-ABC score E0, normalized to the 5-year median age of the analysis",
        "population (Lee 2024 Table 1 footnote: 'beta_age, exponent of power relationship",
        "between normalized age (5 years) and E0')."
      ),
      source_name        = "AGE"
    )
  )

  # Screened during covariate selection but NOT retained in the final model
  # (Lee 2024 Methods 'Pharmacodynamic model': "Potential covariates, including
  # age, weight, sex, and baseline total serine level, were evaluated to determine
  # whether they affected PD parameters."; only age survived forward selection /
  # backward elimination). Documented here so the paper's covariate screen is
  # preserved without carrying an unused-covariate convention warning.
  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened as a PD covariate but not retained in the final model. Study population was 120 male (82.8%) / 25 female (17.2%)."
    ),
    SERINE_BL = list(
      description = "Baseline total serum serine concentration",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Screened as a PD covariate but not retained in the final model. Median (min-max) 128 (74.8-609) umol/L."
    )
  )

  population <- list(
    n_subjects    = 145,
    n_observations = 570,
    n_studies     = 1,
    species       = "human",
    age_range     = "2-11 years",
    weight_range  = "10.7-58.1 kg (median 20.5)",
    sex_female_pct = 17.2,
    race_ethnicity = "Korean",
    disease_state = "Autism spectrum disorder (pediatric)",
    dose_range    = paste(
      "Weight-tiered fixed doses. High-dose arm 2, 4, 7, 10 and 14 g twice daily for",
      "the 10-14, 15-24, 25-37, 38-51 and 52-60 kg weight bands respectively",
      "(~400 mg/kg/day); low-dose arm half of those (~200 mg/kg/day). 12-week primary",
      "period plus a 12-week extension; follow-up at 36 weeks."
    ),
    regions       = "Republic of Korea",
    notes         = paste(
      "Phase II multi-center, randomized, double-blind, placebo-controlled trial",
      "(KCT0007519); 49 placebo / 50 low-dose / 46 high-dose. NO pediatric PK samples",
      "were collected -- the PK layer is fixed from the upstream healthy-adult model",
      "(Lee 2022, 24 healthy Korean male volunteers, 648 concentrations) and",
      "extrapolated by empirical allometry. Only the PD layer was estimated here."
    )
  )

  ini({
    # ==================================================================
    # PK layer -- ALL parameters FIXED from the upstream healthy-adult
    # model, Lee 2022 Table 1 (final estimates). Lee 2024 'Pharmacokinetic
    # model': "All PK parameters were fixed using final model estimates."
    #
    # Time-unit conversion: Lee 2022 reports CL/F and Q/F in L/h and D1 in
    # h; this model's time unit is days (the PD parameters Ke0 and Kprog
    # are per day and the trial runs 12-36 weeks), so the flow terms are
    # multiplied by 24 h/day and D1 divided by 24.
    #
    # The upstream model ALSO carried a zero-order endogenous L-serine
    # production term R1 = 0.287 g/h. It is deliberately ABSENT here:
    # Lee 2024 'Pharmacokinetic model' states "We used the final PK model
    # that excluded endogenous production because endogenous L-serine
    # levels interfered with parameter estimation related to drug effects
    # by externally administered AST-001." Cc is therefore the
    # exogenous-only (baseline-adjusted) L-serine concentration.
    # ==================================================================
    lcl <- fixed(log(22.9 * 24)); label("Apparent clearance CL/F at 70 kg (L/day)")                      # Lee 2022 Table 1: CL/F = 22.9 L/h (RSE 5.8%), x24 h/day = 549.6 L/day
    lvc <- fixed(log(68.4));      label("Apparent central volume of distribution V1/F at 70 kg (L)")     # Lee 2022 Table 1: V1/F = 68.4 L (RSE 6.7%)
    lq  <- fixed(log(16.4 * 24)); label("Apparent inter-compartmental clearance Q/F at 70 kg (L/day)")   # Lee 2022 Table 1: Q/F = 16.4 L/h (RSE 7.3%), x24 h/day = 393.6 L/day
    lvp <- fixed(log(196));       label("Apparent peripheral volume of distribution V2/F at 70 kg (L)")  # Lee 2022 Table 1: V2/F = 196 L (RSE 11.6%)
    ld1 <- fixed(log(1.26 / 24)); label("Duration of zero-order absorption D1 (day)")                    # Lee 2022 Table 1: D1 = 1.26 h (RSE 7.6%), /24 h/day = 0.0525 day

    # Allometric exponents -- fixed at the empirical values, not estimated.
    # Lee 2024 'Pharmacokinetic model': "By applying the empirical allometry
    # coefficient of 0.75 and 1 in clearance and volume of distribution,
    # respectively, we assumed the PK profile of pediatrics." Lee 2022 had
    # already rejected its own estimated exponents as implausible and fixed
    # them for its pediatric simulation.
    e_wt_cl <- fixed(0.75); label("Allometric exponent of body weight on CL/F and Q/F")        # Lee 2024 'Pharmacokinetic model'; Lee 2022 'Model-Based Simulation' Equation 1
    e_wt_vc <- fixed(1);    label("Allometric exponent of body weight on V1/F and V2/F")      # Lee 2024 'Pharmacokinetic model'; Lee 2022 'Model-Based Simulation' Equation 1

    # PK inter-individual variability -- Lee 2022 Table 1, reported as CV%
    # (table footnote b: "Inter-individual variability is presented as
    # coefficient of variation (%)"), converted to log-scale variance via
    # omega^2 = log(1 + CV^2). The V1/F-CL/F covariance is reported
    # untransformed (footnote c) and is used as-is. Fixed, like the rest of
    # the PK layer.
    etalcl + etalvc ~ fixed(log(1 + 0.222^2),
                            0.0483, log(1 + 0.290^2))  # Lee 2022 Table 1: omega_CL/F = 22.2% CV (RSE 23.9%); omega_V1/F = 29.0% CV (RSE 13.1%); cov(V1/F, CL/F) = 0.0483 (RSE 41.8%)
    etald1 ~ fixed(log(1 + 0.380^2))                   # Lee 2022 Table 1: omega_D1 = 38.0% CV (RSE 13.6%)

    # ==================================================================
    # PD layer -- Lee 2024 Table 1 (final parameter estimates of the
    # population pharmacodynamic model of K-VABS-II-ABC score). Estimated
    # in Monolix, so the reported Omega values are random-effect standard
    # deviations, not variances.
    #
    # Structure, from Lee 2024 Figure 1 (which prints the equations):
    #   K-VABS-II-ABC score = E0 + Drug effect + Natural progression
    #   Drug effect         = Deff x Ce
    #   Natural progression = Kprog x time
    # ==================================================================
    lke0 <- log(0.0065); label("Effect-site first-order equilibration rate constant Ke0 (1/day)")  # Lee 2024 Table 1: Ke0 = 0.0065 /day (RSE 7.92%); consistent with the stated ~15-week equilibration half-life (log(2)/0.0065 = 107 days = 15.2 weeks)
    lrbase <- log(48.51); label("Baseline K-VABS-II-ABC score E0 at the reference age of 5 years") # Lee 2024 Table 1: E0 = 48.51 (RSE 1.66%)
    e_age_rbase <- -0.21; label("Power exponent of age (normalized to 5 years) on the baseline K-VABS-II-ABC score E0")  # Lee 2024 Table 1: beta_age = -0.21 (RSE 18.4%); Table 1 footnote defines it as the exponent on age normalized to 5 years

    # ------------------------------------------------------------------
    # UNRESOLVED (operator sidecar oare_PMC11682956 request-002). Deff is
    # the only parameter carrying the drug effect, and the paper's own
    # numbers cannot all be true at once. Three mutually inconsistent
    # statements of the same quantity appear in Lee 2024:
    #
    #   (i)   Table 1:      Deff = 0.0022 "L/ug".  Dimensionally, a
    #         dimensionless score increment = [L/ug] x [Ce] requires Ce in
    #         ug/L -- i.e. 2.2 score per (ug/mL), NOT 0.0022 per (ug/mL).
    #   (ii)  Discussion:   subgroup median Deff "1.45 (mL/ug)" and
    #         "1.65 (mL/ug)". mL/ug x ug/mL is dimensionless, so this is
    #         score per (ug/mL), and it corroborates (i): 0.0022 L/ug
    #         = 2.2 mL/ug, with the two subgroup medians just below the
    #         typical value (as expected for a log-normal Deff).
    #   (iii) Table 2 / Figure 2: the paper's own simulated 12-week
    #         delta-scores (median 1.67 -> 3.74 across the five weight
    #         bands). Reproducing these requires Deff = 0.017 - 0.080
    #         score per (ug/mL) -- about 1/30 of (i) and (ii), and not a
    #         single constant: the required value rises as WT^0.895
    #         across the bands (4.79x, versus a 2.15x rise in Ce), so NO
    #         constant Deff in ANY unit reproduces Table 2.
    #
    # Verified by simulation (see the vignette Errata and the task report):
    # 2.2 score/(ug/mL) gives a median 12-week gain of 27-58 points on a
    # 48-point baseline (target attainment 98-99% vs published 40.5-73%);
    # 0.0022 score/(ug/mL) gives 0.03-0.06 points, i.e. indistinguishable
    # from placebo (target attainment 32-35%). The PK layer and the
    # placebo/progression layer both reproduce the paper independently, so
    # the inconsistency is localized to this one parameter.
    #
    # Pending the ruling, the line below encodes the printed Table 1 value
    # literally against Ce in ug/mL. This is the conservative placeholder,
    # NOT a resolution -- do not treat the drug-effect magnitude of this
    # model as validated until request-002 is answered.
    # ------------------------------------------------------------------
    lslope <- log(0.0022); label("Linear drug-effect slope Deff on effect-site L-serine concentration (K-VABS-II-ABC score per ug/mL)")  # Lee 2024 Table 1: Deff = 0.0022 L/ug (RSE 37.5%); see the UNRESOLVED note above -- Table 1, the Discussion's mL/ug figures, and Table 2 are mutually irreconcilable

    # Kprog is NOT log-transformed. Lee 2024 Results: "The slope of the
    # natural K-VABS-II-ABC score progression (Kprog) was estimated to be
    # 0.015 day-1 with a standard deviation of random effects of 0.018 in a
    # normal eta distribution." A normal eta of SD 0.018 around 0.015 puts
    # ~20% of the population at a negative progression slope (declining
    # adaptive behavior), which a log-normal parameterization cannot
    # represent, so the parameter and its eta are both untransformed.
    dp_slope <- 0.015; label("Linear natural disease-progression slope Kprog on the K-VABS-II-ABC score (score/day)")  # Lee 2024 Table 1: Kprog = 0.015 /day (RSE 12.0%)

    # PD inter-individual variability -- Lee 2024 Table 1 'Random effects'.
    # Monolix reports Omega as the random-effect SD, so the variances below
    # are the squares. E0 and Kprog are correlated (Table 1 'Correlation'
    # row: Correlation E0-Kprog = 0.48, RSE 16.9%):
    #   cov = 0.48 * 0.2 * 0.018 = 0.001728
    etalrbase + etadp_slope ~ c(0.2^2,
                                0.48 * 0.2 * 0.018, 0.018^2)  # Lee 2024 Table 1: Omega_E0 = 0.2 (RSE 6.01%); Omega_Kprog = 0.018 (RSE 8.27%); correlation E0-Kprog = 0.48
    etalke0 ~ 0.21^2      # Lee 2024 Table 1: Omega_Ke0 = 0.21 (RSE 22.3%)
    etalslope ~ 1.4^2     # Lee 2024 Table 1: Omega_Deff = 1.4 (RSE 20.8%)

    # Residual error. Lee 2024 Methods 'Pharmacodynamic model': "The additive
    # error model was used to account for residual unexplained variability in
    # the K-VABS-II-ABC scores."
    addSd <- 1.6; label("Additive residual error on the K-VABS-II-ABC score")  # Lee 2024 Table 1: additive error = 1.6 (RSE 4.56%)

    # Residual error on the PK observation. No pediatric PK samples were
    # collected in Lee 2024, so this comes from the upstream adult model and
    # is fixed with the rest of the PK layer.
    propSd <- fixed(0.183); label("Proportional residual error on plasma L-serine Cc (fraction)")  # Lee 2022 Table 1: proportional residual error = 0.183 (RSE 7.5%), reported as a coefficient of variation (footnote d)
  })

  model({
    # ---- PK layer: allometrically scaled to pediatric size, 70 kg reference
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    q  <- exp(lq)           * (WT / 70)^e_wt_cl
    vp <- exp(lvp)          * (WT / 70)^e_wt_vc
    d1 <- exp(ld1 + etald1)

    # ---- PD individual parameters
    ke0    <- exp(lke0 + etalke0)
    rbase  <- exp(lrbase + etalrbase) * (AGE / 5)^e_age_rbase
    slope  <- exp(lslope + etalslope)
    kprog  <- dp_slope + etadp_slope

    # ---- Two-compartment disposition with zero-order absorption directly
    # into the central compartment (Lee 2022 Figure 1: no depot; the dose is
    # delivered as a zero-order input of duration D1). Endogenous production
    # is excluded, so `central` carries exogenous AST-001 only.
    d/dt(central)     <- -(cl / vc) * central - (q / vc) * central + (q / vp) * peripheral1
    d/dt(peripheral1) <-  (q / vc) * central - (q / vp) * peripheral1
    dur(central) <- d1

    Cc <- central / vc

    # ---- Effect compartment (Lee 2024 Figure 1: Ce equilibrates with the
    # central-compartment concentration at rate Ke0). Starts empty: the
    # exogenous-only concentration is zero before the first dose.
    d/dt(effect) <- ke0 * (Cc - effect)
    effect(0)    <- 0

    # ---- K-VABS-II-ABC score (Lee 2024 Figure 1, verbatim):
    #   score = E0 + Drug effect + Natural progression
    #   Drug effect = Deff x Ce;  Natural progression = Kprog x time
    # E0 cancels out of the change-from-baseline endpoint the paper reports,
    # so the simulated delta-score depends only on `slope`, `effect` and `kprog`.
    KVABS <- rbase + slope * effect + kprog * time

    Cc    ~ prop(propSd)
    KVABS ~ add(addSd)
  })
}
