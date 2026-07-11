Lindauer_2017_lacosamide_dropout <- function() {
  description <- "Repeated time-to-event dropout model (dropout not because of a lack of efficacy) for adult patients newly diagnosed with focal or generalized tonic-clonic epilepsy on lacosamide (LCM) or carbamazepine controlled-release (CBZ-CR) monotherapy, from Lindauer 2017 (SP0993 trial; NCT01243177). The base hazard is a smoothed 4-breakpoint step function h0(t) = k1 + (k2-k1)*S1 + (k3-k2)*S2 + (k4-k3)*S3 + (k5-k4)*S4 where S_i is a logistic sigmoid centred at breakpoint BP_i with steepness GAM=50 (paper Frobel et al. parameterisation). Covariate effects are multiplicative on the hazard: hazard_drop = h0(t) * exp(e_sexf_drop*SEXF + e_conmed_lcm_drop*CONMED_LCM). Sister time-to-seizure model for the same trial: modellib('Lindauer_2017_lacosamide_seizure')."
  reference <- paste(
    "Lindauer A, Laveille C, Stockis A.",
    "Time-to-seizure modeling of lacosamide used in monotherapy in patients",
    "with newly diagnosed epilepsy.",
    "Clin Pharmacokinet. 2017;56(11):1403-1413.",
    "doi:10.1007/s40262-017-0530-8.",
    "Open Access under CC BY-NC 4.0.",
    "Sister model file from the same paper:",
    "modellib('Lindauer_2017_lacosamide_seizure').",
    sep = " "
  )
  vignette <- "Lindauer_2017_lacosamide"
  units <- list(
    time          = "day",
    dosing        = "n/a (no drug-dosing events; treatment arm and demographics enter as time-fixed covariates)",
    concentration = "probability (the model output `sur` is a survival probability, not a drug concentration)"
  )

  covariateData <- list(
    SEXF = list(
      description        = "Subject sex indicator: 1 = female, 0 = male.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Time-fixed per subject. Enters the dropout hazard via multiplicative exponential form exp(e_sexf_drop * SEXF), giving a hazard ratio for female-vs-male of exp(0.238) = 1.27 (90% CI 1.04-1.55, Lindauer 2017 Table 2 row 'HR_SEX').",
      source_name        = "SEX"
    ),
    CONMED_LCM = list(
      description        = "Treatment-arm indicator: 1 = subject assigned to the lacosamide (LCM) monotherapy arm; 0 = subject assigned to the carbamazepine controlled-release (CBZ-CR) monotherapy arm.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (CBZ-CR arm)",
      notes              = "Time-fixed per subject in the SP0993 parallel-group monotherapy design (each subject stays in the arm to which they were randomised). Enters the dropout hazard via multiplicative exponential form exp(e_conmed_lcm_drop * CONMED_LCM), giving a hazard ratio for LCM-vs-CBZ-CR of exp(-0.138) = 0.871 (90% CI 0.714-1.06, Lindauer 2017 Table 2 row 'HR_TYPE'). The paper's source column is `TYPE` (1 = LCM, 0 = CBZ-CR).",
      source_name        = "TYPE"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 883L,
    n_studies      = 1L,
    age_range      = "16-87 years (median 40, IQR 26-55)",
    weight_range   = "not reported in the on-disk trimmed paper text",
    sex_female_pct = 46.3,
    race_ethnicity = NULL,
    disease_state  = "Adult patients (>=16 years) newly diagnosed with focal or generalized tonic-clonic seizures without signs of focal onset, provided they had no history or clinical or electroencephalographic findings suggestive of idiopathic generalized epilepsy (SP0993 inclusion criteria; ClinicalTrials.gov NCT01243177).",
    dose_range     = "LCM target dose levels 200, 400, or 600 mg/day BID (randomisation starting dose 100 mg/day); CBZ-CR target dose levels 400, 800, or 1200 mg/day BID (randomisation starting dose 200 mg/day). Dose escalation to the next level was triggered by seizure occurrence during the 26-week evaluation period at the current dose.",
    regions        = "multinational (SP0993)",
    notes          = "Randomised 883 patients (LCM 443, CBZ-CR 440) analysed for dropout. 205 dropped out for reasons other than lack of efficacy (LCM 111, CBZ-CR 94). Baseline demographics from Lindauer 2017 Table 1. Dropout model was originally developed on the historic N01061 dataset (comparing levetiracetam and CBZ-CR; NCT00150735) and re-estimated on SP0993 as described in Lindauer 2017 Section 2.2 'Modeling Strategy and Software'."
  )

  ini({
    # ----------------------------------------------------------------------
    # Final estimates from Lindauer 2017 Table 2 'Parameter estimates for
    # the dropout model'. Breakpoints are reported as BP1 raw and
    # BP2/BP3/BP4 as increments from the preceding breakpoint, so
    # BP2_abs = BP1 + dBP2, BP3_abs = BP2_abs + dBP3, BP4_abs = BP3_abs + dBP4.
    # Hazards k1-k5 are back-transformed via note (b) 'Parameter back-
    # transformed to normal scale as exp(x)'; on the log scale they enter
    # ini() as log(k_i). Covariate coefficients Coeff_SEX and Coeff_TYPE
    # are on the log-hazard scale (their exp() gives the hazard ratio).
    # ----------------------------------------------------------------------

    bp1_drop  <- log(20);    label("Log breakpoint 1 (day) -- transition from k1 to k2")   # Lindauer 2017 Table 2 BP1 = 20 days (90% CI 18.2-20.9)
    dbp2_drop <- log(9.76);  label("Log delta between BP2 and BP1 (day)")                  # Lindauer 2017 Table 2 D(BP2-BP1) = 9.76 (90% CI 5.29-11.5)
    dbp3_drop <- log(105);   label("Log delta between BP3 and BP2 (day)")                  # Lindauer 2017 Table 2 D(BP3-BP2) = 105 (90% CI 95.7-140)
    dbp4_drop <- log(340);   label("Log delta between BP4 and BP3 (day)")                  # Lindauer 2017 Table 2 D(BP4-BP3) = 340 (90% CI 251-357)

    lk1_drop <- log(0.00205);  label("Log dropout hazard from t=0 to BP1 (1/day)")         # Lindauer 2017 Table 2 k1 = 0.00205 (90% CI 0.0014-0.00267)
    lk2_drop <- log(0.00358);  label("Log dropout hazard from BP1 to BP2 (1/day)")         # Lindauer 2017 Table 2 k2 = 0.00358 (90% CI 0.00254-0.00601)
    lk3_drop <- log(0.000946); label("Log dropout hazard from BP2 to BP3 (1/day)")         # Lindauer 2017 Table 2 k3 = 0.000946 (90% CI 0.000759-0.00123)
    lk4_drop <- log(0.000613); label("Log dropout hazard from BP3 to BP4 (1/day)")         # Lindauer 2017 Table 2 k4 = 0.000613 (90% CI 0.000456-0.00074)
    lk5_drop <- log(0.00221);  label("Log dropout hazard from BP4 to end of trial (1/day)") # Lindauer 2017 Table 2 k5 = 0.00221 (90% CI 0.00134-0.00358)

    e_sexf_drop       <- 0.238;  label("Log-hazard coefficient for female sex (SEXF vs male); HR_SEX = exp(0.238) = 1.27")     # Lindauer 2017 Table 2 Coeff_SEX = 0.238 (90% CI 0.0346-0.437)
    e_conmed_lcm_drop <- -0.138; label("Log-hazard coefficient for lacosamide arm (CONMED_LCM vs CBZ-CR); HR_TYPE = exp(-0.138) = 0.871") # Lindauer 2017 Table 2 Coeff_TYPE = -0.138 (90% CI -0.333-0.0399)

    # ----------------------------------------------------------------------
    # IIV. The dropout model does not include a subject-level random
    # effect (Lindauer 2017 Table 2 lists no IIV row). No eta* parameters
    # are added.
    #
    # Residual error placeholder. The source NONMEM run uses the survival
    # / event-density likelihood (Laplace); no observation-error model.
    # This nlmixr2 translation is intended for forward simulation:
    # `hazard_drop`, `cumhaz_drop`, and `sur` are exposed as derived
    # outputs and a tiny additive residual is attached to `sur` so
    # the nlmixr2 likelihood machinery accepts the model.
    # ----------------------------------------------------------------------
    addSd <- 0.001; label("Placeholder additive residual error on the survival-probability output (unitless); not from the source -- see vignette Assumptions")
  })

  model({
    # Back-transform breakpoints (raw and increments) into absolute days.
    bp1 <- exp(bp1_drop)
    bp2 <- bp1 + exp(dbp2_drop)
    bp3 <- bp2 + exp(dbp3_drop)
    bp4 <- bp3 + exp(dbp4_drop)

    # Back-transform hazard levels.
    k1 <- exp(lk1_drop)
    k2 <- exp(lk2_drop)
    k3 <- exp(lk3_drop)
    k4 <- exp(lk4_drop)
    k5 <- exp(lk5_drop)

    # Frobel-style smoothed step: logistic sigmoids S_i at each breakpoint
    # with steepness gam. Paper Section 2.4 states 'a high value of 50 was
    # chosen for the GAM variable, to mimic a very steep change but
    # avoiding discontinuity'. Sigmoid S_i = 1 / (1 + exp(-gam * (t - bp_i)))
    # transitions from 0 to 1 as t crosses bp_i. The connected form
    # h0(t) = k1 + (k2 - k1)*S1 + (k3 - k2)*S2 + (k4 - k3)*S3 + (k5 - k4)*S4
    # reduces to k_j on the j-th interval (a Riemann-sum expression of the
    # piecewise-constant hazard reported by Table 2).
    gam <- 50
    s1 <- 1.0 / (1.0 + exp(-gam * (t - bp1)))
    s2 <- 1.0 / (1.0 + exp(-gam * (t - bp2)))
    s3 <- 1.0 / (1.0 + exp(-gam * (t - bp3)))
    s4 <- 1.0 / (1.0 + exp(-gam * (t - bp4)))
    h0_drop <- k1 + (k2 - k1) * s1 + (k3 - k2) * s2 + (k4 - k3) * s3 + (k5 - k4) * s4

    # Multiplicative covariate effects (log-linear on hazard).
    hazard_drop <- h0_drop * exp(e_sexf_drop * SEXF + e_conmed_lcm_drop * CONMED_LCM)

    # Cumulative hazard and dropout-survivor probability.
    d/dt(cumhaz_drop) <- hazard_drop
    cumhaz_drop(0)    <- 0
    sur          <- exp(-cumhaz_drop)

    # Forward-simulation observation: placeholder additive residual on the
    # survival-probability output so the nlmixr2 likelihood machinery
    # accepts the model.
    sur ~ add(addSd)
  })
}
