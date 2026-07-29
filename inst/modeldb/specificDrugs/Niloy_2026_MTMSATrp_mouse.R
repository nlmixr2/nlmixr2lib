Niloy_2026_MTMSATrp_mouse <- function() {
  description <- paste(
    "Preclinical (mouse). One-compartment population PK model for",
    "MTMSA-Trp, a novel mithramycin analogue investigated for Ewing",
    "sarcoma, in female athymic nu/nu mice following single IV bolus",
    "doses of 0.3, 1, 3, 5, or 10 mg/kg. First-order elimination",
    "from the central compartment with an empirical power-function",
    "effect of dose on apparent clearance (CL decreases with",
    "increasing dose; reference dose 3 mg/kg, exponent beta = -0.30).",
    "Parameters are expressed in per-kg body-weight units (mL/h/kg",
    "for CL, mL/kg for V) so the dose record carries the per-kg dose",
    "directly (mg/kg) without an explicit body-weight covariate.",
    "Parameter values from Niloy 2026 Table 1 (final model)."
  )
  reference <- paste(
    "Niloy KK, Horn J, Bhuiyan NH, Shaaban KA, Bhosale SS,",
    "Prisinzano T, Thorson JS, Rohr J, Leggas M. (2026).",
    "Nonlinear Mixed-Effects Modeling to Characterize the",
    "Pharmacokinetics of a Novel Mithramycin Analogue for Ewing",
    "Sarcoma in Mice. Research Square preprint.",
    "doi:10.21203/rs.3.rs-9035594/v1.",
    sep = " "
  )
  vignette <- "Niloy_2026_MTMSATrp_mouse"
  units    <- list(time = "h", dosing = "mg/kg", concentration = "ng/mL")

  covariateData <- list(
    COHDOSE = list(
      description        = "Per-subject IV bolus dose cohort, expressed in mg/kg",
      units              = "mg/kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power scaling on apparent clearance with reference value",
        "3 mg/kg and exponent beta = -0.30 (Niloy 2026 Methods,",
        "Model development, Eq. 2 and Table 1). Each mouse received",
        "a single IV bolus at one of 0.3, 1, 3, 5, or 10 mg/kg, so",
        "this is a subject-level (time-fixed) covariate. The dose-CL",
        "effect is empirical (Discussion): it should be interpreted",
        "as a description of non-proportional exposure over 0.3-10",
        "mg/kg, not as a mechanistic identification of saturable",
        "elimination. A Michaelis-Menten elimination model was",
        "explored but did not yield stable estimates and was not",
        "retained (Methods, Model development)."
      ),
      source_name        = "Dose_i"
    )
  )

  population <- list(
    species         = "mouse (female athymic nu/nu, Envigo RMS)",
    n_subjects      = 70L,
    n_observations  = 121L,
    n_studies       = 1L,
    age_range       = "(not tabulated; adult athymic nude mice acclimated >= 1 week prior to dosing)",
    weight_range    = "(not tabulated; parameters expressed in per-kg units)",
    sex_female_pct  = 100,
    disease_state   = "Healthy female athymic nu/nu mice (PK dataset; xenograft-naive)",
    dose_range      = "Single IV bolus at 0.3, 1, 3, 5, or 10 mg/kg",
    regions         = "United States (St. Jude Children's Research Hospital, Memphis TN; University of Kentucky, Lexington KY)",
    notes           = paste(
      "Athymic nu/nu mice from Envigo RMS housed under AAALAC",
      "guidelines (12 h light-dark cycle, ad libitum food and",
      "water). Plasma quantified by validated LC-MS/MS (Niloy 2026",
      "Methods reference [8]). Sparse / destructive sampling",
      "design: some animals contributed a single concentration.",
      "Shrinkage on CL was high (81%) and on V moderate (39.1%);",
      "the source emphasizes simulation-based evaluation (VPC,",
      "NPDE) over EBE-based diagnostics. Niloy 2026 Methods",
      "(Mice, Pharmacokinetic data) and Results."
    )
  )

  ini({
    # ------------------------------------------------------------
    # Structural PK parameters -- Niloy 2026 Table 1 (final model).
    # Reported in per-kg body-weight units so the dose carried in
    # the event record (mg/kg) feeds the model without an explicit
    # body-weight covariate. Reference dose = 3 mg/kg (Methods,
    # Model development, Eq. 2).
    # ------------------------------------------------------------
    lcl <- log(39.18); label("Apparent clearance CL at the reference dose of 3 mg/kg (mL/h/kg)")  # Niloy 2026 Table 1: CL = 39.18 (RSE 4.97%)
    lvc <- log(53.06); label("Apparent volume of distribution V (mL/kg)")                          # Niloy 2026 Table 1: V  = 53.06 (RSE 7.55%)

    # ------------------------------------------------------------
    # Dose-as-covariate-on-CL power coefficient. Niloy 2026 Methods
    # Eq. 2:
    #     CL_i = CL_pop * (Dose_i / Dose_median)^beta
    # with Dose_median = 3 mg/kg. Negative beta means CL decreases
    # with increasing dose (greater-than-proportional AUC).
    # ------------------------------------------------------------
    e_cohdose_cl <- -0.30; label("Dose (mg/kg) exponent on CL (unitless)")  # Niloy 2026 Table 1: beta = -0.30 (RSE 10.9%)

    # ------------------------------------------------------------
    # Inter-individual variability. Niloy 2026 Methods: lognormal
    # (exponential) IIV on CL and V. The %CV column of Table 1 is
    # converted to NONMEM-style omega^2 via the Monolix log-normal
    # convention:
    #     omega^2 = log(1 + CV^2)
    #   omega^2_CL = log(1 + 0.058^2) = 0.003359
    #   omega^2_V  = log(1 + 0.25^2)  = 0.060625
    # ------------------------------------------------------------
    etalcl ~ 0.003359  # Niloy 2026 Table 1: IIV on CL = 5.8% CV (RSE 42%); shrinkage 81%
    etalvc ~ 0.060625  # Niloy 2026 Table 1: IIV on V  = 25%  CV (RSE 31%); shrinkage 39.1%

    # ------------------------------------------------------------
    # Residual error. Niloy 2026 Methods / Results: proportional
    # residual error model. Monolix reports the proportional
    # residual coefficient b directly. Both the estimate (0.51)
    # and bootstrap median (0.5) are decimal fractions (proportional
    # SD = 51% CV), as is conventional for Monolix output; the
    # bootstrap 95% CI 0.41-0.58 confirms the decimal interpretation
    # (the "(%)" tag in the Table 1 row label appears to be carried
    # over from the IIV rows and does not apply literally to b).
    # ------------------------------------------------------------
    propSd <- 0.51; label("Proportional residual error (fraction)")  # Niloy 2026 Table 1: b = 0.51 (RSE 8.13%)
  })

  model({
    # ------------------------------------------------------------
    # 1. Individual PK parameters with dose-power covariate on CL
    #    (Niloy 2026 Eq. 2). Reference dose 3 mg/kg.
    # ------------------------------------------------------------
    cl <- exp(lcl + etalcl) * (COHDOSE / 3)^e_cohdose_cl
    vc <- exp(lvc + etalvc)

    # ------------------------------------------------------------
    # 2. Micro-constant
    # ------------------------------------------------------------
    kel <- cl / vc

    # ------------------------------------------------------------
    # 3. One-compartment IV bolus. Dose lands in central directly.
    # ------------------------------------------------------------
    d/dt(central) <- -kel * central

    # ------------------------------------------------------------
    # 4. Observation. Per-kg parameterisation: dose carried in mg/kg
    #    on the event record, vc in mL/kg, so central/vc is mg/mL;
    #    multiply by 1e6 to express plasma concentration in ng/mL
    #    (typical LC-MS/MS reporting units).
    # ------------------------------------------------------------
    Cc <- central / vc * 1e6

    Cc ~ prop(propSd)
  })
}
