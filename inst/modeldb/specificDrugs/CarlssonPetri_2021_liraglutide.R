CarlssonPetri_2021_liraglutide <- function() {
  description <- "Liraglutide PK model in adolescents (Carlsson Petri 2021)"
  reference <- "Carlsson Petri KC, Hale PM, Hesse D, Rathor N, Mastrandrea LD. Liraglutide pharmacokinetics and exposure-response in adolescents with obesity. Pediatric Obesity. 2021;16(10):e12799. doi:10.1111/ijpo.12799"
  units <- list(time = "hr", dosing = "mg", concentration = "nmol/L") # Carlsson Petri 2021 Table 3 / Methods; LLOQ 0.03 nmol/L

  covariateData <- list(
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on CL/F and V/F; reference weight 100 kg per Carlsson Petri 2021 Table 3 footnote (reference: female, 100 kg, adult).",
      source_name        = "WT"
    ),
    CHILD = list(
      description        = "Indicator for child age group (7-11 years)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (adult, >=18 years)",
      notes              = "Paired with ADOLESCENT; both 0 indicates adult (>=18 y) reference. Age cutoffs per Carlsson Petri 2021 trial design: children 7-11 y (NN2211-4181), adolescents 12-17 y (NN2211-4180 / NN2211-3967), adults >=18 y (NN2211-3630).",
      source_name        = "CHILD"
    ),
    ADOLESCENT = list(
      description        = "Indicator for adolescent age group (12-17 years)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (adult, >=18 years)",
      notes              = "Paired with CHILD; both 0 indicates adult reference. Age cutoffs per Carlsson Petri 2021 trial design: adolescents 12-17 y.",
      source_name        = "ADOLESCENT"
    ),
    SEXF = list(
      description        = "Biological sex indicator, 1 = female, 0 = male",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Carlsson Petri 2021 reports a male/female CL/F ratio of 1.12 (Table 3). Implemented as e_sex_cl^(1 - SEXF), which evaluates to 1.12 for males (SEXF=0) and 1.00 for females (SEXF=1), matching the paper. Source table footnote identifies female as the reference category.",
      source_name        = "SEXM"
    )
  )

  population <- list(
    n_subjects     = 176L, # 121 + 13 + 13 + 29 (Carlsson Petri 2021 Table 1; pooled analysis across 4 trials)
    n_studies      = 4L,   # NN2211-4180, NN2211-3967, NN2211-4181, NN2211-3630 (Carlsson Petri 2021 Methods)
    age_range      = "7-64 years (children 7-11, adolescents 12-17, adults >=18)",
    weight_range   = "Adolescents (trial 4180): 62.1-178.2 kg (mean 99.4, SD 19.7); children, adolescents (3967), adults pooled with different distributions",
    sex_female_pct = 55.4, # Carlsson Petri 2021 Table 1, trial NN2211-4180 adolescent phase 3a (67/121)
    race_ethnicity = c(White = 84.3, Black = 10.7, Asian = 1.7, Other = 3.3), # Carlsson Petri 2021 Table 1 (trial NN2211-4180)
    disease_state  = "Adolescents with obesity (primary population); pooled with children and adults for PK bridging.",
    dose_range     = "0.6 to 3.0 mg once daily SC (weekly escalation 0.6, 1.2, 1.8, 2.4, 3.0). Maintenance 3.0 mg in 107/121 adolescents (Carlsson Petri 2021 Table 1).",
    regions        = "Multi-national (trials listed on ClinicalTrials.gov: NN2211-4180, NN2211-3967, NN2211-4181, NN2211-3630).",
    trials         = c("NN2211-4180", "NN2211-3967", "NN2211-4181", "NN2211-3630"),
    tanner_stage   = "Stages 2-5 (prepubertal excluded) for adolescents (Carlsson Petri 2021 Methods).",
    notes          = "Pooled pop-PK bridging analysis: 121 adolescents (NN2211-4180), 13 adolescents (NN2211-3967), 13 children (NN2211-4181), 29 adults (NN2211-3630)."
  )

  ini({
    lka <- fixed(log(0.0813)); label("Absorption rate (1/hr)") # Carlsson Petri 2021 Table 3 (fixed)
    lcl <- log(1.01); label("Apparent clearance (L/h)") # Carlsson Petri 2021 Table 3 (95% CI 0.922-1.09)
    e_wt_cl <- 0.762; label("Body weight exponent on CL/F") # Carlsson Petri 2021 Table 3 (95% CI 0.565-0.958), reference 100 kg
    e_sex_cl <- 1.12; label("Male/female CL/F ratio") # Carlsson Petri 2021 Table 3 (95% CI 0.993-1.24); applied as ratio^(1 - SEXF)
    e_age_child_cl <- 1.11; label("Child/adult CL/F ratio") # Carlsson Petri 2021 Table 3 (90% CI 0.89-1.34); applied as ratio^CHILD
    e_age_adolescent_cl <- 1.06; label("Adolescent/adult CL/F ratio") # Carlsson Petri 2021 Table 3 (90% CI 0.931-1.19); applied as ratio^ADOLESCENT
    lvc <- fixed(log(13.8)); label("Apparent central volume of distribution (L)") # Carlsson Petri 2021 Table 3 (fixed)
    e_wt_vc <- 0.587; label("Body weight exponent on V/F") # Carlsson Petri 2021 Table 3 (95% CI 0.475-0.700), reference 100 kg

    # IIV reported as CV% in Carlsson Petri 2021 Table 3. Table header/footnote defines
    # CV% = sqrt(exp(omega^2) - 1) * 100 (log-normal convention), so omega^2 = log(1 + CV^2).
    etalcl ~ log(1 + 0.312^2) # Carlsson Petri 2021 Table 3: IIV CL/F = 31.2% CV -> omega^2 = log(1 + 0.312^2)
    etalvc ~ log(1 + 0.317^2) # Carlsson Petri 2021 Table 3: IIV V/F  = 31.7% CV -> omega^2 = log(1 + 0.317^2)
    propSd <- 0.433; label("Proportional residual error (fraction)") # Carlsson Petri 2021 Table 3: 43.3%
  })
  model({
    ka <- exp(lka)
    cl_wt <- (WT / 100)^e_wt_cl # Carlsson Petri 2021 Table 3: (WT/100 kg)^0.762
    cl_sex <- e_sex_cl^(1 - SEXF) # Carlsson Petri 2021: male/female ratio 1.12 applied as 1.12^SEXM; SEXM = 1 - SEXF
    cl_age <- e_age_child_cl^CHILD * e_age_adolescent_cl^ADOLESCENT # Carlsson Petri 2021 Table 3: 1.11^CHILD * 1.06^ADOLESCENT (adult reference)
    cl <- exp(lcl + etalcl) * cl_wt * cl_sex * cl_age # Carlsson Petri 2021: CL/F_i = TVCL * E(WT) * E(sex) * E(age) * exp(eta_CL)
    vc_wt <- (WT / 100)^e_wt_vc # Carlsson Petri 2021 Table 3: (WT/100 kg)^0.587 on V/F (same 100 kg reference as CL)
    vc <- exp(lvc + etalvc) * vc_wt # Carlsson Petri 2021: V/F_i = TVV * E(WT) * exp(eta_V)

    Cc <- linCmt()
    Cc ~ prop(propSd)
  })
}
