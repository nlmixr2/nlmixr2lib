Gong_2025_talquetamab <- function() {
  description <- paste(
    "Two-compartment population pharmacokinetic model for talquetamab, a",
    "GPRC5D x CD3 bispecific antibody, in patients with relapsed/refractory",
    "multiple myeloma (Gong 2025, MonumenTAL-1 phase I/II; 5354 serum",
    "concentrations from 492 patients dosed intravenously or subcutaneously).",
    "Subcutaneous absorption is sequential zero-order (duration D1) into a",
    "depot then first-order (Ka) into the central compartment. Total",
    "clearance is CL(t) = CL0 * (CLPCT + (1 - CLPCT) * exp(-KDES * t)): a",
    "time-independent share CL0 * CLPCT plus a time-dependent share that",
    "decays from the start of treatment with first-order rate KDES. Body",
    "weight acts on CL0 and the central volume; non-IgG myeloma and ISS",
    "stage II or III act on CL0. The paper's exposure-response analyses",
    "for response, survival and adverse events are graphical or",
    "logistic-regression summaries without a reusable model, so only the",
    "population PK model is encoded."
  )
  reference <- paste(
    "Gong J, Zhou J, Yuan D, Ma X, Vishwamitra D, Hilder B, Masterson TJ,",
    "Tolbert J, Renaud T, Heuck C, Kane C, Samtani MN, Girgis S,",
    "Haddish-Berhane N, Berdeja J, Krishnan A, Ouellet D.",
    "Recommended phase II doses of talquetamab in patients with",
    "relapsed/refractory multiple myeloma from MonumenTAL-1: clinical",
    "pharmacology results. Clin Pharmacol Ther. 2025;118(4):954-966.",
    "doi:10.1002/cpt.70004"
  )
  vignette <- "Gong_2025_talquetamab"

  # Gong 2025 Table S5 reports clearances in L/day, volumes in L, Ka and KDES
  # in 1/day and D1 in days; serum concentrations are reported in ng/mL
  # throughout (Table S3, Figures 2-5). Doses are in mg (the mg/kg regimens
  # are converted to mg per patient). central / vc is mg/L = ug/mL, so the
  # observation is multiplied by 1000 to give ng/mL.
  units <- list(time = "day", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(analyte = "talquetamab", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "talquetamab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "talquetamab", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect normalised to 74.8 kg, the median body weight of the population PK dataset (Table S4, overall column; range 40.2-143 kg), with exponent 0.672 on CL0 and 0.670 on V1 (Table S5 and footnotes a, b). Body weight was included a priori on CL0 and V1 because it is an established covariate for IgG-based antibodies (Methods, Pharmacokinetics).",
      source_name = "BWT"
    ),
    MM_NIGG = list(
      description = "Non-IgG-secreting multiple myeloma indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (IgG multiple myeloma) -- Gong 2025 anchors the typical CL0 to IgG myeloma and applies a proportional reduction to non-IgG patients, the same orientation as Miao_2023_teclistamab. Canonical column semantics (1 = non-IgG) are preserved.",
      notes = "Proportional effect -0.547 on CL0 (Table S5 row 'IgG or non-IgG on CL0'), i.e. a factor 1 - 0.547 = 0.453 for non-IgG myeloma, which is the 0.453^(TPMMG = Non-IgG) printed in Table S5 footnote a. Non-IgG patients therefore have lower clearance and about two-fold higher Cavg,4weeks (Results; Figure 4). Cohort split: IgG 247 (50.2%), non-IgG 245 (49.8%) (Table S4).",
      source_name = "TPMMG"
    ),
    ISS_II = list(
      description = "Baseline International Staging System stage II multiple myeloma indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (ISS stage I is the reference level; ISS_II and ISS_III are mutually exclusive)",
      notes = "Gong 2025 estimates ONE proportional effect, +0.318 on CL0, shared by ISS stage II and stage III (Table S5 row 'ISS II or III on CL0'; footnote a prints 1.318^(ISS II/III)). It is encoded as (1 + 0.318 * (ISS_II + ISS_III)), which applies the same factor to either stage, per the register's instruction to supply both level indicators even when a paper pools them. Baseline ISS was derived from serum beta-2-microglobulin and albumin (Table S4 footnote b). Cohort: stage I 217 (44.1%), II 165 (33.5%), III 103 (20.9%), not reported 7 (1.4%).",
      source_name = "ISS"
    ),
    ISS_III = list(
      description = "Baseline International Staging System stage III multiple myeloma indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (ISS stage I is the reference level; ISS_II and ISS_III are mutually exclusive)",
      notes = "Shares the single pooled ISS II-or-III effect on CL0 with ISS_II; see that entry.",
      source_name = "ISS"
    )
  )

  # Screened in the stepwise covariate search (Table S2) but not retained in
  # the final model; Gong 2025 reports no point estimate for any of them.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened as a demographic covariate (Table S2); not retained. Cohort median 65 years (range 33-86) (Table S4)."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      notes = "Screened as a demographic covariate (Table S2); not retained. Cohort 42.7% female (Table S4)."
    ),
    CRCL = list(
      description = "Baseline creatinine clearance",
      units = "mL/min",
      type = "continuous",
      notes = "Screened as an intrinsic factor (Table S2); not retained. Mild or moderate renal impairment had no clinically meaningful effect on exposure (Results; Figure 4)."
    ),
    ALB = list(
      description = "Baseline serum albumin",
      units = "g/dL",
      type = "continuous",
      notes = "Screened (Table S2); not retained in its own right. Albumin enters the model indirectly because baseline ISS stage is derived from serum beta-2-microglobulin and albumin (Table S4 footnote b)."
    ),
    SBCMA = list(
      description = "Baseline soluble B-cell maturation antigen",
      units = "ng/mL",
      type = "continuous",
      notes = "Screened as a disease biomarker (Table S2); not retained."
    ),
    ADA_POS = list(
      description = "Anti-talquetamab antibody positivity",
      units = "(binary)",
      type = "binary",
      notes = "Screened as a patient-level and time-varying covariate (Table S2); not retained. ADA-positive and ADA-negative patients had comparable exposures (Results; Figure S9)."
    ),
    ECOG_GE1 = list(
      description = "Baseline ECOG performance status >= 1 indicator",
      units = "(binary)",
      type = "binary",
      notes = "Screened as a disease characteristic (Table S2); not retained. Cohort: ECOG 0 in 40.7%, 1 in 54.9%, 2 in 4.5% (Table S4)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 492L,
    n_studies = 1L,
    age_range = "33-86 years",
    age_median = "65 years",
    weight_range = "40.2-143 kg",
    weight_median = "74.8 kg",
    sex_female_pct = 42.7,
    race_ethnicity = c(White = 86.2, Black = 8.7, Asian = 2.2, NHPI = 0.4, Other = 2.4),
    disease_state = "Relapsed/refractory multiple myeloma (IMWG criteria). Phase I patients were intolerant of or had progressed on established therapies; phase II patients had received at least three prior lines including an immunomodulatory drug, a proteasome inhibitor and an anti-CD38 antibody, with or without prior T-cell redirection therapy. Baseline ISS stage I 44.1%, II 33.5%, III 20.9%; IgG myeloma 50.2%; triple-refractory 74.2%.",
    dose_range = "Intravenous 0.0015-0.18 mg/kg QW and 0.0005-0.00338 mg/kg Q2W; subcutaneous 0.005-0.8 mg/kg QW, 0.8 and 1.2 mg/kg Q2W and 1.6 mg/kg monthly. Recommended phase II doses 0.4 mg/kg SC weekly (step-up doses 0.01 and 0.06 mg/kg; 0.405 mg/kg in phase I) and 0.8 mg/kg SC every other week (step-up doses 0.01, 0.06 and 0.3 mg/kg).",
    regions = "Multicenter international (USA 152 patients in the Figure 4 subgroup analysis; other regions not tabulated).",
    renal_function = "Renal function >= 90 mL/min in 31.5%, 60 to < 90 in 45.9%, 30 to < 60 in 22.2%, < 30 in 0.4%.",
    hepatic_function = "Normal in 84.1%, impaired (mild, moderate or severe) in 15.4%, not reported in 0.4%.",
    notes = "Baseline demographics from Gong 2025 Table S4, overall column (n = 492): 100 patients dosed intravenously and 392 subcutaneously. MonumenTAL-1 is a single phase I/II trial with two registrations (NCT03399799, NCT04634552). PK data cutoff 22 April 2022."
  )

  ini({
    # Structural PK parameters: Gong 2025 Table S5 ("Parameter estimates for
    # final population PK model"); the covariate equations are footnotes a
    # and b of that table, and the structure is Figure S5. All parameters
    # were estimated; none is fixed.
    #
    # Total clearance (Table S5 footnote a, Figure S5):
    #   CL = CL0 * CLPCT + CL0 * (1 - CLPCT) * exp(-KDES * t)
    #      = CL0 * (1 + (CLPCT - 1) * (1 - exp(-KDES * t)))
    # which is the registered fractional exponential form with
    # cl_exp_famp = CLPCT - 1 (signed; negative = decline).
    lcl <- log(2.08); label("Total clearance at the start of treatment CL0 at 74.8 kg, IgG myeloma, ISS I (L/day)") # Table S5: CL0 = 2.08 L/day, RSE 6.58%
    cl_exp_famp <- 0.512 - 1; label("Fractional change of CL0 at steady state, CLPCT - 1 (unitless)") # Table S5: CLPCT = 0.512, RSE 4.32%
    lcl_exp_kdes <- log(0.0229); label("First-order rate constant for the decline of the time-dependent clearance KDES (1/day)") # Table S5: KDES = 0.0229 1/day, RSE 7.21%
    lvc <- log(4.30); label("Central volume of distribution V1 at 74.8 kg (L)") # Table S5: V1 = 4.30 L, RSE 2.54%
    lq <- log(1.01); label("Intercompartmental clearance Q (L/day)") # Table S5: Q = 1.01 L/day, RSE 6.10%
    lvp <- log(5.78); label("Peripheral volume of distribution V2 (L)") # Table S5: V2 = 5.78 L, RSE 8.59%
    lka <- log(0.138); label("First-order subcutaneous absorption rate constant Ka (1/day)") # Table S5: Ka = 0.138 1/day, RSE 5.47%
    lfdepot <- log(0.619); label("Subcutaneous bioavailability F (fraction)") # Table S5: F = 0.619, RSE 4.58%
    ld1 <- log(0.142); label("Duration of the zero-order subcutaneous input into the depot D1 (day)") # Table S5: D1 = 0.142 day, RSE 19.6%

    # Covariate effects (Table S5; footnotes a and b)
    e_wt_cl <- 0.672; label("Body-weight exponent on CL0, normalised to 74.8 kg (unitless)") # Table S5: BWT on CL0 = 0.672, RSE 20.7%
    e_wt_vc <- 0.670; label("Body-weight exponent on V1, normalised to 74.8 kg (unitless)") # Table S5: BWT on V1 = 0.670, RSE 17.0%
    e_mm_nigg_cl <- -0.547; label("Proportional change in CL0 for non-IgG vs IgG myeloma (fraction)") # Table S5: IgG or non-IgG on CL0 = -0.547, RSE 5.43%; footnote a factor 0.453
    e_iss_ii_iii_cl <- 0.318; label("Proportional change in CL0 for ISS stage II or III vs stage I (fraction)") # Table S5: ISS II or III on CL0 = 0.318, RSE 26.3%; footnote a factor 1.318

    # Interindividual variability. Table S5 reports IIV as CV% only, with no
    # covariances. Methods: random effects were 'modeled as log-normal
    # distribution', and the paper does not state its CV definition, so the
    # exact log-normal relation omega^2 = log(1 + CV^2) is used. The CLPCT
    # eta acts multiplicatively on CLPCT itself (see model()).
    etalcl ~ log(1 + 0.517^2) # Table S5: IIV CL0 51.7%, shrinkage 20.0%
    etacl_exp_famp ~ log(1 + 0.597^2) # Table S5: IIV CLPCT 59.7%, shrinkage 23.3%
    etalvc ~ log(1 + 0.223^2) # Table S5: IIV V1 22.3%, shrinkage 53.2%
    etalcl_exp_kdes ~ log(1 + 0.674^2) # Table S5: IIV KDES 67.4%, shrinkage 52.9%
    etalq ~ log(1 + 0.394^2) # Table S5: IIV Q 39.4%, shrinkage 65.9%
    etalvp ~ log(1 + 0.828^2) # Table S5: IIV V2 82.8%, shrinkage 45.1%
    etalka ~ log(1 + 0.659^2) # Table S5: IIV Ka 65.9%, shrinkage 25.7%
    etalfdepot ~ log(1 + 0.808^2) # Table S5: IIV F 80.8%, shrinkage 47.0%
    etald1 ~ log(1 + 1.453^2) # Table S5: IIV D1 145.3%, shrinkage 60.8%

    # Residual error: Table S5 'Proportional error (CV%)' = 0.220
    propSd <- 0.220; label("Proportional residual error (fraction)") # Table S5: proportional error 0.220
  })

  model({
    # 1. Individual parameters. MM_NIGG, ISS_II and ISS_III enter through the
    # proportional structure the paper describes for categorical covariates
    # (Methods): CL0 * (1 + theta * indicator). ISS_II and ISS_III are
    # mutually exclusive, so (ISS_II + ISS_III) is 1 for stage II or III.
    cl0 <- exp(lcl + etalcl) *
      (WT / 74.8)^e_wt_cl *
      (1 + e_mm_nigg_cl * MM_NIGG) *
      (1 + e_iss_ii_iii_cl * (ISS_II + ISS_III))

    # CLPCT is log-normal around its typical value 1 + cl_exp_famp (Methods:
    # all random effects log-normal); the individual amplitude follows.
    clpct <- (1 + cl_exp_famp) * exp(etacl_exp_famp)
    cl_exp_famp_i <- clpct - 1
    cl_exp_kdes <- exp(lcl_exp_kdes + etalcl_exp_kdes)

    vc <- exp(lvc + etalvc) * (WT / 74.8)^e_wt_vc
    q <- exp(lq + etalq)
    vp <- exp(lvp + etalvp)
    ka <- exp(lka + etalka)
    fdepot <- exp(lfdepot + etalfdepot)
    d1 <- exp(ld1 + etald1)

    # 2. Total clearance. `t` is time since the first talquetamab dose (the
    # paper's 'time (in days)'), so the first dose must be at t = 0.
    cl <- cl0 * (1 + cl_exp_famp_i * (1 - exp(-cl_exp_kdes * t)))

    # 3. Micro-constants
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 4. ODEs. SC doses enter `depot` zero-order over d1 (dose records need
    # rate = -2), then pass first-order into `central`. IV doses go to
    # `central`.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    f(depot) <- fdepot
    dur(depot) <- d1

    # 5. Observation: mg/L -> ng/mL
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
