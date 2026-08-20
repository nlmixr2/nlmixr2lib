VanWart_2025_telavancin <- function() {
  description <- paste(
    "Two-compartment population PK model for intravenous telavancin with a",
    "coupled epithelial lining fluid (ELF) biophase sub-model, fitted to 9,088",
    "plasma concentrations from 1,205 healthy subjects and patients with",
    "complicated skin and skin-structure infection, hospital-acquired or",
    "ventilator-associated bacterial pneumonia, or uncomplicated bacteremia",
    "pooled across 21 Phase 1-4 studies (Van Wart 2025). Total clearance is",
    "the sum of a non-renal intercept and a renal arm driven by BSA-normalized",
    "creatinine clearance through a sigmoidal Hill function, further scaled by",
    "power functions of total body weight and age and by proportional shifts",
    "for infection type. Central and peripheral volumes carry weight, age and",
    "infection-type effects (plus sex on Vc and body mass index on Vp), and",
    "distributional clearance carries weight only. An additive dialysis",
    "clearance of 1.77 L/h is gated on an active intermittent-hemodialysis",
    "session, and central volume gains a fixed additive 1.55 L in",
    "dialysis-dependent subjects sampled more than 48 h after their last",
    "session. Plasma residual variability is stratified by study phase. The",
    "ELF sub-model is a unidirectional biophase compartment whose state is the",
    "ELF concentration, giving a steady-state ELF/plasma ratio of k13/k30 =",
    "6.95% of total drug (73.0% of free drug at 90% protein binding)."
  )
  reference <- paste(
    "Van Wart SA, Safir MC, Bhavnani SM, Lodise TP, Rubino CM. Population",
    "pharmacokinetic analyses for telavancin using data from healthy subjects",
    "and patients with infections. Antimicrob Agents Chemother.",
    "2025;69(7):e01382-24. doi:10.1128/aac.01382-24.",
    "Plasma model parameters from Table 2; ELF sub-model parameters from",
    "Table S4 of the supplement (AAC01382-24-s0001.pdf).",
    sep = " "
  )
  vignette <- "VanWart_2025_telavancin"
  units    <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Plasma residual variability is stratified by study phase (Phases 1 and 4
  # pooled / Phase 2 / Phase 3), so the canonical propSd consumed by the error
  # model is derived inside model() from three phase-specific ini() magnitudes.
  # Same construction as Cammarata_2024_sulbactam_durlobactam (also an ICPD
  # antibacterial popPK paper with a phase-stratified sigma) and
  # vanIersel_2018_posaconazole.
  paper_specific_residual_sds <- c(
    "propSdPhase14", "propSdPhase2", "propSdPhase3"
  )

  compartmentData <- list(
    central     = list(analyte = "telavancin", units = "mg",   specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "telavancin", units = "mg",   specimen = "plasma", verified = TRUE),
    # The ELF biophase state holds a CONCENTRATION, not an amount: Table S4
    # reports only k13 and k30 with no ELF volume, and the resulting
    # steady-state ratio k13/k30 = 0.0695 reproduces the Table 3 median
    # total-drug ELF penetration ratio of 0.0683 directly. Had the state been
    # an amount, the penetration ratio would additionally scale by Vc / V_ELF
    # and no such volume is reported. Being a biophase (effect) compartment it
    # does not drain `central`; the plasma disposition is unchanged by the
    # ELF sub-model, which is required because Van Wart 2025 fitted the ELF
    # data sequentially with the plasma parameters fixed to the individual
    # post hoc values.
    effect      = list(analyte = "telavancin", units = "mg/L", specimen = "epithelial lining fluid", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 1205,
    n_studies      = 21,
    n_observations = 9088,
    age_range      = "18-100 years",
    age_mean       = "47.1 years (SD 18.4)",
    weight_range   = "33.6-227 kg",
    weight_mean    = "79.7 kg (SD 21.5)",
    bmi_range      = "12.3-88.8 kg/m^2",
    bmi_mean       = "27.3 kg/m^2 (SD 6.99)",
    sex_female_pct = 38.3,
    race_ethnicity = c(Caucasian = 76.1, Black = 13.1, Asian = 4.57, Other = 6.39),
    disease_state  = paste(
      "Pooled healthy subjects (33.9%) and patients with complicated skin and",
      "skin-structure infection (46.2%), hospital-acquired or",
      "ventilator-associated bacterial pneumonia (18.3%), or uncomplicated",
      "Staphylococcus aureus bacteremia (1.5%)"
    ),
    renal_function = paste(
      "Full spectrum: normal 47.3% (CLcr >= 90 mL/min/1.73 m^2), mild",
      "impairment 27.8%, moderate 15.0%, severe 9.29%, and chronic kidney",
      "disease stage 5 on intermittent hemodialysis 0.66% (8 subjects);",
      "CLcr 83.7 mL/min/1.73 m^2 (SD 36.2), range 0-203"
    ),
    dose_range     = paste(
      "0.25-15 mg/kg intravenously over 0.5-2 h as single or once-daily",
      "multiple doses (Table S1); the approved regimen is 10 mg/kg q24h",
      "infused over 1 h"
    ),
    regions        = "Not reported by region; includes a Phase 1 study in Japanese and Caucasian subjects",
    notes          = paste(
      "Baseline demographics from Van Wart 2025 Table 1 (PK analysis",
      "population, N = 1,205). The ELF sub-model was informed by only the 20",
      "healthy subjects of Phase 1 study I6424-108a, each contributing a",
      "single bronchoalveolar lavage sample at 4, 8, 12 or 24 h after the Day",
      "3 dose of 10 mg/kg q24h."
    )
  )

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power function on CL, Vc, CLd (= q) and Vp, normalized to the",
        "population mean of 79.7 kg (Table 1). Van Wart 2025 prints no",
        "covariate equations and states no normalization constant anywhere;",
        "see the vignette Errata for the three independent checks that",
        "establish mean-normalization (dimensional analysis of the Table 2",
        "'Coefficient (L)' / '(L/hour)' unit labels, the Figure 1 absolute CL",
        "scale, and the Table 3 steady-state AUC0-24)."
      ),
      source_name        = "TBW"
    ),
    AGE = list(
      description        = "Age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power function on CL, Vc and Vp, normalized to the population mean of",
        "47.1 years (Table 1). Not retained on CLd: age on CLd entered forward",
        "selection at step 13 (Table S3) but was the single relationship",
        "removed during backward elimination (P = 0.05841), and Table 2",
        "reports no CLd-age power."
      ),
      source_name        = "AGE"
    ),
    BMI = list(
      description        = "Body mass index",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Power function on Vp only, normalized to the population mean of 27.3",
        "kg/m^2 (Table 1). Retained alongside the Vp weight effect; Van Wart",
        "2025 Discussion reads the pair as separating body size from relative",
        "obesity. The exponent is negative (-0.308), so at a fixed weight a",
        "more obese subject has a smaller peripheral volume."
      ),
      source_name        = "BMI"
    ),
    CRCL = list(
      description        = "BSA-normalized creatinine clearance",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying. Drives the renal clearance arm through a sigmoidal Hill",
        "function and is the one covariate exempted from backward elimination",
        "(included a priori). Van Wart 2025 Methods: computed with ideal body",
        "weight (total body weight when total body weight was below ideal body",
        "weight), then normalized to a BSA of 1.73 m^2; serum creatinine was",
        "linearly interpolated between measurements so CLcr changes gradually",
        "during treatment. Enters the Hill function on its raw",
        "mL/min/1.73 m^2 scale (not normalized), since CLcr50 = 68.3 is",
        "reported on that scale."
      ),
      source_name        = "CLcr"
    ),
    SEXF = list(
      description        = "Female sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = paste(
        "Proportional shift on Vc only. The Table 2 coefficient is NEGATIVE",
        "(-0.0584) despite the row being labelled 'Vc-proportional increase",
        "for females', so females have 5.84% LOWER central volume than males;",
        "the 90% CI (-0.0898 to -0.0293) excludes zero on the negative side,",
        "confirming the sign is not a typesetting artifact. Encoded as",
        "(1 + e_sexf_vc * SEXF)."
      ),
      source_name        = "SEX"
    ),
    DIS_HABP = list(
      description        = "Hospital-acquired bacterial pneumonia infection-type indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not HABP; the shared all-zero reference across the four infection-type columns is the uninfected healthy subject)",
      notes              = paste(
        "Van Wart 2025 estimated ONE shared coefficient for the pooled",
        "bacteremia / HABP / VABP stratum on each of CL, Vc and Vp, so the",
        "shared coefficient is applied to (DIS_BACTEREMIA + DIS_HABP +",
        "DIS_VABP) inside model(). The three columns are deliberately kept",
        "distinct per the DIS_VABP register discipline: sibling analyses (and",
        "Cammarata 2024, which separates them) may resolve them individually."
      ),
      source_name        = "Infection type = HABP"
    ),
    DIS_VABP = list(
      description        = "Ventilator-associated bacterial pneumonia infection-type indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not VABP; shared all-zero reference is the uninfected healthy subject)",
      notes              = paste(
        "Shares one coefficient with DIS_HABP and DIS_BACTEREMIA on CL, Vc and",
        "Vp. Van Wart 2025 Table 1 pools HABP and VABP into a single 18.3%",
        "demographic stratum and never reports them separately."
      ),
      source_name        = "Infection type = VABP"
    ),
    DIS_BACTEREMIA = list(
      description        = "Uncomplicated bacteremia infection-type indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not bacteremia; shared all-zero reference is the uninfected healthy subject)",
      notes              = paste(
        "Uncomplicated Staphylococcus aureus bacteremia, contributed by the",
        "single Phase 2 study in Table 1 (18 subjects, 1.5%). Shares one",
        "coefficient with DIS_HABP and DIS_VABP on CL, Vc and Vp."
      ),
      source_name        = "Infection type = bacteremia"
    ),
    DIS_CSSSI = list(
      description        = "Complicated skin and skin-structure infection infection-type indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not cSSSI; shared all-zero reference is the uninfected healthy subject)",
      notes              = paste(
        "The largest infected stratum (557 subjects, 46.2%). Carries its own",
        "coefficient on CL, Vc and Vp, separate from the pooled",
        "bacteremia / HABP / VABP coefficient. Distinct from the",
        "severity-WITHIN-cohort indicator DIS_INFECT_CSSSI_SEV."
      ),
      source_name        = "Infection type = cSSSI"
    ),
    RRT_HEMODIAL_ACTIVE = list(
      description        = "Intermittent-hemodialysis-active indicator (time-varying per-session gate)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no hemodialysis session running; also 0 for every subject not on dialysis)",
      notes              = paste(
        "Gates the additive dialysis clearance term. Van Wart 2025:",
        "'CL_DL was estimated only during those periods where intermittent",
        "hemodialysis (IHD) was active and was fixed to a value of zero when",
        "IHD was not operative.' Informed by only 8 CKD5 subjects, which is",
        "why the IIV on CL_DL is imprecise (%SEM 135)."
      ),
      source_name        = "IHD active"
    ),
    RRT_HEMODIAL_STATUS = list(
      description        = "Dialysis-dependent (CKD5) subject indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not dialysis-dependent)",
      notes              = paste(
        "Subject-level flag identifying the 8 chronic kidney disease stage 5",
        "subjects on intermittent hemodialysis (Table 1, 0.66%). Used only to",
        "gate the additive central-volume increase together with",
        "T_POST_HEMODIAL; the dialysis CLEARANCE term is gated by the",
        "time-varying RRT_HEMODIAL_ACTIVE instead."
      ),
      source_name        = "CKD5"
    ),
    T_POST_HEMODIAL = list(
      description        = "Time elapsed since the end of the last intermittent-hemodialysis session",
      units              = "h",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying. Used only through the > 48 h threshold that switches on",
        "the additive central-volume increase in dialysis-dependent subjects.",
        "Van Wart 2025: 'An additive increase in the central volume of",
        "distribution (Vc) was estimated in the model when PK samples were",
        "collected more than 48 hours after the last active IHD session ...",
        "likely due to fluid depletion during IHD and accumulation between",
        "sessions.' Set to 0 for subjects never dialysed; the term is",
        "additionally gated by RRT_HEMODIAL_STATUS so the value is inert",
        "outside the CKD5 stratum."
      ),
      source_name        = "Time since last IHD session"
    ),
    STUDY_TLV_PHASE2 = list(
      description        = "Phase 2 study cohort indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (Phase 1 or Phase 4 study, when STUDY_TLV_PHASE3 is also 0)",
      notes              = paste(
        "Selects the Phase 2 proportional residual-error magnitude. Van Wart",
        "2025 estimated three constant-coefficient-of-variation terms by study",
        "phase, with Phases 1 and 4 combined into a single value; the additive",
        "component is shared across all phases. Paired with STUDY_TLV_PHASE3;",
        "both 0 selects the pooled Phase 1 / Phase 4 term."
      ),
      source_name        = "Study phase = 2"
    ),
    STUDY_TLV_PHASE3 = list(
      description        = "Phase 3 study cohort indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (Phase 1 or Phase 4 study, when STUDY_TLV_PHASE2 is also 0)",
      notes              = paste(
        "Selects the Phase 3 proportional residual-error magnitude. See",
        "STUDY_TLV_PHASE2 for the full phase-stratified residual-error",
        "rationale."
      ),
      source_name        = "Study phase = 3"
    )
  )

  # Screened in the Van Wart 2025 covariate analysis but NOT retained in the
  # final model, so they are documented rather than referenced in model().
  # Methods: "The following continuous covariates were investigated ...: age in
  # years, weight in kg, height in cm, BSA in m^2, BMI in kg/m^2, and CLcr in
  # mL/min/1.73 m^2 ... Categorical covariates that were investigated included
  # sex, race, and infection type." Table S3 lists the 13 relationships that
  # survived forward selection; height, BSA and race appear in none of them.
  covariatesDataExcluded <- list(
    HT = list(
      description = "Height",
      units       = "cm",
      type        = "continuous",
      notes       = paste(
        "Screened as a continuous covariate (Methods) and summarised in Table 1",
        "(171 cm, SD 10.6, range 122-203) but retained on no parameter."
      )
    ),
    BSA = list(
      description = "Body surface area",
      units       = "m^2",
      type        = "continuous",
      notes       = paste(
        "Screened as a continuous covariate and summarised in Table 1 (1.92",
        "m^2, SD 0.259, range 1.22-2.96) but retained on no parameter. BSA",
        "still enters the model indirectly, because CRCL is normalized to",
        "1.73 m^2."
      )
    ),
    RACE_BLACK = list(
      description = "Black race indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Race was screened as a categorical covariate (Methods) and summarised",
        "in Table 1 (Caucasian 76.1%, Black 13.1%, Asian 4.57%, Other 6.39%)",
        "but no race relationship survived forward selection (Table S3)."
      )
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened with race; not retained. See RACE_BLACK."
    )
  )

  ini({
    # =====================================================================
    # Clearance (Van Wart 2025 Table 2, "CL" block).
    #
    #   CL_R = CL_R,max * CLcr^hill / (CLcr50^hill + CLcr^hill)
    #   CL   = (CL_NR + CL_R) * (TBW/79.7)^0.532 * (AGE/47.1)^0.0921
    #                         * (1 + 0.418 * bacteremia/HABP/VABP)
    #                         * (1 + 0.228 * cSSSI)
    #
    # The Hill-on-CRCL idiom and the parameter names lcl_nonren / lcl_renal /
    # crcl50_cl_renal / hill_cl_renal follow Ganesan_2023_tebipenem.R, which
    # the parameter-name register cites as the preferred form for exactly this
    # structure (an asymptote and a half-max on a covariate axis).
    #
    # Results: the Hill form "was chosen over the alternative functional forms
    # tested within NONMEM (linear and power) as it provided the best apparent
    # fit and allowed the relationship between CL_R and CLcr to be relatively
    # flat across the range of normal renal function."
    # =====================================================================
    lcl_nonren <- log(0.407); label("Non-renal clearance CL_NR (L/h)")         # Table 2 CL, CL_NR = 0.407 L/h (%SEM 2.50; bootstrap median 0.402, 90% CI 0.374-0.414)
    lcl_renal  <- log(1.04);  label("Maximum renal clearance CL_R,max (L/h)")  # Table 2 CL, CL_R,max = 1.04 L/h (%SEM 2.59; bootstrap median 1.04, 90% CI 0.982-1.08)

    crcl50_cl_renal <- 68.3; label("Creatinine clearance giving half-maximal renal clearance (mL/min/1.73 m^2)") # Table 2 CL, Baseline CLcr50 = 68.3 (%SEM 3.00; bootstrap median 67.8, 90% CI 61.9-69.5)
    hill_cl_renal   <- 1.88; label("Hill coefficient of the CLcr - CL_R relationship (unitless)")                # Table 2 CL, Hill coefficient = 1.88 (%SEM 3.02; bootstrap median 1.87, 90% CI 1.76-1.98)

    e_wt_cl  <- 0.532;  label("Power exponent: total body weight on CL (unitless)") # Table 2 CL, CL-TBW power = 0.532 (%SEM 6.21; bootstrap median 0.516, 90% CI 0.422-0.598)
    e_age_cl <- 0.0921; label("Power exponent: age on CL (unitless)")               # Table 2 CL, CL-age power = 0.0921 (%SEM 16.3; bootstrap median 0.0912, 90% CI 0.0509-0.133)

    e_bacteremia_habp_vabp_cl <- 0.418; label("Proportional change in CL for bacteremia / HABP / VABP patients (unitless)") # Table 2 CL, CL-proportional increase for bacteremia/HABP/VABP patients = 0.418 (%SEM 5.80; bootstrap median 0.403, 90% CI 0.319-0.500)
    e_csssi_cl                <- 0.228; label("Proportional change in CL for cSSSI patients (unitless)")                   # Table 2 CL, CL-proportional increase for cSSSI patients = 0.228 (%SEM 6.92; bootstrap median 0.226, 90% CI 0.183-0.266)

    # =====================================================================
    # Dialysis clearance (Van Wart 2025 Table 2). An ADDITIVE arm on total
    # clearance, live only while an IHD session is running. Canonical name
    # per the RRT_HEMODIAL_ACTIVE register entry ("Pair with
    # lcl_hemodialysis when the dialysis arm is a primary estimated
    # structural parameter"), as in Veinstein_2013_gentamicin.R.
    # =====================================================================
    lcl_hemodialysis <- log(1.77); label("Dialysis clearance CL_DL during an active IHD session (L/h)") # Table 2 CL_DL = 1.77 L/h (%SEM 11.1; bootstrap median 1.77, 90% CI 1.43-2.18)

    # =====================================================================
    # Central volume (Van Wart 2025 Table 2, "Vc" block).
    #
    #   Vc = 5.75 * (TBW/79.7)^0.469 * (AGE/47.1)^0.188
    #             * (1 - 0.0584 * female)
    #             * (1 + 0.578 * bacteremia/HABP/VABP)
    #             * (1 + 0.313 * cSSSI)
    #      + 1.55 * (CKD5 and > 48 h since the last IHD session)
    # =====================================================================
    lvc <- log(5.75); label("Central volume of distribution Vc (L)") # Table 2 Vc, Coefficient = 5.75 L (%SEM 1.24; bootstrap median 5.71, 90% CI 5.54-5.85)

    e_wt_vc   <- 0.469;   label("Power exponent: total body weight on Vc (unitless)")     # Table 2 Vc, Vc-TBW power = 0.469 (%SEM 11.2; bootstrap median 0.442, 90% CI 0.303-0.584)
    e_age_vc  <- 0.188;   label("Power exponent: age on Vc (unitless)")                   # Table 2 Vc, Vc-age power = 0.188 (%SEM 9.09; bootstrap median 0.183, 90% CI 0.133-0.236)
    e_sexf_vc <- -0.0584; label("Proportional change in Vc for female sex (unitless)")    # Table 2 Vc, Vc-proportional increase for females = -0.0584 (%SEM 19.5; bootstrap median -0.0577, 90% CI -0.0898 to -0.0293). The value is negative: females have LOWER Vc despite the row label saying 'increase'

    e_bacteremia_habp_vabp_vc <- 0.578; label("Proportional change in Vc for bacteremia / HABP / VABP patients (unitless)") # Table 2 Vc, Vc-proportional increase for bacteremia/HABP/VABP patients = 0.578 (%SEM 7.39; bootstrap median 0.599, 90% CI 0.494-0.731)
    e_csssi_vc                <- 0.313; label("Proportional change in Vc for cSSSI patients (unitless)")                   # Table 2 Vc, Vc-proportional increase for cSSSI patients = 0.313 (%SEM 7.20; bootstrap median 0.300, 90% CI 0.236-0.354)

    # Held constant by the authors: "this parameter was later fixed in the
    # model to the final estimate in order to achieve successful minimization"
    # (Results, plasma model development).
    e_t_post_hemodial_vc <- fixed(1.55); label("Additive increase in Vc for dialysis-dependent subjects more than 48 h after their last IHD session (L)") # Table 2 Vc, Vc-increase for CKD5 subjects after >48 hours has elapsed between IHD sessions = 1.55 L (%SEM 25.8; bootstrap median 1.56, 90% CI 0.514-2.54)

    # =====================================================================
    # Distributional clearance CLd -> canonical lq (Van Wart 2025 Table 2,
    # "CLd" block).
    #
    #   CLd = 3.73 * (TBW/79.7)^0.772
    #
    # The Table 2 row is printed as "Vc-CLd-TBW power"; it sits under the CLd
    # heading, Vc already has its own TBW power on the preceding line, Table S3
    # step 9 is "CLd / Weight / Power", and the Abstract states "Only body
    # weight was found to be a significant predictor of the IIV in
    # distributional clearance." The leading "Vc-" is a typesetting artifact
    # carried down from the Vc rows. See the vignette Errata.
    # =====================================================================
    lq <- log(3.73); label("Distributional clearance CLd (L/h)") # Table 2 CLd, Coefficient = 3.73 L/h (%SEM 1.76; bootstrap median 3.75, 90% CI 3.64-3.96)

    e_wt_q <- 0.772; label("Power exponent: total body weight on CLd (unitless)") # Table 2 CLd, "Vc-CLd-TBW power" = 0.772 (%SEM 14.2; bootstrap median 0.764, 90% CI 0.444-1.05)

    # =====================================================================
    # Peripheral volume (Van Wart 2025 Table 2, "Vp" block).
    #
    #   Vp = 5.52 * (TBW/79.7)^0.976 * (AGE/47.1)^0.272 * (BMI/27.3)^-0.308
    #             * (1 + 0.329 * bacteremia/HABP/VABP)
    #             * (1 + 0.118 * cSSSI)
    # =====================================================================
    lvp <- log(5.52); label("Peripheral volume of distribution Vp (L)") # Table 2 Vp, Coefficient = 5.52 L (%SEM 1.30; bootstrap median 5.55, 90% CI 5.43-5.74)

    e_wt_vp  <- 0.976;  label("Power exponent: total body weight on Vp (unitless)") # Table 2 Vp, Vp-TBW power = 0.976 (%SEM 7.15; bootstrap median 0.969, 90% CI 0.824-1.14)
    e_age_vp <- 0.272;  label("Power exponent: age on Vp (unitless)")               # Table 2 Vp, Vp-age power = 0.272 (%SEM 8.28; bootstrap median 0.289, 90% CI 0.240-0.346)
    e_bmi_vp <- -0.308; label("Power exponent: body mass index on Vp (unitless)")   # Table 2 Vp, Vp-BMI power = -0.308 (%SEM 24.0; bootstrap median -0.317, 90% CI -0.468 to -0.212)

    e_bacteremia_habp_vabp_vp <- 0.329; label("Proportional change in Vp for bacteremia / HABP / VABP patients (unitless)") # Table 2 Vp, Vp-proportional increase for bacteremia/HABP/VABP patients = 0.329 (%SEM 14.6; bootstrap median 0.295, 90% CI 0.129-0.423)
    e_csssi_vp                <- 0.118; label("Proportional change in Vp for cSSSI patients (unitless)")                   # Table 2 Vp, Vp-proportional increase for cSSSI patients = 0.118 (%SEM 23.1; bootstrap median 0.125, 90% CI 0.0781-0.181)

    # =====================================================================
    # Infusion duration D1. Van Wart 2025 reports IIV on D1 (Table 2) but no
    # typical value, because D1 came from the recorded infusion start / stop
    # times: "Actual dosing (i.e., infusion start and stop) and PK sampling
    # times were used for the analysis. However, the infusion duration (D1)
    # was fixed to a value of 1 hour for subjects in the Phase 2 or 3
    # clinical studies in which this information was not recorded." One hour
    # is therefore the paper's documented default and the modal duration
    # across the 21 studies (Table S1), so it is anchored here rather than
    # estimated; the variability around it is the estimated omega^2 below.
    # =====================================================================
    ldur <- fixed(log(1)); label("Intravenous infusion duration D1, paper default (h)") # Methods, Data handling: "the infusion duration (D1) was fixed to a value of 1 hour for subjects in the Phase 2 or 3 clinical studies in which this information was not recorded"; no D1 typical value appears in Table 2

    # =====================================================================
    # Epithelial lining fluid biophase sub-model (Van Wart 2025 Table S4).
    #
    #   dCelf/dt = k13 * Cc - k30 * Celf
    #
    # Results: "By fixing the plasma PK parameters to the individual post hoc
    # values, it was possible to predict ELF concentrations using a biophase
    # model with unidirectional first-order transfer from the central
    # compartment to ELF (k13) and first-order elimination from the ELF
    # compartment (k30)."
    #
    # Structural check: k13 / k30 = 0.0107 / 0.154 = 0.0695 is the
    # steady-state ELF / plasma AUC ratio implied by this ODE, and Table 3
    # reports a median total-drug ELF penetration ratio of 0.0683 (mean
    # 0.0730) in the same 20 subjects. The free-drug ratio of 73.0% follows
    # from the 90% protein binding cited in Table 3 footnote d.
    # =====================================================================
    lk13 <- log(0.0107); label("First-order transfer rate constant from plasma to the ELF biophase (1/h)") # Table S4 k13 = 0.0107 1/h (%SEM 42.9)
    lk30 <- log(0.154);  label("First-order elimination rate constant from the ELF biophase (1/h)")        # Table S4 k30 = 0.154 1/h (%SEM 48.2)

    # =====================================================================
    # Inter-individual variability (Van Wart 2025 Table 2 for the plasma
    # terms, Table S4 for k30). The paper reports omega^2 directly; the
    # parenthetical %CV printed in Table 2 is sqrt(omega^2) * 100, not the
    # log-normal sqrt(exp(omega^2) - 1), so the variances below are the
    # published values verbatim. All IIV is exponential: "Interindividual
    # variability (IIV) was estimated for systemic clearance (CL), CL_DL, Vc,
    # distributional clearance ... (CLd), volume of distribution for the
    # peripheral compartment (Vp), and for some individuals, the duration of
    # infusion (D1) using exponential error models."
    #
    # Only two off-diagonals were retained: "A scatterplot matrix of IIV
    # terms (ETAs) ... revealed strong correlations between CL and Vc and
    # between CLd and Vp. Inclusion of these covariance terms ... resulted in
    # a more statistically significant reduction in the ... objective function
    # (441 units, P < 0.00001) than estimating the full omega block."
    # =====================================================================
    etalcl + etalvc ~ c(0.0810,
                        0.0622, 0.0783)
    # Table 2: omega^2 for CL = 0.0810 (28.5% CV, %SEM 3.26, bootstrap 90% CI 0.0688-0.0934);
    # omega^2 for Vc = 0.0783 (28.0% CV, %SEM 4.41, 90% CI 0.0586-0.0967);
    # Covariance between CL and Vc = 0.0622 (%SEM 4.64, 90% CI 0.0486-0.0770)
    # => correlation 0.0622 / sqrt(0.0810 * 0.0783) = 0.781

    etalq + etalvp ~ c(0.128,
                       0.0565, 0.0477)
    # Table 2: omega^2 for CLd = 0.128 (35.8% CV, %SEM 11.1, 90% CI 0.119-0.173);
    # omega^2 for Vp = 0.0477 (21.8% CV, %SEM 7.69, 90% CI 0.0369-0.0634);
    # Covariance between CLd and Vp = 0.0565 (%SEM 9.10, 90% CI 0.0419-0.0798)
    # => correlation 0.0565 / sqrt(0.128 * 0.0477) = 0.723

    etaldur ~ 0.0793
    # Table 2: omega^2 for D1 = 0.0793 (28.2% CV, %SEM 9.10, 90% CI 0.0727-0.0874)

    etalcl_hemodialysis ~ 0.0651
    # Table 2: omega^2 for CL_DL = 0.0651 (25.5% CV, %SEM 135, 90% CI 0.000944-0.129).
    # Retained despite the imprecision, which the Discussion attributes to
    # there being "only eight subjects in the PK database to inform this
    # parameter": "As evidenced by the relative imprecision of the IIV on
    # CL_DL (135% in the final model fit and 58.6% CV in the bootstrap) ...
    # Additional data would be needed to fully inform the variability".

    etalk30 ~ 0.113
    # Table S4: omega^2 for k30 = 0.113 (33.6% CV, %SEM 62.4, shrinkage 21.6).
    # No IIV is reported for k13.

    # =====================================================================
    # Residual variability. Van Wart 2025 reports sigma^2; the SDs below are
    # sqrt(sigma^2), which reproduces the parenthetical mg/L and %CV values
    # printed in Table S4 and so confirms that the Table 2 components are
    # variances too: sqrt(0.316) = 0.562 mg/L and sqrt(0.00985) = 9.92% CV
    # are exactly the parentheticals Table S4 prints.
    #
    # Plasma: one additive component shared across phases plus a
    # phase-specific constant-coefficient-of-variation (proportional) term.
    # "Separate CCV error terms were estimated in the model based on the study
    # phase. The base structural PK model included three separate CCV terms to
    # describe residual variability for each phase of clinical development
    # (Phases 1 and 4 were combined into a single value)."
    # =====================================================================
    addSd <- sqrt(0.326); label("Additive residual SD, plasma, all phases (mg/L)") # Table 2 Residual variability (sigma^2), Additive component = 0.326 (%SEM 3.55; bootstrap median 0.330, 90% CI 0.277-0.453) -> SD 0.571 mg/L

    propSdPhase14 <- sqrt(0.00984); label("Proportional residual SD, plasma, Phase 1 and 4 studies (fraction)") # Table 2 Residual variability, CCV component for Phase 1 and 4 studies = 0.00984 (%SEM 0.86; bootstrap median 0.00983, 90% CI 0.00833-0.0116) -> 9.92% CV
    propSdPhase2  <- sqrt(0.0169);  label("Proportional residual SD, plasma, Phase 2 studies (fraction)")        # Table 2 Residual variability, CCV component for Phase 2 studies = 0.0169 (%SEM 3.25; bootstrap median 0.0170, 90% CI 0.0134-0.0228) -> 13.0% CV
    propSdPhase3  <- sqrt(0.0441);  label("Proportional residual SD, plasma, Phase 3 studies (fraction)")        # Table 2 Residual variability, CCV component for Phase 3 studies = 0.0441 (%SEM 2.04; bootstrap median 0.0439, 90% CI 0.0347-0.0535) -> 21.0% CV

    # ELF residual variability was held at the plasma Phase 1 / 4 values:
    # "In order to allow for the estimation of IIV in k30, since only one ELF
    # measurement was collected per subject, the residual variability of ELF
    # data was fixed to be equal to that previously estimated for the plasma
    # PK data in healthy subjects from Phases 1 and 4 studies." Table S4 marks
    # both components FIXED.
    addSd_Celf  <- fixed(sqrt(0.316));   label("Additive residual SD, ELF (mg/L)")          # Table S4 Additive component for ELF = 0.316 (0.562 mg/L), FIXED
    propSd_Celf <- fixed(sqrt(0.00985)); label("Proportional residual SD, ELF (fraction)")  # Table S4 CCV component for ELF = 0.00985 (9.92% CV), FIXED
  })

  model({
    # ------------------------------------------------------------------
    # 1. Derived covariate terms.
    #
    # The four infection-type indicators are mutually exclusive and share the
    # uninfected healthy subject as their all-zero reference. Van Wart 2025
    # estimated one coefficient for the pooled bacteremia / HABP / VABP
    # stratum and a second for cSSSI, so the pooled coefficient multiplies
    # the SUM of the three pneumonia / bacteremia indicators.
    #
    # Continuous covariates are normalized to their Table 1 population means
    # (TBW 79.7 kg, AGE 47.1 years, BMI 27.3 kg/m^2) so that each Table 2
    # "Coefficient" is the typical value for the mean subject. Van Wart 2025
    # prints no covariate equation and states no normalization constant; see
    # the vignette Errata for the three checks that establish this reading.
    # CRCL is the exception: it enters the Hill function on its raw
    # mL/min/1.73 m^2 scale, because CLcr50 = 68.3 is reported on that scale.
    # ------------------------------------------------------------------
    infect     <- DIS_BACTEREMIA + DIS_HABP + DIS_VABP

    infect_cl  <- 1 + e_bacteremia_habp_vabp_cl * infect + e_csssi_cl * DIS_CSSSI
    infect_vc  <- 1 + e_bacteremia_habp_vabp_vc * infect + e_csssi_vc * DIS_CSSSI
    infect_vp  <- 1 + e_bacteremia_habp_vabp_vp * infect + e_csssi_vp * DIS_CSSSI

    wtn  <- WT  / 79.7
    agen <- AGE / 47.1
    bmin <- BMI / 27.3

    # Dialysis-dependent subject sampled more than 48 h after the last active
    # IHD session. Gated by RRT_HEMODIAL_STATUS so T_POST_HEMODIAL is inert
    # for every subject who is not on dialysis.
    post48_ihd <- RRT_HEMODIAL_STATUS * (T_POST_HEMODIAL > 48)

    # ------------------------------------------------------------------
    # 2. Individual PK parameters.
    # ------------------------------------------------------------------
    cl_nonren <- exp(lcl_nonren)
    cl_renal  <- exp(lcl_renal) *
      CRCL^hill_cl_renal / (crcl50_cl_renal^hill_cl_renal + CRCL^hill_cl_renal)

    cl <- (cl_nonren + cl_renal) *
      wtn^e_wt_cl *
      agen^e_age_cl *
      infect_cl *
      exp(etalcl)

    # Additive dialysis clearance, live only while a session is running.
    cl_hemodialysis <- exp(lcl_hemodialysis + etalcl_hemodialysis) * RRT_HEMODIAL_ACTIVE

    cl_total <- cl + cl_hemodialysis

    vc <- exp(lvc + etalvc) *
      wtn^e_wt_vc *
      agen^e_age_vc *
      (1 + e_sexf_vc * SEXF) *
      infect_vc +
      e_t_post_hemodial_vc * post48_ihd

    q <- exp(lq + etalq) * wtn^e_wt_q

    vp <- exp(lvp + etalvp) *
      wtn^e_wt_vp *
      agen^e_age_vp *
      bmin^e_bmi_vp *
      infect_vp

    dur_inf <- exp(ldur + etaldur)

    k13 <- exp(lk13)
    k30 <- exp(lk30 + etalk30)

    # ------------------------------------------------------------------
    # 3. Micro-constants.
    # ------------------------------------------------------------------
    kel <- cl_total / vc
    k12 <- q / vc
    k21 <- q / vp

    # ------------------------------------------------------------------
    # 4. Two-compartment IV disposition plus the ELF biophase.
    #
    # `effect` holds the ELF CONCENTRATION and is driven by the plasma
    # concentration central / vc. It deliberately does NOT appear in
    # d/dt(central): a biophase compartment carries no mass balance, which is
    # what makes the sequential fit valid (Van Wart 2025 fixed the plasma
    # parameters to their individual post hoc values before fitting the ELF
    # data, so adding the ELF compartment must leave the plasma profile
    # unchanged).
    # ------------------------------------------------------------------
    d/dt(central)     <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <-   k12 * central        - k21 * peripheral1
    d/dt(effect)      <-   k13 * (central / vc) - k30 * effect

    # ------------------------------------------------------------------
    # 5. Zero-order intravenous input over the modelled infusion duration.
    #    Dose rows must be given with rate = -2 so rxode2 uses dur().
    # ------------------------------------------------------------------
    dur(central) <- dur_inf

    # ------------------------------------------------------------------
    # 6. Observations. Dose in mg and volumes in L give mg/L.
    #
    #    Plasma residual variability is switched by study phase; both phase
    #    indicators 0 selects the pooled Phase 1 / Phase 4 term.
    # ------------------------------------------------------------------
    Cc   <- central / vc
    Celf <- effect

    phase14 <- 1 - STUDY_TLV_PHASE2 - STUDY_TLV_PHASE3
    propSd  <- propSdPhase14 * phase14 +
               propSdPhase2  * STUDY_TLV_PHASE2 +
               propSdPhase3  * STUDY_TLV_PHASE3

    Cc   ~ add(addSd)      + prop(propSd)
    Celf ~ add(addSd_Celf) + prop(propSd_Celf)
  })
}
