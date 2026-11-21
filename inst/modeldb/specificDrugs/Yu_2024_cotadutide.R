Yu_2024_cotadutide <- function() {
  description <- "Discrete-time Markov model of the daily nausea state (0 = none, 1 = mild, 2 = moderate/severe) for once-daily subcutaneous cotadutide (MEDI0382, a GLP-1/glucagon receptor dual agonist) in adults with type 2 diabetes, obesity or NASH (Yu 2024). From each current state the logit of leaving it towards the first other state is a baseline plus a linear effect of the daily cotadutide AUC (current daily dose / CL/F) plus a linear tolerance term in the days since the last change of the daily dose level; the logit towards the second other state is offset from it by a second baseline. The first day's state follows pre-treatment probabilities. Sex and indication (obesity, NASH) shift the baselines additively and NASH and obesity scale the 0 -> {1,2} tolerance slope. A one-compartment first-order-absorption PK layer carries the typical CL/F, V/F and ka of the upstream cotadutide popPK. The 9 transition probabilities and 3 initial-state probabilities are model outputs; the Markov chain is advanced day by day outside the solve, as in the vignette."
  reference <- paste(
    "Yu H, Ueckert S, Zhou L, Cheng J, Robertson D, Hansen L, Flor A, Parker V,",
    "Hamren B, Khan AA. Exposure-response modeling for nausea incidence for",
    "cotadutide using a Markov modeling approach. CPT Pharmacometrics Syst",
    "Pharmacol. 2024;13(9):1582-1594. doi:10.1002/psp4.13194.",
    "Typical PK values from the cotadutide popPK cited in its Methods: Yu H, et al.",
    "Clin Pharmacokinet. 2024;63:255. doi:10.1007/s40262-023-01337-0."
  )
  vignette <- "Yu_2024_cotadutide"
  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "ng/mL",
    outcome = "daily nausea-state transition probabilities (fraction)"
  )

  compartmentData <- list(
    depot = list(analyte = "cotadutide", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "cotadutide", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    DOSE_COTADUTIDE_MG = list(
      description = "Current once-daily cotadutide dose level in mg, time-varying along the titration schedule; 0 on a day with no dose.",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "The control stream's DOSE column (ug) is the daily dose level on each record; the nausea model's exposure is AUC = DOSE / CL (Appendix S1 $ERROR: AUC = DOSE/CL, EFF = AUC). This model multiplies the mg column by 1000 so the AUC is in ug*h/L, the scale on which Table 3's alpha slopes were estimated (600 ug QD with CL/F 1.05 L/h gives AUC 571). The column also defines the tolerance clock: T_DOSECHANGE resets whenever its value changes.",
      source_name = "DOSE"
    ),
    T_DOSECHANGE = list(
      description = "Time since the most recent change of the daily cotadutide dose level, in hours (0 on the day of the change).",
      units = "h",
      type = "continuous",
      reference_category = NULL,
      notes = "The paper's d_t = DAY - T_esc in whole days (Appendix S1: T_esc = DAY whenever DOSE differs from the previous record's DOSE). Supply 24 * (whole days since the change); the model divides by 24 so the Table 3 beta slopes keep their per-day units. The first treatment day counts as a change (previous DOSE = 0), as does a skipped day (DOSE 0) and the resumption after it. The Methods text calls the clock the time since the last dose escalation; the readings differ only when the dose goes down or stops (see the vignette).",
      source_name = "d_t (DAY - T_esc)"
    ),
    SEXF = list(
      description = "Biological sex, 1 = female, 0 = male.",
      units = "(binary)",
      type = "binary",
      reference_category = "1 (female); the paper's effect is on males (forest-plot reference 'T2DM female').",
      notes = "Table 3 'Effect of sex on BL1,2|0' = -0.39, applied to males as an additive shift of the 0 -> {1,2} logit: male = 1 - SEXF. The Results report males with a 24% lower nausea rate, which only the additive-logit reading reproduces (see the vignette).",
      source_name = "SEX"
    ),
    DIS_OBESE = list(
      description = "Obesity-indication participant (study D5672C00001, BMI >= 35 kg/m^2 at screening), 1 = yes, 0 = no.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (T2DM participant).",
      notes = "Additive logit shifts of +1.99 on the 1 -> {0,2} and +1.53 on the 2 -> {0,1} transitions, and a factor exp(-0.75) on the 0 -> {1,2} tolerance slope (Table 3). 34 of 675 participants.",
      source_name = "indication = obesity"
    ),
    DIS_NASH = list(
      description = "NASH-indication participant (study D5671C00002: biopsy-confirmed NAFLD/NASH with NAS >= 4, fibrosis F1-F3), 1 = yes, 0 = no.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (T2DM participant).",
      notes = "Additive logit shift of +1.44 on the 2 -> {0,1} transition and a factor exp(-0.58) on the 0 -> {1,2} tolerance slope (Table 3). 49 of 675 participants. NASH is now termed MASH.",
      source_name = "indication = NASH"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 675L,
    n_studies = 4L,
    age_range = "19-82 years (median 57)",
    weight_range = "56-207 kg (median 98)",
    sex_female_pct = 54.4,
    race_ethnicity = c(Caucasian = 92.4, Black = 2.4, Asian = 4.3, Other = 0.9),
    disease_state = "Type 2 diabetes (592), obesity with BMI >= 35 kg/m^2 (34) and biopsy-confirmed NAFLD/NASH (49).",
    dose_range = "Cotadutide 50-600 ug once daily subcutaneously, titrated weekly, biweekly or monthly (steps of 50-150 ug or two-fold) to maintenance doses of 100-600 ug.",
    regions = "Studies D5670C00004, D5671C00002, D5671C00003 (Japanese T2DM) and D5672C00001.",
    biomarkers = "Any-time nausea: mild 23.7%, moderate/severe 9.9% overall (Table 2).",
    notes = "Daily nausea state derived from reported start and end dates; the most severe event of the day was kept. Demographics from Yu 2024 Table 2; titration schemes from Appendix S1."
  )

  ini({
    # Cotadutide PK: typical values of the upstream popPK printed in Yu 2024
    # Methods 'PK model'. The paper drove the nausea model with individual
    # post hoc CL/F; this layer carries the typical value.
    lka <- fixed(log(0.38)); label("Cotadutide absorption rate constant (1/h)") # Methods 'PK model': Ka = 0.38 h-1
    lcl <- fixed(log(1.05)); label("Cotadutide apparent clearance CL/F (L/h)") # Methods 'PK model': CL/F = 1.05 L/h
    lvc <- fixed(log(20.0)); label("Cotadutide apparent volume of distribution V/F (L)") # Methods 'PK model': V/F = 20.0 L

    # Initial (pre-treatment) state, Table 3; Appendix S1 A00 / A01.
    bl_pre0 <- 2.66; label("Logit of state 0 on the first day (unitless logit)") # Table 3 BL0|-1 = 2.66
    bl_pre1 <- 1.22; label("Logit of state 1 among non-0 first-day states (unitless logit)") # Table 3 BL1|-1 = 1.22

    # Transition baselines, Table 3. bl_<from><to>: the second-listed move
    # from each state is an offset on the first (Appendix S1 LP02 = A2 + LP01).
    bl_01 <- -3.79; label("Baseline logit of moving from state 0 to 1 (unitless logit)") # Table 3 BL1|0 = -3.79
    bl_02 <- -0.80; label("Offset of the 0 -> 2 logit on the 0 -> 1 logit (unitless logit)") # Table 3 BL2|0 = -0.80
    bl_10 <- -2.53; label("Baseline logit of moving from state 1 to 0 (unitless logit)") # Table 3 BL0|1 = -2.53
    bl_12 <- -4.23; label("Offset of the 1 -> 2 logit on the 1 -> 0 logit (unitless logit)") # Table 3 BL2|1 = -4.23
    bl_20 <- -2.47; label("Baseline logit of moving from state 2 to 0 (unitless logit)") # Table 3 BL0|2 = -2.47
    bl_21 <- -3.56; label("Offset of the 2 -> 1 logit on the 2 -> 0 logit (unitless logit)") # Table 3 BL1|2 = -3.56

    # Exposure slopes (alpha), per ug*h/L of daily AUC = DOSE / CL.
    slp_auc_0 <- 0.00091; label("AUC slope on the logits out of state 0 (per ug*h/L)") # Table 3 alpha 1,2|0 = 0.00091
    slp_auc_1 <- -0.00091; label("AUC slope on the logits out of state 1 (per ug*h/L)") # Table 3 alpha 0,2|1 = -0.00091
    slp_auc_2 <- -0.0031; label("AUC slope on the logits out of state 2 (per ug*h/L)") # Table 3 alpha 0,1|2 = -0.0031

    # Tolerance slopes (beta), per day since the last dose-level change.
    slp_t_0 <- -0.34; label("Tolerance slope on the logits out of state 0 (per day)") # Table 3 beta 1,2|0 = -0.34
    slp_t_1 <- -0.0089; label("Tolerance slope on the logits out of state 1 (per day)") # Table 3 beta 0,2|1 = -0.0089
    slp_t_2 <- -0.014; label("Tolerance slope on the logits out of state 2 (per day)") # Table 3 beta 0,1|2 = -0.014

    # Covariate effects, Table 3. Baseline effects are additive on the logit;
    # the indication effects on beta 1,2|0 are exponential (see the vignette).
    e_male_bl_01 <- -0.39; label("Additive shift of the 0 -> {1,2} logit for males (unitless logit)") # Table 3 'Effect of sex on BL1,2|0' = -0.39
    e_obese_bl_10 <- 1.99; label("Additive shift of the 1 -> {0,2} logit for obesity participants (unitless logit)") # Table 3 'Effect of Obesity on BL0,2|1' = 1.99
    e_nash_bl_20 <- 1.44; label("Additive shift of the 2 -> {0,1} logit for NASH participants (unitless logit)") # Table 3 'Effect of NASH on BL0,1|2' = 1.44
    e_obese_bl_20 <- 1.53; label("Additive shift of the 2 -> {0,1} logit for obesity participants (unitless logit)") # Table 3 'Effect of Obesity on BL0,1|2' = 1.53
    e_nash_slp_t_0 <- -0.58; label("Log-scale effect of NASH on the 0 -> {1,2} tolerance slope (unitless)") # Table 3 'Effect of NASH on beta1,2|0' = -0.58
    e_obese_slp_t_0 <- -0.75; label("Log-scale effect of obesity on the 0 -> {1,2} tolerance slope (unitless)") # Table 3 'Effect of Obesity on beta1,2|0' = -0.75

    # Between-subject variability, Table 3 %CV with the footnote's
    # CV = 100 * sqrt(exp(variance) - 1), so variance = log(1 + CV^2).
    etaslp_t_0 ~ log(1 + 1.93^2) # Table 3 'BSV beta 1,2|0' 193%CV, exponential on slp_t_0
    etabl_01 ~ log(1 + 1.76^2) # Table 3 'BSV BL1,2|0' 176%CV, additive on the 0 -> {1,2} logits
    etabl_10 ~ log(1 + 3.12^2) # Table 3 'BSV BL0,2|1' 312%CV, additive on the 1 -> {0,2} logits
    etabl_20 ~ log(1 + 1.38^2) # Table 3 'BSV BL0,1|2' 138%CV, additive on the 2 -> {0,1} logits
  })

  model({
    # Cotadutide PK (one compartment, first-order absorption).
    ka <- exp(lka)
    cl <- exp(lcl)
    vc <- exp(lvc)
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - (cl / vc) * central
    # mg / L * 1000 = ng/mL
    Cc <- 1000 * central / vc

    # Nausea-model drivers (Appendix S1 $ERROR): AUC = DOSE / CL with DOSE in
    # ug, and d_t in days since the last change of the daily dose level.
    auc <- 1000 * DOSE_COTADUTIDE_MG / cl
    d_t <- T_DOSECHANGE / 24

    male <- 1 - SEXF
    slp_t0_i <- slp_t_0 *
      exp(e_nash_slp_t_0 * DIS_NASH + e_obese_slp_t_0 * DIS_OBESE + etaslp_t_0)

    # Transition logits (Appendix S1 LP01 ... LP21).
    lp01 <- bl_01 + e_male_bl_01 * male + slp_auc_0 * auc + slp_t0_i * d_t + etabl_01
    lp02 <- bl_02 + lp01
    lp10 <- bl_10 + e_obese_bl_10 * DIS_OBESE + slp_auc_1 * auc + slp_t_1 * d_t + etabl_10
    lp12 <- bl_12 + lp10
    lp20 <- bl_20 + e_nash_bl_20 * DIS_NASH + e_obese_bl_20 * DIS_OBESE +
      slp_auc_2 * auc + slp_t_2 * d_t + etabl_20
    lp21 <- bl_21 + lp20

    # Transition probabilities p<from><to> (Equations 4-6, Appendix S1 Y).
    # The first-listed move takes its probability directly; the second is
    # conditional on not taking the first. Each row sums to 1.
    p01 <- expit(lp01)
    p02 <- expit(lp02) * (1 - p01)
    p00 <- 1 - p01 - p02
    p10 <- expit(lp10)
    p12 <- expit(lp12) * (1 - p10)
    p11 <- 1 - p10 - p12
    p20 <- expit(lp20)
    p21 <- expit(lp21) * (1 - p20)
    p22 <- 1 - p20 - p21

    # First-day state probabilities from the pre-treatment logits.
    p_init0 <- expit(bl_pre0)
    p_init1 <- expit(bl_pre1) * (1 - p_init0)
    p_init2 <- 1 - p_init0 - p_init1
  })
}
