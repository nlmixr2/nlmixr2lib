Khaowroongrueng_2024_sufentanil <- function() {
  description <- paste(
    "Two-compartment IV population PK model for sufentanil in adult Korean",
    "patients undergoing cardiac surgery with cardiopulmonary bypass (CPB).",
    "First-order elimination, no absorption (IV bolus plus continuous",
    "infusion into the central compartment). Two time-varying",
    "intra-operative CPB-phase indicators carry the CPB effect: clearance",
    "increases 2.80-fold during the CPB and rewarming phases, and the",
    "central volume increases 2.74-fold during the CPB phase (hemodilution",
    "from the 1700 mL circuit priming volume). The post-CPB phase collapses",
    "to the pre-CPB reference. Inter-individual variability was estimated on",
    "clearance only (24.2 %CV); Table 2 reports the V1, Q and V2 IIV terms as",
    "0 FIX. Proportional residual error 25.3%. Parameter values from",
    "Khaowroongrueng 2024 Table 2."
  )
  reference <- paste(
    "Khaowroongrueng V, Son KH, Lee S-M, Lee J, Park C-G, Lee SI, Shin D,",
    "Shin K-H. Population pharmacokinetic modeling of sufentanil in adult",
    "Korean patients undergoing cardiopulmonary bypass surgery.",
    "CPT Pharmacometrics Syst Pharmacol 2024;13(10):1682-1692.",
    "doi:10.1002/psp4.13205.",
    sep = " "
  )
  vignette <- "Khaowroongrueng_2024_sufentanil"
  units <- list(
    time          = "h",
    dosing        = "ug",
    concentration = "ng/mL"
  )

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Arterial blood samples were assayed for sufentanil
  # plasma concentration by LC-MS/MS (Methods "Study sample handling"), so
  # the observed matrix is plasma.
  compartmentData <- list(
    central     = list(analyte = "sufentanil", units = "ug", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "sufentanil", units = "ug", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CPB_ON = list(
      description        = "Cardiopulmonary bypass phase indicator (before rewarming begins)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (pre-CPB or post-CPB)",
      notes              = paste(
        "Time-varying within subject. 1 = the record falls in the CPB phase",
        "proper, i.e. from commencement of cardiopulmonary bypass until",
        "rewarming begins; 0 otherwise. Khaowroongrueng 2024 Methods",
        "'CPB-adjusted model' records the phase as a single 0/1/2/3 column",
        "(0 = pre-CPB, 1 = CPB, 2 = warming, 3 = post-CPB); CPB_ON is the",
        "decomposed indicator for that source column's level 1 and",
        "CPB_REWARM for level 2. Both are 0 during the pre-CPB (level 0) and",
        "post-CPB (level 3) phases, so the post-CPB phase collapses to the",
        "pre-CPB reference. CPB_ON alone carries the central-volume effect",
        "(V1 increases 2.74-fold during CPB only, Table 2); the clearance",
        "effect spans both windows and therefore enters as",
        "CPB_ON + CPB_REWARM. The patient is still on the bypass circuit",
        "during rewarming, so 'on CPB in the broad sense' is",
        "CPB_ON + CPB_REWARM == 1, not CPB_ON == 1."
      ),
      source_name        = "CPB phase (0/1/2/3, level 1)"
    ),
    CPB_REWARM = list(
      description        = "Cardiopulmonary bypass rewarming phase indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (pre-CPB, CPB proper, or post-CPB)",
      notes              = paste(
        "Time-varying within subject. 1 = the record falls in the rewarming",
        "phase, i.e. from the start of rewarming until separation from",
        "cardiopulmonary bypass; 0 otherwise. Decoded from level 2 of the",
        "source paper's single 0/1/2/3 CPB-phase column (Methods",
        "'CPB-adjusted model'). Mutually exclusive with CPB_ON. In this study",
        "all patients were rewarmed to a rectal temperature of at least 36 C",
        "before separation from CPB (Methods 'CPB system'), and a 50 ug",
        "sufentanil bolus was given at the onset of rewarming (Methods 'Drug",
        "administration'). Khaowroongrueng 2024 Results found the CPB effect",
        "on clearance during the CPB and rewarming phases to be similar and",
        "so estimated a single typical value spanning both windows; the",
        "central-volume effect was retained for the CPB phase only."
      ),
      source_name        = "CPB phase (0/1/2/3, level 2)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 20L,
    n_studies      = 1L,
    age_range      = "23-77 years",
    age_median     = "66 years",
    weight_range   = "46.7-99.8 kg",
    weight_median  = "66.4 kg",
    height_range   = "150.1-189.2 cm",
    bmi_range      = "16.0-34.3 kg/m^2",
    sex_female_pct = 36.4,
    race_ethnicity = c(Asian = 100),
    disease_state  = paste(
      "Adults (>19 years) undergoing cardiac surgery with cardiopulmonary",
      "bypass support: valvuloplasty or valve replacement (13), coronary",
      "artery bypass graft (6), left ventricular assist device (1), and",
      "combination surgery (2). Patients undergoing renal replacement",
      "therapy or taking medications known to interact with sufentanil were",
      "excluded. Median (range) CPB duration 2.3 (1.4-5.6) h; median total",
      "surgical duration 7.7 (5.3-14.8) h. CPB was run under mild",
      "hypothermia (core temperature >32 C), alpha-stat pH management and",
      "non-pulsatile flow during aortic cross-clamping, with the circuit",
      "primed with 1700 mL Plasma Solution A, 100 mL 20% albumin,",
      "2000 units heparin and 40 mEq bicarbonate."
    ),
    dose_range     = paste(
      "IV bolus of 1 ug/kg (generally 50 ug) at intubation, plus two further",
      "50 ug IV boluses at the onset of CPB and at the onset of rewarming.",
      "Continuous IV infusion at 20 ug/h between incision and CPB",
      "commencement and between CPB rewarming and the end of the procedure,",
      "adjusted at the anaesthetist's discretion. Median (range) total IV",
      "bolus dose 2.3 (1.5-3.3) ug/kg and total IV infusion dose",
      "2.1 (0.9-4.3) ug/kg (Table 1)."
    ),
    regions        = "South Korea (single-centre, Gachon University Gil Medical Center, Incheon)",
    notes          = paste(
      "Prospective PK study conducted May 2021 to June 2022 (IRB",
      "GDIRB 2021-184; CRIS KCT0008742). 22 participants enrolled and",
      "described in Table 1; 110 samples collected. Two participants (12",
      "observations) were excluded for markedly high pre-dose sufentanil",
      "concentrations and two further outlier samples were removed as",
      "suspected sampling errors, leaving 20 participants with 96",
      "observations in the final dataset (Results 'Structural model').",
      "The demographic ranges above are Table 1 values for the 22 enrolled",
      "participants; the paper does not re-tabulate demographics for the 20",
      "analysed. Sex distribution 14 male / 8 female of 22 enrolled.",
      "Arterial samples were drawn pre-dose, 1 h after the first dose,",
      "pre-CPB initiation, 1 h after CPB commencement, and 1, 2 and 3 h",
      "after sufentanil cessation; sufentanil was quantified by validated",
      "LC-MS/MS. Age, height, weight, body mass index, sex, CYP3A4*1G",
      "genotype (rs2242480; 13 GG / 6 GA / 3 AA), fluid replacement volumes,",
      "and the time-varying clinical laboratory parameters and inflammatory",
      "cytokines that changed significantly during CPB (total protein, IL-6,",
      "IL-8, TNF-alpha) were all screened as covariates; none was retained.",
      "NONMEM 7.5 with FOCE-I, PsN 5.3.0, Xpose 3.4.4, Pirana 21.11.1.",
      "Evaluated by 1000-sample bootstrap (99% success rate), prediction-",
      "corrected VPC (n = 1000) and NPDE (n = 1000; mean 0.05403,",
      "variance 1.058)."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Structural PK parameters -- Khaowroongrueng 2024 Table 2 ("Final
    # model", Estimate column). Two-compartment model with first-order
    # elimination, parameterised in CL / V1 / Q / V2 (Methods "Structural
    # model"); NONMEM ADVAN subroutines, FOCE-I.
    # These are the PRE-CPB (reference-phase) typical values -- the CPB
    # multipliers below are applied inside model() via the phase
    # indicators, exactly as in the paper's Equations 1-4.
    # Clearances are in L/h and volumes in L, so with amounts in ug the
    # observed concentration Cc = central / vc comes out in ug/L = ng/mL,
    # the paper's reporting unit.
    # ------------------------------------------------------------------
    lcl <- log(51.8); label("Clearance (L/h)")                              # Table 2: CL = 51.8 L/h (RSE 16.4%)
    lvc <- log(72.5); label("Central volume of distribution (L)")           # Table 2: V1 = 72.5 L (RSE 12.1%)
    lq  <- log(55.2); label("Inter-compartmental clearance (L/h)")          # Table 2: Q = 55.2 L/h (RSE 17.2%)
    lvp <- log(390);  label("Peripheral volume of distribution (L)")        # Table 2: V2 = 390 L (RSE 43.8%)

    # ------------------------------------------------------------------
    # CPB-phase covariate effects -- Khaowroongrueng 2024 Table 2.
    # The paper's Equations 1-4 give, for the jth parameter,
    #     pre-CPB:   P_i,j = P_pop,j * exp(eta_i,j) * 1
    #     CPB:       P_i,j = P_pop,j * exp(eta_i,j) * COV_CPB
    #     warming:   P_i,j = P_pop,j * exp(eta_i,j) * COV_warming
    #     post-CPB:  P_i,j = P_pop,j * exp(eta_i,j) * COV_post-CPB
    # i.e. a plain multiplicative fold-change on the typical value within
    # each phase. Results "CPB-adjusted model": the CL effect during the
    # CPB and warming phases "exhibited similarity; therefore, the CPB
    # effect in these phases was estimated using the same typical value",
    # and the V1 effect was retained for the CPB phase only. Only
    # CL_CPB&warming and V1_CPB appear in Table 2, so COV_post-CPB was NOT
    # retained and the post-CPB phase reverts to the pre-CPB reference.
    # The effects on Q and V2 were "not incorporated due to the
    # considerably high relative standard error" (Results), so neither is
    # phase-dependent here.
    # Applied as fold-change^indicator inside model(), which reproduces the
    # multiplicative form above exactly because the indicators are binary
    # (the register documents this same shape for INTRAOP and
    # SURG_SEV_MAJOR).
    # ------------------------------------------------------------------
    e_cpb_cl <- 2.80; label("Fold-change in CL during the CPB and rewarming phases (unitless)")  # Table 2: CL CPB&warming = 2.80 (RSE 14.4%)
    e_cpb_vc <- 2.74; label("Fold-change in V1 during the CPB phase (unitless)")                 # Table 2: V1 CPB = 2.74 (RSE 18.5%)

    # ------------------------------------------------------------------
    # IIV -- Khaowroongrueng 2024 Table 2 "Random effects".
    # Methods "Structural model": "A log-normal distribution was assumed
    # for the IIV of the pharmacokinetic parameters", and Table 2 reports
    # IIV as a coefficient of variation in percent. For a log-normal
    # random effect the internal variance is therefore
    #     omega^2 = log(CV^2 + 1) = log(0.242^2 + 1) = 0.056921.
    # IIV was estimated on CL only (24.2 %CV, RSE 22.4%, shrinkage 10.5%).
    # Table 2 reports the V1, Q and V2 IIV terms as "0 FIX", i.e. no
    # between-subject variability was estimated on those parameters (the
    # Discussion attributes this to the small sample size limiting "the
    # estimation of IIV for the pharmacokinetic parameters"), so no
    # etalvc / etalq / etalvp term is included here.
    # ------------------------------------------------------------------
    etalcl ~ log(0.242^2 + 1)  # Table 2: IIV CL = 24.2 %CV -> omega^2 = log(CV^2 + 1) = 0.056921

    # ------------------------------------------------------------------
    # Residual error -- Khaowroongrueng 2024 Table 2: RSV sigma_prop =
    # 25.3% (RSE 8.10%, shrinkage 7.20%). Results "Structural model": "The
    # proportional error was selected based on the decrease in the OFV",
    # so the residual model is purely proportional with no additive term.
    # ------------------------------------------------------------------
    propSd <- 0.253; label("Proportional residual SD (fraction)")  # Table 2: RSV sigma_prop = 25.3%
  })

  model({
    # 1. Individual PK parameters.
    # The two binary CPB-phase indicators are mutually exclusive, so
    # CPB_ON + CPB_REWARM is 1 while the patient is on the bypass circuit
    # (either before or during rewarming) and 0 pre- and post-CPB. The
    # fold-change^indicator form therefore multiplies the typical value by
    # e_cpb_cl during CPB and rewarming and by 1 otherwise, matching the
    # paper's Equations 1-4. V1 carries the CPB-phase effect only.
    # Q and V2 are phase-independent (their CPB effects were not retained).
    cl <- exp(lcl + etalcl) * e_cpb_cl^(CPB_ON + CPB_REWARM)
    vc <- exp(lvc)          * e_cpb_vc^CPB_ON
    q  <- exp(lq)
    vp <- exp(lvp)

    # 2. Micro-constants (1/h with cl, q in L/h and vc, vp in L).
    # cl, vc and hence kel and k12 are time-varying across CPB phases.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # 3. ODE system -- two-compartment disposition with IV dosing only
    # (bolus doses at intubation, CPB onset and rewarming onset plus
    # continuous infusion, all into the central compartment; sufentanil
    # was given via the internal jugular vein, Methods "Drug
    # administration"). There is no absorption compartment.
    d/dt(central)     <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <-         k12  * central - k21 * peripheral1

    # 4. Observation and error.
    # central in ug, vc in L -> Cc in ug/L = ng/mL, matching the paper's
    # 0.15-0.7 ng/mL target concentrations (Methods "Simulations").
    Cc <- central / vc

    Cc ~ prop(propSd)
  })
}
