Chen_2026_danicopan <- function() {
  description <- paste(
    "Population PK and PK-PD model for oral danicopan, a complement factor D",
    "inhibitor given as add-on therapy to ravulizumab or eculizumab, in 316",
    "healthy participants and 91 patients with paroxysmal nocturnal",
    "hemoglobinuria (PNH) pooled from 14 studies (Chen 2026). PK is",
    "two-compartment with linear elimination; absorption is a zero-order",
    "release into the depot over D1 followed by first-order absorption (ka)",
    "after a lag time. Body weight scales CL/F and Q/F (exponent 0.645) and",
    "Vc/F and Vp/F (exponent 0.876) relative to 75 kg; female sex lowers",
    "CL/F by 22.1% and severe renal impairment by 43.3%. Formulation",
    "(liquid-filled capsule vs tablet), food status (fasted, standard meal,",
    "high-fat meal) and the per-administration dose act on ka, D1, relative",
    "bioavailability and the lag time. The PK-PD layer is a direct",
    "inhibitory Emax model of serum complement alternative pathway (AP)",
    "activity by hemolysis assay (percent), fitted to 69 patients with PNH",
    "from the phase 3 ALPHA study: AP falls from a baseline of 29.1%",
    "(scaled by (WT/69)^1.29 and 46% lower in East Asian participants)",
    "towards a floor of 3.05% with an IC50 of 12 ng/mL.",
    sep = " "
  )
  reference <- paste(
    "Chen J, Yang F, Li H, Pradhan R, Ortiz S. Population PK-PD Modeling of",
    "Danicopan Add-On Therapy in Participants With Paroxysmal Nocturnal",
    "Hemoglobinuria Treated With Ravulizumab or Eculizumab. CPT",
    "Pharmacometrics Syst Pharmacol. 2026;15:e70209.",
    "doi:10.1002/psp4.70209.",
    sep = " "
  )
  vignette <- "Chen_2026_danicopan"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # The depot carries the zero-order released dose that then absorbs at ka;
  # danicopan was measured in plasma by validated LC (Methods, Data Sources).
  compartmentData <- list(
    depot = list(analyte = "danicopan", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "danicopan", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "danicopan", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters the PK model as power functions centred on 75 kg (CL/F and",
        "Q/F exponent 0.645, Vc/F and Vp/F exponent 0.876; Chen 2026 PK",
        "equations and Table 2) and the PK-PD model as a power function on",
        "baseline AP activity centred on 69 kg (exponent 1.29; AP0 equation",
        "and Table 4). The PK-PD text calls it baseline body weight; the PK",
        "control stream uses the record-level WTKG column."
      ),
      source_name = "WTKG / BWT"
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = paste(
        "1 = female. CL/F is multiplied by (1 - 0.221) in women (Table 2",
        "'Female sex on CL/F' = -0.221; supplement control stream",
        "CLSEXF = 1 + THETA(22) when SEXF = 1)."
      ),
      source_name = "SEXF"
    ),
    RENALIMP_SEV = list(
      description = "Severe renal impairment indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (normal renal function or mild / moderate impairment)",
      notes = paste(
        "1 = severe renal impairment. CL/F is multiplied by (1 - 0.433)",
        "(Table 2 'RI on CL/F' = -0.433; the Table 2 abbreviations define RI",
        "as severe renal impairment, and the control stream applies the",
        "effect when its RI column equals 2). Table 1 classifies renal",
        "function by eGFR: severe 15-29 mL/min/1.73 m^2 (4 participants) and",
        "end-stage < 15 mL/min/1.73 m^2 or dialysis (3 participants), all",
        "healthy-volunteer renal-impairment study participants; the",
        "Discussion refers to 'eight participants with severe renal",
        "impairment'. The paper does not state whether the end-stage",
        "participants were coded into the severe group; set 1 for eGFR",
        "< 30 mL/min/1.73 m^2. The Results report no PK change for mild or",
        "moderate impairment."
      ),
      source_name = "RI"
    ),
    FED = list(
      description = "Fed-state dose-record indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (fasted)",
      notes = paste(
        "Per-dose-record indicator. 1 = dose taken with a meal (standard or",
        "high-fat), 0 = fasted. The paper's three-level food status (fasted",
        "/ standard meal / high-fat meal; control stream FOOD = 0 / 1 / 2)",
        "has the standard meal as reference, so the paper's FAST indicator",
        "is 1 - FED and its HIFAT indicator is FED_HIGHFAT: standard meal is",
        "FED = 1, FED_HIGHFAT = 0; fasted is FED = 0, FED_HIGHFAT = 0;",
        "high-fat meal is FED = 1, FED_HIGHFAT = 1. Fasting acts on ka",
        "(additive +30.9 1/h), D1 (exp(0.0787)), relative bioavailability",
        "(exp(-0.143)) and removes the tablet lag time. All phase 2/3 PNH",
        "dosing used standard meals (Data Sources)."
      ),
      source_name = "FAST (= 1 - FED); FOOD"
    ),
    FED_HIGHFAT = list(
      description = "High-fat-meal-at-dosing indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (standard meal or fasted)",
      notes = paste(
        "1 = dose taken with a high-fat meal (healthy-participant phase 1",
        "food study ACH471-016; 16 participants), 0 otherwise; must be 0",
        "whenever FED = 0. Acts on ka (exp(1.09)), D1 (exp(-0.661)), relative",
        "bioavailability (exp(0.101)) and sets the tablet lag time to 0.293 h.",
        "The ka, D1 and lag-time effects were estimated in the food study",
        "and held at those values in the final model (Table 2 footnote)."
      ),
      source_name = "HIFAT (FOOD = 2)"
    ),
    FORM_DANICOPAN_LFC = list(
      description = "Danicopan liquid-filled-capsule formulation indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (approved tablet)",
      notes = paste(
        "1 = dose given as the liquid-filled capsule (LFC) used in the early",
        "phase 1 studies ACH471-001, -002 and -006 (75 participants); 0 =",
        "the approved tablet. Switches ka (0.801 vs 1.58 1/h) and D1 (1.70 vs",
        "3.20 h), multiplies relative bioavailability by 1.29 and removes the",
        "lag time (the lag time applies to tablet doses only). Control stream",
        "LFC = 1 when DFRM = 15."
      ),
      source_name = "LFC (DFRM)"
    ),
    DOSE = list(
      description = "Danicopan dose of the current administration",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Per-administration danicopan dose (not the daily total). Relative",
        "bioavailability scales as (DOSE / 200)^-0.232 (Table 2 'Dose on F';",
        "control stream TVF1 = TTVF * (DOSE/200)**THETA(10)). Doses in the",
        "analysis ranged from 75 to 1200 mg. Use case (a) of the register",
        "entry: in an rxode2 event table put the DOSE column after amt."
      ),
      source_name = "DOSE"
    ),
    REGION_EASTASIA = list(
      description = "East Asian participant indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (outside East Asia)",
      notes = paste(
        "1 = participant from East Asia. Lowers baseline AP activity by 46%",
        "(Table 4 'East Asian on AP0' = -0.460; AP0 equation factor",
        "1 - 0.460 if East Asian; PK-PD control stream column EASTAS). The",
        "Discussion describes it as the East Asian region; the Methods",
        "simulations describe 'participants of East Asian descent'. In the",
        "ALPHA PK-PD population 27 of 69 (39.1%) were Asian (Table 1). No",
        "effect on danicopan PK."
      ),
      source_name = "EASTAS"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Screened on CL/F (the supplement control stream carries",
        "(AGE/33)^THETA(19) with THETA(19) held at 0) and not retained;",
        "the Results report no impact of age on danicopan PK."
      ),
      source_name = "AGE"
    ),
    ALB = list(
      description = "Baseline serum albumin",
      units = "g/dL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Screened on CL/F (control stream (ALBB/4.6)^THETA(20) with",
        "THETA(20) held at 0) and not retained."
      ),
      source_name = "ALBB"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 407L,
    n_studies = 14L,
    n_observations = 7491L,
    age_range = "16.9-82.0 years",
    age_median = "38.0 years",
    weight_range = "43.4-143.0 kg",
    weight_median = "77.0 kg",
    sex_female_pct = 74.2,
    race_ethnicity = c(
      White = 68.3,
      `Black or African American` = 4.9,
      `American Indian or Alaska Native` = 0.5,
      Asian = 17.0,
      `Native Hawaiian or Pacific Islander` = 2.7,
      Other = 2.9,
      Missing = 2.0
    ),
    disease_state = paste(
      "316 healthy participants from 11 phase 1 studies (including renal-",
      "and hepatic-impairment, food-effect, elderly and Japanese studies)",
      "and 91 adults with PNH on background C5 inhibitor therapy from two",
      "phase 2 studies (NCT03181633, NCT03472885) and the phase 3 ALPHA",
      "study (NCT04469465; 69 patients with clinically significant",
      "extravascular hemolysis on ravulizumab or eculizumab)."
    ),
    dose_range = paste(
      "Single doses of 200-1200 mg and multiple doses of 75-800 mg bid or",
      "75-200 mg tid; ALPHA started at 150 mg tid with optional escalation",
      "to 200 mg tid (Table S1)."
    ),
    regions = "Global",
    renal_function = paste(
      "eGFR >= 90 mL/min/1.73 m^2 58.7%, 60-89 35.1%, 30-59 4.4%, 15-29",
      "1.0%, < 15 or dialysis 0.7% (Table 1)."
    ),
    formulation = "Tablet 81.6%, liquid-filled capsule 18.4% (Table 1).",
    food_status = "Fasted 26.3%, standard meal 69.8%, high-fat meal 3.9% (Table 1).",
    pkpd_population = paste(
      "PK-PD (AP activity) model: 69 ALPHA patients with PNH, 475",
      "time-matched AP activity measurements; 59.4% female, median weight",
      "70.0 kg (44.6-143.0), 39.1% Asian (Table 1)."
    ),
    notes = paste(
      "Baseline demographics from Chen 2026 Table 1. 290 BLQ samples",
      "(3.8%) and 6 outliers were excluded from the PK dataset."
    )
  )

  ini({
    # Structural PK -- Chen 2026 Table 2 and the PK equations. The reference
    # individual is a 75 kg male taking a 200 mg tablet with a standard meal.
    lcl <- log(79.5); label("Apparent clearance CL/F (L/h)") # Table 2 CL/F = 79.5 L/h
    lvc <- log(163); label("Apparent central volume Vc/F (L)") # Table 2 Vc/F = 163 L
    lq <- log(34.9); label("Apparent intercompartmental clearance Q/F (L/h)") # Table 2 Q/F = 34.9 L/h
    lvp <- log(232); label("Apparent peripheral volume Vp/F (L)") # Table 2 Vp/F = 232 L

    # Absorption. The tablet and LFC values are separate thetas, mixed by
    # FORM_DANICOPAN_LFC in model() (Ka and D1 equations).
    lka <- log(1.58); label("First-order absorption rate constant, tablet (1/h)") # Table 2 Ka (tablet) = 1.58 1/h
    lka_lfc <- log(0.801); label("First-order absorption rate constant, liquid-filled capsule (1/h)") # Table 2 Ka (LFC) = 0.801 1/h
    ld1 <- log(3.2); label("Zero-order release duration D1, tablet (h)") # Table 2 D1 (tablet) = 3.2 h
    ld1_lfc <- log(1.7); label("Zero-order release duration D1, liquid-filled capsule (h)") # Table 2 D1 (LFC) = 1.7 h
    ltlag <- fixed(log(0.278)); label("Absorption lag time, tablet with standard meal (h)") # Table 2 'Tlag for tablet with standard meal' = 0.278 h, held at the food-study value
    ltlag_highfat <- fixed(log(0.293)); label("Absorption lag time, tablet with high-fat meal (h)") # Table 2 'High-fat meal on Tlag' = 0.293, held at the food-study value; Tlag equation uses it as the tablet + HIFAT lag

    # Relative bioavailability: 1 for the reference tablet, standard meal and
    # 200 mg dose (F equation and Table 2 abbreviations).
    lfdepot <- fixed(log(1)); label("Relative bioavailability, tablet with standard meal at 200 mg (fraction)") # F equation leading factor 1 (reference)
    e_form_danicopan_lfc_fdepot <- 1.29; label("Relative bioavailability of the liquid-filled capsule vs tablet (ratio)") # Table 2 'F of LFC relative to tablet' = 1.29
    e_dose_fdepot <- -0.232; label("Power exponent of (DOSE/200) on relative bioavailability (unitless)") # Table 2 'Dose on F' = -0.232

    # Food effects (fasted = 1 - FED; high-fat = FED_HIGHFAT).
    e_fasted_fdepot <- -0.143; label("Fasted effect on relative bioavailability, exponential coefficient (unitless)") # Table 2 'Fasted on F' = -0.143
    e_fed_highfat_fdepot <- 0.101; label("High-fat-meal effect on relative bioavailability, exponential coefficient (unitless)") # Table 2 'High-fat meal on F' = 0.101
    e_fasted_d1 <- 0.0787; label("Fasted effect on D1, exponential coefficient (unitless)") # Table 2 'Fasted on D1' = 0.0787
    e_fed_highfat_d1 <- fixed(-0.661); label("High-fat-meal effect on D1, exponential coefficient (unitless)") # Table 2 'High-fat meal on D1' = -0.661, held at the food-study value
    e_fasted_ka <- 30.9; label("Fasted additive increment in ka (1/h)") # Table 2 'Fasted on Ka' = 30.9; Ka equation '+ 30.9 if FAST'
    e_fed_highfat_ka <- fixed(1.09); label("High-fat-meal effect on ka, exponential coefficient (unitless)") # Table 2 'High-fat meal on Ka' = 1.09, held at the food-study value; Ka equation 'x e^(1.09 if HIFAT)'

    # Covariates on disposition.
    e_wt_cl_q <- 0.645; label("Body-weight exponent on CL/F and Q/F (unitless)") # Table 2 'Weight exponent on CL/F and Q/F' = 0.645
    e_wt_vc_vp <- 0.876; label("Body-weight exponent on Vc/F and Vp/F (unitless)") # Table 2 'Weight exponent on Vc/F and Vp/F' = 0.876
    e_sexf_cl <- -0.221; label("Female-sex fractional change in CL/F (fraction)") # Table 2 'Female sex on CL/F' = -0.221
    e_renalimp_sev_cl <- -0.433; label("Severe-renal-impairment fractional change in CL/F (fraction)") # Table 2 'RI on CL/F' = -0.433

    # IIV -- Table 2 variances (BSV = 100 * sqrt(exp(omega^2) - 1)).
    etalcl + etalvc ~ c(0.0634, 0.0774, 0.412) # Table 2 omega2 CL/F 0.0634 (25.6%), covariance CL/F x Vc/F 0.0774, omega2 Vc/F 0.412 (71.4%)
    etalvp ~ 0.141 # Table 2 omega2 Vp/F 0.141 (38.9%)
    etalka ~ 2.1 # Table 2 omega2 Ka 2.1 (268%)

    # Residual error, combined: W = sqrt(F^2 * RUVCV^2 + RUVSD^2) in the
    # supplement $ERROR block.
    propSd <- 0.435; label("Proportional residual error (fraction)") # Table 2 Proportional = 0.435 (control stream RUVCV, initial 0.445)
    addSd <- 0.254; label("Additive residual error (ng/mL)") # Table 2 Additive = 0.254 ng/mL

    # PK-PD (AP activity by hemolysis assay) -- Chen 2026 Table 4 and the
    # AP0 / AP equations.
    e0 <- 29.1; label("Baseline AP activity for a 69 kg non-East-Asian patient (%)") # Table 4 AP0 = 29.1 %
    lic50 <- log(12); label("Danicopan concentration giving half-maximal AP inhibition (ng/mL)") # Table 4 IC50 = 12.0 ng/mL
    lrmax_alternative_pathway_activity <- log(3.05); label("AP activity reached at full inhibition (%)") # Table 4 Imax = 3.05 %; AP equation AP0 - (AP0 - 3.05) * C / (12 + C)
    e_wt_e0 <- 1.29; label("Body-weight exponent on baseline AP activity (unitless)") # Table 4 'Weight on AP0' = 1.29
    e_region_eastasia_e0 <- -0.46; label("East-Asian fractional change in baseline AP activity (fraction)") # Table 4 'East Asian on AP0' = -0.460

    etae0 ~ 91.1 # Table 4 omega2 BASE 91.1, additive IIV (SD 9.5 AP percentage units)
    etalic50 ~ 5.26 # Table 4 omega2 IC50 5.26, exponential IIV

    propSd_alternative_pathway_activity <- 0.671; label("Proportional residual error on AP activity (fraction)") # Table 4 Proportional = 67.1 %
  })

  model({
    # Indicator for the paper's FAST covariate (fasted = not fed).
    fasted <- 1 - FED
    tablet <- 1 - FORM_DANICOPAN_LFC

    # Disposition (PK equations).
    cl <- exp(lcl + etalcl) * (WT / 75)^e_wt_cl_q *
      (1 + e_sexf_cl * SEXF) * (1 + e_renalimp_sev_cl * RENALIMP_SEV)
    vc <- exp(lvc + etalvc) * (WT / 75)^e_wt_vc_vp
    q <- exp(lq) * (WT / 75)^e_wt_cl_q
    vp <- exp(lvp + etalvp) * (WT / 75)^e_wt_vc_vp

    # Absorption (Ka, D1, F and Tlag equations).
    ka <- (exp(lka) * tablet + exp(lka_lfc) * FORM_DANICOPAN_LFC + e_fasted_ka * fasted) *
      exp(e_fed_highfat_ka * FED_HIGHFAT) * exp(etalka)
    d1 <- (exp(ld1) * tablet + exp(ld1_lfc) * FORM_DANICOPAN_LFC) *
      exp(e_fasted_d1 * fasted) * exp(e_fed_highfat_d1 * FED_HIGHFAT)
    fdepot <- exp(lfdepot) * e_form_danicopan_lfc_fdepot^FORM_DANICOPAN_LFC *
      (DOSE / 200)^e_dose_fdepot *
      exp(e_fasted_fdepot * fasted) * exp(e_fed_highfat_fdepot * FED_HIGHFAT)
    tlag <- exp(ltlag) * tablet * FED * (1 - FED_HIGHFAT) +
      exp(ltlag_highfat) * tablet * FED_HIGHFAT

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Zero-order release into the depot: dose records need rate = -2.
    dur(depot) <- d1
    alag(depot) <- tlag
    f(depot) <- fdepot

    # Dose in mg, volume in L -> mg/L; x 1000 gives ng/mL (control stream
    # S2 = V2/1000).
    Cc <- central / vc * 1000

    # Direct inhibitory Emax model of AP activity (%), with an additive
    # random effect on the covariate-scaled baseline.
    ap_base <- e0 * (WT / 69)^e_wt_e0 * (1 + e_region_eastasia_e0 * REGION_EASTASIA) + etae0
    ic50 <- exp(lic50 + etalic50)
    rmax_alternative_pathway_activity <- exp(lrmax_alternative_pathway_activity)
    alternative_pathway_activity <- ap_base -
      (ap_base - rmax_alternative_pathway_activity) * Cc / (ic50 + Cc)

    Cc ~ add(addSd) + prop(propSd)
    alternative_pathway_activity ~ prop(propSd_alternative_pathway_activity)
  })
}
