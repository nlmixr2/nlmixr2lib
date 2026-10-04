Wojciechowski_2022_abrocitinib <- function() {
  description <- paste(
    "Two-compartment population PK model for oral abrocitinib (JAK1 inhibitor)",
    "in healthy adults and in adolescents and adults with psoriasis or atopic",
    "dermatitis (Wojciechowski 2022). Parallel first-order (depot) and",
    "zero-order (central) absorption, with the first-order arm capped at a",
    "fixed amount per dose; an absorption lag for tablets; absolute oral",
    "bioavailability on the logit scale anchored at 0.5977 with additive",
    "shifts for formulation, repeated dosing, CYP DDIs, race, disease, hepatic",
    "impairment, adolescence, sex and the 800 mg dose; clearance that falls",
    "exponentially with time on treatment and depends on the effective",
    "(bioavailable) daily dose; allometric weight scaling; and a study-specific",
    "proportional residual error."
  )
  reference <- paste(
    "Wojciechowski J, Malhotra BK, Wang X, Fostvedt L, Valdez H, Nicholas T.",
    "Population Pharmacokinetics of Abrocitinib in Healthy Individuals and",
    "Patients with Psoriasis or Atopic Dermatitis.",
    "Clin Pharmacokinet. 2022;61:709-723. doi:10.1007/s40262-021-01104-z.",
    "Model equations from the Online Resource 2 NONMEM control stream.",
    sep = " "
  )
  vignette <- "Wojciechowski_2022_abrocitinib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power scaling referenced to 70 kg: exponent 0.453 on CL and Q, 0.52 on Vc and Vp (Table 2). Online Resource 2 uses baseline weight (BWT).",
      source_name = "BWT"
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Additive shift of +0.353 on logit F for females. Source column SEX is 1 = male, 2 = female; SEXF = SEX - 1.",
      source_name = "SEX"
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (White, Black or unknown race)",
      notes = "Shares one additive shift of +0.815 on logit F with RACE_OTHER (source RACE = 3 or 4). Japanese participants are Asian here; a separate Japanese effect was tested and not retained.",
      source_name = "RACE"
    ),
    RACE_OTHER = list(
      description = "Other race indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (White, Black or unknown race)",
      notes = "Shares the +0.815 logit-F shift with RACE_ASIAN (source RACE = 4). Mutually exclusive with RACE_ASIAN.",
      source_name = "RACE"
    ),
    ADOLESCENT = list(
      description = "Adolescent age-group indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (adult)",
      notes = "Source AGEGRP = 1 (adolescents, 12 to 17 years in the phase III atopic-dermatitis studies). Additive shift of -0.589 on logit F.",
      source_name = "AGEGRP"
    ),
    DIS_PSORIASIS = list(
      description = "Plaque psoriasis patient indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (healthy volunteer)",
      notes = "Source PTST = 1. Shares one additive shift of +0.489 on logit F with DIS_ATOPIC_DERMATITIS; the paper could not distinguish the two diseases. Mutually exclusive with DIS_ATOPIC_DERMATITIS, HEPIMP_MILD and HEPIMP_MOD.",
      source_name = "PTST"
    ),
    DIS_ATOPIC_DERMATITIS = list(
      description = "Atopic dermatitis patient indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (healthy volunteer)",
      notes = "Source PTST = 2. Shares the +0.489 logit-F shift with DIS_PSORIASIS.",
      source_name = "PTST"
    ),
    HEPIMP_MILD = list(
      description = "Mild hepatic impairment indicator (Child-Pugh A, score 5-6)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (normal hepatic function)",
      notes = "Source PTST = 3. Shares one additive shift of +1.3 on logit F with HEPIMP_MOD (Child-Pugh classification).",
      source_name = "PTST"
    ),
    HEPIMP_MOD = list(
      description = "Moderate hepatic impairment indicator (Child-Pugh B, score 7-9)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (normal hepatic function)",
      notes = "Source PTST = 4. Shares the +1.3 logit-F shift with HEPIMP_MILD.",
      source_name = "PTST"
    ),
    CONMED_RIFAMPICIN = list(
      description = "Concomitant rifampicin (rifampin) indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no rifampicin)",
      notes = "Rifampin 600 mg once daily for 8 days in B7451019. CL multiplied by (1 + 0.264) and an additive shift of -2.08 on logit F (CYP2C19 / CYP3A4 / CYP2C9 induction).",
      source_name = "RIFDDI"
    ),
    CONMED_FLUCONAZOLE = list(
      description = "Concomitant fluconazole indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no fluconazole)",
      notes = "Fluconazole 400 mg then 200 mg once daily in B7451017. CL multiplied by (1 - 0.541); shares the +1.31 logit-F shift with CONMED_FLUVOXAMINE (the shift applies once if either is present).",
      source_name = "FLZDDI"
    ),
    CONMED_FLUVOXAMINE = list(
      description = "Concomitant fluvoxamine indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no fluvoxamine)",
      notes = "Fluvoxamine 50 mg once daily in B7451017. CL multiplied by (1 - 0.234); shares the +1.31 logit-F shift with CONMED_FLUCONAZOLE.",
      source_name = "FLUDDI"
    ),
    FED_HIGHFAT = list(
      description = "High-fat meal at dosing indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (fasted, or food not controlled)",
      notes = "Source FOOD = 1. Multiplies the first-order amount cap by (1 - 1) = 0, so the whole bioavailable dose is absorbed zero-order. Source FOOD = 2 (not controlled, phase III) is grouped with fasted.",
      source_name = "FOOD"
    ),
    FORM_ABROCITINIB_SUSP = list(
      description = "Abrocitinib oral suspension formulation indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (phase II 100 mg tablet when all three FORM_ABROCITINIB_* indicators are 0)",
      notes = "Source FORMS = 1 (suspension in 0.5% methylcellulose). Multiplies the first-order amount cap by (1 + 1.17); no absorption lag (the 0.183 h lag applies to tablets only, source FORM = 2).",
      source_name = "FORMS"
    ),
    FORM_ABROCITINIB_TAB_PH2B = list(
      description = "Abrocitinib phase IIb 10 mg and 50 mg tablet formulation indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (phase II 100 mg tablet when all three FORM_ABROCITINIB_* indicators are 0)",
      notes = "Source FORMS = 3 (tablets D1500091 / D1500093, phase IIb atopic-dermatitis study B7451006). Additive shift of -1.02 on logit F and the first-order amount cap multiplied by (1 - 0.68).",
      source_name = "FORMS"
    ),
    FORM_ABROCITINIB_TAB_PH3 = list(
      description = "Abrocitinib phase III 100 mg film-coated tablet formulation indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (phase II 100 mg tablet when all three FORM_ABROCITINIB_* indicators are 0)",
      notes = "Source FORMS = 4 (tablet D1700147, the commercial formulation and the paper's reference scenario). Additive shift of -0.766 on logit F.",
      source_name = "FORMS"
    ),
    MULTI_DOSE = list(
      description = "Repeated-dosing indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (single dose, or the first dose of a regimen)",
      notes = "Source MULTI (0 = single dose, 1 = repeated dose). Additive shift of +0.241 on logit F. Because logit F also feeds the effective daily dose on CL, a single-dose subject stays 0 after the dose.",
      source_name = "MULTI"
    ),
    DOSE_ABROCITINIB_MG = list(
      description = "Abrocitinib amount per administration",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "Source DOSE. Sets the first-order fraction min(1, amt_fo / DOSE) and switches on the -0.778 logit-F shift at exactly 800 mg (source DOSE.EQ.800). Must equal the amt of the dose records.",
      source_name = "DOSE"
    ),
    DOSE_ABROCITINIB_MGD = list(
      description = "Randomized total daily abrocitinib dose",
      units = "mg/day",
      type = "continuous",
      reference_category = NULL,
      notes = "Source DOSR. CL scales as (F * DOSR / 200)^-0.169, where F is the individual absolute bioavailability, referenced to 200 mg once daily. Equals the per-administration dose for once-daily regimens and twice it for twice-daily regimens.",
      source_name = "DOSR"
    ),
    STUDY_B7451005 = list(
      description = "Phase II psoriasis study B7451005 (NCT02201524) indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (other studies)",
      notes = "Moderate-variability study: proportional residual SD multiplied by (1 + 0.495).",
      source_name = "PROT"
    ),
    STUDY_B7451006 = list(
      description = "Phase IIb atopic dermatitis study B7451006 (NCT02780167) indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (other studies)",
      notes = "High-variability study: proportional residual SD multiplied by (1 + 1.16).",
      source_name = "PROT"
    ),
    STUDY_B7451012 = list(
      description = "Phase III atopic dermatitis study B7451012 (JADE MONO-1, NCT03349060) indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (other studies)",
      notes = "Moderate-variability study: proportional residual SD multiplied by (1 + 0.495).",
      source_name = "PROT"
    ),
    STUDY_B7451043 = list(
      description = "Phase I probenecid DDI study B7451043 (NCT03937258) indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (other studies)",
      notes = "Moderate-variability study: proportional residual SD multiplied by (1 + 0.495).",
      source_name = "PROT"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "abrocitinib", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "abrocitinib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "abrocitinib", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 995L,
    n_studies = 11L,
    n_observations = 6206L,
    age_range = "12-84 years",
    age_median = "34 years",
    weight_range = "34.4-180 kg",
    weight_median = "76 kg",
    sex_female_pct = 39.2,
    race_ethnicity = c(White = 66.6, Black = 11.0, Asian = 17.9, Other = 4.0, Unknown = 0.5),
    disease_state = "Healthy volunteers (n = 165), moderate-to-severe plaque psoriasis (n = 45), moderate-to-severe atopic dermatitis (n = 769, including 90 adolescents), mild (n = 8) and moderate (n = 8) hepatic impairment",
    dose_range = "3-800 mg/day oral (single doses 3-800 mg; 10-400 mg once daily or 100-200 mg twice daily)",
    regions = "Multinational (Western and Japanese participants)",
    notes = "Seven phase I, two phase II and two phase III studies (Online Resource 1). Demographics from Table 1. Most concentrations came from 200 mg (45.3%) and 100 mg (29.8%) doses."
  )

  ini({
    # Disposition, absolute (not apparent) because F is anchored at its
    # measured absolute value. Reference: 70 kg.
    lcl <- log(22); label("Clearance (L/h)") # Table 2 'CL, L/h' = 22
    lvc <- log(87.8); label("Central volume of distribution (L)") # Table 2 'Vc, L' = 87.8
    lq <- log(1.16); label("Intercompartmental clearance (L/h)") # Table 2 'Q, L/h' = 1.16
    lvp <- log(8.25); label("Peripheral volume of distribution (L)") # Table 2 'Vp, L' = 8.25

    # Absorption
    lr1 <- log(75.3); label("Zero-order input rate into central (mg/h)") # Table 2 'k0, mg/h' = 75.3
    lamt_fo <- log(121); label("Maximum amount of a dose absorbed by the first-order process (mg)") # Table 2 'AK1, mg' = 121
    lka <- log(4.01); label("First-order absorption rate constant (1/h)") # Table 2 'ka, h-1' = 4.01
    ltlag <- log(0.183); label("Absorption lag time for tablet formulations (h)") # Table 2 'Effect of tablet formulations on ALAG1' = 0.183; base ALAG1 = 0 (Online Resource 2)

    # Absolute oral bioavailability, logit scale. Anchored at the 14C
    # microtracer value for the oral solution, assumed equal for the suspension.
    logitfdepot <- fixed(logit(0.5977)); label("Logit of absolute oral bioavailability, suspension and phase II 100 mg tablet (logit)") # Methods 2.4.2 and Online Resource 2 'POPF = 0.5977'

    # Covariate effects on logit F (additive on the logit scale)
    e_rif_f <- -2.08; label("Rifampicin shift on logit F (logit)") # Table 2 'Effect of rifampin on F' = -2.08
    e_cypinh_f <- 1.31; label("Fluconazole or fluvoxamine shift on logit F (logit)") # Table 2 'Effect of fluconazole or fluvoxamine on F' = 1.31
    e_ph2b_f <- -1.02; label("Phase IIb tablet shift on logit F (logit)") # Table 2 'Effect of phase IIb tablets on F' = -1.02
    e_ph3_f <- -0.766; label("Phase III tablet shift on logit F (logit)") # Table 2 'Effect of phase III tablet on F' = -0.766
    e_multidose_f <- 0.241; label("Repeated-dosing shift on logit F (logit)") # Table 2 'Effect of multiple dosing on F' = 0.241
    e_race_f <- 0.815; label("Asian or other race shift on logit F (logit)") # Table 2 'Effect of Asian/other race on F' = 0.815
    e_dis_f <- 0.489; label("Psoriasis or atopic dermatitis shift on logit F (logit)") # Table 2 'Combined effect of psoriasis and AD on F' = 0.489
    e_hepimp_f <- 1.3; label("Mild or moderate hepatic impairment shift on logit F (logit)") # Table 2 'Combined effect of mild and moderate hepatic impairment on F' = 1.3
    e_adol_f <- -0.589; label("Adolescent shift on logit F (logit)") # Table 2 'Effect of adolescent age on F' = -0.589
    e_dose800_f <- -0.778; label("800 mg dose shift on logit F (logit)") # Table 2 'Effect of abrocitinib 800-mg dose on F' = -0.778
    e_sexf_f <- 0.353; label("Female sex shift on logit F (logit)") # Table 2 'Effect of female sex on F' = 0.353

    # Covariate effects on the first-order amount cap (proportional)
    e_highfat_amt_fo <- fixed(-1); label("High-fat meal fractional change in the first-order amount cap (fraction)") # Table 2 'Effect of high-fat meal on Ak1' = -1, held constant (Online Resource 3 run 133)
    e_susp_amt_fo <- 1.17; label("Suspension fractional change in the first-order amount cap (fraction)") # Table 2 'Effect of suspension on Ak1' = 1.17
    e_ph2b_amt_fo <- -0.68; label("Phase IIb tablet fractional change in the first-order amount cap (fraction)") # Table 2 'Effect of phase IIb tablets on Ak1' = -0.68

    # Covariate effects on CL
    e_rif_cl <- 0.264; label("Rifampicin fractional change in CL (fraction)") # Table 2 'Effect of rifampin on CL' = 0.264
    e_fluco_cl <- -0.541; label("Fluconazole fractional change in CL (fraction)") # Table 2 'Effect of fluconazole on CL' = -0.541
    e_fluvox_cl <- -0.234; label("Fluvoxamine fractional change in CL (fraction)") # Table 2 'Effect of fluvoxamine on CL' = -0.234
    e_dose_cl <- -0.169; label("Power exponent of effective daily dose / 200 mg on CL (unitless)") # Table 2 'Effect of effective daily dose on CL' = -0.169

    # Time-dependent CL: cl * (1 + cl_exp_famp * (1 - exp(-cl_exp_kdes * t)))
    cl_exp_famp <- -0.186; label("Maximum fractional change in CL with time on treatment (fraction)") # Table 2 'Maximum change in CL with respect to time (TAFO)' = -0.186 (fraction, printed with a percent unit)
    lcl_exp_kdes <- log(log(2) / 21.6); label("Rate constant of the time-dependent change in CL (1/h)") # Table 2 'Rate of change in CL with respect to time (half-life), h' = 21.6; rate = ln(2) / 21.6

    # Allometric exponents (reference 70 kg)
    e_wt_cl <- 0.453; label("Power exponent of WT / 70 kg on CL and Q (unitless)") # Table 2 'Effect of weight on CL and Q' = 0.453
    e_wt_vc <- 0.52; label("Power exponent of WT / 70 kg on Vc and Vp (unitless)") # Table 2 'Effect of weight on Vc and Vp' = 0.52

    # IIV: Table 2 CV = sqrt(omega^2) x 100 (footnote a), so omega = CV / 100
    etalcl + etalvc ~ c(0.332929, 0.0778712, 0.171396) # Table 2 omega CL 57.7 and omega Vc 41.4 (CV = sqrt(omega^2)); rho CL-Vc 0.326 so cov = 0.326 * 0.577 * 0.414

    # Residual error: W = sqrt(IPRED^2 * PRO^2 + ADD^2) with EPS variance 1 held constant
    propSd <- 0.437; label("Proportional residual SD, reference studies (fraction)") # Table 2 'RUV PRO, SD' = 0.437
    addSd <- 0.509; label("Additive residual SD (ng/mL)") # Table 2 'RUV ADD, SD' = 0.509
    e_studymod_propSd <- 0.495; label("Moderate-variability study fractional change in proportional SD (fraction)") # Table 2 'Effect of moderate-variability studies on RUV PRO' = 0.495
    e_studyhi_propSd <- 1.16; label("High-variability study fractional change in proportional SD (fraction)") # Table 2 'Effect of high-variability studies on RUV PRO' = 1.16
  })

  model({
    # Logit-scale absolute bioavailability (Online Resource 2, FT / FI)
    cypinh <- CONMED_FLUCONAZOLE + CONMED_FLUVOXAMINE
    if (cypinh > 1) cypinh <- 1
    dose800 <- 0
    if (DOSE_ABROCITINIB_MG == 800) dose800 <- 1
    logit_fi <- logitfdepot +
      e_rif_f * CONMED_RIFAMPICIN +
      e_cypinh_f * cypinh +
      e_ph2b_f * FORM_ABROCITINIB_TAB_PH2B +
      e_ph3_f * FORM_ABROCITINIB_TAB_PH3 +
      e_multidose_f * MULTI_DOSE +
      e_race_f * (RACE_ASIAN + RACE_OTHER) +
      e_dis_f * (DIS_PSORIASIS + DIS_ATOPIC_DERMATITIS) +
      e_hepimp_f * (HEPIMP_MILD + HEPIMP_MOD) +
      e_adol_f * ADOLESCENT +
      e_dose800_f * dose800 +
      e_sexf_f * SEXF
    fi <- expit(logit_fi)

    # Split of the bioavailable dose: up to amt_fo mg goes first-order into
    # the depot, the remainder zero-order into central at rate r1.
    amt_fo <- exp(lamt_fo) *
      (1 + e_highfat_amt_fo * FED_HIGHFAT) *
      (1 + e_susp_amt_fo * FORM_ABROCITINIB_SUSP + e_ph2b_amt_fo * FORM_ABROCITINIB_TAB_PH2B)
    ffo <- amt_fo / DOSE_ABROCITINIB_MG
    if (ffo > 1) ffo <- 1

    # Clearance: weight, DDIs, time on treatment and effective daily dose
    cl_exp_kdes <- exp(lcl_exp_kdes)
    cl_time <- 1 + cl_exp_famp * (1 - exp(-cl_exp_kdes * t))
    cl_dose <- (fi * DOSE_ABROCITINIB_MGD / 200)^e_dose_cl
    cl_ddi <- (1 + e_rif_cl * CONMED_RIFAMPICIN) *
      (1 + e_fluco_cl * CONMED_FLUCONAZOLE) *
      (1 + e_fluvox_cl * CONMED_FLUVOXAMINE)
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl * cl_ddi * cl_time * cl_dose
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    q <- exp(lq) * (WT / 70)^e_wt_cl
    vp <- exp(lvp) * (WT / 70)^e_wt_vc
    ka <- exp(lka)
    r1 <- exp(lr1)
    tlag <- exp(ltlag) * (1 - FORM_ABROCITINIB_SUSP)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Each administration is two dose records with the same amt: one into
    # depot (rate 0) and one into central with rate = -1, so rxode2 applies
    # rate(central) and the infusion lasts fi * (1 - ffo) * amt / r1.
    f(depot) <- fi * ffo
    f(central) <- fi * (1 - ffo)
    rate(central) <- r1
    alag(depot) <- tlag
    # The central record carries no drug when the whole dose fits under
    # amt_fo (ffo = 1, e.g. any tablet dose <= 121 mg). Its lag is then
    # meaningless, and a lagged zero-amount rate = -1 record makes rxode2
    # 5.1.8 fail ('Rate is zero/negative') from the third dose on when a
    # covariate such as MULTI_DOSE changes during the regimen.
    alag(central) <- tlag * (ffo < 1)

    # Study-specific proportional residual SD
    propSd_study <- propSd *
      (1 + e_studymod_propSd * (STUDY_B7451005 + STUDY_B7451012 + STUDY_B7451043) +
        e_studyhi_propSd * STUDY_B7451006)

    # mg / L x 1000 = ng/mL
    Cc <- central / vc * 1000
    Cc ~ add(addSd) + prop(propSd_study)
  })
}
