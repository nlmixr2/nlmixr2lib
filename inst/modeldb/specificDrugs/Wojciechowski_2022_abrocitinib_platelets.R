Wojciechowski_2022_abrocitinib_platelets <- function() {
  description <- paste(
    "Population PK-PD model of the platelet time-course after oral abrocitinib",
    "(JAK1 inhibitor) in adolescents and adults with atopic dermatitis or",
    "psoriasis (Wojciechowski 2022, Br J Clin Pharmacol). A Friberg",
    "semi-mechanistic myelosuppression model (proliferating pool, three transit",
    "compartments and circulating platelets, mean transit time 7.21 days,",
    "feedback exponent 0.232) in which abrocitinib concentration inhibits",
    "platelet proliferation by an Emax function. Emax scales with the",
    "individual baseline platelet count (power -1.11) and wanes with time on",
    "treatment (tachyphylaxis half-life 151 days). Baseline platelet count",
    "depends on baseline white blood cell count, haematocrit, age, sex and",
    "psoriasis. The abrocitinib concentrations come from the companion",
    "Wojciechowski 2022 population PK model, carried here unchanged."
  )
  reference <- paste(
    "Wojciechowski J, Malhotra BK, Wang X, Fostvedt L, Valdez H, Nicholas T.",
    "Population pharmacokinetic-pharmacodynamic modelling of platelet",
    "time-courses following administration of abrocitinib.",
    "Br J Clin Pharmacol. 2022;88(8):3856-3871. doi:10.1111/bcp.15334.",
    "PK structure from Wojciechowski J et al. Population Pharmacokinetics of",
    "Abrocitinib in Healthy Individuals and Patients with Psoriasis or Atopic",
    "Dermatitis. Clin Pharmacokinet. 2022;61:709-723.",
    "doi:10.1007/s40262-021-01104-z; see",
    "modellib('Wojciechowski_2022_abrocitinib').",
    sep = " "
  )
  vignette <- "Wojciechowski_2022_abrocitinib_platelets"
  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "ng/mL",
    PLT = "10^9/L (identical to 1000/uL)"
  )

  covariateData <- list(
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power scaling referenced to 70 kg: exponent 0.453 on CL and Q, 0.52 on Vc and Vp (Clin Pharmacokinet Table 2). Clin Pharmacokinet Online Resource 2 uses baseline weight (BWT).",
      source_name = "BWT"
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "PK: additive shift of +0.353 on logit F for females. PD: baseline platelet count multiplied by (1 + 0.0584) for females (Br J Clin Pharmacol Equation 9 and Table 2). Source column SEX is 1 = male, 2 = female; SEXF = SEX - 1.",
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
      reference_category = "0 (healthy volunteer for the PK shift; atopic dermatitis patient for the PD shift)",
      notes = "Source PTST = 1. PK: shares one additive shift of +0.489 on logit F with DIS_ATOPIC_DERMATITIS; the PK paper could not distinguish the two diseases. PD: baseline platelet count multiplied by (1 - 0.118) for psoriasis relative to atopic dermatitis (Br J Clin Pharmacol Table 2); every PD subject had one of the two diseases. Mutually exclusive with DIS_ATOPIC_DERMATITIS, HEPIMP_MILD and HEPIMP_MOD.",
      source_name = "PTST"
    ),
    DIS_ATOPIC_DERMATITIS = list(
      description = "Atopic dermatitis patient indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (healthy volunteer)",
      notes = "Source PTST = 2. Shares the +0.489 logit-F shift with DIS_PSORIASIS. Set to 1 for every atopic dermatitis patient, including when simulating the platelet model, so that the PK reflects patient bioavailability.",
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
    ),
    AGE = list(
      description = "Age at baseline",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on baseline platelet count referenced to the 30-year median, (AGE / 30)^-0.0695 (Br J Clin Pharmacol Table 2 and Equation 10). Enters the PD model only; the PK model uses the ADOLESCENT indicator instead.",
      source_name = "AGE"
    ),
    WBC = list(
      description = "Baseline white blood cell count",
      units = "10^9/L (1000/uL)",
      type = "continuous",
      reference_category = NULL,
      notes = "Baseline value, time-fixed. Power effect on baseline platelet count referenced to 6.9 x 1000/uL, (WBC / 6.9)^0.199 (Br J Clin Pharmacol Table 2).",
      source_name = "WBC"
    ),
    HCT = list(
      description = "Baseline haematocrit",
      units = "%",
      type = "continuous",
      reference_category = NULL,
      notes = "Baseline value, time-fixed, percent scale. Power effect on baseline platelet count referenced to 43%, (HCT / 43)^-0.203 (Br J Clin Pharmacol Table 2).",
      source_name = "HCT"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "abrocitinib", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "abrocitinib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "abrocitinib", units = "mg", specimen = "plasma", verified = TRUE),
    prol = list(analyte = "platelet progenitors", units = "10^9/L", specimen = "not applicable", verified = TRUE),
    transit1 = list(
      analyte = "maturing platelet precursors",
      units = "10^9/L",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit2 = list(
      analyte = "maturing platelet precursors",
      units = "10^9/L",
      specimen = "not applicable",
      verified = TRUE
    ),
    transit3 = list(
      analyte = "maturing platelet precursors",
      units = "10^9/L",
      specimen = "not applicable",
      verified = TRUE
    ),
    circ = list(analyte = "platelets", units = "10^9/L", specimen = "whole blood", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 2333L,
    n_studies = 5L,
    n_observations = 12189L,
    age_range = "12-84 years",
    age_median = "30 years",
    weight_range = "34-204 kg",
    weight_median = "73 kg",
    sex_female_pct = 44.4,
    race_ethnicity = c(White = 71.8, Asian = 17.8, Black = 7.6, Other = 2.3, Unknown = 0.6),
    disease_state = "Moderate-to-severe atopic dermatitis (n = 2274, including 340 adolescents aged 12-17 years) and moderate-to-severe plaque psoriasis (n = 59)",
    dose_range = "Placebo (n = 224) or oral abrocitinib 10, 30, 100, 200 or 400 mg once daily or 200 mg twice daily for 4-12 weeks",
    regions = "Multinational (includes 44 Japanese patients)",
    baseline_platelets = "Median 274 x 1000/uL (range 118-651)",
    notes = "Two phase II studies (B7451005 psoriasis, B7451006 atopic dermatitis) and three phase III atopic dermatitis studies (B7451012 JADE MONO-1, B7451013 JADE MONO-2, B7451014 JADE REGIMEN; Table A1). Platelet counts were sampled up to week 16. Demographics from Table 1. PK parameters were the empirical Bayes estimates of the companion population PK model."
  )

  ini({
    # ------------------------------------------------------------------
    # Abrocitinib PK, carried unchanged from modellib('Wojciechowski_2022_abrocitinib').
    # 'CPK' in the comments of this block is the Clin Pharmacokinet 2022
    # popPK paper; its Table 2 and Online Resource 2 are the sources.
    # ------------------------------------------------------------------
    # Disposition, absolute (not apparent) because F is anchored at its
    # measured absolute value. Reference: 70 kg.
    lcl <- log(22); label("Clearance (L/h)") # CPK Table 2 'CL, L/h' = 22
    lvc <- log(87.8); label("Central volume of distribution (L)") # CPK Table 2 'Vc, L' = 87.8
    lq <- log(1.16); label("Intercompartmental clearance (L/h)") # CPK Table 2 'Q, L/h' = 1.16
    lvp <- log(8.25); label("Peripheral volume of distribution (L)") # CPK Table 2 'Vp, L' = 8.25

    # Absorption
    lr1 <- log(75.3); label("Zero-order input rate into central (mg/h)") # CPK Table 2 'k0, mg/h' = 75.3
    lamt_fo <- log(121); label("Maximum amount of a dose absorbed by the first-order process (mg)") # CPK Table 2 'AK1, mg' = 121
    lka <- log(4.01); label("First-order absorption rate constant (1/h)") # CPK Table 2 'ka, h-1' = 4.01
    ltlag <- log(0.183); label("Absorption lag time for tablet formulations (h)") # CPK Table 2 'Effect of tablet formulations on ALAG1' = 0.183; base ALAG1 = 0 (CPK Online Resource 2)

    # Absolute oral bioavailability, logit scale. Anchored at the 14C
    # microtracer value for the oral solution, assumed equal for the suspension.
    logitfdepot <- fixed(logit(0.5977)); label("Logit of absolute oral bioavailability, suspension and phase II 100 mg tablet (logit)") # CPK Methods 2.4.2 and Online Resource 2 'POPF = 0.5977'

    # Covariate effects on logit F (additive on the logit scale)
    e_rif_f <- -2.08; label("Rifampicin shift on logit F (logit)") # CPK Table 2 'Effect of rifampin on F' = -2.08
    e_cypinh_f <- 1.31; label("Fluconazole or fluvoxamine shift on logit F (logit)") # CPK Table 2 'Effect of fluconazole or fluvoxamine on F' = 1.31
    e_ph2b_f <- -1.02; label("Phase IIb tablet shift on logit F (logit)") # CPK Table 2 'Effect of phase IIb tablets on F' = -1.02
    e_ph3_f <- -0.766; label("Phase III tablet shift on logit F (logit)") # CPK Table 2 'Effect of phase III tablet on F' = -0.766
    e_multidose_f <- 0.241; label("Repeated-dosing shift on logit F (logit)") # CPK Table 2 'Effect of multiple dosing on F' = 0.241
    e_race_f <- 0.815; label("Asian or other race shift on logit F (logit)") # CPK Table 2 'Effect of Asian/other race on F' = 0.815
    e_dis_f <- 0.489; label("Psoriasis or atopic dermatitis shift on logit F (logit)") # CPK Table 2 'Combined effect of psoriasis and AD on F' = 0.489
    e_hepimp_f <- 1.3; label("Mild or moderate hepatic impairment shift on logit F (logit)") # CPK Table 2 'Combined effect of mild and moderate hepatic impairment on F' = 1.3
    e_adol_f <- -0.589; label("Adolescent shift on logit F (logit)") # CPK Table 2 'Effect of adolescent age on F' = -0.589
    e_dose800_f <- -0.778; label("800 mg dose shift on logit F (logit)") # CPK Table 2 'Effect of abrocitinib 800-mg dose on F' = -0.778
    e_sexf_f <- 0.353; label("Female sex shift on logit F (logit)") # CPK Table 2 'Effect of female sex on F' = 0.353

    # Covariate effects on the first-order amount cap (proportional)
    e_highfat_amt_fo <- fixed(-1); label("High-fat meal fractional change in the first-order amount cap (fraction)") # CPK Table 2 'Effect of high-fat meal on Ak1' = -1, held constant (Online Resource 3 run 133)
    e_susp_amt_fo <- 1.17; label("Suspension fractional change in the first-order amount cap (fraction)") # CPK Table 2 'Effect of suspension on Ak1' = 1.17
    e_ph2b_amt_fo <- -0.68; label("Phase IIb tablet fractional change in the first-order amount cap (fraction)") # CPK Table 2 'Effect of phase IIb tablets on Ak1' = -0.68

    # Covariate effects on CL
    e_rif_cl <- 0.264; label("Rifampicin fractional change in CL (fraction)") # CPK Table 2 'Effect of rifampin on CL' = 0.264
    e_fluco_cl <- -0.541; label("Fluconazole fractional change in CL (fraction)") # CPK Table 2 'Effect of fluconazole on CL' = -0.541
    e_fluvox_cl <- -0.234; label("Fluvoxamine fractional change in CL (fraction)") # CPK Table 2 'Effect of fluvoxamine on CL' = -0.234
    e_dose_cl <- -0.169; label("Power exponent of effective daily dose / 200 mg on CL (unitless)") # CPK Table 2 'Effect of effective daily dose on CL' = -0.169

    # Time-dependent CL: cl * (1 + cl_exp_famp * (1 - exp(-cl_exp_kdes * t)))
    cl_exp_famp <- -0.186; label("Maximum fractional change in CL with time on treatment (fraction)") # CPK Table 2 'Maximum change in CL with respect to time (TAFO)' = -0.186 (fraction, printed with a percent unit)
    lcl_exp_kdes <- log(log(2) / 21.6); label("Rate constant of the time-dependent change in CL (1/h)") # CPK Table 2 'Rate of change in CL with respect to time (half-life), h' = 21.6; rate = ln(2) / 21.6

    # Allometric exponents (reference 70 kg)
    e_wt_cl <- 0.453; label("Power exponent of WT / 70 kg on CL and Q (unitless)") # CPK Table 2 'Effect of weight on CL and Q' = 0.453
    e_wt_vc <- 0.52; label("Power exponent of WT / 70 kg on Vc and Vp (unitless)") # CPK Table 2 'Effect of weight on Vc and Vp' = 0.52

    # IIV: CPK Table 2 CV = sqrt(omega^2) x 100 (footnote a), so omega = CV / 100
    etalcl + etalvc ~ c(0.332929, 0.0778712, 0.171396) # CPK Table 2 omega CL 57.7 and omega Vc 41.4 (CV = sqrt(omega^2)); rho CL-Vc 0.326 so cov = 0.326 * 0.577 * 0.414

    # Residual error: W = sqrt(IPRED^2 * PRO^2 + ADD^2) with EPS variance 1 held constant
    propSd <- 0.437; label("Proportional residual SD, reference studies (fraction)") # CPK Table 2 'RUV PRO, SD' = 0.437
    addSd <- 0.509; label("Additive residual SD (ng/mL)") # CPK Table 2 'RUV ADD, SD' = 0.509
    e_studymod_propSd <- 0.495; label("Moderate-variability study fractional change in proportional SD (fraction)") # CPK Table 2 'Effect of moderate-variability studies on RUV PRO' = 0.495
    e_studyhi_propSd <- 1.16; label("High-variability study fractional change in proportional SD (fraction)") # CPK Table 2 'Effect of high-variability studies on RUV PRO' = 1.16

    # ------------------------------------------------------------------
    # Platelet PD (Br J Clin Pharmacol 2022 Table 2 and Equations 1-7).
    # Reference individual for the typical values: male, 30 years, atopic
    # dermatitis, baseline WBC 6.9 x 1000/uL and haematocrit 43%.
    # ------------------------------------------------------------------
    lrbase <- log(270); label("Baseline platelet count CIRC0 (10^9/L)") # Table 2 'Baseline platelet count (CIRC0; x1000/uL)' = 270
    lgamma <- log(0.232); label("Feedback exponent gamma on (CIRC0 / CIRC) (unitless)") # Table 2 'Feedback exponent (GAM; gamma)' = 0.232
    lmtt <- log(7.21 * 24); label("Mean transit time MTT from the proliferating pool to the circulation (h)") # Table 2 'Mean transit time (MTT; days)' = 7.21 days = 173.04 h
    lemax <- log(0.109); label("Maximum fractional inhibition of platelet proliferation Emax at the typical baseline (unitless)") # Table 2 'Maximum drug effect on platelet proliferation (Emax)' = 0.109
    lec50 <- log(55.3); label("Abrocitinib concentration at half-maximal effect EC50 (ng/mL)") # Table 2 'Concentration at 50% of maximum drug effect (EC50; ng/mL)' = 55.3
    e_rbase_emax <- -1.11; label("Power exponent of (CIRC0 / 270) on Emax (unitless)") # Table 2 'Effect of CIRC0 on Emax' = -1.11 (95% CI -1.21 to -1.01); Results 3.2 power model referenced to the typical CIRC0

    # Tachyphylaxis, Equation 7: TCPX = DELTATCPX * (1 - exp(-ln2 / kTCPX * t))
    # and Emax(t) = Emax * (1 - TCPX). Written in the time-varying-parameter
    # form Emax(t) = Emax * (1 + emax_exp_famp * (1 - exp(-emax_exp_kdes * t))),
    # so emax_exp_famp = -DELTATCPX.
    emax_exp_famp <- fixed(-1); label("Maximum fractional change in Emax with time on treatment, -DELTATCPX (fraction)") # Results 3.2: DELTATCPX estimated 0.991, then fixed to 1 in the final model
    lemax_exp_kdes <- log(log(2) / (151 * 24)); label("Rate constant of the loss of Emax with time on treatment (1/h)") # Table 2 'Half-life of tachyphylaxis (kTCPX; days)' = 151 days = 3624 h; rate = ln(2) / 3624

    # Covariate effects on CIRC0 (power for continuous, Equation 10;
    # fractional for categorical, Equation 9)
    e_wbc_rbase <- 0.199; label("Power exponent of (WBC / 6.9 x 1000/uL) on CIRC0 (unitless)") # Table 2 'Effect of baseline WBC on CIRC0 (referenced to 6.9 x1000/uL)' = 0.199
    e_sexf_rbase <- 0.0584; label("Fractional change in CIRC0 for females (fraction)") # Table 2 'Effect of female sex on CIRC0' = 0.0584
    e_age_rbase <- -0.0695; label("Power exponent of (AGE / 30 years) on CIRC0 (unitless)") # Table 2 'Effect of age on CIRC0 (referenced to 30 years)' = -0.0695 (95% CI -0.0854 to -0.0536)
    e_hct_rbase <- -0.203; label("Power exponent of (HCT / 43%) on CIRC0 (unitless)") # Table 2 'Effect of baseline haematocrit on CIRC0 (referenced to 43%)' = -0.203 (95% CI -0.301 to -0.105)
    e_dis_psoriasis_rbase <- -0.118; label("Fractional change in CIRC0 for psoriasis versus atopic dermatitis patients (fraction)") # Table 2 'Effect of psoriasis patients on CIRC0' = -0.118 (95% CI -0.159 to -0.0768)

    # IIV on CIRC0 only. Table 2 footnote: CV = sqrt(omega^2) x 100, so
    # omega^2 = 0.175^2.
    etalrbase ~ 0.030625 # Table 2 'omega CIRC0 (% CV)' = 17.5 -> 0.175^2

    # Residual error: NONMEM W = sqrt(RUVPRO^2 * IPRE^2 + RUVADD^2) with the
    # EPS variance fixed to 1 (Table 2 footnote).
    propSd_PLT <- 0.120; label("Proportional residual SD on platelet count (fraction)") # Table 2 'Proportional residual error (RUV PRO; SD)' = 0.120
    addSd_PLT <- 16.4; label("Additive residual SD on platelet count (10^9/L)") # Table 2 'Additive residual error (RUV ADD; SD)' = 16.4
  })

  model({
    # ---- Abrocitinib PK (Clin Pharmacokinet 2022 popPK model) ----
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

    # ---- Platelet PD (Br J Clin Pharmacol 2022 Equations 1-7) ----
    # Individual baseline platelet count
    rbase <- exp(lrbase + etalrbase) *
      (WBC / 6.9)^e_wbc_rbase *
      (AGE / 30)^e_age_rbase *
      (HCT / 43)^e_hct_rbase *
      (1 + e_sexf_rbase * SEXF) *
      (1 + e_dis_psoriasis_rbase * DIS_PSORIASIS)
    gamma <- exp(lgamma)
    mtt <- exp(lmtt)
    ec50 <- exp(lec50)

    # Equation 6 with n = 3 transit compartments; kPROL = kTR = kCIRC by the
    # pre-dose steady-state assumption.
    ktr <- 4 / mtt

    # Emax scales with the individual baseline relative to the typical 270,
    # and wanes with time after the first dose (Equation 7). The solver time
    # t is the time after first dose, so the first abrocitinib dose must sit
    # at t = 0.
    emax_exp_kdes <- exp(lemax_exp_kdes)
    emax <- exp(lemax) * (rbase / exp(lrbase))^e_rbase_emax *
      (1 + emax_exp_famp * (1 - exp(-emax_exp_kdes * t)))
    edrug <- emax * Cc / (ec50 + Cc)

    d/dt(prol) <- ktr * prol * (1 - edrug) * (rbase / circ)^gamma - ktr * prol
    d/dt(transit1) <- ktr * prol - ktr * transit1
    d/dt(transit2) <- ktr * transit1 - ktr * transit2
    d/dt(transit3) <- ktr * transit2 - ktr * transit3
    d/dt(circ) <- ktr * transit3 - ktr * circ

    prol(0) <- rbase
    transit1(0) <- rbase
    transit2(0) <- rbase
    transit3(0) <- rbase
    circ(0) <- rbase

    PLT <- circ
    PLT ~ add(addSd_PLT) + prop(propSd_PLT)
  })
}
