Usansky_2026_brensocatib <- function() {
  description <- "Two-compartment population PK model for oral brensocatib (DPP1 inhibitor) in healthy adults and in adult and adolescent patients with non-cystic fibrosis bronchiectasis or cystic fibrosis (Usansky 2026), pooled from 11 phase 1-3 trials. Absorption runs through a depot and a chain of three transit compartments that share one rate constant; relative bioavailability rises with dose by a sigmoid Emax function. Covariates: body weight (fixed allometry on all clearances and volumes), sex, Black and Asian race, hepatic impairment and rifampicin, clarithromycin and verapamil coadministration on CL/F; age on Vc/F; fed state and oral-solution formulation on ka; clarithromycin and verapamil on relative bioavailability. Between-occasion variability on ka and bioavailability, and a combined additive and proportional residual error per study phase."
  reference <- "Usansky H, Cammarata A, Li C, Au Yeung S, Rubino CM, Wang L, Teper A, Fernandez C, Mange K. Population Pharmacokinetic Analysis of Brensocatib in Adult and Adolescent Patients With Non-Cystic Fibrosis Bronchiectasis. CPT Pharmacometrics Syst Pharmacol. 2026;15(9):e70340. doi:10.1002/psp4.70340"
  vignette <- "Usansky_2026_brensocatib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric scaling relative to 70 kg with exponents fixed at 0.75 on CL/F and Q/F and 1 on Vc/F and Vp/F (Usansky 2026 Methods and Table 2 footnote a; control stream WTKG). Estimating the exponents worsened the fit (dOFV = +2348).",
      source_name = "WTKG"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on Vc/F centred at 49 years: (AGE/49)^0.314 (Usansky 2026 Vc/F equation and Table 2; control stream VCAGE). Adolescents aged 12 to < 18 years were included.",
      source_name = "AGE"
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "CL/F multiplied by (1 - 0.0957) in females (Usansky 2026 Table 2 'CL/F-female'; control stream CLSEXF).",
      source_name = "SEXF"
    ),
    RACE_BLACK = list(
      description = "Black or African American race indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (any race other than Black or Asian; White is the forest-plot reference)",
      notes = "CL/F multiplied by (1 + 0.492) (Usansky 2026 Table 2 'CL/F-Black'). Source column RACEN = 2; derive RACE_BLACK = as.integer(RACEN == 2).",
      source_name = "RACEN"
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (any race other than Black or Asian)",
      notes = "CL/F multiplied by (1 + 0.229) (Usansky 2026 Table 2 'CL/F-Asian'). Source column RACEN = 3; derive RACE_ASIAN = as.integer(RACEN == 3).",
      source_name = "RACEN"
    ),
    HEPIMP = list(
      description = "Any hepatic impairment (Child-Pugh mild, moderate or severe) indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (normal hepatic function)",
      notes = "CL/F multiplied by (1 + 0.187), fixed from the dedicated hepatic-impairment study NCT05517525 (Usansky 2026 Table 2 footnote b). The printed CL/F equation writes the factor as 1.187^HI with HI the Child-Pugh classification number, but the published control stream applies it as a single indicator, IF(CPCLASSN.NE.0) CLCPCLASSN = (1 + THETA(17)), so mild, moderate and severe share one effect; the control stream form is used.",
      source_name = "CPCLASSN"
    ),
    FED = list(
      description = "Fed state at dosing",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (fasted)",
      notes = "Selects the absorption rate constant: ka = 7.61 1/h fasted and 5.80 1/h fed (Usansky 2026 Table 2; control stream KA = (THETA(5)*FAST + THETA(6)*FED)*KAFORMN with FED = 1 - FAST). Keep the value constant from a dose record through the observations that follow it.",
      source_name = "FAST (FED = 1 - FAST)"
    ),
    FORM_SOLUTION = list(
      description = "Oral-solution formulation indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (phase 2 or phase 3 tablet)",
      notes = "ka multiplied by (1 + 0.809) for the oral solution (Usansky 2026 Table 2 'Ka-oral solution formulation'; control stream KAFORMN, FORMN = 1). The phase 2 and phase 3 tablets share the reference ka.",
      source_name = "FORMN (FORM_SOLUTION = as.integer(FORMN == 1))"
    ),
    CONMED_RIFAMPICIN = list(
      description = "Coadministered rifampicin (rifampin), a strong CYP3A inducer",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no rifampicin)",
      notes = "CL/F multiplied by (1 + 0.383), fixed from the drug-interaction study NCT05826574 (Usansky 2026 Table 2 footnote c; control stream CONMEDN = 1). The source codes coadministered rifampicin, clarithromycin and verapamil as one categorical column, so at most one of the three CONMED_* indicators is 1 on any record.",
      source_name = "CONMEDN (CONMED_RIFAMPICIN = as.integer(CONMEDN == 1))"
    ),
    CONMED_CLARITHROMYCIN = list(
      description = "Coadministered clarithromycin, a strong CYP3A and P-gp inhibitor",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no clarithromycin)",
      notes = "CL/F multiplied by (1 - 0.0779) and relative bioavailability by (1 + 0.388), both fixed from the drug-interaction study NCT05965570 (Usansky 2026 Table 2 footnote d; control stream CONMEDN = 3).",
      source_name = "CONMEDN (CONMED_CLARITHROMYCIN = as.integer(CONMEDN == 3))"
    ),
    CONMED_VERAPAMIL = list(
      description = "Coadministered verapamil, a P-gp inhibitor",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no verapamil)",
      notes = "CL/F multiplied by (1 + 0.212) and relative bioavailability by (1 + 0.657), both fixed from the drug-interaction study NCT02653872 (Usansky 2026 Table 2 footnote e; control stream CONMEDN = 37).",
      source_name = "CONMEDN (CONMED_VERAPAMIL = as.integer(CONMEDN == 37))"
    ),
    DOSE = list(
      description = "Administered brensocatib dose of the current dosing record",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "Use case (a) of the register entry: drives the dose-dependent relative bioavailability F = 0.29 + 0.346 * DOSE^1.54 / (21.6^1.54 + DOSE^1.54) (Usansky 2026 relF equation and Table 2; control stream DOSEMG). Studied range 5-120 mg single dose and 10-40 mg once daily. Put the event columns (id, time, amt, evid, cmt) before DOSE in the event table.",
      source_name = "DOSEMG"
    ),
    STUDY_PHASE2 = list(
      description = "Phase 2 study record indicator (residual-error stratum)",
      units = "(binary)",
      type = "binary",
      reference_category = "0; with STUDY_PHASE3 = 0 selects the phase 1 stratum",
      notes = "Selects the phase 2 combined additive and proportional residual error (control stream $ERROR P2_flag, PHASE = 2). Phase 2 studies: NCT05090904 and WILLOW (NCT03218917). Touches no structural parameter.",
      source_name = "PHASE (STUDY_PHASE2 = as.integer(PHASE == 2))"
    ),
    STUDY_PHASE3 = list(
      description = "Phase 3 study record indicator (residual-error stratum)",
      units = "(binary)",
      type = "binary",
      reference_category = "0; with STUDY_PHASE2 = 0 selects the phase 1 stratum",
      notes = "Selects the phase 3 combined additive and proportional residual error (control stream $ERROR P3_flag, PHASE = 3). The only phase 3 study is ASPEN (NCT04594369). Touches no structural parameter.",
      source_name = "PHASE (STUDY_PHASE3 = as.integer(PHASE == 3))"
    ),
    OCC = list(
      description = "Occasion index for between-occasion variability on ka and relative bioavailability",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = "The control stream defines two occasions, IF(OCC.EQ.1) OCC1 = 1 and IF(OCC.EQ.2) OCC2 = 1, with one shared variance per parameter ($OMEGA BLOCK(1) SAME). Records with any other OCC value (for example 0) carry no between-occasion deviation. The paper does not say how occasions were defined.",
      source_name = "OCC"
    )
  )

  covariatesDataExcluded <- list(
    CRCL = list(
      description = "Creatinine clearance by Cockcroft-Gault, normalized to body surface area",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      notes = "Part of the base model on CL/F but not statistically significant in the final pooled model (Usansky 2026 Results: 'renal functions (e.g., CrCl) ... were not statistically significant PK covariates'). No coefficient is reported."
    ),
    RENALIMP = list(
      description = "Renal impairment indicator (mild, moderate or severe)",
      units = "(binary)",
      type = "binary",
      notes = "Tested and not retained; no coefficient is reported (Usansky 2026 Abstract and Discussion)."
    ),
    DIS_HEALTHY = list(
      description = "Healthy-volunteer indicator (versus bronchiectasis or cystic fibrosis patients)",
      units = "(binary)",
      type = "binary",
      notes = "Disease population was not a statistically significant covariate (Usansky 2026 Results). No coefficient is reported."
    ),
    DIS_CF = list(
      description = "Cystic fibrosis patient indicator",
      units = "(binary)",
      type = "binary",
      notes = "Disease population was not a statistically significant covariate (Usansky 2026 Results). No coefficient is reported."
    )
  )

  compartmentData <- list(
    depot = list(analyte = "brensocatib", units = "mg", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "brensocatib", units = "mg", specimen = "administration site", verified = TRUE),
    transit2 = list(analyte = "brensocatib", units = "mg", specimen = "administration site", verified = TRUE),
    transit3 = list(analyte = "brensocatib", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "brensocatib", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "brensocatib", units = "mg", specimen = "plasma", verified = TRUE)
  )

  # Bioavailability has no estimated log-scale typical value to pair with:
  # its typical value is the dose-dependent Emax function in model(), and the
  # log-normal IIV multiplies it (control stream BIO = TVBIO*EXP(ETA(6)+IOV)).
  paper_specific_etas <- c("etalfdepot")

  population <- list(
    species = "human",
    n_subjects = 1098,
    n_studies = 11,
    n_observations = 11643,
    age_mean = "51 years (CV 43.1%); adolescents 14.4 years (CV 11.4%), adults 52.9 years (CV 36.8%)",
    age_range = "12 years and older (31 adolescents aged 12 to < 18 years)",
    weight_mean = "69.2 kg (CV 23%); adolescents 50.4 kg, adults 69.8 kg",
    sex_female_pct = 55.4,
    race_ethnicity = c(
      White = 71.8,
      Asian = 11.6,
      Black = 6.0,
      AmericanIndianAlaskaNative = 1.2,
      NativeHawaiianPacificIslander = 0.3,
      OtherUnknown = 9.2
    ),
    renal_function = "Mean CrCl 82.3 mL/min/1.73 m^2 (CV 34.3%); 18 participants with mild, moderate or severe renal impairment from a dedicated study",
    hepatic_function = "18 participants with mild, moderate or severe (Child-Pugh) hepatic impairment from a dedicated study",
    disease_state = "Healthy volunteers (26.5%), non-cystic fibrosis bronchiectasis (71.3%, including 31 adolescents) and cystic fibrosis (2.2%)",
    dose_range = "5-120 mg single oral doses; 10-40 mg once daily (ASPEN: 10 or 25 mg once daily)",
    co_medication = "Rifampicin (14), clarithromycin (22), verapamil (15) and esomeprazole (16) participants from drug-interaction studies",
    notes = "Pooled 8 phase 1, 2 phase 2 (NCT05090904, WILLOW) and 1 phase 3 (ASPEN, NCT04594369) trials; intensive sampling in phase 1, sparse in phase 2/3. Baseline demographics in Usansky 2026 Table 1."
  )

  ini({
    # ---- Disposition, typical values for a 70 kg, 49-year-old reference ----
    lcl <- log(3.25)
    label("Apparent clearance CL/F for 70 kg (L/h)") # Table 2 'CL/F' 3.25
    lvc <- log(78.8)
    label("Apparent central volume Vc/F for 70 kg, 49 years (L)") # Table 2 'Vc/F' 78.8
    lq <- log(5.97)
    label("Apparent intercompartmental clearance Q/F for 70 kg (L/h)") # Table 2 'Q/F' 5.97
    lvp <- log(50.4)
    label("Apparent peripheral volume Vp/F for 70 kg (L)") # Table 2 'Vp/F' 50.4

    # ---- Absorption: one rate constant for depot, the 3 transits and entry to central ----
    lka_fasted <- log(7.61)
    label("Absorption and transit rate constant, fasted (1/h)") # Table 2 'Ka-fasted' 7.61
    lka_fed <- log(5.8)
    label("Absorption and transit rate constant, fed (1/h)") # Table 2 'Ka-fed' 5.8
    e_form_solution_ka <- 0.809
    label("Proportional change in ka for the oral solution (fraction)") # Table 2 'Ka-oral solution formulation' 0.809

    # ---- Dose-dependent relative bioavailability (sigmoid Emax in dose) ----
    e0_dose_fdepot <- 0.29
    label("Relative bioavailability extrapolated to zero dose (fraction)") # Table 2 'E0' 0.29
    emax_dose_fdepot <- 0.346
    label("Maximum increase in relative bioavailability with dose (fraction)") # Table 2 'E max' 0.346
    ed50_dose_fdepot <- 21.6
    label("Dose giving half the maximum bioavailability increase (mg)") # Table 2 'EC50, mg' 21.6
    hill_dose_fdepot <- 1.54
    label("Hill coefficient of the dose effect on relative bioavailability (unitless)") # Table 2 'H' 1.54

    # ---- Allometric exponents, fixed ----
    e_wt_cl_q <- fixed(0.75)
    label("Allometric exponent on WT/70 for CL/F and Q/F (unitless)") # Methods; Table 2 footnote a
    e_wt_vc_vp <- fixed(1)
    label("Allometric exponent on WT/70 for Vc/F and Vp/F (unitless)") # Methods; Table 2 footnote a

    # ---- Covariate effects ----
    e_age_vc <- 0.314
    label("Power exponent on AGE/49 for Vc/F (unitless)") # Table 2 'Vc/F-age' 0.314
    e_sexf_cl <- -0.0957
    label("Proportional change in CL/F for females (fraction)") # Table 2 'CL/F-female' -0.0957
    e_race_black_cl <- 0.492
    label("Proportional change in CL/F for Black participants (fraction)") # Table 2 'CL/F-Black' 0.492
    e_race_asian_cl <- 0.229
    label("Proportional change in CL/F for Asian participants (fraction)") # Table 2 'CL/F-Asian' 0.229
    e_hepimp_cl <- fixed(0.187)
    label("Proportional change in CL/F with hepatic impairment (fraction)") # Table 2 'CL/F-hepatic impairment' 0.187; control stream THETA(17) 0.187 FIX
    e_conmed_rifampicin_cl <- fixed(0.383)
    label("Proportional change in CL/F with rifampicin (fraction)") # Table 2 'CL/F-rifampin' 0.383; control stream THETA(12) 0.383 FIX
    e_conmed_clarithromycin_cl <- fixed(-0.0779)
    label("Proportional change in CL/F with clarithromycin (fraction)") # Table 2 'CL/F-clarithromycin' -0.0779; control stream THETA(13) -0.0779 FIX
    e_conmed_verapamil_cl <- fixed(0.212)
    label("Proportional change in CL/F with verapamil (fraction)") # Table 2 'CL/F-verapamil' 0.212; control stream THETA(15) 0.212 FIX
    e_conmed_clarithromycin_fdepot <- fixed(0.388)
    label("Proportional change in relative bioavailability with clarithromycin (fraction)") # Table 2 'relF-clarithromycin' 0.388; control stream THETA(14) 0.388 FIX
    e_conmed_verapamil_fdepot <- fixed(0.657)
    label("Proportional change in relative bioavailability with verapamil (fraction)") # Table 2 'relF-verapamil' 0.657; control stream THETA(16) 0.657 FIX

    # ---- Between-subject variability (log-normal). The paper does not report
    # the final variances; these are the $OMEGA values of the published final
    # control stream (Data S1), which are the starting values of the final run.
    # IIV on Q/F is fixed to zero there ($OMEGA 0 FIX) and is omitted. ----
    etalcl ~ 0.121 # Data S1 $OMEGA IIV_CL 0.121
    etalvc ~ 0.0492 # Data S1 $OMEGA IIV_VC 0.0492
    etalvp ~ 0.0233 # Data S1 $OMEGA IIV_VP 0.0233
    etalka ~ 0.164 # Data S1 $OMEGA IIV_KA 0.164
    etalfdepot ~ 0.0132 # Data S1 $OMEGA IIV_BIO 0.0132

    # ---- Between-occasion variability, two occasions sharing one variance ----
    etaiov_ka_1 ~ 0.184 # Data S1 $OMEGA IOV_KA_OCC1 0.184
    etaiov_ka_2 ~ fixed(0.184) # Data S1 $OMEGA IOV_KA_OCC2 SAME
    etaiov_fdepot_1 ~ 0.0354 # Data S1 $OMEGA IOV_BIO_OCC1 0.0354
    etaiov_fdepot_2 ~ fixed(0.0354) # Data S1 $OMEGA IOV_BIO_OCC2 SAME

    # ---- Residual error per study phase. $SIGMA entries are variances with
    # Y = IPRED*(1 + EPS_prop) + EPS_add, so the SDs are their square roots. ----
    propSdPhase1 <- 0.1562
    label("Proportional residual SD, phase 1 studies (fraction)") # Data S1 $SIGMA RV_prop_P1 0.0244; sqrt = 0.1562
    addSdPhase1 <- 0.3032
    label("Additive residual SD, phase 1 studies (ng/mL)") # Data S1 $SIGMA RV_add_P1 0.0919; sqrt = 0.3032
    propSdPhase2 <- 0.2530
    label("Proportional residual SD, phase 2 studies (fraction)") # Data S1 $SIGMA RV_prop_P2 0.064; sqrt = 0.2530
    addSdPhase2 <- 0.6380
    label("Additive residual SD, phase 2 studies (ng/mL)") # Data S1 $SIGMA RV_add_P2 0.407; sqrt = 0.6380
    propSdPhase3 <- 0.3450
    label("Proportional residual SD, phase 3 study (fraction)") # Data S1 $SIGMA RV_prop_P3 0.119; sqrt = 0.3450
    addSdPhase3 <- 0.6565
    label("Additive residual SD, phase 3 study (ng/mL)") # Data S1 $SIGMA RV_add_P3 0.431; sqrt = 0.6565
  })

  model({
    # 1. Occasion indicators; any OCC other than 1 or 2 carries no IOV
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    iov_ka <- oc1 * etaiov_ka_1 + oc2 * etaiov_ka_2
    iov_fdepot <- oc1 * etaiov_fdepot_1 + oc2 * etaiov_fdepot_2

    # 2. Covariate factors (control stream $PK)
    cov_cl <- (WT / 70)^e_wt_cl_q *
      (1 + e_conmed_rifampicin_cl * CONMED_RIFAMPICIN) *
      (1 + e_conmed_clarithromycin_cl * CONMED_CLARITHROMYCIN) *
      (1 + e_conmed_verapamil_cl * CONMED_VERAPAMIL) *
      (1 + e_hepimp_cl * HEPIMP) *
      (1 + e_race_black_cl * RACE_BLACK) *
      (1 + e_race_asian_cl * RACE_ASIAN) *
      (1 + e_sexf_cl * SEXF)
    cov_vc <- (AGE / 49)^e_age_vc * (WT / 70)^e_wt_vc_vp
    cov_q <- (WT / 70)^e_wt_cl_q
    cov_vp <- (WT / 70)^e_wt_vc_vp

    # 3. Individual parameters
    cl <- exp(lcl + etalcl) * cov_cl
    vc <- exp(lvc + etalvc) * cov_vc
    q <- exp(lq) * cov_q
    vp <- exp(lvp + etalvp) * cov_vp
    ka <- (exp(lka_fasted) * (1 - FED) + exp(lka_fed) * FED) *
      (1 + e_form_solution_ka * FORM_SOLUTION) * exp(etalka + iov_ka)

    frel_dose <- e0_dose_fdepot +
      emax_dose_fdepot * DOSE^hill_dose_fdepot /
        (ed50_dose_fdepot^hill_dose_fdepot + DOSE^hill_dose_fdepot)
    frel_conmed <- (1 + e_conmed_clarithromycin_fdepot * CONMED_CLARITHROMYCIN) *
      (1 + e_conmed_verapamil_fdepot * CONMED_VERAPAMIL)
    fdepot <- frel_dose * frel_conmed * exp(etalfdepot + iov_fdepot)

    # 4. Micro-constants
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 5. ODEs: depot -> transit1 -> transit2 -> transit3 -> central, all at ka
    d/dt(depot) <- -ka * depot
    d/dt(transit1) <- ka * depot - ka * transit1
    d/dt(transit2) <- ka * transit1 - ka * transit2
    d/dt(transit3) <- ka * transit2 - ka * transit3
    d/dt(central) <- ka * transit3 - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # 6. Bioavailability
    f(depot) <- fdepot

    # 7. Observation: dose in mg, volume in L, concentration in ng/mL
    Cc <- 1000 * central / vc
    phase1 <- 1 - STUDY_PHASE2 - STUDY_PHASE3
    propSd <- propSdPhase1 * phase1 + propSdPhase2 * STUDY_PHASE2 + propSdPhase3 * STUDY_PHASE3
    addSd <- addSdPhase1 * phase1 + addSdPhase2 * STUDY_PHASE2 + addSdPhase3 * STUDY_PHASE3
    Cc ~ add(addSd) + prop(propSd)
  })
}
