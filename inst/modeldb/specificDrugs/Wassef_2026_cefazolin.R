Wassef_2026_cefazolin <- function() {
  description <- "Two-compartment population PK model for intravenous bolus cefazolin surgical prophylaxis (2 g below 120 kg, 3 g at or above 120 kg) in adults undergoing sleeve gastrectomy or laparoscopic abdominal surgery, most with obesity (Wassef 2026, n = 87). Total plasma concentration plus a biophase (effect) compartment whose concentration times a partition coefficient gives the total subcutaneous adipose tissue concentration. Final covariates are bioimpedance (BIA) measured body composition: age and fat-free mass on CL, abdominal-trunk fat-free mass and trunk muscle mass on V1, total body water percentage on V2, body weight on ke0 and body fat percentage on the adipose partition coefficient. Log-normal IIV on all six parameters with a CL-V1 correlation; combined additive + proportional plasma error and additive tissue error. Estimated in Monolix 2024R1."
  reference <- "Wassef A, Kosicka-Noworzyn K, Back H, Siemiatkowska A, Yohn C, Sadek R, Kagan L, Brunetti L. Cefazolin surgical prophylaxis in obesity: a body composition-driven population pharmacokinetic approach. Antimicrob Agents Chemother. 2026. doi:10.1128/aac.01677-25. PMCID PMC13321816. Structural ODEs from Fig. 4b; all parameter estimates from Supplementary Table S2; covariates enter uncentred (see model())."
  vignette <- "Wassef_2026_cefazolin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Wassef 2026 Methods: total plasma cefazolin by HPLC-UV (0.5-250 ug/mL) and
  # total subcutaneous adipose tissue cefazolin by LC-MS/MS (0.15-150 ug/g).
  # The biophase compartment carries a concentration (Fig. 4b eq. 3); it holds
  # no drug mass, so it does not deplete `central`.
  compartmentData <- list(
    central = list(analyte = "cefazolin", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "cefazolin", units = "mg", specimen = "plasma", verified = TRUE),
    effect = list(
      analyte = "cefazolin",
      units = "mg/L",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Uncentred. Effect on CL: exp(-0.011 * AGE). Table 1 dose-group means 49 years (n = 57) and 41 years (n = 30).",
      source_name = "AGE"
    ),
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Uncentred. Effect on ke0: exp(0.025 * WT), so ke0 is about 32/h at the Table S1 median of 112.9 kg. Measured on the Tanita MC-780 analyser at the BIA visit.",
      source_name = "WTKG"
    ),
    FFM = list(
      description = "Whole-body fat-free mass measured by segmental bioimpedance analysis",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "MEASURED (Tanita MC-780 segmental BIA, mean of repeated readings), not the Janmahasatian formula estimate the FFM register entry describes; the same biological quantity, so it uses the same canonical. Uncentred. Effect on CL: exp(0.0052 * FFM). Table S1 cohort median 65 kg.",
      source_name = "FFMAMT"
    ),
    FFM_TRUNK = list(
      description = "Abdominal-trunk segmental fat-free mass measured by bioimpedance analysis",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Tanita MC-780 trunk segment. Uncentred. Effect on V1: exp(-0.64 * FFM_TRUNK). Table S1 cohort median 37.1 kg. Near-collinear with MUSCLE_MASS_TRUNK (the two differ by the trunk bone-mineral share, cohort means 37.33 vs 35.67 kg) and the two V1 coefficients are large and opposite, so they must be supplied as a measured pair, never drawn independently.",
      source_name = "MBDFFMAMT"
    ),
    MUSCLE_MASS_TRUNK = list(
      description = "Abdominal-trunk segmental muscle mass measured by bioimpedance analysis",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Tanita MC-780 trunk segment. Uncentred. Effect on V1: exp(0.69 * MUSCLE_MASS_TRUNK). Table S1 cohort mean 35.67 kg (its printed median of 25.4 kg lies outside the 32.7 and 41.2 kg dose-group medians and is a misprint). Together with FFM_TRUNK, V1 scales with exp(0.69 * (MUSCLE_MASS_TRUNK - FFM_TRUNK) + 0.05 * FFM_TRUNK): a 1 kg change in the trunk fat-free-mass minus muscle gap halves or doubles V1.",
      source_name = "MBDMSLAMT"
    ),
    BODYWATER_PCT = list(
      description = "Total body water as a percentage of total body weight, measured by bioimpedance analysis",
      units = "%",
      type = "continuous",
      reference_category = NULL,
      notes = "Tanita MC-780. Uncentred. Effect on V2: exp(0.033 * BODYWATER_PCT). Table S1 cohort median 42.7 %.",
      source_name = "WTRPER"
    ),
    BODYFAT_PCT = list(
      description = "Whole-body fat percentage measured by bioimpedance analysis",
      units = "%",
      type = "continuous",
      reference_category = NULL,
      notes = "Tanita MC-780 segmental BIA. Uncentred. Effect on the adipose partition coefficient: exp(-0.063 * BODYFAT_PCT), giving Kp of about 0.029 at the Table S1 median of 40.2 % and spanning about 12-fold across the observed 18-58 % range and drives the paper's finding that adipose target attainment falls with adiposity.",
      source_name = "FATPER"
    )
  )

  covariatesDataExcluded <- list(
    HT = list(
      description = "Height",
      units = "cm",
      type = "continuous",
      reference_category = NULL,
      notes = "Selected on V2 in the Group 1 (standard anthropometrics) and Group 2 (calculated body composition) covariate searches (Fig. 5) but replaced by BODYWATER_PCT in the final Group 3 (BIA-measured) model."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Selected on Kp in the Group 1 search (Fig. 5) but replaced by BODYFAT_PCT in the final Group 3 model."
    ),
    SEXF = list(
      description = "Sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Selected on Kp in the Group 1 search (Fig. 5) but not retained in the final Group 3 model. 46 of 87 subjects female."
    ),
    DIS_DIAB = list(
      description = "Type 2 diabetes mellitus indicator (1 = present)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (absent)",
      notes = "Selected on V1 in the Group 1 and Group 2 searches (and on Kp in Group 2 search 2, Boer formula) per Fig. 5, but not retained in the final Group 3 model. 17 of 87 subjects (Table 1)."
    ),
    LBM = list(
      description = "Lean body weight calculated by the Boer formula",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Selected on CL in the best calculated-metric search (Group 2 search 2, Boer formula; OFV 4643.3, Fig. 5); the BIA-measured Group 3 model (OFV 4629.67) was retained as final."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 87L,
    n_studies = 1L,
    age_range = "18-65 years (inclusion); dose-group means 49 (SD 13.4) and 41 (SD 14.6) years (Table 1)",
    weight_range = "median 112.9 kg, range spanning 162.5 kg (Table S1); 2 g group mean 95.4 kg, 3 g group mean 146.9 kg (Table 1)",
    sex_female_pct = 53,
    race_ethnicity = "White 76%, African American 13%, Asian 3%, unknown 6%, more than one race 2% (Table 1)",
    disease_state = "Adults undergoing sleeve gastrectomy (n = 57 analysed) or laparoscopic abdominal hernia surgery, receiving cefazolin surgical prophylaxis. BMI median 38 kg/m^2; 36 of 87 had BMI >= 40 and 70 of 87 had an unhealthy body-fat percentage (> 30% male, > 35% female). Excluded: CrCl < 50 mL/min, transaminases > 3 x ULN, beta-lactam allergy.",
    dose_range = "2 g intravenous bolus (total body weight < 120 kg, n = 57) or 3 g (>= 120 kg, n = 30), repeated every 4 h for the duration of surgery; 85 subjects received one dose and 2 received two.",
    regions = "USA (Robert Wood Johnson Barnabas University Hospital, New Jersey; single centre)",
    renal_function = "Cockcroft-Gault CrCl mean 153 mL/min (2 g) and 271 mL/min (3 g); CKD-EPI eGFR mean 104-108 mL/min/1.73 m^2 (Table 1).",
    notes = "Plasma sampled at baseline, 5, 30, 60, 120 and 240 min after the bolus; two subcutaneous abdominal adipose biopsies (cold scalpel) at incision and at closure. 101 enrolled, 87 analysed for the popPK model. Body composition by Tanita MC-780 segmental bioimpedance."
  )

  ini({
    # Structural parameters -- Supplementary Table S2 (final model).
    lcl <- log(3.9); label("Clearance at covariate value 0, no clinical meaning (L/h)") # Table S2 'CLpop' = 3.9 L/h (RSE 20.25%)
    lvc <- log(3.6); label("Central volume of distribution V1 at covariate value 0, no clinical meaning (L)") # Table S2 'V1pop' = 3.6 L (RSE 21.48%)
    lq <- log(11); label("Intercompartmental clearance (L/h)") # Table S2 'Qpop' = 11 L/h (RSE 10.33%)
    lvp <- log(1.1); label("Peripheral volume of distribution V2 at covariate value 0, no clinical meaning (L)") # Table S2 'V2pop' = 1.1 L (RSE 40.53%)
    lke0 <- log(1.91); label("Plasma-to-biophase equilibration rate constant at covariate value 0, no clinical meaning (1/h)") # Table S2 'ke0pop' = 1.91 1/h (RSE 91.49%)
    lkp_adipose <- log(0.37); label("Subcutaneous adipose to plasma partition coefficient at covariate value 0, no clinical meaning (unitless)") # Table S2 'Kppop' = 0.37 (RSE 32.67%)

    # Covariate coefficients -- Table S2. Encoded in the uncentred exponential
    # form P = Ppop * exp(beta * COV); see model() for why.
    e_age_cl <- -0.011; label("Exponential age coefficient on CL (1/year)") # Table S2 'beta CL-AGE' = -0.011 (RSE 20.49%)
    e_ffm_cl <- 0.0052; label("Exponential fat-free-mass coefficient on CL (1/kg)") # Table S2 'beta CL-FFMAMT' = 0.0052 (RSE 39.8%)
    e_ffm_trunk_vc <- -0.64; label("Exponential trunk fat-free-mass coefficient on V1 (1/kg)") # Table S2 'beta V1-MBDFFMAMT' = -0.64 (RSE 41.01%)
    e_muscle_mass_trunk_vc <- 0.69; label("Exponential trunk muscle-mass coefficient on V1 (1/kg)") # Table S2 'beta V1-MBDMSLAMT' = 0.69 (RSE 39.46%)
    e_bodywater_pct_vp <- 0.033; label("Exponential body-water-percentage coefficient on V2 (1/%)") # Table S2 'beta V2-WTRPER' = 0.033 (RSE 24.3%)
    e_wt_ke0 <- 0.025; label("Exponential body-weight coefficient on ke0 (1/kg)") # Table S2 'beta ke0-WTKG' = 0.025 (RSE 22.59%)
    e_bodyfat_pct_kp_adipose <- -0.063; label("Exponential body-fat-percentage coefficient on Kp (1/%)") # Table S2 'beta Kp-FATPER' = -0.063 (RSE 12.31%)

    # IIV -- Table S2 'Standard Deviation of the Random Effects' (log-scale
    # SDs; each row's C.V.(%) column equals sqrt(exp(omega^2) - 1)). Variances
    # are the squared SDs; cov(CL, V1) = 0.59 * 0.27 * 0.25 = 0.039825.
    etalcl + etalvc ~ c(0.0729, 0.039825, 0.0625) # Table S2 Omega CL 0.27 (CV 26.98%), Omega V1 0.25 (CV 25.47%), correlation V1 and CL 0.59
    etalq ~ 0.1681 # Table S2 Omega Q 0.41 (CV 43.21%)
    etalvp ~ 0.1089 # Table S2 Omega V2 0.33 (CV 33.65%)
    etalke0 ~ 0.4624 # Table S2 Omega ke0 0.68 (CV 77.37%)
    etalkp_adipose ~ 0.2601 # Table S2 Omega Kp 0.51 (CV 54.55%)

    # Residual error -- Table S2. Plasma: Monolix 'combined 2'
    # (Results: 'combine 2'), sqrt(a^2 + (b * f)^2). Tissue: additive.
    addSd <- 51.4; label("Additive residual error, plasma (mg/L)") # Table S2 'sigma additive' = 51.4 ug/mL (RSE 21.38%)
    propSd <- 0.095; label("Proportional residual error, plasma (fraction)") # Table S2 'sigma proportional' = 9.5% (RSE 8.38%)
    addSd_Cadipose <- 6.1; label("Additive residual error, subcutaneous adipose tissue (mg/L)") # Table S2 'sigma constant' = 6.1 ug/mL (RSE 17.37%)
  })

  model({
    # Covariates enter uncentred, P = Ppop * exp(beta * COV), so each Ppop is
    # the value at covariate 0 (Monolix's output for an untransformed
    # covariate). The paper states median centring and prints the power form
    # (COV/median)^beta, but neither reproduces its data: the power form moves
    # Kp only 0.39 to 0.36 over 18-58% body fat, and median centring puts
    # adipose concentrations ~15-fold above Fig. 3 with a plasma half-life of
    # ~1 h. Uncentred, the same estimates reproduce Fig. 2, the Table 2 NCA
    # and the Table 3 adipose target attainment (see the vignette).
    cl <- exp(lcl + etalcl + e_age_cl * AGE + e_ffm_cl * FFM)
    vc <- exp(lvc + etalvc + e_ffm_trunk_vc * FFM_TRUNK + e_muscle_mass_trunk_vc * MUSCLE_MASS_TRUNK)
    q <- exp(lq + etalq)
    vp <- exp(lvp + etalvp + e_bodywater_pct_vp * BODYWATER_PCT)
    ke0 <- exp(lke0 + etalke0 + e_wt_ke0 * WT)
    kp_adipose <- exp(lkp_adipose + etalkp_adipose + e_bodyfat_pct_kp_adipose * BODYFAT_PCT)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Fig. 4b equations (1)-(3).
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(effect) <- ke0 * (central / vc - effect)

    # Fig. 4b equations (4)-(5). Tissue concentrations were assayed in ug/g and
    # are reported here in mg/L assuming a tissue density of 1 g/mL.
    Cc <- central / vc
    Cadipose <- kp_adipose * effect

    Cc ~ add(addSd) + prop(propSd) + combined2()
    Cadipose ~ add(addSd_Cadipose)
  })
}
