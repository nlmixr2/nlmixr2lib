Koele_2025_btz043 <- function() {
  description <- "Joint parent-plus-two-metabolite population PK model for the first-in-class benzothiazinone antituberculosis drug BTZ-043 and its metabolites M1 and M2 in adults with drug-susceptible pulmonary tuberculosis: two-compartment BTZ-043 and M2 disposition with one-compartment M1, dual parallel oral absorption (a two-transit-compartment chain feeding an absorption compartment, plus a lag-time secondary route that carries 36% of the bioavailable dose only when the dose is taken with food), allometric scaling on 70 kg with fixed exponents, prandial-state and high-dose effects on relative bioavailability, a prandial-state effect on mean transit time and on the relative fraction metabolised to M1, slower BTZ-043 clearance in Cape-coloured participants, and a step decrease in M2 clearance after 10 days on treatment. Amounts are carried in nmol (the mg dose is converted with the 431.39 g/mol BTZ-043 molecular weight inside the bioavailability term) so that parent and metabolite mass balance is preserved; fractions metabolised and absolute bioavailability were fixed to 1, so every clearance and volume is an apparent value."
  reference <- paste(
    "Koele S. E., Heinrich N., De Jager V. R., Dreisbach J., Phillips P. P. J.,",
    "Gross-Demel P., Dawson R., Narunsky K., Wildner L. M., Mchugh T. D.,",
    "Te Brake L. H. M., Diacon A. H., Aarnoutse R. E., Hoelscher M.,",
    "Svensson E. M. (2025).",
    "Population pharmacokinetics and exposure-response relationship of the",
    "antituberculosis drug BTZ-043.",
    "Journal of Antimicrobial Chemotherapy 80(5):1319-1327.",
    "doi:10.1093/jac/dkaf076.",
    "Structural equations and random-effect variances transcribed from the",
    "final NONMEM control stream in the Supplementary data",
    "('Pharmacokinetic model code'); typical values and residual errors from",
    "Table 2.",
    sep = " "
  )
  vignette <- "Koele_2025_btz043"
  units <- list(time = "h", dosing = "mg", concentration = "nmol/L")

  compartmentData <- list(
    depot1 = list(
      analyte = "BTZ-043", units = "nmol",
      specimen = "administration site", verified = TRUE
    ),
    transit1 = list(
      analyte = "BTZ-043", units = "nmol",
      specimen = "administration site", verified = TRUE
    ),
    transit2 = list(
      analyte = "BTZ-043", units = "nmol",
      specimen = "administration site", verified = TRUE
    ),
    depot2 = list(
      analyte = "BTZ-043", units = "nmol",
      specimen = "administration site", verified = TRUE
    ),
    central = list(
      analyte = "BTZ-043", units = "nmol",
      specimen = "plasma", verified = TRUE
    ),
    peripheral1 = list(
      analyte = "BTZ-043", units = "nmol",
      specimen = "plasma", verified = TRUE
    ),
    central_m1 = list(
      analyte = "BTZ-043 metabolite M1", units = "nmol",
      specimen = "plasma", verified = TRUE
    ),
    central_m2 = list(
      analyte = "BTZ-043 metabolite M2", units = "nmol",
      specimen = "plasma", verified = TRUE
    ),
    peripheral1_m2 = list(
      analyte = "BTZ-043 metabolite M2", units = "nmol",
      specimen = "plasma", verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight, used for allometric scaling of every apparent clearance and volume around a 70 kg reference adult.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Koele 2025 Methods 'PK model development': 'All volumes and clearances were normalized to 70 kg and allometrically scaled, using fixed exponents of 1 and 0.75, respectively.' The supplementary control stream implements this as AlloCL = (WT/70)**0.75 and AlloV = (WT/70)**1, applied to CL, QBTZ-043, CLM1, CLM2, QM2 (0.75) and to V, VpBTZ-043, VM1, VM2, VpM2 (1). Median weight in the combined Stage 1 + Stage 2 analysis population was 54 kg (range 42-81 kg; Koele 2025 Table 1), so 70 kg is well above the observed median and the reference value is a convention rather than a cohort centre.",
      source_name        = "WT"
    ),
    FED = list(
      description        = "1 = the dose was swallowed into a stomach that already contained food (a standard breakfast started about 30 min earlier, or a high-fat breakfast); 0 = the dose was taken before any food (fasted, or immediately prior to the start of a standard breakfast).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "1 (dose taken 30 min after the start of a standard breakfast) is the model reference for relative bioavailability; the parameters gated on FED = 0 quantify the deviation from it.",
      notes              = "This is exactly the WITHFOOD variable of the supplementary control stream, IF(FOOD.EQ.3.OR.FOOD.EQ.1) WITHFOOD=1, where the control stream's FOOD codes are 0 = no food, 1 = high-fat food, 2 = standard food after dose (i.e. dose taken prior to the meal), 3 = standard food before dose (i.e. dose taken 30 min after the start of the meal). FED therefore separates 'food already present' (FOOD 1 and 3) from 'no food present at the moment of dosing' (FOOD 0 and 2). It gates three effects: the secondary lag-time absorption route exists only when FED = 1 (the parallel route was not identifiable without food); mean transit time is multiplied by e_fed_mtt = 0.360 when FED = 0; and the apparent M1 clearance and volume are divided by e_fed_fm_m1 = 1.39 when FED = 0. Per dose record, not per subject: Stage 1 participants dosed fasted on Days 1-12 (FED = 0) and after a high-fat breakfast on Day 14 (FED = 1).",
      source_name        = "FOOD"
    ),
    FED_HIGHFAT = list(
      description        = "1 = the dose was taken with a high-fat breakfast; 0 = any other prandial state.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (standard breakfast started 30 min before the dose, dose taken prior to a standard breakfast, or fasted)",
      notes              = "Control-stream FOOD = 1, gating IF(FOOD.EQ.1) F_highfat = THETA(9) = 1.41. Applies to the Stage 1 Day 14 occasion, on which participants received BTZ-043 with a high-fat breakfast after 12 days of fasted dosing (Koele 2025 Methods 'Clinical study'). Koele 2025 does not print the kcal / percent-fat composition of the high-fat breakfast. Mutually exclusive with the standard-breakfast arms; a FED_HIGHFAT = 1 record also carries FED = 1.",
      source_name        = "FOOD"
    ),
    FASTED_STRICT = list(
      description        = "1 = the dose was taken with no food administered on either side of it (the Stage 1 Days 1-12 fasted occasions); 0 = any less strict prandial state, i.e. the dose was taken immediately before a standard breakfast, or with food already present.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (relaxed fast or fed)",
      notes              = "Control-stream FOOD = 0, gating IF(FOOD.EQ.0) F_nofood = THETA(10) = 0.458. The two protocols this indicator separates, verbatim from Koele 2025 Methods 'Clinical study': FASTED_STRICT = 1 is 'participants received BTZ-043 in a fasted state in Stage 1 from Day 1 to Day 12'; FASTED_STRICT = 0 with FED = 0 is the Stage 2 arm in which 'participants took BTZ-043 either prior to ... the start of intake of a standard breakfast', so food arrives shortly after the dose. Together with FED this spans the register's documented three-level prandial factor: strict fast (FASTED_STRICT = 1, FED = 0), relaxed fast (FASTED_STRICT = 0, FED = 0), fed (FASTED_STRICT = 0, FED = 1); FED_HIGHFAT then refines the fed level. Per dose record.",
      source_name        = "FOOD"
    ),
    DOSE_HIGH = list(
      description        = "1 = the administered daily BTZ-043 dose exceeded 1250 mg (the 1500 and 1750 mg Stage 1 cohorts); 0 = 1250 mg or below.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (daily dose of 1250 mg or less)",
      notes              = "Threshold for this model: DOSE > 1250 mg, per the control stream IF(DOSE.GT.1250)FDose = THETA(8) = 0.710. Koele 2025 Results: 'An apparent plateau in bioavailability was observed for doses over 1250 mg daily. No further decrease in the bioavailability was detected between the 1500 and 1750 mg doses (dOFV = -3.3, df = 1)', attributed in the Discussion to possible saturable absorption. Only Stage 1 escalated above 1250 mg, and those cohorts were dosed fasted, but the control stream applies the effect unconditionally on dose, so the encoding here does the same. Per dose record.",
      source_name        = "DOSE"
    ),
    RACE_COLOURED = list(
      description        = "1 = participant self-identified as Coloured (the paper's term is 'Cape-coloured', the South African population group common in the South-Western Cape); 0 = any other self-identified race.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (Black or White in this cohort)",
      notes              = "Control-stream RACE, documented in its own data dictionary as 'Race [1=Cape-colored, 0=other]', gating IF(RACE.EQ.1)CLeffRace = THETA(22) = 0.762. Koele 2025 Results: 'BTZ-043 clearance was estimated to be 24% (95% CI: 12%-35%) slower in Cape-coloured participants'; the Discussion offers a genetic polymorphism in a BTZ-043-metabolizing enzyme as a hypothesis and states the effect 'should be confirmed in a larger population'. Cohort composition (Koele 2025 Table 1, combined Stage 1 + 2): Black 44/68 (65%), Cape-coloured 23/68 (34%), White 1/68 (1%). Time-fixed per subject.",
      source_name        = "RACE"
    ),
    OCC = list(
      description        = "Intensive-PK sampling-occasion index used for the between-occasion random effects on relative bioavailability, mean transit time and the absorption rate constant.",
      units              = "(index)",
      type               = "categorical",
      reference_category = NULL,
      notes              = "Three occasions, matching the control stream's OCC1 / OCC2 / OCC3 multiplexers: OCC = 1 is Day 1, OCC = 2 is Day 12 of Stage 1, OCC = 3 is Day 14 (Koele 2025 Methods 'PK analysis': intensive sampling on Days 1, 12 and 14 for Stage 1 and on Days 1 and 14 for Stage 2). Koele 2025 Methods 'PK model development': 'Each PK sampling day was treated as a separate occasion to specify the inter-occasion variability (IOV).' For simulation, set OCC to the occasion index of the dosing interval being simulated; a single-occasion simulation may use OCC = 1 throughout.",
      source_name        = "DAY"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 68L,
    n_studies      = 1L,
    age_range      = "18-57 years",
    age_median     = "27 years",
    weight_range   = "42-81 kg",
    weight_median  = "54 kg",
    height_range   = "1.5-1.9 m",
    height_median  = "1.7 m",
    sex_female_pct = 16.2,
    race_ethnicity = c(
      Black            = 64.7,
      `Cape-coloured`  = 33.8,
      White            = 1.5
    ),
    hiv_status     = "HIV-1 negative, 68/68 (100%)",
    disease_state  = "Adults aged 18-64 years with drug-susceptible pulmonary tuberculosis enrolled in the sequential Phase 1b/2a dose-escalation and dose-expansion trial NCT04044001 (Stage 1 Phase 1b dose escalation, Stage 2 Phase 2a randomized dose expansion).",
    dose_range     = "Oral BTZ-043 250-1750 mg once daily for 14 days. Stage 1 escalated through 250, 500, 750, 1000, 1250, 1500 and 1750 mg with three participants per cohort and six in the highest; Stage 2 randomized 54 participants to 250, 500 or 1000 mg daily or to the Rifafour e-275 control regimen in a 3:3:3:2 ratio.",
    prandial_state = "Stage 1: fasted on Days 1-12 and with a high-fat breakfast on Day 14. Stage 2: BTZ-043 taken either prior to, or 30 min after, the start of intake of a standard breakfast. Occasion counts in Koele 2025 Table 1: fasted 24, high-fat 19, standard together with dose 24, standard 30 min prior to dose 20.",
    regions        = "South Africa (TASK, Cape Town; University of Cape Town Lung Institute)",
    notes          = "Baseline demographics from Koele 2025 Table 1, 'Combined Stage 1 + 2' column. Twenty-four Stage 1 and 44 Stage 2 participants were analysed; one Stage 2 participant in the 500 mg group withdrew before any study procedure and contributed no data. The PK dataset held 1808 BTZ-043, 1808 M1 and 1793 M2 observations, of which 600 (33%), 159 (8.8%) and 146 (8.1%) respectively were below the 20 ng/mL limit of quantification. BLQ observations were EXCLUDED from the final analysis: models using the NONMEM M3 method were highly unstable and, when the final model was re-estimated with M3, it significantly underestimated the central tendency of the BTZ-043 and M2 elimination phases (Koele 2025 Results, 'PK model')."
  )

  ini({
    # -----------------------------------------------------------------------
    # Structural parameters, Koele 2025 Table 2 'Structural parameters',
    # cross-checked line by line against the $THETA block of the final
    # control stream in the Supplementary data. Every clearance and volume is
    # an apparent value: Koele 2025 Methods 'PK model development' states
    # 'To keep the model parameters identifiable, the bioavailability and
    # fractions metabolized per path were set to 1, rendering disposition
    # parameters relative to the true values.'
    # -----------------------------------------------------------------------
    lcl     <- log(404)     ; label("Apparent BTZ-043 clearance CL/F at 70 kg (L/h)")                                   # Koele 2025 Table 2 'CL/F 404 (327-479)'; control stream THETA(1)
    lvc     <- log(764)     ; label("Apparent BTZ-043 central volume V/F at 70 kg (L)")                                 # Koele 2025 Table 2 'V/F 764 (303-953)'; control stream THETA(2)
    lka     <- log(1.69)    ; label("First-order absorption rate constant ka from the absorption compartment and from the secondary depot (1/h)") # Koele 2025 Table 2 'ka 1.69 (1.35-2.21)'; control stream THETA(3)
    lmtt    <- log(0.340)   ; label("Mean transit time MTT with food (h)")                                              # Koele 2025 Table 2 'MTT (h) 0.340 (0.270-0.419)'; control stream THETA(4)
    lq      <- log(68.0)    ; label("Apparent BTZ-043 inter-compartmental clearance Q/F at 70 kg (L/h)")                # Koele 2025 Table 2 'QBTZ-043/F 68.0 (55.0-81.3)'; control stream THETA(24)
    lvp     <- log(382)     ; label("Apparent BTZ-043 peripheral volume Vp/F at 70 kg (L)")                             # Koele 2025 Table 2 'VpBTZ-043/F 382 (310-449)'; control stream THETA(25)
    lcl_m1  <- log(36.0)    ; label("Apparent M1 clearance CL_M1/(F*fm_M1) at 70 kg, with food (L/h)")                   # Koele 2025 Table 2 'CLM1/(F.fmM1) 36.0 (29.3-43.2)'; control stream THETA(13)
    lvc_m1  <- log(894)     ; label("Apparent M1 central volume V_M1/(F*fm_M1) at 70 kg, with food (L)")                 # Koele 2025 Table 2 'VM1/(F.fmM1) 894 (699-1100)'; control stream THETA(14)
    lcl_m2  <- log(45.8)    ; label("Apparent M2 clearance CL_M2/(F*fm_M2) at 70 kg, first 10 days of treatment (L/h)")  # Koele 2025 Table 2 'CLM2/(F.fmM2) 45.8 (37.8-54.4)'; control stream THETA(16)
    lvc_m2  <- log(19.0)    ; label("Apparent M2 central volume V_M2/(F*fm_M2) at 70 kg (L)")                            # Koele 2025 Table 2 'VM2/(F.fmM2) 19.0 (13.4-25.8)'; control stream THETA(17)
    lq_m2   <- log(64.2)    ; label("Apparent M2 inter-compartmental clearance Q_M2/(F*fm_M2) at 70 kg (L/h)")           # Koele 2025 Table 2 'QM2/(F.fmM2) 64.2 (38.6-99.5)'; control stream THETA(19)
    lvp_m2  <- log(26.8)    ; label("Apparent M2 peripheral volume Vp_M2/(F*fm_M2) at 70 kg (L)")                        # Koele 2025 Table 2 'VpM2/(F.fmM2) 26.8 (21.4-33.2)'; control stream THETA(20)

    # Parallel-absorption split. TVF1 is the fraction of the bioavailable
    # dose that travels the transit chain; the remaining 1 - TVF1 = 0.356
    # travels the lag-time secondary route. The split applies only when the
    # dose is taken with food; without food TVF1 = 1 and the secondary route
    # carries nothing (control stream: TVF1 = 1; IF(WITHFOOD.EQ.1)TVF1 = THETA(6)).
    lfdepot <- log(0.644)   ; label("Fraction of the bioavailable dose absorbed through the transit chain when the dose is taken with food (unitless)") # Koele 2025 Table 2 'Fraction parallel absorption with food (%) 0.644 (0.606-0.679)'; control stream THETA(6); the complement 0.356 matches Results '36% (95% CI: 32%-39%) of the total bioavailable dose ... absorbed through the secondary delayed absorption route'
    ltlag   <- log(1.83)    ; label("Lag time of the secondary parallel absorption route (h)")                          # Koele 2025 Table 2 'Lag time parallel absorption with food (h) 1.83 (1.76-1.90)'; control stream THETA(7)

    # Allometric exponents, fixed at the theory-based values.
    e_wt_cl <- fixed(0.75)  ; label("Allometric exponent on every apparent clearance (CL/F, Q/F, CL_M1, CL_M2, Q_M2; unitless)") # Koele 2025 Methods 'PK model development': 'allometrically scaled, using fixed exponents of 1 and 0.75'; control stream AlloCL = (WT/70)**0.75
    e_wt_vc <- fixed(1)     ; label("Allometric exponent on every apparent volume (V/F, Vp/F, V_M1, V_M2, Vp_M2; unitless)")     # Koele 2025 Methods 'PK model development'; control stream AlloV = (WT/70)**1

    # -----------------------------------------------------------------------
    # Covariate effects, Koele 2025 Table 2 'Covariates'. All are
    # multiplicative factors expressed as a percentage of the reference in
    # the source table; the reference prandial state is a standard breakfast
    # started about 30 min before the dose.
    # -----------------------------------------------------------------------
    e_dose_high_fdepot      <- 0.710 ; label("Multiplicative factor on relative bioavailability for a daily dose above 1250 mg (unitless)")            # Koele 2025 Table 2 'Dose (>1250 mg) on F (%) 71.0 (51.0-98.0)'; control stream THETA(8)
    e_fed_highfat_fdepot    <- 1.41  ; label("Multiplicative factor on relative bioavailability for a high-fat breakfast (unitless)")                  # Koele 2025 Table 2 'High-fat food on F (%) 141 (105-178)'; control stream THETA(9); Results '41% (95% CI: 5.0%-78%)' higher
    e_fasted_strict_fdepot  <- 0.458 ; label("Multiplicative factor on relative bioavailability when the dose is taken without any food (unitless)")   # Koele 2025 Table 2 'No food on F (%) 45.8 (35.2-56.9)'; control stream THETA(10); Results 'decreased by 54% (95% CI: 43%-65%)'
    e_fasted_relaxed_fdepot <- 0.727 ; label("Multiplicative factor on relative bioavailability when the dose is taken immediately prior to a standard breakfast (unitless)") # Koele 2025 Table 2 'Dose prior to standard food on F (%) 72.7 (58.3-88.3)'; control stream THETA(11); Results '27% (95% CI: 12%-42%) lower'
    e_fed_mtt               <- 0.360 ; label("Multiplicative factor on mean transit time when no food is present at dosing, i.e. FED = 0 (unitless)")  # Koele 2025 Table 2 'Dose prior to standard food/no food on MTT (%) 36.0 (27.6-46.8)'; control stream THETA(5), applied under IF(WITHFOOD.EQ.0)
    e_fed_fm_m1             <- 1.39  ; label("Divisor applied to the apparent M1 clearance and volume when no food is present at dosing, i.e. FED = 0 (unitless)") # Koele 2025 Table 2 'Administration with food on FM1 (%) 139 (134-144)'; control stream THETA(23), applied under IF(WITHFOOD.EQ.0) as CLM1 = TVCLM1/FM1FOOD and VM1 = TVVM1/FM1FOOD
    e_race_coloured_cl      <- 0.762 ; label("Multiplicative factor on apparent BTZ-043 clearance in Cape-coloured participants (unitless)")           # Koele 2025 Table 2 'Cape-coloured race on CL/F (%) 76.2 (65.4-88.5)'; control stream THETA(22)
    e_time_cl_m2            <- 0.734 ; label("Multiplicative factor on apparent M2 clearance after 10 days on treatment (unitless)")                    # Koele 2025 Table 2 'Time effect on CLM2/(F.fmM2) (%) 73.4 (70.8-76.1)'; control stream THETA(21), applied under IF(DAY.GT.10)

    # Breakpoint of the piecewise-constant M2 clearance. The control stream
    # reads a DAY data item and switches under IF(DAY.GT.10); with Day 1
    # spanning 0-24 h, DAY > 10 is 240 h or more after the start of
    # treatment. Carried as a fixed structural constant so the breakpoint is
    # discoverable by name rather than buried in the model block.
    ltclchange_m2 <- fixed(log(240)) ; label("Time from the start of treatment at which apparent M2 clearance steps down (h)") # Koele 2025 Results: 'The clearance of M2 decreased by 27% (95% CI: 24%-29%) after 10 days on BTZ-043 treatment. The effect was implemented as a dichotomous step after 10 days'; control stream IF(DAY.GT.10)

    # -----------------------------------------------------------------------
    # Between-subject variability. Koele 2025 Table 2 reports CV% and the
    # table footnote gives the transform used: 'CV% was calculated as
    # sqrt(e^OM2 - 1)'. The variances below are the $OMEGA values of the
    # supplementary control stream; each was verified to reproduce the
    # printed CV% to the digits given (for example CL:
    # sqrt(exp(0.0715) - 1) = 0.2723 = 27.2%).
    # -----------------------------------------------------------------------
    etalcl    ~ 0.0715  # Koele 2025 Table 2 CL/F CV 27.2% (21.0-35.4); control stream $OMEGA 1
    etalvc    ~ 0.172   # Koele 2025 Table 2 V/F CV 43.3% (33.3-55.8); control stream $OMEGA 2
    etalfdepot ~ 0.0681 # Koele 2025 Table 2 F IIV 26.5% (19.3-35.7); control stream $OMEGA 12
    etalcl_m1 ~ 0.174   # Koele 2025 Table 2 CLM1 CV 43.6% (33.2-56.9); control stream $OMEGA 13
    etalvc_m1 ~ 0.212   # Koele 2025 Table 2 VM1 CV 48.6% (39.0-60.8); control stream $OMEGA 14
    etalcl_m2 ~ 0.0878  # Koele 2025 Table 2 CLM2 CV 30.3% (24.6-36.0); control stream $OMEGA 15
    etalvc_m2 ~ 0.77    # Koele 2025 Table 2 VM2 CV 108% (73.3-180); control stream $OMEGA 16

    # -----------------------------------------------------------------------
    # Between-occasion variability over the three intensive-PK sampling
    # occasions. Relative bioavailability carries IOV on top of its IIV;
    # mean transit time and ka carry IOV ONLY (the control stream has no IIV
    # etas on either), estimated as a 2x2 block whose correlation reproduces
    # the 44.4% printed in Table 2:
    # 0.345 / sqrt(0.621 * 0.972) = 0.444.
    # $OMEGA BLOCK(1) SAME and BLOCK(2) SAME repeats carry the occasion-1
    # values forward to occasions 2 and 3.
    # -----------------------------------------------------------------------
    etaiov_fdepot_1 ~ 0.0543       # Koele 2025 Table 2 F IOV 23.6% (20.0-27.4); control stream $OMEGA BLOCK(1) 3
    etaiov_fdepot_2 ~ fix(0.0543)  # control stream $OMEGA BLOCK(1) SAME (ETA 4)
    etaiov_fdepot_3 ~ fix(0.0543)  # control stream $OMEGA BLOCK(1) SAME (ETA 5)
    etaiov_mtt_1 + etaiov_ka_1 ~ c(0.621,
                                   0.345, 0.972)  # Koele 2025 Table 2 MTT CV 92.8% (76.3-116), ka CV 128% (97.6-173), 'Correlation MTT-ka (%) 44.4 (27.8-54.9)'; control stream $OMEGA BLOCK(2) 6-7
    etaiov_mtt_2 + etaiov_ka_2 ~ c(0.621,
                                   0.345, 0.972)  # control stream $OMEGA BLOCK(2) SAME (ETA 8-9)
    etaiov_mtt_3 + etaiov_ka_3 ~ c(0.621,
                                   0.345, 0.972)  # control stream $OMEGA BLOCK(2) SAME (ETA 10-11)

    # -----------------------------------------------------------------------
    # Residual error. The control stream builds W = err * IPRED with
    # $SIGMA 1 FIX, i.e. a pure proportional error whose SD is the THETA.
    # The values below are the FINAL estimates printed in Koele 2025
    # Table 2; the supplementary control stream's $THETA initials for these
    # three (0.477, 0.329, 0.285) are the only entries in that block that do
    # not equal the published final estimates, so Table 2 is used. See the
    # vignette Errata section.
    # -----------------------------------------------------------------------
    propSd    <- 0.505 ; label("BTZ-043 proportional residual error (fraction)") # Koele 2025 Table 2 'Proportional error BTZ-043 (CV%) 50.5 (47.7-53.4)'
    propSd_m1 <- 0.338 ; label("M1 proportional residual error (fraction)")      # Koele 2025 Table 2 'Proportional error M1 (CV%) 33.8 (32.4-35.1)'
    propSd_m2 <- 0.290 ; label("M2 proportional residual error (fraction)")      # Koele 2025 Table 2 'Proportional error M2 (CV%) 29.0 (28.0-30.1)'
  })

  model({
    # 1. Allometric scaling on 70 kg.
    allcl <- (WT / 70)^e_wt_cl
    allv  <- (WT / 70)^e_wt_vc

    # 2. Occasion multiplexer for the between-occasion random effects
    #    (control stream OCC1 / OCC2 / OCC3).
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    iov_fdepot <- oc1 * etaiov_fdepot_1 + oc2 * etaiov_fdepot_2 +
      oc3 * etaiov_fdepot_3
    iov_mtt    <- oc1 * etaiov_mtt_1 + oc2 * etaiov_mtt_2 + oc3 * etaiov_mtt_3
    iov_ka     <- oc1 * etaiov_ka_1  + oc2 * etaiov_ka_2  + oc3 * etaiov_ka_3

    # 3. Prandial state. FED = 1 is the reference (standard breakfast started
    #    about 30 min before the dose). The relaxed-fast arm - the dose taken
    #    immediately prior to a standard breakfast - is the one prandial state
    #    with no column of its own; it is the complement of the other two
    #    indicators, exactly as in the control stream where FOOD = 2 is the
    #    only code with WITHFOOD = 0 and F_nofood = 1.
    fastrelax <- (1 - FED) * (1 - FASTED_STRICT)

    # 4. Relative bioavailability multiplier and the mg -> nmol conversion.
    #    BTZ-043 molecular weight 431.39 g/mol, so 1 mg = 1e6/431.39 nmol
    #    = 2318 nmol (control stream comment on F1: 'DV in nmol/L. AMT in mg.
    #    Mw = 431.39 g/mol. To get dose in nmol = 10^6/431.39').
    frel <- e_fed_highfat_fdepot^FED_HIGHFAT *
      e_fasted_strict_fdepot^FASTED_STRICT *
      e_fasted_relaxed_fdepot^fastrelax *
      e_dose_high_fdepot^DOSE_HIGH
    nmolpermg <- 1000000 / 431.39

    # 5. Individual parameters. The transit-chain fraction is 0.644 with food
    #    and 1 without (exp(lfdepot * FED) is 0.644^FED).
    fdepot1 <- exp(lfdepot * FED)
    mtt <- exp(lmtt) * e_fed_mtt^(1 - FED) * exp(iov_mtt)
    ktr <- 1 / mtt
    ka  <- exp(lka) * exp(iov_ka)

    cl  <- exp(lcl + etalcl) * allcl * e_race_coloured_cl^RACE_COLOURED
    vc  <- exp(lvc + etalvc) * allv
    q   <- exp(lq)  * allcl
    vp  <- exp(lvp) * allv

    # M1: the food effect divides BOTH the apparent clearance and the
    #     apparent volume, so the M1 half-life is unchanged and only the M1
    #     exposure moves (control stream CLM1 = TVCLM1/FM1FOOD,
    #     VM1 = TVVM1/FM1FOOD, with FM1FOOD = 1.39 when WITHFOOD = 0).
    fm1food <- e_fed_fm_m1^(1 - FED)
    cl_m1 <- exp(lcl_m1 + etalcl_m1) / fm1food * allcl
    vc_m1 <- exp(lvc_m1 + etalvc_m1) / fm1food * allv

    # M2: apparent clearance is piecewise constant, stepping down after 10
    #     days on treatment. The control stream reads a DAY data item and
    #     applies the factor under IF(DAY.GT.10); with Day 1 spanning 0-24 h,
    #     DAY > 10 is t >= 240 h after the start of treatment, which is how the
    #     breakpoint is derived here so that no extra data column is needed.
    #     Written in the canonical tclchange / cl_late two-arm form; note the
    #     paper reports the LATE arm as a ratio to the early one (73.4%), not
    #     as a second absolute clearance, so e_time_cl_m2 carries the printed
    #     value and cl_late_m2 is derived from it.
    tclchange_m2 <- exp(ltclchange_m2)
    cl_early_m2  <- exp(lcl_m2 + etalcl_m2) * allcl
    cl_late_m2   <- cl_early_m2 * e_time_cl_m2
    cl_m2 <- cl_early_m2 * (t < tclchange_m2) + cl_late_m2 * (t >= tclchange_m2)
    vc_m2 <- exp(lvc_m2 + etalvc_m2) * allv
    q_m2  <- exp(lq_m2)  * allcl
    vp_m2 <- exp(lvp_m2) * allv

    # 6. ODE system, transcribed one line per DADT from the supplementary
    #    $DES block. Compartment map: depot1 = DEPOT (1), depot2 = DEPOT2 (2),
    #    central = BTZ043 (3), transit1 = TRANS1 (4), transit2 = TRANS2 (5),
    #    central_m1 = M1 (6), central_m2 = M2 (7),
    #    peripheral1_m2 = M2peripheral (8), peripheral1 = BTZ043peripheral (9).
    #
    #    Note that the FULL parent elimination flux cl * central / vc feeds
    #    BOTH metabolite compartments. That is the direct consequence of
    #    fixing every fraction metabolised to 1: each metabolite's apparent
    #    clearance and volume absorb its own fm, so the metabolite mass
    #    balance is relative rather than absolute. Do not "fix" this by
    #    splitting the flux - it would change every published estimate.
    d/dt(depot1)   <- -ktr * depot1
    d/dt(transit1) <-  ktr * depot1 - ktr * transit1
    d/dt(transit2) <-  ktr * transit1 - ka * transit2
    d/dt(depot2)   <- -ka * depot2
    d/dt(central)  <-  ka * transit2 + ka * depot2 -
      cl * central / vc -
      q * central / vc + q * peripheral1 / vp
    d/dt(peripheral1) <- q * central / vc - q * peripheral1 / vp
    d/dt(central_m1)  <- cl * central / vc - cl_m1 * central_m1 / vc_m1
    d/dt(central_m2)  <- cl * central / vc - cl_m2 * central_m2 / vc_m2 -
      q_m2 * central_m2 / vc_m2 + q_m2 * peripheral1_m2 / vp_m2
    d/dt(peripheral1_m2) <- q_m2 * central_m2 / vc_m2 -
      q_m2 * peripheral1_m2 / vp_m2

    # 7. Bioavailability and lag time. Both depots must receive a dose record
    #    for the same administration; f() then splits the amount between the
    #    transit chain and the delayed secondary route and converts mg to
    #    nmol. Absolute bioavailability was fixed to 1 in the source
    #    analysis, so frel is a RELATIVE factor against the standard-breakfast
    #    reference (control stream F1 and F2).
    f(depot1) <- fdepot1 * nmolpermg * frel *
      exp(iov_fdepot + etalfdepot)
    f(depot2) <- (1 - fdepot1) * nmolpermg * frel *
      exp(iov_fdepot + etalfdepot)
    alag(depot2) <- exp(ltlag)

    # 8. Observations. Amount (nmol) / volume (L) = nmol/L. Cc_total is the
    #    BTZ-043total analyte of the paper: the bioanalytical method converts
    #    all M2 back to BTZ-043 before quantification, so the measured
    #    BTZ-043total is the molar sum of parent and M2. It carries no
    #    residual error of its own here because the paper fitted BTZ-043, M1
    #    and M2 separately; Cc_total is provided because it is the exposure
    #    metric that drives the exposure-response model in
    #    Koele_2025_btz043_bacterialload.
    Cc       <- central / vc
    Cc_m1    <- central_m1 / vc_m1
    Cc_m2    <- central_m2 / vc_m2
    Cc_total <- Cc + Cc_m2

    Cc    ~ prop(propSd)
    Cc_m1 ~ prop(propSd_m1)
    Cc_m2 ~ prop(propSd_m2)
  })
}
