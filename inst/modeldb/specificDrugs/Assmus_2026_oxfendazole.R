Assmus_2026_oxfendazole <- function() {
  description <- "Joint parent + two-metabolite population PK model for oral oxfendazole (immediate-release tablet, fasted) and its metabolites fenbendazole and oxfendazole sulfone in healthy African adults (Assmus 2026). First-order absorption with first-pass metabolism: a logit-scale fraction of the absorbed dose is converted to oxfendazole sulfone before reaching the systemic circulation, the remainder enters the oxfendazole central compartment. One-compartment disposition for each analyte with linear elimination; oxfendazole clearance is split equally (fraction fixed at 0.5 for identifiability) between formation of fenbendazole and formation of oxfendazole sulfone. Dose-limited relative bioavailability F = (dose / 34.7 mg)^-0.541, coefficients fixed from an earlier single-ascending-dose study. Allometric body-weight scaling (0.75 on clearances, 1 on volumes, reference 70 kg) on all six disposition parameters. Inter-occasion variability on F and KA across up to five dosing occasions. The model is on the authors' molar scale (dose in nmol, concentrations in nmol/L); metabolite clearances and volumes are apparent values relative to the unknown true fraction metabolised."
  reference <- paste(
    "Assmus F, Adehin A, Hoglund RM, Nyaulingo G, Mbarak H, Jongo S,",
    "Ackermann E, Reus E, Keiser J, Barreira Da Silva Rocha F, Specht S,",
    "Scandale I, Tarning J.",
    "Repurposing Oxfendazole for Onchocerciasis: Population Pharmacokinetics",
    "of a Tablet Formulation in Healthy African Adults.",
    "CPT Pharmacometrics Syst Pharmacol. 2026. doi:10.1002/psp4.70189.",
    "Structure and parameter values taken from the final NONMEM control",
    "stream (Supplementary Code S1), whose $THETA equals Table 1 of the main",
    "paper; the IIV on KA follows Table 1 (see vignette Errata).",
    sep = " "
  )
  vignette <- "Assmus_2026_oxfendazole"
  units <- list(time = "h", dosing = "nmol", concentration = "nmol/L")

  # The source model is fitted on the molar scale (Code S1 $INPUT: 'AMT ;
  # Dose amount (nmol)', 'DV ; ln of dependent variable (nmol/L)'), so every
  # transfer between analytes is mole-for-mole. verified = TRUE: analyte and
  # specimen are read off the $MODEL block comments of Code S1.
  compartmentData <- list(
    depot = list(analyte = "oxfendazole", units = "nmol", specimen = "administration site", verified = TRUE),
    central = list(analyte = "oxfendazole", units = "nmol", specimen = "plasma", verified = TRUE),
    central_fenbendazole = list(analyte = "fenbendazole", units = "nmol", specimen = "plasma", verified = TRUE),
    central_oxfsfn = list(analyte = "oxfendazole sulfone", units = "nmol", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Baseline. Allometric scaling implemented a priori on all clearances (exponent 0.75) and volumes (exponent 1) of the three analytes, reference 70 kg (Methods 2.4; Code S1 $PK '(WT/70)**0.75' and '(WT/70)**1.00'). Analysis cohort median 58 kg (42-75), Table S2. Source column WT.",
      source_name = "WT"
    ),
    DOSE_OXFENDAZOLE_MG = list(
      description = "Administered oxfendazole dose of the current dosing record, in mg",
      units = "mg",
      type = "continuous",
      reference_category = NULL,
      notes = "Drives the dose-limited relative bioavailability F = (DOSE_OXFENDAZOLE_MG / 34.7)^-0.541 (Eq. 1; Code S1 'COV = (DOSE/34.7)**(-0.541)'). This is the per-administration dose, not the daily dose: the paper's twice-daily simulations rely on splitting the daily dose to reduce the dose-limited absorption (Discussion 4.2). Must not be named DOSE: a covariate column literally named DOSE is consumed by rxode2's event translation and never reaches model(). Because amt is in nmol, supply the mg dose separately in this column (amt nmol = mg * 1e6 / 315.35 for oxfendazole, C15H13N3O3S). Studied doses 100 and 400 mg. Source column DOSE.",
      source_name = "DOSE"
    ),
    OCC = list(
      description = "Integer-valued dosing-occasion indicator for inter-occasion-variability multiplexing.",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = "Each drug administration is a separate occasion (Methods 2.4). Values 1-5: cohorts 1 and 2 (single dose) carry OCC = 1; cohort 3 (400 mg once daily for 5 days) carries OCC = 1..5 for doses 1..5. Decomposed inside model() into oc1 .. oc5, which multiplex the five IOV etas on F and on KA exactly as Code S1 builds IOV1 and IOV2 from ETA(11)-ETA(20). An OCC value outside 1-5 switches every indicator off, so that occasion carries no IOV. Source column OCC.",
      source_name = "OCC"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 24,
    n_studies = 1,
    age_range = "19-41 years",
    age_median = "28 years",
    weight_range = "42-75 kg",
    weight_median = "58 kg",
    height_median = "161.5 cm (145.5-182.5)",
    bmi_median = "21.4 kg/m^2 (18.3-28.9)",
    sex_female_pct = 41.7,
    race_ethnicity = "100% Black African",
    disease_state = "Healthy adults (aged 18-45, BMI 18.0-29.9 kg/m^2, no evidence of parasitic or chronic infection).",
    dose_range = "Oxfendazole immediate-release tablet, fasted: 100 mg single dose, 400 mg single dose, or 400 mg once daily for 5 days (8 subjects per cohort).",
    regions = "Tanzania (Bagamoyo Clinical Trial Facility); HELP-OFZ Phase 1 study, NCT04920292, conducted 2022.",
    notes = "30 enrolled; 6 placebo recipients excluded from the PK analysis (Results 3.1). 984 post-dose samples, 328 per analyte; 8.4% below the LLOQ (2 ng/mL oxfendazole, 1 ng/mL each metabolite) and omitted. Demographics from Table S2 and Results 3.1."
  )

  ini({
    # Structural parameters. Code S1 $THETA equals Table 1's population
    # estimates digit-for-digit, so these are final estimates. Clearances
    # and volumes are apparent (relative to F) and for a 70-kg adult.
    lcl <- log(3.54); label("Apparent oxfendazole clearance CL_OXF/F at WT = 70 kg (L/h)") # Table 1 'CL OXF/F (L/h) 3.54'; $THETA 1
    lvc <- log(67.2); label("Apparent oxfendazole volume of distribution V_OXF/F at WT = 70 kg (L)") # Table 1 'V OXF/F (L) 67.2'; $THETA 2
    lka <- log(0.695); label("First-order absorption rate constant KA (1/h)") # Table 1 'K A (h-1) 0.695'; $THETA 3
    lfdepot <- fixed(log(1)); label("Population relative oral bioavailability F at the 34.7-mg reference dose (unitless)") # Table 1 'F 1 fixed'; $THETA 4 '(1) FIX'
    lcl_fenbendazole <- log(94.4); label("Apparent fenbendazole clearance CL_FEN/F at WT = 70 kg (L/h)") # Table 1 'CL FEN/F (L/h) 94.4'; $THETA 5
    lvc_fenbendazole <- log(3370); label("Apparent fenbendazole volume of distribution V_FEN/F at WT = 70 kg (L)") # Table 1 'V FEN/F (L) 3370'; $THETA 6
    fm_fenbendazole <- fixed(0.5); label("Fraction of systemic oxfendazole clearance forming fenbendazole; the remainder forms oxfendazole sulfone (unitless)") # Table 1 'f M 0.5 fixed'; $THETA 7 '(0.5) FIX'
    lcl_oxfsfn <- log(15.8); label("Apparent oxfendazole sulfone clearance CL_OXF-SO2/F at WT = 70 kg (L/h)") # Table 1 'CL OXF-SO2/F (L/h) 15.8'; $THETA 8
    lvc_oxfsfn <- log(206); label("Apparent oxfendazole sulfone volume of distribution V_OXF-SO2/F at WT = 70 kg (L)") # Table 1 'V OXF-SO2/F (L) 206'; $THETA 9
    logitfm_oxfsfn <- log(0.0549 / (1 - 0.0549)); label("Logit of the fraction of absorbed oxfendazole converted to oxfendazole sulfone by first-pass metabolism, f_FP = 0.0549 (unitless)") # $THETA 10 '(0, 0.0549, 1)'; Table 1 'f FP 0.055'

    # Covariate effects, all fixed.
    e_wt_cl <- fixed(0.75); label("Allometric exponent on all three clearances (unitless)") # Methods 2.4 (fixed a priori); $PK '(WT/70)**0.75'
    e_wt_vc <- fixed(1); label("Allometric exponent on all three volumes (unitless)") # Methods 2.4 (fixed a priori); $PK '(WT/70)**1.00'
    e_dose_fdepot <- fixed(-0.541); label("Power exponent of dose on relative bioavailability (unitless)") # Eq. 1a; Table 1 last row '-0.541 fixed'; $PK 'COV = (DOSE/34.7)**(-0.541)'
    dose_ref_fdepot <- fixed(34.7); label("Reference dose of the bioavailability power function (mg)") # Eq. 1a; Table 1 last row '34.7 fixed'

    # Inter-individual variability. $OMEGA holds variances; Table 1 %CV =
    # 100 * sqrt(exp(omega^2) - 1).
    etalcl ~ 0.0161 # $OMEGA 1 IIV_CL; 12.74% CV = Table 1 'CL OXF/F 12.7'
    etalvc ~ 0.0387 # $OMEGA 2 IIV_V2; 19.86% CV = Table 1 'V OXF/F 19.9'
    etalka ~ 0.0924 # Table 1 'K A IIV 31.1% CV' gives 0.0924; $OMEGA 3 prints 0.924 (123% CV), a dropped leading zero -- see vignette Errata
    etalfdepot ~ 0.424 # $OMEGA 4 IIV_F1; 72.67% CV = Table 1 'F IIV 72.7'
    etalcl_fenbendazole ~ 0.618 # $OMEGA 5 IIV_CL_M1; 92.48% CV = Table 1 'CL FEN/F 92.5'
    etalvc_fenbendazole ~ 0.569 # $OMEGA 6 IIV_V_M1; 87.55% CV = Table 1 'V FEN/F 87.5'
    etalcl_oxfsfn ~ 0.0552 # $OMEGA 8 IIV_CL_M2; 23.82% CV = Table 1 'CL OXF-SO2/F 23.8'
    etalogitfm_oxfsfn ~ 0.192 # $OMEGA 10 IIV_FFP on the logit scale; Table 1 '44.2% CV' is by simulation (Code S2)

    # No IIV on f_M or V_OXF-SO2 ($OMEGA 7 and 9 are '0 FIX'; Table 1 '-').
    # Omitted rather than written as `~ fixed(0)` because a zero-variance
    # diagonal makes OMEGA singular and breaks the Cholesky sampler used by
    # rxSolve.

    # Inter-occasion variability over five dosing occasions. $OMEGA
    # BLOCK(1) with SAME on occasions 2-5; nlmixr2 has no SAME shortcut, so
    # the variance is fixed to the shared value after occasion 1.
    etaiov_fdepot_1 ~ 0.141 # $OMEGA 11 'BLOCK(1) 0.141 ; IOV F1'; 38.91% CV = Table 1 'F IOV 38.9'
    etaiov_fdepot_2 ~ fix(0.141) # $OMEGA 12 'BLOCK(1) SAME'
    etaiov_fdepot_3 ~ fix(0.141) # $OMEGA 13 'BLOCK(1) SAME'
    etaiov_fdepot_4 ~ fix(0.141) # $OMEGA 14 'BLOCK(1) SAME'
    etaiov_fdepot_5 ~ fix(0.141) # $OMEGA 15 'BLOCK(1) SAME'
    etaiov_ka_1 ~ 0.257 # $OMEGA 16 'BLOCK(1) 0.257 ; IOV KA'; 54.10% CV = Table 1 'K A IOV 54.1'
    etaiov_ka_2 ~ fix(0.257) # $OMEGA 17 'BLOCK(1) SAME'
    etaiov_ka_3 ~ fix(0.257) # $OMEGA 18 'BLOCK(1) SAME'
    etaiov_ka_4 ~ fix(0.257) # $OMEGA 19 'BLOCK(1) SAME'
    etaiov_ka_5 ~ fix(0.257) # $OMEGA 20 'BLOCK(1) SAME'

    # Residual error: additive on the natural-log scale for each analyte
    # (nlmixr2 `~ lnorm(expSd)`). $SIGMA holds variances; expSd = sqrt.
    expSd <- sqrt(0.0231); label("Residual SD on the natural-log scale, oxfendazole (log units)") # $SIGMA 1 = 0.0231; Table 1 'sigma OXF 0.152'
    expSd_fenbendazole <- sqrt(0.177); label("Residual SD on the natural-log scale, fenbendazole (log units)") # $SIGMA 2 = 0.177; Table 1 'sigma FEN 0.421'
    expSd_oxfsfn <- sqrt(0.0189); label("Residual SD on the natural-log scale, oxfendazole sulfone (log units)") # $SIGMA 3 = 0.0189; Table 1 'sigma OXF-SO2 0.137'
  })

  model({
    # ---- 1. Occasion indicators and IOV (Code S1 $PK) -----------------
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    oc5 <- (OCC == 5)
    iov_fdepot <- oc1 * etaiov_fdepot_1 + oc2 * etaiov_fdepot_2 +
      oc3 * etaiov_fdepot_3 + oc4 * etaiov_fdepot_4 + oc5 * etaiov_fdepot_5
    iov_ka <- oc1 * etaiov_ka_1 + oc2 * etaiov_ka_2 +
      oc3 * etaiov_ka_3 + oc4 * etaiov_ka_4 + oc5 * etaiov_ka_5

    # ---- 2. Individual parameters -------------------------------------
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    ka <- exp(lka + etalka + iov_ka)
    cl_fenbendazole <- exp(lcl_fenbendazole + etalcl_fenbendazole) * (WT / 70)^e_wt_cl
    vc_fenbendazole <- exp(lvc_fenbendazole + etalvc_fenbendazole) * (WT / 70)^e_wt_vc
    cl_oxfsfn <- exp(lcl_oxfsfn + etalcl_oxfsfn) * (WT / 70)^e_wt_cl
    vc_oxfsfn <- exp(lvc_oxfsfn) * (WT / 70)^e_wt_vc

    # Dose-limited relative bioavailability, Eq. 1a-1b:
    # F_i = F * (Dose / 34.7 mg)^-0.541 * exp(eta_F + IOV_F).
    fdepot <- exp(lfdepot + etalfdepot + iov_fdepot) *
      (DOSE_OXFENDAZOLE_MG / dose_ref_fdepot)^e_dose_fdepot

    # First-pass fraction on the logit scale (Code S1 FFP).
    logitfm_oxfsfn_ind <- logitfm_oxfsfn + etalogitfm_oxfsfn
    fm_oxfsfn <- expit(logitfm_oxfsfn_ind)

    # ---- 3. Micro-rate constants (ADVAN5, 1/h) ------------------------
    k_depot_central <- ka * (1 - fm_oxfsfn) # K12
    k_depot_oxfsfn <- ka * fm_oxfsfn # K14, first-pass formation
    kel <- cl / vc
    k_central_fenbendazole <- fm_fenbendazole * kel # K23
    k_central_oxfsfn <- (1 - fm_fenbendazole) * kel # K24
    kel_fenbendazole <- cl_fenbendazole / vc_fenbendazole # K30
    kel_oxfsfn <- cl_oxfsfn / vc_oxfsfn # K40

    # ---- 4. ODE system ------------------------------------------------
    d/dt(depot) <- -(k_depot_central + k_depot_oxfsfn) * depot
    d/dt(central) <- k_depot_central * depot -
      (k_central_fenbendazole + k_central_oxfsfn) * central
    d/dt(central_fenbendazole) <- k_central_fenbendazole * central -
      kel_fenbendazole * central_fenbendazole
    d/dt(central_oxfsfn) <- k_depot_oxfsfn * depot +
      k_central_oxfsfn * central -
      kel_oxfsfn * central_oxfsfn

    f(depot) <- fdepot

    # ---- 5. Observations and residual error (nmol/L) ------------------
    Cc <- central / vc
    Cc_fenbendazole <- central_fenbendazole / vc_fenbendazole
    Cc_oxfsfn <- central_oxfsfn / vc_oxfsfn

    Cc ~ lnorm(expSd)
    Cc_fenbendazole ~ lnorm(expSd_fenbendazole)
    Cc_oxfsfn ~ lnorm(expSd_oxfsfn)
  })
}
