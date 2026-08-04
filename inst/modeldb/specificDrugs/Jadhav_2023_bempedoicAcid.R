Jadhav_2023_bempedoicAcid <- function() {
  description <- "Population PK model for bempedoic acid (ATP-citrate lyase inhibitor) in healthy subjects and patients with dyslipidemia, renal or hepatic impairment, or type 2 diabetes mellitus (Jadhav 2023): two-compartment disposition with a single transit absorption compartment and linear elimination, pooled across 22 phase 1/2/3 studies. Covariate effects on apparent clearance (sex, body weight, Black race, hyperlipidemia, type 2 diabetes, eGFR, ezetimibe), on apparent central volume (sex, age, body weight, simvastatin), on the absorption rate constant (food), and on relative oral bioavailability (atorvastatin)."
  reference   <- "Jadhav SB, Amore BM, Bockbrader H, Crass RL, Chapel S, Sasiela WJ, Emery MG. Population pharmacokinetic and pharmacokinetic-pharmacodynamic modeling of bempedoic acid and low-density lipoprotein cholesterol in healthy subjects and patients with dyslipidemia. J Pharmacokinet Pharmacodyn. 2023;50(5):351-364. doi:10.1007/s10928-023-09864-w"
  vignette    <- "Jadhav_2023_bempedoicAcid"
  units       <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Derived mechanically; verified = FALSE means it has
  # NOT been checked against the source paper.
  compartmentData <- list(
    depot       = list(analyte = "bempedoicAcid", units = "mg", specimen = "administration site", verified = FALSE),
    transit1    = list(analyte = "bempedoicAcid", units = "mg", specimen = "administration site", verified = FALSE),
    central     = list(analyte = "bempedoicAcid", units = "mg", specimen = "plasma", verified = FALSE),
    peripheral1 = list(analyte = "bempedoicAcid", units = "mg", specimen = "plasma", verified = FALSE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight at baseline",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power covariate on CL/F (exponent 0.61) and on Vc/F (exponent 0.94), Jadhav 2023 Table 2. The paper states the covariate form as theta_TV = theta_REF * (x / x_REF)^theta_x (Methods, PK model covariate analysis) but never prints x_REF. The model file uses the popPK analysis-set median of 83.7 kg (Table 1) -- the standard normalisation for a pooled full-covariate popPK model with ESTIMATED (not fixed-allometric) exponents, and the choice that reconciles the model with the paper's own reported mean steady-state exposure of 12.5 ug/mL at 180 mg/day (see the vignette Assumptions and deviations section for the reconstruction and for the discarded 70 kg alternative). Only the absolute typical values of CL/F and Vc/F depend on this choice -- all published covariate RATIOS are invariant to it.",
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Age at baseline",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power covariate on Vc/F (exponent 0.743), Jadhav 2023 Table 2. Reference age not printed in the paper; the model file uses the popPK analysis-set median of 62.0 years (Table 1; mean 60.5, range 18-89 years), matching the median-normalisation convention applied to WT and CRCL.",
      source_name        = "AGE"
    ),
    SEXF = list(
      description        = "Biological sex indicator, 1 = female, 0 = male",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Proportional-shift covariate on CL/F (theta = -0.127, i.e. 12.7% lower apparent clearance in women) and on Vc/F (theta = -0.0895), Jadhav 2023 Table 2. The popPK analysis set was 40.6% female (Table 1). The paper's Fig. 2 forest plot reports a 1.39-fold (90% CI 1.34, 1.47) AUCss ratio for women vs. men; that marginal ratio also carries the correlated body-weight difference between the sexes and is therefore larger than the isolated 1/(1 - 0.127) = 1.15-fold sex effect encoded here.",
      source_name        = "SEXF"
    ),
    RACE_BLACK = list(
      description        = "Black / African American race indicator, 1 = Black, 0 = other race",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Black; 89.0% of the popPK analysis set was White)",
      notes              = "Proportional-shift covariate on CL/F (theta = -0.143), Jadhav 2023 Table 2, labelled 'Black race'. 9.2% of the popPK analysis set was Black (Table 1). Race in Table 1 is reported as separate White / Black / Asian / Native American / Native Hawaiian / Other groups, so the plain RACE_BLACK dichotomy (not the composite RACE_BLACK_OTH) is the correct canonical.",
      source_name        = "RACE_BLACK"
    ),
    DIS_HYPERLIP = list(
      description        = "Hyperlipidemia diagnosis indicator, 1 = participant has hyperlipidemia, 0 = no hyperlipidemia diagnosis",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no hyperlipidemia diagnosis)",
      notes              = "Proportional-shift covariate on CL/F (theta = -0.0945), Jadhav 2023 Table 2. 89.6% of the popPK analysis set carried a hyperlipidemia diagnosis and 8.2% were healthy subjects (Table 1). Coexists with DIS_DIAB: the popPK dataset included 49 participants with diabetes alone plus 310 with hyperlipidemia and diabetes (Table 1 footnote c), so the two indicators are not mutually exclusive and neither can be derived from the other.",
      source_name        = "Hyperlipidemia"
    ),
    DIS_DIAB = list(
      description        = "Type 2 diabetes mellitus indicator, 1 = participant has T2DM, 0 = no diabetes",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no diabetes)",
      notes              = "Proportional-shift covariate on CL/F (theta = -0.177), Jadhav 2023 Table 2 row 'T2DM'. 16.1% of the popPK analysis set had diabetes (Table 1). The source paper specifies type 2 diabetes; the canonical column pools T1/T2 per the register.",
      source_name        = "T2DM"
    ),
    CRCL = list(
      description        = "MDRD-estimated glomerular filtration rate, expressed in absolute mL/min WITHOUT body-surface-area adjustment",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power covariate on CL/F (exponent 0.574), Jadhav 2023 Table 2 row 'eGFR'. Table 1 footnote a: 'eGFR was calculated using the MDRD formula and expressed in absolute units (mL/min) without body surface area adjustment' -- so the values are NOT in the register's usual mL/min/1.73 m^2 units; use raw mL/min when supplying this column to the model (same non-BSA-normalized usage as Georges_2009_ceftazidime.R and Delattre_2010_amikacin.R). Reference value not printed in the paper; the model file uses the popPK analysis-set median of 89.3 mL/min (Table 1; mean 91.4, range 16.9-286), matching the median-normalisation convention applied to WT and AGE and sitting essentially on the paper's own normal-renal-function threshold (Fig. 2 defines normal as eGFR >= 90 mL/min, mild 60 to < 90, moderate 30 to < 60).",
      source_name        = "eGFR"
    ),
    CONMED_EZE = list(
      description        = "Concomitant ezetimibe indicator, 1 = receiving concomitant ezetimibe, 0 = not receiving ezetimibe",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant ezetimibe)",
      notes              = "Proportional-shift covariate on CL/F (theta = -0.0934), Jadhav 2023 Table 2. Retained by the ad hoc concomitant-medication forward-selection step (Methods). Jadhav 2023 Results states the exposure impact was not clinically meaningful.",
      source_name        = "Ezetimibe"
    ),
    CONMED_SIMVASTATIN = list(
      description        = "Concomitant simvastatin indicator, 1 = receiving concomitant simvastatin, 0 = not receiving simvastatin",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant simvastatin)",
      notes              = "Proportional-shift covariate on Vc/F (theta = -0.154), Jadhav 2023 Table 2. Retained by the ad hoc concomitant-medication forward-selection step, which screened atorvastatin, pravastatin, rosuvastatin, simvastatin, metformin, and ezetimibe against Ka, F1, CL/F, and Vc/F (Methods).",
      source_name        = "Simvastatin"
    ),
    CONMED_ATORVASTATIN = list(
      description        = "Concomitant atorvastatin indicator, 1 = receiving concomitant atorvastatin, 0 = not receiving atorvastatin",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant atorvastatin; the F1 = 1 anchor)",
      notes              = "Proportional-shift covariate on relative oral bioavailability F1 (theta = 0.142), Jadhav 2023 Table 2. Table 2 footnote b: 'F1, relative oral bioavailability was estimated for participants receiving concomitant atorvastatin relative to those not receiving concomitant atorvastatin.' F1 is fixed at 1 in the reference (no-atorvastatin) condition so that CL/F and Vc/F are apparent parameters (Fig. 1 caption).",
      source_name        = "Atorvastatin"
    ),
    FED = list(
      description        = "Fed-state dose-record indicator, 1 = dose administered with food, 0 = fasted",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (fasted)",
      notes              = "Proportional-shift covariate on the absorption rate constant Ka (theta = -0.777), Jadhav 2023 Table 2 row 'Food on Ka'. Evaluated as a structural covariate in the base model. Jadhav 2023 Discussion: 'Food was predicted to decrease the transit rate of bempedoic acid absorption by 78% without affecting the extent of absorption' -- i.e. Ka is multiplied by (1 - 0.777) = 0.223 for fed doses and F1 is unchanged. Per-dose-record indicator.",
      source_name        = "Food"
    ),
    SAMPLE_INTENSIVE = list(
      description        = "Per-observation sampling-intensity indicator, 1 = serial (intensive) PK sampling, 0 = sparse PK sampling",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (sparse PK sampling)",
      notes              = "Record-level indicator that switches the proportional residual-error magnitude: 31.9% for serial PK sampling and 54.3% for sparse PK sampling (Jadhav 2023 Table 2). Methods: 'separate log-additive residual error terms were estimated for serial and sparse sampling conditions'. Set SAMPLE_INTENSIVE = 1 on observations drawn from a rich post-dose profile (the phase 1 and several phase 2 studies in Online Resource 1) and 0 on pre-dose trough observations (the phase 3 studies).",
      source_name        = "SAMPLE_INTENSIVE"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 2232L,
    n_observations = 10347L,
    n_studies      = 22L,
    phases         = "Nine phase 1, nine phase 2, and four phase 3 studies (Online Resource 1)",
    age_range      = "18.0-89.0 years",
    age_median     = "62.0 years (mean 60.5, SD 12.3)",
    weight_range   = "42.5-152 kg",
    weight_median  = "83.7 kg (mean 85.1, SD 17.3)",
    sex_female_pct = 40.6,
    race_ethnicity = c(White = 89.0, Black = 9.2, Asian = 1.0, NativeAmerican = 0.2, NativeHawaiian = 0.2, Other = 0.6),
    ethnicity      = "Hispanic 9.2%, Not Hispanic 90.7%",
    disease_state  = "Healthy subjects (8.2%) and patients with hyperlipidemia (89.6%), type 2 diabetes mellitus (16.1%), or renal or hepatic impairment.",
    renal_function = "eGFR (MDRD, absolute mL/min) mean 91.4 (SD 24.5), median 89.3, range 16.9-286 mL/min; the pool includes a dedicated renal-impairment study (Study 023) spanning normal function through severe impairment.",
    dose_range     = "Oral bempedoic acid 2.5-250 mg, single and once-daily multiple dose; the commercial regimen is 180 mg once daily.",
    regions        = "Multi-regional pool of 22 Esperion-sponsored clinical studies.",
    notes          = "Demographics from Jadhav 2023 Table 1 (popPK column). The PK dataset comprised 11,124 samples from an intended target population of 2261 participants; post-dose samples below the limit of quantification (5.8%) and pre-first-dose samples were excluded, leaving 10,347 quantifiable samples from 2232 participants. Analysis performed in NONMEM 7.3. The bempedoic acid metabolite ESP15228 was evaluated but excluded from the final structure because it circulates at a constant ~20% of parent concentrations, so parent drug serves as the surrogate for total active exposure."
  )

  ini({
    # =========================================================================
    # Structural parameters -- Jadhav 2023 Table 2 (Final PopPK model parameters)
    #
    # Reference participant: male, non-Black, no hyperlipidemia, no T2DM, no
    # concomitant ezetimibe / simvastatin / atorvastatin, fasted dosing,
    # WT = 83.7 kg, AGE = 62 years, eGFR = 89.3 mL/min.
    #
    # Covariate model forms (Jadhav 2023 Methods, "PK model covariate analysis"):
    #   continuous : theta_TV,ij = theta_REF * (x_ij / x_REF)^theta_x
    #   categorical: theta_TV,ij = theta_REF * (1 + theta_x * x_ij)
    # =========================================================================
    lcl  <- log(0.755); label("Apparent systemic clearance CL/F at the reference condition (L/h)")               # Table 2: CL/F = 0.755 L/h (%RSE 2.6; 95% CI 0.716, 0.794)
    lvc  <- log(19.1);  label("Apparent central distribution volume Vc/F at the reference condition (L)")         # Table 2: Vc/F = 19.1 L (%RSE 6.9; 95% CI 16.5, 21.7)
    lka  <- log(1.41);  label("Absorption / transit rate constant Ka under fasted dosing (1/h)")                  # Table 2: Ka = 1.41 1/h (%RSE 5.6; 95% CI 1.25, 1.56)
    lk12 <- log(0.184); label("Central-to-peripheral distribution rate constant K23 (1/h)")                       # Table 2: K23 = 0.184 1/h (%RSE 8.3; 95% CI 0.154, 0.215)
    lk21 <- log(0.156); label("Peripheral-to-central distribution rate constant K32 (1/h)")                       # Table 2: K32 = 0.156 1/h (%RSE 4.3; 95% CI 0.143, 0.169)

    # Relative oral bioavailability. Fig. 1 caption: "F1 oral bioavailability
    # (typical value fixed at 1, such that systemic parameters are apparent)".
    lfdepot <- fixed(log(1)); label("Relative oral bioavailability F1 in the no-atorvastatin reference (fraction)")  # Fig. 1 caption: F1 typical value fixed at 1

    # ---- Covariate effects on CL/F (Table 2, "Covariates of CL/F") ----
    e_wt_cl              <- 0.61;    label("Power exponent of (WT / 83.7 kg) on CL/F (unitless)")                                   # Table 2: Body weight = 0.61 (%RSE 9.8; 95% CI 0.493, 0.727)
    e_crcl_cl            <- 0.574;   label("Power exponent of (CRCL / 89.3 mL/min) on CL/F (unitless)")                             # Table 2: eGFR = 0.574 (%RSE 6.3; 95% CI 0.503, 0.645)
    e_sexf_cl            <- -0.127;  label("Proportional shift of female sex on CL/F: CL/F * (1 + theta * SEXF) (unitless)")      # Table 2: Female sex = -0.127 (%RSE 16.9; 95% CI -0.169, -0.0846)
    e_race_black_cl      <- -0.143;  label("Proportional shift of Black race on CL/F (unitless)")                                 # Table 2: Black race = -0.143 (%RSE 18.0; 95% CI -0.194, -0.0924)
    e_dis_hyperlip_cl    <- -0.0945; label("Proportional shift of hyperlipidemia on CL/F (unitless)")                             # Table 2: Hyperlipidemia = -0.0945 (%RSE 28.4; 95% CI -0.147, -0.0416)
    e_dis_diab_cl        <- -0.177;  label("Proportional shift of type 2 diabetes mellitus on CL/F (unitless)")                   # Table 2: T2DM = -0.177 (%RSE 16.7; 95% CI -0.235, -0.119)
    e_conmed_eze_cl      <- -0.0934; label("Proportional shift of concomitant ezetimibe on CL/F (unitless)")                      # Table 2: Ezetimibe = -0.0934 (%RSE 28.4; 95% CI -0.146, -0.0413)

    # ---- Covariate effects on Vc/F (Table 2, "Covariates of Vc/F") ----
    e_wt_vc              <- 0.94;    label("Power exponent of (WT / 83.7 kg) on Vc/F (unitless)")                                   # Table 2: Body weight = 0.94 (%RSE 20.8; 95% CI 0.557, 1.32)
    e_age_vc             <- 0.743;   label("Power exponent of (AGE / 62 years) on Vc/F (unitless)")                               # Table 2: Age = 0.743 (%RSE 15.9; 95% CI 0.511, 0.976)
    e_sexf_vc            <- -0.0895; label("Proportional shift of female sex on Vc/F (unitless)")                                 # Table 2: Female sex = -0.0895 (%RSE 33.4; 95% CI -0.148, -0.0308)
    e_conmed_simvastatin_vc <- -0.154; label("Proportional shift of concomitant simvastatin on Vc/F (unitless)")                  # Table 2: Simvastatin = -0.154 (%RSE 29.2; 95% CI -0.242, -0.0654)

    # ---- Covariate effect on Ka (Table 2, structural food covariate) ----
    e_fed_ka             <- -0.777;  label("Proportional shift of fed-state dosing on Ka: Ka * (1 + theta * FED) (unitless)")     # Table 2: Food on Ka = -0.777 (%RSE 1.4; 95% CI -0.799, -0.754)

    # ---- Covariate effect on F1 (Table 2, "Covariates of F1") ----
    e_conmed_atorvastatin_fdepot <- 0.142; label("Proportional shift of concomitant atorvastatin on relative oral bioavailability F1 (unitless)")  # Table 2: Atorvastatin = 0.142 (%RSE 19.1; 95% CI 0.0886, 0.195)

    # ---- Inter-individual variability (Table 2, "IIV, %CV") ----
    # Methods: "Interindividual variability was assumed to follow a log-normal
    # distribution", so omega^2 = log(CV^2 + 1). Jadhav 2023 reports only the
    # diagonal %CV values; no OMEGA off-diagonals / correlations are published.
    etalcl ~ 0.084533   # Table 2: IIV CL/F = 29.7 %CV -> log(1 + 0.297^2) = 0.084533
    etalvc ~ 0.693147   # Table 2: IIV Vc/F = 100  %CV -> log(1 + 1.000^2) = 0.693147
    etalka ~ 0.435749   # Table 2: IIV Ka   = 73.9 %CV -> log(1 + 0.739^2) = 0.435749

    # ---- Residual error (Table 2, "Residual error, %") ----
    # Methods: "separate log-additive residual error terms were estimated for
    # serial and sparse sampling conditions". A NONMEM log-additive residual is
    # proportional in nlmixr2's linear space, so both magnitudes enter as prop().
    # Table 2 footnote a: residual error values are reported as positive numbers
    # via sqrt(estimate^2).
    propSdSerial <- 0.319; label("Proportional residual error, serial (intensive) PK sampling (fraction)")   # Table 2: Serial PK sampling = 31.9% (%RSE 1.1; 95% CI 31.2, 32.6)
    propSdSparse <- 0.543; label("Proportional residual error, sparse PK sampling (fraction)")               # Table 2: Sparse PK sampling = 54.3% (%RSE 1.3; 95% CI 52.9, 55.7)
  })

  model({
    # ---- Individual PK parameters (reference: 83.7 kg, 62 y, eGFR 89.3 mL/min,
    #      male, non-Black, no hyperlipidemia / T2DM, no concomitant LMT, fasted) ----
    cl <- exp(lcl + etalcl) *
      (WT / 83.7)^e_wt_cl *
      (CRCL / 89.3)^e_crcl_cl *
      (1 + e_sexf_cl * SEXF) *
      (1 + e_race_black_cl * RACE_BLACK) *
      (1 + e_dis_hyperlip_cl * DIS_HYPERLIP) *
      (1 + e_dis_diab_cl * DIS_DIAB) *
      (1 + e_conmed_eze_cl * CONMED_EZE)

    vc <- exp(lvc + etalvc) *
      (WT / 83.7)^e_wt_vc *
      (AGE / 62)^e_age_vc *
      (1 + e_sexf_vc * SEXF) *
      (1 + e_conmed_simvastatin_vc * CONMED_SIMVASTATIN)

    ka <- exp(lka + etalka) * (1 + e_fed_ka * FED)

    k12 <- exp(lk12)
    k21 <- exp(lk21)
    kel <- cl / vc

    fdepot <- exp(lfdepot) * (1 + e_conmed_atorvastatin_fdepot * CONMED_ATORVASTATIN)

    # ---- ODE system (Jadhav 2023 Fig. 1) ----
    # Absorption (F1) --Ka--> Transit --Ka--> Central (Vc/F) <--K23/K32--> Peripheral,
    # with CL/F eliminating from Central. The single transit compartment shares the
    # same rate constant Ka as the depot-to-transit step (Fig. 1 labels both arrows Ka).
    d/dt(depot)       <- -ka * depot
    d/dt(transit1)    <-  ka * depot - ka * transit1
    d/dt(central)     <-  ka * transit1 - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    f(depot) <- fdepot

    # ---- Observation and residual error ----
    # Dose in mg and vc in L give central / vc in mg/L = ug/mL, the paper's
    # bempedoic acid concentration unit.
    Cc <- central / vc

    # Residual-error magnitude switched per observation by SAMPLE_INTENSIVE:
    #   SAMPLE_INTENSIVE = 1 -> serial (intensive) sampling, 31.9%
    #   SAMPLE_INTENSIVE = 0 -> sparse sampling,             54.3%
    propSd <- propSdSerial * SAMPLE_INTENSIVE + propSdSparse * (1 - SAMPLE_INTENSIVE)
    Cc ~ prop(propSd)
  })
}
