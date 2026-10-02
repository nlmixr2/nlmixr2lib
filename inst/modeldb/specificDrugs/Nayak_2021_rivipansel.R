Nayak_2021_rivipansel <- function() {
  description <- paste(
    "Two-compartment IV population PK model for rivipansel in healthy",
    "subjects, subjects with renal or hepatic impairment, and patients with",
    "sickle cell disease (SCD) with and without an active vaso-occlusive",
    "crisis (VOC), fit jointly to plasma and urine concentrations from",
    "eight studies (Nayak 2021). Total clearance is split into a renal arm",
    "and a small non-renal arm; the renal arm is a power function of",
    "BSA-derived absolute creatinine clearance (CKD-EPI in adults, bedside",
    "Schwartz in children; reference 119 mL/min) with an additional",
    "exponent below 60 mL/min, and is 11.6% higher in the phase II SCD-VOC",
    "study. Central and peripheral volumes share one body-weight exponent",
    "(reference 75 kg). The renally cleared amount accumulates in a urine",
    "compartment that is read as a urinary concentration over the",
    "collected urine volume."
  )
  reference <- paste(
    "Nayak S, Tammara B, Harnisch LO.",
    "Population Pharmacokinetic Analysis of Rivipansel in Healthy Subjects",
    "and Subjects with Sickle Cell Disease.",
    "Drugs R D. 2021;21(2):217-229.",
    "doi:10.1007/s40268-021-00346-3."
  )
  vignette <- "Nayak_2021_rivipansel"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  compartmentData <- list(
    central = list(analyte = "rivipansel", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "rivipansel", units = "mg", specimen = "plasma", verified = TRUE),
    urine = list(analyte = "rivipansel", units = "mg", specimen = "urine", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = paste(
        "Creatinine clearance in ABSOLUTE mL/min (NOT indexed to 1.73 m^2).",
        "Adults: CKD-EPI estimated GFR; children younger than 18 years:",
        "bedside Schwartz estimated GFR; both de-indexed with the Du Bois",
        "body-surface area (eGFR x BSA / 1.73)."
      ),
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Source column CRCL3 ('CRCL calculated by CKD-EPI formula', ESM",
        "control-stream header). Methods 2.3 describe it as 'normalized by",
        "the BSA'; Table 2 settles the direction: children aged 6-11 years",
        "with mean serum creatinine 0.4 mg/dL and mean height 136 cm have a",
        "bedside-Schwartz eGFR of about 140 mL/min/1.73 m^2 but a tabulated",
        "creatinine clearance of 91.6 mL/min, which is that eGFR times their",
        "Du Bois BSA (about 1.1 m^2) over 1.73. Reference 119 mL/min is the",
        "control-stream MDCRCL3 ('median creatinine clearance'). Enters the",
        "renal clearance arm as (CRCL/119)^0.477, with the exponent raised",
        "by 0.413 when CRCL < 60 mL/min (Methods 2.6)."
      ),
      source_name = "CRCL3"
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power covariate on BOTH the central and the peripheral volume with",
        "one shared exponent, reference 75 kg (control stream: MDWT = 75,",
        "SWTB = (WT/MDWT)**THETA(9) applied to V1 and V2). Table 3 labels",
        "the exponent 'Weight factor exponent on Vc'; the stream and Methods",
        "2.3 ('volume in the central and peripheral compartment was modeled",
        "as ... F_WT') apply it to both volumes."
      ),
      source_name = "WT"
    ),
    STUDY_RIV201 = list(
      description = paste(
        "Phase II rivipansel study B5201012 / NCT01119833 indicator (SCD",
        "patients hospitalized for a VOC). 1 = subject from that study;",
        "0 = subject from any of the seven other fitted studies (healthy",
        "volunteers, renal- and hepatic-impairment subjects, elderly",
        "subjects, and adults with stable SCD not in VOC)."
      ),
      units = "(binary)",
      type = "binary",
      reference_category = "0 (all fitted studies other than B5201012)",
      notes = paste(
        "Control stream: IF(STUDY.EQ.1012) SCLF = 1 + THETA(7), and the",
        "same STUDY.EQ.1012 test selects the SCD-VOC residual-error",
        "magnitudes. The paper interprets the 11.6% clearance increment as",
        "SCD hyperfiltration during an active VOC ('Fractional increase in",
        "CL for SCD-VOC subjects', Table 3); B5201012 was the only fitted",
        "study enrolling SCD-VOC patients. As coded, the phase III study",
        "(B5201002, not in the fit) does not match STUDY.EQ.1012; set",
        "STUDY_RIV201 = 1 to simulate SCD patients in VOC with the",
        "hyperfiltration factor, 0 without it."
      ),
      source_name = "STUDY"
    ),
    URINE_VOL_INTERVAL = list(
      description = "Urine volume collected in the urine collection interval containing the current urinary observation",
      units = "mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Control stream: S3 = UVOL, 'Urine is the third (output)",
        "compartment', so the urinary prediction is the amount excreted in",
        "the collection interval divided by the collected volume. Urine",
        "observations were CMT = -3 records, which in NONMEM observe the",
        "output compartment and then reset it; in rxode2 reset the urine",
        "state with an evid = 5, amt = 0 record just after each interval",
        "boundary. The source does not state the UVOL units; stored here in",
        "mL and divided by 1000 in model() so that the mg amount yields",
        "mg/L (= ug/mL). Only the urine observation uses it; it has no",
        "effect on plasma concentrations."
      ),
      source_name = "UVOL"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "The final control stream carries SAGE = (AGE/MDAGE)**THETA(11) with",
        "MDAGE = 32 years on renal clearance, but THETA(11) is '(0 FIX)'",
        "('exponent on AGE for clearance (fixed to 0)'), so the term is 1",
        "for every subject; Table 3 lists no age effect. Age enters only",
        "through the CKD-EPI / Schwartz renal-function estimate."
      ),
      source_name = "AGE"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 217,
    n_studies = 8,
    age_range = "12-80 years in the fitted studies (phase III external validation cohort: 6-56 years)",
    weight_range = "41-145 kg across the adult cohorts of Table 2 (the SCD-VOC summaries pool phase II and phase III); 18-65 kg in the phase III children aged 6-11 years",
    sex_female_pct = 30.4,
    race_ethnicity = "Phase I: 52.3% White, 35.6% Black, 8.0% Asian, 4.0% other; phase II: 95.3% Black, 4.7% other",
    disease_state = paste(
      "Healthy adults (including elderly subjects), adults with mild,",
      "moderate or severe renal impairment, adults with moderate hepatic",
      "impairment, adults with stable SCD not in VOC (B5201011), and",
      "adolescents and adults with SCD hospitalized for a VOC (phase II",
      "B5201012). The phase III RESET study (B5201002; 156 SCD-VOC",
      "patients aged 6 years and older) was held out of the fit and used",
      "for external validation."
    ),
    dose_range = paste(
      "Intravenous infusion. Phase I single doses and multiple doses in",
      "healthy volunteers, single doses in special populations; phase II",
      "20 mg/kg loading + 10 mg/kg q12h (low dose) or 40 mg/kg loading +",
      "20 mg/kg q12h (high dose)."
    ),
    regions = "United States (Pfizer / GlycoMimetics development program)",
    renal_function = "Creatinine clearance 14.8-214.8 mL/min across the fitted cohorts (Table 2)",
    notes = paste(
      "Table 2: 174 phase I subjects (B5201001, B5201004, B5201005,",
      "B5201006, B5201009, B5201010, B5201011) and 43 phase II subjects",
      "(B5201012) were fitted (217 subjects; 66 female). Urine data came",
      "from B5201005, B5201009, B5201010, B5201011 and B5201012. NONMEM",
      "7.4.1, FOCE with interaction (run58). Supersedes the Tammara 2017",
      "three-compartment model (Tammara_2017_rivipansel), which failed to",
      "converge on the enlarged dataset."
    )
  )

  ini({
    # Structural parameters -- Nayak 2021 Table 3 (NONMEM estimates).
    lcl_renal <- log(1.15)
    label("Renal clearance at CRCL = 119 mL/min, non-SCD-VOC (L/h)") # Table 3: CL (renal) = 1.15 L/h (RSE 1.91%)
    lcl_nonren <- log(0.0718)
    label("Non-renal clearance (L/h)") # Table 3: CLn = 0.0718 L/h (RSE 22.2%)
    lvc <- log(6.75)
    label("Central volume at WT = 75 kg (L)") # Table 3: Vc = 6.75 L (RSE 2.56%)
    lq <- log(2.01)
    label("Inter-compartmental clearance (L/h)") # Table 3: Q = 2.01 L/h (RSE 6.60%)
    lvp <- log(4.48)
    label("Peripheral volume at WT = 75 kg (L)") # Table 3: Vp = 4.48 L (RSE 2.80%)

    # Covariate effects -- Table 3; forms from the ESM control stream.
    e_study_riv201_cl <- 0.116
    label("Fractional increase in renal CL in the phase II SCD-VOC study (unitless)") # Table 3: fractional increase in CL for SCD-VOC subjects = 0.116 (RSE 50.4%)
    e_crcl_cl_renal <- 0.477
    label("CRCL exponent on renal CL (unitless)") # Table 3: CrCl factor exponent on CL = 0.477 (RSE 12.8%)
    e_crcl_cl_renal_lt60 <- 0.413
    label("Additional CRCL exponent on renal CL when CRCL < 60 mL/min (unitless)") # Table 3: additional additive exponent at low CrCl (< 60 mL/min) = 0.413 (RSE 17.6%)
    e_wt_vc_vp <- 0.512
    label("Body-weight exponent shared by Vc and Vp (unitless)") # Table 3: weight factor exponent on Vc = 0.512 (RSE 14.4%); stream applies it to V1 and V2

    # IIV -- Table 3 reports sqrt(omega^2) x 100 (footnote), so
    # omega^2 = (CV/100)^2. CL-Vc block: cov = 0.513 * 0.200 * 0.325.
    etalcl_renal + etalvc ~ c(0.0400, 0.0333450, 0.105625) # Table 3: IIV CL 20.0 %, IIV Vc 32.5 %, correlation CL-Vc 0.513
    etalq ~ 0.061504 # Table 3: IIV Q 24.8 %
    etalvp ~ 0.050176 # Table 3: IIV Vp 22.4 %
    etalcl_nonren ~ 0.682276 # Table 3: IIV CLn 82.6 %

    # Residual error. Stream: Y = IPRED + SQRT(IPRED**2 + THETA**2) * ERR,
    # so the proportional SD is sqrt(sigma) (Table 3 'sigma ... (CV%)') and
    # the additive SD is THETA * sqrt(sigma) (Table 3 'Residual error
    # (ng/mL)' is the THETA; ug/mL -- see vignette).
    propSd_study201 <- 0.2149
    label("Proportional residual SD, phase II SCD-VOC study (fraction)") # Table 3: sigma residual error for SCD-VOC = 21.49 CV%
    addSd_study201 <- 1.53 * 0.2149
    label("Additive residual SD, phase II SCD-VOC study (ug/mL)") # Table 3: residual error for SCD-VOC = 1.53 times sqrt(sigma) 0.2149
    propSd_nonstudy201 <- 0.0971
    label("Proportional residual SD, other studies (fraction)") # Table 3: sigma residual error for non-SCD-VOC = 9.71 CV%
    addSd_nonstudy201 <- 0.865 * 0.0971
    label("Additive residual SD, other studies (ug/mL)") # Table 3: residual error for non-SCD-VOC = 0.865 times sqrt(sigma) 0.0971
    propSd_Curine <- 0.311
    label("Proportional residual SD, urine concentrations (fraction)") # Table 3: sigma residual error for urine data = 31.1 CV%
  })

  model({
    # Renal clearance (ESM control stream $PK):
    #   CLR = THETA(1) * EXP(ETA(1)) * SCLF * (CRCL3/119)**THETA(10)
    #   and, IF (CRCL3 .LT. 60), exponent THETA(10) + THETA(13).
    # The weight and age factors on CLR (THETA(8), THETA(11)) are FIX 0.
    crcl_exponent <- e_crcl_cl_renal + e_crcl_cl_renal_lt60 * (CRCL < 60)
    cl_renal <- exp(lcl_renal + etalcl_renal) * (CRCL / 119)^crcl_exponent *
      (1 + e_study_riv201_cl * STUDY_RIV201)
    cl_nonren <- exp(lcl_nonren + etalcl_nonren)
    cl <- cl_renal + cl_nonren

    vc <- exp(lvc + etalvc) * (WT / 75)^e_wt_vc_vp
    q <- exp(lq + etalq)
    vp <- exp(lvp + etalvp) * (WT / 75)^e_wt_vc_vp

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ADVAN3 TRANS4 with urine as the output compartment (F0 = CLR/CL):
    # the renal share of elimination accumulates in `urine`. Reset the
    # state at each collection-interval boundary (evid = 5, amt = 0).
    d/dt(central) <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(urine) <- cl_renal / vc * central

    Cc <- central / vc
    urine_volume <- max(URINE_VOL_INTERVAL / 1000, 0.001)
    Curine <- urine / urine_volume

    addSd <- addSd_study201 * STUDY_RIV201 + addSd_nonstudy201 * (1 - STUDY_RIV201)
    propSd <- propSd_study201 * STUDY_RIV201 + propSd_nonstudy201 * (1 - STUDY_RIV201)

    Cc ~ add(addSd) + prop(propSd)
    Curine ~ prop(propSd_Curine)
  })
}
