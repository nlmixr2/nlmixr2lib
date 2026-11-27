Masters_2026_sasanlimab <- function() {
  description <- paste(
    "Two-compartment population PK model for sasanlimab (PF-06801591), a",
    "subcutaneously administered anti-PD-1 IgG4 monoclonal antibody, in adults",
    "with non-muscle invasive bladder cancer (NMIBC) and other advanced solid",
    "tumors (Masters 2026; pooled studies B8011001, B8011006 [CREST] and",
    "B8011007). Linear clearance; sequential zero-order (duration D1) then",
    "first-order absorption from a subcutaneous depot with bioavailability F;",
    "estimated allometric exponents of baseline body weight on CL and V2; baseline",
    "serum creatinine, albumin, alkaline phosphatase and tumor type (NSCLC,",
    "urothelial carcinoma, other; NMIBC reference) on CL. Correlated log-normal",
    "IIV on CL and V2 and an additive residual error on log-transformed",
    "concentrations.",
    sep = " "
  )
  reference <- paste(
    "Masters JC, Jermain B, Liu A, Cesari R, Calella AM, Michelon E, Wei C,",
    "Liu S, Piscitelli J. Population Pharmacokinetics and Exposure-Response of",
    "Sasanlimab, a Checkpoint Inhibitor, in Patients with Non-Muscle Invasive",
    "Bladder Cancer and Other Solid Tumors. J Clin Pharmacol. 2026;66(9):e70295.",
    "doi:10.1002/jcph.70295.",
    "Final parameter estimates are in Table 3; the covariate equations for CL",
    "and V2 are Supporting Information Equations S1 and S2.",
    sep = " "
  )
  vignette <- "Masters_2026_sasanlimab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed (baseline). Allometric power on CL (estimated exponent 0.980) and V2 (estimated exponent 1.21) normalised to 70 kg (Equations S1, S2). Pooled-cohort median 74 kg, range 38.2-173 kg (Table S1).",
      source_name = "BWT"
    ),
    CREAT = list(
      description = "Baseline serum creatinine",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed (baseline). Power effect on CL normalised to 0.90 mg/dL (Equation S1): (CREAT / 0.90)^-0.129. Pooled-cohort median 0.9 mg/dL (Table S1).",
      source_name = "BSCR"
    ),
    ALB = list(
      description = "Baseline serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed (baseline). Masters 2026 reports albumin in g/dL and normalises to 4.3 g/dL (Equation S1): (BALB / 4.3)^-0.775. The model takes the canonical SI column in g/L and converts inline (alb_gdl = ALB / 10), so the coefficient keeps its published calibration. Pooled-cohort median 4.3 g/dL = 43 g/L (Table S1).",
      source_name = "BALB"
    ),
    ALP = list(
      description = "Baseline serum alkaline phosphatase",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed (baseline). Linear effect on CL centred at 80 U/L (Equation S1): 1 + 0.0005 * (ALP - 80). Pooled-cohort median 79.8 U/L, range 31-1132 U/L (Table S1).",
      source_name = "BALK"
    ),
    TUMTP_NSCLC = list(
      description = "Non-small cell lung cancer tumor-type indicator (1 = NSCLC, 0 = otherwise)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (NMIBC is the reference tumor type when all tumor-type indicators are 0)",
      notes = "Additive term in the (1 + Tumor Type) factor on CL of Equation S1: +33.2% CL relative to NMIBC. 198 of 1007 participants (19.7%, Table S2).",
      source_name = "TTTYPG = NSCLC"
    ),
    TUMTP_BLADDER = list(
      description = "Urothelial carcinoma (advanced / metastatic urothelial cancer cohort of study B8011001) tumor-type indicator (1 = UC, 0 = otherwise)",
      units = "(binary)",
      type = "binary",
      reference_category = "0. Participants with NMIBC carry 0: NMIBC is a separate tumor-type level in Masters 2026 and is the reference",
      notes = "Additive term in the (1 + Tumor Type) factor on CL of Equation S1: +31.0% CL relative to NMIBC. 39 of 1007 participants (3.9%, Table S2), all but one from study B8011001. NMIBC (non-muscle invasive bladder cancer, study B8011006) is also a urothelial cancer of the bladder but is coded 0 here because the source analysis treats it as its own tumor-type level.",
      source_name = "TTTYPG = UC"
    ),
    TUMTP_OTHER = list(
      description = "Other solid tumor types indicator (1 = a tumor type other than NMIBC, NSCLC or UC, 0 = otherwise)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (NMIBC reference)",
      notes = "Additive term in the (1 + Tumor Type) factor on CL of Equation S1: +10.8% CL relative to NMIBC. 64 of 1007 participants (6.4%, Table S2), from the dose-escalation and solid-tumor cohorts of studies B8011001 and B8011007.",
      source_name = "TTTYPG = Others"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "sasanlimab", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "sasanlimab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "sasanlimab", units = "mg", specimen = "serum", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 1007L,
    n_studies = 3L,
    n_observations = 5360L,
    age_range = "26-90 years",
    age_median = "66 years",
    weight_range = "38.2-173 kg",
    weight_median = "74 kg",
    sex_female_pct = 21.8,
    race_ethnicity = c(`Non-Asian` = 65.0, Asian = 33.4, Missing = 1.6),
    disease_state = "BCG-naive high-risk non-muscle invasive bladder cancer (70.1%, study B8011006 CREST Cohort A, sasanlimab combined with BCG in 698 of 706), non-small cell lung cancer (19.7%), urothelial carcinoma (3.9%) and other advanced solid tumors (6.4%).",
    dose_range = "IV infusion over 60 min of 0.5, 1, 3 or 10 mg/kg every 3 weeks (B8011001 dose escalation); SC 300 mg every 4 weeks (all three studies) and SC 600 mg every 6 weeks (B8011007).",
    regions = "United States 13.7%, Canada and Western Europe 17.5%, rest of world 68.8% (Table S2)",
    renal_function = "Creatinine clearance normal 39.4%, mild 41.7%, moderate 17.9%, severe impairment 0.3% (Table S2)",
    notes = "Pooled popPK dataset of 5360 serum concentrations from 1007 participants (B8011001 n = 146, B8011006 n = 706, B8011007 n = 155; Table 1); 12 outlying concentrations (|CWRES| > 6) removed. Baseline demographics in Supporting Information Tables S1 and S2. ADA-positive 5.6%."
  )

  ini({
    # Structural parameters: Table 3, typical values for a 70-kg NMIBC
    # participant at BSCR 0.90 mg/dL, BALB 4.3 g/dL and BALK 80 U/L.
    lcl <- log(0.142); label("Clearance CL (L/day)") # Table 3: thetaCL = 0.142 L/day (RSE 6.8%)
    lvc <- log(3.67); label("Central volume V2 (L)") # Table 3: thetaV2 = 3.67 L (RSE 11%)
    lq <- log(0.576); label("Intercompartmental clearance Q (L/day)") # Table 3: thetaQ = 0.576 L/day (RSE 20%)
    lvp <- log(2.19); label("Peripheral volume V3 (L)") # Table 3: thetaV3 = 2.19 L (RSE 13%)
    lka <- log(0.217); label("First-order SC absorption rate constant ka (1/day)") # Table 3: thetaka = 0.217 per day (RSE 11%)
    lfdepot <- log(0.697); label("Subcutaneous bioavailability F (fraction)") # Table 3: thetaF = 0.697 (RSE 6.7%)
    ld1 <- log(0.116); label("Duration of the zero-order input into the SC depot D1 (day)") # Table 3: thetaD1 = 0.116 days (RSE 25%)

    # Allometric exponents of baseline body weight (reference 70 kg)
    e_wt_cl <- 0.980; label("Allometric exponent of baseline body weight on CL (unitless)") # Table 3: thetaBBW on CL = 0.980 (RSE 5.5%); Equation S1
    e_wt_vc <- 1.21; label("Allometric exponent of baseline body weight on V2 (unitless)") # Table 3: thetaBBW on V2 = 1.21 (RSE 11%); Equation S2

    # Covariate effects on CL (Equation S1)
    e_creat_cl <- -0.129; label("Power exponent of baseline serum creatinine on CL (unitless)") # Table 3: thetaBSCR effect on CL = -0.129 (RSE 38%)
    e_alb_cl <- -0.775; label("Power exponent of baseline albumin on CL (unitless)") # Table 3: thetaBALB effect on CL = -0.775 (RSE 18%)
    e_alp_cl <- 0.0005; label("Linear effect of baseline alkaline phosphatase on CL (per U/L)") # Table 3: thetaBALK effect on CL = 0.0005 (RSE 43%)
    e_tumtp_nsclc_cl <- 0.332; label("Fractional change in CL for NSCLC vs NMIBC (unitless)") # Table 3: thetaChange in CL for NSCLC = 0.332 (RSE 13%)
    e_tumtp_bladder_cl <- 0.310; label("Fractional change in CL for urothelial carcinoma vs NMIBC (unitless)") # Table 3: thetaChange in CL for UC = 0.310 (RSE 29%)
    e_tumtp_other_cl <- 0.108; label("Fractional change in CL for other tumor types vs NMIBC (unitless)") # Table 3: thetaChange in CL for other cancer types = 0.108 (RSE 54%)

    # IIV: Table 3 omega^2 values (CV% column = sqrt(omega^2)); the
    # 'wCL/V2 Correlation' row 0.0625 is the CL-V2 covariance (its CV column
    # 25% = sqrt(0.0625)); implied correlation 0.358.
    etalcl + etalvc ~ c(0.0972, 0.0625, 0.313) # Table 3: wCL2 = 0.0972 (CV 31%), wCL/V2 = 0.0625, wV22 = 0.313 (CV 56%)

    # Residual error: additive on log-transformed concentrations (Methods)
    expSd <- 0.229; label("Additive residual error SD on log-transformed concentration (log ug/mL)") # Table 3: thetaRes. Error = 0.229 (RSE 2.9%)
  })

  model({
    # Albumin column is SI (g/L); Equation S1 is calibrated in g/dL.
    alb_gdl <- ALB / 10

    # Tumor type enters as (1 + Tumor Type) with NMIBC as reference (Equation S1)
    tumtp_cl <- 1 + e_tumtp_nsclc_cl * TUMTP_NSCLC + e_tumtp_bladder_cl * TUMTP_BLADDER + e_tumtp_other_cl * TUMTP_OTHER

    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl * (CREAT / 0.90)^e_creat_cl *
      (alb_gdl / 4.3)^e_alb_cl * (1 + e_alp_cl * (ALP - 80)) * tumtp_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    q <- exp(lq)
    vp <- exp(lvp)
    ka <- exp(lka)
    d1 <- exp(ld1)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # SC doses enter the depot by a zero-order process of duration D1
    # (dose records rate = -2), then are absorbed first-order into central.
    # IV doses go directly to central as a 60-min infusion given in the data.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    f(depot) <- exp(lfdepot)
    dur(depot) <- d1

    # mg / L = ug/mL
    Cc <- central / vc
    Cc ~ lnorm(expSd)
  })
}
