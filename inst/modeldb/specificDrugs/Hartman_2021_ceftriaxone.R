Hartman_2021_ceftriaxone <- function() {
  description <- "Two-compartment population PK model for intravenous ceftriaxone in critically ill children (0-18 years) in two Dutch paediatric intensive care units. Linear first-order elimination of TOTAL ceftriaxone, with clearance scaled by body weight and by the patient's study-period median creatinine-based eGFR (bedside Schwartz) through estimated power exponents, and central volume scaled by body weight. The unbound concentration is derived from the total central concentration by saturable (single-site) protein binding, Cu = fu * Cc with fu from the De Cock quadratic, whose maximum binding capacity scales linearly with serum albumin. Hartman 2021, n = 45 subjects, 205 total and 43 time-matched unbound plasma samples."
  reference <- "Hartman SJF, Upadhyay PJ, Hagedoorn NN, Mathot RAA, Moll HA, van der Flier M, Schreuder MF, Bruggemann RJ, Knibbe CAJ, de Wildt SN. Current Ceftriaxone Dose Recommendations are Adequate for Most Critically Ill Children: Results of a Population Pharmacokinetic Modeling and Simulation Study. Clin Pharmacokinet. 2021;60(10):1361-1372. doi:10.1007/s40262-021-01035-9. PMID 34036552. Model structure and fixed constants from the Electronic Supplementary Material 'Supplementary Model Code' (NONMEM control stream, pages 16-19)."
  vignette <- "Hartman_2021_ceftriaxone"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed per subject: the ESM control stream $INPUT defines WT as 'Bodyweight of the patient at the start of ICU admission in kg'. Power effect on CL (estimated exponent 0.67) and on central volume V1 (estimated exponent 1.28), both centred on the cohort median 14 kg (Table 2 equations; stream 'WT_median = 14'). Cohort median 14 kg (IQR 8.3-32.2), range 3.8-75 kg (Table 1).",
      source_name = "WT"
    ),
    CRCL = list(
      description = "Creatinine-based estimated glomerular filtration rate (bedside Schwartz), taken as the patient's MEDIAN value over the study period",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed per subject: the covariate on CL is the stream column EGFR_KREAT_patmedian ('Median EGFR_KREAT value per patient'), not the time-varying EGFR_KREAT. EGFR_KREAT is described in $INPUT as 'Estimated glomerular filtration rate in ml/min/1.73m2 based on height and serum creatinine using Schwartz' 2012 formula (42.3*(height/creatinine in mg per dl)^0.79)'. Power effect on CL with estimated exponent 0.575, centred on 85.22 (Table 2; stream 'EGFR_KREAT_popmedian = 85.21937', the median of patient medians). Baseline eGFR in the cohort: median 77.9 (IQR 53.4-102.6), range 13.8-269.4 (Table 1). The authors caution that extrapolation above 120 mL/min/1.73 m^2 in patients over 25 kg is poorly supported by the data (Results 3.1).",
      source_name = "EGFR_KREAT_patmedian"
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Scales the maximum protein-binding capacity linearly, BMAX = BMAXpop * (ALB / 27)^1, with the exponent FIXED to 1 (Results 3.1: an estimated exponent 'performed not significantly different from a fixed exponent of 1'; stream '(1) FIX ; COV_ALB_FU'). Centred on the cohort median 27 g/L (stream 'ALBUMIN_popmedian = 27'). Enters only the unbound-concentration prediction, never the total-drug disposition. The source dataset carried albumin as a laboratory value and coded a missing value as 0, which the stream maps to the reference ('IF (ALBUMIN.EQ.0) COV_ALB_FU = 1'); supply the median 27 g/L for a subject with no albumin value. Baseline albumin median 27 (IQR 20-32.3), range 14-45 g/L (Table 1; 36 of 45 patients).",
      source_name = "ALBUMIN"
    )
  )

  compartmentData <- list(
    central = list(
      analyte = "ceftriaxone",
      units = "mg",
      specimen = "plasma",
      verified = TRUE,
      notes = "Total (bound + unbound) ceftriaxone amount in the central compartment; Cc = central / vc is the measured total plasma concentration (stream $DES 'C1 = A(1)/V')."
    ),
    peripheral1 = list(
      analyte = "ceftriaxone",
      units = "mg",
      specimen = "plasma",
      verified = TRUE,
      notes = "Total ceftriaxone amount in the peripheral compartment (stream COMP(PERI))."
    ),
    auc_free = list(
      analyte = "ceftriaxone",
      units = "mg*h/L",
      specimen = "not applicable",
      verified = TRUE,
      notes = "Running integral of the unbound plasma concentration Cu, reproducing the stream's COMP(FREE) with 'DADT(3) = A(1)*FU/S1'. Output bookkeeping only; it does not feed back on any other state."
    )
  )

  paper_specific_compartments <- c("auc_free")

  population <- list(
    species = "human",
    n_subjects = 45,
    n_studies = 2,
    age_range = "0.08-16.67 years",
    age_median = "2.53 years (IQR 0.55-11.01)",
    weight_range = "3.8-75 kg",
    weight_median = "14 kg (IQR 8.3-32.2)",
    sex_female_pct = 46.7,
    disease_state = "Critically ill children in a level-3 paediatric intensive care unit receiving intravenous ceftriaxone (main admission reasons: infection 28.9%, respiratory failure 26.7%, surgery 22.2%, neurological impairment 13.3%). Acute kidney injury in 35.6% and augmented kidney clearance in 13.3% of patients during the study period.",
    dose_range = "100 mg/kg once daily as a 30-min IV infusion (maximum 2000 mg prophylactic, 4000 mg therapeutic); administered dose median 100 (IQR 66.7-100), range 13.3-105.3 mg/kg/day",
    regions = "Netherlands (Radboudumc Nijmegen, POPSICLE study NCT03248349; Erasmus MC-Sophia Rotterdam, PERFORM study NCT03502993)",
    renal_function = "Baseline creatinine-based eGFR median 77.9 (IQR 53.4-102.6), range 13.8-269.4 mL/min/1.73 m^2. Study-period median eGFR < 30: 11.1%; 30-80: 33.3%; 80-120: 42.2%; > 120: 13.3%. Patients on kidney replacement therapy or ECMO were excluded.",
    notes = "Table 1 of Hartman 2021. 26 richly sampled patients (POPSICLE, up to 4 samples/day for 3 days then daily) and 19 sparsely sampled patients (PERFORM, one sample on days 1-3). 205 total ceftriaxone concentrations (1.89-500 mg/L) and 43 unbound concentrations (one per patient, chosen closest to the end of the dosing interval); median unbound fraction 13.6% (range 7.6-70.3%). Baseline albumin median 27 g/L (range 14-45)."
  )

  ini({
    # Structural parameters (Table 2 'Final model estimate'; identical to
    # the $THETA initial values in the ESM control stream). Typical values
    # are for a 14-kg child with eGFR 85.22 mL/min/1.73 m^2.
    lcl <- log(0.708); label("Clearance CLpop for a 14-kg child with eGFR 85.22 (L/h)") # Table 2: CLpop 0.708 L/h (RSE 6.4%); stream THETA(1) TVCL
    lvc <- log(2.8); label("Central volume V1pop for a 14-kg child (L)") # Table 2: V1pop 2.8 L (RSE 25.2%); stream THETA(2) TVV1
    lq <- log(1.77); label("Intercompartmental clearance Qpop (L/h)") # Table 2: Qpop 1.77 L/h (RSE 36%); stream THETA(3) TVQ
    lvp <- log(2.9); label("Peripheral volume V2pop (L)") # Table 2: V2pop 2.9 L (RSE 17.6%); stream THETA(4) TVV2

    # Covariate exponents (Table 2)
    e_wt_cl <- 0.67; label("Power exponent of body weight (WT/14) on CL (unitless)") # Table 2: Theta1 0.67 (RSE 11.9%); stream THETA(5) COV_WT_CL exponent
    e_crcl_cl <- 0.575; label("Power exponent of median eGFR (eGFR/85.22) on CL (unitless)") # Table 2: Theta2 0.575 (RSE 19.1%); stream THETA(7) COV_EGFR_KREAT_EXP_CL
    e_wt_vc <- 1.28; label("Power exponent of body weight (WT/14) on V1 (unitless)") # Table 2: Theta3 1.28 (RSE 15.3%); stream THETA(6) COV_WT_V1 exponent. Results text rounds it to 1.29.

    # Saturable plasma protein binding (Table 2; ESM Equation 1)
    lbmax_pb <- log(223); label("Maximum protein-binding capacity BMAXpop at albumin 27 g/L (mg/L)") # Table 2: BMAXpop 223 mg/L (RSE 22.3%). The ESM stream $THETA(9) prints 228; the Table 2 final estimate is used (see vignette).
    lkd_pb <- log(30.3); label("Dissociation constant Kd of ceftriaxone protein binding (mg/L)") # Table 2: Kd 30.3 mg/L (RSE 26.6%); stream THETA(10) KD
    e_alb_bmax_pb <- fixed(1); label("Power exponent of albumin (ALB/27) on BMAX (unitless)") # Table 2 equation 'BMAX = BMAXpop x (ALBi/27)^1'; stream '(1) FIX ; COV_ALB_FU'

    # IIV: log-normal on CL and V1 with a full 2x2 block (variances).
    # CV check: sqrt(exp(0.135) - 1) = 38%, sqrt(exp(0.427) - 1) = 73%,
    # matching the ESM Supplementary Results '38% and 73%'.
    etalcl + etalvc ~ c(0.135, 0.178, 0.427) # Table 2: 'Cl IIV' 0.135, 'Block matrix' 0.178, 'V1 IIV' 0.427; stream $OMEGA BLOCK(2)

    # Residual error: one proportional epsilon applied to both total and
    # unbound observations in the source (stream 'Y1=IPRED*(1+ERR(1))',
    # 'Y2=IPRED*(1+ERR(1))'). $SIGMA 0.0596 is a variance, so the SD is
    # sqrt(0.0596) = 0.2441. rxode2 cannot share one parameter between two
    # endpoints, so the same value is carried twice.
    propSd <- 0.2441; label("Proportional residual SD on total ceftriaxone Cc (fraction)") # Table 2: proportional error 0.0596 (variance, RSE 14.8%); stream $SIGMA 0.0596 ERR PROP -> SD sqrt(0.0596)
    propSd_Cu <- 0.2441; label("Proportional residual SD on unbound ceftriaxone Cu (fraction)") # Same Table 2 / $SIGMA 0.0596 epsilon as propSd (shared ERR(1) in the source)
  })

  model({
    # 1. Individual PK parameters (stream $PK)
    cl <- exp(lcl + etalcl) * (WT / 14)^e_wt_cl * (CRCL / 85.22)^e_crcl_cl
    vc <- exp(lvc + etalvc) * (WT / 14)^e_wt_vc
    q <- exp(lq)
    vp <- exp(lvp)

    # 2. Saturable protein binding constants (stream $PK)
    bmax_pb <- exp(lbmax_pb) * (ALB / 27)^e_alb_bmax_pb
    kd_pb <- exp(lkd_pb)

    # 3. Micro-constants
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 4. Total central concentration and fraction unbound (ESM Equation 1;
    #    stream $DES 'FU = ((C1-BMAX-KD)+SQRT((C1-BMAX-KD)**2+4*KD*C1))/(2*C1+1E-31)')
    Cc <- central / vc
    fu <- ((Cc - bmax_pb - kd_pb) + sqrt((Cc - bmax_pb - kd_pb)^2 + 4 * kd_pb * Cc)) / (2 * Cc + 1e-31)
    Cu <- fu * Cc

    # 5. Linear two-compartment disposition of TOTAL ceftriaxone; binding
    #    does not enter the mass balance (stream $DES DADT(1), DADT(2)).
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(auc_free) <- Cu

    # 6. Observations: total (FLAG_DV_UNBOUND = 0) and unbound
    #    (FLAG_DV_UNBOUND = 1) plasma ceftriaxone
    Cc ~ prop(propSd)
    Cu ~ prop(propSd_Cu)
  })
}
