Zhao_2014_ciprofloxacin <- function() {
  description <- "Two-compartment population PK model with first-order elimination for intravenous ciprofloxacin in neonates and young infants less than three months of age (Zhao 2014). Central and peripheral volumes (V1, V2) scale allometrically with current body weight (fixed exponent 1, reference 1.955 kg); clearance (CL) and inter-compartmental clearance (Q) scale with current body weight at a fixed exponent of 0.75. CL is further multiplied by a renal-maturation factor in gestational age and postnatal age (F_age), a renal-function factor in serum creatinine (RF = exp((CREAT - 42 umol/L) * theta7)), and a fractional reduction (factor 0.708) when inotropic / vasoactive agents are coadministered. IIV is reported on V1, V2, and CL as %CV on an exponential model. Residual error is proportional. Inter-occasion variability on CL (16.4%CV) reported by Zhao 2014 is not encoded structurally here -- the source paper does not define an operational occasion mapping for the model-library use case; users who need IOV can add an OCC indicator and per-occasion eta downstream."
  reference   <- "Zhao W, Hill H, Le Guellec C, Neal T, Mahoney S, Paulus S, Castellan C, Kassai B, van den Anker JN, Kearns GL, Turner MA, Jacqz-Aigrain E. Population pharmacokinetics of ciprofloxacin in neonates and young infants less than three months of age. Antimicrob Agents Chemother. 2014;58(11):6572-6580. doi:10.1128/AAC.03568-14"
  vignette    <- "Zhao_2014_ciprofloxacin"
  units       <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Current body weight (time-varying).",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying current weight on the day of pharmacokinetic sampling. Reference value is the Zhao 2014 cohort median 1955 g = 1.955 kg (Table 4 footnote). Allometric scaling: exponent 0.75 on CL and Q (per Zhao 2014 Methods, 'Covariate analysis' paragraph 1) and exponent 1 on V1 and V2 (per the same passage; Table 4 lists '(CW/1955)^theta1' for V1 and '(CW/1955)^theta2' for V2 but the surrounding prose makes clear the exponent is the fixed allometric value of 1, not a parameter equal to theta1 or theta2).",
      source_name        = "CW"
    ),
    GA = list(
      description        = "Gestational age at birth (time-fixed per subject).",
      units              = "weeks",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed. Cohort median 27.9 weeks (range 23.3-42.0; Zhao 2014 Table 2). Reference value 27.9 weeks (Table 4 footnote). Enters CL via power form F_age = (GA / 27.9)^theta5 (Table 4) as the antenatal-maturation component.",
      source_name        = "GA"
    ),
    PNA = list(
      description        = "Postnatal age (chronological since birth, time-varying).",
      units              = "months",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Canonical PNA is in MONTHS. Zhao 2014 reports PNA in DAYS (cohort median 27 days, range 5-121; Table 2 and Table 4 footnote). The paper's expression F_age = (PNA_days / 27)^theta6 is reparameterised inside model() as (PNA_months / 0.8870)^theta6 using PNA_months = PNA_days / 30.4375 and reference 0.8870 months = 27 days / 30.4375; both numerator and denominator carry the same units factor so the ratio is unchanged. Users should supply PNA in months in the input dataset.",
      source_name        = "PNA"
    ),
    CREAT = list(
      description        = "Serum creatinine concentration (time-varying).",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying serum creatinine on the day of pharmacokinetic sampling (within 48 h per Zhao 2014 Methods). Cohort median 41 umol/L (range 22-164; Table 2). Reference value 42 umol/L (Table 4 footnote; given to the nearest integer). Enters CL via exponential centered-deviation form RF = exp((CREAT - 42) * theta7), with theta7 = -0.00335 per umol/L so increases in serum creatinine reduce CL.",
      source_name        = "CREA"
    ),
    CONMED_INOTROPE = list(
      description        = "Concomitant inotrope / vasoactive coadministration indicator (per-occasion / time-varying).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (not on any inotrope)",
      notes              = "1 = subject is receiving at least one inotropic / vasoactive agent at the time of the observation; 0 = not on any such agent. In Zhao 2014 the comedication summary (Table 2) reports 22/60 subjects with inotropic-agent coadministration. The paper applies F_inotrope = theta8 = 0.708 directly to CL when the indicator is 1 (CL reduced by 29.2% in patients on inotropes; see Zhao 2014 Table 4 footnote and Discussion paragraph 2, attributing the effect to inotrope-driven hemodynamic instability with reduced GFR, paralleling the Seay et al. 1998 dopamine / vancomycin observation in neonates).",
      source_name        = "inotropic agents"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 60L,
    n_studies      = 1L,
    n_observations = "430 ciprofloxacin plasma concentrations (265 pharmacokinetic + 165 scavenged); pharmacokinetic samples 450-15,976 ng/mL, scavenged samples 52-10,961 ng/mL (Zhao 2014 Results 'Model building' paragraph 1).",
    age_range      = "Postmenstrual age (PMA) 24.9-47.9 weeks; gestational age 23.3-42.0 weeks; postnatal age 5-121 days",
    age_median     = "PMA 36.5 weeks; GA 27.9 weeks; PNA 27 days (Table 2)",
    weight_range   = "Current weight 700-4200 g; birth weight 540-3850 g",
    weight_median  = "Current weight 1955 g; birth weight 1115 g (Table 2)",
    sex_female_pct = 35,
    race_ethnicity = c(White = 88.3, Asian = 8.3, Unknown = 3.3),
    disease_state  = "Neonates and young infants less than three months of age with suspected or documented Gram-negative serious infections receiving ciprofloxacin in two UK neonatal / paediatric intensive care units. Comedication: inotropic agents (22/60), teicoplanin (41/60), diuretics (30/60), caffeine (15/60), amoxicillin-clavulanic acid (12/60), nystatin (12/60), colistin-tobramycin-amphotericin B (10/60). Intrauterine growth retardation flagged in 3/60.",
    dose_range     = "Intravenous ciprofloxacin 5 mg/kg/dose BID (7/60), 10 mg/kg/dose BID (47/60), or 10 mg/kg/dose TID (6/60); 30 or 60 min infusion via syringe pump with microbore tubing (Zhao 2014 Methods 'Dosing regimen and pharmacokinetic sampling'). Median 18.7 mg/dose (range 4.5-40 mg/dose); median 9.1 mg/kg/dose (range 4.4-11 mg/kg/dose).",
    regions        = "United Kingdom (Liverpool Women's Hospital neonatal ICU; Alder Hey Children's Hospital paediatric ICU). Multicentre TINN consortium (EudraCT 2010-019955-23).",
    notes          = "Demographics from Zhao 2014 Table 2. Sixty newborns retained for the population pharmacokinetic analysis after exclusions (2 withdrawn, 1 on dialysis, 1 received ciprofloxacin within 36 h of inclusion). Race recorded as Caucasian (53/60), Asian (5/60), and Unknown (2/60). Comedication categories reported as one-or-the-other counts (a subject may appear in multiple categories). The two PMA strata used in the dose-optimisation simulation (PMA < 34 weeks and PMA >= 34 weeks) are not modelled covariates but post-hoc bins for the dosing-regimen recommendations."
  )

  ini({
    # Final-model estimates from Zhao 2014 Table 4 ("Final estimate" column).
    # The bootstrap medians and 5th-95th percentiles agree with the final
    # estimates and indicate model stability (Table 4 right columns and
    # Results "Model evaluation" paragraph 2).
    # The source paper reports IIV as %CV on an exponential model
    # eta_i ~ N(0, omega^2); the conversion to log-normal variance for
    # nlmixr2 is omega^2 = log(1 + CV^2).

    # Structural parameters -- reference values for the cohort medians
    # CW = 1955 g (= 1.955 kg), GA = 27.9 weeks, PNA = 27 days (= 0.8870
    # months), CREAT = 42 umol/L (Table 4 footnote).
    lvc <- log(1.97);  label("Central volume of distribution at WT = 1.955 kg (V1, L)")             # Zhao 2014 Table 4 theta1 = 1.97 L (RSE 17.7%)
    lvp <- log(1.93);  label("Peripheral volume of distribution at WT = 1.955 kg (V2, L)")          # Zhao 2014 Table 4 theta2 = 1.93 L (RSE 21.9%)
    lq  <- log(2.5);   label("Inter-compartmental clearance at WT = 1.955 kg (Q, L/h)")             # Zhao 2014 Table 4 theta3 = 2.5 L/h (RSE 32.6%)
    lcl <- log(0.366); label("Clearance at cohort-median reference (CL, L/h)")                      # Zhao 2014 Table 4 theta4 = 0.366 L/h (RSE 6.0%)

    # Allometric exponents -- fixed at the standard adult-to-paediatric
    # values per Zhao 2014 Methods "Covariate analysis" paragraph 1
    # ("allometric coefficients of 0.75 for CL and Q and of 1 for V1 and
    # V2"). Table 4 row "V1 = theta1 * (CW/1955)^theta1" appears to be a
    # rendering / OCR artifact in the published table; the surrounding
    # text and the Q row "(CW/1955)^0.75" confirm that V1 and V2 use the
    # fixed exponent 1, not theta1 / theta2.
    allo_cl <- fixed(0.75); label("Allometric exponent on CL (unitless)")   # Zhao 2014 Methods, Covariate analysis
    allo_q  <- fixed(0.75); label("Allometric exponent on Q (unitless)")    # Zhao 2014 Methods, Covariate analysis
    allo_vc <- fixed(1);    label("Allometric exponent on V1 (unitless)")   # Zhao 2014 Methods, Covariate analysis
    allo_vp <- fixed(1);    label("Allometric exponent on V2 (unitless)")   # Zhao 2014 Methods, Covariate analysis

    # Covariate effects on CL.
    e_ga_cl       <-  2.11;     label("Gestational-age power exponent on F_age (unitless; reference 27.9 weeks)")            # Zhao 2014 Table 4 theta5 = 2.11 (RSE 11.9%)
    e_pna_cl      <-  0.494;    label("Postnatal-age power exponent on F_age (unitless; reference 0.8870 months = 27 days)") # Zhao 2014 Table 4 theta6 = 0.494 (RSE 10.8%)
    e_creat_cl    <- -0.00335;  label("Serum-creatinine exponential coefficient on CL (per umol/L from CREAT = 42)")         # Zhao 2014 Table 4 theta7 = -0.00335 (RSE 46.0%)
    e_inotrope_cl <-  0.708;    label("Multiplicative factor on CL when CONMED_INOTROPE = 1 (unitless)")                     # Zhao 2014 Table 4 theta8 = 0.708 (RSE 10.9%)

    # IIV (exponential model). Source paper reports %CV in Table 4;
    # log-normal variance via omega^2 = log(1 + CV^2). No correlation
    # block is reported in the paper (independent etas).
    # V1: log(1 + 0.481^2) = 0.20814
    etalvc ~ 0.20814
    # V2: log(1 + 0.493^2) = 0.21752
    etalvp ~ 0.21752
    # CL: log(1 + 0.332^2) = 0.10448
    etalcl ~ 0.10448

    # Residual error -- proportional only (Zhao 2014 Results "Model
    # building" paragraph 2: "A proportional model best described
    # residual variability"). Table 4 reports 19.3% CV (RSE 28.2%).
    propSd <- 0.193; label("Proportional residual error (fraction)") # Zhao 2014 Table 4 = 19.3 %CV
  })

  model({
    # ----- Derived covariate terms -----
    # F_age = (GA / 27.9)^theta5 * (PNA_months / 0.8870)^theta6
    # Reference 0.8870 months = 27 days / 30.4375 (canonical PNA is in
    # months; see covariateData[[PNA]]$notes).
    pna_ref_months <- 27 / 30.4375
    f_age <- (GA / 27.9) ^ e_ga_cl * (PNA / pna_ref_months) ^ e_pna_cl

    # RF = exp((CREAT - 42) * theta7); centered-deviation exponential.
    f_renal <- exp(e_creat_cl * (CREAT - 42))

    # F_inotrope = theta8^CONMED_INOTROPE; multiplicative selection
    # (= 1 when indicator = 0, = theta8 = 0.708 when indicator = 1).
    f_inotrope <- e_inotrope_cl ^ CONMED_INOTROPE

    # ----- Individual parameters -----
    cl <- exp(lcl + etalcl) * (WT / 1.955) ^ allo_cl * f_age * f_renal * f_inotrope
    vc <- exp(lvc + etalvc) * (WT / 1.955) ^ allo_vc
    vp <- exp(lvp + etalvp) * (WT / 1.955) ^ allo_vp
    q  <- exp(lq)           * (WT / 1.955) ^ allo_q

    # ----- Micro-constants -----
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ----- ODE system (IV infusion to central) -----
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # ----- Observation and error -----
    # Dose in mg, volume in L -> concentration in mg/L.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
