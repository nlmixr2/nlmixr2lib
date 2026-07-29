Przybylowski_2015_propofol <- function() {
  description <- "Three-compartment IV population PK plus effect-compartment sigmoidal Emax PD model for propofol in adult ASA III cancer patients undergoing major lung surgery under propofol-fentanyl total intravenous anesthesia (Przybylowski 2015; N = 23). The PD response is the AAI (A-line ARX-Index) auditory-evoked-potential depth-of-anesthesia index with the maximum effect fixed to 1 and the pretreatment baseline fixed to 87 from a prior study. Inter-individual variability was estimated on Vc, CL, and the deep-compartment intercompartmental clearance Q2 for PK and on Ce50, gamma (Hill), and ke0 for PD; IIV on Vt1, Q1, Vt2 was fixed to 0 (data uninformative). No demographic, biochemical, or hemodynamic covariates were retained in the final model (Results)."
  reference <- paste(
    "Przybylowski K, Tyczka J, Szczesny D, Bienert A, Wiczling P, Kut K,",
    "Plenzler E, Kaliszan R, Grzeskowiak E.",
    "Pharmacokinetics and pharmacodynamics of propofol in cancer patients",
    "undergoing major lung surgery.",
    "J Pharmacokinet Pharmacodyn. 2015;42(3):111-122.",
    "doi:10.1007/s10928-015-9404-6.",
    sep = " "
  )
  vignette <- "Przybylowski_2015_propofol"
  units    <- list(time = "minute", dosing = "mg", concentration = "mg/L")

  covariateData <- list()  # No covariates retained in the final model (Przybylowski 2015 Results: no statistically significant (p < 0.01) relationships identified across body weight, gender, age, blood pressure, heart rate, laboratory blood tests, and stage of lung cancer).

  population <- list(
    species        = "human",
    n_subjects     = 23L,                                        # Przybylowski 2015 Methods + Table 1
    n_studies      = 1L,                                         # single-centre, Poznan + Gdansk
    age_range      = "51-75 years",                              # Przybylowski 2015 Table 1 (median 60)
    age_median     = "60 years",                                 # Przybylowski 2015 Table 1
    weight_range   = "44-125 kg",                                # Przybylowski 2015 Table 1 (median 77)
    weight_median  = "77 kg",                                    # Przybylowski 2015 Table 1
    height_range   = "152-183 cm (median 172)",                  # Przybylowski 2015 Table 1
    lbm_range      = "34.7-77.1 kg (lean body mass; median 56.4)", # Przybylowski 2015 Table 1
    sex_female_pct = 34.8,                                       # Przybylowski 2015 Table 1: 8 F / 15 M (8/23 = 34.8%)
    disease_state  = "ASA III patients scheduled for a major lung surgery due to lung cancer between December 2010 and September 2011. Comorbidities observed in the cohort: hypertension, diabetes, major depression, obesity, chronic obstructive pulmonary disease, renal failure, gastritis, hyperthyroidism, atrial fibrillation, coronary artery disease and post-myocardial infarction. Some patients had hypoalbuminemia and increased leukocytes.",
    dose_range     = "Oral premedication: 7.5 mg midazolam. Induction: fentanyl 3 ug/kg IV + propofol 2 mg/kg IV bolus. Maintenance: propofol continuous IV infusion at 8 mg/kg/h, adjusted to maintain AAI 15-25. Median propofol infusion duration 140 min (range 67-214). Thoracic epidural at T5 (6 mL bolus with 0.1 mg fentanyl + 20 mg bupivacaine, then 0.125 % bupivacaine at 4-6 mL/h). Rocuronium 0.6 mg/kg IV to facilitate intubation.",
    regions        = "Poland (single-centre; surgery in Poznan, modelling at Medical University of Gdansk).",
    n_observations = "423 propofol plasma concentrations and 462 AAI index measurements.",
    notes          = "Plasma propofol assayed by HPLC-fluorescence within 8 weeks; LLOQ 0.01 mg/L. AAI measured with the AEP/2 Monitor (Danmeter, software 1.6); AAI scaled per Vereecke et al. with a hard upper limit of 60. AAI values above 60 were handled in NONMEM using the Beal M3 method with the F-FLAG option (Methods 'Handling the AAI index measurements with upper limit')."
  )

  ini({
    # Final PK estimates (Przybylowski 2015 Table 2). The paper uses NONMEM
    # ADVAN6 with volumes and clearances; CL and Q in L/min, volumes in L,
    # plasma concentration in mg/L. The paper's table names map as
    #   VP  (paper, central compartment)  -> vc  (nlmixr2lib)
    #   VT1 (paper, peripheral 1)         -> vp
    #   VT2 (paper, peripheral 2)         -> vp2
    #   Q1  (paper)                       -> q
    #   Q2  (paper, deep)                 -> q2
    lvc  <- log(5.11);   label("Central volume of distribution Vc (L)")                  # Przybylowski 2015 Table 2: VP = 5.11 (CV% 17.9; 5-95 CI 3.61-6.61)
    lcl  <- log(2.38);   label("Clearance CL (L/min)")                                   # Przybylowski 2015 Table 2: CL = 2.38 (CV% 8.4; 5-95 CI 2.07-2.71)
    lvp  <- log(14.2);   label("Peripheral volume of distribution Vp (L)")               # Przybylowski 2015 Table 2: VT1 = 14.2 (CV% 33.2; 5-95 CI 6.44-21.9)
    lq   <- log(1.17);   label("Inter-compartmental clearance Q (L/min)")                # Przybylowski 2015 Table 2: Q1 = 1.17 (CV% 14.5; 5-95 CI 0.891-1.45)
    lvp2 <- log(189);    label("Second peripheral volume of distribution Vp2 (L)")       # Przybylowski 2015 Table 2: VT2 = 189 (CV% 44.6; 5-95 CI 50.3-327)
    lq2  <- log(0.608);  label("Deep inter-compartmental clearance Q2 (L/min)")          # Przybylowski 2015 Table 2: Q2 = 0.608 (CV% 46.3; 5-95 CI 0.145-1.07)

    # Final PD estimates (Przybylowski 2015 Table 3). Effect-compartment sigmoidal
    # Emax model linked to plasma. Baseline AAI0 fixed at 87 from Vereecke et al.
    # [reference 41 in the paper]; Emax fixed to 1 so AAI -> 0 at very high
    # effect-site concentrations. The paper writes the shape factor as gamma;
    # named lhill here to avoid clashing with the gamma() R primitive and to
    # match the Talke 2018 dexmedetomidine convention.
    laai0  <- fixed(log(87));    label("Baseline (pretreatment) AAI index (unitless)")        # Przybylowski 2015 Table 3: AAI0 = 87 (fixed; footnote a "Fixed based on study [41]")
    lemax  <- fixed(log(1));     label("Maximum AAI suppression Emax (fraction)")             # Przybylowski 2015 Table 3: EMAX = 1 (fixed)
    lce50  <- log(1.40);         label("Effect-site Ce50 (mg/L)")                             # Przybylowski 2015 Table 3: Ce50 = 1.40 (CV% 9.3; 5-95 CI 1.18-1.61)
    lhill  <- log(2.76);         label("Hill / shape factor gamma (unitless)")                # Przybylowski 2015 Table 3: gamma = 2.76 (CV% 14.3; 5-95 CI 2.11-3.41)
    lke0   <- log(0.103);        label("Effect-compartment equilibration rate ke0 (1/min)")   # Przybylowski 2015 Table 3: ke0 = 0.103 (CV% 10.7; 5-95 CI 0.085-0.121); half-life log(2)/ke0 = 6.72 min

    # IIV (log-normal; Przybylowski 2015 Eq. for inter-individual variability).
    # The paper reports IIV magnitudes as CV%, e.g., "73, 22 and 59 %" for Vc,
    # CL, and the deep Q. Internal omega^2 = log(CV^2 + 1).
    # IIV on Vt1, Q1, Vt2 was fixed to 0 in the source (Table 2 footnote a
    # "0 FIX"; shrinkage 100%); those etas are omitted entirely.
    etalvc  ~ 0.429667    # Przybylowski 2015 Table 2: IIV VP   73.3% CV; omega^2 = log(0.733^2 + 1)
    etalcl  ~ 0.046035    # Przybylowski 2015 Table 2: IIV CL   21.7% CV; omega^2 = log(0.217^2 + 1)
    etalq2  ~ 0.301430    # Przybylowski 2015 Table 2: IIV Q2   59.3% CV; omega^2 = log(0.593^2 + 1)
    etalce50 ~ 0.063454   # Przybylowski 2015 Table 3: IIV Ce50 25.6% CV; omega^2 = log(0.256^2 + 1)
    etalhill ~ 0.147685   # Przybylowski 2015 Table 3: IIV gamma 39.9% CV; omega^2 = log(0.399^2 + 1)
    etalke0  ~ 0.172587   # Przybylowski 2015 Table 3: IIV ke0  43.4% CV; omega^2 = log(0.434^2 + 1)

    # Residual error. The paper reports a proportional error model on plasma
    # (Cp,obs = Cp * (1 + eps_prop)) and a combined additive + proportional
    # model on AAI (AAI_obs = AAI * (1 + eps_prop_AAI) + eps_add_AAI). Table 2
    # and Table 3 report the proportional components as CV% (30.0 and 31.8)
    # and the additive component as the variance value 0.553 in AAI^2; the
    # standard deviation passed to add() below is sqrt(0.553) = 0.7436.
    propSd     <- 0.30;             label("Proportional residual error on plasma propofol (fraction)")  # Przybylowski 2015 Table 2: sigma^2 prop;Cp = 30.0% CV (CV% on estimate 6.3; shrinkage 5.6; 5-95 CI 26.9-33.1)
    addSd_aai  <- 0.7436;           label("Additive residual error on AAI (AAI units)")                 # Przybylowski 2015 Table 3: sigma^2 add;AAI = 0.553 AAI^2 (CV% on estimate 48.8; shrinkage 5.2; 5-95 CI 0.11-1.00); SD = sqrt(0.553) = 0.7436
    propSd_aai <- 0.318;            label("Proportional residual error on AAI (fraction)")              # Przybylowski 2015 Table 3: sigma^2 prop;AAI = 31.8% CV (CV% on estimate 6.7; shrinkage 5.2; 5-95 CI 28.3-35.3)
  })

  model({
    # Individual PK parameters
    vc  <- exp(lvc  + etalvc)
    cl  <- exp(lcl  + etalcl)
    vp  <- exp(lvp)
    q   <- exp(lq)
    vp2 <- exp(lvp2)
    q2  <- exp(lq2  + etalq2)

    # Individual PD parameters
    aai0 <- exp(laai0)
    emax <- exp(lemax)
    ce50 <- exp(lce50 + etalce50)
    hill <- exp(lhill + etalhill)
    ke0  <- exp(lke0  + etalke0)

    # Micro-rate constants for the explicit ODE form (matches Przybylowski 2015
    # NONMEM ADVAN6 differential equations for VP, CT1, CT2 with R0 input).
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # Three-compartment IV disposition; drug is administered into `central` via
    # the user data set (induction bolus + continuous maintenance infusion).
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1 - k13 * central + k31 * peripheral2
    d/dt(peripheral1) <-                  k12 * central - k21 * peripheral1
    d/dt(peripheral2) <-                                                      k13 * central - k31 * peripheral2

    # Plasma propofol concentration and effect-compartment equilibration
    Cc <- central / vc                                              # mg/L (= ug/mL)
    d/dt(effect) <- ke0 * (Cc - effect)                             # effect-compartment concentration in mg/L

    # Sigmoidal Emax PD on the AAI index. With emax fixed to 1, AAI -> 0 as
    # effect -> Inf and AAI -> aai0 (= 87) at effect = 0.
    aai <- aai0 * (1 - emax * effect^hill / (ce50^hill + effect^hill))

    # Observation models. Cc carries the bare propSd (parent observation
    # convention); aai is a paper-named secondary output with its own
    # additive + proportional residual.
    Cc  ~ prop(propSd)
    aai ~ add(addSd_aai) + prop(propSd_aai)
  })
}
