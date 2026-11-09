Pokorna_2024_meropenem <- function() {
  description <- "One-compartment IV population PK model for meropenem in 45 critically ill neonates and children aged 0-5 years, 38 of them on extracorporeal membrane oxygenation (ECMO) and 31 on continuous renal replacement therapy (CRRT) (Pokorna 2024). Clearance and volume of distribution both scale linearly with body weight normalised to the cohort median of 7.88 kg (CL 1.09 L/h, V 3.98 L). CRRT running at the time of the sample raises the volume of distribution by a factor of (1 + 1.04), about two-fold; ECMO had no retained effect. Log-normal interindividual variability on CL and V, and a proportional residual error."
  reference <- "Pokorna P, Michalickova D, Tibboel D, Berner J. Meropenem Disposition in Neonatal and Pediatric Extracorporeal Membrane Oxygenation and Continuous Renal Replacement Therapy. Antibiotics (Basel). 2024;13(5):419. doi:10.3390/antibiotics13050419"
  vignette <- "Pokorna_2024_meropenem"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Table 2: 'CL [L/h] = CLp x (BW/7.88)' and 'Vd [L] = Vp x (BW/7.88) x (1 + thetaCRRT)^CRRT_on',",
        "a linear (exponent 1) scaling of both parameters normalised to the cohort median weight of",
        "7.88 kg (Table 1: median 7.88 kg, IQR 3.62-11.97). Time-varying in the source data (Methods 4.4:",
        "'Multiple time-varying measurements were attainable for all continuous covariates')."
      ),
      source_name = "BW"
    ),
    RRT_CRRT_ACTIVE = list(
      description = "Continuous renal replacement therapy running (1) or not (0) at the time of the observation",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "Source column CRRT_on, defined in the Table 2 footnote as 'binary parameter indicating whether",
        "CRRT treatment was on (1) or off (0)'. Time-varying: Results 2.1 counts 91 concentrations",
        "drawn during CRRT and 61 off CRRT, and Table S1 splits concentrations by modality at the time",
        "of sampling, so the flag gates sampling time rather than the patient. CRRT was continuous",
        "veno-venous hemodialysis (CVVHD) or hemodiafiltration (CVVHDF) on a Prismaflex circuit",
        "connected pre-pump to the ECMO circuit when both were running; 19 of 31 CRRT patients were",
        "treated for fluid overload and 12 for acute kidney injury. Enters V as (1 + 1.04)^CRRT_on,",
        "i.e. about a two-fold increase; the authors attribute it to fluid overload (Discussion)."
      ),
      source_name = "CRRT_on"
    )
  )

  covariatesDataExcluded <- list(
    PNA = list(
      description = "Postnatal age",
      units = "days",
      type = "continuous",
      notes = "Maturation covariate tested (Methods 4.4); not retained. Table 1: median 3 days (IQR 0-465)."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "umol/L",
      type = "continuous",
      notes = "Tested (Methods 4.4); not retained. Table 1: median 45 (IQR 29-67) umol/L."
    ),
    CRCL = list(
      description = "Estimated creatinine clearance (Schwartz)",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      notes = "Tested (Methods 4.4); not retained. Table 1: median 63 (IQR 34-113) mL/min/1.73 m^2."
    ),
    BUN = list(
      description = "Serum urea",
      units = "mmol/L",
      type = "continuous",
      notes = "Tested as serum urea (Methods 4.4); not retained. Table 1: median 5.6 (IQR 4.2-9) mmol/L."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      notes = "Tested both over the whole course and at the start of treatment (Methods 4.4); not retained. Table 1: median 26 (IQR 22-30) g/L; at treatment start 23 (22-26) g/L."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units = "umol/L",
      type = "continuous",
      notes = "Tested (Methods 4.4); not retained. Table 1: median 15 (IQR 6-49) umol/L."
    ),
    CRP = list(
      description = "C-reactive protein",
      units = "mg/L",
      type = "continuous",
      notes = "Tested (Methods 4.4); not retained. Table 1: median 49 (IQR 18-112) mg/L."
    ),
    AST = list(
      description = "Aspartate transaminase",
      units = "IU/L",
      type = "continuous",
      notes = "Tested (Methods 4.4); not retained. Table 1 prints median 0.96 (IQR 0.53-2.18) 'IU/L'."
    ),
    ALT = list(
      description = "Alanine transaminase",
      units = "IU/L",
      type = "continuous",
      notes = "Tested (Methods 4.4); not retained. Table 1 prints median 0.66 (IQR 0.31-1.33) 'IU/L'."
    ),
    CONMED_DIURETIC = list(
      description = "Concomitant diuretic",
      units = "(binary)",
      type = "binary",
      notes = "Tested as a categorical covariate (Methods 4.4); not retained. Table 1: 38/45 (84%)."
    ),
    CONMED_INOTROPE = list(
      description = "Concomitant inotropic drug",
      units = "(binary)",
      type = "binary",
      notes = "Tested as a categorical covariate (Methods 4.4); not retained. Table 1: 41/45 (91%)."
    ),
    ECMO_STATUS = list(
      description = "ECMO on (1) or off (0)",
      units = "(binary)",
      type = "binary",
      notes = "On/off ECMO, ECMO modality (veno-venous vs veno-arterial) and ECMO circuit change were tested as categorical covariates (Methods 4.4); none retained. 38/45 patients had ECMO (Table 1); 94 of 152 concentrations were drawn on ECMO (Table S1)."
    ),
    ECMO_PUMP_SPEED = list(
      description = "ECMO pump speed",
      units = "RPM",
      type = "continuous",
      notes = "Tested (Methods 4.4); not retained. Table 1: median 2400 (IQR 2150-2800) revolutions/min."
    ),
    Q_ECMO = list(
      description = "ECMO blood flow",
      units = "L/min",
      type = "continuous",
      notes = "Tested (Methods 4.4); not retained. Table 1 prints median 360 (IQR 300-680) with the unit 'L/min', which reads as mL/min for this cohort."
    ),
    T_ECMO = list(
      description = "Duration of ECMO treatment",
      units = "h",
      type = "continuous",
      notes = "Tested together with time after the start and stop of ECMO (Methods 4.4); not retained. Table 1: ECMO length median 230 (IQR 107-251) h."
    )
  )

  compartmentData <- list(
    central = list(analyte = "meropenem", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 45L,
    n_studies = 1L,
    n_concentrations = 152L,
    age_range = "0-5 years",
    age_median = "3 days (IQR 0-465 days) at the start of therapy",
    weight_median = "7.88 kg (IQR 3.62-11.97 kg)",
    sex_female_pct = 48.9,
    race_ethnicity = "Not reported (Swedish paediatric ICU cohort)",
    disease_state = paste(
      "Critically ill neonates (27, 60%), infants (4, 9%) and children (14, 31%) treated with",
      "meropenem in the paediatric ICU / ECMO centre. 38 on ECMO (31 veno-arterial, 7 veno-venous;",
      "HILITE 2400 / Quadrox 7000 oxygenators, DP3 / Centrimag pumps), 31 on CRRT (CVVHD or",
      "CVVHDF; 19 for fluid overload, 12 for acute kidney injury). ECMO indications: pneumonia, ARDS,",
      "septic shock, respiratory insufficiency, meconium aspiration, congenital diaphragmatic hernia."
    ),
    renal_function = "Serum creatinine median 45 (IQR 29-67) umol/L; Schwartz estimated creatinine clearance median 63 (IQR 34-113) mL/min/1.73 m^2 (Table 1).",
    dose_range = paste(
      "Off ECMO: 20-40 mg/kg IV over 15-30 min 2-4 times a day. On ECMO: 20 mg/kg loading dose",
      "(maximum 2 g) then continuous infusion of 30-80 mg/kg/day (0-2 months) or 60-120 mg/kg/day",
      "(3 months-12 years). Doses adjusted by TDM toward 20-40 mg/L in 33 cases (Methods 4.2)."
    ),
    regions = "Sweden (Astrid Lindgren Children's Hospital, Karolinska University Hospital, Stockholm)",
    notes = paste(
      "Retrospective TDM study, 1 January 2018 - 31 December 2020; demographics from Table 1.",
      "152 plasma concentrations (0.68-75 mg/L) by LC-MS/MS, 94 on ECMO and 91 on CRRT (Table S1).",
      "NONMEM 7.3.0 with PsN 3.4.2. Blood pH and the PIM3 score (median 30.5, IQR 13.4-66.1) were also",
      "screened and not retained."
    )
  )

  ini({
    # Table 2 'Final Model (RSE %)' column. Both typical values refer to a
    # patient of the median 7.88 kg off CRRT.
    lcl <- log(1.09)
    label("Clearance CL at 7.88 kg (L/h)") # Table 2: CLp = 1.09 L/h (RSE 8%; bootstrap 1.09, 0.98-1.21)
    lvc <- log(3.98)
    label("Volume of distribution V at 7.88 kg off CRRT (L)") # Table 2: Vp = 3.98 L (RSE 14%; bootstrap 4.01, 2.22-7.30)
    e_crrt_vc <- 1.04
    label("Fractional increase in V while CRRT is running (unitless)") # Table 2: thetaCRRT = 1.04 (RSE 76%; bootstrap 1.06, 0.14-2.20)

    # Table 2 prints the IIV rows as 'CL (%)' 0.0887 and 'Vd (%)' 0.916 and
    # the residual row as 'Proportional' 0.17. These are carried as NONMEM
    # variances: as SDs (or CV fractions) the CL variability would be ~9%
    # and the residual 17%, which cannot produce the roughly two-fold
    # observed-vs-PRED scatter of Figure 2, whereas the variance reading
    # (CL omega 0.30, residual SD 0.41) does.
    etalcl ~ 0.0887 # Table 2: IIV CL = 0.0887 (RSE 48%; bootstrap 0.0882, 0.0429-0.1510)
    etalvc ~ 0.916 # Table 2: IIV Vd = 0.916 (RSE 68%; bootstrap 0.800, 0.060-1.729)

    propSd <- 0.41231
    label("Proportional residual error (fraction)") # Table 2: proportional sigma^2 = 0.17 (RSE 17%); sqrt(0.17) = 0.41231
  })
  model({
    # Table 2: CL = CLp x (BW/7.88); Vd = Vp x (BW/7.88) x (1 + thetaCRRT)^CRRT_on.
    # With a 0/1 indicator (1 + theta)^CRRT_on equals 1 + theta x CRRT_on.
    cl <- exp(lcl + etalcl) * (WT / 7.88)
    vc <- exp(lvc + etalvc) * (WT / 7.88) * (1 + e_crrt_vc * RRT_CRRT_ACTIVE)
    kel <- cl / vc

    d/dt(central) <- -kel * central

    # Dose in mg, V in L -> mg/L.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
