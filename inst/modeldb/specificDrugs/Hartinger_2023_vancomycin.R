Hartinger_2023_vancomycin <- function() {
  description <- paste(
    "Two-compartment population PK model for intraperitoneally administered",
    "vancomycin in adults with end-stage renal disease on continuous",
    "ambulatory peritoneal dialysis (CAPD) treated for PD-associated",
    "peritonitis (Hartinger 2023). The first compartment is the peritoneal",
    "cavity, whose volume V1 is not estimated but fixed per subject to the",
    "actual instilled dialysate volume (1-2 L); the second is the central",
    "(systemic) compartment V2. Vancomycin exchanges bidirectionally between",
    "the two at intercompartmental clearance Q = 0.544 L/h and is eliminated",
    "first-order from central at CL. Clearance is 0.192 L/h in oliguric",
    "patients and is multiplied by (1 + 1.26) * (eGFR / 6.76) in patients with",
    "preserved residual diuresis (> 500 mL/day), so residual renal function",
    "raises CL 2.26-fold at the cohort-median eGFR. Central volume is linear",
    "in body weight on the natural scale: V2 = 23.6 + 50.9 * (BW / 75) L,",
    "i.e. 74.5 L at the median 75 kg. Interindividual variability is",
    "log-normal on CL (34.2% CV), V2 (30.3% CV) and Q (50.7% CV), with",
    "additional inter-occasion variability on Q (31% CV) across peritonitis",
    "episodes -- peritoneal membrane permeability changes with the severity of",
    "each inflammatory episode. Residual error is proportional on peritoneal",
    "dialysate concentrations and combined additive + proportional on plasma.",
    "There is no drainage term in the ODEs: the end-of-dwell drain and the",
    "next instillation are dosing EVENTS on the peritoneum compartment (see",
    "the validation vignette), not model parameters."
  )
  reference <- paste(
    "Hartinger JM, Michalickova D, Dvorackova E, Hronova K, Krekels EHJ,",
    "Szonowska B, Bednarova V, Benakova H, Kroneislova G, Zavora J,",
    "Tesar V, Slanar O.",
    "Intraperitoneally Administered Vancomycin in Patients with Peritoneal",
    "Dialysis-Associated Peritonitis: Population Pharmacokinetics and Dosing",
    "Implications.",
    "Pharmaceutics. 2023;15(5):1394.",
    "doi:10.3390/pharmaceutics15051394"
  )
  vignette <- "Hartinger_2023_vancomycin"

  # The peritoneal cavity is a paper-mechanistic state that is not in the
  # canonical compartment register. Precedent: Royer_2011_cisplatin.R declares
  # the same state (also called `peritoneum`) as paper-specific for
  # intraperitoneal perioperative chemotherapy. Hartinger 2023 is the second
  # independent paper in the library to carry an intraperitoneal-cavity state,
  # which per the operator's standing ruling is the trigger to consider
  # promoting `peritoneum` to a canonical compartment; that promotion is a
  # register decision and is proposed in the PR rather than taken unilaterally.
  paper_specific_compartments <- c("peritoneum")

  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. `peritoneum` is both the intraperitoneal dosing site and
  # a sampled matrix (drained dialysate is assayed), but "peritoneal dialysate"
  # is not in conventions$specimenVocabulary; "administration site" follows the
  # Royer_2011_cisplatin.R precedent for the same anatomical state.
  compartmentData <- list(
    peritoneum = list(analyte = "vancomycin", units = "mg", specimen = "administration site", verified = TRUE),
    central    = list(analyte = "vancomycin", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Actual body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Cohort median 75 kg (IQR 70-84) per Hartinger 2023 Table 1. Enters the central volume as a linear-on-the-natural-scale ratio term, V2i = V2p + theta_BWV * (BW / 75) (Hartinger 2023 Table 2 row 'V2i [L] = V2p + theta BWV x (BW/75)'). Note this is a through-origin RATIO, not a centred difference: at the median 75 kg the typical V2 is 23.6 + 50.9 = 74.5 L, not 23.6 L. The Results text quotes '23.6 L (78%)' as the typical V2 for a 75 kg patient, but that number is the Table 2 population intercept V2p together with its RSE, and the additive-ratio equation is the one that reproduces the paper's own simulations (see the vignette source-trace: the day-1 plasma level of ~16 mg/L after a 25 mg/kg intraperitoneal loading dose in Figure 3E/3F requires V2 near 75 L, and would be ~50 mg/L at 23.6 L). The blow-up of the V2p RSE from 9% to 78% on adding this covariate, which the paper reports and retains deliberately, is itself the signature of an intercept that carries only a third of the typical volume. Simulated range 50-100 kg (Hartinger 2023 Figures 3, 4, S6, S7).",
      source_name        = "BW"
    ),
    CRCL = list(
      description        = "Residual estimated glomerular filtration rate, CKD-EPI 2009 equation",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Cohort median 6.76 mL/min/1.73 m^2 (IQR 5.07-7.92) per Hartinger 2023 Table 1, which states the value is 'Calculated according to CKD-EPI 2009 formula, only in patients with preserved diuresis.' Enters clearance as the through-origin ratio (CRCL / 6.76) INSIDE the preserved-diuresis bracket, so it acts only when URINE_VOL_24H > 500 mL/day (Hartinger 2023 Table 2: CLi = CLp * ((1 + theta_RESDIU) * (CRCL/6.76))^(RESDIU>500)); this is the categorical-gates-a-continuous interaction the paper describes as 'eGFR was included as a covariate on CL in a linear relationship FOR PATIENTS WITH PRESERVED DIURESIS'. Because the exponent is 0 in oliguric subjects the column is unused there and any positive placeholder (e.g. the median 6.76) reproduces the model exactly. Hartinger 2023 Table 2 writes the column 'CRCL' while the abbreviation list and Methods both define it as CKD-EPI eGFR; the eGFR reading is authoritative.",
      source_name        = "CRCL"
    ),
    URINE_VOL_24H = list(
      description        = "Residual diuresis, total 24-hour urine volume",
      units              = "mL/24h",
      type               = "continuous",
      reference_category = "n/a -- used as the binary preserved-diuresis gate (URINE_VOL_24H > 500)",
      notes              = "Hartinger 2023 Methods 'Covariate analysis' defines the categorical covariate as 'Preserved diuresis (yes = over 500 mL urine daily or no = less than 500 mL urine daily)'; Table 2 calls the resulting indicator RESDIU>500. Cohort: 31/41 (75.6%) with preserved diuresis > 500 mL/day, 10/41 (24.4%) oliguric < 500 mL/day; median residual diuresis 1000 mL/day (IQR 350-1450) per Table 1. The indicator is the EXPONENT on the whole bracket ((1 + theta_RESDIU) * (CRCL/6.76)), so it simultaneously switches on the 2.26-fold clearance increase and the eGFR ratio. The 500 mL/24h threshold is the preserved-diuresis cutoff already named in this register entry's reference-category note; Huppe 2023 uses the same column at the 100 mL/24h anuria cutoff instead.",
      source_name        = "RESDIU"
    ),
    RRT_PERIT_DIAL_FILL_VOLUME = list(
      description        = "Volume of peritoneal dialysate instilled at the start of the dwell",
      units              = "L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Hartinger 2023 Table 2 fixes V1, the volume of the peritoneal compartment, to 'the actual volume of peritoneal solution used (range 1-2 L)' rather than estimating it, so this column IS the peritoneal compartment volume and appears directly in model() as v1. Observed exchange volumes (Table 1): 2000 mL in 23/41 patients (56%), 1500 mL in 16/41 (39%), 1200 mL and 1000 mL in one patient each (2% each). The paper's simulations use 2 L for 75 kg and 100 kg patients and 1.5 L for 50 kg patients 'due to the smaller peritoneal cavity' (Results 3.2). Confirmed against the paper's own figure axes: a 20 mg/kg intraperitoneal loading dose in a 75 kg patient (1500 mg) reads exactly 750 mg/L on the Figure 4 peritoneal-concentration axis, i.e. dose / 2 L.",
      source_name        = "V1"
    ),
    OCC = list(
      description        = "Integer-valued occasion indicator for inter-occasion-variability multiplexing",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = "Hartinger 2023 Methods defines an occasion as 'a particular vancomycin treatment course that ended up with vancomycin withdrawal due to the cure of the infection, switch to a different antibiotic when ineffective, or the removal of the PD catheter' -- i.e. one peritonitis episode. The data comprise 57 treatment occasions in 41 patients. The paper reports one shared IOV magnitude on Q (31% CV) and never states a maximum occasion count, so three occasions are encoded here (as with Ding 2026 vancomycin and Stoschus 2025 phenobarbital, where the source likewise gives the IOV magnitude without an occasion count); occasions 2 and 3 are fix()'d to the occasion-1 variance to reproduce NONMEM $OMEGA BLOCK(1) SAME. A single-episode simulation uses OCC = 1 throughout.",
      source_name        = "OCC"
    )
  )

  covariatesDataExcluded <- list(
    LBM = list(
      description = "Lean body mass",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened on both V2 and CL in the stepwise covariate search (Hartinger 2023 Methods 'Covariate analysis'); not retained. Cohort median 54.31 kg (IQR 46.18-59.89) per Table 1."
    ),
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened on both V2 and CL; not retained. Cohort median 68 years (IQR 53-74)."
    ),
    ALB = list(
      description = "Serum albumin at the start of vancomycin treatment",
      units       = "g/L",
      type        = "continuous",
      notes       = "Screened on both V2 and CL; not retained. Cohort median 28 g/L (IQR 26-31)."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened on V2, CL and Q; not retained. Cohort 17 female / 24 male."
    ),
    UREA = list(
      description = "Serum urea at the start of vancomycin treatment",
      units       = "mmol/L",
      type        = "continuous",
      notes       = "Screened on CL; not retained. Cohort median 18 mmol/L (IQR 14.9-22.3)."
    ),
    SCR = list(
      description = "Serum creatinine at the start of vancomycin treatment",
      units       = "umol/L",
      type        = "continuous",
      notes       = "Screened on CL; not retained. Cohort median 694 umol/L (IQR 564-849). Used upstream to compute CRCL via CKD-EPI 2009."
    ),
    CRP = list(
      description = "C-reactive protein at the start of the peritonitis treatment",
      units       = "mg/L",
      type        = "continuous",
      notes       = "Screened on Q as a marker of peritoneal inflammation; not retained. Cohort median 31.7 mg/L (IQR 10.9-96.5). Hartinger 2023 Discussion notes that a more sensitive inflammation marker such as procalcitonin might have correlated with intercompartmental clearance where CRP did not."
    ),
    POTASSIUM = list(
      description = "Serum potassium at the start of vancomycin treatment",
      units       = "mmol/L",
      type        = "continuous",
      notes       = "Screened on Q; not retained. Cohort median 4.1 mmol/L (IQR 3.8-4.6)."
    ),
    PERIT_DIAL_SOLUTION = list(
      description = "Type of peritoneal dialysis solution (low / medium / high glucose content, or icodextrin-based)",
      units       = "(categorical)",
      type        = "categorical",
      notes       = "Screened on Q as a categorical covariate; not retained. Not registered as a canonical column because no model in the library retains it; the solution type does matter analytically -- glucose-based solutions depressed the nephelometric vancomycin assay by ~20%, and all dialysate concentrations were multiplied by 1.2885 before modelling (Hartinger 2023 Methods 2.2 and Figure S1)."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 41L,
    n_studies      = 1L,
    n_occasions    = "57 treatment occasions (hospitalizations), including recurrent and relapsing peritonitis episodes",
    n_observations = "373 vancomycin concentrations (132 peritoneal dialysate + 241 plasma)",
    age_range      = "median 68 years (IQR 53-74)",
    weight_range   = "median 75 kg (IQR 70-84); BMI median 26.67 kg/m^2 (IQR 21.6-28.34); lean body weight median 54.31 kg (IQR 46.18-59.89)",
    sex_female_pct = 41.5,
    race_ethnicity = "All but one patient were Caucasian; the paper states ethnic origin could not be tested as a covariate for this reason.",
    disease_state  = "End-stage renal disease treated by peritoneal dialysis, with suspected or culture-confirmed PD-associated peritonitis. Cultivations from peritoneal dialysate: G+ 36 (63%), culture-negative 13 (23%), G- 4 (7%), mixed 3 (5%), G+ with candida 1 (2%). Concomitant exit-site infection in 27 (47%). Median time on PD 22.5 months (IQR 9-46.75).",
    renal_function = "Residual diuresis > 500 mL/day (preserved) in 31/41 (75.6%) and < 500 mL/day (oliguria) in 10/41 (24.4%); median residual diuresis 1000 mL/day (IQR 350-1450). Residual eGFR (CKD-EPI 2009, computed only in patients with preserved diuresis) median 6.76 mL/min/1.73 m^2 (IQR 5.07-7.92). Serum creatinine median 694 umol/L (IQR 564-849).",
    dose_range     = "Intraperitoneal vancomycin per ISPD guidance, predominantly a 15-30 mg/kg loading dose followed by 25 mg per litre of instilled dialysate in every subsequent exchange, with maintenance doses adjusted by therapeutic drug monitoring to a plasma AUC24 of 400-600 mg*h/L. Intravenous dosing was combined with intraperitoneal dosing when systemic infection was present; 3 patients received intravenous vancomycin only and contributed peritoneal concentrations. On admission all patients were switched to CAPD with 4-5 manual exchanges per day; exchange volumes 2000 mL (56%), 1500 mL (39%), 1200 mL (2%), 1000 mL (2%).",
    regions        = "Czech Republic (two nephrology departments of the General University Hospital in Prague), June 2016 - August 2022",
    notes          = "Open-label retrospective observational therapeutic-drug-monitoring study. Vancomycin assayed by nephelometry (Beckman Coulter); the plasma analytical range was 3.5-40 mg/L with dilution above it. The method was validated for dialysate over 10-250 mg/L: icodextrin did not interfere but glucose-based solutions depressed the measured value by ~20% independent of glucose or vancomycin concentration, so all dialysate concentrations were multiplied by 1.2885 before modelling (Methods 2.2, Figure S1). NONMEM 7.4.0 with PsN 3.4.2 under Pirana 2.9.0, FOCE-I. Model validated by a 500-replicate nonparametric bootstrap (all medians except V2p within 10% of the final estimates) and by normalized prediction distribution errors from 1000 simulated datasets (plasma NPDE mean 0.1027 / variance 0.848; peritoneal 0.1228 / 0.8438; neither significantly different from 0 and 1)."
  )

  ini({
    # ------------------------------------------------------------------
    # Structural parameters -- Hartinger 2023 Table 2, 'Final Model (RSE %)'
    # column. The peritoneal compartment volume V1 is NOT estimated: Table 2
    # records it as 'FIXED *' with the footnote 'Fixed to the actual volume of
    # peritoneal solution used (range 1-2 L)', so it enters model() directly
    # from the RRT_PERIT_DIAL_FILL_VOLUME data column and has no ini() entry.
    # ------------------------------------------------------------------
    lcl <- log(0.192) ; label("Population clearance CLp from the central compartment (L/h; oliguric reference)") # Hartinger 2023 Table 2: CLp = 0.192 L/h (RSE 17%); bootstrap median 0.186 (0.144-0.267)
    lvc <- log(23.6)  ; label("Population central volume intercept V2p (L)")                                      # Hartinger 2023 Table 2: V2p = 23.6 L (RSE 78%); bootstrap median 27.2 (1.68-65.3)
    lq  <- log(0.544) ; label("Intercompartmental clearance Q between peritoneum and central (L/h)")              # Hartinger 2023 Table 2: Q = 0.544 L/h (RSE 16%); bootstrap median 0.53 (0.40-0.76)

    # ------------------------------------------------------------------
    # Covariate effects -- Hartinger 2023 Table 2.
    #
    # CL: CLi = CLp * ((1 + theta_RESDIU) * (CRCL / 6.76))^(RESDIU>500). The
    # preserved-diuresis indicator is the EXPONENT on the whole bracket, so in
    # oliguric patients (exponent 0) CL is exactly CLp and the eGFR term drops
    # out entirely. At the cohort-median eGFR of 6.76 the bracket is
    # (1 + 1.26) = 2.26, i.e. CL = 0.434 L/h with preserved diuresis -- close
    # to the 0.443 L/h non-dialysis vancomycin clearance the Discussion cites
    # for high-flux haemodialysis patients.
    #
    # V2: V2i = V2p + theta_BWV * (BW / 75) -- an additive term in a
    # through-origin weight RATIO, so both parameters contribute at the median
    # weight (23.6 + 50.9 = 74.5 L at 75 kg). See the WT covariateData note for
    # why the printed equation, not the Results sentence quoting '23.6 L', is
    # the one encoded.
    # ------------------------------------------------------------------
    e_urine_vol_24h_cl <- 1.26 ; label("Fractional clearance increase with preserved diuresis, theta_RESDIU (unitless)") # Hartinger 2023 Table 2: theta_RESDIU = 1.26 (RSE 31%); bootstrap median 1.36 (0.50-2.25)
    e_wt_vc            <- 50.9 ; label("Additive effect of the weight ratio BW/75 on central volume, theta_BWV (L)")     # Hartinger 2023 Table 2: theta_BWV = 50.9 (RSE 38%); bootstrap median 45.9 (9.6-76.0)

    # ------------------------------------------------------------------
    # Interindividual variability. Hartinger 2023 Table 2 reports IIV as %CV
    # under 'Inter-individual variability'; the log-normal log-scale variance
    # is omega^2 = log(1 + CV^2). No IIV was estimated on the peritoneal volume
    # V1 (Methods: IIV was tested on each PK parameter 'with the exception of
    # the volume of the peritoneal compartment (V1)'), which is consistent with
    # V1 being fixed to a measured per-subject quantity.
    # ------------------------------------------------------------------
    etalcl ~ 0.1106143 # log(1 + 0.342^2); Hartinger 2023 Table 2: IIV CL = 34.2% CV (RSE 22%); bootstrap 32.2% (16.0-45.8%)
    etalvc ~ 0.0878360 # log(1 + 0.303^2); Hartinger 2023 Table 2: IIV V2 = 30.3% CV (RSE 21%); bootstrap 29.2% (13.8-40.4%)
    etalq  ~ 0.2287669 # log(1 + 0.507^2); Hartinger 2023 Table 2: IIV Q  = 50.7% CV (RSE 23%); bootstrap 45.3% (20.1-79.6%)

    # ------------------------------------------------------------------
    # Inter-occasion variability on Q. Hartinger 2023 Table 2 row 'IOC on Q(%)'
    # = 31% CV (RSE 37%), bootstrap 30.1% (17.3-51.7%). An occasion is one
    # peritonitis treatment course; the paper adds IOV on Q because peritoneal
    # membrane permeability changes with the severity of each inflammatory
    # episode (Results, Figure 2). One shared magnitude is reported, so
    # occasions 2 and 3 are fix()'d to the occasion-1 variance to reproduce
    # NONMEM $OMEGA BLOCK(1) SAME (nlmixr2 has no SAME shortcut).
    # ------------------------------------------------------------------
    etaiov_q_1 ~ 0.0917584      # log(1 + 0.31^2); Hartinger 2023 Table 2: IOC on Q = 31% CV (estimated, occasion 1)
    etaiov_q_2 ~ fix(0.0917584) # occasion-2 variance equal to occasion-1 ($OMEGA BLOCK(1) SAME pattern)
    etaiov_q_3 ~ fix(0.0917584) # occasion-3 variance equal to occasion-1 ($OMEGA BLOCK(1) SAME pattern)

    # ------------------------------------------------------------------
    # Residual error -- Hartinger 2023 Table 2, 'Residual variability' block.
    # The paper reports proportional error for peritoneal concentrations and a
    # combined additive + proportional error for plasma (Results 3.1). The
    # three values are printed as bare numbers while every IIV row in the same
    # table is printed as a percentage, and the analysis is NONMEM 7.4.0, whose
    # $SIGMA reports VARIANCES -- so the tabulated values are read as variances
    # and converted to the SD scale nlmixr2 expects. See the vignette
    # Assumptions and deviations section for the alternative reading and why it
    # is rejected.
    # ------------------------------------------------------------------
    propSd_Cip <- 0.301662 ; label("Proportional residual error, peritoneal dialysate concentration (fraction)") # sqrt(0.091);   Hartinger 2023 Table 2: proportional error, peritoneal concentration = 0.091 (RSE 15%); bootstrap 0.089 (0.064-0.146)
    addSd      <- 2.172556 ; label("Additive residual error SD, plasma concentration (mg/L)")                    # sqrt(4.72);    Hartinger 2023 Table 2: additive error, plasma concentration = 4.72 (RSE 31%); bootstrap 4.74 (1.38-7.11)
    propSd     <- 0.077717 ; label("Proportional residual error, plasma concentration (fraction)")               # sqrt(0.00604); Hartinger 2023 Table 2: proportional error, plasma concentration = 0.00604 (RSE 54%); bootstrap 0.00649 (0.000450-0.01489)
  })

  model({
    # Decompose the integer occasion column into binary indicators that
    # multiplex the three IOV etas onto log-Q. An occasion is one peritonitis
    # treatment course (Hartinger 2023 Methods 'Development of structural and
    # statistical model').
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    iov_q <- oc1 * etaiov_q_1 + oc2 * etaiov_q_2 + oc3 * etaiov_q_3

    # Preserved-diuresis indicator, RESDIU>500 in Hartinger 2023 Table 2:
    # 1 when residual diuresis exceeds 500 mL/day, 0 otherwise (Methods
    # 'Covariate analysis').
    resdiu <- (URINE_VOL_24H > 500)

    # Individual parameters.
    #
    # Clearance: the preserved-diuresis indicator is the exponent on the whole
    # bracket, so the eGFR ratio acts only in patients with residual renal
    # function. 6.76 mL/min/1.73 m^2 is the cohort-median eGFR the ratio is
    # normalised to (Hartinger 2023 Table 1 and Table 2).
    cl <- exp(lcl + etalcl) * ((1 + e_urine_vol_24h_cl) * (CRCL / 6.76))^resdiu

    # Central volume: additive linear-in-weight-ratio term on the natural
    # scale, Hartinger 2023 Table 2 'V2i [L] = V2p + theta BWV x (BW/75)'.
    vc <- (exp(lvc) + e_wt_vc * (WT / 75)) * exp(etalvc)

    # Intercompartmental clearance between the peritoneal cavity and the
    # central compartment, carrying both IIV and IOV.
    q <- exp(lq + etalq + iov_q)

    # Peritoneal compartment volume. Not estimated: Hartinger 2023 Table 2
    # fixes V1 to the actual instilled dialysate volume for each subject, so it
    # is read straight from the data column.
    v1 <- RRT_PERIT_DIAL_FILL_VOLUME

    # Concentrations. Cip is the vancomycin concentration in the peritoneal
    # dialysate (the quantity assayed in the drained-out fluid); Cc is the
    # plasma concentration.
    Cip <- peritoneum / v1
    Cc  <- central / vc

    # ODE system -- Hartinger 2023 Figure 1. Vancomycin exchanges between the
    # peritoneal cavity and the central compartment down its concentration
    # gradient at intercompartmental clearance Q (so drug diffuses out of the
    # blood into fresh dialysate whenever plasma exceeds peritoneal
    # concentrations, as the Introduction describes), and is cleared
    # first-order from central at CL. The instillation and end-of-dwell
    # drainage of dialysate (VOL IN / VOL OUT in Figure 1) are handled as
    # dosing and replacement EVENTS on the peritoneum compartment in the data,
    # not as terms here.
    d/dt(peritoneum) <- q * (Cc - Cip)
    d/dt(central)    <- q * (Cip - Cc) - cl * Cc

    # Observation model. Proportional residual error on peritoneal dialysate
    # concentrations, combined additive + proportional on plasma
    # (Hartinger 2023 Results 3.1 and Table 2).
    Cip ~ prop(propSd_Cip)
    Cc  ~ add(addSd) + prop(propSd)
  })
}
