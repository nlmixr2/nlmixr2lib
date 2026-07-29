Fattinger_1991_quinidine <- function() {
  description <- "Two-compartment population PK model for oral quinidine in adults treated for supraventricular or ventricular arrhythmias (Fattinger 1991). Zero-order absorption from the gastrointestinal tract with formulation-specific absorption duration: immediate-release quinidine sulphate (Chinidin sulfuricum) with a typical absorption duration of 1.37 h, and slow-release quinidine bisulphate (Kinidin duriles) with a typical absorption duration of 6.0 h and a 1.36-fold higher relative bioavailability versus quinidine sulphate. Apparent total clearance is the sum of a renal arm proportional to creatinine clearance (proportionality 0.0566 L/h per mL/min) and a non-renal arm of 12.6 L/h that is halved to 6.8 L/h in patients with severe heart failure or severe liver failure. Apparent central volume is 161 L. Inter-compartmental clearance Q is 12.6 L/h and peripheral volume V2 is 66.7 L. Inter-individual variability is assigned to total clearance (40.2% CV), central volume (75.6% CV), and the quinidine sulphate absorption duration (49.4% CV); residual variability is proportional (22% CV)."
  reference <- paste(
    "Fattinger K, Vozeh S, Ha HR, Borner M, Follath F (1991).",
    "Population pharmacokinetics of quinidine.",
    "Br J Clin Pharmacol 31(3):279-286.",
    "doi:10.1111/j.1365-2125.1991.tb05531.x.",
    sep = " "
  )
  vignette <- "Fattinger_1991_quinidine"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    CRCL = list(
      description        = "Creatinine clearance, raw modified Cockcroft-Gault (NOT BSA-normalized)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Estimated from serum creatinine via the modified Cockcroft-Gault equation reported on page 281: CLcr = (150 - age) * Body weight / Scr, with +10% for males and -10% for females (Dettli 1983). Used as a linear effect on the renal-clearance arm of total apparent CL: CL_renal = e_crcl_cl_renal * CRCL (L/h per mL/min). Source column 'CLcr' renamed to canonical 'CRCL' on input. The published cohort had median CRCL 62.5 mL/min (range 17 to >100 mL/min, page 281).",
      source_name        = "CLcr"
    ),
    DIS_HF_OR_LF_SEV = list(
      description        = "Severe heart failure OR severe liver failure pooled indicator (1 = either present; 0 = neither)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no severe heart failure and no severe liver failure)",
      notes              = "Time-fixed. Severe heart failure (n = 2) was defined as low cardiac output or pulmonary oedema; severe liver failure (n = 3) was serum bilirubin > 30 umol/L AND prothrombin time < 60% of normal (page 281). The two conditions were pooled into a single covariate because their estimated effects on non-renal CL were of similar magnitude and patient counts were small (page 282). Mild and moderate heart or liver dysfunction did not show an effect and are excluded from this indicator. Drives a multiplicative reduction of the non-renal CL arm from 12.6 L/h to 6.8 L/h (a factor of 0.5397) per Table 1.",
      source_name        = "(pooled severe HF or severe LF)"
    ),
    FORM_QUIN_SR = list(
      description        = "Slow-release quinidine bisulphate (Kinidin duriles, Astra) vs immediate-release quinidine sulphate (Chinidin sulfuricum, Siegfried) formulation indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (immediate-release quinidine sulphate; the typical-value absorption-duration and bioavailability reference)",
      notes              = "Per-dose-occasion indicator: in the Fattinger 1991 cohort 45/60 patients received both formulations across the dataset and 15/60 received only one (page 280 Figure 2). When a subject receives both formulations, FORM_QUIN_SR is set on each dose record. Drives the structural switch between (a) durations of zero-order absorption (1.37 h for QS reference; 6.0 h for QBS) and (b) relative bioavailability (1 for QS reference; 1.36 for QBS). The per-subject IIV on the QS absorption duration applies only when FORM_QUIN_SR = 0; the QBS duration carries no IIV per Methods page 282. Doses must be entered in mg of quinidine base (apply the paper's stoichiometric factors of 0.829 mg base per mg quinidine sulphate and 0.663 mg base per mg quinidine bisulphate, Windholz 1983, before passing to the model; the model's bioavailability term then captures only the formulation-driven absorption difference).",
      source_name        = "(formulation: QS vs QBS)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 60L,
    n_studies      = 1L,
    n_observations = "260 serum drug concentration measurements in the model-development cohort; 30 separate patients (single concentration each) in the external-validation cohort",
    age_range      = "28-82 years",
    age_median     = "65.5 years",
    weight_range   = "45-105 kg",
    weight_median  = "70.5 kg",
    sex_female_pct = 23.3,
    race_ethnicity = NA_character_,
    disease_state  = "Adults treated for supraventricular or ventricular arrhythmias. Cohort breakdown by hepatic and cardiac status (page 281): severe liver failure n = 3, moderate liver dysfunction n = 22; mild heart failure n = 19, moderate heart failure n = 20, severe heart failure n = 2. Median creatinine clearance 62.5 mL/min (range 17 to >100 mL/min). Nine patients received concomitant nifedipine (5 at 30-40 mg/day, 4 at 60 mg/day); nifedipine had no effect on quinidine PK in this cohort.",
    dose_range     = "Usual regimen: first dose 400 or 600 mg quinidine sulphate orally, followed 3 h later by 500 mg quinidine bisulphate twice or three times daily; 56 of 260 samples obtained under steady state (>= 96 h on unchanged regimen).",
    regions        = "Switzerland (University Hospital, Basel)",
    renal_function = "Median creatinine clearance 62.5 mL/min (range 17 to >100 mL/min); covers mild to severe but not end-stage renal impairment. Apply with caution outside the 17 to >100 mL/min range observed (Discussion page 283).",
    notes          = "Demographics from Methods 'Study population' page 280-281. The model was fit by NONMEM using a stepwise covariate-selection procedure (Methods page 281). Bioavailability of oral quinidine sulphate was assumed to be 100% per the abstract item 3 (so CL and V in this file are apparent CL/F and V/F at the quinidine-sulphate reference). External validation in a separate cohort of 30 consecutive patients (median CrCl 48.7 mL/min, range 19.5 to >100 mL/min, no severe HF/LF) showed 28/30 measured concentrations within the 90% prediction interval (page 282)."
  )

  ini({
    # Final population PK parameters from Table 1 of Fattinger 1991 (page
    # 282). Interindividual variability is reported as coefficient of
    # variation; converted to log-scale variance via omega^2 = log(1 + CV^2)
    # per the lognormal IIV convention. Residual error is reported as 22%
    # CV with the paper stating "A log-additive error distribution was
    # assumed for description of the interindividual variability ... and
    # of the intraindividual variability" (page 281), which is a log-
    # normal proportional residual error in linear space (mapped to
    # nlmixr2's propSd).

    # Renal-clearance arm: proportionality between creatinine clearance and
    # apparent renal clearance of quinidine, in (L/h) per (mL/min).
    # Encoded as a linear-additive covariate effect with the typical value
    # of an intercept-free renal-clearance arm.
    e_crcl_cl_renal <- log(0.0566); label("Log proportionality constant of CRCL on apparent renal CL of quinidine (log of L/h per mL/min)")  # Table 1 'CL_renal' = 0.0566 L/h per mL/min, s.e. 0.0242

    # Non-renal clearance arm in patients without severe heart or liver
    # failure. The Table 1 effect of severe HF or LF is encoded as a
    # log-multiplicative covariate effect that reduces this value by a
    # factor of 0.5397 (= 6.8 / 12.6) when DIS_HF_OR_LF_SEV = 1.
    lcl_nonren            <- log(12.6); label("Apparent non-renal clearance in patients without severe HF/LF (L/h)")            # Table 1 'CL_nonrenal for patients without severe HF and LF' = 12.6 L/h, s.e. 1.8
    e_dis_hf_or_lf_sev_cl_nonren <- log(6.8 / 12.6); label("Log-multiplicative effect of severe HF or LF on apparent non-renal CL (unitless)")  # Derived from Table 1 'CL_nonrenal for patients with severe HF or LF' / 'CL_nonrenal for patients without severe HF/LF' = 6.8 / 12.6 = 0.5397

    # Central volume of distribution.
    lvc <- log(161); label("Apparent central volume of distribution V1 (L)")  # Table 1 'V1' = 161 L, s.e. 14

    # Inter-compartmental clearance and peripheral volume.
    lq  <- log(12.6); label("Apparent inter-compartmental clearance Q (L/h)")          # Table 1 'Q' = 12.6 L/h, s.e. 5.8 (footnote e identifies Q as inter-compartmental clearance)
    lvp <- log(66.7); label("Apparent peripheral volume of distribution V2 (L)")       # Table 1 'V2' = 66.7 L, s.e. 16.4

    # Zero-order absorption durations -- one per formulation.
    # FORM_QUIN_SR = 0 (QS, immediate-release) uses ldur_qs;
    # FORM_QUIN_SR = 1 (QBS, slow-release Kinidin duriles) uses ldur_qbs.
    ldur_qs  <- log(1.37); label("Zero-order absorption duration for quinidine sulphate, immediate-release (h)")  # Table 1 't_max,QS' = 1.37 h, s.e. 0.04
    ldur_qbs <- log(6.00); label("Zero-order absorption duration for quinidine bisulphate, slow-release (h)")     # Table 1 't_max,QBS' = 6.00 h, s.e. 0.25

    # Relative bioavailability of quinidine bisulphate compared with
    # quinidine sulphate (the abstract item 3 and Table 1 footnote g
    # state that QS bioavailability was assumed = 1, so this parameter
    # is the relative scalar applied when FORM_QUIN_SR = 1).
    lfdepot <- log(1.36); label("Log relative bioavailability of quinidine bisulphate vs quinidine sulphate (unitless; QS = 1 by convention)")  # Table 1 'F' = 1.36, s.e. 0.12

    # Interindividual variability. Paper Methods page 281 / Results page
    # 282 assign IIV to three parameters: total clearance, central
    # volume, and the QS absorption duration. No IIV on Q, V2, QBS
    # duration, or F.
    #
    # The single IIV on "clearance" applies multiplicatively to the
    # covariate-driven typical-value total CL (renal + non-renal arms),
    # matching the convention in NONMEM stepwise covariate analyses
    # (Methods page 281: "Interindividual variability was assigned to
    # the following parameters: clearance ...").
    etalcl     ~ log(1 + 0.402^2)   # Table 1 IIV on total CL  = 40.2% CV, s.e. 55%
    etalvc     ~ log(1 + 0.756^2)   # Table 1 IIV on V1        = 75.6% CV, s.e. 53%
    etaldur_qs ~ log(1 + 0.494^2)   # Table 1 IIV on t_max,QS  = 49.4% CV, s.e. 65%

    # Proportional residual error (linear-space CV; the paper's "log-
    # additive error distribution" -- Methods page 281 -- maps to
    # nlmixr2's propSd because additive on log-scale is proportional on
    # linear scale at the small-CV limit and is conventionally encoded
    # as propSd in linear-space proportional residuals).
    propSd <- 0.22; label("Proportional residual error (fraction)")  # Table 1 'sigma' = 22% CV, s.e. 39%
  })

  model({
    # Typical-value clearance arms and individual total CL.
    # The single IIV (etalcl) acts multiplicatively on the combined
    # (renal + non-renal) typical clearance, matching the paper's
    # NONMEM parameterisation in which a log-additive eta is attached
    # to total CL.
    cl_renal_tv  <- exp(e_crcl_cl_renal) * CRCL
    cl_nonren_tv <- exp(lcl_nonren + e_dis_hf_or_lf_sev_cl_nonren * DIS_HF_OR_LF_SEV)
    cl <- (cl_renal_tv + cl_nonren_tv) * exp(etalcl)

    # Central and peripheral volumes; inter-compartmental clearance.
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp)
    q  <- exp(lq)

    # Formulation-dependent zero-order absorption duration. The IIV on
    # QS absorption duration applies only to QS doses (Methods page
    # 282: no improvement from allowing IIV on QBS duration). For
    # subjects who receive both formulations, the dur(central) switches
    # between formulations on a per-dose-record basis driven by
    # FORM_QUIN_SR; the eta is identifiable from the QS-dosing records.
    dur_qs_i <- exp(ldur_qs + etaldur_qs)
    dur_qbs  <- exp(ldur_qbs)
    dur_central <- (1 - FORM_QUIN_SR) * dur_qs_i + FORM_QUIN_SR * dur_qbs

    # Formulation-dependent bioavailability (QS = 1 reference, QBS =
    # exp(lfdepot) = 1.36). No IIV on F per Methods page 282.
    fdepot <- exp(lfdepot)
    f_central <- (1 - FORM_QUIN_SR) + FORM_QUIN_SR * fdepot

    # Two-compartment ODE system. Doses target the central compartment
    # as a zero-order infusion of duration dur_central with
    # bioavailability f_central.
    d/dt(central)     <- -cl / vc * central - q / vc * central + q / vp * peripheral1
    d/dt(peripheral1) <-  q / vc * central - q / vp * peripheral1

    dur(central) <- dur_central
    f(central)   <- f_central

    # Plasma quinidine concentration (mg/L). Dose must be supplied in
    # mg of quinidine BASE (apply the Windholz 1983 stoichiometric
    # factors of 0.829 mg base per mg quinidine sulphate and 0.663 mg
    # base per mg quinidine bisulphate before passing to the model);
    # the f(central) term then captures only the formulation-driven
    # relative absorption difference.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
