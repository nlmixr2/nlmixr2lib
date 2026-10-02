Yu_2020_inotersen <- function() {
  description <- "Two-compartment population PK and indirect-response PD model for the 2'-MOE antisense oligonucleotide inotersen targeting transthyretin (TTR) mRNA, fit to pooled data from a phase 1 study in healthy volunteers and the phase 2/3 NEURO-TTR study plus its open-label extension in patients with hereditary transthyretin amyloidosis with polyneuropathy (Yu 2020). First-order SC absorption, linear power-form scaling of CL/F, Q/F, Vc/F, and Vp/F on lean body mass (exponents fixed at 1), proportional disease-status effects on CL/F and Vc/F, and an indirect-response model with inotersen-driven inhibition of TTR production carrying a proportional disease-status effect on baseline TTR."
  reference <- "Yu RZ, Collins JW, Hall S, Ackermann EJ, Geary RS, Monia BP, Henry SP, Wang Y. Population Pharmacokinetic-Pharmacodynamic Modeling of Inotersen, an Antisense Oligonucleotide for Treatment of Patients with Hereditary Transthyretin Amyloidosis. Nucleic Acid Ther. 2020;30(3):153-163. doi:10.1089/nat.2019.0822"
  vignette <- "Yu_2020_inotersen"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot = list(analyte = "inotersen", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "inotersen", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "inotersen", units = "mg", specimen = "tissue", verified = TRUE),
    effect = list(analyte = "transthyretin", units = "mg/dL", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    LBM = list(
      description = "Lean body mass (baseline)",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-fixed baseline. Power exponent fixed at 1 on CL/F and Q/F (Yu 2020 Table 2 row 'Lean body mass * CL/Q') and on Vc/F and Vp/F (row 'Lean body mass * Vc/Vp'), each 'Power-centered on median'; reference 51.6 kg = cohort median LBM in Table 1 (range 31.3-80.3 kg). The paper does not name the LBM formula; Table 1 reports lean body mass (median 51.6 kg) separately from lean body weight (median 54.9 kg), and the model uses the former.",
      source_name = "LBM"
    ),
    DIS_HEALTHY = list(
      description = "Healthy-participant indicator (1 = healthy volunteer, 0 = patient with hATTR-PN)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (patient with hereditary transthyretin amyloidosis with polyneuropathy)",
      notes = "Time-fixed per subject. Healthy stratum = the 51 phase 1 healthy volunteers; reference complement = the 151 hATTR-PN patients of the phase 2/3 study and its open-label extension. Enters as proportional (1 + theta * DIS_HEALTHY) effects on CL/F (+0.111) and Vc/F (-0.284) (Yu 2020 Table 2) and on baseline TTR (+0.269) (Table 3). Table 3 footnote defines the baseline effect as 'the proportional shift in baseline TTR levels of healthy volunteers compared to hATTR patients', which fixes the orientation of the PD term; the PK terms use the same orientation, which reproduces the paper's Table 4 hATTR-patient exposures (Dose / CL = 300 / 3.4 = 88.2 ug*h/mL vs the simulated 89.9) and the Discussion's statement that post hoc CL/F in hATTR patients (3.18 L/h) is close to the typical 3.40 L/h.",
      source_name = "Disease (Table 2 rows 'Disease * CL' and 'Disease * Vc'; Table 3 row 'Disease effect on baseline (healthy volunteers)')"
    )
  )

  covariatesDataExcluded <- list(
    DOSE = list(
      description = "Administered SC dose, normalized to 300 mg",
      units = "mg",
      type = "continuous",
      notes = "Yu 2020 Table 2 row 'Dose effect * CL' (covariate model 'Exponential', estimate 3.97, 95% CI 2.01-5.93) captures a higher apparent clearance at the lowest phase 1 dose (50 mg). The Methods describe it only as 'a simplified exponential model with dose administered being normalized to 300 mg dose'; the literal reading exp(3.97 * (1 - DOSE / 300)) would raise CL/F 7.3-fold at 150 mg, which contradicts the paper's own Table 4 simulations (exposure exactly proportional to dose between 150 and 300 mg). The exact functional form cannot be recovered from the publication, so the term is omitted; the paper states the nonlinearity 'will not be clinically relevant because the projected clinical dose levels are in the range of 150-300 mg'. Not appropriate for doses below 150 mg."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 202L,
    n_subjects_pk = 197L,
    n_observations_pk = 3602L,
    n_studies = 3L,
    n_healthy = 51L,
    n_patients_hattr = 151L,
    age_range = "25 - 81 years",
    age_median = "57 years",
    weight_range = "37.0 - 140 kg",
    weight_median = "71.8 kg",
    lbm_range = "31.3 - 80.3 kg",
    lbm_median = "51.6 kg",
    sex_female_pct = 30.2,
    race_ethnicity = c(Caucasian = 87.6, `African American` = 7.43, Asian = 3.47),
    disease_state = "Healthy volunteers (phase 1) pooled with patients with hereditary transthyretin amyloidosis with polyneuropathy (hATTR-PN; phase 2/3 NEURO-TTR and its open-label extension).",
    dose_range = "Phase 1: SC inotersen 50, 100, 200, or 400 mg single dose, or 50, 100, 200, 300, or 400 mg on days 1, 3, 5, 8, 15, and 22. Phase 2/3: 300 mg SC on days 1, 3, and 5 then once weekly for 65 weeks. OLE: 300 mg SC once weekly for up to 260 weeks.",
    baseline_ttr = "22.0 (5.95) mg/dL mean (SD), range 3.80-39.7 mg/dL (Table 1)",
    notes = "Yu 2020 Table 1 and Clinical studies section. Five hATTR-PN patients had no PK data at the time of analysis, so the PK model used 197 subjects. PK samples collected after the onset of antidrug antibodies were excluded. Plasma inotersen was quantified by hybridization ELISA (LLOQ 1.00 ng/mL); serum TTR by immunoturbidimetry (LLOQ 3 mg/dL). Sequential PK then PD fit (post hoc PK drove the PD model) in NONMEM 7.2, FOCE-I."
  )

  ini({
    # ---- Structural PK parameters (Yu 2020 Table 2) ----
    # Reference subject: hATTR-PN patient (DIS_HEALTHY = 0) with lean body mass
    # 51.6 kg (the cohort median, Table 1). Parameters are apparent (CL/F etc.)
    # because dosing is SC only.
    lka <- log(0.261); label("First-order SC absorption rate constant (1/h)")                     # Yu 2020 Table 2 'Absorption rate constant, Ka (1/h)' = 0.261
    lcl <- log(3.4);   label("Apparent clearance CL/F at reference LBM, hATTR patient (L/h)")      # Yu 2020 Table 2 'Clearance, CL (L/h)' = 3.4
    lvc <- log(20.7);  label("Apparent central volume Vc/F at reference LBM, hATTR patient (L)")   # Yu 2020 Table 2 'Central volume, Vc (L)' = 20.7
    lq  <- log(0.266); label("Apparent intercompartmental clearance Q/F at reference LBM (L/h)")   # Yu 2020 Table 2 'Intercompartmental clearance, Q (L/h)' = 0.266
    lvp <- log(230);   label("Apparent peripheral volume Vp/F at reference LBM (L)")               # Yu 2020 Table 2 'Peripheral volume, Vp (L)' = 230

    # ---- Covariate effects on PK (Yu 2020 Table 2) ----
    e_lbm_cl <- fixed(1); label("Power exponent of LBM on CL/F (unitless)")  # Yu 2020 Table 2 'Lean body mass * CL/Q' Power-centered on median = 1 FIXED
    e_lbm_q  <- fixed(1); label("Power exponent of LBM on Q/F (unitless)")   # Yu 2020 Table 2 'Lean body mass * CL/Q' Power-centered on median = 1 FIXED
    e_lbm_vc <- fixed(1); label("Power exponent of LBM on Vc/F (unitless)")  # Yu 2020 Table 2 'Lean body mass * Vc/Vp' Power-centered on median = 1 FIXED
    e_lbm_vp <- fixed(1); label("Power exponent of LBM on Vp/F (unitless)")  # Yu 2020 Table 2 'Lean body mass * Vc/Vp' Power-centered on median = 1 FIXED
    e_dis_healthy_cl <- 0.111;  label("Proportional effect of healthy status on CL/F (fraction)") # Yu 2020 Table 2 'Disease * CL' Proportional = 0.111
    e_dis_healthy_vc <- -0.284; label("Proportional effect of healthy status on Vc/F (fraction)") # Yu 2020 Table 2 'Disease * Vc' Proportional = -0.284

    # ---- Structural PD parameters (Yu 2020 Table 3) ----
    # Indirect response with inhibition of TTR production (Methods, PD equation):
    #   dTTR/dt = kin * (1 - Imax * Cp / (IC50 + Cp)) - kout * TTR
    lrbase <- log(20.4);   label("Baseline serum TTR in hATTR patients (mg/dL)")                          # Yu 2020 Table 3 'Estimated baseline' = 20.4 mg/dL
    e_dis_healthy_rbase <- 0.269; label("Proportional effect of healthy status on baseline TTR (fraction)") # Yu 2020 Table 3 'Disease effect on baseline (healthy volunteers)' = 0.269
    lkout  <- log(0.00308); label("First-order TTR loss rate constant (1/h)")                              # Yu 2020 Table 3 'kout' = 0.00308 1/h
    imax   <- 0.913;        label("Maximum fractional inhibition of TTR production (unitless)")           # Yu 2020 Table 3 'Imax' = 0.913
    lic50  <- log(9.07);    label("Plasma inotersen concentration giving half-maximal inhibition (ng/mL)") # Yu 2020 Table 3 'IC50' = 9.07 ng/mL

    # ---- IIV on PK (Yu 2020 Table 2; omega^2 on the log scale) ----
    # Order etalcl, etalvc, etalvp. Table 2 reports covariances Vc*CL and CL*Vp
    # only; the Vc*Vp covariance is not reported and is set to 0.
    etalcl + etalvc + etalvp ~ c(
      0.071,
      0.049, 0.335,
      -0.145, 0, 0.677
    ) # Yu 2020 Table 2: omega2 CL = 0.071, Covariance Vc*CL = 0.049, omega2 Vc = 0.335, Covariance CL*Vp = -0.145, omega2 Vp = 0.677

    # ---- IIV on PD (Yu 2020 Table 3) ----
    etalic50 + etalkout ~ c(
      0.953,
      -0.182, 0.288
    ) # Yu 2020 Table 3: omega2 IC50 = 0.953, Covariance IC50*kout = -0.182, omega2 kout = 0.288
    etalrbase ~ 0.0537 # Yu 2020 Table 3 'omega2 Baseline' = 0.0537

    # ---- Residual error ----
    # PK: Table 2 'Residual variability, sigma2' = 0.168 with 'Log additive
    # error' on log-transformed concentrations -> lnorm SD = sqrt(0.168).
    expSd <- 0.4099; label("Log-scale residual SD on plasma inotersen Cc") # Yu 2020 Table 2 sigma2 = 0.168 (log additive) -> sqrt(0.168) = 0.4099
    # PD: Table 3 'Residual variability (log space), (SD)' reports a
    # proportional and an additive component for log-transformed TTR; encoded
    # as combined proportional + additive error on linear-scale TTR.
    propSd_ttr <- 0.13; label("Proportional residual SD on serum TTR (fraction)") # Yu 2020 Table 3 'Proportional error' = 0.13
    addSd_ttr  <- 1.35; label("Additive residual SD on serum TTR (mg/dL)")       # Yu 2020 Table 3 'Additive error' = 1.35
  })

  model({
    # ---- 1. Individual PK parameters ----
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (LBM / 51.6)^e_lbm_cl * (1 + e_dis_healthy_cl * DIS_HEALTHY)
    vc <- exp(lvc + etalvc) * (LBM / 51.6)^e_lbm_vc * (1 + e_dis_healthy_vc * DIS_HEALTHY)
    q  <- exp(lq) * (LBM / 51.6)^e_lbm_q
    vp <- exp(lvp + etalvp) * (LBM / 51.6)^e_lbm_vp

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ---- 2. Individual PD parameters ----
    rbase <- exp(lrbase + etalrbase) * (1 + e_dis_healthy_rbase * DIS_HEALTHY)
    kout  <- exp(lkout + etalkout)
    ic50  <- exp(lic50 + etalic50)
    kin   <- rbase * kout

    # ---- 3. PK ODEs (amounts in mg; volumes in L) ----
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                                k12 * central - k21 * peripheral1

    # central / vc is mg/L = ug/mL; x 1000 gives ng/mL, the units of IC50.
    Cc <- central / vc * 1000

    # ---- 4. PD: indirect response, inhibition of TTR production ----
    inh          <- 1 - imax * Cc / (ic50 + Cc)
    d/dt(effect) <- kin * inh - kout * effect
    effect(0)    <- rbase

    # Serum TTR (mg/dL), the paper-named alias of the `effect` state.
    ttr <- effect

    # ---- 5. Observations ----
    Cc  ~ lnorm(expSd)
    ttr ~ add(addSd_ttr) + prop(propSd_ttr)
  })
}
