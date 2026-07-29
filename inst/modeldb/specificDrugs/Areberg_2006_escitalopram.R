Areberg_2006_escitalopram <- function() {
  description <- "Two-compartment population PK model with first-order absorption and lag time for escitalopram in healthy and hepatic-impaired adults (Areberg 2006)"
  reference <- "Areberg J, Christophersen JS, Poulsen MN, Larsen F, Molz K-H. The Pharmacokinetics of Escitalopram in Patients With Hepatic Impairment. AAPS J. 2006;8(1):E14-E19 (Article 2). doi:10.1208/aapsj080102"
  vignette <- "Areberg_2006_escitalopram"
  units <- list(time = "hour", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject in the single-dose Areberg 2006 study. Linear-additive deviation effect on apparent central volume V/F centered on the pooled cohort mean (~79 kg from Table 1).",
      source_name        = "WT"
    ),
    CYP2C19 = list(
      description        = "CYP2C19 metabolic-activity proxy (urinary S/R-mephenytoin ratio after a single 100 mg dose of racemic mephenytoin; Methods 'CYP2C19 Phenotyping')",
      units              = "(unitless ratio)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per subject (one-time probe-substrate assay). Higher S/R ratio indicates LOWER CYP2C19 activity because high CYP2C19 activity selectively metabolizes S-mephenytoin and lowers the residual S enantiomer in urine; the model uses a linear-additive deviation form centered on the population mean ratio of 0.769, and a NEGATIVE coefficient on the (CYP2C19 - 0.769) deviation so that subjects with higher S/R ratio have lower apparent escitalopram clearance. This orientation is the OPPOSITE of the dextromethorphan-probe CYP2D6 / CYP3A4 columns (where higher value = higher activity); see covariate-columns.md CYP2C19 entry Notes.",
      source_name        = "CYP2C19"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 24,
    n_studies       = 1,
    age_range       = "43-69 years",
    age_median      = "58 years (group-mean range 57.6-59.0; Table 1)",
    weight_range    = "63-101 kg",
    weight_median   = "~79 kg (pooled-cohort mean; group-mean range 76.3-83.9; Table 1)",
    sex_female_pct  = 25,
    race_ethnicity  = c(White = 100),
    disease_state   = "Stratified Child-Pugh classification: 8 with normal hepatic function, 8 with mild hepatic impairment (Child-Pugh 5-6), 8 with moderate hepatic impairment (Child-Pugh 7-9; alcohol-induced cirrhosis).",
    dose_range      = "20 mg single oral dose",
    regions         = "Single site, Munich, Germany",
    notes           = "All Caucasian; no clinically meaningful differences between groups in age, body weight, or creatinine clearance (all CrCL > 80 mL/min). Concomitant CYP2C19 / CYP2D6 inhibitors and inducers excluded for 6 weeks prior to dosing. Baseline demographics per Areberg 2006 Table 1."
  )

  ini({
    # Structural parameters: typical values for the 'average' subject
    # (WT = 79 kg pooled mean and CYP2C19 = 0.769 S/R-mephenytoin ratio).
    # 'The "average" subject would have values for V/F and CL/F of 8.0*10^2 L
    # and 19 L/h, respectively.' (Areberg 2006 Results, p E17).
    lka      <- fixed(log(3.6)); label("Absorption rate constant (Ka, 1/h; fixed)")                      # Areberg 2006 Table 3 (theta_5, 'fixed'; Results: 'k_a had to be fixed in order to run the model successfully ... set to 3.6 h^-1 ... based on modeling the data from each subject separately')
    lcl      <- log(19);         label("Apparent oral clearance at the cohort-mean CYP2C19 ratio (CL/F, L/h)")  # Areberg 2006 Table 3 (theta_2, 19 L/h, RSE 7%) / Results 'The average subject would have values for V/F and CL/F of ... 19 L/h'
    lvc      <- log(800);        label("Apparent central volume of distribution at the cohort-mean weight (V/F, L)")  # Areberg 2006 Table 3 (theta_1, 8.0*10^2 L, RSE 6%) / Results 'The average subject would have values for V/F ... of 8.0*10^2 L'
    lvp      <- log(470);        label("Apparent peripheral volume of distribution (V2/F, L)")            # Areberg 2006 Table 3 (theta_3, 4.7*10^2 L, RSE 8%)
    lq       <- log(98);         label("Apparent inter-compartmental clearance (CLD2/F, L/h)")            # Areberg 2006 Table 3 (theta_4, 98 L/h, RSE 13%)
    ltlag    <- log(0.93);       label("Absorption lag time (t_lag, h)")                                  # Areberg 2006 Table 3 (theta_6, 0.93 h, RSE 1%)

    # Covariate-effect coefficients (linear-additive deviation form;
    # 'Mean normalized (centered) covariate values were used.', Methods p E16).
    e_wt_vc       <-  12;  label("Weight effect on V/F (linear-additive coefficient, L per kg of (WT - 79))")     # Areberg 2006 Table 3 (theta_7, 12 L/kg, RSE 26%) / Eq. 1
    e_cyp2c19_cl <- -17;  label("CYP2C19 effect on CL/F (linear-additive coefficient, L/h per unit of (CYP2C19 - 0.769); negative because higher S/R ratio = lower CYP2C19 activity)")  # Areberg 2006 Table 3 (theta_8 magnitude 17 L/h, RSE 23%) / Eq. 2; sign biologically inferred -- higher S/R-mephenytoin ratio indicates LOWER CYP2C19 activity, so a positive deviation in CYP2C19 lowers escitalopram clearance.

    # Inter-individual variability (omega^2 = log(CV^2 + 1)).
    # Final-model CV%: V/F = 17, CL/F = 34, ka = 148; diagonal covariance
    # (Methods 'A diagonal covariance matrix was used; that is, covariances
    # between the structural model parameters were assumed to be negligible.').
    etalvc ~ 0.02849   # log(0.17^2 + 1)
    etalcl ~ 0.10940   # log(0.34^2 + 1)
    etalka ~ 1.16025   # log(1.48^2 + 1); 31% RSE on the omega per Table 3 -- the typical-value Ka is fixed at 3.6 1/h but eta on Ka is estimated, per Discussion 'k_a values between subjects were allowed to vary because of interindividual variability'.

    # Residual error (proportional).
    propSd <- 0.096; label("Proportional residual error (fraction)")  # Areberg 2006 Table 3 (epsilon_1, 9.6%, RSE 23%) / Results 'A proportional model was used for the residual error.'
  })
  model({
    # Linear-additive covariate parameterisation, Eq. 1 and Eq. 2:
    #   V/F  = 800 + 12 * (WT     - 79)
    #   CL/F = 19  - 17 * (CYP2C19 - 0.769)
    # Multiplicative exponential IIV applied to the covariate-adjusted typical value.
    ka   <- exp(lka + etalka)
    vc   <- (exp(lvc) + e_wt_vc       * (WT      - 79))    * exp(etalvc)
    cl   <- (exp(lcl) + e_cyp2c19_cl * (CYP2C19 - 0.769)) * exp(etalcl)
    vp   <- exp(lvp)
    q    <- exp(lq)
    tlag <- exp(ltlag)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    alag(depot) <- tlag

    # Observation: dose in mg, volume in L -> mg/L == ug/mL.
    # The source paper reports serum concentrations in nmol/L. The validation
    # vignette converts via the escitalopram free-base molecular weight
    # 324.39 g/mol (Cc [nmol/L] = Cc [ug/mL] * 1e6 / 324.39, since
    # 1 ug/mL == 1 mg/L == 1e6/MW nmol/L).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
