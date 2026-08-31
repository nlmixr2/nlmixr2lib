Dohmann_2025_piperacillin <- function() {
  description <- "One-compartment population PK model for intravenous piperacillin (administered as piperacillin/tazobactam) in adult non-intensive-care ward patients with chronic kidney disease, including patients on thrice-weekly intermittent haemodialysis (n=49; 135 total-piperacillin serum samples; Dohmann 2025). Elimination is linear and is decomposed into a non-dialysis clearance carrying a power effect of BSA-individualized MDRD eGFR, plus an additive haemodialysis-arm clearance gated on/off by the time-varying RRT_HEMODIAL_ACTIVE regressor. Central volume carries a power effect of body surface area. Residual error is proportional; there is no interindividual variability on the haemodialysis clearance. The model was built in Monolix 2023R1 and was used for Monte Carlo probability-of-target-attainment simulations against a conservative (60% fT > MIC) and an aggressive (100% fT > 4 x MIC) PK/PD target for Pseudomonas aeruginosa."
  reference <- "Dohmann E, Hagel S, Kurlbaum M, Schellong P, Scherf-Clavel O, Surat G. Probability of pharmacokinetic/pharmacodynamic target attainment for different piperacillin/tazobactam dosing regimens in renally impaired patients in a non-intensive care unit setting. Br J Clin Pharmacol. 2025;91(11):3070-3081. doi:10.1002/bcp.70153"
  vignette <- "Dohmann_2025_piperacillin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Dohmann 2025 Methods 2.4 and 2.6:
  # piperacillin was quantified in serum and the model simulates TOTAL (not
  # unbound) piperacillin concentrations.
  compartmentData <- list(
    central = list(analyte = "piperacillin", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    BSA = list(
      description        = "Body surface area (Mosteller formula)",
      units              = "m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Dohmann 2025 Methods 2.5: BSA computed with the Mosteller formula and retained as the only covariate on Vd, entering as the power term (BSA / 1.88 m^2)^beta_BSA. IMPORTANT -- the reference value encoded here (1.88 m^2) is NOT the value printed in the paper's typeset Vd equation on p.3074, which reads 1.66 m^2. The printed 1.66 m^2 is treated as a typographical error and 1.88 m^2 is encoded instead, per operator decision (sidecar request-001 / response-001, question q1, answer A). Five independent lines of evidence inside the paper itself support 1.88 m^2: (i) Methods 2.6 states the simulation cohort drew BSA from a lognormal distribution with MEAN 1.88 m^2 and SD 0.12 m^2, i.e. the authors' own cohort-typical BSA; (ii) Table 1 mean height 170.9 cm and mean weight 76.1 kg give a Mosteller BSA of sqrt(170.9 x 76.1 / 3600) = 1.90 m^2, so a cohort median of 1.66 m^2 is not attainable; (iii) the Discussion states 'the population estimate for volume of distribution of the herein presented model is 15.4 L for an adult patient', which only holds if the 15.42 L typical value corresponds to a cohort-typical rather than a small adult; (iv) the companion eGFR reference in the SAME pair of equations (21.8 mL/min) is demonstrably the cohort central value on the individualized scale (Table 1 mean 19.9 mL/min/1.73 m^2 x 1.90 / 1.73 = 21.85 mL/min), so the same convention applied to BSA gives 1.88-1.90 m^2; (v) reimplementing the paper's own Monte Carlo simulation reproduces the six mid-range cells of Table 3 to within about 1 percentage point with 1.88 m^2 but over-predicts every cell with 1.66 m^2. See the vignette 'Assumptions, deviations and errata' section, which reports both readings side by side.",
      source_name        = "BSA"
    ),
    CRCL = list(
      description        = "Estimated glomerular filtration rate, MDRD formula individualized using body surface area (i.e. de-normalized to absolute mL/min, NOT mL/min/1.73 m^2)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Dohmann 2025 Methods 2.5 lists the available covariate as 'eGFR (Modification of Diet in Renal Disease formula individualized using BSA)', and the typeset CL equation on p.3074 divides by '21.8 mL/min' -- absolute mL/min units, not mL/min/1.73 m^2. This is the raw un-normalized variant of the CRCL canonical (same convention as Georges_2009_ceftazidime.R and Delattre_2010_amikacin.R); supplying a BSA-normalized value here would silently rescale the renal term. Consistency check on the reference value: Table 1 reports a cohort mean eGFR of 19.9 mL/min/1.73 m^2, which individualized to the cohort-typical BSA gives 19.9 x 1.90 / 1.73 = 21.85 mL/min, matching the printed 21.8 mL/min. Study eGFR was computed with MDRD at Wuerzburg and CKD-EPI at Jena (Table 1 note). Enters as the power term (CRCL / 21.8)^beta_eGFR on the non-dialysis clearance only; the haemodialysis arm is renal-function independent. In the paper's own Monte Carlo simulations (Methods 2.6) eGFR was 'set exactly to the value according to the group', i.e. the group labels 40 / 30 / 20 / 10 mL/min are covariate values on this absolute scale.",
      source_name        = "eGFR"
    ),
    RRT_HEMODIAL_ACTIVE = list(
      description        = "Haemodialysis-active indicator (1 while an intermittent haemodialysis session is running, 0 otherwise)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (interdialytic / no dialysis running)",
      notes              = "Time-varying within subject. Dohmann 2025 Methods 2.5: 'Patients with iHD were modelled by including a second clearance process during the times of dialysis. The dialysis clearance was turned on or off by a regressor variable (theta_HD = 0 or 1) according to the documented dialysis times.' The paper's own symbol is theta_HD; it is a data-supplied 0/1 regressor, not an estimated parameter, so it is registered here as a covariate column. Group 3 patients received thrice-weekly intermittent haemodialysis with a mean session length of 258.8 +/- 24.8 min (Table 1); the Monte Carlo simulations (Methods 2.6) used a 4-h session once daily beginning directly after the end of the piperacillin infusion. Gates the additive cl_hemodialysis arm exactly as in Veinstein_2013_gentamicin.R and VanWart_2025_telavancin.R.",
      source_name        = "theta_HD"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age", units = "years", type = "continuous",
      notes = "Dohmann 2025 Methods 2.5 lists age among the available covariates for model building; Results 3.2 reports that 'in the stepwise testing of covariates, only BSA and eGFR were significant covariates', so age was screened and not retained. No point estimate is reported (see Table S1 / Figure S4, not on disk)."
    ),
    BMI = list(
      description = "Body mass index", units = "kg/m^2", type = "continuous",
      notes = "Dohmann 2025 Methods 2.5: BMI was among the covariates considered on Vd; screened and not retained in the final model (Results 3.2)."
    ),
    HT = list(
      description = "Height", units = "cm", type = "continuous",
      notes = "Dohmann 2025 Methods 2.5: height was among the covariates considered on Vd; screened and not retained in the final model (Results 3.2). Height is nevertheless an input to the Mosteller BSA that IS retained."
    ),
    WT = list(
      description = "Body weight", units = "kg", type = "continuous",
      notes = "Dohmann 2025 Methods 2.5: weight was among the covariates considered on Vd; screened and not retained in the final model (Results 3.2). Weight is nevertheless an input to the Mosteller BSA that IS retained."
    ),
    SEXF = list(
      description = "Female sex indicator (1 = female, 0 = male)", units = "(binary)", type = "binary",
      notes = "Dohmann 2025 Methods 2.5 gives the categorical-covariate form explicitly -- Vd_ind = Vd_pop * exp(beta_sex * theta_sex) * exp(eta_Vd) with theta_sex = 1 for female and 0 for male -- and tested sex on both CL and Vd. Results 3.2: only BSA and eGFR were significant, so no sex effect is present in the final model and no beta_sex estimate is reported."
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 49L,
    n_studies        = 1L,
    age_range        = "35-88 years (mean 69.3 +/- 11.1; Table 1 and Discussion, strengths (i))",
    weight_range     = "mean 76.1 +/- 17.1 kg; BMI mean 26.1 +/- 5.6 kg/m^2 spanning <18.5 to >40 kg/m^2 (Table 1 and Discussion, strengths (i))",
    sex_female_pct   = 36.7,
    race_ethnicity   = "Not reported (two German university hospitals)",
    disease_state    = "Adults on general (non-intensive-care) wards with chronic kidney disease receiving intermittent intravenous piperacillin/tazobactam. Three prespecified renal-function groups: group 1 eGFR 20-40 mL/min/1.73 m^2 (n=20), group 2 eGFR <20 mL/min/1.73 m^2 (n=19), group 3 thrice-weekly intermittent haemodialysis (n=10). Indications included empirical therapy for nosocomial infection of unclear focus, and targeted therapy for hospital-acquired pneumonia or severe skin infection with or without bacteraemia (Methods 2.3). Exclusion criteria: pregnancy or breastfeeding, age under 18 years, beta-lactam hypersensitivity, participation in another clinical trial.",
    dose_range       = "Piperacillin/tazobactam 4.5 g (4000 mg piperacillin / 500 mg tazobactam) as a 30-min intermittent intravenous infusion. Group 1 q8h, groups 2 and 3 q12h, per the Summary of Product Characteristics; group 3 additionally received 2.25 g after haemodialysis (Table 1). A minority of patients deviated from the recommended interval (six patients in group 1 received q12h, three in group 2 received q8h, one in group 3 received q24h). Only piperacillin was quantified.",
    regions          = "Germany (Wuerzburg University Hospital n=30, Jena University Hospital n=19)",
    renal_function   = "Chronic kidney disease with eGFR <40 mL/min/1.73 m^2, or thrice-weekly intermittent haemodialysis. Cohort mean eGFR 19.9 +/- 9.3 mL/min/1.73 m^2 (group 1 29.2 +/- 6.6, group 2 14.6 +/- 3.6, group 3 11.2 +/- 3.3). eGFR was calculated with MDRD at Wuerzburg and CKD-EPI at Jena; creatinine clearance was not measured because it is not routinely determined in non-ICU patients (Discussion, limitations (i)).",
    notes            = "Prospective multicentre study conducted 2 March 2022 to 17 March 2023 (ethics approvals 110/21 and 2021-2399). Sampling at steady state after 3 days of therapy at three time points for groups 1 and 2 (within 30 min pre-dose, within 30 min after end of infusion, and 4 or 6 h post-infusion) and immediately before and after a dialysis session for group 3. A total of 135 samples were analysed (80 from 30 Wuerzburg patients and 55 from 19 Jena patients); one sample at each site was excluded for an implausibly high concentration. Bioanalysis by HPLC-MS/MS (Sciex QTRAP 4500MD with Agilent 1290 UHPLC), analytical range 0.5-190.0 mg/L, validated to EMA guidance. Concentrations are reported and modelled as TOTAL piperacillin; the paper converts to unbound by multiplying by 0.8 (an 80% unbound fraction taken from the literature), so the free target of 16 mg/L corresponds to a total concentration of 19.8-20 mg/L (Methods 2.4 and 2.6). Model estimated in Monolix 2023R1; simulations in Simulx 2023R1. Supplementary Table S1 and Figures S1-S4 (covariate screening and goodness-of-fit diagnostics) are not on disk; they contain no final parameter values, all of which are in Table 2."
  )

  ini({
    # ---- Structural parameters (Dohmann 2025 Table 2, final model) --------
    # The typeset covariate equations are given on p.3074 (Methods 2.5):
    #   CL_ind = CL_pop * (eGFR / 21.8 mL/min)^beta_eGFR * exp(eta_CL)
    #            + (CL_HD * theta_HD)
    #   Vd_ind = Vd_pop * (BSA  / 1.66 m^2  )^beta_BSA  * exp(eta_Vd)
    # NOTE on the BSA reference: 1.88 m^2 is encoded in model() below rather
    # than the printed 1.66 m^2. See covariateData$BSA$notes and the vignette
    # errata section for the full evidence and the operator decision.
    lcl              <- log(3.52);  label("Non-dialysis clearance CL_pop at eGFR 21.8 mL/min (L/h)")  # Table 2: CL_pop = 3.52 L/h (SE 0.2, RSE 5.68%)
    lvc              <- log(15.42); label("Central volume of distribution Vd_pop at BSA 1.88 m^2 (L)") # Table 2: Vd_pop = 15.42 L (SE 1.08, RSE 7.02%)
    lcl_hemodialysis <- log(3.96);  label("Additive haemodialysis-arm clearance CL_HD (L/h)")          # Table 2: CL_HD = 3.96 L/h (SE 1.36, RSE 34.4%)

    # ---- Covariate effects (Dohmann 2025 Table 2) ------------------------
    e_crcl_cl <- 0.54; label("Power exponent of (CRCL / 21.8 mL/min) on non-dialysis CL (unitless)")  # Table 2: beta_eGFR = 0.54 (SE 0.11, RSE 20.0%)
    e_bsa_vc  <- 1.46; label("Power exponent of (BSA / 1.88 m^2) on Vd (unitless)")                   # Table 2: beta_BSA  = 1.46 (SE 0.42, RSE 28.7%)

    # ---- Interindividual variability (Dohmann 2025 Table 2) --------------
    # Monolix lognormal random effects; the Table 2 note defines omega_CL and
    # omega_Vd as the "standard deviation of random effects", so the nlmixr2
    # ini() VARIANCES are omega^2:
    #   omega_CL = 0.34 -> 0.34^2 = 0.1156
    #   omega_Vd = 0.24 -> 0.24^2 = 0.0576
    # Cross-check that these are SDs and not variances: for a variance
    # estimated from n = 49 subjects the RSE cannot fall much below
    # sqrt(2/49) = 20.2%, yet omega_CL is reported with an RSE of 12.4%, which
    # is only attainable on the SD scale.
    etalcl ~ 0.1156  # Table 2: omega_CL = 0.34 (SE 0.042, RSE 12.4%)
    etalvc ~ 0.0576  # Table 2: omega_Vd = 0.24 (SE 0.056, RSE 23.2%)

    # No interindividual variability on the haemodialysis clearance. Results
    # 3.2: "Interindividual variability on haemodialysis clearance was not
    # included in the model, since the estimate was unreliable due to high
    # standard error, which was also probably because of the small patient
    # number."

    # ---- Residual error --------------------------------------------------
    # Results 3.2: "Different residual error models were evaluated, and the
    # proportional error model was chosen, as it led to no obvious pattern in
    # the residuals".
    propSd <- 0.26; label("Proportional residual error (fraction)")  # Table 2: proportional error = 0.26 (SE 0.024, RSE 9.38%)
  })

  model({
    # Non-dialysis (intrinsic body) clearance with the power effect of the
    # BSA-individualized MDRD eGFR, reference 21.8 mL/min as printed in the
    # Methods 2.5 CL equation.
    cl <- exp(lcl + etalcl) * (CRCL / 21.8)^e_crcl_cl  # L/h

    # Central volume with the power effect of body surface area. The
    # reference divisor is 1.88 m^2 (the authors' own simulation-cohort mean
    # BSA, Methods 2.6) rather than the 1.66 m^2 printed in the Methods 2.5
    # Vd equation; see covariateData$BSA$notes for the evidence and the
    # vignette errata section for the side-by-side comparison.
    vc <- exp(lvc + etalvc) * (BSA / 1.88)^e_bsa_vc    # L

    # Haemodialysis-arm clearance: an additive term with no interindividual
    # variability, gated on only while a session is running.
    cl_hemodialysis <- exp(lcl_hemodialysis)           # L/h

    # Methods 2.5: "Patients with iHD were modelled by including a second
    # clearance process during the times of dialysis. The dialysis clearance
    # was turned on or off by a regressor variable (theta_HD = 0 or 1)".
    cl_total <- cl + RRT_HEMODIAL_ACTIVE * cl_hemodialysis  # L/h

    kel <- cl_total / vc

    # One-compartment model with linear elimination and intravenous
    # (zero-order infusion) input; no absorption compartment.
    d/dt(central) <- -kel * central

    # Dose in mg of piperacillin, vc in L -> mg/L. Concentrations are TOTAL
    # piperacillin (Methods 2.4); multiply by 0.8 to obtain unbound
    # concentrations for PK/PD target evaluation (Methods 2.6).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
