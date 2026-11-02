Samuels_2025_infliximab <- function() {
  description <- "Two-compartment population PK model of intravenous infliximab in children and young adults with Crohn's disease, extended with dose-number-driven trajectories of the four clearance biomarkers to forecast exposure during induction and early maintenance (Samuels 2025). The structural model, covariate exponents and variability are the Xiong 2021 real-world paediatric model (REFINE cohort), all held at the Xiong estimates: clearance scales with body weight, serum albumin, erythrocyte sedimentation rate (ESR), neutrophil CD64 activity ratio and antibody-to-infliximab concentration; V1 with weight; Q and V2 with weight at exponents 1 and 0.586. Samuels 2025 adds an Emax-in-dose-number function that predicts each biomarker from its pre-treatment value as baseline * exp(Emax * (n - 1) / (D50 - 1 + n - 1)) at infusion n, with D50 held at 1.1 and Emax estimated for weight (+0.046), albumin (+0.110), ESR (-0.978) and nCD64 (-0.362), so that the model anticipates the improvement in clearance biomarkers that accompanies a response to induction. Holding the infusion counter at 1 throughout reproduces the baseline-covariates-only Xiong prediction. Inter-individual variability on clearance only; residual error additive on log-transformed concentrations."
  reference <- paste(
    "Samuels A, Irie K, Mizuno T, Reifenberg J, Punt N, Vinks AA, Minar P.",
    "Integrating early response biomarkers in pharmacokinetic models: A novel",
    "method to individualize the initial infliximab dose in patients with",
    "Crohn's disease. Clin Transl Sci. 2025;18(2):e70086.",
    "doi:10.1111/cts.70086. Structural model and parameter values from: Xiong Y,",
    "Mizuno T, Colman R, Hyams J, Noe JD, Boyle B, et al. Real-world infliximab",
    "pharmacokinetic study informs an electronic health record-embedded",
    "dashboard to guide precision dosing in children with Crohn's disease.",
    "Clin Pharmacol Ther. 2021;109(6):1639-1647. doi:10.1002/cpt.2148.",
    sep = " "
  )
  vignette <- "Samuels_2025_infliximab"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  # What each ODE state holds. Infliximab is given as an intravenous infusion
  # into the central compartment (Data S1 code: ADVAN3 TRANS4, S1 = V1); serum
  # infliximab was measured by the Esoterix drug-tolerant ECLIA (Methods,
  # 'Trough concentrations').
  compartmentData <- list(
    central = list(analyte = "infliximab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "infliximab", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight at the first infusion (pre-treatment baseline)",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Enter the PRE-TREATMENT value and hold it constant: the model predicts the on-treatment weight itself as WT * exp(wt_dosenum_max * f(CYCLE)), so supplying an observed on-treatment weight would apply the predicted gain twice. The predicted weight enters CL as (WT_pred/65)^0.594 and V1 as (WT_pred/65)^0.55. Q and V2 are scaled by the UNPREDICTED baseline weight, (WT/65)^1 and (WT/65)^0.586, exactly as in the Data S1 code (Q=TVQ*EXP(ETA(3))*(WT/65), V2=TVV2*EXP(ETA(4))*(WT/65)**0.586). Reference weight 65 kg (Xiong 2021 Table 2 'CL, L/h/65 kg'). Samuels 2025 Table 1 cohort mean 48.3 (SD 20.5) kg at the first infusion; the Emax fit (Table S1) gives E0 = 48.8 kg. A missing weight in the source code (WT = -99) was handled by assuming the reference value; impute 65 kg to reproduce that.",
      source_name = "WT"
    ),
    ALB = list(
      description = "Serum albumin at the first infusion (pre-treatment baseline)",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Canonical SI g/L; the model converts to the g/dL in which the coefficients were estimated (alb_gdl <- ALB * 0.1). Enter the PRE-TREATMENT value and hold it constant: the on-treatment albumin is predicted as ALB * exp(alb_dosenum_max * f(CYCLE)). Power effect on CL (ALB_pred/3.5 g/dL)^-1.07 (Xiong 2021 Table 2 'ALB for CL' -1.07; Data S1 THETA(5) = -1.07 and (ALB*E1)/3.5). Samuels 2025 Table 1 cohort mean 3.3 (SD 0.61) g/dL; Emax-fit E0 = 3.3 g/dL. A missing albumin (ALB = -99 in the Data S1 code) contributes E1^THETA(5), which is what imputing the 3.5 g/dL reference (35 g/L) reproduces.",
      source_name = "ALB"
    ),
    ESR = list(
      description = "Erythrocyte sedimentation rate at the first infusion (pre-treatment baseline)",
      units = "mm/h",
      type = "continuous",
      reference_category = NULL,
      notes = "Enter the PRE-TREATMENT value and hold it constant: the on-treatment ESR is predicted as ESR * exp(esr_dosenum_max * f(CYCLE)). Power effect on CL (ESR_pred/9)^0.101 (Xiong 2021 Table 2 'ESR for CL' 0.101; Data S1 THETA(8) and (ESR*E3)/9). Samuels 2025 Table 1 cohort mean 26 (SD 24) mm/h; Emax-fit E0 = 27.6 mm/h. A missing ESR in the Data S1 code contributes E3^THETA(8), reproduced by imputing the 9 mm/h reference. The power form is undefined at 0 mm/h; floor a reported 0 at the assay's lowest reportable value.",
      source_name = "ESR"
    ),
    NEUTROPHIL_CD64_RATIO = list(
      description = "Neutrophil CD64 activity ratio (nCD64) at the first infusion (pre-treatment baseline)",
      units = "(ratio)",
      type = "continuous",
      reference_category = NULL,
      notes = "Flow-cytometric ratio of granulocyte to lymphocyte CD64 mean fluorescence intensity (Xiong 2021 Methods, 'Biologic assays'; Xiong 2021 Figure S1 gives the conversion from the older CD64 index, nCD64 = 5.37 * index + 0.9147). Enter the PRE-TREATMENT value and hold it constant: the on-treatment ratio is predicted as NEUTROPHIL_CD64_RATIO * exp(neutrophil_cd64_ratio_dosenum_max * f(CYCLE)). Power effect on CL (nCD64_pred/4.6)^0.168 (Xiong 2021 Table 2 'nCD64 for CL' 0.168; Data S1 THETA(7) and (NCAR_n*E4)/4.6). Samuels 2025 Table 1 cohort mean 7.2 (SD 4.6); Emax-fit E0 = 7.7. A missing value in the Data S1 code contributes E4^THETA(7), reproduced by imputing the 4.6 reference.",
      source_name = "NCAR_n"
    ),
    CONC_ADA_NGML = list(
      description = "Antibody-to-infliximab (ATI) concentration",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Drug-tolerant Esoterix ECLIA (lower limit of detection 22 ng/mL; Samuels 2025 Methods and Xiong 2021 Methods). Thresholded power effect on CL: 1 when ATI <= 22 ng/mL, otherwise (ATI/22)^0.134 (Data S1 code IF(ATI.LE.22) CLATI=1 ELSE CLATI=(ATI/22)**THETA(6); Xiong 2021 Table 2 'ATI for CL' 0.134). Unlike the four biomarkers above, ATI has no predicted trajectory and is used as observed (time-varying). Samuels 2025 Discussion states that ATI was not part of the trajectory extension because antibodies were scarce in both cohorts; for a forward simulation of an antibody-negative patient set it to any value <= 22 (the vignette uses 0).",
      source_name = "ATI"
    ),
    CYCLE = list(
      description = "Infliximab infusion number (1 = first infusion), the driver of the biomarker trajectories",
      units = "(count)",
      type = "count",
      reference_category = NULL,
      notes = "Integer infusion counter, 1 from the first infusion until the second, 2 from the second until the third, and so on (Samuels 2025 Methods equation 'Dose Number'; Figure 1 x-axis). The Data S1 code carries the zero-based count INFUSION = dose number - 1, which is what the D50 - 1 = 0.1 constant in E1-E4 = EXP(Emax*INFUSION/(0.1+INFUSION)) implies; the model forms that count as CYCLE - 1. At CYCLE = 1 every biomarker equals its baseline, so holding CYCLE at 1 for the whole simulation gives the 'baseline covariates alone' prediction of Samuels 2025 Figure 2 and Table 3 scenario a. Records before the first infusion should carry CYCLE = 1; values below 1 are treated as 1.",
      source_name = "INFUSION (= CYCLE - 1)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 239L,
    n_studies = 2L,
    age_range = "2-22 years (APPDASH inclusion criterion; REFINE enrolled children and young adults under 22)",
    age_mean = "13.6 (SD 3.7) years at the first infusion; 92.1% under 18 years",
    weight_mean = "48.3 (SD 20.5) kg at the first infusion",
    sex_female_pct = 40.6,
    race_ethnicity = "86.6% White",
    disease_state = "Crohn's disease starting infliximab (anti-TNF naive in REFINE); Paris location L1:L2:L3 = 31:46:162, behaviour B1:B2:B3:B4 = 186:26:18:9, perianal disease 12.6%",
    dose_range = "intravenous infliximab, mean starting dose 7.7 (SD 2.3) mg/kg (REFINE 6.5, APPDASH 8.2 mg/kg); induction at weeks 0, 2 and 6, first maintenance infusion near week 14",
    regions = "United States (Cincinnati Children's Hospital Medical Center, Nationwide Children's Hospital, Medical College of Wisconsin, Connecticut Children's Medical Center)",
    notes = "Two cohorts: REFINE (prospective multicentre, n = 78, 2014-2019) and APPDASH (retrospective single-centre, n = 161, 2019-2021) (Samuels 2025 Methods and Table 1). The biomarker-trajectory Emax parameters were estimated on the combined 239 patients. The PK parameters are not re-estimated: they are the Xiong 2021 final model, developed on 671 peak and trough concentrations from the 78 REFINE patients, and are held at those values in the Samuels 2025 Data S1 code. Baseline biomarker means (Table 1): albumin 3.3 g/dL, ESR 26 mm/h, nCD64 7.2."
  )

  ini({
    # Structural parameters, Xiong 2021 final model, held at those estimates by
    # Samuels 2025 (Data S1 $THETA, every entry marked FIX). Reference subject
    # is 65 kg, albumin 3.5 g/dL, ESR 9 mm/h, nCD64 4.6, ATI <= 22 ng/mL.
    lcl <- fixed(log(0.0138)); label("Clearance for the reference subject (L/h)")                           # Data S1 THETA(1) = 0.0138; Xiong 2021 Table 2 'CL, L/h/65 kg' 0.0138 (RSE 5.7%)
    lvc <- fixed(log(2.97)); label("Central volume of distribution for a 65 kg subject (L)")               # Data S1 THETA(2) = 2.97; Xiong 2021 Table 2 'V1, L/65kg' 2.97 (RSE 4.5%)
    lq <- fixed(log(0.00953)); label("Inter-compartmental clearance for a 65 kg subject (L/h)")            # Data S1 THETA(3) = 0.00953; Xiong 2021 Table 2 'Q, L/h/65 kg' 0.0095, itself carried from Fasanmade 2011
    lvp <- fixed(log(2.84)); label("Peripheral volume of distribution for a 65 kg subject (L)")            # Data S1 THETA(4) = 2.84; Xiong 2021 Table 2 'V2, L/65 kg' 2.84 (RSE 3.6%)

    # Covariate exponents on clearance and the volumes (power form on
    # reference-normalised covariates).
    e_wt_cl <- fixed(0.594); label("Power exponent of body weight on CL (unitless)")                        # Data S1 THETA(9) 'CLWT1' = 0.594; Xiong 2021 Table 2 'Weight for CL' 0.594 (RSE 10.2%). See vignette: the Data S1 code's non-missing branch cites THETA(10) instead
    e_alb_cl <- fixed(-1.07); label("Power exponent of serum albumin on CL (unitless)")                     # Data S1 THETA(5) 'CLALB1' = -1.07; Xiong 2021 Table 2 'ALB for CL' -1.07 (RSE 9.6%)
    e_conc_ada_cl <- fixed(0.134); label("Power exponent of antibody-to-infliximab concentration above 22 ng/mL on CL (unitless)") # Data S1 THETA(6) 'CLATI1' = 0.134; Xiong 2021 Table 2 'ATI for CL' 0.134 (RSE 9.6%)
    e_neutrophil_cd64_ratio_cl <- fixed(0.168); label("Power exponent of neutrophil CD64 activity ratio on CL (unitless)") # Data S1 THETA(7) 'CLNCAR_n' = 0.168; Xiong 2021 Table 2 'nCD64 for CL' 0.168 (RSE 16%)
    e_esr_cl <- fixed(0.101); label("Power exponent of erythrocyte sedimentation rate on CL (unitless)")     # Data S1 THETA(8) 'CLESR1' = 0.101; Xiong 2021 Table 2 'ESR for CL' 0.101 (RSE 27.7%)
    e_wt_vc <- fixed(0.55); label("Power exponent of body weight on V1 (unitless)")                         # Data S1 THETA(10) 'V1WT1' = 0.55; Xiong 2021 Table 2 'Weight for V1' 0.55 (RSE 20.2%)
    e_wt_q <- fixed(1); label("Power exponent of body weight on Q (unitless)")                              # Data S1 code Q=TVQ*EXP(ETA(3))*(WT/65); Xiong 2021 Table 2 'Weight for Q' 1, not estimated
    e_wt_vp <- fixed(0.586); label("Power exponent of body weight on V2 (unitless)")                        # Data S1 code V2=TVV2*EXP(ETA(4))*(WT/65)**0.586; Xiong 2021 Table 2 'Weight for V2' 0.586, not estimated

    # Biomarker trajectories (Samuels 2025 Methods equation and Table S1):
    # biomarker(n) = baseline * exp(Emax * (n - 1) / (D50 - 1 + n - 1)), n the
    # infusion number. Emax is the change on the log scale; exp(Emax) - 1 is
    # the asymptotic fractional change quoted in Results.
    wt_dosenum_max <- 0.046; label("Maximal log-scale change in body weight over successive infusions (unitless)")    # Table S1 WT Emax = 0.046 (bootstrap median 0.049, 95% CI 0.028-0.063); Data S1 E2 = EXP(0.046*INFUSION/(0.1+INFUSION))
    alb_dosenum_max <- 0.110; label("Maximal log-scale change in serum albumin over successive infusions (unitless)") # Table S1 ALB Emax = 0.110 (bootstrap median 0.110, 95% CI 0.094-0.127); Data S1 E1 = EXP(0.11*INFUSION/(0.1+INFUSION))
    esr_dosenum_max <- -0.978; label("Maximal log-scale change in erythrocyte sedimentation rate over successive infusions (unitless)") # Table S1 ESR Emax = -0.978 (bootstrap median -0.979); Data S1 E3 = EXP((-0.978)*INFUSION/(0.1+INFUSION))
    neutrophil_cd64_ratio_dosenum_max <- -0.362; label("Maximal log-scale change in neutrophil CD64 activity ratio over successive infusions (unitless)") # Table S1 nCD64 Emax = -0.362 (bootstrap median -0.404); Data S1 E4 = EXP((-0.362)*INFUSION/(0.1+INFUSION)); Figure 1d also prints -0.362
    dosenum50 <- fixed(1.1); label("Infusion number at half of the maximal biomarker change (count)")                 # Methods: 'D50 was fixed to 1.1'; Data S1 denominator 0.1 = D50 - 1

    # Inter-individual variability on CL only; Data S1 OMEGA(2)-(4) for V1, Q
    # and V2 are 0 and held there, matching Xiong 2021 Table 2, which reports
    # IIV on CL alone.
    etalcl ~ fixed(0.0629) # Data S1 OMEGA(1,1) = 0.0629; Xiong 2021 Table 2 omega CL 25.1% CV

    # Residual error additive on log-transformed concentrations (Data S1
    # $ERROR Y = LOG(F) + ERR(1)); SIGMA variance 0.189, so SD = sqrt(0.189).
    expSd <- fixed(0.4347); label("Residual error SD on the log scale (unitless)")                       # Data S1 $SIGMA 0.189; sqrt(0.189) = 0.4347; Xiong 2021 Table 2 sigma add 0.435
  })

  model({
    # 1. Biomarker trajectories. ninf is the zero-based count of previous
    #    infusions (the Data S1 INFUSION column); ftraj runs from 0 at the
    #    first infusion to 1 asymptotically, reaching 0.5 at infusion
    #    dosenum50.
    ninf <- CYCLE - 1
    if (ninf < 0) {
      ninf <- 0
    }
    ftraj <- ninf / ((dosenum50 - 1) + ninf)
    wt_pred <- WT * exp(wt_dosenum_max * ftraj)
    alb_pred <- ALB * 0.1 * exp(alb_dosenum_max * ftraj) # SI g/L -> g/dL, the scale of the Xiong 3.5 g/dL reference
    esr_pred <- ESR * exp(esr_dosenum_max * ftraj)
    ncd64_pred <- NEUTROPHIL_CD64_RATIO * exp(neutrophil_cd64_ratio_dosenum_max * ftraj)

    # 2. Antibody-to-infliximab effect: none at or below the 22 ng/mL limit of
    #    detection.
    cl_ada <- 1
    if (CONC_ADA_NGML > 22) {
      cl_ada <- (CONC_ADA_NGML / 22)^e_conc_ada_cl
    }

    # 3. Individual parameters (Data S1 $PK). Q and V2 use the baseline WT,
    #    not the predicted weight.
    cl <- exp(lcl + etalcl) *
      (wt_pred / 65)^e_wt_cl *
      (alb_pred / 3.5)^e_alb_cl *
      (esr_pred / 9)^e_esr_cl *
      (ncd64_pred / 4.6)^e_neutrophil_cd64_ratio_cl *
      cl_ada
    vc <- exp(lvc) * (wt_pred / 65)^e_wt_vc
    q <- exp(lq) * (WT / 65)^e_wt_q
    vp <- exp(lvp) * (WT / 65)^e_wt_vp

    # 4. Micro-constants and ODEs (ADVAN3 TRANS4).
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # 5. Observation: mg / L = ug/mL.
    Cc <- central / vc
    Cc ~ lnorm(expSd)
  })
}
