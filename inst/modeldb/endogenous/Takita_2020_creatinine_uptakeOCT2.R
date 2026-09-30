Takita_2020_creatinine_uptakeOCT2 <- function() {
  description <- paste(
    "PBPK (semi-mechanistic kidney model, MATLAB/Simulink). Physiologically-based",
    "creatinine model extended to chronic kidney disease (CKD), coupled to",
    "operational PK models of three renal-transporter inhibitors (trimethoprim,",
    "cimetidine, famotidine) to predict creatinine-drug interactions (Takita 2020),",
    "uptake-OCT2 variant. Creatinine is synthesised at a rate set by age, sex and",
    "weight, distributes into a central compartment, and is handled by the",
    "proximal tubule through three states (peritubular blood/interstitium, tubular",
    "cell, tubular filtrate): glomerular filtration, basolateral OAT2 uptake,",
    "basolateral OCT2 uptake of the cationic species only, apical MATE1 and MATE2-K",
    "efflux into the filtrate, and transcellular and paracellular passive",
    "permeability, with a fixed fraction of distal filtrate reabsorbed. CKD scales",
    "GFR from the covariates, proximal-tubule volumes and passive permeability in",
    "proportion to GFR (intact nephron hypothesis), and renal blood flow by CKD",
    "stage. Transporter activity follows the final non-INH scenario: OAT2 declines",
    "more than GFR and OCT2/MATE decline less than GFR (Eqs S9-S11); setting",
    "coeff_ckd_tp = 1, fx_oat2_slope = 0 and fx_oat2_int = 1 gives the INH",
    "scenario. Each inhibitor's unbound plasma concentration inhibits the",
    "transporters competitively (Eq S13). Serum creatinine starts at its analytic",
    "steady state.",
    sep = " "
  )
  reference <- paste(
    "Takita H, Scotcher D, Chinnadurai R, Kalra PA, Galetin A.",
    "Physiologically-Based Pharmacokinetic Modelling of Creatinine-Drug",
    "Interactions in the Chronic Kidney Disease Population.",
    "CPT Pharmacometrics Syst Pharmacol. 2020;9(12):695-706.",
    "doi:10.1002/psp4.12566.",
    "Model structure, system parameters and trimethoprim inputs are taken from the",
    "deposited MATLAB/Simulink code (Supporting Information PSP4-9-695-s002.zip:",
    "Creatinine_Inh_model_uptake.slx, DefineSystem.m, ReabIVIVE_Freab.m,",
    "Run_simulation.m, Drug_info_S.mat, CKDdata.xlsx) and the Supplementary",
    "Material (PSP4-9-695-s001.docx: Eqs S1-S16, Tables S6-S8, Text S1). The",
    "healthy-subject creatinine model it extends is Scotcher D et al.",
    "CPT Pharmacometrics Syst Pharmacol. 2020;9(5):310-321 (doi:10.1002/psp4.12509).",
    sep = " "
  )
  vignette <- "Takita_2020_creatinine"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    CRCL = list(
      description = "BSA-normalised estimated GFR (CKD-EPI equation from baseline serum creatinine, age, sex and race; Eq S1), or measured GFR from an exogenous marker when available",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Multiplied by BSA / 1.73 to give absolute GFR (Run_simulation.m column",
        "'eGFR_BSA'). Also sets the CKD stage for the renal-blood-flow reduction",
        "(>= 60: none; 30-59: 27%; < 30: 42%) and is the regressor of the",
        "inhibitor renal-clearance regressions used for Kp,uu,filtrate (Figures S13-S15)."
      ),
      source_name = "eGFR"
    ),
    BSA = list(
      description = "Body surface area (Du Bois formula, Eq S2)",
      units = "m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Converts the BSA-normalised CRCL into absolute GFR: GFR = CRCL * BSA / 1.73.",
      source_name = "BSA"
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Scales the creatinine synthesis rate (Eq S5) and the trimethoprim CL and V (Table S6, reference 70 kg).",
      source_name = "WT"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Creatinine synthesis rate regression (Eq S5, Bjornsson 1979).",
      source_name = "AGE"
    ),
    SEXF = list(
      description = "Sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      reference_category = "male (0)",
      notes = "Selects the sex-specific coefficients of the creatinine synthesis regression (Eq S5). Source column 'SEX(M1F0)' is 1 = male, so SEXF = 1 - SEX.",
      source_name = "SEX(M1F0)"
    )
  )

  # The three proximal-tubule states are specific to this model and have no
  # registered canonical. The inhibitor compartments carry drug-name suffixes
  # that are not in the register either. `central_creatinine` and
  # `urine_creatinine` use the registered `creatinine` co-analyte suffix.
  paper_specific_compartments <- c(
    "pt_blood",
    "pt_cell",
    "pt_filtrate",
    "depot_trimethoprim",
    "central_trimethoprim",
    "depot_cimetidine",
    "central_cimetidine",
    "peripheral1_cimetidine",
    "depot_famotidine",
    "central_famotidine",
    "peripheral1_famotidine"
  )
  paper_specific_etas <- c(
    "etalcl_trimethoprim",
    "etalvc_trimethoprim",
    "etalka_trimethoprim"
  )
  paper_specific_residual_sds <- c("propSd_trimethoprim")

  compartmentData <- list(
    central_creatinine = list(analyte = "creatinine", units = "mg", specimen = "serum", verified = TRUE),
    pt_blood = list(analyte = "creatinine", units = "mg", specimen = "tissue", verified = TRUE),
    pt_cell = list(analyte = "creatinine", units = "mg", specimen = "tissue", verified = TRUE),
    pt_filtrate = list(analyte = "creatinine", units = "mg", specimen = "urine", verified = TRUE),
    urine_creatinine = list(analyte = "creatinine", units = "mg", specimen = "urine", verified = TRUE),
    depot_trimethoprim = list(
      analyte = "trimethoprim",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central_trimethoprim = list(analyte = "trimethoprim", units = "mg", specimen = "plasma", verified = TRUE),
    depot_cimetidine = list(analyte = "cimetidine", units = "mg", specimen = "administration site", verified = TRUE),
    central_cimetidine = list(analyte = "cimetidine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_cimetidine = list(analyte = "cimetidine", units = "mg", specimen = "tissue", verified = TRUE),
    depot_famotidine = list(analyte = "famotidine", units = "mg", specimen = "administration site", verified = TRUE),
    central_famotidine = list(analyte = "famotidine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1_famotidine = list(analyte = "famotidine", units = "mg", specimen = "tissue", verified = TRUE)
  )

  dosing <- c(
    "depot_trimethoprim",
    "depot_cimetidine",
    "central_cimetidine",
    "depot_famotidine",
    "central_famotidine"
  )

  population <- list(
    species = "human",
    n_subjects = 64L,
    n_studies = 8L,
    age_range = "22-88 years",
    weight_range = "44-96 kg",
    sex_female_pct = 51.6,
    disease_state = "chronic kidney disease, stages G3-G4 (eGFR 15-59 mL/min/1.73 m^2), non-dialysis",
    dose_range = paste(
      "Trimethoprim 100-200 mg/day (Salford Kidney Study) and 100-320 mg/day oral",
      "(literature); cimetidine 200-1200 mg/day oral or 300 mg IV; famotidine",
      "20 mg twice daily oral or 10 mg IV"
    ),
    regions = "United Kingdom (Salford Kidney Study) and published literature",
    notes = paste(
      "Creatinine model optimised against baseline serum creatinine in 64 CKD",
      "patients from 8 clinical studies (35 G3, 29 G4; 31 men, 33 women;",
      "Table 1). Verified against 42 independent CKD patients (Table S8).",
      "Creatinine-drug interactions evaluated in 12 studies (90 patients, G3-4).",
      "Trimethoprim popPK fitted to 9 patients (Rieder 1974, 160 mg oral single",
      "dose; Text S1); cimetidine and famotidine models fitted to mean",
      "profiles by naive pooling (Table 3)."
    )
  )

  ini({
    # ---------------------------------------------------------------
    # Creatinine system parameters, healthy reference values carried from the
    # Scotcher 2020 healthy-subject model (Table 2 'Value in healthy' column;
    # exact values from the deposited DefineSystem.m and Run_simulation.m).
    gfr_healthy <- fixed(7.5); label("GFR of a healthy reference subject (L/h; 125 mL/min)") # Table 2 GFR 7.5 L/h; Eq S6 GFRhealthy 125 mL/min
    vd_creatinine <- fixed(43.7); label("Total creatinine volume of distribution (L)") # DefineSystem.m 'Vold' 43.7; unchanged in CKD (Methods)
    cl_nr_creatinine <- fixed(0.17); label("Non-renal creatinine clearance (L/h)") # DefineSystem.m 'CLnr' 0.17
    q_pt_healthy <- fixed(58.449); label("Blood flow to the proximal tubule, healthy (L/h)") # DefineSystem.m 'Qrpt' 58.449; Table 2 Q_PT,blood 58
    f_qr_g3 <- fixed(0.73); label("Renal blood flow fraction of healthy, CKD G3 (unitless)") # Methods 'decreased by 27% ... in CKD G3'; CKDdata.xlsx 'Qrenal(% normal)' 0.73
    f_qr_g4 <- fixed(0.58); label("Renal blood flow fraction of healthy, CKD G4 (unitless)") # Methods 'decreased by ... 42% in CKD G3 and G4'; CKDdata.xlsx 0.58
    v_ptbi_healthy <- fixed(0.08176); label("Blood and interstitial water volume of the cortex, healthy (L)") # DefineSystem.m 'Volptbi' 0.08176; Table 2 V_PT,bi 0.082
    v_ptc_healthy <- fixed(0.066076); label("Proximal tubule cell water volume, healthy (L)") # DefineSystem.m 'Volptc' 0.066076; Table 2 V_PT,cell 0.066
    v_ptfilt_healthy <- fixed(0.053529); label("Proximal tubule filtrate volume, healthy (L)") # DefineSystem.m 'Volptfilt' 0.053529; Table 2 V_PT,filt 0.054
    f_reab_water <- fixed(0.64); label("Fraction of filtered water reabsorbed in the proximal tubule (unitless)") # DefineSystem.m 'PT_H20_reab' 0.64; Section S7 F_reab,water 0.64
    pka_creatinine <- fixed(4.74); label("Creatinine pKa (base)") # Run_simulation.m 'creatinine_pKa' 4.74; Section S4
    ph_blood <- fixed(7.4); label("pH of the peritubular blood compartment") # DefineSystem.m 'pH_b' 7.4

    # Passive permeability of creatinine: optimised in vitro apparent
    # permeability scaled by the tubular surface area (ReabIVIVE_Freab.m).
    papp_creatinine <- fixed(28.87); label("Creatinine apparent permeability, uptake-OCT2 model (1e-6 cm/s)") # Run_simulation.m 'papp' (uptake) 28.87; reproduces Table 2 CL_PD,trans 0.89 and CL_PD,para 5.9 L/h
    sa_pt <- fixed(61072.5611857856); label("Proximal tubule surface area (cm^2)") # ReabIVIVE_Freab.m
    f_transcellular <- fixed(0.07); label("Fraction of proximal tubule permeability through the transcellular route (unitless)") # ReabIVIVE_Freab.m 'CLpd_pt*0.07*2' and 'CLpd_pt*0.93'
    sa_loh <- fixed(1628.60163); label("Loop of Henle surface area (cm^2)") # ReabIVIVE_Freab.m
    sa_dt <- fixed(2073.45115); label("Distal tubule surface area (cm^2)") # ReabIVIVE_Freab.m
    sa_cd <- fixed(453.48631793); label("Collecting duct surface area (cm^2)") # ReabIVIVE_Freab.m
    q_loh <- fixed(35); label("Filtrate flow in the loop of Henle, healthy (mL/min)") # ReabIVIVE_Freab.m
    q_dt <- fixed(17.8); label("Filtrate flow in the distal tubule, healthy (mL/min)") # ReabIVIVE_Freab.m
    q_cd <- fixed(6.3); label("Filtrate flow in the collecting duct, healthy (mL/min)") # ReabIVIVE_Freab.m

    # Creatinine synthesis rate regression (Bjornsson 1979): (C0 - C1 * Age) * WT / 24
    rsyn_c0_male <- fixed(27); label("Creatinine synthesis intercept, men (mg/kg/day)") # Eq S5 C0 = 27 for male
    rsyn_c1_male <- fixed(0.173); label("Creatinine synthesis age slope, men (mg/kg/day/year)") # Eq S5 C1 = 0.173 for male
    rsyn_c0_female <- fixed(25); label("Creatinine synthesis intercept, women (mg/kg/day)") # Eq S5 C0 = 25 for female
    rsyn_c1_female <- fixed(0.175); label("Creatinine synthesis age slope, women (mg/kg/day/year)") # Eq S5 C1 = 0.175 for female

    # Healthy transporter intrinsic clearances of the uptake-OCT2 model
    # (Run_simulation.m 'Optim_list.CLint(1,:)'). clint_oct2_healthy acts on
    # the cationic fraction only; times that fraction at pH 7.4 (0.002183) it
    # is the 23.9 L/h printed in Table 2.
    clint_oat2_healthy <- fixed(20.771); label("OAT2 intrinsic clearance, healthy (L/h)") # Run_simulation.m 20.771; Table 2 CL_int,OAT2 (Uptake) 20.8
    clint_oct2_healthy <- fixed(10950.1960990373); label("OCT2 intrinsic clearance of cationic creatinine, healthy (L/h)") # Run_simulation.m 10950.196; x cationic fraction = Table 2 CL_int,OCT2 (Uptake) 23.9
    clint_mate1_healthy <- fixed(0.157); label("MATE1 intrinsic clearance, healthy (L/h)") # Run_simulation.m 0.157; Table 2 CL_int,MATE1 (Uptake) 0.16
    clint_mate2k_healthy <- fixed(0.510); label("MATE2-K intrinsic clearance, healthy (L/h)") # Run_simulation.m 0.510; Table 2 CL_int,MATE2-K (Uptake) 0.51

    # CKD transporter scaling, non-INH scenario
    fx_oat2_slope <- 0.4973; label("Slope of the additional OAT2 deterioration Fx_OAT2 on GFR_CKD/GFR_healthy (unitless)") # Eq S9 0.4973 (deposited code 0.4972)
    fx_oat2_int <- 0.5027; label("Intercept of the additional OAT2 deterioration Fx_OAT2 (unitless)") # Eq S9 0.5027
    coeff_ckd_tp <- 0.8763; label("Coeff_CKD,TP: slope of relative OCT2/MATE clearance on relative GFR, uptake-OCT2 model (unitless)") # Run_simulation.m 0.8763; Results 'Coeff_CKD,TP ... 0.88'

    # Inhibitory concentration at the apical MATE transporters: 0 = unbound
    # plasma concentration (primary analysis, Figure 4, Table S9); 1 = unbound
    # proximal-tubular-filtrate concentration C_PT,filt (Figure S17, Table S10).
    mate_ic_filtrate <- fixed(0); label("Use C_PT,filt instead of C_p,u for MATE inhibition (0 or 1)") # Run_simulation.m 'Kpuuapply=0'; Section S7

    # ---------------------------------------------------------------
    # Trimethoprim: one-compartment popPK with first-order absorption, fitted
    # in NONMEM to 9 CKD patients (Table S7, Text S1).
    lcl_trimethoprim <- log(5.07); label("Trimethoprim clearance at 70 kg (L/h)") # Table S7 CL 5.07 L/h; Text S1 THETA(1)
    lvc_trimethoprim <- log(105); label("Trimethoprim volume of distribution at 70 kg (L)") # Table S7 Vd 105 L; Text S1 THETA(2)
    lka_trimethoprim <- log(1.24); label("Trimethoprim absorption rate constant (1/h)") # Table S7 ka 1.24 1/h (Text S1 $THETA 1.23)
    e_wt_cl_trimethoprim <- 0.822; label("Power exponent of WT/70 on trimethoprim CL (unitless)") # Table S7 'WT on CL' 0.822; Text S1 THETA(5)
    e_wt_vc_trimethoprim <- 0.532; label("Power exponent of WT/70 on trimethoprim V (unitless)") # Table S7 'WT on Vd' 0.532; Text S1 THETA(4)
    fu_trimethoprim <- fixed(0.51); label("Trimethoprim fraction unbound in plasma (unitless)") # Table S6 f_u,p 0.51
    mw_trimethoprim <- fixed(290.32); label("Trimethoprim molecular weight (g/mol)") # Drug_info_S.mat 'Mr' 290.32
    ic50_oat2_trimethoprim <- fixed(1000); label("Trimethoprim IC50 on OAT2 (uM)") # Table S8 1000 (no inhibition observed; set to 1000 uM)
    ic50_oct2_trimethoprim <- fixed(25.8); label("Trimethoprim IC50 on OCT2 (uM)") # Table S8 25.8
    ic50_mate1_trimethoprim <- fixed(1.62); label("Trimethoprim IC50 on MATE1 (uM)") # Table S8 1.62
    ic50_mate2k_trimethoprim <- fixed(0.58); label("Trimethoprim IC50 on MATE2-K (uM)") # Table S8 0.58
    clr_slope_trimethoprim <- fixed(0.6525); label("Slope of trimethoprim renal clearance on GFR (unitless)") # Figure S13A 'Renal clearance = 0.6525*GFR + 18.12'
    clr_int_trimethoprim <- fixed(18.12); label("Intercept of trimethoprim renal clearance on GFR (mL/min)") # Figure S13A

    # ---------------------------------------------------------------
    # Cimetidine: two-compartment model, naive-pooled fit to mean profiles.
    lcl_cimetidine <- log(27.30); label("Cimetidine clearance (L/h)") # Table S6 CL 27.30
    lvc_cimetidine <- log(35.60); label("Cimetidine central volume (L)") # Table S6 V2 35.60
    lka_cimetidine <- log(0.23); label("Cimetidine absorption rate constant (1/h)") # Table S6 ka 0.23
    lq_cimetidine <- log(38.80); label("Cimetidine intercompartmental clearance (L/h)") # Table S6 Q 38.80
    lvp_cimetidine <- log(39.30); label("Cimetidine peripheral volume (L)") # Table S6 V3 39.30
    lfdepot_cimetidine <- fixed(log(1)); label("Cimetidine oral bioavailability (unitless)") # Table S6 F 1
    fu_cimetidine <- fixed(0.84); label("Cimetidine fraction unbound in plasma, CKD (unitless)") # Table S6 f_u,p 0.84
    mw_cimetidine <- fixed(252.34); label("Cimetidine molecular weight (g/mol)") # computed from the molecular formula C10H16N6S; not printed in the source
    ic50_oat2_cimetidine <- fixed(102.3); label("Cimetidine IC50 on OAT2 (uM)") # Table S8 102.3
    ic50_oct2_cimetidine <- fixed(36.3); label("Cimetidine IC50 on OCT2 (uM)") # Table S8 36.3
    ic50_mate1_cimetidine <- fixed(3.78); label("Cimetidine IC50 on MATE1 (uM)") # Table S8 3.78
    ic50_mate2k_cimetidine <- fixed(23.7); label("Cimetidine IC50 on MATE2-K (uM)") # Table S8 23.7
    clr_slope_cimetidine <- fixed(1.5327); label("Slope of cimetidine renal clearance on GFR (unitless)") # Figure S14A 'Renal clearance = 1.5327*GFR + 84.495'
    clr_int_cimetidine <- fixed(84.495); label("Intercept of cimetidine renal clearance on GFR (mL/min)") # Figure S14A

    # ---------------------------------------------------------------
    # Famotidine: two-compartment model, naive-pooled fit to mean profiles.
    lcl_famotidine <- log(2.61); label("Famotidine clearance (L/h)") # Table S6 CL 2.61
    lvc_famotidine <- log(24.19); label("Famotidine central volume (L)") # Table S6 V2 24.19
    lka_famotidine <- log(0.44); label("Famotidine absorption rate constant (1/h)") # Table S6 ka 0.44
    lq_famotidine <- log(74.24); label("Famotidine intercompartmental clearance (L/h)") # Table S6 Q 74.24
    lvp_famotidine <- log(46.67); label("Famotidine peripheral volume (L)") # Table S6 V3 46.67
    lfdepot_famotidine <- log(0.47); label("Famotidine oral bioavailability (unitless)") # Table S6 F 0.47
    fu_famotidine <- fixed(0.72); label("Famotidine fraction unbound in plasma (unitless)") # Table S6 f_u,p 0.72
    mw_famotidine <- fixed(337.44); label("Famotidine molecular weight (g/mol)") # computed from the molecular formula C8H15N7O2S3; not printed in the source
    ic50_oat2_famotidine <- fixed(184); label("Famotidine IC50 on OAT2 (uM)") # Table S8 184
    ic50_oct2_famotidine <- fixed(27.9); label("Famotidine IC50 on OCT2 (uM)") # Table S8 27.9
    ic50_mate1_famotidine <- fixed(0.27); label("Famotidine IC50 on MATE1 (uM)") # Table S8 0.27
    ic50_mate2k_famotidine <- fixed(7.3); label("Famotidine IC50 on MATE2-K (uM)") # Table S8 7.3
    clr_slope_famotidine <- fixed(2.8172); label("Slope of famotidine renal clearance on creatinine clearance (unitless)") # Figure S15A 'Renal clearance = 2.8172*C_Cr - 7.1038'
    clr_int_famotidine <- fixed(-7.1038); label("Intercept of famotidine renal clearance on creatinine clearance (mL/min)") # Figure S15A

    # ---------------------------------------------------------------
    # Trimethoprim between-subject variability and residual error
    etalcl_trimethoprim ~ 0.0552 # Text S1 $OMEGA 0.0552 (Table S7 IIV 23.5%)
    etalvc_trimethoprim ~ 0.0157 # Text S1 $OMEGA 0.0157 (Table S7 IIV 12.5%)
    etalka_trimethoprim ~ 0.0817 # Text S1 $OMEGA 0.0817 (Table S7 IIV 28.6%)
    propSd_trimethoprim <- 0.08497; label("Trimethoprim proportional residual error (fraction)") # Table S7 Sigma 0.00722 (variance); sqrt = 0.08497
  })

  model({
    # --- CKD-dependent system parameters (Section S4) ---
    gfr <- CRCL * BSA / 1.73 * 60 / 1000
    gfr_rel <- gfr / gfr_healthy
    rsyn <- ((1 - SEXF) * (rsyn_c0_male - rsyn_c1_male * AGE) +
      SEXF * (rsyn_c0_female - rsyn_c1_female * AGE)) * WT / 24
    f_qr <- 1
    if (CRCL < 60) {
      f_qr <- f_qr_g3
    }
    if (CRCL < 30) {
      f_qr <- f_qr_g4
    }
    q_pt <- q_pt_healthy * f_qr
    v_ptbi <- v_ptbi_healthy * gfr_rel
    v_ptc <- v_ptc_healthy * gfr_rel
    v_ptfilt <- v_ptfilt_healthy * gfr_rel
    vc_creatinine <- vd_creatinine - v_ptbi - v_ptc - v_ptfilt
    q_ptdt <- gfr * (1 - f_reab_water)

    # --- Passive permeability and distal reabsorption (ReabIVIVE_Freab.m) ---
    clpd_pt <- papp_creatinine / 1e6 * 60 * sa_pt * 60 / 1000 * gfr_rel
    clpd_mem <- clpd_pt * f_transcellular * 2
    clpd_para <- clpd_pt * (1 - f_transcellular)
    clint_loh <- papp_creatinine / 1e6 * 60 * sa_loh
    clint_dt <- papp_creatinine / 1e6 * 60 * sa_dt
    clint_cd <- papp_creatinine / 1e6 * 60 * sa_cd
    f_reab_dt <- 1 - (1 - clint_loh / (clint_loh + q_loh)) *
      (1 - clint_dt / (clint_dt + q_dt)) *
      (1 - clint_cd / (clint_cd + q_cd))

    # --- Transporter clearances in CKD, non-INH scenario (Eqs S9-S11) ---
    clint_oat2 <- clint_oat2_healthy * gfr_rel * (fx_oat2_slope * gfr_rel + fx_oat2_int)
    f_tp <- coeff_ckd_tp * (gfr_rel - 1) + 1
    clint_oct2 <- clint_oct2_healthy * f_tp
    clint_mate1 <- clint_mate1_healthy * f_tp
    clint_mate2k <- clint_mate2k_healthy * f_tp

    # --- OCT2: uptake of the cationic fraction of creatinine only ---
    f_cation_blood <- 10^(pka_creatinine - ph_blood) / (10^(pka_creatinine - ph_blood) + 1)
    oct2_in <- clint_oct2 * f_cation_blood
    oct2_out <- 0

    # --- Analytic steady state of the uninhibited creatinine system ---
    # Blood, cell and filtrate concentrations per unit central concentration
    # (Cramer's rule), then the central concentration from the whole-body
    # balance rsyn = cl_nr * C + urinary excretion rate.
    ss_m11 <- -(clpd_para + q_pt + clpd_mem + clint_oat2 + oct2_in)
    ss_m12 <- clpd_mem + oct2_out
    ss_m13 <- clpd_para
    ss_m21 <- clint_oat2 + oct2_in + clpd_mem
    ss_m22 <- -(2 * clpd_mem + clint_mate1 + clint_mate2k + oct2_out)
    ss_m23 <- clpd_mem
    ss_m31 <- clpd_para
    ss_m32 <- clint_mate1 + clint_mate2k + clpd_mem
    ss_m33 <- -(q_ptdt + clpd_mem + clpd_para)
    ss_r1 <- -q_pt
    ss_r3 <- -gfr
    ss_det <- ss_m11 * (ss_m22 * ss_m33 - ss_m23 * ss_m32) -
      ss_m12 * (ss_m21 * ss_m33 - ss_m23 * ss_m31) +
      ss_m13 * (ss_m21 * ss_m32 - ss_m22 * ss_m31)
    ss_xb <- (ss_r1 * (ss_m22 * ss_m33 - ss_m23 * ss_m32) -
      ss_m12 * (-ss_m23 * ss_r3) +
      ss_m13 * (-ss_m22 * ss_r3)) / ss_det
    ss_xt <- (ss_m11 * (-ss_m23 * ss_r3) -
      ss_r1 * (ss_m21 * ss_m33 - ss_m23 * ss_m31) +
      ss_m13 * (ss_m21 * ss_r3)) / ss_det
    ss_xf <- (ss_m11 * (ss_m22 * ss_r3) -
      ss_m12 * (ss_m21 * ss_r3) +
      ss_r1 * (ss_m21 * ss_m32 - ss_m22 * ss_m31)) / ss_det
    css_creatinine <- rsyn / (cl_nr_creatinine + q_ptdt * (1 - f_reab_dt) * ss_xf)

    # --- Inhibitor PK ---
    cl_trimethoprim <- exp(lcl_trimethoprim + etalcl_trimethoprim) * (WT / 70)^e_wt_cl_trimethoprim
    vc_trimethoprim <- exp(lvc_trimethoprim + etalvc_trimethoprim) * (WT / 70)^e_wt_vc_trimethoprim
    ka_trimethoprim <- exp(lka_trimethoprim + etalka_trimethoprim)
    cl_cimetidine <- exp(lcl_cimetidine)
    vc_cimetidine <- exp(lvc_cimetidine)
    ka_cimetidine <- exp(lka_cimetidine)
    q_cimetidine <- exp(lq_cimetidine)
    vp_cimetidine <- exp(lvp_cimetidine)
    cl_famotidine <- exp(lcl_famotidine)
    vc_famotidine <- exp(lvc_famotidine)
    ka_famotidine <- exp(lka_famotidine)
    q_famotidine <- exp(lq_famotidine)
    vp_famotidine <- exp(lvp_famotidine)

    Cc_trimethoprim <- central_trimethoprim / vc_trimethoprim
    Cc_cimetidine <- central_cimetidine / vc_cimetidine
    Cc_famotidine <- central_famotidine / vc_famotidine
    # Unbound plasma concentrations in uM (mg/L / g/mol * 1000)
    cu_trimethoprim <- fu_trimethoprim * Cc_trimethoprim / mw_trimethoprim * 1000
    cu_cimetidine <- fu_cimetidine * Cc_cimetidine / mw_cimetidine * 1000
    cu_famotidine <- fu_famotidine * Cc_famotidine / mw_famotidine * 1000

    # Unbound filtrate-to-plasma ratio Kp,uu,filtrate (Eqs S14-S16) for the
    # optional C_PT,filt inhibition of MATE1 / MATE2-K
    gfr_mlmin <- CRCL * BSA / 1.73
    kpuu_trimethoprim <- (clr_slope_trimethoprim * CRCL + clr_int_trimethoprim) /
      (gfr_mlmin * (1 - f_reab_water) * fu_trimethoprim)
    kpuu_cimetidine <- (clr_slope_cimetidine * CRCL + clr_int_cimetidine) /
      (gfr_mlmin * (1 - f_reab_water) * fu_cimetidine)
    kpuu_famotidine <- (clr_slope_famotidine * CRCL + clr_int_famotidine) /
      (gfr_mlmin * (1 - f_reab_water) * fu_famotidine)
    kmate_trimethoprim <- 1 + mate_ic_filtrate * (kpuu_trimethoprim - 1)
    kmate_cimetidine <- 1 + mate_ic_filtrate * (kpuu_cimetidine - 1)
    kmate_famotidine <- 1 + mate_ic_filtrate * (kpuu_famotidine - 1)

    # Fraction of transporter activity remaining (Eq S13)
    finh_oat2 <- 1 / (1 + cu_trimethoprim / ic50_oat2_trimethoprim +
      cu_cimetidine / ic50_oat2_cimetidine + cu_famotidine / ic50_oat2_famotidine)
    finh_oct2 <- 1 / (1 + cu_trimethoprim / ic50_oct2_trimethoprim +
      cu_cimetidine / ic50_oct2_cimetidine + cu_famotidine / ic50_oct2_famotidine)
    finh_mate1 <- 1 / (1 + kmate_trimethoprim * cu_trimethoprim / ic50_mate1_trimethoprim +
      kmate_cimetidine * cu_cimetidine / ic50_mate1_cimetidine +
      kmate_famotidine * cu_famotidine / ic50_mate1_famotidine)
    finh_mate2k <- 1 / (1 + kmate_trimethoprim * cu_trimethoprim / ic50_mate2k_trimethoprim +
      kmate_cimetidine * cu_cimetidine / ic50_mate2k_cimetidine +
      kmate_famotidine * cu_famotidine / ic50_mate2k_famotidine)

    # --- Creatinine fluxes (mg/h), Creatinine_Inh_model_uptake.slx ---
    c_cent <- central_creatinine / vc_creatinine
    c_blood <- pt_blood / v_ptbi
    c_cell <- pt_cell / v_ptc
    c_filt <- pt_filtrate / v_ptfilt
    j_oat2 <- clint_oat2 * finh_oat2 * c_blood
    j_oct2 <- finh_oct2 * (oct2_in * c_blood - oct2_out * c_cell)
    j_mate <- (clint_mate1 * finh_mate1 + clint_mate2k * finh_mate2k) * c_cell
    j_out_pt <- q_ptdt * c_filt

    d/dt(central_creatinine) <- rsyn + j_out_pt * f_reab_dt + q_pt * c_blood -
      (cl_nr_creatinine + gfr + q_pt) * c_cent
    d/dt(pt_blood) <- q_pt * c_cent + clpd_mem * c_cell + clpd_para * c_filt -
      (clpd_para + q_pt + clpd_mem) * c_blood - j_oct2 - j_oat2
    d/dt(pt_cell) <- j_oat2 + j_oct2 + clpd_mem * c_blood + clpd_mem * c_filt -
      2 * clpd_mem * c_cell - j_mate
    d/dt(pt_filtrate) <- j_mate + clpd_mem * c_cell + gfr * c_cent + clpd_para * c_blood -
      (q_ptdt + clpd_mem + clpd_para) * c_filt
    d/dt(urine_creatinine) <- j_out_pt * (1 - f_reab_dt)

    d/dt(depot_trimethoprim) <- -ka_trimethoprim * depot_trimethoprim
    d/dt(central_trimethoprim) <- ka_trimethoprim * depot_trimethoprim -
      cl_trimethoprim / vc_trimethoprim * central_trimethoprim
    d/dt(depot_cimetidine) <- -ka_cimetidine * depot_cimetidine
    d/dt(central_cimetidine) <- ka_cimetidine * depot_cimetidine -
      (cl_cimetidine + q_cimetidine) / vc_cimetidine * central_cimetidine +
      q_cimetidine / vp_cimetidine * peripheral1_cimetidine
    d/dt(peripheral1_cimetidine) <- q_cimetidine / vc_cimetidine * central_cimetidine -
      q_cimetidine / vp_cimetidine * peripheral1_cimetidine
    d/dt(depot_famotidine) <- -ka_famotidine * depot_famotidine
    d/dt(central_famotidine) <- ka_famotidine * depot_famotidine -
      (cl_famotidine + q_famotidine) / vc_famotidine * central_famotidine +
      q_famotidine / vp_famotidine * peripheral1_famotidine
    d/dt(peripheral1_famotidine) <- q_famotidine / vc_famotidine * central_famotidine -
      q_famotidine / vp_famotidine * peripheral1_famotidine

    f(depot_cimetidine) <- exp(lfdepot_cimetidine)
    f(depot_famotidine) <- exp(lfdepot_famotidine)

    # Start the creatinine system at its steady state
    central_creatinine(0) <- css_creatinine * vc_creatinine
    pt_blood(0) <- ss_xb * css_creatinine * v_ptbi
    pt_cell(0) <- ss_xt * css_creatinine * v_ptc
    pt_filtrate(0) <- ss_xf * css_creatinine * v_ptfilt

    # --- Outputs ---
    Cc_creatinine <- c_cent
    scr <- Cc_creatinine / 10
    scr_baseline <- css_creatinine / 10
    pct_change_scr <- 100 * (Cc_creatinine - css_creatinine) / css_creatinine
    clr_creatinine <- j_out_pt * (1 - f_reab_dt) / c_cent
    ccr_gfr_ratio <- clr_creatinine / gfr

    Cc_trimethoprim ~ prop(propSd_trimethoprim)
  })
}
