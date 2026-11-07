Yan_2023_tafolecimab <- function() {
  description <- "Sequential population PK/PD model for tafolecimab (anti-PCSK9 IgG2 monoclonal antibody) and serum LDL cholesterol in Chinese healthy volunteers and patients with non-familial, heterozygous familial or homozygous familial hypercholesterolemia (Yan 2023). Two-compartment PK with first-order SC absorption, fixed bioavailability, and parallel linear plus Michaelis-Menten (target-mediated approximation) elimination from the central compartment; a type IV (stimulation-of-loss) indirect response model links free tafolecimab concentration to LDL-C elimination. Body weight, sex and baseline unbound PCSK9 act on the PK, and baseline LDL-C, sex and disease type act on the LDL-C parameters."
  reference <- "Yan S, Zhao X, Xie Q, Du W, Ma Q, Zhu T, Deng H, Qian L, Zheng S, Cui Y. Pharmacokinetic/LDL-C and exposure-response analysis of tafolecimab in Chinese hypercholesterolemia patients: Results from phase I, II, and III studies. Clin Transl Sci. 2023;16:2791-2803. doi:10.1111/cts.13674"
  vignette <- "Yan_2023_tafolecimab"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # The three PK states hold amounts. `ldl` is the LDL-C turnover pool and is
  # carried as a serum CONCENTRATION (mg/dL), because Equation 2 is written
  # directly on the measured LDL-C concentration with no volume term.
  compartmentData <- list(
    depot = list(analyte = "tafolecimab", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "tafolecimab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "tafolecimab", units = "mg", specimen = "serum", verified = TRUE),
    ldl = list(analyte = "low-density lipoprotein cholesterol", units = "mg/dL", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Yan 2023 calls this BBWT. Time-fixed (baseline). Power effects (WT/70.7)^beta on ka (-0.74), central volume (1.15) and linear clearance (1.4), Table 2; 70.7 kg is the median of the pooled data set (Table 2 note, Table S1; range 41.2-123 kg).",
      source_name = "BBWT"
    ),
    SEXF = list(
      description = "Biological sex indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Yan 2023 codes SEX as M / F with male as the reference (Equation 5); SEXF = 1 for SEX = F. Exponential shifts exp(beta * SEXF): ka -0.73, central volume +0.24, linear clearance +0.3, LDL-C Emax -0.15 (Table 2). 414/1153 (35.9%) female (Table S1).",
      source_name = "SEX"
    ),
    PCSK9 = list(
      description = "Baseline unbound (free) serum PCSK9 concentration",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Yan 2023 calls this BPCSK9, measured by a ligand-binding assay that captures unbound PCSK9 on immobilized tafolecimab (Bioanalysis section). Power effect (PCSK9/775.762)^0.22 on the central volume (Table 2); 775.762 ng/mL is the pooled median (Table 2 note; Table S1 median 775.76, range 64.71-2454.06 ng/mL). Missing values were imputed with the median of the same study population.",
      source_name = "BPCSK9"
    ),
    LDLC = list(
      description = "Baseline serum LDL cholesterol concentration",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Yan 2023 calls this BLDLC. Time-fixed. Power effects (LDLC/112.1285)^beta on the indirect-response baseline R0 (0.94) and on kout (-1.63), Table 2; 112.1285 mg/dL is the pooled median (Table 2 note; Table S1 median 112.13, range 57.22-586.93 mg/dL). Dual role: it is the covariate that sets the individual baseline rbase, which is in turn the initial condition of the `ldl` state and anchors kin = kout * rbase.",
      source_name = "BLDLC"
    ),
    DIS_HEALTHY = list(
      description = "Healthy-volunteer indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-familial hypercholesterolemia patient, when DIS_HEFH = 0 and DIS_HOFH = 0)",
      notes = "One of three indicators derived from the four-level Yan 2023 DISEASE covariate (HV, non_FH, HeFH, HoFH); non_FH is the reference level, so at most one of DIS_HEALTHY, DIS_HEFH and DIS_HOFH is 1. Exponential shift exp(-0.081) on Emax (Table 2), giving Emax 2.19 in male HVs (Results). 44/1153 (3.8%), all in the phase I study (Table S1).",
      source_name = "DISEASE"
    ),
    DIS_HEFH = list(
      description = "Heterozygous familial hypercholesterolemia indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-familial hypercholesterolemia patient, when DIS_HEALTHY = 0 and DIS_HOFH = 0)",
      notes = "DISEASE = HeFH in Yan 2023. Exponential shift exp(-0.37) on Emax (Table 2), giving Emax 1.64 in male HeFH patients (Results). 187/1153 (16.2%) (Table S1).",
      source_name = "DISEASE"
    ),
    DIS_HOFH = list(
      description = "Homozygous familial hypercholesterolemia indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-familial hypercholesterolemia patient, when DIS_HEALTHY = 0 and DIS_HEFH = 0)",
      notes = "DISEASE = HoFH in Yan 2023. Exponential shift exp(-3.07) on Emax (Table 2), giving Emax 0.11 in male HoFH patients (Results). 18/1153 (1.6%), all in study NCT04031742 (Table S1).",
      source_name = "DISEASE"
    )
  )

  # Covariates screened (Covariate screening section, Table S1) but not
  # retained in the final model.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened; not retained. Median 57 (18-76) years (Table S1)."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      notes = "Screened; not retained. Median 45.91 g/L (Table S1)."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Screened; not retained. Median 22 U/L (Table S1)."
    ),
    BILI = list(
      description = "Total bilirubin",
      units = "umol/L",
      type = "continuous",
      notes = "Screened as TBIL; not retained. Median 13 umol/L (Table S1)."
    ),
    CRCL = list(
      description = "Creatinine clearance",
      units = "mL/min",
      type = "continuous",
      notes = "Screened as CLCR; not retained. Median 101.84 mL/min (Table S1)."
    ),
    ADA_POS = list(
      description = "Treatment-induced anti-drug antibody positive",
      units = "(binary)",
      type = "binary",
      notes = "Screened as ADA; not retained. 118/1153 (10.2%) positive (Table S1)."
    ),
    CONMED_STATIN = list(
      description = "Concomitant statin",
      units = "(binary)",
      type = "binary",
      notes = "Screened as STATIN (0 none, 2 moderate, 3 high intensity); not retained because almost all patients received a statin and BPCSK9 absorbs the statin-induced PCSK9 up-regulation (Discussion)."
    ),
    CONMED_EZE = list(
      description = "Concomitant ezetimibe",
      units = "(binary)",
      type = "binary",
      notes = "Screened as EZETIMIBE; not retained. 142/1153 (12.3%) (Table S1)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 1153L,
    n_subjects_ldlc = 1149L,
    n_observations = "9527 tafolecimab concentrations; 12,007 LDL-C values",
    n_studies = 6L,
    phases = "Phase I (NCT03366688), phase IIa (NCT03815812) and four phase III studies (NCT04031742, NCT04709536, NCT04289285, NCT04179669)",
    age_range = "18-76 years",
    age_median = "57 years",
    weight_range = "41.2-123 kg",
    weight_median = "70.7 kg",
    sex_female_pct = 35.9,
    race_ethnicity = "Chinese (all studies conducted in China)",
    disease_state = "904 (78.4%) non-familial hypercholesterolemia, 187 (16.2%) heterozygous familial hypercholesterolemia, 18 (1.6%) homozygous familial hypercholesterolemia, 44 (3.8%) healthy volunteers.",
    dose_range = "SC 25-600 mg single dose and IV 75 or 450 mg single dose (phase I); SC 75-600 mg every 2, 4, 6 or 8 weeks (phase II/III); pivotal regimens 150 mg Q2W, 450 mg Q4W and 600 mg Q6W.",
    regions = "China",
    co_medication = "Statin 1102 (95.6%; moderate intensity 89.7%, high intensity 5.9%); ezetimibe 142 (12.3%).",
    baseline_ldlc = "Median 112.13 mg/dL (57.22-586.93)",
    baseline_pcsk9 = "Unbound PCSK9 median 775.76 ng/mL (64.71-2454.06)",
    notes = "Baseline characteristics from Table S1 (median (min, max)). Concentrations below the LLOQ (0.4 ug/mL) were handled by the M4 method. The analyte is free (unbound) tafolecimab."
  )

  ini({
    # PK layer - Yan 2023 Table 2, 'PopPK final model'. Rates in the paper are
    # per hour, volumes in L, Vm in mg/h and Km in mg/L; serum concentration is
    # therefore mg/L (= ug/mL).
    lka <- log(0.01)
    label("Absorption rate constant ka for a 70.7 kg male (1/h)") # Table 2 'Absorption rate constant Ka (1/h)' 0.01
    lvc <- log(3.31)
    label("Central volume V for a 70.7 kg male at 775.762 ng/mL baseline PCSK9 (L)") # Table 2 'Central volume V (L)' 3.31
    lcl <- log(0.0061)
    label("Linear clearance CL for a 70.7 kg male (L/h)") # Table 2 'Clearance CL (L/h)' 0.0061
    lvmax <- log(0.045)
    label("Michaelis-Menten maximum elimination rate Vm (mg/h)") # Table 2 'Michaelis-Menten Vm (mg/h)' 0.045
    lkm <- fixed(log(0.2))
    label("Michaelis-Menten constant Km (mg/L)") # Table 2 'Michaelis-Menten Km (mg/L)' 0.2 (fixed)
    lq <- fixed(log(0.011))
    label("Intercompartmental clearance Q (L/h)") # Table 2 'Peripheral clearance Q (L/h)' 0.011 (fixed)
    lvp <- fixed(log(1.74))
    label("Peripheral volume V2 (L)") # Table 2 'Peripheral volume V2 (L)' 1.74 (fixed)
    lfdepot <- fixed(log(0.58))
    label("SC bioavailability F (fraction)") # Table 2 'Bioavailability F (%)' 0.58 (fixed)

    e_wt_ka <- -0.74
    label("Power exponent of WT/70.7 kg on ka (unitless)") # Table 2 'Weight effect beta WT on Ka'
    e_sexf_ka <- -0.73
    label("Exponential effect of female sex on ka (unitless)") # Table 2 'Sex effect beta SEX_F on ka'
    e_wt_vc <- 1.15
    label("Power exponent of WT/70.7 kg on V (unitless)") # Table 2 'Weight effect beta WT on V'
    e_pcsk9_vc <- 0.22
    label("Power exponent of baseline PCSK9/775.762 ng/mL on V (unitless)") # Table 2 'BPCSK9 effect beta BPCSK9 on V'
    e_sexf_vc <- 0.24
    label("Exponential effect of female sex on V (unitless)") # Table 2 'Sex effect beta SEX_F on V'
    e_wt_cl <- 1.4
    label("Power exponent of WT/70.7 kg on CL (unitless)") # Table 2 'Weight effect beta WT on CL'
    e_sexf_cl <- 0.3
    label("Exponential effect of female sex on CL (unitless)") # Table 2 'Sex effect beta SEX_F on CL'

    # LDL-C layer - Yan 2023 Table 2, 'PopPK/LDL-C final model'. Fitted with
    # the PK parameters fixed (two-step procedure).
    lrbase <- log(109.8)
    label("Typical baseline LDL-C R0 at 112.1285 mg/dL observed baseline (mg/dL)") # Table 2 'Baseline LDL-C R0 (mg/dL)' 109.8
    lkout <- log(0.006)
    label("First-order LDL-C loss rate kout at 112.1285 mg/dL baseline (1/h)") # Table 2 'Response loss rate kout (1/h)' 0.006
    lemax <- log(2.38)
    label("Maximum stimulation of kout, Emax, male non-FH (unitless)") # Table 2 'Maximum drug-induced effect Emax' 2.38
    lec50 <- fixed(log(0.54))
    label("Free tafolecimab concentration at half-maximal effect EC50 (mg/L)") # Table 2 'Typical value of EC50 (mg/L)' 0.54 (fixed)
    lhill <- fixed(log(1))
    label("Hill coefficient gamma (unitless)") # Table 2 'Hill coefficient gamma' 1 (fixed)

    e_ldlc_rbase <- 0.94
    label("Power exponent of LDLC/112.1285 mg/dL on R0 (unitless)") # Table 2 'Effect of BLDLC beta BLDLC on R0'
    e_ldlc_kout <- -1.63
    label("Power exponent of LDLC/112.1285 mg/dL on kout (unitless)") # Table 2 'Effect of BLDLC beta BLDLC on kout'
    e_dis_healthy_emax <- -0.081
    label("Exponential effect of healthy-volunteer status on Emax (unitless)") # Table 2 'Effect of DISEASE beta HV on Emax'
    e_dis_hefh_emax <- -0.37
    label("Exponential effect of HeFH on Emax (unitless)") # Table 2 'Effect of DISEASE beta HeFH on Emax'
    e_dis_hofh_emax <- -3.07
    label("Exponential effect of HoFH on Emax (unitless)") # Table 2 'Effect of DISEASE beta HoFH on Emax'
    e_sexf_emax <- -0.15
    label("Exponential effect of female sex on Emax (unitless)") # Table 2 'Effect of SEX beta F on Emax'

    # IIV. Table 2 reports the Monolix 'Standard deviation of the random
    # effects'; the variances below are those SDs squared. The Results text
    # confirms the SD scale ('omega ... on Ka ... reduced from 63% to 55%').
    etalka ~ 0.3025 # Table 2 'IIV_ka' SD 0.55
    etalvc ~ 0.1156 # Table 2 'IIV_V' SD 0.34
    etalvmax ~ 0.04 # Table 2 'IIV_Vm' SD 0.2
    etalcl ~ 0.1089 # Table 2 'IIV_Cl' SD 0.33
    etalrbase ~ 0.002116 # Table 2 'IIV_R0' SD 0.046
    etalkout ~ 0.5476 # Table 2 'IIV_kout' SD 0.74
    etalemax ~ 0.2601 # Table 2 'IIV_Emax' SD 0.51

    propSd <- 0.35
    label("Proportional residual error on tafolecimab concentration (fraction)") # Table 2 PK 'Proportional residual error' 0.35
    propSd_ldl <- 0.2
    label("Proportional residual error on LDL-C (fraction)") # Table 2 LDL-C 'Proportional residual error' 0.2
  })
  model({
    # Individual PK parameters (Results, 'Final PopPK model' equations and the
    # Table 2 note). Sex enters as exp(beta) for SEX = F (Equation 5).
    ka <- exp(lka + etalka) * (WT / 70.7)^e_wt_ka * exp(e_sexf_ka * SEXF)
    vc <- exp(lvc + etalvc) * (WT / 70.7)^e_wt_vc * (PCSK9 / 775.762)^e_pcsk9_vc * exp(e_sexf_vc * SEXF)
    cl <- exp(lcl + etalcl) * (WT / 70.7)^e_wt_cl * exp(e_sexf_cl * SEXF)
    vmax <- exp(lvmax + etalvmax)
    km <- exp(lkm)
    q <- exp(lq)
    vp <- exp(lvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Individual LDL-C parameters (Table 2 note).
    rbase <- exp(lrbase + etalrbase) * (LDLC / 112.1285)^e_ldlc_rbase
    kout <- exp(lkout + etalkout) * (LDLC / 112.1285)^e_ldlc_kout
    emax <- exp(lemax + etalemax) * exp(e_sexf_emax * SEXF + e_dis_healthy_emax * DIS_HEALTHY + e_dis_hefh_emax * DIS_HEFH + e_dis_hofh_emax * DIS_HOFH)
    ec50 <- exp(lec50)
    hill <- exp(lhill)
    # Drug-free steady state of Equation 2: R0 = kin / kout.
    kin <- kout * rbase

    # Free tafolecimab concentration L = A/V (Equation 1).
    Cc <- central / vc

    # Equation 1: two-compartment disposition, first-order SC absorption and
    # parallel linear plus Michaelis-Menten elimination (Vm in mg/h).
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - vmax * Cc / (km + Cc) - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Equation 2: type IV indirect response, drug stimulating the loss of LDL-C.
    d/dt(ldl) <- kin - kout * (1 + emax * Cc^hill / (Cc^hill + ec50^hill)) * ldl
    ldl(0) <- rbase

    # F applies to SC doses; an IV dose goes to cmt = central.
    f(depot) <- exp(lfdepot)

    Cc ~ prop(propSd)
    ldl ~ prop(propSd_ldl)
  })
}
