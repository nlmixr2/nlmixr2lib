Kim_2024_meropenem <- function() {
  description <- paste(
    "Two-compartment intravenous population PK model for meropenem in",
    "healthy Korean adults with normal renal function (Kim 2024; n = 12,",
    "84 plasma samples after a single 500 mg 30-min IV infusion).",
    "Total clearance depends on serum creatinine through a power model",
    "centred at the cohort median 0.86 mg/dL:",
    "CL = 12.4 * (CREAT/0.86)^-0.392 L/h, so CL falls from 14.3 L/h at",
    "CREAT = 0.6 mg/dL to 11.7 L/h at CREAT = 1.0 mg/dL. Because the",
    "estimated correlation between the random effects on CL and Vc was",
    "0.99, the paper codes a SINGLE eta shared by both parameters, with",
    "the Vc random effect constructed as 1.53 * eta_CL (Table 2 'V1'",
    "row under interindividual variability). Interindividual variability",
    "on Q and Vp is fixed, not estimated. Residual error is",
    "proportional. The unbound concentration Cu = fu * Cc (fu = 0.98)",
    "is the driver of the fT>MIC PK/PD targets (40% fT>MIC, 40%",
    "fT>4MIC, 100% fT>MIC, 100% fT>4MIC) used in the paper's Monte",
    "Carlo probability-of-target-attainment simulations for",
    "intermittent, 3-h extended, and continuous infusion regimens."
  )
  reference <- paste(
    "Kim YK, Kang G, Zang DY, Lee DH. (2024). Precision Dosing of",
    "Meropenem in Adults with Normal Renal Function: Insights from a",
    "Population Pharmacokinetic and Monte Carlo Simulation Study.",
    "Antibiotics 13(9):849. doi:10.3390/antibiotics13090849.",
    sep = " "
  )
  vignette <- "Kim_2024_meropenem"
  units    <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. Verified against Kim 2024 Methods 4.2 (plasma
  # meropenem measured by LC-MS/MS) and 4.3 (ADVAN3 TRANS4 two-compartment
  # model with central and one peripheral compartment).
  compartmentData <- list(
    central     = list(analyte = "meropenem", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "meropenem", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CREAT = list(
      description        = "Serum creatinine concentration",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "The only covariate retained in the Kim 2024 final model.",
        "Enters CL as a power term centred at 0.86 mg/dL, the cohort",
        "median (Table 1: mean 0.863, CV 19.0%; median 0.860, IQR",
        "0.738-1.02 mg/dL). The structural-model row of Table 2 gives",
        "CL = theta1 * (CR/0.86)^theta2 with theta1 = 12.4 L/h and",
        "theta2 = -0.392, so clearance DECREASES as serum creatinine",
        "rises. The Discussion confirms both anchor points: CL = 14.3",
        "L/h at CREAT = 0.6 mg/dL and 11.7 L/h at CREAT = 1.0 mg/dL.",
        "Note the units are mg/dL, not the umol/L form used by several",
        "other models in this registry (1 mg/dL = 88.42 umol/L).",
        "Time-fixed per subject in the source dataset (a single",
        "screening chemistry panel per healthy volunteer)."
      ),
      source_name        = "CR"
    )
  )

  # Kim 2024 Methods 4.3 screened a wide covariate set by stepwise forward
  # selection (p < 0.01) and backward elimination (p < 0.001); only serum
  # creatinine survived. The Discussion attributes this to the small sample
  # (n = 12) and the narrow renal-function range of a healthy cohort. These
  # covariates are recorded for provenance but are NOT in the final model,
  # and no point estimates are published for them (Table S2 holds only the
  # stepwise OFV trace).
  covariatesDataExcluded <- list(
    SEXF = list(
      description        = "Female sex indicator (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Screened as 'gender' in the Kim 2024 stepwise covariate search (Methods 4.3); not retained. Cohort was 4 female / 8 male (Results 2.1).",
      source_name        = "SEX"
    ),
    AGE = list(
      description        = "Age at screening",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened in the Kim 2024 stepwise covariate search; not retained. Table 1: mean 36.8 years (CV 19.9%), median 36.0 (IQR 31.5-39.3); inclusion criterion 19-55 years.",
      source_name        = "AGE"
    ),
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened in the Kim 2024 stepwise covariate search; not retained, so NO allometric scaling is present in this model. Table 1: mean 65.7 kg (CV 20.8%), median 61.7 (IQR 56.7-73.3). Table 3 reports per-kilogram post-hoc parameters only as descriptive statistics, not as a structural weight relationship.",
      source_name        = "WT"
    ),
    HT = list(
      description        = "Body height",
      units              = "cm",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened in the Kim 2024 stepwise covariate search; not retained. Table 1: mean 168 cm (CV 4.29%), median 168 (IQR 163-173).",
      source_name        = "HT"
    ),
    BMI = list(
      description        = "Body mass index",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened in the Kim 2024 stepwise covariate search; not retained. Table 1: mean 23.0 kg/m^2 (CV 13.9%), median 21.5 (IQR 21.4-24.0).",
      source_name        = "BMI"
    ),
    BSA = list(
      description        = "Body surface area",
      units              = "m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened in the Kim 2024 stepwise covariate search; not retained. Table 1: mean 1.74 m^2 (CV 11.6%), median 1.71 (IQR 1.61-1.88). Also used to normalize creatinine clearance to mL/min/1.73 m^2.",
      source_name        = "BSA"
    ),
    ALB = list(
      description        = "Serum albumin concentration",
      units              = "g/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened in the Kim 2024 stepwise covariate search; not retained. Table 1: mean 4.88 g/dL (CV 4.35%), median 4.80 (IQR 4.78-5.03).",
      source_name        = "ALB"
    ),
    CYSC = list(
      description        = "Serum cystatin C concentration",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened in the Kim 2024 stepwise covariate search; not retained. Table 1 reports it in mg/dL (mean 0.790, CV 15.9%; median 0.765, IQR 0.705-0.873), which is a factor of 10 below the mg/L units used elsewhere in the registry for this analyte.",
      source_name        = "CYSC"
    ),
    CRCL = list(
      description        = "Renal function estimated as creatinine clearance or eGFR",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Kim 2024 screened Cockcroft-Gault creatinine clearance, MDRD eGFR, and three CKD-EPI eGFR variants (creatinine, cystatin C, combined), each in both raw and BSA-normalized form (Table 1 footnotes a-g); none was retained in the final model, and only the raw serum creatinine concentration reached significance. The Discussion attributes this to the narrow renal-function range of a healthy cohort (normalized Cockcroft-Gault CLCR mean 93.8 mL/min/1.73 m^2, CV 16.4%).",
      source_name        = "CLCR"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 12,
    n_studies      = 1,
    n_observations = 84,
    age_range      = "19-55 years by protocol; observed mean 36.8 years (CV 19.9%), median 36.0 (IQR 31.5-39.3) (Table 1).",
    weight_range   = "Mean 65.7 kg (CV 20.8%), median 61.7 kg (IQR 56.7-73.3) (Table 1).",
    sex_female_pct = "33% (4 of 12 participants female).",
    race_ethnicity = "Korean. The paper states this is the first population PK study of meropenem in healthy Korean adults.",
    disease_state  = "Healthy volunteers with normal renal and hepatic function; no clinically significant laboratory abnormalities and no adverse events reported. Model was then applied by Monte Carlo simulation to patients with normal renal function (CLCR > 50 mL/min) and infections caused by Pseudomonas aeruginosa.",
    renal_function = "Normal. Cockcroft-Gault CLCR mean 105 mL/min (CV 21.2%); BSA-normalized CLCR mean 93.8 mL/min/1.73 m^2 (CV 16.4%); CKD-EPI creatinine eGFR mean 111 mL/min/1.73 m^2 (Table 1).",
    dose_range     = "Estimation data: a single 500 mg dose of meropenem in 100 mL saline infused intravenously over 30 min. Simulated regimens: 0.5 / 1 / 1.5 / 2 g at 0.5-h and 3-h infusion durations every 6, 8 or 12 h, and continuous infusions of 2, 4, 6 and 8 g/day.",
    regions        = "Republic of Korea (Clinical Trial Center, Hallym University Sacred Heart Hospital, Anyang), January 2023.",
    notes          = paste(
      "14 volunteers consented; 2 were excluded for positive allergic",
      "skin tests, leaving 12 analysed. Sampling was pre-dose (0 h) and",
      "at 0.5, 0.75, 1, 2, 3 and 6 h from the start of the 30-min",
      "infusion; plasma meropenem was measured by validated LC-MS/MS.",
      "Estimation used NONMEM 7.5 FOCE with interaction, subroutine",
      "ADVAN3 TRANS4; model evaluation used PsN 5.3.1 (VPC with 1000",
      "simulations, 2000-sample nonparametric bootstrap). A",
      "three-compartment model gave a lower OFV (-160.7 vs -150.2) but",
      "was rejected as unphysiological (CL 0.0138 L/h, second",
      "peripheral volume 4404 L, RSEs not estimable); the",
      "two-compartment model was retained. Final-model OFV -175.610."
    )
  )

  ini({
    # ---- Structural parameters (Kim 2024 Table 2, final model) ----
    # Reported as linear-scale typical values; log-transformed here.
    lcl <- log(12.4); label("Clearance at the reference serum creatinine of 0.86 mg/dL (L/h)")  # Table 2 theta1 = 12.4 L/h (RSE 7.87%; bootstrap 12.3, 95% CI 10.8-14.7)
    lvc <- log(8.26); label("Central volume of distribution (L)")                               # Table 2 V1 = 8.26 L (RSE 12.5%; bootstrap 8.31, 95% CI 6.56-11.0)
    lq  <- log(5.22); label("Intercompartmental clearance (L/h)")                                # Table 2 Q = 5.22 L/h (RSE 16.1%; bootstrap 5.05, 95% CI 3.34-7.33)
    lvp <- log(4.06); label("Peripheral volume of distribution (L)")                             # Table 2 V2 = 4.06 L (RSE 11.1%; bootstrap 4.01, 95% CI 3.00-5.07)

    # ---- Covariate effect (Kim 2024 Table 2 structural-model row) ----
    # CL = theta1 * (CR / 0.86)^theta2, CR in mg/dL, reference 0.86 = the
    # cohort median serum creatinine (Table 1). theta2 is NEGATIVE, so
    # clearance falls as serum creatinine rises. Both anchor points printed
    # in the Discussion reproduce exactly:
    #   CREAT = 0.6 -> 12.4 * (0.6/0.86)^-0.392 = 14.28 L/h  (paper: 14.3)
    #   CREAT = 1.0 -> 12.4 * (1.0/0.86)^-0.392 = 11.69 L/h  (paper: 11.7)
    e_creat_cl <- -0.392; label("Power exponent of serum creatinine on CL (unitless)")           # Table 2 theta2 = -0.392 (RSE 19.2%; bootstrap -0.378, 95% CI -0.579 to -0.115)

    # ---- Plasma protein binding ----
    # Kim 2024 Methods 4.5 "Dosage Simulation": "The fraction of unbound
    # drug (f) was fixed at 0.98." Not estimated from the PK data, hence
    # fixed(). fu converts the total plasma concentration Cc to the free
    # concentration Cu that drives the fT>MIC PK/PD targets.
    fu <- fixed(0.98); label("Fraction unbound in plasma (unitless)")                            # Methods 4.5, Dosage Simulation

    # ---- Inter-individual variability (Kim 2024 Table 2) ----
    # Methods 4.3: theta_i = theta * exp(eta_i) with eta ~ N(0, omega^2).
    # Table 2 reports the IIV magnitudes as percentages, converted to the
    # internal log-scale variance via omega^2 = log(CV^2 + 1).
    #
    # Results 2.2: the correlation between the CL and V1 random effects was
    # 0.99, so the paper codes ONE eta shared by both parameters:
    #   CL = TVCL * EXP[ETA(1)]  and  V1 = TVV1 * EXP[THETA(1) * ETA(1)]
    # with THETA(1) = 1.53 in the final model. The Table 2 footnote states
    # "the estimate indicates that the interindividual variability (IIV) of
    # V1 is 1.53 times the IIV of CL". This is encoded as scale_etalvc
    # below and applied in model(); no independent eta on Vc exists.
    etalcl ~ 0.0663906  # Table 2: IIV on CL = 26.2% (RSE 30.4%, shrinkage 1.82%; bootstrap 25.4%, 95% CI 7.8-38.3) -> log(1 + 0.262^2)
    etalq  ~ fixed(0.0205239)  # Table 2: IIV on Q = 14.4% (footnote b: fixed; shrinkage 49.1%) -> log(1 + 0.144^2)
    etalvp ~ fixed(0.0315384)  # Table 2: IIV on V2 = 17.9% (footnote b: fixed; shrinkage 11.2%) -> log(1 + 0.179^2)

    scale_etalvc <- 1.53; label("Scaling of the CL inter-individual eta onto Vc (unitless)")     # Table 2 IIV block, V1 row = 1.53 (RSE 4.80%; bootstrap 1.53, 95% CI 1.07-2.43); NONMEM THETA(1)

    # ---- Residual error (Kim 2024 Table 2) ----
    propSd <- 0.109; label("Proportional residual error (fraction)")                             # Table 2: proportional error 10.9% (RSE 20.4%; bootstrap 10.4%, 95% CI 6.50-14.0)
  })
  model({
    # 1. Individual PK parameters. Serum creatinine acts on CL only. The Vc
    #    random effect is the CL eta scaled by scale_etalvc (perfect
    #    correlation, matching the paper's shared-eta coding); Q and Vp
    #    carry their own etas whose variances were fixed.
    cl <- exp(lcl + etalcl) * (CREAT / 0.86)^e_creat_cl
    vc <- exp(lvc + scale_etalvc * etalcl)
    q  <- exp(lq  + etalq)
    vp <- exp(lvp + etalvp)

    # 2. Micro-constants
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # 3. ODE system. Meropenem is given intravenously only (30-min infusion
    #    in the estimation study; 0.5-h / 3-h intermittent or continuous
    #    infusions in the simulations), so doses go straight to `central`.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # 4. Observations. Cc is total plasma meropenem (mg/L; dose in mg,
    #    volume in L). Cu is the free concentration compared against the MIC
    #    for the fT>MIC targets; it is a deterministic transform of Cc and is
    #    not separately measured, so it carries no residual error.
    Cc <- central / vc
    Cu <- fu * Cc

    Cc ~ prop(propSd)
  })
}
