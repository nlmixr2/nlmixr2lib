Troisi_2024_meropenem <- function() {
  description <- paste(
    "Sequential population PK/PD model of continuous-infusion meropenem and plasma C-reactive",
    "protein (C-RP) in critically ill adults with documented Gram-negative hospital-acquired or",
    "ventilator-associated pneumonia. One-compartment intravenous PK with the volume fixed at 20 L",
    "and clearance rising exponentially (uncentered) with CKD-EPI eGFR. Meropenem concentration",
    "inhibits C-RP production in an indirect-response (turnover) model with FULL inhibition",
    "(Imax = 1); the C-RP baseline R0 is estimated and kin = R0 * kout. The IC50 rises",
    "exponentially (uncentered) with the meropenem MIC of the subject's isolate and with",
    "empirical anti-MRSA co-therapy.",
    sep = " "
  )
  reference <- paste(
    "Troisi C, Cojutti PG, Rinaldi M, Tonetti T, Siniscalchi A, van Hasselt C, Viale P, Pea F.",
    "Impact of Continuous Infusion Meropenem PK/PD Target Attainment on C-Reactive Protein",
    "Dynamics in Critically Ill Patients With Documented Gram-Negative Hospital-Acquired or",
    "Ventilator-Associated Pneumonia. Clin Pharmacokinet. 2024;63(11):1573-1583.",
    "doi:10.1007/s40262-024-01436-6",
    sep = " "
  )
  vignette <- "Troisi_2024_meropenem"

  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "meropenem", units = "mg", specimen = "plasma", verified = TRUE),
    # The C-RP turnover state holds a CONCENTRATION (mg/dL): the paper's
    # equation is written directly on C-RP with no volume term, so kin carries
    # units of mg/dL per h.
    crp = list(analyte = "C-reactive protein", units = "mg/dL", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = "CKD-EPI estimated glomerular filtration rate, BSA-normalized",
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Source column eGFR, estimated with the CKD-EPI formula (Patients and Methods) and",
        "reported in mL/min/1.73 m^2 (Table 1 median 70, IQR 39-100). Table 3 labels the",
        "coefficient 'CLCR on CLMER' but Results, the ESM covariate equation and the Discussion",
        "all name the covariate eGFR. The effect is EXPONENTIAL and UNCENTERED: Results print",
        "'CL (L/h) = 2.72 x e^(0.011 x eGFR)', so exp(lcl) is the clearance extrapolated to",
        "eGFR = 0. Closed by the Discussion: 4.66 L/h at eGFR 49 and 7.44 L/h at eGFR 91.5."
      ),
      source_name = "eGFR"
    ),
    MIC = list(
      description = "Meropenem minimum inhibitory concentration against the subject's Gram-negative clinical isolate",
      units = "mg/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-fixed per subject. Enters the IC50 as an UNCENTERED exponential,",
        "IC50 = 0.79 x e^(0.44 x MIC) x e^(2.46 x DRUG) (Results; ESM covariate equation 2).",
        "Table 2: overall median 0.25 mg/L (IQR 0.12-4.00, range 0.12-64) across 47 isolates",
        "(Pseudomonas aeruginosa 13, Klebsiella pneumoniae 13, Acinetobacter baumannii 7, other",
        "Enterobacterales 14). Many isolates sit at 0.12 mg/L, the lowest dilution reported. The",
        "uncentered exponential makes the IC50 very large at high MIC (about 1100-fold at 16 mg/L),",
        "which is how the model expresses a near-absent C-RP response to a resistant isolate."
      ),
      source_name = "MIC"
    ),
    CONMED_ANTIMRSA = list(
      description = "Empirical anti-MRSA (anti-Gram-positive) agent co-administered with meropenem",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = meropenem monotherapy (Fig. 2 Group A)",
      notes = paste(
        "Source column DRUG (ESM covariate equation 2; Table 3 'Group on IC50'). 1 = Fig. 2",
        "Group B, 'meropenem + an anti-Gram positive agent' given EMPIRICALLY without a documented",
        "Gram-positive pathogen (23/47 of the PD cohort, Table 1); 0 = Group A, meropenem",
        "monotherapy (24/47). The paper does not name the anti-MRSA agents. Patients on a second",
        "anti-Gram-negative agent (Group C) or on TARGETED anti-MRSA therapy (Group D) were in the",
        "PK fit but excluded from the PD fit, so the model is not defined for them. Effect:",
        "exp(2.46) = 11.7-fold higher IC50."
      ),
      source_name = "DRUG"
    )
  )

  # Screened on the PK parameters (Patients and Methods: "age, sex, weight,
  # height, serum creatinine and eGFR") and not retained; Results keep only eGFR
  # on CL. Documentation only -- none is referenced in model().
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened on the PK parameters but not retained. Table 1 median 65.5 years (IQR 57.0-74.0)."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "binary",
      notes = "Screened on the PK parameters but not retained. Table 1 M/F 36/28 (56.25%/43.75%)."
    ),
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened on the PK parameters but not retained; no allometric scaling. Table 1 median 75.0 kg (IQR 65.0-85.0)."
    ),
    HT = list(
      description = "Height",
      units = "cm",
      type = "continuous",
      notes = "Screened on the PK parameters but not retained. Table 1 median 1.70 m (IQR 1.68-1.75)."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "mg/dL",
      type = "continuous",
      notes = "Screened on the PK parameters but not retained; renal function entered through CKD-EPI eGFR. Table 1 median 0.96 mg/dL (IQR 0.62-1.31)."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 64L,
    n_studies = 1L,
    n_observations = "211 meropenem steady-state plasma concentrations (PK, 64 patients) and 328 C-RP concentrations (PD, 47 patients)",
    age_median = "65.5 years (IQR 57.0-74.0)",
    weight_median = "75.0 kg (IQR 65.0-85.0)",
    sex_female_pct = 43.75,
    race_ethnicity = "Not reported; single-centre Italian cohort.",
    disease_state = paste(
      "Critically ill adults with documented Gram-negative hospital-acquired (HAP) or",
      "ventilator-associated pneumonia (VAP) treated with continuous-infusion meropenem under",
      "real-time TDM. PD cohort (n = 47, Table 2): VAP 35, HAP 6, VAP + bloodstream infection 5,",
      "HAP + bloodstream infection 1; isolates P. aeruginosa 13, K. pneumoniae 13, A. baumannii 7,",
      "other Enterobacterales 14; median SOFA 10, APACHE 17, mechanical ventilation 95.7%, shock",
      "53.2% (Table 1). Isolate type (P. aeruginosa vs other) was screened on the PD parameters and",
      "not retained."
    ),
    renal_function = "CKD-EPI eGFR median 70.0 mL/min/1.73 m^2 (IQR 39.0-100.0), Table 1.",
    baseline_crp = "Observed baseline C-RP median 15.7 mg/dL (IQR 9.2-26.1) in the PD cohort (Table 1); estimated typical C-RP0 22.58 mg/dL (Table 3).",
    dose_range = paste(
      "Loading dose 2 g over 2 h, then continuous infusion started at 1 g q6h (eGFR >= 60) or",
      "0.5 g q6h (eGFR < 60), each bag infused over 6-8 h, then TDM-adjusted to a Css/MIC target",
      "of 4-8. Median daily dose 2 g (IQR 1.5-4); median treatment 9 days (IQR 6.5-14); median",
      "Css 15.9 mg/L (IQR 9.8-27.6)."
    ),
    regions = "Italy (IRCCS Azienda Ospedaliero Universitaria di Bologna)",
    notes = paste(
      "Retrospective single-centre study, December 2020 to August 2023. Treatment groups (Fig. 2):",
      "A meropenem monotherapy (n = 24), B meropenem + empirical anti-MRSA agent (n = 23), C",
      "meropenem + another anti-Gram-negative agent, D meropenem + targeted anti-MRSA therapy. The",
      "PK model used all 64 patients; the PD model only groups A and B (n = 47). Two-step",
      "sequential fit in Monolix 2023R1 (SAEM): the PK parameters were estimated first and then",
      "fixed during the PD fit. Meropenem was measured by LC-MS/MS."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # PK (Troisi 2024 Table 3, 'PK model'). Estimated in the first step of
    # the sequential fit and held at these values while the PD part was fit.
    # V was fixed from the group's earlier meropenem model (ref 16,
    # modellib('Cojutti_2021_meropenem')) because continuous infusion made
    # it inestimable.
    # ------------------------------------------------------------------
    lcl <- log(2.72); label("Meropenem clearance intercept at CRCL = 0 (L/h)")                       # Table 3: CLMER = 2.72 (RSE 16.7%; bootstrap 2.32-3.13)
    lvc <- fixed(log(20)); label("Meropenem volume of distribution (L)")                            # Table 3: V = 20 L (fixed); Methods: fixed from ref 16
    e_crcl_cl <- 0.011; label("Exponential coefficient of CRCL on CL (per mL/min/1.73 m^2; uncentered)") # Table 3: 'CLCR on CLMER' = 0.011 (RSE 20.0%); Results: CL = 2.72 x e^(0.011 x eGFR)

    # ------------------------------------------------------------------
    # PD: indirect response, full inhibition of C-RP production
    # (Table 3, 'PD model'; equations in Methods and the ESM).
    # ------------------------------------------------------------------
    lrbase <- log(22.58); label("Baseline plasma C-RP, C-RP0 (mg/dL)")                              # Table 3: C-RP0 = 22.58 mg/dL (RSE 9.55%)
    lic50 <- log(0.79); label("Meropenem concentration for half-maximal inhibition of C-RP production at MIC = 0, no anti-MRSA co-therapy (mg/L)") # Table 3: IC50 = 0.79 mg/L (RSE 53.4%)
    lkout <- log(0.012); label("First-order C-RP degradation rate constant (1/h)")                  # Table 3: kout = 0.012 1/h (RSE 27.2%)
    e_mic_ic50 <- 0.44; label("Exponential coefficient of MIC on IC50 (per mg/L; uncentered)")     # Table 3: 'MIC on IC50' = 0.44 (RSE 15.9%)
    e_conmed_antimrsa_ic50 <- 2.46; label("Exponential coefficient of empirical anti-MRSA co-therapy on IC50 (unitless)") # Table 3: 'Group on IC50' = 2.46 (RSE 31.9%)

    # ------------------------------------------------------------------
    # IIV. Table 3 footnote: 'omega standard deviation of inter-individual
    # variability' (Monolix omegas, SD of the log-scale random effect), so
    # the nlmixr2 variance is omega^2. No correlations were retained.
    # ------------------------------------------------------------------
    etalcl ~ 0.3481    # 0.59^2; Table 3 omega CLMER = 0.59 (RSE 10.4%)
    etalrbase ~ 0.2601 # 0.51^2; Table 3 omega C-RP0 = 0.51 (RSE 15.8%)
    etalic50 ~ 2.5281  # 1.59^2; Table 3 omega IC50 = 1.59 (RSE 22.5%)
    etalkout ~ 1.1664  # 1.08^2; Table 3 omega kout = 1.08 (RSE 21.6%)

    # ------------------------------------------------------------------
    # Residual error. Table 3 footnote: observation = prediction x
    # (1 + proportional error), i.e. proportional on each output.
    # ------------------------------------------------------------------
    propSd <- 0.37; label("Proportional residual error on meropenem plasma concentration (fraction)") # Table 3 PK: proportional error = 0.37 (RSE 6.53%)
    propSd_crp <- 0.35; label("Proportional residual error on plasma C-RP (fraction)")               # Table 3 PD: proportional error = 0.35 (RSE 4.87%)
  })

  model({
    # 1. Individual PK parameters. eGFR effect is exponential and uncentered
    #    (Results: CL = 2.72 x e^(0.011 x eGFR)); V carries no IIV.
    cl <- exp(lcl + etalcl) * exp(e_crcl_cl * CRCL)
    vc <- exp(lvc)
    kel <- cl / vc

    # 2. Individual PD parameters. IC50 covariate equation (Results; ESM 2):
    #    IC50 = IC50pop x e^(beta_MIC x MIC) x e^(beta_DRUG x DRUG).
    rbase <- exp(lrbase + etalrbase)
    kout <- exp(lkout + etalkout)
    ic50 <- exp(lic50 + etalic50) * exp(e_mic_ic50 * MIC) *
      exp(e_conmed_antimrsa_ic50 * CONMED_ANTIMRSA)
    kin <- rbase * kout

    # 3. Meropenem plasma concentration (mg/L; dose mg, volume L). This
    #    drives C-RP production inhibition directly (no effect compartment).
    Cc <- central / vc

    # 4. ODEs. C-RP equation as in Methods / ESM:
    #    dCRP/dt = kin x (1 - C / (C + IC50)) - kout x CRP, CRP(0) = C-RP0.
    d/dt(central) <- -kel * central
    d/dt(crp) <- kin * (1 - Cc / (Cc + ic50)) - kout * crp
    crp(0) <- rbase

    # 5. Observations: meropenem (mg/L) and C-RP (mg/dL).
    Cc ~ prop(propSd)
    crp ~ prop(propSd_crp)
  })
}
