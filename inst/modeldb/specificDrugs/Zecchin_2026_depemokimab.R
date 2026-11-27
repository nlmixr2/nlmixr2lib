Zecchin_2026_depemokimab <- function() {
  description <- "One-compartment PK with first-order subcutaneous absorption coupled to an indirect-response model of blood eosinophil count (inhibition of production, sigmoid Imax) for depemokimab, used to bridge the mepolizumab population PK/PD meta-analysis to depemokimab in severe asthma, CRSwNP, EGPA and HES (Zecchin 2026, Table 2). Baseline BEC (KRO) depends on measured baseline BEC, disease and mepolizumab Study 200622; the residual-production fraction Imax depends on measured baseline BEC; IC50 and the between-patient and residual variability were re-estimated on the depemokimab Phase I data. Zecchin 2026 does not print its depemokimab PK estimates, so the PK layer carries the published Thorsted 2026 depemokimab population PK model evaluated for the same first-time-in-human asthma study population."
  reference <- paste(
    "Zecchin C, Schalkwijk S, Pouliquen IJ, Berges A, Bird N, Follows R, Austin D. Straight to Phase III: Model-Informed Approach Speeds Depemokimab Clinical Development in Interleukin-5-Driven Diseases. Clin Pharmacol Ther. 2026;119(4):1005. doi:10.1002/cpt.70183",
    "PK layer: Thorsted A, Tanneau L, Lavalley-Morelle A, Follows R, Jacques L, Bird N, Gevaert P, Pavord I, Ribbing J, Gupta A, Howarth P, Schalkwijk S. Population Pharmacokinetics and Pharmacodynamics of Depemokimab in People with Asthma and Chronic Rhinosinusitis with Nasal Polyps. J Clin Pharmacol. 2026;66(8):e70187. doi:10.1002/jcph.70187",
    sep = " "
  )
  vignette <- "Zecchin_2026_depemokimab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  compartmentData <- list(
    depot = list(analyte = "depemokimab", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "depemokimab", units = "mg", specimen = "plasma", verified = TRUE),
    blood_eosinophils = list(analyte = "eosinophils", units = "10^9 cells/L", specimen = "whole blood", verified = TRUE)
  )

  covariateData <- list(
    EOS = list(
      description = "Measured baseline blood eosinophil count (time-fixed per subject)",
      units = "cells/uL",
      type = "continuous",
      reference_category = NULL,
      notes = "Table 2 BEOS (GI/L), entering as LOG(BEOS/0.26), i.e. log(EOS / 260) with EOS in cells/uL. Covariate of KRO (BEOS1, exponent 0.563) and of Imax (BEOS2, exponent -0.380).",
      source_name = "BEOS"
    ),
    DIS_SEVERE_ASTHMA = list(
      description = "Severe asthma participant (1 = yes, 0 = no)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Table 2 'hasSA'; log-scale additive effect 0.438 on KRO. The reference (all disease indicators 0) is the remaining meta-analysis population (healthy participants, mild-to-moderate asthma including the depemokimab first-time-in-human cohort, and other eosinophilic conditions); Zecchin 2026 does not list it explicitly.",
      source_name = "hasSA"
    ),
    DIS_CHRONIC_RHINOSINUSITIS_NASAL_POLYPS = list(
      description = "Chronic rhinosinusitis with nasal polyps (CRSwNP) participant (1 = yes, 0 = no)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Table 2 'hasCRSwNP'; log-scale additive effect 0.132 on KRO (fixed in the source, SE 0).",
      source_name = "hasCRSwNP"
    ),
    DIS_EOSINOPHILIC_GRANULOMATOSIS_POLYANGIITIS = list(
      description = "Eosinophilic granulomatosis with polyangiitis (EGPA) participant (1 = yes, 0 = no)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Table 2 'hasEGPA'; log-scale additive effect 0.208 on KRO (fixed in the source, SE 0).",
      source_name = "hasEGPA"
    ),
    DIS_HYPEREOSINOPHILIC_SYNDROME = list(
      description = "Hypereosinophilic syndrome (HES) participant (1 = yes, 0 = no)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Table 2 'hasHES'; log-scale additive effect 1.27 on KRO.",
      source_name = "hasHES"
    ),
    STUDY_200622 = list(
      description = "Participant of the mepolizumab Phase III HES Study 200622 (1 = yes, 0 = no)",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = "Table 2 'Mepolizumab (Study 200622) HES KRO' = 0.5 (fixed), multiplying the whole KRO expression ('KROi = EXP(...)*theta200622'). The Figure 4 HES depemokimab predictions are reproduced only with this multiplier switched on (GSK 200622 is the HES reference trial), see the vignette.",
      source_name = "theta200622 (study flag)"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 2564L,
    n_studies = 17L,
    age_range = "Phase I depemokimab cohort >= 18 years; mepolizumab meta-analysis included adolescents and children (Study MEE103219, 2-17 years)",
    weight_range = NA_character_,
    sex_female_pct = NA_real_,
    race_ethnicity = NA_character_,
    disease_state = "Mepolizumab meta-analysis: healthy participants and patients with asthma (mild-to-moderate and severe), CRSwNP, EGPA, HES and eosinophilic oesophagitis; depemokimab re-estimation: mild-to-moderate asthma with BEC >= 200 cells/uL",
    dose_range = "Depemokimab single SC 2-300 mg (Phase I); mepolizumab IV and SC regimens across 16 studies",
    regions = "Multinational",
    notes = "The mepolizumab population PK/PD meta-analysis used 29,575 BEC measurements from 2516 participants in 16 studies (including 205 CRSwNP, 135 EGPA and 192 HES patients); the depemokimab layer used 675 BEC measurements from 48 patients in the first-time-in-human Study 205722 (Zecchin 2026 Results). Only IC50, the two between-patient variances and the residual SD were re-estimated on the depemokimab data; every other PD parameter, the covariate effects and the model structure were fixed to the mepolizumab estimates (Table 2 grey rows; Supplementary Methods). The PK layer (Thorsted 2026, Table S3) was fitted to 961 participants from three Phase I and four Phase III depemokimab studies; it is evaluated here for the first-time-in-human asthma study population at the Thorsted reference covariates (76 kg, 50 years, non-Asian, albumin 45 g/L, eGFR 94 mL/min/1.73 m2, upper-arm injection)."
  )

  ini({
    # ---- PK layer: Thorsted 2026 Table S3, FTIH Asthma study population ----
    # Zecchin 2026 re-estimated a one-compartment first-order-absorption PK
    # model on Study 205722 but does not report the estimates. Thorsted 2026
    # categorical covariates are fractional changes, TV = theta * (1 + effect)
    # (Thorsted 2026 Supplementary Methods); evaluated at the FTIH Asthma study
    # population and the Thorsted reference subject. Held fixed: no estimate
    # from Zecchin 2026 exists to replace them.
    lcl <- fixed(log(0.0920 * (1 + 0.155))); label("Apparent clearance CL/F, FTIH asthma population (L/day)") # Thorsted 2026 Table S3: CL/F 0.0920 L/day x (1 + 'FTIH Asthma study population on CL' 0.155) = 0.10626
    lvc <- fixed(log(6.29)); label("Apparent volume V/F (L)") # Thorsted 2026 Table S3: V/F 6.29 L (non-Asian patient reference)
    lka <- fixed(log(0.212 * (1 + 0.627))); label("First-order SC absorption rate ka, FTIH asthma population (1/day)") # Thorsted 2026 Table S3: ka 0.212 /day x (1 + 'FTIH Asthma study population on ka' 0.627) = 0.34492
    lfdepot <- fixed(log(1 - 0.173)); label("Relative SC bioavailability Frel, FTIH asthma population (fraction)") # Thorsted 2026 Table S3: Frel 1.00 FIX x (1 + 'FTIH Asthma study population on Frel' -0.173) = 0.827

    # ---- PD layer: Zecchin 2026 Table 2 ----
    # Thetas are on the log scale except GAMA and the covariate effects.
    lrbase <- fixed(-1.73); label("Baseline BEC KRO at BEOS = 0.26 GI/L, reference disease (GI/L)") # Table 2 'KRO (GI/L)' theta -1.73 (SE 0.0529) -> 0.177 GI/L; grey row (mepolizumab estimate, fixed)
    lkout <- fixed(-4.35 + log(24)); label("First-order BEC elimination rate KOUT (1/day)") # Table 2 'KOUT (1/h)' theta -4.35 (SE 0.0921) -> 0.0129 /h; KOUTi = 24*EXP(theta) = 0.310 /day; grey row (fixed)
    lic50 <- -2.34; label("Depemokimab concentration at half-maximal inhibition IC50 (ug/mL)") # Table 2 'IC50 (ug/mL)' theta -2.34 (SE 0.0789) -> 0.0960 ug/mL; estimated on depemokimab Phase I data
    limax <- fixed(-1.77); label("Residual fraction of BEC production at maximal inhibition Imax, BEOS = 0.26 GI/L (fraction)") # Table 2 'I max' theta -1.77 (SE 0.0265) -> 0.170; grey row (fixed)
    lhill <- fixed(log(1.61)); label("Hill coefficient of the inhibition GAMA (unitless)") # Table 2 'GAMA' theta 1.61 (SE 0.0559), untransformed; grey row and Supplementary Methods ('fixed to the value estimated for mepolizumab')
    e_eos_imax <- fixed(-0.380); label("Exponent of BEOS/0.26 on Imax, BEOS2 (unitless)") # Table 2 'BEOS2 on I max' -0.380 (SE 0.0339); grey row (fixed)
    e_eos_rbase <- fixed(0.563); label("Exponent of BEOS/0.26 on KRO, BEOS1 (unitless)") # Table 2 'BEOS1 on KRO' 0.563 (SE 0.0226); grey row (fixed)
    e_dis_severe_asthma_rbase <- fixed(0.438); label("Log-scale additive effect of severe asthma on KRO (unitless)") # Table 2 'SA KRO' 0.438 (SE 0.0546); grey row (fixed)
    e_dis_chronic_rhinosinusitis_nasal_polyps_rbase <- fixed(0.132); label("Log-scale additive effect of CRSwNP on KRO (unitless)") # Table 2 'CRSwNP KRO' 0.132 (SE 0); grey row (fixed)
    e_dis_eosinophilic_granulomatosis_polyangiitis_rbase <- fixed(0.208); label("Log-scale additive effect of EGPA on KRO (unitless)") # Table 2 'EGPA KRO' 0.208 (SE 0); grey row (fixed)
    e_dis_hypereosinophilic_syndrome_rbase <- fixed(1.27); label("Log-scale additive effect of HES on KRO (unitless)") # Table 2 'HES KRO' 1.27 (SE 0.103); grey row (fixed)
    e_study_200622_rbase <- fixed(0.5); label("Multiplier on KRO for mepolizumab Study 200622 participants (unitless)") # Table 2 'Mepolizumab (Study 200622) HES KRO' 0.5 (SE 0); grey row (fixed)

    # Exponential between-patient variability (Table 2 footnote c: %CV = sqrt(exp(OMEGA) - 1)).
    etalrbase ~ 0.0853 # Table 2 'BPV KRO' OMEGA2 0.0853 (SE 0.0161), 29.8% CV
    etalimax ~ 0.294 # Table 2 'BPV I max' OMEGA2 0.294 (SE 0.0866), 58.5% CV

    # Y = IPRED * exp(eps), eps ~ N(0, SD^2): additive on the log scale.
    expSd <- 0.350; label("Residual SD of log BEC (log GI/L)") # Table 2 'Residual error' SD 0.350 (SE 0.0203)
  })

  model({
    # PK (Thorsted 2026 structure; FTIH asthma population)
    ka <- exp(lka)
    cl <- exp(lcl)
    vc <- exp(lvc)
    kel <- cl / vc

    # PD (Zecchin 2026 Table 2 'PD parameter estimation')
    leos_rel <- log(EOS / 260)
    study_mult <- 1 + (e_study_200622_rbase - 1) * STUDY_200622
    rbase <- exp(lrbase + e_eos_rbase * leos_rel +
      e_dis_severe_asthma_rbase * DIS_SEVERE_ASTHMA +
      e_dis_hypereosinophilic_syndrome_rbase * DIS_HYPEREOSINOPHILIC_SYNDROME +
      e_dis_eosinophilic_granulomatosis_polyangiitis_rbase * DIS_EOSINOPHILIC_GRANULOMATOSIS_POLYANGIITIS +
      e_dis_chronic_rhinosinusitis_nasal_polyps_rbase * DIS_CHRONIC_RHINOSINUSITIS_NASAL_POLYPS +
      etalrbase) * study_mult
    kout <- exp(lkout)
    ic50 <- exp(lic50)
    imax <- exp(limax + e_eos_imax * leos_rel + etalimax)
    hill <- exp(lhill)

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central
    f(depot) <- exp(lfdepot)

    Cc <- central / vc

    # Zecchin 2026 equation: deos/dt = KOUT * (KRO * (1 + (Imax - 1) * C^GAMA / (C^GAMA + IC50^GAMA)) - eos)
    inh <- Cc^hill / (Cc^hill + ic50^hill)
    d/dt(blood_eosinophils) <- kout * (rbase * (1 + (imax - 1) * inh) - blood_eosinophils)
    blood_eosinophils(0) <- rbase

    blood_eosinophils ~ lnorm(expSd)
  })
}
