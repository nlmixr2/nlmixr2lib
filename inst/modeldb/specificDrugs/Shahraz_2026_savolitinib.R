Shahraz_2026_savolitinib <- function() {
  description <- paste(
    "One-compartment population PK model for oral savolitinib in adults",
    "with locally advanced or metastatic solid tumours, given alone or with",
    "osimertinib (Shahraz 2026, pooled analysis of 10 studies). Absorption",
    "runs through a depot and two transit compartments that share the",
    "first-order rate constant KA; elimination is first order. Covariates:",
    "baseline body weight and alkaline phosphatase on V/F; baseline total",
    "bilirubin, alkaline phosphatase and osimertinib co-administration on",
    "CL/F; proton pump inhibitor use on V/F and KA; female sex on CL/F and",
    "V/F; and White race (against an Asian, Black and Other reference) on",
    "CL/F, V/F and KA."
  )
  reference <- paste(
    "Shahraz A, Derippe T, Kode KK, Allmond J, Nelson S, Igwegbe I,",
    "Sawant-Basak A, Vishwanathan K, Miah K. Population Pharmacokinetics and",
    "Exposure-Response Analyses of Savolitinib Alone or in Combination With",
    "Osimertinib in Patients With Locally Advanced or Metastatic Solid",
    "Tumors. Clin Transl Sci. 2026;19:e70700. doi:10.1111/cts.70700.",
    "PMCID PMC13504015.",
    sep = " "
  )
  vignette <- "Shahraz_2026_savolitinib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description = "Baseline body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effect on V/F, reference 62.2 kg (the median of the pooled dataset, Table S2); Data S1 equation for Vd/F cont,cov.",
      source_name = "WT"
    ),
    ALP = list(
      description = "Baseline serum alkaline phosphatase",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Baseline value (source column BALP), time-fixed. Power effects on V/F and CL/F, reference 89 U/L (the pooled median, Table S2); Data S1 equations. The source abbreviation list expands BALP as 'baseline of alanine alkaline phosphatase'; Table S2 and Figure S7 call it baseline alkaline phosphatase.",
      source_name = "BALP"
    ),
    TBILI = list(
      description = "Baseline total bilirubin",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Baseline value (source column BBILI), time-fixed, in umol/L (Table S2). Power effect on CL/F, reference 8.5 umol/L; Data S1 equation for CL/F cont,cov.",
      source_name = "BBILI"
    ),
    SEXF = list(
      description = "Female sex indicator (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Source column SEX with SEX = 1 the effect level. Data S1 does not name the level; Figure 2 labels it: females have an AUC24,ss ratio of 1.24 = 1 / (1 - 0.197), so SEX = 1 is female and SEXF = SEX.",
      source_name = "SEX"
    ),
    RACE_WHITE = list(
      description = "White race indicator (1 = White, 0 = Asian, Black or African American, or Other)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (Asian, Black or African American, Other; source RACE = 2, 3, 6)",
      notes = "Source column RACE with RACE = 1 the effect level and RACE = 3, 6, 2 the reference (Data S1 equations). The Results text states that Asian is the reference group, and Figure 2 shows White as the only race with a ratio different from 1 (AUC24,ss 0.827 = 1 / 1.212), with Asian, Black and Other at 1.00; Table 1's footnote 'RACE (1), Asian' contradicts both and is treated as an error.",
      source_name = "RACE"
    ),
    CONMED_OSIMERTINIB = list(
      description = "Osimertinib co-administration indicator (1 = savolitinib given with osimertinib 80 mg once daily, 0 = savolitinib monotherapy or with placebo)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (savolitinib monotherapy)",
      notes = "Source column COMB (Data S1 equation for CL/F cat,cov; Figure 2 labels COMB = 1 'W_Osimertinib'). Time-fixed per subject in the source dataset (Table S3: 590 of 998 patients on the combination).",
      source_name = "COMB"
    ),
    CONMED_PPI = list(
      description = "Concomitant proton pump inhibitor use (1 = yes, 0 = no)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no concomitant PPI)",
      notes = "Source column PPI, subject-level in the source dataset (Table S3: 379 of 998 patients with concomitant PPI). Effects on V/F and KA only, so a PPI leaves AUC24,ss unchanged (Figure 2a).",
      source_name = "PPI"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Tested in the forward-addition / backward-elimination covariate search (Table S2) and not retained."
    ),
    ALB = list(
      description = "Baseline serum albumin",
      units = "g/L",
      type = "continuous",
      notes = "Source column BALB; tested (Table S2) and not retained."
    ),
    AST = list(
      description = "Baseline aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Source column BAST; tested (Table S2) and not retained."
    ),
    ALT = list(
      description = "Baseline alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      notes = "Source column BALT; tested (Table S2) and not retained."
    ),
    CREAT = list(
      description = "Baseline serum creatinine",
      units = "umol/L",
      type = "continuous",
      notes = "Source column BSCR; tested (Table S2) and not retained. Renal impairment category (RIST), hepatic impairment category (HIST), liver metastases (METS), baseline ECOG performance status (BECOG) and smoking status (SMK) were also tested (Table S3) and not retained."
    )
  )

  compartmentData <- list(
    depot = list(analyte = "savolitinib", units = "mg", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "savolitinib", units = "mg", specimen = "administration site", verified = TRUE),
    transit2 = list(analyte = "savolitinib", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "savolitinib", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 998L,
    n_studies = 10L,
    n_observations = 9706L,
    age_range = "23-92 years",
    age_median = "63 years",
    weight_range = "36.0-177 kg",
    weight_median = "62.2 kg",
    sex_female_pct = 48.5,
    race_ethnicity = c(White = 40.5, `Black or African American` = 1.4, Asian = 55.7, Other = 2.4),
    disease_state = "Locally advanced or metastatic solid tumours: EGFR-mutated MET-overexpressed and/or amplified non-small cell lung cancer after osimertinib (SAVANNAH, TATTON, and others), papillary renal cell carcinoma, pulmonary sarcomatoid carcinoma and NSCLC with MET exon 14 alterations, and other advanced solid tumours",
    dose_range = "Savolitinib 100-1000 mg once daily or 100-600 mg twice daily orally (Table S4), alone (40.9%) or with osimertinib 80 mg once daily (59.1%)",
    regions = "Multinational; 55.7% Asian",
    co_medication = "Concomitant proton pump inhibitor in 38.0%",
    notes = "Ten phase I-III studies (Table S1). Demographics from Tables S2 (continuous) and S3 (categorical), pooled column. 190 samples (1.92%) below the 1 ng/mL LLOQ were excluded. Median (IQR) baseline total bilirubin 8.45 (6.00-10.7) umol/L; baseline alkaline phosphatase 89.0 (67.3-128) U/L."
  )

  ini({
    lka <- log(2.64)
    label("First-order rate constant shared by depot -> transit1 -> transit2 -> central (1/h)") # Table 1, theta1 KA 2.64 1/h
    lvc <- log(173)
    label("Apparent volume of distribution V/F (L)") # Table 1, theta2 Vd/F 173 L
    lcl <- log(31.7)
    label("Apparent clearance CL/F (L/h)") # Table 1, theta3 CL/F 31.7 L/h

    e_wt_vc <- 0.364
    label("Power exponent of body weight / 62.2 kg on V/F (unitless)") # Table 1, theta6 'Vd/F ~ weight' 0.364; Data S1 Vd/F cont,cov equation
    e_alp_vc <- -0.0866
    label("Power exponent of baseline alkaline phosphatase / 89 U/L on V/F (unitless)") # Table 1, theta23 'Vd/F ~ BALP' -0.0866; Data S1
    e_alp_cl <- -0.0932
    label("Power exponent of baseline alkaline phosphatase / 89 U/L on CL/F (unitless)") # Table 1, theta24 'CL/F ~ BALP' -0.0932; Data S1
    e_tbili_cl <- -0.118
    label("Power exponent of baseline total bilirubin / 8.5 umol/L on CL/F (unitless)") # Table 1, theta8 'CL/F ~ BBILI' -0.118; Data S1
    e_conmed_osimertinib_cl <- 0.119
    label("Fractional change in CL/F with osimertinib co-administration (unitless)") # Table 1, theta10 'CL/F ~ COMB' 0.119; Data S1 (1 + 0.119 for COMB = 1)
    e_sexf_vc <- -0.160
    label("Fractional change in V/F for females (unitless)") # Table 1, theta11 'Vd/F ~ SEX' -0.160; Data S1 (1 - 0.160 for SEX = 1)
    e_sexf_cl <- -0.197
    label("Fractional change in CL/F for females (unitless)") # Table 1, theta12 'CL/F ~ SEX' -0.197; Data S1 (1 - 0.197 for SEX = 1)
    e_race_white_vc <- 0.198
    label("Fractional change in V/F for White patients (unitless)") # Table 1, theta13 'Vd/F ~ RACE (1)' 0.198; Data S1 (1 + 0.198 for RACE = 1)
    e_race_white_cl <- 0.212
    label("Fractional change in CL/F for White patients (unitless)") # Table 1, theta16 'CL/F ~ RACE (1)' 0.212; Data S1 (1 + 0.212 for RACE = 1)
    e_race_white_ka <- 0.326
    label("Fractional change in KA for White patients (unitless)") # Table 1, theta19 'KA ~ RACE (1)' 0.326; Data S1 (1 + 0.326 for RACE = 1)
    e_conmed_ppi_vc <- 0.0981
    label("Fractional change in V/F with a concomitant proton pump inhibitor (unitless)") # Table 1, theta9 'Vd/F ~ PPI' 0.0981; Data S1 (1 + 0.0981 for PPI = 1)
    e_conmed_ppi_ka <- -0.164
    label("Fractional change in KA with a concomitant proton pump inhibitor (unitless)") # Table 1, theta22 'KA ~ PPI' -0.164; Data S1 (1 - 0.164 for PPI = 1)

    etalka ~ 0.435 # Table 1, Omega(1,1) IIV-KA 0.435 (CV 73.8%)
    etalvc ~ 0.0560 # Table 1, Omega(2,2) IIV-Vd/F 0.0560 (CV 24.0%)
    etalcl ~ 0.135 # Table 1, Omega(3,3) IIV-CL/F 0.135 (CV 38.0%)

    addSd <- 1.07
    label("Additive residual error SD (ng/mL)") # Table 1, theta4 err_add 1.07
    propSd <- 0.495
    label("Proportional residual error SD (fraction)") # Table 1, theta5 err_prop 0.495
  })

  model({
    # Covariate multipliers (Data S1): continuous power terms and
    # categorical fractional-change terms, multiplied together.
    cov_ka <- (1 + e_conmed_ppi_ka * CONMED_PPI) *
      (1 + e_race_white_ka * RACE_WHITE)
    cov_vc <- (WT / 62.2)^e_wt_vc *
      (ALP / 89)^e_alp_vc *
      (1 + e_conmed_ppi_vc * CONMED_PPI) *
      (1 + e_sexf_vc * SEXF) *
      (1 + e_race_white_vc * RACE_WHITE)
    cov_cl <- (ALP / 89)^e_alp_cl *
      (TBILI / 8.5)^e_tbili_cl *
      (1 + e_sexf_cl * SEXF) *
      (1 + e_race_white_cl * RACE_WHITE) *
      (1 + e_conmed_osimertinib_cl * CONMED_OSIMERTINIB)

    ka <- exp(lka + etalka) * cov_ka
    vc <- exp(lvc + etalvc) * cov_vc
    cl <- exp(lcl + etalcl) * cov_cl

    kel <- cl / vc

    # Figure S3: Tr0 (dosing) -> Tr1 -> Tr2 -> central, each step at KA.
    d / dt(depot) <- -ka * depot
    d / dt(transit1) <- ka * depot - ka * transit1
    d / dt(transit2) <- ka * transit1 - ka * transit2
    d / dt(central) <- ka * transit2 - kel * central

    # Dose in mg and V/F in L give mg/L; x 1000 for ng/mL.
    Cc <- 1000 * central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
