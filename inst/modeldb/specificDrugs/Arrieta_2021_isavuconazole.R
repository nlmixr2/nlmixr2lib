Arrieta_2021_isavuconazole <- function() {
  description <- "Three-compartment population PK model for isavuconazole (administered as the prodrug isavuconazonium sulfate) in immunocompromised children aged 1 to <18 years pooled with adults from a phase 1 intravenous study (Arrieta 2021). Zero-order (1-h intravenous infusion) and first-order (oral capsule) input into the central compartment, linear elimination, and allometric body-weight scaling (reference 70 kg) on all clearances and volumes; no other covariates retained."
  reference <- "Arrieta AC, Neely M, Day JC, Rheingold SR, Sue PK, Muller WJ, Danziger-Isakov LA, Chu J, Yildirim I, McComsey GA, Frangoul HA, Chen TK, Statler VA, Steinbach WJ, Yin DE, Hamed K, Jones ME, Lademacher C, Desai A, Micklus K, Phillips DL, Kovanda LL, Walsh TJ. Safety, Tolerability, and Population Pharmacokinetics of Intravenous and Oral Isavuconazonium Sulfate in Pediatric Patients. Antimicrobial Agents and Chemotherapy. 2021;65(8):e00290-21. doi:10.1128/AAC.00290-21"
  vignette <- "Arrieta_2021_isavuconazole"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    depot = list(analyte = "isavuconazole", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "isavuconazole", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "isavuconazole", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral2 = list(analyte = "isavuconazole", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric scaling P = P_adult * (WT/70)^x on CL, Q3 and Q4 (x = 0.75) and on V2, V3 and V4 (x = 1), reference 70 kg (Arrieta 2021 Methods, 'Population PK model', allometric equation). The paper names x only as 'the relevant allometric component' and does not print its value; no exponent appears among the estimated parameters in Table S2, so the conventional fixed values 0.75 / 1 are used. Pediatric cohort weight 10.9-103.5 kg (Table 1).",
      source_name = "WT"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened by stepwise covariate modeling (PsN forward inclusion / backward elimination) on clearance and volume; not statistically significant for any PK parameter (Arrieta 2021 Results, 'Population PK model analysis')."
    ),
    SEXF = list(
      description = "Female sex indicator (1 = female)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male)",
      notes = "Sex screened and not retained (Arrieta 2021 Results, 'Population PK model analysis'); the paper's own coding of sex is not stated."
    ),
    RACE_WHITE = list(
      description = "White race indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-White)",
      notes = "Race screened and not retained (Arrieta 2021 Results); the paper's race coding is not stated. Discussion notes too few non-White patients to analyse race formally."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened and not retained (Arrieta 2021 Results, 'Population PK model analysis')."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units = "(unit not stated)",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened and not retained (Arrieta 2021 Results, 'Population PK model analysis')."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened and not retained (Arrieta 2021 Results, 'Population PK model analysis')."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened and not retained (Arrieta 2021 Results, 'Population PK model analysis')."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units = "(unit not stated)",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened and not retained (Arrieta 2021 Results, 'Population PK model analysis')."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "(unit not stated)",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened and not retained (Arrieta 2021 Results, 'Population PK model analysis')."
    ),
    ALP = list(
      description = "Alkaline phosphatase",
      units = "U/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Screened and not retained (Arrieta 2021 Results, 'Population PK model analysis')."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 45,
    n_studies = 2,
    age_range = "1-17 years (pediatric cohorts); adults from a phase 1 intravenous study (demographics not reported in this paper)",
    age_median = "10.0 years (i.v. cohort), 12.0 years (oral cohort)",
    weight_range = "10.9-103.5 kg (pediatric)",
    weight_median = "33.7 kg (i.v. cohort), 42.6 kg (oral cohort)",
    sex_female_pct = 37.0,
    race_ethnicity = c(
      White = 71.7,
      Black = 10.9,
      Asian = 6.5,
      `American Indian or Alaska native` = 2.2,
      `Pacific Islander` = 2.2,
      Other = 6.5
    ),
    disease_state = "Immunocompromised children at risk for invasive mycoses (mainly acute myeloid / lymphoblastic leukemia, neuroblastoma and other solid tumors, aplastic anemia; 15% prior hematopoietic stem cell transplant) receiving antifungal prophylaxis.",
    dose_range = "Isavuconazonium sulfate 10 mg/kg (maximum 372 mg; = 5.38 mg/kg or maximum 200 mg isavuconazole) q8h for 6 doses then once daily for up to 26 days; i.v. as a 1-h infusion (1 to <18 years; 372 mg above 37 kg) or orally as 74.5-mg capsules (6 to <18 years; 372 mg above 32 kg).",
    regions = "United States (11 i.v. and 12 oral centers)",
    n_observations = "551 pediatric plasma samples (333 i.v., 218 oral) plus adult i.v. data",
    notes = "Pediatric PK analysis set: 26 i.v. (9 aged 1 to <6, 8 aged 6 to <12, 9 aged 12 to <18 years) and 19 oral (9 aged 6 to <12, 10 aged 12 to <18 years) patients (Arrieta 2021 Table 1, Table 2; NCT03241550). Pooled with i.v. data from the adult phase 1 renal-impairment study of Townsend 2017 (Eur J Clin Pharmacol 73:669-678); adult numbers are not reported in this paper. Sex and race percentages are over the 46-patient safety analysis set (Table 1). Estimation in NONMEM with PsN 4.7.0; LLOQ 100 ng/mL."
  )

  ini({
    # Structural parameters for a 70-kg adult (Arrieta 2021 Table S2, 'Parameter
    # estimates for the best population pharmacokinetic model'). V2 is the
    # central volume, Q3/V3 the first and Q4/V4 the second (deep) peripheral
    # compartment.
    lcl <- log(2.55); label("Clearance CL for a 70-kg subject (L/h)") # Table S2: CL = 2.55 L/h
    lvc <- log(17.8); label("Central volume V2 for a 70-kg subject (L)") # Table S2: V2 = 17.80 L
    lq <- log(30.3); label("Intercompartmental clearance Q3 to the first peripheral compartment for a 70-kg subject (L/h)") # Table S2: Q3 = 30.30 L/h
    lvp <- log(26.0); label("First peripheral volume V3 for a 70-kg subject (L)") # Table S2: V3 = 26.00 L
    lq2 <- log(24.5); label("Intercompartmental clearance Q4 to the second peripheral compartment for a 70-kg subject (L/h)") # Table S2: Q4 = 24.50 L/h
    lvp2 <- log(254); label("Second peripheral volume V4 for a 70-kg subject (L)") # Table S2: V4 = 254.0 L
    lka <- log(0.162); label("First-order oral absorption rate constant Ka (1/h)") # Table S2: Ka = 0.162 (unit printed as 'h'; a rate constant, so 1/h)
    lfdepot <- log(0.95); label("Oral bioavailability F1 of the pediatric capsule (fraction)") # Table S2: F1 = 0.95

    # Allometric exponents. Arrieta 2021 Methods writes
    # P_pediatric = P_adults * (WT/70)^x without printing x; Table S2 lists no
    # exponent among the estimates, so the conventional fixed values are used.
    e_wt_cl_q <- fixed(0.75); label("Allometric exponent of (WT/70) on CL, Q3 and Q4 (unitless)") # Methods 'Population PK model' allometric equation; value not printed, conventional 0.75 assumed
    e_wt_vc_vp <- fixed(1); label("Allometric exponent of (WT/70) on V2, V3 and V4 (unitless)") # Methods 'Population PK model' allometric equation; value not printed, conventional 1 assumed

    # Inter-individual variability. Table S2 'Variability (%)' reports each
    # IIV as a CV%. Its SE column is on the variance scale and matches
    # omega^2 = (CV/100)^2 (CL: 0.034 / 0.4516^2 = 17% RSE, printed 17; V4:
    # 0.128 / 0.6811^2 = 28%, printed 28; Q4: 0.052 / 0.4571^2 = 25%, printed
    # 25), so the variances are the squared CVs.
    etalcl ~ 0.20394 # Table S2: IIV CL = 45.16% -> 0.4516^2
    etalvp ~ 0.18697 # Table S2: IIV V3 = 43.24% -> 0.4324^2
    etalq2 ~ 0.20894 # Table S2: IIV Q4 = 45.71% -> 0.4571^2
    etalvp2 ~ 0.46390 # Table S2: IIV V4 = 68.11% -> 0.6811^2

    # Residual error. Table S2 reports 'Residual error sigma^2' = 40.86 under
    # the Variability (%) heading, with SE 0.00934 and 6% RSE: 0.00934 /
    # 0.4086^2 = 5.6%, so 40.86 is the SD in percent. Figure S2 goodness-of-fit
    # plots are on log concentration (log-transform-both-sides), i.e. an
    # additive error on the log scale, encoded as proportional in linear space.
    propSd <- 0.4086; label("Proportional residual error (fraction)") # Table S2: Residual error = 40.86%
  })

  model({
    # Allometric size scaling, reference 70 kg (Methods, allometric equation).
    wt_cl_q <- (WT / 70)^e_wt_cl_q
    wt_vc_vp <- (WT / 70)^e_wt_vc_vp

    cl <- exp(lcl + etalcl) * wt_cl_q
    vc <- exp(lvc) * wt_vc_vp
    q <- exp(lq) * wt_cl_q
    vp <- exp(lvp + etalvp) * wt_vc_vp
    q2 <- exp(lq2 + etalq2) * wt_cl_q
    vp2 <- exp(lvp2 + etalvp2) * wt_vc_vp
    ka <- exp(lka)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # Oral doses enter the depot (first-order input); i.v. doses are 1-h
    # zero-order infusions into central (rate or dur on the dosing record).
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1 - k13 * central + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2

    f(depot) <- exp(lfdepot)

    # Dose in mg isavuconazole (372 mg isavuconazonium sulfate = 200 mg
    # isavuconazole) / volume in L = mg/L (= ug/mL).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
