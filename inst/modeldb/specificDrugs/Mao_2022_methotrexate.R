Mao_2022_methotrexate <- function() {
  description <- "Two-compartment population PK model with linear elimination for high-dose intravenous methotrexate in adults with primary central nervous system lymphoma (Mao 2022). Apparent clearance carries three covariates: a power function of Cockcroft-Gault creatinine clearance normalized to 98 mL/min (exponent 0.49), a power function of serum albumin normalized to 40 g/L (exponent 0.35), and a multiplicative factor of 0.89 for patients older than 60 years. Between-subject variability on CL, Vc, Q and Vp, a single shared between-occasion (per-course) variability on CL, and a proportional residual error."
  reference <- paste(
    "Mao J, Li Q, Li P, Qin W, Chen B, Zhong M. Evaluation and Application",
    "of Population Pharmacokinetic Models for Identifying Delayed",
    "Methotrexate Elimination in Patients With Primary Central Nervous",
    "System Lymphoma. Front Pharmacol. 2022;13:817673.",
    "doi:10.3389/fphar.2022.817673.",
    sep = " "
  )
  vignette <- "Mao_2022_methotrexate"
  units <- list(time = "h", dosing = "mg", concentration = "umol/L")

  # Doses are entered in mg and converted to umol/L inside model() with the
  # methotrexate molecular weight (see the Cc line). Mao 2022 fitted the
  # model to molar doses, so the state is carried in mg and only the
  # observation is molar.
  compartmentData <- list(
    central = list(analyte = "methotrexate", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "methotrexate", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = "Creatinine clearance by the Cockcroft-Gault equation (raw mL/min, NOT BSA-normalized)",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Table 1 footnote c: 'CrCL = [(140-Age (year)) x WT (kg)] / (0.818 x",
        "Scr (umol/L)) x (0.85 for female)', i.e. the Cockcroft-Gault",
        "equation on total body weight with serum creatinine in umol/L; the",
        "result is absolute mL/min and is not divided by BSA. Power form on",
        "CL/F normalized to 98 mL/min, the cohort median (Table 1 median 98,",
        "range 15.1-326.5; mean 104.2 +/- 34.2). Serum creatinine was",
        "recorded before each methotrexate infusion and alongside the",
        "plasma samples (Methods 2.1), so CRCL may be supplied per course or",
        "time-varying."
      ),
      source_name = "CrCL"
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power form on CL/F normalized to 40 g/L (Results 3.2.2 final-model",
        "equation). The cohort median is 39.0 g/L (Table 1, range",
        "24.0-50.0), so the reference is a rounded value rather than the",
        "median. Albumin was preferred to hematocrit as the protein-binding",
        "covariate (delta OFV -59.9 vs -14.8; Supplementary Table S5)."
      ),
      source_name = "ALB"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Enters only as a dichotomy: CL/F is multiplied by 0.89 when AGE > 60",
        "years (Results 3.2.2: 'The CL/F of elderly patients (age > 60) was",
        "11.0% lower than that of the younger patients'; final-model",
        "equation '... x 0.89, if age > 60'). The indicator is derived",
        "inside model() from the continuous AGE column. Supplementary Table",
        "S7 labels the simulated groups '< 60' and '>= 60'; the strict",
        "inequality of the printed equation is encoded. Cohort median 56",
        "years (range 28-76); 24.6% of patients were older than 60."
      ),
      source_name = "AGE"
    ),
    OCC = list(
      description = "Integer-valued methotrexate course (occasion) indicator for between-occasion variability on CL",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = paste(
        "Methods 2.2.2: 'IOV was assumed to be the same for all occasions', so",
        "a single between-occasion variance (24.7%, Table 3) applies to every",
        "course. Seventeen slots are provided because patients received",
        "1-17 courses (Table 1, Occasions median 4, range 1-17). Decomposed",
        "inside model() into binary indicators oc1..oc17; OCC = 0 or any",
        "value outside 1..17 zeros every indicator and gives CL with",
        "between-subject variability only. Pass OCC = 1 for a single-course",
        "simulation."
      ),
      source_name = "OCC"
    )
  )

  # Covariates that Mao 2022 screened but did not retain (Supplementary
  # Table S5). Documentation only -- none is referenced in model().
  covariatesDataExcluded <- list(
    HCT = list(
      description = "Hematocrit",
      units = "%",
      type = "continuous",
      notes = paste(
        "Tested on CL/F as the alternative protein-binding covariate to",
        "albumin (delta OFV -14.8 vs -59.9 for albumin; Supplementary Table",
        "S5 model 4) and not retained. Table 1 median 36.1% (range",
        "15.7-48.4)."
      )
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = paste(
        "Tested allometrically on the PK parameters (Supplementary Table S5",
        "model 6, delta OFV -3.8) and not retained. Table 1 median 69.0 kg",
        "(range 41.0-94.0)."
      )
    ),
    BSA = list(
      description = "Body surface area",
      units = "m^2",
      type = "continuous",
      notes = paste(
        "Tested allometrically on the PK parameters (Supplementary Table S5",
        "model 7, delta OFV -3.73) and not retained. Table 1 median 1.61",
        "m^2 (range 0.85-2.32)."
      )
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = paste(
        "Tested allometrically on the PK parameters (Supplementary Table S5",
        "model 8, OFV rose by 8.58) and not retained. Only two patients had",
        "BMI >= 30 (Discussion)."
      )
    )
  )

  population <- list(
    species = "human",
    n_subjects = 77L,
    n_studies = 1L,
    age_range = "28 to 76 years",
    age_median = "56 years (mean 54.6 +/- 9.2); 24.6% older than 60 years",
    weight_range = "41.0 to 94.0 kg",
    weight_median = "69.0 kg (mean 67.8 +/- 10.7)",
    height_range = "150 to 185 cm (median 170)",
    bsa_range = "0.85 to 2.32 m^2 (median 1.61)",
    sex_female_pct = 100 * 28 / 77,
    race_ethnicity = "Chinese (single centre in Shanghai); not tabulated",
    disease_state = "Primary central nervous system lymphoma treated with high-dose methotrexate (> 1 g/m^2).",
    renal_function = paste(
      "Serum creatinine median 66.0 umol/L (range 22.0-480.0); Cockcroft-Gault",
      "creatinine clearance median 98 mL/min (range 15.1-326.5). Hydration",
      "and urine alkalinization (urine pH > 7) began 12 h before each",
      "infusion."
    ),
    albumin = "Serum albumin median 39.0 g/L (range 24.0-50.0)",
    dose_range = paste(
      "Methotrexate 2.0-15.8 g per course (median 4.0 g; median 2.8 g/m^2,",
      "range 1.1-10.2 g/m^2) as an intravenous infusion of median 3 h (range",
      "1-28.25 h); 92.3% of courses were 1-4 h infusions. 377 courses in",
      "total (1-17 per patient)."
    ),
    regions = "China (Huashan Hospital, Fudan University, Shanghai).",
    co_medication = paste(
      "Leucovorin rescue every 6 h from 24 h after the start of infusion until",
      "methotrexate fell to 0.2 umol/L or below. Proton-pump inhibitors",
      "(lansoprazole in 969 samples, omeprazole 214, pantoprazole 58,",
      "esomeprazole 19) and dexamethasone (1,122 samples) were common; more",
      "than 85.4% of samples had a drug interaction score of 2 or more."
    ),
    notes = paste(
      "Retrospective therapeutic-drug-monitoring data collected June",
      "2011 to November 2016: 1,458 plasma methotrexate concentrations",
      "(median 16 per patient, range 3-67) sampled at 24, 48 and 72 h and",
      "until 0.2 umol/L or below, assayed by EMIT (Siemens Viva-Emit 2000;",
      "limit 0.3 umol/L). 567 concentrations were below the limit and were",
      "handled with Beal's M6 method (first BLQ value of each run set to",
      "LLOQ/2, later ones dropped). NONMEM 7.4, FOCE-I; ADVAN3 TRANS4.",
      "Baseline demographics from Table 1; final estimates from Table 3;",
      "covariate model from the Results 3.2.2 equation. A drug interaction",
      "score >= 2 was tested on CL/F (delta OFV -3.6) and not retained."
    )
  )

  ini({
    # Structural parameters -- Mao 2022 Table 3 'Final model' Estimate
    # column, applied in the final-model equation printed in Results 3.2.2:
    #   CL/F = 4.91 x (CrCL/98)^0.49 x (ALB/40)^0.35 x 0.89, if age > 60
    # Methotrexate was given intravenously; Table 3 footnote defines F as 'the
    # bioavailability relative to 1', so CL/F, V/F etc. are CL, V etc.
    lcl <- log(4.91); label("Clearance at CRCL 98 mL/min, ALB 40 g/L and age <= 60 years (L/h)") # Table 3 CL/F = 4.91 L/h (RSE 3.7%; bootstrap median 4.97, 95% CI 4.37-5.44)
    lvc <- log(18.4); label("Central volume of distribution (L)") # Table 3 Vc/F = 18.4 L (RSE 3.8%; bootstrap median 18.0, 95% CI 16.5-20.3)
    lq <- log(0.063); label("Intercompartmental clearance (L/h)") # Table 3 Q/F = 0.063 L/h (RSE 9.8%; bootstrap median 0.073, 95% CI 0.022-0.10)
    lvp <- log(2.18); label("Peripheral volume of distribution (L)") # Table 3 Vp/F = 2.18 L (RSE 14.1%; bootstrap median 2.17, 95% CI 1.59-2.77)

    # Covariate effects on CL (Table 3 'Covariate effect on CL/F').
    e_crcl_cl <- 0.49; label("Power exponent of creatinine clearance (/98 mL/min) on CL (unitless)") # Table 3 CrCL = 0.49 (RSE 22.6%; bootstrap median 0.50, 95% CI 0.29-0.69)
    e_alb_cl <- 0.35; label("Power exponent of serum albumin (/40 g/L) on CL (unitless)") # Table 3 ALB = 0.35 (RSE 50.6%; bootstrap median 0.35, 95% CI 0.031-0.72)
    e_age_gt60_cl <- 0.89; label("CL ratio for age > 60 vs <= 60 years (unitless)") # Table 3 AGE = 0.89 (RSE 9.4%; bootstrap median 0.90, 95% CI 0.72-1.06); Results 3.2.2 '0.89, if age > 60'

    # Between-subject variability. Table 3 reports %CV; converted to log-scale
    # variance as log(1 + CV^2). The paper reports no BSV correlations
    # (Supplementary Figure S1 shows eta-eta scatter only), so the block is
    # diagonal.
    etalcl ~ 0.042754 # Table 3 BSV CL/F 20.9% (RSE 25.5%; shrinkage 25.5%) -> log(1 + 0.209^2)
    etalvc ~ 0.037696 # Table 3 BSV Vc/F 19.6% (RSE 26.9%; shrinkage 20.9%) -> log(1 + 0.196^2)
    etalq ~ 0.152580 # Table 3 BSV Q/F 40.6% (RSE 18.1%; shrinkage 19.4%) -> log(1 + 0.406^2)
    etalvp ~ 0.088392 # Table 3 BSV Vp/F 30.4% (RSE 17.5%; shrinkage 47.7%) -> log(1 + 0.304^2)

    # Between-occasion (per-course) variability on log-CL. Methods 2.2.2: 'IOV
    # was assumed to be the same for all occasions', so occasions 2..17 each
    # carry their own eta with the variance fixed equal to the occasion-1
    # estimate (the NONMEM BLOCK SAME construct; Olivo_2024_methotrexate.R
    # precedent).
    etaiov_cl_1 ~ 0.059220 # Table 3 IOV on CL 24.7% (RSE 20.4%; shrinkage 43.8%) -> log(1 + 0.247^2)
    etaiov_cl_2 ~ fixed(0.059220) # same shared IOV variance
    etaiov_cl_3 ~ fixed(0.059220)
    etaiov_cl_4 ~ fixed(0.059220)
    etaiov_cl_5 ~ fixed(0.059220)
    etaiov_cl_6 ~ fixed(0.059220)
    etaiov_cl_7 ~ fixed(0.059220)
    etaiov_cl_8 ~ fixed(0.059220)
    etaiov_cl_9 ~ fixed(0.059220)
    etaiov_cl_10 ~ fixed(0.059220)
    etaiov_cl_11 ~ fixed(0.059220)
    etaiov_cl_12 ~ fixed(0.059220)
    etaiov_cl_13 ~ fixed(0.059220)
    etaiov_cl_14 ~ fixed(0.059220)
    etaiov_cl_15 ~ fixed(0.059220)
    etaiov_cl_16 ~ fixed(0.059220)
    etaiov_cl_17 ~ fixed(0.059220)

    # Residual error. Results 3.2.2: 'The exponential model provided the best
    # result for the residual variability'; an exponential error on the
    # untransformed scale is proportional to first order.
    propSd <- 0.401; label("Proportional residual error (fraction)") # Table 3 Residual variability Proportional 40.1% (RSE 6.3%; bootstrap median 39.6, 95% CI 34.9-44.9)
  })

  model({
    # Decompose the integer course column into binary occasion indicators that
    # multiplex the shared-variance IOV etas onto log-CL.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    oc5 <- (OCC == 5)
    oc6 <- (OCC == 6)
    oc7 <- (OCC == 7)
    oc8 <- (OCC == 8)
    oc9 <- (OCC == 9)
    oc10 <- (OCC == 10)
    oc11 <- (OCC == 11)
    oc12 <- (OCC == 12)
    oc13 <- (OCC == 13)
    oc14 <- (OCC == 14)
    oc15 <- (OCC == 15)
    oc16 <- (OCC == 16)
    oc17 <- (OCC == 17)
    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2 + oc3 * etaiov_cl_3 +
      oc4 * etaiov_cl_4 + oc5 * etaiov_cl_5 + oc6 * etaiov_cl_6 +
      oc7 * etaiov_cl_7 + oc8 * etaiov_cl_8 + oc9 * etaiov_cl_9 +
      oc10 * etaiov_cl_10 + oc11 * etaiov_cl_11 + oc12 * etaiov_cl_12 +
      oc13 * etaiov_cl_13 + oc14 * etaiov_cl_14 + oc15 * etaiov_cl_15 +
      oc16 * etaiov_cl_16 + oc17 * etaiov_cl_17

    # Age dichotomy of the final-model equation ('0.89, if age > 60').
    age_gt60 <- (AGE > 60)

    # Individual PK parameters (Results 3.2.2 final-model equation).
    cl <- exp(lcl + etalcl + iov_cl) * (CRCL / 98)^e_crcl_cl * (ALB / 40)^e_alb_cl *
      e_age_gt60_cl^age_gt60
    vc <- exp(lvc + etalvc)
    q <- exp(lq + etalq)
    vp <- exp(lvp + etalvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Two-compartment disposition; intravenous infusion into central.
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Plasma methotrexate in umol/L: central / vc is mg/L, and dividing by the
    # molecular weight 454.44 g/mol (mg/mmol) then x 1000 gives umol/L. Results
    # 3.1 prints the molecular weight used to convert doses to molar units as
    # '222 g mol-1', which is not the molecular weight of methotrexate
    # (C20H22N8O5, 454.44 g/mol). The paper's own Monte Carlo simulation
    # (Supplementary Table S7) is reproduced with 454.44 and not with 222 (see
    # the vignette), so 454.44 is used here.
    Cc <- central / vc / 454.44 * 1000

    Cc ~ prop(propSd)
  })
}
