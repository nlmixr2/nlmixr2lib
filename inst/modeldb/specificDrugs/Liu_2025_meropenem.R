Liu_2025_meropenem <- function() {
  description <- paste(
    "Two-compartment population PK model for intravenous meropenem in 101",
    "critically ill children (202 therapeutic-drug-monitoring concentrations)",
    "from two Chinese paediatric intensive care units (Liu 2025). The cohort is",
    "dominated by the very young -- 28.7% neonates and 47.5% under three months",
    "-- and by augmented renal clearance: the median bedside-Schwartz eGFR is",
    "123.4 mL/min/1.73 m2, well above the age-specific healthy reference in",
    "every stratum. Table 2 labels every structural row per kilogram of body",
    "weight (CL 0.24 L/h/kg, V1 1.53 L/kg, Q 0.014 L/h/kg, V2 6.06 L/kg); those",
    "are read here as typical ABSOLUTE values quoted per median kilogram, so",
    "each is multiplied by the cohort median weight of 7.5 kg rather than by the",
    "individual weight. Multiplying by the individual weight instead would count",
    "body weight twice and is contradicted by the paper's own Figure 4 -- see the",
    "vignette's Assumptions and deviations section. Clearance then carries a",
    "power term on body weight (exponent 0.43) and on eGFR (exponent 0.96, i.e.",
    "very close to proportional, as expected for a drug cleared almost entirely",
    "by the kidney), the central volume a power term on body weight (0.37) and",
    "the inter-compartmental clearance a power term on body weight (1.54); all",
    "are normalised to the cohort medians of 7.5 kg and 123.4 mL/min/1.73 m2.",
    "Q is small relative to CL and V2 is large, so the peripheral compartment",
    "acts as a slowly-equilibrating deep sink that carries only a few percent of",
    "the elimination flux over an 8 h dosing interval. Log-normal IIV is placed",
    "on clearance and on the central volume, and residual variability is",
    "proportional. The paper's companion whole-body PBPK model was built in",
    "PK-Sim and is NOT encoded here: its physiology (Rodgers and Rowland tissue",
    "partitioning, 'PK-Sim Standard' cellular permeabilities, organ volumes and",
    "blood flows, and the OAT3 / tubular-secretion ontogeny functions) lives in",
    "the platform database rather than in the paper or its supplement, so the",
    "PBPK ODEs cannot be reproduced from on-disk sources."
  )
  reference <- paste(
    "Liu Y, He H, Zhang SS, Zhou J, Zhu JW, Xu J, Miao HJ, Chen JH, Hao K.",
    "PopPK and PBPK Models Guide Meropenem Dosing in Critically Ill Children",
    "with Augmented Renal Clearance.",
    "Pharmaceutics. 2025;17(12):1544.",
    "doi:10.3390/pharmaceutics17121544"
  )
  vignette <- "Liu_2025_meropenem"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central     = list(analyte = "meropenem", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "meropenem", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Body weight enters ONCE, through the power terms of Equations 4, 5 and 6,",
        "normalised to BWmed and encoded here as (WT/7.5)^beta. BWmed is the Table 1",
        "cohort median body weight, 7.5 kg (IQR 3.4-12.1). Table 2 gives",
        "beta_CL,WT = 0.43 (RSE 10.3%), beta_V1,WT = 0.37 (RSE 19.7%) and",
        "beta_Q,WT = 1.54 (RSE 12.1%).",
        "Table 2 additionally LABELS each structural row per kilogram ('CL (L/h/kg)',",
        "'V1 (L/kg)', 'Q (L/h/kg)', 'V2 (L/kg)'). Those labels are read here as typical",
        "absolute values quoted per median kilogram, i.e. the absolute clearance is",
        "0.24 * 7.5 = 1.8 L/h at the cohort medians, NOT 0.24 * WT. Multiplying by the",
        "individual weight as well as applying the (WT/7.5)^0.43 term would count body",
        "weight twice and give a net exponent of 1.43 on absolute clearance.",
        "The paper's own Figure 4 rules that out. Its four PopPK panels all give",
        "20 mg/kg and report an essentially flat steady-state AUC across strata",
        "(roughly 160, 140, 158 and 158 mg.h/L for the non-ARC subgroups) even though",
        "median eGFR rises from 88 to 175.8 and median weight from about 3.4 to 12 kg.",
        "A flat mg/kg-normalised AUC under a rising eGFR requires a net weight exponent",
        "on absolute clearance near 0.26-0.46 depending on the assumed stratum weights,",
        "which brackets beta_CL,WT = 0.43 and is more than a full unit away from the",
        "per-kg reading's 1.43. That reading predicts AUC falling roughly threefold from",
        "neonates to children over three months, and Cmax falling with age where Figure 4",
        "shows it rising. Regressing log AUC(tau) on log WT and log eGFR over the",
        "vignette's simulated cohort confirms it a third time: the weight slope comes",
        "back at 0.61 against the 0.57 this encoding predicts and the -0.43 the per-kg",
        "reading predicts. What does NOT reproduce is the exposure LEVEL, which sits",
        "about twofold below Figure 4 uniformly across panels; see the vignette's",
        "Assumptions and deviations section.",
        "No V2 covariate was retained; the paper prints no equation for V2."
      ),
      source_name        = "BW"
    ),
    CRCL = list(
      description        = "Estimated glomerular filtration rate by the bedside Schwartz equation",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Methods Equation 1: eGFR = k * H / Scr, with k = 0.33 for preterm, 0.45 for",
        "full-term and 0.41 for other paediatric patients, H the height in cm and Scr the",
        "serum creatinine in mg/dL; the result is expressed in mL/min/1.73 m2. This is a",
        "creatinine-based, body-surface-area-normalised estimate and so is the CRCL",
        "canonical, not an absolute mL/min clearance -- do NOT de-normalise it before use.",
        "The reference value used here is the Table 1 cohort median of",
        "123.4 mL/min/1.73 m2 (IQR 80.9-180.1), which is the natural referent for the",
        "'eGFRmed' the paper's covariate model normalises to.",
        "Renal function in this cohort is markedly supranormal: Table 1 gives stratum",
        "medians of 56.7 for neonates, 109.7 for 28 d to 3 months, 175.8 for 3 months to",
        "2 years, 164.7 for 2-6 years and 177.9 for 6-15 years, and the Discussion notes",
        "that term neonates ran at 88 against a healthy reference of 59 at four weeks.",
        "Table 2 gives beta_CL,eGFR = 0.96 (RSE 9.19%, bootstrap 95% CI 0.72-1.29), an",
        "exponent statistically indistinguishable from 1, i.e. clearance is essentially",
        "proportional to eGFR. The paper stratifies its dosing recommendations (Table 5) by",
        "eGFR expressed as a percentage of the age-expected level: 50% renal impairment,",
        "100% no ARC, 150% moderate ARC and 200% severe ARC."
      ),
      source_name        = "eGFR"
    )
  )

  covariatesDataExcluded <- list(
    HT = list(
      description = "Height",
      units       = "cm",
      type        = "continuous",
      notes       = paste(
        "Methods 2.2 lists height among the candidate continuous covariates screened by",
        "stepwise regression. It was not retained in the final model. It does enter the",
        "model indirectly, as an input to the Schwartz eGFR of Equation 1."
      )
    ),
    WBC = list(
      description = "White blood cell count",
      units       = "10^9/L",
      type        = "continuous",
      notes       = "Screened as a candidate covariate (Methods 2.2); not retained in the final model."
    ),
    NEUTPCT = list(
      description = "Neutrophil percentage",
      units       = "%",
      type        = "continuous",
      notes       = "Screened as a candidate covariate (Methods 2.2); not retained in the final model."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened as a candidate covariate (Methods 2.2); not retained in the final model."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened as a candidate covariate (Methods 2.2); not retained in the final model."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units       = "umol/L",
      type        = "continuous",
      notes       = paste(
        "Screened as a candidate covariate (Methods 2.2); not retained in the final model",
        "in its own right. It enters indirectly as the denominator of the Schwartz eGFR of",
        "Equation 1. Table 1 median 22.1 umol/L (IQR 15.5-33.9)."
      )
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/L",
      type        = "continuous",
      notes       = paste(
        "Screened as a candidate covariate (Methods 2.2); not retained in the final model.",
        "Table 1 prints 33.8 +/- 5.35 under a 'g/dL' header, but the magnitude is a g/L",
        "value; the unit label in that table is wrong (the same applies to total protein)."
      )
    ),
    TP = list(
      description = "Total protein",
      units       = "g/L",
      type        = "continuous",
      notes       = paste(
        "Screened as a candidate covariate (Methods 2.2); not retained in the final model.",
        "Table 1 prints 56.5 +/- 10.28 under a 'g/dL' header; as for albumin the magnitude",
        "is a g/L value and the printed unit is wrong."
      )
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "unitless",
      type        = "categorical",
      notes       = paste(
        "Gender was one of the two categorical covariates screened (Methods 2.2); not",
        "retained in the final model. Table 1: 41 female (40.6%), 60 male (59.4%).",
        "Figure S2E also found no relationship between gender and PBPK prediction error."
      )
    ),
    PRETERM = list(
      description = "Preterm birth indicator",
      units       = "unitless",
      type        = "categorical",
      notes       = paste(
        "Preterm status was the second categorical covariate screened (Methods 2.2); not",
        "retained in the final model. Of the 29 neonates, 13 were preterm and 16 term.",
        "Preterm status does enter indirectly through the Schwartz coefficient k of",
        "Equation 1 (0.33 preterm versus 0.45 full-term)."
      )
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 101L,
    n_studies      = 2L,
    n_observations = 202L,
    age_range      = "preterm neonates to <14 years; Table 1 median 20.3 weeks (IQR 3.2-114.6)",
    age_median     = "20.3 weeks",
    weight_range   = "Table 1 median 7.5 kg (IQR 3.4-12.1)",
    weight_median  = "7.5 kg",
    sex_female_pct = 40.6,
    disease_state  = paste(
      "Critically ill children in the paediatric intensive care unit.",
      "Infection-related diagnoses: sepsis or septic shock 78 (77.2%), pneumonia",
      "58 (57.4%), bacterial meningitis 32 (31.7%), intra-abdominal infection 16,",
      "skin and skin-structure infection 3."
    ),
    renal_function = paste(
      "Augmented renal clearance is the defining feature of the cohort. Bedside-Schwartz",
      "eGFR median 123.4 mL/min/1.73 m2 (IQR 80.9-180.1); by stratum, neonates 56.7",
      "(preterm 33.1, term 88.0), 28 d to 3 months 109.7, 3 months to 2 years 175.8,",
      "2-6 years 164.7, 6-15 years 177.9. Serum creatinine median 22.1 umol/L",
      "(IQR 15.5-33.9)."
    ),
    dose_range     = paste(
      "Intravenous meropenem, daily dose 22-157 mg/kg/day (mean 125.2 +/- 28.7),",
      "dosing intervals from every 24 h to every 6 h. Simulated regimens were",
      "20 mg/kg q12h or q8h given as 60 min infusions."
    ),
    regions        = "China (Shanghai and Nanjing)",
    notes          = paste(
      "Retrospective therapeutic-drug-monitoring data collected January 2020 to",
      "December 2023 in the paediatric intensive care units of Xinhua Hospital",
      "(Shanghai Jiao Tong University School of Medicine) and the Children's Hospital",
      "of Nanjing Medical University. Baseline characteristics are Table 1. Of the 202",
      "measured concentrations, 14 were below the limit of detection. Age strata:",
      "neonates <28 d 29 (28.7%; 13 preterm, 16 term), 28 d to <3 months 19 (18.8%),",
      "3 months to <2 years 21 (20.8%), 2 to <6 years 20 (19.8%), 6 to <14 years",
      "12 (11.9%). Fitted with SAEM in Monolix 2023."
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # Structural parameters: Liu 2025 Table 2, 'Stochastic Approximation /
    # Estimate' column. Every structural row of Table 2 is LABELLED per
    # kilogram; each is read here as a typical ABSOLUTE value quoted per
    # median kilogram, so each theta below is multiplied by BWmed = 7.5 kg
    # in model() to give the absolute L/h or L used in the ODEs. See the
    # vignette's Assumptions and deviations section for the Figure 4
    # evidence that rules out multiplying by the individual weight instead.
    # ---------------------------------------------------------------------

    # Table 2 row 'CL (L/h/kg)' = 0.24 (RSE 11.1%; bootstrap median 0.23,
    # 95% CI 0.16-0.34). Absolute clearance at the cohort medians of 7.5 kg
    # and 123.4 mL/min/1.73 m2 is 0.24 * 7.5 = 1.8 L/h, because both
    # covariate terms equal exactly 1 there.
    lcl <- log(0.24); label("Clearance per median kilogram at the median covariates (L/h/kg)")
    # Table 2 row 'V1 (L/kg)' = 1.53 (RSE 15.8%; bootstrap median 1.64,
    # 95% CI 0.90-2.79); absolute V1 at the median weight is 11.5 L.
    lvc <- log(1.53); label("Central volume of distribution per median kilogram (L/kg)")
    # Table 2 row 'Q (L/h/kg)' = 0.014 (RSE 42.4%; bootstrap median 0.014,
    # 95% CI 0.0002-0.038). Poorly identified, as the RSE and the bootstrap
    # interval that nearly touches zero both show.
    lq <- log(0.014); label("Inter-compartmental clearance per median kilogram (L/h/kg)")
    # Table 2 row 'V2 (L/kg)' = 6.06 (RSE 23.3%). The bootstrap columns for
    # this row are internally inconsistent -- the printed 2.5%ile of 9.4
    # exceeds the printed median of 6.08 -- so only the point estimate is
    # used. No covariate and no IIV were retained on V2, and the paper prints
    # no equation for it.
    lvp <- log(6.06); label("Peripheral volume of distribution per median kilogram (L/kg)")

    # ---------------------------------------------------------------------
    # Covariate effects. Equations 4-6, normalised to the Table 1 cohort
    # medians BWmed = 7.5 kg and eGFRmed = 123.4 mL/min/1.73 m2. See the
    # vignette's Assumptions and deviations section for why the printed
    # exp(beta * BW/BWmed) form is read as the power form (BW/BWmed)^beta.
    # ---------------------------------------------------------------------

    # Table 2 row 'beta CL,WT' = 0.43 (RSE 10.3%; bootstrap median 0.42,
    # 95% CI 0.33-0.65). This is the NET exponent on absolute clearance,
    # which is why it sits below the 0.75 of classical allometry: eGFR is
    # already in the model and absorbs much of the size and maturation
    # signal in a cohort spanning preterm neonates to 14-year-olds.
    e_wt_cl <- 0.43; label("Power exponent of body weight on clearance (unitless)")
    # Table 2 row 'beta CL,eGFR' = 0.96 (RSE 9.19%; bootstrap median 0.97,
    # 95% CI 0.72-1.29). Indistinguishable from 1: clearance is essentially
    # proportional to eGFR.
    e_crcl_cl <- 0.96; label("Power exponent of eGFR on clearance (unitless)")
    # Table 2 row 'beta V1,WT' = 0.37 (RSE 19.7%; bootstrap median 0.35,
    # 95% CI 0.09-0.56)
    e_wt_vc <- 0.37; label("Power exponent of body weight on central volume (unitless)")
    # Table 2 row 'beta Q,WT' = 1.54 (RSE 12.1%; bootstrap median 1.62,
    # 95% CI 1.04-4.27)
    e_wt_q <- 1.54; label("Power exponent of body weight on inter-compartmental clearance (unitless)")

    # ---------------------------------------------------------------------
    # Inter-individual variability. Results 3.2: 'incorporating log-normal
    # IIV in clearance (CL) and central volume of distribution (V1)'.
    # Monolix reports the omega rows of a parameter table as the STANDARD
    # DEVIATION of the random effect on the log scale, so the variance
    # encoded here is the square of the printed value.
    # ---------------------------------------------------------------------

    # Table 2 row 'omega CL' = 0.41 (RSE 13.1%; bootstrap median 0.41,
    # 95% CI 0.27-0.58); 0.41^2 = 0.1681, i.e. 42.8% CV.
    etalcl ~ 0.1681
    # Table 2 row 'omega V1' = 0.45 (RSE 21.1%; bootstrap median 0.36,
    # 95% CI 0.11-0.91); 0.45^2 = 0.2025, i.e. 47.4% CV.
    # Encoded DIAGONAL. Section 3.5 describes the simulations as using
    # 'log-normal IIV on CL and V1 with covariance', and Methods 2.2 says a
    # correlation was retained when it exceeded 0.30, but Table 2 prints no
    # correlation or covariance row and no such value appears anywhere in the
    # paper or its supplement. The off-diagonal is therefore omitted rather
    # than invented; see the vignette's Assumptions and deviations section.
    etalvc ~ 0.2025

    # ---------------------------------------------------------------------
    # Residual error. Results 3.2: 'Residual variability was effectively
    # captured by a proportional error model.'
    # ---------------------------------------------------------------------

    # Table 2 row 'Proportion error' = 0.34 (RSE 11.6%; bootstrap median 0.34,
    # 95% CI 0.24-0.41). Monolix reports the proportional coefficient b of a
    # proportional error model directly on the SD scale, so 0.34 is a 34%
    # proportional error and is used as-is.
    propSd <- 0.34; label("Proportional residual error (fraction)")
  })

  model({
    # Equation 4. The Table 2 value is a typical ABSOLUTE clearance quoted
    # per median kilogram, so it is multiplied by BWmed = 7.5 kg (a
    # constant), not by the individual weight. Both covariate terms are
    # normalised to the Table 1 cohort medians, so a 7.5 kg subject with an
    # eGFR of 123.4 mL/min/1.73 m2 has cl = 0.24 * 7.5 = 1.8 L/h exactly.
    cl <- exp(lcl + etalcl) * 7.5 * (WT / 7.5)^e_wt_cl * (CRCL / 123.4)^e_crcl_cl

    # Equation 5, same convention: 1.53 * 7.5 = 11.5 L at the median weight.
    # Note that the paper's own text under Equation 5 calls the coefficient
    # 'beta_CL,BW'; that is a typo for beta_V1,BW, which Table 2 lists
    # separately as 0.37.
    vc <- exp(lvc + etalvc) * 7.5 * (WT / 7.5)^e_wt_vc

    # Equation 6. Q carries no IIV.
    q <- exp(lq) * 7.5 * (WT / 7.5)^e_wt_q

    # V2 has no covariate and no IIV, so it is constant across subjects at
    # the Table 2 'V2 (L/kg)' row times the median weight: 6.06 * 7.5 = 45.5 L.
    vp <- exp(lvp) * 7.5

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Two-compartment model with linear elimination from the central
    # compartment (Results 3.2). Q is small relative to CL and V2 is large,
    # so peripheral1 behaves as a slowly-equilibrating deep sink: at the
    # cohort medians k12 = 0.0092/h and k21 = 0.0023/h against kel = 0.157/h.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    # Dose in mg and vc in L give mg/L, which is the unit of the measured
    # plasma concentrations. The paper's pharmacodynamic target is written on
    # UNBOUND concentration (100% fT > MIC); meropenem is only about 2% bound
    # (supplementary Table S3 gives fu = 0.98), and the paper applies no
    # unbound-fraction correction of its own, so Cc is compared against MIC
    # targets directly.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
