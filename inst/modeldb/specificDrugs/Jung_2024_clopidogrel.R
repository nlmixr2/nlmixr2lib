Jung_2024_clopidogrel <- function() {
  description <- paste(
    "Joint parent-plus-two-metabolite population PK model for oral",
    "clopidogrel with a sequentially-fitted platelet-reactivity (PRU)",
    "turnover PD submodel, in healthy Korean male adults stratified by",
    "CYP2C19 metabolizer phenotype (Jung 2024). A hepatic first-pass",
    "compartment (`liver`) sits between the depot and the central",
    "compartment; all absorbed parent drug is assumed to be metabolized in",
    "the liver, so the parent's only elimination pathway is the hepatic",
    "metabolic clearance CLc. Clopidogrel itself is described by a",
    "two-compartment model (`central` + `peripheral1`) exchanging with the",
    "liver at Qc and with the periphery at Qp; the hepatic volume was set",
    "equal to the central volume (VH = Vc) because the two were not",
    "separately identifiable. The hepatic metabolic flux CLc * (liver / VH)",
    "is split three ways by two nested logit-scale fractions: fm1 * fm2 to",
    "the active thiol metabolite clopidogrel H4 (one compartment,",
    "`central_h4`), (1 - fm1) * fm2 to the inactive clopidogrel carboxylic",
    "acid (two compartments, `central_cloca` + `peripheral1_cloca`),",
    "and 1 - fm2 to unmeasured other metabolites. Each metabolite formation",
    "flux carries a molar-mass ratio scaling factor (1.106 for H4, 0.956",
    "for the carboxylic acid) so that a mass-unit parent flux produces the",
    "correct mass-unit metabolite amount. CYP2C19 phenotype is the only",
    "retained covariate and acts as an additive shift on the logit scale of",
    "both fm1 and fm2, reducing the active-metabolite fraction from 0.120",
    "in extensive metabolizers to 0.071 in intermediate and 0.034 in poor",
    "metabolizers. The PD endpoint is the P2Y12 reaction unit (PRU) from",
    "the VerifyNow P2Y12 assay, modeled as a turnover pool in which",
    "clopidogrel H4 stimulates the fractional turnover rate Kout through a",
    "sigmoid Emax function; the drug-free baseline is Kin / Kout = 212.67",
    "PRU. Inter-individual variability is reported on Vc, CLc, Qc (with a",
    "correlation between Vc and CLc), fm1, fm2, Vp2, CLm2, Kin, and Emax.",
    "Residual error is proportional for all three analytes and additive for",
    "PRU. PK and PD were fitted sequentially: all PK parameters were fixed",
    "to their final PK-model estimates before the PD parameters were",
    "estimated.")
  reference <- "Jung YS, Jin BH, Park MS, Kim CO, Chae D. Population pharmacokinetic-pharmacodynamic modeling of clopidogrel for dose regimen optimization based on CYP2C19 phenotypes: A proof of concept study. CPT Pharmacometrics Syst Pharmacol. 2024;13(1):29-40. doi:10.1002/psp4.13053"
  vignette <- "Jung_2024_clopidogrel"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot                = list(analyte = "clopidogrel",                  units = "mg",  specimen = "administration site", verified = TRUE),
    liver                = list(analyte = "clopidogrel",                  units = "mg",  specimen = "tissue",              verified = TRUE),
    central              = list(analyte = "clopidogrel",                  units = "mg",  specimen = "plasma",              verified = TRUE),
    peripheral1          = list(analyte = "clopidogrel",                  units = "mg",  specimen = "plasma",              verified = TRUE),
    central_h4           = list(analyte = "clopidogrel H4 active thiol",  units = "mg",  specimen = "plasma",              verified = TRUE),
    central_cloca        = list(analyte = "clopidogrel carboxylic acid",  units = "mg",  specimen = "plasma",              verified = TRUE),
    peripheral1_cloca    = list(analyte = "clopidogrel carboxylic acid",  units = "mg",  specimen = "plasma",              verified = TRUE),
    pru                  = list(analyte = "P2Y12 reaction unit",          units = "PRU", specimen = "whole blood",         verified = TRUE)
  )

  covariateData <- list(
    CYP2C19_IM = list(
      description        = "CYP2C19 intermediate-metabolizer phenotype indicator: 1 = subject is a CYP2C19 intermediate metabolizer; 0 = otherwise.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 together with CYP2C19_PM = 0, i.e. the extensive-metabolizer (EM) pool.",
      notes              = "Time-fixed germline-genotype-derived phenotype determined from a blood sample taken on the day of admission (Jung 2024 Methods 'Plasma assay for PK, PD, and CYP phenotyping'). Jung 2024 Results 'Subjects and data' gives the genotype-to-phenotype map explicitly: *1/*1 and *1/*17 = EM; *1/*2, *1/*3, *2/*17 and *3/*17 = IM; *2/*2, *2/*3 and *3/*3 = PM. Cohort counts EM 17 / IM 15 / PM 4. Enters as an ADDITIVE shift on the logit scale of both nested metabolized fractions: logit(fm1) = logit(0.125) + e_cyp2c19_im_logitfm1 * CYP2C19_IM + e_cyp2c19_pm_logitfm1 * CYP2C19_PM, and likewise for fm2. The paper does not distinguish CYP2C19 ultrarapid metabolizers (UM) from EM -- Jung 2024 Discussion limitation 1 states *1/*17 subjects were pooled into EM -- so the reference category is an EM/UM pool exactly as in the register entry.",
      source_name        = "CYP2C19 phenotype (IM)"
    ),
    CYP2C19_PM = list(
      description        = "CYP2C19 poor-metabolizer phenotype indicator: 1 = subject is a CYP2C19 poor metabolizer; 0 = otherwise.",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 together with CYP2C19_IM = 0, i.e. the extensive-metabolizer (EM) pool.",
      notes              = "Time-fixed germline-genotype-derived phenotype; see the CYP2C19_IM entry for the genotype-to-phenotype map and the logit-additive covariate form. n = 4 PM subjects (*2/*2, *2/*3, *3/*3) out of 36. The PM shift is applied to both fm1 (e_cyp2c19_pm_logitfm1 = -0.996) and fm2 (e_cyp2c19_pm_logitfm2 = -2.432), reproducing Jung 2024 Table 3 exactly: fm1 = 0.050, fm2 = 0.678, fmH4 = 0.034, fmcarbo = 0.644, fmothers = 0.322.",
      source_name        = "CYP2C19 phenotype (PM)"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened by stepwise covariate modeling (forward p < 0.01, backward p < 0.001) against the individual etas but not retained. Jung 2024 Results 'PK submodel': 'The only significant covariate was the CYP2C19 phenotype in fm1 and fm2.' Cohort weight 71.32 +/- 8.04 (EM), 72.80 +/- 8.69 (IM), 74.93 +/- 7.68 (PM) kg -- Table 1."
    ),
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened but not retained (see the WT entry). Cohort age 32.35 +/- 6.99 (EM), 31.27 +/- 4.51 (IM), 27.00 +/- 4.90 (PM) years -- Jung 2024 Table 1."
    ),
    HT = list(
      description = "Height",
      units       = "cm",
      type        = "continuous",
      notes       = "Screened but not retained (see the WT entry). Cohort height 172.96 +/- 3.89 (EM), 175.31 +/- 5.19 (IM), 179.00 +/- 8.21 (PM) cm -- Jung 2024 Table 1."
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/dL",
      type        = "continuous",
      notes       = "Screened but not retained (see the WT entry). Cohort albumin 4.62 +/- 0.21 (EM), 4.73 +/- 0.32 (IM), 4.58 +/- 0.35 (PM) g/dL -- Jung 2024 Table 1."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units       = "mg/dL",
      type        = "continuous",
      notes       = "Screened but not retained (see the WT entry). Cohort creatinine 0.91 +/- 0.14 (EM), 0.92 +/- 0.15 (IM), 0.95 +/- 0.07 (PM) mg/dL -- Jung 2024 Table 1."
    ),
    AST = list(
      description = "Aspartate transaminase",
      units       = "IU/L",
      type        = "continuous",
      notes       = "Screened but not retained (see the WT entry). Cohort AST 18.53 +/- 4.90 (EM), 17.93 +/- 3.08 (IM), 17.50 +/- 6.95 (PM) IU/L -- Jung 2024 Table 1."
    ),
    ALT = list(
      description = "Alanine transaminase",
      units       = "IU/L",
      type        = "continuous",
      notes       = "Screened but not retained (see the WT entry). Cohort ALT 20.29 +/- 12.62 (EM), 15.40 +/- 5.64 (IM), 15.00 +/- 7.75 (PM) IU/L -- Jung 2024 Table 1."
    ),
    ALKPHOS = list(
      description = "Alkaline phosphatase",
      units       = "IU/L",
      type        = "continuous",
      notes       = "Screened but not retained (see the WT entry). Cohort ALP 69.06 +/- 14.74 (EM), 64.40 +/- 12.88 (IM), 72.00 +/- 8.87 (PM) IU/L -- Jung 2024 Table 1."
    ),
    GGT = list(
      description = "Gamma-glutamyl transferase",
      units       = "IU/L",
      type        = "continuous",
      notes       = "Screened but not retained (see the WT entry). Cohort gamma-GT 23.12 +/- 13.79 (EM), 22.20 +/- 12.19 (IM), 22.50 +/- 10.97 (PM) IU/L -- Jung 2024 Table 1."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units       = "mg/dL",
      type        = "continuous",
      notes       = "Screened but not retained (see the WT entry). Cohort total bilirubin 0.69 +/- 0.11 (EM), 0.85 +/- 0.38 (IM), 0.80 +/- 0.50 (PM) mg/dL -- Jung 2024 Table 1."
    ),
    BUN = list(
      description = "Blood urea nitrogen",
      units       = "mg/dL",
      type        = "continuous",
      notes       = "Screened but not retained (see the WT entry). Cohort BUN 13.11 +/- 3.60 (EM), 12.87 +/- 3.83 (IM), 12.68 +/- 3.17 (PM) mg/dL -- Jung 2024 Table 1."
    ),
    TPROT = list(
      description = "Total protein",
      units       = "g/dL",
      type        = "continuous",
      notes       = "Screened but not retained (see the WT entry). Cohort total protein 7.02 +/- 0.35 (EM), 7.19 +/- 0.41 (IM), 6.98 +/- 0.30 (PM) g/dL -- Jung 2024 Table 1."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 36L,
    n_studies      = 1L,
    n_observations = "503 clopidogrel, 239 clopidogrel H4, 504 clopidogrel carboxylic acid, and 280 PRU observations (Jung 2024 Results 'Subjects and data'). Proportions below the quantification limit were 27.0% (clopidogrel), 11.7% (H4) and 3.8% (carboxylic acid); no PRU value was BQL. BQL records were handled with the M4 method (simultaneous categorical modeling of BQL observations constrained to (0, LLOQ)).",
    age_range      = "19-55 years by protocol; observed means 32.35 +/- 6.99 (EM), 31.27 +/- 4.51 (IM), 27.00 +/- 4.90 (PM) years (Table 1).",
    weight_range   = "55-90 kg by protocol; observed means 71.32 +/- 8.04 (EM), 72.80 +/- 8.69 (IM), 74.93 +/- 7.68 (PM) kg (Table 1).",
    sex_female_pct = 0,
    race_ethnicity = c(Asian = 100),
    disease_state  = "Healthy male volunteers. Body mass index 18.5-27.0 kg/m^2. Excluded: clinically significant pulmonary, cardiovascular, hepatobiliary, neurological, endocrine or immune disease; current smokers; gastrointestinal disease or surgery affecting absorption; clinically significant bleeding history; and screening PRU outside +/- 10% of the normal-range limits.",
    dose_range     = "Clopidogrel 75 mg tablet orally once daily for 7 days, with no loading dose (aligned with non-acute-coronary-syndrome practice). PK samples at predose and 0.33, 0.67, 1, 1.5, 2, 3, 4, 6, 8, 12 and 24 h; H4 sampled only at predose and 0.33, 0.67, 1, 2, 4 and 6 h. PRU measured predose on days 1, 3, 5 and 7 and at 1, 4, 6 and 24 h after the day-7 dose, in duplicate with the arithmetic mean used as the endpoint.",
    regions        = "Single centre, Severance Hospital, Seoul, South Korea (October 2019 to April 2020).",
    genotype_distribution = "CYP2C19 phenotype: EM n = 17 (*1/*1, *1/*17), IM n = 15 (*1/*2, *1/*3, *2/*17, *3/*17), PM n = 4 (*2/*2, *2/*3, *3/*3). Ultrarapid metabolizers were not distinguished from extensive metabolizers (Jung 2024 Discussion limitation 1).",
    notes          = "Prospective phase I trial, IRB 4-2019-0740, NCT04171687. Estimation in MonolixSuite 2021R2 with the SAEM algorithm; simulations in Simulx 2021R2. Assay LLOQs 0.09 ng/mL (clopidogrel), 0.3 ng/mL (H4) and 50 ng/mL (carboxylic acid); calibration ranges 0.09-10, 0.3-50 and 50-5000 ng/mL respectively. PRU by the VerifyNow P2Y12 assay (Accumetrics) on citrated whole blood analysed within 4 h of collection. Baseline PRU was 207.41 +/- 23.60 (EM), 199.83 +/- 28.60 (IM) and 207.50 +/- 33.48 (PM). Covariates screened and rejected are recorded in `covariatesDataExcluded`."
  )

  ini({
    # ==================================================================
    # CLOPIDOGREL (PARENT) STRUCTURAL PARAMETERS
    # ==================================================================
    # All values are the population means of Jung 2024 Table 2
    # ("Parameter estimates of the final PK-PD model"); the parenthesised
    # number in that table is the %RSE, not an IIV. Volumes are in L and
    # clearances in L/h. The parent's apparent volumes and clearances are
    # large because clopidogrel undergoes near-complete first-pass
    # metabolism -- only ~8% of the hepatic flux reaches systemic plasma
    # as parent -- and because no absolute bioavailability term is
    # identifiable from oral-only data.

    ltlag <- log(0.196)
    label("Absorption lag time (Tlag, h)")                                     # Jung 2024 Table 2: T_lag = 0.196 h (%RSE 4.59)

    lka <- log(19.64)
    label("First-order absorption rate constant depot -> liver (ka, 1/h)")     # Jung 2024 Table 2: k_a = 19.64 1/h (%RSE 28.0); Results 'PK submodel' confirms "corresponding to an absorption half-life of 2.13 min" (log(2)/19.64 h = 2.12 min)

    lvc <- log(1463.92)
    label("Clopidogrel central volume (Vc, L; the hepatic volume VH is set equal to it)")  # Jung 2024 Table 2: V_c (= V_H) = 1463.92 L (%RSE 7.56). Methods 'PK model': "Because of parameter unidentifiability, we set the hepatic volume to be the same as the central volume".

    lvp <- log(2823.98)
    label("Clopidogrel peripheral volume (Vp, L)")                             # Jung 2024 Table 2: V_p = 2823.98 L (%RSE 7.72)

    lcl <- log(9257.28)
    label("Clopidogrel hepatic metabolic clearance out of the liver compartment (CLc, L/h)")  # Jung 2024 Table 2: CL_c = 9257.28 L/h (%RSE 6.92). This is the parent's ONLY elimination pathway -- Methods 'PK model': "We further assumed that all parent drugs ultimately underwent metabolism in the hepatic compartment."

    lq_liver <- log(845.70)
    label("Inter-compartmental clearance liver <-> central (Qc, L/h)")          # Jung 2024 Table 2: Q_c = 845.70 L/h (%RSE 7.04); Table 2 footnote defines Q_c as "intercompartmental clearance between hepatic and central compartment"

    lq <- log(587.93)
    label("Inter-compartmental clearance central <-> peripheral1 (Qp, L/h)")    # Jung 2024 Table 2: Q_p = 587.93 L/h (%RSE 9.23); Table 2 footnote defines Q_p as "intercompartmental clearance between central and peripheral compartment"

    # ==================================================================
    # NESTED METABOLIZED FRACTIONS (LOGIT SCALE)
    # ==================================================================
    # Jung 2024 Methods 'PK model' defines two nested fractions with
    #   fmH4     = fm1 * fm2
    #   fmcarbo  = (1 - fm1) * fm2
    #   fmothers = 1 - fm2
    # i.e. fm2 is the fraction of the hepatic metabolic flux that goes to
    # the two MEASURED metabolites and fm1 is the share of that flux
    # routed to the active H4 thiol. Both typical values were FIXED, not
    # estimated -- Table 2 flags them "0.125 FIX" and "0.960 FIX", and
    # Methods 'PK model' explains why: "Because of parameter
    # unidentifiability, we ... fixed the typical fm1 and fm2 values to
    # 0.125 and 0.96, respectively, based on the known fractions
    # metabolized to clopidogrel H4 and clopidogrel carboxylic acid,
    # which were ~10-15% and 85%, respectively."
    #
    # Both are carried on the LOGIT scale because Methods 'PK-PD model'
    # states: "For fm1 and fm2, whose values were constrained between 0
    # and 1, IIVs were applied using a logit function", with the
    # transformation written out as log(Pi/(1-Pi)) = log(TVP/(1-TVP)) +
    # eta_i. The CYP2C19 shifts below are therefore additive on this
    # same logit scale.
    #   logit(0.125) = log(0.125/0.875) = -1.945910
    #   logit(0.960) = log(0.960/0.040) = +3.178054

    logitfm1 <- fixed(log(0.125 / (1 - 0.125)))
    label("Logit of the typical fraction of the measured metabolic flux routed to clopidogrel H4 (logit(fm1), unitless)")  # Jung 2024 Table 2: fm1 = 0.125 FIX -> logit(0.125) = -1.94591

    logitfm2 <- fixed(log(0.960 / (1 - 0.960)))
    label("Logit of the typical fraction of the hepatic metabolic flux routed to the two measured metabolites (logit(fm2), unitless)")  # Jung 2024 Table 2: fm2 = 0.960 FIX -> logit(0.960) = 3.178054

    # ==================================================================
    # CYP2C19 PHENOTYPE COVARIATE EFFECTS (ADDITIVE ON THE LOGIT SCALE)
    # ==================================================================
    # Jung 2024 Table 2 reports these as the rows "fm1 ~ IM", "fm1 ~ PM",
    # "fm2 ~ IM" and "fm2 ~ PM". The EM phenotype is the reference (both
    # indicators zero). Arithmetic check against Table 3, which tabulates
    # the resulting fractions:
    #   fm1 IM: logit(0.125) - 0.450 = -2.39591 -> 0.0835 -> Table 3 0.083
    #   fm1 PM: logit(0.125) - 0.996 = -2.94191 -> 0.0501 -> Table 3 0.050
    #   fm2 IM: logit(0.960) - 1.428 = +1.75005 -> 0.8520 -> Table 3 0.852
    #   fm2 PM: logit(0.960) - 2.432 = +0.74605 -> 0.6783 -> Table 3 0.678
    # and the derived fractions then reproduce Table 3 exactly:
    #   EM fmH4 0.120 / fmcarbo 0.840 / fmothers 0.040
    #   IM fmH4 0.071 / fmcarbo 0.781 / fmothers 0.148
    #   PM fmH4 0.034 / fmcarbo 0.644 / fmothers 0.322

    e_cyp2c19_im_logitfm1 <- -0.450
    label("Additive logit-scale shift of fm1 for CYP2C19 intermediate metabolizers (unitless)")  # Jung 2024 Table 2: fm1 ~ IM = -0.450 (%RSE 29.8)

    e_cyp2c19_pm_logitfm1 <- -0.996
    label("Additive logit-scale shift of fm1 for CYP2C19 poor metabolizers (unitless)")          # Jung 2024 Table 2: fm1 ~ PM = -0.996 (%RSE 17.6)

    e_cyp2c19_im_logitfm2 <- -1.428
    label("Additive logit-scale shift of fm2 for CYP2C19 intermediate metabolizers (unitless)")  # Jung 2024 Table 2: fm2 ~ IM = -1.428 (%RSE 29.8)

    e_cyp2c19_pm_logitfm2 <- -2.432
    label("Additive logit-scale shift of fm2 for CYP2C19 poor metabolizers (unitless)")          # Jung 2024 Table 2: fm2 ~ PM = -2.432 (%RSE 24.2)

    # ==================================================================
    # CLOPIDOGREL H4 (ACTIVE THIOL METABOLITE) -- ONE COMPARTMENT
    # ==================================================================

    lvc_h4 <- log(51.45)
    label("Clopidogrel H4 central volume (Vm1, L)")                            # Jung 2024 Table 2: V_m1 = 51.45 L (%RSE 9.81)

    lcl_h4 <- log(74.25)
    label("Clopidogrel H4 clearance (CLm1, L/h)")                              # Jung 2024 Table 2: CL_m1 = 74.25 L/h (%RSE 8.87)

    # ==================================================================
    # CLOPIDOGREL CARBOXYLIC ACID (INACTIVE METABOLITE) -- TWO COMPARTMENTS
    # ==================================================================

    lvc_cloca <- log(17.34)
    label("Clopidogrel carboxylic acid central volume (Vm2, L)")                # Jung 2024 Table 2: V_m2 = 17.34 L (%RSE 3.54)

    lvp_cloca <- log(51.89)
    label("Clopidogrel carboxylic acid peripheral volume (Vp2, L)")             # Jung 2024 Table 2: V_p2 = 51.89 L (%RSE 6.48)

    lcl_cloca <- log(7.248)
    label("Clopidogrel carboxylic acid clearance (CLm2, L/h)")                  # Jung 2024 Table 2: CL_m2 = 7.248 L/h (%RSE 3.53)

    lq_cloca <- log(4.476)
    label("Clopidogrel carboxylic acid inter-compartmental clearance (Qm2, L/h)")  # Jung 2024 Table 2: Q_m2 = 4.476 L/h (%RSE 3.59)

    # ==================================================================
    # PRU TURNOVER PD SUBMODEL
    # ==================================================================
    # Sequentially fitted: Jung 2024 Methods 'PK-PD model' states "we
    # fixed all PK parameters to the estimates obtained from the final PK
    # model" before estimating the PD parameters. The PK values above are
    # nonetheless left estimable here because they ARE the final PK-model
    # estimates -- the fixing was a sequential-estimation device, not a
    # statement that the PK values were externally imposed.
    #
    # Kout is entered as 0.00576 1/h, not the rounded "0.006" printed in
    # Table 2. Jung 2024 Results 'PD submodel' gives the unrounded value
    # explicitly ("Kout was estimated to be 0.00576 h-1 (~0.006)") and
    # only the unrounded value reproduces the two derived quantities the
    # same paragraph reports:
    #   baseline PRU = Kin / Kout = 1.225 / 0.00576 = 212.67  (paper: 212.67)
    #   turnover half-life = log(2) / 0.00576 = 120.3 h = 5.01 days (paper: "a half-life of 5 days")
    # Using 0.006 would give 204.17 and 4.81 days, matching neither.

    lkin <- log(1.225)
    label("Zero-order PRU production rate (Kin, PRU/h)")                       # Jung 2024 Table 2: K_in = 1.225 (%RSE 3.81)

    lkout <- log(0.00576)
    label("First-order fractional PRU turnover rate (Kout, 1/h)")              # Jung 2024 Results 'PD submodel': "Kout was estimated to be 0.00576 h-1 (~0.006)" -- Table 2 prints the rounded 0.006

    lemax <- log(57.84)
    label("Maximum fractional stimulation of PRU turnover by clopidogrel H4 (Emax, unitless multiplier on Kout)")  # Jung 2024 Table 2: E_max = 57.84 (%RSE 12.7). Results 'PD submodel': "Emax was estimated to be 57.84, implying a maximum reduction in PRU at a steady state of 1.7% (1/(1 + 57.84))" -- i.e. the steady-state PRU floor is baseline / (1 + Emax) = 212.67 / 58.84 = 3.61 PRU.

    lec50 <- log(67.32)
    label("Clopidogrel H4 concentration giving half-maximal stimulation of PRU turnover (EC50, ng/mL)")  # Jung 2024 Table 2: EC_50 = 67.32 ng/mL (%RSE 12.0)

    lhill <- log(1.851)
    label("Hill coefficient of the H4-on-PRU sigmoid (unitless)")              # Jung 2024 Table 2: Hill = 1.851 (%RSE 7.46)

    # ==================================================================
    # INTER-INDIVIDUAL VARIABILITY
    # ==================================================================
    # Jung 2024 Table 2's second column is headed "Standard deviation",
    # so the tabulated values are omega (SD), not omega^2. The paper
    # confirms this reading in prose for Emax -- Results 'PD submodel':
    # "The standard deviation (omega) associated with the IIV of Emax was
    # estimated to be 0.501". nlmixr2's `eta ~ value` takes a VARIANCE, so
    # every value below is the squared Table 2 entry.
    #
    # IIV is exponential on every parameter except fm1 and fm2, which are
    # logit-additive (Methods 'PK-PD model': "Pi = TVP * exp(eta_i) ... or
    # log(Pi/(1-Pi)) = log(TVP/(1-TVP)) + eta_i for fm1 and fm2").
    #
    # Results 'PK submodel' lists exactly which parameters carry IIV:
    # "IIV was incorporated into Vc(=VH), CLc, Qc, fm1, fm2, Vp2, and
    # CLm2", and Results 'PD submodel' adds "IIV was incorporated into
    # Kin and Emax".
    #
    # Vc and CLc are correlated. Table 2 reports Cor(Vc, CLc) = 0.702
    # (%RSE 21.2) as a CORRELATION, so the covariance must be rebuilt:
    #   cov = 0.702 * 0.331 * 0.343 = 0.0797  (see the block comment below)

    etalvc + etalcl ~ c(0.109561,
                        0.079700, 0.117649)                                    # Jung 2024 Table 2: SD(Vc) = 0.331 -> var 0.331^2 = 0.109561; SD(CLc) = 0.343 -> var 0.343^2 = 0.117649; Cor(Vc, CLc) = 0.702 -> cov = 0.702 * 0.331 * 0.343 = 0.0797

    etalq_liver ~ 0.094249                                                     # Jung 2024 Table 2: SD(Qc) = 0.307 -> var 0.307^2 = 0.094249
    etalogitfm1 ~ 0.065025                                                     # Jung 2024 Table 2: SD(fm1) = 0.255 (on the logit scale) -> var 0.255^2 = 0.065025
    etalogitfm2 ~ 1.276900                                                     # Jung 2024 Table 2: SD(fm2) = 1.130 (on the logit scale) -> var 1.130^2 = 1.2769
    etalvp_cloca ~ 0.081225                                                 # Jung 2024 Table 2: SD(Vp2) = 0.285 -> var 0.285^2 = 0.081225
    etalcl_cloca ~ 0.015129                                                 # Jung 2024 Table 2: SD(CLm2) = 0.123 -> var 0.123^2 = 0.015129
    etalkin ~ 0.014400                                                         # Jung 2024 Table 2: SD(Kin) = 0.120 -> var 0.120^2 = 0.0144
    etalemax ~ 0.251001                                                        # Jung 2024 Table 2: SD(Emax) = 0.501 -> var 0.501^2 = 0.251001; Results 'PD submodel' calls 0.501 the standard deviation (omega), confirming the column is an SD

    # ==================================================================
    # RESIDUAL ERROR
    # ==================================================================
    # Jung 2024 Methods 'PK-PD model' writes the combined form
    # Y_ij = F_ij * (1 + eps_Pro,ij) + eps_Add,ij, but Results 'PK
    # submodel' states which component was actually retained per
    # endpoint: "A proportional residual error model was selected for all
    # end points, with estimated magnitudes of 35.7%, 36.3%, and 20.9%
    # for clopidogrel, clopidogrel H4, and clopidogrel carboxylic acid",
    # and Results 'PD submodel' reports for PRU that "The additive
    # residual error model yielded an estimated error magnitude of
    # 13.69". So the three PK endpoints are purely proportional and the
    # PD endpoint is purely additive; no endpoint carries both.

    propSd <- 0.357
    label("Clopidogrel proportional residual SD (fraction)")                   # Jung 2024 Table 2: sigma_clopidogrel = 0.357 (%RSE 2.68); Results 'PK submodel' quotes it as 35.7%

    propSd_h4 <- 0.363
    label("Clopidogrel H4 proportional residual SD (fraction)")                # Jung 2024 Table 2: sigma_H4 = 0.363 (%RSE 3.91); Results 'PK submodel' quotes it as 36.3%

    propSd_cloca <- 0.209
    label("Clopidogrel carboxylic acid proportional residual SD (fraction)")   # Jung 2024 Table 2: sigma_carbo = 0.209 (%RSE 2.11); Results 'PK submodel' quotes it as 20.9%

    addSd_PRU <- 13.69
    label("PRU additive residual SD (PRU)")                                    # Jung 2024 Table 2: sigma_PRU = 13.69 (%RSE 2.66); Table 2 footnote: "sigma_PRU, additive error about PRU"
  })

  model({
    # ------------------------------------------------------------------
    # Individual parent PK parameters. VH = Vc exactly (Jung 2024
    # Methods 'PK model'), so the hepatic volume is not a separate
    # parameter and the individual Vc eta propagates to both.
    # ------------------------------------------------------------------
    tlag    <- exp(ltlag)
    ka      <- exp(lka)
    vc      <- exp(lvc + etalvc)
    vh      <- vc
    vp      <- exp(lvp)
    cl      <- exp(lcl + etalcl)
    q_liver <- exp(lq_liver + etalq_liver)
    q       <- exp(lq)

    # ------------------------------------------------------------------
    # Nested metabolized fractions. CYP2C19 phenotype enters additively
    # on the logit scale alongside the logit-additive eta, so every
    # individual fm1 and fm2 is guaranteed to stay inside (0, 1)
    # regardless of the phenotype / eta combination.
    # ------------------------------------------------------------------
    # Linear predictor first, inverse logit second, following the
    # `logitfm_ind` / `fm` form of the founding logitfm model
    # (Mitra_2026_ziftomenib.R).
    logitfm1_ind <- logitfm1 + e_cyp2c19_im_logitfm1 * CYP2C19_IM +
                    e_cyp2c19_pm_logitfm1 * CYP2C19_PM + etalogitfm1
    logitfm2_ind <- logitfm2 + e_cyp2c19_im_logitfm2 * CYP2C19_IM +
                    e_cyp2c19_pm_logitfm2 * CYP2C19_PM + etalogitfm2
    fm1 <- 1 / (1 + exp(-logitfm1_ind))
    fm2 <- 1 / (1 + exp(-logitfm2_ind))

    # ------------------------------------------------------------------
    # Metabolite PK parameters
    # ------------------------------------------------------------------
    vc_h4       <- exp(lvc_h4)
    cl_h4       <- exp(lcl_h4)
    vc_cloca <- exp(lvc_cloca)
    vp_cloca <- exp(lvp_cloca + etalvp_cloca)
    cl_cloca <- exp(lcl_cloca + etalcl_cloca)
    q_cloca  <- exp(lq_cloca)

    # Molar-mass ratio scaling factors applied to the two metabolite
    # formation fluxes. Jung 2024 Methods 'PK model': "Given the
    # molecular weights of clopidogrel, clopidogrel H4, and clopidogrel
    # carboxylic acid, viz. 321.82, 355.83, and 307.80 g/mol,
    # respectively, we applied scaling factors of 1.106 and 0.956 to the
    # formation rates of clopidogrel H4 and clopidogrel carboxylic acid".
    # Check: 355.83 / 321.82 = 1.1057 and 307.80 / 321.82 = 0.9564, so
    # each factor converts a parent-mass flux into the equivalent
    # metabolite-mass flux at 1:1 molar stoichiometry.
    mpr_h4       <- 1.106
    mpr_cloca <- 0.956

    # ------------------------------------------------------------------
    # PD parameters
    # ------------------------------------------------------------------
    kin  <- exp(lkin + etalkin)
    kout <- exp(lkout)
    emax <- exp(lemax + etalemax)
    ec50 <- exp(lec50)
    hill <- exp(lhill)

    # ------------------------------------------------------------------
    # ODE system -- Jung 2024 Results 'PK-PD modeling', page 34, the
    # eight printed differential equations, transcribed term by term.
    # Source symbols map as: Aa -> depot, AH -> liver, Ac -> central,
    # Ap -> peripheral1, Am1 -> central_h4, Am2 -> central_cloca,
    # Ap2 -> peripheral1_cloca.
    #
    #   dAa/dt  = -ka * Aa
    #   dAH/dt  =  ka * Aa - (CLc + Qc)/VH * AH + Qc/Vc * Ac
    #   dAc/dt  =  Qc/VH * AH - (Qc + Qp)/Vc * Ac + Qp/Vp * Ap
    #   dAp/dt  =  Qp/Vc * Ac - Qp/Vp * Ap
    #   dAm1/dt =  1.106 * fm1 * fm2 * CLc/VH * AH - CLm1/Vm1 * Am1
    #   dAm2/dt =  0.956 * (1 - fm1) * fm2 * CLc/VH * AH
    #              - (CLm2 + Qm2)/Vm2 * Am2 + Qm2/Vp2 * Ap2
    #   dAp2/dt =  Qm2/Vm2 * Am2 - Qm2/Vp2 * Ap2
    #
    # Note that the fm1 * fm2 and (1 - fm1) * fm2 shares are taken out of
    # the SAME hepatic metabolic flux CLc/VH * AH that already leaves the
    # liver state; the remaining 1 - fm2 share goes to unmeasured "other"
    # metabolites and is therefore not carried as a state.
    # ------------------------------------------------------------------
    d/dt(depot)       <- -ka * depot
    d/dt(liver)       <-  ka * depot -
                          (cl + q_liver) / vh * liver +
                          q_liver / vc * central
    d/dt(central)     <-  q_liver / vh * liver -
                          (q_liver + q) / vc * central +
                          q / vp * peripheral1
    d/dt(peripheral1) <-  q / vc * central - q / vp * peripheral1

    d/dt(central_h4)  <-  mpr_h4 * fm1 * fm2 * cl / vh * liver -
                          cl_h4 / vc_h4 * central_h4

    d/dt(central_cloca)     <- mpr_cloca * (1 - fm1) * fm2 * cl / vh * liver -
                                  (cl_cloca + q_cloca) / vc_cloca * central_cloca +
                                  q_cloca / vp_cloca * peripheral1_cloca
    d/dt(peripheral1_cloca) <- q_cloca / vc_cloca * central_cloca -
                                  q_cloca / vp_cloca * peripheral1_cloca

    alag(depot) <- tlag

    # ------------------------------------------------------------------
    # PRU turnover. Jung 2024 Results 'PK-PD modeling', page 35:
    #   dPRU/dt = Kin - Kout * (1 + Emax * Cm1^Hill / (EC50^Hill +
    #             Cm1^Hill)) * PRU,  PRU(0) = Kin / Kout
    # The ODE state is the lowercase `pru` (library convention: ODE
    # states are lowercase) and the uppercase `PRU` below is the
    # observation output, following the registered `anc` / `ANC` pairing.
    #
    # The H4 concentration driving the sigmoid is written INLINE as
    # (central_h4 / vc_h4 * 1000) rather than routed through the named
    # `Cc_h4` intermediate defined below: in the nlmixr2 model-function
    # form a named intermediate that references an ODE state can silently
    # evaluate to zero inside a d/dt() expression, which would delete the
    # entire drug effect without any error. The *1000 is the mg/L ->
    # ng/mL conversion described at the observation block below, and it
    # matters here because EC50 is reported in ng/mL.
    # ------------------------------------------------------------------
    d/dt(pru) <- kin -
                 kout * (1 + emax * (central_h4 / vc_h4 * 1000)^hill /
                             (ec50^hill + (central_h4 / vc_h4 * 1000)^hill)) * pru
    pru(0)    <- kin / kout

    # ------------------------------------------------------------------
    # Observations. Doses are in mg and volumes in L, so amount / volume
    # returns mg/L; the * 1000 expresses each concentration in the source
    # paper's ng/mL. (Jung 2024's own dataset carried amounts in ug so
    # that ug/L = ng/mL fell out directly; scaling by 1000 here keeps the
    # library convention of mg dosing while reproducing the paper's
    # concentration scale, on which EC50 = 67.32 ng/mL is defined.)
    # ------------------------------------------------------------------
    Cc          <- central       / vc       * 1000
    Cc_h4       <- central_h4    / vc_h4    * 1000
    Cc_cloca    <- central_cloca / vc_cloca * 1000

    # PRU needs no unit conversion -- the state is already in the assay's
    # own P2Y12-reaction-unit scale. The uppercase alias is the
    # observation output name (see the ODE block above).
    PRU <- pru

    Cc          ~ prop(propSd)
    Cc_h4       ~ prop(propSd_h4)
    Cc_cloca    ~ prop(propSd_cloca)
    PRU         ~ add(addSd_PRU)
  })
}
