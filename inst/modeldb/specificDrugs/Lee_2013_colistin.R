Lee_2013_colistin <- function() {
  description <- "One-compartment population PK model of colistin in adult burn-ICU patients receiving colistimethate sodium (CMS) as a 30-minute IV infusion every 12 hours, with first-order CMS-to-colistin conversion (Lee 2013). Apparent CL and Vc of colistin are scaled inversely by the relative fraction of CMS converted to colistin (RFM = 1 - theta4 * (CRCL/128)); the CMS-to-colistin turnover rate constant TR is reduced in patients with clinically-evident peripheral edema (TR = theta3 - theta5 * DIS_EDEMA)."
  reference <- "Lee J, Han S, Jeon S, Hong T, Song W, Woo H, Yim D-S. Population pharmacokinetic analysis of colistin in burn patients. Antimicrob Agents Chemother. 2013;57(5):2421-2427. doi:10.1128/AAC.00271-13"
  vignette <- "Lee_2013_colistin"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    CRCL = list(
      description        = "Cockcroft-Gault creatinine clearance (raw, not BSA-normalized)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source column CLCR. Computed by the Cockcroft-Gault equation from age, body weight, and serum creatinine, expressed in raw mL/min (NOT BSA-normalized to mL/min/1.73 m^2). Stored under the canonical CRCL column per inst/references/covariate-columns.md (CRCL accepts raw mL/min when the source paper does not apply BSA normalization, with the per-model description recording the assay form). Reference value 128 mL/min (cohort median, Lee 2013 Table 1; range 22.6-309 mL/min). The CRCL effect enters the relative fraction of CMS converted to colistin: RFM = 1 - theta4 * (CRCL/128); the apparent CL and Vc of colistin are then divided by RFM (effective CL = (CL/fm*) / RFM; effective Vc = (Vc/fm*) / RFM).",
      source_name        = "CLCR"
    ),
    DIS_EDEMA = list(
      description        = "Clinically-evident peripheral edema indicator (puffy face and pitting edema in the legs, per clinical exam on Day 1 of CMS administration)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no clinical edema)",
      notes              = "Source column EDEMA. Stored under canonical DIS_EDEMA per inst/references/covariate-columns.md (DIS_EDEMA entry added in this PR alongside the Lee 2013 extraction; operator-resolved sidecar request-001 / response-001 on the canonical name). Lee 2013 records edema as a single clinical diagnosis on Day 1 of CMS administration (Table 1 footnote d). 18 of 50 enrolled patients were edematous at baseline. Time-fixed per subject in this paper (no serial reassessment). The covariate enters as an additive linear deviation on the CMS-to-colistin turnover rate constant TR: TR = theta3 - theta5 * DIS_EDEMA, reducing TR from 0.796 (non-edematous) to 0.371 h^-1 (edematous).",
      source_name        = "EDEMA"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 50L,
    n_studies      = 1L,
    age_range      = "26-80 years",
    age_median     = "48 years (mean, SD 13)",
    weight_range   = "50-98 kg",
    weight_median  = "65.8 kg (mean, SD 10.3)",
    sex_female_pct = 22,
    race_ethnicity = "Not reported (single-center South Korean burn-ICU cohort)",
    disease_state  = "Adult burn patients with 4-85% total body surface area affected (median TBSA 50.5%, SD 21.8); treated for nosocomial multidrug-resistant Gram-negative bacterial infections (Acinetobacter baumannii, Pseudomonas aeruginosa, etc.) in a burn intensive care unit.",
    dose_range     = "150 mg CMS as colistin base activity (CBA) by IV infusion over 30 min, every 12 h. PK sampling started >=3 days after the first dose (steady-state assumption).",
    regions        = "South Korea (single center: Burn Intensive Care Unit, Hangang Sacred Heart Hospital, Hallym University Medical Center)",
    abs_i          = "ABSI 9.82 (SD 2.34, range 5-14)",
    days_post_burn = "15.5 days median (SD 10.4, range 3-58) from burn injury to first CMS dose",
    renal_function = "Cockcroft-Gault CRCL 128 mL/min (SD 75.2, range 22.6-309); raw, not BSA-normalized",
    crrt_pct       = "17/50 (34%) on continuous renal replacement therapy at baseline",
    sepsis_pct     = "29/50 (58%) with sepsis at baseline",
    edema_pct      = "18/50 (36%) clinically edematous at baseline",
    albumin        = "Serum albumin 2.5 g/dL (SD 0.3, range 1.9-3.1); low across cohort",
    notes          = "Baseline demographics per Lee 2013 Table 1. 50 burn-ICU adults enrolled June 2010 - May 2011. Burn etiology: 39 flame, 5 electrical, 3 scalding, 2 chemical, 1 contact. All but one patient had TBSA 11-85% (the one outlier had a 4% electrical burn). Patients < 18 years, pregnant / breastfeeding, or allergic to CMS / colistin were excluded. The cohort spans a wide renal-function range, which drives the CRCL covariate effect on RFM in the final model. CRRT was tested as a covariate but not retained -- individual CL and RFM in CRRT vs non-CRRT patients did not differ in t tests (Lee 2013 Discussion paragraph 6)."
  )

  ini({
    # Structural parameters at the typical non-edematous patient with median renal
    # function (CRCL = 128 mL/min) -- Lee 2013 Table 2 final-model column.
    # The reported CL/fm* and Vc/fm* are the apparent CL and Vc scaled by the
    # theoretical fraction of CMS converted to colistin at CRCL = 0 (RFM = 1);
    # the effective CL and Vc are obtained by dividing by RFM (see model() block).
    # At CRCL = 128 mL/min, RFM = 0.787 and effective CL = 8.49 / 0.787 = 10.8 L/h
    # and effective Vc = 81.1 / 0.787 = 103 L (Lee 2013 Results paragraph 4).
    lcl <- log(8.49); label("Apparent CL/fm* at RFM=1 (L/h)")                     # Lee 2013 Table 2: CL/fm* (theta1) = 8.49 L/h
    lvc <- log(81.1); label("Apparent Vc/fm* at RFM=1 (L)")                       # Lee 2013 Table 2: Vc/fm* (theta2) = 81.1 L
    lka <- log(0.796); label("CMS->colistin turnover rate TR (non-edema) (1/h)")  # Lee 2013 Table 2: TR (theta3) = 0.796 1/h for non-edematous

    # Covariate effects -- Lee 2013 Table 2.
    # e_dis_edema_ka is an additive linear deviation on TR (theta5):
    #   TR = exp(lka) - e_dis_edema_ka * DIS_EDEMA
    # e_crcl_cl_vc is the slope of (1 - RFM) per (CRCL/128) (theta4); the same
    # numeric value enters both CL and Vc divisively via 1/RFM:
    #   RFM = 1 - e_crcl_cl_vc * (CRCL/128)
    #   effective CL = (CL/fm*) / RFM, effective Vc = (Vc/fm*) / RFM
    e_dis_edema_ka <- 0.425; label("TR decrease with edema (1/h)")                # Lee 2013 Table 2: theta5 = 0.425 1/h
    e_crcl_cl_vc   <- 0.213; label("Slope of RFM decrease per (CRCL/128)")        # Lee 2013 Table 2: theta4 = 0.213

    # Inter-individual variability (Lee 2013 Table 2 final-model CV%);
    # omega^2 = log(CV^2 + 1) for log-normal etas (Pij = theta_j * exp(eta_ij)
    # per Lee 2013 Methods, "Population PK model development" paragraph 3).
    etalcl ~ 0.131604  # log(0.375^2 + 1); 37.5% CV on CL/fm*
    etalvc ~ 0.057380  # log(0.243^2 + 1); 24.3% CV on Vc/fm*
    etalka ~ 0.372681  # log(0.672^2 + 1); 67.2% CV on TR

    # Combined residual error (Lee 2013 Table 2).
    # The paper reports sigma_add = 99.2 ng/mL; converted to mg/L for unit
    # consistency with units$concentration = "mg/L" -> addSd = 0.0992 mg/L.
    addSd  <- 0.0992; label("Additive residual error (mg/L)")             # Lee 2013 Table 2: sigma_add = 99.2 ng/mL
    propSd <- 0.0672; label("Proportional residual error (fraction)")     # Lee 2013 Table 2: sigma_prop = 6.72%
  })
  model({
    # 1. Derived terms: relative fraction of CMS converted to colistin.
    # RFM = 1 - theta4 * (CRCL / 128); at the cohort median CRCL = 128 mL/min,
    # RFM = 0.787; at CRCL = 0 (anephric), RFM = 1; at CRCL = 309 mL/min (cohort
    # max), RFM ~= 0.485. RFM stays positive across the observed CRCL range.
    rfm <- 1 - e_crcl_cl_vc * (CRCL / 128)

    # 2. Individual PK parameters.
    # CL and Vc are scaled inversely by RFM so the effective CL = (CL/fm*) / RFM
    # and effective Vc = (Vc/fm*) / RFM. TR is reduced by an additive amount
    # when edema is present: TR = theta3 - theta5 * DIS_EDEMA, then log-normal
    # IIV is applied multiplicatively (Pij = theta_j * exp(eta_ij)).
    cl <- exp(lcl + etalcl) / rfm
    vc <- exp(lvc + etalvc) / rfm
    ka <- (exp(lka) - e_dis_edema_ka * DIS_EDEMA) * exp(etalka)

    kel <- cl / vc

    # 3. ODE system: CMS depot -> colistin central -> first-order elimination.
    # CMS was assumed to distribute in a single compartment with first-order
    # rate TR to colistin (CMS concentrations were not measured) per Lee 2013
    # Methods, "Population PK model development" paragraph 2. NONMEM subroutine
    # ADVAN2 TRANS2 was used (CL, V, F1, Ka parameterization).
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Dose in mg, volumes in L -> central/vc has units mg/L = ug/mL.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
