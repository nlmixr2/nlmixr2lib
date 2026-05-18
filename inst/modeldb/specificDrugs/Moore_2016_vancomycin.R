Moore_2016_vancomycin <- function() {
  description <- "Two-compartment IV population PK model for vancomycin in adult patients on extracorporeal membrane oxygenation (ECMO) therapy (Moore 2016). Linear (additive) covariate effects on CL (Cockcroft-Gault creatinine clearance), Vc, and Vp (body weight), each centered on the cohort median (CRCL 84 mL/min; WT 95 kg). Proportional residual error; IIV on CL and Vc only (Q and Vp had no IIV)."
  reference <- "Moore JN, Healy JR, Thoma BN, Peahota MM, Ahamadi M, Schmidt L, Cavarocchi NC, Kraft WK. A population pharmacokinetic model for vancomycin in adult patients receiving extracorporeal membrane oxygenation therapy. CPT Pharmacometrics Syst Pharmacol. 2016;5(9):495-502. doi:10.1002/psp4.12112"
  vignette <- "Moore_2016_vancomycin"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight (baseline)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Moore 2016 Table 1: mean 95 kg, SD 27 (n = 14 adult ECMO patients). Centering value 95 kg confirmed in Figure 4 caption (dose simulation 'for the patient of median weight (95 kg) and creatinine clearance (84 ml/min)'). Used with an additive linear effect on the typical values of Vc and Vp: Vc_typ = 24.2 + 0.00638 * (WT - 95); Vp_typ = 32.3 + 0.0169 * (WT - 95) (Moore 2016 Table 2 V1WT and V2WT rows).",
      source_name        = "WT"
    ),
    CRCL = list(
      description        = "Cockcroft-Gault creatinine clearance (raw, not BSA-normalized)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Moore 2016 Table 1: mean 84 mL/min, SD 37 (n = 14). Centering value 84 mL/min confirmed in Figure 4 caption. Cockcroft-Gault method per Moore 2016 Table 1 footnote ('Renal impairment definition based on creatinine clearance calculated using the Cockcroft-Gault equation'). Stored under the canonical CRCL column with units mL/min (raw Cockcroft-Gault, not BSA-normalized), matching the existing Delattre_2010_amikacin.R and Goti_2018_vancomycin.R precedents in inst/references/covariate-columns.md. Used with an additive linear effect on the typical value of CL: CL_typ = 2.83 + 0.0154 * (CRCL - 84) (Moore 2016 Table 2 CLCRCL row).",
      source_name        = "CRCL"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 14L,
    n_studies        = 1L,
    age_range        = "19-72 years",
    age_mean         = "47 years (SD 16)",
    weight_mean      = "95 kg (SD 27)",
    sex_female_pct   = 21,
    race_ethnicity   = "Not reported",
    disease_state    = "Critically ill adults (>= 18 years) receiving extracorporeal membrane oxygenation (ECMO) therapy at Thomas Jefferson University Hospital. 86% venoarterial (VA) ECMO; 14% venovenous (VV) ECMO. Renal impairment by Cockcroft-Gault: 50% any (28.6% mild [60-89], 14.3% moderate [30-59], 7.15% severe [15-29]); 0% on renal replacement therapy.",
    dose_range       = "Intravenous vancomycin per standard-of-care dosing; concentration sampling at 30, 60, 120, 240, and 360 minutes after the first infusion plus routine trough monitoring.",
    regions          = "United States (single-center, Thomas Jefferson University Hospital, Philadelphia, PA)",
    renal_function   = "Cockcroft-Gault CRCL mean 84 mL/min (SD 37); 50% had some renal impairment; no patient on renal replacement therapy at enrollment",
    ecmo_type        = "86% venoarterial (VA), 14% venovenous (VV)",
    n_concentrations = 65L,
    notes            = "Demographics from Moore 2016 Table 1. ECMO circuit: ROTAFLOW centrifugal pump + CARDIOHELP system (Maquet) with poly-methyl-pentene QUADROX-D oxygenator; circuit primed with ~600 mL normal saline; cannulas Fem-Flex II (arterial) or Femtrak (venous). Assay: Roche Cobas C501 enzyme immunoassay, LLOQ 1.7 ug/mL. Modeling in NONMEM 7.3 with full-covariate-model approach (no stepwise selection). Final OFV minimization successful with condition number 345; ETA shrinkage 6.9% (V1) and 9.6% (CL); epsilon shrinkage 18%. Nonparametric bootstrap (n = 1000) used for 95% CIs."
  )

  ini({
    # Structural parameters at the cohort-median covariate values
    # (WT = 95 kg, CRCL = 84 mL/min). Moore 2016 Table 2 final-model
    # column "Estimate".
    lcl <- log(2.83); label("Typical CL at CRCL = 84 mL/min (L/h)")           # Moore 2016 Table 2: CL = 2.83 L/h, RSE 33.5%
    lvc <- log(24.2); label("Typical central volume Vc at WT = 95 kg (L)")    # Moore 2016 Table 2: V1 = 24.2 L, RSE 14.5%
    lq  <- log(11.2); label("Intercompartmental clearance Q (L/h)")           # Moore 2016 Table 2: Q  = 11.2 L/h, RSE 15%
    lvp <- log(32.3); label("Typical peripheral volume Vp at WT = 95 kg (L)") # Moore 2016 Table 2: V2 = 32.3 L, RSE 11.8%

    # Additive linear covariate effects centered on the cohort median:
    #   CL_typ = exp(lcl) + e_crcl_cl * (CRCL - 84)
    #   Vc_typ = exp(lvc) + e_wt_vc   * (WT   - 95)
    #   Vp_typ = exp(lvp) + e_wt_vp   * (WT   - 95)
    # Q has no covariates (Moore 2016 Table 2; not listed).
    e_crcl_cl <- 0.0154 ; label("Additive slope of CL on (CRCL - 84) (L/h per mL/min)")  # Moore 2016 Table 2: CLCRCL = 0.0154, RSE 21.3%
    e_wt_vc   <- 0.00638; label("Additive slope of Vc on (WT - 95) (L per kg)")          # Moore 2016 Table 2: V1WT   = 0.00638, RSE 98%
    e_wt_vp   <- 0.0169 ; label("Additive slope of Vp on (WT - 95) (L per kg)")          # Moore 2016 Table 2: V2WT   = 0.0169, RSE 14.6%

    # Inter-individual variability (Moore 2016 Table 2 "Intersubject
    # variability" column reports apparent %CV computed as
    # sqrt(exp(OMEGA) - 1) * 100% per Moore 2016 Methods; therefore
    # omega^2 = log(CV^2 + 1). IIV on Q and Vp was tested and removed
    # because removal did not increase the OFV significantly (Moore 2016
    # Results, "Removal of the interindividual variability terms for
    # intercompartmental clearance (Q) and peripheral volume of
    # distribution (V2) did not increase the objective value function
    # significantly").
    etalcl ~ 0.4655 # log(0.77^2 + 1); 77% CV on CL (Moore 2016 Table 2)
    etalvc ~ 0.1094 # log(0.34^2 + 1); 34% CV on Vc (Moore 2016 Table 2)

    # Proportional residual error. Moore 2016 Table 2 reports the
    # proportional residual variance (sigma^2) as 0.0067; the
    # corresponding SD is sqrt(0.0067) = 0.08185 (~ 8.2% CV).
    propSd <- 0.08185; label("Proportional residual error (fraction)") # Moore 2016 Table 2: proportional error (r^2) = 0.0067, RSE 46.9%
  })
  model({
    # Individual PK parameters. The paper uses an additive linear
    # covariate effect on the typical value, then multiplies by the
    # exponential IIV term -- Moore 2016 Methods: "ISV was modeled using
    # exponential functions" with covariates "described using a linear
    # approach centered around the median covariate."
    cl <- (exp(lcl) + e_crcl_cl * (CRCL - 84)) * exp(etalcl)
    vc <- (exp(lvc) + e_wt_vc   * (WT   - 95)) * exp(etalvc)
    vp <-  exp(lvp) + e_wt_vp   * (WT   - 95)
    q  <-  exp(lq)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # Dose in mg, volumes in L -> central/vc has units mg/L.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
