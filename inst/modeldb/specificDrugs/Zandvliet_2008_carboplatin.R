Zandvliet_2008_carboplatin <- function() {
  description <- "Two-compartment population PK model for free (ultrafilterable) carboplatin in adult cancer patients receiving combination chemotherapy with indisulam (Zandvliet 2008). Clearance is modelled as a renal + non-renal split: a Cockcroft-Gault creatinine-clearance-proportional renal component (theta1 = 0.76) plus a fixed non-renal component (theta2 = 1.5 L/h, fixed at the Calvert 1989 estimate)."
  reference   <- "Zandvliet AS, Schellens JHM, Dittrich C, Wanders J, Beijnen JH, Huitema ADR. Population pharmacokinetic and pharmacodynamic analysis to support treatment optimization of combination chemotherapy with indisulam and carboplatin. Br J Clin Pharmacol. 2008;66(4):485-497. doi:10.1111/j.1365-2125.2008.03230.x"
  vignette    <- "Zandvliet_2008_carboplatin"
  units       <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    AGE = list(
      description        = "Subject age at baseline",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Used as 'age' in the Cockcroft-Gault creatinine-clearance term in Equation 1 of Zandvliet 2008 (CLcr_CG = (140 - age) * weight * [0.85 if female] * 0.074 / SCr_uM). Table 1 reports median 63 years (range 19-81) across the n = 16 oncology cohort. Time-fixed at baseline.",
      source_name        = "age"
    ),
    WT = list(
      description        = "Body weight at baseline",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Used as 'weight' in the Cockcroft-Gault creatinine-clearance term in Equation 1 of Zandvliet 2008. Table 1 reports median 66 kg (range 43-116). The Discussion notes the indisulam patient population simulated in the dosing-regimen evaluation used sex-stratified geometric means (males 77.1 kg +/- 21%, females 63.7 kg +/- 15%). Time-fixed at baseline.",
      source_name        = "weight"
    ),
    SEXF = list(
      description        = "Biological sex indicator (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Equation 1 of Zandvliet 2008 applies a 0.85 multiplier on Cockcroft-Gault CLcr for female subjects, mapped here as (1 - 0.15 * SEXF). The canonical SEXF = 1 for female matches the paper's '0.85 (for female)' rule. Table 1 reports 11 male / 5 female across the n = 16 cohort.",
      source_name        = "sex"
    ),
    CREAT = list(
      description        = "Serum creatinine concentration at baseline",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Used as 'serum creatinine' in the denominator of the Cockcroft-Gault CLcr term in Equation 1 of Zandvliet 2008; the equation's leading 0.074 constant assumes serum creatinine in micromolar units. Table 1 reports median 74 umol/L (range 33-132). The Discussion notes that in this study serum creatinine was measured by an enzymatic method, which under-reads relative to the alkaline-picrate Jaffe method used in the original Cockcroft and Gault 1976 derivation -- one of the reasons the renal-clearance scaling theta1 = 0.76 (< 1) needed to be estimated.",
      source_name        = "serum creatinine"
    )
  )

  population <- list(
    species              = "human",
    n_subjects           = 16L,
    n_studies            = 1L,
    n_observations_carbo = 111L,
    age_range            = "19-81 years",
    age_median           = "63 years",
    weight_range         = "43-116 kg",
    weight_median        = "66 kg",
    bsa_range            = "1.36-2.22 m^2",
    bsa_median           = "1.73 m^2",
    sex_female_pct       = 31.25,
    serum_creatinine_range_uM    = "33-132 (median 74)",
    race_ethnicity       = "Caucasian (16 / 16)",
    disease_state        = "Adults with solid tumours refractory to standard chemotherapy enrolled in a Phase I dose-escalation study of indisulam + carboplatin combination chemotherapy.",
    dose_range           = "Carboplatin doses were calculated with the Calvert formula to target an AUC of 5-6 mg.min/mL; administered as a 30 min IV infusion on day 2 of each 3- or 4-weekly cycle, the day after a 1 h indisulam infusion (350-600 mg/m^2).",
    regions              = "Netherlands (Antoni van Leeuwenhoek Hospital / Slotervaart Hospital, Amsterdam) and Austria (KFJ-Spital, Vienna).",
    notes                = "Carboplatin PK was fit independently of indisulam exposure as the first step of the sequential PK / PK-PD analysis (Methods 'Pharmacokinetic model of carboplatin'). The PD myelosuppression layer of this paper is NOT extracted here; per operator sidecar (frompeople-847 request-001 / response-001 selecting option B), it is deferred because the indisulam PK structural model the combination PD relies on lives in upstream refs [20] / [21] (Zandvliet 2006 J Pharmacokinet Pharmacodyn 33:543-570; Zandvliet 2007 Clin Cancer Res 13:2970-2976) which are not on disk."
  )

  ini({
    # Carboplatin clearance is modelled per Equation 2 of Zandvliet 2008 as
    #   CL_i (L/h) = [CLcr_CG * theta1 + theta2] * exp(eta_CL)
    # where theta1 = 0.76 is the renal-clearance fraction and theta2 = 1.5 L/h
    # is the non-renal clearance, fixed at the Calvert 1989 estimate (paper
    # ref [16]). Encoding follows the Tod_1998_amikacin precedent
    # (non-renal CL on log scale; IIV applied multiplicatively to the sum).
    cl_renal_fraction <- 0.76;       label("Renal-fraction multiplier on CLcr (Zandvliet theta1, unitless)") # Zandvliet 2008 Table 2 row 1 (theta1 = 0.76, RSE 0.05; Equation 2)
    lcl_nonrenal      <- fixed(log(1.5)); label("Log non-renal CL (Zandvliet theta2, L/h; fixed from Calvert 1989)") # Zandvliet 2008 Table 2 footnote ** (theta2 = 1.5 L/h fixed); Calvert 1989 ref [16]

    lvc <- log(15.5); label("Log central volume V1 (L)")               # Zandvliet 2008 Table 2 row 2 (V_central = 15.5 L, RSE 0.19)
    lq  <- log(3.46); label("Log inter-compartmental clearance Q (L/h)") # Zandvliet 2008 Table 2 row 3 (Q = 3.46 L/h, RSE 0.18)
    lvp <- log(9.86); label("Log peripheral volume V2 (L)")             # Zandvliet 2008 Table 2 row 4 (V_peripheral = 9.86 L, RSE 0.11)

    # Inter-individual variability. Table 2 reports IIV magnitudes as CV%;
    # omega^2 = log(CV^2 + 1) for the log-normal etas. Per Equation 2 the
    # CL eta is applied multiplicatively to the sum (renal + non-renal),
    # i.e. on log(CL_total) rather than on log(theta1) or log(theta2).
    etalcl ~ 0.01679 # log(0.13^2 + 1); IIV CL = 13% CV (Zandvliet 2008 Table 2 row 1 IIV)
    etalvc ~ 0.26236 # log(0.54^2 + 1); IIV V1 = 54% CV (Zandvliet 2008 Table 2 row 2 IIV)
    etalq  ~ 0.19121 # log(0.46^2 + 1); IIV Q  = 46% CV (Zandvliet 2008 Table 2 row 3 IIV)
    etalvp ~ 0.09159 # log(0.31^2 + 1); IIV V2 = 31% CV (Zandvliet 2008 Table 2 row 4 IIV)

    # Residual error. NONMEM was run with FOCE-INTER on log-transformed data,
    # ln(OBS) = ln(PRED) + eps, which is equivalent to a proportional error
    # in linear space with SD = sigma_log (8.2% here). The 'Residual error (%)'
    # row of Table 2 reports the percentage directly.
    propSd <- 0.082; label("Proportional residual error (fraction)") # Zandvliet 2008 Table 2 row 5 (residual error = 8.2%, RSE 0.11)
  })

  model({
    # Cockcroft-Gault creatinine clearance per Equation 1 of Zandvliet 2008.
    # CLcr_CG (L/h) = (140 - AGE) * WT * (1 - 0.15 * SEXF) * 0.074 / CREAT
    # where SEXF in {0 male, 1 female} mirrors the '0.85 (for female)' factor.
    # The 0.074 constant is the unit-conversion combination (mg/dL <-> umol/L
    # via 88.4 and mL/min <-> L/h via 0.06) of the original Cockcroft 1976
    # formula, yielding CLcr in L/h directly when CREAT is in umol/L.
    crcl_cg <- (140 - AGE) * WT * (1 - 0.15 * SEXF) * 0.074 / CREAT

    # Total carboplatin CL = renal + non-renal contributions, log-normal IIV
    # applied multiplicatively to the sum per Equation 2 of Zandvliet 2008.
    cl <- (cl_renal_fraction * crcl_cg + exp(lcl_nonrenal)) * exp(etalcl)
    vc <- exp(lvc + etalvc)
    q  <- exp(lq  + etalq)
    vp <- exp(lvp + etalvp)

    # Two-compartment IV micro-constants.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # Dose in mg, volumes in L -> central / vc has units mg/L. The paper
    # reports concentrations as platinum in umol/L (LLQ 0.24 umol/L); convert
    # via the carboplatin / platinum molar mass (371.25 g/mol),
    # 1 mg/L = 2.694 umol/L.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
