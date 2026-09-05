Chen_2024_biapenem <- function() {
  description <- "Two-compartment IV population PK model for biapenem in adult patients with sepsis, with a linear centered creatinine-clearance effect on clearance and a linear centered blood-urea-nitrogen effect on intercompartmental clearance (Chen 2024)"
  reference <- "Chen D, Wu X, Zhang H, Yao H, Jin L, Luo X, Liu J, Wu Z, Li Y, Xu W, Ge W, Chen X, Zhu H. Population pharmacokinetics, dosing optimization and clinical outcomes of biapenem in patients with sepsis. Front Pharmacol. 2024;15:1388150. doi:10.3389/fphar.2024.1388150"
  vignette <- "Chen_2024_biapenem"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central     = list(analyte = "biapenem", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "biapenem", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description        = "Cockcroft-Gault creatinine clearance (raw, not BSA-normalized)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source column CLCr. Computed by the Cockcroft-Gault equation (Chen 2024 Methods 2.1) in raw mL/min, NOT BSA-normalized to mL/min/1.73 m^2; the paper reports a separate MDRD eGFR column in mL/min/1.73 m^2 that was screened but did not enter the final model. Stored under the canonical CRCL column per inst/references/covariate-columns.md, which accepts raw mL/min when the source paper does not apply BSA normalization. Enters CL as a linear centered term TVCL = 8.33 * (1 + 0.0046 * (CRCL - 78.2)) with reference 78.2 mL/min (Chen 2024 Table 2; the Discussion confirms 'the mean CL of BPM in patients with sepsis was 8.33 L/h for CLCr of 78.2 mL/min'). Modeling-cohort CLCr median 84.92 mL/min (range 3.5-295.5). The covariate factor stays positive across all attainable CLCr: it would reach zero only at CLCr = -139 mL/min, so no guard is needed.",
      source_name        = "CLCr"
    ),
    BUN = list(
      description        = "Blood urea nitrogen",
      units              = "mmol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Enters intercompartmental clearance as a linear centered term TVQ = 3.75 * (1 + 0.112 * (BUN - 6.8)) with reference 6.8 mmol/L (Chen 2024 Table 2). Reported in SI units (mmol/L), not mg/dL; 1 mmol/L urea ~= 2.80 mg/dL BUN. Modeling-cohort BUN median 6.2 mmol/L (range 0.4-66.9). This is the first registered model to place BUN on Q rather than on a clearance or absorption parameter; the authors interpret it (Discussion) as a marker of sepsis-associated catabolism and neurohormonal activation rather than of filtration per se. The covariate factor stays positive across all attainable BUN: it would reach zero only at BUN = -2.1 mmol/L, so no guard is needed.",
      source_name        = "BUN"
    )
  )

  covariatesDataExcluded <- list(
    ALB = list(
      description = "Serum albumin",
      units       = "g/L",
      type        = "continuous",
      notes       = "Reached significance on CL during forward inclusion but was removed during backward elimination and is NOT in the final model (Chen 2024 Results 3.2). No point estimate for its effect is reported anywhere in the paper, so it cannot be encoded. Retained here to preserve the provenance of the covariate screen; hypoalbuminemia is separately reported as a risk factor for clinical failure (Table 5, OR 0.33)."
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 245L,
    n_studies        = 1L,
    age_range        = "18-97 years",
    age_median       = "63 years",
    weight_range     = "36.8-100 kg",
    weight_median    = "62 kg",
    sex_female_pct   = 36.3,
    race_ethnicity   = "Not reported (single-center Chinese cohort, Nanjing Drum Tower Hospital)",
    disease_state    = "Adults with sepsis as defined by the Third International Consensus Definitions (Sepsis-3). Septic shock in 26.9%, immunocompromised in 13.9%. Primary infection site respiratory 44.1%, intra-abdominal 48.6%, other 7.3%. Patients receiving renal replacement therapy or extracorporeal membrane oxygenation during biapenem therapy were excluded, as were patients treated for less than 48 h and those without therapeutic drug monitoring.",
    dose_range       = "Biapenem 300-600 mg per administration, 2-4 times daily (product-label maximum 1.2 g/day), given as a 1-hour intravenous infusion in 100 mL sodium chloride or glucose.",
    regions          = "China (single center: Nanjing Drum Tower Hospital, Nanjing). Admissions January 2018 to May 2022.",
    renal_function   = "Cockcroft-Gault creatinine clearance median 84.92 mL/min (range 3.5-295.5), raw mL/min and not BSA-normalized. MDRD eGFR median 108.1 mL/min/1.73 m^2 (range 2.6-332.9). Serum creatinine median 66 umol/L (range 28-1655). Renal replacement therapy was an exclusion criterion.",
    n_concentrations = 351L,
    notes            = "Retrospective single-center therapeutic-drug-monitoring study. 466 biapenem measurements from 317 adults were split 3:1 into a modeling cohort (351 samples from 245 patients, the data this model was fit to) and an external evaluation cohort (115 samples from 72 patients). Sparse sampling: on average 1.5 samples per participant, usually drawn after the third dosing interval; median time after the last dose 6 h (range 0.5-20.73). Median trough concentration 2.1 mg/L (modeling) and 1.8 mg/L (external). Assay HPLC-UV at 300 nm, linear 0.3-30 mg/L with a 0.3 mg/L lower limit of quantitation. Estimation in NONMEM 7.3 with the first-order (FO) method; final-model OFV 550.508 versus 683.499 for the base model. Validated by 1000-replicate nonparametric bootstrap (Table 2) and prediction-corrected VPC (Figure 2); external evaluation gave MDPE 6.75%, MAPE 26.25%, F20 42.61%, F30 53.04%. Baseline demographics per Chen 2024 Table 1 (modeling cohort column)."
  )

  ini({
    # Structural parameters at the reference patient (CRCL 78.2 mL/min, BUN 6.8 mmol/L).
    # Chen 2024 Table 2, "Final model / Estimate" column.
    lcl <- log(8.33); label("Clearance (L/h)")                          # Chen 2024 Table 2: theta1 = 8.33 L/h (RSE 6.4%; bootstrap median 8.331, 5th-95th 7.48-9.92)
    lvc <- log(13.4); label("Central volume of distribution V1 (L)")    # Chen 2024 Table 2: theta2 = 13.4 L (RSE 15%; bootstrap median 13.23, 5th-95th 10.94-23.55)
    lq  <- log(3.75); label("Intercompartmental clearance Q (L/h)")     # Chen 2024 Table 2: theta3 = 3.75 L/h (RSE 13.6%; bootstrap median 3.58, 5th-95th 2.26-5.73)
    lvp <- log(60.4); label("Peripheral volume of distribution V2 (L)") # Chen 2024 Table 2: theta4 = 60.4 L (RSE 17.5%; bootstrap median 64.72, 5th-95th 23.38-169.07)

    # Covariate effects, both linear on a subtractively centered covariate:
    #   CL = theta1 * (1 + theta5 * (CLCr - 78.2))
    #   Q  = theta3 * (1 + theta6 * (BUN  - 6.8))
    # NOTE: the running text of Results 3.2 prints theta5 as "0.046", a factor of
    # ten larger than the 0.0046 given in Table 2. Table 2 is correct and the
    # in-text equation has lost a zero -- see the vignette Errata section for the
    # falsification (0.046 drives CL negative for any CLCr below 56.5 mL/min,
    # which covers a large part of this cohort, and the bootstrap 5th-95th
    # interval of 0.0034-0.0064 brackets 0.0046 and excludes 0.046).
    e_crcl_cl <- 0.0046; label("Linear coefficient for centered CRCL on CL (per mL/min)")  # Chen 2024 Table 2: theta5 = 0.0046 (RSE 11.9%; bootstrap median 0.0049, 5th-95th 0.0034-0.0064)
    e_bun_q   <- 0.112;  label("Linear coefficient for centered BUN on Q (per mmol/L)")    # Chen 2024 Table 2: theta6 = 0.112 (RSE 1.9% as printed; bootstrap median 0.114, 5th-95th 0.034-0.15)

    # Inter-individual variability, exponential (Chen 2024 Eq. 1). Only CL and Q
    # carried IIV in the final model; V1 and V2 did not. The Table 2 values are
    # the NONMEM $OMEGA diagonal elements, i.e. VARIANCES on the log scale, so
    # they are used here unchanged -- see the vignette Errata for the
    # adjudication against the paper's own Monte Carlo target-attainment results.
    etalcl ~ 0.0591 # Chen 2024 Table 2: omega_CL = 0.0591 (RSE 19.8%; bootstrap median 0.057, 5th-95th 0.031-0.088); 24.7% CV
    etalq  ~ 1.12   # Chen 2024 Table 2: omega_Q  = 1.12   (RSE 25.2%; bootstrap median 1.10, 5th-95th 0.26-2.94); 144% CV

    # Residual error: additive (Chen 2024 Eq. 2; Results 3.2 "the residual
    # unexplained variability was best described by an additive residual error
    # model"). Table 2 reports the NONMEM $SIGMA element 0.591, a variance, so
    # the standard deviation is sqrt(0.591) = 0.769 mg/L.
    addSd <- sqrt(0.591); label("Additive residual error (mg/L)") # Chen 2024 Table 2: sigma = 0.591 (RSE 17.1%; bootstrap median 0.59, 5th-95th 0.39-0.85); sqrt keeps the published variance literal, = 0.769 mg/L
  })
  model({
    # Individual PK parameters. Both covariate effects are linear multiplicative
    # factors on a subtractively centered covariate, with exponential IIV on the
    # typical value (Chen 2024 Table 2 parameter equations).
    cl <- exp(lcl + etalcl) * (1 + e_crcl_cl * (CRCL - 78.2))
    vc <- exp(lvc)
    q  <- exp(lq + etalq) * (1 + e_bun_q * (BUN - 6.8))
    vp <- exp(lvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    Cc <- central / vc
    Cc ~ add(addSd)
  })
}
