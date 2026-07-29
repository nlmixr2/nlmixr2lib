Feng_2006_enoxaparin <- function() {
  description <- "Two-compartment population PK model for enoxaparin in adult inpatients receiving continuous intravenous infusion (CII) or subcutaneous (SC) dosing (Feng 2006)"
  reference <- "Feng Y, Green B, Duffull SB, Kane-Gill SL, Bobek MB, Bies RR. Development of a dosage strategy in patients receiving enoxaparin by continuous intravenous infusion using modelling and simulation. Br J Clin Pharmacol. 2006;62(2):165-176. doi:10.1111/j.1365-2125.2006.02650.x"
  vignette <- "Feng_2006_enoxaparin"
  units <- list(time = "hour", dosing = "IU", concentration = "IU/L")

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed baseline weight. Used as a linear (no allometric exponent) effect on V2 with reference 70 kg: V2 = theta_V2 * (WT / 70) per Feng 2006 Results equation. Pooled cohort range 16-108 kg (CII + SC studies).",
      source_name        = "WT"
    ),
    CRCL = list(
      description        = "Cockcroft-Gault creatinine clearance (raw, not BSA-normalized)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-fixed per dosing record. Computed by the Cockcroft-Gault equation using ideal body weight for patients with stable serum creatinine, and by the Brater equation when two SCr concentrations measured over 12 h apart differed by > 0.2 mg/dL (Feng 2006 Methods, p. 167). Raw mL/min (NOT BSA-normalized to mL/min/1.73 m^2); stored under the canonical CRCL column per inst/references/covariate-columns.md (CRCL accepts raw mL/min when the source paper does not apply BSA normalization, with the per-model description recording the assay form). The effect on CL is encoded as an additive linear covariate term TVCL = theta_NR + theta_CrCL * (CRCL_L_per_h / 4.8), with 4.8 L/h = 80 mL/min the literature 'normal' renal-function cutoff (Feng 2006 Results equation; cutoff cited from refs 42-43). When CRCL is missing the paper used TVCL = 0.972 L/h (Feng 2006 Table 2 'CLmissing'); this branch is documented in the vignette but the model file encodes only the with-CRCL form, so users must supply CRCL.",
      source_name        = "CRCL"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 83L,
    n_studies      = 2L,
    age_range      = "16-90 years (CII study) / 44-86 years (SC study)",
    age_median     = "CII general medical unit 60.9 y, CII ICU 59.3 y, SC general medical unit 75.1 y; combined 66.6 y",
    weight_range   = "16-108 kg (combined CII + SC)",
    weight_median  = "CII general medical unit 74.1 kg, CII ICU 73.3 kg, SC general medical unit 67.7 kg; combined 71.0 kg",
    sex_female_pct = 51.8,
    race_ethnicity = "Not reported",
    disease_state  = "Adult inpatients receiving therapeutic enoxaparin. CII cohort (n = 48): 29 general medical unit + 19 intensive care unit patients receiving continuous IV infusion (mean infusion 138 +/- 158 h, rate range 100-1600 IU/h, mean 500 +/- 210 IU/h) for treatment of acute thromboembolic disease (Cleveland Clinic Foundation, Jan 1997 - Dec 1998). SC cohort (n = 35): general medical unit patients in the Green 2003 study, included to stabilize PK estimates.",
    dose_range     = "Initial CII regimen 100 IU/kg per 12 h (8.3 IU/kg/h); infusion rates 100-1600 IU/h (mean 500 IU/h). SC cohort dosed per Green 2003 protocol.",
    regions        = "United States (Cleveland Clinic Foundation, OH) and the Green 2003 SC cohort.",
    renal_function = "CrCL median 45.0 mL/min combined; subgroup medians: SC general 39.2 mL/min (range 14.9-95.7), CII general 63.5 mL/min (range 31.1-128.3), CII ICU 26.8 mL/min (range 7.6-49.6). 27 of 48 CII patients had missing SCr concentrations.",
    notes          = "Baseline demographics per Feng 2006 Table 1 (combined CII + SC dataset, n = 83). The pooled analysis used 363 anti-Xa observations from the CII study and 309 from the SC study (672 total). Anti-Xa concentrations were measured by chromogenic assay (LMWH activity). Patients with unstable serum creatinine (defined as two values > 0.2 mg/dL apart over 12 h) had CrCL computed by the Brater equation; otherwise Cockcroft-Gault on ideal body weight. The model also reports a missing-CrCL branch CL = 0.972 L/h (Table 2 CLmissing) for the 27 CII patients without SCr; this branch is not encoded in the structural model file (see vignette Assumptions and deviations)."
  )

  ini({
    # Structural PK parameters at reference weight 70 kg and reference CrCL 80 mL/min
    # (4.8 L/h); Feng 2006 Table 2 Final-model column.
    lka      <- log(0.476);          label("Subcutaneous absorption rate (Ka, 1/h)")                            # Feng 2006 Table 2 Final: Ka = 0.476 1/h (27.3% SE)
    lcl_nr   <- fixed(log(0.229));   label("Non-renal clearance component theta_NR (L/h, FIXED)")               # Feng 2006 Table 2 Final: theta_NR = 0.229 L/h (FIXED from Green 2003 ref [20])
    e_crcl_cl <- 0.744;              label("Renal CL slope per (CRCL_L_h / 4.8) i.e. per CRCL/80 mL/min (L/h)") # Feng 2006 Table 2 Final: theta_CrCL = 0.744 L/h (18.7% SE)
    lvc      <- log(6.78);           label("Central volume V2 at reference WT 70 kg (L per 70 kg)")             # Feng 2006 Table 2 Final: V2 = 6.78 L/70 kg (19.2% SE)
    lvp      <- log(6.19);           label("Peripheral volume V3 (L)")                                          # Feng 2006 Table 2 Final: V3 = 6.19 L (24.9% SE)
    lq       <- log(0.429);          label("Intercompartmental clearance Q (L/h)")                              # Feng 2006 Table 2 Final: Q = 0.429 L/h (24.7% SE)
    lfdepot  <- log(0.94);           label("Subcutaneous bioavailability F1 (fraction)")                        # Feng 2006 Table 2 Final: F1 = 0.94 (9.7% SE)

    # Inter-individual variability on log-normal PK parameters; omega^2 = log(CV^2 + 1).
    # Feng 2006 Table 2 Final CV%.
    etalcl ~ 0.15328 # log(0.407^2 + 1) = 0.15328; Feng 2006 Table 2 Final: omega_CL = 40.7%
    etalvc ~ 0.08290 # log(0.294^2 + 1) = 0.08290; Feng 2006 Table 2 Final: omega_V2 = 29.4%

    # Residual error (combined additive + proportional) from the general medical unit
    # subgroup (Feng 2006 Table 2 Final). The ICU subgroup had a proportional-only
    # residual error with higher CV (sigma3 = 44.0%) which the model file does not encode
    # separately; see vignette Assumptions and deviations for the deviation rationale.
    propSd   <- 0.121; label("Proportional residual SD on anti-Xa Cc (general medical unit, fraction)") # Feng 2006 Table 2 Final: sigma1 = 12.1% (100% SE; small subgroup precision)
    addSd    <- 132;   label("Additive residual SD on anti-Xa Cc (general medical unit, IU/L)")         # Feng 2006 Table 2 Final: sigma2 = 132 IU/L (44.7% SE)
  })

  model({
    # 1. Typical-value clearance has an additive linear covariate term on CRCL with
    # divisive normalization to the literature "normal" reference (4.8 L/h = 80 mL/min,
    # Feng 2006 refs 42-43). theta_NR is the intercept (non-renal CL); theta_CrCL is
    # the slope. The log-normal eta wraps the entire (intercept + slope * CRCL) sum.
    # Convert raw CRCL (mL/min) to L/h via x60/1000 = x0.06 implicitly inside the ratio:
    #   (CRCL_L_h / 4.8) = (0.06 * CRCL_mL_min / 4.8) = (CRCL / 80).
    cl <- (exp(lcl_nr) + e_crcl_cl * (CRCL / 80)) * exp(etalcl)

    # 2. Volumes -- V2 scales linearly with WT/70 (no allometric exponent, per Feng
    # 2006 Results equation V2 = theta_2 * (weight/70) * exp(eta_V2)).
    vc <- exp(lvc + etalvc) * (WT / 70)
    vp <- exp(lvp)
    q  <- exp(lq)

    # 3. Absorption (SC route) -- first-order from depot.
    ka <- exp(lka)

    # 4. Micro-constants.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # 5. ODE system. CII (continuous IV infusion) is administered into central via
    # the user's dosing record (cmt = central, rate > 0). SC doses enter depot.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-               k12 * central - k21 * peripheral1

    # 6. Subcutaneous bioavailability applies to the depot compartment only; IV
    # doses delivered directly to central are assumed F = 1.
    f(depot) <- exp(lfdepot)

    # 7. Anti-Xa concentration: dose in IU, volume in L -> Cc in IU/L.
    # Therapeutic range 0.5-1.2 IU/mL = 500-1200 IU/L (Feng 2006 Methods).
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
