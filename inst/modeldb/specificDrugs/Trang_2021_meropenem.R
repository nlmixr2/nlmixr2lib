Trang_2021_meropenem <- function() {
  description <- "Two-compartment population PK model for intravenous meropenem (given with vaborbactam) in noninfected adults and adults with complicated urinary tract, bloodstream or other serious infections, with a sigmoidal Hill relationship between renal clearance and MDRD eGFR, a lower nonrenal clearance at eGFR <= 30 mL/min/1.73 m^2, a power effect of age on clearance, fixed allometric weight scaling, and a cumulative urine compartment for the urinary concentrations"
  reference <- paste(
    "Trang M, Griffith DC, Bhavnani SM, Loutit JS, Dudley MN, Ambrose PG,",
    "Rubino CM. 2021. Population pharmacokinetics of meropenem and",
    "vaborbactam based on data from noninfected subjects and infected",
    "patients. Antimicrob Agents Chemother 65:e02606-20.",
    "doi:10.1128/AAC.02606-20. Covariate-equation forms and reference",
    "values from the FDA Clinical Pharmacology review of NDA 209776",
    "(Vabomere, 2017), section 4.2, Equations (1)-(3).",
    sep = " "
  )
  vignette <- "Trang_2021_meropenem_vaborbactam"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "meropenem", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "meropenem", units = "mg", specimen = "plasma", verified = TRUE),
    urine = list(analyte = "meropenem", units = "mg", specimen = "urine", verified = TRUE)
  )

  covariateData <- list(
    CRCL = list(
      description = paste0(
        "Estimated glomerular filtration rate by the Modification of Diet in ",
        "Renal Disease equation, BSA-normalised (mL/min/1.73 m^2); time-varying"
      ),
      units = "mL/min/1.73 m^2",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Materials and Methods, Demographics: eGFR was calculated from serum ",
        "creatinine, age and sex by the MDRD equation at every serum ",
        "creatinine measurement and treated as time-varying, with serum ",
        "creatinine capped at a lower bound of 0.5 mg/dL. Drives the renal ",
        "clearance arm through the sigmoidal Hill term (no normalising ",
        "value), and also defines the renal group: eGFR <= 30 mL/min/1.73 m^2 ",
        "takes the reduced nonrenal clearance (Results, Covariate analyses; ",
        "RGRP indicator of the FDA review Equation 1). Pooled median 90.1, ",
        "range 4.50-338 mL/min/1.73 m^2 (Table 1)."
      ),
      source_name = "eGFR"
    ),
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Enters CL and CLd as (WT / 80)^0.75 and Vc and Vp as (WT / 80)^1, ",
        "all exponents fixed (Table 2). The final paper does not print the ",
        "normalising weight; 80 kg is the reference printed in the FDA review ",
        "Equations (2)-(3) for the same analysts' initial model of this ",
        "dataset, which the final model updates (see the vignette). Pooled ",
        "median 75.0, range 40.0-177 kg (Table 1)."
      ),
      source_name = "WTKG"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Enters total CL as (AGE / 58)^-0.526 (Table 2). The 58-year reference ",
        "is printed in the FDA review Equation 1 for the initial model and is ",
        "the age the Discussion uses for its typical-patient clearance ",
        "('about 10 liters/h in a 58-year-old patient with an eGFR of 100'). ",
        "Pooled median 53.0, range 18.0-92.0 years (Table 1)."
      ),
      source_name = "AGE"
    ),
    URINE_VOL_INTERVAL = list(
      description = "Urine volume collected in the urine collection interval containing the current urinary observation",
      units = "mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste0(
        "Observation-equation divisor only (Curine = urine amount / collected ",
        "volume). Urine was collected in the phase 1 studies over 0-4, 4-8, ",
        "8-12, 12-24, 24-48 and 48-72 h intervals (supplement Table S3), so ",
        "the urine state must be reset to zero at each interval boundary ",
        "(evid = 5, amt = 0, cmt = 'urine'). Not needed for plasma-only ",
        "simulation; any positive value may then be supplied."
      ),
      source_name = "urine volume"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 413L,
    n_studies = 4L,
    age_range = "18-92 years",
    age_median = "53 years",
    weight_range = "40.0-177 kg",
    weight_median = "75.0 kg",
    sex_female_pct = 55.5,
    race_ethnicity = "Not tabulated; race was screened and not retained (Figure S7).",
    disease_state = paste(
      "Pooled phase 1 noninfected adults (study 501, healthy volunteers; study",
      "504, normal renal function to end-stage renal disease) and phase 3",
      "infected patients (study 505 TANGO I, complicated urinary tract",
      "infection or acute pyelonephritis; study 506 TANGO II, infections due",
      "to confirmed or suspected carbapenem-resistant Enterobacterales,",
      "including bacteremia, HABP/VABP, cIAI and cUTI)."
    ),
    dose_range = paste(
      "Meropenem 1 or 2 g (0.5 g in severe renal impairment in study 506) as",
      "3-h intravenous infusions (1 h in study 501 group 6), single dose or",
      "every 8 h, co-administered with vaborbactam; renal dose adjustments per",
      "supplement Table S3."
    ),
    regions = "Multinational (phase 3 TANGO I and TANGO II).",
    renal_function = "eGFR (MDRD) 4.50-338 mL/min/1.73 m^2; median 90.1.",
    notes = paste(
      "Final data set: 4,264 plasma concentrations from 91 noninfected subjects",
      "and 322 infected patients plus 834 urine concentrations from 84",
      "noninfected subjects (Results). Table 1 demographics are for all 431",
      "enrolled subjects/patients (239 female, 55%). Assay LC-MS/MS, 0.02-100",
      "mg/L; BLQ handled with Beal M3. NONMEM 7.2, FOCE-I."
    )
  )

  ini({
    # Structural parameters -- Trang 2021 Table 2 (final meropenem model).
    # CL = (CLNR * (1 + shift * RGRP) + CLR,max * eGFR^h / (eGFR50^h + eGFR^h))
    #      * (WT/80)^0.75 * (AGE/58)^theta_age
    # (functional form from the FDA review Equation 1; see covariateData).
    lcl_nonren <- log(3.85)
    label("Nonrenal clearance at eGFR > 30 mL/min/1.73 m^2, 80 kg, 58 years (L/h)") # Table 2, 'CL NR' = 3.85 L/h (%SEM 0.70)
    lcl_renal_max <- log(6.58)
    label("Maximal renal clearance, 80 kg, 58 years (L/h)") # Table 2, 'CL R,max' = 6.58 L/h (%SEM 3.20)
    lcrcl50 <- log(40.0)
    label("eGFR at half-maximal renal clearance (mL/min/1.73 m^2)") # Table 2, 'eGFR50' = 40.0 (%SEM 0.10)
    lhill <- log(1.95)
    label("Hill coefficient of the renal clearance-eGFR relationship (unitless)") # Table 2, 'Hill coefficient' = 1.95 (%SEM 13.9)
    lvc <- log(17.0)
    label("Central volume of distribution, 80 kg (L)") # Table 2, 'V c' = 17.0 L (%SEM 22.9)
    lq <- log(1.36)
    label("Distributional clearance, 80 kg (L/h)") # Table 2, 'CL d' = 1.36 L/h (%SEM 0.10)
    lvp <- log(2.32)
    label("Peripheral volume of distribution, 80 kg (L)") # Table 2, 'V p' = 2.32 L (%SEM 0.30)

    # Allometric weight exponents, fixed (Table 2 %SEM column 'Fixed').
    e_wt_cl <- fixed(0.75)
    label("Allometric exponent of weight on CL (unitless)") # Table 2, 'Power coefficient of WTKG on CL' = 0.75 Fixed
    e_wt_vc <- fixed(1.00)
    label("Allometric exponent of weight on Vc (unitless)") # Table 2, 'Power coefficient of WTKG on V c' = 1.00 Fixed
    e_wt_q <- fixed(0.75)
    label("Allometric exponent of weight on CLd (unitless)") # Table 2, 'Power coefficient of WTKG on CL d' = 0.75 Fixed
    e_wt_vp <- fixed(1.00)
    label("Allometric exponent of weight on Vp (unitless)") # Table 2, 'Power coefficient of WTKG on V p' = 1.00 Fixed

    # Covariate effects.
    e_age_cl <- -0.526
    label("Power exponent of age on CL (unitless)") # Table 2, 'Power coefficient of age on CL' = -0.526 (%SEM 16.7)
    e_rgrp_cl_nonren <- -0.650
    label("Proportional shift in nonrenal CL for eGFR <= 30 mL/min/1.73 m^2 (fraction)") # Table 2, 'Proportional shift with renal group on CL NR' = -0.650 (%SEM 10.1)

    # IIV. Table 2 reports %CV; omega^2 = log(CV^2 + 1). The final model
    # fitted a full covariance matrix whose off-diagonal elements are not
    # published, so the etas are entered as independent (vignette errata).
    etalcl ~ 0.18067 # Table 2, IIV 44.5 %CV (printed on the 'CL R,max' row; applied to total CL, see vignette)
    etalvc ~ 0.21047 # Table 2, 'V c' IIV 48.4 %CV
    etalq ~ 0.23688 # Table 2, 'CL d' IIV 51.7 %CV
    etalvp ~ 0.13289 # Table 2, 'V p' IIV 37.7 %CV

    # Residual error. Table 2 prints unlabelled sigma values; they are
    # variances (sigma^2): read as SDs the plasma proportional error would be
    # 4.2%, inconsistent with the +/-2 IWRES spread at 25-75 mg/L and the
    # observed-versus-IPRED scatter of supplement Figure S1 (vignette).
    propSd <- sqrt(0.0423)
    label("Plasma proportional residual SD (fraction)") # Table 2, 'Plasma proportional error' = 0.0423 (variance); SD 0.2057
    addSd <- sqrt(0.0204)
    label("Plasma additive residual SD (mg/L)") # Table 2, 'Plasma additive error' = 0.0204 (variance); SD 0.1428 mg/L
    propSd_Curine <- sqrt(0.207)
    label("Urine proportional residual SD (fraction)") # Table 2, 'Urine proportional error' = 0.207 (variance); SD 0.4550
    addSd_Curine <- sqrt(0.0511)
    label("Urine additive residual SD (mg/L)") # Table 2, 'Urine additive error' = 0.0511 (variance); SD 0.2261 mg/L
  })

  model({
    # 1. Derived covariate terms. Renal group (FDA review Equation 1,
    #    RGRP): 1 when eGFR <= 30 mL/min/1.73 m^2, evaluated on the
    #    time-varying eGFR.
    rgrp <- (CRCL <= 30)
    hill <- exp(lhill)
    crcl50 <- exp(lcrcl50)

    # Nonrenal and renal clearance arms at the reference weight and age.
    cl_nonren <- exp(lcl_nonren) * (1 + e_rgrp_cl_nonren * rgrp)
    cl_renal <- exp(lcl_renal_max) * CRCL^hill / (crcl50^hill + CRCL^hill)
    frac_renal <- cl_renal / (cl_nonren + cl_renal)

    # 2. Individual PK parameters. Weight and age scale total CL; the
    #    renal fraction is carried unchanged into the urine arm.
    cl <- (cl_nonren + cl_renal) * (WT / 80)^e_wt_cl * (AGE / 58)^e_age_cl *
      exp(etalcl)
    vc <- exp(lvc + etalvc) * (WT / 80)^e_wt_vc
    q <- exp(lq + etalq) * (WT / 80)^e_wt_q
    vp <- exp(lvp + etalvp) * (WT / 80)^e_wt_vp

    # 3. Micro-constants.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 4. ODE system. Zero-order intravenous infusion into central via the
    #    event-table rate; renally cleared drug accumulates in urine.
    d/dt(central) <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(urine) <- frac_renal * kel * central

    # 5. Observations. Urine concentration = amount excreted in the
    #    collection interval / collected volume (mL converted to L; guarded
    #    against a zero volume).
    Cc <- central / vc
    urine_volume <- max(URINE_VOL_INTERVAL / 1000, 0.001)
    Curine <- urine / urine_volume

    Cc ~ add(addSd) + prop(propSd)
    Curine ~ add(addSd_Curine) + prop(propSd_Curine)
  })
}
