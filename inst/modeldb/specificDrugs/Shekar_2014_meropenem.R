Shekar_2014_meropenem <- function() {
  description <- "Two-compartment IV population PK model for meropenem in critically ill adult patients on extracorporeal membrane oxygenation (ECMO) and historical critically ill control patients with sepsis, with a piecewise covariate on clearance that switches between a fixed RRT-cohort CL and a Cockcroft-Gault-CrCL-driven non-RRT CL (Shekar 2014)"
  reference   <- "Shekar K, Fraser JF, Taccone FS, Welch S, Wallis SC, Mullany DV, Lipman J, Roberts JA; ASAP ECMO Study Investigators. The combined effects of extracorporeal membrane oxygenation and renal replacement therapy on meropenem pharmacokinetics: a matched cohort study. Crit Care. 2014;18(6):565. doi:10.1186/s13054-014-0565-2"
  vignette    <- "Shekar_2014_meropenem"
  units       <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    CRCL = list(
      description        = "Cockcroft-Gault creatinine clearance (raw, not BSA-normalized)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source column CrCL. Computed by the Cockcroft-Gault equation in raw mL/min (NOT BSA-normalized to mL/min/1.73 m^2). Stored under the canonical CRCL column per inst/references/covariate-columns.md (CRCL accepts raw mL/min when the source paper does not apply BSA normalization, with the per-model description recording the assay form; precedent: Delattre 2010 amikacin). Cockcroft-Gault CrCL is conventionally not defined for RRT-dependent subjects; Shekar 2014 reports CrCL only for non-RRT subjects (Table 1: 106 mL/min IQR 98-127 in controls, 108 mL/min IQR 65-183 in ECMO) and the model's CL formula switches off the CrCL-driven term when RRT_CRRT_STATUS = 1. Inside model() the column is converted from mL/min to L/h via the factor 60/1000 = 0.06 to match the L/h units of the structural CL parameter; the published clearance coefficient (CL_CRCL = 1.89) is then dimensionless.",
      source_name        = "CrCL"
    ),
    RRT_CRRT_STATUS = list(
      description        = "Subject-level binary indicator for continuous or extended-session renal replacement therapy active during the PK sampling period",
      units              = "(binary)",
      type               = "binary",
      reference_category = 0,
      notes              = "Source column RRT. 1 = subject was on continuous venovenous hemofiltration (CVVHF, control cohort) or extended daily diafiltration (EDD-f, ECMO cohort); 0 = no RRT. The Shekar 2014 cohort mixes CVVH and EDD-f and the published model treats them identically as a single binary covariate. Stored under the canonical RRT_CRRT_STATUS column per inst/references/covariate-columns.md (distinct from HEMODIAL = intermittent-hemodialysis-only and from HEMODIALYSIS = per-time-point session gate). 5/11 ECMO and 5/10 controls were on RRT. Treated as time-fixed at the subject level (all RRT patients were continuously / daily on RRT during sampling; the indicator does not resolve session timing).",
      source_name        = "RRT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 21L,
    n_studies      = 1L,
    age_range      = "16-66 years (Table 1 IQRs)",
    age_median     = "approx 47 years (medians 29-56 across the four subgroups)",
    weight_range   = "60-100 kg (Table 1 IQRs)",
    weight_median  = "70 kg",
    sex_female_pct = 52,
    race_ethnicity = "Not reported (single-centre Australian university-affiliated tertiary referral ICU)",
    disease_state  = "Critically ill adults requiring intravenous meropenem. ECMO arm (n = 11): venovenous (n = 6) or peripheral venoarterial (n = 5) ECMO for severe cardiorespiratory failure; indications included pneumonia / septic shock (n = 7), cardiogenic shock (n = 2), sickle-cell crisis (n = 1), and primary graft dysfunction post lung transplant (n = 1); 5 of 11 ECMO patients were also on EDD-f RRT. Control arm (n = 10): historical critically ill adults with sepsis (Roberts 2009 [reference 15] without renal dysfunction; Bilgrami 2010 [reference 14] on high-volume CVVHF); 5 of 10 controls were on CVVHF.",
    dose_range     = "Meropenem IV bolus q8h. ECMO: 1 g q8h (n = 8, with 1 g loading), 1 g q8h (n = 2, with 1.5 g loading), 1 g q8h (n = 1, with 2 g loading). Controls without renal dysfunction (Roberts 2009): 1.5 g first dose then 1 g q8h. Controls on CVVHF (Bilgrami 2010): 1 g q8h. None of the RRT-dependent patients received post-RRT supplemental doses.",
    regions        = "Australia (single-centre 650-bed university-affiliated tertiary referral hospital, 27-bed mixed ICU, predominantly cardiothoracic cohort; ethics approval HREC/11/QPCH/121 from Prince Charles Hospital, Brisbane QLD)",
    sofa_score     = "Median Day 1 SOFA: 3 (3-4) controls no RRT; 15 (14-16) controls RRT; 9 (7-14) ECMO no RRT; 16 (13-17) ECMO RRT (Table 1)",
    renal_function = "Cockcroft-Gault CrCL median 106 mL/min (98-127) in controls without RRT and 108 mL/min (65-183) in ECMO without RRT (Table 1). RRT-dependent subjects: CrCL not defined / not reported.",
    rrt_modalities = "Controls RRT: continuous venovenous hemofiltration (CVVHF) using Nephral ST500 AN69 hollow-fibre filters (surface area 2.15 m^2), ultrafiltrate rate 66-100 mL/min, target blood flow 250 mL/min, initiated at least 8 hours prior to sampling. ECMO RRT: extended daily diafiltration (EDD-f) using Fresenius 4008s ARrT plus haemodialysis with AV600S filters connected to the post-oxygenator site of the ECMO circuit, blood flow 200-300 mL/min, dialysate flow 200 mL/min, 6-8 hour sessions.",
    notes          = "Baseline demographics per Shekar 2014 Table 1 (medians with IQR). Median time to PK sampling in ECMO patients was 2 days (range 1-7). The 10 historical controls are drawn from Roberts 2009 (n = 5 without renal dysfunction; J Antimicrob Chemother 64:142-150) and Bilgrami 2010 (n = 5 on high-volume CVVHF; Antimicrob Agents Chemother 54:2974-2978). Sample handling: HPLC with UV detection at 304 nm, ertapenem internal standard, LLOQ 1.0 mg/L."
  )

  ini({
    # Structural fixed-effect parameters from Shekar 2014 Table 2 (final covariate model column).
    # The published covariate model on CL is piecewise on RRT_CRRT_STATUS:
    #   TVCL = theta_CL * I(RRT) + theta_CL_CRCL * CrCL_in_Lh * I(not RRT)
    # The paper's printed equation on page 5 (TVCL = theta_1 * (CL_RRT) + theta_1 * (CL_NORRT * CrCL))
    # collapses both thetas onto a single subscript, but the parameter table reports two distinct
    # fixed effects (CL = 5.1 L/h and CL_CRCL = 1.89) and the Table 3 Monte-Carlo trough simulations
    # at CrCL = 20 / 50 / 80 / 120 / 180 mL/min only reproduce when the second theta is the CL_CRCL
    # slope of 1.89, not theta_1 = 5.1. See vignette Assumptions and deviations for the cross-check.
    lcl       <- log(5.1);  label("Typical CL on RRT (L/h)")              # Shekar 2014 Table 2: CL = 5.1 L/h (model mean) -- applies when RRT_CRRT_STATUS = 1
    e_crcl_cl <- 1.89;      label("CrCL slope on CL off RRT (unitless)")  # Shekar 2014 Table 2: CL_CRCL = 1.89 -- slope applied to (CrCL in L/h) when RRT_CRRT_STATUS = 0
    lvc       <- log(18.7); label("Central volume Vc (L)")                # Shekar 2014 Table 2: Vc = 18.7 L
    lvp       <- log(13.2); label("Peripheral volume Vp (L)")             # Shekar 2014 Table 2: Vp = 13.2 L
    lq        <- log(21.0); label("Intercompartmental clearance Q (L/h)") # Shekar 2014 Table 2: Q = 21.0 L/h

    # Between-subject variability from Shekar 2014 Table 2 (Random effects BSV % CV column).
    # omega^2 = log(CV^2 + 1) for log-normal etas. No IIV reported on Q.
    etalcl ~ 0.23560 # log(0.516^2 + 1); Shekar 2014 Table 2: BSV on CL = 51.6 %CV
    etalvc ~ 0.19073 # log(0.458^2 + 1); Shekar 2014 Table 2: BSV on Vc = 45.8 %CV
    etalvp ~ 0.07916 # log(0.287^2 + 1); Shekar 2014 Table 2: BSV on Vp = 28.7 %CV

    # Combined residual error (Shekar 2014 Methods: "combined exponential and additive random error model";
    # Table 2 final-model row).
    propSd <- 0.137; label("Proportional residual error (fraction)") # Shekar 2014 Table 2: RUV = 13.7 %CV
    addSd  <- 2.3;   label("Additive residual error (mg/L)")         # Shekar 2014 Table 2: RUV = 2.3 mg/L SD
  })
  model({
    # Cockcroft-Gault CrCL (raw mL/min, not BSA-normalized) converted to L/h so the
    # published slope CL_CRCL = 1.89 is applied as a dimensionless multiplier.
    crcl_Lh <- CRCL * 0.06

    # Piecewise typical CL: fixed RRT-cohort value when RRT_CRRT_STATUS = 1; CrCL-driven value when 0.
    tvcl <- exp(lcl) * RRT_CRRT_STATUS + e_crcl_cl * crcl_Lh * (1 - RRT_CRRT_STATUS)

    # Individual PK parameters with log-normal IIV.
    cl <- tvcl * exp(etalcl)
    vc <- exp(lvc + etalvc)
    vp <- exp(lvp + etalvp)
    q  <- exp(lq)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # Dose in mg, volumes in L -> central / vc has units mg/L.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
