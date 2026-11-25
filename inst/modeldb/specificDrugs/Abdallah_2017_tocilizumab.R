Abdallah_2017_tocilizumab <- function() {
  description <- "Two-compartment population PK model for intravenous and subcutaneous tocilizumab in adults with rheumatoid arthritis (Abdallah 2017; pooled SUMMACTA and BREVACTA), with first-order SC absorption, SC bioavailability, parallel linear and Michaelis-Menten elimination from the central compartment, and power body-weight effects on clearances and volumes. Typical-value model: the source reports no IIV or residual-error estimates."
  reference <- "Abdallah H, Hsu JC, Lu P, Fettner S, Zhang X, Douglass W, Bao M, Rowell L, Burmester GR, Kivitz A. Pharmacokinetic and Pharmacodynamic Analysis of Subcutaneous Tocilizumab in Patients With Rheumatoid Arthritis From 2 Randomized, Controlled Trials: SUMMACTA and BREVACTA. J Clin Pharmacol. 2017;57(4):459-468. doi:10.1002/jcph.826"
  vignette <- "Abdallah_2017_tocilizumab"
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  compartmentData <- list(
    depot = list(analyte = "tocilizumab", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "tocilizumab", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "tocilizumab", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power effects (WT/70)^0.51 on the clearance terms (linear CL and inter-compartmental Q) and (WT/70)^0.68 on the volume terms (V2 and V3) (Results, 'Effect of Body Weight on Tocilizumab Exposure and PK Parameters'). The paper does not state whether baseline or time-varying weight was used; Table 2 stratifies by baseline weight.",
      source_name = "WT"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Biological sex indicator, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      notes = "Tested, not retained as a strong covariate (Methods; Results state that weight on CL and V was the only strong covariate dependence). About 83% of the safety population was female (Supplementary Table S1)."
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Tested, not retained (Methods). Mean age about 52 years in every arm (Supplementary Table S1)."
    ),
    CRCL = list(
      description = "Creatinine clearance",
      units = "mL/min",
      type = "continuous",
      notes = "Tested, not retained (Methods). The estimating equation and normalisation are not stated."
    ),
    RHEUMATOID_FACTOR = list(
      description = "Rheumatoid factor",
      units = "U/mL",
      type = "continuous",
      notes = "Tested, not retained (Methods)."
    ),
    ADA_POS = list(
      description = "Anti-tocilizumab antibody positive indicator (1 = positive)",
      units = "(binary)",
      type = "binary",
      notes = "Not identified as a covariate influencing tocilizumab PK (Results, Immunogenicity); 0.8-1.6% of tocilizumab-treated patients were ADA-positive by week 24."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 1759L,
    n_observations = 13642L,
    n_studies = 2L,
    age_mean = "52-53 years by arm (SD 11.5-12.5)",
    weight_mean = "70.3-74.6 kg by arm (SD 16.6-19.1)",
    weight_range = "about 40-140 kg (Results: '40 and 140 kg (values close to the extreme values of weight in the study population)')",
    sex_female_pct = 83,
    disease_state = "Adults with moderate-to-severe rheumatoid arthritis of at least 6 months' duration and an inadequate response to DMARDs (up to 20% anti-TNF inadequate responders), on stable background non-biologic DMARDs.",
    dose_range = "Tocilizumab 162 mg SC every week (SUMMACTA), 162 mg SC every 2 weeks (BREVACTA), or 8 mg/kg IV every 4 weeks (SUMMACTA; flat 800 mg for patients >= 100 kg), for 24 weeks.",
    regions = "International (SUMMACTA NCT01194414 and BREVACTA NCT01232569).",
    notes = paste(
      "PK dataset of 1759 patients pooled from the 24-week double-blind periods of SUMMACTA (n = 1262 in the PK/PD population) and BREVACTA (n = 437); 13 642 serum samples (Methods; Results).",
      "Mostly predose samples plus intensive sampling in PK substudies.",
      "Baseline demographics from Supplementary Table S1 (safety populations): SC QW n = 631, IV Q4W n = 631, SC Q2W n = 437, placebo n = 218.",
      "BREVACTA patients who escaped to rescue therapy were included up to the time of escape."
    )
  )

  ini({
    # Structural parameters - Abdallah 2017 Table 1A (final PK model). The table
    # reports CL and Q in mL/day; they are converted to L/day here so that
    # concentrations come out in mg/L = ug/mL with dose in mg and volumes in L.
    # Reference subject: 70 kg.
    lcl <- log(0.216); label("Linear clearance CL (L/day)") # Table 1A 'CL (mL/day)' theta1 = 216
    lvc <- log(4.5); label("Central volume of distribution V2 (L)") # Table 1A 'V2 (L)' theta2 = 4.5
    lq <- log(0.274); label("Inter-compartmental clearance Q (L/day)") # Table 1A 'Q (mL/day)' theta3 = 274
    lvp <- log(2.8); label("Peripheral volume of distribution V3 (L)") # Table 1A 'V3 (L)' theta4 = 2.8

    # Michaelis-Menten elimination. VM is reported as a concentration rate
    # (mg/L/day); the amount rate leaving the central compartment is
    # VM * Vc * C / (KM + C).
    lvmax <- log(1.9); label("Maximum Michaelis-Menten elimination rate VM (mg/L/day)") # Table 1A 'VM (mg/L/day)' theta5 = 1.9
    lkm <- log(0.34); label("Michaelis-Menten constant KM (ug/mL)") # Table 1A 'KM (ug/mL)' theta6 = 0.34

    # Subcutaneous absorption.
    lka <- log(0.23); label("First-order SC absorption rate ka (1/day)") # Table 1A 'ka (1/day)' theta7 = 0.23
    lfdepot <- log(0.795); label("Absolute SC bioavailability FSC (fraction)") # Results text 'estimated to be 79.5%' (Table 1A 'FSC' theta8 = 0.80, rounded)

    # Body-weight power effects - Results, 'Effect of Body Weight on Tocilizumab
    # Exposure and PK Parameters': CL ~ (WT/70)^CLWT and V ~ (WT/70)^VWT, applied
    # to the clearances (CL, Q) and the volumes (V2, V3).
    e_wt_cl <- 0.51; label("Power exponent of WT/70 on CL and Q (unitless)") # Results 'power coefficients estimated at 0.51 (95%CI, 0.47-0.56)'
    e_wt_vc <- 0.68; label("Power exponent of WT/70 on V2 and V3 (unitless)") # Results '0.68 (95%CI, 0.63-0.74)'

    # Residual error: not reported anywhere in the article or its supplement.
    # Fixed to 0 so the model returns typical-value (noise-free) predictions.
    propSd <- fixed(0); label("Proportional residual error (fraction; 0 - not reported by the source)") # not reported
  })

  model({
    logwt <- log(WT / 70)

    cl <- exp(lcl + e_wt_cl * logwt)
    q <- exp(lq + e_wt_cl * logwt)
    vc <- exp(lvc + e_wt_vc * logwt)
    vp <- exp(lvp + e_wt_vc * logwt)
    vmax <- exp(lvmax)
    km <- exp(lkm)
    ka <- exp(lka)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    Cc <- central / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - vmax * vc * Cc / (km + Cc) - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    f(depot) <- exp(lfdepot)

    Cc ~ prop(propSd)
  })
}
