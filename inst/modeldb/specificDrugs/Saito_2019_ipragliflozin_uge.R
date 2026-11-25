Saito_2019_ipragliflozin_uge <- function() {
  description <- paste(
    "Exposure-response (PD-only) Emax model for the change from baseline in",
    "24 h urinary glucose excretion (delta UGE24h) as a function of the",
    "ipragliflozin AUC over the same 24 h interval, in Japanese healthy",
    "adults and adults with type 2 diabetes mellitus (Saito 2019). Emax",
    "is a power function of fasting plasma glucose and absolute GFR while",
    "their product FAC = FPG x GFR is at most 18000, and a constant",
    "107 g/24 h above it. No IIV. Exposure enters as the AUC_IPRA",
    "covariate; the companion popPK model Saito_2019_ipragliflozin gives",
    "the steady-state AUC as 1000 x dose / CL."
  )
  reference <- paste(
    "Saito M, Kaibara A, Kadokura T, Toyoshima J, Yoshida S, Kazuta K,",
    "Ueyama E. Pharmacokinetic and pharmacodynamic modelling for renal",
    "function dependent urinary glucose excretion effect of ipragliflozin,",
    "a selective sodium-glucose cotransporter 2 inhibitor, both in healthy",
    "subjects and patients with type 2 diabetes mellitus.",
    "Br J Clin Pharmacol. 2019;85(8):1808-1819. doi:10.1111/bcp.13972"
  )
  vignette <- "Saito_2019_ipragliflozin"
  units <- list(
    time = "h",
    dosing = "n/a (no dosing events; ipragliflozin exposure enters as the AUC_IPRA covariate)",
    concentration = "g/24 h (change from baseline in 24 h urinary glucose excretion; observation output, not a drug concentration)",
    AUC_IPRA = "ng*h/mL"
  )

  covariateData <- list(
    AUC_IPRA = list(
      description = "Ipragliflozin plasma AUC over the 24 h interval matching the urine collection (AUC24h)",
      units = "ng*h/mL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Saito 2019 fitted the Emax model to noncompartmental AUC0-24h",
        "(Phoenix WinNonlin) after the first and the last dose, paired with",
        "the UGE24h of the same collection interval. For simulation the",
        "paper used the steady-state AUC24h = dose / individual post-hoc",
        "CL/F from the popPK model (Saito_2019_ipragliflozin): in these",
        "units AUC_IPRA = 1000 x dose (mg) / CL (L/h). AUC_IPRA = 0 gives",
        "zero drug effect."
      ),
      source_name = "AUC24h"
    ),
    FPG = list(
      description = "Baseline fasting plasma glucose",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Reference 100 mg/dL (Equation 2). Enters Emax as (FPG / 100)^1.37 and the FAC threshold FPG x CRCL <= 18000.",
      source_name = "FPG"
    ),
    CRCL = list(
      description = "Absolute glomerular filtration rate: Japanese-modified MDRD eGFR de-normalised by the individual Du Bois body surface area",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Baseline GFR in mL/min, NOT BSA-normalised: Saito 2019 Results",
        "'Demographics' -- 'GFR (mL/min) corrected by individual BSA was",
        "used for modelling'. Supply eGFR x BSA / 1.73. Reference 90 mL/min",
        "(Equation 2). Enters Emax as (CRCL / 90)^0.623 and the FAC",
        "threshold FPG x CRCL <= 18000."
      ),
      source_name = "GFR"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 111L,
    n_studies = 3L,
    studies = paste(
      "Study A (CL-0101, NCT01121198; phase I single and multiple dose in",
      "healthy Japanese men), study B (CL-0070, NCT01023945; 14-day",
      "multiple dose in T2DM) and study C (CL-0073, NCT01097681; single",
      "50 mg dose in T2DM with normal, mild or moderately impaired renal",
      "function). 155 delta UGE24h records from 65 healthy subjects and",
      "46 T2DM patients; the 1 and 3 mg doses were excluded."
    ),
    age_range = "20-41 years (healthy) and 34-75 years (T2DM) (Table 1)",
    weight_range = "51.4-80.1 kg (healthy) and 45.6-100.8 kg (T2DM) (Table 1)",
    sex_female_pct = 100 * 16 / 137,
    race_ethnicity = "Japanese",
    disease_state = "Healthy adult men and adults with type 2 diabetes mellitus (eGFR >= 30 mL/min/1.73 m^2)",
    dose_range = "Ipragliflozin 10-300 mg single dose, 20-100 mg once daily (healthy); 50-100 mg once daily or 50 mg single dose (T2DM)",
    regions = "Japan",
    fpg_range = "83-110 mg/dL (healthy) and 84-255 mg/dL (T2DM) (Table 1)",
    renal_function = "GFR 72.2-153.0 mL/min (healthy) and 29.8-169.8 mL/min (T2DM) (Table 1)",
    notes = "Table 1 counts (84 healthy, 53 T2DM incl. placebo) describe the enrolled cohorts; 111 subjects contributed to the PK/PD fit. Female percentage computed over the 137 enrolled subjects."
  )

  ini({
    lec50 <- log(1590)
    label("AUC24h giving half-maximal delta UGE24h, EC50 (ng*h/mL)") # Table 2 'EC50 (ng.h/mL)' = 1590
    lemax <- log(72.3)
    label("Emax at FPG 100 mg/dL and GFR 90 mL/min when FAC <= 18000 (g/24 h)") # Table 2 'Emax (g/24 h): FAC <= 18 000' = 72.3; Equation 2
    lemax_hifac <- log(107)
    label("Constant Emax when FAC = FPG x GFR > 18000 (g/24 h)") # Table 2 'Emax (g/24 h): FAC > 18 000' = 107; Equation 3
    e_fpg_emax <- 1.37
    label("Power exponent of FPG (mg/dL, ref 100) on Emax (unitless)") # Table 2 'FPG effect on Emax' = 1.37; Equation 2
    e_crcl_emax <- 0.623
    label("Power exponent of GFR (mL/min, ref 90) on Emax (unitless)") # Table 2 'GFR effect on Emax' = 0.623; Equation 2

    addSd <- 0.352
    label("Additive residual error (g/24 h)") # Table 2 'Additive error' = 352 (mg; Results '+/-352 mg') -> 0.352 g
    propSd <- 0.198
    label("Proportional residual error (fraction)") # Table 2 'Proportional error' = 0.198
  })
  model({
    # FAC threshold 18000 (mg/dL x mL/min) preset from the exploratory plots
    # (Saito 2019 Results 'Exploratory assessment of PD'), not estimated.
    fac <- FPG * CRCL
    emax_lofac <- exp(lemax) * (FPG / 100)^e_fpg_emax * (CRCL / 90)^e_crcl_emax
    emax <- emax_lofac
    if (fac > 18000) {
      emax <- exp(lemax_hifac)
    }
    ec50 <- exp(lec50)

    # Equation 1
    ugecfb <- emax * AUC_IPRA / (ec50 + AUC_IPRA)
    ugecfb ~ add(addSd) + prop(propSd)
  })
}
