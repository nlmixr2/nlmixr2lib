Hirt_2008_ibuprofen <- function() {
  description <- "One-compartment population PK model with linear elimination for intravenous ibuprofen-lysine (15-min infusion) administered for closure of patent ductus arteriosus in preterm neonates (Hirt 2008). Total-body clearance increases with postnatal age via a power function (CL = 9.49 mL/h x (PNA / 96.3 h)^1.49) anchored at the cohort median PNA of 96.3 h; the apparent volume of distribution is not influenced by postnatal age, gestational age, body weight, Apgar score, or baseline serum sodium / creatinine / albumin / urine output. Exponential inter-individual variability on CL and V; proportional residual error. The PK-PD link reported by the authors (AUC1D > 600 mg L^-1 h or AUC3D > 900 mg L^-1 h associated with >= 91% PDA closure) is illustrated in the validation vignette rather than carried in this model file."
  reference   <- "Hirt D, Van Overmeire B, Treluyer JM, Langhendries JP, Marguglio A, Eisinger MJ, Schepens P, Urien S. An optimized ibuprofen dosing scheme for preterm neonates with patent ductus arteriosus, based on a population pharmacokinetic and pharmacodynamic study. Br J Clin Pharmacol. 2008;65(5):629-636. doi:10.1111/j.1365-2125.2008.03118.x"
  vignette    <- "Hirt_2008_ibuprofen"
  units       <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    PNA = list(
      description        = "Postnatal age (chronological since birth). Time-fixed at the first dose in the source population fit.",
      units              = "months",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Hirt 2008 reports postnatal age in HOURS (cohort median 96.3 h, range 14-262 h) and parameterises the covariate on CL as `CL = 9.49 x (PNA_hours / 96.3)^1.49`. The canonical PNA column is in MONTHS, so the reference is rescaled inside model() via 96.3 hours / 730.5 hours-per-month = 0.13183 months. The dynamic relationship is unchanged because both numerator and denominator carry the same units factor and cancel. Users supply PNA in months in the dataset (per Hirt 2008 Patients section, PNA was time-fixed at the first dose; subsequent doses still use the entry value).",
      source_name        = "PNA"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description        = "Current body weight.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Tested as a covariate on CL and V via generalised additive modelling on the base model but did not meet the inclusion criteria (Hirt 2008 Methods 'Population pharmacokinetic modelling of ibuprofen' and Results 'Population pharmacokinetics'). Cohort median bodyweight 1015 g (range 490-1986 g) at the first dose. PNA was the only retained covariate in the final model.",
      source_name        = "BW"
    ),
    GA = list(
      description        = "Gestational age at birth.",
      units              = "weeks",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Tested on CL and V but not retained in the final model (Hirt 2008 Results 'Population pharmacokinetics'). Cohort median 28 weeks (range 25-34 weeks). GA was significantly higher in the PDA-closure group than the non-closure group in the downstream PD analysis (28.6 vs 27 weeks, p = 0.04) but did not improve the structural PK fit.",
      source_name        = "GA"
    ),
    APGAR1 = list(
      description        = "Apgar score at 1 minute.",
      units              = "score (0-10)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Tested on CL and V but not retained (Hirt 2008 Methods).",
      source_name        = "APGAR1"
    ),
    APGAR5 = list(
      description        = "Apgar score at 5 minutes.",
      units              = "score (0-10)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Tested on CL and V but not retained (Hirt 2008 Methods).",
      source_name        = "APGAR5"
    ),
    APGAR10 = list(
      description        = "Apgar score at 10 minutes.",
      units              = "score (0-10)",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Tested on CL and V but not retained (Hirt 2008 Methods).",
      source_name        = "APGAR10"
    ),
    CREAT = list(
      description        = "Serum creatinine before treatment.",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Tested on CL and V but not retained (Hirt 2008 Methods). Cohort median 1.0 mg/dL (range 0.6-1.3).",
      source_name        = "serum creatinine"
    ),
    ALB = list(
      description        = "Serum albumin before treatment (mass concentration; SI canonical g/L per the 2026-06-19 register standardization audit).",
      units              = "g/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Source paper reports serum albumin in g/dL (cohort median 1.9 g/dL, range 1.2-3.0; = 19 g/L, range 12-30 g/L in SI). Canonical column is now SI g/L per the 2026-06-19 register standardization audit (1 g/dL = 10 g/L). Excluded covariate: tested on CL and V but not retained (Hirt 2008 Methods). Because ALB is not retained in the final model it is not referenced in model(); no inline g/L -> g/dL conversion is required.",
      source_name        = "albumin (g/dL)"
    ),
    NAS = list(
      description        = "Serum sodium before treatment.",
      units              = "mEq/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Tested on CL and V but not retained (Hirt 2008 Methods). Cohort median 141 mEq/L (range 124-155).",
      source_name        = "serum sodium"
    ),
    UR = list(
      description        = "Urine output before treatment.",
      units              = "mL/kg/h",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Tested on CL and V but not retained (Hirt 2008 Methods). Cohort median 3.4 mL/kg/h (range 0.9-7.7).",
      source_name        = "urine output"
    ),
    OCC = list(
      description        = "Occasion / day-of-administration indicator.",
      units              = "(categorical)",
      type               = "categorical",
      reference_category = "day 1 of treatment",
      notes              = "Tested as occasion (DD = 0 on day 1, DD = 1 on days 2-3) to capture any auto-induction or auto-inhibition independent of age and weight; not retained -- clearance and volume of distribution on days 2-3 were not significantly different from day 1 (Hirt 2008 Results 'Population pharmacokinetics' last sentence).",
      source_name        = "occasion (DAY)"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 66L,
    n_studies       = 1L,
    n_observations  = 129L,
    age_range       = "Postnatal age 14-262 h (median 96.3 h, i.e. 1-11 days) at the first dose; gestational age 25-34 weeks (median 28 weeks)",
    age_median      = "Postnatal age 96.3 h",
    weight_range    = "490-1986 g at the first dose",
    weight_median   = "1015 g",
    sex_female_pct  = NA_real_,
    race_ethnicity  = NULL,
    disease_state   = "Preterm neonates with haemodynamically significant patent ductus arteriosus (PDA) and respiratory distress syndrome, ductal patency confirmed by colour Doppler echography. Exclusions: congenital malformations, hydrops fetalis, intraventricular haemorrhage, renal insufficiency, documented infection.",
    co_medication   = "Standard NICU supportive care; no other ductal-closure agent reported.",
    dose_range      = "Three IV infusions of ibuprofen-lysine at 24 h intervals: day 1 median 10 mg/kg (range 5.9-50 mg/kg; one accidental overdose at 50 mg/kg was retained in the modelling dataset), day 2 median 5 mg/kg (range 2-10.7, n = 63), day 3 median 5 mg/kg (range 2-6.4, n = 49). Each dose calculated from birthweight rounded to nearest 100 g and administered as a 15-min infusion via a peripheral vein with a saline flush afterwards.",
    sampling_window = "Arterial-line plasma samples on days 1, 2, and 3; median sampling time 13.2 h after the preceding infusion (range 8 min - 24 h 50 min). 49 infants contributed 3 samples, 14 infants 2, and 3 infants 1; 129 ibuprofen concentrations total.",
    regions         = "Belgium (Antwerp University Hospital; CHC St Vincent Hospital, Liege).",
    notes           = "Demographics from Hirt 2008 Table 1. Sex distribution and race/ethnicity were not tabulated in the publication; cohort phenotype is preterm neonates (gestational age 25-34 weeks) treated within the first 11 days of life. Assay LLOQ 1 mg/L; assay imprecision 5% CV. Modelling software NONMEM V level 1.1, ADVAN1 TRANS2, FOCE estimation. Final-model evaluation by 1000-replicate bootstrap and visual predictive check (Hirt 2008 Figures 2 and 4)."
  )

  ini({
    # Final-model fixed-effect estimates from Hirt 2008 Table 2 ("Final model
    # Original dataset" column). RSE values shown in parentheses in Table 2.
    # Bootstrap means (1000 replicates) are similar to the original-dataset
    # estimates, confirming model stability.
    #
    # Structural parameters -- ibuprofen central one-compartment disposition.
    # The paper reports CL and V as total quantities (mL h^-1 and mL) without
    # body-weight normalisation in the structural model; PNA is the only
    # retained covariate on CL.
    lcl     <- log(0.00949);   label("Total clearance CL at the cohort reference PNA = 96.3 h (L/h)")  # Hirt 2008 Table 2 Final: CL = 9.49 mL/h (RSE 10%); converted to 0.00949 L/h
    lvc     <- log(0.375);     label("Total central volume of distribution V (L)")                     # Hirt 2008 Table 2 Final: V = 375 mL (RSE 10%); converted to 0.375 L

    # Covariate effect: PNA power exponent on CL.
    # Hirt 2008 final-covariate equation (Results 'Population pharmacokinetics'):
    # CL = 9.49 * (PNA / 96.3)^1.49 with PNA in hours.
    e_pna_cl <- 1.49;          label("PNA power exponent on CL (unitless); PNA reference 96.3 h = 0.13183 months")  # Hirt 2008 Table 2 Final: CL, qAGE = 1.49 (RSE 17%)

    # Inter-individual variability (exponential model on CL and V).
    # Hirt 2008 Table 2 Final reports w(CL) = 62% CV and w(V) = 75% CV;
    # converted to log-normal variance via omega^2 = log(1 + CV^2).
    etalcl ~ 0.32506   # log(1 + 0.62^2) = 0.32506; Hirt 2008 Table 2 Final: w(CL) = 62% CV (RSE 28%)
    etalvc ~ 0.44629   # log(1 + 0.75^2) = 0.44629; Hirt 2008 Table 2 Final: w(V) = 75% CV (RSE 44%)

    # Residual error. Hirt 2008 Methods 'Population pharmacokinetic modelling'
    # selected an exponential model for IIV and a proportional model for
    # residual variability; Table 2 Final reports s = 18% CV.
    propSd <- 0.18;             label("Proportional residual error (fraction)")  # Hirt 2008 Table 2 Final: sigma = 18% CV (RSE 27%)
  })

  model({
    # Convert canonical postnatal age (months) to the source-paper reference
    # (96.3 hours = 96.3 / (30.4375 * 24) months = 0.13183 months) so the
    # ratio reproduces Hirt 2008's `CL = 9.49 * (PNA_hours / 96.3)^1.49`
    # term verbatim.
    pna_ref_months <- 96.3 / (30.4375 * 24)
    f_pna_cl       <- (PNA / pna_ref_months) ^ e_pna_cl

    # Individual PK parameters.
    cl <- exp(lcl + etalcl) * f_pna_cl
    vc <- exp(lvc + etalvc)

    kel <- cl / vc

    # One-compartment with intravenous infusion. The dose record targets
    # `central` with `rate` matching the 15-min infusion (rate = dose / 0.25 h);
    # no depot compartment.
    d/dt(central) <- -kel * central

    # Observation: dose in mg, vc in L -> mg/L (the unit used throughout
    # Hirt 2008 for plasma ibuprofen).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
