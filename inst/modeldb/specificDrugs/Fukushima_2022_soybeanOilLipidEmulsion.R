Fukushima_2022_soybeanOilLipidEmulsion <- function() {
  description <- "One-compartment plasma triglyceride (TG) kinetic model for a soybean oil-based intravenous lipid emulsion (Intralipos 20%) in adult Japanese inpatients receiving parenteral nutrition: zero-order endogenous TG production (EndVLDL) plus the zero-order emulsion infusion into a plasma volume fixed at 0.422 dL/kg, with first-order elimination whose rate constant falls with baseline TG and body weight (Fukushima 2022)."
  reference <- paste(
    "Fukushima K, Omura K, Goshi S, Okada A, Tanaka M, Tsujimoto T, Iriyama K, Sugioka N.",
    "Individualization of the infusion rate of a soybean oil-based intravenous lipid emulsion",
    "for inpatients, based on baseline triglyceride concentrations: A population",
    "pharmacokinetic approach. JPEN J Parenter Enteral Nutr. 2022;46:104-113.",
    "doi:10.1002/jpen.2111.",
    sep = " "
  )
  vignette <- "Fukushima_2022_soybeanOilLipidEmulsion"
  units <- list(time = "h", dosing = "mg", concentration = "mg/dL")

  covariateData <- list(
    TRIG = list(
      description = "Baseline (pre-infusion) serum triglyceride concentration",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Measured at the ~6:00 AM routine blood draw on the day of administration,",
        "before the infusion; held constant per subject (baseline, not time-varying).",
        "Power effect on Kel, (TRIG / 94)^e_trig_kel, centred on the cohort median",
        "94 mg/dL (Table 4 note). Also sets the initial plasma TG amount,",
        "central(0) = TRIG * vc (see the vignette for how this was established).",
        "Patients with baseline TG > 300 mg/dL were excluded (Methods, Patients)."
      ),
      source_name = "TGbase"
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Baseline body weight. Scales the fixed plasma volume (0.422 dL/kg) and enters",
        "Kel as a power effect (WT / 55.8)^e_wt_kel centred on the cohort median",
        "55.8 kg (Table 4 note)."
      ),
      source_name = "body weight"
    )
  )

  covariatesDataExcluded <- list(
    APOC3 = list(
      description = "Serum apolipoprotein C-III concentration",
      units = "mg/dL",
      type = "continuous",
      notes = "Significant on Kel in the univariate screen (Table 2, dOFV -24.5) but its added value after baseline TG (Table 3 model 2, dOFV -7.6) was far below body weight's (model 4, dOFV -26.3), so it was not carried into the full model; the authors attribute the screen signal to collinearity with baseline TG (Discussion)."
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/dL",
      type = "continuous",
      notes = "Significant on Kel in the univariate screen (Table 2, dOFV -12.8) but not retained after baseline TG (Table 3 model 5, dOFV -6.3)."
    )
  )

  compartmentData <- list(
    central = list(
      analyte = "triglyceride",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 83L,
    n_studies = 1L,
    age_range = "20-95 years",
    age_median = "74 years",
    weight_range = "31.2-86.1 kg",
    weight_median = "55.8 kg",
    sex_female_pct = 47.0,
    race_ethnicity = c(Asian = 100),
    disease_state = paste(
      "Adult inpatients receiving peripheral (n = 75) or central venous (n = 8)",
      "parenteral nutrition, mostly with acute inflammatory gastrointestinal disease",
      "(colon diverticulitis, acute cholecystitis, peptic ulcer, ischemic colitis,",
      "colorectal cancer, ...); low serum albumin (median 3.52 g/dL) and raised CRP",
      "(median 1.56 mg/dL). Baseline TG median 94 mg/dL (range 30-291)."
    ),
    dose_range = paste(
      "Soybean oil-based 20% lipid emulsion (Intralipos 20%), single infusion of",
      "20-120 g fat (median 50 g) at 0.053-0.363 g/kg/h (median 0.217) over",
      "0.75-19.8 h (median 5.2 h)."
    ),
    regions = "Japan (2 hospitals: Ageo Central General Hospital, Joetsu General Hospital)",
    notes = paste(
      "Prospective observational study, October 2016 - March 2018. 238 TG",
      "concentrations: 2 samples during the infusion per patient (first within",
      "0-1, 1-2 or 2-3 h by randomised group, second at the end of the infusion)",
      "plus 72 next-morning routine samples. 44 male / 39 female. Demographics",
      "from Table 1."
    )
  )

  ini({
    # Plasma volume fixed from the literature: V (dL) = 0.422 x body weight
    lvc <- fixed(log(0.422)); label("Plasma volume per kg body weight (dL/kg)") # Table 4 'V(dl) = 0.422 x body weight'; Methods 'fixed with the previously reported mean value of healthy individuals'

    lkin <- log(5709); label("Zero-order endogenous TG production rate EndVLDL (mg/h)") # Table 4 theta1 = 5709 mg/h (RSE 11.4%)
    lkel <- log(2.51); label("First-order TG elimination rate constant at median TRIG and WT (1/h)") # Table 4 theta2 = 2.51 1/h (RSE 10.3%)

    e_trig_kel <- -0.833; label("Power exponent of baseline TG (TRIG/94) on kel (unitless)") # Table 4 theta3 = -0.833 (RSE 14.9%)
    e_wt_kel <- -1.27; label("Power exponent of body weight (WT/55.8) on kel (unitless)") # Table 4 theta4 = -1.27 (RSE 20.0%)

    etalkel ~ 0.187489 # Table 4 omega Kel = 43.3%, read as the SD of eta (0.433^2); see vignette

    propSd <- 0.318; label("Proportional residual error (fraction)") # Table 4 delta = 31.8% (RSE 6.4%), Eq 3 CObs = C x (1 + eps)
  })
  model({
    # Individual parameters (Table 4 final model; no IIV on EndVLDL or V,
    # Table 3 model 9 and Supplementary Model Selection)
    vc <- exp(lvc) * WT
    kin <- exp(lkin)
    kel <- exp(lkel + etalkel) * (TRIG / 94)^e_trig_kel * (WT / 55.8)^e_wt_kel

    # Eq 1: dA/dt = EndVLDL - A x Kel + InfRate. InfRate is the lipid
    # emulsion infusion, given as a zero-order infusion (rate) into central
    # with amt in mg of triglyceride (1 g fat = 1000 mg TG).
    d/dt(central) <- kin - kel * central
    # Initial plasma TG amount = measured baseline TG x plasma volume
    central(0) <- TRIG * vc

    # C = A / V
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
