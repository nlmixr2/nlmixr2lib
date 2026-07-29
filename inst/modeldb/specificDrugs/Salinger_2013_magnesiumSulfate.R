Salinger_2013_magnesiumSulfate <- function() {
  description <- "One-compartment population PK model of magnesium sulphate (MgSO4-7H2O) with first-order intramuscular absorption, IV dosing into the central compartment, and an endogenous baseline magnesium term added to the administered drug, in pregnant women with pre-eclampsia (Salinger 2013)."
  reference <- "Salinger DH, Mundle S, Regi A, Bracken H, Winikoff B, Vicini P, Easterling T. Magnesium sulphate for prevention of eclampsia: are intramuscular and intravenous regimens equivalent? A population pharmacokinetic study. BJOG 2013;120:894-900. doi:10.1111/1471-0528.12222"
  vignette  <- "Salinger_2013_magnesiumSulfate"
  units     <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Maternal body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power effect on V; normalised as (WT/55)^theta_1 per Salinger 2013 Table 2 footnote (reference 55 kg, the median maternal weight in the Indian cohort).",
      source_name        = "WT"
    ),
    CREAT = list(
      description        = "Maternal serum creatinine concentration",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Inverse power effect on CL: CL_i = CL * (0.8 / CREAT_i)^theta_2 per Salinger 2013 Table 2 footnote (reference 0.8 mg/dL, equivalent to ~61 umol/L; the paper denotes the column 'CrConc' for serum creatinine concentration in mg/dL).",
      source_name        = "CrConc"
    )
  )

  population <- list(
    n_subjects     = 258L,
    n_studies      = 1L,
    age_range      = "18-41 years (mean ~24-25)",
    age_median     = "~24 years",
    weight_range   = "39-94 kg (mean 55.6-57.3)",
    weight_median  = "~55 kg",
    sex_female_pct = 100,
    race_ethnicity = "Indian (women enrolled at low-resource obstetric hospitals in Nagpur and Vellore, India)",
    disease_state  = "Pre-eclampsia (gestational age 23-41 weeks, mean ~34) at risk for eclampsia",
    dose_range     = "IV arm: 4 g MgSO4-7H2O IV loading over 20 minutes followed by 1 g/h IV maintenance for 12 h (16 g total). IM arm: 4 g MgSO4-7H2O IV loading over 20 minutes, then 10 g IM (5 g per buttock), then 5 g IM every 4 h for 12 h (24 g total).",
    regions        = "India",
    notes          = "Pharmacokinetic study within a randomised trial comparing IV and IM maintenance regimens (NCT00666133). 147 IV + 153 IM enrolled; 126 IV + 132 IM women contributed PK data after exclusions for missing samples and mislabelled times. A single magnesium concentration was drawn per woman; the population analysis pooled 258 sparse data points and did not estimate between-subject variability. Baseline demographics in Salinger 2013 Table 1. Doses are administered as MgSO4-7H2O (heptahydrate); serum concentrations are reported as elemental Mg. The PK parameters in this file are expressed in dose units of mg of elemental Mg and concentration units of mg/L of total magnesium (administered + endogenous baseline); convert administered MgSO4-7H2O grams to mg Mg by multiplying by 24.305/246.47 = 0.0986 (e.g., 4 g MgSO4-7H2O = 394 mg Mg)."
  )

  ini({
    # Structural parameters from Salinger 2013 Table 2 (Final PK Model). Reference
    # subject: 55 kg maternal body weight, 0.8 mg/dL serum creatinine.
    # Dose units used here: mg of elemental Mg (multiply g of MgSO4-7H2O by 0.0986).
    # Concentration units: mg/L of total magnesium (administered + endogenous BL).
    lcl     <- log(4.81);  label("Clearance for the reference subject (CL, L/h)")             # Salinger 2013 Table 2 (48.1 dL/h = 80.1 mL/min)
    lvc     <- log(15.6);  label("Volume of distribution for the reference subject (V, L)")    # Salinger 2013 Table 2 (156 dL)
    lka     <- log(0.317); label("Intramuscular first-order absorption rate (KA, 1/h)")        # Salinger 2013 Table 2
    lfdepot <- log(0.862); label("Intramuscular bioavailability (F, fraction)")                # Salinger 2013 Table 2
    lrbase     <- log(20.8);  label("Endogenous steady-state baseline magnesium (BL, mg/L)")      # Salinger 2013 Table 2 (0.85 mmol/L = 2.08 mg/dL = 20.8 mg/L)

    # Covariate effects from Salinger 2013 Table 2.
    e_wt_vc    <- 0.692; label("Power exponent of body weight on V (theta_1, unitless)")              # Salinger 2013 Table 2
    e_creat_cl <- 1.48;  label("Power exponent of (0.8 / serum creatinine) on CL (theta_2, unitless)") # Salinger 2013 Table 2

    # No inter-individual variability is included. Salinger 2013 Methods state:
    # "because we had only a single data point per woman, we did not attempt to
    # model the between-subject variability." All variability is captured by the
    # proportional residual error term below.

    # Residual error from Salinger 2013 Table 2: proportional, 22.9% CV.
    # The paper reports residual variability on the *total* observed magnesium
    # concentration (administered + BL), consistent with the Methods statement
    # that "the administered magnesium was modelled as additive to BL."
    propSd <- 0.229; label("Proportional residual error (SD, fraction)")  # Salinger 2013 Table 2
  })
  model({
    # Individual PK parameters. Reference subject: WT = 55 kg, CREAT = 0.8 mg/dL.
    # Power covariate forms per the Salinger 2013 Table 2 footnote:
    #   V_i  = V  * (WT_i  / 55)^theta_1
    #   CL_i = CL * (0.8 / CREAT_i)^theta_2  (inverse normalisation: higher SCR -> lower CL)
    cl <- exp(lcl) * (0.8 / CREAT)^e_creat_cl
    vc <- exp(lvc) * (WT  / 55)^e_wt_vc
    ka <- exp(lka)
    rbase <- exp(lrbase)

    # One-compartment model with first-order IM absorption. The depot compartment
    # is dosed for IM administration; IV doses go directly into the central
    # compartment (set cmt = central on IV dose rows in the events file).
    kel <- cl / vc

    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Intramuscular bioavailability applied to depot doses only (IV doses to
    # central are unaffected, matching the Salinger 2013 Table 2 footnote
    # "*Not appropriate for intravenous administration").
    f(depot) <- exp(lfdepot)

    # Observation: administered magnesium concentration (central / vc) plus the
    # endogenous baseline (rbase). Concentration units are mg/L of elemental Mg.
    Cc <- central / vc + rbase

    Cc ~ prop(propSd)
  })
}
