Frymoyer_2013_gentamicin <- function() {
  description <- "One-compartment IV population PK model for gentamicin in 29 term neonates with hypoxic ischemic encephalopathy (HIE) receiving therapeutic hypothermia (Frymoyer 2013), with fixed allometric birth-weight scaling (exponent 0.75 on CL, 1 on Vc, reference 3.3 kg) and a power effect of serum creatinine on PNA day 1 on CL (reference 0.9 mg/dL)."
  reference <- paste(
    "Frymoyer A, Meng L, Bonifacio SL, Verotta D, Guglielmo BJ.",
    "Gentamicin pharmacokinetics and dosing in neonates with hypoxic ischemic",
    "encephalopathy receiving hypothermia.",
    "Pharmacotherapy. 2013;33(7):718-726.",
    "doi:10.1002/phar.1263.",
    sep = " "
  )
  vignette <- "Frymoyer_2013_gentamicin"
  units    <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT_BIRTH = list(
      description        = "Birth weight (time-fixed per subject).",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometric scaling on CL (fixed exponent 0.75) and on Vc (fixed",
        "exponent 1). Reference birth weight 3.3 kg corresponds to the typical",
        "neonate as discussed in the Abstract and Discussion (cohort median 3.32",
        "kg per Table 1, IQR 2.97-3.50 kg). The covariate is birth weight",
        "(Table 1 row 'Birthweight, kg'); Methods 'Population Pharmacokinetic",
        "Analysis' states 'the effect of birthweight on clearance and volume of",
        "distribution was implemented using an allometric model with the",
        "exponent defining the relationship fixed to 0.75 and 1, respectively.'"
      ),
      source_name        = "BW"
    ),
    CREAT = list(
      description        = "Serum creatinine on postnatal age (PNA) day 1 (second day of life).",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-fixed per subject (a single SCr value drawn on PNA day 1).",
        "Reference 0.9 mg/dL is the cohort median (Table 1 IQR 0.8-1.2 mg/dL).",
        "Enters CL as a power effect (CREAT/0.9)^e_creat_cl with",
        "e_creat_cl = -0.566 so that higher SCr lowers CL. Cross-check from",
        "Results: 'a typical neonate with SCr 1.4 mg/dL had a clearance 27%",
        "lower than a neonate with SCr 0.8 mg/dL' = (1.4/0.8)^-0.566 = 0.729,",
        "i.e. 27.1% lower (Results 'Population Pharmacokinetic Analysis')."
      ),
      source_name        = "SCr"
    )
  )

  covariatesDataExcluded <- list(
    GA = list(
      description = "Gestational age at birth (weeks). Screened during covariate model building.",
      units       = "weeks",
      type        = "continuous",
      notes       = paste(
        "Screened; not retained. Cohort median 40.0 weeks (IQR 37.6-40.7;",
        "Table 1). Methods note inclusion required GA >= 36 weeks; Discussion:",
        "'Gestational age did not affect clearance in our study and this was",
        "likely due to the narrow range of gestational ages intrinsic to our",
        "institution's eligibility criteria for hypothermia (>=36 weeks).'"
      )
    ),
    PNA = list(
      description = "Postnatal age (days). Screened during covariate model building.",
      units       = "days",
      type        = "continuous",
      notes       = paste(
        "Screened; not retained. Discussion: 'Since all neonates were of the",
        "same postnatal age +/- 1 day, we were unable to examine this effect.'",
        "TDM was performed at the third or fourth dose, typically PNA days 3-5."
      )
    ),
    PHA = list(
      description = "First arterial or capillary pH. Screened during covariate model building.",
      units       = "pH units",
      type        = "continuous",
      notes       = paste(
        "Screened as a marker of asphyxia severity; not retained. Cohort",
        "median 7.0, IQR 6.9-7.1 (Table 1). Methods 'Population Pharmacokinetic",
        "Analysis' lists 'first blood pH' among the continuous covariates",
        "tested on CL; Results: 'No other covariates examined significantly",
        "affected gentamicin clearance.'"
      )
    ),
    APGAR10 = list(
      description = "10-minute APGAR <= 5 indicator. Screened during covariate model building.",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened as a categorical marker of asphyxia severity; not retained.",
        "Cohort median 10-minute APGAR 5, IQR 3-6 (Table 1). Methods lists",
        "'APGAR <= 5 at 10 minutes' among the categorical covariates tested",
        "on CL; Results: 'No other covariates examined significantly affected",
        "gentamicin clearance.'"
      )
    ),
    CONMED_DOPAMINE = list(
      description = "Concomitant dopamine indicator. Screened during covariate model building.",
      units       = "(binary)",
      type        = "binary",
      notes       = paste(
        "Screened as a marker of cardiovascular instability that may affect",
        "renal perfusion; not retained. 62% of cohort received dopamine",
        "(Table 1, n = 18/29). Methods lists 'concomitant dopamine' among",
        "the categorical covariates tested on CL; Results: not significant."
      )
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 29L,
    n_studies       = 1L,
    n_observations  = 47L,
    age_range       = "GA 37.6-40.7 weeks (median 40.0) at birth; PNA approximately 2-5 days at TDM",
    age_median      = "GA 40.0 weeks; PNA day 1 (SCr) and PNA days 2-5 (gentamicin observations)",
    weight_range    = "2.97-3.50 kg (IQR; median 3.32 kg)",
    weight_median   = "3.32 kg",
    sex_female_pct  = 62,
    race_ethnicity  = "Not reported (single-centre US tertiary NICU).",
    disease_state   = paste(
      "Term neonates (>= 36 weeks GA) with presumed hypoxic ischemic",
      "encephalopathy (HIE) treated with whole-body therapeutic hypothermia",
      "at 33.5 degrees C in a Level III NICU. Asphyxia severity: median",
      "first arterial/capillary pH 7.0 (IQR 6.9-7.1), base deficit -19 mmol/L",
      "(IQR -15 to -27), 10-minute APGAR median 5 (IQR 3-6). 83% required",
      "assisted ventilation; 55% had seizures; 62% received dopamine.",
      "Mortality before discharge: 21%. All neonates were undergoing",
      "hypothermia for the entire gentamicin course; the paper could not",
      "evaluate hypothermia as an independent covariate because there were",
      "no normothermic controls (Discussion: 'We were unable to evaluate",
      "the independent impact of hypothermia on gentamicin pharmacokinetics",
      "since all neonates in our study were cooled.')."
    ),
    dose_range      = paste(
      "Empiric gentamicin: 5 mg/kg every 24 h (25 neonates), 4.5 mg/kg",
      "every 24 h (1 neonate), 4 mg/kg every 24 h (3 neonates). Each dose",
      "administered as a 30-minute IV infusion. Median treatment duration",
      "6 days (IQR 4-7). TDM at the third or fourth dose."
    ),
    regions         = "Single centre (University of California San Francisco NICU, USA). Data collected November 2007 to March 2010.",
    creatinine_summary = "Serum creatinine on PNA day 1: median 0.9 mg/dL, IQR 0.8-1.2 mg/dL (Table 1).",
    notes           = paste(
      "Retrospective chart-review cohort, n = 29 term neonates (5 excluded:",
      "3 ECMO, 1 cardiomyopathy, 1 incomplete dosing record). 47 gentamicin",
      "concentrations available: 18 patients had paired peak + trough, 11",
      "had trough only. Initial trough median 1.8 mg/L (IQR 1.4-2.6 mg/L);",
      "initial peak median 11.3 mg/L (IQR 9.9-12.6 mg/L). All initial peaks",
      ">= 6 mg/L; 38% had initial trough >= 2 mg/L. Baseline demographics",
      "from Table 1; PK estimates from Table 2 (final model). Fit in NONMEM",
      "VII with FOCE-I."
    )
  )

  ini({
    # Structural PK parameters (Frymoyer 2013 Table 2 'Final Model' column),
    # reported at the typical neonate body weight 3.3 kg and SCr 0.9 mg/dL.
    # All log-transformed.
    lcl <- log(0.111); label("Clearance at reference BW = 3.3 kg, SCr = 0.9 mg/dL (L/h)") # Frymoyer 2013 Table 2 final-model 'Typical Value of CL' = 0.111 L/h (RSE 4.2%)
    lvc <- log(1.56);  label("Volume of distribution at reference BW = 3.3 kg (L)")        # Frymoyer 2013 Table 2 final-model 'V' = 1.56 L (RSE 4.7%)

    # Allometric birth-weight exponents. Frymoyer 2013 Methods 'Population
    # Pharmacokinetic Analysis': 'the effect of birthweight on clearance and
    # volume of distribution was implemented using an allometric model with
    # the exponent defining the relationship fixed to 0.75 and 1,
    # respectively.' Reported without uncertainty (canonical fixed exponents).
    e_wt_birth_cl <- fixed(0.75); label("Allometric birth-weight exponent on CL (unitless)") # Frymoyer 2013 Methods: BW exponent on CL fixed at 0.75
    e_wt_birth_vc <- fixed(1.0);  label("Allometric birth-weight exponent on Vc (unitless)")  # Frymoyer 2013 Methods: BW exponent on V fixed at 1

    # Power effect of serum creatinine (PNA day 1) on CL, centred on the
    # cohort median 0.9 mg/dL. Table 2 prints the exponent magnitude (0.566)
    # while the direction is conveyed by the Results narrative: 'a typical
    # neonate with SCr 1.4 mg/dL had a clearance 27% lower than a neonate
    # with SCr 0.8 mg/dL'. (1.4/0.8)^(-0.566) = 0.729 -> 27.1% lower, so
    # the model exponent applied to (CREAT/0.9) is -0.566.
    e_creat_cl <- -0.566; label("Power exponent of (CREAT/0.9 mg/dL) on CL (unitless)")    # Frymoyer 2013 Table 2 final-model 'Exponent accounting for SCr effect on CL' = 0.566 (RSE 21.9%); sign assigned per Results narrative

    # Inter-individual variability (Frymoyer 2013 Table 2 final-model
    # 'Inter-individual variability of CL, %CV' = 16.1%). Reported as
    # %CV; for log-normal eta the variance on the log scale is
    # omega^2 = log(1 + CV^2) = log(1 + 0.161^2) = 0.02559. IIV on V
    # could not be estimated with precision and was removed (Results
    # 'Population Pharmacokinetic Analysis' paragraph 1).
    etalcl ~ 0.02559 # Frymoyer 2013 Table 2: IIV CL 16.1% CV -> log(1 + 0.161^2) = 0.02559

    # Residual error. Frymoyer 2013 Table 2 final-model 'Residual
    # variability, %CV' = 16.2% (RSE 24.1%). A proportional residual
    # error model was selected (Results 'Population Pharmacokinetic
    # Analysis': 'The residual variability ... was best described by a
    # proportional error model.').
    propSd <- 0.162; label("Proportional residual error (fraction)") # Frymoyer 2013 Table 2 final-model proportional residual = 16.2% CV
  })

  model({
    # ----- Derived covariate terms (centred on cohort medians) -----
    # BW reference 3.3 kg (typical neonate per Abstract and Discussion;
    # cohort median 3.32 kg per Table 1). CREAT reference 0.9 mg/dL
    # (cohort median per Table 1 IQR 0.8-1.2).
    f_creat_cl <- (CREAT / 0.9) ^ e_creat_cl

    # ----- Individual PK parameters -----
    cl <- exp(lcl + etalcl) * (WT_BIRTH / 3.3) ^ e_wt_birth_cl * f_creat_cl
    vc <- exp(lvc)          * (WT_BIRTH / 3.3) ^ e_wt_birth_vc

    kel <- cl / vc

    # ----- ODE system: one-compartment IV PK with first-order elimination -----
    # Gentamicin was administered as a 30-minute IV infusion in the source
    # paper; the model does not hard-code the infusion duration so users
    # specify rate / dur per dose in their event table.
    d/dt(central) <- -kel * central

    # ----- Observation and error -----
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
