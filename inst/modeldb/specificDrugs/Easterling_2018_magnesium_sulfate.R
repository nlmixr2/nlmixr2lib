Easterling_2018_magnesium_sulfate <- function() {
  description <- "One-compartment population PK model of magnesium sulfate (MgSO4-7H2O) with intravenous administration and an endogenous baseline magnesium term added to the administered drug, in pregnant women with severe preeclampsia comparing continuous IV infusion vs serial IV bolus dosing (Easterling 2018)."
  reference <- "Easterling T, Hebert M, Bracken H, Darwish E, Ramadan MC, Shaarawy S, Charles D, Abdel-Aziz T, Nasr AS, Safwal SM, Winikoff B. A randomized trial comparing the pharmacology of magnesium sulfate when used to treat severe preeclampsia with serial intravenous boluses versus a continuous intravenous infusion. BMC Pregnancy Childbirth 2018;18:290. doi:10.1186/s12884-018-1919-6"
  vignette  <- "Easterling_2018_magnesium_sulfate"
  units     <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Maternal body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Exponential effect on Vc centered on the study mean of 90.54 kg per Easterling 2018 Methods ('Weight adjustments were normalized for the mean body weight in the study (90.54 kg)') and Table 4 (Weight effect on volume exponential = -0.426). Encoded here as vc = exp(lvc) * exp(e_wt_vc * (WT/90.54 - 1)).",
      source_name        = "WT"
    ),
    CREAT = list(
      description        = "Maternal serum creatinine at the start of the PK study",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Linear effect on CL centered on the study mean of 66.3 umol/L (= 0.75 mg/dL) per Easterling 2018 Table 1 (initial serum creatinine) and Table 4 (linear coefficient = 0.553). Encoded here as cl = exp(lcl) * (1 + e_creat_cl * (CREAT/66.3 - 1)). The Easterling 2018 text identifies the covariate as 'serum creatinine at the start of the study (linear)' on CL; Table 4 reports the coefficient verbatim as +0.553. With this direct-ratio linear form the as-reported positive sign implies CL increases with serum creatinine; the prior Salinger 2013 paper from the same author group used the inverse-ratio power form (0.8/CREAT)^theta that yields the conventional renal-function direction (higher SCR -> lower CL). See the vignette's Assumptions and deviations section for the discussion. CREAT values in the paper are reported in both umol/L (mean 66.3) and mg/dL (mean 0.75); 1 mg/dL = 88.4 umol/L.",
      source_name        = "CREAT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 200L,
    n_studies      = 1L,
    age_range      = "mean 29 +/- 6 years (both arms)",
    age_median     = "~29 years",
    weight_range   = "mean 91.4 +/- 17.2 kg (serial IV bolus arm) and 89.6 +/- 14.8 kg (continuous infusion arm); pooled study mean 90.54 kg",
    weight_median  = "~90 kg",
    sex_female_pct = 100,
    race_ethnicity = "Egyptian (women enrolled at El Galaa Teaching Hospital in Cairo and Shatby Maternity Hospital in Alexandria)",
    disease_state  = "Severe preeclampsia (systolic BP >= 160 mmHg or diastolic BP >= 110 mmHg with >= 1+ proteinuria); gestational age mean ~35.7 +/- 2.8 weeks (serial IV bolus) and ~35.2 +/- 3.3 weeks (continuous infusion); exclusion of subjects with serum creatinine > 106 umol/L (1.2 mg/dL); 0-2% postpartum at enrolment; 4-5% multiple pregnancy; 24-46% primigravid",
    dose_range     = "Continuous IV infusion arm: 4 g MgSO4-7H2O loading dose IV over ~20 min, then 1 g/h IV continuous infusion for ~12 h. Serial IV bolus arm: 6 g MgSO4-7H2O loading dose IV over ~30 min, then 2 g IV bolus over 10 min every 2 h via Springfusor pump for ~12 h.",
    regions        = "Egypt (Cairo and Alexandria)",
    notes          = "Open-label randomised trial (NCT02091401) of two MgSO4 administration regimens for preeclampsia, enrolled January 2015 - February 2016. 200 women randomised 1:1; six sparse PK samples per subject (plus baseline) over a 12-h treatment period; actual administration and sampling times used in the analysis. 1347 magnesium concentrations from 200 women; 16 excluded for contamination or mislabelling, leaving 1331 for the pop-PK fit. Estimation in NONMEM v7.3.0. Baseline demographics in Table 1 (this script); pop-PK estimates in Table 4. Doses are administered as MgSO4-7H2O (heptahydrate, MW 246.47); serum concentrations are reported as elemental Mg. The PK parameters in this file are expressed in dose units of mg of elemental Mg and concentration units of mg/L of total magnesium (administered + endogenous baseline); convert administered MgSO4-7H2O grams to mg Mg by multiplying by 24.305/246.47 = 0.0986 (e.g., 4 g MgSO4-7H2O = 394.4 mg Mg, 1 g MgSO4-7H2O/h = 98.6 mg Mg/h infusion rate, 2 g bolus = 197.2 mg Mg). The Easterling 2018 dataset has only intravenous dosing (continuous infusion or short serial bolus); both are administered into central via the rxode2 events file (cmt = central, evid = 1, rate set per dosing record)."
  )

  ini({
    # Structural parameters from Easterling 2018 Table 4 (Final PK Model). Reference
    # subject: 90.54 kg maternal body weight, 66.3 umol/L (= 0.75 mg/dL) serum creatinine.
    # Dose units used here: mg of elemental Mg (multiply g of MgSO4-7H2O by 98.6).
    # Concentration units: mg/L of total magnesium (administered + endogenous BL).
    # Conversions from Table 4 (paper reports CL in mL/min and BL in mmol/L):
    #   CL  = 56.3 mL/min  = 56.3 * 60 / 1000 = 3.378 L/h
    #   BL  = 0.925 mmol/L * 24.305 mg/mmol   = 22.48 mg/L
    lcl <- log(3.378); label("Clearance for the reference subject (CL, L/h)")           # Easterling 2018 Table 4 (56.3 mL/min)
    lvc <- log(65.3);  label("Volume of distribution for the reference subject (V, L)") # Easterling 2018 Table 4
    lbl <- log(22.48); label("Endogenous steady-state baseline magnesium (BL, mg/L)")   # Easterling 2018 Table 4 (0.925 mmol/L = 2.25 mg/dL)

    # Covariate effects from Easterling 2018 Table 4. Centering references are
    # taken from Methods text (WT_REF = 90.54 kg study mean) and Table 1 (CREAT_REF
    # = 66.3 umol/L study mean of initial serum creatinine).
    e_wt_vc    <- -0.426; label("Exponential coefficient of body weight on Vc (unitless)")  # Easterling 2018 Table 4 (-0.426 +/- 0.196)
    e_creat_cl <-  0.553; label("Linear coefficient of serum creatinine on CL (unitless)")  # Easterling 2018 Table 4 (+0.553 +/- 0.138)

    # Inter-individual variability was estimated and used by the paper's Monte
    # Carlo simulations ('The simulations utilized the final population model
    # parameter estimates of fixed effects as well as inter-individual and
    # residual random variability'), but the variance estimates are not reported
    # in Table 4 nor anywhere else in the published manuscript. Following the
    # convention of the matched Salinger 2013 model in this library (same author
    # group, also lacking reported IIV estimates), the file omits IIV terms;
    # simulations using this model produce the typical-value response. The
    # vignette's Assumptions and deviations section documents this gap.

    # Residual error from Easterling 2018 Table 4. The paper reports 'Residual
    # variability: 0.0387 +/- 0.0013' with no explicit form statement, but the
    # small relative SE (0.0013 / 0.0387 ~= 3.4%) is consistent with a NONMEM
    # $SIGMA variance estimate (asymptotic SE/var ~= sqrt(2/N) for N ~= 1331
    # observations gives ~3.9%, in good agreement). Treated here as a
    # proportional residual variance: SD = sqrt(0.0387) = 0.197 (~ 19.7% CV).
    propSd <- 0.197; label("Proportional residual error (SD, fraction of Cc)") # Easterling 2018 Table 4 (variance 0.0387)
  })
  model({
    # Individual PK parameters. Reference subject: WT = 90.54 kg, CREAT = 66.3 umol/L.
    # Covariate forms per Easterling 2018 Methods + Table 4:
    #   V_i  = V  * exp(e_wt_vc    * (WT_i    / 90.54 - 1))   (exponential, centered+normalised)
    #   CL_i = CL * (1 + e_creat_cl * (CREAT_i / 66.3  - 1))   (linear, centered+normalised)
    cl <- exp(lcl) * (1 + e_creat_cl * (CREAT / 66.3 - 1))
    vc <- exp(lvc) * exp(e_wt_vc * (WT / 90.54 - 1))
    bl <- exp(lbl)

    # One-compartment model with first-order elimination. Easterling 2018 dosed
    # only intravenously (no IM/oral route in this study); all doses go directly
    # into central via cmt = central, evid = 1 (instantaneous IV bolus) or a
    # non-zero rate (zero-order IV infusion) on the dose record.
    kel <- cl / vc

    d/dt(central) <- -kel * central

    # Observation: administered magnesium concentration (central / vc) plus the
    # endogenous baseline (bl). Concentration units are mg/L of elemental Mg.
    # Per Methods: 'The administered magnesium was modeled as additive to BL.'
    Cc <- central / vc + bl

    Cc ~ prop(propSd)
  })
}
