Takechi_2018_propranolol <- function() {
  description <- "One-compartment first-order absorption population PK model for oral propranolol in Japanese infants with infantile hemangioma (35-150 days postnatal age), with fixed allometric body-weight scaling and a power effect of postnatal age on apparent oral clearance; the companion logistic-regression PD model relating exposure (AUC), treatment duration, and gestational age to treatment-success probability is reproduced in the validation vignette."
  reference <- "Takechi T, Kumokawa T, Kato R, Higuchi T, Kaneko T, Ieiri I. Population Pharmacokinetics and Pharmacodynamics of Oral Propranolol in Pediatric Patients With Infantile Hemangioma. J Clin Pharmacol. 2018;58(10):1361-1370. doi:10.1002/jcph.1149"
  vignette <- "Takechi_2018_propranolol"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight (kg), time-varying. Used for fixed-exponent allometric scaling on CL (0.75) and V (1) normalized to the study median 6.115 kg.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Median 6.115 kg (range 3.15-8.71) at study start (Table 1). The paper centres the allometric power function on this median.",
      source_name        = "body weight"
    ),
    PNA = list(
      description        = "Postnatal age (days), time-varying. Power effect on apparent oral clearance normalized to the study median 113 days, with fixed exponent 1 (Table 2 reports 'PowerAGE for CL/F = 1 fixed').",
      units              = "days",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Median 113 days (range 53-150) at study start (Table 1). Interpreted to capture postnatal maturation of CYP1A2 / CYP2D6 / CYP2C19 / UGT pathways relevant to propranolol metabolism. The paper says 'The effect of postnatal age on CL/F was fixed at 1 through a backward elimination step' (Results, page 1365).",
      source_name        = "postnatal age"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 32L,
    n_studies      = 1L,
    age_range      = "53-150 days postnatal (35-150 days at enrollment per inclusion criteria)",
    age_median     = "113 days postnatal",
    weight_range   = "3.15-8.71 kg",
    weight_median  = "6.115 kg",
    sex_female_pct = 71.9,
    disease_state  = "infantile hemangioma (proliferating target lesion minimum diameter 1.5 cm)",
    dose_range     = "3 mg/kg/day oral propranolol solution (3.75 mg/mL base) divided into 2 administrations, after a 1-1-1 mg/kg titration during the first days (titration to 3 mg/kg in 1-mg/kg increments every other day); treatment 24 weeks",
    regions        = "Japan (multicenter open-label phase 3 across 13 sites)",
    n_observations = 63L,
    notes          = "PK dataset: 63 plasma propranolol concentration-time records from 32 patients during the first 12 weeks (sparse: first sample 1/2/3/4/6 h after day-1 maintenance morning dose; second sample 2 h after morning intake at week 12). PD dataset: 64 success/failure assessments at weeks 12 and 24 (success rates 43.8% at week 12, 78.1% at week 24). Gestational age at birth (median 272 days, range 213-293; 4/32 preterm) enters the companion PD logistic-regression model reproduced in the validation vignette but is not used by the PK structural model. Demographics from Table 1."
  )

  ini({
    # ---- Structural PK (Takechi 2018 Table 2, regression equations on page 1364) ----
    lcl <- log(9.34); label("Apparent oral clearance at WT = 6.115 kg and PNA = 113 days (L/h)") # Table 2 CL/F = 9.34 L/h
    lvc <- log(146);  label("Apparent volume of distribution at WT = 6.115 kg (L)")              # Table 2 V/F  = 146 L
    lka <- fixed(log(1.03)); label("First-order oral absorption rate (1/h)")                     # Table 2 ka  = 1.03 fixed (calculated from tmax equation using upstream parameters)

    # ---- Fixed structural covariate exponents (Takechi 2018 Table 2) ----
    e_wt_cl  <- fixed(0.75); label("Allometric exponent of body weight on apparent oral clearance (unitless)") # Table 2 Power_WT for CL/F = 0.75 fixed (theoretical allometric)
    e_wt_vc  <- fixed(1);    label("Allometric exponent of body weight on apparent volume of distribution (unitless)") # Table 2 Power_WT for V/F = 1 fixed (theoretical allometric)
    e_pna_cl <- fixed(1);    label("Postnatal-age power exponent on apparent oral clearance (unitless)") # Table 2 Power_AGE for CL/F = 1 fixed via backward elimination

    # ---- IIV (Takechi 2018 Table 2, "IIV CL/F (%) = 76.9 (29.3)") ----
    # NONMEM exponential variance: omega^2 = log(CV^2 + 1) for lognormal CV%
    # CV = 0.769 -> omega^2 = log(1 + 0.769^2) = log(1.5914) = 0.4646
    etalcl ~ 0.4646  # Table 2 IIV CL/F = 76.9%

    # ---- Residual error (Takechi 2018 Table 2, "Proportional error (%) = 34.0") ----
    propSd <- 0.34; label("Proportional residual error on propranolol concentration (fraction)") # Table 2 proportional error = 34.0%
  })

  model({
    # ---- Individual PK parameters (Takechi 2018 page 1364 regression equations) ----
    # CL/F (L/h) = 9.34 * (body weight / 6.115)^0.75 * (postnatal age / 113) * exp(eta_CL)
    # V/F  (L)   = 146  * (body weight / 6.115)
    # ka   (1/h) = 1.03
    cl <- exp(lcl + etalcl) * (WT / 6.115)^e_wt_cl * (PNA / 113)^e_pna_cl
    vc <- exp(lvc)          * (WT / 6.115)^e_wt_vc
    ka <- exp(lka)

    # ---- Micro-constants ----
    kel <- cl / vc

    # ---- ODE system (1-compartment with first-order absorption, Supplemental Figure S2) ----
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # ---- Observation ----
    Cc <- central / vc * 1000  # mg/L -> ng/mL (dose in mg, vc in L)
    Cc ~ prop(propSd)
  })
}
