Rambiritch_2016_glibenclamide_dose_meanglucose <- function() {
  description <- "Dose-response model (Model 3 of Rambiritch 2016) for the effect of once-daily oral glibenclamide on the mean blood glucose concentration over the 8-hour sampling day in South African adults with poorly controlled type 2 diabetes. Inhibitory Emax model driven directly by the daily dose: Gmean = E0 * (1 - DOSE * Emax / (DOSE + ED50)). No PK layer; the dose-escalation study was sampled at steady state on each dose level. Sibling models from the same paper: Rambiritch_2016_glibenclamide_cpss_meanglucose (steady-state concentration as the driver) and the fasting-glucose and full-glucose-profile models; the companion population PK model is Rambiritch_2016_glibenclamide."
  reference <- paste(
    "Rambiritch V, Naidoo P, Pillai G. Glibenclamide population",
    "pharmacokinetic/pharmacodynamic modeling in South African type 2 diabetic",
    "subjects. Clin Pharmacol. 2016;8:141-150. doi:10.2147/CPAA.S102674.",
    "Parameter values are the 'Dose-mean glucose model' column of Table 2;",
    "the structural equation is Equation 1 with the mean glucose concentration as the response.",
    sep = " "
  )
  vignette <- "Rambiritch_2016_glibenclamide_glucose_response"
  units <- list(time = "h", dosing = "mg", concentration = "mmol/L")

  covariateData <- list(
    DOSE = list(
      description = "Current once-daily glibenclamide dose",
      units = "mg/day",
      type = "continuous",
      reference_category = NULL,
      notes = "TIME-VARYING within subject: the study escalated each subject through 0, 2.5, 5, 10 and 20 mg/day at 2-week intervals and sampled glucose at steady state on each level. Set to 0 for the drug-free (placebo) level. Enters directly as the driving force of the inhibitory Emax model (Rambiritch 2016 Equation 1).",
      source_name = "Dose"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 24L,
    n_studies = 1L,
    age_range = "39-73 years (mean 54); from the companion population PK paper on the same cohort",
    weight_range = "42.0-107.8 kg (mean 71.1); from the companion population PK paper on the same cohort",
    sex_female_pct = 91,
    race_ethnicity = "South African; race not reported",
    disease_state = "Poorly controlled type 2 diabetes mellitus requiring oral antidiabetic therapy; mean fasting blood glucose at entry 15.4 mmol/L (range 9.9-21.8), mean HbA1c 12.2%",
    dose_range = "Glibenclamide 0, 2.5, 5, 10 and 20 mg once daily, escalated at 2-week intervals",
    regions = "South Africa (Durban / Chatsworth, KwaZulu-Natal)",
    notes = "Rambiritch 2016 Patients and methods: 24 subjects; subject 16 did not escalate beyond 5 mg/day (hypoglycaemia) and subjects 20, 21 and 24 had incomplete data sets, but all 24 were included in the PKPD data set. The paper does not reprint baseline demographics; the age and weight ranges are those of the same study cohort as printed in the companion population PK paper (Rambiritch 2016, Clin Pharmacol 2016;8:83-92, Table 1, which characterises 22 subjects at entry: 20 female, 2 male). Glucose was measured with a glucometer over an 8-hour sampling day (breakfast at time 0, lunch about 4 hours later) at each dose level; the mean glucose concentration is the area under the glucose-time curve over the sampling day divided by its duration (Rambiritch 2016, Model 3 description)."
  )

  ini({
    le0 <- log(16.70); label("Baseline (drug-free) mean blood glucose over the sampling day E0 (mmol/L)") # Rambiritch 2016 Table 2, Dose-mean glucose model: E0 = 16.70 mmol/L (RSE 5.74%)
    lemax <- log(0.34); label("Maximum fractional inhibition of mean blood glucose Emax (unitless)") # Rambiritch 2016 Table 2, Dose-mean glucose model: Emax = 0.34 (RSE 14.32%)
    led50 <- log(1.85); label("Daily dose producing half-maximal inhibition ED50 (mg/day)") # Rambiritch 2016 Table 2, Dose-mean glucose model: ED50 = 1.85 mg (RSE 40.49%)

    # BSV (%CV) column converted as omega^2 = (CV/100)^2, the same sqrt(omega^2)
    # reading the paper applies to every BSV row (see the vignette).
    etale0 ~ 0.0732785 # Rambiritch 2016 Table 2: BSV E0 27.07% CV -> 0.2707^2
    etalemax ~ 0.259998 # Rambiritch 2016 Table 2: BSV Emax 50.99% CV -> 0.5099^2
    etaled50 ~ 1.4500976 # Rambiritch 2016 Table 2: BSV ED50 120.42% CV -> 1.2042^2

    propSd <- 0.0875; label("Proportional residual error on mean blood glucose (fraction)") # Rambiritch 2016 Table 2: residual variance 0.008 (8.75% CV); sqrt(0.008) = 0.0894 agrees with the printed CV to rounding of the variance
  })

  model({
    e0 <- exp(le0 + etale0)
    emax <- exp(lemax + etalemax)
    ed50 <- exp(led50 + etaled50)

    # Rambiritch 2016 Equation 1 (Model 3: Effect is the mean glucose concentration): Effect = E0 x (1 - [Dose x Emax]/[Dose + ED50])
    glucose_mean <- e0 * (1 - DOSE * emax / (DOSE + ed50))

    # Glucose was log-transformed and fitted with an additive log-scale error;
    # Equation 7 prints the equivalent proportional form Cij = C*ij (1 + eps_ij).
    glucose_mean ~ prop(propSd)
  })
}
