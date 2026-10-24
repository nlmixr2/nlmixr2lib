Rambiritch_2016_glibenclamide_dose_fpg <- function() {
  description <- "Dose-response model (Model 1 of Rambiritch 2016) for the effect of once-daily oral glibenclamide on fasting blood glucose in South African adults with poorly controlled type 2 diabetes. Inhibitory Emax model driven directly by the daily dose: FBG = E0 * (1 - DOSE * Emax / (DOSE + ED50)). No PK layer; the dose-escalation study was sampled at steady state on each dose level. Sibling models from the same paper: Rambiritch_2016_glibenclamide_cpss_fpg (steady-state concentration as the driver) and the mean-glucose and full-glucose-profile models; the companion population PK model is Rambiritch_2016_glibenclamide."
  reference <- paste(
    "Rambiritch V, Naidoo P, Pillai G. Glibenclamide population",
    "pharmacokinetic/pharmacodynamic modeling in South African type 2 diabetic",
    "subjects. Clin Pharmacol. 2016;8:141-150. doi:10.2147/CPAA.S102674.",
    "Parameter values are the 'Dose-fasting glucose model' column of Table 1;",
    "the structural equation is Equation 1.",
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
    notes = "Rambiritch 2016 Patients and methods: 24 subjects; subject 16 did not escalate beyond 5 mg/day (hypoglycaemia) and subjects 20, 21 and 24 had incomplete data sets, but all 24 were included in the PKPD data set. The paper does not reprint baseline demographics; the age and weight ranges are those of the same study cohort as printed in the companion population PK paper (Rambiritch 2016, Clin Pharmacol 2016;8:83-92, Table 1, which characterises 22 subjects at entry: 20 female, 2 male). Glucose was measured with a glucometer over an 8-hour sampling day (breakfast at time 0, lunch about 4 hours later) at each dose level; FBG is the time-0 fasting sample."
  )

  ini({
    le0 <- log(14.10); label("Baseline (drug-free) fasting blood glucose E0 (mmol/L)") # Rambiritch 2016 Table 1, Dose-fasting glucose model: E0 = 14.10 mmol/L (RSE 6.35%)
    lemax <- log(0.37); label("Maximum fractional inhibition of fasting blood glucose Emax (unitless)") # Rambiritch 2016 Table 1, Dose-fasting glucose model: Emax = 0.37 (RSE 37.47%)
    led50 <- log(4.56); label("Daily dose producing half-maximal inhibition ED50 (mg/day)") # Rambiritch 2016 Table 1, Dose-fasting glucose model: ED50 = 4.56 mg (RSE 77.19%)

    # BSV (%CV) column converted as omega^2 = (CV/100)^2, the same sqrt(omega^2)
    # reading the paper applies to every BSV row (see the vignette).
    etale0 ~ 0.0843903 # Rambiritch 2016 Table 1: BSV E0 29.05% CV -> 0.2905^2
    etalemax ~ 0.2369742 # Rambiritch 2016 Table 1: BSV Emax 48.68% CV -> 0.4868^2
    etaled50 ~ 1e-08 # Rambiritch 2016 Table 1: BSV ED50 0.01% CV -> 0.0001^2 (effectively zero)

    propSd <- 0.15; label("Proportional residual error on fasting blood glucose (fraction)") # Rambiritch 2016 Table 1: residual variance 0.02 (15% CV); sqrt(0.02) = 0.141 agrees with the printed CV to rounding of the variance
  })

  model({
    e0 <- exp(le0 + etale0)
    emax <- exp(lemax + etalemax)
    ed50 <- exp(led50 + etaled50)

    # Rambiritch 2016 Equation 1: Effect = E0 x (1 - [Dose x Emax]/[Dose + ED50])
    fpg <- e0 * (1 - DOSE * emax / (DOSE + ed50))

    # Glucose was log-transformed and fitted with an additive log-scale error;
    # Equation 7 prints the equivalent proportional form Cij = C*ij (1 + eps_ij).
    fpg ~ prop(propSd)
  })
}
