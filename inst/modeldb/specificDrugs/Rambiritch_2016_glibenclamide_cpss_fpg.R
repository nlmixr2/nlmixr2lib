Rambiritch_2016_glibenclamide_cpss_fpg <- function() {
  description <- "Exposure-response model (Model 2 of Rambiritch 2016) for the effect of once-daily oral glibenclamide on fasting blood glucose in South African adults with poorly controlled type 2 diabetes. Inhibitory Emax model driven by the average steady-state plasma concentration Cpss = DOSE / (CL/F x 24), with CL/F fixed from the companion population PK model (Rambiritch_2016_glibenclamide). Sibling models from the same paper: Rambiritch_2016_glibenclamide_dose_fpg (dose as the driver) and the mean-glucose and full-glucose-profile models."
  reference <- paste(
    "Rambiritch V, Naidoo P, Pillai G. Glibenclamide population",
    "pharmacokinetic/pharmacodynamic modeling in South African type 2 diabetic",
    "subjects. Clin Pharmacol. 2016;8:141-150. doi:10.2147/CPAA.S102674.",
    "Parameter values are the 'Cpss-fasting glucose model' column of Table 1;",
    "the structural equations are Equations 2 and 3. CL/F and its",
    "between-subject variability are from Rambiritch V, Naidoo P, Maharaj B,",
    "Pillai G. Population pharmacokinetic modeling of glibenclamide in poorly",
    "controlled South African type 2 diabetic subjects. Clin Pharmacol.",
    "2016;8:83-92. doi:10.2147/CPAA.S102676 (Table 3).",
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
      notes = "TIME-VARYING within subject: the study escalated each subject through 0, 2.5, 5, 10 and 20 mg/day at 2-week intervals and sampled glucose at steady state on each level. Set to 0 for the drug-free (placebo) level. Converted to the average steady-state concentration Cpss = DOSE / (CL/F x 24) (Rambiritch 2016 Equation 3), which drives the inhibitory Emax model.",
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
    # CL/F is not re-estimated here: Equation 3 uses the CL/F of the companion
    # population PK model (Rambiritch 2016 reference 9), which the Table 1
    # footnote prints as 2.16 L/h.
    lcl <- fixed(log(2.16)); label("Apparent clearance CL/F of glibenclamide (L/h)") # Rambiritch 2016 Table 1 footnote: CL/f = 2.16 L/h (from the companion population PK model, its Table 3)

    le0 <- log(14.30); label("Baseline (drug-free) fasting blood glucose E0 (mmol/L)") # Rambiritch 2016 Table 1, Cpss-fasting glucose model: E0 = 14.30 mmol/L (RSE 6.04%)
    lemax <- log(0.39); label("Maximum fractional inhibition of fasting blood glucose Emax (unitless)") # Rambiritch 2016 Table 1, Cpss-fasting glucose model: Emax = 0.39 (RSE 22.64%)
    lec50 <- log(85.20); label("Average steady-state concentration producing half-maximal inhibition EC50 (ng/mL)") # Rambiritch 2016 Table 1, Cpss-fasting glucose model: EC50 = 85.20 (RSE 26.29%); printed header 'mg/L' is ug/L = ng/mL, see the vignette

    # The companion PK model's BSV on CL/F carries the between-subject PK
    # variability that the Cpss driver was introduced to separate from PD
    # variability; held fixed because it was not estimated in this analysis.
    etalcl ~ fixed(0.1149888) # Companion popPK paper (doi:10.2147/CPAA.S102676) Table 3: BSV CL/F 33.91% CV -> 0.3391^2
    # BSV (%CV) column converted as omega^2 = (CV/100)^2, the same sqrt(omega^2)
    # reading the paper applies to every BSV row (see the vignette).
    etale0 ~ 0.075076 # Rambiritch 2016 Table 1: BSV E0 27.40% CV -> 0.2740^2
    etalemax ~ 0.308025 # Rambiritch 2016 Table 1: BSV Emax 55.50% CV -> 0.5550^2
    etalec50 ~ 1.6e-07 # Rambiritch 2016 Table 1: BSV EC50 0.04% CV -> 0.0004^2 (effectively zero)

    propSd <- 0.14; label("Proportional residual error on fasting blood glucose (fraction)") # Rambiritch 2016 Table 1: residual variance 0.02 (14% CV); sqrt(0.02) = 0.141 agrees with the printed CV to rounding of the variance
  })

  model({
    cl <- exp(lcl + etalcl)
    e0 <- exp(le0 + etale0)
    emax <- exp(lemax + etalemax)
    ec50 <- exp(lec50 + etalec50)

    # Rambiritch 2016 Equation 3: Cpss = Dose / (CL/f x 24). DOSE in mg/day and
    # CL/F in L/h give mg/L; x 1000 converts to ng/mL, the units of EC50.
    cpss <- 1000 * DOSE / (cl * 24)

    # Rambiritch 2016 Equation 2: Effect = E0 x (1 - [Cpss x Emax]/[Cpss + EC50])
    fpg <- e0 * (1 - cpss * emax / (cpss + ec50))

    # Glucose was log-transformed and fitted with an additive log-scale error;
    # Equation 7 prints the equivalent proportional form Cij = C*ij (1 + eps_ij).
    fpg ~ prop(propSd)
  })
}
