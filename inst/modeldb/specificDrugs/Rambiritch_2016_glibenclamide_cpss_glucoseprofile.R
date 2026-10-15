Rambiritch_2016_glibenclamide_cpss_glucoseprofile <- function() {
  description <- "Exposure-response model (Model 6 of Rambiritch 2016) for the effect of once-daily oral glibenclamide on the full blood glucose versus time profile over an 8-hour sampling day in South African adults with poorly controlled type 2 diabetes. The drug-free baseline is an empirical two-harmonic Fourier series with an 8-hour period (Krzyzanski 2000 method) capturing the breakfast and lunch glucose excursions; an inhibitory Emax model driven by the average steady-state plasma concentration Cpss = DOSE / (CL/F x 24) scales the whole profile: glucose = E0(t) * (1 - Cpss * Emax / (Cpss + EC50)), with CL/F fixed from the companion population PK model (Rambiritch_2016_glibenclamide). Sibling models from the same paper: Rambiritch_2016_glibenclamide_dose_glucoseprofile (dose as the driver) and the fasting-glucose and mean-glucose models."
  reference <- paste(
    "Rambiritch V, Naidoo P, Pillai G. Glibenclamide population",
    "pharmacokinetic/pharmacodynamic modeling in South African type 2 diabetic",
    "subjects. Clin Pharmacol. 2016;8:141-150. doi:10.2147/CPAA.S102674.",
    "Parameter values are the 'Cpss-full glucose profile model' column of",
    "Table 3; the structural equations are Equations 2-7. CL/F and its",
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
    notes = "Rambiritch 2016 Patients and methods: 24 subjects; subject 16 did not escalate beyond 5 mg/day (hypoglycaemia) and subjects 20, 21 and 24 had incomplete data sets, but all 24 were included in the PKPD data set. The paper does not reprint baseline demographics; the age and weight ranges are those of the same study cohort as printed in the companion population PK paper (Rambiritch 2016, Clin Pharmacol 2016;8:83-92, Table 1, which characterises 22 subjects at entry: 20 female, 2 male). Glucose was measured with a glucometer over an 8-hour sampling day at each dose level, with breakfast at time 0 and lunch about 4 hours later; the profile shows two glucose peaks at about 2 and 6 hours."
  )

  ini({
    # CL/F is not re-estimated here: Equation 3 uses the CL/F of the companion
    # population PK model (Rambiritch 2016 reference 9), which the Table 3
    # footnote prints as 2.16 L/h.
    lcl <- fixed(log(2.16)); label("Apparent clearance CL/F of glibenclamide (L/h)") # Rambiritch 2016 Table 3 footnote: CL/f = 2.16 L/h (from the companion population PK model, its Table 3)

    le0 <- log(15.40); label("Constant term A0 of the drug-free glucose baseline (mmol/L)") # Rambiritch 2016 Table 3, Cpss-full glucose profile model: A0 = 15.40 mmol/L (RSE 5.40%)

    # Fourier coefficients of Equation 4. Signed and on the linear scale, with
    # additive between-subject variability (Equation 6). The paper writes A for
    # the cosine and B for the sine coefficients; harmonic 2 uses 4*pi*t/P.
    e0_glucose_fourier_cos1 <- 0.97; label("Cosine coefficient of the first harmonic of the glucose baseline, A1 (mmol/L)") # Rambiritch 2016 Table 3, Cpss-full glucose profile model: A1 = 0.97 (RSE 16.36%)
    e0_glucose_fourier_sin1 <- 1.45; label("Sine coefficient of the first harmonic of the glucose baseline, B1 (mmol/L)") # Rambiritch 2016 Table 3, Cpss-full glucose profile model: B1 = 1.45 (RSE 34.21%)
    e0_glucose_fourier_cos2 <- -1.86; label("Cosine coefficient of the second harmonic of the glucose baseline, A2 (mmol/L)") # Rambiritch 2016 Table 3, Cpss-full glucose profile model: A2 = -1.86 (RSE -15.38%)
    e0_glucose_fourier_sin2 <- 2.57; label("Sine coefficient of the second harmonic of the glucose baseline, B2 (mmol/L)") # Rambiritch 2016 Table 3, Cpss-full glucose profile model: B2 = 2.57 (RSE 8.05%)
    e0_glucose_fourier_period <- fixed(8); label("Period of the Fourier glucose baseline (h)") # Rambiritch 2016 Equation 4 text: 'Period was 8 hours, ie, the time during which glucose concentrations were measured'

    lemax <- log(0.31); label("Maximum fractional inhibition of blood glucose Emax (unitless)") # Rambiritch 2016 Table 3, Cpss-full glucose profile model: Emax = 0.31 (RSE 9.20%)
    lec50 <- log(43.60); label("Average steady-state concentration producing half-maximal inhibition EC50 (ng/mL)") # Rambiritch 2016 Table 3, Cpss-full glucose profile model: EC50 = 43.60 (RSE 21.06%); printed header 'mg/L' is ug/L = ng/mL, see the vignette

    # The companion PK model's BSV on CL/F carries the between-subject PK
    # variability that the Cpss driver was introduced to separate from PD
    # variability; held fixed because it was not estimated in this analysis.
    etalcl ~ fixed(0.1149888) # Companion popPK paper (doi:10.2147/CPAA.S102676) Table 3: BSV CL/F 33.91% CV -> 0.3391^2
    # Exponential BSV (Equation 5); BSV (%CV) column converted as
    # omega^2 = (CV/100)^2, the sqrt(omega^2) reading the paper applies to
    # every BSV row (see the vignette).
    etale0 ~ 0.1069944 # Rambiritch 2016 Table 3: BSV A0 32.71% CV -> 0.3271^2
    # Additive BSV on the four Fourier coefficients (Equation 6) with one
    # variance shared by all four ('the omega matrix for all the coefficients
    # was constrained to the same value'). The 7.05 printed on every row is
    # sqrt(omega^2) x 100 of that additive eta, i.e. an SD of 0.0705 mmol/L.
    etae0_glucose_fourier_cos1 ~ 0.0049702 # Rambiritch 2016 Table 3: BSV A1 7.05 -> additive SD 0.0705 mmol/L, squared
    etae0_glucose_fourier_sin1 ~ 0.0049702 # Rambiritch 2016 Table 3: BSV B1 7.05 -> additive SD 0.0705 mmol/L, squared
    etae0_glucose_fourier_cos2 ~ 0.0049702 # Rambiritch 2016 Table 3: BSV A2 7.05 -> additive SD 0.0705 mmol/L, squared
    etae0_glucose_fourier_sin2 ~ 0.0049702 # Rambiritch 2016 Table 3: BSV B2 7.05 -> additive SD 0.0705 mmol/L, squared
    etalemax ~ 0.2060252 # Rambiritch 2016 Table 3: BSV Emax 45.39% CV -> 0.4539^2
    etalec50 ~ 4.8801228 # Rambiritch 2016 Table 3: BSV EC50 220.91% CV -> 2.2091^2

    propSd <- 0.13; label("Proportional residual error on blood glucose (fraction)") # Rambiritch 2016 Table 3: residual variance 0.02 (13% CV)
  })

  model({
    cl <- exp(lcl + etalcl)
    a0 <- exp(le0 + etale0)
    fcos1 <- e0_glucose_fourier_cos1 + etae0_glucose_fourier_cos1
    fsin1 <- e0_glucose_fourier_sin1 + etae0_glucose_fourier_sin1
    fcos2 <- e0_glucose_fourier_cos2 + etae0_glucose_fourier_cos2
    fsin2 <- e0_glucose_fourier_sin2 + etae0_glucose_fourier_sin2
    emax <- exp(lemax + etalemax)
    ec50 <- exp(lec50 + etalec50)

    # Rambiritch 2016 Equation 4. Time is hours from the start of the 8-hour
    # sampling day (the time-0 fasting sample, breakfast). Because 24 h is a
    # whole number of 8-hour periods, a time axis running across days stays in
    # phase as long as each day's sampling starts at the same clock time; the
    # series is only meaningful inside each day's 0-8 h window.
    ang <- 2 * pi * t / e0_glucose_fourier_period
    e0 <- a0 + fcos1 * cos(ang) + fsin1 * sin(ang) + fcos2 * cos(2 * ang) + fsin2 * sin(2 * ang)

    # Rambiritch 2016 Equation 3: Cpss = Dose / (CL/f x 24). DOSE in mg/day and
    # CL/F in L/h give mg/L; x 1000 converts to ng/mL, the units of EC50.
    cpss <- 1000 * DOSE / (cl * 24)

    # Rambiritch 2016 Equation 2 with the Equation 4 baseline (Model 6)
    glucose <- e0 * (1 - cpss * emax / (cpss + ec50))

    # Glucose was log-transformed and fitted with an additive log-scale error;
    # Equation 7 prints the equivalent proportional form Cij = C*ij (1 + eps_ij).
    glucose ~ prop(propSd)
  })
}
