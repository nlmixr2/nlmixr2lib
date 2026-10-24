Rambiritch_2016_glibenclamide_dose_glucoseprofile <- function() {
  description <- "Dose-response model (Model 5 of Rambiritch 2016) for the effect of once-daily oral glibenclamide on the full blood glucose versus time profile over an 8-hour sampling day in South African adults with poorly controlled type 2 diabetes. The drug-free baseline is an empirical two-harmonic Fourier series with an 8-hour period (Krzyzanski 2000 method) capturing the breakfast and lunch glucose excursions; an inhibitory Emax model driven directly by the daily dose scales the whole profile: glucose = E0(t) * (1 - DOSE * Emax / (DOSE + ED50)). No PK layer. Sibling models from the same paper: Rambiritch_2016_glibenclamide_cpss_glucoseprofile (steady-state concentration as the driver) and the fasting-glucose and mean-glucose models; the companion population PK model is Rambiritch_2016_glibenclamide."
  reference <- paste(
    "Rambiritch V, Naidoo P, Pillai G. Glibenclamide population",
    "pharmacokinetic/pharmacodynamic modeling in South African type 2 diabetic",
    "subjects. Clin Pharmacol. 2016;8:141-150. doi:10.2147/CPAA.S102674.",
    "Parameter values are the 'Dose-full glucose profile model' column of",
    "Table 3; the structural equations are Equations 1 and 4-7.",
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
    notes = "Rambiritch 2016 Patients and methods: 24 subjects; subject 16 did not escalate beyond 5 mg/day (hypoglycaemia) and subjects 20, 21 and 24 had incomplete data sets, but all 24 were included in the PKPD data set. The paper does not reprint baseline demographics; the age and weight ranges are those of the same study cohort as printed in the companion population PK paper (Rambiritch 2016, Clin Pharmacol 2016;8:83-92, Table 1, which characterises 22 subjects at entry: 20 female, 2 male). Glucose was measured with a glucometer over an 8-hour sampling day at each dose level, with breakfast at time 0 and lunch about 4 hours later; the profile shows two glucose peaks at about 2 and 6 hours."
  )

  ini({
    le0 <- log(15.90); label("Constant term A0 of the drug-free glucose baseline (mmol/L)") # Rambiritch 2016 Table 3, Dose-full glucose profile model: A0 = 15.90 mmol/L (RSE 4.71%)

    # Fourier coefficients of Equation 4. Signed and on the linear scale, with
    # additive between-subject variability (Equation 6). The paper writes A for
    # the cosine and B for the sine coefficients; harmonic 2 uses 4*pi*t/P.
    e0_glucose_fourier_cos1 <- 0.95; label("Cosine coefficient of the first harmonic of the glucose baseline, A1 (mmol/L)") # Rambiritch 2016 Table 3, Dose-full glucose profile model: A1 = 0.95 (RSE 17.44%)
    e0_glucose_fourier_sin1 <- 1.45; label("Sine coefficient of the first harmonic of the glucose baseline, B1 (mmol/L)") # Rambiritch 2016 Table 3, Dose-full glucose profile model: B1 = 1.45 (RSE 36.34%)
    e0_glucose_fourier_cos2 <- -1.95; label("Cosine coefficient of the second harmonic of the glucose baseline, A2 (mmol/L)") # Rambiritch 2016 Table 3, Dose-full glucose profile model: A2 = -1.95 (RSE -14.26%)
    e0_glucose_fourier_sin2 <- 2.63; label("Sine coefficient of the second harmonic of the glucose baseline, B2 (mmol/L)") # Rambiritch 2016 Table 3, Dose-full glucose profile model: B2 = 2.63 (RSE 8.75%)
    e0_glucose_fourier_period <- fixed(8); label("Period of the Fourier glucose baseline (h)") # Rambiritch 2016 Equation 4 text: 'Period was 8 hours, ie, the time during which glucose concentrations were measured'

    lemax <- log(0.34); label("Maximum fractional inhibition of blood glucose Emax (unitless)") # Rambiritch 2016 Table 3, Dose-full glucose profile model: Emax = 0.34 (RSE 14.97%)
    led50 <- log(2.21); label("Daily dose producing half-maximal inhibition ED50 (mg/day)") # Rambiritch 2016 Table 3, Dose-full glucose profile model: ED50 = 2.21 mg (RSE 29.46%)

    # Exponential BSV (Equation 5); BSV (%CV) column converted as
    # omega^2 = (CV/100)^2, the sqrt(omega^2) reading the paper applies to
    # every BSV row (see the vignette).
    etale0 ~ 0.0824264 # Rambiritch 2016 Table 3: BSV A0 28.71% CV -> 0.2871^2
    # Additive BSV on the four Fourier coefficients (Equation 6) with one
    # variance shared by all four ('the omega matrix for all the coefficients
    # was constrained to the same value'). The 120.83 printed on every row is
    # sqrt(omega^2) x 100 of that additive eta, i.e. an SD of 1.2083 mmol/L.
    etae0_glucose_fourier_cos1 ~ 1.4599889 # Rambiritch 2016 Table 3: BSV A1 120.83 -> additive SD 1.2083 mmol/L, squared
    etae0_glucose_fourier_sin1 ~ 1.4599889 # Rambiritch 2016 Table 3: BSV B1 120.83 -> additive SD 1.2083 mmol/L, squared
    etae0_glucose_fourier_cos2 ~ 1.4599889 # Rambiritch 2016 Table 3: BSV A2 120.83 -> additive SD 1.2083 mmol/L, squared
    etae0_glucose_fourier_sin2 ~ 1.4599889 # Rambiritch 2016 Table 3: BSV B2 120.83 -> additive SD 1.2083 mmol/L, squared
    etalemax ~ 0.0033063 # Rambiritch 2016 Table 3: BSV Emax 5.75% CV -> 0.0575^2
    etaled50 ~ 0.0228917 # Rambiritch 2016 Table 3: BSV ED50 15.13% CV -> 0.1513^2

    propSd <- 0.13; label("Proportional residual error on blood glucose (fraction)") # Rambiritch 2016 Table 3: residual variance 0.02 (13% CV)
  })

  model({
    a0 <- exp(le0 + etale0)
    fcos1 <- e0_glucose_fourier_cos1 + etae0_glucose_fourier_cos1
    fsin1 <- e0_glucose_fourier_sin1 + etae0_glucose_fourier_sin1
    fcos2 <- e0_glucose_fourier_cos2 + etae0_glucose_fourier_cos2
    fsin2 <- e0_glucose_fourier_sin2 + etae0_glucose_fourier_sin2
    emax <- exp(lemax + etalemax)
    ed50 <- exp(led50 + etaled50)

    # Rambiritch 2016 Equation 4. Time is hours from the start of the 8-hour
    # sampling day (the time-0 fasting sample, breakfast). Because 24 h is a
    # whole number of 8-hour periods, a time axis running across days stays in
    # phase as long as each day's sampling starts at the same clock time; the
    # series is only meaningful inside each day's 0-8 h window.
    ang <- 2 * pi * t / e0_glucose_fourier_period
    e0 <- a0 + fcos1 * cos(ang) + fsin1 * sin(ang) + fcos2 * cos(2 * ang) + fsin2 * sin(2 * ang)

    # Rambiritch 2016 Equation 1 with the Equation 4 baseline (Model 5)
    glucose <- e0 * (1 - DOSE * emax / (DOSE + ed50))

    # Glucose was log-transformed and fitted with an additive log-scale error;
    # Equation 7 prints the equivalent proportional form Cij = C*ij (1 + eps_ij).
    glucose ~ prop(propSd)
  })
}
