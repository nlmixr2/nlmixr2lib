Bienert_2010_propofol <- function() {
  description <- paste(
    "Two-compartment population PK model of propofol with a direct-effect",
    "inhibitory Emax model for the Bispectral Index (BIS), in adult intensive",
    "care unit patients sedated with a prolonged constant intravenous",
    "infusion (Bienert 2010). Metabolic clearance falls linearly with the",
    "hourly-recorded systolic blood pressure, centred at 135 mmHg. BIS is",
    "driven by the plasma concentration with no effect-compartment delay;",
    "baseline BIS and Emax are both fixed to 100, so EC50 is the only",
    "estimated PD parameter. Twenty-four-hour cosine rhythms on clearance and",
    "on EC50 were tested and not retained, so the model carries no circadian",
    "term. Time in minutes."
  )
  reference <- paste(
    "Bienert A, Kusza K, Wawrzyniak K, Grzeskowiak E, Kokot ZJ, Matysiak J,",
    "Grabowski T, Wolc A, Wiczling P, Regulski M (2010). Assessing circadian",
    "rhythms in propofol PK and PD during prolonged infusion in ICU patients.",
    "J Pharmacokinet Pharmacodyn 37(3):289-304. doi:10.1007/s10928-010-9161-5.",
    sep = " "
  )
  vignette <- "Bienert_2010_propofol_circadian"
  units <- list(time = "min", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    SBP = list(
      description = "Systolic blood pressure",
      units = "mmHg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-varying. Recorded every hour at the blood-collection time",
        "points during the second and third day of infusion (Bienert 2010",
        "Methods, 'Patients') and entered as a time-dependent covariate on",
        "metabolic clearance in the linear-deviation form",
        "CL = 2.65 * (1 - 0.00714 * (SBP - 135)) (Table 2, Results). 135 mmHg",
        "is the centring value printed in the equation and the cohort mean of",
        "Table 1 (135.0 +/- 19.4, range 85-191). Over the observed range the",
        "factor runs from 1.36 (85 mmHg) to 0.60 (191 mmHg); clearance would",
        "reach zero at 275 mmHg, so do not extrapolate far beyond the observed",
        "range."
      ),
      source_name = "SBP"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = paste(
        "Screened graphically against the individual PK/PD parameter",
        "estimates (Methods, 'Covariate search and rhythm assessment') and",
        "not retained."
      )
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Screened graphically and not retained."
    ),
    HR = list(
      description = "Heart rate",
      units = "beats/min",
      type = "continuous",
      notes = paste(
        "Time-dependent covariate screened against the weighted residuals",
        "and not retained; systolic blood pressure was the only covariate",
        "kept (Results)."
      )
    ),
    BODYTEMP = list(
      description = "Body temperature",
      units = "degC",
      type = "continuous",
      notes = "Time-dependent covariate screened and not retained."
    )
  )

  compartmentData <- list(
    central = list(
      analyte = "propofol",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "propofol",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 10L,
    n_studies = 1L,
    age_range = "53-85 years",
    age_median = "63 years (mean, SD 11)",
    weight_range = "50-95 kg",
    weight_median = "72.6 kg (mean, SD 22.8)",
    sex_female_pct = 27.3,
    disease_state = paste(
      "Critically ill adult intensive care unit patients sedated with",
      "propofol for more than 24 h. Excluded: neuromuscular blockade, head",
      "trauma, and neurological disorders that could alter the EEG",
      "(metabolic or hypoxic encephalopathy, meningoencephalitis, stroke,",
      "tumours). No vasopressors during the study period."
    ),
    dose_range = paste(
      "Initial dose 1.5 mg/kg followed by a constant maintenance infusion of",
      "0.7-0.9 mg/kg/h (one patient required rate changes), plus occasional",
      "extra boluses. Mean total infusion duration 98.4 +/- 84.5 h; sampling",
      "hourly on days 2-3 (38.2 +/- 11.0 h) and 3-240 min after stopping the",
      "infusion. Mean steady-state concentration 0.48 +/- 0.21 mg/L."
    ),
    regions = "Poland (Poznan)",
    notes = paste(
      "Table 1 demographics are for the 11 enrolled patients (8 male, 3",
      "female); patient 3 had no blood samples and was excluded, so 10",
      "patients contributed PK and BIS data. Mean BMI 25.5 (17.3-33.2).",
      "Plasma propofol by HPLC with fluorescence detection, LOQ 0.01 mg/L.",
      "BIS recorded hourly at the blood-collection times (Philips IntelliVue",
      "MP60)."
    )
  )

  ini({
    # Final PK model, Table 2 'PK model with covariate' column (NONMEM 6,
    # ADVAN3 TRANS4, FOCE-I). Parenthesised values in Table 2 are RSEs.
    lcl <- log(2.65)
    label("Metabolic clearance at SBP = 135 mmHg, CL (L/min)") # Table 2: theta Cl = 2.65 (RSE 15%)
    lvc <- log(27.7)
    label("Central volume of distribution, VC (L)") # Table 2: theta VC = 27.7 (RSE 27%)
    lq <- log(2.73)
    label("Inter-compartmental clearance, CLD (L/min)") # Table 2: theta CLD = 2.73 (RSE 28%)
    lvp <- log(801)
    label("Peripheral volume of distribution, VT (L)") # Table 2: theta VT = 801 (RSE 33%)

    # Linear-deviation effect of systolic blood pressure on clearance,
    # CL = 2.65 * (1 - 0.00714 * (SBP - 135)) (Table 2, Results, Abstract).
    # The minus sign lives in the model equation, as printed.
    e_sbp_cl <- 0.00714
    label("Fractional decrease in clearance per mmHg of SBP above 135 (1/mmHg)") # Table 2: 0.00714 (RSE 66%)

    # Direct-effect inhibitory Emax model on BIS, Eq. 2:
    #   BIS = BIS0 - Emax * Cp / (EC50 + Cp)
    # BIS0 and Emax were both fixed to 100 (Methods, 'PK/PD model of
    # propofol'; Table 3 'theta EMAX 100 (Fixed)'). No Hill coefficient.
    le0 <- fixed(log(100))
    label("Baseline (fully awake) BIS, BIS0 (BIS units)") # Methods: BIS0 fixed to 100
    limax <- fixed(log(100))
    label("Maximum decrease in BIS, Emax (BIS units)") # Table 3: theta EMAX = 100 (Fixed)
    lec50 <- log(1.98)
    label("Plasma propofol concentration giving half the maximal BIS decrease, EC50 (mg/L)") # Table 3: theta EC50 = 1.98 (RSE 20%)

    # Log-normal IIV, Eq. 3. Tables 2-3 print the omega rows in percent and
    # the text calls them 'IIV ... 48.6%' and 'IIV of 45.1%': coefficients of
    # variation, CV = sqrt(omega^2), so omega^2 = CV^2.
    etalcl ~ 0.236196 # Table 2: 'omega Cl2 (%)' = 48.6 (RSE 50%); 0.486^2
    etalec50 ~ 0.203401 # Table 3: 'omega EC502 [%]' = 45.1 (RSE 37%); 0.451^2

    # Proportional residual errors, Eqs. 4-5. Printed in percent.
    propSd <- 0.418
    label("Proportional residual error on plasma propofol (fraction)") # Table 2: 'sigma prop,Cp2 (%)' = 41.8 (RSE 17%)
    propSd_BIS <- 0.191
    label("Proportional residual error on BIS (fraction)") # Table 3: residual row (printed 'sigma prop,Cp2') = 19.1 (RSE 28%)
  })

  model({
    # Individual parameters. Only clearance and EC50 carry IIV ('We failed
    # to estimate IIV for any other parameters', Results).
    cl <- exp(lcl + etalcl) * (1 - e_sbp_cl * (SBP - 135))
    vc <- exp(lvc)
    q <- exp(lq)
    vp <- exp(lvp)

    e0 <- exp(le0)
    imax <- exp(limax)
    ec50 <- exp(lec50 + etalec50)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Eq. 1, written on amounts: infusion and boluses enter central.
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # mg / L = mg/L
    Cc <- central / vc

    # Eq. 2: direct effect of the plasma concentration on BIS.
    BIS <- e0 - imax * Cc / (ec50 + Cc)

    Cc ~ prop(propSd)
    BIS ~ prop(propSd_BIS)
  })
}
