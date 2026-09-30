Lim_2017_ethanol <- function() {
  description <- "Population PK model for orally ingested ethanol (alcohol) in 24 healthy Korean adult males (Lim 2017, the 'original alcohol PK model' used for a forensic Bayesian back-estimation). One-compartment model with first-order absorption and Michaelis-Menten elimination; body weight scales the volume and Vmax linearly (reference 70 kg). Concentrations are blood alcohol concentrations in percent weight/volume (g/dL) estimated by breath alcohol test."
  reference <- "Lim HS, Soung JH, Bae KS. Forensic science meets clinical pharmacology: pharmacokinetic model based estimation of alcohol concentration of a defendant as requested by a local prosecutor's office. Transl Clin Pharmacol. 2017;25(1):5-9. doi:10.12793/tcp.2017.25.1.5"
  vignette <- "Lim_2017_ethanol"
  # The paper prints the dose as '32 mg', V in 'L', Vmax in '%/hour' and Km in
  # '%'. Those labels are not mutually consistent. The numeric model maps the
  # dose to the concentration as Cc = amount / V, and it reproduces the
  # paper's own Table 3 profile with dose = 32 per subject (see vignette).
  # Reading the dose as grams of ethanol (4 cups x 8 g) and Cc as g/dL (percent
  # w/v) makes V dimensionally dL (typical 372 dL = 37.2 L, 0.53 L/kg) and Vmax
  # g/h (amount per hour), which are physiologic values for ethanol.
  units <- list(time = "h", dosing = "g", concentration = "g/dL")

  compartmentData <- list(
    depot = list(analyte = "ethanol", units = "g", specimen = "administration site", verified = TRUE),
    central = list(analyte = "ethanol", units = "g", specimen = "whole blood", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Linear scaling of V and Vmax, V = V(70) * (WT/70) and Vmax = Vmax(70) * (WT/70), per the Table 1 footnote.",
      source_name = "WT"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 24L,
    n_observations = 178L,
    n_studies = 1L,
    age_range = "adults; not reported",
    weight_range = "not reported; reference weight 70 kg",
    sex_female_pct = 0,
    race_ethnicity = "Korean",
    disease_state = "Healthy adult male volunteers (phase 1 clinical trial)",
    dose_range = "Single oral ethanol dose printed as '32 mg' (Table 1 title); numerically 32 dose units, interpreted as 32 g ethanol",
    regions = "Republic of Korea",
    notes = "The original alcohol PK model (OAPKM) of Lim 2017 Methods and Table 1: 178 blood alcohol concentrations (percent) estimated by breath alcohol test in 24 healthy Korean adult males. Structure follows Bruno et al. 1983 (linear absorption, saturable elimination). The paper's application is a MAP Bayesian estimation, over 1,000 bootstrap replicates, of one defendant's individual parameters (Table 2) and blood alcohol profile (Table 3, Figure 2)."
  )

  ini({
    # Structural parameters: Lim 2017 Table 1, 'Estimate' column (typical values for a 70 kg subject).
    lka <- log(5.62); label("Absorption rate constant ka (1/h)") # Table 1 row 'Ka, 1/hour' = 5.62 (95% CI 2.39 to 8.85)
    lvc <- log(372); label("Volume of distribution V at 70 kg (dL; printed as L)") # Table 1 row 'V, L' = 372.00 (342.60 to 401.40)
    lvmax <- log(72.4); label("Maximum elimination rate Vmax at 70 kg (g/h)") # Table 1 row 'Vmax, %/hour' = 72.40 (-68.52 to 213.32)
    lkm <- log(0.47); label("Michaelis-Menten constant Km (g/dL)") # Table 1 row 'Km, %' = 0.47 (-0.56 to 1.50)

    # IIV: Table 1 prints omega^2 (log-normal variance). The Ka and Km CV% entries
    # equal sqrt(exp(omega^2) - 1) and the bootstrap V entry (0.023, 15.4) does
    # too; the single-run V and Vmax CV% entries (35.7, 40.1) do not match any
    # transform of the printed variances and are treated as typographical.
    etalka ~ 1.420 # Table 1 row 'IIVKa (CV %)' = 1.420 (177.1)
    etalvc ~ 0.026 # Table 1 row 'IIVV (CV %)' = 0.026 (35.7)
    etalvmax ~ 1.090 # Table 1 row 'IIVVmax (CV %)' = 1.090 (40.1)
    etalkm ~ 0.416 # Table 1 row 'IIVKm (CV %)' = 0.416 (71.8)

    # Residual error: Table 1 footnote states epsilon is reported as a standard deviation.
    addSd <- 0.005; label("Additive residual SD (g/dL)") # Table 1 row 'epsilon (additive), percent' = 0.005 (0.004 to 0.005)
    propSd <- 0.041; label("Proportional residual SD (fraction)") # Table 1 row 'epsilon (proportional)' = 0.041 (0.006 to 0.076)
  })

  model({
    # Linear weight scaling of V and Vmax (Table 1 footnote).
    ka <- exp(lka + etalka)
    vc <- exp(lvc + etalvc) * (WT / 70)
    vmax <- exp(lvmax + etalvmax) * (WT / 70)
    km <- exp(lkm + etalkm)

    # One-compartment model, first-order absorption, Michaelis-Menten
    # elimination of amount (Methods; Bruno et al. 1983 structure).
    Cc <- central / vc

    d / dt(depot) <- -ka * depot
    d / dt(central) <- ka * depot - vmax * Cc / (km + Cc)

    Cc ~ add(addSd) + prop(propSd)
  })
}
