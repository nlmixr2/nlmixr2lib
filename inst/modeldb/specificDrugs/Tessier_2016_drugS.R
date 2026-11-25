Tessier_2016_drugS <- function() {
  description <- paste(
    "Two-compartment population PK model of the anonymised oral 'drug S'",
    "(Institut de Recherches Servier) used by Tessier 2016 to simulate phase I",
    "and phase II pharmacogenetic studies. Each dose is absorbed by two",
    "parallel routes whose split and total depend on the dose: a lagged",
    "zero-order input straight into the central compartment (Tlag1 0.401 h,",
    "duration Tk0 1.59 h) carries (1 - FRAC) of the bioavailable dose, and a",
    "late first-order input (Tlag2 22.7 h, ka 0.203 /h) carries FRAC, with",
    "FRAC = 0.45 * DOSE / (18.6 + DOSE). Bioavailability is 1 up to 20 dose",
    "units and above that falls as F = 1 - 0.8 * (DOSE - 20) / (41.7 + DOSE - 20).",
    "Linear elimination (CL 94.9), V1 1520, Q 147, V2 2130; log-normal IIV on",
    "F, Tk0, Tlag1, Q, V2 and CL; 20% proportional residual error. The",
    "parameters are the estimates from the 78-subject phase I data of drug S;",
    "the drug's identity and its dose and volume units are not disclosed."
  )
  reference <- paste(
    "Tessier A, Bertrand J, Chenel M, Comets E. Combined Analysis of Phase I",
    "and Phase II Data to Enhance the Power of Pharmacogenetic Tests.",
    "CPT Pharmacometrics Syst Pharmacol. 2016;5(3):123.",
    "doi:10.1002/psp4.12054.",
    "Parameter values from Table 1 and the dose equations printed beneath it;",
    "model structure and the scale of the simulated concentrations from the",
    "authors' simulation R code (Supporting Information, Data S1) and",
    "Supplementary Figure S4.",
    sep = " "
  )
  vignette <- "Tessier_2016_drugS"

  units <- list(
    time = "h",
    dosing = "units (the paper's undisclosed dose unit of drug S)",
    concentration = "milliunits/L (1000 x dose units per volume unit; ng/mL if the dose unit is mg and volumes are in L)"
  )

  covariateData <- list(
    DOSE_DRUGS_UNITS = list(
      description = "Drug S amount per administration, in the paper's dose units",
      units = "units",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Sets both dose-dependent terms: the first-order fraction",
        "ffo = emax_ffo * DOSE / (d50_ffo + DOSE) and the bioavailability,",
        "1 up to 20 units and 1 - imax_fdepot * (DOSE - 20) / (d50_fdepot + DOSE - 20)",
        "above. Must equal the amt of the two dose records of each",
        "administration. The authors' Monolix dataset carried it as the",
        "regressor DMG = DOSE / 1000, where DOSE was the dose record amount in",
        "milliunits."
      ),
      source_name = "DMG"
    )
  )

  compartmentData <- list(
    depot = list(
      analyte = "drug S (anonymised)",
      units = "units",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "drug S (anonymised)",
      units = "units",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "drug S (anonymised)",
      units = "units",
      specimen = "plasma",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 78L,
    n_studies = 3L,
    disease_state = "Phase I clinical study participants (health status not reported)",
    dose_range = paste(
      "Single oral doses of 5, 10, 20, 50, 100, 200, 400 or 800 units",
      "(6, 6, 24, 12, 12, 6, 6 and 6 subjects)"
    ),
    sampling = "16 samples per subject at 0.5-192 h after the dose",
    genotyping = "All subjects genotyped at baseline for 176 SNPs on a PK-focused DNA microarray",
    notes = paste(
      "Tessier 2016 Methods 'Simulation study' and 'Simulated designs':",
      "drug S (IRIS, Institut de Recherches Servier) PK data from 78 subjects",
      "in three phase I studies; the population PK analysis itself is",
      "described in Tessier 2013 PAGE abstract 2836. No demographics are",
      "reported. The paper's simulated phase II studies (306 subjects, 20,",
      "50 or 100 units once daily at steady state) are designs, not data."
    )
  )

  ini({
    # Disposition -- Table 1 'mu' column; volume and clearance units are not printed
    lcl <- log(94.9); label("Linear elimination clearance CL (volume units/h)") # Table 1 'CL' mu = 94.9
    lvc <- log(1520); label("Central volume V1 (volume units)") # Table 1 'V1' mu = 1520
    lq <- log(147); label("Intercompartmental clearance Q (volume units/h)") # Table 1 'Q' mu = 147
    lvp <- log(2130); label("Peripheral volume V2 (volume units)") # Table 1 'V2' mu = 2130

    # Zero-order route into central
    ld1 <- log(1.59); label("Zero-order absorption duration Tk0 (h)") # Table 1 'Tk0' mu = 1.59
    ltlag <- log(0.401); label("Lag time of the zero-order absorption Tlag1 (h)") # Table 1 'Tlag1' mu = 0.401

    # First-order route through depot
    lka <- log(0.203); label("First-order absorption rate constant Ka (1/h)") # Table 1 'Ka' mu = 0.203
    ltlag2 <- log(22.7); label("Lag time of the first-order absorption Tlag2 (h)") # Table 1 'Tlag2' mu = 22.7

    # Dose-dependent first-order fraction FRAC = EmaxFRAC * DOSE / (D50FRAC + DOSE)
    lemax_ffo <- log(0.45); label("Maximum fraction of the bioavailable dose absorbed first-order EmaxFRAC (fraction)") # Table 1 'FRAC EmaxFRAC' mu = 0.45
    ld50_ffo <- log(18.6); label("Dose giving half the maximum first-order fraction D50FRAC (dose units)") # Table 1 'FRAC D50FRAC' mu = 18.6

    # Dose-dependent bioavailability: F = 1 up to 20 units, F = 1 - ImaxF * (DOSE - 20) / (D50F + DOSE - 20) above
    lfdepot <- fixed(log(1)); label("Bioavailability at doses up to 20 units (fraction)") # Table 1 footnote 'For doses<20 units F=1'; Data S1 'Fpop<-ifelse(DOSE/1000<=20,1,...)'
    limax_fdepot <- log(0.8); label("Maximum fractional decrease of bioavailability ImaxF (fraction)") # Table 1 'F ImaxF' mu = 0.8
    ld50_fdepot <- log(41.7); label("Dose increment above 20 units giving half the maximum decrease of bioavailability D50F (dose units)") # Table 1 'F D50F' mu = 41.7

    # IIV: Table 1 'omega (%)' is the SD of the log-normal random effect
    # (Data S1 code: rnorm(sd = OMEGA)), so the variance is its square.
    etalfdepot ~ 0.108241 # Table 1 'F' omega = 32.9%; 0.329^2
    etald1 ~ 0.099856 # Table 1 'Tk0' omega = 31.6%; 0.316^2
    etaltlag ~ 0.123201 # Table 1 'Tlag1' omega = 35.1%; 0.351^2
    etalq ~ 0.808201 # Table 1 'Q' omega = 89.9%; 0.899^2
    etalvp ~ 0.195364 # Table 1 'V2' omega = 44.2%; 0.442^2
    etalcl ~ 0.063001 # Table 1 'CL' omega = 25.1%; 0.251^2 (Methods: omega2 CL = 0.06)

    propSd <- 0.2; label("Proportional residual error (fraction)") # Table 1 'sigma slope (%)' = 20; Data S1 SIG_inter = 0, SIG_slope = 0.2
  })

  model({
    # Dose-dependent first-order fraction and bioavailability (Table 1 footnote equations)
    emax_ffo <- exp(lemax_ffo)
    d50_ffo <- exp(ld50_ffo)
    ffo <- emax_ffo * DOSE_DRUGS_UNITS / (d50_ffo + DOSE_DRUGS_UNITS)

    imax_fdepot <- exp(limax_fdepot)
    d50_fdepot <- exp(ld50_fdepot)
    # Bioavailability does not fall below this dose (Table 1 footnote: 'For doses < 20 units F = 1')
    dose_fdepot_thr <- 20
    dose_above <- max(DOSE_DRUGS_UNITS - dose_fdepot_thr, 0)
    fdepot <- exp(lfdepot + etalfdepot) * (1 - imax_fdepot * dose_above / (d50_fdepot + dose_above))

    # Individual PK parameters; V1, Ka, Tlag2, EmaxFRAC and D50FRAC carry no IIV (Table 1)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc)
    q <- exp(lq + etalq)
    vp <- exp(lvp + etalvp)
    ka <- exp(lka)
    d1 <- exp(ld1 + etald1)
    tlag <- exp(ltlag + etaltlag)
    tlag2 <- exp(ltlag2)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d / dt(depot) <- -ka * depot
    d / dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d / dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Each administration is TWO dose records with the same amt: one to depot
    # (first-order route) and one to central with rate = -2 (zero-order route;
    # rxode2 applies dur() only to rate = -2 records).
    f(depot) <- fdepot * ffo
    alag(depot) <- tlag2
    f(central) <- fdepot * (1 - ffo)
    alag(central) <- tlag
    dur(central) <- d1

    # Data S1 simulates concentrations as (1000 x dose units) / volume
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
