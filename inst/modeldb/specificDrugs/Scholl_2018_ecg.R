Scholl_2018_ecg <- function() {
  description <- "Two-compartment population PK model for epicatechin gallate (ECG) after oral green tea extract in healthy adults (Scholl 2018), with sequential zero-order dissolution into the depot and lagged first-order absorption, no covariates, and combined proportional + additive residual error."
  reference <- paste(
    "Scholl C, Lepper A, Lehr T, Hanke N, Schneider KL, Brockmoller J, Seufferlein T, Stingl JC. (2018).",
    "Population nutrikinetics of green tea extract.",
    "PLoS One 13(2):e0193074. doi:10.1371/journal.pone.0193074"
  )
  vignette <- "Scholl_2018_green_tea_catechins"
  units <- list(time = "h", dosing = "mg", concentration = "ug/L")

  compartmentData <- list(
    depot = list(analyte = "epicatechin gallate", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "epicatechin gallate", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "epicatechin gallate", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list()

  population <- list(
    species = "human",
    n_subjects = 84L,
    n_studies = 1L,
    age_range = "19-49 years",
    age_mean = "24.5 years (SD 3.9)",
    weight_mean = "68.0 kg (SD 11.7)",
    sex_female_pct = 64.3,
    race_ethnicity = "Western European ancestry (all participants)",
    disease_state = "Healthy volunteers (BMI 18-30 kg/m2; no acute or chronic disease; no regular medication except oral contraceptives).",
    dose_range = "Green tea extract capsules containing 35 mg ECG (standardised to 150 mg EGCG) taken orally twice daily with a meal for 4 days, then one capsule after an overnight fast on day 5; PK sampled on day 5 pre-dose and 0.5-9 h post-dose.",
    regions = "Germany (University of Ulm)",
    notes = "Scholl 2018 Table 1 (84 analysed of 89 enrolled; 30 male, 54 female). NCT01360320. The three catechin models (EGCG, EGC, ECG) were developed independently on the same subjects; no covariate was retained for ECG."
  )

  ini({
    # Structural parameters -- Scholl 2018 Table 3, ECG column (final
    # estimates). Structure from the S6 File (Modeling Supplement) ECG NONMEM
    # control stream: ADVAN6 with explicit $DES, zero-order input of duration
    # D1 into the depot (RATE = -2 in the data), first-order K12 from depot to
    # central, lag time and F1 on the depot.
    lka <- log(1.78); label("First-order absorption rate constant ka, depot to central (1/h)") # Table 3 ECG ka = 1.78 1/h (S6 File K12 = THETA(1))
    lcl <- log(196); label("Apparent clearance CL/F (L/h)") # Table 3 ECG CL/F = 196 L/h
    lvc <- log(395); label("Apparent central volume Vcentral/F (L)") # Table 3 ECG Vcentral/F = 395 L
    lq <- log(56.2); label("Apparent intercompartmental clearance Q/F (L/h)") # Table 3 ECG Q/F = 56.2 L/h
    lvp <- log(1300); label("Apparent peripheral volume Vperipheral/F (L)") # Table 3 ECG Vperipheral/F = 1300 L
    ld1 <- log(0.422); label("Zero-order dissolution duration D1 into the depot (h)") # Table 3 ECG D1 = 0.422 h
    ltlag <- log(0.319); label("Absorption lag time (h)") # Table 3 ECG lag time = 0.319 h
    lfdepot <- fixed(log(1)); label("Relative bioavailability F1 (unitless)") # Table 3 ECG F = 1.00 fixed; S6 File $THETA (1 FIX) ;F

    # IIV -- Table 3 ECG IIV rows (%CV), omega^2 = (CV/100)^2 (see the
    # vignette for the convention). The S6 File ECG stream estimates a full
    # $OMEGA BLOCK(5) over K12, CL, D1, F1 and ALAG, but the paper reports
    # only the variances; the unreported covariances are set to zero.
    etalka ~ 0.6084 # Table 3 ECG IIV ka 78 %CV -> 0.78^2
    etalcl ~ 0.1521 # Table 3 ECG IIV CL/F 39 %CV -> 0.39^2
    etald1 ~ 0.49 # Table 3 ECG IIV D1 70 %CV -> 0.70^2
    etalfdepot ~ 0.3249 # Table 3 ECG IIV F 57 %CV -> 0.57^2
    etaltlag ~ 0.1849 # Table 3 ECG IIV lag time 43 %CV -> 0.43^2

    # Residual error -- S6 File $ERROR Y = IPRED*(1+EPS(1)) + EPS(2).
    propSd <- 0.18; label("Proportional residual error (fraction)") # Table 3 ECG prop. residual variability 18 %
    addSd <- 1; label("Additive residual error (ug/L)") # Table 3 ECG add. residual variability 1 ug/L
  })

  model({
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc)
    q <- exp(lq)
    vp <- exp(lvp)
    d1 <- exp(ld1 + etald1)
    tlag <- exp(ltlag + etaltlag)
    fdepot <- exp(lfdepot + etalfdepot)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Doses enter the depot as a zero-order infusion of duration D1 (dose
    # records need rate = -2), delayed by the lag time.
    dur(depot) <- d1
    alag(depot) <- tlag
    f(depot) <- fdepot

    # Amounts in mg, volumes in L; S2 = V2/1000 gives ug/L.
    Cc <- central / vc * 1000
    Cc ~ prop(propSd) + add(addSd)
  })
}
