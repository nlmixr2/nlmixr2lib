Perlstein_2022_risperidone_tv46000 <- function() {
  description <- paste(
    "Two-compartment population PK model for the risperidone total active",
    "moiety (TAM = risperidone + 9-hydroxyrisperidone * 410/426) following",
    "subcutaneous TV-46000, a copolymer-based long-acting risperidone",
    "suspension, in 97 adults with schizophrenia or schizoaffective disorder",
    "from one phase 1 single- and multiple-dose study (Perlstein 2022).",
    "Absorption is the paper's convolution-based prescribed input: a double",
    "Weibull release profile that splits the dose between a first release",
    "process (shape below 1, front-loaded) and a second, sigmoidal release",
    "process, implemented as two parallel depot compartments each emptying",
    "with its own Weibull hazard. Allometric body weight on CL/F (exponent",
    "0.75) and V/F (exponent 1) fixed at the 70 kg reference. Exponential",
    "inter-individual variability on every structural parameter and a",
    "combined additive + proportional residual error. Carries the published",
    "Emax dopamine D2 receptor occupancy layer (Kd 10.1 ng/mL) the authors",
    "used in their simulations as an algebraic observable. This is the",
    "phase 1 TAM model that selected the phase 3 doses; the later pooled",
    "parent-metabolite TV-46000 model is",
    "modellib('Perlstein_2025_risperidone_tv46000').",
    sep = " "
  )
  reference <- paste(
    "Perlstein I, Merenlender Wagner A, Gomeni R, Lamson M, Harary E,",
    "Spiegelstein O, Kalmanczhelyi A, Tiver R, Loupe P, Levi M, Elgart A (2022).",
    "Population Pharmacokinetic Modeling and Simulation of TV-46000: A",
    "Long-Acting Injectable Formulation of Risperidone.",
    "Clin Pharmacol Drug Dev 11(7):865-877. doi:10.1002/cpdd.1078.",
    sep = " "
  )
  vignette <- "Perlstein_2022_risperidone_tv46000"
  units <- list(time = "week", dosing = "mg", concentration = "ng/mL")
  # Declared explicitly because buildModelDb()'s fallback heuristic recognises
  # only the literal names "depot" and "central"; every TV-46000 injection must
  # be recorded on BOTH depots (the release fraction is applied as a
  # bioavailability split across the two parallel Weibull depots), and central
  # is never dosed by this subcutaneous-only model.
  dosing <- c("depot", "depot2")

  compartmentData <- list(
    depot = list(
      analyte = "risperidone active moiety (risperidone + 9-OH-risperidone)",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    depot2 = list(
      analyte = "risperidone active moiety (risperidone + 9-OH-risperidone)",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "risperidone active moiety (risperidone + 9-OH-risperidone)",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "risperidone active moiety (risperidone + 9-OH-risperidone)",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Allometric scaling on CL/F and V/F with the exponents FIXED at 0.75 and 1",
        "(Methods, Population Pharmacokinetic Analysis, Equation 5: 'the allometric",
        "coefficients of 0.75 for clearance and 1 for volume were fixed in the model",
        "according to Anderson and Holford'), centred on WT = 70 kg. The study mean",
        "weight was 87 kg (Results, Population Pharmacokinetic Model Development)."
      ),
      source_name = "WT"
    )
  )

  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      notes = "Formally tested on CL/F (Table 1 run 4, change in -2LL +12.69) and on V/F (run 5, -2.132) but NOT retained; the final model is the reference run 1. Cohort mean 44.4 years (SD 8.8)."
    ),
    INJSITE_ARM = list(
      description = "Upper-arm (vs abdomen) subcutaneous injection site indicator (1 = upper arm, 0 = abdomen)",
      units = "(binary)",
      type = "binary",
      notes = "Formally tested on td (Table 1 run 2, +20.85) and on V/F (run 3, +10.85) but NOT retained. Only cohort 8 (225 mg single dose) was injected in the upper arm."
    ),
    BMI = list(
      description = "Body mass index",
      units = "kg/m^2",
      type = "continuous",
      notes = "Screened graphically against the empirical Bayes estimates (Results) but not formally tested; no coefficient reported. Cohort mean 28.6 kg/m^2 (SD 4.4)."
    ),
    CRCL_BASE = list(
      description = "Baseline creatinine clearance (Cockcroft-Gault, not BSA-normalized)",
      units = "mL/min",
      type = "continuous",
      notes = "Screened graphically (Methods Equation 6) but not formally tested; no coefficient reported. Cohort mean 115.76 mL/min (SD 21.9)."
    ),
    SEXF = list(
      description = "Sex (1 = female, 0 = male)",
      units = "(binary)",
      type = "binary",
      notes = "Screened graphically but not formally tested; the Discussion calls the evaluation inconclusive because the cohort was 80% male."
    ),
    RACE_BLACK = list(
      description = "Black or African American race indicator (1 = yes, 0 = no)",
      units = "(binary)",
      type = "binary",
      notes = "Screened graphically but not formally tested; the Discussion calls the evaluation inconclusive because the cohort was 95.9% Black or African American."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 97,
    n_studies = 1,
    age_mean = "44.4 years (SD 8.8)",
    weight_mean = "87 kg",
    bmi_mean = "28.6 kg/m^2 (SD 4.4)",
    crcl_mean = "115.76 mL/min (SD 21.9)",
    sex_female_pct = 20,
    race_ethnicity = c(Black = 95.9),
    disease_state = "Clinically stable adults with a diagnosis of schizophrenia or schizoaffective disorder, not on an antipsychotic other than oral risperidone.",
    dose_range = paste(
      "Single subcutaneous TV-46000 injections of 50, 75, 100, 150 or 225 mg (abdomen; cohorts 1-5),",
      "225 mg in the upper arm (cohort 8), or three once-monthly injections of 75 or 150 mg",
      "(abdomen; cohorts 6-7). Each cohort first received 7 days of oral risperidone 2-6 mg/day",
      "followed by a 7-day washout."
    ),
    regions = "United States",
    notes = paste(
      "Phase 1 study TV-46000-SAD-10055, open-label, adaptive, single- and multiple-dose.",
      "99 patients enrolled; two cohort-4 patients were excluded for self-medicating with oral",
      "risperidone, leaving 97 in the analysis (Methods, General Considerations for Data",
      "Management). Concentrations below the LLQ were treated as missing. Estimation by SAEM",
      "in NONMEM 7.4."
    )
  )

  ini({
    # ---------------------------------------------------------------------
    # Double Weibull release input function (Methods Equation 3, Figure 2):
    #   r(t) = ff*exp(-(t/td)^ss) + (1 - ff)*exp(-(t/td1)^ss1)
    # r(t) is the fraction of the dose still UNRELEASED (the Results check it:
    # 50% released by 2.6 weeks, 79% by week 4, 93% by week 8), so the input
    # into the central compartment is Dose*(-dr/dt). Written with the
    # registered Weibull stems ra / gam1 and ra2 / gam2:
    #   r(t) = frel*exp(-(ra*t)^gam1) + (1 - frel)*exp(-(ra2*t)^gam2)
    # so ra = 1/td and ra2 = 1/td1. The paper's 'ss' cannot be used as a name
    # here in any case: ss is the rxode2 steady-state event flag.
    # ---------------------------------------------------------------------
    lra <- log(1 / 2.69); label("Weibull rate-scaling parameter of the first release process (ra = 1/td, 1/week)") # Table 2: td = 2.69 week (RSE 16.2%); ra = 1/2.69 = 0.372 1/week
    lgam1 <- log(0.616); label("Weibull shape (sigmoidicity) of the first release process (unitless)") # Table 2: ss = 0.616 (RSE 3.3%)
    lra2 <- log(1 / 3.32); label("Weibull rate-scaling parameter of the second release process (ra2 = 1/td1, 1/week)") # Table 2: td1 = 3.32 week (RSE 9.5%); ra2 = 1/3.32 = 0.301 1/week
    lgam2 <- log(3.66); label("Weibull shape (sigmoidicity) of the second release process (unitless)") # Table 2: ss1 = 3.66 (RSE 8%)

    # Fraction of the dose released by the first process, held on the logit
    # scale so it stays inside (0, 1) for every eta draw. Table 2 prints it as
    # a fraction ('ff,%' = 0.511; the Results text's 'ff = 0.55%' is a typo).
    # logit(0.511) = 0.04400.
    logitfrel <- log(0.511 / (1 - 0.511)); label("Logit of the fraction of the dose released by the first process (unitless)") # Table 2: ff = 0.511 (RSE 5.7%)

    # Disposition (Table 2). Time unit is the week, as printed.
    lcl <- log(354); label("Apparent clearance CL/F at the 70 kg reference weight (L/week)") # Table 2: CL/F = 354 L/wk (RSE 3.8%)
    lvc <- log(374); label("Apparent central volume of distribution V/F at the 70 kg reference weight (L)") # Table 2: V/F = 374 L (RSE 9.4%)
    lk12 <- log(1.48); label("First-order transfer rate constant central to peripheral1 (1/week)") # Table 2: k23 = 1.48 1/week (RSE 30.4%)
    lk21 <- log(2.34); label("First-order transfer rate constant peripheral1 to central (1/week)") # Table 2: k32 = 2.34 1/week (RSE 32.2%)

    # Allometric exponents fixed per Methods Equation 5 (Anderson and Holford).
    e_wt_cl <- fixed(0.75); label("Body-weight allometric exponent on CL/F (unitless)") # Methods Eq. 5: 0.75 fixed
    e_wt_vc <- fixed(1); label("Body-weight allometric exponent on V/F (unitless)") # Methods Eq. 5: 1 fixed

    # Published Emax D2 receptor occupancy layer (Methods, PK/D2RO Model,
    # Equation 7: RO = ROmax * Cp / (Kd + Cp)). Neither value was estimated
    # in this study, so both are fixed.
    emax <- fixed(100); label("Maximal attainable dopamine D2 receptor occupancy ROmax (%)") # Methods PK/D2RO Model: ROmax fixed to 100%
    lec50 <- fixed(log(10.1)); label("TAM concentration giving 50% D2 receptor occupancy, Kd (ng/mL)") # Methods PK/D2RO Model: Kd = 10.1 ng/mL

    # ---------------------------------------------------------------------
    # Inter-individual variability, Table 2 Random effect block. The rows are
    # NONMEM OMEGA VARIANCES of log-normal etas (Methods: 'It was assumed that
    # the IIV of the model parameters was log-normally distributed'). The SE
    # column is printed as SE*100 (e.g. td 0.163 * 45.3% = 0.0738 -> '7.38').
    # var(log ra) = var(-log td) = var(log td), so the omega carries onto the
    # reciprocal parameterisation unchanged.
    # ---------------------------------------------------------------------
    etalra ~ 0.163 # Table 2 Random effect td = 0.163 (RSE 45.3%, shrinkage 54.6%)
    etalra2 ~ 0.733 # Table 2 Random effect td1 = 0.733 (RSE 14.1%, shrinkage 6.49%)
    etalgam1 ~ 0.0247 # Table 2 Random effect ss = 0.0247 (RSE 59.9%, shrinkage 24.7%)
    etalgam2 ~ 0.185 # Table 2 Random effect ss1 = 0.185 (RSE 30.4%, shrinkage 20.5%)
    # The paper puts a log-normal eta on ff itself: 0.511*exp(eta) with
    # variance 0.115 exceeds 1 for about 2.4% of subjects, which would make the
    # second fraction negative. Re-expressed on the logit scale by the delta
    # method, var(logit ff) = var(log ff) / (1 - ff)^2 = 0.115 / 0.489^2 = 0.481.
    etalogitfrel ~ 0.481 # Table 2 Random effect ff = 0.115 (RSE 26%, shrinkage 3.67%) on the log scale; delta-method logit-scale variance, see Errata
    etalcl ~ 0.117 # Table 2 Random effect CL/F = 0.117 (RSE 27.9%, shrinkage 27.8%)
    etalvc ~ 0.147 # Table 2 Random effect V/F = 0.147 (RSE 48.1%, shrinkage 19.0%)
    etalk12 ~ 1.41 # Table 2 Random effect k23 = 1.41 (RSE 25%, shrinkage 29.1%)
    etalk21 ~ 1.9 # Table 2 Random effect k32 = 1.9 (RSE 71.6%, shrinkage 23.4%)

    # Residual error, Table 2 Residual effect block. Read as standard deviations
    # of a single-epsilon combined error model: the table reports one epsilon
    # shrinkage for the pair, the layout NONMEM gives when the two components
    # are THETAs scaling one EPS.
    propSd <- 0.152; label("Proportional residual error (fraction)") # Table 2: Error prop. = 0.152 (RSE 7.8%)
    addSd <- 1.58; label("Additive residual error (ng/mL)") # Table 2: Error_add. = 1.58 (RSE 18.3%)
  })

  model({
    # 1. Release input function parameters.
    ra <- exp(lra + etalra)
    gam1 <- exp(lgam1 + etalgam1)
    ra2 <- exp(lra2 + etalra2)
    gam2 <- exp(lgam2 + etalgam2)
    # The fixed effect and the eta are collected on their own line so the
    # term stays mu-referenced.
    logitfrel_ind <- logitfrel + etalogitfrel
    frel <- expit(logitfrel_ind)

    # 2. Disposition, allometrically scaled to the 70 kg reference weight
    #    (Methods Equation 5).
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    k12 <- exp(lk12 + etalk12)
    k21 <- exp(lk21 + etalk21)
    kel <- cl / vc

    # 3. Weibull release hazards. Each depot empties with the hazard of its
    #    own Weibull survival function S(t) = exp(-(rate*t)^shape), which is
    #    h(t) = shape*rate*(rate*t)^(shape - 1); a depot that empties this way
    #    releases exactly Dose*(-dS/dt), the paper's convolution input. Time
    #    after dose is floored at a small positive number because
    #    (rate*t)^(shape - 1) is not finite at t = 0 for the first process's
    #    shape below 1. tad0() returns 0 rather than NA before the first dose,
    #    so records ahead of the first injection solve cleanly. A new injection
    #    restarts the hazard clock for any amount still in the depot -- see the
    #    vignette for the size of that multiple-dose approximation.
    #
    #    The time argument is also capped where the cumulative hazard
    #    (rate*t)^shape reaches 40, i.e. once the depot holds less than
    #    exp(-40) = 4e-18 of its dose. Beyond that point the hazard of a
    #    shape > 1 process keeps growing as a power of time although the depot
    #    is empty; with the large second-process IIV a subject can draw a
    #    shape near 15, the hazard reaches ~1e10 per week within one dosing
    #    interval and the solver fails. Freezing the hazard there changes no
    #    released amount by more than exp(-40) of the dose.
    tr1 <- min(max(tad0(depot), 1e-6), 40^(1 / gam1) / ra)
    tr2 <- min(max(tad0(depot2), 1e-6), 40^(1 / gam2) / ra2)
    h1 <- gam1 * ra * (ra * tr1)^(gam1 - 1)
    h2 <- gam2 * ra2 * (ra2 * tr2)^(gam2 - 1)

    # 4. ODE system (Methods Equation 4: A1 = depot, A2 = central,
    #    A3 = peripheral1; the paper's k23 / k32 are k12 / k21 here).
    d/dt(depot) <- -h1 * depot
    d/dt(depot2) <- -h2 * depot2
    d/dt(central) <- h1 * depot + h2 * depot2 - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # 5. Dose split across the two parallel release processes. Each injection
    #    is recorded with the same amt on both depots.
    f(depot) <- frel
    f(depot2) <- 1 - frel

    # 6. Observation: TAM in risperidone-equivalent ng/mL. Dose in mg over
    #    volume in L gives mg/L; the factor of 1000 converts to ng/mL.
    Cc <- 1000 * central / vc

    # D2 receptor occupancy (percent), a deterministic transform of the TAM
    # concentration with no residual error (Methods Equation 7).
    ec50 <- exp(lec50)
    D2RO <- emax * Cc / (ec50 + Cc)

    Cc ~ add(addSd) + prop(propSd)
  })
}
