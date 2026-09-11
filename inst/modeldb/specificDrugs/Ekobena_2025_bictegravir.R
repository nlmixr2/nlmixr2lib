Ekobena_2025_bictegravir <- function() {
  description <- "One-compartment population PK model for bictegravir in real-world people with HIV followed by therapeutic drug monitoring, with first-order absorption fixed at 0.64 1/h, apparent oral clearance and volume, a power body-weight effect and a median-centred exponential-linear age effect on apparent clearance, between-subject variability on clearance only, and a proportional residual error"
  reference <- paste(
    "Ekobena P, Briki M, Dao K, Marzolini C, Andre P, Buclin T, Cavassini M,",
    "Guidi M, Thoueille P; Swiss HIV Cohort Study. Population pharmacokinetics",
    "of bictegravir in real-world people with HIV. J Antimicrob Chemother.",
    "2025;80(11):2782-2789. doi:10.1093/jac/dkaf297.",
    sep = " "
  )
  vignette <- "Ekobena_2025_bictegravir"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Enters apparent clearance as a power (the paper's 'log-linear') term centred on BWRef = 70 kg, the median reference value named in the Results text immediately below the final covariate equation. Note that 70 kg is a rounded reference, NOT the cohort median: Table 1 gives a cohort median body weight of 74 kg (range 37-135 kg). The exponent 0.35 is estimated, not fixed to an allometric 0.75. Self-check against the Results sentence: exp(0.35 * log(100/70)) = 1.133, i.e. the '13% increased CL' quoted for a 100 kg individual aged 51 years. Baseline weight; the paper draws it from the routine Swiss HIV Cohort Study visits (every 3-6 months) and does not describe it as time-varying in the model.",
      source_name        = "BW"
    ),
    AGE = list(
      description        = "Subject age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Enters apparent clearance as the paper's 'linear model' form, which is linear in the FRACTIONAL deviation from the cohort median and sits inside the same exponential as the weight term: exp(beta_Age * (AGE - AgeM) / AgeM) with AgeM = 51 years, the median age of the study population (Table 1: 51 years, range 19-81). The division by AgeM is load-bearing and easy to drop -- without it the coefficient -0.20 applied to a raw 29-year deviation would give exp(-5.8), an 0.3% clearance. Self-check against the Results sentence: exp(-0.20 * (80 - 51) / 51) = 0.8925, i.e. the '11% decreased CL' quoted for an 80-year-old weighing 70 kg. Note the SIGN: older age lowers clearance and therefore RAISES exposure, which is the direction the paper's Table 3 simulations show (65-80 year olds have ~20% higher Ctrough than 20-65 year olds).",
      source_name        = "age"
    )
  )

  covariatesDataExcluded <- list(
    BMI = list(
      description = "Body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Screened and statistically significant on CL in the univariate step (Ekobena 2025 Results, dBICc = -3.2), but discarded at the multivariate / backward-deletion step in favour of body weight (dBICc = -21.0), so it does not appear in the final model. Cohort median 24.5 kg/m^2, range 12.5-44.5 (Table 1). No coefficient is reported for it, so it cannot be encoded."
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/L",
      type        = "continuous",
      notes       = "Considered a priori because bictegravir is >99% protein bound, but not retained. The paper states why in the Discussion limitations: albumin measured within +/- 30 days of a bictegravir sample was available for only 34% of drug levels, and the available values were tightly clustered (median 45 g/L, range 28-55), 'limiting the statistical power to detect a potential effect'. No coefficient is reported."
    )
  )

  compartmentData <- list(
    depot   = list(analyte = "bictegravir", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "bictegravir", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 572L,
    n_studies      = 1L,
    n_observations = 708L,
    age_range      = "19-81 years",
    age_median     = "51 years",
    weight_range   = "37-135 kg",
    weight_median  = "74 kg",
    sex_female_pct = 31,
    race_ethnicity = c(White = 65, Black = 26, `Hispanic American` = 4, Asian = 4, `Other/missing` = 1),
    disease_state  = "HIV infection; unselected real-world outpatients enrolled in the Swiss HIV Cohort Study and receiving bictegravir as part of routine care",
    dose_range     = "bictegravir 50 mg once daily orally, available only as the fixed-dose bictegravir/emtricitabine/tenofovir alafenamide combination",
    regions        = "Switzerland (therapeutic drug monitoring performed at Lausanne University Hospital, July 2019 to July 2024)",
    notes          = paste(
      "Baseline characteristics from Ekobena 2025 Table 1 (n = 572). This is a",
      "sparse therapeutic-drug-monitoring dataset, not a rich PK study: 815",
      "concentrations were collected and 107 excluded (42 undetectable, most",
      "likely non-adherence; 40 missing time information; 25 with uncertainties",
      "about reported information), leaving 708 concentrations. A median of ONE",
      "sample was collected per individual (range 1-4) and sampling times after",
      "the last dose spanned 0.25-32 h. Steady state was assumed for every",
      "sample, on the strength of the long follow-up and a reference half-life",
      "of 17.3 h. Assay LLOQ 15 ng/mL by UHPLC-MS/MS. Additional Table 1",
      "characteristics not used as covariates: height median 173 cm (147-198),",
      "BMI median 24.5 kg/m^2 (12.5-44.5), albumin median 45 g/L (28-55), CD4",
      ">= 500 cells/mm^3 in 33% (54% missing) and plasma HIV RNA < 50 copies/mL",
      "in 50% (38% missing). Three participants were pregnant. Co-medication",
      "with metabolism inhibitors was recorded in 15 (2%) and with inducers in",
      "4 (<1%) participants; darunavir and ritonavir were the most frequent",
      "concomitant drugs. No relationship between bictegravir concentrations and",
      "inhibitor or inducer co-administration was detected, which the authors",
      "attribute to lack of power rather than absence of an interaction",
      "(Discussion), so this model must NOT be used to predict bictegravir",
      "exposure under strong CYP3A4/UGT1A1 inhibition or induction. Sex,",
      "gender and ethnicity were extracted and tested but no categorical",
      "covariate was retained. No external validation of the model was",
      "performed (Discussion limitations).",
      sep = " "
    )
  )

  ini({
    # =========================================================================
    # Structural parameters. CL and V are APPARENT (oral) values -- the study
    # has no intravenous reference arm and bioavailability was neither
    # estimated nor reported, so F is subsumed into both. All values are the
    # 'Final model Estimate (RSE, %)' column of Ekobena 2025 Table 2; the
    # bootstrap medians and 95% CIs quoted alongside come from the same table
    # (2000 replicates) and agree with the point estimates to within 3%, as
    # the Results state.
    # =========================================================================
    lka <- fixed(log(0.64))
    label("First-order absorption rate constant ka (1/h)")
    # Ekobena 2025 Table 2 row 'ka (h-1)' = '0.64 FIX' (bootstrap column also
    # '0.64 FIX'). Fixed rather than estimated, and the paper says why in
    # Methods: 0.64 was 'the value initially estimated in our study with a
    # relative standard error (RSE) of 38%, due to the limited availability of
    # data collected shortly after drug intake. This decision was upheld
    # following a bootstrap analysis, which revealed instability in the
    # parameter's estimation.' Consequence for users: the absorption phase of
    # this model is NOT informed by the data and should not be relied on for
    # Cmax or Tmax predictions -- the Discussion says as much ('The few
    # observations collected in the early phase after drug intake prevented to
    # describe precisely bictegravir absorption'). No absorption lag time was
    # included, unlike the manufacturer's unpublished analysis.

    lcl <- log(0.46)
    label("Apparent oral clearance CL/F for a 70 kg, 51-year-old individual (L/h)")
    # Ekobena 2025 Table 2 row 'CL (L/h)' = 0.46, RSE 1% (bootstrap median
    # 0.46, 95% CI 0.44-0.47). This is the typical value at the covariate
    # reference point BWRef = 70 kg and AgeM = 51 years.

    lvc <- log(10.9)
    label("Apparent volume of distribution V/F (L)")
    # Ekobena 2025 Table 2 row 'V (L)' = 10.9, RSE 5% (bootstrap median 11.2,
    # 95% CI 10.1-12.4). No covariate was retained on V and no between-subject
    # variability was supportable on it (see the etalcl block below).
    # Cross-check on the pair: kel = 0.46 / 10.9 = 0.04220 1/h gives a terminal
    # half-life of 16.4 h and, with ka = 0.64, a Tmax of 4.55 h -- matching the
    # 't1/2 of 16.3 h' and 'Tmax of 4.5 h' the Results quote for the base
    # model.

    # =========================================================================
    # Covariate effects on apparent clearance. The final covariate model is
    # printed as a single equation in Ekobena 2025 Results:
    #
    #   CLi = CL * exp( beta_Bodyweight * log(BW / BWRef)
    #                 + beta_Age * (Age - AgeM) / AgeM
    #                 + eta_i )
    #
    # with BWRef = 70 kg and AgeM = 51 years named in the sentence directly
    # below it. Both coefficients live inside ONE exponential together with
    # the IIV term; they are the two generic forms defined in Methods (the
    # 'log-linear model' for BW and the 'linear model' for age), and both of
    # those generic forms are themselves exponential -- see the Methods
    # equations Pi = Ppop * exp(beta_cov * log(cov_i / covRef) + eta_i) and
    # Pi = Ppop * exp(beta_cov * (cov_i - covRef) / covRef + eta_i).
    # =========================================================================
    e_wt_cl <- 0.35
    label("Body-weight power exponent on apparent clearance (unitless)")
    # Ekobena 2025 Table 2 row 'beta Bodyweight' = 0.35, RSE 19% (bootstrap
    # median 0.36, 95% CI 0.22-0.49). Estimated, NOT fixed at an allometric
    # 0.75. Because the term is written exp(beta * log(BW / 70)), it is
    # algebraically the power form (BW / 70)^0.35.

    e_age_cl <- -0.20
    label("Age coefficient on apparent clearance, per unit fractional deviation from 51 years (unitless)")
    # Ekobena 2025 Table 2 row 'beta Age' = -0.20, RSE 28% (bootstrap median
    # -0.19, 95% CI -0.31 to -0.08). The coefficient multiplies the
    # FRACTIONAL deviation (AGE - 51) / 51, not the raw deviation in years --
    # see the covariateData$AGE note. Negative, so clearance falls and
    # exposure rises with age.

    # =========================================================================
    # Between-subject variability. Carried on apparent clearance ONLY:
    # Ekobena 2025 Results, 'IIV associated with CL, adequately described
    # bictegravir PK. Assignment of IIV on V did not improve data description
    # (dBICc = +8.4)'. The Discussion adds that the absence of variability on
    # V is consistent with the manufacturer's unpublished analysis (11%) and
    # reflects the limited role of the distribution phase.
    #
    # Table 2 heads the row 'omega CL (CV%)' = 22, so the printed 22 is a
    # coefficient of variation and is converted to the log-normal variance
    # exactly: omega^2 = log(1 + 0.22^2) = 0.04726519. If the authors instead
    # meant omega itself (the sqrt(omega^2) approximation that many Monolix
    # tables print under a '(CV%)' header), the variance would be 0.0484 --
    # a 1.2% difference in the standard deviation, immaterial to every
    # simulation gate in the vignette. The exact reading is used here because
    # the header says CV%, and it matches the convention used by the sibling
    # models in this register (e.g. Punyawudho_2025_dolutegravir.R). PK
    # parameters were assumed log-normally distributed (Methods).
    # =========================================================================
    etalcl ~ 0.04726519
    # Ekobena 2025 Table 2 row 'omega CL (CV%)' = 22, RSE 6% (bootstrap median
    # 22, 95% CI 18-25); omega^2 = log(1 + 0.22^2) = 0.04726519.

    # =========================================================================
    # Residual unexplained variability. Proportional only: Ekobena 2025
    # Results, 'The RUV was best captured by a proportional error model'
    # (Methods tested proportional, additive and combined). For a proportional
    # error model the tabulated CV% is the standard deviation of the
    # proportional term directly, so no log-normal conversion applies here.
    # =========================================================================
    propSd <- 0.27
    label("Proportional residual error (fraction)")
    # Ekobena 2025 Table 2 row 'sigma prop (CV%)' = 27, RSE 4% (bootstrap
    # median 28, 95% CI 25-30).
  })

  model({
    # 1. Individual parameters. The clearance line reproduces the single final
    #    covariate equation printed in Ekobena 2025 Results, written in the
    #    paper's own exp(beta * log(BW / BWRef) + beta * (Age - AgeM) / AgeM)
    #    form rather than the algebraically identical power form, so that the
    #    source trace is a literal one. BWRef = 70 kg and AgeM = 51 years.
    cl <- exp(lcl + etalcl) *
      exp(e_wt_cl * log(WT / 70) + e_age_cl * (AGE - 51) / 51)

    # No between-subject variability on V (dBICc = +8.4 when it was added) and
    # no covariate was retained on it, so vc is the typical value for every
    # individual.
    vc <- exp(lvc)

    # ka is fixed (Table 2, '0.64 FIX') and carries no IIV; no lag time was
    # estimated.
    ka <- exp(lka)

    # 2. Micro-constant.
    kel <- cl / vc

    # 3. One-compartment disposition with first-order absorption from an oral
    #    depot and first-order elimination from the central compartment
    #    (Ekobena 2025 Results: 'A one-compartment model with linear
    #    elimination, parameterized in terms CL, V, ka fixed to 0.64 h-1').
    #    A second compartment was tested and rejected (dBICc = +18.6).
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # 4. Bioavailability is NOT applied: CL and V above are already apparent
    #    (CL/F, V/F) values, so introducing an f(depot) term would double-count
    #    the unknown bioavailability.

    # 5. Observation and residual error. A 50 mg dose with vc in L gives Cc in
    #    mg/L; the paper reports concentrations in ng/mL, so multiply Cc by
    #    1000 to compare against its 2781 ng/mL typical Ctrough, its 760 ng/mL
    #    minimum concentration recommended for efficacy and its 15 ng/mL assay
    #    LLOQ.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
