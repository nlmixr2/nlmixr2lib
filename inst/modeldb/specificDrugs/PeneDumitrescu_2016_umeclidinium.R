PeneDumitrescu_2016_umeclidinium <- function() {
  description <- "Two-compartment population PK model for umeclidinium after a single topical application of [14C]umeclidinium 1.85% (w/w) solution to the occluded axilla of healthy men, pooled with 30-min intravenous infusion data, with sequential zero-order (direct to central) then lagged first-order dermal absorption, a very low absolute dermal bioavailability and linear elimination (Pene Dumitrescu 2016)"
  reference <- paste(
    "Pene Dumitrescu T, Santos LL, Hughes SC, Pereira AI, Young GC, Hussey E,",
    "Charlton P, Baptiste-Brown S, Stuart JS, Vincent V, van Marle SP, Schmith VD.",
    "A Novel Method for Studying the Pharmacokinetics of [14C]Umeclidinium After",
    "Application to the Axilla or Palm of Healthy Male Subjects. Clin Transl Sci.",
    "2016;9(4):183-191. doi:10.1111/cts.12406",
    sep = " "
  )
  vignette <- "PeneDumitrescu_2016_umeclidinium"
  units <- list(time = "h", dosing = "ug", concentration = "pg/mL")
  # Unit note: doses are in ug and volumes in L, so central / vc is in ug/L
  # (= ng/mL); the factor 1000 in the observation line converts to pg/mL, the
  # unit of every concentration and of the AMS assay LOQ (0.348 pg/mL) in
  # Pene Dumitrescu 2016 (Methods 'Analytical methods', Table 2, Figure 3).

  covariateData <- list(
    ROUTE_IV = list(
      description = "Indicator for intravenous administration",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (topical application to the occluded axilla)",
      notes = paste(
        "Per-dose-record indicator (1 = 30-min intravenous infusion, 0 = topical dose to the",
        "occluded axilla). The pooled dataset combined the six occluded-axilla subjects with",
        "nine subjects given 65 ug UMEC as a 30-min i.v. infusion in a prior study",
        "(NCT01110018). Figure 2 routes the i.v. dose and the zero-order part of the dermal",
        "dose (F*F1*dose) into the same central compartment, but only the dermal dose carries",
        "the absolute bioavailability F. rxode2 applies one bioavailability per compartment,",
        "so ROUTE_IV selects between them: f(central) collapses to 1 for the i.v. infusion and",
        "to F*F1 for the zero-order dermal input. An i.v. record is an infusion into",
        "cmt = 'central' with an explicit rate (which makes rxode2 ignore dur(central)). A",
        "dermal administration needs TWO records at the same time with the same amt: one",
        "into cmt = 'depot' and one into cmt = 'central' carrying rate = -2 so the modelled",
        "zero-order duration is applied. ROUTE_IV = 0 on both dermal records.",
        sep = " "
      ),
      source_name = "administration route (IV vs dermal)"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "umeclidinium", units = "ug", specimen = "administration site", verified = TRUE),
    central = list(analyte = "umeclidinium", units = "ug", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "umeclidinium", units = "ug", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 15L,
    n_studies = 2L,
    age_range = "30-55 years (occluded-axilla cohort, Table 1)",
    age_mean = "47.8 years (occluded-axilla cohort, Table 1)",
    weight_mean = "82.9 kg (occluded-axilla cohort, Table 1)",
    sex_female_pct = 0,
    race_ethnicity = c(White = 100),
    disease_state = "healthy volunteers",
    dose_range = "single topical dose of 165 mg of a 1.85% (w/w) [14C]umeclidinium solution (about 3.06 mg umeclidinium, ~666 kBq) over 40 cm2 of one occluded axilla, removed after 8 h; plus a single 65 ug i.v. infusion over 30 min (prior study)",
    regions = "The Netherlands (dermal study NCT01934153, GSK 117157)",
    notes = paste(
      "The model was fit to the 6 healthy men of the occluded-axilla (OA) cohort of the",
      "single-dose [14C]UMEC study (BMI 18-27 kg/m2 by protocol; mean BMI 26.0, height",
      "178.5 cm; Table 1) pooled with 9 subjects given a 30-min i.v. infusion of 65 ug",
      "UMEC in a prior study (NCT01110018, demographics not reported). The unoccluded",
      "axilla and occluded palm cohorts (6 men each) were not modelled because of too many",
      "nonquantifiable concentrations. Plasma [14C]UMEC was measured by accelerator mass",
      "spectrometry (LOQ 0.348 pg/mL). No covariate effects were investigated.",
      sep = " "
    )
  )

  ini({
    # Structural parameters: Pene Dumitrescu 2016 Table 3, column 'Model
    # parameter estimates' (Phoenix NLME 1.2, FOCE-I).
    lcl <- log(55.0); label("Clearance CL (L/h)") # Table 3 'CL (L/h)' 55.0 (bootstrap median 47.5, RSE 23.0%)
    lvc <- log(7.56); label("Central volume of distribution V1 (L)") # Table 3 'V1 (L)' 7.56 (bootstrap median 7.61, RSE 15.4%)
    lq <- log(40.9); label("Intercompartmental clearance CL2 (L/h)") # Table 3 'CL2 (L/h)' 40.9 (bootstrap median 41.6, RSE 30.8%)
    lvp <- log(194.7); label("Peripheral volume of distribution V2 (L)") # Table 3 'V2 (L)' 194.7 (bootstrap median 255.1, RSE 57.7%)
    lka <- log(0.0934); label("First-order dermal absorption rate constant Ka (1/h)") # Table 3 'Ka (1/h)' 0.0934 (bootstrap median 0.0910, RSE 45.4%)
    # Tlag is both the lag of the first-order input and the duration of the
    # zero-order input: the absorption processes are 'sequential' (Abstract,
    # Results) and Table 3 has no separate duration parameter.
    ltlag <- log(10.6); label("Lag time of the first-order input = duration of the zero-order input, Tlag (h)") # Table 3 'Tlag (h)' 10.6 (bootstrap median 4.71, RSE 49.1%)
    lfdepot <- log(0.00275); label("Absolute bioavailability of the dermal dose to the occluded axilla, F (fraction)") # Table 3 'F' 0.00275 (bootstrap median 0.00217, RSE 53.4%)
    # Table 3 reports F1, the fraction of the bioavailable dermal dose absorbed
    # by the zero-order process. The canonical logitffo is the complementary
    # first-order fraction: logitffo = logit(1 - F1).
    logitffo <- qlogis(1 - 0.312); label("Logit of the fraction of the bioavailable dermal dose absorbed first-order, 1 - F1 (unitless)") # Table 3 'F1' 0.312 (bootstrap median 0.171, RSE 76.6%)

    # Interindividual variability: Table 3 reports CV% of exponential
    # (log-normal) random effects; omega^2 = log(1 + CV^2).
    etalcl ~ 0.67915 # Table 3 'CL' IIV CV 98.6%
    etalvc ~ 0.22475 # Table 3 'V1' IIV CV 50.2%
    etalq ~ 0.73314 # Table 3 'CL2' IIV CV 104%
    etalka ~ 1.07572 # Table 3 'Ka' IIV CV 139%
    etaltlag ~ 1.09464 # Table 3 'Tlag' IIV CV 141%

    # Residual error: proportional (Methods 'Population PK modeling').
    propSd <- 0.406; label("Proportional residual error (fraction)") # Table 3 'Residual error, CV%' 40.6%
  })

  model({
    # Individual parameters (Methods: Theta_i = Theta_TV * exp(eta_i)).
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q <- exp(lq + etalq)
    vp <- exp(lvp)
    ka <- exp(lka + etalka)
    tlag <- exp(ltlag + etaltlag)
    fdepot <- exp(lfdepot)
    ffo <- expit(logitffo)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Figure 2: the dermal dose splits into F*F1 entering central zero-order and
    # F*(1-F1) entering the absorption pool, absorbed first-order after Tlag.
    # The i.v. infusion enters central with no bioavailability term.
    f(central) <- ROUTE_IV + (1 - ROUTE_IV) * fdepot * (1 - ffo)
    dur(central) <- tlag
    f(depot) <- fdepot * ffo
    alag(depot) <- tlag

    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
