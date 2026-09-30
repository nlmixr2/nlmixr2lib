Gupta_2021_brigatinib <- function() {
  description <- paste(
    "Three-compartment population PK model for oral brigatinib in healthy",
    "volunteers and patients with cancer (mostly ALK-positive non-small cell",
    "lung cancer) (Gupta 2021). Absorption follows the Savic transit",
    "compartment model (a non-integer number of transit compartments with",
    "inter-individual variability, and a mean transit time) feeding the",
    "central compartment directly, with no separate first-order absorption",
    "step. Elimination is linear from the central compartment. Apparent",
    "clearance carries a power effect of serum albumin centred on 38 g/L;",
    "inter-individual variability is on CL/F and V1/F (correlated), on the",
    "first peripheral volume, and on both transit parameters, with a",
    "proportional residual error."
  )
  reference <- paste(
    "Gupta N, Wang X, Offman E, Prohn M, Narasimhan N, Kerstein D,",
    "Hanley MJ, Venkatakrishnan K. (2021).",
    "Population Pharmacokinetics of Brigatinib in Healthy Volunteers and",
    "Patients With Cancer.",
    "Clin Pharmacokinet 60:235-247.",
    "doi:10.1007/s40262-020-00929-4.",
    sep = " "
  )
  vignette <- "Gupta_2021_brigatinib"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    ALB = list(
      description = "Serum albumin concentration",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power effect on CL/F, (ALB / 38)^0.661 (Table 3; Figure 1a).",
        "Table 2 and the Figure 4 legend label the median of 38 as g/dL,",
        "which is physiologically impossible; Section 2.4 gives the",
        "simulated-population albumin as 36 (20-47) g/L and the successor",
        "ALTA-1L analysis (Gupta 2022, doi:10.1111/cts.13231) tabulates",
        "albumin in g/L with the same model, so 38 is g/L. Baseline value."
      ),
      source_name = "Albumin"
    )
  )

  compartmentData <- list(
    depot = list(
      analyte = "brigatinib",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte = "brigatinib",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "brigatinib",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    peripheral2 = list(
      analyte = "brigatinib",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 442L,
    n_studies = 5L,
    age_range = "18-83 years",
    age_median = "52 years",
    weight_range = "41-172 kg",
    weight_median = "73 kg",
    sex_female_pct = 48.6,
    race_ethnicity = c(White = 68.8, Asian = 23.3, Black = 5.9, Other = 2.0),
    disease_state = paste(
      "105 healthy volunteers and 337 patients with cancer (201 with",
      "ALK-positive NSCLC from the phase II ALTA trial; 136 from a phase I/II",
      "dose-escalation study in advanced malignancies, mostly NSCLC)"
    ),
    dose_range = paste(
      "Single oral doses of 90, 120 or 180 mg (healthy volunteers);",
      "30-300 mg qd or 60-120 mg bid (phase I/II); 90 mg qd or 180 mg qd",
      "after a 7-day 90 mg qd lead-in (ALTA)"
    ),
    albumin = "38 (20-56) g/L, median (range)",
    egfr = "85.6 (32.7-277.5) mL/min/1.73 m^2, median (range)",
    notes = paste(
      "Demographics from Gupta 2021 Table 2; studies from Table 1.",
      "6086 PK samples, of which 247 below the LLOQ, 80 with suspected",
      "recording errors and 28 with |CWRES| > 4 were excluded.",
      "Healthy-volunteer food-effect and drug-interaction studies",
      "contributed only fasted, brigatinib-alone periods."
    )
  )

  ini({
    lcl <- log(10.6); label("Apparent clearance CL/F at albumin 38 g/L (L/h)") # Table 3, CL/F = 10.6 (RSE 2.7%)
    lvc <- log(207); label("Apparent central volume V1/F (L)") # Table 3, V1/F = 207 (RSE 3.8%)
    lq <- log(12.6); label("Apparent intercompartmental clearance to the first peripheral compartment Q1/F (L/h)") # Table 3, Q1/F = 12.6 (RSE 9.2%)
    lvp <- log(114); label("Apparent first peripheral volume V2/F (L)") # Table 3, V2/F = 114 (RSE 11.2%)
    lq2 <- log(2.7); label("Apparent intercompartmental clearance to the second peripheral compartment Q2/F (L/h)") # Table 3, Q2/F = 2.7 (RSE 18.2%)
    lvp2 <- log(78.5); label("Apparent second peripheral volume V3/F (L)") # Table 3, V3/F = 78.5 (RSE 8.6%)
    lntr <- fixed(log(2.35)); label("Number of Savic transit compartments (unitless)") # Table 3, 2.35; footnote c: fixed at the population estimate in the final model
    lmtt <- log(0.9); label("Mean transit time (h)") # Table 3, mean transit time = 0.9 (RSE 3.6%)

    e_alb_cl <- 0.661; label("Power exponent of albumin (ALB/38 g/L) on CL/F (unitless)") # Table 3, 'x (albumin/38)^0.661' (RSE 10.7%); Figure 1a

    # IIV: Table 3 reports sqrt(omega^2) x 100% (footnote a), so omega^2 = (CV/100)^2.
    # The CL/F - V1/F covariance 0.228 is reported on the variance scale (correlation 0.846).
    etalcl + etalvc ~ c(0.234256, 0.228, 0.309136) # Table 3, IIV CL/F 48.4%, V1/F 55.6%, Covariance (CL/F, V1/F) = 0.228
    etalvp ~ 0.906304 # Table 3, IIV V2/F 95.2%
    etalntr ~ 1.0609 # Table 3, IIV number of transit compartments 103%
    etalmtt ~ 0.3481 # Table 3, IIV mean transit time 59%

    propSd <- 0.269; label("Proportional residual error (fraction)") # Table 3, proportional error = 0.269 (SD scale; see vignette)
  })
  model({
    # Individual parameters (Section 2.2 Equation 1 exponential IIV; Figure 1a).
    cl <- exp(lcl + etalcl) * (ALB / 38)^e_alb_cl
    vc <- exp(lvc + etalvc)
    q <- exp(lq)
    vp <- exp(lvp + etalvp)
    q2 <- exp(lq2)
    vp2 <- exp(lvp2)
    ntr <- exp(lntr + etalntr)
    mtt <- exp(lmtt + etalmtt)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k13 <- q2 / vc
    k31 <- q2 / vp2

    # Savic transit absorption (Figure 1a): the dose passes through the
    # transit chain at rate ktr straight into the central compartment; there
    # is no separate first-order absorption compartment and no ka in Table 3.
    # The chain is collapsed into its closed form: the rate leaving it is the
    # Savic gamma density ktr * Dose * (ktr t)^ntr * exp(-ktr t) / Gamma(ntr + 1)
    # with ktr = (ntr + 1) / mtt, which allows a non-integer ntr. The depot
    # holds the whole transit chain, so it receives the dose as a normal bolus
    # and drains by exactly that density (depot + central + peripherals then
    # conserve mass). As in the NONMEM Savic implementation, the density is
    # driven by the most recent dose only; with a mean transit time near 1 h
    # the carry-over between doses given 12 or 24 h apart is negligible.
    tdos <- tad(depot)
    ktr <- (ntr + 1) / mtt
    ktt <- ktr * tdos
    trin <- 0
    if (ktt > 0) {
      trin <- exp(log(podo(depot) * ktr) - lgamma(ntr + 1) + ntr * log(ktt) - ktt)
    }

    d/dt(depot) <- -trin
    d/dt(central) <- trin - kel * central - k12 * central + k21 * peripheral1 - k13 * central + k31 * peripheral2
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(peripheral2) <- k13 * central - k31 * peripheral2

    # Dose in mg, volume in L: mg/L x 1000 = ng/mL.
    Cc <- central / vc * 1000
    Cc ~ prop(propSd)
  })
}
