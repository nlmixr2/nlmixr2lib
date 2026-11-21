Dickinson_2026_dolutegravir <- function() {
  description <- "Two-compartment population PK model for oral dolutegravir with lagged first-order absorption in healthy volunteers dosed to steady state and followed through drug washout, with allometric weight scaling and a lower absorption rate constant in men; developed to compare dose-reduction strategies for fixed-dose tenofovir-lamivudine-dolutegravir"
  reference <- paste(
    "Dickinson L, Else L, Venter WDF, Boffito M, Brand R, Waitt C, Hill A,",
    "Khoo S. A model-based comparison of dose reduction strategies for",
    "fixed-dose dolutegravir-containing regimens. J Antimicrob Chemother.",
    "2026;81(1):dkaf445. doi:10.1093/jac/dkaf445.",
    "Structure adapted from Kawuma AN, Wasmann RE, Dooley KE, Boffito M,",
    "Maartens G, Denti P. Population pharmacokinetic model and alternative",
    "dosing regimens for dolutegravir coadministered with rifampicin.",
    "Antimicrob Agents Chemother. 2022;66(6):e00215-22.",
    "doi:10.1128/aac.00215-22; every parameter VALUE below is from the",
    "2026 paper's own Table 2.",
    sep = " "
  )
  vignette <- "Dickinson_2026_dolutegravir"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Allometric scaling included a priori on all four disposition parameters, exponents fixed at 0.75 for CL/F and Q/F and 1 for Vc/F and Vp/F (Methods 'Pharmacokinetic modelling'; Table 2 footnote b). The paper does not print the reference weight; 70 kg is assumed, as in the Kawuma model it adapted (the cohort median was 64.4 kg, Table 1).",
      source_name = "Weight"
    ),
    SEXF = list(
      description = "Biological sex indicator, 1 = female, 0 = male",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (male) in the canonical register; the SOURCE paper's reference category is FEMALE",
      notes = "Table 2 reports 'Male sex on ka' = 0.274, a multiplicative factor against the female reference (Results: 'Males had a 72.6% lower ka than females (typical values 0.71 versus 2.58 h-1)'). Encoded as ka * (SEXF + (1 - SEXF) * e_sexf_ka) so the printed typical ka of 2.58 /h is the female value and the printed factor is kept verbatim.",
      source_name = "Male sex"
    ),
    OCC = list(
      description = "Integer-valued occasion indicator for between-occasion-variability multiplexing",
      units = "(count)",
      type = "categorical",
      reference_category = NULL,
      notes = "The paper retains interoccasion variability on relative bioavailability, with occasions defined as in the Kawuma model it adapted ('with the pre-dose being treated as a separate occasion', Methods). One variance is reported and shared by every occasion, so occasion 1 carries the estimated variance and occasions 2-4 fix it to the same value (the NONMEM $OMEGA BLOCK(1) SAME idiom; see Kawuma_2023_dolutegravir.R). Extending the chain is a mechanical copy of the fix() lines. For single-occasion records pass OCC = 1.",
      source_name = "OCC"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "dolutegravir", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "dolutegravir", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "dolutegravir", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 73L,
    n_studies = 3L,
    n_observations = 825L,
    age_median = "33.3 years",
    age_range = "20.0-67.0 years",
    weight_median = "64.4 kg",
    weight_range = "41.5-97.5 kg",
    sex_female_pct = 60.3,
    race_ethnicity = c(White = 46.6, Black = 49.3, Asian = 4.1),
    disease_state = "healthy HIV-negative volunteers",
    dose_range = "dolutegravir 50 mg once daily to steady state (6-15 days), alone or with emtricitabine/tenofovir alafenamide or lamivudine/tenofovir disoproxil fumarate; two studies followed concentrations for 10 and 14 days after the last dose",
    regions = "United States (Pittsburgh), United Kingdom (London) and Uganda (Kampala)",
    notes = "Pooled from three studies (Table 1): Pittsburgh NCT04302896 (n = 30, 14-day washout after day-15 final dose), SSAT061 NCT02219217 (n = 17, 9-day washout after 10 days of dosing) and DolACT NCT02242799 (n = 26, dolutegravir-alone period only, 24-h steady-state profile). Doses were taken with food. Concentrations below the assay LLQ (10 ng/mL Pittsburgh and DolACT; 0.75 ng/mL SSAT061) were fitted with the M3 method."
  )

  ini({
    # =========================================================================
    # Structural parameters. Apparent (oral) values for the typical 70 kg
    # (assumed reference) woman, from Table 2 'Estimate (95% CI)'; the
    # confidence intervals are bootstrap (n = 1000, Table 2 footnote a).
    # =========================================================================
    lcl <- log(0.787)
    label("Apparent clearance CL/F (L/h)")
    # Table 2 row 'CL/F (L/h)' = 0.787 (0.751-0.824)

    lvc <- log(15.0)
    label("Apparent central volume of distribution Vc/F (L)")
    # Table 2 row 'Vc/F (L)' = 15.0 (14.2-15.2)

    lq <- log(0.0093)
    label("Apparent inter-compartmental clearance Q/F (L/h)")
    # Table 2 row 'Q/F (L/h)' = 0.0093 (0.0072-0.011)

    lvp <- log(0.671)
    label("Apparent peripheral volume of distribution Vp/F (L)")
    # Table 2 row 'Vp/F (L)' = 0.671 (0.586-0.756)

    lka <- log(2.58)
    label("Absorption rate constant ka in women (1/h)")
    # Table 2 row 'ka (h-1)' = 2.58 (0.246-4.92); Results 'typical values 0.71
    # versus 2.58 h-1' (male versus female), so 2.58 is the female value.

    ltlag <- log(0.836)
    label("Absorption lag time, fed dosing (h)")
    # Table 2 row 'Lag time (h)' = 0.836 (0.776-0.897). All three studies dosed
    # with food (Methods 'Analysis population'), so this is a fed lag time.

    lfdepot <- fixed(log(1))
    label("Relative bioavailability (unitless)")
    # Table 2 row 'Relative bioavailability (F)' = 1 fixed

    # =========================================================================
    # Covariate effects.
    # =========================================================================
    e_wt_cl <- fixed(0.75)
    label("Allometric exponent on CL/F and Q/F (unitless)")
    # Table 2 row 'Weight on CL/F & Q/F' = 0.75 fixed (footnote b: allometric scaling)

    e_wt_vc <- fixed(1)
    label("Allometric exponent on Vc/F and Vp/F (unitless)")
    # Table 2 row 'Weight on Vc/F & Vp/F' = 1 fixed (footnote b: allometric scaling)

    e_sexf_ka <- 0.274
    label("Ratio of male to female ka, applied to (1 - SEXF) (unitless)")
    # Table 2 row 'Male sex on ka' = 0.274 (0.122-0.426); Results 'Males had a
    # 72.6% lower ka than females (typical values 0.71 versus 2.58 h-1)'.

    # =========================================================================
    # Between-subject variability (Table 2 'Random effects', %). The paper does
    # not state how the percentages were derived from the variances, so the
    # exponential-model default omega^2 = log(1 + CV^2) is used.
    # =========================================================================
    etalcl ~ log(1 + 0.215^2) # Table 2 row 'IIV CL/F (%)' = 21.5 (17.6-24.8)
    etalka ~ log(1 + 1.36^2) # Table 2 row 'IIV ka (%)' = 136 (111-157)

    # =========================================================================
    # Between-occasion variability on relative bioavailability. One variance is
    # reported and shared across occasions; occasion 1 carries the estimate and
    # occasions 2-4 fix it to the same value.
    # =========================================================================
    etaiov_fdepot_1 ~ log(1 + 0.303^2) # Table 2 row 'IOV F (%)' = 30.3 (24.4-35.2)
    etaiov_fdepot_2 ~ fix(log(1 + 0.303^2)) # SAME-equivalent: equal to the occasion-1 variance
    etaiov_fdepot_3 ~ fix(log(1 + 0.303^2)) # SAME-equivalent: equal to the occasion-1 variance
    etaiov_fdepot_4 ~ fix(log(1 + 0.303^2)) # SAME-equivalent: equal to the occasion-1 variance

    # =========================================================================
    # Residual unexplained variability (Table 2 'Residual variability'):
    # combined proportional plus additive on the linear mg/L scale, the
    # additive part restricted to 20% of the assay LLQ (Methods).
    # =========================================================================
    propSd <- 0.223
    label("Proportional residual error (fraction)")
    # Table 2 row 'Proportional (%)' = 22.3 (19.8-24.7)

    addSd <- 0.0013
    label("Additive residual error (mg/L)")
    # Table 2 row 'Additive (mg/L)' = 0.0013 (0.00058-0.0019)
  })

  model({
    # 1. Occasion multiplexing of the between-occasion eta on bioavailability.
    #    For single-occasion records pass OCC = 1.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    iov_fdepot <- oc1 * etaiov_fdepot_1 + oc2 * etaiov_fdepot_2 +
      oc3 * etaiov_fdepot_3 + oc4 * etaiov_fdepot_4

    # 2. Individual parameters. Allometric scaling with a 70 kg reference
    #    (assumed); men absorb with 0.274 x the female ka.
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    vc <- exp(lvc) * (WT / 70)^e_wt_vc
    q <- exp(lq) * (WT / 70)^e_wt_cl
    vp <- exp(lvp) * (WT / 70)^e_wt_vc
    ka <- exp(lka + etalka) * (SEXF + (1 - SEXF) * e_sexf_ka)
    tlag <- exp(ltlag)
    fdepot <- exp(lfdepot + iov_fdepot)

    # 3. Micro-constants.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 4. Two-compartment disposition with lagged first-order absorption from an
    #    oral depot (Figure S1).
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central -
      k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    f(depot) <- fdepot
    alag(depot) <- tlag

    # 5. Observation: dose in mg over volume in L gives mg/L, the unit of the
    #    additive error; the EC90 of 320 ng/mL is 0.32 mg/L.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
