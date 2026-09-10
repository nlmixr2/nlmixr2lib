Prohn_2021_letermovir_hsct <- function() {
  description <- "Two-compartment steady-state population PK model for letermovir given orally or intravenously to allogeneic hematopoietic stem cell transplant (HSCT) recipients, with first-order absorption, an absorption lag time and linear elimination. Concomitant cyclosporine carries separate estimates of both clearance (3.38 vs 4.84 L/h) and oral bioavailability (0.849 vs 0.346); healthy phase I participants pooled into the same fit carry a bioavailability fixed to 1 and an absorption rate roughly 8-fold faster than HSCT recipients. Asian participants have a 39.1 percent lower peripheral volume. Between-subject variability is carried on clearance, bioavailability, peripheral volume and the absorption rate, with between-occasion variability on bioavailability."
  reference <- paste(
    "Prohn M, Viberg A, Zhang D, Dykstra K, Davis C, Macha S, Sabato P,",
    "de Alwis D, Iwamoto M, Fancourt C, Cho CR (2021). Population",
    "pharmacokinetics of letermovir following oral and intravenous",
    "administration in healthy participants and allogeneic hematopoietic",
    "cell transplantation recipients. CPT Pharmacometrics Syst Pharmacol",
    "10(3):255-267. doi:10.1002/psp4.12593.",
    "This file encodes the HSCT recipient (phase III) model of Table 3 and",
    "Figure 1b; the healthy participant (phase I) model of Table 2 and",
    "Figure 1a is a separate, independently fitted model and is encoded in",
    "Prohn_2021_letermovir_healthy.R.",
    sep = " "
  )
  vignette <- "Prohn_2021_letermovir"
  # Prohn works throughout in ng/mL (Figure 1b defines Cp = A_central / Vc *
  # 1000, and the additive residual error is tabulated in ng/mL), so the
  # observation is scaled to ng/mL rather than the mg/L that mg and L give
  # directly.
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  compartmentData <- list(
    depot       = list(analyte = "letermovir", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "letermovir", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "letermovir", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    CONMED_CSA = list(
      description        = "Concomitant cyclosporine (CsA) coadministration indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no concomitant cyclosporine)",
      notes              = paste(
        "Prohn 2021 estimated cyclosporine cotreatment as two SEPARATE",
        "parameters rather than as a fractional shift: CL is a distinct THETA",
        "with (4.84 L/h) and without (3.38 L/h) CsA, and oral bioavailability",
        "likewise (0.849 with, 0.346 without). Table 3 reports the four values",
        "independently. They are re-expressed here as a reference value plus a",
        "log-ratio covariate effect, which is algebraically identical and keeps",
        "the file in the registry's covariate-effect idiom; the trailing",
        "comments carry the four original Table 3 numbers. The Results text",
        "confirms the direction and size of the clearance term:",
        "'approximately 30% lower CL with CSA' (exp(-0.3591) = 0.698, a 30.2%",
        "reduction). Cyclosporine inhibits OATP1B1-mediated hepatic uptake of",
        "letermovir; the Discussion hypothesises that the same mechanism also",
        "reduces first-pass metabolism, which is why apparent bioavailability",
        "rises at the same time clearance falls. Per protocol, HSCT recipients",
        "on cyclosporine received 240 mg/day rather than 480 mg/day.",
        sep = " "
      ),
      source_name        = "CSA"
    ),
    TX_HCT = list(
      description        = "Allogeneic hematopoietic stem cell transplant (HSCT) recipient indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (healthy phase I participant)",
      notes              = paste(
        "1 = HSCT recipient (the phase IIb and phase III cohorts, 363 of the",
        "399 subjects); 0 = healthy participant contributing steady-state",
        "phase I data (36 subjects) that were pooled in to anchor the",
        "structural model. The indicator switches TWO parameters, both of",
        "which the Results text calls out explicitly ('Healthy participants",
        "had very high bioavailability and an almost 10-fold faster absorption",
        "rate than HSCT recipients'): bioavailability, fixed to 1 in healthy",
        "participants because models estimating it were numerically unstable",
        "and returned values close to 100%; and the absorption rate, 1.26 1/h",
        "in healthy participants against 0.150 1/h in HSCT recipients. No",
        "healthy participant in this dataset received cyclosporine, so the",
        "CONMED_CSA arm of the bioavailability expression is only ever",
        "exercised with TX_HCT = 1.",
        sep = " "
      ),
      source_name        = "HP (healthy participant flag, complemented)"
    ),
    RACE_ASIAN = list(
      description        = "Asian race indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (non-Asian: White, Black, Hispanic and Other)",
      notes              = paste(
        "The only covariate retained by the stepwise search for this model.",
        "Multiplicative on the peripheral volume: Vp is 0.609 times the",
        "non-Asian value in Asian participants (Table 3, 'Asian effect",
        "peripheral volume'), i.e. 39.1% lower. Results: 'the stepwise",
        "covariate search identified a lower V2 for Asian subjects'.",
        "Figure S7 contrasts predicted profiles for non-Asian (White, Black,",
        "Hispanic, Other) versus Asian (Japanese plus Asian participants from",
        "other countries) populations, which is the grouping this indicator",
        "encodes. Adding body weight on CL, either alongside or in place of",
        "the Asian effect on Vp, did not improve the model, so no allometric",
        "term is carried.",
        sep = " "
      ),
      source_name        = "ASIAN"
    ),
    OCC = list(
      description        = "Dosing-occasion index for between-occasion variability on bioavailability",
      units              = "(index)",
      type               = "categorical",
      reference_category = "no reference category; 0 or any value outside 1-8 zeroes every occasion indicator",
      notes              = paste(
        "Prohn 2021 carried interoccasion variability on bioavailability",
        "(Table 3, 'IOV, bioavailability'), added because repeated trough",
        "samples in the same participant spanned a 100-fold range: 'IOV was",
        "included on bioavailability to account for the 100-fold difference in",
        "minimum concentration (Ctrough)'. The source does not state a maximum",
        "occasion count; phase III participants contributed sparse pre-dose",
        "samples at weeks 2, 4, 6, 8, 10, 12 and 14, so eight occasion slots",
        "are provided here, matching the idiom already used for the same drug",
        "in Royston_2025_letermovir.R. Occasion 1 carries the estimated",
        "variance and occasions 2-8 repeat it as a fixed value, which is the",
        "equivalent of NONMEM's $OMEGA BLOCK(1) SAME. Setting OCC to 0 (or",
        "omitting per-occasion structure) collapses the model to IIV only,",
        "which is what Prohn's own AUCss simulations did: exposure was",
        "computed 'without interoccasion variability'.",
        sep = " "
      ),
      source_name        = "OCC"
    )
  )

  covariatesDataExcluded <- list(
    WT = list(
      description = "Baseline body weight",
      units       = "kg",
      type        = "continuous",
      notes       = paste(
        "Pre-specified and tested on CL and bioavailability, and on V1 and V2",
        "(Supplementary Information, 'HSCT recipient (phase III model)",
        "covariate analysis'), but not retained: 'adding the effect of body",
        "weight on CL without the Asian effect on V2 or replacing the Asian",
        "covariate on V2 with body weight did not improve the model'. Body",
        "weight IS retained in the companion healthy participant (phase I)",
        "model, on both CLmax and Vd. Median 75 kg, range 35-142 kg (Table 1).",
        sep = " "
      ),
      source_name = "WT"
    ),
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Pre-specified and tested on CL, bioavailability, V1 and V2; not retained. Median 51 years, range 18-75 (Table 1).",
      source_name = "AGE"
    ),
    CRCL = list(
      description = "Creatinine clearance",
      units       = "mL/min",
      type        = "continuous",
      notes       = "Pre-specified and tested on CL, bioavailability, V1 and V2; not retained. No point estimate is reported.",
      source_name = "CrCl"
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Pre-specified and tested on CL, bioavailability, V1 and V2; not retained. 190 of 399 subjects (48%) were female (Table 1).",
      source_name = "sex"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 399,
    n_studies      = 5,
    age_range      = "median 51 years, range 18-75 (Table 1)",
    weight_range   = "median 75 kg, range 35-142 (Table 1)",
    sex_female_pct = 48,
    race_ethnicity = "Reported only as the Asian / non-Asian contrast retained in the model. Figure S7 defines non-Asian as White, Black, Hispanic and Other, and Asian as Japanese plus Asian participants from other countries. Counts by race are not tabulated.",
    disease_state  = "Cytomegalovirus-seropositive allogeneic hematopoietic stem cell transplant recipients receiving letermovir prophylaxis (363 of 399 subjects, 91%), pooled with 36 healthy participants (9%) who contributed steady-state phase I data at the same dosing schedule to anchor the structural model.",
    dose_range     = "240-480 mg once daily (480 mg/day alone, 240 mg/day with concomitant cyclosporine), orally or as a 1-hour intravenous infusion",
    regions        = "Not reported.",
    notes          = paste(
      "Data sources (Table 1): one phase III trial NCT02137772 (n = 350), one",
      "phase IIb trial NCT01063829 (n = 13) and three phase I trials (n = 36).",
      "2888 concentration observations, 2566 (89%) after oral and 322 (11%)",
      "after intravenous dosing; all at steady state, defined as at least one",
      "week of dosing and less than 72 h after the last dose. 53 observations",
      "(1.8%) below the 1 ng/mL lower limit of quantification were excluded.",
      "74 phase III participants contributed rich profiles and the remaining",
      "275 contributed sparse pre-dose samples only, which is why the source",
      "estimates separate residual error magnitudes for the two sampling",
      "schemes. Fitted in NONMEM 7.3.",
      sep = " "
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Structural parameters -- Prohn 2021 Table 3, "Estimate" column.
    # Table 3 reports arithmetic values; they are carried here as logs so
    # the log-normal IIV terms add on the natural-log scale.
    # ------------------------------------------------------------------
    lcl <- log(4.84)
    label("Clearance without concomitant cyclosporine (L/h)")  # Table 3 'CL non-CSA treatment, L/h' 4.84 (95% CI 4.3-5.45)
    lvc <- log(19.7)
    label("Central volume of distribution (L)")  # Table 3 'Central volume, L' 19.7 (95% CI 17.6-22.1)
    lvp <- log(25.8)
    label("Peripheral volume of distribution in non-Asian participants (L)")  # Table 3 'Peripheral volume, L' 25.8 (95% CI 19.1-34.9)
    lq <- log(1.54)
    label("Intercompartmental clearance (L/h)")  # Table 3 'Intercompartment CL, L/h' 1.54 (95% CI 1.17-2.04)
    lka <- log(0.150)
    label("First-order absorption rate in HSCT recipients (1/h)")  # Table 3 'Absorption rate, 1/h' 0.150 (95% CI 0.104-0.215)
    ltlag <- log(0.674)
    label("Absorption lag time (h)")  # Table 3 'Absorption lag, h' 0.674 (95% CI 0.59-0.769)
    lfdepot <- log(0.346)
    label("Oral bioavailability in HSCT recipients without cyclosporine (fraction)")  # Table 3 'Bioavailability without CSA' 0.346 (95% CI 0.278-0.42)

    # Healthy-participant branch. Bioavailability was FIXED to 1: "models in
    # which the bioavailability for healthy participants was estimated were
    # numerically unstable and resulted in bioavailability estimates close to
    # 100%. Therefore, to stabilize the model, bioavailability was fixed to
    # 100%." Table 3 accordingly gives no confidence interval for this row.
    lfdepot_hp <- fixed(log(1.00))
    label("Oral bioavailability in healthy participants (fraction)")  # Table 3 'Bioavailability HP' 1.00, no CI reported
    lka_hp <- log(1.26)
    label("First-order absorption rate in healthy participants (1/h)")  # Table 3 'Absorption rate HP, 1/h' 1.26 (95% CI 0.933-1.71)

    # ------------------------------------------------------------------
    # Covariate effects. Both cyclosporine terms are log-ratios of two
    # independently estimated Table 3 values (see covariateData$CONMED_CSA).
    # ------------------------------------------------------------------
    e_csa_cl <- log(3.38 / 4.84)
    label("Log-effect of concomitant cyclosporine on CL (unitless)")  # Table 3 'CL CSA treatment' 3.38 / 'CL non-CSA treatment' 4.84; = -0.3591, exp = 0.698, the "approximately 30% lower CL with CSA" of the Results
    e_csa_fdepot <- log(0.849 / 0.346)
    label("Log-effect of concomitant cyclosporine on oral bioavailability (unitless)")  # Table 3 'Bioavailability with CSA' 0.849 / 'Bioavailability without CSA' 0.346; = 0.8977, exp = 2.454
    e_asian_vp <- log(0.609)
    label("Log-effect of Asian race on peripheral volume (unitless)")  # Table 3 'Asian effect peripheral volume' 0.609 (95% CI 0.53-0.7), i.e. 39.1% lower Vp

    # ------------------------------------------------------------------
    # Inter-individual variability. Table 3's IIV / IOV column reports the
    # NONMEM $OMEGA elements, i.e. VARIANCES on the log scale, and is carried
    # verbatim below. The column is unlabelled and unitless, so this was
    # settled against the source's own reported prediction intervals rather
    # than assumed:
    #
    #   Intravenous 480 mg/day. AUCss = Dose / CL, so only the CL variance
    #   enters. Reported median 100,000 ng*h/mL, 90% PI 65,300-148,000
    #   ("Exposure predictions"). ln(148000/100000)/1.645 = 0.238 and
    #   ln(100000/65300)/1.645 = 0.259, i.e. an observed log-scale SD of
    #   about 0.249. sqrt(0.0605) = 0.246 reproduces it; reading 0.0605 as an
    #   SD instead would predict a 90% PI of 90,500-110,500, which is a
    #   four-fold narrower log-width than reported.
    #
    #   Oral 480 mg/day. AUCss = F * Dose / CL, so the CL and F variances add.
    #   Reported median 34,400, 90% PI 16,900-73,700. sqrt(0.0605 + 0.137) =
    #   0.444 predicts 16,600-71,400; the SD reading predicts 26,900-44,000.
    #
    # Both arms agree to within 3%, and the intravenous arm isolates CL, so
    # the variance reading is not in doubt. Note that two secondary sources
    # have re-used these numbers as SDs -- Fromage 2025 (doi:10.1371/
    # journal.pone.0321180) squares them again in its mrgsolve [OMEGA] block,
    # and Royston 2025 (doi:10.1128/aac.00697-25) copies 0.719 through as
    # "0.72" under a footnote declaring its table to hold SDs -- so a
    # comparison against either of those papers will show materially less
    # between-subject spread than this file produces.
    #
    # There is no IIV on Vc: "the model included between-subject variability
    # in CL and V2 but not in V1", which the Discussion gives as the reason
    # Cmax after intravenous dosing is not well predicted.
    # ------------------------------------------------------------------
    etalcl ~ 0.0605
    label("IIV on clearance")  # Table 3 'IIV, CL' 0.0605 (95% CI 0.0403-0.0807); log-scale variance, CV 24.9%
    etalfdepot ~ 0.137
    label("IIV on bioavailability")  # Table 3 'IIV, bioavailability' 0.137 (95% CI 0.0714-0.203); CV 38.4%
    etalvp ~ 0.229
    label("IIV on peripheral volume")  # Table 3 'IIV, peripheral volume' 0.229 (95% CI 0.0643-0.393); CV 50.7%
    etalka ~ 0.719
    label("IIV on absorption rate")  # Table 3 'IIV, absorption rate' 0.719 (95% CI 0.136-1.3); CV 104%

    # Between-occasion variability on bioavailability. Occasion 1 carries the
    # estimated variance; occasions 2-8 repeat it as a fixed value, the
    # equivalent of NONMEM's $OMEGA BLOCK(1) SAME.
    etaiov_fdepot_1 ~ 0.197
    label("IOV on bioavailability, occasion 1")  # Table 3 'IOV, bioavailability' 0.197 (95% CI 0.136-0.259)
    etaiov_fdepot_2 ~ fixed(0.197)  # Table 3 'IOV, bioavailability'; shared variance
    etaiov_fdepot_3 ~ fixed(0.197)  # Table 3 'IOV, bioavailability'; shared variance
    etaiov_fdepot_4 ~ fixed(0.197)  # Table 3 'IOV, bioavailability'; shared variance
    etaiov_fdepot_5 ~ fixed(0.197)  # Table 3 'IOV, bioavailability'; shared variance
    etaiov_fdepot_6 ~ fixed(0.197)  # Table 3 'IOV, bioavailability'; shared variance
    etaiov_fdepot_7 ~ fixed(0.197)  # Table 3 'IOV, bioavailability'; shared variance
    etaiov_fdepot_8 ~ fixed(0.197)  # Table 3 'IOV, bioavailability'; shared variance

    # ------------------------------------------------------------------
    # Residual error. Table 3 tabulates FIVE residual terms, because the
    # source estimated separate magnitudes for healthy participants, for
    # HSCT recipients with rich sampling, and for HSCT recipients with sparse
    # (trough-only) sampling:
    #
    #   Proportional residual variability                0.517  <- carried here
    #   Additive residual variability, ng/mL             383    <- carried here
    #   Proportional residual variability phase I        0.244
    #   Proportional residual variability sparse data    0.612
    #   Additive residual variability sparse data, ng/mL 267
    #
    # This file carries the HSCT-recipient rich-sampling pair, which is the
    # population the model file describes; the other three are recorded above
    # for anyone re-fitting. Unlike the IIV column these are on the SD scale:
    # the additive rows carry an explicit ng/mL unit, and the corresponding
    # row of the companion phase I model's Table 2 is headed "Proportional
    # residual variability, %" with the value 0.283, which is a 28.3%
    # proportional error rather than a 53% one.
    # ------------------------------------------------------------------
    propSd <- 0.517
    label("Proportional residual error, HSCT recipients with rich sampling (fraction)")  # Table 3 'Proportional residual variability' 0.517 (RSE 12.2%, 95% CI 0.394-0.641)
    addSd <- 383
    label("Additive residual error, HSCT recipients with rich sampling (ng/mL)")  # Table 3 'Additive residual variability, ng/ml' 383.0 (RSE 22.2%, 95% CI 216-550)
  })

  model({
    # Decompose the integer-valued OCC column into binary occasion indicators
    # for IOV multiplexing on bioavailability. OCC = 1..8 selects the matching
    # per-occasion eta; OCC = 0 or any value outside 1..8 zeroes every
    # indicator and leaves bioavailability at its IIV-only value, which is the
    # configuration the source used for its own exposure simulations.
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    oc5 <- (OCC == 5)
    oc6 <- (OCC == 6)
    oc7 <- (OCC == 7)
    oc8 <- (OCC == 8)
    iov_fdepot <- oc1 * etaiov_fdepot_1 + oc2 * etaiov_fdepot_2 +
      oc3 * etaiov_fdepot_3 + oc4 * etaiov_fdepot_4 +
      oc5 * etaiov_fdepot_5 + oc6 * etaiov_fdepot_6 +
      oc7 * etaiov_fdepot_7 + oc8 * etaiov_fdepot_8

    # --- Individual parameters (Prohn 2021 Table 3, Figure 1b) -----------
    cl <- exp(lcl + e_csa_cl * CONMED_CSA + etalcl)
    vc <- exp(lvc)
    vp <- exp(lvp + e_asian_vp * RACE_ASIAN + etalvp)
    q  <- exp(lq)

    # Absorption rate and bioavailability both switch on transplant status.
    # The typical values are switched on the log scale and the single IIV eta
    # is then applied to whichever branch is active, which is how NONMEM
    # multiplexes a two-THETA parameter carrying one ETA.
    ka <- exp(lka * TX_HCT + lka_hp * (1 - TX_HCT) + etalka)
    fdepot <- exp(
      (lfdepot + e_csa_fdepot * CONMED_CSA) * TX_HCT +
        lfdepot_hp * (1 - TX_HCT) +
        etalfdepot + iov_fdepot
    )

    tlag <- exp(ltlag)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # --- ODE system (Figure 1b) ------------------------------------------
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    f(depot)    <- fdepot
    alag(depot) <- tlag

    # Figure 1b: Cp = (A_central / Vc) * 1000. Amounts are in mg and volumes
    # in L, so central/vc is mg/L and the factor of 1000 converts to ng/mL,
    # the unit in which Table 3's additive residual error is reported.
    Cc <- central / vc * 1000
    Cc ~ prop(propSd) + add(addSd)
  })
}
