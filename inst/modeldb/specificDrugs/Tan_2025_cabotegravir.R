Tan_2025_cabotegravir <- function() {
  description <- "One-compartment population PK model with first-order absorption for long-acting intramuscular cabotegravir in adults with HIV-1 followed in routine outpatient care (Tan 2025, JABS-PKInSITE). Absorption is far slower than elimination (ka = 0.000972 1/h, absorption half-life 30 days, against an elimination half-life near 47 h), so the profile is flip-flop and the 8-weekly trough is absorption-limited. The one retained covariate is the ultrasound-determined location of the injected depot: when a gluteal injection deposited primarily into subcutaneous tissue rather than muscle, ka fell by 56.3%, which RAISES rather than lowers the steady-state trough. Body weight enters CL/F and V/F allometrically at the 70 kg reference printed in the paper's own parameter-table headers. Sex, body mass index and skin-to-muscle thickness were screened and not retained (see covariatesDataExcluded); the authors report that sex and BMI lost significance once depot location entered the model, so the long-observed female and high-BMI cabotegravir effects are re-expressed here as the probability that an intended intramuscular injection actually lands intramuscularly. CL/F carries between-subject variability, and ka carries both between-subject and inter-occasion variability across the three study injections."
  reference <- paste(
    "Tan B, John M, Castley A, Williams L, Joyce D, Nolan D, O'Halloran S,",
    "Salman S.",
    "Exploring the interaction between injection site and biological sex on",
    "the real-world population pharmacokinetics of long-acting cabotegravir",
    "and rilpivirine in people with HIV.",
    "Open Forum Infect Dis. 2025;12(10):ofaf614. doi:10.1093/ofid/ofaf614.",
    "Parameter estimates are from Tan 2025 Table 2 ('Final Population",
    "Pharmacokinetic Estimates and Bootstrap Results for Cabotegravir and",
    "Rilpivirine in Patients With HIV Receiving Long-acting Injections'),",
    "cabotegravir block. The companion rilpivirine model from the same",
    "analysis is modellib('Tan_2025_rilpivirine').",
    sep = " "
  )
  vignette <- "Tan_2025_cabotegravir_rilpivirine"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric power scaling on CL/F and V/F at a 70 kg reference. The reference weight is stated by Tan 2025 Table 2 itself, whose structural rows are headed 'CL/F (liters/h/70 kg)' and 'V/F (liters/70 kg)' for both analytes; a per-70-kg volume unit is meaningful only if the parameter is weight-normalised. The paper never prints the exponents, so the theory-based values 0.75 (CL/F) and 1.0 (V/F) are used and are encoded as fixed() -- see the model file's ini() comment and the vignette Errata for the arithmetic that supports this reading (at the cohort median 84 kg both drugs' predicted typical troughs land within a few percent of the observed medians of Table 1, whereas dropping the weight term leaves both 17-23% high in the same direction). Cohort weight 52-114 kg, median 84 kg (Tan 2025 Table 1). Note that the paper's separate statement that 'BMI and weight did not show any correlation with Ctrough' concerns the observed trough correlation analysis and the covariate SEARCH, not this a-priori structural normalisation.",
      source_name        = "Weight (kg)"
    ),
    ROUTE_SC = list(
      description        = "Ultrasound-determined subcutaneous location of the depot laid down by the dose being absorbed: 1 = the injectate was found primarily in subcutaneous tissue, 0 = intramuscular or mixed",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (intramuscular, and also the 'mixed' depots that spanned both muscle and subcutaneous tissue -- see notes)",
      notes              = "Per-dose-record covariate. Every injection in this study was INTENDED to be an intramuscular ventrogluteal injection with the 23 gauge / 38 mm needle supplied in the Cabenuva pack; ROUTE_SC records where the injectate was actually found by a post-injection ultrasound performed within 15 minutes, not a route the clinician chose. Of 134 imaged injections, 75 (56%) were intramuscular, 40 (30%) subcutaneous and 19 (14%) mixed, and all injections were subcutaneous once skin-to-muscle thickness exceeded 30 mm (Tan 2025 Results and Table 1). Tan 2025 Methods screened the site as a three-level categorical covariate (intramuscular, mixed, subcutaneous) but the final model reports a SINGLE coefficient, 'Effect of subcutaneous location on k a (%)', so the retained relationship is binary and the 19 mixed depots fall in the reference category together with the intramuscular ones; that grouping is an assumption of this encoding and is recorded in the vignette Errata. Tan 2025 Results describes the effect as attaching to the 'subcutaneous location of prior dose', which is the natural per-dose reading for an 8-weekly trough: the depot whose absorption governs a pre-injection sample is the one laid down by the preceding injection, and in an event table that covariate rides on the dose record it modifies. Because ka enters below the elimination rate constant the drug is flip-flop, so the 56.3% slower subcutaneous absorption INCREASES the trough; Tan 2025 Discussion makes this point explicitly and uses it to argue that ectopic subcutaneous deposition does not threaten efficacy.",
      source_name        = "Injection depot disposition = Subcutaneous (SC)"
    ),
    OCC = list(
      description        = "Integer injection-occasion index driving the inter-occasion variability on ka",
      units              = "(count)",
      type               = "categorical",
      reference_category = "n/a -- decomposed inside model() into three mutually exclusive binary indicators multiplied against the per-occasion etaiov_lka_<k> slots",
      notes              = "Values 1, 2 or 3. Tan 2025 Table 2 reports one inter-occasion variability magnitude for ka but never states an occasion count or how an occasion was delimited, and no control stream was deposited. Three occasions are encoded because Tan 2025 Methods fixes the observation window exactly: 'Participants were observed over 16 weeks, corresponding to 3 clinic appointments (8-week intervals) with corresponding injection administration', so each participant contributed at most three injections to the analysis. All three slots share the single estimated magnitude, the standard NONMEM $OMEGA BLOCK(1) SAME idiom, so occasions 2 and 3 are fixed() to the occasion-1 variance. Setting OCC = 0 on every record zeroes all three indicators and switches inter-occasion variability off.",
      source_name        = "OCC"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Sex at birth indicator, 1 = female, 0 = male",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened and NOT retained. Tan 2025 Methods lists gender among the covariates collected and tested, and the paper's premise is the well-replicated observation that cabotegravir absorption is about 50% slower in women. Tan 2025 Discussion reports the outcome of the search: 'When depot location and skin to muscle depth were incorporated into our population pharmacokinetic models, sex, and BMI lost their significance.' No sex coefficient appears in Table 2, so there is no point estimate to encode. The mechanism the paper offers is mediation rather than confounding: women in this cohort had 17.7 mm greater median skin-to-muscle thickness at the same weight or BMI and 72.4% of their injections landed subcutaneously versus 18.1% of men's, so ROUTE_SC carries the sex effect. Cohort: 8 of 31 participants (26%) were female.",
      source_name = "Gender"
    ),
    BMI = list(
      description = "Body mass index at baseline",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Screened and NOT retained, for the same reason as SEXF and with the same Discussion sentence as the source. Tan 2025 Results additionally reports that 'BMI and weight did not show any correlation with Ctrough for either of the compounds' in the observed-data analysis. Skin-to-muscle thickness correlated positively with BMI in both sexes (r^2 = 0.260 in men, 0.272 in women), which is the path by which BMI acts on absorption in this model. Cohort median 26.3 kg/m^2 (range 19.0-38.1); 8 of 31 participants (25.8%) had BMI > 30 kg/m^2.",
      source_name = "BMI (kg/m2)"
    ),
    INJDEPTH = list(
      description = "Ultrasound-measured skin-to-muscle thickness at the ventrogluteal injection site",
      units       = "mm",
      type        = "continuous",
      notes       = "Screened, SIGNIFICANT, and still not retained. Tan 2025 Results reports that skin-to-muscle thickness dichotomised at a 20 mm cutoff was associated with slower absorption (dOFV = -6.824, P < .01) but that subcutaneous depot location gave the larger drop (dOFV = -7.237, P < .01), so 'Given a greater reduction in OFV, and with consideration of biological plausibility, subcutaneous location of prior dose was included in the model'. The two are near-collinear -- every injection with a skin-to-muscle thickness above 30 mm was subcutaneous -- so only one could enter. Table 2 therefore contains no thickness coefficient and none is encoded here. Cohort median 18.5 mm (range 3.7-50.6; men 16.9, women 34.6). NOTE ON THE NAME: INJDEPTH has no entry in inst/references/covariate-columns.md. It is documentation-only here, is never referenced in model(), and would require register ratification before any model actually used it; it is spelled to sit alongside the existing INJSITE_* family.",
      source_name = "Skin-to-muscle thickness (mm)"
    ),
    AGE = list(
      description = "Age at enrolment",
      units       = "years",
      type        = "continuous",
      notes       = "Screened and NOT retained. Tan 2025 Methods lists age among the collected covariates and Results states that 'No further significant covariate relationships were identified in the modeling process' beyond depot location. Cohort median 45 years (range 23-72).",
      source_name = "Age (years)"
    ),
    HT = list(
      description = "Body height at baseline",
      units       = "cm",
      type        = "continuous",
      notes       = "Screened and NOT retained; collected only as an input to BMI. Tan 2025 Methods lists height among the collected covariates; no height coefficient appears in Table 2. Cohort median 173 cm (range 159-187).",
      source_name = "Height (cm)"
    )
  )

  compartmentData <- list(
    depot   = list(analyte = "cabotegravir", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "cabotegravir", units = "mg", specimen = "serum",              verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 31L,
    n_studies      = 1L,
    age_range      = "23-72 years",
    age_median     = "45 years",
    weight_range   = "52-114 kg",
    weight_median  = "84 kg",
    sex_female_pct = 25.8,
    race_ethnicity = "Not reported",
    disease_state  = "Adults with virologically suppressed HIV-1 (at least one recent HIV-1 RNA < 40 copies/mL) already established on long-acting injectable cabotegravir plus rilpivirine as routine care; median duration on injectable therapy at the trough samples was 56 weeks",
    dose_range     = "Not stated in the paper. Participants received the fixed-dose co-packaged cabotegravir/rilpivirine nanosuspension by ventrogluteal intramuscular injection at 8-week intervals; the paper names the Cabenuva pack and its 23 gauge / 38 mm needle but never prints a milligram dose. The vignette therefore simulates the labelled every-2-months maintenance dose of 600 mg cabotegravir and flags it as a non-paper-derived input.",
    regions        = "Australia (Royal Perth Hospital immunology outpatient clinic, Perth, Western Australia)",
    notes          = "JABS-PKInSITE, a prospective non-interventional observational study. 31 participants recruited 4 October 2023 to 21 March 2024 and followed over 16 weeks spanning 3 clinic visits and 3 injections. 141 serum samples were assayed, of which 78 were pre-injection troughs at approximately 8 weeks; additional samples were sought at 24 h and at 2, 4 and 6 weeks post injection (12, 19, 17 and 15 samples respectively). 134 injection sites were imaged by post-injection ultrasound (< 15 minutes after administration, single operator, Clarius L15HD3 5-15 MHz linear probe) from 67 visits: 75 (56%) intramuscular, 40 (30%) subcutaneous, 19 (14%) mixed. Observed cabotegravir trough concentrations had a median of 1390 ng/mL (range 340-4423), 1292 ng/mL in men and 1727 ng/mL in women (P < .05). Serum was assayed by validated LC-MS/MS with a 100-10000 ng/mL quantification range for cabotegravir. Estimation used NONMEM 7.5.1; evaluation used a 1000-sample bootstrap and prediction- and variability-corrected visual predictive checks. No virological failure or rebound occurred during the study."
  )

  ini({
    # =====================================================================
    # All point estimates are the 'Mean' column of Tan 2025 Table 2,
    # cabotegravir block. The bootstrap median and 95% empirical confidence
    # interval from the same table's second column are quoted alongside each
    # value as the published uncertainty. No NONMEM control stream was
    # deposited: the PMC supplement (ofaf614_supplementary_data.zip) holds
    # only Supplementary Figures 1-4 and their captions, so Table 2 is the
    # sole source for every number below.
    #
    # VARIABILITY SCALE. Tan 2025 prints its own back-transform in the
    # footnote above Table 2: "Variability parameters are presented as
    # 100% x sqrt(variability estimate)". Every percentage in the
    # "Variable model parameters" block is therefore 100 * omega on the
    # STANDARD-DEVIATION scale, and the variance nlmixr2 wants is
    # (value/100)^2. This settles the usual SD-versus-variance ambiguity
    # arithmetically rather than by convention.
    # =====================================================================

    # --- Absorption ------------------------------------------------------
    # ka sits an order of magnitude BELOW kel = CL/V = 0.0147 1/h, so
    # cabotegravir is flip-flop: ln(2)/ka = 713 h = 29.7 days of absorption
    # half-life against ln(2)/kel = 47 h of elimination half-life. Tan 2025
    # Discussion cross-checks this ka against the two published comparators,
    # 0.00102 1/h in the Swiss real-world cohort and 0.000642 1/h in the
    # registrational analysis.
    lka <- log(0.000972) ; label("Absorption rate constant (ka, 1/h)")  # Table 2 'k a (h-1) 0.000972'; bootstrap median 0.000978 [0.00069-0.00118]

    # --- Disposition -----------------------------------------------------
    lcl <- log(0.136) ; label("Apparent clearance at 70 kg (CL/F, L/h)")                     # Table 2 'CL/F (liters/h/70 kg) 0.136'; bootstrap median 0.1361 [0.1167-0.1628]
    lvc <- log(9.28)  ; label("Apparent central volume of distribution at 70 kg (V/F, L)")   # Table 2 'V/F (liters/70 kg) 9.28'; bootstrap median 9.16 [5.95-17.77]. Tan 2025 Discussion notes the registrational two-compartment central and peripheral volumes "approximately summed to our single compartment estimate"

    # --- Allometric weight scaling ---------------------------------------
    # NOT PRINTED BY THE PAPER. Tan 2025 states the reference weight in the
    # Table 2 row headers themselves ("liters/h/70 kg", "liters/70 kg") but
    # never gives the exponents, and no control stream was deposited. The
    # theory-based values are used and held fixed. Supporting arithmetic,
    # laid out in full in the vignette Errata: at the cohort median weight
    # of 84 kg these exponents put the typical steady-state trough at
    # 1493 ng/mL for an intramuscular depot (observed cohort median
    # 1390 ng/mL over a 70/30 mix of non-subcutaneous and subcutaneous
    # injections) and, for the companion rilpivirine model, 56.9 ng/mL
    # against an observed 56.0 ng/mL. Deleting the weight term entirely
    # leaves BOTH analytes 17-23% high in the same direction (1706 ng/mL
    # and 65.3 ng/mL), which is the
    # signature of a missing normalisation from the 70 kg reference to the
    # cohort's heavier median. A user who prefers the unscaled reading can
    # set both exponents to 0.
    e_wt_cl <- fixed(0.75) ; label("Power exponent of body weight on CL/F, reference 70 kg (unitless)")  # NOT in Tan 2025; theory-based allometric exponent, fixed. Reference 70 kg is from the Table 2 row header 'CL/F (liters/h/70 kg)'
    e_wt_vc <- fixed(1)    ; label("Power exponent of body weight on V/F, reference 70 kg (unitless)")   # NOT in Tan 2025; theory-based allometric exponent, fixed. Reference 70 kg is from the Table 2 row header 'V/F (liters/70 kg)'

    # --- Depot-location effect on absorption -----------------------------
    # Entered as a fractional change on the typical ka, matching the way the
    # paper reports it: a percentage reduction rather than a log-shift or a
    # ratio. ka_SC = ka * (1 - 0.563) = 0.000425 1/h, an absorption half-life
    # of 68 days.
    e_route_sc_ka <- -0.563 ; label("Fractional change in ka for a subcutaneous depot versus intramuscular or mixed (unitless)")  # Table 2 'Effect of subcutaneous location on k a (%) -56.3'; bootstrap median -55.3 [-32.1 to -71.8]. Tan 2025 Results: "subcutaneous location of prior dose was included in the model, associated with a 56.3% reduction in the rate of absorption"

    # --- Between-subject variability -------------------------------------
    etalka ~ 0.1089  # Table 2 'IIV in k a  33 [35]'; (33/100)^2 = 0.1089. Bootstrap median 33 [11-48]. Shrinkage 35%
    etalcl ~ 0.1089  # Table 2 'IIV in CL/F 33 [16]'; (33/100)^2 = 0.1089. Bootstrap median 32 [19-44]. Shrinkage 16%

    # --- Inter-occasion variability on ka --------------------------------
    # One magnitude shared by three occasions (the three 8-weekly study
    # injections; see the OCC covariateData entry for why three). Occasions 2
    # and 3 are fixed to the occasion-1 variance, the NONMEM
    # $OMEGA BLOCK(1) SAME idiom.
    etaiov_lka_1 ~ 0.0625         # Table 2 'IOV in k a  25 [33]'; (25/100)^2 = 0.0625. Bootstrap median 25 [0.2-60]. Shrinkage 33%
    etaiov_lka_2 ~ fixed(0.0625)  # same magnitude, second of three occasions
    etaiov_lka_3 ~ fixed(0.0625)  # same magnitude, third of three occasions

    # --- Residual error --------------------------------------------------
    # Tan 2025 does not name the residual-error form. It is read as
    # proportional because Table 2 reports the residual row 'RV' on the same
    # 100 * sqrt(estimate) percentage scale as the log-normal IIV and IOV
    # rows rather than in ng/mL, and because the assayed concentrations span
    # more than a decade (340-4423 ng/mL). See vignette Errata.
    propSd <- 0.17 ; label("Proportional residual standard deviation (fraction)")  # Table 2 'RV 17 [31]'; 17/100 = 0.17. Bootstrap median 16 [12-21]. Shrinkage 31%
  })

  model({
    # =====================================================================
    # 1. Inter-occasion variability on ka, multiplexed by the occasion
    #    indicator. Exactly one indicator is non-zero on a record with
    #    OCC in 1:3, so iov_ka is that occasion's eta; OCC = 0 switches
    #    inter-occasion variability off entirely.
    # =====================================================================
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    iov_ka <- oc1 * etaiov_lka_1 + oc2 * etaiov_lka_2 + oc3 * etaiov_lka_3

    # =====================================================================
    # 2. Individual parameters. The depot-location effect is a fractional
    #    change on the typical ka, so a subcutaneous depot multiplies ka by
    #    1 - 0.563 = 0.437.
    # =====================================================================
    ka <- exp(lka + etalka + iov_ka) * (1 + e_route_sc_ka * ROUTE_SC)
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    vc <- exp(lvc) * (WT / 70)^e_wt_vc

    kel <- cl / vc

    # =====================================================================
    # 3. One-compartment disposition with first-order absorption from the
    #    injected depot. F is not identifiable from injection-only data, so
    #    CL and V are apparent (CL/F, V/F) and no f() is applied.
    # =====================================================================
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # =====================================================================
    # 4. Observation. Doses are in mg and vc is in L, so central / vc is
    #    mg/L; the factor 1000 converts to the ng/mL in which Tan 2025
    #    reports every concentration.
    # =====================================================================
    Cc <- 1000 * central / vc
    Cc ~ prop(propSd)
  })
}
