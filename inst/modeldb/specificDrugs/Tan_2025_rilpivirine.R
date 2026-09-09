Tan_2025_rilpivirine <- function() {
  description <- "One-compartment population PK model with first-order absorption for long-acting intramuscular rilpivirine in adults with HIV-1 followed in routine outpatient care (Tan 2025, JABS-PKInSITE). Absorption is roughly half as fast as its cabotegravir counterpart from the same co-packaged injection (ka = 0.000500 1/h, absorption half-life 58 days) and two orders of magnitude below elimination, so the profile is strongly flip-flop and unusually flat across the 8-week interval. NO covariate was retained: unlike cabotegravir, rilpivirine absorption was insensitive to whether the injectate actually landed in muscle or in subcutaneous fat, which the authors attribute to differing physicochemical properties of the two nanosuspensions. Sex, body mass index, age, height, skin-to-muscle thickness and depot location were all screened and rejected (see covariatesDataExcluded). Body weight enters CL/F and V/F allometrically at the 70 kg reference printed in the paper's own parameter-table headers. ka carries between-subject variability and CL/F carries both between-subject and inter-occasion variability across the three study injections."
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
    "rilpivirine block. The companion cabotegravir model from the same",
    "analysis is modellib('Tan_2025_cabotegravir').",
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
      notes              = "Allometric power scaling on CL/F and V/F at a 70 kg reference. The reference weight is stated by Tan 2025 Table 2 itself, whose structural rows are headed 'CL/F (liters/h/70 kg)' and 'V/F (liters/70 kg)' for both analytes; a per-70-kg volume unit is meaningful only if the parameter is weight-normalised. The paper never prints the exponents, so the theory-based values 0.75 (CL/F) and 1.0 (V/F) are used and are encoded as fixed(). Rilpivirine gives the cleaner of the two arithmetic checks that support this reading: at the cohort median 84 kg the predicted typical steady-state trough is 56.9 ng/mL against an observed median of 56.0 ng/mL (Tan 2025 Table 1), whereas deleting the weight term puts it at 65.3 ng/mL, 17% high -- and the companion cabotegravir model is biased high in the same direction by the same omission. See the vignette Errata for the full calculation. Cohort weight 52-114 kg, median 84 kg (Tan 2025 Table 1). The paper's separate statement that 'BMI and weight did not show any correlation with Ctrough' concerns the observed-trough correlation analysis and the covariate SEARCH, not this a-priori structural normalisation.",
      source_name        = "Weight (kg)"
    ),
    OCC = list(
      description        = "Integer injection-occasion index driving the inter-occasion variability on CL/F",
      units              = "(count)",
      type               = "categorical",
      reference_category = "n/a -- decomposed inside model() into three mutually exclusive binary indicators multiplied against the per-occasion etaiov_lcl_<k> slots",
      notes              = "Values 1, 2 or 3. Tan 2025 Table 2 reports one inter-occasion variability magnitude for CL/F but never states an occasion count or how an occasion was delimited, and no control stream was deposited. Three occasions are encoded because Tan 2025 Methods fixes the observation window exactly: 'Participants were observed over 16 weeks, corresponding to 3 clinic appointments (8-week intervals) with corresponding injection administration', so each participant contributed at most three injections to the analysis. All three slots share the single estimated magnitude, the standard NONMEM $OMEGA BLOCK(1) SAME idiom, so occasions 2 and 3 are fixed() to the occasion-1 variance. Setting OCC = 0 on every record zeroes all three indicators and switches inter-occasion variability off. Note that the parameter carrying inter-occasion variability differs between the two analytes of this paper: ka for cabotegravir, CL/F for rilpivirine.",
      source_name        = "OCC"
    )
  )

  covariatesDataExcluded <- list(
    ROUTE_SC = list(
      description = "Ultrasound-determined subcutaneous location of the depot laid down by the dose being absorbed: 1 = the injectate was found primarily in subcutaneous tissue, 0 = intramuscular or mixed",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened and NOT retained -- this is the paper's headline negative result, and the contrast with the companion cabotegravir model is the point of the study. The identical ultrasound covariate (75 intramuscular / 40 subcutaneous / 19 mixed depots over 134 imaged injections) was offered to both analytes from the same co-packaged injection; it cut cabotegravir ka by 56.3% and did nothing measurable to rilpivirine. Tan 2025 Results: 'None of the tested covariates met criteria for inclusion in the model.' Tan 2025 Discussion: 'Subcutaneous deposition slowed the absorption of CAB by 56.3% but did not slow absorption of RPV', and 'The differing behavior of CAB and RPV serves to reflect the impact that the physicochemical properties (eg, pH, particle size and charge, solubility) of drugs may have on absorption.' Table 2's rilpivirine block has no covariate row, so there is no point estimate to encode. To explore the counterfactual, apply the cabotegravir coefficient from modellib('Tan_2025_cabotegravir').",
      source_name = "Injection depot disposition = Subcutaneous (SC)"
    ),
    SEXF = list(
      description = "Sex at birth indicator, 1 = female, 0 = male",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened and NOT retained. Tan 2025 Results found women's observed rilpivirine troughs higher than men's (62.5 vs 53.0 ng/mL, P = .02) but the effect did not survive the model: 'Although our exploratory analysis suggested a potential small sex difference with RPV concentrations (P < .05), no injection site covariate, nor other covariate, was identified in the modelling.' The Discussion adds that the one prior report of a rilpivirine sex effect 'was small and only affected 1 of the 2 components of absorption with limited effect on trough concentrations'. No coefficient appears in Table 2. Cohort: 8 of 31 participants (26%) were female.",
      source_name = "Gender"
    ),
    BMI = list(
      description = "Body mass index at baseline",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Screened and NOT retained. Tan 2025 Results: 'BMI and weight did not show any correlation with Ctrough for either of the compounds', and 'None of the tested covariates met criteria for inclusion in the model.' Cohort median 26.3 kg/m^2 (range 19.0-38.1); 8 of 31 participants (25.8%) had BMI > 30 kg/m^2.",
      source_name = "BMI (kg/m2)"
    ),
    INJDEPTH = list(
      description = "Ultrasound-measured skin-to-muscle thickness at the ventrogluteal injection site",
      units       = "mm",
      type        = "continuous",
      notes       = "Screened and NOT retained. Tan 2025 Methods tested skin-to-muscle thickness 'as both continuous and categorical covariates' for both analytes; for cabotegravir the 20 mm dichotomy reached significance but lost to depot location, and for rilpivirine no covariate met the forward-selection criterion at all. Cohort median 18.5 mm (range 3.7-50.6; men 16.9, women 34.6). NOTE ON THE NAME: INJDEPTH has no entry in inst/references/covariate-columns.md. It is documentation-only here, is never referenced in model(), and would require register ratification before any model actually used it; it is spelled to sit alongside the existing INJSITE_* family.",
      source_name = "Skin-to-muscle thickness (mm)"
    ),
    AGE = list(
      description = "Age at enrolment",
      units       = "years",
      type        = "continuous",
      notes       = "Screened and NOT retained. Tan 2025 Methods lists age among the collected covariates; Results reports that no tested covariate entered the rilpivirine model. Cohort median 45 years (range 23-72).",
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
    depot   = list(analyte = "rilpivirine", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "rilpivirine", units = "mg", specimen = "serum",              verified = TRUE)
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
    dose_range     = "Not stated in the paper. Participants received the fixed-dose co-packaged cabotegravir/rilpivirine nanosuspension by ventrogluteal intramuscular injection at 8-week intervals; the paper names the Cabenuva pack and its 23 gauge / 38 mm needle but never prints a milligram dose. The vignette therefore simulates the labelled every-2-months maintenance dose of 900 mg rilpivirine and flags it as a non-paper-derived input.",
    regions        = "Australia (Royal Perth Hospital immunology outpatient clinic, Perth, Western Australia)",
    notes          = "JABS-PKInSITE, a prospective non-interventional observational study. 31 participants recruited 4 October 2023 to 21 March 2024 and followed over 16 weeks spanning 3 clinic visits and 3 injections. 141 serum samples were assayed, of which 78 were pre-injection troughs at approximately 8 weeks; additional samples were sought at 24 h and at 2, 4 and 6 weeks post injection (12, 19, 17 and 15 samples respectively). 134 injection sites were imaged by post-injection ultrasound from 67 visits: 75 (56%) intramuscular, 40 (30%) subcutaneous, 19 (14%) mixed. Observed rilpivirine trough concentrations had a median of 56.0 ng/mL (range 20-257), 53.0 ng/mL in men and 62.5 ng/mL in women (P = .02). Serum was assayed by validated LC-MS/MS with a 5-500 ng/mL quantification range for rilpivirine. Estimation used NONMEM 7.5.1; evaluation used a 1000-sample bootstrap and prediction- and variability-corrected visual predictive checks. Tan 2025 Discussion reports that these real-world troughs were 23% lower than those of the phase 3 ATLAS-2M study and that 38.5% of rilpivirine measurements sat below 4x the protein-adjusted 90% inhibitory concentration."
  )

  ini({
    # =====================================================================
    # All point estimates are the 'Mean' column of Tan 2025 Table 2,
    # rilpivirine block. The bootstrap median and 95% empirical confidence
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
    # (value/100)^2.
    # =====================================================================

    # --- Absorption ------------------------------------------------------
    # ka is two orders of magnitude below kel = CL/V = 0.0837 1/h:
    # ln(2)/ka = 1386 h = 57.8 days of absorption half-life against
    # ln(2)/kel = 8.3 h of elimination half-life. Tan 2025 Results:
    # "Compared to CAB the absorption rate of RPV was approximately half
    # (0.000500 vs 0.000972 hours-1), consistent with a flatter
    # pharmacokinetic profile."
    lka <- log(0.000500) ; label("Absorption rate constant (ka, 1/h)")  # Table 2 'k a (h-1) 0.000500'; bootstrap median 0.000507 [0.000402-0.000624]

    # --- Disposition -----------------------------------------------------
    lcl <- log(7.24) ; label("Apparent clearance at 70 kg (CL/F, L/h)")                     # Table 2 'CL/F (liters/h/70 kg) 7.24'; bootstrap median 7.16 [6.28-8.16]
    lvc <- log(86.5) ; label("Apparent central volume of distribution at 70 kg (V/F, L)")   # Table 2 'V/F (liters/70 kg) 86.5'; bootstrap median 86.0 [33.8-183]

    # --- Allometric weight scaling ---------------------------------------
    # NOT PRINTED BY THE PAPER. Tan 2025 states the reference weight in the
    # Table 2 row headers themselves ("liters/h/70 kg", "liters/70 kg") but
    # never gives the exponents, and no control stream was deposited. The
    # theory-based values are used and held fixed; see the WT covariateData
    # entry and the vignette Errata for the trough arithmetic that supports
    # the reading. A user who prefers the unscaled reading can set both
    # exponents to 0.
    e_wt_cl <- fixed(0.75) ; label("Power exponent of body weight on CL/F, reference 70 kg (unitless)")  # NOT in Tan 2025; theory-based allometric exponent, fixed. Reference 70 kg is from the Table 2 row header 'CL/F (liters/h/70 kg)'
    e_wt_vc <- fixed(1)    ; label("Power exponent of body weight on V/F, reference 70 kg (unitless)")   # NOT in Tan 2025; theory-based allometric exponent, fixed. Reference 70 kg is from the Table 2 row header 'V/F (liters/70 kg)'

    # --- Between-subject variability -------------------------------------
    etalka ~ 0.1024  # Table 2 'IIV in k a  32 [15]'; (32/100)^2 = 0.1024. Bootstrap median 30 [13-41]. Shrinkage 15%
    etalcl ~ 0.1089  # Table 2 'IIV in CL/F 33 [45]'; (33/100)^2 = 0.1089. Bootstrap median 34 [13-53]. Shrinkage 45%

    # --- Inter-occasion variability on CL/F ------------------------------
    # One magnitude shared by three occasions (the three 8-weekly study
    # injections; see the OCC covariateData entry for why three). Occasions 2
    # and 3 are fixed to the occasion-1 variance, the NONMEM
    # $OMEGA BLOCK(1) SAME idiom.
    etaiov_lcl_1 ~ 0.0256         # Table 2 'IOV in CL/F 16 [47]'; (16/100)^2 = 0.0256. Bootstrap median 18 [5-28]. Shrinkage 47%
    etaiov_lcl_2 ~ fixed(0.0256)  # same magnitude, second of three occasions
    etaiov_lcl_3 ~ fixed(0.0256)  # same magnitude, third of three occasions

    # --- Residual error --------------------------------------------------
    # Tan 2025 does not name the residual-error form. It is read as
    # proportional because Table 2 reports the residual row 'RV' on the same
    # 100 * sqrt(estimate) percentage scale as the log-normal IIV and IOV
    # rows rather than in ng/mL, and because the assayed concentrations span
    # more than a decade (20-257 ng/mL). See vignette Errata.
    propSd <- 0.25 ; label("Proportional residual standard deviation (fraction)")  # Table 2 'RV 25 [20]'; 25/100 = 0.25. Bootstrap median 24 [17-29]. Shrinkage 20%
  })

  model({
    # =====================================================================
    # 1. Inter-occasion variability on CL/F, multiplexed by the occasion
    #    indicator. Exactly one indicator is non-zero on a record with
    #    OCC in 1:3, so iov_cl is that occasion's eta; OCC = 0 switches
    #    inter-occasion variability off entirely.
    # =====================================================================
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    iov_cl <- oc1 * etaiov_lcl_1 + oc2 * etaiov_lcl_2 + oc3 * etaiov_lcl_3

    # =====================================================================
    # 2. Individual parameters. No covariate was retained on any of them
    #    beyond the a-priori allometric weight normalisation.
    # =====================================================================
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl + iov_cl) * (WT / 70)^e_wt_cl
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
