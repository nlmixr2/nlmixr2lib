Tsirizani_2025_ritonavir <- function() {
  description <- paste(
    "Two-compartment population PK model for low-dose oral ritonavir used as a",
    "pharmacokinetic booster of lopinavir, atazanavir or darunavir in African",
    "children with HIV failing first-line ART (CHAPAS-4 trial, ISRCTN22964075;",
    "170 children aged 3.16-15.6 years and 14.2-64.2 kg from Zambia, Uganda and",
    "Zimbabwe). Absorption is an absorption lag followed by sequential",
    "zero-order (duration D1) then first-order (ka) input. Clearance and both",
    "volumes are allometrically scaled to fat-free mass with exponents fixed at",
    "0.75 and 1, referenced to FFM 21.0 kg (the cohort median, corresponding to",
    "a child weighing 26 kg). Relative bioavailability is fixed to 1 in the",
    "darunavir reference arm, so all clearances and volumes are apparent",
    "(CL/F, V/F) on that reference. Companion protease inhibitor is the only",
    "retained covariate on disposition: atazanavir raises relative",
    "bioavailability by 137% and clearance by 20.7%, and lopinavir lowers",
    "relative bioavailability by 23.4%. In the twice-daily lopinavir/ritonavir",
    "arm the evening dose absorbs 3.61-fold more slowly (fold-change on the",
    "absorption lag). Between-subject variability was retained only on",
    "clearance; between-occasion variability is carried on relative",
    "bioavailability, absorption lag, zero-order duration and ka, with the",
    "bioavailability BOV standard deviation inflated 2.14-fold on the",
    "unwitnessed dose preceding the sampling window. No effect of the NRTI",
    "backbone (including tenofovir alafenamide) or of age was found. The",
    "authors note the model should not be used to simulate ritonavir doses",
    "above 100 mg, because clearance saturation reported at higher doses could",
    "not be characterised from boosting-dose data alone."
  )
  reference <- paste(
    "Tsirizani L, Waalewijn H, Szubert A, Mulenga V, Chabala C,",
    "Bwakura-Dangarembizi M, Chitsamatanga M, Rutebarika DA, Musiime V,",
    "Kasozi M, Lugemwa A, McIlleron HM, Burger DM, Gibb DM, Colbers A,",
    "Denti P, Wasmann RE, the CHAPAS-4 trial team (2025).",
    "Population pharmacokinetics of ritonavir as a booster of lopinavir,",
    "atazanavir, or darunavir in African children with HIV.",
    "Antimicrob Agents Chemother. doi:10.1128/aac.00771-25.",
    "Parameter values from Table 3; model structure and the fat-free-mass",
    "reference value from the NONMEM control stream in supplemental Data S1",
    "(AAC00771-25-s0001.docx).",
    sep = " "
  )
  vignette <- "Tsirizani_2025_ritonavir"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    FFM = list(
      description        = "Fat-free mass, the body-size descriptor for allometric scaling of clearance and volume",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Reference value 21.0 kg, taken from the supplemental NONMEM control",
        "stream (Data S1, $PK block 'TVFFM = 21.0 ;MEDIAN'). Table 3 footnote",
        "b states the typical values 'refer to a child weighing 26 kg', so",
        "FFM 21.0 kg is the cohort-median fat-free mass of a 26 kg child in",
        "this population (FFM/WT = 0.808). Fat-free mass was preferred over",
        "total body weight and over fat mass as the size descriptor (Results,",
        "'Population pharmacokinetic analysis': dOFV = -5.0 versus weight).",
        "IMPORTANT: the paper does NOT report the equation used to predict FFM",
        "from weight, height, age and sex. FFM entered the analysis as a",
        "pre-computed data column ($INPUT 'FFM') and the Methods cite only",
        "Holford & Anderson 2017 (reference 21) for the allometric-size theory,",
        "not a specific FFM prediction equation. A downstream user must supply",
        "FFM directly, or compute it with an equation of their choice; the",
        "validation vignette scales FFM from weight using this paper's own",
        "21.0 kg / 26 kg anchor rather than importing an uncited equation.",
        "Exponents are fixed at 0.75 on CL and Q and 1 on Vc and Vp.",
        sep = " "
      ),
      source_name        = "FFM"
    ),
    CONMED_ATAZANAVIR = list(
      description        = "Concomitant atazanavir, i.e. the child is in the once-daily atazanavir/ritonavir arm",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (darunavir/ritonavir reference arm)",
      notes              = paste(
        "1 = atazanavir/ritonavir arm (N = 60), 0 = not on atazanavir. Time-fixed:",
        "the companion protease inhibitor was randomised at trial entry and stable",
        "through the week-6 intensive PK day. Carries TWO effects: +137% on relative",
        "bioavailability and +20.7% on clearance (Table 3). Together with",
        "CONMED_LOPINAVIR this encodes the paper's three-level companion-PI",
        "covariate; the darunavir arm (N = 59) is the reference and is represented",
        "by both indicators being 0, so no CONMED_DARUNAVIR column is needed.",
        "Note the direction of effect is the opposite of most existing",
        "CONMED_ATAZANAVIR models: here ritonavir is the analyte and atazanavir the",
        "perpetrator, whereas Arab-Alameddine 2012, von Hentig 2009 and Bukkems",
        "2021 all use the indicator on a different victim drug.",
        sep = " "
      ),
      source_name        = "PI_BCK_BONE == 2"
    ),
    CONMED_LOPINAVIR = list(
      description        = "Concomitant lopinavir, i.e. the child is in the twice-daily lopinavir/ritonavir arm",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (darunavir/ritonavir reference arm)",
      notes              = paste(
        "1 = lopinavir/ritonavir arm (N = 51), 0 = not on lopinavir. Time-fixed.",
        "Carries -23.4% on relative bioavailability (Table 3). This is also the",
        "only arm dosed twice daily, so it is the only arm in which the evening",
        "dose effect on the absorption lag is identified: the supplemental control",
        "stream gates that effect on 'PI_BCK_BONE.EQ.3.AND.OCC.EQ.1', i.e. the",
        "product CONMED_LOPINAVIR * (OCC == 1).",
        sep = " "
      ),
      source_name        = "PI_BCK_BONE == 3"
    ),
    OCC = list(
      description        = "Dosing-occasion indicator distinguishing the unwitnessed dose preceding the sampling window from the witnessed dose",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = paste(
        "Two occasions, per the supplemental control stream's OCC-gated ETA",
        "assignments. OCC = 1 is the last dose taken before the intensive PK",
        "sampling window, which was NOT taken under direct observation; in the",
        "twice-daily lopinavir/ritonavir arm this is the previous evening's dose.",
        "OCC = 2 is the dose given under direct observation on the PK day, in the",
        "morning with a 5% fat, ~250 kCal breakfast (Methods, 'Procedures').",
        "OCC is decomposed inside model() into binary indicators oc1 / oc2 that",
        "multiplex the per-occasion BOV etas, following the Chen 2023 nemonoxacin",
        "and Jonsson 2011 ethambutol precedent. OCC = 1 additionally (a) inflates",
        "the bioavailability BOV standard deviation 2.14-fold, because the dose was",
        "unwitnessed and both its amount and its timing are uncertain, and (b) in",
        "the lopinavir arm only, multiplies the absorption lag by 3.61.",
        "For a single witnessed morning dose matching the Table 3 reference",
        "condition, pass OCC = 2.",
        sep = " "
      ),
      source_name        = "OCC"
    )
  )

  compartmentData <- list(
    depot       = list(analyte = "ritonavir", units = "mg", specimen = "administration site", verified = TRUE),
    central     = list(analyte = "ritonavir", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "ritonavir", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 170,
    n_studies      = 1,
    age_range      = "3.16-15.6 years",
    age_median     = "10.5 years",
    weight_range   = "14.2-64.2 kg",
    weight_median  = "26.0 kg",
    height_range   = "97.0-169 cm",
    height_median  = "131 cm",
    sex_female_pct = 51.2,
    race_ethnicity = "not reported by category; all participants enrolled in Zambia, Uganda and Zimbabwe",
    disease_state  = "HIV-1 infection failing first-line antiretroviral therapy by WHO virological, CD4 or clinical criteria, starting second-line ritonavir-boosted protease-inhibitor ART",
    dose_range     = "ritonavir 50-200 mg total daily dose (median 100 mg; 1.56-6.90 mg/kg/day) by WHO weight band, as 200/50 mg lopinavir/ritonavir twice daily, 25 mg or 100 mg ritonavir or co-formulated 300/100 mg atazanavir/ritonavir once daily, or 100 mg ritonavir once daily with darunavir",
    regions        = "Zambia, Uganda, Zimbabwe",
    co_medication  = "two NRTIs: tenofovir alafenamide/emtricitabine (54.7%), abacavir/lamivudine (25.3%) or zidovudine/lamivudine (20.0%); no NRTI effect on ritonavir PK was found",
    notes          = paste(
      "Baseline characteristics in Table 1, stratified by boosted protease",
      "inhibitor arm (lopinavir N = 51, atazanavir N = 60, darunavir N = 59).",
      "Nested PK sub-study of the CHAPAS-4 trial; intensive sampling after week",
      "6 of study treatment at pre-dose, 0.5 h (TAF/FTC arms only), 1, 2, 4, 6,",
      "8, 12 and 24 h post-dose. 1,254 ritonavir concentrations, 6.9% below the",
      "0.045 mg/L LLOQ and 3.7% undetectable. Median weight-for-age Z-score",
      "-1.4 (-4.5 to 1.7) and height-for-age Z-score -1.3 (-4.4 to 3.6), so the",
      "cohort is substantially stunted and underweight relative to WHO",
      "references. One profile was excluded for implausibly low concentrations",
      "throughout, and all 24 h samples in the lopinavir arm were excluded",
      "because the 12 h dosing times were insufficiently documented.",
      sep = " "
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Structural disposition and absorption. Table 3 typical values apply
    # to the darunavir/ritonavir reference arm at FFM 21.0 kg (a 26 kg
    # child; Table 3 footnote b) on a witnessed morning dose (OCC = 2).
    #
    # Relative bioavailability is FIXED to 1 on that reference, so every
    # clearance and volume below is APPARENT (CL/F, Vc/F, Q/F, Vp/F).
    # ------------------------------------------------------------------
    lcl <- log(10.5);    label("Apparent oral clearance CL/F at FFM 21.0 kg, darunavir arm (L/h)")                       # Tsirizani 2025 Table 3 'Clearance (L/h) 10.5 (9.22-11.9)'
    lvc <- log(54.5);    label("Apparent central volume of distribution Vc/F at FFM 21.0 kg (L)")                        # Tsirizani 2025 Table 3 'Central volume of distribution (L) 54.5 (48.9-60.4)'
    lq  <- log(1.19);    label("Apparent intercompartmental clearance Q/F at FFM 21.0 kg (L/h)")                         # Tsirizani 2025 Table 3 'Intercompartmental clearance (L/h) 1.19 (0.851-1.58)'
    lvp <- log(134);     label("Apparent peripheral volume of distribution Vp/F at FFM 21.0 kg (L)")                     # Tsirizani 2025 Table 3 'Peripheral volume of distribution (L) 134 (38.8-330)'

    # Bioavailability anchor. Not estimable without a ritonavir-alone
    # control arm (Discussion, 'Our study had some strengths and
    # weaknesses'), so it is fixed to 1 on the darunavir reference and the
    # companion-PI effects below are read as RELATIVE bioavailability.
    lfdepot <- fixed(log(1)); label("Relative bioavailability on the darunavir reference arm (fraction)")               # Tsirizani 2025 Table 3 'Bioavailability 1 FIXED'; Data S1 '$THETA 1 FIX ; 4 BIO'

    ltlag <- log(0.981); label("Absorption lag time on a daytime dose (h)")                                             # Tsirizani 2025 Table 3 'Lag time (h) 0.981 (0.934-1.06)'
    ld1   <- log(2.45);  label("Duration of the zero-order absorption phase (h)")                                       # Tsirizani 2025 Table 3 'Zero-order absorption duration (h) 2.45 (2.05-2.80)'
    lka   <- log(1.16);  label("First-order absorption rate constant (1/h)")                                            # Tsirizani 2025 Table 3 'First-order absorption rate constant (1/h) 1.16 (0.893-1.89)'

    # ------------------------------------------------------------------
    # Allometric exponents, both FIXED at the theory-based values rather
    # than estimated.
    # ------------------------------------------------------------------
    e_ffm_cl <- fixed(0.75); label("Allometric exponent of fat-free mass on CL/F and Q/F (unitless)")                    # Tsirizani 2025 Methods 'Population pharmacokinetic analysis': "allometric scaling of clearance and volume parameters with a fixed exponent of 0.75 and 1, respectively"; Data S1 'ALLMCL_FFM = (FFM/TVFFM)**0.75'
    e_ffm_vc <- fixed(1);    label("Allometric exponent of fat-free mass on Vc/F and Vp/F (unitless)")                   # Tsirizani 2025 Methods, same sentence; Data S1 'ALLMV_FFM = (FFM/TVFFM)'

    # ------------------------------------------------------------------
    # Companion protease inhibitor effects, referenced to darunavir.
    # Table 3 reports them as percent changes; the supplemental control
    # stream carries them as multiplicative selectors, and the two agree:
    # 1 + 1.37 = 2.37 vs '$THETA (0,2.32224,10) ; 12 ATV_BIO';
    # 1 - 0.234 = 0.766 vs '$THETA (0,0.779626,5) ; 11 LPV_BIO';
    # 1 + 0.207 = 1.207 vs '$THETA (0,1.20875,5) ; 15 ATV_CL'
    # (the $THETA values are initial estimates, Table 3 the finals).
    # ------------------------------------------------------------------
    e_atazanavir_fdepot <- 1.37;   label("Fractional change in relative bioavailability with concomitant atazanavir (unitless)")   # Tsirizani 2025 Table 3 'Atazanavir on relative bioavailability (%) 137 (107-190)'
    e_lopinavir_fdepot  <- -0.234; label("Fractional change in relative bioavailability with concomitant lopinavir (unitless)")    # Tsirizani 2025 Table 3 'Lopinavir on relative bioavailability (%) -23.4 (-8.20 to -34.4)'
    e_atazanavir_cl     <- 0.207;  label("Fractional change in CL/F with concomitant atazanavir (unitless)")                      # Tsirizani 2025 Table 3 'Atazanavir on clearance (%) +20.7 (+11.3 to +31.3)'

    # Evening-dose effect on the absorption lag, identified only in the
    # twice-daily lopinavir arm and only on the unwitnessed OCC = 1 dose.
    # Applied as a power-of-binary multiplier so it collapses to 1 when
    # either indicator is 0.
    e_evening_tlag <- 3.61; label("Fold change in absorption lag time on the evening lopinavir/ritonavir dose (fold)")             # Tsirizani 2025 Table 3 'Night dose additional lag time (Fold) 3.61 (2.57-4.56)'; Data S1 'IF(PI_BCK_BONE.EQ.3.AND.OCC.EQ.1)LPV_NIGHT_LAG = THETA(10)'

    # Inflation of the bioavailability BOV standard deviation on the
    # unwitnessed dose. Data S1 applies it INSIDE the exponent, scaling the
    # eta itself: 'IF(OCC.EQ.1)BOVBIO = ETA(9)*EBOV'. Same encoding shape as
    # sd_ratio_cl_m5 in Keunecke_2020_regorafenib_phase3.R.
    sd_ratio_fdepot_occ1 <- 2.14; label("Ratio of the occasion-1 to occasion-2 bioavailability BOV standard deviation (unitless)") # Tsirizani 2025 Table 3 'Extra variability for unobserved doses (fold change) 2.14 (1.53-2.42)' with footnote c 'This parameter was on between-occasion variability in bioavailability'; Data S1 '$THETA (0,2.10111,5) ; 13 EBOV'

    # ------------------------------------------------------------------
    # Between-subject variability. Table 3 reports variability as a
    # percentage that equals the raw log-scale omega standard deviation
    # times 100, NOT the log-normal CV sqrt(exp(omega^2) - 1). The
    # supplemental $OMEGA initial estimates settle this: the BOV on the
    # first-order absorption rate constant is 114% in Table 3, and
    # sqrt(1.25373) = 1.120 while sqrt(exp(1.25373) - 1) = 1.582. Every
    # other row agrees on the same reading to within a few percent.
    #
    # Only CL carried BSV in the final model; the supplement fixes BSV on
    # Vc, ka, F, Vp and Q, and BOV on CL, to zero ('$OMEGA BLOCK(1) FIX 0'),
    # and Table 3 lists no such rows. Those five zero-variance etas are
    # omitted here rather than encoded as ~ fixed(0), which would make
    # OMEGA singular.
    # ------------------------------------------------------------------
    etalcl ~ 0.017161  # Tsirizani 2025 Table 3 'Between-subject variability / Clearance (%) 13.1 (10.2-15.8)' -> omega^2 = 0.131^2; Data S1 '$OMEGA BLOCK(1) 0.0154914 ; 1 BSVCL' (initial)

    # ------------------------------------------------------------------
    # Between-occasion variability. Each parameter has one variance shared
    # across both occasions, encoded as a free occasion-1 eta plus an
    # occasion-2 eta at ~ fixed(<same value>) to reproduce NONMEM's
    # '$OMEGA BLOCK(1) SAME'.
    # ------------------------------------------------------------------
    etaiov_fdepot_1 ~ 0.160801         # Tsirizani 2025 Table 3 'Between-occasion variability / Bioavailability (%) 40.1 (35.2-46.0)' -> 0.401^2; Data S1 '$OMEGA BLOCK(1) 0.159089 ; 9 BOVBIO' (initial)
    etaiov_fdepot_2 ~ fixed(0.160801)  # Data S1 '$OMEGA BLOCK(1) SAME' following BOVBIO
    etaiov_tlag_1   ~ 0.384400         # Tsirizani 2025 Table 3 'Lag time (%) 62.0 (54.5-70.2)' -> 0.620^2; Data S1 '$OMEGA BLOCK(1) 0.372121 ; 13 BOVLAG' (initial)
    etaiov_tlag_2   ~ fixed(0.384400)  # Data S1 '$OMEGA BLOCK(1) SAME' following BOVLAG
    etaiov_d1_1     ~ 0.275625         # Tsirizani 2025 Table 3 'Zero-order rate of absorption (%) 52.5 (42.7-63.0)' -> 0.525^2; Data S1 '$OMEGA BLOCK(1) 0.289698 ; 15 BOVD1' (initial)
    etaiov_d1_2     ~ fixed(0.275625)  # Data S1 '$OMEGA BLOCK(1) SAME' following BOVD1
    etaiov_ka_1     ~ 1.299600         # Tsirizani 2025 Table 3 'First-order absorption rate constant (%) 114 (90.0-141)' -> 1.14^2; Data S1 '$OMEGA BLOCK(1) 1.25373 ; 11 BOVKA' (initial)
    etaiov_ka_2     ~ fixed(1.299600)  # Data S1 '$OMEGA BLOCK(1) SAME' following BOVKA

    # ------------------------------------------------------------------
    # Combined additive-and-proportional residual error, per Data S1
    # '$ERROR': W = SQRT(ADD**2 + PROP**2), Y = IPRED + W*ERR(1) with
    # '$SIGMA 1 FIX'.
    # ------------------------------------------------------------------
    propSd <- 0.147;  label("Proportional residual error (fraction)")     # Tsirizani 2025 Table 3 'Proportional error (%) 14.7 (13.4-16.1)'; Data S1 'PROP = IPRED*THETA(5)'
    addSd  <- 0.0152; label("Additive residual error (mg/L)")             # Tsirizani 2025 Table 3 'Additive error (mg/L) 0.0152 (0.0121-0.0184)'; Data S1 'ADD = THETA(6)+(LLOQ*0.2)' with LLOQ = 0.045, so this is the TOTAL additive SD at an uncensored record -- see the file header note
  })

  model({
    # ------------------------------------------------------------------
    # 1. Occasion indicators and the per-occasion BOV terms.
    #
    #    OCC = 1 is the unwitnessed dose preceding the sampling window
    #    (the previous evening's dose in the twice-daily lopinavir arm);
    #    OCC = 2 is the witnessed morning dose taken with breakfast.
    #
    #    The occasion-1 bioavailability eta is scaled by
    #    sd_ratio_fdepot_occ1 INSIDE the exponent, matching Data S1's
    #    'IF(OCC.EQ.1)BOVBIO = ETA(9)*EBOV', so the effective occasion-1
    #    BOV standard deviation is 0.401 * 2.14 = 0.858 on the log scale.
    # ------------------------------------------------------------------
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)

    iov_fdepot <- oc1 * sd_ratio_fdepot_occ1 * etaiov_fdepot_1 + oc2 * etaiov_fdepot_2
    iov_tlag   <- oc1 * etaiov_tlag_1 + oc2 * etaiov_tlag_2
    iov_d1     <- oc1 * etaiov_d1_1   + oc2 * etaiov_d1_2
    iov_ka     <- oc1 * etaiov_ka_1   + oc2 * etaiov_ka_2

    # ------------------------------------------------------------------
    # 2. Derived covariate multipliers. All collapse to 1 at the Table 3
    #    reference covariate vector: FFM 21.0 kg, darunavir arm
    #    (CONMED_ATAZANAVIR = CONMED_LOPINAVIR = 0), witnessed morning
    #    dose (OCC = 2).
    # ------------------------------------------------------------------
    ffm_cl <- (FFM / 21.0)^e_ffm_cl
    ffm_v  <- (FFM / 21.0)^e_ffm_vc

    # Evening-dose lag multiplier: 3.61 only when the child is in the
    # twice-daily lopinavir arm AND the dose is the unwitnessed OCC = 1
    # (evening) dose; 1 otherwise.
    evening_tlag <- e_evening_tlag^(CONMED_LOPINAVIR * oc1)

    # ------------------------------------------------------------------
    # 3. Individual parameters.
    # ------------------------------------------------------------------
    cl <- exp(lcl + etalcl) * ffm_cl * (1 + e_atazanavir_cl * CONMED_ATAZANAVIR)
    vc <- exp(lvc) * ffm_v
    q  <- exp(lq)  * ffm_cl
    vp <- exp(lvp) * ffm_v

    ka     <- exp(lka + iov_ka)
    d1     <- exp(ld1 + iov_d1)
    tlag   <- exp(ltlag + iov_tlag) * evening_tlag
    fdepot <- exp(lfdepot + iov_fdepot) *
      (1 + e_atazanavir_fdepot * CONMED_ATAZANAVIR) *
      (1 + e_lopinavir_fdepot  * CONMED_LOPINAVIR)

    # ------------------------------------------------------------------
    # 4. Micro-constants for the two-compartment system (Data S1 uses
    #    ADVAN4 TRANS1 with K = CL/V, K23 = Q/V, K32 = Q/V3).
    # ------------------------------------------------------------------
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ------------------------------------------------------------------
    # 5. Two-compartment disposition with an absorption lag followed by
    #    sequential zero-order then first-order absorption (Results,
    #    'Population pharmacokinetic analysis': "a lag time in absorption
    #    followed by sequential zero- and first-order absorption",
    #    dOFV = -488 against no absorption delay). The dose enters `depot`
    #    at a constant rate over the window of duration d1 beginning tlag
    #    after the dose record, and `depot` then drains into `central`
    #    first-order at ka.
    #
    #    Because d1 is a MODELLED duration, dose records must carry
    #    rate = -2 so rxode2 honours dur(depot); a plain bolus collapses
    #    the zero-order phase and biases Cmax upward.
    # ------------------------------------------------------------------
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1

    f(depot)    <- fdepot
    dur(depot)  <- d1
    alag(depot) <- tlag

    # ------------------------------------------------------------------
    # 6. Observation. Dose is in mg and vc in L, so central / vc is mg/L,
    #    the unit of Table 3's additive error (0.0152 mg/L), of the
    #    0.045 mg/L LLOQ, and of the Table 4 Cmax values.
    # ------------------------------------------------------------------
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
