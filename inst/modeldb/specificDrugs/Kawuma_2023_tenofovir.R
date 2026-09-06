Kawuma_2023_tenofovir <- function() {
  description <- "Joint semimechanistic two-compartment population PK model for plasma tenofovir when administered as either tenofovir disoproxil fumarate (TDF) or tenofovir alafenamide (TAF) in South African adults living with HIV (Kawuma 2023). Both prodrugs share one tenofovir disposition model; TDF appears in plasma by a single first-order process, whereas a TAF dose splits into a fast fraction absorbed first-order into the systemic circulation and a slow fraction sequestered in an intracellular reservoir (interpreted as PBMCs) that releases tenofovir with a fixed 6.83-day half-life. Doses are expressed as tenofovir-equivalent amounts."
  reference <- "Kawuma AN, Wasmann RE, Sinxadi P, Sokhela SM, Chandiwana N, Venter WDF, Wiesner L, Maartens G, Denti P. Population pharmacokinetics of tenofovir given as either tenofovir disoproxil fumarate or tenofovir alafenamide in an African population. CPT Pharmacometrics Syst Pharmacol. 2023;12(6):821-830. doi:10.1002/psp4.12955"
  vignette <- "Kawuma_2023_tenofovir"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # The TAF dose enters two separate input compartments (depot2, depot3) and the
  # TDF dose a third (depot), so buildModelDb()'s two-name depot/central
  # heuristic cannot infer the dosing targets and would otherwise report
  # `dosing = central`.
  dosing <- c("depot", "depot2", "depot3")

  covariateData <- list(
    WT = list(
      description        = "Total body weight, used for allometric scaling of all clearance and volume parameters with exponents fixed to 0.75 and 1 relative to a 70 kg reference individual.",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Kawuma 2023 Methods: 'Allometry with either total body weight or fat-free mass was tested in the model and allometric exponents for clearance and volume were fixed to 0.75 and 1, respectively.' Results: allometric scaling with weight was retained (dOFV = -18 vs no allometry) and applied to CL, Q, Vc and Vp; fat-free mass did not improve the fit (dOFV = -3.80). Table 2 footnote c gives the 70 kg reference. Cohort median (IQR) weight 73.1 (67.2-85.2) kg (Table 1).",
      source_name        = "WT"
    ),
    OCC = list(
      description        = "Integer occasion index for between-occasion variability on the absorption parameters and on TDF bioavailability. 1 = the unobserved dose taken the day before the pharmacokinetic visit (which produces the predose sample); 2 = the observed dose taken on the day of the pharmacokinetic visit (with samples at 1, 2, 4, 6, 8 and 24 h).",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = "Kawuma 2023 Methods: \"an 'occasion' was defined as a dose with its proceeding sample. For example, the unobserved dose the day before the pharmacokinetic visit (leading to the predose sample) was considered a separate occasion from the observed dose on the day of the pharmacokinetic visit\". The rich-sampling design therefore has exactly two occasions, encoded here as OCC = 1 and OCC = 2 with a shared variance (the NONMEM $OMEGA BLOCK(1) SAME idiom). Records outside those two occasions should carry OCC = 0, which zeroes every indicator and removes the between-occasion contribution; Kawuma 2023 Table 3 footnote a states the published AUC simulations were themselves run 'excluding between-occasion variability on absorption parameters'.",
      source_name        = "OCC"
    )
  )

  # Covariates screened by Kawuma 2023 but NOT retained in the final model:
  # "Adding creatinine clearance (collected at screening) as a covariate on
  # clearance, did not improve the model significantly and neither did age."
  # No point estimate is published for either effect, so neither can be encoded.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age at screening.",
      units       = "years",
      type        = "continuous",
      notes       = "Tested on tenofovir pharmacokinetic parameters by stepwise inclusion (dOFV > 3.84) and not retained; Kawuma 2023 Results. Cohort median (IQR) 31.0 (29.0-36.0) years, Table 1. No coefficient is reported."
    ),
    CRCL_BASE = list(
      description = "Baseline creatinine clearance at screening, calculated by Cockcroft and Gault (not BSA-normalized).",
      units       = "mL/min",
      type        = "continuous",
      notes       = "Tested as a covariate on tenofovir clearance and not retained; Kawuma 2023 Results and Discussion attribute the negative finding to the narrow distribution in the cohort (median 120, IQR 96.0-140 mL/min, Table 1) and to creatinine not being measured at the time of drug sampling. No coefficient is reported."
    )
  )

  compartmentData <- list(
    depot = list(
      analyte  = "tenofovir",
      units    = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    depot2 = list(
      analyte  = "tenofovir",
      units    = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    depot3 = list(
      analyte  = "tenofovir",
      units    = "mg",
      specimen = "blood cell",
      verified = FALSE
    ),
    central = list(
      analyte  = "tenofovir",
      units    = "mg",
      specimen = "plasma",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte  = "tenofovir",
      units    = "mg",
      specimen = "plasma",
      verified = TRUE
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 41,
    n_studies      = 1,
    n_observations = 279,
    age_median     = "31.0 years (IQR 29.0-36.0)",
    weight_median  = "73.1 kg (IQR 67.2-85.2)",
    height_median  = "167 cm (IQR 161-174)",
    sex_female_pct = 65.9,
    race_ethnicity = "South African adults; the source paper reports no race or ethnicity breakdown.",
    renal_function = "Creatinine clearance at screening (Cockcroft-Gault) median 120 mL/min (IQR 96.0-140); no participant had impaired renal function.",
    disease_state  = "Treatment-naive adults living with HIV, sampled after at least 48 weeks of antiretroviral therapy.",
    dose_range     = "TDF 300 mg once daily (= 136 mg tenofovir) in 21 participants, or TAF 25 mg once daily (= 15 mg tenofovir) in 20 participants; both co-formulated with emtricitabine 200 mg and given with dolutegravir 50 mg.",
    co_medication  = "Dolutegravir 50 mg and emtricitabine 200 mg once daily in both arms.",
    regions        = "South Africa",
    notes          = "Pharmacokinetic substudy nested within the ADVANCE trial (NCT03122262), an open-label phase III randomized noninferiority trial. Rich sampling predose and at 1, 2, 4, 6, 8 and 24 h postdose. Baseline demographics are in Kawuma 2023 Table 1. None of the 279 samples was below the 0.0005 mg/L limit of quantification."
  )

  ini({
    # =========================================================================
    # Tenofovir disposition -- shared by both prodrugs. Kawuma 2023 Table 2,
    # typical values for a 70 kg individual (Table 2 footnote c).
    # =========================================================================
    lcl <- log(44.7)
    label("Tenofovir clearance for a 70 kg individual (CL, L/h)")            # Table 2: CL 44.7 (40.2-49.5) L/h
    lvc <- log(378)
    label("Tenofovir central volume of distribution for a 70 kg individual (Vc, L)")   # Table 2: Vc 378 (319-459) L
    lq <- log(157)
    label("Tenofovir intercompartmental clearance for a 70 kg individual (Q, L/h)")    # Table 2: Q 157 (103-233) L/h
    lvp <- log(356)
    label("Tenofovir peripheral volume of distribution for a 70 kg individual (Vp, L)")  # Table 2: Vp 356 (298-438) L

    # Allometric exponents, fixed at the canonical values (Methods: "allometric
    # exponents for clearance and volume were fixed to 0.75 and 1,
    # respectively"). Applied to all clearance and all volume parameters
    # (Results: "applied to all clearance and volume parameters").
    e_wt_cl_q <- fixed(0.75)
    label("Allometric exponent on WT/70 shared by CL and Q (unitless)")      # Methods, Table 2 footnote c
    e_wt_vc_vp <- fixed(1)
    label("Allometric exponent on WT/70 shared by Vc and Vp (unitless)")     # Methods, Table 2 footnote c

    # =========================================================================
    # TDF absorption -- a single first-order process into the systemic
    # circulation. F_TDF is the structural bioavailability anchor.
    # =========================================================================
    lka_tdf <- log(3.04)
    label("First-order rate constant for appearance of tenofovir in plasma after a TDF dose (Ka_TDF, 1/h)")  # Table 2: Ka_TDF 3.04 (2.11-3.88) 1/h
    lfdepot_tdf <- fixed(log(1))
    label("Bioavailability of tenofovir given as TDF (F_TDF, fraction); the reference against which F_TAF is relative")  # Table 2: F_TDF 1-Fixed

    # =========================================================================
    # TAF absorption -- two parallel pathways (Kawuma 2023 Figure 2 legend):
    # a fast fraction Frac_TAF-Fast absorbed first-order into the systemic
    # circulation, and the remaining slow fraction (1 - Frac_TAF-Fast) first
    # taken up into an intracellular reservoir and then released to plasma by a
    # first-order process with half-life t1/2_TAF-Slow.
    # =========================================================================
    lfdepot_taf <- log(0.822)
    label("Bioavailability of tenofovir given as TAF, relative to TDF (F_TAF, fraction)")  # Table 2: F_TAF 0.822 (0.723-0.939); Results: "relative bioavailability of tenofovir when given as TAF, was estimated to be 82.2% (95% CI, 72.3-93.9)"
    logitffo <- log(0.324 / (1 - 0.324))
    label("Logit of the fraction of the bioavailable TAF dose entering the fast (first-order) pathway (Frac_TAF-Fast, fraction)")  # Table 2: Frac_TAF-Fast 32.4% (27.0-37.7)
    lka_taf_fast <- log(1.45)
    label("First-order rate constant for the fast TAF pathway (Ka_TAF, 1/h)")  # Table 2: Ka_TAF 1.45 (0.924-2.60) 1/h
    # The slow pathway is reported as a terminal half-life in DAYS; the
    # equivalent first-order release rate constant in 1/h is log(2) / (6.83 d *
    # 24 h/d) = 0.004229 1/h. Fixed by the authors after a likelihood-profiling
    # exercise showed values of 5-60 days were nearly equally supported; 6.83
    # days is the previously reported intracellular TFV-DP decay half-life
    # (Kawuma 2023 Results, citing reference 16).
    lka_taf_slow <- fixed(log(log(2) / (6.83 * 24)))
    label("First-order rate constant for release of tenofovir from the intracellular reservoir to plasma (1/h)")  # Table 2: t1/2_TAF-Slow 6.83 days-Fixed

    # =========================================================================
    # Random effects. Kawuma 2023 Table 2 footnote b: "Between-subject
    # variability (BSV) and between-occasion variability (BOV) were assumed to
    # be log-normally distributed and calculated by CV% = sqrt(omega^2) * 100".
    # The reported percentages are therefore omega on the SD scale directly, so
    # omega^2 = (CV%/100)^2 -- NOT the exp-based log(1 + CV^2) conversion.
    # =========================================================================
    etalcl ~ 0.040401      # Table 2 CL row, footnote d (BSV) 20.1% (16.1-24.7); omega^2 = 0.201^2

    # Between-occasion variability across the two occasions defined in Methods.
    # rxode2 has no NONMEM-style `| occasion` level, so each occasion gets its
    # own eta and the indicators in model() select among them. The second
    # occasion's variance is fixed to the first, which is the $OMEGA BLOCK(1)
    # SAME idiom: one magnitude shared by every occasion.
    etaiov_fdepot_tdf_1 ~ 0.057121        # Table 2 F_TDF row, footnote e (BOV) 23.9% (18.2-30.3); omega^2 = 0.239^2
    etaiov_fdepot_tdf_2 ~ fixed(0.057121) # shared magnitude, $OMEGA BLOCK(1) SAME equivalent
    etaiov_ka_tdf_1 ~ 1.311025            # Table 2 Ka_TDF row, footnote e (BOV) 114.5% (68.4-162); omega^2 = 1.145^2
    etaiov_ka_tdf_2 ~ fixed(1.311025)     # shared magnitude, $OMEGA BLOCK(1) SAME equivalent
    etaiov_ka_taf_1 ~ 0.439569            # Table 2 Ka_TAF row, footnote e (BOV) 66.3% (31.0-91.6); omega^2 = 0.663^2
    etaiov_ka_taf_2 ~ fixed(0.439569)     # shared magnitude, $OMEGA BLOCK(1) SAME equivalent

    # No variability on Vc, Q, Vp, F_TAF or Frac_TAF-Fast: Table 2 leaves those
    # variability cells blank, and Results states "Within the TAF arm, the model
    # did not support the estimation of any variability (BOV or BSV) on the
    # bioavailability parameter."

    # =========================================================================
    # Residual error (Kawuma 2023 Table 2). The additive component was
    # constrained to be at least 20% of the LLOQ (Methods) and is reported as
    # fixed at that bound; the LLOQ was 0.0005 mg/L (Table 2 footnote f), so
    # the additive SD is 0.2 * 0.0005 = 0.0001 mg/L.
    # =========================================================================
    propSd <- 0.119
    label("Proportional residual error (fraction)")   # Table 2: Proportional error 11.9% (10.8-13.4)
    addSd <- fixed(0.0001)
    label("Additive residual error (mg/L)")           # Table 2: Additive error = 20% of LLOQ-Fixed, LLOQ = 0.0005 mg/L
  })

  model({
    # ---- 1. Occasion indicators and between-occasion variability ------------
    # OCC = 0 (or any value outside 1:2) zeroes every indicator, which
    # reproduces the paper's own simulation setting of "excluding
    # between-occasion variability on absorption parameters" (Table 3 note a).
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    iov_fdepot_tdf <- oc1 * etaiov_fdepot_tdf_1 + oc2 * etaiov_fdepot_tdf_2
    iov_ka_tdf <- oc1 * etaiov_ka_tdf_1 + oc2 * etaiov_ka_tdf_2
    iov_ka_taf <- oc1 * etaiov_ka_taf_1 + oc2 * etaiov_ka_taf_2

    # ---- 2. Allometric size terms, reference 70 kg -------------------------
    wt_cl_q <- (WT / 70)^e_wt_cl_q
    wt_vc_vp <- (WT / 70)^e_wt_vc_vp

    # ---- 3. Individual parameters ------------------------------------------
    cl <- exp(lcl + etalcl) * wt_cl_q
    vc <- exp(lvc) * wt_vc_vp
    q <- exp(lq) * wt_cl_q
    vp <- exp(lvp) * wt_vc_vp

    ka_tdf <- exp(lka_tdf + iov_ka_tdf)
    ka_taf_fast <- exp(lka_taf_fast + iov_ka_taf)
    ka_taf_slow <- exp(lka_taf_slow)

    fdepot_tdf <- exp(lfdepot_tdf + iov_fdepot_tdf)
    fdepot_taf <- exp(lfdepot_taf)
    ffo <- expit(logitffo)

    # ---- 4. Micro-constants -------------------------------------------------
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # ---- 5. ODE system ------------------------------------------------------
    # depot  : tenofovir appearing in plasma from a TDF dose (first-order).
    # depot2 : the fast TAF pathway -- tenofovir immediately available for
    #          absorption into the systemic circulation.
    # depot3 : the slow TAF pathway -- tenofovir first sequestered
    #          intracellularly (interpreted by the authors as the PBMC pool,
    #          where it is phosphorylated to TFV-DP) and then released back to
    #          plasma once the diphosphate decays.
    d/dt(depot) <- -ka_tdf * depot
    d/dt(depot2) <- -ka_taf_fast * depot2
    d/dt(depot3) <- -ka_taf_slow * depot3
    d/dt(central) <- ka_tdf * depot + ka_taf_fast * depot2 + ka_taf_slow * depot3 -
      kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # ---- 6. Bioavailability -------------------------------------------------
    # A TDF regimen doses `depot` only. A TAF regimen doses BOTH `depot2` and
    # `depot3` with the same tenofovir-equivalent amount; the f() terms split
    # that amount between the fast and slow pathways so the two records
    # together deliver F_TAF of one dose. Doses are tenofovir-equivalent
    # amounts (Methods: 300 mg TDF provides 136 mg tenofovir; 25 mg TAF
    # provides 15 mg tenofovir).
    f(depot) <- fdepot_tdf
    f(depot2) <- fdepot_taf * ffo
    f(depot3) <- fdepot_taf * (1 - ffo)

    # ---- 7. Observation and residual error ----------------------------------
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
