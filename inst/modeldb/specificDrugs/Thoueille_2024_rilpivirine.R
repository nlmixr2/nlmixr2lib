Thoueille_2024_rilpivirine <- function() {
  description <- "Two-compartment population PK model for rilpivirine covering both the oral 25 mg daily lead-in and the long-acting intramuscular nanosuspension in people with HIV followed in routine clinical care in the Swiss HIV Cohort Study (Thoueille 2024). Oral doses enter the central compartment as a zero-order input of duration Doral fixed to 4 h and carry a relative bioavailability Foral = 65.4% with between-subject variability; the intramuscular injection is assumed completely bioavailable and is split between two parallel first-order release pathways, a fast pathway taking a fraction Fi.m.fast = 27.6% of the dose with kafast = 0.00214 1/h and a slow pathway taking the remainder with kaslow = 0.000229 1/h. Because both absorption rate constants are far below the elimination rate constant, long-acting rilpivirine shows flip-flop kinetics with an apparent half-life of 18 weeks driven by kaslow. Clearance carries both between-subject and inter-occasion variability across up to six injection occasions, and the residual error switches by route: additive after oral dosing and proportional after intramuscular injection. No covariate was retained in the final validated model; a female-sex effect reducing Fi.m.fast by 45.6% was estimated but judged not clinically relevant and left out (see covariatesDataExcluded)."
  reference <- paste(
    "Thoueille P, Saldanha SA, Schaller F, Choong E, Veuve F, Munting A,",
    "Cavassini M, Braun D, Gunthard HF, Duran Ramirez JJ, Surial B, Furrer H,",
    "Rauch A, Ustero P, Calmy A, Stockle M, Di Benedetto C, Bernasconi E,",
    "Schmid P, Marzolini C, Girardin FR, Buclin T, Decosterd LA, Guidi M;",
    "for the Swiss HIV Cohort Study.",
    "Population pharmacokinetics of rilpivirine following oral administration",
    "and long-acting intramuscular injection in real-world people with HIV.",
    "Front Pharmacol. 2024;15:1437400. doi:10.3389/fphar.2024.1437400.",
    sep = " "
  )
  vignette <- "Thoueille_2024_rilpivirine"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Three dose targets: an intramuscular injection is entered as simultaneous
  # records into `depot` (fast release) and `depot2` (slow release), and an oral
  # dose is a zero-order input into `central`. buildModelDb()'s dosing heuristic
  # only looks for compartments literally named `depot` and `central`, so
  # without this field the registry would report `depot,central` -- true but one
  # short. Unrelated to units$dosing above, which is a unit spelling.
  dosing <- c("depot", "depot2", "central")

  # Residual-error SDs are route-specific (Thoueille 2024 Supplementary
  # NONMEM code, $ERROR: Y0 = IPRED0 + ERR(1) for oral records and
  # Y1 = IPRED1 * (1 + ERR(2)) for intramuscular records), so neither
  # carries the bare canonical `addSd` / `propSd` name.
  paper_specific_residual_sds <- c("addSdOral", "propSdIm")

  covariateData <- list(
    ROUTE_ORAL = list(
      description        = "1 = the record belongs to the oral rilpivirine phase (25 mg tablet once daily), 0 = the record belongs to the long-acting intramuscular nanosuspension phase",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (long-acting intramuscular gluteal injection; there is no third route in this analysis)",
      notes              = "Row-level indicator carried on BOTH dose and observation records. Thoueille 2024 supplementary NONMEM control stream ships this column as LAI with the opposite polarity (LAI = 1 for the intramuscular records, LAI = 0 for the oral records), so ROUTE_ORAL = 1 - LAI. In the source model the column does two jobs. (i) In $PK it selects the absorption block: IF(LAI.EQ.0) the zero-order oral duration D3 = THETA(7) applies, ELSE the two intramuscular first-order rate constants KA1 and KA2 apply. That selection is structural rather than a covariate effect and is reproduced here by the dose record's target compartment instead -- an oral dose is placed in `central` as a zero-order input and never reads the depot parameters, while an intramuscular dose is placed in `depot` and `depot2` and never reads the oral parameters -- so the model body does not need ROUTE_ORAL for it. (ii) In $ERROR it selects the residual-error structure, which cannot be expressed by the target compartment, so ROUTE_ORAL is read in model() to switch between the additive oral SD and the proportional intramuscular SD. Set ROUTE_ORAL = 1 on every oral observation and 0 on every intramuscular observation.",
      source_name        = "LAI"
    ),
    OCC = list(
      description        = "Integer injection-occasion index used for the inter-occasion variability on clearance",
      units              = "(count)",
      type               = "categorical",
      reference_category = "n/a -- decomposed inside model() into six mutually exclusive binary indicators multiplied against the per-occasion etaiov_cl_<k> slots",
      notes              = "Thoueille 2024 Methods 'Model building and selection': 'occasions were coded to be consistent with the duration of the follow-up by including an occasion variable constructed with an incremental number within each subject for a maximum of six occasions/injections'. The supplementary control stream implements exactly six slots via $ABBR REPLACE ETA(OCC)=ETA(5,6,7,8,9,10), backed by $OMEGA BLOCK(1) 0.0167 followed by five BLOCK(1) SAME lines, so all six occasions share one variance. Values run 1 to 6; records outside the six injection occasions (in particular the oral lead-in records, which precede the first injection) take OCC = 1 so that exactly one indicator is active on every row. Thoueille 2024 Results states that the inter-occasion variability could be supported on clearance but not on the intramuscular absorption parameters.",
      source_name        = "OCC"
    )
  )

  covariatesDataExcluded <- list(
    SEXF = list(
      description = "Sex at birth indicator, 1 = female, 0 = male",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened, estimated, and then deliberately NOT retained. Thoueille 2024 Results 3.1 reports that univariate analysis found an effect of female sex on the fast intramuscular absorption fraction (dOFV = -11, p < 0.001) which survived forward insertion and backward deletion, entering the logit as TEMP = ln(theta_Fi.m.fast * (1 + theta_Female) / (1 - theta_Fi.m.fast * (1 + theta_Female))) with theta_Female = -0.456, i.e. females had an Fi.m.fast 45.6% lower than males. The authors nevertheless excluded it: Results 3.2 states 'although statistically significant, the effect of sex on long-acting rilpivirine Ctrough was not considered clinically relevant, and this model was not validated', because the sex covariate explained only 11% of the between-subject variability on Fi.m.fast and moved trough concentrations by no more than 15% over 48 weeks. Table 2 ('Final population PK parameter estimates') contains no theta_Female row, and in the supplementary control stream both the covariate TEMP line and the -0.456 THETA are commented out, which is what fixes the final model as the covariate-free one encoded here. A user who wants the sex model can reinstate it by replacing logitfdepot with log(0.276 * (1 - 0.456 * SEXF) / (1 - 0.276 * (1 - 0.456 * SEXF))); the vignette shows the resulting trough comparison against Thoueille 2024 Supplementary Table S2."
    ),
    AGE = list(
      description = "Age at the last recorded value",
      units       = "years",
      type        = "continuous",
      notes       = "Screened with a linear function on the base-model parameters (Thoueille 2024 Methods 'Model building and selection') and not retained. Cohort median 46 years, range 20-79 (Table 1). Thoueille 2024 Discussion notes that the small number of older participants may have masked an age effect reported by a PBPK analysis."
    ),
    WT = list(
      description = "Body weight at the last recorded value",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened both as a linear function and as an allometric scaling relationship with the exponent estimated (Thoueille 2024 Methods 'Model building and selection') and not retained. Cohort median 78 kg, range 50-126 (Table 1)."
    ),
    BMI = list(
      description = "Body mass index at the last recorded value",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Screened as a linear function and as an allometric scaling relationship and reached univariate significance on the fast intramuscular absorption fraction (dOFV = -6, p < 0.05; Thoueille 2024 Results 3.1) but was dropped at backward deletion (p < 0.01) and does not appear in the final model. Cohort median 25.4 kg/m^2, range 18.2-43.3 (Table 1). Thoueille 2024 Discussion notes that no morbidly obese participants were enrolled, so a reported PBPK obesity effect could not be tested."
    ),
    RACE_BLACK = list(
      description = "Self-reported Black ethnicity indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Ethnicity was screened as a covariate (Thoueille 2024 Methods 'Model building and selection') and not retained. Table 1 reports White 133 (56%), Black 36 (15%), Hispanic American 19 (8%), Asian 11 (5%), Other/Missing 39 (16%). The paper reports no point estimate for any ethnicity contrast, so no coefficient can be recovered."
    ),
    EGFR = list(
      description = "Estimated glomerular filtration rate, CKD-EPI (Levey 2009)",
      units       = "mL/min/1.73m^2",
      type        = "continuous",
      notes       = "Screened as CKD-EPI eGFR categories (Thoueille 2024 Methods 'Model building and selection') and not retained. Table 1 reports G1 (>= 90) 158 (66%), G2 (60-89) 76 (32%), G3 (30-59) 4 (2%). Rilpivirine is cleared hepatically, so a renal effect was not expected."
    )
  )

  compartmentData <- list(
    depot        = list(analyte = "rilpivirine", units = "mg", specimen = "administration site", verified = TRUE),
    depot2       = list(analyte = "rilpivirine", units = "mg", specimen = "administration site", verified = TRUE),
    central      = list(analyte = "rilpivirine", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1  = list(analyte = "rilpivirine", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 238L,
    n_studies      = 1L,
    n_observations = "1038 rilpivirine plasma concentrations: 186 after oral administration from 176 people and 852 after intramuscular injection from 222 people, with detailed within-interval sampling in 28 people (Thoueille 2024 Results, first paragraph). Median 4 samples per person (range 1-15); median follow-up 26 weeks (range 3-196). Only 10 people had concentrations assumed to be at intramuscular steady state (from week 96).",
    age_range      = "20-79 years (median 46; Thoueille 2024 Table 1)",
    age_median     = "46 years",
    weight_range   = "50-126 kg (median 78; Thoueille 2024 Table 1)",
    weight_median  = "78 kg",
    height_range   = "151-198 cm (median 176; Thoueille 2024 Table 1)",
    bmi_range      = "18.2-43.3 kg/m^2 (median 25.4; Thoueille 2024 Table 1). BMI < 25 in 104 (44%), 25-30 in 103 (43%), > 30 in 31 (13%). No morbidly obese participants were enrolled.",
    sex_female_pct = 20.2,
    race_ethnicity = c(White = 56, Black = 15, `Hispanic American` = 8, Asian = 5, `Other/Missing` = 16),
    disease_state  = "People with HIV-1 on suppressive antiretroviral therapy switching to, or established on, long-acting cabotegravir plus rilpivirine. Plasma HIV RNA < 50 copies/mL in 233 (98%); CD4 >= 500 cells/mm^3 in 186 (78%). Liver cirrhosis in 2 (1%, both Child-Pugh class A).",
    renal_function = "CKD-EPI eGFR category G1 (>= 90 mL/min/1.73m^2) in 158 (66%), G2 (60-89) in 76 (32%), G3 (30-59) in 4 (2%)",
    dose_range     = "Oral rilpivirine 25 mg once daily during the lead-in; long-acting intramuscular gluteal rilpivirine 900 mg with cabotegravir 600 mg every 2 months, and 600 mg with cabotegravir 400 mg every 4 weeks in two people treated for compassionate use before Swiss market authorisation",
    regions        = "Switzerland (Lausanne, Zurich, Bern, Geneva, Basel, Lugano, St Gallen)",
    co_medication  = "Cabotegravir is co-administered in every intramuscular injection. Thoueille 2024 Methods state that no clinically relevant interacting comedication, such as a potent CYP3A4 inducer, was encountered in the study population.",
    notes          = "Real-world therapeutic-drug-monitoring cohort nested in the Swiss HIV Cohort Study, sampled mostly sparsely at the discretion of physicians between March 2022 and June 2023, with a richer within-interval sampling substudy (pre-dose and 1, 2, 4 and 8 weeks after injection) offered to consenting participants in Lausanne and Geneva. This is an observational cohort rather than a registrational trial, which is the stated contrast with the phase III popPK analyses of Neyens 2021 and Benaboud 2023."
  )

  ini({
    # =====================================================================
    # All point estimates below are from Thoueille 2024 Table 2 ("Final
    # population PK parameter estimates of rilpivirine with their bootstrap
    # evaluations"), cross-checked against the $THETA / $OMEGA / $SIGMA
    # blocks of the supplementary NONMEM control stream ("NONMEM CODE FOR
    # FINAL MODEL", Supplementary Information p. 7-9). Where Table 2
    # reports a rounded CV% and the control stream reports the underlying
    # variance, the control-stream variance is used because it is the
    # unrounded quantity that was actually estimated.
    #
    # NONMEM compartment map (control stream $MODEL):
    #   COMP(DEPOT1)  -> depot        (fast intramuscular release)
    #   COMP(DEPOT2)  -> depot2       (slow intramuscular release)
    #   COMP(CENTRAL) -> central
    #   COMP(PERIPH)  -> peripheral1
    # =====================================================================

    # --- Disposition -----------------------------------------------------
    lcl <- log(6.74)  ; label("Apparent clearance (CL, L/h)")                                 # Table 2 'CL (L/h) 6.74 (RSE 3)', bootstrap median 6.68 [5.41-7.37]; control stream $THETA(1) 6.74
    lvc <- log(277)   ; label("Apparent central volume of distribution (V3, L)")              # Table 2 'V3 (L) 277 (RSE 25)', bootstrap median 274 [184-433]; control stream $THETA(2) 277
    lq  <- log(4.08)  ; label("Apparent inter-compartmental clearance (Q, L/h)")              # Table 2 'Q (L/h) 4.08 (RSE 40)', bootstrap median 4.03 [1.75-9.30]; control stream $THETA(3) 4.08
    lvp <- log(839)   ; label("Apparent peripheral volume of distribution (V4, L)")           # Table 2 'V4 (L) 839 (RSE 11)', bootstrap median 853 [407-1365]; control stream $THETA(4) 839

    # --- Oral absorption -------------------------------------------------
    # The oral dose is a zero-order input directly into the central
    # compartment: the control stream sets D3 (the duration for NONMEM
    # compartment 3 = CENTRAL) when LAI = 0, and scales that dose by
    # F3 = THETA(5) * EXP(ETA(2)). Intramuscular administration is taken as
    # completely bioavailable, so Foral is a RELATIVE bioavailability of the
    # oral route versus the intramuscular route (Thoueille 2024 Results 3.1).
    lfdepot_oral <- log(0.654) ; label("Relative bioavailability of oral versus intramuscular rilpivirine (Foral, fraction)")  # Table 2 'F oral (%) 65.4 (RSE 5)', bootstrap median 64.8 [52.6-73.3]; control stream $THETA(5) 0.654
    ld1_oral <- fixed(log(4))  ; label("Zero-order absorption duration for the oral dose (Doral, h)")                          # Table 2 'D oral (h) 4 FIX'; control stream $THETA(7) '4 FIX'. Thoueille 2024 Methods: fixed to 4 h from the Edurant label and preliminary estimation because too few samples were drawn just after oral intake

    # --- Long-acting intramuscular absorption ----------------------------
    # The injected dose is split between two parallel first-order pathways.
    # The fast fraction is carried on the logit scale exactly as in the
    # paper's eq. for TEMP: TEMP = ln(theta / (1 - theta)) and
    # Fi.m.fast_i = exp(TEMP + eta) / (1 + exp(TEMP + eta)), which keeps the
    # fraction inside (0, 1) for every eta draw. The slow fraction is the
    # remainder, Fi.m.slow = 1 - Fi.m.fast.
    logitfdepot <- log(0.276 / (1 - 0.276)) ; label("Logit of the fraction of the intramuscular dose released via the fast absorption pathway (Fi.m.fast, fraction)")  # Table 2 'F i.m.fast (%) 27.6 (RSE 9)', bootstrap median 27.5 [22.5-32.4]; control stream $THETA(6) 0.276; logit(0.276) = -0.9642
    lka  <- log(0.00214)  ; label("First-order absorption rate constant of the fast intramuscular pathway (kafast, 1/h)")  # Table 2 'k a fast (h-1) 0.00214 (RSE 11)', bootstrap median 0.00211 [0.00167-0.00266]; control stream $THETA(8) 0.00214
    lka2 <- log(0.000229) ; label("First-order absorption rate constant of the slow intramuscular pathway (kaslow, 1/h)")  # Table 2 'k a slow (h-1) 0.000229 (RSE 11)', bootstrap median 0.000225 [0.000108-0.000292]; control stream $THETA(9) 0.000229. ln(2)/kaslow = 3027 h = 18.0 weeks, the reported apparent half-life

    # --- Between-subject variability -------------------------------------
    # Values are the $OMEGA variances of the supplementary control stream.
    # Each reproduces the CV% printed in Table 2 under that row's stated
    # transformation, which is the check that the scale is right:
    #   log-normal rows use CV = sqrt(exp(omega^2) - 1)
    #   the logit row uses the paper's Table 2 footnote c approximation,
    #   CV_Fi.m.fast = theta * (1 - theta) * omega
    etalcl ~ 0.0649          # control stream $OMEGA 'IIV CL' 0.0649; sqrt(exp(0.0649) - 1) = 0.2589 = Table 2 'omega CL (CV%) 25.9 (RSE 9)'
    etalfdepot_oral ~ 0.129  # control stream $OMEGA 'IIV F3' 0.129; sqrt(exp(0.129) - 1) = 0.3711 = Table 2 'omega F oral (CV%) 37.1 (RSE 11)'
    etalogitfdepot ~ 0.708   # control stream $OMEGA 'IIV F1' 0.708; 0.276 * (1 - 0.276) * sqrt(0.708) = 0.1681 = Table 2 'omega F i.m.fast (CV%) 16.8 (RSE 15)'. Table 2 also notes a modest 40% shrinkage on this eta
    etalka2 ~ 0.521          # control stream $OMEGA 'IIV KA2 LAI' 0.521; sqrt(exp(0.521) - 1) = 0.8269 = Table 2 'omega k a slow (CV%) 82.7 (RSE 12)'

    # --- Inter-occasion variability on clearance -------------------------
    # $ABBR REPLACE ETA(OCC)=ETA(5,6,7,8,9,10) with $OMEGA BLOCK(1) 0.0167
    # followed by five BLOCK(1) SAME lines: six occasion slots sharing a
    # single variance. Occasions 2-6 are therefore fixed to the estimated
    # occasion-1 variance rather than estimated separately.
    etaiov_cl_1 ~ 0.0167         # control stream $OMEGA BLOCK(1) 'IOV CL' 0.0167; sqrt(exp(0.0167) - 1) = 0.1298 = Table 2 'omega IOV (CV%) 13.0 (RSE 10)'
    etaiov_cl_2 ~ fixed(0.0167)  # $OMEGA BLOCK(1) SAME
    etaiov_cl_3 ~ fixed(0.0167)  # $OMEGA BLOCK(1) SAME
    etaiov_cl_4 ~ fixed(0.0167)  # $OMEGA BLOCK(1) SAME
    etaiov_cl_5 ~ fixed(0.0167)  # $OMEGA BLOCK(1) SAME
    etaiov_cl_6 ~ fixed(0.0167)  # $OMEGA BLOCK(1) SAME

    # --- Residual error --------------------------------------------------
    # Thoueille 2024 Results 3.1: "a common mixed error model for both
    # routes of administration failed to estimate both components of the
    # error model. An additive error model best described rilpivirine RUV
    # after oral administration, while a proportional error model was
    # retained for rilpivirine RUV when administered i.m." The control
    # stream's $SIGMA holds variances, so both are square-rooted here.
    addSdOral <- sqrt(317) ; label("Additive residual standard deviation for oral observations (ng/mL)")        # control stream $SIGMA 317 'Add PO'; sqrt(317) = 17.80 = Table 2 'sigma add-oral (ng/mL) 18 (RSE 22)', bootstrap median 17.7 [10.5-28.1]
    propSdIm  <- sqrt(0.031) ; label("Proportional residual standard deviation for intramuscular observations (fraction)")  # control stream $SIGMA 0.031 'Prop LAI'; sqrt(0.031) = 0.1761 = Table 2 'sigma prop-LA (CV%) 18 (RSE 5)', bootstrap median 17.5 [15.2-20.1]
  })

  model({
    # =====================================================================
    # 1. Occasion decomposition for the inter-occasion variability on CL.
    #    Reproduces NONMEM's ETA(OCC) indirection: exactly one indicator is
    #    non-zero on each record, so iov_cl equals the eta of the current
    #    occasion. Records with OCC outside 1-6 would contribute no IOV, so
    #    the oral lead-in records are given OCC = 1 (see covariateData).
    # =====================================================================
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)
    oc5 <- (OCC == 5)
    oc6 <- (OCC == 6)
    iov_cl <-
      oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2 + oc3 * etaiov_cl_3 +
      oc4 * etaiov_cl_4 + oc5 * etaiov_cl_5 + oc6 * etaiov_cl_6

    # =====================================================================
    # 2. Individual parameters. Control stream $PK:
    #      CL = TVCL * EXP(ETA(1) + ETA(OCC))
    #      F3 = TVF3 * EXP(ETA(2))
    #      F1 = EXP(TEMP + ETA(3)) / (1 + EXP(TEMP + ETA(3)))
    #      KA1 = THETA(8)                  (no eta)
    #      KA2 = TVKA2 * EXP(ETA(4))
    #    V3, Q, V4 and D3 carry no random effect.
    # =====================================================================
    cl <- exp(lcl + etalcl + iov_cl)
    vc <- exp(lvc)
    q <- exp(lq)
    vp <- exp(lvp)

    ka <- exp(lka)
    ka2 <- exp(lka2 + etalka2)

    fdepot_oral <- exp(lfdepot_oral + etalfdepot_oral)
    d1_oral <- exp(ld1_oral)
    fdepot <- expit(logitfdepot + etalogitfdepot)

    # =====================================================================
    # 3. Micro-constants (control stream $PK: K34 = Q/V3, K43 = Q/V4,
    #    K30 = CL/V3).
    # =====================================================================
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # =====================================================================
    # 4. ODE system, transcribed from the control stream $DES:
    #      DADT(1) = -KA1*A(1)
    #      DADT(2) = -KA2*A(2)
    #      DADT(3) = KA1*A(1) + KA2*A(2) - K34*A(3) + K43*A(4) - K30*A(3)
    #      DADT(4) = K34*A(3) - K43*A(4)
    # =====================================================================
    d/dt(depot) <- -ka * depot
    d/dt(depot2) <- -ka2 * depot2
    d/dt(central) <- ka * depot + ka2 * depot2 + k21 * peripheral1 -
      k12 * central - kel * central
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # =====================================================================
    # 5. Dose routing.
    #
    #    An INTRAMUSCULAR injection is entered as TWO simultaneous dose
    #    records, each carrying the full injected amount (900 mg for the
    #    two-monthly regimen), one into `depot` and one into `depot2`; the
    #    f() multipliers below split it into the fast and slow pathways.
    #    Intramuscular administration is assumed 100% bioavailable
    #    (Thoueille 2024 Results 3.1), so the two fractions sum to 1.
    #
    #    An ORAL dose is entered as a single record into `central` with
    #    rate = -2, which makes rxode2 read the modelled duration below and
    #    deliver Foral * amt as a 4 h zero-order input. Note that f(central)
    #    applies to ANY dose placed in `central`; a user simulating an
    #    intravenous dose would have to account for the 0.654 multiplier.
    # =====================================================================
    f(depot) <- fdepot
    f(depot2) <- 1 - fdepot
    f(central) <- fdepot_oral
    dur(central) <- d1_oral

    # =====================================================================
    # 6. Observation. Doses are in mg and vc is in L, so central / vc is
    #    mg/L; multiplying by 1000 gives ng/mL, matching the control
    #    stream's S3 = V3/1000.
    #
    #    The residual error switches by route ($ERROR selects Y0 or Y1 on
    #    LAI). Writing it as a combined additive-plus-proportional model
    #    whose two magnitudes are gated by ROUTE_ORAL reproduces both
    #    branches exactly: on an oral record the proportional term is zero
    #    so the SD is addSdOral, and on an intramuscular record the additive
    #    term is zero so the SD is propSdIm * Cc.
    # =====================================================================
    Cc <- 1000 * central / vc
    ruvAdd <- addSdOral * ROUTE_ORAL
    ruvProp <- propSdIm * (1 - ROUTE_ORAL)
    Cc ~ add(ruvAdd) + prop(ruvProp)
  })
}
