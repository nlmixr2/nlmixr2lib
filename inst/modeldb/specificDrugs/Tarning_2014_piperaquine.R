# Population pharmacokinetic model of oral piperaquine in adults with
# uncomplicated Plasmodium falciparum malaria from a randomized study of
# concomitant food (low-fat chocolate milk, 6.4 g fat) vs fasting in
# Thailand (Tarning 2014, Antimicrob Agents Chemother 58(4):2052-2058;
# doi:10.1128/AAC.02318-13).

Tarning_2014_piperaquine <- function() {
  description <- paste(
    "Population PK model for oral piperaquine in adults with uncomplicated",
    "Plasmodium falciparum malaria in Thailand (Tarning 2014; n = 30, fed",
    "vs fasting parallel design). Three-transit-compartment absorption",
    "(ka = ktr) feeding a three-compartment disposition model. Allometric",
    "body-weight scaling on all clearances (fixed exponent 0.75) and",
    "volumes (fixed exponent 1.0); 70 kg reference. Linear dose-occasion",
    "effect on relative bioavailability (+25.3% per consecutive dose, OCC",
    "= 1, 2, 3). Linear age effect on the first peripheral volume of",
    "distribution (+4.10% per year of age). Relative bioavailability",
    "anchored at 1 with between-dose-occasion variability (no BSV in the",
    "final model). Concomitant low-fat food was tested as a covariate but",
    "was not retained in the final model.",
    sep = " "
  )
  reference <- paste(
    "Tarning J, Lindegardh N, Lwin KM, Annerberg A, Kiricharoen L, Ashley",
    "E, White NJ, Nosten F, Day NPJ (2014).",
    "Population Pharmacokinetic Assessment of the Effect of Food on",
    "Piperaquine Bioavailability in Patients with Uncomplicated Malaria.",
    "Antimicrobial Agents and Chemotherapy 58(4):2052-2058.",
    "doi:10.1128/AAC.02318-13.",
    sep = " "
  )
  vignette <- "Tarning_2014_piperaquine"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight at admission",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometric scaling with fixed exponents 0.75 on apparent",
        "clearances (CL/F, Q1/F, Q2/F) and 1.0 on apparent volumes",
        "(Vc/F, Vp1/F, Vp2/F). Tarning 2014 Methods page 2053: 'Body",
        "weight was investigated as a covariate by using an allometric",
        "function with power values of 3/4 for clearance parameters and 1",
        "for volume parameters'. The paper does not state the allometric",
        "reference weight explicitly; 70 kg (standard adult anchor) is",
        "used here because it reproduces the published post-hoc per-kg",
        "estimates (CL/F = 1.02 L/h/kg, V/F = 516 L/kg in Table 3) to",
        "within ~3% for a cohort-median patient, whereas using the",
        "cohort-median weight (~51 kg) as the reference would imply per-",
        "kg CL of ~1.33 L/h/kg and contradict Table 3. See vignette",
        "Assumptions and deviations. Treated as time-fixed at admission.",
        sep = " "
      ),
      source_name        = "WT"
    ),
    AGE = list(
      description        = "Subject age at admission",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Linear effect on the first peripheral volume of distribution.",
        "Tarning 2014 Results page 2054: 'This covariate relationship",
        "resulted in a linear increase in the peripheral volume of",
        "distribution of 4.10% for each year of age increase'. The paper",
        "does not state the reference age; the cohort-pooled median",
        "(~33 years, computed as the midpoint of the fasting median 38",
        "and fed median 28 reported in Table 1) is used here as the",
        "centring constant so that the published typical Vp1/F = 6240 L",
        "corresponds to a median-age patient. See vignette Assumptions",
        "and deviations. Time-fixed at admission.",
        sep = " "
      ),
      source_name        = "AGE"
    ),
    OCC = list(
      description        = "Dose-occasion counter (1 = first daily dose, 2 = second, 3 = third)",
      units              = "(count)",
      type               = "count",
      reference_category = NULL,
      notes              = paste(
        "Linear dose-occasion effect on relative bioavailability.",
        "Tarning 2014 Results page 2054 and Table 2 'Dose effect on F",
        "(%)': 'this could be simplified to a linear covariate",
        "relationship (i.e., 25.3% increase in bioavailability per dose)",
        "without a significant reduction in model fit'. Encoded",
        "additively in OCC: F_OCC = 1 + 0.253 * (OCC - 1), so OCC = 1",
        "-> F = 1.000, OCC = 2 -> F = 1.253, OCC = 3 -> F = 1.506. This",
        "linear form reproduces the bootstrap categorical estimates",
        "reported in the same paragraph (21.9% at dose 2 and 50.9% at",
        "dose 3 relative to dose 1) to within ~3%. OCC is carried on",
        "each dose event row; observation rows inherit the most recent",
        "dose's OCC value but the value is only used by f(depot) at",
        "dose-administration time.",
        sep = " "
      ),
      source_name        = "OCC"
    )
  )

  population <- list(
    species         = "human",
    n_subjects      = 30L,
    n_studies       = 1L,
    n_fasting       = 15L,
    n_fed           = 15L,
    age_range       = "18-55 years (Table 1; fasting 18-55, fed 19-45)",
    age_median      = "38 years fasting; 28 years fed (Table 1)",
    weight_range    = "39-73 kg (Table 1; fasting 39-62, fed 45-73)",
    weight_median   = "50 kg fasting; 53 kg fed (Table 1)",
    sex_female_pct  = 13.3,
    disease_state   = paste(
      "Uncomplicated Plasmodium falciparum malaria in Thai adults",
      "presenting at the Shoklo Malaria Research Unit, Mae Sot, Thailand.",
      "Inclusion: microscopy-confirmed asexual P. falciparum or mixed",
      "infection, no signs of severe malaria, hematocrit >= 30%, ages",
      "16-65 years (Methods page 2053)."
    ),
    dose_range      = paste(
      "Dihydroartemisinin-piperaquine fixed-dose combination (DuoCotecxin;",
      "40 mg dihydroartemisinin + 320 mg piperaquine phosphate per",
      "tablet), once daily for 3 days, weight-based dosing to achieve a",
      "daily target of 18 mg piperaquine phosphate/kg (Methods page",
      "2053). Median delivered dose was 17.2 mg/kg/day (fasting) and",
      "17.5 mg/kg/day (fed) per Table 1."
    ),
    regions         = "Thailand (Shoklo Malaria Research Unit, Mae Sot)",
    food_arms       = paste(
      "15 fasting and 15 fed (chocolate milk 200 mL = 6.4 g fat per",
      "dose). Concomitant food was tested as a covariate but was not",
      "retained in the final model (Results page 2054)."
    ),
    notes           = paste(
      "Demographics from Tarning 2014 Table 1. PK analysis used 985",
      "post-dose piperaquine samples (1,016 collected, 91 BLQ below",
      "1.2 ng/mL omitted; 82.6% of BLQ samples were in the long terminal",
      "phase > 70 days). Intensive sampling pre-dose and at 0.5-12 h",
      "after each daily dose; sparse follow-up venous samples at days",
      "4, 5, 7, 14, 21, 28, 42, 56, 70, 84, 98, 112, and 126. NONMEM",
      "v7.2 with ADVAN5 and FOCE-I."
    )
  )

  ini({
    # Structural parameters from Tarning 2014 Table 2 ("Value" column).
    # Apparent values (relative to F = 1) reported on the linear scale;
    # log() applied here for the nlmixr2 internal log scale. Typical
    # values correspond to the standard adult allometric reference weight
    # 70 kg (the paper does not state the reference; see WT covariate
    # notes for the inference from the post-hoc per-kg estimates in Table
    # 3).
    lcl  <- log(67.6)
    label("Apparent piperaquine elimination clearance CL/F at WT = 70 kg (L/h)")
    # Tarning 2014 Table 2: CL/F = 67.6 L/h (RSE 11.6%; 95% CI 54.0-85.5)

    lvc  <- log(3030)
    label("Apparent central volume of distribution Vc/F at WT = 70 kg (L)")
    # Tarning 2014 Table 2: Vc/F = 3030 L (RSE 16.4%; 95% CI 2160-4180)

    lq   <- log(408)
    label("Apparent inter-compartmental clearance to first peripheral compartment Q1/F at WT = 70 kg (L/h)")
    # Tarning 2014 Table 2: Q1/F = 408 L/h (RSE 15.0%; 95% CI 309-557)

    lvp  <- log(6240)
    label("Apparent first peripheral volume of distribution Vp1/F at WT = 70 kg and AGE = 33 y (L)")
    # Tarning 2014 Table 2: Vp1/F = 6240 L (RSE 14.6%; 95% CI 4890-8530)

    lq2  <- log(109)
    label("Apparent inter-compartmental clearance to second peripheral compartment Q2/F at WT = 70 kg (L/h)")
    # Tarning 2014 Table 2: Q2/F = 109 L/h (RSE 13.6%; 95% CI 83.3-143)

    lvp2 <- log(24400)
    label("Apparent second peripheral volume of distribution Vp2/F at WT = 70 kg (L)")
    # Tarning 2014 Table 2: Vp2/F = 24400 L (RSE 10.1%; 95% CI 20000-29500)

    lmtt <- log(2.04)
    label("Mean transit time of the 3-compartment transit-absorption chain MTT (h)")
    # Tarning 2014 Table 2: MTT = 2.04 h (RSE 7.50%; 95% CI 1.80-2.41)

    # Relative bioavailability fixed at 1 (Tarning 2014 Table 2: 'F (%)
    # = 100 (fixed)'). The paper retained between-dose-occasion
    # variability (IOV) on F (48.8% CV) but dropped between-subject
    # variability after IOV implementation (Results page 2054:
    # 'However, the between-subject variability could be removed after
    # implementation of between-dose occasion variability without a
    # major impact on the model fit'). Because nlmixr2lib does not
    # natively encode IOV, the published F IOV is captured here as a
    # forward-simulation IIV term (same idiom as the IOV terms in
    # Hoglund_2012_piperaquine.R; see vignette Assumptions and
    # deviations).
    lfdepot <- fixed(log(1))
    label("Relative bioavailability F at OCC = 1 (unitless, fixed at 1)")
    # Tarning 2014 Table 2: F (%) = 100 (fixed)

    # Allometric exponents fixed by the source paper (strong biological
    # prior; the exponents are written into the model rather than
    # estimated). Tarning 2014 Methods page 2053: 'Body weight was
    # investigated as a covariate by using an allometric function with
    # power values of 3/4 for clearance parameters and 1 for volume
    # parameters'.
    e_wt_cl <- fixed(0.75)
    label("Allometric WT exponent on all apparent clearance parameters (CL/F, Q1/F, Q2/F)")
    # Tarning 2014 Methods page 2053: clearance power = 3/4 (fixed)

    e_wt_vc <- fixed(1.00)
    label("Allometric WT exponent on all apparent volume parameters (Vc/F, Vp1/F, Vp2/F)")
    # Tarning 2014 Methods page 2053: volume power = 1 (fixed)

    # Age effect on first peripheral volume of distribution. Tarning
    # 2014 Table 2 'Age effect on Vp1 (%)' = 4.10 (RSE 18.1%; 95% CI
    # 2.38-5.32) and Results page 2054: 'linear increase in the
    # peripheral volume of distribution of 4.10% for each year of age
    # increase'. Centred on the cohort-pooled median age 33 y (paper
    # does not state the reference; see AGE covariate notes).
    e_age_vp <- 0.0410
    label("Linear age effect on Vp1/F (fractional change per year of age)")
    # Tarning 2014 Table 2: Age effect on Vp1 (%) = 4.10 (RSE 18.1%; 95% CI 2.38-5.32)

    # Dose-occasion effect on relative bioavailability. Tarning 2014
    # Table 2 'Dose effect on F (%)' = 25.3 (RSE 34.4%; 95% CI 9.82-
    # 53.2) and Results page 2054: '25.3% increase in bioavailability
    # per dose'. Encoded additively in OCC: F_OCC = 1 + 0.253 *
    # (OCC - 1).
    e_doseocc_f <- 0.253
    label("Increment in relative bioavailability per additional dose occasion (fraction)")
    # Tarning 2014 Table 2: Dose effect on F (%) = 25.3 (RSE 34.4%; 95% CI 9.82-53.2)

    # IIV. Coefficient of variation reported in Tarning 2014 Table 2
    # ("IIV (% CV)" column). The same NONMEM convention used by the
    # Tarning group (cf. Hoglund 2017 footnote: '100 * (e^variance -
    # 1)^(1/2)') gives the internal log-scale variance as
    # omega^2 = log(CV^2 + 1).
    #
    #   CL   IIV 24.4%   -> omega^2 = log(0.244^2 + 1) = 0.057846
    #   Vc   IIV 51.6%   -> omega^2 = log(0.516^2 + 1) = 0.234710
    #   Vp1  IIV 45.6%   -> omega^2 = log(0.456^2 + 1) = 0.188103
    #   Q2   IIV 25.8%   -> omega^2 = log(0.258^2 + 1) = 0.064375
    #   MTT  IIV 24.1%   -> omega^2 = log(0.241^2 + 1) = 0.056498  (BSV)
    #   F    IOV 48.8%   -> omega^2 = log(0.488^2 + 1) = 0.212074  (IOV treated as forward-simulation IIV; see lfdepot comment)
    #
    # No IIV on Q1 or Vp2 (Table 2 reports no IIV value for either).
    # MTT additionally carries a 39.4% IOV (Table 2) that is not encoded
    # here as a separate slot; only the published BSV is retained on
    # MTT in the forward-simulation model.
    etalcl     ~ 0.057846
    # Tarning 2014 Table 2: IIV on CL/F = 24.4% CV (RSE 26.0%; 95% CI 17.4-29.6)

    etalvc     ~ 0.234710
    # Tarning 2014 Table 2: IIV on Vc/F = 51.6% CV (RSE 32.3%; 95% CI 31.2-68.1)

    etalvp     ~ 0.188103
    # Tarning 2014 Table 2: IIV on Vp1/F = 45.6% CV (RSE 48.8%; 95% CI 18.8-68.4)

    etalq2     ~ 0.064375
    # Tarning 2014 Table 2: IIV on Q2/F = 25.8% CV (RSE 48.2%; 95% CI 6.67-37.9)

    etalmtt    ~ 0.056498
    # Tarning 2014 Table 2: BSV on MTT = 24.1% CV (RSE 52.7%; 95% CI 8.77-35.6); separate IOV 39.4% CV not encoded

    etalfdepot ~ 0.212074
    # Tarning 2014 Table 2: IOV on F = 48.8% CV (RSE 16.6%; 95% CI 38.3-56.0); BSV removed after IOV implementation. Encoded as forward-simulation IIV; see vignette Errata.

    # Residual error. Tarning 2014 Methods page 2053: 'Piperaquine
    # plasma concentrations were transformed into their natural
    # logarithms ... The random residual variability was assumed to be
    # additive, since data were modeled as natural logarithms (i.e.,
    # essentially equivalent to an exponential error model for
    # untransformed data)'. This NONMEM additive-on-log-scale residual
    # maps to a proportional residual in linear concentration space (see
    # references/parameter-names.md Residual error and the matching
    # comment in Hoglund_2012_piperaquine.R). Table 2 reports
    # sigma = 30.7% (% CV); for the log-additive parameterisation this
    # value is the SD on the log scale, which equals the proportional
    # CV to first order.
    propSd <- 0.307
    label("Proportional residual SD for piperaquine plasma concentration (SD on log scale)")
    # Tarning 2014 Table 2: sigma (% CV) = 30.7 (RSE 4.42%; 95% CI 27.8-33.5)
  })

  model({
    # Age-dependent first peripheral volume scaling. Centred on the
    # cohort-pooled median age 33 y (Tarning 2014 Table 1 medians; paper
    # does not state the reference age, see AGE covariate notes).
    age_vp <- 1 + e_age_vp * (AGE - 33)

    # Individual PK parameters. Allometric scaling on clearances
    # (exponent 0.75) and volumes (exponent 1) centred on the standard
    # adult anchor WT = 70 kg (paper does not state the reference; see
    # WT covariate notes). The age effect multiplies Vp1/F. No IIV on
    # Q1/F or Vp2/F (Tarning 2014 Table 2).
    cl  <- exp(lcl  + etalcl)  * (WT / 70)^e_wt_cl
    vc  <- exp(lvc  + etalvc)  * (WT / 70)^e_wt_vc
    q   <- exp(lq)             * (WT / 70)^e_wt_cl
    vp  <- exp(lvp  + etalvp)  * (WT / 70)^e_wt_vc * age_vp
    q2  <- exp(lq2  + etalq2)  * (WT / 70)^e_wt_cl
    vp2 <- exp(lvp2)           * (WT / 70)^e_wt_vc

    # Mean transit time and chain rate constant. Tarning 2014 Methods
    # page 2053 and Results page 2054: 'a more flexible transit
    # compartment absorption with a fixed number of transit compartments
    # for the population' and 'A transit compartment (n = 3) absorption
    # model described the absorption phase well ... The absorption rate
    # from the last transit compartment could be set as identical to the
    # rate constant between transit compartments without a significant
    # impact on the model'. With NN = 3 transit compartments and ka =
    # ktr, the absorption chain depot -> transit1 -> transit2 ->
    # transit3 -> central has (NN + 1) = 4 equal-rate transitions; the
    # mean transit time is MTT = (NN + 1) / ktr, so ktr = 4 / MTT
    # (Savic & Karlsson 2007 convention; same idiom as
    # Hoglund_2012_piperaquine.R).
    mtt <- exp(lmtt + etalmtt)
    ktr <- 4 / mtt

    # Three-compartment disposition micro-constants.
    kel <- cl  / vc
    k12 <- q   / vc
    k21 <- q   / vp
    k13 <- q2  / vc
    k31 <- q2  / vp2

    # ODE system: 3-transit-compartment absorption chain feeding into a
    # 3-compartment disposition model (Tarning 2014 Figure 1). The same
    # ktr propagates the dose through the depot and three transit
    # compartments into central.
    d/dt(depot)       <- -ktr * depot
    d/dt(transit1)    <-  ktr * depot    - ktr * transit1
    d/dt(transit2)    <-  ktr * transit1 - ktr * transit2
    d/dt(transit3)    <-  ktr * transit2 - ktr * transit3
    d/dt(central)     <-  ktr * transit3 - kel * central -
                          k12 * central + k21 * peripheral1 -
                          k13 * central + k31 * peripheral2
    d/dt(peripheral1) <-  k12 * central  - k21 * peripheral1
    d/dt(peripheral2) <-  k13 * central  - k31 * peripheral2

    # Relative bioavailability applied to the depot. Anchor lfdepot =
    # log(1) is fixed (Tarning 2014 final model). Dose-occasion-
    # dependent multiplier F_OCC = 1 + e_doseocc_f * (OCC - 1)
    # increments F by 25.3% per successive dose (Tarning 2014 Table 2
    # Dose effect on F = 25.3%). Variability captured by etalfdepot
    # (the published IOV, encoded here as forward-simulation IIV).
    f_occ <- 1 + e_doseocc_f * (OCC - 1)
    f(depot) <- f_occ * exp(lfdepot + etalfdepot)

    # Piperaquine plasma concentration. The model's typical CL/F,
    # volumes, and bioavailability are reported on a piperaquine-base
    # basis (the LC-MS assay measures piperaquine, not the dosed
    # phosphate salt), so input doses must be supplied in mg of
    # piperaquine base. Users with phosphate-salt dose records should
    # apply a salt-to-base conversion before passing the dose to the
    # model (the Tarning group convention; cf.
    # Hoglund_2017_piperaquine.R Methods page 4 -- piperaquine
    # phosphate is converted to piperaquine base by a 57.7% scale
    # factor for the tetra-phosphate anhydrous salt). Vc in L gives
    # central / vc in mg/L; the 1000 factor converts to ng/mL for
    # direct comparison against Tarning 2014 Table 3 (Cmax, day-7
    # concentrations all in ng/mL).
    Cc <- 1000 * central / vc

    # Proportional residual error on the linear-concentration scale
    # (NONMEM additive-on-log-scale -> nlmixr2 proportional; see
    # ini() comment on propSd).
    Cc ~ prop(propSd)
  })
}
