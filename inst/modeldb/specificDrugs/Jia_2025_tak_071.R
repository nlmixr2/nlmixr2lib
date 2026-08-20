Jia_2025_tak_071 <- function() {
  description <- "One-compartment population PK model for TAK-071 (muscarinic M1 positive allosteric modulator) in healthy adults and participants with Parkinson disease, with a fixed absorption lag, a 3-step transit chain, and parallel slow/quick first-order absorption; dose effects on relative bioavailability and absorption rate, formulation (tablet vs drug-in-capsule) effects on the slow-absorption fraction and rate, age on clearance, and body weight on central volume"
  reference <- paste(
    "Jia H, Facius A, Jennings R, Hang Y, Padmanabhan J, Shanbhag NM,",
    "Harel BT, Simen A, Yin W. Population Pharmacokinetic and",
    "Exposure-Response Analysis of the Cognitive Effects of TAK-071 in",
    "Participants With Parkinson Disease and Cognitive Impairment.",
    "Clin Pharmacol Drug Dev. 2025;14(11):856-868. doi:10.1002/cpdd.1579.",
    "Structural equations from the Supplemental Code (NONMEM control stream",
    "of the final model, Supplemental Information page S15-S17); parameter",
    "values from Table 1.",
    sep = " "
  )
  vignette <- "Jia_2025_tak_071"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    AGE = list(
      description        = "Baseline age",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power covariate on clearance, normalized to 70 years (Methods, Population PK Model Development: 'age to 70 years'; Table 1 footnote c). Control-stream line TVCL = TVCL * (AGE/70)**THETA(17). Observed range 18-83 years (Table S1).",
      source_name        = "AGE"
    ),
    WT = list(
      description        = "Baseline body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power covariate on central volume, normalized to 80 kg (Methods, Population PK Model Development: 'body weight was normalized to 80 kg'; Table 1 footnote d). Control-stream line TVVC = TVVC * (WEIGHT/80)**THETA(16). Observed range 47.3-122 kg (Table S1).",
      source_name        = "WEIGHT"
    ),
    DOSE_TAK071_MG = list(
      description        = "Administered TAK-071 dose level for the current dose record",
      units              = "mg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Per-dose-record dose level in mg, normalized to 5 mg (Methods: 'dose to 5 mg'; Table 1 footnote b). Enters twice in the control stream: a power effect on relative bioavailability (TVFREL = THETA(1) * (DOSE/5)**THETA(14)) and a power effect on the slow absorption rate (TVKA = TVKA * (DOSE/5)**THETA(12)). Relative bioavailability and absorption rate both decrease as dose increases. Studied dose levels: 1, 3, 9, 20, 40, 80, 120, 160 mg single dose and 3, 5, 7.5, 9, 15 mg once daily.",
      source_name        = "DOSE"
    ),
    FORM_TABLET = list(
      description        = "Oral formulation indicator: 1 = tablet, 0 = drug-in-capsule (DIC)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (drug-in-capsule, DIC)",
      notes              = "Per-dose-record formulation flag. NOTE: the comparator here is a CAPSULE (drug-in-capsule), not the non-tablet oral liquid named as this canonical's default reference category -- same orientation as Wada_2023_sparsentan.R. Source column FORM, where FORM == 2 is the tablet; the control stream applies both tablet effects under IF(FORM == 2). The tablet carries (a) a multiplicative factor on the LOGIT of the slow-absorption fraction (IF(FORM == 2) TVFRAC_ = TVFRAC_ * THETA(15)) and (b) a multiplicative factor on the slow absorption rate (IF(FORM == 2) TVKA = TVKA * THETA(13)). Consistent with the observed median tmax of 2.00 h for the tablet vs 4.98 h for the capsule after a single 10 mg dose (Discussion).",
      source_name        = "FORM"
    ),
    SAMPLE_INTENSIVE = list(
      description        = "Sampling-design indicator: 1 = dense (intensive) PK sampling, 0 = sparse PK sampling",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (sparse sampling)",
      notes              = "Source column SAMPLING, where SAMPLING == 2 is the sparse subset (participants with Parkinson disease in the Phase 2 main cohort, TAK-071-2002). Unlike the Macpherson 2015 founding example this flag does NOT vary within a subject in this analysis -- it is a per-subject design attribute, so set it once per participant. It switches three things in the control stream: the between-subject random effect on the logit slow-absorption fraction is suppressed (IF(SAMPLING == 2) FRAC_ = TVFRAC_), the between-subject random effect on the slow absorption rate is suppressed (IF(SAMPLING == 2) KA = TVKA), and the residual error is inflated (IF(SAMPLING == 2) W = W * THETA(11), ratio 1.21). Rationale (Results): 'The sparse PK data collected in patients with PD did not allow for reliable estimation of individual absorption rates and fraction parameters, which is why no random effects for ka and frac were estimated for this subset. In addition, the model could quantify an inflated residual variability in this subset of approximately 21%.'",
      source_name        = "SAMPLING"
    )
  )

  compartmentData <- list(
    depot    = list(analyte = "TAK-071", units = "mg", specimen = "administration site", verified = TRUE),
    transit1 = list(analyte = "TAK-071", units = "mg", specimen = "administration site", verified = TRUE),
    transit2 = list(analyte = "TAK-071", units = "mg", specimen = "administration site", verified = TRUE),
    depot2   = list(analyte = "TAK-071", units = "mg", specimen = "administration site", verified = TRUE),
    depot3   = list(analyte = "TAK-071", units = "mg", specimen = "administration site", verified = TRUE),
    central  = list(analyte = "TAK-071", units = "mg", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 165L,
    n_studies      = 2L,
    age_range      = "18-83 years",
    age_median     = "44 years",
    weight_range   = "47.3-122 kg",
    weight_median  = "78.7 kg",
    sex_female_pct = 11.8,
    race_ethnicity = c(
      White = 65.4, `Black or African American` = 15.0, Asian = 12.4,
      Multiple = 3.3, `Native Hawaiian or other Pacific Islander` = 2.0,
      `American Indian or Alaska Native` = 0.7, `Not reported` = 1.3
    ),
    disease_state  = "healthy adults, and participants with Parkinson disease with cognitive impairment and an elevated risk of falls",
    dose_range     = "1-160 mg single oral dose and 3-15 mg once daily (Phase 1, TAK-071-1001); 7.5 mg single dose and 5 or 7.5 mg once daily (Phase 2, TAK-071-2002)",
    renal_function = "creatinine clearance 56.9-180 mL/min/1.73m2; eGFR 57.1-143 mL/min (Table S1); neither was a significant covariate on clearance",
    notes          = "Two studies: TAK-071-1001 (Phase 1 single- and multiple-ascending-dose with a 3-way crossover food-effect and relative-bioavailability arm, healthy adults 18-55 years, 96% men, 16% Japanese, NCT02769065) and TAK-071-2002 (Phase 2 randomized double-blind placebo-controlled 2-period crossover, sentinel cohort of healthy participants older than 55 years plus a main cohort of participants with Parkinson disease aged 40-85 years, no Japanese participants, NCT04334317). Counts differ between the paper's own sources: the Abstract and the Results Data Set section give 165 participants (112 healthy = 104 from TAK-071-1001 plus 8 sentinel, and 53 with Parkinson disease = 37 at 5 mg plus 16 at 7.5 mg), whereas Table S1 reports 153 individuals (92 healthy, 61 with Parkinson disease). The demographic summaries recorded here (median/range, sex, race) are Table S1 values and are therefore computed on the n = 153 denominator."
  )

  ini({
    # ---- Structural parameters -------------------------------------------
    # Absorption lag on the dosing (gut) compartment.
    ltlag <- fixed(log(0.2)); label("Absorption lag time (log h)")  # Supplemental Code $THETA(2): 0.2 FIX; Results: "including a fixed tlag of 0.2 hours"

    # Mean transit time. Table 1 reports MTT = 1.54 h INCLUDING the fixed
    # 0.2 h lag ("an MTT (including a fixed tlag of 0.2 hours), estimated at
    # 1.54 hours"). The control stream subtracts the lag before forming the
    # transit rate: TVMTT = THETA(3) - ALAG1, so the transit-chain MTT is
    # 1.54 - 0.2 = 1.34 h.
    lmtt <- log(1.54 - 0.2); label("Mean transit time, excluding the lag (log h)")  # Table 1 MTT TV 1.54 h (RSE 2.3%, 95% CI 1.47-1.61) minus the fixed 0.2 h lag; Supplemental Code $PK line TVMTT = THETA(3) - ALAG1

    # Fraction of the dose absorbed by the slow route, on the logit scale.
    logitfslow <- log(0.804 / (1 - 0.804)); label("Fraction absorbed by the slow route (logit)")  # Table 1 Fraction absorbed slowly TV 80.4% (RSE 0.5%, 95% CI 79.7-81.2)

    # Parallel first-order absorption rates. The control stream estimates the
    # QUICK rate as a ratio to the slow rate (KAQ = THETA(6); K5T6 = KA * KAQ),
    # so the two rates share the dose effect, the tablet effect, and the
    # between-subject random effect. Encoding lka_quick as
    # log(slow TV * ratio) reproduces that exactly:
    # 0.883 * 71.0 = 62.7 /h, the value printed in the Results text.
    lka_slow  <- log(0.883);        label("Slow first-order absorption rate at 5 mg, drug-in-capsule (log 1/h)")  # Table 1 Absorption rate TV 0.883 1/h (RSE 1.9%, 95% CI 0.851-0.916), footnote b: typical value for a 5 mg dose with formulation DIC
    lka_quick <- log(0.883 * 71.0); label("Quick first-order absorption rate at 5 mg, drug-in-capsule (log 1/h)")  # Table 1 Absorption ratio quick/slow TV 71.0 (RSE 1.9%, 95% CI 68.3-73.7) times the 0.883 1/h slow TV; Results text: "the faster absorption route had a ka,quick of 62.7"

    lcl <- log(0.566); label("Clearance at 70 years of age (log L/h)")   # Table 1 Clearance TV 0.566 L/h (RSE 2.2%, 95% CI 0.542-0.591), footnote c: typical value for a person aged 70 years
    lvc <- log(47.5);  label("Central volume of distribution at 80 kg (log L)")  # Table 1 Volume TV 47.5 L (RSE 2.7%, 95% CI 45.0-50.0), footnote d: typical value for a person with 80 kg body weight

    # Relative bioavailability anchor: THETA(1) was fixed at 100%, so the
    # 5 mg reference dose has Frel = 1 by construction and only the dose
    # effect below is estimated.
    lfdepot <- fixed(log(1)); label("Relative bioavailability at the 5 mg reference dose (log fraction)")  # Supplemental Code $THETA(1): 100 FIX (FREL, %); F1 = FREL/100

    # ---- Covariate effects -----------------------------------------------
    e_dose_ka_slow        <- -0.885;    label("Power exponent on (DOSE/5) for the slow absorption rate (unitless)")  # Table 1 Absorption rate dose-effect (power) -0.885 (RSE 5.7%, 95% CI -0.984 to -0.787)
    e_form_tablet_ka_slow <- log(27.0); label("Log ratio of tablet to drug-in-capsule slow absorption rate (unitless)")  # Table 1 Absorption rate tablet-effect (ratio) 27.0 (RSE 7.6%, 95% CI 23.0-31.0); control stream IF(FORM == 2) TVKA = TVKA * THETA(13)
    # NOTE: the tablet effect on the slow-absorption fraction multiplies the
    # LOGIT of the fraction, not the fraction itself -- IF(FORM == 2)
    # TVFRAC_ = TVFRAC_ * THETA(15) where TVFRAC_ = logit(TVFRAC). It is a
    # ratio on the logit scale, so a tablet fraction of
    # expit(logit(0.804) * 0.419) = 64.4% is slow-absorbed, not 0.804 * 0.419.
    e_form_tablet_fslow   <- 0.419;     label("Ratio applied to the logit slow-absorption fraction for tablets (unitless)")  # Table 1 Fraction absorbed slowly tablet-effect (ratio) 0.419 (RSE 5.8%, 95% CI 0.372-0.467); control stream IF(FORM == 2) TVFRAC_ = TVFRAC_ * THETA(15)
    e_age_cl              <- -0.171;    label("Power exponent on (AGE/70) for clearance (unitless)")  # Table 1 Clearance age-effect (power) -0.171 (RSE 37.9%, 95% CI -0.298 to -0.0439)
    e_wt_vc               <- 0.913;     label("Power exponent on (WT/80) for central volume (unitless)")  # Table 1 Volume body weight-effect (power) 0.913 (RSE 13.8%, 95% CI 0.666-1.16)
    e_dose_fdepot         <- -0.0965;   label("Power exponent on (DOSE/5) for relative bioavailability (unitless)")  # Table 1 Relative bioavailability dose-effect (power) -0.0965 (RSE 13.8%, 95% CI -0.123 to -0.0705)

    # ---- Between-subject variability -------------------------------------
    # Table 1 footnote a: "Results for random-effect parameters are shown as
    # standard deviations/correlations as reported by the [NONMEM] program",
    # so the tabulated BSV values are omega SDs and are squared here to give
    # variances. Confirmed by Supplemental Equation 2 (CV = sqrt(exp(omega^2) - 1)):
    # sqrt(exp(0.630^2) - 1) = 69.8% (Results, MTT), sqrt(exp(0.373^2) - 1) = 38.6%
    # (Results, CL) and sqrt(exp(0.233^2) - 1) = 23.6% (Results, Vc), each
    # matching the CV% quoted in the Results text.
    etalmtt       ~ 0.630^2;  # Table 1 MTT BSV 0.630 (RSE 2.8%, 95% CI 0.596-0.665) -> 69.8% CV
    etalogitfslow ~ 1.13^2;   # Table 1 Fraction absorbed slowly BSV 1.13 (RSE 3.8%, 95% CI 1.05-1.21), additive on the logit scale
    etalka_slow   ~ 0.788^2;  # Table 1 Absorption rate BSV 0.788 (RSE 4.0%, 95% CI 0.727-0.849)
    etalcl        ~ 0.373^2;  # Table 1 Clearance BSV 0.373 (RSE 3.2%, 95% CI 0.349-0.396) -> 38.6% CV
    etalvc        ~ 0.233^2;  # Table 1 Volume BSV 0.233 (RSE 2.6%, 95% CI 0.221-0.245) -> 23.6% CV
    # No between-subject variability on the quick/slow absorption ratio:
    # Supplemental Code $OMEGA row 4 is "0 FIX" for KA.Q.

    # ---- Residual unexplained variability ---------------------------------
    # Combined proportional + additive error. The control stream applies a
    # single inflation ratio to BOTH components for the sparse subset:
    # W = sqrt((THETA(9)/100 * IPRED)^2 + THETA(10)^2) and then
    # IF(SAMPLING == 2) W = W * THETA(11), which is algebraically the same as
    # scaling the proportional and additive SDs by 1.21 each.
    propSdIntensive <- 0.115;        label("Proportional residual error, dense sampling (fraction)")   # Table 1 RUV Prop. 11.5% (RSE 0.6%, 95% CI 11.3-11.6)
    propSdSparse    <- 0.115 * 1.21; label("Proportional residual error, sparse sampling (fraction)")  # Table 1 RUV Prop. 11.5% inflated by the Sparse (ratio) 1.21 (RSE 3.4%, 95% CI 1.13-1.29)
    addSdIntensive  <- 5.29;         label("Additive residual error, dense sampling (ng/mL)")          # Table 1 RUV Add. 5.29 ng/mL (RSE 1.5%, 95% CI 5.14-5.44)
    addSdSparse     <- 5.29 * 1.21;  label("Additive residual error, sparse sampling (ng/mL)")         # Table 1 RUV Add. 5.29 ng/mL inflated by the Sparse (ratio) 1.21
  })

  model({
    # ---- 1. Individual parameters ----------------------------------------
    tlag <- exp(ltlag)
    mtt  <- exp(lmtt + etalmtt)

    # Transit rate for a chain of NT = 3 transfers: KTR = (NT + 1)/MTT
    # (Supplemental Code $PK: NT = 3; KTR = (NT+1)/MTT).
    ktr <- 4 / mtt

    # Slow-absorption fraction. The tablet effect scales the logit; the
    # between-subject random effect is additive on the logit and is switched
    # off for the sparsely-sampled subset (control stream:
    # IF(SAMPLING == 2) FRAC_ = TVFRAC_).
    logitfslowTv <- logitfslow * ((1 - FORM_TABLET) + e_form_tablet_fslow * FORM_TABLET)
    fslow        <- 1 / (1 + exp(-(logitfslowTv + etalogitfslow * SAMPLE_INTENSIVE)))

    # Parallel absorption rates. Both carry the dose power effect and the
    # tablet ratio, and both share the single random effect on the slow rate,
    # which is switched off for the sparsely-sampled subset (control stream:
    # IF(SAMPLING == 2) KA = TVKA). The quick rate is the slow rate times the
    # estimated quick/slow ratio, which is folded into lka_quick.
    kaCov    <- e_dose_ka_slow * log(DOSE_TAK071_MG / 5) + e_form_tablet_ka_slow * FORM_TABLET +
      etalka_slow * SAMPLE_INTENSIVE
    ka_slow  <- exp(lka_slow  + kaCov)
    ka_quick <- exp(lka_quick + kaCov)

    cl <- exp(lcl + etalcl) * (AGE / 70)^e_age_cl
    vc <- exp(lvc + etalvc) * (WT / 80)^e_wt_vc

    # ---- 2. Micro-constants ----------------------------------------------
    kel <- cl / vc

    # ---- 3. ODE system ----------------------------------------------------
    # Supplemental Code $MODEL / $PK transfer rates, in order:
    #   GUT -> TR1 -> TR2 (K1T2 = K2T3 = KTR), then TR2 splits into the slow
    #   (ABS1, K3T4 = KTR * FRAC) and quick (ABS2, K3T5 = KTR * (1 - FRAC))
    #   absorption compartments, which feed CENTRAL at K4T6 = KA and
    #   K5T6 = KA * KAQ. Elimination is K6T0 = CL/VC.
    # Compartment map: depot = GUT, transit1 = TR1, transit2 = TR2,
    #   depot2 = ABS1 (slow), depot3 = ABS2 (quick), central = CENTRAL.
    d/dt(depot)    <- -ktr * depot
    d/dt(transit1) <-  ktr * depot    - ktr * transit1
    d/dt(transit2) <-  ktr * transit1 - ktr * transit2
    d/dt(depot2)   <-  ktr * fslow * transit2         - ka_slow  * depot2
    d/dt(depot3)   <-  ktr * (1 - fslow) * transit2   - ka_quick * depot3
    d/dt(central)  <-  ka_slow * depot2 + ka_quick * depot3 - kel * central

    # ---- 4. Lag time and relative bioavailability -------------------------
    alag(depot) <- tlag
    f(depot)    <- exp(lfdepot) * (DOSE_TAK071_MG / 5)^e_dose_fdepot

    # ---- 5. Observation and residual error --------------------------------
    # Amounts are in mg and vc is in L, so central/vc is mg/L = ug/mL; the
    # factor of 1000 converts to ng/mL, matching the control stream's
    # S6 = VC/1000 scaling.
    Cc <- 1000 * central / vc

    propSd <- propSdIntensive * SAMPLE_INTENSIVE + propSdSparse * (1 - SAMPLE_INTENSIVE)
    addSd  <- addSdIntensive  * SAMPLE_INTENSIVE + addSdSparse  * (1 - SAMPLE_INTENSIVE)
    Cc ~ prop(propSd) + add(addSd)
  })
}
