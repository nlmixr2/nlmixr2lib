Valle_2005_exemestane <- function() {
  description <- "Three-compartment population PK with first-order absorption + lag time, coupled to an indirect-response PD model on plasma estrone sulphate (E1S), for oral exemestane (25 mg single dose) in healthy postmenopausal women. Crossover study comparing a sugar-coated tablet (SCT) under fasting versus an extemporaneous tablet-suspended-in-water suspension under fasting versus a SCT taken after a standard high-fat breakfast. Disposition is independent of formulation and food; absorption rate ka and apparent bioavailability F depend on formulation (suspension: ka 7.6 vs SCT 2.35 1/h, F 1.2x) and on the high-fat meal (ka 1.13 1/h, F 1.6x). Exemestane inhibits E1S synthesis via a sigmoid Imax function with IC50 22.1 pg/mL and Hill coefficient 1.73."
  reference   <- "Valle M, Di Salle E, Jannuzzo MG, Poggesi I, Rocchetti M, Spinelli R, Verotta D. A predictive model for exemestane pharmacokinetics/pharmacodynamics incorporating the effect of food and formulation. Br J Clin Pharmacol. 2005;59(3):355-364. doi:10.1111/j.1365-2125.2005.02335.x"
  vignette    <- "Valle_2005_exemestane"
  units       <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  paper_specific_compartments <- c("e1s")

  covariateData <- list(
    FORM_SUSPENSION = list(
      description        = "Extemporaneous tablet-suspended-in-water formulation indicator (1 = exemestane delivered as an extemporaneously-prepared suspension administered with 180 mL of vehicle, 0 = sugar-coated tablet swallowed whole with 180 mL of water).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (sugar-coated tablet (SCT); the typical-value reference in Valle 2005 Tables 2-3).",
      notes              = "Per-dose-occasion indicator: each subject received the SCT (fasting and after high-fat breakfast) and the suspension (fasting) in a randomized three-period 3x3 Latin-square crossover. Multiplicative effect on absorption rate ka (suspension/SCT-fasting ratio 7.6/2.35 = 3.234x; Valle 2005 Table 2) and on apparent bioavailability F (suspension/SCT-fasting ratio = 1.2x; Valle 2005 Results). V is intrinsic and shared across formulations -- the apparent V/F differences in Table 2 (1360 vs 1120 vs 844 L for SCT-fasting / suspension-fasting / SCT-food) collapse to V multiplied by 1/F at the treatment-specific F.",
      source_name        = "Treatment 2 (suspension after fasting) vs Treatment 1 (SCT after fasting)"
    ),
    FED_HIGHFAT = list(
      description        = "High-fat-meal-at-dosing indicator (1 = sugar-coated tablet given 15 min after a standard high-fat breakfast, 0 = dose taken under overnight-fasting conditions).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (overnight-fasting condition; the typical-value reference in Valle 2005 Tables 2-3).",
      notes              = "Per-dose-occasion indicator. Valle 2005 Treatment 3 = 25 mg SCT taken 15 min after a 'standard high-fat breakfast'; Treatments 1 (SCT) and 2 (suspension) were given under overnight-fasting conditions. Multiplicative effect on ka (SCT-food/SCT-fasting ratio 1.13/2.35 = 0.481x; Valle 2005 Table 2) and on apparent bioavailability F (SCT-food/SCT-fasting ratio = 1.6x; Valle 2005 Results). The food effect was only estimated on the SCT arm (no suspension-after-food condition was run); the model applies the food effect generically per dose record so SCT-after-food is the documented use case.",
      source_name        = "Treatment 3 (SCT after high-fat breakfast) vs Treatment 1 (SCT after fasting)"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 12,
    n_studies      = 1,
    age_range      = "45-68 years",
    age_median     = "55 years (mean)",
    weight_range   = "46-66 kg",
    weight_median  = "54.9 kg (mean)",
    sex_female_pct = 100,
    race_ethnicity = "Not reported in Valle 2005 Methods; trial conducted at Pontchaillou Hospital (Rennes, France).",
    disease_state  = "Healthy postmenopausal women",
    dose_range     = "25 mg single oral dose per period; three treatment periods per subject in a 3x3 Latin-square crossover; 4-5 week washout between periods.",
    regions        = "France (Pontchaillou Hospital, Rennes)",
    notes          = "Open, three-period, randomized, 3x3 Latin-square crossover. Treatments: (1) sugar-coated tablet 25 mg under fasting, (2) extemporaneous suspension 25 mg under fasting, (3) sugar-coated tablet 25 mg 15 min after a standard high-fat breakfast. Exemestane plasma sampling at 0.25, 0.5, 1, 1.5, 2, 4, 6, 8, 12, 16, 24, 48, 72, 120, and 168 h post-dose (15 points). E1S sampling at 0, 24, 48, 72, 120, 168, and 336 h post-dose. Estimation by NONMEM V FOCE INTERACTION; population PK fitted first with empirical-Bayes individual estimates, then population PD fitted conditional on the individual PK estimates."
  )

  ini({
    # ------------------ PK structural parameters (reference = SCT under fasting) --
    # 3-compartment mammillary disposition (depot -> central <-> peripheral1;
    # central <-> peripheral2) with first-order absorption and absorption lag.
    # Parameterised as rate constants (kel, k12, k21, k13, k31) + central
    # volume V to preserve the paper's reported IIV structure (CV on V/F, ka,
    # and k independently) without cross-parameter covariance translation.
    # Precedent: rate-constant parameterisation with lkel + lvc is used by
    # Reijers_2016_trastuzumab.R and Krause_2017_selexipag.R.
    lka   <- log(2.35);   label("First-order absorption rate constant ka under SCT/fasting reference (1/h)")           # Valle 2005 Table 2: ka1 = 2.35 1/h (SCT/fasting), RSE 25%
    lvc   <- log(1360);   label("Apparent central volume V/F under SCT/fasting reference (L)")                          # Valle 2005 Table 2: V/F = 1360 L (SCT/fasting), RSE 10%; treated as intrinsic V (constant across treatments) with F-treatment effects below
    ltlag <- log(0.21);   label("Absorption lag time t_lag (h)")                                                        # Valle 2005 Table 2: Lag = 0.21 h (common to all three treatments), RSE 10%
    lkel  <- log(0.738);  label("First-order elimination rate constant k from central (1/h)")                           # Valle 2005 Table 3: k = 0.738 1/h, RSE 7%
    lk12  <- log(0.454);  label("Central -> peripheral1 transfer rate k12 (1/h)")                                       # Valle 2005 Table 3: k12 = 0.454 1/h, RSE 19%
    lk21  <- log(0.158);  label("Peripheral1 -> central transfer rate k21 (1/h)")                                       # Valle 2005 Table 3: k21 = 0.158 1/h, RSE 6%
    lk13  <- log(0.174);  label("Central -> peripheral2 transfer rate k13 (1/h)")                                       # Valle 2005 Table 3: k13 = 0.174 1/h, RSE 9%
    lk31  <- log(0.016);  label("Peripheral2 -> central transfer rate k31 (1/h)")                                       # Valle 2005 Table 3: k31 = 0.016 1/h, RSE 10%

    # ------------------ Formulation and food effects (log-scale shifts) ----------
    # All four effects are multiplicative on the linear scale (log-additive on
    # the log scale). Reference treatment is SCT under fasting (FORM_SUSPENSION
    # = 0, FED_HIGHFAT = 0); the other two treatments shift ka and F as
    # reported in Table 2 / Results.
    e_form_suspension_ka     <- log(7.6 / 2.35);   label("Effect of suspension formulation on ka (log scale)")           # Valle 2005 Table 2: ka1 suspension/fasting = 7.6 (RSE 21%) vs SCT/fasting reference 2.35 -> log(7.6/2.35) = log(3.234)
    e_fed_highfat_ka         <- log(1.13 / 2.35);  label("Effect of high-fat meal on ka (log scale)")                    # Valle 2005 Table 2: ka1 SCT/food = 1.13 (RSE 35%) vs SCT/fasting reference 2.35 -> log(1.13/2.35) = log(0.481)
    e_form_suspension_fdepot <- log(1.2);          label("Effect of suspension formulation on apparent F (log scale)")  # Valle 2005 Results: suspension formulation showed a 1.2x higher F than the SCT after fasting -> log(1.2)
    e_fed_highfat_fdepot     <- log(1.6);          label("Effect of high-fat meal on apparent F (log scale)")            # Valle 2005 Results: food intake increased F by a factor of 1.6 -> log(1.6)

    # ------------------ PD structural parameters (indirect-response on E1S) ------
    # The PD model is an indirect-response model in which E1S plasma
    # concentration is determined by a zero-order synthesis rate ks (derived
    # as ks = baseline * kout in model()) and a first-order elimination rate
    # constant kout, and where exemestane inhibits E1S synthesis with a
    # sigmoid Imax function: inh = Cc^hill / (ic50^hill + Cc^hill). The
    # alternative mechanism-based model (equations 5-7) was not retained
    # because the OFV decrease was less than 3.9 points (Valle 2005 Results).
    lrbase <- log(203);   label("Baseline plasma E1S concentration (pg/mL)")                                            # Valle 2005 Table 4: Baseline = 203 pg/mL, RSE 8% (matches ks/kout = 6.5 / 0.032 reported in the Abstract)
    lkout  <- log(0.032); label("First-order elimination rate constant of E1S, kout (1/h)")                             # Valle 2005 Table 4: ko = 0.032 1/h, RSE 9% (terminal half-life ~21.7 h)
    lic50  <- log(22.1);  label("Exemestane concentration giving 50% inhibition of E1S synthesis IC50 (pg/mL)")          # Valle 2005 Table 4: C50 = 22.1 pg/mL, RSE 10%
    lhill  <- log(1.73);  label("Hill coefficient on inhibition of E1S synthesis (unitless)")                            # Valle 2005 Table 4: g = 1.73, RSE 12% (Hill exponent; IIV not estimated -- N.E. in Table 4)

    # ------------------ IIV (variances on log scale) -----------------------------
    # omega^2 = log(CV^2 + 1) for log-normal IIV. The paper reports IIV CV%
    # for ka, V/F, and k in the PK Results paragraph and for Baseline, ko,
    # C50 in Table 4 (with RSE on each omega in parenthesis). No IIV on
    # k12, k21, k13, k31, lag, or Hill per Valle 2005 Results paragraph
    # ("No significant interindividual variability in the other parameters
    # describing the plasma disposition of the drug was supported by the
    # data") and Table 4 N.E. for g.
    etalka    ~ 0.61308  # log(1 + 0.92^2); Valle 2005 Results: ka CV 92%
    etalvc    ~ 0.04727  # log(1 + 0.22^2); Valle 2005 Results: V/F CV 22%
    etalkel   ~ 0.01203  # log(1 + 0.11^2); Valle 2005 Results: k CV 11%
    etalrbase ~ 0.27290  # log(1 + 0.56^2); Valle 2005 Table 4: Baseline CV 56% (RSE 19%)
    etalkout  ~ 0.15533  # log(1 + 0.41^2); Valle 2005 Table 4: ko CV 41% (RSE 47%)
    etalic50  ~ 0.19960  # log(1 + 0.47^2); Valle 2005 Table 4: C50 CV 47% (RSE 32%)

    # ------------------ Residual error -------------------------------------------
    # Exemestane: combined error model with proportional 4% estimated and
    # additive part FIXED to the SD of the analytical technique (Valle 2005
    # Results). The paper does not report the numeric SD value; we
    # approximate it as the assay's near-LLOQ SD = LLOQ x inter-assay CV =
    # 13.5 pg/mL x 0.177 = 2.39 pg/mL ~= 0.0024 ng/mL. This is a documented
    # approximation -- see the vignette's Errata / Assumptions section.
    # E1S: proportional 35% (Valle 2005 Results, PD residual variability).
    propSd     <- 0.04;          label("Proportional residual SD on exemestane Cc (fraction)")                          # Valle 2005 Results: proportional part of exemestane residual error = 4%
    addSd      <- fixed(0.0024); label("Additive residual SD on exemestane Cc (ng/mL) -- assay-SD approximation")        # Valle 2005 Results: fixed to SD of analytical technique; numeric SD not reported, approximated as LLOQ (13.5 pg/mL) x inter-assay CV (17.7%) ~= 0.0024 ng/mL (see vignette Errata)
    propSd_e1s <- 0.35;          label("Proportional residual SD on E1S (fraction)")                                    # Valle 2005 Results: PD residual variability = 35%
  })

  model({
    # 1. Individual PK parameters (rate-constant parameterisation).
    # Formulation and food effects are log-additive on ka and on the F factor
    # applied to the depot dose; V is intrinsic (no covariate effects).
    ka     <- exp(lka  + etalka  + e_form_suspension_ka * FORM_SUSPENSION + e_fed_highfat_ka * FED_HIGHFAT)
    vc     <- exp(lvc  + etalvc)
    kel    <- exp(lkel + etalkel)
    k12    <- exp(lk12)
    k21    <- exp(lk21)
    k13    <- exp(lk13)
    k31    <- exp(lk31)
    tlag   <- exp(ltlag)
    fdepot <- exp(e_form_suspension_fdepot * FORM_SUSPENSION + e_fed_highfat_fdepot * FED_HIGHFAT)

    # 2. Individual PD parameters and derived synthesis rate.
    rbase <- exp(lrbase + etalrbase)
    kout  <- exp(lkout  + etalkout)
    ic50  <- exp(lic50  + etalic50)
    hill  <- exp(lhill)
    kin   <- rbase * kout  # zero-order E1S synthesis rate (pg/mL/h); ks = baseline * ko per Valle 2005

    # 3. ODE system: 3-compartment mammillary PK + indirect-response PD on E1S.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1 - k13 * central + k31 * peripheral2
    d/dt(peripheral1) <-  k12 * central - k21 * peripheral1
    d/dt(peripheral2) <-  k13 * central - k31 * peripheral2

    # Exemestane plasma concentration: central is in mg, vc is in L, so
    # central/vc is in mg/L; scale to ng/mL by *1000.
    Cc      <- central / vc * 1000
    # Convert Cc to pg/mL for the PD inhibition equation so it shares units
    # with IC50.
    Cc_pgml <- Cc * 1000

    # E1S indirect-response equation. E1S starts at its endogenous steady-state
    # baseline rbase = ks/kout and returns there as exemestane is cleared.
    inh         <- Cc_pgml^hill / (ic50^hill + Cc_pgml^hill)
    d/dt(e1s)   <- kin * (1 - inh) - kout * e1s
    e1s(0)      <- rbase

    # 4. Absorption lag and bioavailability on the depot dose.
    alag(depot) <- tlag
    f(depot)    <- fdepot

    # 5. Observations and residual error.
    Cc  ~ add(addSd) + prop(propSd)
    e1s ~ prop(propSd_e1s)
  })
}
