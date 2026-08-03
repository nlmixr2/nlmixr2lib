Gibbs_2011_amg221 <- function() {
  description <- paste(
    "Semi-mechanistic population PK/PD model of subcutaneous adipose",
    "11-beta-hydroxysteroid dehydrogenase type 1 (11-beta-HSD1) activity",
    "after single oral doses of AMG 221, a selective 11-beta-HSD1 inhibitor,",
    "in healthy obese adults (Gibbs 2011). Two-compartment PK with",
    "first-order absorption and lag time; dose-dependent bioavailability",
    "reduced at 3 mg vs 30/100 mg; adipose-tissue effect compartment linked",
    "to plasma by a first-order equilibration rate constant (keo) and to",
    "the tissue-concentration observation by an adipose-plasma density",
    "correction factor (kpp); direct-effect Imax inhibitory model on",
    "11-beta-HSD1 activity driven by the effect-site plasma-equivalent",
    "AMG 221 concentration."
  )
  reference <- paste(
    "Gibbs JP, Emery MG, McCaffery I, Smith B, Gibbs MA, Akrami A, Rossi J,",
    "Paweletz K, Gastonguay MR, Bautista E, Wang M, Perfetti R, Daniels O.",
    "Population pharmacokinetic/pharmacodynamic model of subcutaneous",
    "adipose 11-beta-hydroxysteroid dehydrogenase type 1 (11-beta-HSD1)",
    "activity after oral administration of AMG 221, a selective 11-beta-HSD1",
    "inhibitor. J Clin Pharmacol. 2011;51(6):830-841.",
    "doi:10.1177/0091270010374470"
  )
  vignette <- "Gibbs_2011_amg221"
  units <- list(time = "h", dosing = "mg", concentration = "ng/mL")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix. analyte/specimen proposed by a local model from the
  # model description; units derived from the units block. verified = FALSE
  # means NOT checked against the source paper.
  compartmentData <- list(
    depot       = list(analyte = "AMG 221", units = "mg", specimen = "administration site", verified = FALSE),
    central     = list(analyte = "AMG 221", units = "mg", specimen = "plasma", verified = FALSE),
    peripheral1 = list(analyte = "AMG 221", units = "mg", specimen = "plasma", verified = FALSE),
    effect      = list(analyte = "AMG 221", units = "mg", specimen = "not applicable", verified = FALSE)
  )

  covariateData <- list(
    DOSE_LOW_AMG221 = list(
      description        = "Low-dose AMG 221 indicator (1 = 3 mg oral AMG 221 dose record; 0 = 30 or 100 mg oral AMG 221 dose record).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (30 or 100 mg oral AMG 221, where F1 = 1)",
      notes              = paste(
        "Per-dose-record indicator that scales relative bioavailability F1",
        "downward at the 3 mg dose level (Gibbs 2011 Table III: F1 = 0.546",
        "at 3 mg vs F1 = 1 at 30/100 mg; footnote a). Gibbs 2011 Discussion",
        "attributes the reduced 3 mg bioavailability to a possible",
        "high-affinity intestinal-metabolism / transport process saturating",
        "at higher doses (paragraphs on Caco-2 permeability and CYP3A",
        "metabolism with Km > 100 uM); the paper flags an alternative",
        "Michaelis-Menten dose-F structure that was not fit because only",
        "three discrete dose levels were tested. Apply at the dose record;",
        "observation rows inherit the indicator from the preceding dose."
      ),
      source_name        = "derived from administered dose"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 55L,
    n_studies      = 1L,
    age_range      = "18-45 years",
    age_median     = "30 years (mean 30 +/- 7 SD)",
    weight_range   = "means 104 kg (placebo) and 106 kg (AMG 221) with SD +/- 13-14 kg",
    weight_median  = "approximately 105 kg (mean +/- 13 SD)",
    sex_female_pct = 18,
    race_ethnicity = c(
      "White or Caucasian"           = 49,
      "Black or African American"    = 15,
      "Hispanic or Latino"           = 33,
      "American Indian or Alaska Native" = 2,
      "Other"                        = 2
    ),
    disease_state  = "Healthy obese adults (BMI 29.1-42.8 kg/m^2, mean 34.0 +/- 3.2 SD).",
    dose_range     = "Single oral 3, 30, or 100 mg AMG 221 as a suspension (or matching placebo).",
    regions        = "Not stated; Amgen phase 1 sponsor",
    notes          = paste(
      "n = 55 total: 44 AMG 221 (3 mg n=20, 30 mg n=12, 100 mg n=12) + 11",
      "placebo. Phase 1, randomized, placebo-controlled, double-blind,",
      "exploratory design (Gibbs 2011 Methods). Rich plasma sampling: 13",
      "samples per subject over 24 h (predose, 0.25, 0.5, 1, 1.5, 2, 3, 4,",
      "6, 8, 12, 16, 24 h). Adipose biopsies at baseline (day -1) and on",
      "day 1 at 1.5, 4, 8, or 24 h post dose (24 h only for the 3 mg cohort).",
      "11-beta-HSD1 activity measured ex vivo using prednisone->prednisolone",
      "conversion (Gibbs 2011 Analytical Methods). Baseline demographics",
      "from Gibbs 2011 Table I."
    )
  )

  ini({
    # =================================================================
    # Structural PK parameters (Gibbs 2011 Table III, p. 205)
    # Two-compartment open model with first-order absorption + lag.
    # Time unit: h. Dose unit: mg. Concentration unit: ng/mL.
    # =================================================================
    lka   <- log(2.90);   label("Absorption rate constant ka (1/h)")            # Table III: ka = 2.90 1/h, SE 0.489 (17% RSE)
    lcl   <- log(16.8);   label("Apparent oral clearance CLp/F (L/h)")          # Table III: CLp = 16.8 L/h, SE 1.22 (7.3% RSE)
    lvc   <- log(71.8);   label("Apparent central volume V2/F (L)")             # Table III: V2 = 71.8 L, SE 5.54 (7.7% RSE)
    lvp   <- log(83.4);   label("Apparent peripheral volume V3/F (L)")          # Table III: V3 = 83.4 L, SE 8.08 (9.7% RSE)
    lq    <- log(16.8);   label("Apparent inter-compartmental clearance Qp/F (L/h)")  # Table III: Qp = 16.8 L/h, SE 1.25 (7.4% RSE)
    ltlag <- log(0.231);  label("Absorption lag time ALAG1 (h)")                 # Table III: ALAG1 = 0.231 h, SE 0.00903 (3.9% RSE)

    # Bioavailability: F1 = 1 for 30 and 100 mg (fixed anchor); F1 = 0.546
    # at 3 mg (estimated). Encoded via the DOSE_LOW_AMG221 indicator.
    lfdepot                  <- fixed(log(1));  label("Relative bioavailability F1 at 30/100 mg reference (unitless)")
    e_dose_low_amg221_fdepot <- log(0.546);     label("log F1(3 mg) / F1(30 or 100 mg) ratio (unitless)")   # Table III + footnote a: F1(3 mg) = 0.546, SE 0.0415 (7.6% RSE); F1(30/100 mg) fixed at 1

    # PK IIV (exponential structure per Gibbs 2011 Equation 6; omega^2
    # reported directly in Table III). Partial covariance block: CLp-V2;
    # ka is independent; ALAG1 is independent (Table III row lists an
    # ALAG1 IIV of 0.0704 without SE -- see vignette Assumptions).
    # F1, V3, and Qp had IIV fixed at 0 (Table III footnote b).
    etalcl + etalvc ~ c(0.125, 0.0557, 0.108)   # Table III: omega^2(CLp) = 0.125 (~35% CV), cov(CLp,V2) = 0.0557 (SE 0.0263, 47% RSE), omega^2(V2) = 0.108 (~33% CV)
    etalka  ~ 1.09                              # Table III: omega^2(ka)   = 1.09  (~104% CV), SE 0.235 (22% RSE)
    etaltlag ~ 0.0704                           # Table III: omega^2(ALAG1) = 0.0704 (~27% CV; SE column blank in the Table III ALAG1 row)

    # Residual error for plasma AMG 221. Gibbs 2011 Equation 7: exponential
    # residual error (C = Cpred * exp(eps), eps ~ N(0, sigma^2)); encoded
    # here as ~ lnorm(expSd) where expSd = sqrt(sigma^2) is the SD on the
    # log scale (approximately the CV for small sigma^2).
    expSd <- sqrt(0.0398); label("Log-scale residual SD on plasma AMG 221 Cc (unitless)")  # Table III: sigma^2(plasma AMG 221) = 0.0398 (~20% CV)

    # =================================================================
    # PK/PD parameters (Gibbs 2011 Table IV, p. 206)
    # Sequential PK/PD: post-hoc PK parameters used as a forcing function
    # for the adipose effect-compartment + Imax 11-beta-HSD1 model.
    # =================================================================
    lke0    <- log(0.220);   label("Effect-compartment equilibration rate constant keo (1/h)")            # Table IV: keo = 0.220 1/h, SE 0.0209 (9.5% RSE)
    lkpp    <- log(1.36);    label("Adipose-plasma density correction factor kpp (mL/g)")                  # Table IV: kpp = 1.36 mL/g, SE 0.0657 (4.8% RSE); Gibbs 2011 Discussion notes 1/kpp = 0.735 g/mL approximates the density of fat
    lrbase  <- log(755);     label("Baseline 11-beta-HSD1 enzyme activity E0 (pmol/mg)")                    # Table IV: E0 = 755 pmol/mg, SE 61.2 (8.1% RSE)
    limax   <- log(0.975);   label("Maximal fractional inhibition of 11-beta-HSD1 activity Imax (unitless, 0-1)")   # Table IV: Imax = 0.975, SE 0.00259 (0.3% RSE)
    lic50   <- log(1.19);    label("Plasma-equivalent AMG 221 concentration at half-max inhibition IC50 (ng/mL)")   # Table IV: IC50 = 1.19 ng/mL, SE 0.121 (10% RSE)

    # PK/PD IIV (Table IV). kpp, Imax, and IC50 had IIV fixed at 0
    # (Table IV footnote a).
    etalke0    ~ 0.158    # Table IV: omega^2(keo) = 0.158 (~40% CV), SE 0.0566 (36% RSE)
    etalrbase  ~ 0.120    # Table IV: omega^2(E0)  = 0.120 (~35% CV), SE 0.0946 (79% RSE)

    # Residual errors for adipose AMG 221 tissue concentration and 11-beta-HSD1
    # activity (Gibbs 2011 Equation 7 exponential form; encoded as lnorm).
    expSd_Ccadipose <- sqrt(0.0271); label("Log-scale residual SD on adipose AMG 221 Ccadipose (unitless)")   # Table IV: sigma^2(AMG 221 tissue) = 0.0271 (~16% CV), SE 0.0153 (57% RSE)
    expSd_activity  <- sqrt(0.0984); label("Log-scale residual SD on 11-beta-HSD1 activity (unitless)")       # Table IV: sigma^2(11-beta-HSD1)   = 0.0984 (~31% CV), SE 0.0353 (36% RSE)
  })

  model({
    # =================================================================
    # Individual parameters
    # =================================================================
    ka    <- exp(lka  + etalka)
    cl    <- exp(lcl  + etalcl)
    vc    <- exp(lvc  + etalvc)
    vp    <- exp(lvp)
    q     <- exp(lq)
    tlag  <- exp(ltlag + etaltlag)

    # Dose-dependent bioavailability: F1 = 1 for reference (30 or 100 mg
    # oral AMG 221; DOSE_LOW_AMG221 = 0) and F1 = 0.546 for 3 mg
    # (DOSE_LOW_AMG221 = 1). Gibbs 2011 Table III + footnote a.
    fdepot <- exp(lfdepot + e_dose_low_amg221_fdepot * DOSE_LOW_AMG221)

    keo   <- exp(lke0   + etalke0)
    kpp   <- exp(lkpp)
    rbase <- exp(lrbase + etalrbase)
    imax  <- exp(limax)
    ic50  <- exp(lic50)

    # =================================================================
    # Two-compartment PK disposition with effect compartment
    # =================================================================
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                                k12 * central - k21 * peripheral1

    # Adipose effect compartment. First-order equilibration with the
    # central compartment at rate keo (Gibbs 2011 Table IV keo). The
    # implicit effect-compartment volume equals the central Vc so that
    # the effect-site plasma-equivalent concentration is A(effect) / vc.
    d/dt(effect)      <-  keo * central - keo * effect

    f(depot)    <- fdepot
    alag(depot) <- tlag

    # =================================================================
    # Observations
    # =================================================================
    # Plasma AMG 221 concentration (ng/mL): dose is in mg, vc is in L, so
    # central/vc has units of mg/L = 1000 ng/mL.
    Cc <- central / vc * 1000

    # Effect-site plasma-equivalent AMG 221 concentration (ng/mL), used as
    # the driver for the Imax inhibition and for the tissue-concentration
    # observation via the kpp density correction.
    Ce_plasma_eq <- effect / vc * 1000

    # Adipose tissue AMG 221 concentration (ng/g): plasma-equivalent
    # concentration times the adipose-plasma density correction factor.
    # kpp has units of mL/g so (ng/mL) * (mL/g) = ng/g.
    Ccadipose <- Ce_plasma_eq * kpp

    # 11-beta-HSD1 activity (pmol/mg) with direct-effect Imax inhibition
    # driven by the effect-site plasma-equivalent concentration (Gibbs
    # 2011 Table IV; IC50 is on the plasma concentration scale).
    activity  <- rbase * (1 - imax * Ce_plasma_eq / (ic50 + Ce_plasma_eq))

    Cc        ~ lnorm(expSd)
    Ccadipose ~ lnorm(expSd_Ccadipose)
    activity  ~ lnorm(expSd_activity)
  })
}
