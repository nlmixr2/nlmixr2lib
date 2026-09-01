Zhang_2016_burosumab <- function() {
  description <- paste(
    "One-compartment population PK model with first-order absorption for",
    "subcutaneous burosumab (development code KRN23, a human anti-FGF23 IgG1",
    "monoclonal antibody) in adults with X-linked hypophosphatemia, coupled to",
    "an Emax pharmacodynamic model for the change from baseline in serum",
    "phosphate (Pi) in which the potency EC50 drifts upward over treatment",
    "time. Apparent clearance and apparent central volume carry fixed",
    "allometric exponents on body weight (0.75 and 1) referenced to 70 kg, and",
    "CL/F is 1.15-fold higher at the lowest studied dose level (0.05 mg/kg)",
    "than over the linear 0.1 to 1.0 mg/kg range. The PD potency EC50,t rises",
    "sigmoidally with time since the first dose, from 1799.6 ng/mL at time 0",
    "toward an asymptote 4605.5 ng/mL higher, with a half-maximal rise at a",
    "structurally fixed 32 weeks; this captures the loss of phosphate response",
    "observed after roughly 280 days of monthly dosing."
  )
  reference <- paste(
    "Zhang X, Peyret T, Gosselin NH, Marier JF, Imel EA, Carpenter TO.",
    "Population Pharmacokinetic and Pharmacodynamic Analyses From a 4-Month",
    "Intradose Escalation and Its Subsequent 12-Month Dose Titration Studies",
    "for a Human Monoclonal Anti-FGF23 Antibody (KRN23) in Adults With",
    "X-Linked Hypophosphatemia.",
    "J Clin Pharmacol. 2016;56(4):429-438. doi:10.1002/jcph.611",
    sep = " "
  )
  vignette <- "Zhang_2016_burosumab"
  units <- list(time = "day", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight at baseline",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometric scaling referenced to 70 kg, with exponents fixed at 0.75",
        "on CL/F and 1 on Vc/F (Zhang 2016 Table 2). Body weight was the only",
        "covariate retained in the final PK model; Zhang 2016 Discussion notes",
        "that model-predicted steady-state AUC at 0.1 mg/kg is independent of",
        "body weight (Supplementary Figure S3), i.e. the mg/kg dosing plus the",
        "allometric terms together flatten the exposure-weight relationship."
      ),
      source_name        = "WT"
    ),
    DOSE_LOW = list(
      description        = "Lowest studied burosumab dose level indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (0.1, 0.3, 0.6 or 1.0 mg/kg subcutaneous burosumab)",
      notes              = paste(
        "1 = the dose record is the 0.05 mg/kg dose level, which in study",
        "KRN23-INT-001 is the first step of the four-step intradose escalation",
        "(0.05 -> 0.1 -> 0.3 -> 0.6 mg/kg); 0 = any of the 0.1, 0.3, 0.6 or",
        "1.0 mg/kg levels. Time-varying within a subject, not subject-level:",
        "a KRN23-INT-001 participant carries DOSE_LOW = 1 only during the",
        "first 28-day dosing interval and 0 thereafter. Zhang 2016 Results",
        "'KRN23 Population Pharmacokinetics' fit a separate typical CL/F per",
        "dose level (0.320 L/day at 0.05 mg/kg vs 0.271-0.294 L/day from 0.1",
        "to 1.0 mg/kg, Supplementary Table S2), concluded the PK was linear",
        "from 0.1 to 1.0 mg/kg, and therefore pooled the higher levels against",
        "a single low-dose multiplier (dMOF 37.117). Zhang 2016 Discussion",
        "attributes the higher low-dose CL/F to free circulating intact FGF23",
        "acting as a target sink that is proportionally more important at low",
        "burosumab concentrations."
      ),
      source_name        = "dose level"
    )
  )

  # Covariates Zhang 2016 screened graphically but did not retain in the final
  # model (Results 'KRN23 Population Pharmacokinetics': "Overall, there were no
  # relevant trends that were observed for the covariate effect with the
  # exception of the dose levels on CL/F"). No point estimates are reported for
  # any of them, so they are documented here rather than in covariateData.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened graphically against BSV on CL/F and Vc/F; not retained (Zhang 2016 Methods 'Population PK Methodology' and Results)."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened via box plots; not retained (Zhang 2016 Methods 'Population PK Methodology' and Results)."
    ),
    RACE_BLACK = list(
      description = "African American race indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened; not retained. The cohort was 60/62 white with a single African American participant per multiple-dose study (Zhang 2016 Table 1), so race was not estimable."
    ),
    FGF23 = list(
      description = "Baseline intact fibroblast growth factor 23",
      units       = "pg/mL",
      type        = "continuous",
      notes       = "Screened as an intrinsic factor on CL/F and Vc/F; not retained as a covariate, although Zhang 2016 Discussion invokes free intact FGF23 to explain the higher CL/F at 0.05 mg/kg that is instead encoded by DOSE_LOW."
    ),
    BALP = list(
      description = "Baseline bone alkaline phosphatase",
      units       = "ug/L",
      type        = "continuous",
      notes       = "Screened as an intrinsic factor; not retained (Zhang 2016 Methods 'Population PK Methodology')."
    )
  )

  compartmentData <- list(
    depot = list(
      analyte  = "burosumab",
      units    = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(
      analyte  = "burosumab",
      units    = "mg",
      specimen = "serum",
      verified = TRUE
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 40,
    n_studies      = 3,
    age_range      = "25-68 years (KRN23-US-02); lower bound 20 years in KRN23-INT-001 / KRN23-INT-002 (the printed upper bound is not legible in the source PDF)",
    age_median     = "48 +/- 11 years (KRN23-US-02, mean +/- SD); 42 +/- 14 (KRN23-INT-001); 42 +/- 15 (KRN23-INT-002)",
    weight_range   = "48-103 kg (KRN23-US-02); 46-122 kg (KRN23-INT-001); 51-124 kg (KRN23-INT-002)",
    weight_median  = "73 +/- 16 kg (KRN23-US-02, mean +/- SD); 70 kg median (KRN23-INT-001); 75.3 kg median (KRN23-INT-002)",
    sex_female_pct = 62.5,
    race_ethnicity = c(White = 97.5, Black = 2.5),
    disease_state  = "X-linked hypophosphatemia (XLH) in adults; baseline serum Pi 1.6-1.9 mg/dL and TmP/GFR 1.6-1.9 mg/dL, both below the normal ranges, with median intact FGF23 above the upper limit of normal",
    dose_range     = "0.05-1.0 mg/kg subcutaneous every 28 days; 0.1-1.0 mg/kg subcutaneous single dose",
    regions        = "United States and international sites (NCT00830674, NCT01340482, NCT01571596)",
    notes          = paste(
      "Baseline demographics and laboratory values are in Zhang 2016 Table 1.",
      "PK analysis set: 40 subjects (12 from the single-dose study KRN23-US-02",
      "plus 28 from the dose-escalation study KRN23-INT-001) contributing 1192",
      "measurable KRN23 concentrations. PD analysis set: the 28 multiple-dose",
      "subjects contributing 1621 serum Pi observations; 22 of them continued",
      "into the 12-month extension KRN23-INT-002 and 19 received all 16 doses.",
      "IV single-dose data were excluded because subcutaneous bioavailability",
      "was essentially complete. Sex percentages are pooled across the three",
      "study columns of Table 1 (24 female of 62 study entries is 38.7% male /",
      "61.3% female; the value here uses the 40-subject PK analysis set:",
      "6 + 19 = 25 female of 40)."
    )
  )

  ini({
    # ---- Structural PK, Zhang 2016 Table 2 (typical values for a 70-kg adult) ----
    lka <- log(0.349)
    label("First-order absorption rate constant (1/day)")                       # Zhang 2016 Table 2: Ka = 0.349 /day (90%CI 0.300-0.398, RSE 7.2%); no absorption lag time (Results 'KRN23 Population Pharmacokinetics')
    lcl <- log(0.279)
    label("Apparent clearance CL/F at 70 kg, doses >= 0.1 mg/kg (L/day)")       # Zhang 2016 Table 2: CL/F = 0.279 * (WT/70)^0.75 L/day (90%CI 0.247-0.311, RSE 5.1%)
    lvc <- log(7.17)
    label("Apparent central volume Vc/F at 70 kg (L)")                          # Zhang 2016 Table 2: Vc/F = 7.17 * (WT/70)^1 L (90%CI 6.45-7.89, RSE 5.1%)

    # ---- Allometry: exponents printed as literals with no RSE / CI, i.e. held
    # fixed at the theoretical monoclonal-antibody values rather than estimated
    # (Zhang 2016 Methods: allometric functions "as proposed for monoclonal
    # antibodies", reference 19). ----
    e_wt_cl <- fixed(0.75)
    label("Allometric exponent of body weight on CL/F (unitless)")              # Zhang 2016 Table 2: (WT/70)^0.75 on CL/F
    e_wt_vc <- fixed(1)
    label("Allometric exponent of body weight on Vc/F (unitless)")              # Zhang 2016 Table 2: (WT/70)^1 on Vc/F

    # ---- Dose-level effect on CL/F ----
    e_dose_low_cl <- 1.15
    label("Fold-increase in CL/F at the 0.05 mg/kg dose level (unitless)")      # Zhang 2016 Table 2: CL/F = 0.321 * (WT/70)^0.75 at 0.05 mg/kg = 1.15 * typical (90%CI 1.109-1.194, RSE 1.9%); 0.279 * 1.15 = 0.321

    # ---- PK between-subject variability ----
    # Reported as BSV %CV on log-normally distributed parameters (Zhang 2016
    # equation 3, theta_i = theta_TV * exp(eta_i)), so omega^2 = log(CV^2 + 1).
    # No CL/F-Vc/F correlation: the final structural model is the uncorrelated
    # FIX1cpt_mixederror form of Supplementary Table S1.
    etalka ~ 0.155378
    # Zhang 2016 Table 2: BSV on Ka = 41.0% (RSE 7.2%, shrinkage 18.8%) -> log(0.410^2 + 1) = 0.155378
    etalcl ~ 0.127006
    # Zhang 2016 Table 2: BSV on CL/F = 36.8% (RSE 5.1%, shrinkage 4.6%) -> log(0.368^2 + 1) = 0.127006
    etalvc ~ 0.088949
    # Zhang 2016 Table 2: BSV on Vc/F = 30.5% (RSE 5.1%, shrinkage 7.4%) -> log(0.305^2 + 1) = 0.088949

    # ---- PK residual error (combined additive + proportional) ----
    propSd <- 0.218
    label("Proportional residual error on serum burosumab (fraction)")          # Zhang 2016 Table 2: Error prop = 21.8% (90%CI 7.7%-35.8%, RSE 32.9%)
    addSd <- 0.099
    label("Additive residual error on serum burosumab (ng/mL)")                 # Zhang 2016 Table 2: Error additive = 0.099 ng/mL (90%CI 0.054-0.144, RSE 23.3%)

    # ---- PK-PD, Zhang 2016 Table 3 ----
    # Structure (Zhang 2016 equations 9 and 10, page 435):
    #   dPi     = simC_KRN23 * Emax / (EC50,t + simC_KRN23)
    #   EC50,t  = tvEC50 + a * t^g / (32^g + t^g),  t in weeks after first dose
    lemax <- fixed(log(1.5))
    label("Maximum effect on serum Pi change from baseline (mg/dL)")            # Zhang 2016 Table 3: Emax = 1.5 mg/dL (FIX). Note to Table 3: "Emax value was fixed to achieve the covariance step and obtain precision of the parameters"; the fixed value came from data exploration and previous run results (Results 'Population PK-PD Modeling')
    ltvec50 <- log(1799.6)
    label("EC50,t at time zero, tvEC50 (ng/mL)")                                # Zhang 2016 Table 3: tvEC50 = 1799.6 ng/mL (RSE 15.9%)
    lssec50 <- log(4605.5)
    label("Maximum increase of EC50,t above tvEC50 (ng/mL)")                    # Zhang 2016 Table 3: a = 4605.5 ng/mL/week (RSE 16.8%), labelled "maximum rate of increase of EC50,t"; in equation 10 it multiplies a dimensionless saturating function of time, so it acts as the asymptotic increment of EC50,t (see the model file note below)
    lhill_ec50t <- log(2.88)
    label("Hill coefficient of the EC50,t rise over time (unitless)")           # Zhang 2016 Table 3: g = 2.88 (RSE 17.1%), the exponent on t in equation 10
    lt50_ec50t <- fixed(log(32))
    label("Time to half-maximal rise of EC50,t (weeks)")                        # Zhang 2016 equation 10: the denominator is (32^g + t^g), i.e. a structural constant of 32 weeks printed only inside the equation and absent from Table 3; see the model file note below

    # ---- PK-PD between-subject variability ----
    # Zhang 2016 Table 3 places the single 73.8% BSV on the composite row
    # "EC50,t", with tvEC50 / a / g as its indented components, so the random
    # effect scales the whole EC50,t(t) curve rather than the intercept alone.
    etaltvec50 ~ 0.434793
    # Zhang 2016 Table 3: BSV on EC50,t = 73.8% (RSE 14.8%) -> log(0.738^2 + 1) = 0.434793

    # ---- PK-PD residual error ----
    addSd_dPi <- 0.310
    label("Additive residual error on serum Pi change from baseline (mg/dL)")   # Zhang 2016 Table 3: Residual additive error = 0.310 mg/dL (RSE 1.8%)
  })

  model({
    # ---- 1. Individual PK parameters ----
    # DOSE_LOW is a per-dose-record indicator (1 only for 0.05 mg/kg records).
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl * e_dose_low_cl^DOSE_LOW
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc

    kel <- cl / vc

    # ---- 2. One-compartment disposition with first-order SC absorption ----
    # Subcutaneous bioavailability was essentially complete versus IV
    # (Zhang 2016 Introduction and Discussion), so F is 1 and CL/F, Vc/F are
    # reported as apparent parameters; no bioavailability term is estimated.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # ---- 3. Serum burosumab concentration ----
    # central is in mg and vc in L, so central/vc is mg/L = ug/mL; the assay and
    # every reported PK / PD parameter are in ng/mL, hence the factor of 1000.
    Cc <- 1000 * central / vc

    # ---- 4. Time-varying potency, Zhang 2016 equation 10 ----
    # Equation 10 reads EC50,t = tvEC50 + a * t^g / (32^g + t^g) with t in
    # WEEKS after the first dose; the model time unit is days, hence t / 7.
    #
    # Two features of equation 10 are worth flagging because neither is
    # recoverable from Table 3 alone:
    #   (a) the 32 in the denominator is a structural constant printed only
    #       inside the equation. It is carried here as the fixed parameter
    #       t50_ec50t so it is visible and source-traced rather than buried as
    #       a magic number.
    #   (b) "a" is tabulated with units of ng/mL/week and described as the
    #       "maximum rate of increase", but in equation 10 it multiplies a
    #       dimensionless saturating function of time, so it is the asymptotic
    #       INCREMENT of EC50,t in ng/mL. The paper's own published EC50,t
    #       values confirm this reading: with tvEC50 = 1799.6, a = 4605.5,
    #       g = 2.88 and the 32-week constant, equation 10 returns 4102 ng/mL
    #       at week 32, 5999 ng/mL at week 72 and 6098 ng/mL at 560 days,
    #       matching the three values quoted in Results and Discussion.
    #
    # The 73.8% BSV attaches to the composite EC50,t (Table 3 row hierarchy),
    # so the same eta scales both components: this is algebraically identical
    # to EC50,t(t) * exp(eta) while keeping the parameters mu-referenced.
    tweek <- t / 7
    tvec50 <- exp(ltvec50 + etaltvec50)
    ssec50 <- exp(lssec50 + etaltvec50)
    hill_ec50t <- exp(lhill_ec50t)
    t50_ec50t <- exp(lt50_ec50t)

    ec50t <- tvec50 + ssec50 * tweek^hill_ec50t /
      (t50_ec50t^hill_ec50t + tweek^hill_ec50t)

    # ---- 5. Emax effect on serum phosphate, Zhang 2016 equation 9 ----
    # dPi is the change from baseline in serum Pi, with baseline defined as the
    # first predose Pi level (day 0 of KRN23-INT-001). The model predicts the
    # change, not the absolute concentration; add the observed baseline to
    # recover absolute serum Pi.
    emax <- exp(lemax)
    dPi <- emax * Cc / (ec50t + Cc)

    # ---- 6. Observations ----
    Cc ~ add(addSd) + prop(propSd)
    dPi ~ add(addSd_dPi)
  })
}
