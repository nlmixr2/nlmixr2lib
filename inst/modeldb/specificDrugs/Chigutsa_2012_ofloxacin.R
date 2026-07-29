Chigutsa_2012_ofloxacin <- function() {
  description <- paste(
    "Two-compartment population PK model for oral ofloxacin in South",
    "African adults with multidrug-resistant tuberculosis (MDR-TB)",
    "(Chigutsa 2012; n = 65; pooled Cape Town and Durban cohorts).",
    "Savic 2007 transit-compartment absorption chain (number of",
    "transit compartments NN = 6 estimated). Total apparent oral",
    "clearance is an additive sum of two routes: a glomerular-filtration",
    "component scaling linearly with creatinine clearance (CrCl computed",
    "by a lean-body-weight modification of the Cockcroft-Gault equation;",
    "reference 68 mL/min), and an extraglomerular component (active",
    "tubular secretion + minor biliary excretion) allometrically scaled",
    "to total body weight (exponent 0.75 fixed, reference 70 kg).",
    "Central volume is allometrically scaled to lean body mass (exponent",
    "1 fixed, reference 46 kg LBM); peripheral volume to total body",
    "weight (exponent 1, reference 70 kg); intercompartmental clearance",
    "to total body weight (exponent 0.75, reference 70 kg). Mean transit",
    "time is 2.4-fold longer when ofloxacin is administered after a meal",
    "(Cape Town cohort, FED = 1) than fasted (Durban cohort, FED = 0).",
    "F is fixed at 1; residual error is combined additive (0.6 mg/L) and",
    "proportional (9.6%)."
  )
  reference <- paste(
    "Chigutsa E, Meredith S, Wiesner L, Padayatchi N, Harding J,",
    "Moodley P, Mac Kenzie WR, Weiner M, McIlleron H, Kirkpatrick CMJ",
    "(2012). Population pharmacokinetics and pharmacodynamics of",
    "ofloxacin in South African patients with multidrug-resistant",
    "tuberculosis. Antimicrobial Agents and Chemotherapy 56(7):3857-3863.",
    "doi:10.1128/AAC.00048-12.",
    sep = " "
  )
  vignette <- "Chigutsa_2012_ofloxacin"
  units    <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Baseline body weight (kg). Drives the allometric scaling of",
        "the extraglomerular clearance component (cl_nonren, exponent",
        "0.75), the peripheral volume of distribution Vp (exponent 1),",
        "and the intercompartmental clearance Q (exponent 0.75). All",
        "three scalings use 70 kg as the reference total body weight.",
        "Population median 55 kg (range 39-80 kg; Chigutsa 2012 Table 1)."
      ),
      source_name        = "WT"
    ),
    LBM = list(
      description        = "Lean body mass (lean body weight)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Baseline lean body mass (kg). Drives the allometric scaling of",
        "the central volume of distribution Vc (exponent 1, reference",
        "46 kg LBM) and is the size descriptor in the modified",
        "Cockcroft-Gault equation used to compute the CRCL covariate.",
        "The paper writes LBW (lean body weight) -- same biological",
        "quantity. Computed externally per Janmahasatian et al. 2005",
        "(Chigutsa 2012 Methods cites reference 15: Janmahasatian et",
        "al., Clin Pharmacokinet 2005; 44:1051-1065). Population median",
        "46 kg (range 32-54 kg; Chigutsa 2012 Table 1)."
      ),
      source_name        = "LBW"
    ),
    CRCL = list(
      description        = "Cockcroft-Gault creatinine clearance computed with lean body mass substituted for total body weight (raw mL/min, NOT BSA-normalized)",
      units              = "mL/min",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Stored under the canonical CRCL column per",
        "inst/references/covariate-columns.md (CRCL accepts raw mL/min",
        "when the source paper does not apply BSA normalization, with",
        "the per-model description recording the assay form -- see the",
        "CLCR alias for the Delattre 2010 raw-mL/min precedent).",
        "Chigutsa 2012 Methods reports the modified Cockcroft-Gault",
        "equation",
        "  CrCl (mL/min) = ((140 - AGE) * LBW * K) / sCr",
        "where K = 1.04 for women and 1.23 for men, AGE in years, LBW",
        "in kg, sCr in umol/L. Substituting lean body weight for total",
        "body weight in Cockcroft-Gault improved the model fit by 8",
        "OFV points (p < 0.01) compared with the standard total-body-",
        "weight formulation (Chigutsa 2012 Discussion paragraph 2).",
        "Drives the glomerular-filtration clearance component (cl_renal)",
        "linearly normalised to the population median 68 mL/min",
        "(Chigutsa 2012 Equations describing the final model; Table 1",
        "population median CrCl 109 mL/min is the unmodified total-body-",
        "weight Cockcroft-Gault and is reported only for context -- the",
        "68 mL/min reference used inside ini() is the lean-body-weight",
        "version that the final model is parameterised against)."
      ),
      source_name        = "CrCl_LBW"
    ),
    FED = list(
      description        = "Fed-vs-fasted dose-record indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (fasted; Durban cohort)",
      notes              = paste(
        "Per-dose-record fed indicator. Durban patients (FED = 0)",
        "received the daily ofloxacin dose on an empty stomach; Cape",
        "Town patients (FED = 1) received the dose after a standardised",
        "breakfast of oatmeal porridge, bread, and a cup of tea.",
        "Encoded as a multiplicative effect on mean transit time:",
        "(1 + e_fed_mtt * FED). Reference MTT (fasted) = 0.74 h",
        "(Chigutsa 2012 Table 3 'Durban mean transit time'); fed MTT =",
        "1.76 h (Chigutsa 2012 Table 3 'Cape Town mean transit time'),",
        "implying e_fed_mtt = (1.76 - 0.74) / 0.74 = 1.378 (a 2.4-fold",
        "increase under fed conditions, as reported in Chigutsa 2012",
        "Results 'Ofloxacin pharmacokinetics' final paragraph). The",
        "paper notes that the food-vs-site interpretation is partially",
        "confounded by the differing pharmacokinetic sampling schedules",
        "between the two sites, but a stochastic simulation-estimation",
        "experiment using the Durban absorption parameters with the",
        "Cape Town sampling schedule confirmed that the difference in",
        "MTT is not an artefact of the sampling schedule (bias +0.9%,",
        "precision 13%)."
      ),
      source_name        = "FOOD"
    )
  )

  covariatesDataExcluded <- list(
    HIV_POS = list(
      description        = "HIV-1 infection status indicator (1 = HIV-positive, 0 = HIV-negative)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (HIV-negative)",
      notes              = paste(
        "Screened in the Chigutsa 2012 stepwise covariate analysis but",
        "not retained: 'HIV infection was not a significant covariate",
        "on ofloxacin pharmacokinetics' (Results 'Ofloxacin",
        "pharmacokinetics' paragraph; Discussion paragraph 2). 35 / 65",
        "patients (54%) were HIV-positive, 29 of whom were on",
        "efavirenz-based antiretroviral therapy with 2 NRTIs (Results",
        "first paragraph). Documented here for provenance of the",
        "covariate screen; not referenced inside model()."
      ),
      source_name        = "HIV"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 65L,
    n_studies      = 1L,
    age_range      = "20-63 years (2.5th-97.5th percentile)",
    age_median     = "34 years",
    weight_range   = "39-80 kg (2.5th-97.5th percentile total body weight); 32-54 kg lean body weight",
    weight_median  = "55 kg total / 46 kg lean",
    sex_female_pct = 100 * 13 / 65,
    race_ethnicity = "Not explicitly reported; pooled cohort from two South African MDR-TB referral hospitals (Cape Town and Durban). All 13 female patients were from the Durban site.",
    disease_state  = paste(
      "Adults newly diagnosed with multidrug-resistant tuberculosis",
      "(MDR-TB; rifampicin- and isoniazid-resistant) at two referral",
      "hospitals in South Africa. 35 / 65 (54%) were HIV-1 co-infected,",
      "29 of whom were on antiretroviral therapy (efavirenz with 2",
      "NRTIs)."
    ),
    dose_range     = paste(
      "Once-daily oral ofloxacin 800 mg, taken in addition to other",
      "antitubercular drugs (weight-based kanamycin or amikacin,",
      "pyrazinamide, terizidone, ethionamide). 12 / 27 Durban patients",
      "additionally received 600 mg daily linezolid as part of a",
      "concurrent clinical trial (TBTC Study 30)."
    ),
    sampling       = paste(
      "Single pharmacokinetic sampling occasion at steady state (>= 1",
      "week after starting therapy). Cape Town schedule: 0.5, 3.5,",
      "5.5, 7.5, and 12 h post-dose (after breakfast). Durban",
      "schedule: 0, 1, 2, 4, 8, 11, and 24 h post-dose (fasted)."
    ),
    regions        = "South Africa (Cape Town, Western Cape Province; Durban, KwaZulu-Natal Province).",
    notes          = paste(
      "Baseline demographics per Chigutsa 2012 Table 1 (median and",
      "2.5th / 97.5th percentile). 38 patients enrolled at the Cape",
      "Town site (DP Marais Hospital) and 27 at the Durban site",
      "(CAPRISA / TBTC Study 30). All 13 female patients (20% of the",
      "cohort) were enrolled at Durban. Final parameter estimates were",
      "obtained by FOCE-I in NONMEM 7.1.2 via PsN. Covariates",
      "investigated and not retained in the final model: HIV status,",
      "total body weight on the GFR clearance component, lean body",
      "weight on the extraglomerular clearance component, and a sex",
      "effect on Vc (initially significant under a total-body-weight",
      "scaling of Vc but no longer significant under the final LBW",
      "scaling). The 2.4-fold fed-vs-fasted increase in mean transit",
      "time was independently confirmed by a simulation-estimation",
      "experiment ruling out the differing sampling schedule as the",
      "driver."
    )
  )

  ini({
    # Structural PK parameters (Chigutsa 2012 Table 3 final-model column).
    # Reference subject: total body weight 70 kg, lean body mass 46 kg,
    # creatinine clearance 68 mL/min (computed by the lean-body-weight
    # modification of Cockcroft-Gault per Chigutsa 2012 Equation set in
    # Methods 'Pharmacokinetic analysis' and Results equations describing
    # the final model). Fasted dosing (FED = 0, Durban reference cohort).

    # Clearance is decomposed into two additive routes (Chigutsa 2012 final-
    # model equation):
    #   CL/F = theta_GFR * (CrCl_i / 68) + theta_nonGFR * (WT_i / 70)^0.75
    # Both routes share a single eta on the SUM (etalcl) -- the paper's
    # Table 3 footnote b notes: 'Variability was put on the overall
    # clearance, which was the sum of the two different pathways.' This is
    # the renal-vs-non-renal-function-driven additive-CL pattern; the
    # canonical cl_renal / cl_nonren pair is used here for the
    # CrCl-driven (glomerular) versus weight-driven (extraglomerular)
    # components per parameter-names.md guidance. The 'extraglomerular'
    # route is, anatomically, predominantly active tubular secretion
    # (still renal) plus a small biliary fraction -- documented per the
    # canonical's intent (functional, not strictly anatomic). See the
    # detailed label comments below.
    lcl_renal  <- log(3.7);  label("CL_GFR (glomerular filtration) at CrCl 68 mL/min (L/h)")                       # Chigutsa 2012 Table 3: Glomerular filtration (liters/h/68 mL/min CrCl) = 3.7 (RSE 30%)
    lcl_nonren <- log(4.7);  label("CL_nonGFR (extraglomerular: tubular secretion + minor biliary) at WT 70 kg (L/h)") # Chigutsa 2012 Table 3: Extraglomerular excretion (liters/h/70 kg) = 4.7 (RSE 28%)
    lvc        <- log(52);   label("Central volume of distribution Vc at LBM 46 kg (L)")                           # Chigutsa 2012 Table 3: Central vol (liters/46 kg LBW) = 52 (RSE 20%)
    lvp        <- log(40);   label("Peripheral volume of distribution Vp at WT 70 kg (L)")                         # Chigutsa 2012 Table 3: Peripheral vol (liters/70 kg) = 40 (RSE 25%)
    lq         <- log(59);   label("Intercompartmental clearance Q at WT 70 kg (L/h)")                             # Chigutsa 2012 Table 3: Intercompartmental clearance (liters/h/70 kg) = 59 (RSE 44%)
    lmtt       <- log(0.74); label("Mean absorption transit time MTT (h, fasted reference)")                       # Chigutsa 2012 Table 3: Durban mean transit time (h) = 0.74 (RSE 18%; PPV shared with Cape Town MTT)
    lnn        <- log(6);    label("Number of absorption transit compartments NN (continuous, dimensionless)")     # Chigutsa 2012 Table 3: Number of absorption transit compartments = 6 (RSE 15%)
    lfdepot    <- fixed(log(1)); label("Oral bioavailability F (anchored at 1; not estimated)")                    # Apparent oral CL/F and V/F parameterisation; F was not separately estimated by Chigutsa 2012

    # Allometric exponents (fixed per Chigutsa 2012 Methods 'Pharmacokinetic
    # analysis' paragraph 4, citing reference (1) Anderson & Holford 2008
    # for allometric scaling).
    allo_cl_q <- fixed(0.75); label("Allometric exponent on CL_nonGFR and Q (unitless)")                           # Chigutsa 2012 Methods: 'introduced using allometric scaling (1)' -- standard 0.75 for clearances
    allo_v    <- fixed(1);    label("Allometric exponent on Vc and Vp (unitless)")                                 # Chigutsa 2012 Methods: 'allometrically scaled' -- standard 1.0 for volumes

    # Food / site effect on mean transit time.
    # Chigutsa 2012 Table 3 reports two separately estimated MTT typical
    # values (0.74 h fasted Durban and 1.76 h fed Cape Town) which Results
    # narrate as 'a 2.4-fold increase'. Encoded here as a single fasted-
    # reference MTT (lmtt) with a fractional fed-effect coefficient so the
    # food covariate is structurally explicit:
    #   e_fed_mtt = (1.76 - 0.74) / 0.74 = 1.378
    # No separate RSE is reported for this derived ratio; the underlying
    # MTT estimates carry RSE 11-18%.
    e_fed_mtt <- 1.378; label("Fractional increase in MTT under fed (FED = 1) vs fasted (FED = 0)")               # Chigutsa 2012 Table 3 ratio (Cape Town MTT 1.76) / (Durban MTT 0.74) - 1 = 1.378 (~ 2.4-fold)

    # Inter-individual variability (Chigutsa 2012 Table 3 PPV column,
    # reported as % CV; converted to log-normal variance via
    # omega^2 = log(1 + CV^2)). The eta on CL is applied to the SUM of the
    # two clearance routes per Table 3 footnote b (single shared eta on
    # overall CL/F). The correlation between eta_CL and eta_Vc is
    # reported in Table 3 as 'Covariance between random effects of
    # clearance and central vol of distribution: 0.56 (25)'. A direct
    # NONMEM OMEGA(2,1) of 0.56 is mathematically impossible given the
    # log-scale diagonal variances 0.0654 / 0.0862 (which bound any
    # covariance by sqrt(product) ~ 0.075), so the reported 0.56 is
    # interpreted here as the correlation coefficient r between the two
    # etas (the most common alternative reporting convention). The off-
    # diagonal covariance is then
    #   cov = r * sqrt(var_CL * var_Vc) = 0.56 * sqrt(0.0654 * 0.0862) = 0.0421
    # matching verification-checklist.md item C 'Correlated IIV'.
    # var(eta_CL)   = log(1 + 0.26^2)  -- Table 3 PPV CL = 26% CV (footnote b: variance on overall clearance)
    # cov(eta_CL, eta_Vc) = 0.56 * sqrt(0.06539 * 0.08618)  -- Table 3 'Covariance ... = 0.56' interpreted as correlation r
    # var(eta_Vc)   = log(1 + 0.30^2)  -- Table 3 PPV Vc = 30% CV
    etalcl + etalvc ~ c(0.06539, 0.04203, 0.08618)
    # var(eta_MTT) = log(1 + 0.54^2)  -- Table 3 PPV MTT = 54% CV (shared across fasted Durban and fed Cape Town MTT)
    etalmtt ~ 0.25596

    # Combined additive + proportional residual error (Chigutsa 2012
    # Table 3 'Additive error' and 'Proportional error'). Concentrations
    # in mg/L (= ug/mL).
    addSd  <- 0.6;   label("Additive residual error (mg/L)")               # Chigutsa 2012 Table 3: Additive error (mg/liters) = 0.6 (RSE 6%)
    propSd <- 0.096; label("Proportional residual error (fraction, CV)")   # Chigutsa 2012 Table 3: Proportional error (%) = 9.6 (RSE 9.4%)
  })

  model({
    # 1. Derived covariate terms.
    fed_mtt <- 1 + e_fed_mtt * FED                       # multiplicative food/site effect on MTT (1.0 fasted, 2.378 fed)

    # 2. Individual PK parameters. Total apparent CL is the sum of the
    #    CrCl-driven (glomerular) and weight-driven (extraglomerular)
    #    routes; a single eta is applied to the typical-value sum per
    #    Chigutsa 2012 Table 3 footnote b.
    cl_renal_typ  <- exp(lcl_renal)  * (CRCL / 68)              # linear normalisation to population-median CrCl_LBW 68 mL/min
    cl_nonren_typ <- exp(lcl_nonren) * (WT / 70)^allo_cl_q      # allometric scaling, reference 70 kg total body weight
    cl_typ        <- cl_renal_typ + cl_nonren_typ
    cl            <- cl_typ * exp(etalcl)

    vc <- exp(lvc + etalvc) * (LBM / 46)^allo_v                 # allometric scaling, reference 46 kg lean body mass
    vp <- exp(lvp)          * (WT / 70)^allo_v                  # allometric scaling, reference 70 kg total body weight
    q  <- exp(lq)           * (WT / 70)^allo_cl_q               # allometric scaling, reference 70 kg total body weight

    mtt    <- exp(lmtt + etalmtt) * fed_mtt
    nn     <- exp(lnn)                                          # NN is estimated but carries no PPV in Table 3
    fdepot <- exp(lfdepot)

    # 3. Micro-constants.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # 4. Absorption: Savic 2007 transit-compartment chain via the rxode2
    #    transit() function (analytical gamma-PDF input rate for non-
    #    integer NN). Following the Vinnard 2017 rifampicin precedent
    #    (which also uses transit() to feed a virtual depot with a fast
    #    ka collapsing the depot exponential tail), set ka >> KTR =
    #    (NN + 1)/MTT so transit() delivers the dose directly into
    #    central. f(depot) = 0 suppresses the dose-event bolus into
    #    depot; transit() reads the raw dose amount from podo(depot)
    #    regardless of f(depot), and fdepot enters the absorption
    #    process through the transit() bio argument.
    ka <- 60                                                    # >> (NN + 1) / MTT ~ 9.5 / h fasted, 4 / h fed

    d/dt(depot)       <- transit(nn, mtt, fdepot) - ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                k12 * central - k21 * peripheral1

    f(depot) <- 0

    # 5. Plasma ofloxacin concentration. Dose mg / volume L = mg/L.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
