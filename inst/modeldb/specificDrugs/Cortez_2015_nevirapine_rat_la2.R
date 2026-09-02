Cortez_2015_nevirapine_rat_la2 <- function() {
  description <- paste(
    "Preclinical (rat).",
    "One-compartment population PK model with first-order absorption for the",
    "NVP LA-2 long-acting injectable nevirapine microparticle suspension",
    "(recrystallised nevirapine cuboids, 239 um median diameter,",
    "fluidised-bed coated with poly(D,L-lactide) to a 2 um membrane) after a",
    "single subcutaneous dose to Sprague-Dawley rats (study 2, JNN00030).",
    "Fitted by naive pooling in Phoenix NLME; no between-subject variability",
    "was estimated. IMPORTANT: Cortez 2015 describes this arm's elimination",
    "as SATURABLE, but Table 2 reports only the intrinsic clearance Vm/Km =",
    "48.2 L/day and never reports Vm or Km individually, so the",
    "Michaelis-Menten term cannot be written down from the publication. The",
    "elimination is therefore encoded as first-order with CL/F equal to that",
    "reported intrinsic clearance -- the exact C << Km limit of the",
    "Michaelis-Menten model, in which intrinsic clearance IS the clearance.",
    "No value is invented, but the saturation the authors describe is not",
    "reproduced; see the vignette Assumptions and deviations. Vc and CL are",
    "apparent (/F) formulation-dependent quantities distorted by flip-flop",
    "kinetics, and the authors state that the study 1 and study 2 parameters",
    "cannot be compared directly.",
    sep = " ")
  reference <- paste(
    "Cortez JM Jr, Quintero R, Moss JA, Beliveau M, Smith TJ, Baum MM.",
    "Pharmacokinetics of injectable, long-acting nevirapine for HIV",
    "prophylaxis in breastfeeding infants.",
    "Antimicrob Agents Chemother. 2015;59(1):59-66.",
    "doi:10.1128/AAC.03906-14.",
    sep = " ")
  vignette <- "Cortez_2015_nevirapine"
  units <- list(
    time          = "day",
    dosing        = "mg",
    concentration = "ug/mL"
  )

  compartmentData <- list(
    depot   = list(analyte = "nevirapine", units = "mg", specimen = "administration site", verified = TRUE),
    central = list(analyte = "nevirapine", units = "mg", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Allometric size covariate, reference 0.21 kg (the mean weight of",
        "the n = 10 study-2 rats, Cortez 2015 Methods 'Animals' and Table 2",
        "column header). Cortez 2015 Methods 'Structural PK model buildup'",
        "gives the covariate model theta_i = theta_typical *",
        "(Cov_i / Cov_reference)^theta_eff with nominal exponents -0.25 on",
        "absorption rate, 1 on volumes and 0.75 on clearances; the exponents",
        "were imposed, not estimated, so they are fixed() here."
      ),
      source_name        = "WT"
    )
  )

  population <- list(
    species        = "rat (Sprague-Dawley)",
    n_subjects     = 10L,
    n_studies      = 1L,
    age_range      = "Adult (age not reported)",
    weight_range   = "Mean 0.21 kg (individual weights not reported)",
    weight_median  = "0.21 kg (reported as the mean)",
    disease_state  = "Healthy",
    dose_range     = paste(
      "Single subcutaneous injection of NVP LA-2 at a weight-adjusted",
      "nominal dose of 36 mg/kg (dose concentration 12 mg/mL in 1% w/w",
      "carboxymethylcellulose) through an 18-gauge needle. Cortez 2015",
      "Methods 'Population PK data set construction' assumes 80% of the",
      "intended dose was received, i.e. 6.7 mg (29 mg/kg)."
    ),
    regions        = "United States (Charles River Laboratories, Shrewsbury, MA)",
    notes          = paste(
      "Study 2 (JNN00030, May 2009). Jugular-vein sampling on days -7, 1",
      "(1, 3, 7 and 12 h post-dose), 2, 3, 7, 14, 21 and 28. Plasma",
      "nevirapine by LC-MS/MS (API 4000, TurboIonSpray+, m/z 267.1 ->",
      "226.1; reserpine internal standard) with an LLOQ of 1 ng/mL;",
      "concentrations below the LLOQ were set to missing rather than",
      "imputed, so the plotted late medians in Figure 3B are taken over",
      "only the animals still quantifiable and are biased upward.",
      "Subcutaneous bioavailability was assumed to be 100% (Cortez 2015",
      "Discussion, citing reference 38). The authors note that parameter",
      "precision and residual error were WORSE for this formulation than",
      "for study 1. See Cortez 2015 Table 2 (study 2 column, whose header",
      "misprints the formulation as 'NVP LA-1') and Figure 3B."
    )
  )

  ini({
    # ------------------------------------------------------------------
    # Structural parameters -- Cortez 2015 Table 2, second estimate
    # column ("Study 2 ... n = 10, mean wt 0.21 kg"; the header misprints
    # the formulation name as NVP LA-1, but Methods, Results and the
    # Discussion all identify study 2 with NVP LA-2). Values are apparent
    # (/F) quantities for a 0.21 kg rat.
    # ------------------------------------------------------------------
    lka <- log(15.4) ; label("Absorption rate constant ka (1/day) for a 0.21 kg rat")   # Cortez 2015 Table 2 study 2: Ka = 15.4 /day (RSE 41.8%); consistent with the tabulated t1/2,a = 0.0450 day
    lvc <- log(74.0) ; label("Apparent central volume Vc/F (L) for a 0.21 kg rat")      # Cortez 2015 Table 2 study 2: Vc/F = 74.0 L (RSE 16.2%)

    # ------------------------------------------------------------------
    # Elimination. Cortez 2015 Discussion: "a one-compartment model with
    # linear absorption and saturable elimination adequately described the
    # study 2 formulation." Table 2 reports one number for it,
    #   CL/F = 48.2 L/day (RSE 47.2%),
    # with footnote d: "Represents intrinsic clearance (Vm/Km), where Vm
    # is the maximal metabolic rate and Km is the Michaelis-Menten
    # constant."
    #
    # Vm and Km are NOT reported separately anywhere in the paper, and the
    # Michaelis-Menten term Vm*C/(Km+C) needs both. What the paper does
    # report -- the ratio Vm/Km -- is exactly the limiting clearance of
    # that term when C << Km, so encoding a first-order CL/F of
    # 48.2 L/day uses the published number without inventing anything.
    # The cost is that the saturation itself is not represented; the
    # vignette Assumptions and deviations quantifies the consequence
    # (the linear reading tracks Figure 3B only through roughly day 3).
    # ------------------------------------------------------------------
    lcl <- log(48.2) ; label("Apparent clearance CL/F (L/day) for a 0.21 kg rat, from the reported intrinsic clearance Vm/Km")  # Cortez 2015 Table 2 study 2: CL/F = 48.2 L/day (RSE 47.2%), footnote d: intrinsic clearance Vm/Km

    lfdepot <- fixed(log(1)) ; label("Subcutaneous bioavailability F (fraction)")  # Cortez 2015 Discussion 'Pharmacokinetics of sustained-release NVP in rats': "The subcutaneous bioavailability, F, was assumed to be 100% for both formulations based on previous studies in rats (38)."

    # ------------------------------------------------------------------
    # Allometric exponents -- Cortez 2015 Methods 'Structural PK model
    # buildup'; nominal (imposed) values with no reported uncertainty.
    # ------------------------------------------------------------------
    e_wt_ka <- fixed(-0.25) ; label("Allometric exponent on ka (unitless)")   # Cortez 2015 Methods 'Structural PK model buildup'
    e_wt_vc <- fixed(1)     ; label("Allometric exponent on Vc (unitless)")   # Cortez 2015 Methods 'Structural PK model buildup'
    e_wt_cl <- fixed(0.75)  ; label("Allometric exponent on CL (unitless)")   # Cortez 2015 Methods 'Structural PK model buildup'

    # ------------------------------------------------------------------
    # No IIV -- Cortez 2015 Methods 'Population PK model development':
    # "Due to the low number of rats (n = 4 or n = 10) available for each
    # formulation, no between-subject variability (BSV) was added in the
    # PK model".
    #
    # Residual error -- Cortez 2015 equation 2 allows a combined
    # proportional + additive model, but Table 2 reports only a
    # proportional term for study 2 (+/-79.0%, RSE 9.4%), with no
    # additive row; per equation 2 that is the eps2 = 0 simplification
    # ("This model can be simplified to an additive only model (eps1 = 0)
    # or proportional only model (eps2 = 0)").
    # ------------------------------------------------------------------
    propSd <- 0.790 ; label("Proportional residual error SD (fraction)")   # Cortez 2015 Table 2 study 2 'Error' row: +/-79.0 (9.4) %
  })

  model({
    # ------------------------------------------------------------------
    # 1. Individual parameters, allometrically scaled from the 0.21 kg
    #    study-2 reference weight.
    # ------------------------------------------------------------------
    ka <- exp(lka) * (WT / 0.21)^e_wt_ka
    vc <- exp(lvc) * (WT / 0.21)^e_wt_vc
    cl <- exp(lcl) * (WT / 0.21)^e_wt_cl

    # ------------------------------------------------------------------
    # 2. Micro-constant.
    # ------------------------------------------------------------------
    kel <- cl / vc

    # ------------------------------------------------------------------
    # 3. ODE system -- depot (subcutaneous injection site) feeding a
    #    one-compartment disposition model.
    # ------------------------------------------------------------------
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # ------------------------------------------------------------------
    # 4. Bioavailability (assumed 100%).
    # ------------------------------------------------------------------
    f(depot) <- exp(lfdepot)

    # ------------------------------------------------------------------
    # 5. Observation and error.
    # ------------------------------------------------------------------
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
