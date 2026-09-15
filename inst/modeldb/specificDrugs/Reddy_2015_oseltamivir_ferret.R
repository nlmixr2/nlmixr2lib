Reddy_2015_oseltamivir_ferret <- function() {
  description <- paste(
    "Preclinical (ferret).",
    "Population PK model for oseltamivir carboxylate (OC), the active",
    "metabolite of oseltamivir, after oral dosing of oseltamivir phosphate or",
    "oseltamivir free base: an absorptive compartment feeding two transit",
    "compartments that model the delayed appearance of OC, followed by",
    "two-compartment OC disposition with first-order elimination. Clearance",
    "and volume terms are apparent values conditioned on oral bioavailability",
    "(F) and the fraction of parent converted to metabolite (Fm), and are",
    "normalised per kg of ferret body weight (dose is entered as ug of",
    "oseltamivir free base per kg). Parameter values are the pooled influenza",
    "A + B (ketamine-anaesthetised) Monte Carlo parameter set the authors used",
    "for all reported simulations."
  )
  reference <- paste(
    "Reddy MB, Yang K-H, Rao G, Rayner CR, Nie J, Pamulapati C, Marathe BM,",
    "Forrest A, Govorkova EA (2015). Oseltamivir Population Pharmacokinetics",
    "in the Ferret: Model Application for Pharmacokinetic/Pharmacodynamic",
    "Study Design. PLoS ONE 10(10):e0138069.",
    "doi:10.1371/journal.pone.0138069."
  )
  vignette <- "Reddy_2015_oseltamivir_ferret"
  units <- list(time = "h", dosing = "ug", concentration = "ng/mL")

  # No covariate effect is carried in model(): the paper screened anaesthesia
  # and influenza-inoculation strain by Kruskal-Wallis tests on the individual
  # post-hoc parameter estimates (Tables 4 and 5) rather than by building
  # covariate coefficients into the structural model. The screened covariates
  # are documented below so the provenance of that screen is not lost.
  covariatesDataExcluded <- list(
    CONMED_KETAMINE = list(
      description = paste(
        "1 = ferret received intramuscular ketamine (10 mg/kg in Study 1,",
        "25 mg/kg in Study 2) before blood sampling and dosing; 0 = no",
        "anaesthesia (Study 3). Screened in Reddy 2015 Table 4 by",
        "Kruskal-Wallis test on the individual post-hoc parameter estimates:",
        "Kt 1.27 vs 4.18 1/h (p < 0.001), Ka 0.463 vs 0.335 1/h (p = 0.01),",
        "CLt 1.52 vs 0.919 L/h (p = 0.001), Vc 0.157 vs 1.00 L (p = 0.002),",
        "Vp 5.59 vs 2.08 L (p < 0.001), CLd 0.585 vs 0.878 L/h (p = 0.174).",
        "The differences are statistically significant but the authors did",
        "not estimate a ketamine covariate coefficient; they reported the two",
        "stratified parameter summaries side by side and carried only the",
        "ketamine (Studies 1 + 2) set forward into the simulations."
      ),
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no anaesthesia)",
      notes              = paste(
        "Screened and reported as significant, but NOT retained as a",
        "structural covariate effect in this model. The no-anaesthesia",
        "parameter set (Reddy 2015 Table 4, Study 3 column, n = 8 ferrets)",
        "is reported without a covariance matrix and was never used for",
        "simulation; it is reproduced in the validation vignette rather than",
        "as a second model file. A third anaesthetic, Saffan (Study 4), was",
        "excluded from model development entirely."
      )
    )
  )

  compartmentData <- list(
    depot = list(
      analyte = "oseltamivir", units = "ug",
      specimen = "administration site", verified = TRUE
    ),
    transit1 = list(
      analyte = "oseltamivir carboxylate", units = "ug",
      specimen = "administration site", verified = TRUE
    ),
    transit2 = list(
      analyte = "oseltamivir carboxylate", units = "ug",
      specimen = "administration site", verified = TRUE
    ),
    central = list(
      analyte = "oseltamivir carboxylate", units = "ug",
      specimen = "plasma", verified = TRUE
    ),
    peripheral1 = list(
      analyte = "oseltamivir carboxylate", units = "ug",
      specimen = "plasma", verified = TRUE
    )
  )

  population <- list(
    species        = "ferret (Mustela putorius furo)",
    n_subjects     = 65,
    n_studies      = 3,
    n_observations = 430,
    age_range      = "3-5 months (young adult)",
    weight_range   = "0.6-0.92 kg",
    sex_female_pct = 12.3,
    disease_state  = paste(
      "uninfected (n = 17) or inoculated with influenza A/Shenzheng/406H/2006",
      "(H5N1, n = 18), influenza A/Hong Kong/433581/2009 (H3N2, n = 12) or",
      "influenza B/Yamagata/16/1988 (n = 18); all inoculations produced only",
      "mild, essentially subclinical illness"
    ),
    dose_range     = paste(
      "0.76-25 mg/kg oseltamivir free base (equivalently 1.0-32.9 mg/kg",
      "oseltamivir phosphate) orally, as single doses or every 12 h for 5 days"
    ),
    regions        = "Beijing (China), Memphis (USA), London (UK)",
    notes          = paste(
      "Reddy 2015 Methods 'PK Studies in a Ferret Model' and 'OC PK Model'.",
      "Studies 1 and 2 used male ferrets, Study 3 female ferrets, so the",
      "female fraction is 8/65. The 65 animals contributed 430 OC",
      "concentrations above the 10 ng/mL limit of quantitation. Study 4",
      "(Saffan anaesthesia, Tmax 7 h for both OP and OC) was excluded from",
      "model development because its NCA profile differed from Studies 1-3."
    )
  )

  ini({
    # =======================================================================
    # Structural parameters. Reddy 2015 Table 6 'Diagonals / Mean' column --
    # the pooled influenza A + B (Studies 1 and 2, ketamine anaesthesia)
    # iterative-two-stage parameter set the authors carried into every
    # reported Monte Carlo simulation. The same means appear in Table 4,
    # 'Studies 1 and 2, ketamine' column.
    #
    # Footnote a of Tables 4 and 6: the CL and V parameters are conditioned on
    # oral bioavailability and fraction of parent changed to metabolite, i.e.
    # they are apparent (CL/(F*Fm), V/(F*Fm)) values.
    #
    # Unit system: the modelling dataset (S1 Table) carries dose as ug of
    # oseltamivir free base per kg and OC concentration as ng/mL, with no
    # body-weight column, so the amounts are per-kg amounts in ug and the
    # 'L' / 'L/h' the paper prints on Vc, Vp, CLd and CLt are per-kg values.
    # Cc = central/vc then returns ug/L = ng/mL, matching the dataset.
    # =======================================================================
    lktr <- log(1.27)  ; label("Transit transfer rate constant Kt (1/h)")                                   # Table 6 Kt mean = 1.27 1/h
    lka  <- log(0.463) ; label("Absorption rate constant Ka, second transit compartment to central (1/h)")  # Table 6 Ka mean = 0.463 1/h
    lq   <- log(0.585) ; label("Apparent distribution clearance CLd/(F*Fm) (L/h per kg)")                   # Table 6 CLd mean = 0.585 L/h
    lcl  <- log(1.52)  ; label("Apparent elimination clearance CLt/(F*Fm) (L/h per kg)")                    # Table 6 CLt mean = 1.52 L/h
    lvc  <- log(0.157) ; label("Apparent central volume Vc/(F*Fm) (L per kg)")                              # Table 6 Vc mean = 0.157 L
    lvp  <- log(5.59)  ; label("Apparent peripheral volume Vp/(F*Fm) (L per kg)")                           # Table 6 Vp mean = 5.59 L

    # =======================================================================
    # Between-animal variability. Reddy 2015 Table 6 reports variances and
    # covariances on the NATURAL parameter scale while Methods 'OC PK Model'
    # states 'We assumed parameters to be log-normally distributed'. The
    # log-scale entries below are the moment-matched transforms
    #     omega^2(i)   = log(1 + var(i) / mean(i)^2)
    #     cov_log(i,j) = log(1 + cov(i,j) / (mean(i) * mean(j)))
    # so that exp(l<param>) is the median of the log-normal and the natural-
    # scale CV% reproduces Table 6 exactly. Check of the natural-scale CV%
    # implied by Table 6: Kt 49.9, Ka 66.1, CLd 50.0, CLt 55.2, Vc 49.9,
    # Vp 63.5 -- reproducing the 50% values the authors state they imposed on
    # Kt, CLd and Vc, and the Table 4 ketamine CV% for Ka (66.1), CLt (55.1)
    # and Vp (63.6) which were left as estimated.
    # =======================================================================

    # Kt, CLd and Vc: Methods 'Simulations' -- 'CV% for Vc, Cld, and Kt were
    # empirically reduced to 50%, and all related covariance terms were fixed
    # to zero'. Imposed by the authors rather than estimated, hence fixed();
    # the estimated CV% for these three were 79.8, 96.6 and 165 (Table 4).
    # Being outside the correlated block encodes the zeroed covariances.
    etalktr ~ fixed(0.2225358)   # Table 6 Kt variance 0.402, mean 1.27 -> log(1 + 0.402/1.27^2)
    etalq   ~ fixed(0.2232458)   # Table 6 CLd variance 0.0856, mean 0.585 -> log(1 + 0.0856/0.585^2)
    etalvc  ~ fixed(0.2227459)   # Table 6 Vc variance 0.00615, mean 0.157 -> log(1 + 0.00615/0.157^2)

    # Ka, CLt and Vp: estimated variances with the three off-diagonal
    # covariances Table 6 retained (all other covariances were fixed to zero
    # because their correlation was below 0.15). Natural-scale correlations
    # are CLt-Ka 0.405, Vp-Ka -0.585, Vp-CLt -0.290; the log-scale block
    # below is positive definite (eigenvalues 0.702, 0.199, 0.066).
    etalka + etalcl + etalvp ~ c(
       0.3626250,
       0.1378276,  0.2659801,
      -0.2820095, -0.1071039,  0.3387728
    ) # Table 6: variances Ka 0.0937, CLt 0.704, Vp 12.6; covariances 'CLt, Ka' 0.104, 'Vp, Ka' -0.636, 'Vp, CLt' -0.863

    # =======================================================================
    # Residual error. Reddy 2015 Results 'Population PK Model', final
    # paragraph: 'The intercept for the additive error model was 5 ng/mL,
    # roughly equal to half of the LOQ. The slope was 0.15, roughly equal to
    # the coefficient of variation (CV%) of the assay performance.' Methods
    # 'OC PK Model': 'The residual error model included both additive and
    # proportional error.' The S1 Table records below-limit values as
    # 'BLQ<10.0', confirming the 10 ng/mL LOQ that 5 ng/mL is half of.
    # =======================================================================
    addSd  <- 5    ; label("Additive residual error (ng/mL)")            # Results, 'The intercept for the additive error model was 5 ng/mL'
    propSd <- 0.15 ; label("Proportional residual error (fraction)")     # Results, 'The slope was 0.15'
  })

  model({
    # 1. Individual parameters (log-normal, per Methods 'OC PK Model').
    ktr <- exp(lktr + etalktr)
    ka  <- exp(lka  + etalka)
    q   <- exp(lq   + etalq)
    cl  <- exp(lcl  + etalcl)
    vc  <- exp(lvc  + etalvc)
    vp  <- exp(lvp  + etalvp)

    # 2. ODE system, Reddy 2015 eqs (1)-(5) and Fig 1. Canonical names map to
    #    the paper's states as depot = X1 (absorptive compartment holding the
    #    administered oseltamivir), transit1 = X2 and transit2 = X3 (the two
    #    transit compartments holding OC), central = X4 (Vc) and
    #    peripheral1 = X5 (Vp). Kt moves depot -> transit1 -> transit2 and Ka
    #    moves transit2 -> central, exactly as drawn in Fig 1.
    #      eq (1) dX1/dt = -Kt*X1
    #      eq (2) dX2/dt =  Kt*X1 - Kt*X2
    #      eq (3) dX3/dt =  Kt*X2 - Ka*X3
    #      eq (4) dX4/dt =  Ka*X3 - CLd*X4/Vc + CLd*X5/Vp - CLt*X4/Vc
    #      eq (5) dX5/dt =  CLd*X4/Vc - CLd*X5/Vp
    d/dt(depot)       <- -ktr * depot
    d/dt(transit1)    <-  ktr * depot    - ktr * transit1
    d/dt(transit2)    <-  ktr * transit1 - ka * transit2
    d/dt(central)     <-  ka * transit2 - (q / vc) * central + (q / vp) * peripheral1 - (cl / vc) * central
    d/dt(peripheral1) <-  (q / vc) * central - (q / vp) * peripheral1

    # 3. Observation. Methods 'OC PK Model': 'The concentration of OC in the
    #    plasma is calculated as X4/Vc.' No separate bioavailability term is
    #    applied because F (and Fm) are already folded into the apparent CL
    #    and V parameters (Table 6 footnote a).
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
