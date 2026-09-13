Darwish_2026_remlifanserin <- function() {
  description <- "Population PK model for oral remlifanserin (ACP-204), a selective 5-HT2A receptor inverse agonist (Darwish 2026): one-compartment with first-order absorption, an absorption lag time and linear elimination, pooled across seven Phase 1 studies in healthy young and older adults. Fed status shifts the absorption lag time and female sex lowers the apparent central volume."
  reference <- "Darwish M, Lin N, Dirks B, Jaworowicz D, Maxwell K, Pathak S. Population pharmacokinetics of remlifanserin (ACP-204), a serotonin 2A receptor inverse agonist. Alzheimer's & Dementia: Translational Research & Clinical Interventions. 2026;12(1):e70254. doi:10.1002/trc2.70254"
  vignette <- "Darwish_2026_remlifanserin"

  # The model is parameterised in DAYS throughout, exactly as Darwish 2026
  # Table 1 and the Supplemental Methods NONMEM control stream report it
  # (CL/F 849 L/day, Ka 16.6 1/day, ALAG1 0.0300 day). The control stream's
  # scaling line is `S2=V ;Dose amount = ug; Volume = L; Concentration =
  # ng/mL`, so DOSE AMOUNTS ARE IN MICROGRAMS: the 60 mg target dose is
  # entered as 60000. Darwish 2026 Section 3.4 restates the same estimates in
  # hours (CL/F 35.4 L/h, Ka 0.692 1/h, ALAG1 0.720 h fasted / 0.929 h fed);
  # those are the identical values, not a second parameterisation.
  units <- list(time = "day", dosing = "ug", concentration = "ng/mL")

  covariateData <- list(
    FED = list(
      description        = "Fed-versus-fasted state at the time of dosing",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (fasted)",
      notes              = "Time-varying per dose record, not a subject-level attribute: Darwish 2026 Table 2 footnote states 'ALAG1 is nonstationary per participant. Typical ALAG1 value varies by fed status (time-varying). Some participants in the food effect study had > 1 distinct ALAG1 value.' Proportional shift on the absorption lag time only; the lag rises 29.0% from 0.0300 day (0.720 h) fasted to 0.0387 day (0.929 h) fed. Darwish 2026 Figure 4 shows the resulting steady-state Cmax and AUC0-24 geometric mean ratios fall well inside the 0.8-1.25 clinical-relevance bounds, so the effect is statistically significant (P < 0.001) but not clinically relevant. Records were 71% fasted and 29% fed overall (Supplemental Table 3); the crossover food-effect study 010 contributed both states within subject.",
      source_name        = "FED"
    ),
    SEXF = list(
      description        = "Biological sex indicator",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Additive (not proportional) shift on the apparent central volume, in litres: V/F is 930 L for a typical male and 930 - 131 = 799 L for a typical female (Darwish 2026 Section 3.4). The source column is already named SEXF with the canonical 1 = female orientation (Supplemental Methods $INPUT), so no value transformation is needed. Darwish 2026 Section 3.2 notes that female sex was correlated with lower body weight but that sex was the stronger predictor of V/F in the stepwise analysis, so weight is not in the final model. The analysis population was 30.6% female (64/209, Supplemental Table 2).",
      source_name        = "SEXF"
    )
  )

  # Covariates that Darwish 2026 Section 2.4 screened but did not retain in the
  # final model. Several are still computed in the Supplemental Methods control
  # stream's $PK block (RACEB, GFRCKDEP2, RFCAT2) as leftovers of the stepwise
  # covariate search, but none enters a parameter equation.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened on all PK parameters; not retained. Analysis population 19-75 years, median 37 (Supplemental Table 2). Only 18 participants were 65-75 years old, which Darwish 2026 Section 4.3 flags as a limitation."
    ),
    WT = list(
      description = "Body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened; not retained. Range 47-108 kg, median 77.6 (Supplemental Table 2). Darwish 2026 Section 3.2 reports that weight was correlated with sex and that sex was the stronger predictor of V/F."
    ),
    BMI = list(
      description = "Body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Screened; not retained. Range 18.5-31.9 kg/m^2, median 26.3 (Supplemental Table 2)."
    ),
    BSA = list(
      description = "Body surface area",
      units       = "m^2",
      type        = "continuous",
      notes       = "Carried in the analysis dataset ($INPUT column BSA in the Supplemental Methods control stream) but not retained in the final model."
    ),
    RACE_BLACK = list(
      description = "Black / African American race indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened; not retained. The Supplemental Methods $PK block still derives RACEB = 1 when RACEN == 2, regrouping Asian and Other with White, but never uses it. Analysis population 59.3% White, 33.0% Black, 2.87% Asian, 4.78% Other (Supplemental Table 2)."
    ),
    CRCL = list(
      description = "BSA-normalized estimated glomerular filtration rate (CKD-EPI)",
      units       = "mL/min/1.73 m^2",
      type        = "continuous",
      notes       = "Screened; not retained. Source column GFRCKDEP. Range 58.4-138 mL/min/1.73 m^2, median 103 (Supplemental Table 2). The control stream imputes 9 missing values to the median of 103 but never uses the covariate."
    ),
    RENALIMP_MILD = list(
      description = "Mild renal impairment indicator (baseline renal function category)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened; not retained. Source column RFCAT; the control stream derives RFCAT2 = 1 when RFCAT == 2 and imputes missing to normal, but never uses it. Nearly all participants had normal renal function (Darwish 2026 Section 4.3)."
    ),
    HEPIMP = list(
      description = "Hepatic impairment indicator (NCI ODWG baseline liver function category)",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened; not retained. Source column NCILIV. Nearly all participants had normal liver function (Darwish 2026 Section 4.3)."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened as part of the liver-function assessment; not retained. Range 6-74 U/L, median 20; 10 values missing (Supplemental Table 2)."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened as part of the liver-function assessment; not retained. Range 12-36 U/L, median 19; 10 values missing (Supplemental Table 2)."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units       = "mg/dL",
      type        = "continuous",
      notes       = "Screened as part of the liver-function assessment; not retained. Range 0.2-1.6 mg/dL, median 0.5 (Supplemental Table 2). Reported in US conventional units by Darwish 2026; the register's canonical unit for TBILI is umol/L (multiply by 17.104)."
    )
  )

  compartmentData <- list(
    depot   = list(analyte = "remlifanserin", units = "ug", specimen = "administration site", verified = TRUE),
    central = list(analyte = "remlifanserin", units = "ug", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 209,
    n_studies      = 7,
    n_observations = 3935,
    age_range      = "19-75 years",
    age_median     = "37 years",
    weight_range   = "47-108 kg",
    weight_median  = "77.6 kg",
    bmi_range      = "18.5-31.9 kg/m^2",
    bmi_median     = "26.3 kg/m^2",
    sex_female_pct = 30.6,
    race_ethnicity = c(White = 59.3, Black = 33.0, Asian = 2.87, Other = 4.78),
    disease_state  = "healthy volunteers; 18-55 year old adults plus healthy older adults aged 65-75 years (n = 18)",
    dose_range     = "10-180 mg oral; single dose and once-daily multiple dose for 10 days",
    renal_function = "eGFR (CKD-EPI) 58.4-138 mL/min/1.73 m^2, median 103; predominantly normal",
    prandial_state = "2832 of 3987 dose-associated records fasted (71%), 1155 fed (29%)",
    notes          = "Pooled from seven Phase 1 studies of remlifanserin (ACP-204): a single-ascending-dose study with a food arm (001), a multiple-ascending-dose study including an elderly cohort (002), a PET receptor-occupancy study (003), itraconazole and carbamazepine drug-drug-interaction studies (004, 005; only period 1 retained via the control stream's IGNORE=(DELFN.EQ.52)), a crossover food-effect study (010), and a [14C] mass-balance study (011). Study designs are in Supplemental Table 1 and baseline demographics in Supplemental Table 2. Concentrations were quantified by validated LC-MS/MS with an LLOQ of 0.10 ng/mL; the 1.3% of samples below the LLOQ were treated as missing and excluded."
  )

  ini({
    # Structural parameters. Darwish 2026 Table 1, 'Final parameter estimate /
    # Population mean' column. The Supplemental Methods control stream is the
    # INITIAL-estimate version (THETA 845, 965, 16.7, 0.03, 0.289, -222), so
    # every value below is taken from Table 1, not from $THETA.
    lcl <- log(849); label("Apparent oral clearance (L/day)") # Table 1, 'CL/F (L/day)' = 849 (RSE 2.84%; bootstrap 853, 95% CI 806-899)
    lvc <- log(930); label("Apparent central volume of distribution in a male (L)") # Table 1, 'V/F (L)' = 930 (RSE 2.43%; bootstrap 936, 95% CI 892-978)
    lka <- log(16.6); label("First-order absorption rate constant (1/day)") # Table 1, 'Ka (1/day)' = 16.6 (RSE 4.65%; bootstrap 16.7, 95% CI 15.4-18.1)
    ltlag <- log(0.0300); label("Absorption lag time in the fasted state (day)") # Table 1, 'ALAG1 (day)' = 0.0300 (RSE 1.91%; bootstrap 0.0314, 95% CI 0.0306-0.0321)

    # Covariate effects. Both are the paper's own printed coefficients; the
    # equations they enter are quoted verbatim in Darwish 2026 Section 3.2 and
    # reproduced in model() below.
    e_sexf_vc <- -131; label("Additive shift in apparent central volume for females (L)") # Table 1, 'Additive shift in V/F for females' = -131 (RSE 16.8%; bootstrap -131, 95% CI -174 to -90.4)
    e_fed_tlag <- 0.290; label("Proportional shift in absorption lag time when dosed fed (unitless)") # Table 1, 'Proportional shift in ALAG1 for fed status' = 0.290 (RSE 10.4%)

    # Interindividual variability. Darwish 2026 Table 1 prints IIV only as a
    # percent CV, so each variance below is back-transformed with the
    # log-normal relation omega^2 = log(1 + CV^2). Two independent checks pin
    # that scale rather than the omega = CV approximation:
    #   (1) Table 1 footnote gives r = 0.852 and r^2 = 0.726 for the
    #       CL/F-V/F covariance of 0.109. log(1 + CV^2) reproduces
    #       0.109 / sqrt(0.16676 * 0.09807) = 0.8523; omega = CV would give
    #       0.109 / (0.426 * 0.321) = 0.7971, which contradicts the footnote.
    #   (2) The Supplemental Methods $OMEGA initial estimates (0.16 for CL/F,
    #       0.09 for V/F, 0.10 covariance, 0.30 for Ka, 0.01 for ALAG1) sit
    #       beside the back-transformed finals, confirming these are NONMEM
    #       variances on the log scale.
    etalcl + etalvc ~ c(0.1668,
                        0.1090, 0.0981) # Table 1: CL/F 42.6 %CV -> log(1+0.426^2)=0.1668; covariance 'Covariance (IIV in V/F, IIV in CL/F)' = 0.109; V/F 32.1 %CV -> log(1+0.321^2)=0.0981
    etalka ~ 0.3754 # Table 1, 'Ka' 67.5 %CV -> log(1+0.675^2) = 0.3754
    etaltlag ~ 0.01138 # Table 1, 'ALAG1' 10.7 %CV -> log(1+0.107^2) = 0.01138

    # Residual error. NONMEM $SIGMA is a variance and the $ERROR block is
    # Y = IPRED + IPRED*EPS(1), i.e. purely proportional; sqrt(0.0388) = 0.1970
    # reproduces the 19.7 %CV that Table 1 prints alongside it.
    propSd <- 0.197; label("Proportional residual error (fraction)") # Table 1, 'Residual variability' = 0.0388 (variance; RSE 6.38%), reported as 19.7 %CV
  })

  model({
    # Individual parameters. Structure follows the Supplemental Methods $PK
    # block exactly: the sex shift is additive on the TYPICAL volume, inside
    # the exponential eta (TVV = TVVI + COV2; V = TVV*EXP(ETA(2))), and the
    # food shift is proportional on the typical lag time
    # (TVALAG1 = THETA(4)*(1+THETA(5)*FED); ALAG1 = TVALAG1*EXP(ETA(4))).
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl)
    # Darwish 2026 Section 3.2: V/F_i = 930 L - 131 * female_i
    vc <- (exp(lvc) + e_sexf_vc * SEXF) * exp(etalvc)
    # Darwish 2026 Section 3.2: ALAG1_i = 0.0300 d * (1 + 0.290 * fed_i)
    tlag <- exp(ltlag) * (1 + e_fed_tlag * FED) * exp(etaltlag)

    kel <- cl / vc

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    alag(depot) <- tlag

    # NONMEM ADVAN2 TRANS2 with S2 = V. Only oral data were available, so CL
    # and V are apparent values (CL/F, V/F) and the dose enters the depot
    # unscaled; F is not separately identifiable (Darwish 2026 Section 2.3).
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
