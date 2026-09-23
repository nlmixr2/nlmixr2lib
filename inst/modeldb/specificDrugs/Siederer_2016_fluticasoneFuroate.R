Siederer_2016_fluticasoneFuroate <- function() {
  description <- "Two-compartment population PK model with first-order absorption for inhaled fluticasone furoate in subjects with COPD and healthy volunteers, with a composite race effect on apparent inhaled clearance"
  reference <- "Siederer S, Allen A, Yang S. Population Pharmacokinetics of Inhaled Fluticasone Furoate and Vilanterol in Subjects with Chronic Obstructive Pulmonary Disease. Eur J Drug Metab Pharmacokinet. 2016;41(6):743-758. doi:10.1007/s13318-015-0303-4"
  vignette <- "Siederer_2016_fluticasoneFuroate_vilanterol"
  units <- list(time = "h", dosing = "ug", concentration = "ng/mL")
  # Unit note: doses are entered in ug and volumes are in L, so `Cc` is in
  # ug/L == ng/mL. Siederer 2016 reports fluticasone furoate concentrations and
  # exposures in pg/mL and pg*h/mL (assay LLQ 10 pg/mL, Sect. 2.2); multiply
  # `Cc` by 1000 to compare against the published values. No scale factor is
  # applied inside the model so that dose / volume / clearance stay mutually
  # consistent.

  covariateData <- list(
    RACE_ASIAN_EAST_SE = list(
      description = "Composite East Asian / Japanese / South East Asian heritage indicator (Siederer 2016 RACE1 = 2)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = White/Caucasian (RACE1 = 1), the model reference category",
      notes = "1 if the subject's race is East Asian, Japanese, or South East Asian heritage; 0 otherwise. 14% of the fluticasone furoate dataset (Sect. 2.3.1). Mutually exclusive with RACE_BLACK and RACE_ASIAN_CENTRAL_ARABIC_AMIND_OTH; all three = 0 selects the White/Caucasian reference. Explicitly EXCLUDES Central / South Asian heritage, which Siederer 2016 pools into RACE1 = 4 instead.",
      source_name = "RACE1 = 2"
    ),
    RACE_BLACK = list(
      description = "Black / African American race indicator (Siederer 2016 RACE1 = 3)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = White/Caucasian (RACE1 = 1), the model reference category",
      notes = "3% of the fluticasone furoate dataset (Sect. 2.3.1). The coefficient is very imprecisely estimated (%RSE 199, 95% CI spans 1) and Sect. 3.1.2 warns it should be interpreted with caution.",
      source_name = "RACE1 = 3"
    ),
    RACE_ASIAN_CENTRAL_ARABIC_AMIND_OTH = list(
      description = "Composite Central/South Asian, White-Arabic/North African, American Indian/Alaska Native and 'other' heritage indicator (Siederer 2016 RACE1 = 4)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = White/Caucasian (RACE1 = 1), the model reference category",
      notes = "1 if the subject's race is Central/South Asian, White-Arabic/North African, American Indian/Native Alaskan, or 'other'; 0 otherwise. 2% of the fluticasone furoate dataset (Sect. 2.3.1). Mutually exclusive with RACE_ASIAN_EAST_SE and RACE_BLACK.",
      source_name = "RACE1 = 4"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "fluticasone furoate", units = "ug", specimen = "administration site", verified = TRUE),
    central = list(analyte = "fluticasone furoate", units = "ug", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "fluticasone furoate", units = "ug", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 1307,
    n_studies = 4,
    age_range = "18-85 years",
    age_median = "61.0 years",
    weight_range = "35.0-174.6 kg",
    weight_mean = "75.7 kg",
    sex_female_pct = 67,
    race_ethnicity = c(
      `White/Caucasian/European` = 82,
      `African American/African` = 3,
      `Asian - East Asian` = 5,
      `Asian - Japanese` = 4,
      `Asian - South East Asian` = 5,
      `Asian - Central/South Asian` = 1,
      `American Indian/Native Alaskan` = 1,
      `White - Arabic/North African` = 1,
      Other = 1
    ),
    disease_state = "chronic obstructive pulmonary disease (94% of subjects) pooled with healthy volunteers (6%)",
    dose_range = "fluticasone furoate 50, 100, 200, 400 or 800 ug once daily by oral inhalation, alone or as the fluticasone furoate/vilanterol combination",
    regions = "global (Argentina, Chile, Czech Republic, Estonia, Germany, Japan, Korea, Mexico, Norway, Philippines, Poland, Russia, Sweden, USA and others)",
    notes = "Three Phase III studies (HZC112206, HZC112207, HZC110946) in COPD plus one Phase I study (HZA102936) in healthy volunteers; the Phase II study HZC111348 was excluded from the fluticasone furoate analysis because it provided only 0-4 h post-dose data (Sect. 2.1). Demographics from Online Resource Table S2; mean BMI 26.3 kg/m2, mean height 169 cm, mean percent-predicted FEV1 48.1% in the COPD population. 11,789 observations, 39% below the 10 pg/mL assay LLQ, handled by the NONMEM M3 likelihood method."
  )

  ini({
    # -----------------------------------------------------------------------
    # Structural parameters. Siederer 2016 Table 1 reports each THETA twice:
    # as the estimated log-scale value ('Ln estimate') and as its exponential
    # ('Estimate'). The untransformed column is used here because it carries
    # more significant figures; log() of it reproduces the Ln column
    # (e.g. log(230) = 5.438 vs the printed 5.44).
    # -----------------------------------------------------------------------
    lka <- log(0.0523); label("Apparent first-order absorption rate constant after oral inhalation (1/h)") # Table 1 'ka (h-1)' = 0.0523 (95% CI 0.0493, 0.0556; RSE 1.06%); Ln estimate -2.95
    lcl <- log(230); label("Apparent inhaled clearance CL/F for a White/Caucasian subject (L/h)") # Table 1 'CL/F (L/h)' = 230 (95% CI 219, 242; RSE 0.47%); Ln estimate 5.44; Sect. 3.1.2 'The typical value of CL/F was 230 L/h for a white Caucasian subject with COPD'
    lvc <- fixed(log(1.36)); label("Apparent central volume of distribution V2/F (L)") # Table 1 'V2/F (L)' = 1.36 FIXED; Ln estimate 0.31 FIXED; Sect. 3.1.2 'the volume of the central compartment (V2/F) was fixed to a value appropriate for central V2 (1.36 L) following evaluation of a range of values (unpublished data, GSK, UK, 2012)'
    lq <- log(268); label("Apparent intercompartmental clearance Q/F (L/h)") # Table 1 'Q/F (L/h)' = 268 (95% CI 221, 324; RSE 1.73%); Ln estimate 5.59
    lvp <- log(111); label("Apparent peripheral volume of distribution V3/F (L)") # Table 1 'V3/F (L)' = 111 (95% CI 90.9, 136; RSE 2.21%); Ln estimate 4.71

    # -----------------------------------------------------------------------
    # Race effect on CL/F. Siederer 2016 Eq. (a) in Sect. 3.1.2 is
    #   Ln(theta) = theta_1 + COV
    # i.e. the category coefficient is added directly to the log-scale
    # population estimate with no multiplying theta, and COV = 0 for the
    # RACE1 = 1 (White/Caucasian) reference. The coefficients below are the
    # 'Ln estimate' column of Table 1 and reproduce the untransformed
    # multipliers and the typical CL/F values quoted in Sect. 3.1.2:
    #   exp(5.44 - 0.211) = 187 L/h  ('186 L/h' for RACE1 = 2)
    #   exp(5.44 + 0.0602) = 245 L/h ('244 L/h' for RACE1 = 3)
    #   exp(5.44 - 0.265)  = 177 L/h ('176 L/h' for RACE1 = 4)
    # -----------------------------------------------------------------------
    e_race_asian_east_se_cl <- -0.211; label("Log-scale effect of East Asian/Japanese/South East Asian heritage on CL/F (unitless)") # Table 1 'RACE1 = 2 on CL/F' Ln estimate -0.211 (95% CI -0.329, -0.0930), multiplier 0.810 (95% CI 0.720, 0.911; RSE 28.5%)
    e_race_black_cl <- 0.0602; label("Log-scale effect of Black/African American race on CL/F (unitless)") # Table 1 'RACE1 = 3 on CL/F' Ln estimate 0.0602 (95% CI -0.175, 0.295), multiplier 1.062 (95% CI 0.839, 1.343; RSE 199.0%)
    e_race_asian_central_arabic_amind_oth_cl <- -0.265; label("Log-scale effect of Central Asian/White-Arabic/American Indian/other heritage on CL/F (unitless)") # Table 1 'RACE1 = 4 on CL/F' Ln estimate -0.265 (95% CI -0.528, -0.002), multiplier 0.767 (95% CI 0.590, 0.998; RSE 50.6%)

    # -----------------------------------------------------------------------
    # Inter-individual variability. Sect. 3.1.2 states only that
    # 'Inter-individual variances (exponential model) were estimated with
    # reasonable precision (%RSE <= 36 %)' -- Table 1 lists no OMEGA row, and
    # neither the paper nor the Online Resource reports which parameters
    # carried an eta or how large any variance was. Only IIV on CL/F is
    # positively evidenced, because Sect. 2.5 derives 'individual post hoc
    # estimates of CL/F' per subject. It is therefore declared here at
    # fixed(0) so the structure is preserved without inventing a variance;
    # see the vignette 'Assumptions and deviations' section.
    # -----------------------------------------------------------------------
    etalcl ~ fixed(0) # Sect. 3.1.2 'Inter-individual variances (exponential model) were estimated with reasonable precision' - magnitude not reported anywhere in the paper or Online Resource

    # -----------------------------------------------------------------------
    # Residual error. Sect. 3.1.2: 'An additive error model described the
    # residual variability.' Fig. 1 plots observed FF concentration (pg/mL)
    # against population/individual prediction (pg/mL) on LINEAR axes, so the
    # additive error is additive on the untransformed concentration scale
    # (contrast the vilanterol model, whose Fig. 3 axes are log-Value /
    # log-Prediction). The SIGMA magnitude is not reported, so it is declared
    # at fixed(0).
    # -----------------------------------------------------------------------
    addSd <- fixed(0); label("Additive residual error (ng/mL)") # Sect. 3.1.2 'An additive error model described the residual variability' - magnitude not reported
  })

  model({
    # Apparent inhaled clearance; race enters log-additively per Sect. 3.1.2
    # Eq. (a). All three indicators are 0 for the White/Caucasian reference.
    cl <- exp(
      lcl + etalcl +
        e_race_asian_east_se_cl * RACE_ASIAN_EAST_SE +
        e_race_black_cl * RACE_BLACK +
        e_race_asian_central_arabic_amind_oth_cl * RACE_ASIAN_CENTRAL_ARABIC_AMIND_OTH
    )
    ka <- exp(lka)
    vc <- exp(lvc)
    q <- exp(lq)
    vp <- exp(lvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    Cc <- central / vc
    Cc ~ add(addSd)
  })
}
