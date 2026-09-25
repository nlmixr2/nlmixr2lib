Mehta_2020_fluticasoneFuroate <- function() {
  description <- "Two-compartment population PK model with first-order absorption for inhaled fluticasone furoate in adults with COPD receiving FF/UMEC/VI single-inhaler triple therapy, FF/VI + UMEC or FF/VI, with Japanese-heritage and FF/VI-arm effects on apparent inhaled clearance"
  reference <- "Mehta R, Farrell C, Hayes S, Birk R, Okour M, Lipson DA. Population Pharmacokinetic Analysis of Fluticasone Furoate/Umeclidinium Bromide/Vilanterol in Patients with Chronic Obstructive Pulmonary Disease. Clin Pharmacokinet. 2020;59(1):67-79. doi:10.1007/s40262-019-00794-w"
  vignette <- "Mehta_2020_fluticasoneFuroate_umeclidinium_vilanterol"
  units <- list(time = "h", dosing = "ug", concentration = "ng/mL")
  # Unit note: doses are entered in ug and volumes are in L, so `Cc` is in
  # ug/L == ng/mL. Mehta 2020 reports concentrations and exposures in pg/mL
  # and pg*h/mL (assay LLOQ 10 pg/mL, Sect. 2.2); multiply `Cc` by 1000 to
  # compare against the published values. No scale factor is applied inside
  # the model so that dose / volume / clearance stay mutually consistent
  # (same convention as the Siederer 2016 fluticasone furoate model of the
  # same lineage).

  covariateData <- list(
    RACE_JAPANESE = list(
      description = "Japanese-heritage indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = not of Japanese heritage (White, African American/African, East Asian and other)",
      notes = "1 = Asian-Japanese heritage (13% of the FF dataset, Table 2). East Asian heritage was tested separately in the full model and dropped (Sect. 3.2), so East Asian subjects stay in the reference group.",
      source_name = "Japanese heritage"
    ),
    CONMED_UMECLIDINIUM = list(
      description = "Umeclidinium co-administration indicator (treatment-arm assignment)",
      units = "(binary)",
      type = "binary",
      reference_category = "1 = umeclidinium co-administered (FF/UMEC/VI single inhaler or FF/VI + UMEC two inhalers), the model's reference",
      notes = "0 = FF/VI dual-therapy arm (IMPACT study only), which carries the 1.42-fold CL/F multiplier (Table 3 'FF/VI on CL/F'). The paper's categorical indicator is CAT = 1 for FF/VI (Eq. 2), i.e. CAT = 1 - CONMED_UMECLIDINIUM. The paper attributes no mechanism to the difference; the indicator records arm assignment only.",
      source_name = "treatment (FF/VI)"
    )
  )

  compartmentData <- list(
    depot = list(analyte = "fluticasone furoate", units = "ug", specimen = "administration site", verified = TRUE),
    central = list(analyte = "fluticasone furoate", units = "ug", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "fluticasone furoate", units = "ug", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 714,
    n_studies = 3,
    n_observations = 2948,
    age_range = "41-88 years",
    age_median = "66 years",
    weight_range = "35.4-154 kg",
    weight_median = "72.0 kg",
    bmi_median = "25.3 kg/m2 (range 14.4-49.2)",
    sex_female_pct = 29,
    race_ethnicity = c(
      `White/Caucasian/European` = 67,
      `Asian - East Asian` = 16,
      `Asian - Japanese` = 13,
      `African American/African` = 3,
      `White - Arabic/North African` = 0.3
    ),
    smoker_pct = 40,
    fev1_pct_predicted_median = "39.7% (range 6.65-78.2)",
    disease_state = "chronic obstructive pulmonary disease (symptomatic, at risk of exacerbation)",
    dose_range = "fluticasone furoate 100 ug once daily by oral inhalation as FF/UMEC/VI 100/62.5/25 ug (single inhaler), FF/VI 100/25 ug + UMEC 62.5 ug, or FF/VI 100/25 ug",
    regions = "multinational; all East Asian and Japanese subjects were resident in China, Japan or Korea",
    notes = "Pooled Phase III studies FULFIL (CTT116853), IMPACT (CTT116855) and 200812 (Table 1). Demographics from Table 2 (FF dataset). 41% of FF observations were below the 10 pg/mL LLOQ and were handled with the NONMEM M3 method; 79% of trough (> 20 h post-dose) samples were BQL."
  )

  ini({
    # Table 3 reports each THETA on the log scale ('Ln estimate') and
    # untransformed ('Estimate'). The untransformed value is used; log() of it
    # reproduces the Ln column (log(513) = 6.240, log(1.36) = 0.307 vs 0.310,
    # log(268) = 5.591, log(111) = 4.710, log(0.0821) = -2.500).
    lcl <- log(513); label("Apparent inhaled clearance CL/F, non-Japanese, umeclidinium co-administered (L/h)") # Table 3 'CL/F (L/h)' Estimate 513 [493, 534], Ln 6.24, RSE 0.385%; Sect. 3.2 'typical value of FF CL/F was 513 L/h for a subject with COPD of non-Japanese heritage'
    lvc <- fixed(log(1.36)); label("Apparent central volume V2/F (L)") # Table 3 'V2/F (L)' 1.36 fixed (Ln 0.310 fixed); Sect. 3.2 'V2/F, Q/F, and V3/F were fixed to values estimated from the previous FF population pharmacokinetic model'
    lq <- fixed(log(268)); label("Apparent intercompartmental clearance Q/F (L/h)") # Table 3 'Q/F (L/h)' 268 fixed (Ln 5.59 fixed)
    lvp <- fixed(log(111)); label("Apparent peripheral volume V3/F (L)") # Table 3 'V3/F (L)' 111 fixed (Ln 4.71 fixed)
    lka <- log(0.0821); label("Absorption rate constant KA (1/h)") # Table 3 'KA (h-1)' 0.0821 [0.0805, 0.0837], Ln -2.50, RSE 0.604%

    # Categorical covariates enter as TVPK = theta_pop * theta^CAT (Eq. 2),
    # i.e. log-additive; the coefficient below is log(theta).
    e_race_japanese_cl <- -0.436; label("Log-scale effect of Japanese heritage on CL/F (unitless)") # Table 3 'Japanese heritage on CL/F' Ln -0.436 [-0.466, -0.406], Estimate 0.647, RSE 11.8%; Sect. 3.2 '35% lower'
    e_ffvi_cl <- 0.351; label("Log-scale effect of the FF/VI dual-therapy arm (no umeclidinium) on CL/F (unitless)") # Table 3 'FF/VI on CL/F' Ln 0.351 [0.321, 0.381], Estimate 1.42, RSE 11.0%; Sect. 3.2 '42% higher'

    # IIV: Table 3 'IIV, CV%' converted with omega^2 = log(1 + CV^2)
    # (log-normal exponential IIV). No OMEGA variances or correlations are
    # printed; etas are taken as independent. IIV was estimated on V2/F, Q/F
    # and V3/F although their THETAs were fixed (Table 3 lists a CV% for
    # each).
    etalcl ~ 0.3913 # Table 3 CL/F IIV 69.2 CV% -> log(1 + 0.692^2)
    etalvc ~ 2.8190 # Table 3 V2/F IIV 397 CV% -> log(1 + 3.97^2)
    etalq ~ 0.4685 # Table 3 Q/F IIV 77.3 CV% -> log(1 + 0.773^2)
    etalvp ~ 0.3782 # Table 3 V3/F IIV 67.8 CV% -> log(1 + 0.678^2)
    etalka ~ 0.4044 # Table 3 KA IIV 70.6 CV% -> log(1 + 0.706^2)

    # Residual error: Mehta 2020 does not report the residual-error model or
    # its magnitude anywhere (Sect. 2.3 states only SAEM with interaction and
    # the M3 method for BQL data). Declared at fixed(0); see the vignette
    # 'Assumptions and deviations'.
    propSd <- fixed(0); label("Proportional residual SD (fraction; 0 -- not reported in the source)")
  })

  model({
    cl <- exp(lcl + e_race_japanese_cl * RACE_JAPANESE + e_ffvi_cl * (1 - CONMED_UMECLIDINIUM) + etalcl)
    vc <- exp(lvc + etalvc)
    q <- exp(lq + etalq)
    vp <- exp(lvp + etalvp)
    ka <- exp(lka + etalka)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
