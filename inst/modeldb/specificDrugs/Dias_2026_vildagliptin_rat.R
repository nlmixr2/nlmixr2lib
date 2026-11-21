Dias_2026_vildagliptin_rat <- function() {
  description <- paste(
    "Preclinical (rat, male Wistar). Four-compartment population PK model",
    "for intravenous vildagliptin in healthy rats and in rats with",
    "alloxan-induced diabetes, fit jointly to plasma and to unbound skeletal",
    "muscle and liver interstitial concentrations measured by",
    "microdialysis (Dias 2026). The ODE states carry UNBOUND drug amounts",
    "per kg body weight: a central (plasma) compartment with linear",
    "elimination, a peripheral compartment, and muscle and liver",
    "compartments each linked to central by separate influx and efflux",
    "clearances. Diabetes has its own central-peripheral clearance Q1 and",
    "liver efflux clearance. Total plasma is unbound plasma divided by the",
    "unbound fraction 0.91. IIV on CL, Q1, the muscle and liver volumes and",
    "the liver efflux clearance; log-additive residual error per matrix.",
    sep = " "
  )
  reference <- paste(
    "Dias BB, Olivo LB, de Araujo BV. Preclinical Modeling and Simulation to",
    "Explore the Tissue/Plasma Exposure and Pharmacodynamic Effect of",
    "Vildagliptin in Diabetes Treatment. CPT Pharmacometrics Syst Pharmacol.",
    "2026;15(2):e70165. doi:10.1002/psp4.70165. PMCID: PMC12896384.",
    "Parameter estimates: Table 1. Structural ODEs: Equations 4-7 and",
    "Figure 2. Group sizes: Supporting Information Table S1.",
    sep = " "
  )
  vignette <- "Dias_2026_vildagliptin"
  units <- list(time = "h", dosing = "mg/kg", concentration = "mg/L")

  compartmentData <- list(
    central = list(analyte = "vildagliptin (unbound)", units = "mg/kg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "vildagliptin (unbound)", units = "mg/kg", specimen = "tissue", verified = TRUE),
    muscle = list(
      analyte = "vildagliptin (unbound, skeletal muscle interstitial fluid)",
      units = "mg/kg",
      specimen = "tissue",
      verified = TRUE
    ),
    liver = list(
      analyte = "vildagliptin (unbound, liver interstitial fluid)",
      units = "mg/kg",
      specimen = "tissue",
      verified = TRUE
    )
  )

  covariateData <- list(
    DIS_DIAB = list(
      description = "Diabetes status: 1 = rat with alloxan-induced diabetes; 0 = healthy rat.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (healthy)",
      source_name = "Status",
      notes = paste(
        "Experimental diabetes: alloxan 40 mg/kg i.v., animals diabetic when",
        "blood glucose > 300 mg/dL at 72 h, studied 4 weeks after induction",
        "(Methods 'Animal Model of Diabetic'). Methods 'Covariate Analysis':",
        "'Status was defined as 0 for healthy and 1 for diabetic animals'.",
        "Results: Q1 and Qout,liver were estimated 'separately for diabetic",
        "and healthy animals', so the indicator selects the diabetic Q1",
        "(2.61 L/h/kg) and liver efflux clearance (19.7 L/h/kg) in place of",
        "the healthy values (0.338 and 218 L/h/kg); each pair shares one eta."
      )
    )
  )

  population <- list(
    species = "rat (Wistar, male)",
    n_subjects = 37L,
    n_studies = 1L,
    age_range = "8 weeks at purchase",
    weight_range = "250-350 g",
    sex_female_pct = 0,
    disease_state = paste(
      "Healthy rats and rats with alloxan-induced diabetes",
      "(blood glucose > 300 mg/dL; studied 4 weeks after induction)."
    ),
    dose_range = "Single 50 mg/kg intravenous bolus.",
    regions = "Brazil (Federal University of Rio Grande do Sul, Porto Alegre)",
    notes = paste(
      "Separate animals per matrix (Table S1, animals/observations):",
      "plasma healthy 7/62 and diabetic 5/45; liver microdialysis healthy",
      "6/42 and diabetic 6/42; muscle microdialysis healthy 7/23 and",
      "diabetic 6/40; 37 animals and 254 observations in total. Plasma by",
      "carotid sampling at 0.083-6 h; microdialysate collected every 40 min",
      "up to 6 h at 1 uL/min and corrected for the in vivo relative recovery",
      "(31.2% +/- 5.0%). Plasma protein binding 9%. NONMEM 7.4, FOCE-I,",
      "ADVAN13."
    )
  )

  ini({
    # Structural parameters: Table 1 (estimate, %RSE), per kg body weight.
    lcl <- log(2.72); label("Clearance CL (L/h/kg)") # Table 1: CL = 2.72 L/h/kg (RSE 5%)
    lvc <- log(1.11); label("Central volume V1, unbound (L/kg)") # Table 1: V1 = 1.11 L/kg (RSE 9%)
    lq <- log(0.338); label("Central-peripheral clearance Q1, healthy rats (L/h/kg)") # Table 1: Q1,healthy = 0.338 L/h/kg (RSE 21%)
    lq_diabetic <- log(2.61); label("Central-peripheral clearance Q1, diabetic rats (L/h/kg)") # Table 1: Q1,diabetic = 2.61 L/h/kg (RSE 21%)
    lvp <- log(1.85); label("Peripheral volume V2 (L/kg)") # Table 1: V2 = 1.85 L/kg (RSE 12%)
    lclin_muscle <- log(2.09); label("Influx clearance central to muscle Qin,muscle (L/h/kg)") # Table 1: Qin,muscle = 2.09 L/h/kg (RSE 22%)
    lclef_muscle <- log(14.3); label("Efflux clearance muscle to central Qout,muscle (L/h/kg)") # Table 1: Qout,muscle = 14.3 L/h/kg (RSE 23%)
    lv_muscle <- log(2.55); label("Muscle compartment volume V3 (L/kg)") # Table 1: V3 = 2.55 L/kg (RSE 26%)
    lclin_liver <- log(2.5); label("Influx clearance central to liver Qin,liver (L/h/kg)") # Table 1: Qin,liver = 2.5 L/h/kg (RSE 31%)
    lclef_liver <- log(218); label("Efflux clearance liver to central Qout,liver, healthy rats (L/h/kg)") # Table 1: Qout,liver,healthy = 218 L/h/kg (RSE 53%)
    lclef_liver_diabetic <- log(19.7); label("Efflux clearance liver to central Qout,liver, diabetic rats (L/h/kg)") # Table 1: Qout,liver,diabetic = 19.7 L/h/kg (RSE 49%)
    lv_liver <- log(0.856); label("Liver compartment volume V4 (L/kg)") # Table 1: V4 = 0.856 L/kg (RSE 49%)
    fu <- fixed(0.91); label("Unbound fraction of vildagliptin in rat plasma (fraction)") # Methods 'Model Building': rat plasma protein binding of 9%

    # IIV: Table 1 footnote gives IIV as CV% = omega x 100 with omega the
    # SD of the random effect, so the variance is (CV/100)^2. The healthy
    # and diabetic Q1 (and Qout,liver) share one eta (same IIV row).
    etalcl ~ 0.036864 # Table 1: IIV CL = 19.2% (RSE 19%), shrinkage 16%
    etalq ~ 0.427716 # Table 1: IIV Q1 = 65.4% (RSE 20%), shrinkage 30%
    etalv_muscle ~ 0.586756 # Table 1: IIV V3 = 76.6% (RSE 29%), shrinkage 30%
    etalclef_liver ~ 0.5184 # Table 1: IIV Qout,liver = 72.0% (RSE 21%), shrinkage 40%
    etalv_liver ~ 3.290596 # Table 1: IIV V4 = 181.4% (RSE 26%), shrinkage 43%

    # Residual error: log-additive, separate for plasma, muscle and liver
    # (Methods 'Model Building'); encoded as SDs on the log scale.
    expSd <- 0.213; label("Log-additive residual SD, plasma (log scale)") # Table 1: residual error plasma = 0.213 (RSE 21%)
    expSd_Cmuscle <- 0.217; label("Log-additive residual SD, muscle microdialysate (log scale)") # Table 1: residual error muscle = 0.217 (RSE 10%)
    expSd_Cliver <- 0.159; label("Log-additive residual SD, liver microdialysate (log scale)") # Table 1: residual error liver = 0.159 (RSE 15%)
  })

  model({
    # Individual parameters (exponential IIV). Q1 and Qout,liver take the
    # diabetic estimate when DIS_DIAB = 1 (Results: 'The estimation of these
    # two parameters separately for diabetic and healthy animals').
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc)
    q <- exp(lq * (1 - DIS_DIAB) + lq_diabetic * DIS_DIAB + etalq)
    vp <- exp(lvp)
    clin_muscle <- exp(lclin_muscle)
    clef_muscle <- exp(lclef_muscle)
    v_muscle <- exp(lv_muscle + etalv_muscle)
    clin_liver <- exp(lclin_liver)
    clef_liver <- exp(lclef_liver * (1 - DIS_DIAB) + lclef_liver_diabetic * DIS_DIAB + etalclef_liver)
    v_liver <- exp(lv_liver + etalv_liver)

    # Equations 4-7: A1 = central, A2 = peripheral1, A3 = muscle,
    # A4 = liver (unbound amounts, mg/kg).
    d/dt(central) <- q / vp * peripheral1 + clef_muscle / v_muscle * muscle + clef_liver / v_liver * liver - (cl / vc + q / vc + clin_muscle / vc + clin_liver / vc) * central
    d/dt(peripheral1) <- q / vc * central - q / vp * peripheral1
    d/dt(muscle) <- clin_muscle / vc * central - clef_muscle / v_muscle * muscle
    d/dt(liver) <- clin_liver / vc * central - clef_liver / v_liver * liver

    # Concentrations: amounts in mg/kg over volumes in L/kg give mg/L, the
    # unit of Figures 1 and 4.
    Cu <- central / vc
    Cmuscle <- muscle / v_muscle
    Cliver <- liver / v_liver
    # Total plasma: the model was fit to unbound plasma concentrations
    # derived from the measured totals with 9% protein binding (Methods
    # 'Model Building'); under log-additive error the residual SD is the
    # same on either scale.
    Cc <- Cu / fu

    Cc ~ lnorm(expSd)
    Cmuscle ~ lnorm(expSd_Cmuscle)
    Cliver ~ lnorm(expSd_Cliver)
  })
}
