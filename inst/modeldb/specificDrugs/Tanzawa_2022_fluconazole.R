Tanzawa_2022_fluconazole <- function() {
  description <- "One-compartment population PK model for fluconazole formed by first-order conversion of intravenous fosfluconazole (dosed as fluconazole equivalent), with allometric weight scaling and power effects of postmenstrual age, serum creatinine and alkaline phosphatase on clearance, in extremely low-birth-weight infants (Tanzawa 2022)"
  reference <- "Tanzawa A, Saito J, Shoji K, Kojo Y, Funaki T, Maruyama H, Isayama T, Ito Y, Nakamura H, Yamatani A. Fluconazole Population Pharmacokinetics after Fosfluconazole Administration and Dosing Optimization in Extremely Low-Birth-Weight Infants. Microbiol Spectr. 2022;10(2):e01952-21. doi:10.1128/spectrum.01952-21"
  vignette <- "Tanzawa_2022_fluconazole"
  units <- list(time = "h", dosing = "mg", concentration = "ug/mL")

  # Fosfluconazole (a phosphate-ester prodrug) is infused intravenously and is
  # converted to fluconazole by alkaline phosphatase at the first-order rate kc
  # (Results 'one-compartment analysis with first-order conversion'). Fosfluconazole
  # was not measured, so the prodrug is carried as an input compartment `depot`
  # (dose in mg of fluconazole equivalent) draining into `central`, following the
  # `Hornik_2019_methylprednisolone_il6.R` precedent for an unmeasured
  # intravenous prodrug. Give infusions into `depot`, not `central`.
  dosing <- "depot"

  compartmentData <- list(
    depot = list(
      analyte = "fosfluconazole (prodrug, as fluconazole equivalent)",
      units = "mg",
      specimen = "administration site",
      verified = TRUE
    ),
    central = list(analyte = "fluconazole", units = "mg", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Current body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying. Unnormalized allometric scaling with fixed exponents (Methods: CL scaled by WT^0.75, V by WT^1.0), so CL is reported in L/h/kg^0.75 and V in L/kg. Table 1 reports weight in grams (median 750 g, IQR 580-923 g); convert to kg.",
      source_name = "WT"
    ),
    PAGE = list(
      description = "Postmenstrual age",
      units = "weeks",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying. Power effect on CL referenced to 29 weeks: (PAGE / 29)^1.52 (Table 2 / Table 3 footnote a). Carried in WEEKS, as the paper states it, per the PAGE register note on neonatal models. Modelled range 22.9-37.2 weeks (Table 5 footnote b).",
      source_name = "PMA"
    ),
    CREAT = list(
      description = "Serum creatinine (enzymatic assay)",
      units = "mg/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying. Power effect on CL referenced to 0.64 mg/dL (the study-period median, Table 1): (CREAT / 0.64)^-0.17. Modelled range 0.24-2.3 mg/dL (Table 5 footnote a). Missing values were imputed from the closest available measurement (carry-forward or back-fill).",
      source_name = "SCr"
    ),
    ALP = list(
      description = "Serum alkaline phosphatase",
      units = "IU/L",
      type = "continuous",
      reference_category = NULL,
      notes = "Time-varying. Power effect on CL referenced to 958 IU/L (Table 1 median): (ALP / 958)^0.10. Table 1 prints the unit as 'IU/mL', a slip for IU/L (the Discussion quotes the preterm range as 100 to 2,000 IU/L).",
      source_name = "ALP"
    )
  )

  covariatesDataExcluded <- list(
    PNA = list(
      description = "Postnatal age",
      units = "days",
      type = "continuous",
      notes = "Significant on CL in the univariate screen (Table 2, reference 6 days, dOFV -234.3) but excluded because it was strongly correlated with PMA (Results). Also screened on V (dOFV -7.6), not retained.",
      source_name = "PNA"
    ),
    HT = list(
      description = "Current height",
      units = "cm",
      type = "continuous",
      notes = "Significant on CL in the univariate screen (Table 2, reference 31 cm, dOFV -162.9) but removed in the multivariate full model for its limited OFV impact (Results). No effect on V (dOFV 0). Current (time-varying) height, not baseline height.",
      source_name = "HT"
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/dL",
      type = "continuous",
      notes = "Screened on CL (Table 2, reference 2.7 g/dL, dOFV -2.3, not significant) and on V (dOFV -12.8) but not retained on V because of a 66.6% RSE and a 95% CI spanning the null (Results). SCr (dOFV -21.8) and ALP (dOFV -13.5) on V were rejected for the same reason.",
      source_name = "Alb"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 18,
    n_studies = 1,
    n_samples = "442 scavenged serum samples assayed; 64 (14.4%) below the 0.0031 ug/mL LLOQ were excluded",
    age_range = "GA 23.2 weeks (IQR 23.1-26.4); PMA 28.5 weeks (IQR 26.3-31.7) over the study, modelled range 22.9-37.2 weeks; dosing started at PNA 0 days",
    weight_range = "Birth weight 748 g (IQR 521-866); current weight 750 g (IQR 580-923)",
    sex_female_pct = 38.9,
    race_ethnicity = "Japanese (single centre in Tokyo; race not tabulated)",
    disease_state = "Extremely low-birth-weight infants (< 1,000 g) with central vascular access receiving fosfluconazole prophylaxis against invasive fungal infection in a NICU",
    dose_range = "Fosfluconazole IV 3 mg/kg as fluconazole equivalent, every 72 h in weeks 1-2 of life, every 48 h in weeks 3-4, every 24 h from week 5; no loading dose",
    regions = "Japan (National Center for Child Health and Development, Tokyo)",
    renal_function = "Serum creatinine 0.64 mg/dL (IQR 0.48-1.07) over the study; modelled range 0.24-2.3 mg/dL",
    co_medication = "Amikacin 94.4%, indomethacin 44.4%, vancomycin 16.7% (Results)",
    notes = "Table 1 baseline demographics. Prospective single-centre observational study, March-August 2020, opportunistic (scavenged residual-sample) sampling; median 42.2 h (IQR 18.5-61.5) after the latest dose. Fitted with Phoenix NLME 8.2."
  )

  ini({
    # Structural parameters (Table 3, final model). CL and V are per allometric
    # weight: CL in L/h/kg^0.75 and V in L/kg, with WT in kg (unnormalized).
    lcl <- log(0.011); label("Clearance per kg^0.75 at PMA 29 weeks, SCr 0.64 mg/dL, ALP 958 IU/L (L/h/kg^0.75)") # Table 3 final model theta_CL = 0.011 (RSE 10.59%)
    lvc <- log(0.95); label("Volume of distribution per kg (L/kg)") # Table 3 final model theta_V = 0.95 (RSE 5.40%)
    lka <- fixed(log(0.43)); label("Fosfluconazole-to-fluconazole conversion rate constant kc (1/h)") # Table 3 theta_kc = 0.43 in base and final model; bootstrap 95% CI 0.43 to 0.43 (zero width) -> held constant

    # Allometric exponents, fixed by design (Methods: 'CL was scaled using
    # allometric weight (WT^0.75) and V was scaled using weight (WT^1.0) as a
    # default structure model').
    e_wt_cl <- fixed(0.75); label("Allometric exponent of weight on CL (unitless)") # Methods, PK analysis
    e_wt_vc <- fixed(1.0); label("Allometric exponent of weight on V (unitless)") # Methods, PK analysis

    # Covariate power exponents on CL (Table 3 final model; Table 3 footnote a:
    # CL = theta_CL * (PMA/29)^theta_PMA * (SCr/0.64)^theta_SCr * (ALP/958)^theta_ALP).
    e_page_cl <- 1.52; label("Power exponent of postmenstrual age on CL (unitless)") # Table 3 theta_PMA = 1.52 (bootstrap median 1.53)
    e_creat_cl <- -0.17; label("Power exponent of serum creatinine on CL (unitless)") # Table 3 theta_SCr = -0.17 (bootstrap 95% CI -0.31 to -0.08)
    e_alp_cl <- 0.10; label("Power exponent of alkaline phosphatase on CL (unitless)") # Table 3 theta_ALP = 0.10 (bootstrap median 0.10)

    # IIV (exponential). Table 3 reports CV% only; omega^2 = log(1 + CV^2).
    # No IIV on kc (Results: removed because its shrinkage was 0.9).
    etalcl ~ 0.07600 # Table 3 final model IIV CL = 28.1 CV% -> log(1 + 0.281^2)
    etalvc ~ 0.02466 # Table 3 final model IIV V = 15.8 CV% -> log(1 + 0.158^2)

    # Combined proportional + additive residual error (Results: 'The model
    # combining proportional and additive errors best described the residual
    # variability').
    propSd <- 0.14; label("Proportional residual error (fraction)") # Table 3 final model proportional = 14.0 CV%
    addSd <- 0.068; label("Additive residual error (ug/mL)") # Table 3 final model additive = 0.068 ug/mL
  })
  model({
    # Covariate ratios, referenced to the study medians (Table 2; Table 3 footnote a)
    page_ratio <- PAGE / 29
    creat_ratio <- CREAT / 0.64
    alp_ratio <- ALP / 958

    # Individual parameters: CL in L/h and V in L after allometric weight scaling
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * WT^e_wt_cl * page_ratio^e_page_cl * creat_ratio^e_creat_cl * alp_ratio^e_alp_cl
    vc <- exp(lvc + etalvc) * WT^e_wt_vc

    kel <- cl / vc

    # depot: fosfluconazole as fluconazole equivalent (mg); central: fluconazole (mg)
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # Dose in mg and volume in L give mg/L = ug/mL
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
