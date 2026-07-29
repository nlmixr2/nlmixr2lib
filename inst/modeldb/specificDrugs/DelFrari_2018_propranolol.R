DelFrari_2018_propranolol <- function() {
  description <- paste(
    "One-compartment population PK model for oral propranolol in infants",
    "(aged 50-243 days, 3.6-9.7 kg) with proliferative Infantile Hemangiomas",
    "(Del Frari 2018). First-order absorption and first-order elimination;",
    "apparent oral clearance CL/F scales with body weight using a fixed",
    "allometric exponent of 0.75 and a reference weight of 6.3 kg (the",
    "median weight pooled across visits D1-D84). Apparent volume V/F has",
    "no covariate effect (the paper tested but did not retain weight on",
    "V/F). Between-subject variability is retained on CL/F and Ka only;",
    "BSV on V/F was dropped from the final model (large 95% CI including 0",
    "and 62.3% eta-shrinkage). Residual error is proportional."
  )
  reference <- paste(
    "Del Frari L, Leaute-Labreze C, Guibaud L, Barbarot S, Lacour J-P,",
    "Chaumont C, Delarue A, Voisard J-J, Brunner V. Propranolol",
    "pharmacokinetics in infants treated for Infantile Hemangiomas requiring",
    "systemic therapy: Modeling and dosing regimen recommendations.",
    "Pharmacol Res Perspect. 2018;6(3):e00399. doi:10.1002/prp2.399.",
    sep = " "
  )
  vignette <- "DelFrari_2018_propranolol"
  units    <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Time-varying body weight in kg. Drives fixed-exponent allometric",
        "scaling on CL/F (exponent 0.75) with reference weight 6.3 kg (the",
        "median weight observed across the pooled cohort from D1 to D84,",
        "stated explicitly in Del Frari 2018 Results section 3.4 below",
        "Equation 5). The paper tested weight on V/F as well but did not",
        "retain it; V/F is therefore weight-independent here.",
        sep = " "
      ),
      source_name        = "WGT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 22L,
    n_studies      = 1L,
    age_range      = "50-243 days (postnatal age at PK assessment); 35-150 days at study inclusion",
    age_median     = "104 days at inclusion (range 50-151)",
    weight_range   = "3.6-9.7 kg across visits D1-D84",
    weight_median  = "6.3 kg (pooled median across visits D1-D84; reference weight for allometric scaling)",
    sex_female_pct = 27.3,
    disease_state  = paste(
      "Infants with proliferating Infantile Hemangiomas (IH) requiring",
      "systemic therapy. Patients were stratified into two groups by age at",
      "inclusion (Group 1: 35-90 days, n = 10, PK at week 4; Group 2:",
      "91-150 days, n = 12, PK at week 12). EudraCT 2009-018102-22.",
      sep = " "
    ),
    dose_range     = paste(
      "Oral solution of propranolol hydrochloride dosed twice daily (BID).",
      "Titration: 1 mg/kg/day in week 1, 2 mg/kg/day in week 2, then target",
      "3 mg/kg/day BID (1.5 mg/kg per dose) for weeks 3-12.",
      "PK assessments performed at steady state at the target dose with a",
      "regular 12-hour dosing interval on the assessment day.",
      sep = " "
    ),
    regions        = "France (4 hospitals: Bordeaux, Lyon, Nantes, Nice)",
    notes          = paste(
      "23 infants enrolled; 22 were evaluable for popPK (6 females, 16 males).",
      "167 plasma concentrations contributed to the analysis (21 trough at D7",
      "+ 18 trough at D14 + 60 + 68 across full PK profiles at D28 / D84).",
      "One D84 outlier concentration (448 ng/mL, weighted residual +4.5) was",
      "excluded after graphical evaluation. Bioanalytical method: LC-MS/MS,",
      "quantification range 0.50-250 ng/mL (Del Frari 2018 Methods 2.2).",
      sep = " "
    )
  )

  ini({
    # Structural parameters from Del Frari 2018 Table 5 (final model).
    # CL/F is the apparent oral clearance at the reference weight 6.3 kg
    # (median weight pooled across visits D1-D84; cited explicitly below
    # Equation 5 in section 3.4). V/F is weight-independent in the final
    # model. Ka is first-order absorption.
    lcl <- log(19.3) ; label("Apparent oral clearance CL/F at WT = 6.3 kg (L/h)")  # Del Frari 2018 Table 5 (CL/F = 19.3 L/h; 95% CI 15.4-23.2, RSE 10.3%)
    lvc <- log(122)  ; label("Apparent volume of distribution V/F (L)")            # Del Frari 2018 Table 5 (V/F = 122 L; 95% CI 104-140, RSE 7.32%)
    lka <- log(0.993); label("First-order absorption rate constant Ka (1/h)")      # Del Frari 2018 Table 5 (Ka = 0.993 1/h; 95% CI 0.558-1.43, RSE 22.4%)

    # Fixed allometric exponent on body weight for CL/F. The paper compared
    # an estimated vs fixed-at-0.75 exponent and retained the fixed-0.75
    # form on AIC grounds (Del Frari 2018 Results section 3.4 paragraph 2).
    e_wt_cl <- fixed(0.75); label("Allometric exponent of (WT/6.3) on CL/F (unitless; fixed)")  # Del Frari 2018 Results 3.4 + Eq 5

    # IIV from Table 5 final-model column. NONMEM exponential IIV with
    # variance reported directly as omega^2; the CV% column equals
    # sqrt(omega^2) * 100 (44.2% CV on CL, 132% CV on Ka). BSV on V/F was
    # tested in the base model but dropped from the final model because
    # the 95% CI of omega^2 included 0 and the eta-shrinkage was 62.3%
    # (Del Frari 2018 Results 3.4 paragraph 4 and Table 4 vs Table 5).
    etalcl ~ 0.195   # Del Frari 2018 Table 5 (omega^2 on CL = 0.195; 44.2% CV; 95% CI 0.0492-0.341, RSE 38.2%)
    etalka ~ 1.75    # Del Frari 2018 Table 5 (omega^2 on Ka = 1.75;  132% CV; 95% CI 0.370-3.13, RSE 40.2%)

    # Proportional residual error. NONMEM reports variance sigma^2; the
    # proportional SD on the linear concentration scale is sqrt(sigma^2).
    propSd <- sqrt(0.0953); label("Proportional residual error (fraction)")  # Del Frari 2018 Table 5 (sigma^2 = 0.0953; 30.9% CV; 95% CI 0.0675-0.123, RSE 14.9%) -> propSd = sqrt(0.0953) = 0.309
  })

  model({
    # Individual PK parameters. Apparent CL/F scales allometrically with
    # body weight (reference 6.3 kg, fixed exponent 0.75); V/F has no
    # covariate effect; Ka has IIV but no covariate. NONMEM ADVAN2 TRANS2
    # with apparent parameters (oral bioavailability is not identifiable
    # from oral-only data; CL/F and V/F absorb the F term).
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) * (WT / 6.3)^e_wt_cl
    vc <- exp(lvc)

    # Micro-constant for first-order elimination.
    kel <- cl / vc

    # One-compartment ODEs with first-order oral absorption.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Plasma concentration. With dose in mg and vc in L, central/vc has
    # units mg/L; multiply by 1000 to express concentration in ng/mL
    # matching the bioanalytical and reported units in the paper.
    Cc <- central / vc * 1000

    Cc ~ prop(propSd)
  })
}
