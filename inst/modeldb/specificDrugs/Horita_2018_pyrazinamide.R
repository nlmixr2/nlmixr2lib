Horita_2018_pyrazinamide <- function() {
  description <- "One-compartment population pharmacokinetic model with three-compartment transit absorption followed by first-order absorption and first-order elimination for oral pyrazinamide in Ghanaian children with active tuberculosis (Horita 2018); allometric weight scaling on V/F (estimated exponent 0.677) and CL/F (estimated exponent 0.735) normalised to the cohort median 14.3 kg."
  reference <- "Horita Y, Alsultan A, Kwara A, Antwi S, Enimil A, Ortsin A, Dompreh A, Yang H, Wiesner L, Peloquin CA. Evaluation of the Adequacy of WHO Revised Dosages of the First-Line Antituberculosis Drugs in Children with Tuberculosis Using Population Pharmacokinetic Modeling and Simulations. Antimicrob Agents Chemother. 2018;62(9):e00008-18. doi:10.1128/AAC.00008-18"
  vignette <- "Horita_2018_pyrazinamide"
  units <- list(time = "hour", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric scaling on CL/F (estimated exponent 0.735) and V/F (estimated exponent 0.677) per Horita 2018 Table 4. Reference weight is the cohort median 14.3 kg (Table 1); see the vignette Errata for the reference-weight derivation.",
      source_name        = "WT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 113L,
    n_studies      = 1L,
    age_range      = "3 months to 14 years (median 5.00 years, IQR 2.17 to 8.25)",
    age_median     = "5.00 years",
    weight_range   = "5-30 kg (median 14.3, IQR 9.70 to 20.1)",
    weight_median  = "14.3 kg",
    sex_female_pct = 44.2,
    hiv_positive_pct = 52.2,
    disease_state  = "Ghanaian children with active tuberculosis (HIV-positive and HIV-negative). 21.2% under 2 years of age.",
    dose_range     = "Pyrazinamide 30-40 mg/kg orally daily (median 24.7 mg/kg, IQR 22.6-29.7; note: median was below the WHO 2010 recommended range). Administered as part of standard four-drug anti-TB regimen.",
    regions        = "Ghana (Komfo Anokye Teaching Hospital, Kumasi).",
    notes          = "Patients enrolled October 2012-August 2015. PK sampling after at least 4 weeks of anti-TB treatment (steady state). Blood samples at 0, 1, 2, 4, 8 h postdose. PZA concentrations 0.20-80 ug/mL by LC-MS/MS. Three children with only-BLQ values (malabsorption group) and one child with no observed peak were excluded from model building. ClinicalTrials.gov NCT01687504. Demographics from Horita 2018 Table 1; structural model and parameters from Table 4."
  )

  ini({
    # Structural PK parameters -- Horita 2018 Table 4 final population pharmacokinetic
    # model for pyrazinamide. Typical values are at the cohort median weight 14.3 kg.
    lktr <- log(12.5);  label("Transit-chain rate constant Ktr (1/h)")                       # Table 4: Ktr = 12.5 1/h (RSE 11%)
    lka  <- log(5.28);  label("First-order absorption rate constant ka from last transit to central (1/h)")  # Table 4: ka = 5.28 1/h (RSE 7%)
    lvc  <- log(13.1);  label("Apparent central volume V/F at WT = 14.3 kg (L)")             # Table 4: V/F = 13.1 L (RSE 4%)
    lcl  <- log(1.6);   label("Apparent oral clearance CL/F at WT = 14.3 kg (L/h)")          # Table 4: CL/F = 1.6 L/h (RSE 4%)

    # Allometric exponents on body weight -- ESTIMATED (paper reports RSE), not fixed.
    # The Table 4 IIV column lists 'Fixed' for these rows, which means no IIV is
    # estimated on the exponent (typical for allometric exponents) -- not that the
    # exponent value itself is held fixed during fitting.
    e_wt_cl <- 0.735; label("Allometric exponent on CL/F (estimated, unitless)")    # Table 4: Exponent (BW on CL/F) = 0.735 (RSE 10%)
    e_wt_vc <- 0.677; label("Allometric exponent on V/F  (estimated, unitless)")    # Table 4: Exponent (BW on V/F)  = 0.677 (RSE 11%)

    # Inter-individual variability. Table 4 IIV column reports 'omega (CV%)' on the
    # log scale (variance = omega^2). The paper additionally reports an MTT IIV
    # which is implicit in the Ktr IIV when N is fixed (see vignette Errata).
    etalktr ~ 0.207   # Table 4: 0.455 (48.0% CV)  -- 0.455^2 = 0.207 (Ktr; absorbs the paper's MTT-IIV in the fixed-N implementation)
    etalka  ~ 4.20    # Table 4: 2.05 (811.5% CV)  -- 2.05^2  = 4.20  (ka; very large IIV consistent with the paper's report)
    etalvc  ~ 0.137   # Table 4: 0.370 (38.3% CV)  -- 0.370^2 = 0.137 (V/F)
    etalcl  ~ 0.148   # Table 4: 0.385 (40.0% CV)  -- 0.385^2 = 0.148 (CL/F)

    # Combined residual error. Table 4: 'Constant a' = 0.268 (RSE 23%), 'Slope b'
    # = 0.112 (RSE 7%).
    addSd  <- 0.268; label("Additive residual SD (ug/mL)")                          # Table 4: constant a = 0.268
    propSd <- 0.112; label("Proportional residual SD (fraction)")                   # Table 4: slope b    = 0.112
  })

  model({
    # Individual PK parameters with allometric weight scaling (reference 14.3 kg).
    ktr <- exp(lktr + etalktr)
    ka  <- exp(lka  + etalka)
    cl  <- exp(lcl  + etalcl) * (WT / 14.3)^e_wt_cl
    vc  <- exp(lvc  + etalvc) * (WT / 14.3)^e_wt_vc

    kel <- cl / vc

    # Three transit-compartment chain at rate Ktr followed by first-order absorption
    # at rate ka into central (Horita 2018 Results 'PZA' paragraph 1: 'a one-compartment
    # model with transit compartment absorption and first-order elimination... the
    # number of transit compartments was about 3'). Follows the Bienczak 2016 nevirapine
    # pattern: depot -> transit1 -> transit2 -> transit3 -> central with the depot-to-
    # transit1 and inter-transit transitions at rate Ktr, and the final transit3-to-
    # central transition at rate ka.
    d/dt(depot)    <- -ktr * depot
    d/dt(transit1) <-  ktr * depot    - ktr * transit1
    d/dt(transit2) <-  ktr * transit1 - ktr * transit2
    d/dt(transit3) <-  ktr * transit2 - ka  * transit3
    d/dt(central)  <-  ka  * transit3 - kel * central

    # Concentration: dose mg / V/F L -> mg/L = ug/mL.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
