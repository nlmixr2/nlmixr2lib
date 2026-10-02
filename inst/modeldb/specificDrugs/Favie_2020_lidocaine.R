Favie_2020_lidocaine <- function() {
  description <- paste(
    "One-compartment population PK model for intravenous lidocaine with a",
    "sequential one-compartment model for its metabolite monoethylglycinexylidide",
    "(MEGX) in preterm and (near-)term neonates treated for seizures, with and",
    "without therapeutic hypothermia (Favie 2020). All lidocaine elimination",
    "feeds the MEGX compartment; because the fraction converted to MEGX was",
    "unknown, the MEGX clearance and volume are apparent values relative to that",
    "fraction. Allometric body-weight scaling (fixed exponents 0.75 on clearance",
    "and 1 on volume, reference 3.5 kg) on all four disposition parameters, a",
    "linear postmenstrual-age effect on both clearances (reference 280 days =",
    "40 weeks), and a linear time-varying body-temperature effect on lidocaine",
    "clearance only (reference 36.5 degC)."
  )
  reference <- paste(
    "Favie LMA, Huitema ADR, van den Broek MPH, Rademaker CMA, de Haan TR,",
    "van Straaten HLM, Simons SHP, Rijken M, Nuytemans DHGM, Egberts TCG,",
    "Groenendaal F; PharmaCool study group. Lidocaine as treatment for neonatal",
    "seizures: Evaluation of previously developed population pharmacokinetic",
    "models and dosing regimen. Br J Clin Pharmacol. 2020;86(1):75-84.",
    "doi:10.1111/bcp.14136.",
    sep = " "
  )
  vignette <- "Favie_2020_lidocaine"
  units <- list(
    time = "h",
    dosing = "mg",
    concentration = "mg/L"
  )

  covariateData <- list(
    WT = list(
      description = "Body weight. Allometric scaling of lidocaine and MEGX clearance (exponent 0.75) and volume (exponent 1) with reference 3.5 kg.",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Table 4 footnote c states the typical values are 'for neonate with a birth weight of 3.5 kg and PMA 40 weeks', so the body-size descriptor is birth weight; the Methods call it 'BW'. Supply the birth weight (a constant per neonate over the few days of lidocaine treatment).",
      source_name = "BW"
    ),
    PAGE = list(
      description = "Postmenstrual age (gestational age + postnatal age). Linear effect on lidocaine clearance (0.69%/day) and MEGX clearance (0.35%/day) about 280 days (40 weeks).",
      units = "weeks",
      type = "continuous",
      reference_category = NULL,
      notes = "Table 4 final-model equations write the effect as (1 + theta * (PMA - 280)) with PMA in DAYS (280 days = the 'PMA 40 weeks' of footnote c). This model takes PAGE in WEEKS, the scale the neonatal PAGE precedents use, and converts inside model() as PAGE * 7 - 280. Time-varying in principle; over a ~28 h lidocaine course it changes by < 0.2 weeks. Observed range 25 - 42.7 weeks (Discussion); the linear form reaches zero clearance at 19.3 weeks and must not be extrapolated below the observed range.",
      source_name = "PMA"
    ),
    BODYTEMP = list(
      description = "Body temperature. Linear effect on lidocaine clearance (7.26%/degC) about 36.5 degC; no effect on MEGX clearance.",
      units = "degC",
      type = "continuous",
      reference_category = NULL,
      notes = "TIME-VARYING. Table 4 footnote d: 'In neonates treated with TH, TEMP was set to 33.5 degC during TH with rewarming at 0.4 degC/h. Normothermia for all neonates was set to 36.5 degC.' So the covariate is the protocolised (not measured) temperature: 36.5 degC for normothermic neonates; 33.5 degC during the 72 h of therapeutic hypothermia, rising linearly at 0.4 degC/h back to 36.5 degC (7.5 h). Methods: 'Body temperature was tested as a continuous variable using a dynamic model as described previously' (ref 36).",
      source_name = "TEMP"
    )
  )

  compartmentData <- list(
    central = list(analyte = "lidocaine", units = "mg", specimen = "plasma", verified = TRUE),
    central_megx = list(
      analyte = "monoethylglycinexylidide (MEGX)",
      units = "mg",
      specimen = "plasma",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 159L,
    n_studies = 4L,
    age_range = "neonates; gestational age mean 37.0 (SD 4.84) weeks; postmenstrual age 25 - 42.7 weeks",
    weight_range = "not reported as a range; mean 2.89 (SD 1.05) kg",
    sex_female_pct = 45.9,
    disease_state = "Neonatal seizures refractory to midazolam and/or phenobarbital, treated with lidocaine as second- or third-line antiepileptic drug; 50 (31.4%) preterm (GA < 36 weeks), 109 (near-)term, of whom 49 (30.8% of all) received therapeutic hypothermia (33.5 degC for 72 h) for hypoxic-ischaemic encephalopathy.",
    dose_range = "Continuous intravenous lidocaine infusion per local protocols, including the weight-banded regimen of Table 1 (2 mg/kg bolus over 10 min, then a 4 h [3.5 h under hypothermia] loading phase of 5-7 mg/kg/h and two 12 h maintenance phases at one half and one quarter of the loading rate). Lidocaine hydrochloride doses were converted to lidocaine base for the analysis.",
    regions = "The Netherlands (multicentre)",
    notes = "Table 2/Table 3. Pooled from clinical care cohort 1 (2004-2008, n = 46, Utrecht), the SHIVER study (2008-2010, n = 21, all hypothermia), the PharmaCool study (2010-2014, n = 22, all hypothermia) and clinical care cohort 2 (2010-2018, n = 70, 6 hypothermia). 444 lidocaine and/or MEGX samples; LC-MS/MS with LLQ 0.2 mg/L for both analytes; a value below LLQ for one compound was set to LLQ/2 (0.1 mg/L). NONMEM 7.3; precision by sampling importance resampling. Male 86 (54.1%)."
  )

  ini({
    # ------------------------------------------------------------------
    # Lidocaine disposition. Table 4 'Lidocaine' column; typical values
    # for a 3.5 kg neonate at PMA 40 weeks (footnote c) at 36.5 degC.
    lcl <- log(1.77)
    label("Log typical lidocaine clearance CL at WT = 3.5 kg, PMA = 40 weeks, BODYTEMP = 36.5 degC (L/h)") # Table 4 row 'Cl, l/h' Lidocaine = 1.77 (SIR 95% CI 1.63 - 2.03)
    lvc <- log(9.32)
    label("Log typical lidocaine volume of distribution V at WT = 3.5 kg (L)") # Table 4 row 'V, l' Lidocaine = 9.32 (SIR 95% CI 8.49 - 9.63)

    # ------------------------------------------------------------------
    # MEGX disposition. Table 4 'MEGX' column; footnote a: 'MEGX
    # estimates are relative to formation fraction F', i.e. CL_MEGX/F and
    # V_MEGX/F (Table 4 final-model equations).
    lcl_megx <- log(1.51)
    label("Log typical apparent MEGX clearance CL_MEGX/F at WT = 3.5 kg, PMA = 40 weeks (L/h)") # Table 4 row 'Cl, l/h' MEGX = 1.51 (SIR 95% CI 1.37 - 1.73)
    lvc_megx <- log(15.8)
    label("Log typical apparent MEGX volume of distribution V_MEGX/F at WT = 3.5 kg (L)") # Table 4 row 'V, l' MEGX = 15.8 (SIR 95% CI 13.6 - 18.7)

    # ------------------------------------------------------------------
    # Allometric exponents. Methods 2.4: 'BW was used as a descriptor for
    # body size and was related to PK parameters using allometric
    # relationships with an exponent of 0.75 on clearance and an exponent
    # of 1 on volume of distribution' -- fixed, not estimated.
    e_wt_cl <- fixed(0.75)
    label("Allometric exponent of body weight on lidocaine clearance (unitless)") # Methods 2.4; Table 4 final-model equation (BW/3.5)^0.75
    e_wt_vc <- fixed(1)
    label("Allometric exponent of body weight on lidocaine volume (unitless)") # Methods 2.4; Table 4 final-model equation (BW/3.5)^1
    e_wt_cl_megx <- fixed(0.75)
    label("Allometric exponent of body weight on apparent MEGX clearance (unitless)") # Methods 2.4; Table 4 final-model equation (BW/3.5)^0.75
    e_wt_vc_megx <- fixed(1)
    label("Allometric exponent of body weight on apparent MEGX volume (unitless)") # Methods 2.4; Table 4 final-model equation (BW/3.5)^1

    # ------------------------------------------------------------------
    # Linear covariate effects on clearance.
    e_page_cl <- 0.0069
    label("Linear postmenstrual-age effect on lidocaine clearance (fraction per day of PMA)") # Table 4 row 'PMA on Cl, %/d' Lidocaine = 0.690 (SIR 95% CI 0.581 - 0.837); equation 0.0069 * (PMA - 280)
    e_page_cl_megx <- 0.0035
    label("Linear postmenstrual-age effect on apparent MEGX clearance (fraction per day of PMA)") # Table 4 row 'PMA on Cl, %/d' MEGX = 0.350 (SIR 95% CI 0.114 - 0.805); equation 0.0035 * (PMA - 280)
    e_bodytemp_cl <- 0.0726
    label("Linear body-temperature effect on lidocaine clearance (fraction per degC)") # Table 4 row 'TEMP on Cl, %/C' Lidocaine = 7.26 (SIR 95% CI 1.63 - 11.2); equation 0.0726 * (TEMP - 36.5); MEGX 'NA'

    # ------------------------------------------------------------------
    # Inter-individual variability. Methods 2.4: 'Interindividual
    # variability was modelled using a proportional model'; Table 4
    # reports variances with RSD = sqrt(variance) (sqrt(0.231) = 48.1%),
    # i.e. exponential etas on the log scale. No covariances reported.
    etalcl ~ 0.231 # Table 4 row 'Cl, variance (RSD)' Lidocaine = 0.231 (48.1%)
    etalvc ~ 0.0673 # Table 4 row 'V, variance (RSD)' Lidocaine = 0.0673 (25.9%)
    etalcl_megx ~ 0.237 # Table 4 row 'Cl, variance (RSD)' MEGX = 0.237 (48.7%)
    etalvc_megx ~ 0.478 # Table 4 row 'V, variance (RSD)' MEGX = 0.478 (69.1%)

    # ------------------------------------------------------------------
    # Residual error. Methods 2.4: 'both proportional and additive error
    # models were used ... in which the additive error was fixed on
    # LLQ/2'; separate error models for lidocaine and MEGX.
    addSd <- fixed(0.1)
    label("Additive residual error SD on lidocaine concentration (mg/L)") # Table 4 row 'Additional, mg/L' Lidocaine = 0.1 (fixed) = LLQ/2
    propSd <- sqrt(0.0379)
    label("Proportional residual error SD on lidocaine concentration (fraction)") # Table 4 row 'Proportional, variance (RSD)' Lidocaine = 0.0379 (19.5%)
    addSd_megx <- fixed(0.1)
    label("Additive residual error SD on MEGX concentration (mg/L)") # Table 4 row 'Additional, mg/L' MEGX = 0.1 (fixed) = LLQ/2
    propSd_megx <- sqrt(0.0550)
    label("Proportional residual error SD on MEGX concentration (fraction)") # Table 4 row 'Proportional, variance (RSD)' MEGX = 0.0550 (23.5%)
  })

  model({
    # Postmenstrual age in days, the scale of the Table 4 equations
    # (280 days = 40 weeks).
    pma_days <- PAGE * 7

    # Table 4 final-model equations.
    cl <- exp(lcl + etalcl) *
      (WT / 3.5)^e_wt_cl *
      (1 + e_page_cl * (pma_days - 280)) *
      (1 + e_bodytemp_cl * (BODYTEMP - 36.5))
    vc <- exp(lvc + etalvc) * (WT / 3.5)^e_wt_vc
    cl_megx <- exp(lcl_megx + etalcl_megx) *
      (WT / 3.5)^e_wt_cl_megx *
      (1 + e_page_cl_megx * (pma_days - 280))
    vc_megx <- exp(lvc_megx + etalvc_megx) * (WT / 3.5)^e_wt_vc_megx

    kel <- cl / vc
    kel_megx <- cl_megx / vc_megx

    # The analysis ran in umol (Methods 2.4), so the formation flux is
    # molar: every umol of lidocaine eliminated forms F umol of MEGX, and
    # F is absorbed into the apparent MEGX parameters. Here doses are in
    # mg lidocaine base and MEGX amounts in mg MEGX, so the flux carries
    # the molar-mass ratio MEGX (C12H18N2O, 206.29 g/mol) / lidocaine
    # (C14H22N2O, 234.34 g/mol). Molar masses are chemistry constants, not
    # printed in the paper.
    mw_ratio_megx <- 206.29 / 234.34

    d / dt(central) <- -kel * central
    d / dt(central_megx) <- kel * central * mw_ratio_megx - kel_megx * central_megx

    Cc <- central / vc
    Cc_megx <- central_megx / vc_megx

    Cc ~ add(addSd) + prop(propSd)
    Cc_megx ~ add(addSd_megx) + prop(propSd_megx)
  })
}
