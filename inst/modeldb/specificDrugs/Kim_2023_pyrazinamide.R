Kim_2023_pyrazinamide <- function() {
  description <- "One-compartment population pharmacokinetic model with first-order absorption and first-order elimination for oral pyrazinamide in Korean adults with drug-susceptible tuberculosis (Kim 2023); lean body mass is an allometric covariate on CL/F and V/F (fixed exponents 0.75 and 1) and geriatric diabetes mellitus (age >= 70 years with diabetes mellitus) increases CL/F by 32%"
  reference <- "Kim R, Jayanti RP, Lee H, Kim H-K, Kang J, Park I-N, Kim J, Oh JY, Kim HW, Lee H, Ghim J-L, Ahn S, Long NP, Cho Y-S, Shin J-G; on behalf of the cPMTb. Development of a population pharmacokinetic model of pyrazinamide to guide personalized therapy: impacts of geriatric and diabetes mellitus on clearance. Front Pharmacol. 2023;14:1116226. doi:10.3389/fphar.2023.1116226"
  vignette <- "Kim_2023_pyrazinamide"
  units <- list(time = "hour", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    LBM = list(
      description        = "Lean body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric size descriptor on both CL/F and V/F with fixed exponents 0.75 and 1 respectively, normalised to a 48 kg reference subject (Supplementary File S1 $PK: 'TVCL = THETA (1) * (LBW/48)**(THETA(6))' and 'TVV = THETA (2) * (LBW/48)**(THETA(7))'). The 48 kg reference is the rounded median lean body weight of the 488-patient training set (48.1 kg, range 23.1-63.79 kg; Kim 2023 Table 1). Lean body weight outperformed total body weight as the size descriptor: adding allometric scaling on lean body weight dropped the objective function value by 112.3 points versus 88.5 points for total body weight (Kim 2023 Results, 'Development of a population PK model').",
      source_name        = "LBW"
    ),
    AGE_GE70 = list(
      description        = "Advanced-age indicator (1 = aged 70 years or older at baseline, 0 = younger than 70 years).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (younger than 70 years)",
      notes              = "Kim 2023 tested advanced age at three thresholds (>= 60, >= 65 and >= 70 years) both alone and crossed with comorbidities. Age alone -- whether continuous or as any of the three threshold indicators -- did not improve the objective function value; only the interaction with diabetes mellitus did, and the >= 70 year threshold gave the largest drop (151.3 points, p < 0.001). This column therefore enters the model ONLY through the product AGE_GE70 * DIS_DIAB, never as a main effect. The underlying continuous AGE is listed in covariatesDataExcluded.",
      source_name        = "OLD"
    ),
    DIS_DIAB = list(
      description        = "Diabetes-mellitus comorbidity indicator (1 = diabetes mellitus, 0 = no diabetes mellitus).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (no diabetes mellitus)",
      notes              = "Diabetes mellitus recorded as a comorbidity at study entry; 55 of the 488 training-set patients (11.27%) were diabetic (Kim 2023 Table 1). Type 1 versus type 2 is not distinguished in the source paper. Like AGE_GE70, diabetes mellitus alone did not improve the objective function value and enters the model ONLY through the product AGE_GE70 * DIS_DIAB. Of the 488 training patients, 23 were both >= 70 years old and diabetic (Supplementary File S2 subgroup table). Kim 2023 attributes the higher apparent clearance in this subgroup to diabetes-induced elevation of xanthine oxidase, the enzyme that converts pyrazinamide to 5-hydroxypyrazinoic acid (Discussion).",
      source_name        = "DM"
    )
  )

  # Covariates that Kim 2023 screened but did not retain in the final model.
  # Documentation only -- these are not referenced in model().
  covariatesDataExcluded <- list(
    WT = list(
      description = "Total body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened as an allometric size descriptor but rejected in favour of lean body weight, which produced the larger objective-function drop (112.3 vs 88.5 points; Kim 2023 Results). Training-set median 60.8 kg (range 28.8-95.3 kg)."
    ),
    AGE = list(
      description = "Subject age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened as a continuous covariate and as threshold indicators at >= 60, >= 65 and >= 70 years; none improved the objective function value on its own. Retained only as the >= 70 year threshold crossed with diabetes mellitus (see AGE_GE70 in covariateData). Training-set median 54.5 years (range 19-96 years)."
    ),
    HT = list(
      description = "Body height at baseline",
      units       = "cm",
      type        = "continuous",
      notes       = "Screened; did not improve the objective function value (Kim 2023 Results). Also an input to the lean-body-weight derivation."
    ),
    SEXF = list(
      description = "Biological sex indicator (1 = female, 0 = male).",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened; did not improve the objective function value. Training set was 33.4% female (Kim 2023 Table 1)."
    ),
    FED = list(
      description = "Fed-versus-fasted status at the dose preceding the sample (1 = fed, 0 = fasted).",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened as 'fasting or food intake status'; did not improve the objective function value. 75.6% of training-set records were fasted (Kim 2023 Table 1)."
    ),
    ALB = list(
      description = "Serum albumin",
      units       = "g/dL (as reported by Kim 2023 Table 1; the canonical SI unit is g/L)",
      type        = "continuous",
      notes       = "Screened; did not improve the objective function value. Training-set median 4.2 g/dL (range 2.0-14.1 g/dL)."
    ),
    TBILI = list(
      description = "Total bilirubin",
      units       = "not reported by Kim 2023 (values are not tabulated)",
      type        = "continuous",
      notes       = "Screened as a continuous covariate; did not improve the objective function value. Present as the NONMEM input column TBIL in Supplementary File S1 but never used in $PK."
    ),
    AST = list(
      description = "Aspartate aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened; did not improve the objective function value. Training-set median 27 U/L (range 4.2-191 U/L)."
    ),
    ALT = list(
      description = "Alanine aminotransferase",
      units       = "U/L",
      type        = "continuous",
      notes       = "Screened; did not improve the objective function value. Training-set median 19 U/L (range 0.8-233 U/L)."
    ),
    CRCL = list(
      description = "Estimated glomerular filtration rate, used by Kim 2023 to define the renal-disease stratum.",
      units       = "mL/min/1.73 m^2",
      type        = "continuous",
      notes       = "Kim 2023 screened 'renal disease', defined as eGFR <= 60 mL/min/1.73 m^2, as a binary covariate and also crossed with advanced age; neither improved the objective function value. 19 of 488 training-set patients (3.9%) were below the threshold (Kim 2023 Table 1). Liver disease (a chart diagnosis, 2 of the 110 patients >= 70 years old) was screened the same way and likewise not retained; it has no canonical covariate column because it was not retained by any model in the library."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 488L,
    n_subjects_enrolled = 613L,
    n_subjects_test     = 125L,
    n_studies      = 1L,
    n_sites        = 18L,
    age_range      = "19-96 years",
    age_median     = "54.5 years",
    weight_range   = "28.8-95.3 kg",
    weight_median  = "60.8 kg",
    lbm_range      = "23.1-63.79 kg",
    lbm_median     = "48.1 kg",
    sex_female_pct = 33.4,
    race_ethnicity = c(Korean = 100),
    disease_state  = "Adults (> 18 years) with bacteriologically confirmed drug-susceptible tuberculosis receiving a pyrazinamide-containing regimen for at least 2 weeks. Comorbidities in the 488-patient training set: 55 (11.27%) diabetes mellitus, 19 (3.9%) renal disease (eGFR <= 60 mL/min/1.73 m^2). 110 of the 613 enrolled patients were older than 70 years, of whom 23 had diabetes mellitus, 11 renal disease and 2 liver disease.",
    dose_range     = "Oral pyrazinamide once daily at 20-30 mg/kg, rounded to the nearest tablet size as prescribed: 500 mg (0.6%), 1,000 mg (11.3%), 1,200 mg (4.1%), 1,250 mg (4.1%), 1,500 mg (67.2%), 1,600 mg (10.0%) and 2,000 mg (1.4%) in the training set. Sampling was at steady state; non-adherent and non-steady-state patients were excluded.",
    regions        = "Republic of Korea (18 hospitals, cPMTb multicentre prospective tuberculosis cohort).",
    renal_function = "Predominantly normal: 96.1% of the training set had eGFR > 60 mL/min/1.73 m^2.",
    co_medication  = "Companion anti-tuberculosis regimens, chiefly RHZE (85.8% of the training set), plus RHZL, RZEM, RZEL and HZEL. Interactions among the anti-tuberculosis drugs were not modelled (Kim 2023 Discussion, limitation 3).",
    notes          = "Data were split at random into a 488-patient training set (model development) and a 125-patient test set (external validation) in a 4:1 ratio. Blood samples were drawn at random times 0-24 h after the last observed dose: typically one sample per outpatient, and at least two samples (from pre-dose, 1, 2 and 5 h post-dose) per inpatient. Pyrazinamide was quantified by validated HPLC-ESI-MS/MS over 2.0-80.0 ug/mL (LLOQ 2.0 ug/mL); concentrations below the LLOQ were imputed as 1 ug/mL per Beal's M5 method. Baseline demographics from Kim 2023 Table 1 (training-data column); final parameter estimates from Kim 2023 Table 2; model structure and reference values from Supplementary File S1 (NONMEM control stream); geriatric-diabetes subgroup sizes from Supplementary File S2."
  )

  ini({
    # Structural parameters. Kim 2023 fits pyrazinamide with a
    # one-compartment model with first-order absorption and first-order
    # elimination (NONMEM ADVAN2 TRANS2; Supplementary File S1 $SUBR).
    # Typical values reference a 48 kg lean-body-weight patient who is not
    # simultaneously >= 70 years old and diabetic. Interindividual
    # variability was log-normal: CL = TVCL * EXP(ETA(1)) etc.
    # (Supplementary File S1 $PK).
    # Absorption models with a lag time, sequential zero- then first-order
    # input, and transit compartments were all evaluated and none improved
    # model performance (Kim 2023 Results), so plain first-order absorption
    # is the final structure.
    # NOTE ON THETA NUMBERING: Kim 2023 Table 2 and Supplementary File S1 use
    # DIFFERENT theta indices for the same quantities. Table 2 numbers them
    # theta1 = CL/F, theta2 = geriatric-DM effect, theta3 = Vd/F, theta4 = Ka.
    # The S1 control stream numbers them THETA(1) = CL, THETA(2) = V,
    # THETA(3) = KA, THETA(4) = additive SD, THETA(5) = proportional (fixed 0),
    # THETA(6) = allometric exponent on CL, THETA(7) = allometric exponent on V,
    # THETA(8) = geriatric-DM effect. Comments below name the source explicitly
    # to keep the two schemes apart.
    lka <- log(1.49); label("Absorption rate constant (1/h)")                                              # Kim 2023 Table 2 'Ka (1/h) = theta4'; = S1 THETA(3)
    lcl <- log(4.49); label("Apparent oral clearance for a 48 kg lean-body-weight non-geriatric-diabetic patient (L/h)")  # Kim 2023 Table 2 'CL/F (L/h); = theta1'; = S1 THETA(1)
    lvc <- log(44.2); label("Apparent volume of distribution for a 48 kg lean-body-weight patient (L)")     # Kim 2023 Table 2 'Vd/F (L) = theta3'; = S1 THETA(2)

    # Allometric exponents on lean body weight. Both are FIXED to the
    # canonical Anderson and Holford (2008) values, not estimated: "Allometric
    # scaling of either body weight or lean body mass was applied to apparent
    # clearance (CL/F) and central volume of distribution (Vd/F) using fixed
    # exponents of 0.75 and 1, respectively" (Kim 2023 Methods). Confirmed by
    # the FIX flags in Supplementary File S1 $THETA: '(0.75 FIX)' = THETA(6)
    # and '(1 FIX)' = THETA(7). Neither is given an %RSE in Table 2.
    e_lbm_cl <- fixed(0.75); label("Allometric exponent of lean body weight on CL/F (fixed, unitless)")     # Kim 2023 Methods; Suppl. S1 $THETA (0.75 FIX) = S1 THETA(6)
    e_lbm_vc <- fixed(1.0);  label("Allometric exponent of lean body weight on V/F (fixed, unitless)")      # Kim 2023 Methods; Suppl. S1 $THETA (1 FIX) = S1 THETA(7)

    # Geriatric-diabetes effect on CL/F. Kim 2023 encodes the categorical
    # effect as a fractional increase from 1: Table 2 gives the covariate row
    # as "CL/F; elder patient (>= 70 years old) with DM = theta1 x (1 + theta2)",
    # matching Supplementary File S1 $PK
    # 'TVCL = THETA (1) * (LBW/48)**(THETA(6)) * (1 + OLDDM1 * THETA(8))'.
    # theta2 = 0.32 therefore raises CL/F by 32%, giving
    # 4.49 * 1.32 = 5.93 L/h, which matches the 5.9 L/h quoted in the Results
    # text for geriatric diabetic patients.
    e_age_ge70_diab_cl <- 0.32; label("Fractional increase in CL/F for patients >= 70 years old with diabetes mellitus (unitless)")  # Kim 2023 Table 2 theta2; = S1 THETA(8)

    # Interindividual variability. Table 2 reports the IIV rows as
    # "omega^2; <parameter> (%)", i.e. the variance multiplied by 100. The
    # scale is pinned by the fixed Vd/F row: Table 2 prints "3 (FIX)" where
    # Supplementary File S1 $OMEGA carries the literal '(0.03 FIX)', so
    # table value = 100 * omega^2.
    #   CL/F: 1   -> omega^2 = 0.01 (%RSE 60, shrinkage 18.0%)
    #   Vd/F: 3   -> omega^2 = 0.03 FIXED (shrinkage 28.4%)
    #   Ka:   100 -> omega^2 = 1.00 (shrinkage 32.7%)
    # The IIV on Vd/F was fixed to stabilise the model and achieve successful
    # minimisation, at the value obtained from a prior run before fixing
    # (Kim 2023 Results); Supplementary File S1 $OMEGA confirms the FIX flag.
    # The large Ka variance and its 32.7% shrinkage reflect sparse sampling in
    # the absorption phase (Kim 2023 Results).
    etalcl ~ 0.01        # Kim 2023 Table 2 IIV CL/F, omega^2 = 1%   -> 0.01
    etalvc ~ fixed(0.03) # Kim 2023 Table 2 IIV Vd/F, omega^2 = 3% (FIX); Suppl. S1 $OMEGA (0.03 FIX)
    etalka ~ 1           # Kim 2023 Table 2 IIV Ka,   omega^2 = 100% -> 1.00

    # Residual error. The final model is additive only: Kim 2023 Table 2 lists
    # a single "Additive" residual row (3.41 ug/mL) and the Results state that
    # additive residual error described the data adequately. Supplementary
    # File S1 codes a combined scaffold,
    # 'W = SQRT (THETA (4)**2 + THETA(5)**2 * IPRED**2)' with Y = IPRED + W*EPS(1)
    # and $SIGMA 1 FIX, but the proportional coefficient THETA(5) is declared
    # '(0 FIX)' -- fixed to zero -- so only the additive term is active. The
    # proportional term is retained here as fixed(0) to preserve that
    # provenance.
    addSd  <- 3.41;     label("Additive residual standard deviation (ug/mL)")                      # Kim 2023 Table 2 'Residual variability / Additive' = 3.41; = S1 THETA(4)
    propSd <- fixed(0); label("Proportional residual standard deviation (fraction; fixed to zero)")  # Suppl. S1 $THETA (0 FIX) = S1 THETA(5); not reported in Table 2
  })

  model({
    # Composite geriatric-diabetes indicator. Supplementary File S1 $PK
    # derives it as a single data column:
    #   OLDDM1 = 0
    #   IF (OLDDM.EQ.1) OLDDM1 = 1
    # where OLDDM flags patients who are simultaneously >= 70 years old and
    # diabetic. It is reconstructed here as the product of the two canonical
    # covariate columns so that neither has to be pre-crossed in the data.
    age_ge70_diab <- AGE_GE70 * DIS_DIAB

    # Individual PK parameters, reproducing Supplementary File S1 $PK with
    # lean body weight normalised to 48 kg:
    #   CL/F = 4.49 * (LBW/48)^0.75 * (1 + 0.32 * OLDDM1)
    #   V/F  = 44.2 * (LBW/48)^1
    #   Ka   = 1.49
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl) * (LBM / 48)^e_lbm_cl * (1 + e_age_ge70_diab_cl * age_ge70_diab)
    vc <- exp(lvc + etalvc) * (LBM / 48)^e_lbm_vc

    kel <- cl / vc

    # ODEs: one-compartment model with first-order absorption, equivalent to
    # NONMEM ADVAN2 TRANS2 with S2 = V and K = CL/V.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # Observation. CL/F and V/F are apparent (oral) parameters, so
    # bioavailability is absorbed into them and no f(depot) term is needed.
    # Dose in mg divided by volume in L gives mg/L = ug/mL, matching the
    # assay units of the source paper.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
