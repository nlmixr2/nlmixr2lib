Green_2003_enoxaparin <- function() {
  description <- "Two-compartment first-order-input population PK model for subcutaneous enoxaparin in adults treated at the Royal Brisbane Hospital for acute coronary syndrome, deep vein thrombosis, pulmonary embolism, or DVT prophylaxis (Green & Duffull 2003). Anti-Xa activity is the observation; lean body weight (LBW; James 1976 formula) is the size descriptor on clearance and total body weight is the size descriptor on the central volume."
  reference <- "Green B, Duffull SB. Development of a dosing strategy for enoxaparin in obese patients. Br J Clin Pharmacol. 2003 Jul;56(1):96-103. doi:10.1046/j.1365-2125.2003.01849.x. PMID:12848781."
  vignette <- "Green_2003_enoxaparin"
  units <- list(time = "hour", dosing = "IU", concentration = "IU/L")

  covariateData <- list(
    WT = list(
      description        = "Total body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Linear size descriptor on the central volume of distribution Vc with reference 70 kg (Green 2003 Table 3 footnote: V2 reported in units of L 70 kg^-1 (WT)). Cohort weight range 41-160 kg (Table 1).",
      source_name        = "WT"
    ),
    LBM = list(
      description        = "Lean body weight (canonical column LBM; source paper uses LBW)",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Linear size descriptor on clearance CL with reference 70 kg",
        "(Green 2003 Table 3 footnote: CL reported in units of L h^-1 70 kg^-1 (LBW)).",
        "LBW computed via the James (1976) formula reproduced in the paper's Methods:",
        "LBW_male (kg)   = 1.1  * WT - 120 * (WT/HT)^2;",
        "LBW_female (kg) = 1.07 * WT - 148 * (WT/HT)^2",
        "(WT in kg, HT in cm). Stored under canonical LBM (lean body mass) per",
        "inst/references/covariate-columns.md; LBW and LBM refer to the same quantity.",
        "Sex is implicit in the LBM derivation; the paper notes that LBW alone outperformed",
        "sex-only and (sex + size) covariate models on CL."
      ),
      source_name        = "LBW"
    )
  )

  covariatesDataExcluded <- list(
    HT = list(
      description = "Body height",
      units       = "cm",
      type        = "continuous",
      notes       = "Required to derive LBM (LBW) via the James (1976) formula reproduced in Methods, but not directly referenced inside model(). Documented here so simulators have the dependency visible."
    ),
    SEXF = list(
      description       = "Sex (1 = female, 0 = male)",
      units             = "(binary)",
      type              = "binary",
      reference_category = "Male (SEXF = 0)",
      notes             = "Required to choose the male vs female James (1976) LBW formula. Not referenced inside model() directly; documented here so simulators have the dependency visible. Sex was screened during covariate analysis but LBW alone outperformed sex-only and (sex + size) models on CL.",
      source_name        = "SEX"
    ),
    AGE = list(
      description = "Subject age",
      units       = "years",
      type        = "continuous",
      notes       = "Age was a covariate in the post-hoc binomial bruising logistic-regression model (Results, 'Logistic regression model'), not in the structural popPK model. Documented here so simulators reproducing the bruising analysis have the column declared."
    ),
    BMI = list(
      description = "Body mass index",
      units       = "kg/m^2",
      type        = "continuous",
      notes       = "Screened as a size descriptor during covariate analysis but not retained in the final popPK model (Methods, 'Population analysis')."
    ),
    IBW = list(
      description = "Ideal body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened (Devine variant: IBW = 45.4 + 0.89 * (HT - 152.4) + 4.5 if male) but not retained (Methods, 'Population analysis')."
    ),
    ABW = list(
      description = "Adjusted body weight",
      units       = "kg",
      type        = "continuous",
      notes       = "Screened (ABW = IBW + 0.4 * (WT - IBW)) but not retained (Methods, 'Population analysis')."
    ),
    BSA = list(
      description = "Body surface area",
      units       = "m^2",
      type        = "continuous",
      notes       = "Screened (BSA = sqrt(HT * WT / 3600)) but not retained (Methods, 'Population analysis')."
    ),
    CRCL = list(
      description = "Creatinine clearance",
      units       = "mL/min",
      type        = "continuous",
      notes       = "Reported in Table 1 demographics. Enrollment required CrCl >= 72 mL/min (Cockcroft-Gault using LBW). Did not improve bruising-regression fit alone or with age (Results, 'Logistic regression model'); not retained in the popPK model either."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 96L,
    n_studies      = 1L,
    age_range      = "mean 56.3 years (SD 16.9); not reported as min-max",
    weight_range   = "41-160 kg (mean 85.0, SD 20.5)",
    height_range   = "mean 173 cm (SD 9.62); not reported as min-max",
    bmi_range      = "15.0-44.9 kg/m^2 (mean 28.1, SD 6.27); enriched for obesity by design (1/3 BMI < 25, 1/3 BMI 25-29.99, 1/3 BMI > 30)",
    sex_female_pct = 26L,
    race_ethnicity = NULL,
    disease_state  = "Mixed indications: acute coronary syndrome (33%), DVT (14%), pulmonary embolism (5%), DVT prophylaxis (48%)",
    dose_range     = "Subcutaneous enoxaparin per normal clinical care; doses ranged from 40 mg once daily (prophylaxis) to weight-based therapeutic doses (~100 IU/kg = 1 mg/kg twice daily for ACS/DVT/PE)",
    regions        = "Single-centre Royal Brisbane Hospital, Australia",
    notes          = paste(
      "Prospective single-centre study at the Royal Brisbane Hospital with approximately 3 anti-Xa samples per patient.",
      "Enrollment required normal hepatic enzymes (<= 2x ULN), normal bilirubin / albumin,",
      "and estimated GFR >= 72 mL/min by Cockcroft-Gault using LBW.",
      "Table 1 (Demographic data) is the source for the demographic distribution.",
      "Race / ethnicity was not tabulated."
    )
  )

  ini({
    # Structural parameters - typical values for the paper's reference patient
    # (LBW = 70 kg, WT = 70 kg). Units: CL in L/h, Vc / Vp in L, Q in L/h, Ka in 1/h.
    # Source: Green 2003 Table 3 ("Final parameter estimates for covariate model").
    # The paper fits CL and Vc as P_pop = P_70 * (size / 70) per the
    # "L h^-1 70 kg^-1 (LBW)" and "L 70 kg^-1 (WT)" unit conventions in Table 3.
    lcl <- log(1.03);  label("Clearance for a 70 kg LBM (lean body weight) patient (CL, L/h)")               # Green 2003 Table 3: CL = 1.03 L/h per 70 kg LBW
    lvc <- log(3.67);  label("Central volume of distribution for a 70 kg WT patient (V2, L)")                # Green 2003 Table 3: V2 = 3.67 L per 70 kg WT
    lka <- log(0.195); label("First-order absorption rate constant from the SC depot (Ka, 1/h)")             # Green 2003 Table 3: Ka = 0.195 1/h
    lvp <- log(13.1);  label("Peripheral volume of distribution (V3, L)")                                    # Green 2003 Table 3: V3 = 13.1 L
    lq  <- log(0.363); label("Inter-compartmental clearance between central and peripheral (Q, L/h)")        # Green 2003 Table 3: Q = 0.363 L/h

    # Linear size descriptors enter as ratios with reference 70 kg, so no
    # explicit allometric exponent is estimated (paper reports CL and V2 in
    # units of L h^-1 70 kg^-1 and L 70 kg^-1, i.e. the exponent is fixed to 1
    # by construction; no per-parameter uncertainty is reported for the
    # exponent). Wrap in fixed() to make the structural assumption explicit.
    e_lbw_cl <- fixed(1); label("Linear exponent of LBM on CL (unitless; fixed by Table 3 unit convention)") # Green 2003 Table 3 unit convention (CL per 70 kg LBW)
    e_wt_vc  <- fixed(1); label("Linear exponent of WT on Vc (unitless; fixed by Table 3 unit convention)")  # Green 2003 Table 3 unit convention (V2 per 70 kg WT)

    # Inter-individual variability. Green 2003 reports BSV as %CV with a
    # log-normal model (P_j = P_pop * exp(eta_j)); for log-normal, the variance
    # on the eta scale is omega^2 = log((CV/100)^2 + 1).
    #   CV(CL) = 35.6% -> omega^2 = log(1 + 0.356^2) = 0.11933
    #   CV(V2) = 58.0% -> omega^2 = log(1 + 0.58^2)  = 0.28991
    # IIV was not reported on Ka, Q, or V3 in the final covariate model.
    etalcl ~ 0.11933 # Green 2003 Table 3: wCL = 35.6% CV
    etalvc ~ 0.28991 # Green 2003 Table 3: wV2 = 58.0% CV

    # Additive residual unexplained variability on anti-Xa concentration.
    # Green 2003 reports sigma^2 = 6430 (IU/L)^2 -> sigma = sqrt(6430) = 80.19 IU/L.
    # The paper modeled additive (not proportional) residual variance.
    addSd <- 80.19; label("Additive residual error SD on anti-Xa concentration (IU/L)") # Green 2003 Table 3: sigma^2 = 6430 (IU/L)^2 -> SD = sqrt(6430) = 80.19 IU/L
  })

  model({
    # Size scaling: linear in (LBM/70) for CL, linear in (WT/70) for Vc.
    # The exponents are fixed to 1 by the paper's unit convention (Table 3),
    # so they appear in fixed() form in ini() and are kept symbolic here for
    # source-trace transparency.
    cl <- exp(lcl + etalcl) * (LBM / 70)^e_lbw_cl
    vc <- exp(lvc + etalvc) * (WT  / 70)^e_wt_vc
    ka <- exp(lka)
    vp <- exp(lvp)
    q  <- exp(lq)

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Two-compartment first-order-input SC model. Subcutaneous enoxaparin enters
    # the depot and absorbs into the central compartment with rate constant Ka.
    # Bioavailability F is not estimated and is taken as 1 (the paper reports
    # parameters as apparent clearance / volume relative to administered IU,
    # so F is absorbed into CL/F and Vc/F; modelling F = 1 is consistent with
    # how the parameters were estimated).
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                                k12 * central - k21 * peripheral1

    # Anti-Xa activity (IU/L) = central amount (IU) / Vc (L).
    Cc <- central / vc
    Cc ~ add(addSd)
  })
}
