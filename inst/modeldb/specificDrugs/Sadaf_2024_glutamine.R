Sadaf_2024_glutamine <- function() {
  description <- paste(
    "One-compartment population PK model with first-order absorption for oral",
    "L-glutamine (Endari) in children and adults with sickle cell disease and",
    "in healthy adult volunteers (Sadaf 2024). The model is fitted to",
    "baseline-subtracted plasma L-glutamine, i.e. the exogenous increment",
    "above each participant's own endogenous pre-dose concentration, so the",
    "predicted Cc is zero before the first dose. Apparent oral clearance and",
    "apparent central volume are allometrically scaled to body weight with the",
    "theoretical exponents 0.75 and 1 (70 kg reference). Apparent clearance",
    "falls as the pre-dose endogenous L-glutamine baseline rises (power -0.96,",
    "reference 683 umol/L), interpreted as competition between endogenous and",
    "exogenous glutamine for elimination. Apparent volume rises with the",
    "administered dose per kg (power 0.27, reference 0.1 g/kg), which is how",
    "the paper captures the capacity-limited, less-than-dose-proportional rise",
    "in exposure seen at 0.3 and 0.6 g/kg after a Michaelis-Menten base model",
    "proved unstable. Inter-occasion variability is carried on clearance across",
    "the four study visits; inter-individual variability on clearance and on",
    "the absorption rate constant, and the additive residual error, were all",
    "held at zero in the published final model.",
    sep = " "
  )
  reference <- paste(
    "Sadaf A, Dong M, Pfeiffer A, Latham T, Kalfa T, Vinks AA, Ware RE,",
    "Quinn CT. A Population Pharmacokinetic Analysis of L-Glutamine Exposure",
    "in Patients with Sickle Cell Disease: Evaluation of Dose and Food",
    "Effects. Clin Pharmacokinet. 2024;63(3):357-365.",
    "doi:10.1007/s40262-024-01349-4.",
    sep = " "
  )
  vignette <- "Sadaf_2024_glutamine"
  units <- list(time = "hour", dosing = "umol", concentration = "umol/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Sadaf 2024 Eq. 3 and Table 2. Allometric scaling to a 70 kg",
        "reference with the theoretical exponents held at 0.75 for CL/F and",
        "1 for V/F. The paper also estimated the exponents freely and",
        "obtained 0.73 (CL/F) and 0.78 (V/F), 'similar to the theoretical",
        "values', but the final model reported in Table 2 uses the fixed",
        "theoretical exponents. Cohort weights 30.3-119.1 kg (Table 1).",
        sep = " "
      ),
      source_name        = "WT"
    ),
    GLN_BL = list(
      description        = "Pre-dose (baseline) endogenous plasma L-glutamine concentration",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Sadaf 2024 Table 2: CLi = theta1 * (WT/70)^0.75 *",
        "(Glu_BSL / Glu_BSLstandard)^theta2 with the standard baseline",
        "Glu_BSLstandard = 683 umol/L given in the Table 2 footnote and",
        "theta2 = -0.96, i.e. higher endogenous baselines go with lower",
        "apparent clearance. Per-subject and time-fixed in the packaged",
        "model; the paper verified that baseline did not drift across the",
        "four visits (Results 3.2, Fig. 2). This same per-subject baseline is",
        "the quantity subtracted from every measured concentration before",
        "fitting (Methods 2.6.1), so it plays a dual role: it defines the",
        "zero of the modelled concentration scale AND it enters the CL/F",
        "covariate equation. Group means 571 (pediatric SCD), 625 (adult",
        "SCD), and 531 (healthy volunteers) umol/L; overall range",
        "364-758 umol/L (Table 1).",
        sep = " "
      ),
      source_name        = "Glu_BSL"
    ),
    DOSE_GLN_GKG = list(
      description        = "Administered L-glutamine dose per kg body weight for the current dosing occasion",
      units              = "g/kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = paste(
        "Sadaf 2024 Table 2: Vi = theta3 * (WT/70) * (DOSE per kg / 0.1)^theta4",
        "with theta4 = 0.27 and the 0.1 g/kg lowest dose level as the",
        "reference. Takes the values 0.1, 0.3, and 0.6 g/kg in the source",
        "trial. Per-dosing-occasion, not per-subject: every participant",
        "received all three levels in ascending order over three weeks. This",
        "covariate is the paper's surrogate for capacity-limited absorption -",
        "a Michaelis-Menten disposition base model was tested first and was",
        "unstable (Results 3.3), so the less-than-proportional exposure",
        "increase was absorbed into an apparent volume that grows with dose.",
        "Because it inflates V/F only, it lowers peak concentration without",
        "changing AUC, which is the paper's stated reading ('the higher the",
        "dose, the lower the relative bioavailability').",
        sep = " "
      ),
      source_name        = "DOSE per kg"
    ),
    OCC = list(
      description        = "Study-visit occasion index for inter-occasion variability on CL/F",
      units              = "(count)",
      type               = "categorical",
      reference_category = NULL,
      notes              = paste(
        "Sadaf 2024 Methods 2.6.1: 'For the full dataset, four occasions were",
        "defined', corresponding to study visits 1-4. Decomposed inside",
        "model() into binary indicators oc1..oc4 that multiplex four",
        "per-occasion etas on log-CL/F, all sharing the single published IOV",
        "variance of 0.081 (Table 2) in the NONMEM $OMEGA BLOCK(1) SAME",
        "idiom. Adding IOV on clearance dropped the objective function by",
        "34.5 units (Results 3.3).",
        sep = " "
      ),
      source_name        = "OCC"
    )
  )

  # Covariates that Sadaf 2024 screened by stepwise selection but did NOT
  # retain in the final model (Results 3.3: "No effects of age, sex, Hb
  # genotype, serum creatinine levels, or baseline l-glutamate concentrations
  # on l-glutamine pharmacokinetics were observed"). No point estimates are
  # published for any of them, so they are documented rather than encoded.
  covariatesDataExcluded <- list(
    AGE = list(
      description = "Age",
      units       = "years",
      type        = "continuous",
      notes       = "Screened (6.8-43.8 years, Table 1); not retained. Sadaf 2024 Results 3.3 and Discussion report no age effect once weight was scaled allometrically."
    ),
    HT = list(
      description = "Height",
      units       = "cm",
      type        = "continuous",
      notes       = "Screened as a body-size descriptor (Methods 2.6.2); not retained - allometric body weight was the size descriptor selected."
    ),
    BSA = list(
      description = "Body surface area",
      units       = "m^2",
      type        = "continuous",
      notes       = "Screened as a body-size descriptor (Methods 2.6.2); not retained - allometric body weight was the size descriptor selected."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units       = "(binary)",
      type        = "binary",
      notes       = "Screened (8 of 12 participants female, Table 1); not retained (Results 3.3)."
    ),
    CREAT = list(
      description = "Serum creatinine",
      units       = "mg/dL",
      type        = "continuous",
      notes       = "Screened (Methods 2.6.2); not retained (Results 3.3). Cohort creatinine was unchanged by treatment (0.63 vs 0.63 mg/dL, p = 0.887, Results 3.4)."
    ),
    FED = list(
      description = "Meal-intake status of the dosing occasion",
      units       = "(0 / 0.5 / 1 / 1.5)",
      type        = "continuous",
      notes       = paste(
        "Screened on CL/F and not retained - this is the paper's headline",
        "negative finding ('Food intake did not alter glutamine clearance,",
        "thus l-glutamine can be taken with or without food'). Note that",
        "Sadaf 2024 Methods 2.6.2 did NOT code food as a plain binary: the",
        "meal covariate took 0 for fasting, 0.5 for a snack, 1 for a full",
        "meal, and 1.5 when a participant had both a snack and a full meal",
        "between two consecutive measurements. It is recorded here under the",
        "canonical FED name for discoverability, but a future extraction that",
        "retains this graded coding should not assume FED is binary.",
        sep = " "
      )
    )
  )

  compartmentData <- list(
    depot   = list(analyte = "L-glutamine", units = "umol", specimen = "administration site", verified = TRUE),
    central = list(analyte = "L-glutamine", units = "umol", specimen = "plasma", verified = TRUE)
  )

  population <- list(
    species        = "human",
    n_subjects     = 12,
    n_studies      = 1,
    age_range      = "6.8-43.8 years",
    weight_range   = "30.3-119.1 kg",
    sex_female_pct = 66.7,
    disease_state  = paste(
      "Sickle cell disease (N = 8: 4 homozygous sickle cell anemia HbSS,",
      "4 sickle-hemoglobin C disease HbSC) plus healthy adult volunteers",
      "with a normal hemoglobin profile (N = 4, HbAA). Concomitant",
      "hydroxyurea was allowed if the dose had been stable for 3 months",
      "and stayed unchanged during the trial; the paper is internally",
      "inconsistent about how many took it (Table 1 gives 1 of 4 pediatric",
      "and 2 of 4 adult, i.e. 3 of 8, whereas Results 3.1 states 'four yes,",
      "four no'). Hydroxyurea was not among the covariates screened",
      "(Methods 2.6.2), so the discrepancy does not affect the model.",
      sep = " "
    ),
    dose_range     = "0.1 g/kg twice daily, then 0.3 g/kg twice daily, then 0.6 g/kg once daily, oral, one dose level per week for three weeks",
    regions        = "United States (single center: Cincinnati Children's Hospital Medical Center)",
    notes          = paste(
      "Open-label, dose-ascending phase IV trial (NCT04684381), January to",
      "June 2021. Baseline demographics are Sadaf 2024 Table 1, reported by",
      "the three allocation groups (pediatric SCD N = 4, adult SCD N = 4,",
      "healthy volunteers N = 4) rather than pooled; the ranges above span",
      "all three groups. 400 plasma L-glutamine concentrations were",
      "available. 15 participants were enrolled and 12 completed at least",
      "the third visit; the 3 withdrawn participants were excluded from the",
      "analysis. Sampling was pre-dose and at 30, 60, 120, 180, and 240 min",
      "after the morning dose at visits 1-4, and after the afternoon dose at",
      "visits 1-2. Doses were given fasting (morning, visits 1-3) or after a",
      "standardised meal (lunch at visits 1-2, breakfast at visit 4). Hb",
      "genotype was screened as a covariate and not retained, so it is not a",
      "model input; the packaged model applies equally to the SCD and",
      "healthy-volunteer cohorts, which is the paper's finding.",
      sep = " "
    )
  )

  ini({
    # Structural parameters. Reference subject: 70 kg, pre-dose L-glutamine
    # baseline 683 umol/L, dose level 0.1 g/kg (Sadaf 2024 Table 2).
    lka <- log(0.91) ; label("Absorption rate constant (1/h)")                     # Table 2: KA = 0.91 1/h (RSE 6.4%; bootstrap 95% CI 0.78-0.99)
    lcl <- log(78.5) ; label("Apparent oral clearance CL/F (L/h per 70 kg)")       # Table 2: theta1 = 78.5 L/h/70kg (RSE 7.6%; bootstrap 95% CI 67.1-91.2)
    lvc <- log(52.3) ; label("Apparent central volume V/F (L per 70 kg)")          # Table 2: theta3 = 52.3 L/70kg (RSE 19.5%; bootstrap 95% CI 30.9-69.9)

    # Allometric exponents. Sadaf 2024 Eq. 3 sets these at the theoretical
    # values rather than estimating them; Table 2 prints them as literal
    # exponents (0.75 / 1) with no estimate, RSE, or bootstrap interval,
    # unlike every estimated theta in the same table.
    e_wt_cl <- fixed(0.75) ; label("Allometric exponent of body weight on CL/F (unitless)")  # Eq. 3 and Table 2 CLi equation
    e_wt_vc <- fixed(1)    ; label("Allometric exponent of body weight on V/F (unitless)")   # Eq. 3 and Table 2 Vi equation

    # Covariate effects, both in the power form of Eq. 5.
    e_gln_bl_cl       <- -0.96 ; label("Power of pre-dose L-glutamine baseline on CL/F (unitless)")  # Table 2: theta2 = -0.96 (RSE 19.8%; bootstrap 95% CI -1.3 to -0.4)
    e_dose_gln_gkg_vc <-  0.27 ; label("Power of dose per kg on V/F (unitless)")                     # Table 2: theta4 = 0.27 (RSE 45.0%; bootstrap 95% CI 0.09-0.60)

    # Inter-individual variability. Table 2 reports variances directly (the
    # block is headed "Inter-patient variability (variance)"), so no CV%-to-
    # variance conversion is applied. IIV on CL and on Ka were held at zero
    # in the final model.
    etalcl ~ fixed(0)  # Table 2: IIV CL variance zero in the final model
    etalvc ~ 0.078     # Table 2: IIV V variance = 0.078 (RSE 19.7%; bootstrap 95% CI 0.004-0.226), i.e. 28% CV
    etalka ~ fixed(0)  # Table 2: IIV Ka variance zero in the final model

    # Inter-occasion variability on log-CL/F across the four study visits.
    # Table 2 reports one IOV variance, so all four per-occasion etas share
    # it (NONMEM $OMEGA BLOCK(1) SAME idiom); the first is left estimable and
    # the remaining three are fixed equal to it.
    etaiov_cl_1 ~ 0.081        # Table 2: IOV CL variance = 0.081 (RSE 36.4%; bootstrap 95% CI 0.028-0.144), i.e. 28% CV (occasion 1)
    etaiov_cl_2 ~ fixed(0.081) # IOV on log-CL/F, occasion 2; variance equal to occasion 1
    etaiov_cl_3 ~ fixed(0.081) # IOV on log-CL/F, occasion 3; variance equal to occasion 1
    etaiov_cl_4 ~ fixed(0.081) # IOV on log-CL/F, occasion 4; variance equal to occasion 1

    # Residual error. Eq. 2 as printed is purely proportional,
    # Y = Cpred * (1 + eps_prop), and Table 2 reports the additive part as
    # "Fixed to 0", so no additive term is carried. Table 2 gives the
    # proportional term as a VARIANCE of 0.49 (the block is headed "Residual
    # variability (variance)"), so the SD entered here is sqrt(0.49) = 0.7.
    propSd <- 0.7 ; label("Proportional residual error (fraction)")  # Table 2: Prop.Err. variance = 0.49 (RSE 10.9%; bootstrap 95% CI 0.39-0.59) -> SD = sqrt(0.49)
  })

  model({
    # 1. Occasion indicators - binary decomposition of the OCC column
    # (study visits 1-4, Sadaf 2024 Methods 2.6.1).
    oc1 <- (OCC == 1)
    oc2 <- (OCC == 2)
    oc3 <- (OCC == 3)
    oc4 <- (OCC == 4)

    # 2. Inter-occasion variability on log-CL/F: exactly one eta is active
    # per record, all four sharing the single published IOV variance.
    iov_cl <- oc1 * etaiov_cl_1 + oc2 * etaiov_cl_2 + oc3 * etaiov_cl_3 + oc4 * etaiov_cl_4

    # 3. Individual PK parameters (Sadaf 2024 Table 2 covariate equations).
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl + iov_cl) * (WT / 70)^e_wt_cl * (GLN_BL / 683)^e_gln_bl_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc * (DOSE_GLN_GKG / 0.1)^e_dose_gln_gkg_vc

    # 4. Micro-constants.
    kel <- cl / vc

    # 5. One-compartment disposition with first-order absorption.
    d/dt(depot)   <- -ka * depot
    d/dt(central) <-  ka * depot - kel * central

    # 6. Observation. Cc is the exogenous, baseline-subtracted plasma
    # L-glutamine concentration - the endogenous pre-dose level GLN_BL was
    # subtracted from every measurement before fitting (Methods 2.6.1). To
    # compare against a raw measured concentration, add GLN_BL back.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
