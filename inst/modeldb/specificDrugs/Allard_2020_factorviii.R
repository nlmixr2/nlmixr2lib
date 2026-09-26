Allard_2020_factorviii <- function() {
  description <- "Two-compartment population PK model for factor VIII activity (FVIII:C, one-stage clotting assay) pooled across eight standard and extended half-life FVIII concentrates in children, adolescents and adults with mainly severe haemophilia A (Allard 2020; weight and age on CL, weight on V1, Elocta (rFVIIIFc) indicator on CL)"
  reference <- "Allard Q, Djerada Z, Pouplard C, Repesse Y, Desprez D, Galinat H, Frotscher B, Berger C, Harroche A, Ryman A, Flaujac C, Chamouni P, Guillet B, Volot F, Szymezak J, Nguyen P, Cazaubon Y. Real Life Population Pharmacokinetics Modelling of Eight Factors VIII in Patients with Severe Haemophilia A: Is It Always Relevant to Switch to an Extended Half-Life? Pharmaceutics. 2020;12(4):380. doi:10.3390/pharmaceutics12040380"
  vignette <- "Allard_2020_factorviii"
  units <- list(time = "h", dosing = "IU", concentration = "IU/mL")

  compartmentData <- list(
    central = list(analyte = "factorviii", units = "IU", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "factorviii", units = "IU", specimen = "plasma", verified = TRUE)
  )

  covariateData <- list(
    WT = list(
      description = "Total body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Power scaling on CL (exponent fixed to 0.75) and on V1 (estimated exponent 0.827), centred on the study median of 64 kg (Allard 2020 Table 2 footnote).",
      source_name = "weight (TBW)"
    ),
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = "Power scaling on CL (exponent -0.214), centred on the study median of 30 years (Allard 2020 Table 2 footnote). Older patients have lower CL.",
      source_name = "Age"
    ),
    FORM_FVIII_FC = list(
      description = "Extended half-life product indicator: 1 = Elocta (efmoroctocog alfa, rFVIIIFc), 0 = any of the seven standard half-life FVIII concentrates in the dataset (Advate, Afstyla, Factane, Kogenate, Kovaltry, Novoeight, Refacto)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (standard half-life FVIII)",
      notes = "The paper's covariate 'EHL (elocta:1 versus others:0)'. Elocta was the only extended half-life product in the dataset, so the fitted effect exp(-0.394) = 0.674 on CL describes the Fc-fusion product specifically, not extended half-life FVIII as a class. Set per dose record for patients who switched products.",
      source_name = "EHL"
    )
  )

  covariatesDataExcluded <- list(
    VWF = list(
      description = "Plasma von Willebrand factor antigen",
      units = "IU/dL",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested but not retained because 130 patients had no vWF value (Allard 2020 Section 2.2.2); effect estimates are in Table S2, which is not reproduced here."
    ),
    BLOOD_GROUP_O = list(
      description = "ABO blood group O indicator",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (non-O blood group)",
      notes = "ABO blood group tested but not retained because 51 patients had no ABO value (Allard 2020 Section 2.2.2); Table S2."
    ),
    FFM = list(
      description = "Fat-free mass",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = "Tested but not retained because 61 patients had no height to derive FFM (Allard 2020 Section 2.2.2); Table S2."
    )
  )

  population <- list(
    species = "human",
    n_subjects = 258L,
    n_studies = 1L,
    age_range = "3-77 years (87 children/adolescents, 171 adults)",
    age_median = "30 years",
    weight_range = "15.1-130 kg",
    weight_median = "64 kg",
    sex_female_pct = NA_real_,
    race_ethnicity = "not reported",
    disease_state = "Haemophilia A without inhibitors: severe n = 244, moderate n = 11, minor n = 3",
    dose_range = "Intravenous FVIII concentrate, median 2750 IU (500-5000 IU) per PK dose",
    regions = "France (13 haemophilia treatment centres)",
    products = "Elocta n = 136, Advate n = 44, Kogenate n = 34, Refacto n = 18, Factane n = 8, Kovaltry n = 7, Novoeight n = 6, Afstyla n = 5 (PK profiles; 44 patients switched from a standard half-life product to Elocta and contribute one profile on each)",
    notes = "Retrospective real-life routine data collected 2012-2019 (Allard 2020 Section 2.1 and Table 1). 935 FVIII:C observations, 1-11 per patient (median 4), 10.6% below the LLOQ (0.004 or 0.01 IU/mL). Sex is not reported."
  )

  ini({
    # Structural parameters at the reference covariates (WT = 64 kg, AGE = 30 years,
    # standard half-life product). Allard 2020 Table 2.
    lcl <- log(204); label("Clearance CL at WT 64 kg, AGE 30 years, standard half-life product (mL/h)") # Allard 2020 Table 2: Cl = 204 mL/h (RSE 3.25%)
    lvc <- log(2640); label("Central volume V1 at WT 64 kg (mL)") # Allard 2020 Table 2: V1 = 2640 mL (RSE 2.04%)
    lq <- log(135); label("Intercompartmental clearance Q (mL/h)") # Allard 2020 Table 2: Q = 135 mL/h (RSE 20.5%)
    lvp <- log(339); label("Peripheral volume V2 (mL)") # Allard 2020 Table 2: V2 = 339 mL (RSE 10.7%)

    # Covariate effects. Signs of the age and Elocta coefficients are those of the
    # Table 2 footnote equation and bootstrap intervals (CLi = CLpop x (Age/30)^-0.214
    # x (weight/64)^0.75 x e^-0.394).
    e_wt_cl <- fixed(0.75); label("Power exponent of WT on CL (unitless)") # Allard 2020 Table 2: beta Weight on Cl = 0.75 FIX
    e_age_cl <- -0.214; label("Power exponent of AGE on CL (unitless)") # Allard 2020 Table 2: beta Age = -0.214 (RSE 21.1%)
    e_form_fviii_fc_cl <- -0.394; label("Log-scale shift in CL for Elocta vs standard half-life FVIII (unitless)") # Allard 2020 Table 2: beta EHL = -0.394 (RSE 10.1%)
    e_wt_vc <- 0.827; label("Power exponent of WT on V1 (unitless)") # Allard 2020 Table 2: beta Weight on V1 = 0.827 (RSE 4.75%)

    # Between-subject variability. Monolix reports omega as the SD of the random
    # effect; variances are the squares. Q carries no random effect.
    # cov(CL, V1) = 0.599 * 0.349 * 0.232 = 0.048500
    etalcl + etalvc ~ c(
      0.121801,
      0.048500, 0.053824
    ) # Allard 2020 Table 2: omega Cl = 34.9%, omega V1 = 23.2%, Corr(V1, Cl) = 0.599
    etalvp ~ 0.219024 # Allard 2020 Table 2: omega V2 = 46.8%

    # Residual error, Monolix combined1 form y = Cc + (a + b*Cc)*e (Allard 2020 Eq. 3)
    addSd <- 0.00868; label("Additive residual error a (IU/mL)") # Allard 2020 Table 2: a = 0.00868 IU/mL (RSE 8.88%)
    propSd <- 0.108; label("Proportional residual error b (fraction)") # Allard 2020 Table 2: b = 10.8% (RSE 6.6%)
  })
  model({
    # Individual parameters (Allard 2020 Eqs. 1-2 and Table 2 footnote)
    cl <- exp(lcl + etalcl + e_form_fviii_fc_cl * FORM_FVIII_FC) *
      (WT / 64)^e_wt_cl * (AGE / 30)^e_age_cl
    vc <- exp(lvc + etalvc) * (WT / 64)^e_wt_vc
    q <- exp(lq)
    vp <- exp(lvp + etalvp)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # Intravenous FVIII doses enter the central compartment.
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # FVIII:C in IU/mL: dose in IU and V1 in mL.
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd) + combined1()
  })
}
