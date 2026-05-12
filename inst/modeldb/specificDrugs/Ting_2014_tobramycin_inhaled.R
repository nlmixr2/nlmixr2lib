Ting_2014_tobramycin_inhaled <- function() {
  description <- "Two-compartment population PK model for inhaled tobramycin powder (TIP / TOBI Podhaler) in cystic fibrosis patients (Ting 2014), with first-order absorption from a depot compartment and apparent (post-bioavailability) clearance and volumes. Body mass index (BMI) and baseline FEV1 percent-predicted are power-form covariates on apparent central volume of distribution (reference 18.8 kg/m^2 and 62.1 % respectively)."
  reference <- paste(
    "Ting L, Aksenov S, Bhansali SG, Ramakrishna R, Tang P, Geller DE. (2014).",
    "Population pharmacokinetics of inhaled tobramycin powder in cystic fibrosis patients.",
    "CPT Pharmacometrics Syst Pharmacol 3(9):e99. doi:10.1038/psp.2013.76"
  )
  vignette <- "Ting_2014_tobramycin_inhaled"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    BMI = list(
      description        = "Body mass index at baseline",
      units              = "kg/m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form effect on apparent central volume Vd/F: (BMI / 18.8)^0.624 (Ting 2014 Results, equation between Tables 1 and 2). Reference 18.8 kg/m^2 is the population median across the combined three-study cohort (Table 1). Baseline / time-fixed.",
      source_name        = "BMI"
    ),
    FEV1_PCTPRED = list(
      description        = "Baseline forced expiratory volume in 1 second as percent of the predicted value",
      units              = "% predicted",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power-form effect on apparent central volume Vd/F: (FEV1_PCTPRED / 62.1)^-0.303 (Ting 2014 Results, equation between Tables 1 and 2). Reference 62.1 % predicted is the population median across the combined three-study cohort (Table 1). Source paper uses the unsubscripted notation 'FEV1% predicted'; the canonical column FEV1_PCTPRED carries the percent-predicted value as a number (e.g. 62.1, not 0.621). Baseline / time-fixed.",
      source_name        = "FEV1% predicted"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 139L,
    n_studies      = 3L,
    n_observations = 662L,
    age_range      = "6-58 years",
    age_median     = "17 years (SD 10.8) across the combined cohort; sub-study medians 21 (TPI001), 14 (C2301), 31 (C2302)",
    age_groups     = c("6-11 yr" = 20.9, "12-17 yr" = 31.7, ">=18 yr" = 47.5),
    weight_range   = "16.2-100.9 kg",
    weight_median  = "49.5 kg (SD 17.7)",
    bmi_range      = "11.4-31 kg/m^2",
    bmi_median     = "18.8 kg/m^2 (SD 4.1)",
    crcl_range     = "63.9-222.5 mL/min",
    crcl_median    = "112.5 mL/min (SD 29.5)",
    fev1_pctpred_range  = "24.1-119.7 % predicted",
    fev1_pctpred_median = "62.1 % predicted (SD 20.6)",
    sex_female_pct = 53.2,
    race_ethnicity = c(White = 86.3, Black = 2.2, Hispanic = 8.6, Other = 2.9),
    disease_state  = "Cystic fibrosis patients with Pseudomonas aeruginosa airway infection (ages 6+), pooled from one phase I (TPI001) and two phase III (C2301, C2302) studies.",
    dose_range     = "Single doses of 28, 56, 84, or 112 mg TIP (phase I dose-escalation TPI001); 112 mg b.i.d. multiple doses for 28 days/cycle in phase III studies C2301 and C2302.",
    administration = "Inhaled dry tobramycin powder delivered via the TOBI Podhaler (Novartis) dry-powder inhaler. Dose entered the model as the prescribed inhaled mass; absolute bioavailability is unidentifiable and absorbed into the apparent parameters CL/F, Vd/F, Q/F, V2/F.",
    regions        = "Multinational (study locations not enumerated in the paper).",
    notes          = "Baseline demographics from Ting 2014 Table 1 ('Combined' column). Three pooled studies: TPI001 (phase I dose-escalation, n=64), C2301 (phase III, n=62), C2302 (phase III, n=13)."
  )

  ini({
    # Structural parameters (Ting 2014 Table 2, "Final model" column). All
    # values are apparent post-bioavailability quantities (CL/F, Vd/F, Q/F,
    # V2/F); the inhaled dose enters the depot as the prescribed mg of TIP
    # and absolute F is folded into the apparent parameters.
    lka <- log(2.39); label("First-order absorption rate ka (1/h)")             # Ting 2014 Table 2: ka = 2.39 1/h
    lcl <- log(14.5); label("Apparent clearance CL/F (L/h)")                    # Ting 2014 Table 2: CL/F = 14.5 L/h
    lvc <- log(85.1); label("Apparent central volume Vd/F at reference covariates (L)")  # Ting 2014 Table 2: Vd/F = 85.1 L
    lq  <- log(6.43); label("Apparent intercompartmental clearance Q/F (L/h)")  # Ting 2014 Table 2: Q/F = 6.43 L/h
    lvp <- log(210);  label("Apparent peripheral volume V2/F (L)")              # Ting 2014 Table 2: V2/F = 210 L

    # Covariate effects on apparent central volume Vd/F (Ting 2014 Results,
    # power-form equation just before Table 2):
    #   Vd/F = 85.1 * (BMI / 18.8)^0.624 * (FEV1% / 62.1)^-0.303
    e_bmi_vc  <-  0.624; label("Power exponent of (BMI/18.8) on Vd/F (unitless)")        # Ting 2014 Table 2: BMI on Vd/F = 0.624
    e_fev1_vc <- -0.303; label("Power exponent of (FEV1_PCTPRED/62.1) on Vd/F (unitless)")  # Ting 2014 Table 2: FEV1% on Vd/F = -0.303

    # Inter-individual variability (Ting 2014 Table 2 reports IIV as
    # variance on log-transformed parameters). CL and Vc are correlated;
    # the paper reports their covariance separately. ka is independent.
    # Inter-occasion variability of CL/F (0.078 variance, identical across
    # all four occasions per Table 2) is not encoded in the library model;
    # see Assumptions and deviations in the vignette.
    etalcl + etalvc ~ c(
      0.164,
      0.123, 0.152
    )                                       # Ting 2014 Table 2: IIV CL/F variance 0.164, IIV Vd/F variance 0.152, COV(CL/F,Vd/F) 0.123
    etalka ~ 0.129                          # Ting 2014 Table 2: IIV ka variance 0.129

    # Residual error (Ting 2014 Table 2). The paper reports two proportional
    # residual SDs: 0.073 for the phase I study TPI001 with prospective
    # sampling, and 0.308 for the phase III studies C2301/C2302 where many
    # sampling times had to be imputed (paper Results, second paragraph
    # under "Population pharmacokinetic model"). The TPI001 value reflects
    # true assay precision and is used here as the default residual error
    # for clean-sampling simulation; see vignette Assumptions and deviations.
    propSd <- 0.073; label("Proportional residual SD (fraction; TPI001 clean-sampling estimate)")  # Ting 2014 Table 2: proportional RUV (TPI001) = 0.073
    addSd  <- 0.007; label("Additive residual SD (mg/L)")                                          # Ting 2014 Table 2: additive RUV = 0.007 ug/mL = 0.007 mg/L
  })
  model({
    # Individual parameters. Vd/F carries the BMI and FEV1% power covariate
    # effects (Ting 2014 Results equation). CL/F, Q/F, V2/F, and ka have no
    # statistically significant covariates in the final model.
    ka <- exp(lka + etalka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc) * (BMI / 18.8)^e_bmi_vc * (FEV1_PCTPRED / 62.1)^e_fev1_vc
    q  <- exp(lq)
    vp <- exp(lvp)

    # Two-compartment micro-constants with first-order absorption.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                                k12 * central - k21 * peripheral1

    # Dose in mg into depot, vc in L -> central / vc has units mg/L (= ug/mL,
    # the unit used throughout Ting 2014 for serum tobramycin concentration).
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
