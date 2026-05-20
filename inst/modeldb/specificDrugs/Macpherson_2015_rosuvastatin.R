Macpherson_2015_rosuvastatin <- function() {
  description <- "Two-compartment population PK model with first-order oral absorption for rosuvastatin in pediatric patients (aged 6 to <18 years) with heterozygous familial hypercholesterolemia (Macpherson 2015 Eur J Clin Pharmacol). Apparent clearance scales with body weight (estimated power exponent 0.352, reference 42 kg) and is 1.41-fold higher in males than females. Residual error is proportional and switches between intensive and sparse PK sampling phases."
  reference <- paste(
    "Macpherson M, Hamren B, Braamskamp MJAM, Kastelein JJP, Lundstrom T, Martin PD.",
    "Population pharmacokinetics of rosuvastatin in pediatric patients with",
    "heterozygous familial hypercholesterolemia.",
    "Eur J Clin Pharmacol. 2016;72(1):19-27.",
    "doi:10.1007/s00228-015-1946-4.",
    sep = " "
  )
  vignette <- "Macpherson_2015_rosuvastatin"
  units <- list(time = "hour", dosing = "mg", concentration = "ng/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight (time-varying through the 2-year CHARON observation window).",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Power covariate on apparent clearance (estimated exponent 0.352, reference 42 kg = CHARON baseline median weight). No weight effect on Vc, Vp, or Q (the paper rejected the fixed allometric form with exponents 0.75 on CL and 1 on V because it gave prediction bias). Time-varying: actual observed weights through the 2-year CHARON study were used so that growth-driven changes in CL/F are captured.",
      source_name        = "WT"
    ),
    SEXF = list(
      description        = "Biological sex indicator (1 = female, 0 = male).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Multiplicative covariate on apparent clearance. The source paper reports a male-vs-female fold-effect (Male CL/F x 1.41; female-children as the implicit reference), so the model applies the factor as 1.41^(1 - SEXF) to keep Table 2's published CL/F = 129 L/h interpretable as the female-at-42-kg typical value. Male children have CL/F ~1.41x higher than female children of the same weight. Worked Table 2 examples reproduced: a 20 kg female has CL/F = 99 L/h, a 20 kg male has CL/F = 140 L/h, a 99 kg male has CL/F = 246 L/h, a 111 kg female has CL/F = 182 L/h.",
      source_name        = "SEXF"
    ),
    SAMPLE_INTENSIVE = list(
      description        = "Per-observation indicator of sampling intensity: 1 = intensive (rich post-dose profile), 0 = sparse (pre-dose / steady-state trough).",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (sparse pre-dose sampling)",
      notes              = "Record-level indicator that switches the proportional residual-error magnitude per observation: intensive sampling 39.4% CV, sparse sampling 59.5% CV (Macpherson 2015 Table 2 final model). In CHARON (n=196) the 12 PK-pilot subjects had intensive 24-h post-dose sampling on Day 0 followed by sparse pre-dose sampling over the 2-year follow-up, while the other 184 subjects had sparse pre-dose sampling only. In Study 4522IL/0086 (n=18) every observation is intensive (rich profiles after single 10/40/80 mg doses and after 7 days of 80 mg once-daily dosing). Set SAMPLE_INTENSIVE = 1 on observations within a rich post-dose profile, 0 otherwise.",
      source_name        = "SAMPLE_INTENSIVE"
    )
  )

  population <- list(
    species          = "human",
    n_subjects       = 214L,
    n_observations   = 2029L,
    n_studies        = 2L,
    age_range        = "6-17 years",
    age_median       = "11 years (CHARON), 14 years (4522IL/0086)",
    weight_range     = "20-116 kg",
    weight_median    = "42 kg (CHARON), 63.2 kg (4522IL/0086)",
    sex_female_pct   = 56,
    race_ethnicity   = c(Caucasian = 89, Black = 3, Asian = 7, Hispanic = 0.5, Other = 1.5),
    disease_state    = "Heterozygous familial hypercholesterolemia (HeFH)",
    dose_range       = "5-80 mg oral once daily; CHARON: 5/10/20 mg once-daily titration; 4522IL/0086: 10/40/80 mg single dose and 80 mg once daily for 7 days",
    regions          = "Multicenter pediatric HeFH centers (CHARON, NCT01078675) and Study 4522IL/0086",
    notes            = "Pooled population PK analysis from two AstraZeneca pediatric studies. Demographics from Macpherson 2015 Table 1. Race was not formally analyzed as a covariate (89% Caucasian, frequency of any other race category <7%). The dataset included 1,735 concentrations from CHARON (median 9 sparse trough samples per subject; 12 PK-pilot subjects had an additional intensive 24-h profile) and 294 concentrations from 4522IL/0086 (median 12 intensive samples per subject)."
  )

  ini({
    # Structural parameters -- Macpherson 2015 Table 2 final covariate model
    # (column "Estimate (relative standard error, %)", final-model column).
    # Reference for CL/F: a female child at WT = 42 kg.
    lka <- log(0.183); label("Absorption rate constant Ka (1/h)")                                  # Table 2 final model: Ka = 0.183 1/h (RSE 10.5%)
    lcl <- log(129);   label("Apparent clearance CL/F at WT = 42 kg, female (L/h)")                # Table 2 final model: CL/F = 129 L/h (RSE 5.70%)
    lvc <- log(303);   label("Apparent central volume Vc/F (L)")                                   # Table 2 final model: Vc/F = 303 L (RSE 28.0%)
    lq  <- log(89.9);  label("Apparent intercompartmental clearance Q/F (L/h)")                    # Table 2 final model: Q/F = 89.9 L/h (RSE 15.6%)
    lvp <- log(5153);  label("Apparent peripheral volume Vp/F (L)")                                # Table 2 final model: Vp/F = 5,153 L (RSE 23.9%)

    # Covariate effects on CL/F (Macpherson 2015 Table 2 final model).
    # Body weight enters as a power function with estimated exponent 0.352;
    # the paper rejected the fixed allometric form (exponents 0.75 on CL, 1 on V)
    # because it produced prediction bias (Results, Covariate analysis).
    e_wt_cl   <- 0.352; label("Power exponent of body weight on CL/F (unitless, reference 42 kg)") # Table 2 final model: Weight CL/F x (WT/42)^theta = 0.352 (RSE 23.9%)
    # Male-vs-female multiplicative effect on CL/F. The paper reports
    # Male CL/F x 1.41 (RSE 6.33%) with female-children as the implicit reference;
    # encoded as 1.41^(1 - SEXF) in model() so that SEXF = 0 (male) yields the
    # 1.41 factor and SEXF = 1 (female) yields 1.0. This preserves Table 2's
    # published CL/F = 129 L/h as the female-at-42-kg typical clearance value
    # and the table's worked examples (140 L/h for a 20 kg male, 99 L/h for a
    # 20 kg female).
    e_sexf_cl <- 1.41;  label("Male fold-effect on CL/F vs female reference (unitless)")           # Table 2 final model: Male CL/F x theta = 1.41 (RSE 6.33%)

    # Inter-individual variability (log-normal). The paper reports CV%; the
    # log-scale variance is omega^2 = log(1 + CV^2):
    #   IIV (%) CL/F = 40.0 -> omega^2 = log(1 + 0.40^2)  = 0.14842
    #   IIV (%) Vc/F = 105  -> omega^2 = log(1 + 1.05^2)  = 0.74295
    #   IIV (%) Q/F  = 64.8 -> omega^2 = log(1 + 0.648^2) = 0.35067
    # The paper estimated IIV on CL/F, Vc/F, and Q/F only; no IIV was reported
    # on Ka or Vp/F.
    etalcl ~ 0.14842                                                                                # Table 2 final model: IIV CL/F = 40.0% CV (RSE 8.09%)
    etalvc ~ 0.74295                                                                                # Table 2 final model: IIV Vc/F = 105% CV (RSE 15.5%)
    etalq  ~ 0.35067                                                                                # Table 2 final model: IIV Q/F  = 64.8% CV (RSE 13.3%)

    # Residual error -- proportional in linear space (the source paper used
    # "proportional in nature (additive on the log scale)"; the log-additive
    # form in NONMEM corresponds to nlmixr2's prop() in linear space). Two
    # magnitudes were estimated, switched by SAMPLE_INTENSIVE in model().
    propSdSparse    <- 0.595; label("Proportional residual error, sparse sampling (fraction)")     # Table 2 final model: Residual error sparse sampling = 59.5% CV (RSE 4.40%)
    propSdIntensive <- 0.394; label("Proportional residual error, intensive sampling (fraction)")  # Table 2 final model: Residual error intensive sampling = 39.4% CV (RSE 7.89%)
  })

  model({
    # Individual PK parameters with body-weight power effect and sex effect on CL/F.
    # CL/F = 129 L/h * (WT/42)^0.352 * 1.41^(1 - SEXF). Worked Table 2 examples:
    #   20 kg female: 129 * (20/42)^0.352 * 1^1 = 99.2 L/h  (paper: 99)
    #   20 kg male:   129 * (20/42)^0.352 * 1.41 = 140.0 L/h (paper: 140)
    #   99 kg male:   129 * (99/42)^0.352 * 1.41 = 246.1 L/h (paper: 246)
    #   111 kg female:129 * (111/42)^0.352      = 180.8 L/h (paper: 182)
    ka <- exp(lka)
    cl <- exp(lcl + etalcl) * (WT / 42)^e_wt_cl * e_sexf_cl^(1 - SEXF)
    vc <- exp(lvc + etalvc)
    q  <- exp(lq  + etalq)
    vp <- exp(lvp)

    # Micro-constants for the two-compartment disposition
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # ODE system: first-order oral absorption from depot, two-compartment central + peripheral.
    d/dt(depot)       <- -ka * depot
    d/dt(central)     <-  ka * depot - kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                               k12 * central - k21 * peripheral1

    # Rosuvastatin plasma concentration. Dose is in mg, volume in L, so
    # central/vc has units of mg/L = ug/mL; multiply by 1000 to convert to
    # ng/mL (the bioanalytical units in both source studies: CHARON LLOQ
    # 0.02 ng/mL, range 0.02-20 ng/mL; 4522IL/0086 LLOQ 0.1 ng/mL, range
    # 0.1-30 ng/mL; Methods, Pharmacokinetic assessments and analysis).
    Cc <- (central / vc) * 1000

    # Residual-error magnitude switched by SAMPLE_INTENSIVE:
    #   SAMPLE_INTENSIVE = 1 -> intensive sampling, 39.4% CV
    #   SAMPLE_INTENSIVE = 0 -> sparse sampling,    59.5% CV
    propSd <- propSdIntensive * SAMPLE_INTENSIVE + propSdSparse * (1 - SAMPLE_INTENSIVE)
    Cc ~ prop(propSd)
  })
}
