Tod_1998_amikacin <- function() {
  description <- "Two-compartment intravenous population PK model for amikacin in febrile, severely neutropenic adults with hematological malignancies (Tod 1998); clearance modeled as the sum of a non-renal intercept and a Cockcroft-Gault-like renal component with sex-stratified slope coefficient (males theta_1, females theta_2), age-correction factor (theta_3 - AGE/100), and Cockcroft-Gault-like renal-function ratio (WT / CREAT). Power-variance residual-error model."
  reference <- "Tod M, Lortholary O, Seytre D, Semaoun R, Uzzan B, Guillevin L, Casassus P, Petitjean O. Population pharmacokinetic study of amikacin administered once or twice daily to febrile, severely neutropenic adults. Antimicrob Agents Chemother. 1998;42(4):849-856. doi:10.1128/aac.42.4.849"
  vignette <- "Tod_1998_amikacin"
  units <- list(time = "hour", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    WT = list(
      description        = "Body weight at baseline",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Used as bw in the Cockcroft-Gault-like renal-CL covariate term (WT / CREAT). Tod 1998 Table 1: mean (SD) 65.8 (13.0) kg in the o.d. group and 68.1 (13.1) kg in the b.i.d. group; range 44-106 kg. Time-fixed at baseline in the published model.",
      source_name        = "bw"
    ),
    AGE = list(
      description        = "Subject age at baseline",
      units              = "years",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Used as age in the renal-CL age-correction factor (theta_3 - AGE/100). Tod 1998 Table 1: mean (SD) 50.2 (16.8) years in the o.d. group and 51.3 (16.0) years in the b.i.d. group; range 18-85 years. At AGE = 0 the age-correction factor equals theta_3 = 0.985; the linear extrapolation reaches zero at AGE = 100 x theta_3 = 98.5 years, so the structural CL is only physiologically meaningful for adult ages below that threshold.",
      source_name        = "age"
    ),
    CREAT = list(
      description        = "Serum creatinine concentration",
      units              = "umol/L",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Used as S_CR (the patient's measured serum creatinine) in the Cockcroft-Gault-like renal-CL covariate term (WT / CREAT). Tod 1998 Discussion (p. 853) reports the per-patient serum creatinine assumed log-normally distributed with mean (SD) 85 (37) umol/L (used for the population simulation of Table 5). The full per-patient SCR vector underlies the CLCR mean (SD) 91 (36) mL/min in the b.i.d. group and 104 (39) mL/min in the o.d. group (Table 1). Time-fixed at baseline in the published model.",
      source_name        = "SCR"
    ),
    SEXF = list(
      description        = "Biological sex indicator (1 = female, 0 = male)",
      units              = "(binary)",
      type               = "binary",
      reference_category = "0 (male)",
      notes              = "Tod 1998 Table 4 reports two distinct renal-CL slope coefficients (theta_1 = 0.797 for males, theta_2 = 0.640 for females). The structural model selects theta_i based on SEXF; the canonical SEXF (1 = female) maps directly to the paper's 'i = 2 for women'. Removing the sex covariate produced a significantly poorer fit (Table 3 step 6).",
      source_name        = "sex"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 57L,
    n_studies      = 1L,
    age_range      = "18-85 years",
    age_median     = "o.d. 51 years; b.i.d. 50 years",
    weight_range   = "44-106 kg",
    weight_median  = "o.d. 64 kg; b.i.d. 66 kg",
    sex_female_pct = 39,
    race_ethnicity = "Not reported (single-centre French cohort, Hopital Avicenne, Bobigny)",
    disease_state  = "Febrile severe neutropenia (neutrophil count <500/mm^3) under treatment for a primary hematological disorder. Diagnoses: acute myeloblastic leukemia (n=18), non-Hodgkin's lymphoma (n=21), acute lymphoblastic leukemia (n=8), myeloma (n=7), Hodgkin's lymphoma (n=1), agranulocytosis (n=1), aplastic anemia (n=1).",
    dose_range     = "7.5 mg/kg b.i.d. (n=29) or 20 mg/kg o.d. (n=28) intravenous infusion over 0.5 h; both regimens given alongside piperacillin (b.i.d. cohort) or piperacillin-tazobactam (o.d. cohort) t.i.d.",
    regions        = "France (single centre, Hopital Avicenne)",
    n_observations = 278L,
    renal_function = "Cockcroft-Gault CLCR mean (SD) 91 (36) mL/min in the b.i.d. group and 104 (39) mL/min in the o.d. group; range 20-213 mL/min (Tod 1998 Table 1).",
    notes          = "Demographics from Tod 1998 Table 1. 278 serum amikacin samples (93 peak, 117 trough, 68 intermediate); median 4 samples per subject (range 1-14). Concentrations were measured by enzyme-multiplied immunoassay (EMIT, Cobas Roche); LOQ 2.5 mg/L. Pregnant women and HIV-infected patients were excluded. The dosing-regimen covariate was not retained in the final model (Table 3 steps 7-14)."
  )

  ini({
    # Structural CL = theta_i x 6 x (theta_3 - AGE/100) x (WT / CREAT) + theta_4,
    # with theta_i = theta_1 for males (SEXF = 0) or theta_2 for females
    # (SEXF = 1). Tod 1998 Table 4 final-model column. The factor 6 is a
    # fitting constant fitted alongside theta_1..theta_4; it does NOT equal
    # the L/h <-> mL/min unit conversion (0.06). CREAT in umol/L and WT in kg
    # so that theta_i x 6 x (theta_3 - AGE/100) x (WT/CREAT) carries units
    # L/h when combined with the dimensional absorption of theta_i.
    cl_renal_male   <- 0.797 ; label("Male renal-CL slope coefficient (Tod theta_1)")     # Tod 1998 Table 4: theta_1 = 0.797 (SE 0.181)
    cl_renal_female <- 0.640 ; label("Female renal-CL slope coefficient (Tod theta_2)")   # Tod 1998 Table 4: theta_2 = 0.640 (SE 0.162)
    e_age_cl_renal  <- 0.985 ; label("Age-correction intercept on renal CL (Tod theta_3; unitless)") # Tod 1998 Table 4: theta_3 = 0.985 (SE 0.106)
    lcl             <- log(1.66) ; label("Log non-renal CL intercept (L/h) (Tod theta_4); the matching IIV etalcl applies multiplicatively to the sum (non-renal + renal) per the Delattre 2010 precedent") # Tod 1998 Table 4: theta_4 = 1.66 L/h (SE 0.34)

    # Two-compartment structural volumes and inter-compartmental clearance
    lvc <- log(8.92) ; label("Log central volume V1 (L)")                       # Tod 1998 Table 4: V1 = 8.92 L (SE 1.17)
    lq  <- log(4.43) ; label("Log inter-compartmental clearance CL_D (L/h)")    # Tod 1998 Table 4: CL_D = 4.43 L/h (SE 0.76)
    lvp <- log(11.4) ; label("Log peripheral volume V_t (L)")                   # Tod 1998 Table 4: V_t = 11.4 L (SE 1.3)

    # Inter-individual variability (Tod 1998 Table 4 final-model CV%);
    # omega^2 = log(CV^2 + 1) for log-normal etas. The paper notes
    # 'Allowing for covariance between eta's did not improve the fit',
    # so all four IIVs are independent.
    etalcl ~ 0.04316  # log(0.21^2 + 1); 21% CV on total CL
    etalvc ~ 0.02225  # log(0.15^2 + 1); 15% CV on V1
    etalq  ~ 0.08618  # log(0.30^2 + 1); 30% CV on CL_D
    etalvp ~ 0.06062  # log(0.25^2 + 1); 25% CV on V_t

    # Power-variance residual-error model. Tod 1998 reports
    #   C_obs = C_pred + eps * C_pred^b  with var(eps) = sigma_eps^2 = 0.189
    #   and b = 0.939
    # encoded in rxode2 / nlmixr2 as 'Cc ~ pow(propSd, powExp)' where the
    # standard deviation of (C_obs - C_pred) equals propSd x C_pred^powExp.
    # So propSd = sigma_eps = sqrt(0.189) = 0.4347 and powExp = b = 0.939.
    propSd <- 0.4347 ; label("Power-error SD coefficient (Tod sigma_eps = sqrt(0.189))") # Tod 1998 Table 4: sigma_eps^2 = 0.189 (SE 0.083)
    powExp <- 0.939  ; label("Power-error exponent (Tod b; unitless)")                   # Tod 1998 Table 4: b = 0.939 (SE 0.060)
  })

  model({
    # Sex-specific renal-CL slope coefficient (Tod theta_i selector).
    cl_renal_coef <- (1 - SEXF) * cl_renal_male + SEXF * cl_renal_female

    # Renal CL contribution: theta_i x 6 x (theta_3 - AGE/100) x (WT/CREAT).
    cl_renal <- cl_renal_coef * 6 * (e_age_cl_renal - AGE / 100) * (WT / CREAT)

    # Total CL = non-renal intercept + renal component; log-normal IIV is
    # applied multiplicatively to the sum (Tod 1998 'CL_j = CL * exp(eta_CL,j)').
    cl <- (exp(lcl) + cl_renal) * exp(etalcl)
    vc <- exp(lvc + etalvc)
    q  <- exp(lq  + etalq)
    vp <- exp(lvp + etalvp)

    # Two-compartment IV micro-constants.
    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                   k12 * central - k21 * peripheral1

    # Dose in mg, volumes in L -> central / vc has units mg/L.
    Cc <- central / vc
    Cc ~ pow(propSd, powExp)
  })
}
