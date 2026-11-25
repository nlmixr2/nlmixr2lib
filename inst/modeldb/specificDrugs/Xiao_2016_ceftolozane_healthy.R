Xiao_2016_ceftolozane_healthy <- function() {
  description <- paste(
    "Two-compartment population PK model for ceftolozane with an added",
    "epithelial lining fluid (ELF) compartment, fitted to plasma and",
    "bronchoalveolar-lavage ELF concentrations from 25 healthy adults given",
    "ceftolozane/tazobactam 1000/500 mg as a 60-minute intravenous infusion",
    "every 8 hours (Xiao 2016, Table 1). Plasma-to-ELF transfer and",
    "ELF-to-plasma return are separate, asymmetric clearances (Qce, Qec)",
    "into an ELF compartment whose volume is fixed at 1 L, so the parametric",
    "plasma-to-ELF penetration ratio is Qce / Qec (0.51). No covariates were",
    "retained. Companion models: Xiao_2016_tazobactam_healthy (same",
    "structure, fitted separately) and Xiao_2016_ceftolozane_pneumonia (the",
    "cIAI-patient plasma PK combined with these ELF kinetics for the",
    "nosocomial-pneumonia dose simulations)."
  )
  reference <- paste(
    "Xiao AJ, Miller BW, Huntington JA, Nicolau DP.",
    "Ceftolozane/tazobactam pharmacokinetic/pharmacodynamic-derived dose",
    "justification for phase 3 studies in patients with nosocomial",
    "pneumonia. J Clin Pharmacol. 2016;56(1):56-66. doi:10.1002/jcph.566.",
    "Structure: Results 'Plasma-to-ELF Kinetics' mass-balance equations and",
    "Figure 2. Parameter estimates: Table 1 (ceftolozane panel). The",
    "underlying ELF study is Chandorkar 2012, J Antimicrob Chemother",
    "67:2463-2469 (doi:10.1093/jac/dks246)."
  )
  vignette <- "Xiao_2016_ceftolozane_tazobactam"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list()

  compartmentData <- list(
    central = list(analyte = "ceftolozane", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "ceftolozane", units = "mg", specimen = "tissue", verified = TRUE),
    elf = list(analyte = "ceftolozane", units = "mg", specimen = "epithelial lining fluid", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 25L,
    n_studies = 1L,
    n_concentrations = "150 plasma + 25 ELF (one bronchoalveolar lavage per subject)",
    age_range = "21-47 years",
    age_median = "31 years",
    weight_range = "56.5-94.5 kg",
    weight_median = "78 kg",
    sex_female_pct = 44,
    disease_state = "Healthy adult volunteers",
    dose_range = paste(
      "Ceftolozane/tazobactam 1000/500 mg (1.5 g) as a 60-minute intravenous",
      "infusion every 8 hours for 3 doses; plasma and bronchoalveolar-lavage",
      "samples taken 0-8 hours after the third dose."
    ),
    renal_function = "Creatinine clearance median 119 mL/min (range 77.7-155.2)",
    notes = paste(
      "Methods, 'Plasma-to-ELF Penetration Kinetics': 25 healthy subjects,",
      "14 men and 11 women; median [min, max] age 31 [21, 47] years, body",
      "weight 78 [56.5, 94.5] kg and CrCL 119 [77.7, 155.2] mL/min. Each",
      "subject contributed a single ELF sample. Estimation in NONMEM 7.2.",
      "No covariate effect was identified for ceftolozane in this small",
      "population (Results)."
    )
  )

  ini({
    lcl <- log(6.3); label("Ceftolozane clearance (L/h)") # Table 1 ceftolozane: CL = 6.3 L/h (SE 0.187)
    lvc <- log(10.3); label("Ceftolozane central volume of distribution (L)") # Table 1 ceftolozane: Vc = 10.3 L (SE 0.57)
    lq <- log(3.39); label("Ceftolozane central-peripheral inter-compartmental clearance Q2 (L/h)") # Table 1 ceftolozane: Q2 = 3.39 L/h (SE 0.711)
    lvp <- log(3.76); label("Ceftolozane peripheral volume of distribution V2 (L)") # Table 1 ceftolozane: V2 = 3.76 L (SE 0.392)
    lq_c_elf <- log(0.85); label("Ceftolozane one-way central-to-ELF clearance Qce (L/h)") # Table 1 ceftolozane: Qce = 0.85 L/h (SE 0.362)
    lq_elf_c <- log(1.66); label("Ceftolozane one-way ELF-to-central clearance Qec (L/h)") # Table 1 ceftolozane: Qec = 1.66 L/h (SE 0.677)
    lv_elf <- fixed(log(1)); label("ELF compartment volume Ve (L)") # Results equations: 'Ve (fixed as 1)'

    etalcl ~ 0.0179 # Table 1 ceftolozane: 'Variance of IIV on CL' = 0.0179 (SE 0.00529)
    etalvc ~ 0.0498 # Table 1 ceftolozane: 'Variance of IIV on Vc' = 0.0498 (SE 0.0174)
    etalq_c_elf ~ 0.174 # Table 1 ceftolozane: 'Variance of IIV on Qce' = 0.174 (SE 0.0454); 42% CV on PR (Results)

    propSd <- 0.07765307; label("Proportional residual error, plasma (fraction)") # Table 1 ceftolozane: 'Variance of RVP' = 0.00603 -> SD sqrt(0.00603)
    addSd <- 0.5449771; label("Additive residual error, plasma (mg/L)") # Table 1 ceftolozane: 'Variance of RVA' = 0.297 -> SD sqrt(0.297)
    propSd_Celf <- 0.07765307; label("Proportional residual error, ELF (fraction)") # Table 1 ceftolozane: one RVP for both matrices, sqrt(0.00603)
    addSd_Celf <- 0.5449771; label("Additive residual error, ELF (mg/L)") # Table 1 ceftolozane: one RVA for both matrices, sqrt(0.297)
  })

  model({
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q <- exp(lq)
    vp <- exp(lvp)
    q_c_elf <- exp(lq_c_elf + etalq_c_elf)
    q_elf_c <- exp(lq_elf_c)
    v_elf <- exp(lv_elf)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k_central_elf <- q_c_elf / vc
    k_elf_central <- q_elf_c / v_elf

    d/dt(central) <- -(kel + k12 + k_central_elf) * central + k21 * peripheral1 + k_elf_central * elf
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(elf) <- k_central_elf * central - k_elf_central * elf

    Cc <- central / vc
    Celf <- elf / v_elf

    Cc ~ add(addSd) + prop(propSd)
    Celf ~ add(addSd_Celf) + prop(propSd_Celf)
  })
}
