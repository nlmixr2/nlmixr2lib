Xiao_2016_tazobactam_healthy <- function() {
  description <- paste(
    "Two-compartment population PK model for tazobactam with an added",
    "epithelial lining fluid (ELF) compartment, fitted to plasma and",
    "bronchoalveolar-lavage ELF concentrations from 25 healthy adults given",
    "ceftolozane/tazobactam 1000/500 mg as a 60-minute intravenous infusion",
    "every 8 hours (Xiao 2016, Table 1). Plasma-to-ELF transfer and",
    "ELF-to-plasma return are separate, asymmetric clearances (Qce, Qec)",
    "into an ELF compartment whose volume is fixed at 1 L, so the parametric",
    "plasma-to-ELF penetration ratio is Qce / Qec (0.46). Creatinine",
    "clearance acts on CL through a power function and body weight",
    "proportionally on central volume. Companion models:",
    "Xiao_2016_ceftolozane_healthy (same structure, fitted separately) and",
    "Xiao_2016_tazobactam_pneumonia (the cIAI-patient plasma PK combined",
    "with these ELF kinetics for the nosocomial-pneumonia dose simulations)."
  )
  reference <- paste(
    "Xiao AJ, Miller BW, Huntington JA, Nicolau DP.",
    "Ceftolozane/tazobactam pharmacokinetic/pharmacodynamic-derived dose",
    "justification for phase 3 studies in patients with nosocomial",
    "pneumonia. J Clin Pharmacol. 2016;56(1):56-66. doi:10.1002/jcph.566.",
    "Structure: Results 'Plasma-to-ELF Kinetics' mass-balance equations and",
    "Figure 2. Parameter estimates: Table 1 (tazobactam panel). The",
    "underlying ELF study is Chandorkar 2012, J Antimicrob Chemother",
    "67:2463-2469 (doi:10.1093/jac/dks246)."
  )
  vignette <- "Xiao_2016_ceftolozane_tazobactam"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  covariateData <- list(
    CRCL = list(
      description = "Baseline creatinine clearance",
      units = "mL/min",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Table 1 tazobactam: 'Exponent of (CrCL/119) on CL' = 0.438. The",
        "reference 119 mL/min is the cohort median (Methods). The paper does",
        "not restate the estimating formula; the companion plasma model it",
        "builds on (Chandorkar 2015) used Cockcroft-Gault on actual body",
        "weight, not normalised to body surface area."
      ),
      source_name = "CrCL"
    ),
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Table 1 tazobactam: 'Vc, L, at 70-kg body weight' = 15.0. Results:",
        "'CrCL and body weight were identified as significant covariates of",
        "clearance and volume of distribution, respectively'. Table 1 lists",
        "no weight exponent among the estimated parameters, so the exponent",
        "was not estimated; it is encoded as a fixed proportional (exponent",
        "1) scaling on WT/70."
      ),
      source_name = "body weight"
    )
  )

  compartmentData <- list(
    central = list(analyte = "tazobactam", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "tazobactam", units = "mg", specimen = "tissue", verified = TRUE),
    elf = list(analyte = "tazobactam", units = "mg", specimen = "epithelial lining fluid", verified = TRUE)
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
      "samples taken 0-8 hours after the third dose. The tazobactam component",
      "is 500 mg."
    ),
    renal_function = "Creatinine clearance median 119 mL/min (range 77.7-155.2)",
    notes = paste(
      "Methods, 'Plasma-to-ELF Penetration Kinetics': 25 healthy subjects,",
      "14 men and 11 women; median [min, max] age 31 [21, 47] years, body",
      "weight 78 [56.5, 94.5] kg and CrCL 119 [77.7, 155.2] mL/min. Each",
      "subject contributed a single ELF sample. Estimation in NONMEM 7.2."
    )
  )

  ini({
    lcl <- log(24.5); label("Tazobactam clearance at CrCL = 119 mL/min (L/h)") # Table 1 tazobactam: CL = 24.5 L/h (SE 0.578)
    lvc <- log(15.0); label("Tazobactam central volume of distribution at 70 kg (L)") # Table 1 tazobactam: 'Vc, L, at 70-kg body weight' = 15.0 (SE 0.557)
    lq <- log(3.91); label("Tazobactam central-peripheral inter-compartmental clearance Q2 (L/h)") # Table 1 tazobactam: Q2 = 3.91 L/h (SE 0.823)
    lvp <- log(4.42); label("Tazobactam peripheral volume of distribution V2 (L)") # Table 1 tazobactam: V2 = 4.42 L (SE 0.483)
    lq_c_elf <- log(0.792); label("Tazobactam one-way central-to-ELF clearance Qce (L/h)") # Table 1 tazobactam: Qce = 0.792 L/h (SE 0.14)
    lq_elf_c <- log(1.72); label("Tazobactam one-way ELF-to-central clearance Qec (L/h)") # Table 1 tazobactam: Qec = 1.72 L/h (SE 0.339)
    lv_elf <- fixed(log(1)); label("ELF compartment volume Ve (L)") # Results equations: 'Ve (fixed as 1)'

    e_crcl_cl <- 0.438; label("Power exponent of CrCL on tazobactam CL, CRCL/119 (unitless)") # Table 1 tazobactam: 'Exponent of (CrCL/119) on CL' = 0.438 (SE 0.112)
    e_wt_vc <- fixed(1); label("Power exponent of body weight on tazobactam Vc, WT/70 (unitless)") # Table 1 tazobactam: Vc 'at 70-kg body weight'; no exponent estimated

    etalcl ~ 0.00412 # Table 1 tazobactam: 'Variance of IIV on CL' = 0.00412 (SE 0.00228)
    etalq_elf_c ~ 0.229 # Table 1 tazobactam: 'Variance of IIV on Qec' = 0.229 (SE 0.0689); 48% CV on PR (Results)

    propSd <- 0.1526434; label("Proportional residual error, plasma (fraction)") # Table 1 tazobactam: 'Variance of RVP' = 0.0233 -> SD sqrt(0.0233)
    propSd_Celf <- 0.1526434; label("Proportional residual error, ELF (fraction)") # Table 1 tazobactam: one RVP for both matrices, sqrt(0.0233)
  })

  model({
    cl <- exp(lcl + etalcl) * (CRCL / 119)^e_crcl_cl
    vc <- exp(lvc) * (WT / 70)^e_wt_vc
    q <- exp(lq)
    vp <- exp(lvp)
    q_c_elf <- exp(lq_c_elf)
    q_elf_c <- exp(lq_elf_c + etalq_elf_c)
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

    Cc ~ prop(propSd)
    Celf ~ prop(propSd_Celf)
  })
}
