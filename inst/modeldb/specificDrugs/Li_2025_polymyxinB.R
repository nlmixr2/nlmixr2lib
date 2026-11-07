Li_2025_polymyxinB <- function() {
  description <- paste(
    "Coupled plasma + epithelial-lining-fluid (ELF) population PK model for polymyxin B given",
    "intravenously and/or by nebulization to critically ill adults with multidrug-resistant",
    "Gram-negative pneumonia (Li 2025). Two two-compartment models with linear elimination:",
    "an IV infusion enters the plasma central compartment (central <-> peripheral1, clearance",
    "from central), a nebulized dose enters the ELF central compartment (elf <-> lung, clearance",
    "from elf) with a fraction F = 0.557 reaching the lung, and the two central compartments",
    "exchange through a single inter-compartmental clearance Q = 0.03 L/h. Age is a power",
    "covariate on the ELF central volume and serum albumin a power covariate on the ELF",
    "clearance, both normalised to the cohort medians (66 years, 32.1 g/L). Exponential IIV on",
    "the plasma and ELF central clearances and volumes; one proportional residual error (43%)",
    "estimated jointly for the plasma and ELF concentrations. Fit in NONMEM 7.5 (FOCE-I)."
  )
  reference <- paste(
    "Li X, Zhou L, Wang D, Wu Q, Huang X, Zhang H, Wu W, Liu M, Wu X, Qiu H, Cheng Y.",
    "Population pharmacokinetics study on nebulized and intravenous administration of",
    "polymyxin B in patients with pneumonia caused by multidrug-resistant gram-negative",
    "bacteria.",
    "Antimicrob Agents Chemother. 2025;69(5):e00044-25.",
    "doi:10.1128/aac.00044-25. PMID 40237505. PMCID PMC12057357."
  )
  vignette <- "Li_2025_polymyxinB"
  # Declared explicitly because the registry's automatic detection only
  # recognises `depot` and `central`: a nebulized dose targets cmt = "elf" (the
  # ELF central compartment, Figure 2 compartment 3) and an IV infusion
  # cmt = "central" (Figure 2 compartment 1).
  dosing <- c("central", "elf")
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # What each ODE state holds, in what amount units, in what biological
  # matrix. Li 2025 Figure 2 and Methods: polymyxin B1 + B2 summed as total
  # polymyxin B in plasma and in urea-corrected bronchoalveolar lavage (ELF).
  compartmentData <- list(
    central = list(analyte = "polymyxinB", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "polymyxinB", units = "mg", specimen = "tissue", verified = TRUE),
    elf = list(analyte = "polymyxinB", units = "mg", specimen = "epithelial lining fluid", verified = TRUE),
    lung = list(analyte = "polymyxinB", units = "mg", specimen = "tissue", verified = TRUE)
  )

  covariateData <- list(
    AGE = list(
      description = "Age",
      units = "years",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power covariate on the ELF central volume V3, V3 = 0.896 * (AGE / 66)^1.28. Table 2",
        "prints only 'AGE on V3 = 1.28' (RSE 23%); the paper gives neither the functional form",
        "nor the normalising value. A power model on the Table 1 cohort median (66 years) is",
        "the only reading under which an estimate of 1.28 is plausible (a linear or",
        "exponential per-year coefficient of 1.28 would change V3 by orders of magnitude",
        "across the observed age range). Baseline (collected at study onset)."
      ),
      source_name = "AGE"
    ),
    ALB = list(
      description = "Serum albumin",
      units = "g/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Power covariate on the ELF clearance CL3, CL3 = 0.186 * (ALB / 32.1)^-0.893. Table 2",
        "prints only 'ALB on CL3 = -0.893' (RSE 17%) without the equation or the normalising",
        "value; the Table 1 cohort median (32.1 g/L) is used as the reference. Under this",
        "reading a lower albumin gives a higher ELF clearance, which is the mechanism the",
        "Discussion proposes. Baseline (collected at study onset)."
      ),
      source_name = "ALB"
    )
  )

  covariatesDataExcluded <- list(
    OTHER_SCREENED = list(
      description = "Remaining candidate covariates in the stepwise search",
      units = "(various)",
      type = "continuous",
      notes = paste(
        "Sex, weight, serum urea, Cockcroft-Gault CrCL, APACHE II, SOFA, total bilirubin, AST,",
        "ALT, GGT, ALP, IL-6, CRP, sepsis, duration of mechanical ventilation, nebulizer type",
        "(jet vs vibrating mesh) and immunosuppressive therapy were screened (Methods,",
        "forward p < 0.05 / backward p < 0.01) and none was retained (Table S1). Grouped",
        "because the paper reports the screen collectively."
      )
    )
  )

  population <- list(
    species = "human",
    n_subjects = 76L,
    n_studies = 1L,
    age_range = "median 66.0 years (IQR 55.0-76.3)",
    weight_range = "median 65.0 kg (IQR 54.38-70.0)",
    sex_female_pct = 22.4,
    race_ethnicity = "Chinese (single-centre cohort; race not otherwise reported)",
    disease_state = paste(
      "Critically ill, mechanically ventilated adults with pneumonia caused by",
      "microbiologically proven multidrug-resistant Gram-negative bacteria (Acinetobacter",
      "baumannii 19.7%, Klebsiella pneumoniae 15.8%, Pseudomonas aeruginosa 18.4%, two or",
      "more pathogens 46.1%); sepsis 55.3%; APACHE II median 21, SOFA median 7. Patients on",
      "continuous renal replacement therapy were excluded."
    ),
    dose_range = paste(
      "Nebulized polymyxin B 25 mg in 5 mL every 12 h over 15 min (vibrating-mesh or jet",
      "nebulizer) for 3 days, alone or with IV polymyxin B (loading 100-150 mg, maintenance",
      "40-75 mg every 12 h as ~1 h infusions). Median daily IV dose 1.82 mg/kg, daily",
      "nebulized dose 0.77 mg/kg."
    ),
    regions = "China (Fujian Medical University Union Hospital, Fuzhou)",
    albumin = "median 32.1 g/L (IQR 28.5-37.5)",
    renal_function = "Cockcroft-Gault CrCL median 55.7 mL/min (IQR 35.4-101.3)",
    n_observations = paste(
      "156 ELF concentrations (urea-corrected bronchoalveolar lavage) from all 76 patients and",
      "188 plasma concentrations from 65 patients, at steady state (after at least four doses),",
      "pre-dose and 1, 2, 4, 6, 8 and 12 h after administration."
    ),
    notes = paste(
      "Prospective single-centre study, January 2022 - October 2024. Demographics in Li 2025",
      "Table 1. Polymyxin B1 + B2 summed, LC-MS/MS. NONMEM 7.5, FOCE with interaction;",
      "1,000-sample bootstrap and VPC (Figure S1)."
    )
  )

  ini({
    # Plasma two-compartment model -- Li 2025 Table 2 (final model). Figure 2:
    # compartment 1 = plasma central, compartment 2 = plasma peripheral; CL2 is
    # the inter-compartmental clearance (Table 2 footnote).
    lcl <- log(1.39); label("Plasma clearance CL1 (L/h)") # Table 2: CL1 = 1.39 L/h (RSE 6%)
    lvc <- log(27.6); label("Plasma central volume V1 (L)") # Table 2: V1 = 27.6 L (RSE 12%)
    lq <- log(13.11); label("Plasma inter-compartmental clearance CL2 (L/h)") # Table 2: CL2 = 13.11 L/h (RSE 18%)
    lvp <- log(17.1); label("Plasma peripheral volume V2 (L)") # Table 2: V2 = 17.1 L (RSE 11%)

    # ELF two-compartment model -- Li 2025 Table 2. Figure 2: compartment 3 =
    # ELF central (receives the nebulized dose), compartment 4 = ELF
    # peripheral; CL4 is the inter-compartmental clearance (Table 2 footnote).
    lcl_elf <- log(0.186); label("ELF central clearance CL3 at ALB = 32.1 g/L (L/h)") # Table 2: CL3 = 0.186 L/h (RSE 21%)
    lv_elf <- log(0.896); label("ELF central volume V3 at AGE = 66 years (L)") # Table 2: V3 = 0.896 L (RSE 19%)
    lq_elf_lung <- log(3.44); label("ELF central-to-peripheral inter-compartmental clearance CL4 (L/h)") # Table 2: CL4 = 3.44 L/h (RSE 14%)
    lv_lung <- log(4.35); label("ELF peripheral volume V4 (L)") # Table 2: V4 = 4.35 L (RSE 10%)

    # Plasma <-> ELF central exchange and nebulization bioavailability.
    lq_elf <- log(0.03); label("Plasma-to-ELF inter-compartmental clearance Q (L/h)") # Table 2: Q = 0.03 L/h (RSE 10%)
    lf_elf <- log(0.557); label("Fraction of a nebulized dose reaching the ELF, F (fraction)") # Table 2: F = 0.557 (RSE 21%)

    # Covariate effects -- Table 2 rows 'AGE on V3' and 'ALB on CL3'; power form
    # on the Table 1 medians (see covariateData notes and the vignette).
    e_age_v_elf <- 1.28; label("Power exponent of (AGE / 66 years) on V3 (unitless)") # Table 2: AGE on V3 = 1.28 (RSE 23%)
    e_alb_cl_elf <- -0.893; label("Power exponent of (ALB / 32.1 g/L) on CL3 (unitless)") # Table 2: ALB on CL3 = -0.893 (RSE 17%)

    # Inter-individual variability -- Table 2 prints omega^2 (variances);
    # shrinkage in parentheses. No IIV on CL2, V2, CL4, V4, Q or F (Discussion).
    etalcl ~ 0.56 # Table 2: omega2 CL1 = 0.56 (shrinkage 22%)
    etalvc ~ 0.31 # Table 2: omega2 V1 = 0.31 (shrinkage 18%)
    etalcl_elf ~ 0.45 # Table 2: omega2 CL3 = 0.45 (shrinkage 25%)
    etalv_elf ~ 0.53 # Table 2: omega2 V3 = 0.53 (shrinkage 14%)

    # Residual variability -- one proportional error for both matrices (Results;
    # Table 2 has a single row). nlmixr2 needs one parameter per endpoint, so
    # the single estimate is carried twice with the same value.
    propSd <- 0.43; label("Proportional residual error, plasma (fraction; shared with ELF in source)") # Table 2: Proportional error = 0.43 (RSE 13%)
    propSd_Celf <- 0.43; label("Proportional residual error, ELF (fraction; shared with plasma in source)") # Table 2: Proportional error = 0.43 (RSE 13%)
  })
  model({
    # Individual parameters.
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q <- exp(lq)
    vp <- exp(lvp)
    cl_elf <- exp(lcl_elf + etalcl_elf) * (ALB / 32.1)^e_alb_cl_elf
    v_elf <- exp(lv_elf + etalv_elf) * (AGE / 66)^e_age_v_elf
    q_elf_lung <- exp(lq_elf_lung)
    v_lung <- exp(lv_lung)
    q_elf <- exp(lq_elf)
    f_elf <- exp(lf_elf)

    # Micro-constants, Figure 2: k10 = CL1/V1, k12 = CL2/V1, k21 = CL2/V2,
    # k13 = Q/V1, k31 = Q/V3, k30 = CL3/V3, k34 = CL4/V3, k43 = CL4/V4.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k_central_elf <- q_elf / vc
    k_elf_central <- q_elf / v_elf
    kel_elf <- cl_elf / v_elf
    k_elf_lung <- q_elf_lung / v_elf
    k_lung_elf <- q_elf_lung / v_lung

    d/dt(central) <- -(kel + k12 + k_central_elf) * central + k21 * peripheral1 + k_elf_central * elf
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(elf) <- -(kel_elf + k_elf_lung + k_elf_central) * elf + k_lung_elf * lung + k_central_elf * central
    d/dt(lung) <- k_elf_lung * elf - k_lung_elf * lung

    # Nebulized dose: only the fraction F reaches the ELF central compartment.
    f(elf) <- f_elf

    Cc <- central / vc
    Celf <- elf / v_elf
    Cc ~ prop(propSd)
    Celf ~ prop(propSd_Celf)
  })
}
