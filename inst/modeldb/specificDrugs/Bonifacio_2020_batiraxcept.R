Bonifacio_2020_batiraxcept <- function() {
  description <- paste(
    "Preclinical (cynomolgus monkey), allometrically scaled to human.",
    "Two-compartment population PK model for batiraxcept (AVB-S6-500, an",
    "AXL-ectodomain / IgG1 Fc fusion protein that neutralises GAS6) with",
    "parallel linear and Michaelis-Menten (target-mediated) elimination from",
    "the central compartment, fit to cynomolgus monkey serum concentrations",
    "(0.1-150 mg/kg) with parameters centred on a 70 kg subject. CL, Q and",
    "Vmax scale with (WT/70)^0.75 and VC, VP with (WT/70)^1. A direct",
    "(no-hysteresis) inhibitory sigmoid Emax relationship drives serum free",
    "GAS6 suppression from the drug concentration (not scaled between",
    "species). The paper used the typical-value model to select the",
    "first-in-human doses (1, 2.5, 5, 10 mg/kg IV) in healthy volunteers."
  )
  reference <- paste(
    "Bonifacio L, Dodds M, Prohaska D, Moss A, Giaccia A, Tabibiazar R,",
    "McIntyre G. (2020). Target-Mediated Drug Disposition",
    "Pharmacokinetic/Pharmacodynamic Model-Informed Dose Selection for the",
    "First-in-Human Study of AVB-S6-500. Clinical and Translational Science",
    "13(1), 204-211. doi:10.1111/cts.12706. Model equations and the PD",
    "parameter table are in Supplementary Material CTS-13-204-s001.docx.",
    sep = " "
  )
  vignette <- "Bonifacio_2020_batiraxcept"
  units <- list(
    time = "day",
    dosing = "mg",
    concentration = "ug/mL"
  )

  covariateData <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Individual body weight. Allometric scaling centred on 70 kg:",
        "exponent 0.75 on CL, Q and Vmax, exponent 1.0 on VC and VP",
        "(Bonifacio 2020 Methods 'Model predictions of clinical PK/PD' and",
        "Supplement 'Scaling to Human Pharmacokinetics'). The same scaling",
        "was used within the monkey fit (individual animal weight was in the",
        "analysis dataset) and for the human projection; human simulations",
        "in the paper used a 75 kg subject."
      ),
      source_name = "WT"
    )
  )

  compartmentData <- list(
    central = list(analyte = "batiraxcept", units = "mg", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "batiraxcept", units = "mg", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "cynomolgus monkey (fit); allometrically scaled to human",
    n_subjects = NA_integer_,
    n_studies = 5L,
    disease_state = "healthy (nonclinical toxicology / PK studies)",
    dose_range = paste(
      "0.1, 0.5, 1 mg/kg (2 animals each), 5 mg/kg (4 animals, 2 studies),",
      "30-100 mg/kg (12 animals per dose level, 3 studies), 150 mg/kg",
      "(18 animals)"
    ),
    regions = "Canada (CiToxLAB) and USA (Shin Nippon Biomedical Laboratories)",
    notes = paste(
      "Five nonclinical cynomolgus monkey studies (Bonifacio 2020 Methods",
      "'Nonclinical studies used for PK/PD model development'). Total",
      "animal count not stated (the per-dose counts above are from the",
      "Methods). Animal body weights not reported. The model was validated",
      "externally against a first-in-human study (NCT03401528) in 31",
      "healthy volunteers receiving batiraxcept (single 1, 2.5, 5, 10 mg/kg",
      "or 5 mg/kg weekly x 4, 60-min IV infusion; median weight 75.5 kg;",
      "Table 1); those human data were not used in the fit."
    )
  )

  ini({
    # Structural PK parameters, centred on a 70 kg subject
    # (Bonifacio 2020 Table 2 = Supplement Table 1).
    lcl <- log(1.10); label("Linear clearance at 70 kg (L/day)") # Table 2: CL 1.10 L/day (RSE 7.87)
    lvc <- log(2.97); label("Central volume at 70 kg (L)") # Table 2: VC 2.97 L (RSE 2.90)
    lvp <- log(2.94); label("Peripheral volume at 70 kg (L)") # Table 2: VP 2.94 L (RSE 6.29)
    lq <- log(2.04); label("Intercompartmental clearance at 70 kg (L/day)") # Table 2: Q 2.04 L/day (RSE 12.0)
    lvmax <- log(3.35); label("Maximal target-mediated elimination rate at 70 kg (mg/day)") # Table 2: Vmax 3.35 mg/day (RSE 3.77)
    lkm <- log(0.102); label("Concentration at half-maximal target-mediated elimination (ug/mL)") # Table 2: KM 102 ng/mL = 0.102 ug/mL (RSE 7.01)

    # Allometric exponents held at the stated values (not estimated)
    e_wt_cl <- fixed(0.75); label("Allometric exponent on CL, Q and Vmax (unitless)") # Methods 'Model predictions of clinical PK/PD'; Supplement 'Scaling to Human Pharmacokinetics'
    e_wt_vc <- fixed(1); label("Allometric exponent on VC and VP (unitless)") # Methods 'Model predictions of clinical PK/PD'; Supplement 'Scaling to Human Pharmacokinetics'

    # PD: direct inhibitory sigmoid Emax on serum GAS6 (Supplement Table 2,
    # fit to time-matched monkey AVB-S6-500 / GAS6 pairs; not scaled to human)
    lrbase <- log(23.8); label("Baseline serum GAS6 (E0) (ng/mL)") # Supplement Table 2: E0 23.800 (SE 0.749)
    lec50 <- log(0.0214); label("Batiraxcept concentration for half-maximal GAS6 suppression (ug/mL)") # Supplement Table 2: EC50 21.400 ng/mL = 0.0214 ug/mL (SE 10.300)
    lhill <- log(0.796); label("Hill coefficient of GAS6 suppression (unitless)") # Supplement Table 2: H 0.796 (SE 0.250)

    # Between-animal variability, log-normal; reported as CV and converted
    # with omega^2 = log(CV^2 + 1). Covariances were estimated but not
    # reported, so a diagonal omega is used.
    etalcl ~ 0.08453 # Table 2: BSV CL 29.7 CV (RSE 1.45)
    etalvc ~ 0.01941 # Table 2: BSV VC 14.0 CV (RSE 16.1)
    etalvp ~ 0.02024 # Table 2: BSV VP 14.3 CV (RSE 15.7)
    etalq ~ 0.22314 # Table 2: BSV Q 50.0 CV (RSE 27.4)
    etalvmax ~ 0.67909 # Table 2: BSV Vmax 98.6 CV (RSE 35.1)
    etalkm ~ 0.72313 # Table 2: BSV KM 103 CV (RSE 0.874)

    # Residual error: proportional per the Supplement, magnitude not reported
    propSd <- fixed(0); label("Proportional residual error (fraction)") # Supplement 'Non-Human Primate Pharmacokinetic Analysis': proportional; value not reported
  })

  model({
    # Individual parameters with allometric scaling centred on 70 kg
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl
    q <- exp(lq + etalq) * (WT / 70)^e_wt_cl
    vmax <- exp(lvmax + etalvmax) * (WT / 70)^e_wt_cl
    vc <- exp(lvc + etalvc) * (WT / 70)^e_wt_vc
    vp <- exp(lvp + etalvp) * (WT / 70)^e_wt_vc
    km <- exp(lkm + etalkm)

    rbase <- exp(lrbase)
    ec50 <- exp(lec50)
    hill <- exp(lhill)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    Cc <- central / vc

    # Supplement Methods: dAp/dt and dAt/dt with parallel linear and
    # Michaelis-Menten elimination of the central amount
    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1 - vmax * Cc / (km + Cc)
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # Supplement 'Non-Human Primate Pharmacodynamic Analysis':
    # GAS6 = E0 * (1 - C^H / (EC50^H + C^H)); serum GAS6 in ng/mL
    GAS6 <- rbase * (1 - Cc^hill / (ec50^hill + Cc^hill))

    Cc ~ prop(propSd)
  })
}
