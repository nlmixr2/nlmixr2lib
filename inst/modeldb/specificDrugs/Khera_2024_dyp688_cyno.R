Khera_2024_dyp688_cyno <- function() {
  description <- "Preclinical (cynomolgus monkey). Population PK model of DYP688, a PMEL17-targeted antibody-drug conjugate carrying two FR900359 (Galpha-q/11 inhibitor) payloads that inactivate while still attached to the antibody (Khera 2024). Same structure as the mouse model: a two-compartment linear-elimination model for three ADC species distinguished by their number of ACTIVE payloads (dar2 / dar1 / dar0), each active payload inactivating non-cooperatively at rate kdec in both compartments, with 1% of the dose inactive in the vial. Fitted to total antibody (serum LBA) and conjugated active and inactive payload (plasma papain-release LC-MS assay, which measures payload directly, so no assay cross-reactivity terms). Dosing: give the SAME nmol amount as three simultaneous records into dar2_central, dar1_central and dar0_central; f() splits it into 1 - fr, (1 - xfr) * fr and xfr * fr, reproducing the authors' three Monolix depot() targets."
  reference <- paste(
    "Khera E, Dharmarajan L, Hainzl D, Engelhardt V, Vostiarova H, Davis J,",
    "Ebel N, Wuersch K, Romanet V, Sharaby S, Kearns JD. (2024). QSP modeling",
    "of a transiently inactivating antibody-drug conjugate highlights benefit",
    "of short antibody half life. J Pharmacokinet Pharmacodyn 52(1):7.",
    "doi:10.1007/s10928-024-09956-1."
  )
  vignette <- "Khera_2024_dyp688"
  units <- list(
    time = "h",
    dosing = "nmol (antibody; one record per ADC species, see description)",
    concentration = "nM (nmol / L)"
  )

  compartmentData <- list(
    dar2_central = list(
      analyte = "DYP688 with two active FR900359 payloads",
      units = "nmol",
      specimen = "plasma",
      verified = TRUE
    ),
    dar1_central = list(
      analyte = "DYP688 with one active FR900359 payload",
      units = "nmol",
      specimen = "plasma",
      verified = TRUE
    ),
    dar0_central = list(
      analyte = "DYP688 with no active FR900359 payload",
      units = "nmol",
      specimen = "plasma",
      verified = TRUE
    ),
    dar2_peripheral1 = list(
      analyte = "DYP688 with two active FR900359 payloads",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    dar1_peripheral1 = list(
      analyte = "DYP688 with one active FR900359 payload",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    dar0_peripheral1 = list(
      analyte = "DYP688 with no active FR900359 payload",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    )
  )

  population <- list(
    species = "cynomolgus monkey",
    n_subjects = 20L,
    n_studies = 1L,
    weight_median = "3 kg (assumed body weight used for the allometric scaling to human)",
    disease_state = "Healthy animals in a GLP toxicology study",
    dose_range = "30 and 75 mg/kg slow IV bolus every two weeks up to Day 71 (six administrations)",
    regions = "Novartis (Switzerland)",
    notes = "Khera 2024 Methods 'Preclinical in vivo study design' (the text describes male animals and also compares males and females). Intensive sampling after the first (Day 1) and sixth (Day 71) doses: pre-dose and 1, 6, 24, 72, 168 and 336 h. Monolix 2020R1 output in the Supplementary Material: 20 animals, 280 observations each of conjugated active payload, conjugated inactive payload and total antibody, 120 doses. Total antibody by sequential sandwich ELISA in serum (LLOQ 0.144 ug/mL); conjugated active and inactive payload by papain-release LC-MS/MS assay in plasma."
  )

  ini({
    # Structural parameters: Khera 2024 Table 2, as printed to full
    # precision in the Supplementary Material Monolix output ('NHP PK model
    # results summary'). Log-normally distributed (Monolix default).
    lkel <- fixed(log(0.0129))
    label("Antibody / ADC elimination rate constant kel (1/h)") # Supplement Monolix output k_pop = 0.0129 (no SE); Table 2 0.013 (RSE '-')
    lk12 <- log(0.0637)
    label("Central-to-peripheral distribution rate constant k12 (1/h)") # Supplement Monolix output k12_pop = 0.0637; Table 2 0.064 (RSE 0.288%)
    lk21 <- log(0.0386)
    label("Peripheral-to-central distribution rate constant k21 (1/h)") # Supplement Monolix output k21_pop = 0.0386; Table 2 0.039 (RSE 11.8%)
    lkdec <- log(0.0114)
    label("On-antibody inactivation rate constant of each active payload ki.a. (1/h)") # Supplement Monolix output kia1_pop = 0.0114; Table 2 0.0114 (RSE 1.62%)
    lvc <- log(0.0829)
    label("Plasma volume of distribution Vc (L)") # Supplement Monolix output V_pop = 0.0829; Table 2 0.083 (RSE 4.83%)

    fr <- fixed(0.01)
    label("Fraction of inactive ADC in the dose fr (fraction)") # Table 2 footnote: 'Fraction of inactive ADC in dose (fr) vial assumed to 1%'; Monolix output fr_pop = 0.01 (no SE)
    xfr <- fixed(0.99)
    label("Fraction of the inactive dose fraction carried as ADC with no active payload xfr (fraction)") # Supplement Monolix output xfr_pop = 0.99 (no SE)

    # IIV (Monolix SD of the log-scale random effect; variances = omega^2).
    # omega_k21 and omega_V carry no standard error (fixed).
    etalk21 ~ fixed(0.2304) # Supplement Monolix output omega_k21 = 0.48 (no SE); Table 2 IISD 0.48
    etalkdec ~ 0.00389376 # Supplement Monolix output omega_kia1 = 0.0624 (RSE 22.1%); Table 2 IISD 0.06
    etalvc ~ fixed(0.042436) # Supplement Monolix output omega_V = 0.206 (no SE); Table 2 IISD 0.21

    # Residual error: Monolix proportional error b1, b2, b3 on the three
    # outputs in OUTPUT-block order {Total_mAb, Active_P, Inactive_P}.
    propSd <- fixed(0.37)
    label("Proportional residual SD of total antibody (fraction)") # Supplement Monolix output b1 = 0.37 (no SE); Table 2 0.37
    propSd_Cc_active_payload <- 0.267
    label("Proportional residual SD of conjugated active payload (fraction)") # Supplement Monolix output b2 = 0.267 (RSE 5.51%); Table 2 0.27
    propSd_Cc_inactive_payload <- 0.725
    label("Proportional residual SD of conjugated inactive payload (fraction)") # Supplement Monolix output b3 = 0.725 (RSE 4.26%); Table 2 0.73
  })

  model({
    kel <- exp(lkel)
    k12 <- exp(lk12)
    k21 <- exp(lk21 + etalk21)
    kdec <- exp(lkdec + etalkdec)
    vc <- exp(lvc + etalvc)

    # Supplement 'NHP PK model structure file' (Monolix EQUATION block),
    # written on amounts; identical to the mouse structure.
    d/dt(dar2_central) <- -k12 * dar2_central + k21 * dar2_peripheral1 - 2 * kdec * dar2_central - kel * dar2_central
    d/dt(dar1_central) <- -k12 * dar1_central + k21 * dar1_peripheral1 + 2 * kdec * dar2_central - kdec * dar1_central - kel * dar1_central
    d/dt(dar0_central) <- -k12 * dar0_central + k21 * dar0_peripheral1 + kdec * dar1_central - kel * dar0_central
    d/dt(dar2_peripheral1) <- k12 * dar2_central - k21 * dar2_peripheral1 - 2 * kdec * dar2_peripheral1
    d/dt(dar1_peripheral1) <- k12 * dar1_central - k21 * dar1_peripheral1 + 2 * kdec * dar2_peripheral1 - kdec * dar1_peripheral1
    d/dt(dar0_peripheral1) <- k12 * dar0_central - k21 * dar0_peripheral1 + kdec * dar1_peripheral1

    f(dar2_central) <- 1 - fr
    f(dar1_central) <- (1 - xfr) * fr
    f(dar0_central) <- xfr * fr

    # Concentrations in nM: states in nmol, vc in L.
    c_dar2 <- dar2_central / vc
    c_dar1 <- dar1_central / vc
    c_dar0 <- dar0_central / vc

    # Observables (Monolix OUTPUT block): total antibody, conjugated active
    # payload (2 * ADC2P + ADC1P) and conjugated inactive payload
    # (total payload - active payload), all nM.
    Cc <- c_dar2 + c_dar1 + c_dar0
    Cc_active_payload <- 2 * c_dar2 + c_dar1
    Cc_inactive_payload <- 2 * Cc - Cc_active_payload

    Cc ~ prop(propSd)
    Cc_active_payload ~ prop(propSd_Cc_active_payload)
    Cc_inactive_payload ~ prop(propSd_Cc_inactive_payload)
  })
}
