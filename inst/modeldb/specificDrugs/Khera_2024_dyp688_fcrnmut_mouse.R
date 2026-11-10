Khera_2024_dyp688_fcrnmut_mouse <- function() {
  description <- "Preclinical (mouse). Population PK model of DYP688 FcRn-mut, a variant of the PMEL17-targeted antibody-drug conjugate DYP688 with an Fc mutation that abolishes FcRn binding and so clears rapidly (Khera 2024 Supplementary Table 2). Same structure as the DYP688 mouse model: a two-compartment linear-elimination model for three ADC species distinguished by their number of ACTIVE FR900359 payloads (dar2 / dar1 / dar0), each active payload inactivating non-cooperatively at rate kdec in both compartments, with a fraction fr of the dose inactive in the vial and ligand-binding-assay cross-reactivity terms pinact / pinact2. Distribution and inactivation rate constants are fixed to the DYP688 values; volume, elimination rate and fr were re-estimated. Used by the authors with the DYP688 tumour growth inhibition model to predict FcRn-mut efficacy (Figure S7B). Dosing: give the SAME nmol amount as three simultaneous records into dar2_central, dar1_central and dar0_central; f() splits it into 1 - fr, (1 - xfr) * fr and xfr * fr."
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
    concentration = "nM (antibody concentration; 1000 * nmol / mL)"
  )

  compartmentData <- list(
    dar2_central = list(
      analyte = "DYP688 FcRn-mut with two active FR900359 payloads",
      units = "nmol",
      specimen = "serum",
      verified = TRUE
    ),
    dar1_central = list(
      analyte = "DYP688 FcRn-mut with one active FR900359 payload",
      units = "nmol",
      specimen = "serum",
      verified = TRUE
    ),
    dar0_central = list(
      analyte = "DYP688 FcRn-mut with no active FR900359 payload",
      units = "nmol",
      specimen = "serum",
      verified = TRUE
    ),
    dar2_peripheral1 = list(
      analyte = "DYP688 FcRn-mut with two active FR900359 payloads",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    dar1_peripheral1 = list(
      analyte = "DYP688 FcRn-mut with one active FR900359 payload",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    ),
    dar0_peripheral1 = list(
      analyte = "DYP688 FcRn-mut with no active FR900359 payload",
      units = "nmol",
      specimen = "tissue",
      verified = TRUE
    )
  )

  population <- list(
    species = "mouse (female nude Crl:NU(NCr)-Foxn1nu, 92.1 uveal melanoma xenograft)",
    n_subjects = NA_integer_,
    sex_female_pct = 100,
    disease_state = "Subcutaneous 92.1 uveal melanoma cell-derived xenograft",
    dose_range = "IV single doses (3, 7.5 and 12.5 mg/kg FcRn-mut appear in the efficacy comparisons, Figure 3B and Figure S9B)",
    notes = "Khera 2024 Results 'PK-TGI modeling of mouse tumors' and Supplementary Figure S6 / Supplementary Table 2. The number of animals in the FcRn-mut PK fit is not reported."
  )

  ini({
    # Supplementary Table 2 'FcRn-mut PK fitting parameters'. Values printed
    # with an RSE were estimated; those without one equal the DYP688 mouse
    # values and were held constant.
    lvc <- log(1.04)
    label("Plasma volume of distribution Vc (mL)") # Supplementary Table 2, Vc = 1.04 mL (RSE 14.3%)
    lkel <- log(0.2709)
    label("Antibody / ADC elimination rate constant kel (1/h)") # Supplementary Table 2, kel = 0.2709 1/h (RSE 12.5%)
    lk12 <- fixed(log(0.035))
    label("Central-to-peripheral distribution rate constant k12 (1/h)") # Supplementary Table 2, k12 = 0.035 (RSE '-')
    lk21 <- fixed(log(0.036))
    label("Peripheral-to-central distribution rate constant k21 (1/h)") # Supplementary Table 2, k21 = 0.036 (RSE '-')
    lkdec <- fixed(log(0.013))
    label("On-antibody inactivation rate constant of each active payload ki.a. (1/h)") # Supplementary Table 2, kia = 0.013 (RSE '-')
    fr <- 0.0082
    label("Fraction of inactive ADC in the dose fr (fraction)") # Supplementary Table 2, fr = 0.0082 (RSE 1.67%)
    xfr <- fixed(0.99)
    label("Fraction of the inactive dose fraction carried as ADC with no active payload xfr (fraction)") # Supplementary Table 2, xfr = 0.99 (RSE '-')
    pinact <- fixed(0.83)
    label("Fraction of ADC with one active payload detected by the active-ADC LBA pinact (fraction)") # Supplementary Table 2, pinact = 0.83 (RSE '-')
    pinact2 <- fixed(0.18)
    label("Fraction of ADC with one active payload detected by the inactive-ADC LBA pinact2 (fraction)") # Supplementary Table 2, pinact2 = 0.18 (RSE '-')

    # IIV: Supplementary Table 2 'IISD' column (SD of the log-scale random
    # effect); variances = IISD^2. Entries without an RSE are the DYP688
    # mouse values held constant.
    etalvc ~ 0.0016 # Supplementary Table 2, Vc IISD = 0.04 (RSE 251%)
    etalkel ~ 0.0441 # Supplementary Table 2, kel IISD = 0.21 (RSE 32%)
    etalk12 ~ fixed(0.0225) # Supplementary Table 2, k12 IISD = 0.15 (RSE '-')
    etalk21 ~ fixed(0.0676) # Supplementary Table 2, k21 IISD = 0.26 (RSE '-')
    etalkdec ~ fixed(0.005625) # Supplementary Table 2, kia IISD = 0.075 (RSE '-')

    # Residual error is not reported for the FcRn-mut fit; the DYP688 mouse
    # values (same model, same assays) are carried over unchanged.
    expSd <- fixed(0.518)
    label("Log-scale residual SD of total antibody (unitless)") # DYP688 mouse model, Supplement Monolix output aTOTAL_MAB = 0.518
    expSd_Cc_active <- fixed(0.409)
    label("Log-scale residual SD of active ADC (unitless)") # DYP688 mouse model, Supplement Monolix output aTOTAL_ADC_ACTIVE_HOMOGENOUS = 0.409
    expSd_Cc_inactive <- fixed(0.244)
    label("Log-scale residual SD of inactive ADC (unitless)") # DYP688 mouse model, Supplement Monolix output aTOTAL_ADC_INACTIVE_SEQUENTIAL = 0.244
  })

  model({
    vc <- exp(lvc + etalvc)
    kel <- exp(lkel + etalkel)
    k12 <- exp(lk12 + etalk12)
    k21 <- exp(lk21 + etalk21)
    kdec <- exp(lkdec + etalkdec)

    # Same structure as the DYP688 mouse model (Supplement 'Mouse PK model
    # structure file').
    d/dt(dar2_central) <- -k12 * dar2_central + k21 * dar2_peripheral1 - 2 * kdec * dar2_central - kel * dar2_central
    d/dt(dar1_central) <- -k12 * dar1_central + k21 * dar1_peripheral1 + 2 * kdec * dar2_central - kdec * dar1_central - kel * dar1_central
    d/dt(dar0_central) <- -k12 * dar0_central + k21 * dar0_peripheral1 + kdec * dar1_central - kel * dar0_central
    d/dt(dar2_peripheral1) <- k12 * dar2_central - k21 * dar2_peripheral1 - 2 * kdec * dar2_peripheral1
    d/dt(dar1_peripheral1) <- k12 * dar1_central - k21 * dar1_peripheral1 + 2 * kdec * dar2_peripheral1 - kdec * dar1_peripheral1
    d/dt(dar0_peripheral1) <- k12 * dar0_central - k21 * dar0_peripheral1 + kdec * dar1_peripheral1

    f(dar2_central) <- 1 - fr
    f(dar1_central) <- (1 - xfr) * fr
    f(dar0_central) <- xfr * fr

    c_dar2 <- 1000 * dar2_central / vc
    c_dar1 <- 1000 * dar1_central / vc
    c_dar0 <- 1000 * dar0_central / vc

    # Conjugated active payload (2 * ADC2P + ADC1P), nM: the exposure that
    # drives the tumour growth inhibition model.
    active_payload <- 2 * c_dar2 + c_dar1

    Cc <- c_dar2 + c_dar1 + c_dar0
    Cc_active <- c_dar2 + pinact * c_dar1
    Cc_inactive <- c_dar0 + pinact2 * c_dar1

    Cc ~ lnorm(expSd)
    Cc_active ~ lnorm(expSd_Cc_active)
    Cc_inactive ~ lnorm(expSd_Cc_inactive)
  })
}
