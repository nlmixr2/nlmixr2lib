Khera_2024_dyp688_human <- function() {
  description <- "Human PK projection of DYP688, a PMEL17-targeted antibody-drug conjugate carrying two FR900359 (Galpha-q/11 inhibitor) payloads that inactivate while still attached to the antibody (Khera 2024), scaled allometrically from the cynomolgus monkey PK model (single species, body-weight exponent 1, 3 kg monkey to 70 kg human; Khera 2024 Table 5). Two-compartment linear-elimination model for three ADC species distinguished by their number of ACTIVE payloads (dar2 / dar1 / dar0), each active payload inactivating non-cooperatively at rate kdec in both compartments; the inactivation rate, distribution rate constants and 1% inactive dose fraction are carried over unchanged from the monkey. Predicts a typical 70 kg patient; the authors compared it with interim Phase I PK at 8 mg/kg Q2W (Khera 2024 Figure 6). Dosing: give the SAME nmol amount as three simultaneous records into dar2_central, dar1_central and dar0_central; f() splits it into 1 - fr, (1 - xfr) * fr and xfr * fr. For an infusion use dur (not rate) so the split scales the infused amount rather than the duration."
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
    species = "human (projected from cynomolgus monkey)",
    n_subjects = NA_integer_,
    weight_median = "70 kg (assumed for the allometric scaling)",
    disease_state = "Metastatic uveal melanoma and other GNAQ/11-mutant melanomas (target population of the first-in-human study NCT05415072)",
    dose_range = "Projection evaluated at 8 mg/kg IV over 1.5 h every 2 weeks (Figure 6); the Phase I study dosed 4, 8, 12 or 16 mg/kg Q2W in 4-week cycles",
    notes = "No human data were used to build this model: every parameter is scaled from the cynomolgus monkey fit (Khera 2024 Table 5). The authors overlaid the projection on interim clinical PK from 12 patients at 8 mg/kg (data cut-off 8 Sept 2023; Figure 6)."
  )

  ini({
    # Khera 2024 Table 5 'Predicted Human PK Median', from the cynomolgus
    # medians by V_human = V_cyno * (BW_human / BW_cyno)^w and
    # CL_human = CL_cyno * (BW_human / BW_cyno)^w with w = 1, BW_human = 70 kg
    # and BW_cyno = 3 kg. With the same exponent on CL and V the elimination
    # rate constant CL / V is unchanged by the scaling.
    lvc <- fixed(log(1.93))
    label("Plasma volume of distribution Vc (L)") # Table 5, predicted human Vc = 1.93 L (cyno 0.083 L * 70 / 3)
    lkel <- fixed(log(0.0251 / 1.93))
    label("Antibody / ADC elimination rate constant kel = CL / Vc (1/h)") # Table 5, predicted human CL = 0.0251 L/h divided by Vc = 1.93 L
    lk12 <- fixed(log(0.064))
    label("Central-to-peripheral distribution rate constant k12 (1/h)") # Table 5, k12 = 0.064 'Assumed same in cyno and human'
    lk21 <- fixed(log(0.039))
    label("Peripheral-to-central distribution rate constant k21 (1/h)") # Table 5, k21 = 0.039 'Assumed same in cyno and human'
    lkdec <- fixed(log(log(2) / 61))
    label("On-antibody inactivation rate constant of each active payload ki.a. (1/h)") # Table 5, t1/2,inactivation = 61 h 'Assumed same in cyno and human'; ki.a. = ln 2 / 61
    fr <- fixed(0.01)
    label("Fraction of inactive ADC in the dose fr (fraction)") # Table 5, initial fraction of inactive ADC = 1%
    xfr <- fixed(0.99)
    label("Fraction of the inactive dose fraction carried as ADC with no active payload xfr (fraction)") # cyno model (Supplement Monolix output xfr_pop = 0.99), carried over unchanged

    # Variability is carried over from the cynomolgus model: Table 5 gives
    # the human IQRs as the cyno IQRs scaled by the same factor, i.e. the
    # monkey between-animal variability on V, k21 and ki.a. (CL inherits the
    # V variability because kel has none).
    etalvc ~ fixed(0.042436) # cyno model omega_V = 0.206 (Supplement Monolix output); Table 5 IQR 1.75-2.21 L
    etalk21 ~ fixed(0.2304) # cyno model omega_k21 = 0.48 (Supplement Monolix output); Table 5 IQR 0.03-0.046 1/h
    etalkdec ~ fixed(0.00389376) # cyno model omega_kia1 = 0.0624 (Supplement Monolix output); Table 5 IQR 59-63 h

    # Residual error carried over from the cynomolgus model (same assay
    # types: LBA total antibody, papain-release LC-MS conjugated payload).
    propSd <- fixed(0.37)
    label("Proportional residual SD of total antibody (fraction)") # cyno model b1 = 0.37 (Supplement Monolix output; Table 2)
    propSd_Cc_active_payload <- fixed(0.267)
    label("Proportional residual SD of conjugated active payload (fraction)") # cyno model b2 = 0.267 (Supplement Monolix output; Table 2)
    propSd_Cc_inactive_payload <- fixed(0.725)
    label("Proportional residual SD of conjugated inactive payload (fraction)") # cyno model b3 = 0.725 (Supplement Monolix output; Table 2)
  })

  model({
    vc <- exp(lvc + etalvc)
    kel <- exp(lkel)
    k12 <- exp(lk12)
    k21 <- exp(lk21 + etalk21)
    kdec <- exp(lkdec + etalkdec)

    # Same structure as the cynomolgus model (Supplement 'NHP PK model
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

    c_dar2 <- dar2_central / vc
    c_dar1 <- dar1_central / vc
    c_dar0 <- dar0_central / vc

    Cc <- c_dar2 + c_dar1 + c_dar0
    Cc_active_payload <- 2 * c_dar2 + c_dar1
    Cc_inactive_payload <- 2 * Cc - Cc_active_payload

    Cc ~ prop(propSd)
    Cc_active_payload ~ prop(propSd_Cc_active_payload)
    Cc_inactive_payload ~ prop(propSd_Cc_inactive_payload)
  })
}
