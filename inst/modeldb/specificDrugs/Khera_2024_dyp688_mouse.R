Khera_2024_dyp688_mouse <- function() {
  description <- "Preclinical (mouse). Population PK model of DYP688, a PMEL17-targeted antibody-drug conjugate carrying two FR900359 (Galpha-q/11 inhibitor) payloads that inactivate while still attached to the antibody (Khera 2024). Two-compartment linear-elimination model for three ADC species distinguished by their number of ACTIVE payloads (dar2 = two active, dar1 = one active, dar0 = none active); each active payload inactivates non-cooperatively at rate kdec in both the central and peripheral compartments (dar2 -> dar1 at 2 * kdec, dar1 -> dar0 at kdec). A small fraction fr of the dose is inactive ADC in the vial. Fitted to total antibody, active ADC and inactive ADC (ligand binding assays, LBA) and percent active conjugated payload in mouse serum. Dosing: give the SAME nmol amount as three simultaneous records into dar2_central, dar1_central and dar0_central; f() splits it into 1 - fr, (1 - xfr) * fr and xfr * fr, reproducing the authors' three Monolix depot() targets."
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
      analyte = "DYP688 with two active FR900359 payloads",
      units = "nmol",
      specimen = "serum",
      verified = TRUE
    ),
    dar1_central = list(
      analyte = "DYP688 with one active FR900359 payload",
      units = "nmol",
      specimen = "serum",
      verified = TRUE
    ),
    dar0_central = list(
      analyte = "DYP688 with no active FR900359 payload",
      units = "nmol",
      specimen = "serum",
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
    species = "mouse (female nude Crl:NU(NCr)-Foxn1nu, 92.1-luc uveal melanoma xenograft)",
    n_subjects = 16L,
    n_studies = 1L,
    age_range = "6-8 weeks at study start",
    sex_female_pct = 100,
    disease_state = "Subcutaneous 92.1-luc uveal melanoma cell-derived xenograft",
    dose_range = "IV 3 and 6 mg/kg single dose; dose fractionation 1.5 mg/kg Q2W vs 3 mg/kg Q4W and 3 mg/kg Q2W vs 6 mg/kg Q4W (six treatment groups)",
    regions = "Novartis, Basel (Switzerland)",
    notes = "Khera 2024 Methods 'Preclinical in vivo study design' and Results 'Modeling with payload inactivation captures mouse and NHP PK'. Monolix 2020R1 output in the Supplementary Material: 16 animals, 108 active-ADC, 126 inactive-ADC, 134 total-antibody and 131 percent-conjugated-payload observations, 34 doses. Total antibody and total inactive ADC by sequential sandwich ELISA, total active ADC by homogeneous sandwich ELISA."
  )

  ini({
    # Structural parameters: Khera 2024 Table 1, final estimates as printed
    # to full precision in the Supplementary Material Monolix output
    # ('Mouse PK model results summary'). All are log-normally distributed
    # (Monolix default), so each is log-transformed here.
    lkel <- log(0.0161)
    label("Antibody / ADC elimination rate constant kel (1/h)") # Supplement Monolix output kel_pop = 0.0161; Table 1 0.016 (RSE 12.1%)
    lk12 <- log(0.0355)
    label("Central-to-peripheral distribution rate constant k12 (1/h)") # Supplement Monolix output k12_pop = 0.0355; Table 1 0.035 (RSE 11.4%)
    lk21 <- log(0.0356)
    label("Peripheral-to-central distribution rate constant k21 (1/h)") # Supplement Monolix output k21_pop = 0.0356; Table 1 0.036 (RSE 14.1%)
    lkdec <- log(0.0125)
    label("On-antibody inactivation rate constant of each active payload ki.a. (1/h)") # Supplement Monolix output kia1_pop = 0.0125; Table 1 0.013 (RSE 2.77%)
    lvc <- log(1.53)
    label("Plasma volume of distribution Vc (mL)") # Supplement Monolix output Vc_pop = 1.53; Table 1 1.53 (RSE 2.5%)

    # Dose composition and assay cross-reactivity (no IIV).
    fr <- 0.0274
    label("Fraction of inactive ADC in the dose fr (fraction)") # Supplement Monolix output fr_pop = 0.0274; Table 1 0.027 (RSE 0.46%)
    xfr <- fixed(0.99)
    label("Fraction of the inactive dose fraction carried as ADC with no active payload xfr (fraction)") # Supplement Monolix output xfr_pop = 0.99 (no SE); Methods: '99% of the inactive ADC in the dose was assumed to be in the ADC0P form'
    pinact <- 0.793
    label("Fraction of ADC with one active payload detected by the active-ADC LBA pinact (fraction)") # Supplement Monolix output pinact_pop = 0.793; Table 1 0.79 (RSE 9.33%)
    pinact2 <- 0.161
    label("Fraction of ADC with one active payload detected by the inactive-ADC LBA pinact2 (fraction)") # Supplement Monolix output pinact2_pop = 0.161; Table 1 0.16 (RSE 24.5%)

    # IIV: Monolix reports the SD of the random effect on the log scale;
    # the variances below are omega^2. omega_k12, omega_k21 and omega_Vc
    # carry no standard error in the Monolix output (fixed).
    etalkel ~ 0.201601 # Supplement Monolix output omega_kel = 0.449 (RSE 18.2%); Table 1 IISD 0.45
    etalk12 ~ fixed(0.0225) # Supplement Monolix output omega_k12 = 0.15 (no SE)
    etalk21 ~ fixed(0.0676) # Supplement Monolix output omega_k21 = 0.26 (no SE)
    etalkdec ~ 0.00567009 # Supplement Monolix output omega_kia1 = 0.0753 (RSE 28.6%); Table 1 IISD 0.075
    etalvc ~ fixed(0.0009) # Supplement Monolix output omega_Vc = 0.03 (no SE)

    # Residual error. The three ELISA outputs use a Monolix constant error on
    # the log scale (lognormal). The assignment of the three values follows
    # the observation-ID names in the Monolix output; Table 1 prints the same
    # three numbers against permuted row labels (see vignette).
    expSd <- 0.518
    label("Log-scale residual SD of total antibody (unitless)") # Supplement Monolix output aTOTAL_MAB = 0.518 (RSE 8.24%)
    expSd_Cc_active <- 0.409
    label("Log-scale residual SD of active ADC (unitless)") # Supplement Monolix output aTOTAL_ADC_ACTIVE_HOMOGENOUS = 0.409 (RSE 8.67%)
    expSd_Cc_inactive <- 0.244
    label("Log-scale residual SD of inactive ADC (unitless)") # Supplement Monolix output aTOTAL_ADC_INACTIVE_SEQUENTIAL = 0.244 (RSE 7.73%)
    addSd_active_payload_pct <- 0.654
    label("Additive residual SD of percent active conjugated payload (percentage points)") # Supplement Monolix output ap_ACTIVE_CONJUGATED_PAYLOAD_HOMOGENOUS = 0.654; Table 1 0.65 (RSE 28.5%)
    propSd_active_payload_pct <- 0.209
    label("Proportional residual SD of percent active conjugated payload (fraction)") # Supplement Monolix output bp_ACTIVE_CONJUGATED_PAYLOAD_HOMOGENOUS = 0.209; Table 1 0.21 (RSE 10.3%)
  })

  model({
    kel <- exp(lkel + etalkel)
    k12 <- exp(lk12 + etalk12)
    k21 <- exp(lk21 + etalk21)
    kdec <- exp(lkdec + etalkdec)
    vc <- exp(lvc + etalvc)

    # Supplement 'Mouse PK model structure file' (Monolix EQUATION block),
    # written on amounts. Each active payload inactivates independently, so
    # an ADC with two active payloads loses one at 2 * kdec (kia2 = 2 *
    # kia1). Elimination is from the central compartment only and is the
    # same for all three species (kel2 = kel).
    d/dt(dar2_central) <- -k12 * dar2_central + k21 * dar2_peripheral1 - 2 * kdec * dar2_central - kel * dar2_central
    d/dt(dar1_central) <- -k12 * dar1_central + k21 * dar1_peripheral1 + 2 * kdec * dar2_central - kdec * dar1_central - kel * dar1_central
    d/dt(dar0_central) <- -k12 * dar0_central + k21 * dar0_peripheral1 + kdec * dar1_central - kel * dar0_central
    d/dt(dar2_peripheral1) <- k12 * dar2_central - k21 * dar2_peripheral1 - 2 * kdec * dar2_peripheral1
    d/dt(dar1_peripheral1) <- k12 * dar1_central - k21 * dar1_peripheral1 + 2 * kdec * dar2_peripheral1 - kdec * dar1_peripheral1
    d/dt(dar0_peripheral1) <- k12 * dar0_central - k21 * dar0_peripheral1 + kdec * dar1_peripheral1

    # Dose composition (Monolix PK block: depot(target = ADC_2P_C, p = 1 - fr),
    # depot(target = ADC_1P_C, p = (1 - xfr) * fr), depot(target = ADC_0P_C,
    # p = xfr * fr)).
    f(dar2_central) <- 1 - fr
    f(dar1_central) <- (1 - xfr) * fr
    f(dar0_central) <- xfr * fr

    # Concentrations in nM: states in nmol, vc in mL.
    c_dar2 <- 1000 * dar2_central / vc
    c_dar1 <- 1000 * dar1_central / vc
    c_dar0 <- 1000 * dar0_central / vc
    # Conjugated active payload (2 * ADC2P + ADC1P) and total payload
    # (2 * total ADC), nM.
    active_payload <- 2 * c_dar2 + c_dar1
    total_payload <- 2 * (c_dar2 + c_dar1 + c_dar0)

    # Observables (Monolix OUTPUT block). The active-ADC LBA detects ADC2P
    # plus a fraction pinact of ADC1P; the inactive-ADC LBA detects ADC0P
    # plus a fraction pinact2 of ADC1P.
    Cc <- c_dar2 + c_dar1 + c_dar0
    Cc_active <- c_dar2 + pinact * c_dar1
    Cc_inactive <- c_dar0 + pinact2 * c_dar1
    active_payload_pct <- 0
    if (total_payload > 0) active_payload_pct <- 100 * active_payload / total_payload

    Cc ~ lnorm(expSd)
    Cc_active ~ lnorm(expSd_Cc_active)
    Cc_inactive ~ lnorm(expSd_Cc_inactive)
    active_payload_pct ~ add(addSd_active_payload_pct) + prop(propSd_active_payload_pct) + combined1()
  })
}
