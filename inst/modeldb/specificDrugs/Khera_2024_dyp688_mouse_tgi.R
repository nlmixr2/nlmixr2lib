Khera_2024_dyp688_mouse_tgi <- function() {
  description <- "Preclinical (mouse, 92.1 uveal melanoma xenograft). Population PK / tumour growth inhibition (TGI) model of DYP688, a PMEL17-targeted antibody-drug conjugate carrying two FR900359 (Galpha-q/11 inhibitor) payloads that inactivate while still attached to the antibody (Khera 2024). The DYP688 mouse PK model (three ADC species by number of active payloads, two compartments, non-cooperative on-antibody inactivation; parameters and their variability held at the PK-model estimates) drives Gompertz tumour growth. Plasma conjugated active payload (2 * ADC2P + ADC1P, nM) sets a linear kill rate kkill * C that reaches the tumour through one signal-distribution transit compartment (transit time tau); kkill has a categorical dose-level covariate (1.5 mg/kg reference, 3 and 6 mg/kg). Initial tumour volume is the TUM_VOL covariate. Dosing: give the SAME nmol amount as three simultaneous records into dar2_central, dar1_central and dar0_central; f() splits it into 1 - fr, (1 - xfr) * fr and xfr * fr."
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
    concentration = "nM for plasma ADC and conjugated payload; mm^3 for tumour volume"
  )

  covariateData <- list(
    TUM_VOL = list(
      description = "Tumour volume at the start of the modelled period; initial condition of the tumour-size state (Monolix TS_0 = TS0).",
      units = "mm^3",
      type = "continuous",
      reference_category = NULL,
      notes = "Monolix regressor TS0 ('TS0 = {use = regressor}'), the observed tumour volume at randomisation (mean approximately 150-400 mm^3, Khera 2024 Methods). Tumour volume = L * W^2 * pi / 6.",
      source_name = "TS0"
    ),
    DOSE_DYP688_MGKG = list(
      description = "DYP688 dose per administration of the animal's treatment group (mg/kg); categorical covariate on the kill rate constant kkill.",
      units = "mg/kg",
      type = "continuous",
      reference_category = "1.5 mg/kg (any value other than 3 or 6 uses the reference kkill)",
      notes = "Monolix categorical covariate Dose with levels 1.5 (reference), 3 and 6 mg/kg (beta_kkill_Dose_3, beta_kkill_Dose_6). The 3 mg/kg Q2W and Q4W groups share the 3 mg/kg level. Doses other than 3 or 6 mg/kg take the reference (smallest) kkill, which is the conservative choice the authors made for all subsequent simulations (Khera 2024 Results, before Table 4).",
      source_name = "Dose"
    )
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
    ),
    effect = list(
      analyte = "tumour kill-rate signal after one transit compartment (K1)",
      units = "1/h",
      specimen = "not applicable",
      verified = TRUE
    ),
    tumor_size = list(
      analyte = "92.1 xenograft tumour volume",
      units = "mm^3",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "mouse (female nude Crl:NU(NCr)-Foxn1nu, 92.1-luc uveal melanoma xenograft)",
    n_subjects = 24L,
    n_studies = 1L,
    age_range = "6-8 weeks at cell inoculation",
    sex_female_pct = 100,
    disease_state = "Subcutaneous 92.1-luc uveal melanoma cell-derived xenograft (Class 1 uveal melanoma, GNAQ-mutant)",
    dose_range = "IV 1.5 mg/kg Q2W, 3 mg/kg Q4W, 3 mg/kg Q2W and 6 mg/kg Q4W, plus saline vehicle",
    regions = "Novartis, Basel (Switzerland)",
    notes = "Khera 2024 Methods 'Preclinical in vivo study design' and Results 'PK-TGI modeling of mouse tumors'. Groups of 4-6 animals randomised at mean tumour volumes of approximately 150-400 mm^3. Vehicle-only Gompertz fit: 13 animals, 110 tumour-volume observations; PK-TGI fit: 24 animals, 616 observations, 71 doses (Supplementary Material Monolix output). Flat, unmeasurable tumours at 3 and 6 mg/kg were censored as BLOQ."
  )

  ini({
    # PK layer: values held constant in the PK-TGI Monolix run (Supplement
    # 'Mouse PK-TGI model results summary', 04_TS0CatCov_...mlxtran), i.e.
    # the DYP688 mouse PK estimates (Table 1) rounded to two significant
    # figures.
    lkel <- fixed(log(0.016))
    label("Antibody / ADC elimination rate constant kel (1/h)") # Supplement PK-TGI Monolix output kel_pop = 0.016 (no SE); Table 1
    lk12 <- fixed(log(0.035))
    label("Central-to-peripheral distribution rate constant k12 (1/h)") # Supplement PK-TGI Monolix output k12_pop = 0.035 (no SE); Table 1
    lk21 <- fixed(log(0.036))
    label("Peripheral-to-central distribution rate constant k21 (1/h)") # Supplement PK-TGI Monolix output k21_pop = 0.036 (no SE); Table 1
    lkdec <- fixed(log(0.013))
    label("On-antibody inactivation rate constant of each active payload ki.a. (1/h)") # Supplement PK-TGI Monolix output kia1_pop = 0.013 (no SE); Table 1
    lvc <- fixed(log(1.53))
    label("Plasma volume of distribution Vc (mL)") # Supplement PK-TGI Monolix output Vc_pop = 1.53 (no SE); Table 1
    fr <- fixed(0.027)
    label("Fraction of inactive ADC in the dose fr (fraction)") # Supplement PK-TGI Monolix output fr_pop = 0.027 (no SE); Table 1
    xfr <- fixed(0.99)
    label("Fraction of the inactive dose fraction carried as ADC with no active payload xfr (fraction)") # Supplement PK-TGI Monolix output xfr_pop = 0.99 (no SE)

    # Tumour growth (Gompertz), from the vehicle-only fit and held constant
    # in the PK-TGI run (Table 3, upper block).
    lbeta <- fixed(log(0.0059))
    label("Gompertz tumour growth rate constant beta (1/h)") # Table 3, beta = 0.0059 1/h (RSE 17.2%, vehicle fit); PK-TGI Monolix output beta_pop = 0.0059 (no SE)
    ltsmax <- fixed(log(929))
    label("Maximum tumour volume TSmax (mm^3)") # PK-TGI Monolix output TSmax_pop = 929 (no SE); Table 3 928.39 mm^3 (RSE 8.48%, vehicle fit)

    # Drug effect (Table 3, lower block; full precision from the PK-TGI
    # Monolix output).
    lkkill <- log(0.000165)
    label("Tumour kill rate constant per nM of conjugated active payload kkill, 1.5 mg/kg reference (1/(nM*h))") # PK-TGI Monolix output kkill_pop = 0.000165 (RSE 47.4%); Table 3 0.00016
    e_dose3_kkill <- 0.719
    label("Log-scale shift in kkill for the 3 mg/kg dose level (unitless)") # PK-TGI Monolix output beta_kkill_Dose_3 = 0.719 (RSE 75.3%); Table 3 0.72
    e_dose6_kkill <- 1.18
    label("Log-scale shift in kkill for the 6 mg/kg dose level (unitless)") # PK-TGI Monolix output beta_kkill_Dose_6 = 1.18 (RSE 51.7%); Table 3 1.18
    ltau <- log(699)
    label("Signal-distribution transit time tau (h)") # PK-TGI Monolix output tau_pop = 699 (RSE 0.0182%); Table 3 698.87 h

    # IIV (Monolix SD of the log-scale random effect; variances = omega^2).
    # The PK omegas and omega_tau carry no standard error in the PK-TGI run.
    # beta and TSmax had IIV in the vehicle fit (Table 3) but none in the
    # PK-TGI run, whose individual beta and TSmax are all equal to the
    # population values (Supplement, conditional-distribution summary).
    etalkel ~ fixed(0.2025) # PK-TGI Monolix output omega_kel = 0.45 (no SE)
    etalk12 ~ fixed(0.0225) # PK-TGI Monolix output omega_k12 = 0.15 (no SE)
    etalk21 ~ fixed(0.0676) # PK-TGI Monolix output omega_k21 = 0.26 (no SE)
    etalkdec ~ fixed(0.005625) # PK-TGI Monolix output omega_kia1 = 0.075 (no SE)
    etalvc ~ fixed(0.0009) # PK-TGI Monolix output omega_Vc = 0.03 (no SE)
    etalkkill ~ 0.804609 # PK-TGI Monolix output omega_kkill = 0.897 (RSE 19.4%); Table 3 IISD 0.9
    etaltau ~ fixed(0.01) # PK-TGI Monolix output omega_tau = 0.1 (no SE); Table 3 IISD 0.1

    # Residual error: Monolix constant error on the log scale (lognormal).
    expSd <- 0.734
    label("Log-scale residual SD of tumour volume (unitless)") # PK-TGI Monolix output a4 = 0.734 (RSE 3.35%); Table 3 0.73
  })

  model({
    kel <- exp(lkel + etalkel)
    k12 <- exp(lk12 + etalk12)
    k21 <- exp(lk21 + etalk21)
    kdec <- exp(lkdec + etalkdec)
    vc <- exp(lvc + etalvc)
    beta <- exp(lbeta)
    tsmax <- exp(ltsmax)
    tau <- exp(ltau + etaltau)

    # Categorical dose-level covariate on log(kkill) (Monolix: lognormal
    # kkill with categorical covariate Dose, reference 1.5 mg/kg).
    dose3 <- 0
    if (DOSE_DYP688_MGKG == 3) dose3 <- 1
    dose6 <- 0
    if (DOSE_DYP688_MGKG == 6) dose6 <- 1
    kkill <- exp(lkkill + e_dose3_kkill * dose3 + e_dose6_kkill * dose6 + etalkkill)

    # PK (Supplement 'Mouse PK-TGI model structure file', TGI model).
    d/dt(dar2_central) <- -k12 * dar2_central + k21 * dar2_peripheral1 - 2 * kdec * dar2_central - kel * dar2_central
    d/dt(dar1_central) <- -k12 * dar1_central + k21 * dar1_peripheral1 + 2 * kdec * dar2_central - kdec * dar1_central - kel * dar1_central
    d/dt(dar0_central) <- -k12 * dar0_central + k21 * dar0_peripheral1 + kdec * dar1_central - kel * dar0_central
    d/dt(dar2_peripheral1) <- k12 * dar2_central - k21 * dar2_peripheral1 - 2 * kdec * dar2_peripheral1
    d/dt(dar1_peripheral1) <- k12 * dar1_central - k21 * dar1_peripheral1 + 2 * kdec * dar2_peripheral1 - kdec * dar1_peripheral1
    d/dt(dar0_peripheral1) <- k12 * dar0_central - k21 * dar0_peripheral1 + kdec * dar1_peripheral1

    f(dar2_central) <- 1 - fr
    f(dar1_central) <- (1 - xfr) * fr
    f(dar0_central) <- xfr * fr

    # Conjugated active payload in plasma, nM (EXPOSURE = Active_P =
    # 2 * ADC_2P_C / Vc + ADC_1P_C / Vc; states in nmol, vc in mL).
    active_payload <- 1000 * (2 * dar2_central + dar1_central) / vc
    Cc <- 1000 * (dar2_central + dar1_central + dar0_central) / vc

    # Signal distribution: one transit compartment from the kill rate
    # kkill * max(0, EXPOSURE) to the effective kill rate K1 (K1_0 = 0).
    d/dt(effect) <- (kkill * max(0, active_payload) - effect) / tau

    # Gompertz growth with first-order kill. The authors bound log() away
    # from zero with max(1e-6, TS) and stop growth above 1e12 mm^3.
    tumor_size(0) <- TUM_VOL
    ts_rate <- tumor_size * beta * log(tsmax / max(1e-6, tumor_size)) - effect * tumor_size
    if (tumor_size > 1e12) ts_rate <- 0
    d/dt(tumor_size) <- ts_rate

    tumor_size ~ lnorm(expSd)
  })
}
