SanzCodina_2024_cefiderocol <- function() {
  description <- "Two-compartment population PK model for intravenous cefiderocol in healthy adult male volunteers, fitted to total and unbound plasma concentrations and to microdialysis-derived unbound interstitial concentrations in subcutaneous adipose and skeletal muscle tissue (Sanz-Codina 2024). Clearance and volumes are referenced to the UNBOUND concentration; total plasma cefiderocol is Cunbound / fu with a constant estimated fraction unbound, and each tissue concentration is the peripheral-compartment concentration multiplied by a tissue scaling factor."
  reference <- paste(
    "Sanz-Codina M, van Os W, Pham AD, Jorda A, Wolf-Duchek M, Bergmann F,",
    "Lackner E, Lier C, van Hasselt JGC, Minichmayr IK, Dorn C, Zeitlinger M,",
    "al Jalali V. Target-site cefiderocol pharmacokinetics in soft tissues of",
    "healthy volunteers. J Antimicrob Chemother. 2024;79(12):3281.",
    "doi:10.1093/jac/dkae359. PMID 39373642; PMCID PMC11638555.",
    sep = " "
  )
  vignette <- "SanzCodina_2024_cefiderocol"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # The final model has no covariates: 'The small and homogeneous healthy
  # volunteer population precluded a systematic covariate analysis' (Methods,
  # 'Population PK analysis').
  covariateData <- list()

  compartmentData <- list(
    central = list(analyte = "cefiderocol", units = "mg", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "cefiderocol", units = "mg", specimen = "tissue", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 8L,
    n_studies = 1L,
    age_range = "18-55 years (inclusion criterion); mean 34.9 +/- 8.9 years",
    bmi_range = "19-30 kg/m^2 (inclusion criterion); mean 24.6 +/- 1.1 kg/m^2",
    sex_female_pct = 0,
    race_ethnicity = "not reported in the source paper",
    disease_state = "Healthy male volunteers; non-smokers, no regular medication, creatinine clearance > 90 mL/min (Cockcroft-Gault)",
    dose_range = "Single 2 g cefiderocol (Fetcroja) intravenous infusion over 3 h",
    regions = "Austria (Medical University of Vienna)",
    renal_function = "Normal; creatinine clearance <= 90 mL/min was an exclusion criterion",
    notes = paste(
      "Plasma sampled pre-dose and at 1, 2, 3 (end of infusion), 3.5, 4, 5, 6, 7 and 8 h;",
      "ultrafiltration (unbound plasma) at 3, 4 and 7 h. Unbound interstitial concentrations",
      "in subcutaneous adipose and skeletal muscle tissue of the same thigh by microdialysis",
      "over the collection intervals 0-1, 1-2, 2-3, 3-3.5, 3.5-4, 4-5, 5-6, 6-7 and 7-8 h,",
      "corrected for the per-catheter in vivo relative recovery (retrodialysis). Two",
      "implausible subcutaneous microdialysate observations were excluded from the model.",
      "NONMEM 7.4 FOCE-I (Methods, 'Population PK analysis'; Results, 'Study population and safety')."
    )
  )

  ini({
    # Structural parameters. Source: Table 2, 'Fixed effects'. Table 2 footnote:
    # 'Estimates are based on unbound cefiderocol concentrations', so CL, Vc, Q
    # and Vp are referenced to the unbound plasma concentration.
    lcl <- log(11.0); label("Clearance referenced to unbound cefiderocol (L/h)") # Table 2, row 'CL (L/h)' = 11.0 (RSE 6.7%)
    lvc <- log(11.3); label("Central volume referenced to unbound cefiderocol (L)") # Table 2, row 'Vc (L)' = 11.3 (RSE 16.3%)
    lq <- log(25.3); label("Intercompartmental clearance referenced to unbound cefiderocol (L/h)") # Table 2, row 'Q12 (L/h)' = 25.3 (RSE 23.0%)
    lvp <- log(15.1); label("Peripheral volume referenced to unbound cefiderocol (L)") # Table 2, row 'Vp (L)' = 15.1 (RSE 14.2%)

    # Constant plasma protein binding; the unbound fraction was estimated from
    # the paired total and ultrafiltrate (unbound) plasma samples at 3, 4 and
    # 7 h (Methods, 'Population PK analysis'). No IIV on fu (Results).
    lfu <- log(0.639); label("Fraction of cefiderocol unbound in plasma (unitless)") # Table 2, row 'fu' = 0.639 (RSE 1.5%)

    # Tissue scaling factors. Interstitial tissue concentrations are the
    # peripheral-compartment concentration scaled by TF (Results, 'Population
    # PK model and PTA analysis': 'tissue concentrations were described by
    # scaling the tissue concentrations to concentrations in the peripheral
    # compartment of the plasma model'). The peripheral compartment carries
    # unbound-referenced drug, so at distribution equilibrium TF is the
    # unbound tissue-to-unbound plasma concentration ratio.
    lkp_adipose <- log(0.943); label("Scaling factor, unbound subcutaneous adipose interstitial to peripheral-compartment concentration (unitless)") # Table 2, row 'TFsub' = 0.943 (RSE 10.5%)
    lkp_muscle <- log(0.901); label("Scaling factor, unbound skeletal muscle interstitial to peripheral-compartment concentration (unitless)") # Table 2, row 'TFmus' = 0.901 (RSE 11.2%)

    # IIV. Exponential model (Methods). Table 2 reports CV% calculated as
    # sqrt(exp(omega^2) - 1) * 100 (Table 2 footnote), so omega^2 = log(CV^2 + 1).
    etalcl ~ 0.0298267 # Table 2, IIV row 'CL (CV%)' = 17.4 -> log(0.174^2 + 1)
    etalvc ~ 0.0678689 # Table 2, IIV row 'Vc (CV%)' = 26.5 -> log(0.265^2 + 1)
    etalkp_adipose ~ 0.0754785 # Table 2, IIV row 'TFsub (CV%)' = 28.0 -> log(0.280^2 + 1)
    etalkp_muscle ~ 0.0615691 # Table 2, IIV row 'TFmus (CV%)' = 25.2 -> log(0.252^2 + 1)

    # Residual error. One proportional error for plasma (Table 2 'Propplasma')
    # applied to both the total and the unbound plasma observations, and ONE
    # combined proportional + additive error shared by the two tissues (Results:
    # 'The residual unexplained variability for tissue PK was described using
    # one combined error model'). nlmixr2 needs a parameter per endpoint, so
    # each shared estimate is carried on every endpoint it applies to.
    propSd <- 0.0585; label("Proportional residual error, total plasma cefiderocol (fraction)") # Table 2, row 'Propplasma (CV%)' = 5.85 (RSE 6.9%)
    propSd_Cunbound <- 0.0585; label("Proportional residual error, unbound plasma cefiderocol (fraction)") # Table 2, row 'Propplasma (CV%)' = 5.85; same plasma estimate as propSd
    propSd_Cadipose <- 0.303; label("Proportional residual error, subcutaneous adipose tissue (fraction)") # Table 2, row 'Proptissue (CV%)' = 30.3 (RSE 13.5%)
    addSd_Cadipose <- 2.3; label("Additive residual error, subcutaneous adipose tissue (mg/L)") # Table 2, row 'Addtissue (mg/L)' = 2.3 (RSE 5.6%)
    propSd_Cmuscle <- 0.303; label("Proportional residual error, skeletal muscle tissue (fraction)") # Table 2, row 'Proptissue (CV%)' = 30.3; same tissue estimate as propSd_Cadipose
    addSd_Cmuscle <- 2.3; label("Additive residual error, skeletal muscle tissue (mg/L)") # Table 2, row 'Addtissue (mg/L)' = 2.3; same tissue estimate as addSd_Cadipose
  })

  model({
    # 1. Individual parameters (exponential IIV on CL, Vc, TFsub and TFmus).
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc + etalvc)
    q <- exp(lq)
    vp <- exp(lvp)
    fu <- exp(lfu)
    kp_adipose <- exp(lkp_adipose + etalkp_adipose)
    kp_muscle <- exp(lkp_muscle + etalkp_muscle)

    # 2. Micro-constants.
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # 3. Two-compartment disposition with first-order elimination (Results).
    #    Doses are IV infusions into central.
    d/dt(central) <- -(kel + k12) * central + k21 * peripheral1
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1

    # 4. Observations. The volumes are unbound-referenced, so central / vc is
    #    the unbound plasma concentration and total plasma cefiderocol follows
    #    from the constant binding model, C_total = C_unbound / fu. The tissue
    #    outputs are the unbound interstitial concentrations measured by
    #    microdialysis. The paper fitted microdialysate as the integral of the
    #    tissue curve over each collection interval; these outputs are the
    #    instantaneous tissue concentrations.
    Cunbound <- central / vc
    Cc <- Cunbound / fu
    Cadipose <- kp_adipose * peripheral1 / vp
    Cmuscle <- kp_muscle * peripheral1 / vp

    Cc ~ prop(propSd)
    Cunbound ~ prop(propSd_Cunbound)
    Cadipose ~ add(addSd_Cadipose) + prop(propSd_Cadipose)
    Cmuscle ~ add(addSd_Cmuscle) + prop(propSd_Cmuscle)
  })
}
