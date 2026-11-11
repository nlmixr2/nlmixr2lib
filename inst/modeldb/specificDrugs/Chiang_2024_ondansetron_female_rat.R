Chiang_2024_ondansetron_female_rat <- function() {
  description <- paste(
    "Preclinical (rat, female). Semi-physiological five-compartment PK model of",
    "intravenous ondansetron in female Sprague-Dawley rats (Chiang 2024): central",
    "and peripheral systemic disposition compartments with first-order",
    "elimination from central, and brain, spinal cord and cerebrospinal-fluid",
    "(CSF) compartments exchanging bidirectionally with central and with each",
    "other, plus first-order P-glycoprotein (Pgp) efflux from brain and spinal",
    "cord back to central. Fitted jointly to mean profiles of wild-type rats",
    "(ondansetron alone), Pgp-knockout (Mdr1a -/-) rats (Pgp efflux set to",
    "zero) and wild-type rats pretreated with the Pgp inhibitor tariquidar",
    "(separate central volume, elimination rate constant and Pgp efflux rate).",
    "CAUTION: with the published Table III values the model reproduces the",
    "paper's concentrations only for an input dose of about 0.7-0.8 mg/kg,",
    "1/12 to 1/14 of the stated 10 mg/kg (see the vignette).",
    sep = " "
  )
  reference <- paste(
    "Chiang M, Back H, Lee JB, Oh S, Guo T, Girgis S, Park C, Haroutounian S,",
    "Kagan L. Pharmacokinetic Modeling of the Effect of Tariquidar on",
    "Ondansetron Disposition into the Central Nervous System.",
    "Pharm Res. 2024;41(7):1401-1411. doi:10.1007/s11095-024-03739-6.",
    sep = " "
  )
  vignette <- "Chiang_2024_ondansetron"
  units <- list(time = "h", dosing = "ug", concentration = "ug/mL")

  covariateData <- list(
    CONMED_TARIQUIDAR = list(
      description = "Tariquidar co-administration indicator (1 = tariquidar pretreatment, 0 = ondansetron alone)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = ondansetron alone",
      notes = paste(
        "Chiang 2024 OT groups: tariquidar 7.5 mg/kg IV (5 mg/mL in 2% aqueous",
        "dextran) 30 min before ondansetron 10 mg/kg IV. Tariquidar plasma and",
        "tissue concentrations were nearly constant over the ondansetron",
        "experiment (Fig. 5), so the paper treats tariquidar as a categorical",
        "covariate: it switches V1 to V1,TRQ, kel to kel,TRQ and kPgp to",
        "kPgp,TRQ. Time-fixed per subject. Do not combine with STRAIN_ABCB1A_KO = 1;",
        "the paper studied tariquidar only in wild-type rats.",
        sep = " "
      ),
      source_name = "OT group (ondansetron + tariquidar)"
    ),
    STRAIN_ABCB1A_KO = list(
      description = "P-glycoprotein knockout rat indicator (1 = Mdr1a(-/-) SD-Abcb1a tm1sage, 0 = wild type)",
      units = "(binary)",
      type = "binary",
      reference_category = "0 = wild-type Sprague-Dawley",
      notes = paste(
        "Chiang 2024 Methods (Pharmacokinetic Modeling): Pgp efflux at the BBB",
        "'was included for wild type animals only; the term was set to zero in",
        "KO rats'. The KO strain is described in the upstream Chiang 2020 paper",
        "(Pharm Res 37:205) as Mdr1a(-/-), SD-Abcb1a tm1sage (Horizon",
        "Discovery / Sage Labs). Time-fixed (germline).",
        sep = " "
      ),
      source_name = "KO group"
    )
  )

  compartmentData <- list(
    central = list(analyte = "ondansetron", units = "ug", specimen = "plasma", verified = TRUE),
    peripheral1 = list(analyte = "ondansetron", units = "ug", specimen = "plasma", verified = TRUE),
    brain = list(analyte = "ondansetron", units = "ug", specimen = "tissue", verified = TRUE),
    spinal_cord = list(analyte = "ondansetron", units = "ug", specimen = "tissue", verified = TRUE),
    csf = list(analyte = "ondansetron", units = "ug", specimen = "CSF", verified = TRUE)
  )

  population <- list(
    species = "rat (Sprague-Dawley, wild type and Mdr1a(-/-) Pgp knockout)",
    sex_female_pct = 100,
    n_studies = 2,
    age_range = "12-14 weeks",
    weight_range = "250-280 g",
    disease_state = "healthy",
    dose_range = "ondansetron 10 mg/kg IV bolus (jugular cannula or tail vein), alone or 30 min after tariquidar 7.5 mg/kg IV",
    regions = "United States (Rutgers University)",
    notes = paste(
      "Mean (naive-pooled) concentration profiles fitted with the Ubiquity R",
      "package. Data: OT-F (ondansetron + tariquidar, this study; serial",
      "plasma sampling to 4 h plus terminal plasma, brain, spinal cord and",
      "cisterna-magna CSF at 0.16-2.5 h, n = 4-6 per time point) and WT-F /",
      "KO-F (ondansetron alone) from Chiang 2020 (Pharm Res 37:205,",
      "doi:10.1007/s11095-020-02929-2). Male rats were fitted separately; see",
      "Chiang_2024_ondansetron_male_rat.",
      sep = " "
    )
  )

  ini({
    # Chiang 2024 Table III, 'Female' column (estimate, CV%). All structural
    # values are rate constants (1/h) or volumes (mL); amounts are in ug so
    # concentrations are in ug/mL. Paper state numbering: A1 central,
    # A2 peripheral, A3 brain, A4 spinal cord, A5 CSF (Eqs. 1-5).
    lvc <- log(13.3); label("Central volume of distribution, ondansetron alone, V1 (mL)") # Table III V1 = 13.3 mL (CV 7%)
    lkel <- log(5.06); label("Systemic elimination rate constant, ondansetron alone, kel (1/h)") # Table III kel = 5.06 1/h (CV 11%)
    lk12 <- log(4.37); label("Central-to-peripheral rate constant k12 (1/h)") # Table III k12 = 4.37 1/h (CV 19%)
    lk21 <- log(2.32); label("Peripheral-to-central rate constant k21 (1/h)") # Table III k21 = 2.32 1/h (CV 10%)
    lk_central_brain <- log(13.0); label("Central-to-brain rate constant k13 (1/h)") # Table III k13 = 13.0 1/h (CV 14%)
    lk_brain_central <- log(25.2); label("Brain-to-central passive rate constant k31 (1/h)") # Table III k31 = 25.2 1/h (CV 18%)
    lk_central_spinal_cord <- log(3.60); label("Central-to-spinal-cord rate constant k14 (1/h)") # Table III k14 = 3.60 1/h (CV 17%)
    lk_spinal_cord_central <- log(31.5); label("Spinal-cord-to-central passive rate constant k41 (1/h)") # Table III k41 = 31.5 1/h (CV 22%)
    lk_central_csf <- log(0.518); label("Central-to-CSF rate constant k15 (1/h)") # Table III k15 = 0.518 1/h (CV 17%)
    lk_csf_central <- log(47.3); label("CSF-to-central rate constant k51 (1/h)") # Table III k51 = 47.3 1/h (CV 25%)
    lk_cns_csf <- log(0.914); label("Brain-to-CSF and spinal-cord-to-CSF rate constant k35 = k45 (1/h)") # Table III k35=45 = 0.914 1/h (CV 18%)
    lk_csf_cns <- log(79.9); label("CSF-to-brain and CSF-to-spinal-cord rate constant k53 = k54 (1/h)") # Table III k53=54 = 79.9 1/h (CV 15%)
    lk_cns_central_pgp <- log(101); label("Pgp-mediated efflux rate constant from brain and spinal cord to central, wild type, kPgp (1/h)") # Table III kPgp = 101 1/h (CV 19%)

    lv_brain <- fixed(log(1.8)); label("Brain volume Vbrain (mL)") # Table III Vbrain = 1.8 mL, fixed
    lv_spinal_cord <- fixed(log(0.6)); label("Spinal cord volume Vspinal (mL)") # Table III Vspinal = 0.6 mL, fixed
    lvcsf <- fixed(log(0.25)); label("CSF volume VCSF (mL)") # Table III VCSF = 0.25 mL, fixed

    # Tariquidar (OT) group: Table III estimates separate V1,TRQ, kel,TRQ and
    # kPgp,TRQ; each is encoded as the log-ratio to the ondansetron-alone
    # value so that exp(l<param> + e * CONMED_TARIQUIDAR) returns the printed
    # OT estimate exactly.
    e_conmed_tariquidar_vc <- log(11.3 / 13.3); label("Log-ratio of V1,TRQ to V1 (unitless)") # Table III V1,TRQ = 11.3 mL (CV 7%)
    e_conmed_tariquidar_kel <- log(3.02 / 5.06); label("Log-ratio of kel,TRQ to kel (unitless)") # Table III kel,TRQ = 3.02 1/h (CV 11%)
    e_conmed_tariquidar_k_cns_central_pgp <- log(0.0001 / 101); label("Log-ratio of kPgp,TRQ to kPgp (unitless)") # Table III kPgp,TRQ = 0.0001 1/h (CV 59%)

    # Residual error: Table III var_P = 0.03 (CV 19%). Chiang 2020 Methods
    # defines the Ubiquity variance model as VAR_i = sigma1 * Yhat_i^2, i.e.
    # proportional with SD = sqrt(sigma1); a single sigma1 is shared by every
    # observed matrix, so all four outputs carry the same value.
    propSd <- sqrt(0.03); label("Proportional residual SD, plasma (fraction)") # Table III var_P = 0.03
    propSd_Cbrain <- sqrt(0.03); label("Proportional residual SD, brain (fraction)") # Table III var_P = 0.03 (shared)
    propSd_Cspinal_cord <- sqrt(0.03); label("Proportional residual SD, spinal cord (fraction)") # Table III var_P = 0.03 (shared)
    propSd_Ccsf <- sqrt(0.03); label("Proportional residual SD, CSF (fraction)") # Table III var_P = 0.03 (shared)
  })
  model({
    vc <- exp(lvc + e_conmed_tariquidar_vc * CONMED_TARIQUIDAR)
    kel <- exp(lkel + e_conmed_tariquidar_kel * CONMED_TARIQUIDAR)
    k12 <- exp(lk12)
    k21 <- exp(lk21)
    k_central_brain <- exp(lk_central_brain)
    k_brain_central <- exp(lk_brain_central)
    k_central_spinal_cord <- exp(lk_central_spinal_cord)
    k_spinal_cord_central <- exp(lk_spinal_cord_central)
    k_central_csf <- exp(lk_central_csf)
    k_csf_central <- exp(lk_csf_central)
    k_cns_csf <- exp(lk_cns_csf)
    k_csf_cns <- exp(lk_csf_cns)
    # Pgp efflux is absent in knockout rats (Methods: 'set to zero in KO rats').
    k_cns_central_pgp <- exp(lk_cns_central_pgp + e_conmed_tariquidar_k_cns_central_pgp * CONMED_TARIQUIDAR) * (1 - STRAIN_ABCB1A_KO)
    v_brain <- exp(lv_brain)
    v_spinal_cord <- exp(lv_spinal_cord)
    vcsf <- exp(lvcsf)

    # Chiang 2024 Eqs. 1-5 (A1 central, A2 peripheral1, A3 brain,
    # A4 spinal_cord, A5 csf). Pgp efflux returns drug to central.
    d/dt(central) <- k21 * peripheral1 + (k_brain_central + k_cns_central_pgp) * brain +
      (k_spinal_cord_central + k_cns_central_pgp) * spinal_cord + k_csf_central * csf -
      (kel + k12 + k_central_brain + k_central_spinal_cord + k_central_csf) * central
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(brain) <- k_central_brain * central + k_csf_cns * csf -
      (k_brain_central + k_cns_csf + k_cns_central_pgp) * brain
    d/dt(spinal_cord) <- k_central_spinal_cord * central + k_csf_cns * csf -
      (k_spinal_cord_central + k_cns_csf + k_cns_central_pgp) * spinal_cord
    d/dt(csf) <- k_central_csf * central + k_cns_csf * (brain + spinal_cord) -
      (k_csf_central + 2 * k_csf_cns) * csf

    Cc <- central / vc
    Cbrain <- brain / v_brain
    Cspinal_cord <- spinal_cord / v_spinal_cord
    Ccsf <- csf / vcsf

    Cc ~ prop(propSd)
    Cbrain ~ prop(propSd_Cbrain)
    Cspinal_cord ~ prop(propSd_Cspinal_cord)
    Ccsf ~ prop(propSd_Ccsf)
  })
}
