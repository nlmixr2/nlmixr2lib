Shinha_2020_irinotecan_invitro <- function() {
  description <- paste(
    "In vitro (multi-organ-on-a-chip: HepG2 liver part + A549 lung-cancer part).",
    "Deterministic parent-metabolite PK-PD model for the prodrug irinotecan",
    "(CPT-11) and its active metabolite SN-38 in the recirculating culture",
    "medium of a microfluidic chip. Both species share one well-mixed medium",
    "volume (the microchannel volume, vc); CPT-11 is converted to SN-38 by the",
    "liver part at a first-order rate q_liver * eh / vc (flow rate times",
    "extraction ratio, Shinha 2020 Eq. 5) and SN-38 is eliminated by the same",
    "liver part at q_liver * eh_sn38 / vc. Cancer-cell density, as a percentage",
    "of the untreated control, is a log-linear function of the cumulative SN-38",
    "AUC (Eq. 7). Concomitant simvastatin halves the CPT-11 extraction ratio",
    "(CES2 inhibition, Table II). The flow rate and volume default to the",
    "physiological-flow-ratio chip (with bypass channel); the chip without the",
    "bypass channel is simulated by overriding lq_liver and lvc (Table I). The",
    "culture medium is exchanged every 24 h, which must be encoded in the event",
    "table as replacement events (see the vignette).",
    sep = " "
  )

  reference <- paste(
    "Shinha K, Nihei W, Ono T, Nakazato R, Kimura H.",
    "A pharmacokinetic-pharmacodynamic model based on multi-organ-on-a-chip for",
    "drug-drug interaction studies.",
    "Biomicrofluidics. 2020;14(4):044108. doi:10.1063/5.0011545.",
    "Model equations from Section II.D (Eqs. 1-8); chip parameters from",
    "Table I (estimation experiments) and Table II (DDI simulations);",
    "extraction-ratio estimates from Section III.B.",
    sep = " "
  )

  vignette <- "Shinha_2020_irinotecan_invitro"

  # auc_sn38 is a running integral of the SN-38 medium concentration that
  # drives the PD (Eq. 8), not a biological compartment. Same idiom as
  # auc_central in Beguin_2024_carboplatin_thrombocytopenia_dog.R.
  paper_specific_compartments <- c("auc_sn38")

  units <- list(
    time = "h",
    dosing = "ng",
    concentration = "ng/uL"
  )

  covariateData <- list(
    CONMED_SIMVASTATIN = list(
      description = "Concomitant simvastatin (CES2 inhibitor) in the culture medium, 1 = present, 0 = absent.",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "Shinha 2020 Section II.F: simvastatin at 1 uM, a concentration reported to reduce CES2 expression to 50 percent (ref. 20, Fukami 2010).",
        "Section III.C and Table II: the CPT-11 extraction ratio was therefore set to 0.2 percent with simvastatin versus 0.4 percent without, i.e. eh is multiplied by (1 - 0.5).",
        "In this in-vitro model the column flags co-incubation in the chip medium rather than a patient co-medication; the quantity (presence of the co-administered drug) is the same."
      ),
      source_name = "w/ SV"
    )
  )

  covariatesDataExcluded <- list(
    CONMED_RTV = list(
      description = "Concomitant ritonavir (CYP3A4 inhibitor) in the culture medium, 1 = present, 0 = absent.",
      units = "(binary)",
      type = "binary",
      reference_category = 0,
      notes = paste(
        "Ritonavir 10 uM was tested (Section II.F), but the model assigns it NO effect: Table II carries the same CPT-11 extraction ratio (0.004) with and without ritonavir because HepG2 cells express very little CYP3A4 and CPT-11 is metabolised only by CES2 in this chip (fm = 1).",
        "The measured cell density ratio with ritonavir (20.3 percent) was lower than without (29.5 percent, not significant); the authors attribute this to UGT1A1 inhibition of SN-38 glucuronidation, which the model does not include."
      ),
      source_name = "w/ RTV"
    )
  )

  compartmentData <- list(
    central = list(analyte = "irinotecan", units = "ng", specimen = "administration site", verified = TRUE),
    central_sn38 = list(analyte = "SN-38", units = "ng", specimen = "administration site", verified = TRUE),
    auc_sn38 = list(analyte = "SN-38", units = "h*ng/uL", specimen = "not applicable", verified = TRUE)
  )

  population <- list(
    species = "in vitro (HepG2 liver model + A549 lung cancer cells on a multi-organ-on-a-chip)",
    n_subjects = NA_integer_,
    n_studies = 1L,
    age_range = NA_character_,
    weight_range = NA_character_,
    sex_female_pct = NA_real_,
    race_ethnicity = NA_character_,
    disease_state = "Not applicable -- A549 human lung adenocarcinoma cells as the drug-target part, HepG2 hepatocellular carcinoma cells as the metabolising liver part.",
    dose_range = "CPT-11 15 uM (9.35 ng/uL) in the culture medium, replaced every 24 h for 72 h; simvastatin 1 uM or ritonavir 10 uM in the DDI experiments.",
    regions = "Japan (Tokai University)",
    notes = paste(
      "Polydimethylsiloxane multi-organ-on-a-chip with a stirrer-based micropump at 2800 rpm (Sections II.A-B, III.A).",
      "HepG2 seeded at 1.7-2.0 x 10^5 cells/cm^2 in the liver chamber and, 48 h later, A549 at 1.7-2.0 x 10^4 cells/cm^2 in the lung-cancer chamber (Section II.C).",
      "Two chip designs: without the bypass channel (liver:lung-cancer flow ratio 1.0:1.0, Q = 89.60 uL/h, Vd = 70.19 uL) and with the bypass channel (physiological 1.0:3.3 ratio, Q = 26.45 uL/h, Vd = 85.01 uL) (Table I).",
      "Endpoint: A549 nuclear density (Hoechst 33342) after 72 h relative to CPT-11-free control chips; mean +/- SD of n = 3-6 chips per condition (Figs. 3-4).",
      "The PD relationship (Eq. 7) was built from published SN-38 cytotoxicity on A549 cells (ref. 17, Mijatovic 2006); its coefficients are inputs to this study, not estimated from the chip data."
    )
  )

  ini({
    # Chip design constants (Table II, 'w/o inhibitors' column = Table I 'w/
    # bypass channel' column). Determined from the chip design, not estimated
    # (Section II.E: 'Vd, Q, and X0 were determined from the MOoC design and
    # the experimental conditions').
    lvc <- fixed(log(85.01)); label("Distribution volume = microchannel medium volume, chip with bypass channel (uL)") # Table I/II 'Vd' w/ bypass = 85.01 uL (w/o bypass: 70.19 uL)
    lq_liver <- fixed(log(26.45)); label("Medium flow rate through the liver part, chip with bypass channel (uL/h)") # Table I/II 'Q' w/ bypass = 26.45 uL/h (w/o bypass: 89.60 uL/h)

    # Drug-specific parameters, estimated from the chip experiments (Section
    # III.B: 'The extraction ratios of CPT-11 and SN-38 in the liver part were
    # estimated to be 0.4% and 8.4%'). No uncertainty was reported.
    eh <- 0.004; label("Liver-part extraction ratio of CPT-11 (irinotecan) (fraction)") # Section III.B and Table II 'Ep' = 0.004
    eh_sn38 <- 0.084; label("Liver-part extraction ratio of SN-38 (fraction)") # Section III.B and Table II 'Em' = 0.084
    fm <- fixed(1); label("Fraction of CPT-11 metabolised to SN-38 by CES2 (fraction)") # Table I/II 'fm' = 1; Section III.B 'fmCES2 ... was set to 1'

    # Simvastatin: CES2 expression halved, so Ep 0.004 -> 0.002 (Table II
    # 'w/ SV'; Section III.C 'set to 0.2% because the CES2 expression would be
    # half', citing ref. 20).
    e_conmed_simvastatin_eh <- fixed(0.5); label("Fractional reduction of the CPT-11 extraction ratio with concomitant simvastatin (fraction)") # Table II 'Ep' w/ SV 0.002 vs 0.004 = 1 - 0.5

    # PD: Eq. 7, cell density / control cell density = -0.086 * ln(AUC) + 0.512,
    # AUC of SN-38 in h*ng/uL. Built from the literature (ref. 17) before the
    # chip experiments and held constant here.
    viability_ref <- fixed(0.512); label("Cell density ratio at an SN-38 AUC of 1 h*ng/uL (fraction of control)") # Eq. 7 intercept 0.512
    e_auc_sn38 <- fixed(-0.086); label("Change in cell density ratio per unit natural-log SN-38 AUC (fraction of control)") # Eq. 7 slope -0.086
  })

  model({
    vc <- exp(lvc)
    q_liver <- exp(lq_liver)

    # Simvastatin inhibits CES2-mediated conversion of CPT-11 (Table II)
    eh_i <- eh * (1 - e_conmed_simvastatin_eh * CONMED_SIMVASTATIN)

    # Eq. 5: k = CL / Vd = Q * E / Vd
    kel <- q_liver * eh_i / vc
    kel_sn38 <- q_liver * eh_sn38 / vc

    # Eqs. 1 and 3, written on amounts in the shared medium volume
    d/dt(central) <- -kel * central
    d/dt(central_sn38) <- fm * kel * central - kel_sn38 * central_sn38
    # Eq. 8: the PD is driven by the cumulative SN-38 AUC in the medium
    d/dt(auc_sn38) <- central_sn38 / vc

    Cc <- central / vc
    Cc_sn38 <- central_sn38 / vc

    # Eq. 7, expressed in percent of control as plotted in Figs. 3-4. The
    # relationship is an empirical log-linear fit and is only meaningful for
    # SN-38 AUC values near the 72-h experimental range (about 3-20 h*ng/uL);
    # it diverges as auc_sn38 -> 0.
    viability <- 100 * (viability_ref + e_auc_sn38 * log(auc_sn38))
  })
}
