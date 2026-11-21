Soeorg_2026_meropenem_colistin_polymyxinB <- function() {
  description <- "In vitro (Acinetobacter baumannii, 60 isolates of 53 strains pooled from 21 published time-kill and hollow-fibre studies). Semi-mechanistic PK/PD model of meropenem and colistin or polymyxin B, alone and in combination, fitted to digitised log10 CFU/mL time-kill data. A single bacterial population grows logistically and is killed by each drug through an Emax term whose EC50 is a multiple of the isolate MIC (with interisolate variability). Each drug's kill is lost over time through a two-state attenuation switched on while the drug is present. The meropenem Emax depends on the starting inoculum (power model) and on an acquired carbapenemase; the colistin/polymyxin B Emax is reduced at an 8 log10 CFU/mL inoculum. The two kills combine under Bliss independence, with a general pharmacodynamic interaction (GPDI) through which colistin/polymyxin B raises the meropenem Emax and lowers its EC50. Colistin and polymyxin B share one set of parameters. Drug concentrations, MICs, inoculum and carbapenemase status are covariates (no dosing events); a time-varying drug concentration can be supplied to drive the model from a PK profile."
  reference <- "Soeorg H, Lyu S, Kloprogge F, Standing JF. Pharmacokinetic-pharmacodynamic modeling using combined time-kill data for meropenem and colistin or polymyxin B combination effects in Acinetobacter baumannii. Antimicrob Agents Chemother. 2026. doi:10.1128/aac.00274-26. Structural equations: main-text Eq 1-13; final estimates: Table 1; model-building steps: Supplementary Table S2."
  vignette <- "Soeorg_2026_meropenem_colistin_polymyxinB"
  units <- list(time = "h", dosing = "mg/L", concentration = "mg/L")

  compartmentData <- list(
    cfu = list(
      analyte = "Acinetobacter baumannii (total viable bacteria)",
      units = "CFU/mL",
      specimen = "not applicable",
      verified = TRUE
    ),
    aron_mero = list(
      analyte = "fraction of bacteria in the reduced-killing (attenuated) state for meropenem",
      units = "fraction",
      specimen = "not applicable",
      verified = TRUE
    ),
    aroff_mero = list(
      analyte = "fraction of bacteria in the maximum-killing state for meropenem",
      units = "fraction",
      specimen = "not applicable",
      verified = TRUE
    ),
    aron_col = list(
      analyte = "fraction of bacteria in the reduced-killing (attenuated) state for colistin/polymyxin B",
      units = "fraction",
      specimen = "not applicable",
      verified = TRUE
    ),
    aroff_col = list(
      analyte = "fraction of bacteria in the maximum-killing state for colistin/polymyxin B",
      units = "fraction",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list(
    CONC_MEM_MGL = list(
      description = "Meropenem concentration in the time-kill medium. Drives the meropenem Emax kill and switches on its attenuation while above zero.",
      units = "mg/L",
      type = "continuous",
      reference_category = "n/a -- 0 mg/L for growth controls and colistin/polymyxin B monotherapy",
      notes = "Static in the static time-kill experiments (0 to 125 x MIC across the pooled studies, Table S1) and time-varying in the hollow-fibre experiments. In the paper's clinical simulations it is the predicted unbound epithelial-lining-fluid meropenem concentration (plasma x 0.625 unbound x 0.25 lung penetration).",
      source_name = "C_MERO"
    ),
    CONC_CST_MGL = list(
      description = "Colistin concentration in the time-kill medium. Colistin and polymyxin B share one set of parameters; the model uses CONC_CST_MGL + CONC_PMB_MGL as the polymyxin concentration (C_COLPMB), so supply 0 for whichever drug was not used.",
      units = "mg/L",
      type = "continuous",
      reference_category = "n/a -- 0 mg/L when colistin was not used",
      notes = "The fit found no difference between colistin and polymyxin B potency (Table S2, model P7 vs P6), so the paper pools them as C_COLPMB. Switches on the colistin/polymyxin B attenuation while the pooled concentration is above zero.",
      source_name = "C_COLPMB"
    ),
    CONC_PMB_MGL = list(
      description = "Polymyxin B concentration in the time-kill medium. Summed with CONC_CST_MGL into the pooled colistin/polymyxin B concentration (C_COLPMB); supply 0 when polymyxin B was not used.",
      units = "mg/L",
      type = "continuous",
      reference_category = "n/a -- 0 mg/L when polymyxin B was not used",
      notes = "Total concentration as reported by the source time-kill studies; see CONC_CST_MGL.",
      source_name = "C_COLPMB"
    ),
    MIC_MEM = list(
      description = "Meropenem MIC of the isolate. The meropenem EC50 is EC50_MERO x MIC_MEM (EC50_MERO fixed to 1).",
      units = "mg/L",
      type = "continuous",
      reference_category = "n/a -- continuous; EC50 scales linearly with the MIC",
      notes = "Reported by each source study (broth microdilution in 15 of 21 studies). 57 of the 60 isolates were meropenem-resistant (EUCAST MIC > 8 mg/L; Fig 1A).",
      source_name = "MIC_MERO"
    ),
    MIC_CST = list(
      description = "Colistin MIC of the isolate. Used as MIC_COLPMB (scaling the colistin/polymyxin B EC50 and the interaction potency) when the polymyxin in use is colistin, i.e. when CONC_PMB_MGL is 0.",
      units = "mg/L",
      type = "continuous",
      reference_category = "n/a -- continuous; EC50 scales linearly with the MIC",
      notes = "14 of the 60 isolates were colistin/polymyxin B-resistant (EUCAST MIC > 2 mg/L; Fig 1A). The paper does not distinguish colistin and polymyxin B MICs.",
      source_name = "MIC_COLPMB"
    ),
    MIC_PMB = list(
      description = "Polymyxin B MIC of the isolate. Used as MIC_COLPMB when the polymyxin in use is polymyxin B, i.e. when CONC_PMB_MGL is above 0.",
      units = "mg/L",
      type = "continuous",
      reference_category = "n/a -- continuous; EC50 scales linearly with the MIC",
      notes = "See MIC_CST.",
      source_name = "MIC_COLPMB"
    ),
    CARBAPENEMASE_POS = list(
      description = "1 = the isolate carries an acquired meropenem-degrading beta-lactamase (OXA-23, OXA-24/40, OXA-58, OXA-72, VIM-1 or NDM-1); 0 = none reported. Multiplies the meropenem Emax by e_carbapenemase_emax_mero.",
      units = "(binary)",
      type = "binary",
      reference_category = "0 (no acquired carbapenemase)",
      notes = "The intrinsic A. baumannii OXA-51-like enzymes (OXA-51, OXA-69) do not count (Methods, Data extraction). 30 of the 60 isolates carried an acquired carbapenemase (Fig 1B). The paper writes BL for this indicator.",
      source_name = "BL"
    ),
    INOC_LOG10CFU = list(
      description = "Starting inoculum of the time-kill experiment, log10 CFU/mL. Sets the initial bacterial count (10^INOC_LOG10CFU), scales the meropenem Emax as (INOC_LOG10CFU/6)^e_inoc_emax_mero, and at 8 log10 CFU/mL multiplies the colistin/polymyxin B Emax by e_inoc8_emax_col.",
      units = "log10 CFU/mL",
      type = "continuous",
      reference_category = "6 log10 CFU/mL (meropenem Emax reference inoculum, Eq 7)",
      notes = "The intended inocula of the pooled studies were 4, 5, 6, 7 and 8 log10 CFU/mL and served as the typical inoculum (Methods). The 8-log10 indicator I(INOC = 8) of Eq 8 is coded as INOC_LOG10CFU >= 7.5 so that it switches on at the nominal 8 log10 inoculum. The paper's clinical simulations drew the inoculum from a normal distribution with mean 6 and SD 0.25 log10 CFU/mL.",
      source_name = "INOC"
    )
  )

  population <- list(
    species = "in vitro (Acinetobacter baumannii clinical and reference strains)",
    n_subjects = 60L,
    n_studies = 21L,
    disease_state = "53 A. baumannii strains (60 strain-paper instances, called isolates) from 21 published time-kill studies; 57 isolates meropenem-resistant (MIC > 8 mg/L), 14 colistin/polymyxin B-resistant (MIC > 2 mg/L), and 30 carrying an acquired carbapenemase (Fig 1).",
    model_system = "369 experiments: 70 growth controls, 104 meropenem, 46 colistin, 37 polymyxin B, 53 meropenem + colistin and 59 meropenem + polymyxin B. 18 papers ran static time-kill experiments and 5 ran hollow-fibre infection models (Table S1). Mostly cation-adjusted Mueller-Hinton broth at 35-37 degrees C. Counts digitised from the published figures; values below each study's limit of quantification were handled with the M3 method.",
    inoculum_range = "Intended 4 to 8 log10 CFU/mL",
    dose_range = "Static meropenem 0 to 125 x MIC and colistin/polymyxin B 0 to 4 x MIC; dynamic meropenem 0.5 to 8 g q8h regimens and colistin/polymyxin B loading/maintenance regimens (Table S1)",
    regions = "Pooled published studies (Asia, Europe, North America)",
    notes = paste(
      "Fitted in NONMEM 7.6 in three steps (growth controls, then each monotherapy, then the combination), fixing each step's estimates before the next.",
      "Interisolate variability (IIV) is retained on the carrying capacity and on both EC50s; no IIV or carbapenemase effect on the interaction was supported (Table S2).",
      "The residual error is additive on the log10 CFU/mL scale."
    )
  )

  ini({
    # --- Growth (Table 1, Eq 1-2) ---
    # Table 1 prints Kgrowth = 0.33. Read as the log-scale estimate, so that
    # Kgrowth = exp(0.33) = 1.39 1/h: a net growth rate of 0.33 1/h cannot
    # produce the growth seen in the paper's own simulations (Fig 2 controls,
    # Fig 3 colistin-alone regrowth, 0% colistin-alone PTA); 1.39 1/h does.
    lkg <- 0.33; label("Net bacterial growth rate constant, Kgrowth, log scale (log(1/h))") # Table 1: Kgrowth 0.33 (RSE 30.0%), read as log(Kgrowth)
    log10bmax <- 8.41; label("Carrying capacity, Bmax (log10 CFU/mL)") # Table 1: Bmax 8.41 log10 cfu/mL (RSE 0.7%)

    # --- Meropenem kill (Table 1, Eq 3-4, 7, 9) ---
    emax_mero <- 3.47; label("Maximum meropenem kill rate at inoculum 6 log10 CFU/mL, no colistin/polymyxin B (1/h)") # Table 1: EmaxMERO 3.47 (RSE 12.7%)
    kon_mero <- 0.033; label("Attenuation rate of meropenem killing, Katten,MERO (1/h)") # Table 1: Katten,MERO 0.033 (RSE 19.1%)
    e_inoc_emax_mero <- -0.63; label("Power exponent of the inoculum effect on the meropenem Emax, gammaMERO (unitless)") # Table 1: gammaMERO -0.63 (RSE 24.0%)
    e_carbapenemase_emax_mero <- 0.83; label("Multiplicative change in the meropenem Emax with an acquired carbapenemase, thetaBL (unitless)") # Table 1: thetaBL 0.83 (RSE 9.2%)
    lec50_mero <- fixed(log(1)); label("Meropenem EC50 as a multiple of the meropenem MIC, log scale (log(x MIC))") # Table 1: EC50MERO 'Fixed to MIC, if no colistin/polymyxin B'

    # --- Colistin/polymyxin B kill (Table 1, Eq 5-6, 8, 10) ---
    emax_col <- 3.62; label("Maximum colistin/polymyxin B kill rate at inoculum below 8 log10 CFU/mL (1/h)") # Table 1: EmaxCOLPMB 3.62 (RSE 2.7%)
    kon_col <- 0.075; label("Attenuation rate of colistin/polymyxin B killing, Katten,COLPMB (1/h)") # Table 1: Katten,COLPMB 0.075 (RSE 6.0%)
    e_inoc8_emax_col <- 0.005; label("Multiplicative change in the colistin/polymyxin B Emax at inoculum 8 log10 CFU/mL, gammaCOLPMB (unitless)") # Table 1: gammaCOLPMB 0.005 (RSE 11.7%)
    lec50_col <- log(0.33); label("Colistin/polymyxin B EC50 as a multiple of the polymyxin MIC, log scale (log(x MIC))") # Table 1: EC50COLPMB 0.33 'If MIC = 1 mg/L' (RSE 13.2%); linear EC50 = a x MIC (Table S2, model P6)

    # --- General pharmacodynamic interaction (Table 1, Eq 11) ---
    int_emax_col_mero <- 0.29; label("Maximum fractional increase of the meropenem kill by colistin/polymyxin B, Delta_Emax (fraction)") # Table 1: Delta_Emax 0.29 (RSE 27.1%)
    int_ec50_col_mero <- fixed(1); label("Maximum fractional decrease of the meropenem EC50 by colistin/polymyxin B, Delta_EC50 (fraction)") # Table 1: Delta_EC50 1 (fix)
    # Not reported in Table 1 or the supplement (the Table 1 footnote cites a
    # 'Supplementary text' that the published supplement does not contain).
    # Assumed fixed to the MIC: Methods models each EC50 as a x MIC^b with a
    # and b 'estimated or fixed to 1', and Table 1 likewise omits a value for
    # EC50_MERO, which is fixed to the MIC.
    ec50int_col_mero <- fixed(1); label("Colistin/polymyxin B potency for the interaction as a multiple of the polymyxin MIC, EC50_INTER (x MIC)") # Eq 11 EC50_INTER: NOT REPORTED; assumed a = b = 1 (Methods, EC50 = a x MIC^b)

    # --- Interisolate variability (Table 1; exponential, SD reported) ---
    etalog10bmax ~ 0.0121 # Table 1: Bmax interisolate variability SD 0.11; variance 0.11^2
    etalec50_mero ~ 1.2769 # Table 1: EC50MERO interisolate variability SD 1.13; variance 1.13^2
    etalec50_col ~ 1.8225 # Table 1: EC50COLPMB interisolate variability SD 1.35; variance 1.35^2

    # --- Residual error (Table 1) ---
    addSd <- 1.21; label("Additive residual SD on log10 CFU/mL (log10 CFU/mL)") # Table 1: residual variability SD 1.21 (RSE 9.2%), printed as 'Proportional residual error'; encoded additive on log10 CFU/mL
  })

  model({
    # Carrying capacity with exponential interisolate variability on the
    # log10-scale parameter used in Eq 1 (B / 10^BMAX).
    kg <- exp(lkg)
    bmax_i <- log10bmax * exp(etalog10bmax)

    # Drug concentrations. Colistin and polymyxin B share one parameter set;
    # only one of the two is present in an experiment.
    c_mero <- CONC_MEM_MGL
    c_col <- CONC_CST_MGL + CONC_PMB_MGL
    if (CONC_PMB_MGL > 0) {
      mic_col <- MIC_PMB
    } else {
      mic_col <- MIC_CST
    }

    # Inoculum effects (Eq 7-8).
    emax_mero_inoc <- emax_mero * (INOC_LOG10CFU / 6)^e_inoc_emax_mero
    inoc8 <- 0
    if (INOC_LOG10CFU >= 7.5) {
      inoc8 <- 1
    }
    emax_col_inoc <- emax_col * e_inoc8_emax_col^inoc8

    # EC50 as a multiple of the MIC with interisolate variability.
    ec50_mero_i <- exp(lec50_mero + etalec50_mero) * MIC_MEM
    ec50_col_i <- exp(lec50_col + etalec50_col) * mic_col

    # GPDI: saturable functions of the colistin/polymyxin B concentration (Eq 11).
    fint_col <- c_col / (ec50int_col_mero * mic_col + c_col)
    int_emax <- int_emax_col_mero * fint_col
    int_ec50 <- int_ec50_col_mero * fint_col

    # Single-drug kill rates (Eq 9-10); the carbapenemase effect multiplies the
    # meropenem Emax when CARBAPENEMASE_POS = 1.
    k_mero <- emax_mero_inoc * e_carbapenemase_emax_mero^CARBAPENEMASE_POS *
      (1 + int_emax) * (1 - aron_mero) * c_mero /
      (ec50_mero_i * (1 - int_ec50) + c_mero)
    k_col <- emax_col_inoc * (1 - aron_col) * c_col / (ec50_col_i + c_col)

    # Bliss independence rescaled to the larger single-drug Emax (Eq 12-13).
    eref <- max(emax_mero_inoc, emax_col_inoc)
    kkill <- k_mero + k_col - (k_mero / emax_mero_inoc) * (k_col / emax_col_inoc) * eref

    # Attenuation switches on only while the drug is present (Eq 3-6).
    kon_mero_i <- kon_mero * (c_mero > 0)
    kon_col_i <- kon_col * (c_col > 0)

    # Logistic growth with killing on the linear CFU scale (Eq 1).
    d/dt(cfu) <- (kg * (1 - cfu / 10^bmax_i) - kkill) * cfu
    d/dt(aron_mero) <- kon_mero_i * aroff_mero
    d/dt(aroff_mero) <- -kon_mero_i * aroff_mero
    d/dt(aron_col) <- kon_col_i * aroff_col
    d/dt(aroff_col) <- -kon_col_i * aroff_col

    # Initial conditions: inoculum (Eq 2); all bacteria start in the
    # maximum-killing state.
    cfu(0) <- 10^INOC_LOG10CFU
    aron_mero(0) <- 0
    aroff_mero(0) <- 1
    aron_col(0) <- 0
    aroff_col(0) <- 1

    # Observation on log10 CFU/mL. The floor (-30 log10 CFU/mL, far below any
    # quantifiable count) keeps the logarithm finite when solver round-off
    # takes a fully killed population to or below zero.
    log_cfu <- log10(max(cfu, 1e-30))
    log_cfu ~ add(addSd)
  })
}
