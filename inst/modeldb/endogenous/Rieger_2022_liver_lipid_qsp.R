Rieger_2022_liver_lipid_qsp <- function() {
  description <- paste(
    "QSP. Human hepatocyte lipid metabolism model for non-alcoholic fatty",
    "liver disease (NAFLD) (Rieger 2022): five ODEs for cytosolic fatty acid",
    "(FA), cytosolic triglyceride (TG), endoplasmic-reticulum FA and TG, and",
    "lumped plasma TG, on a 24-h-average time scale (day). Fluxes are NEFA",
    "uptake from adipose with optional cytosolic-FA feedback, de novo",
    "lipogenesis (DNL), cubic mass-action esterification of FA to TG in the",
    "cytosol and ER, cytosolic lipolysis, beta-oxidation with DNL feedback,",
    "saturable VLDL export from the ER, chylomicron input, hepatic uptake of",
    "plasma TG and peripheral lipase clearance. Outputs are liver fat",
    "(volume %) and plasma TG (mM). Deterministic: the authors build virtual",
    "patients by Metropolis-Hastings / acceptance-rejection sampling against",
    "a joint log-normal of liver fat and plasma TG; no IIV or residual error",
    "is reported. Default ini() values are the basal (typical) parameter set,",
    "which is an exact steady state. Five scale_* multipliers (default 1)",
    "carry the interventions: pioglitazone = scale_nefa_uptake 0.725; diet =",
    "scale_chylo 0.8, scale_dnl 0.448, scale_nefa_uptake 0.89.",
    sep = " "
  )
  reference <- paste(
    "Rieger TR, Allen RJ, Musante CJ (2022).",
    "A Quantitative Systems Pharmacology Model of Liver Lipid Metabolism for",
    "Investigation of Non-Alcoholic Fatty Liver Disease.",
    "Front Pharmacol 13:910789. doi:10.3389/fphar.2022.910789.",
    "Equations and parameter values are from the authors' deposited source",
    "code (Rieger TR, Allen RJ, Musante CJ (2022) Source Code for a",
    "Quantitative Pharmacology Model of Liver Lipid Metabolism, v1.1,",
    "doi:10.5281/zenodo.6621096: dxdt.jl, parameters_pluto.csv, util.jl,",
    "derived_parameters.jl, main.jl); basal values agree with Supplementary",
    "Table S1.",
    sep = " "
  )
  vignette <- "Rieger_2022_liver_lipid_qsp"

  paper_specific_compartments <- c("fa_cy", "tg_cy", "fa_er", "tg_er", "tg_p")

  units <- list(
    time = "day",
    dosing = "mmol",
    concentration = "mM"
  )

  compartmentData <- list(
    fa_cy = list(analyte = "fatty acids (hepatocyte cytosol)", units = "mM", specimen = "tissue", verified = TRUE),
    tg_cy = list(analyte = "triglyceride (hepatocyte cytosol)", units = "mM", specimen = "tissue", verified = TRUE),
    fa_er = list(
      analyte = "fatty acids (hepatocyte endoplasmic reticulum)",
      units = "mM",
      specimen = "tissue",
      verified = TRUE
    ),
    tg_er = list(
      analyte = "triglyceride (hepatocyte endoplasmic reticulum)",
      units = "mM",
      specimen = "tissue",
      verified = TRUE
    ),
    tg_p = list(
      analyte = "triglyceride (all circulating lipoproteins)",
      units = "mM",
      specimen = "plasma",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 1900,
    disease_state = paste(
      "General adult population spanning normal liver fat, hyperlipidaemia,",
      "simple NAFLD and NAFLD with hyperlipidaemia; the NAFLD cohort is the",
      "900 virtual patients with liver fat > 5% (Rieger 2022 Sections 2.4,",
      "3.1).",
      sep = " "
    ),
    dose_range = paste(
      "No drug PK. Pioglitazone 45 mg QD for 24 weeks (Belfort 2006) is",
      "represented as a 27.5% reduction of NEFA uptake; a 26-week diet",
      "(Haufe 2013) as -20% chylomicron, -55% DNL and -11% NEFA flux.",
      sep = " "
    ),
    notes = paste(
      "The 1,900 virtual patients are model constructs selected from 500,000",
      "plausible patients so that steady-state liver fat and plasma TG match",
      "a joint log-normal fitted to the digitised individual data of",
      "Kotronen 2007 (simultaneous liver fat % and fasting serum TG; TG",
      "divided by 0.8 to approximate a 24-h mean). Basal parameters are",
      "derived from literature fluxes and steady-state constraints for a",
      "70 kg adult with a 1.5 kg liver eating 2400 kcal/day (35% fat).",
      sep = " "
    )
  )

  ini({
    # All parameters are calibrated, not estimated: the values below are the
    # basal (typical-value, TV) column of parameters_pluto.csv in the
    # deposited source code (Zenodo 10.5281/zenodo.6621096), matching the
    # two-significant-figure Basal Value column of Supplementary Table S1.
    # The plausible-patient search varies parameters 1-14 (klipase_clear to
    # nefa_uptake_flux) between the LV/HV bounds of the same file.

    # --- Plasma TG disposition ---------------------------------------------
    lklipase_clear <- fixed(log(16.801971423166567)); label("TG clearance from plasma by peripheral lipases (1/day)") # parameters_pluto.csv row 1; Table S1 klipase_clear 1.7E+01
    lkuptake_liver_tg <- fixed(log(0.5734207183549842)); label("Hepatic uptake rate constant of plasma TG (1/day)") # parameters_pluto.csv row 5; Table S1 kuptake_liver_tg 5.7E-01
    lvd_tg_p <- fixed(log(4.515)); label("Distribution volume of plasma TG (L)") # parameters_pluto.csv row 15

    # --- Feedback sensitivities (free parameters in the virtual population) -
    sens_nefa_uptake <- fixed(0.001); label("Exponent of cytosolic-FA feedback on NEFA uptake (unitless)") # parameters_pluto.csv row 2; Table S1 sens_nefa_uptake 1.0E-03
    sens_betaox_dnl <- fixed(0.1); label("Exponent of DNL feedback on beta-oxidation (unitless)") # parameters_pluto.csv row 3; Table S1 sens_betaox_dnl 1.0E-01

    # --- Hepatocyte FA / TG handling -----------------------------------------
    lkuptake_er <- fixed(log(929.808692197991)); label("Rate constant of FA transfer from cytosol to ER (1/day)") # parameters_pluto.csv row 4; Table S1 kuptake_er 9.3E+02
    lksynth_cy_tg <- fixed(log(14463.690767524306)); label("Cytosolic esterification rate constant, 3 FA to 1 TG (1/(mM^2*day))") # parameters_pluto.csv row 6; Table S1 ksynth_cy_tg 1.4E+04
    lklipo_cy_tg <- fixed(log(0.27914923181321905)); label("Cytosolic lipolysis rate constant, 1 TG to 3 FA (1/day)") # parameters_pluto.csv row 7; Table S1 klipo_cy_tg 2.8E-01
    lksynth_er_tg <- fixed(log(129140.09613860988)); label("ER esterification rate constant, 3 FA to 1 TG (1/(mM^2*day))") # parameters_pluto.csv row 8; Table S1 ksynth_er_tg 1.3E+05
    lkbetaox <- fixed(log(1667.225980228621)); label("Beta-oxidation rate constant of cytosolic FA (1/day)") # parameters_pluto.csv row 9; Table S1 kbetaox 1.7E+03
    lemax_vldl_prod <- fixed(log(33.62629534)); label("Maximum VLDL-TG export from the ER (mmol/day)") # parameters_pluto.csv row 10; Table S1 emax_vldl_prod 3.4E+01
    lec50_vldl_prod <- fixed(log(26.67240932642487)); label("ER TG concentration at half-maximal VLDL export (mM)") # parameters_pluto.csv row 11; Table S1 ec50_vldl_prod 2.7E+01

    # --- Input fluxes ------------------------------------------------------
    lchylo_basal_flux <- fixed(log(97.62521588946458)); label("Basal chylomicron TG appearance in plasma (mmol/day)") # parameters_pluto.csv row 12; Table S1 chylo_basal_flux 9.8E+01
    ldnl_basal_flux <- fixed(log(9.665448866236932)); label("Basal de novo lipogenesis FA flux (mmol/day)") # parameters_pluto.csv row 13; Table S1 dnl_basal_flux 9.7E+00
    lnefa_uptake_flux <- fixed(log(171.69402414356475)); label("Basal NEFA uptake flux from plasma into hepatocyte cytosol (mmol/day)") # parameters_pluto.csv row 14; Table S1 nefa_uptake_flux 1.7E+02

    # --- Hepatocyte volumes ------------------------------------------------
    lvd_cyt <- fixed(log(0.4962299999999999)); label("Total hepatocyte cytosol volume (L)") # parameters_pluto.csv row 16; derived_parameters.jl 0.7 * 3.4e-9 cm^3 * 139e6 cells/g * 1500 g
    lvd_er <- fixed(log(0.15879359999999998)); label("Total hepatocyte endoplasmic reticulum volume (L)") # parameters_pluto.csv row 17; derived_parameters.jl vd_cyt * 0.16/0.5

    # --- Intervention multipliers (1 = untreated) -----------------------------
    scale_chylo <- fixed(1); label("Multiplier on chylomicron input (unitless)") # parameters_pluto.csv row 18; diet 0.8 (main.jl, Section 2.6)
    scale_dnl <- fixed(1); label("Multiplier on DNL flux (unitless)") # parameters_pluto.csv row 19; diet 0.448 (main.jl, Section 2.6 reports -55%)
    scale_nefa_uptake <- fixed(1); label("Multiplier on NEFA uptake flux (unitless)") # parameters_pluto.csv row 20; pioglitazone 0.725 (main.jl; Section 2.5 reports -28%), diet 0.89
    scale_tg_ester <- fixed(1); label("Multiplier on cytosolic esterification (unitless)") # parameters_pluto.csv row 21; sensitivity analysis only (Section 2.7)
    scale_vldl_prod <- fixed(1); label("Multiplier on VLDL export (unitless)") # parameters_pluto.csv row 22; sensitivity analysis only (Section 2.7)

    # --- Basal concentrations: initial conditions and the FA feedback reference
    lfa_cy_basal <- fixed(log(0.15)); label("Basal cytosolic FA concentration (mM)") # parameters_pluto.csv row 23; Holzhutter and Berndt 2021
    ltg_cy_basal <- fixed(log(58.29015544041452)); label("Basal cytosolic TG concentration (mM)") # parameters_pluto.csv row 24; derived_parameters.jl 5% volume fraction * 0.9 g/mL / 772 g/mol
    lfa_er_basal <- fixed(log(0.15)); label("Basal ER FA concentration (mM)") # parameters_pluto.csv row 25; set equal to cytosol
    ltg_er_basal <- fixed(log(58.29015544041452)); label("Basal ER TG concentration (mM)") # parameters_pluto.csv row 26; set equal to cytosol
    ltg_p_basal <- fixed(log(1.5385)); label("Basal 24-h mean plasma TG concentration (mM)") # parameters_pluto.csv row 27
  })

  model({
    klipase_clear <- exp(lklipase_clear)
    kuptake_liver_tg <- exp(lkuptake_liver_tg)
    vd_tg_p <- exp(lvd_tg_p)
    kuptake_er <- exp(lkuptake_er)
    ksynth_cy_tg <- exp(lksynth_cy_tg)
    klipo_cy_tg <- exp(lklipo_cy_tg)
    ksynth_er_tg <- exp(lksynth_er_tg)
    kbetaox <- exp(lkbetaox)
    emax_vldl_prod <- exp(lemax_vldl_prod)
    ec50_vldl_prod <- exp(lec50_vldl_prod)
    chylo_basal_flux <- exp(lchylo_basal_flux)
    dnl_basal_flux <- exp(ldnl_basal_flux)
    nefa_uptake_flux <- exp(lnefa_uptake_flux)
    vd_cyt <- exp(lvd_cyt)
    vd_er <- exp(lvd_er)
    fa_cy_basal <- exp(lfa_cy_basal)
    tg_cy_basal <- exp(ltg_cy_basal)
    fa_er_basal <- exp(lfa_er_basal)
    tg_er_basal <- exp(ltg_er_basal)
    tg_p_basal <- exp(ltg_p_basal)

    # Cytosolic-FA feedback on NEFA uptake (dxdt.jl). The denominator is a
    # smooth max(fa_cy, fa_cy_basal / 10): a logistic switch with steepness
    # k = 20 that keeps the base of the power term away from zero.
    fa_floor <- fa_cy_basal / 10
    h_switch <- expit(2 * 20 * (fa_cy - fa_floor)) # dxdt.jl h(k, x, minx) = 1 / (1 + exp(-2k(x - minx)))
    fa_denom <- (fa_cy - fa_floor) * h_switch + fa_floor
    nefa_uptake_feedback <- (fa_cy_basal / fa_denom)^sens_nefa_uptake

    # Fluxes (dxdt.jl). Amount fluxes are mmol/day; rate-constant fluxes are
    # concentration rates in the compartment named.
    nefa_uptake <- nefa_uptake_flux * nefa_uptake_feedback * scale_nefa_uptake # mmol FA/day
    er_uptake <- kuptake_er * fa_cy # mM FA (cytosol)/day
    tg_cy_ester <- ksynth_cy_tg * fa_cy^3 * scale_tg_ester # mM FA (cytosol)/day
    tg_cy_lipo <- klipo_cy_tg * tg_cy # mM TG (cytosol)/day
    tg_p_liver_uptake <- kuptake_liver_tg * tg_p # mM TG (plasma)/day
    tg_er_ester <- ksynth_er_tg * fa_er^3 # mM FA (ER)/day
    # ER lipolysis is zero in the source: ATGL is treated as cytosolic.
    tg_er_vldl_prod <- emax_vldl_prod * tg_er / (ec50_vldl_prod + tg_er) * scale_vldl_prod # mmol TG/day
    tg_p_clear <- klipase_clear * tg_p # mM TG (plasma)/day
    chylo_prod <- chylo_basal_flux * scale_chylo # mmol TG/day
    dnl <- dnl_basal_flux * scale_dnl # mmol FA/day
    beta_ox_feedback <- (dnl_basal_flux / dnl)^sens_betaox_dnl
    beta_ox <- kbetaox * fa_cy * beta_ox_feedback # mM FA (cytosol)/day

    d/dt(fa_cy) <- nefa_uptake / vd_cyt - er_uptake - tg_cy_ester + 3 * tg_cy_lipo -
      beta_ox + 3 * tg_p_liver_uptake * vd_tg_p / vd_cyt + dnl / vd_cyt
    d/dt(tg_cy) <- tg_cy_ester / 3 - tg_cy_lipo
    d/dt(fa_er) <- er_uptake * vd_cyt / vd_er - tg_er_ester
    d/dt(tg_er) <- tg_er_ester / 3 - tg_er_vldl_prod / vd_er
    d/dt(tg_p) <- tg_er_vldl_prod / vd_tg_p - tg_p_liver_uptake - tg_p_clear + chylo_prod / vd_tg_p

    fa_cy(0) <- fa_cy_basal
    tg_cy(0) <- tg_cy_basal
    fa_er(0) <- fa_er_basal
    tg_er(0) <- tg_er_basal
    tg_p(0) <- tg_p_basal

    # Liver fat as a volume percentage (util.jl getLTGprct, vector method):
    # TG mass (g) = mM * L * 772 g/mol / 1000; TG volume (L) = mass / 900 g/L;
    # non-TG hepatocyte volume = 0.95 * 0.7089 L (a 1.5 kg liver of
    # 139e6 hepatocytes/g at 3.4e-9 mL each, 5% of it TG at basal).
    liver_tg_vol <- (tg_cy * vd_cyt + tg_er * vd_er) * 772 / 1000 / 900
    liver_fat <- 100 * liver_tg_vol / (liver_tg_vol + 0.95 * 0.7089) # volume %
    tg_plasma <- tg_p # mM
  })
}
