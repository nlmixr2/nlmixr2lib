Hong_2019_evofosfamide_cellular <- function() {
  description <- paste(
    "QSP. In vitro / preclinical (HCT116 and H460 human carcinoma cells;",
    "nude-mouse xenograft tissue parameterization). Hong 2019 cellular",
    "pharmacokinetic/pharmacodynamic model for the hypoxia-activated prodrug",
    "evofosfamide (TH-302). Ten ODEs carrying four chemical species --",
    "the prodrug TH-302, its cytotoxic bromo-isophosphoramide mustard",
    "metabolite Br-IPM, a single notional intermediate INT standing for the",
    "three chloro-substitution intermediates, and the dichloro product IPM --",
    "each in an extracellular and an intracellular compartment coupled by",
    "first-order membrane influx and efflux. Intracellular bioreductive",
    "activation of TH-302 to Br-IPM is oxygen-inhibited through a hyperbolic",
    "term with KO2 = 0.27 uM, so activation is near-maximal under anoxia and",
    "essentially switched off in well-oxygenated cells. Two alternative",
    "cell-kill readouts are computed side by side from intracellular AUC:",
    "a 'bystander' model driven by Br-IPM + IPM exposure and a 'no bystander'",
    "model driven by TH-302 exposure and its oxygen-dependent rate of",
    "reduction; the paper's conclusion is that the two are nearly equivalent,",
    "so a bystander effect is not needed to explain TH-302 monotherapy",
    "activity. Deterministic: no IIV and no residual error are reported.",
    "SCOPE -- this is the WELL-MIXED (no-gradient) limit of the paper's",
    "equations (1) and (2). The published headline model is spatially",
    "resolved: the same equations are solved with a Laplacian diffusion term",
    "by Green's function methods over digitized R3230Ac and FaDu",
    "microvascular networks, which rxode2 cannot express and whose network",
    "geometry was not available when this model was built. Diffusion coefficients are therefore NOT",
    "parameters of this file. See the vignette for what this scope does and",
    "does not reproduce.",
    sep = " "
  )
  reference <- paste(
    "Hong CR, Wilson WR, Hicks KO (2019). An intratumor",
    "pharmacokinetic/pharmacodynamic model for the hypoxia-activated prodrug",
    "evofosfamide (TH-302): monotherapy activity is not dependent on a",
    "bystander effect. Neoplasia 21(2):159-171.",
    "doi:10.1016/j.neo.2018.11.009. PMCID: PMC6314220.",
    "Model equations from main-text equations (1)-(5); parameter values from",
    "Supplementary Table S1 and the Supplementary Methods section 'Kinetics",
    "of Br-IPM, intermediate and IPM formation'.",
    sep = " "
  )
  vignette <- "Hong_2019_evofosfamide"

  # Every state is a chemical species of the TH-302 activation cascade held in
  # one of two subcellular spaces. None maps onto a canonical PK compartment
  # role: `central` / `peripheral1` describe a body-level disposition model,
  # whereas these are the extracellular medium and the cytosol of one cell
  # population. The canonical `int_tumor` / `is_tumor` PBPK sub-compartment
  # pair carries only ONE species per organ and so cannot hold four.
  # Suffixes: `_ec` extracellular, `_ic` intracellular.
  paper_specific_compartments <- c(
    "th302_ec",
    "th302_ic",
    "bripm_ec",
    "bripm_ic",
    "intm_ec",
    "intm_ic",
    "ipm_ec",
    "ipm_ic",
    "auc_th302_ic",
    "auc_metab_ic"
  )

  units <- list(
    time = "h",
    dosing = paste(
      "uM (the states are CONCENTRATIONS, not amounts; an rxode2 dose record",
      "into th302_ec sets the initial extracellular prodrug concentration, so",
      "amt = 30 means 30 uM TH-302 in the medium as in Figure 2E)",
      sep = " "
    ),
    concentration = "uM"
  )

  covariateData <- list(
    STIM_OXYGEN_UM = list(
      description = paste(
        "Oxygen concentration in the medium / tissue surrounding the cells.",
        "Sets the degree of inhibition of bioreductive TH-302 activation",
        "through equation (3) and, in the 'no bystander' cell-kill model,",
        "scales exposure to cell kill through equation (5).",
        sep = " "
      ),
      units = "uM",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Per-record covariate; constant within an experimental arm. Anchors",
        "used by the source: 0 uM = anoxia (95% N2 gas phase); 180 uM =",
        "humidified 20% O2 in the gas phase (Supplementary Methods,",
        "'Estimation of KO2 for TH-302'); < 1 uM is the source's definition of",
        "the tumor hypoxic fraction. The source reports KO2 in both units",
        "(0.2 mmHg = 0.27 uM); this model works in uM throughout.",
        sep = " "
      ),
      source_name = "[O2]"
    )
  )

  # Issue #482. Every species state holds a CONCENTRATION (uM) rather than an
  # amount, because equations (1) and (2) are written per unit volume of the
  # respective phase. `specimen` is "tumor" for the species states: the cells
  # are HCT116 / H460 carcinoma cells, whether grown as monolayers, as
  # multicellular layers, or as a xenograft microregion. The two AUC states are
  # bookkeeping integrators.
  compartmentData <- list(
    th302_ec = list(analyte = "evofosfamide (TH-302)", units = "uM", specimen = "tumor", verified = TRUE),
    th302_ic = list(analyte = "evofosfamide (TH-302)", units = "uM", specimen = "tumor", verified = TRUE),
    bripm_ec = list(
      analyte = "bromo-isophosphoramide mustard (Br-IPM)",
      units = "uM",
      specimen = "tumor",
      verified = TRUE
    ),
    bripm_ic = list(
      analyte = "bromo-isophosphoramide mustard (Br-IPM)",
      units = "uM",
      specimen = "tumor",
      verified = TRUE
    ),
    intm_ec = list(
      analyte = "notional chloro-substitution intermediate (INT) between Br-IPM and IPM",
      units = "uM",
      specimen = "tumor",
      verified = TRUE
    ),
    intm_ic = list(
      analyte = "notional chloro-substitution intermediate (INT) between Br-IPM and IPM",
      units = "uM",
      specimen = "tumor",
      verified = TRUE
    ),
    ipm_ec = list(
      analyte = "isophosphoramide mustard (IPM)",
      units = "uM",
      specimen = "tumor",
      verified = TRUE
    ),
    ipm_ic = list(
      analyte = "isophosphoramide mustard (IPM)",
      units = "uM",
      specimen = "tumor",
      verified = TRUE
    ),
    auc_th302_ic = list(
      analyte = "cumulative intracellular evofosfamide (TH-302) exposure",
      units = "uM*h",
      specimen = "not applicable",
      verified = TRUE
    ),
    auc_metab_ic = list(
      analyte = "cumulative intracellular Br-IPM + IPM exposure",
      units = "uM*h",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  population <- list(
    species = "in vitro (HCT116 human colon carcinoma cell line; H460 NSCLC cells cross-checked) + mouse (nude-mouse HCT116 and H460 xenograft tissue parameterization)",
    n_subjects = 3L,
    n_studies = 1L,
    disease_state = "HCT116 colon carcinoma and H460 non-small-cell lung carcinoma; hypoxic tumor microenvironment",
    dose_range = "30 uM TH-302 or 100 uM Br-IPM applied to the donor compartment in vitro; 50 mg/kg i.v. TH-302 in nude mice, giving a plasma AUC of 25 uM*h used as the tumor inflow",
    notes = paste(
      "Not a patient population. Transport parameters were measured in",
      "multicellular layers (MCLs) grown from HCT116 cells and in HCT116",
      "monolayer cultures (10^6 cells / 0.5 mL); mean and SE are from 3 MCLs",
      "per condition, with mean MCL thickness 139 +/- 2 um (Figure 2 legend).",
      "The PD parameter AUC10 was estimated from clonogenic survival of",
      "anoxic HCT116 cells exposed to 0.01-0.5 uM TH-302 for 1 h",
      "(Figure S3 legend). The tumor intracellular volume fraction",
      "phicell = 0.45 (Table S1) applies to tumors and MCLs; the monolayer",
      "experiments of Figure 2E have a far smaller cell volume fraction --",
      "see the vignette.",
      sep = " "
    )
  )

  ini({
    # =====================================================================
    # Hong 2019 Supplementary Table S1 lists every parameter. The source
    # reports all rate constants in SECONDS^-1; this file keeps the printed
    # s^-1 values here so each number is a literal source trace, and converts
    # to h^-1 with an explicit * 3600 in model(). AUC10 is printed in uM*h
    # and the model time unit is h, so the PD equations need no conversion.
    #
    # Every parameter is fixed(): none is estimated by this file. The source
    # fitted most of them (to MCL flux, monolayer steady state, or clonogenic
    # survival) and fixed the rest from prior work; either way they enter here
    # as published point estimates of a deterministic mechanistic model.
    #
    # NOT included: the tissue diffusion coefficients D (TH-302 1.82e-7,
    # Br-IPM and IPM 1.33e-7 cm^2/s) and the support-membrane / medium
    # diffusion coefficients. They multiply the Laplacian term of equation (1),
    # which is identically zero in this well-mixed scope, so carrying them
    # would imply a spatial capability this file does not have.
    # =====================================================================

    # --- Volume fraction -------------------------------------------------
    phicell <- fixed(0.45); label("Intracellular volume fraction of the tissue or culture (unitless)") # Table S1 'phi_i' = 0.45, cited to Foehrenbacher 2013

    # --- Oxygen-dependent bioreductive activation of the prodrug ---------
    #     Equation (3): kmet = ko2 / (ko2 + [O2]) * kmet0
    kmet0 <- fixed(0.0115); label("Maximum rate constant for intracellular bioreductive metabolism of TH-302, attained under anoxia (1/s)") # Table S1 'k_met,0' = 0.0115 +/- 0.05 s^-1, from anoxic HCT116 MCL flux
    ko2 <- fixed(0.27); label("Oxygen concentration giving half-maximal inhibition of TH-302 bioreduction (uM)") # Suppl Methods 'Estimation of KO2': 0.27 uM (Table S1 prints the same value as 0.2 mmHg)

    # --- Membrane transfer rate constants (source symbols kin / kout) ----
    #     Named kmemin_/kmemout_ rather than kin/kout because the canonical
    #     nlmixr2lib kin/kout are indirect-response turnover constants
    #     (production into, and loss from, a turnover pool), which is a
    #     different quantity from a plasma-membrane permeability rate.
    kmemin_th302 <- fixed(0.15); label("TH-302 plasma-membrane influx rate constant (1/s)") # Table S1 'k_in' TH-302 = 0.15 s^-1
    kmemout_th302 <- fixed(0.05); label("TH-302 plasma-membrane efflux rate constant (1/s)") # Table S1 'k_out' TH-302 = 0.05 s^-1
    kmemin_bripm <- fixed(0.001); label("Br-IPM plasma-membrane influx rate constant (1/s)") # Table S1 'k_in' Br-IPM = 0.001 s^-1
    kmemout_bripm <- fixed(0.001); label("Br-IPM plasma-membrane efflux rate constant (1/s)") # Table S1 'k_out' Br-IPM = 0.001 s^-1
    kmemin_ipm <- fixed(0.0005); label("IPM plasma-membrane influx rate constant (1/s)") # Table S1 'k_in' IPM = 0.0005 s^-1
    kmemout_ipm <- fixed(0.0005); label("IPM plasma-membrane efflux rate constant (1/s)") # Table S1 'k_out' IPM = 0.0005 s^-1
    kmemin_intm <- fixed(0.0005); label("INT plasma-membrane influx rate constant, assumed equal to IPM (1/s)") # Suppl Methods: INT assumed to share D, k_in and k_out with IPM
    kmemout_intm <- fixed(0.0005); label("INT plasma-membrane efflux rate constant, assumed equal to IPM (1/s)") # Suppl Methods: INT assumed to share D, k_in and k_out with IPM

    # --- Extracellular (medium) chemical conversion, source symbol re -----
    #     Each rate constant is the LOSS of the named species; the product is
    #     the next species in the cascade (TH-302 -> Br-IPM -> INT -> IPM ->
    #     untracked downstream products).
    rec_th302 <- fixed(6.67e-7); label("Extracellular chemical reduction of TH-302 to Br-IPM, anoxia only (1/s)") # Table S1 're' TH-302 = 6.67e-7 s^-1; Fig S1 legend, fitted in anoxic medium without cells
    rec_bripm <- fixed(0.0033); label("Extracellular conversion of Br-IPM to INT (1/s)") # Table S1 're' Br-IPM = 0.0033 s^-1 (see the vignette Errata: Suppl Methods prose misprints this as 0.033)
    rec_intm <- fixed(0.00017); label("Extracellular conversion of INT to IPM (1/s)") # Suppl Methods 'Kinetics of Br-IPM...': r_e,INT = 0.00017 s^-1, from the Fig 2C,D and Fig 3A,B flux data
    rec_ipm <- fixed(6.67e-6); label("Extracellular loss of IPM to downstream products (1/s)") # Table S1 're' IPM = 6.67e-6 s^-1

    # --- Intracellular chemical conversion, source symbol ri --------------
    ric_bripm <- fixed(0.015); label("Intracellular conversion of Br-IPM to INT (1/s)") # Table S1 'ri' Br-IPM = 0.015 s^-1, fitted to monolayer intracellular concentrations
    ric_intm <- fixed(0.01); label("Intracellular conversion of INT to IPM (1/s)") # Suppl Methods 'Kinetics of Br-IPM...': r_i,INT = 0.01 s^-1
    ric_ipm <- fixed(0.015); label("Intracellular loss of IPM to downstream products (1/s)") # Table S1 'ri' IPM = 0.015 s^-1

    # --- Pharmacodynamics: intracellular AUC giving 10% clonogenic survival
    auc10_metab <- fixed(0.771); label("Intracellular Br-IPM + IPM AUC giving 10% clonogenic survival (uM*h)") # Table S1 'AUC10' Br-IPM and IPM (equation 4) = 0.771 uM.h, from Figure S3A
    auc10_th302 <- fixed(0.695); label("Intracellular TH-302 AUC giving 10% clonogenic survival (uM*h)") # Table S1 'AUC10' TH-302 (equation 5) = 0.695 uM.h, from Figure S3B
  })

  model({
    # ===================================================================
    # 0. Unit conversion and labels.
    #    All source rate constants are per SECOND; the model time unit is h.
    # ===================================================================
    sPerH <- 3600 # s/h, converts every published rate constant to the model time unit

    kmemInP <- kmemin_th302 * sPerH
    kmemOutP <- kmemout_th302 * sPerH
    kmemInB <- kmemin_bripm * sPerH
    kmemOutB <- kmemout_bripm * sPerH
    kmemInN <- kmemin_intm * sPerH
    kmemOutN <- kmemout_intm * sPerH
    kmemInI <- kmemin_ipm * sPerH
    kmemOutI <- kmemout_ipm * sPerH

    recP <- rec_th302 * sPerH
    recB <- rec_bripm * sPerH
    recN <- rec_intm * sPerH
    recI <- rec_ipm * sPerH

    ricB <- ric_bripm * sPerH
    ricN <- ric_intm * sPerH
    ricI <- ric_ipm * sPerH

    # ===================================================================
    # 1. Oxygen dependence of bioreductive activation. Equation (3).
    # ===================================================================
    kmet <- ko2 / (ko2 + STIM_OXYGEN_UM) * kmet0 * sPerH

    # Chemical (non-enzymatic) reduction of TH-302 in the medium is, per the
    # Supplementary Methods, "assumed to be zero except under anoxia". Anoxia
    # is the source's 95% N2 / 0% O2 gas phase, so the switch is exact at
    # [O2] = 0 and needs no invented threshold. The constant is tiny
    # (t1/2 ~ 12 days) and matters only over the multi-hour MCL experiments.
    anoxic <- (STIM_OXYGEN_UM <= 0)

    # ===================================================================
    # 2. Phase-volume ratio. Equations (1) and (2) are written per unit
    #    volume of each phase and the membrane flux term is multiplied by
    #    phi_i in BOTH, so dividing equation (1) through by phi_e leaves the
    #    extracellular side scaled by phi_i / phi_e while the intracellular
    #    side is unscaled. That asymmetry is what conserves mass.
    # ===================================================================
    phiRatio <- phicell / (1 - phicell)

    # ===================================================================
    # 3. The cascade. Equation (1) governs each extracellular species and
    #    equation (2) each intracellular species, with the Laplacian term
    #    D * del^2(Ce) dropped (well-mixed scope, see description).
    #
    #    TH-302 --kmet--> Br-IPM --ri--> INT --ri--> IPM --ri--> (sink)   [cells]
    #    TH-302 --re----> Br-IPM --re--> INT --re--> IPM --re--> (sink)   [medium]
    #    with the medium arm of the first step active only under anoxia.
    # ===================================================================
    d/dt(th302_ec) <- -phiRatio * (kmemInP * th302_ec - kmemOutP * th302_ic) -
      anoxic * recP * th302_ec
    d/dt(th302_ic) <- kmemInP * th302_ec - kmemOutP * th302_ic -
      kmet * th302_ic

    d/dt(bripm_ec) <- -phiRatio * (kmemInB * bripm_ec - kmemOutB * bripm_ic) -
      recB * bripm_ec + anoxic * recP * th302_ec
    d/dt(bripm_ic) <- kmemInB * bripm_ec - kmemOutB * bripm_ic -
      ricB * bripm_ic + kmet * th302_ic

    d/dt(intm_ec) <- -phiRatio * (kmemInN * intm_ec - kmemOutN * intm_ic) -
      recN * intm_ec + recB * bripm_ec
    d/dt(intm_ic) <- kmemInN * intm_ec - kmemOutN * intm_ic -
      ricN * intm_ic + ricB * bripm_ic

    d/dt(ipm_ec) <- -phiRatio * (kmemInI * ipm_ec - kmemOutI * ipm_ic) -
      recI * ipm_ec + recN * intm_ec
    d/dt(ipm_ic) <- kmemInI * ipm_ec - kmemOutI * ipm_ic -
      ricI * ipm_ic + ricN * intm_ic

    # ===================================================================
    # 4. Intracellular exposure integrators feeding the PD models.
    #    Equation (4) is driven by Br-IPM + IPM; INT is excluded because the
    #    source fitted AUC10 to "the intracellular concentrations of
    #    Br-IPM + IPM" (Suppl Methods, 'Monolayer metabolism model').
    # ===================================================================
    d/dt(auc_th302_ic) <- th302_ic
    d/dt(auc_metab_ic) <- bripm_ic + ipm_ic

    # ===================================================================
    # 5. Cell kill. Equations (4) and (5). Log cell kill is -log10(SF), so
    #    an intracellular AUC equal to AUC10 gives exactly 1 log of kill,
    #    i.e. 10% clonogenic survival -- which is the definition of AUC10.
    # ===================================================================
    lck_bystander <- auc_metab_ic / auc10_metab # equation (4), 'bystander' model
    lck_nobystander <- ko2 / (ko2 + STIM_OXYGEN_UM) * auc_th302_ic / auc10_th302 # equation (5), 'no bystander' model

    sf_bystander <- 10^(-lck_bystander) # surviving fraction, 'bystander' model
    sf_nobystander <- 10^(-lck_nobystander) # surviving fraction, 'no bystander' model

    # No residual-error model and no IIV. Hong 2019 fits deterministic
    # reaction-diffusion solutions by least squares (MatLab nlinfit) and
    # tabulates no residual-error magnitude or between-experiment variance for
    # any output; the only uncertainties reported are SEs across 3 MCLs.
  })
}
