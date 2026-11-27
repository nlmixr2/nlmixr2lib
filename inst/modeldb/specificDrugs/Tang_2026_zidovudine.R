Tang_2026_zidovudine <- function() {
  description <- paste(
    "Intracellular transport and phosphorylation model for zidovudine",
    "(AZT, ZDV) in peripheral blood mononuclear cells (PBMCs) of healthy",
    "Chinese adults after a single 600 mg oral dose (Tang 2026). Plasma",
    "AZT is supplied as the time-varying covariate CP_ZIDOVUDINE_UM",
    "(umol/L) and drives first-order uptake into the PBMC AZT pool;",
    "intracellular AZT effluxes back to plasma and is phosphorylated",
    "stepwise and REVERSIBLY to AZT monophosphate (MP), diphosphate (DP)",
    "and triphosphate (TP), each step with its own forward and backward",
    "first-order rate constant. Efflux of parent AZT is the only exit from",
    "the cell. All four intracellular species are carried in ng per 10^6",
    "cells. The paper describes plasma AZT with a two-compartment model",
    "but prints none of its parameters (the intracellular layer was fitted",
    "with each subject's observed plasma profile as the input function),",
    "so no plasma PK layer is encoded here and users supply the plasma",
    "trajectory. Inter-individual variability is the %CV of the six",
    "individual ADAPT II MAP-Bayesian estimates; residual error is not",
    "reported and is fixed to zero.",
    sep = " "
  )
  reference <- paste(
    "Tang F, Ji S, Lai X, Yan T, Peng J. (2026).",
    "The study on the transport and phosphorylation kinetics of",
    "zidovudine in the Chinese population and its significance for",
    "individualized antiretroviral therapy.",
    "Sci Rep 16:18986. doi:10.1038/s41598-026-49905-1.",
    sep = " "
  )
  vignette <- "Tang_2026_zidovudine"
  units <- list(
    time = "h",
    dosing = "(none; intracellular model driven by the plasma zidovudine covariate CP_ZIDOVUDINE_UM, umol/L)",
    concentration = "ng/10^6 cells"
  )

  covariateData <- list(
    CP_ZIDOVUDINE_UM = list(
      description = "Instantaneous total zidovudine plasma concentration, supplied as a time-varying covariate (observed plasma samples or an upstream PK simulation).",
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-varying per event row; interpolate linearly between plasma",
        "samples (rxSolve(..., covsInterpolation = 'linear')). Drives the",
        "uptake term kinf_pbmc * CP_ZIDOVUDINE_UM, whose coefficient Tang",
        "2026 Table 2 reports in ng/10^6 cells/uM/h, so the column must be",
        "in umol/L. Zidovudine MW = 267.2 g/mol (the value used by",
        "Fauchet 2013), so a plasma concentration in ng/mL is divided by",
        "267.2 before being supplied here. Tang 2026 Figure 4 mean Cmax",
        "after 600 mg orally is about 3000 ng/mL (11 umol/L). Set to 0",
        "for drug-free periods."
      ),
      source_name = "X1 / C1 (plasma AZT concentration, Tang 2026 Eqs. 1-3)"
    )
  )

  compartmentData <- list(
    pbmc = list(analyte = "zidovudine", units = "ng/10^6 cells", specimen = "blood cell", verified = TRUE),
    pbmc_zdvmp = list(
      analyte = "zidovudine monophosphate",
      units = "ng/10^6 cells",
      specimen = "blood cell",
      verified = TRUE
    ),
    pbmc_zdvdp = list(
      analyte = "zidovudine diphosphate",
      units = "ng/10^6 cells",
      specimen = "blood cell",
      verified = TRUE
    ),
    pbmc_zdvtp = list(
      analyte = "zidovudine triphosphate",
      units = "ng/10^6 cells",
      specimen = "blood cell",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = 6L,
    n_studies = 1L,
    age_range = "22-26 years",
    weight_range = "49-80 kg",
    height_range = "158-180 cm",
    sex_female_pct = 66.7,
    race_ethnicity = "Chinese (100%)",
    disease_state = "Healthy adult volunteers (HIV-negative).",
    dose_range = "Single 600 mg oral zidovudine (two 300 mg tablets) after a 12-hour overnight fast.",
    regions = "China (Nanhua Hospital, University of South China, Hengyang; enrolled 2024)",
    notes = paste(
      "Methods 'Volunteer recruitment': six healthy volunteers (two male,",
      "four female) aged 22-26 years; Results 'Transportation and",
      "phosphorylation of AZT clinical study': body weight 49-80 kg,",
      "height 158-180 cm. Plasma and PBMC samples pre-dose and at 0.5, 1,",
      "1.5, 2, 2.5, 3, 3.5, 4, 5, 6, 8 and 10 h post-dose; AZT, AZT-MP,",
      "AZT-DP and AZT-TP quantified in PBMCs by HPLC-MS/MS. Individual",
      "intracellular profiles were fitted in ADAPT II by maximum-likelihood",
      "estimation with Bayesian priors, using each subject's plasma AZT",
      "data as the input function. The same paper reports in-vitro Molt-4",
      "and HUVEC uptake / phosphorylation experiments that were not fitted",
      "to a model; they are described in the vignette only."
    )
  )

  ini({
    # Table 2 'Parameter estimates for AZT transportation and phosphorylation
    # in PBMCs': Mean and %CV over the six individual ADAPT II estimates.
    # Paper symbols map to canonical names as kin -> kinf_pbmc, kout ->
    # keff_pbmc, kO-MP -> kmet_zdvmp, kMP-O -> kicv, kMP-DP -> kmet_zdvdp,
    # kDP-MP -> kicv_zdvmp, kDP-TP -> kmet_zdvtp, kTP-DP -> kicv_zdvdp.
    lkinf_pbmc <- log(1.22); label("Plasma-to-PBMC zidovudine uptake coefficient kin (ng/10^6 cells per umol/L per h)") # Table 2: kin = 1.22 ng/10^6 cells/uM/h, %CV 39
    lkeff_pbmc <- log(3.76); label("PBMC-to-plasma zidovudine efflux rate constant kout (1/h)") # Table 2: kout = 3.76 1/h, %CV 72
    lkmet_zdvmp <- log(0.478); label("Phosphorylation rate constant AZT -> AZT-MP, kO-MP (1/h)") # Table 2: kO-MP = 0.478 1/h, %CV 45
    lkicv <- log(0.792); label("Dephosphorylation rate constant AZT-MP -> AZT, kMP-O (1/h)") # Table 2: kMP-O = 0.792 1/h, %CV 62
    lkmet_zdvdp <- log(0.659); label("Phosphorylation rate constant AZT-MP -> AZT-DP, kMP-DP (1/h)") # Table 2: kMP-DP = 0.659 1/h, %CV 81
    lkicv_zdvmp <- log(1.72); label("Dephosphorylation rate constant AZT-DP -> AZT-MP, kDP-MP (1/h)") # Table 2: kDP-MP = 1.72 1/h, %CV 76
    lkmet_zdvtp <- log(3.18); label("Phosphorylation rate constant AZT-DP -> AZT-TP, kDP-TP (1/h)") # Table 2: kDP-TP = 3.18 1/h, %CV 79
    lkicv_zdvdp <- log(4.62); label("Dephosphorylation rate constant AZT-TP -> AZT-DP, kTP-DP (1/h)") # Table 2: kTP-DP = 4.62 1/h, %CV 88

    # IIV: Table 2 %CV converted to log-normal variance, omega^2 = log(1 + CV^2).
    etalkinf_pbmc ~ log(1 + 0.39^2) # Table 2 kin %CV 39
    etalkeff_pbmc ~ log(1 + 0.72^2) # Table 2 kout %CV 72
    etalkmet_zdvmp ~ log(1 + 0.45^2) # Table 2 kO-MP %CV 45
    etalkicv ~ log(1 + 0.62^2) # Table 2 kMP-O %CV 62
    etalkmet_zdvdp ~ log(1 + 0.81^2) # Table 2 kMP-DP %CV 81
    etalkicv_zdvmp ~ log(1 + 0.76^2) # Table 2 kDP-MP %CV 76
    etalkmet_zdvtp ~ log(1 + 0.79^2) # Table 2 kDP-TP %CV 79
    etalkicv_zdvdp ~ log(1 + 0.88^2) # Table 2 kTP-DP %CV 88

    # Residual error: the Results paragraph before Table 2 says the table
    # presents residual error variances, but Table 2 prints none and the
    # ADAPT II variance model is not described.
    propSd_Cpbmc <- fixed(0); label("Proportional residual SD, PBMC zidovudine (fraction); not reported") # not reported in Table 2 or text
    propSd_Cpbmc_zdvmp <- fixed(0); label("Proportional residual SD, PBMC AZT-MP (fraction); not reported") # not reported in Table 2 or text
    propSd_Cpbmc_zdvdp <- fixed(0); label("Proportional residual SD, PBMC AZT-DP (fraction); not reported") # not reported in Table 2 or text
    propSd_Cpbmc_zdvtp <- fixed(0); label("Proportional residual SD, PBMC AZT-TP (fraction); not reported") # not reported in Table 2 or text
  })

  model({
    kinf_pbmc <- exp(lkinf_pbmc + etalkinf_pbmc)
    keff_pbmc <- exp(lkeff_pbmc + etalkeff_pbmc)
    kmet_zdvmp <- exp(lkmet_zdvmp + etalkmet_zdvmp)
    kicv <- exp(lkicv + etalkicv)
    kmet_zdvdp <- exp(lkmet_zdvdp + etalkmet_zdvdp)
    kicv_zdvmp <- exp(lkicv_zdvmp + etalkicv_zdvmp)
    kmet_zdvtp <- exp(lkmet_zdvtp + etalkmet_zdvtp)
    kicv_zdvdp <- exp(lkicv_zdvdp + etalkicv_zdvdp)

    # Figure 1 scheme with first-order transfers throughout (Methods
    # 'Transportation and phosphorylation pathway model'). Uptake does not
    # deplete plasma (the plasma profile is an input). Eq. 2 writes the
    # AZT balance as dX2/dt = kab*X1 - (kba + k0)*X2; the MP -> AZT return
    # term below is the kMP-O arrow of Figure 1 and Table 2.
    d/dt(pbmc) <- kinf_pbmc * CP_ZIDOVUDINE_UM - (keff_pbmc + kmet_zdvmp) * pbmc + kicv * pbmc_zdvmp
    d/dt(pbmc_zdvmp) <- kmet_zdvmp * pbmc - (kicv + kmet_zdvdp) * pbmc_zdvmp + kicv_zdvmp * pbmc_zdvdp
    d/dt(pbmc_zdvdp) <- kmet_zdvdp * pbmc_zdvmp - (kicv_zdvmp + kmet_zdvtp) * pbmc_zdvdp + kicv_zdvdp * pbmc_zdvtp
    d/dt(pbmc_zdvtp) <- kmet_zdvtp * pbmc_zdvdp - kicv_zdvdp * pbmc_zdvtp

    # States are already in ng/10^6 cells (Table 1 units), so each
    # observation equals its state.
    Cpbmc <- pbmc
    Cpbmc_zdvmp <- pbmc_zdvmp
    Cpbmc_zdvdp <- pbmc_zdvdp
    Cpbmc_zdvtp <- pbmc_zdvtp

    Cpbmc ~ prop(propSd_Cpbmc)
    Cpbmc_zdvmp ~ prop(propSd_Cpbmc_zdvmp)
    Cpbmc_zdvdp ~ prop(propSd_Cpbmc_zdvdp)
    Cpbmc_zdvtp ~ prop(propSd_Cpbmc_zdvtp)
  })
}
