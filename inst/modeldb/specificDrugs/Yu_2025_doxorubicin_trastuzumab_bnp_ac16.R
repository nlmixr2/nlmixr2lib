Yu_2025_doxorubicin_trastuzumab_bnp_ac16 <- function() {
  description <- "In vitro (AC16 human cardiomyocyte cell line). QSP / quantitative systems toxicology cardiomyocyte-injury biomarker model for B-type natriuretic peptide (BNP) released into the culture medium after doxorubicin (DOX), trastuzumab (TmAb) or DOX + TmAb exposure: an indirect-response model with two transit compartments (BNP1 -> BNP2 -> BNP3, turnover kbnp) whose production is stimulated ADDITIVELY by the DOX medium concentration (first-order degradation, kdeg) and by a hypothetical TmAb stimulus compartment that starts at the applied TmAb concentration and decays first-order (ks), capturing the attenuation of the BNP response over time. Output is BNP3, the relative BNP level in the medium normalised to its t = 0 value. No inter-group variability was estimated. Sibling models: Yu_2025_doxorubicin_trastuzumab_ac16 (caspase signalling and cell viability) and Yu_2025_doxorubicin_trastuzumab_bnp_human (clinical scale-up of this layer)."
  reference <- "Yu S, Mody H, Vaidya TR, Kagan L, Ait-Oudhia S. Mitigating Trastuzumab-Doxorubicin Cardiotoxicity With Multiscale Quantitative Systems Toxicology and PBPK-Toxicodynamic Predictive Modeling Framework. CPT Pharmacometrics Syst Pharmacol. 2025;14(10):1625-1636. doi:10.1002/psp4.70087"
  vignette <- "Yu_2025_trastuzumab_doxorubicin_cardiotoxicity"

  units <- list(
    time = "h",
    dosing = "nM (doxorubicin concentration applied to the culture medium at t = 0)",
    concentration = "nM (medium doxorubicin concentration); relative BNP is a dimensionless ratio to its t = 0 value"
  )

  # bnp1-bnp3 are the paper's Figure 1 BNP compartments; bnp_tmab is the
  # paper's hypothetical BNPTmAb stimulus compartment (Eq 12); conc_dox
  # follows the in-vitro medium-concentration precedent of
  # Mody_2023_doxorubicin_dexrazoxane_jimt1.R from the same group.
  paper_specific_compartments <- c("conc_dox", "bnp_tmab", "bnp1", "bnp2", "bnp3")

  compartmentData <- list(
    conc_dox = list(analyte = "doxorubicin", units = "nM", specimen = "administration site", verified = TRUE),
    bnp_tmab = list(
      analyte = "hypothetical time-attenuated trastuzumab stimulus of BNP release",
      units = "nM",
      specimen = "not applicable",
      verified = TRUE
    ),
    bnp1 = list(analyte = "BNP transit 1", units = "ratio to t = 0", specimen = "not applicable", verified = TRUE),
    bnp2 = list(analyte = "BNP transit 2", units = "ratio to t = 0", specimen = "not applicable", verified = TRUE),
    bnp3 = list(
      analyte = "BNP (cell-culture supernatant)",
      units = "ratio to t = 0",
      specimen = "not applicable",
      verified = TRUE
    )
  )

  covariateData <- list(
    CONC_TRASTUZUMAB_NM = list(
      description = "Static trastuzumab concentration applied to the AC16 culture medium",
      units = "nM",
      type = "continuous",
      reference_category = NULL,
      notes = "Read once, as the t = 0 initial condition of the decaying bnp_tmab stimulus compartment (Eq 12, BNPTmAb(0) = cTmAb_0). 0 for control and DOX-only wells; 23,000 nM (23 uM) in the TmAb-only and DOX + TmAb BNP experiments.",
      source_name = "CTmabI"
    )
  )

  population <- list(
    species = "in vitro (AC16 human cardiomyocyte cell line)",
    n_subjects = NA_integer_,
    n_studies = 1L,
    disease_state = "Immortalised HER2-positive human cardiomyocytes (AC16) used as an in vitro model of anthracycline / anti-HER2 cardiotoxicity.",
    dose_range = "Control, DOX 25 nM, TmAb 23 uM and DOX 25 nM + TmAb 23 uM; supernatant BNP sampled over 72 h.",
    cell_line = "AC16 (MilliporeSigma-EMD Millipore); HER2-positive immortalised human cardiomyocytes.",
    culture = "DMEM with 10% FBS, 1% sodium bicarbonate, 1% MEM non-essential amino acids and 1% penicillin/streptomycin; 37 C, humidified 5% CO2.",
    notes = "BNP measured by ELISA in the culture supernatant and normalised to baseline. Fitted in Monolix 2021R1."
  )

  ini({
    lkdeg_dox <- fixed(log(0.022)); label("First-order degradation rate constant of DOX in medium kdeg (1/h)") # Table 1 'k deg' 0.022 (fix), from Vaidya 2021 (ref 27); Supplement code Kdeg = 0.022
    lks <- log(1); label("First-order decay rate constant of the TmAb stimulus of BNP release ks (1/h)") # Table 1 'k s (h-1)' 1 (RSE 0.1%)
    lkbnp <- fixed(log(0.202)); label("BNP turnover rate constant kbnp (1/h)") # Table 1 'k bnp (h-1)' 0.20 (fix) from ref 27; Supplement code Kbnp = 0.202
    lsbnp1 <- log(0.086); label("Slope for DOX stimulation of BNP release Sbnp1 (1/nM)") # Table 1 'Sbnp1 (nM-1)' 0.086 (RSE 23%)
    lsbnp2 <- log(0.005); label("Slope for TmAb stimulation of BNP release Sbnp2 (1/nM)") # Table 1 'Sbnp2 (nM-1)' 0.005 (RSE 17%)

    propSd_bnp3 <- 0.44; label("Proportional residual error on relative BNP (fraction)") # Table 1 BNP model 'B' 0.44
  })

  model({
    kdeg_dox <- exp(lkdeg_dox)
    ks <- exp(lks)
    kbnp <- exp(lkbnp)
    sbnp1 <- exp(lsbnp1)
    sbnp2 <- exp(lsbnp2)

    # Eq 1 (DOX in medium): bolus "dose" in nM at t = 0
    d/dt(conc_dox) <- -kdeg_dox * conc_dox

    # Eq 12: hypothetical TmAb stimulus compartment, starts at the applied
    # TmAb concentration and decays first-order
    d/dt(bnp_tmab) <- -ks * bnp_tmab
    bnp_tmab(0) <- CONC_TRASTUZUMAB_NM

    # Eqs 13-15: indirect response with two transit compartments; additive
    # DOX and TmAb stimulation of production
    d/dt(bnp1) <- kbnp * (1 + sbnp2 * bnp_tmab + sbnp1 * conc_dox) - kbnp * bnp1
    d/dt(bnp2) <- kbnp * (bnp1 - bnp2)
    d/dt(bnp3) <- kbnp * (bnp2 - bnp3)

    bnp1(0) <- 1
    bnp2(0) <- 1
    bnp3(0) <- 1

    bnp3 ~ prop(propSd_bnp3)
  })
}
