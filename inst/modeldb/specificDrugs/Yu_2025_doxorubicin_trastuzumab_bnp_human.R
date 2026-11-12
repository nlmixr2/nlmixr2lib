Yu_2025_doxorubicin_trastuzumab_bnp_human <- function() {
  description <- "QSP / PBPK-toxicodynamic clinical scale-up of the AC16 cardiomyocyte-injury biomarker model: relative serum B-type natriuretic peptide (BNP, converted from NT-proBNP) in HER2-positive breast cancer patients on trastuzumab (TmAb), as an indirect-response model with two transit compartments (BNP1 -> BNP2 -> BNP3, turnover kbnp) whose production is stimulated additively by the TmAb concentration in cardiac interstitial fluid and by the doxorubicin (DOX) concentration. All in vitro parameters (kbnp, Sbnp1, Sbnp2) are fixed to the AC16 estimates; only a scaling factor fSct on the TmAb slope was estimated against clinical data. The PBPK exposure layer is NOT included: as in the authors' deposited Monolix code, both drug concentrations are external time-varying inputs (CP_TRASTUZUMAB_NM, CP_DOXORUBICIN_NGML). The paper generated the TmAb input with the Glassman and Balthasar 2016 mAb PBPK model, whose cardiac interstitial profile is not reported, so users must supply their own driver. Output is BNP3, relative BNP normalised to the pre-dose value. No inter-individual variability was estimated (mean data). Sibling models: Yu_2025_doxorubicin_trastuzumab_bnp_ac16 (the in vitro layer this model scales) and Yu_2025_doxorubicin_trastuzumab_ac16."
  reference <- "Yu S, Mody H, Vaidya TR, Kagan L, Ait-Oudhia S. Mitigating Trastuzumab-Doxorubicin Cardiotoxicity With Multiscale Quantitative Systems Toxicology and PBPK-Toxicodynamic Predictive Modeling Framework. CPT Pharmacometrics Syst Pharmacol. 2025;14(10):1625-1636. doi:10.1002/psp4.70087. TmAb cardiac interstitial exposure in the source from Glassman PM, Balthasar JP. J Pharmacokinet Pharmacodyn. 2016;43(4):427-446 (not encoded here)."
  vignette <- "Yu_2025_trastuzumab_doxorubicin_cardiotoxicity"

  units <- list(
    time = "h",
    dosing = "none (both drug concentrations are supplied as time-varying covariates)",
    concentration = "nM (driver concentrations); relative BNP is a dimensionless ratio to its pre-dose value"
  )

  paper_specific_compartments <- c("bnp1", "bnp2", "bnp3")

  compartmentData <- list(
    bnp1 = list(analyte = "BNP transit 1", units = "ratio to pre-dose", specimen = "not applicable", verified = TRUE),
    bnp2 = list(analyte = "BNP transit 2", units = "ratio to pre-dose", specimen = "not applicable", verified = TRUE),
    bnp3 = list(analyte = "BNP", units = "ratio to pre-dose", specimen = "serum", verified = TRUE)
  )

  covariateData <- list(
    CP_TRASTUZUMAB_NM = list(
      description = "Instantaneous trastuzumab concentration in CARDIAC INTERSTITIAL FLUID, supplied as a time-varying driver",
      units = "nM",
      type = "continuous",
      reference_category = NULL,
      notes = "Cardiac interstitial fluid, NOT plasma: the source simulated it with the Glassman and Balthasar 2016 human mAb PBPK model (Methods, Clinical BNP Model Development). Supplying a plasma concentration overstates the driver, because interstitial IgG concentrations sit below plasma. 0 before the first infusion.",
      source_name = "CTmab"
    ),
    CP_DOXORUBICIN_NGML = list(
      description = "Instantaneous doxorubicin concentration, supplied as a time-varying driver",
      units = "ng/mL",
      type = "continuous",
      reference_category = NULL,
      notes = "Carried in the deposited Monolix human model as the regressor Cdox in nM; converted here to nM with the doxorubicin molar mass (543.52 g/mol) because Sbnp1 is in 1/nM. The source does not name the matrix and its clinical datasets (Advani 2017 NT-proBNP; Jacquinot 2018 LVEF) are trastuzumab-only, so it was 0 in every published fit and simulation. Set to 0 for trastuzumab monotherapy.",
      source_name = "Cdox"
    )
  )

  population <- list(
    species = "human",
    n_subjects = 11L,
    n_studies = 1L,
    disease_state = "Breast cancer (HER2-positive) treated with trastuzumab.",
    dose_range = "Trastuzumab 8 mg/kg IV over 90 min, then 6 mg/kg IV over 30-90 min every 3 weeks (two cycles in the development data; 12 months in the qualification data).",
    notes = "Development data: mean NT-proBNP digitised from 11 breast cancer patients treated with trastuzumab for two cycles (Advani 2017, ref 32), sampled at baseline and on Days 1, 2, 3 and 7 of each cycle, converted to BNP with the Taki 2018 linear relationship (ref 33) and normalised to baseline. External qualification: mean LVEF over 12 months of adjuvant trastuzumab and 6 months of follow-up in 163 HER2-positive early breast cancer patients (Jacquinot 2018, ref 35), converted to relative BNP through an NT-proBNP-LVEF regression (Krittayaphong 2008, ref 36). No individual-level data, so no IIV."
  )

  ini({
    lkbnp <- fixed(log(0.202)); label("BNP turnover rate constant kbnp (1/h)") # Table 1 'k bnp (h-1)' 0.20 (fix); Supplement human-model code Kbnp = 0.202
    lsbnp1 <- fixed(log(0.086)); label("Slope for DOX stimulation of BNP release Sbnp1 (1/nM)") # Table 1 'Sbnp1 (nM-1)' 0.086, fixed from the in vitro fit (Supplement: 'Parameters from the in vitro BNP model are fixed')
    lsbnp2 <- fixed(log(0.005)); label("Slope for TmAb stimulation of BNP release Sbnp2 (1/nM)") # Table 1 'Sbnp2 (nM-1)' 0.005, fixed from the in vitro fit
    lfsct <- log(0.86); label("In vitro-to-human scaling factor on the TmAb slope fSct (unitless)") # Table 1 'f Sct' 0.86 (RSE 31%)

    propSd_bnp3 <- 0.63; label("Proportional residual error on relative human BNP (fraction)") # Table 1 human model 'B' 0.63
  })

  model({
    kbnp <- exp(lkbnp)
    sbnp1 <- exp(lsbnp1)
    sbnp2 <- exp(lsbnp2)
    fsct <- exp(lfsct)

    # Doxorubicin molar mass 543.52 g/mol (C27H29NO11; not stated in the
    # source): ng/mL x 1000 / 543.52 = nM, the unit of the Cdox regressor.
    cdox_nm <- CP_DOXORUBICIN_NGML * 1000 / 543.52

    # Supplement human BNP model code: in vitro Eqs 13-15 with the TmAb
    # slope scaled by fSct and both drug concentrations as regressors
    d/dt(bnp1) <- kbnp * (1 + fsct * sbnp2 * CP_TRASTUZUMAB_NM + sbnp1 * cdox_nm) - kbnp * bnp1
    d/dt(bnp2) <- kbnp * (bnp1 - bnp2)
    d/dt(bnp3) <- kbnp * (bnp2 - bnp3)

    bnp1(0) <- 1
    bnp2(0) <- 1
    bnp3(0) <- 1

    bnp3 ~ prop(propSd_bnp3)
  })
}
