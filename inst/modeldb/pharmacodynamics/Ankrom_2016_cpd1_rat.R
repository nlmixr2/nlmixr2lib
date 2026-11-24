Ankrom_2016_cpd1_rat <- function() {
  description <- paste(
    "Preclinical (rat). Translational PK/PD exposure-response model for CPD1,",
    "a selective small-molecule coagulation factor IXa inhibitor (compound 1),",
    "in a male Sprague-Dawley rat arteriovenous-shunt (AVS) thrombosis and",
    "cuticle-bleeding-time (CBT) model (Ankrom 2016 Table 2, CPD1 columns).",
    "PD-only: the steady-state CPD1 plasma concentration is supplied as a",
    "time-varying covariate CP_CPD1_UM (umol/L) rather than from a coupled PK",
    "model, because the authors fitted clot weight and bleeding time directly",
    "against measured plasma concentration. Clot weight (CW, mg) follows a",
    "vehicle-proportional Imax model CW = CWvehicle * (1 - Imax * C / (IC50 +",
    "C)); cuticle bleeding time (CBT, s) follows a linear model CBT =",
    "CBTvehicle + slope * C. A derived fIXa enzyme-occupancy observable EO =",
    "100 * C / (Ki/fu + C) (%) reproduces the paper's translational",
    "target-engagement axis. The sibling Ankrom_2016_apixaban_rat carries the",
    "comparator factor Xa inhibitor apixaban."
  )
  reference <- paste(
    "Ankrom W, Wood HB, Xu J, Geissler W, Bateman T, Chatterjee MS, Feng K,",
    "Metzger JM, Strapps WR, Tadin-Strapps M, Seiffert D, Andre P.",
    "Preclinical and translational evaluation of coagulation factor IXa",
    "as a novel therapeutic target.",
    "Pharmacol Res Perspect. 2016;4(1):e00207.",
    "doi:10.1002/prp2.207"
  )
  vignette <- "Ankrom_2016_factor_ixa_anticoagulation"
  units <- list(
    time = "h",
    dosing = "(none; PD-only model fed by an external CPD1 plasma-concentration covariate)",
    concentration = paste(
      "driving covariate CP_CPD1_UM is in umol/L; outputs are CW (clot",
      "weight, mg) and CBT (cuticle bleeding time, s); derived EO is fIXa",
      "enzyme occupancy in percent"
    )
  )

  covariateData <- list(
    CP_CPD1_UM = list(
      description = paste(
        "Instantaneous (steady-state) plasma concentration of CPD1 (compound",
        "1, a selective factor IXa small-molecule inhibitor) at the time of",
        "each PD observation, supplied as a time-varying covariate. In Ankrom",
        "2016 a bolus-plus-infusion regimen held plasma CPD1 nearly constant",
        "across the AVS/CBT assessment window, and a terminal sample set the",
        "single concentration paired with each CW/CBT measure. Drives the Imax",
        "clot-weight term and the linear bleeding-time term."
      ),
      units = "umol/L",
      type = "continuous",
      reference_category = NULL,
      notes = paste(
        "Time-varying per event row; set to 0 outside the drug-exposure window",
        "(the drug effect then collapses to zero). Reference values observed:",
        "the five CPD1 regimens (C1-C5) achieved total steady-state plasma",
        "concentrations from 0.28 to 153 umol/L (calculated fIXa EO 7-98%)",
        "(Ankrom 2016 Results)."
      ),
      source_name = "C"
    )
  )

  population <- list(
    species = "rat (Sprague-Dawley, male)",
    n_subjects = 50L,
    n_studies = 1L,
    age_range = NA_character_,
    weight_range = "350-400 g at the time of the AVS/CBT procedures",
    sex_female_pct = 0,
    disease_state = "healthy anesthetized rats in an arteriovenous-shunt thrombosis and cuticle-bleeding model",
    dose_range = paste(
      "Five CPD1 bolus-plus-infusion regimens (bolus / infusion): C1",
      "0.06 / 0.09 mg/kg, C2 0.2 / 0.3 mg/kg, C3 0.6 / 0.9 mg/kg, C4",
      "2 / 3 mg/kg, C5 6.5 / 9.8 mg/kg, each >= 10 rats per group",
      "(Ankrom 2016 Methods)."
    ),
    regions = NA_character_,
    notes = paste(
      "Male Sprague-Dawley rats (Charles River) anesthetized with",
      "thiobutabarbital; compounds given as an i.v. bolus followed by an",
      "infusion selected to hold plasma concentrations nearly constant",
      "(within 10% of the mean) across the AVS/CBT window. Fifteen minutes",
      "after dosing an extracorporeal shunt was connected for 15 min and clot",
      "weight recorded; five minutes after shunt placement two cuticles were",
      "cut and bleeding time recorded to a 600 s assay ceiling. CBT was right",
      "censored at 600 s and fit using the likelihood of the datum exceeding",
      "600 s. Models fit with MATLAB R2013a and Monolix; Table 2 reports",
      "typical values and %RSE only (no IIV or residual-error variance)."
    )
  )

  ini({
    # Ankrom 2016 Table 2, CPD1 columns. Clot-weight Imax model and
    # cuticle-bleeding-time linear model; typical values only (the paper
    # reports %RSE of each estimate but no IIV or residual-error variance).

    # ---- Clot weight CW (mg): vehicle-proportional Imax ------------------
    lrbase_CW <- log(121)
    label("Vehicle (baseline) clot weight CWvehicle (mg)") # Table 2 CPD1 CW: CWvehicle = 121 mg (RSE 17%)
    limax <- log(0.936)
    label("Maximal fractional inhibition of clot weight Imax (unitless, 0-1)") # Table 2 CPD1 CW: Imax = 0.936 (RSE 2%)
    lic50 <- log(1.69)
    label("CPD1 concentration giving half-maximal clot-weight inhibition IC50 (umol/L)") # Table 2 CPD1 CW: IC50 = 1.69 umol/L (RSE 17%)

    # ---- Cuticle bleeding time CBT (s): linear ---------------------------
    lrbase_CBT <- log(122)
    label("Vehicle (baseline) cuticle bleeding time CBTvehicle (s)") # Table 2 CPD1 CBT: CBTvehicle = 122 s (RSE 15%)
    lslope_CBT <- log(5.93)
    label("Linear CPD1 slope on bleeding time (s per umol/L)") # Table 2 CPD1 CBT: Slope = 5.93 (RSE 34%)

    # ---- fIXa enzyme-occupancy constants (Table 1, rat CPD1) -------------
    # EO = 100 * C / (Ki/fu + C). Ki is the rat fIXa inhibition constant and
    # fu the rat plasma fraction unbound; both are fixed measured constants.
    lki <- fixed(log(0.022))
    label("Rat fIXa inhibition constant Ki (umol/L)") # Table 1 CPD1 rat fIXa Ki = 22 nmol/L = 0.022 umol/L
    lfu <- fixed(log(0.006))
    label("Rat plasma fraction unbound fu (unitless)") # Table 1 CPD1 % unbound (rat) = 0.6% = 0.006

    # ---- Residual error (NOT reported by Ankrom 2016) --------------------
    # Table 2 tabulates typical values and %RSE only. The paper states model
    # residuals were approximately normally distributed around zero, so an
    # additive error is the natural structure; the magnitudes below are
    # placeholders (documented as non-paper-derived in the vignette Errata).
    # Typical-value predictions do not depend on them.
    addSd_CW <- 10
    label("Additive residual SD on clot weight (mg; placeholder, not reported by Ankrom 2016)")
    addSd_CBT <- 40
    label("Additive residual SD on bleeding time (s; placeholder, not reported by Ankrom 2016)")
  })

  model({
    rbase_CW <- exp(lrbase_CW)
    imax <- exp(limax)
    ic50 <- exp(lic50)
    rbase_CBT <- exp(lrbase_CBT)
    slope_CBT <- exp(lslope_CBT)
    ki <- exp(lki)
    fu <- exp(lfu)

    # fIXa enzyme occupancy (%), Ankrom 2016 Methods / Figure 4 equation:
    # EO = 100 * C / (Ki/fu + C). Derived output; not a fitted endpoint.
    EO <- 100 * CP_CPD1_UM / (ki / fu + CP_CPD1_UM)

    # Clot weight: vehicle-proportional Imax (Ankrom 2016 Results). The
    # printed equation 'CW = CWvehicle * (1 - Imax * C/IC50 + C)' is a
    # transcription typo for the standard fractional-Imax form below (see
    # vignette Errata).
    CW <- rbase_CW * (1 - imax * CP_CPD1_UM / (ic50 + CP_CPD1_UM))

    # Cuticle bleeding time: linear in concentration (Ankrom 2016 Results).
    CBT <- rbase_CBT + slope_CBT * CP_CPD1_UM

    CW ~ add(addSd_CW)
    CBT ~ add(addSd_CBT)
  })
}
