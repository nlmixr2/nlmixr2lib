Horita_2018_ethambutol <- function() {
  description <- "Two-compartment population pharmacokinetic model with zero-order absorption (lag time + zero-order duration) and first-order elimination for oral ethambutol in Ghanaian children with active tuberculosis (Horita 2018); allometric weight scaling on CL/F, Q/F, V1/F, V2/F with non-canonical estimated exponents (0.382, 0.474, 0.228, 0.858) normalised to the cohort median 14.3 kg."
  reference <- "Horita Y, Alsultan A, Kwara A, Antwi S, Enimil A, Ortsin A, Dompreh A, Yang H, Wiesner L, Peloquin CA. Evaluation of the Adequacy of WHO Revised Dosages of the First-Line Antituberculosis Drugs in Children with Tuberculosis Using Population Pharmacokinetic Modeling and Simulations. Antimicrob Agents Chemother. 2018;62(9):e00008-18. doi:10.1128/AAC.00008-18"
  vignette <- "Horita_2018_ethambutol"
  units <- list(time = "hour", dosing = "mg", concentration = "ug/mL")

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric scaling with non-canonical estimated exponents on each PK parameter: 0.382 on CL/F (RSE 24%), 0.474 on Q/F (RSE 44%), 0.228 on V1/F (RSE 84%), 0.858 on V2/F (RSE 60%). Reference weight is the cohort median 14.3 kg (Table 1); see the vignette Errata for the reference-weight derivation. The Horita 2018 EMB exponents differ markedly from the canonical 0.75 / 1.0 because the cohort exhibits notably low EMB exposure in younger children with otherwise stable CL/F across age bands (Methods 'EMB' / Discussion paragraph 5); the paper's note 'Fixed allometric scaling exponents of 0.382 for CL/F, 0.474 for Q/F, 0.228 for V1/F, and 0.858 for V2/F improved the goodness-of-fit plots and the distribution of PK parameters' refers to fixing the exponents during covariate model selection -- they are point estimates with reported RSEs.",
      source_name        = "WT"
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 113L,
    n_studies      = 1L,
    age_range      = "3 months to 14 years (median 5.00 years, IQR 2.17 to 8.25)",
    age_median     = "5.00 years",
    weight_range   = "5-30 kg (median 14.3, IQR 9.70 to 20.1)",
    weight_median  = "14.3 kg",
    sex_female_pct = 44.2,
    hiv_positive_pct = 52.2,
    disease_state  = "Ghanaian children with active tuberculosis (HIV-positive and HIV-negative). 21.2% under 2 years of age.",
    dose_range     = "Ethambutol 15-25 mg/kg orally daily (median 16.8 mg/kg, IQR 15.0-20.3). Administered as part of standard four-drug anti-TB regimen during the initial 2-month intensive phase.",
    regions        = "Ghana (Komfo Anokye Teaching Hospital, Kumasi).",
    notes          = "Patients enrolled October 2012-August 2015. PK sampling after at least 4 weeks of anti-TB treatment (steady state). Blood samples at 0, 1, 2, 4, 8 h postdose. EMB concentrations 0.0844-5.46 ug/mL by LC-MS/MS. Three children with only-BLQ values (malabsorption group) were excluded from model building; two of those were under 2 years old. ClinicalTrials.gov NCT01687504. Demographics from Horita 2018 Table 1; structural model and parameters from Table 4."
  )

  ini({
    # Structural PK parameters -- Horita 2018 Table 4 final population pharmacokinetic
    # model for ethambutol. Typical values are at the cohort median weight 14.3 kg.
    ltlag <- log(0.723); label("Absorption lag time Tlag (h)")                                # Table 4: Tlag = 0.723 h (RSE 7%)
    ltk0  <- log(0.9);   label("Duration of zero-order absorption Tk0 (h)")                   # Table 4: Tk0  = 0.9   h (RSE 16%)
    lcl   <- log(32.5);  label("Apparent oral clearance CL/F at WT = 14.3 kg (L/h)")          # Table 4: CL/F = 32.5 L/h (RSE 5%)
    lvc   <- log(112);   label("Apparent central volume V1/F at WT = 14.3 kg (L)")            # Table 4: V1/F = 112 L (RSE 7%)
    lq    <- log(15.4);  label("Apparent inter-compartmental clearance Q/F at WT = 14.3 kg (L/h)")  # Table 4: Q/F  = 15.4 L/h (RSE 10%)
    lvp   <- log(97.8);  label("Apparent peripheral volume V2/F at WT = 14.3 kg (L)")         # Table 4: V2/F = 97.8 L (RSE 8%)

    # Allometric exponents on body weight -- ESTIMATED (paper reports RSE), not the
    # canonical 0.75 / 1.0 theoretical values. The Table 4 IIV column lists 'Fixed'
    # for these rows, which means no IIV is estimated on the exponent -- not that the
    # exponent value itself is held fixed during fitting.
    e_wt_cl <- 0.382; label("Allometric exponent on CL/F (estimated, unitless)")   # Table 4: Exponent (BW on CL/F) = 0.382 (RSE 24%)
    e_wt_q  <- 0.474; label("Allometric exponent on Q/F  (estimated, unitless)")   # Table 4: Exponent (BW on Q/F)  = 0.474 (RSE 44%)
    e_wt_vc <- 0.228; label("Allometric exponent on V1/F (estimated, unitless)")   # Table 4: Exponent (BW on V1/F) = 0.228 (RSE 84%)
    e_wt_vp <- 0.858; label("Allometric exponent on V2/F (estimated, unitless)")   # Table 4: Exponent (BW on V2/F) = 0.858 (RSE 60%)

    # Inter-individual variability. Table 4 IIV column reports 'omega (CV%)' on the
    # log scale (variance = omega^2).
    etaltlag ~ 0.0303   # Table 4: 0.174 (17.5% CV) -- 0.174^2 = 0.0303 (Tlag)
    etaltk0  ~ 0.366    # Table 4: 0.605 (66.5% CV) -- 0.605^2 = 0.366  (Tk0)
    etalcl   ~ 0.210    # Table 4: 0.458 (48.3% CV) -- 0.458^2 = 0.210  (CL/F)
    etalvc   ~ 0.368    # Table 4: 0.607 (66.7% CV) -- 0.607^2 = 0.368  (V1/F)
    etalq    ~ 0.0751   # Table 4: 0.274 (27.9% CV) -- 0.274^2 = 0.0751 (Q/F)
    etalvp   ~ 0.0961   # Table 4: 0.310 (31.8% CV) -- 0.310^2 = 0.0961 (V2/F)

    # Proportional residual error only. Table 4: 'Slope b' = 0.272 (RSE 6%); no
    # additive component (Results 'EMB' paragraph 1: 'the proportional residual
    # error model was selected based on the OFV').
    propSd <- 0.272; label("Proportional residual SD (fraction)")    # Table 4: slope b = 0.272
  })

  model({
    # Individual PK parameters with allometric weight scaling (reference 14.3 kg).
    # Each parameter has its own estimated exponent per Table 4.
    tlag <- exp(ltlag + etaltlag)
    tk0  <- exp(ltk0  + etaltk0)
    cl   <- exp(lcl   + etalcl) * (WT / 14.3)^e_wt_cl
    vc   <- exp(lvc   + etalvc) * (WT / 14.3)^e_wt_vc
    q    <- exp(lq    + etalq)  * (WT / 14.3)^e_wt_q
    vp   <- exp(lvp   + etalvp) * (WT / 14.3)^e_wt_vp

    kel <- cl / vc
    k12 <- q  / vc
    k21 <- q  / vp

    # Two-compartment PK with zero-order absorption directly into central. The
    # dose event must specify cmt = central with rate = -2 so rxode2 uses the
    # model's dur() for the zero-order infusion duration. The alag() delays the
    # start of absorption by Tlag.
    d/dt(central)     <- -kel * central - k12 * central + k21 * peripheral1
    d/dt(peripheral1) <-                  k12 * central - k21 * peripheral1

    dur(central)  <- tk0
    alag(central) <- tlag

    # Concentration: dose mg / V1 L -> mg/L = ug/mL.
    Cc <- central / vc
    Cc ~ prop(propSd)
  })
}
