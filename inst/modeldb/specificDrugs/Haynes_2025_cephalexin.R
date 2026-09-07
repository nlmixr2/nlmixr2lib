Haynes_2025_cephalexin <- function() {
  description <- "One-compartment oral population PK model for cephalexin in neonates and young infants 7-60 days old (Haynes 2025; NCT04916951, Children's Hospital Colorado). First-order absorption with an absorption lag time and linear elimination. Apparent clearance (CL/F) and apparent volume of distribution (Vd/F) are scaled to body weight with exponents fixed to 0.75 and 1.0 respectively (reference 70 kg), with an estimated power effect of postmenstrual age on CL/F (reference 41.06 weeks) and of postnatal age on the absorption rate constant (reference 29.14 days). Inter-individual variability is estimated on lag time, absorption rate constant and CL/F but not on Vd/F; residual error is combined additive plus proportional."
  reference <- "Haynes AS, Wei Z, Scheetz MH, Gonzalez D, Messacar K, Tang Girdwood S, Peloquin CA, Fish DN, Anderson P. Oral Cephalexin Population Pharmacokinetics and Target Attainment Analysis in Infants 7-60 Days Old. J Pediatric Infect Dis Soc. 2025;14(10):piaf088. doi:10.1093/jpids/piaf088"
  vignette <- "Haynes_2025_cephalexin"
  units <- list(time = "h", dosing = "mg", concentration = "mg/L")

  # Issue #482: what each ODE state holds, in what amount units, in what
  # biological matrix.
  compartmentData <- list(
    depot = list(
      analyte = "cephalexin", units = "mg",
      specimen = "administration site", verified = TRUE
    ),
    central = list(
      analyte = "cephalexin", units = "mg",
      specimen = "plasma", verified = TRUE
    )
  )

  covariateData <- list(
    WT = list(
      description        = "Body weight",
      units              = "kg",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Allometric scaling on CL/F with the exponent fixed to 0.75 and linear scaling on Vd/F with the exponent fixed to 1.0, both normalized to a 70 kg reference weight (Haynes 2025 Results 'Final Model' and Table 2 rows 'beta WT_CL' / 'beta WT_V'). The study cohort weighed 2.20-5.39 kg (median 3.36 kg; Table 1), so the 70 kg reference is far outside the observed range and is a normalization constant only.",
      source_name        = "WT"
    ),
    PNA = list(
      description        = "Postnatal age (chronological age since birth)",
      units              = "months",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying. Drives an estimated power effect on the absorption rate constant: ka = ka_pop * (PNA_days / 29.14)^0.92 (Haynes 2025 Results 'Final Model' equation 2). The canonical PNA column is carried in MONTHS per inst/references/covariate-columns.md, whereas the source equation is written in days, so model() recovers days as PNA * 30.4375 before forming the age ratio. This is the same reparameterisation used by Zhao_2018_omeprazole.R (days) and Bardhi_2026_ampicillin_foal.R (hours). The 29.14-day reference corresponds to 0.9574 months.",
      source_name        = "PNA"
    ),
    PAGE = list(
      description        = "Postmenstrual age (gestational age at birth plus postnatal age)",
      units              = "weeks",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Time-varying, carried in WEEKS rather than the register-default months because the source equation and its 41.06-week reference constant are written in weeks (permitted by the PAGE register entry; same convention as Germovsek_2018_meropenem.R and Riccobene_2017_ceftaroline.R). Drives an estimated power effect on CL/F: CL/F = CL/F_pop * (WT/70)^0.75 * (PMA / 41.06)^2.92 (Haynes 2025 Results 'Final Model' equation 4). PMA was selected over eGFR for the final model because it is the more physiologically appropriate marker of renal maturation and eGFR added no predictive value beyond PMA.",
      source_name        = "PMA"
    )
  )

  # Screened during covariate selection (stepwise forward addition p < 0.05,
  # backward elimination p < 0.01; Haynes 2025 Methods 'Model Development' and
  # Table S3) but NOT retained in the final model. Documented so downstream
  # users know what was tested. Also screened and not retained, but without a
  # canonical covariate column: IV beta-lactam co-administration, feeding
  # category (continuous vs intermittent feeds), route of administration
  # (oral vs gastric vs post-pyloric tube) and blood sampling source
  # (capillary vs venous vs arterial).
  covariatesDataExcluded <- list(
    GA = list(
      description        = "Gestational age at birth",
      units              = "weeks",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened as a covariate and not retained in the final model; postmenstrual age (PAGE) carries the maturational signal instead. GA is still required to CONSTRUCT PAGE (PAGE = GA + postnatal age in weeks), and the dosing simulations stratified the virtual population by GA (30-34 weeks vs >= 35 weeks; Haynes 2025 Methods 'Dosing Simulations' and Table 3). Cohort range 29 3/7 to 40 6/7 weeks, median 37 2/7 weeks (Table 1)."
    ),
    CRCL = list(
      description        = "Estimated glomerular filtration rate (BSA-normalized)",
      units              = "mL/min/1.73 m^2",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened on CL/F and SIGNIFICANT: eGFR and PMA gave similar reductions in objective function value, but PMA was selected for the final model as the more physiologically appropriate marker of renal maturation, and eGFR provided no additional predictive value beyond PMA (Haynes 2025 Results 'Final Model'). No eGFR coefficient is reported for the final model, so none is encoded here. Cohort median 68.8 mL/min/1.73 m^2 (range 37.8-126.1; Table 1), available within 48 h of dosing for 25 of 33 subjects."
    ),
    CREAT = list(
      description        = "Serum creatinine",
      units              = "mg/dL",
      type               = "continuous",
      reference_category = NULL,
      notes              = "Screened as a covariate and not retained; the eGFR derived from it was also screened and not retained (see CRCL). Only clinically obtained creatinine results were used, so 8 of 33 subjects had no creatinine within 48 h of dosing (Haynes 2025 Methods 'Study design' and Table 1). Units not stated in the source; mg/dL is the conventional US reporting unit for the study site."
    )
  )

  population <- list(
    species        = "human",
    n_subjects     = 33L,
    n_studies      = 1L,
    n_observations = 144L,
    age_range      = "9.49-56.54 days postnatal",
    age_median     = "31.16 days postnatal",
    ga_range       = "29 3/7 to 40 6/7 weeks gestational age at birth",
    ga_median      = "37 2/7 weeks gestational age at birth",
    weight_range   = "2.20-5.39 kg",
    weight_median  = "3.36 kg",
    sex_female_pct = 30,
    race_ethnicity = c(
      White = 70, Black = 0, Asian = 9,
      `Native Hawaiian/Other Pacific Islander` = 0,
      `Unknown or not reported` = 21,
      Hispanic = 30, `Non-Hispanic` = 58,
      `Ethnicity unknown or not reported` = 12
    ),
    disease_state  = "Hospitalized neonates and young infants receiving antibiotics, enrolled either while receiving enteral cephalexin as standard of care (9 of 33) or while receiving IV antibiotics and given a single 25 mg/kg enteral research dose of cephalexin (24 of 33). Common indications include bacteremic pyelonephritis and other Enterobacterales or MSSA infections.",
    renal_function = "eGFR median 68.8 mL/min/1.73 m^2 (range 37.8-126.1); available within 48 h of dosing for 25 of 33 subjects",
    dose_range     = "12.0-29.5 mg/kg per dose as an oral suspension (median 24.7 mg/kg); one subject received 12.0 mg/kg and all others 22.5-29.5 mg/kg. Administered by mouth (10 of 33), via nasogastric / orogastric / gastric tube (20 of 33) or via post-pyloric tube (3 of 33).",
    regions        = "United States (single center: Children's Hospital Colorado, Aurora, CO)",
    notes          = "Baseline demographics, dosing, sampling and laboratory availability are in Haynes 2025 Table 1 (with Tables S1-S2). 33 subjects contributed 144 plasma concentrations after data cleaning. Sampling was 3-5 samples per subject across 1-5 dosing intervals. NOTE: the earlier IDWeek 2024 conference abstract of this work (Open Forum Infect Dis 2025;12(Suppl 1):S781, abstract P-1222, doi:10.1093/ofid/ofae631.1404) reported a preliminary analysis of 27 subjects / 114 concentrations with eGFR (not PMA) on CL/F and 15% protein binding; the values encoded here are the FINAL peer-reviewed estimates."
  )

  ini({
    # Structural parameters. Vd/F and CL/F are reported normalized to a 70 kg
    # reference weight, so exp(lvc) and exp(lcl) are the 70 kg values and the
    # allometric terms in model() scale them down to infant weights.
    ltlag <- log(0.61); label("Absorption lag time (h)")  # Haynes 2025 Table 2, 'Lag time (Tlag, h)' = 0.61 (RSE 21.0%, 95% CI 0.41-0.91)
    lka <- log(1.32); label("Absorption rate constant at the reference postnatal age (1/h)")  # Haynes 2025 Table 2, 'Absorption rate constant (Ka, 1/h)' = 1.32 (RSE 21.2%, 95% CI 0.88-1.97)
    lvc <- log(30.63); label("Apparent volume of distribution Vd/F (L/70 kg)")  # Haynes 2025 Table 2, 'Apparent volume of distribution (Vd/F, L/70kg)' = 30.63 (RSE 5.64%, 95% CI 27.43-34.20)
    lcl <- log(6.33); label("Apparent clearance CL/F (L/h/70 kg)")  # Haynes 2025 Table 2, 'Apparent clearance (Cl/F, L/h/70kg)' = 6.33 (RSE 4.49%, 95% CI 5.80-6.91)

    # Covariate effects. Both weight exponents were held fixed by the authors
    # (Results 'Final Model': "exponent fixed to 0.75" / "linear scaling ...
    # exponent fixed to 1.0"), and Table 2 reports them with no RSE or CI,
    # confirming they were not estimated.
    e_pna_ka <- 0.92; label("Power exponent of postnatal age on ka (unitless)")  # Haynes 2025 Table 2, 'beta PNA_Ka' = 0.92 (RSE 33.0%, 95% CI 0.32-1.52)
    e_page_cl <- 2.92; label("Power exponent of postmenstrual age on CL/F (unitless)")  # Haynes 2025 Table 2, 'beta PMA_CL' = 2.92 (RSE 19.7%, 95% CI 1.79-4.05)
    e_wt_cl <- fixed(0.75); label("Allometric exponent of body weight on CL/F (unitless)")  # Haynes 2025 Table 2, 'beta WT_CL' = 0.75, no RSE reported; Results 'Final Model' states the exponent was fixed
    e_wt_vc <- fixed(1.0); label("Allometric exponent of body weight on Vd/F (unitless)")  # Haynes 2025 Table 2, 'beta WT_V' = 1.00, no RSE reported; Results 'Final Model' states linear scaling with the exponent fixed to 1.0

    # Inter-individual variability. Haynes 2025 Table 2 reports these under
    # "Standard Deviation of the Random Effects" as a 'Value' column alongside
    # a 'C.V.(%)' column; the Value column is the log-scale omega SD, which the
    # CV column confirms via CV = sqrt(exp(omega^2) - 1):
    #   Tlag  omega 0.74 -> CV 85.4% (printed 85.35%)
    #   Ka    omega 0.63 -> CV 70.2% (printed 70.19%)
    #   Cl/F  omega 0.20 -> CV 20.2% (printed 20.66%)
    # nlmixr2 expects VARIANCES on the diagonal, so each value below is the
    # printed SD squared. No IIV is encoded on Vd/F: the authors did not
    # estimate one "due to high eta shrinkage and imprecise eta estimation"
    # (Results 'Final Model').
    etaltlag ~ 0.5476; label("IIV on absorption lag time (variance, log scale)")  # Haynes 2025 Table 2, IIV for Tlag SD 0.74 (RSE 20.9%, 95% CI 0.50-1.10); 0.74^2 = 0.5476
    etalka ~ 0.3969; label("IIV on absorption rate constant (variance, log scale)")  # Haynes 2025 Table 2, IIV for Ka SD 0.63 (RSE 24.2%, 95% CI 0.40-0.99); 0.63^2 = 0.3969
    etalcl ~ 0.04; label("IIV on apparent clearance (variance, log scale)")  # Haynes 2025 Table 2, IIV for Cl/F SD 0.20 (RSE 18.5%, 95% CI 0.14-0.29); 0.20^2 = 0.04

    # Residual error: combined constant plus proportional (Results 'Final Model').
    addSd <- 0.26; label("Additive residual error (mg/L)")  # Haynes 2025 Table 2, 'constant (a, mg/L)' = 0.26 (RSE 41.8%, 95% CI 0.12-0.54)
    propSd <- 0.22; label("Proportional residual error (fraction)")  # Haynes 2025 Table 2, 'proportional (b)' = 0.22 (RSE 13.0%, 95% CI 0.17-0.28)
  })

  model({
    # ------------------------------------------------------------------
    # Individual parameters. Haynes 2025 Results 'Final Model' prints the
    # four covariate-adjusted equations reproduced here one-for-one:
    #   Tlag_i = Tlag_pop * exp(eta_Tlag)
    #   Ka_i   = Ka_pop * (PNA / 29.14)^beta_PNA * exp(eta_Ka)
    #   Vd_i/F = Vd_pop/F * (WT / 70)^1
    #   CL_i/F = CL_pop/F * (WT / 70)^0.75 * (PMA / 41.06)^beta_PMA * exp(eta_CL)
    # ------------------------------------------------------------------
    tlag <- exp(ltlag + etaltlag)

    # The canonical PNA column is carried in months; the source equation and
    # its 29.14 reference are in DAYS, so recover days here (1 month =
    # 30.4375 days) rather than restating the paper's constant on a different
    # scale. Matches Zhao_2018_omeprazole.R / Bardhi_2026_ampicillin_foal.R.
    pnaDays <- PNA * 30.4375
    ka <- exp(lka + etalka) * (pnaDays / 29.14)^e_pna_ka

    # PAGE is carried in weeks for this model, matching the paper's 41.06-week
    # reference (see covariateData$PAGE$notes).
    vc <- exp(lvc) * (WT / 70)^e_wt_vc
    cl <- exp(lcl + etalcl) * (WT / 70)^e_wt_cl * (PAGE / 41.06)^e_page_cl

    kel <- cl / vc

    # ------------------------------------------------------------------
    # One-compartment disposition with first-order oral absorption.
    # ------------------------------------------------------------------
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - kel * central

    # Absorption lag applied to the depot: dosing starts tlag hours after the
    # dose record.
    alag(depot) <- tlag

    # ------------------------------------------------------------------
    # Observation. Dose in mg and vc in L give mg/L. CL/F and Vd/F are
    # APPARENT (oral) parameters, so no separate bioavailability term is
    # identifiable and none is encoded.
    # ------------------------------------------------------------------
    Cc <- central / vc
    Cc ~ add(addSd) + prop(propSd)
  })
}
