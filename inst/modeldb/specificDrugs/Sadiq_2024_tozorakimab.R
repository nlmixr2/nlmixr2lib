Sadiq_2024_tozorakimab <- function() {
  description <- "Two-compartment population PK / target-engagement model for tozorakimab (anti-IL-33 IgG1 monoclonal antibody) with explicit central-compartment binding of IL-33 to tozorakimab and to the soluble decoy receptor sST2, in healthy adults and patients with mild COPD"
  reference <- "Sadiq MW, Yu H, Astrand M, et al. Population pharmacokinetic/target engagement modelling of tozorakimab in healthy volunteers and patients with chronic obstructive pulmonary disease. Br J Clin Pharmacol. 2024;90(12):3286-3295. doi:10.1111/bcp.16195"
  vignette <- "Sadiq_2024_tozorakimab"
  units <- list(time = "day", dosing = "nmol", concentration = "nM")

  # The model is molar throughout: every target-engagement parameter in Table 2
  # carries molar units (Kon in 1/nM/day, Kd and the baselines in nM), so the
  # free-tozorakimab state that multiplies Kon must also be nM. The paper
  # reports doses in mg and concentrations in ug/mL (tozorakimab, Figure 2) or
  # pg/mL (both complexes, Figures 3 and 4) but never reports the molar masses
  # needed to convert; see the vignette "Assumptions and deviations" section.

  # `sst2` (free soluble ST2 decoy receptor) and `complex_sst2` (the IL-33/sST2
  # complex) have no canonical role in compartment-names.md. They are declared
  # here per the documented paper-specific-compartment mechanism.
  paper_specific_compartments <- c("sst2", "complex_sst2")

  covariateData <- list()

  # Sex, age, race, body weight and the study population (healthy volunteer vs
  # COPD) were all screened by stepwise covariate modelling in PsN, and none was
  # retained: "During the stepwise covariate modelling analysis, no significant
  # effects of covariates on the PK parameters were identified" (Section 3.3).
  # They are documented here rather than in covariateData because the final
  # model references none of them.
  covariatesDataExcluded <- list(
    WT = list(
      description = "Body weight",
      units = "kg",
      type = "continuous",
      notes = "Screened on CL and Vc by stepwise covariate modelling (Section 2.2); not retained in the final model (Section 3.3). Baseline means by cohort are in Table 1."
    ),
    AGE = list(
      description = "Subject age",
      units = "years",
      type = "continuous",
      notes = "Screened on CL and Vc (Section 2.2); not retained (Section 3.3). Table 1 does not tabulate age."
    ),
    SEXF = list(
      description = "Female sex indicator",
      units = "(binary)",
      type = "categorical",
      reference_category = "0 (male)",
      notes = "Screened on CL and Vc (Section 2.2); not retained (Section 3.3). Table 1 reports male sex percentages by cohort."
    ),
    RACE_ASIAN = list(
      description = "Asian race indicator",
      units = "(binary)",
      type = "categorical",
      reference_category = "0 (not Asian)",
      notes = "Race was screened on CL and Vc (Section 2.2); not retained (Section 3.3). Table 1 reports Asian / Black or African American / White / Other."
    ),
    RACE_BLACK = list(
      description = "Black or African American race indicator",
      units = "(binary)",
      type = "categorical",
      reference_category = "0 (not Black or African American)",
      notes = "Race was screened on CL and Vc (Section 2.2); not retained (Section 3.3)."
    ),
    DIS_COPD = list(
      description = "Chronic obstructive pulmonary disease indicator (GOLD grade I-II COPD patient in the MAD cohorts vs healthy volunteer in the SAD / Japanese cohorts)",
      units = "(binary)",
      type = "categorical",
      reference_category = "0 (healthy volunteer)",
      notes = "The study population was screened as a covariate on CL and Vc (Sections 2.2 and 3.3) and was not retained. Documentation only -- no canonical covariate column has been ratified for this concept, and the final model does not reference it."
    )
  )

  compartmentData <- list(
    depot = list(analyte = "tozorakimab", units = "nmol", specimen = "administration site", verified = TRUE),
    central = list(analyte = "tozorakimab", units = "nmol", specimen = "serum", verified = TRUE),
    peripheral1 = list(analyte = "tozorakimab", units = "nmol", specimen = "tissue", verified = TRUE),
    target = list(analyte = "interleukin-33 (free)", units = "nM", specimen = "serum", verified = TRUE),
    sst2 = list(analyte = "soluble ST2 (free)", units = "nM", specimen = "serum", verified = TRUE),
    complex = list(analyte = "IL-33/tozorakimab complex", units = "nM", specimen = "serum", verified = TRUE),
    complex_sst2 = list(analyte = "IL-33/sST2 complex", units = "nM", specimen = "serum", verified = TRUE)
  )

  population <- list(
    species = "human",
    n_subjects = 60,
    n_studies = 1,
    disease_state = "healthy adults with mild atopy and house dust mite sensitivity (SAD and Japanese cohorts); adults with GOLD grade I-II chronic obstructive pulmonary disease (MAD cohorts)",
    weight_median = "cohort means 77.4 kg (SAD), 78.4 kg (MAD), 70.4 kg (Japanese)",
    sex_female_pct = 13.6,
    race_ethnicity = c(White = 72.7, Asian = 18.2, Other = 7.6, Black = 1.5),
    dose_range = "SAD: 1, 3, 10, 30, 100 or 300 mg SC single dose, or 300 mg IV single dose; MAD: 30 or 300 mg SC once every 2 weeks on days 1, 15 and 29; Japanese cohort: 300 mg IV single dose",
    regions = "single first-in-human study (NCT03096795) approved by the London Hampstead Research Ethics Committee; includes a dedicated healthy Japanese cohort",
    notes = paste(
      "Baseline demographics are Table 1, which tabulates 66 tozorakimab-treated participants (42 SAD, 18 MAD, 6 Japanese).",
      "The analysis data set is 60 tozorakimab-treated participants (Abstract and Section 2.1) because the 100 mg MAD cohort was excluded:",
      "its PK levels were 'consistently lower than expected (both over time and between individuals)' and no root cause was found (Section 3.1).",
      "Figures 2-4 accordingly show 30 mg and 300 mg MAD panels but no 100 mg MAD panel.",
      "sex_female_pct and race_ethnicity are computed from Table 1 over all 66 tozorakimab-treated participants; the paper does not tabulate them for the 60-subject analysis set.",
      "23.7% (216/912) of SAD/MAD PK samples were below the 0.01 ug/mL LLOQ and were imputed at half the LLOQ (Sections 2.2 and 3.3)."
    )
  )

  ini({
    # -- Tozorakimab PK, Table 2 --------------------------------------------
    lka <- log(0.48); label("Absorption rate constant (1/day)") # Sadiq 2024 Table 2: Ka = 0.48 1/day (RSE 10.11%; bootstrap 0.41-0.57)
    lcl <- log(0.87); label("Clearance (L/day)") # Sadiq 2024 Table 2: CL = 0.87 L/day (RSE 14.89%; bootstrap 0.74-1.09)
    lvc <- log(12.64); label("Central volume of distribution (L)") # Sadiq 2024 Table 2: Vc = 12.64 L (RSE 19.91%; bootstrap 10.43-16.70)
    lvp <- fixed(log(2.61)); label("Peripheral volume of distribution (L)") # Sadiq 2024 Table 2: Vp = 2.61 L (fixed)
    lq <- fixed(log(0.21)); label("Inter-compartmental clearance (L/day)") # Sadiq 2024 Table 2: Q = 0.21 L/day (fixed)
    lfdepot <- log(0.45); label("Subcutaneous bioavailability (fraction)") # Sadiq 2024 Table 2: Fsc = 0.45 (RSE 10.08%; bootstrap 0.41-0.48)

    # -- Target engagement, Table 2 ------------------------------------------
    lkdeg <- log(4.83); label("Elimination rate constant of free IL-33 (1/day)") # Sadiq 2024 Table 2: Kel_IL-33 = 4.83 1/day (RSE 45.94%; bootstrap 2.23-9.66)
    lkdeg_sst2 <- log(0.069); label("Elimination rate constant of free sST2 and of the IL-33/sST2 complex (1/day)") # Sadiq 2024 Table 2: Kel_sST2 = 0.069 1/day (RSE 22.54%). The printed IL-33/sST2 complex ODE eliminates the complex with this same constant.
    lkon <- log(68.85); label("Association rate constant of IL-33 with tozorakimab (1/(nM*day))") # Sadiq 2024 Table 2: Kon_IL-33_tozo = 68.85 1/nM/day (RSE 29.07%; bootstrap 49.45-96.60)
    lkon_sst2 <- log(91.38); label("Association rate constant of IL-33 with sST2 (1/(nM*day))") # Sadiq 2024 Table 2: Kon_IL-33_sST2 = 91.38 1/nM/day (RSE 23.57%; bootstrap 64.89-125.62)
    lkd <- fixed(log(0.00003)); label("Equilibrium dissociation constant of the IL-33/tozorakimab pair (nM)") # Sadiq 2024 Table 2: Kd,IL-33_tozo = 0.00003 (fixed). Equilibrium constant, not an off-rate: Koff = Kd * Kon (see model() and vignette Errata).
    lkd_sst2 <- fixed(log(0.0006)); label("Equilibrium dissociation constant of the IL-33/sST2 pair (nM)") # Sadiq 2024 Table 2: Kd,IL-33_sST2 = 0.0006 (fixed). 20-fold weaker than the tozorakimab pair, matching the Discussion's "preferential binding of IL-33red to tozorakimab over sST2".
    lbl_sst2 <- fixed(log(0.20)); label("Baseline free sST2 concentration (nM)") # Sadiq 2024 Table 2: BL_sST2 = 0.20 nM (fixed)
    lbl_complex_sst2 <- log(0.00047); label("Baseline IL-33/sST2 complex concentration (nM)") # Sadiq 2024 Table 2: BL_IL-33/sST2 = 0.00047 nM (RSE 5.29%; bootstrap 0.00044-0.00051)

    # -- Inter-individual variability, Table 2 -------------------------------
    # Table 2 reports IIV as CV%; the table note gives the conversion used,
    # CV% = 100 * (exp(variance) - 1)^(1/2), so variance = log(CV^2 + 1).
    etalcl ~ 0.0423144 # Sadiq 2024 Table 2: IIV CL = 20.79 CV% -> log(0.2079^2 + 1)
    etalkon ~ 0.4224992 # Sadiq 2024 Table 2: IIV Kon_IL-33_tozo = 72.51 CV% -> log(0.7251^2 + 1)
    etalkdeg_sst2 ~ 0.2350058 # Sadiq 2024 Table 2: IIV Kel_sST2 = 51.47 CV% -> log(0.5147^2 + 1)
    etalkon_sst2 ~ 0.1898199 # Sadiq 2024 Table 2: IIV Kon_IL-33_sST2 = 45.72 CV% -> log(0.4572^2 + 1)
    etalbl_complex_sst2 ~ 0.0862878 # Sadiq 2024 Table 2: IIV IL-33/sST2 = 30.02 CV% -> log(0.3002^2 + 1)
    etalfdepot ~ 0.0988288 # Sadiq 2024 Table 2: IIV Fsc = 32.23 CV% -> log(0.3223^2 + 1)

    # -- Residual error, Table 2 ---------------------------------------------
    # Table 2's six residual rows are all labelled "(%)" and are NONMEM $SIGMA
    # variances: read as standard deviations the additive terms are physically
    # impossible (482.11 against an IL-33/tozorakimab complex that never exceeds
    # ~40 pg/mL in Figure 4; 31.51 against an IL-33/sST2 baseline of ~26 pg/mL in
    # Figure 3). Taking them as variances gives additive SDs of 22.0 and 5.6
    # pg/mL, which match the observed scatter in those figures.
    # Only the proportional terms are transferable here: they are scale-free,
    # so propSd = sqrt(reported% / 100). The additive terms carry the paper's
    # ug/mL and pg/mL reporting units, which cannot be converted to this model's
    # nM states without molar masses the paper never reports, so they are NOT
    # encoded. See the vignette "Assumptions and deviations" section.
    propSd <- 0.315753; label("Proportional residual error on tozorakimab (fraction)") # Sadiq 2024 Table 2: Prop error, tozorakimab = 9.97% (RSE 5.86%) -> sqrt(0.0997)
    propSd_complex <- 0.226716; label("Proportional residual error on the IL-33/tozorakimab complex (fraction)") # Sadiq 2024 Table 2: Prop error, IL-33/tozorakimab = 5.14% (RSE 14.65%) -> sqrt(0.0514)
    propSd_Cc_complex_sst2 <- 0.102956; label("Proportional residual error on the IL-33/sST2 complex (fraction)") # Sadiq 2024 Table 2: Prop error, IL-33/sST2 = 1.06% (RSE 17.82%) -> sqrt(0.0106)
  })

  model({
    # -- 1. Individual parameters -------------------------------------------
    ka <- exp(lka)
    cl <- exp(lcl + etalcl)
    vc <- exp(lvc)
    vp <- exp(lvp)
    q <- exp(lq)

    kdeg <- exp(lkdeg)
    kdeg_sst2 <- exp(lkdeg_sst2 + etalkdeg_sst2)
    kon <- exp(lkon + etalkon)
    kon_sst2 <- exp(lkon_sst2 + etalkon_sst2)

    # Table 2 fixes the equilibrium dissociation constants, so the off-rates
    # follow the individual on-rates and hold Kd constant within a subject:
    # Koff = Kd * Kon. Typical values: 0.00003 * 68.85 = 0.00207 1/day and
    # 0.0006 * 91.38 = 0.05483 1/day.
    koff <- exp(lkd) * kon
    koff_sst2 <- exp(lkd_sst2) * kon_sst2

    bl_sst2 <- exp(lbl_sst2)
    bl_complex_sst2 <- exp(lbl_complex_sst2 + etalbl_complex_sst2)

    # -- 2. Baseline free IL-33 ---------------------------------------------
    # Not tabulated. Setting the printed IL-33/sST2 complex ODE to zero at
    # baseline gives BL_IL-33 = (Koff_sST2 + Kel_sST2) * BL_IL-33/sST2 /
    # (Kon_IL-33/sST2 * BL_sST2); typical value 3.184e-06 nM.
    bl_target <- (koff_sst2 + kdeg_sst2) * bl_complex_sst2 / (kon_sst2 * bl_sst2)

    # -- 3. Zero-order production rates -------------------------------------
    # Kin,IL-33 is printed in the paper (Section 2.2). Typical value
    # 4.781e-05 nM/day.
    kin_target <- kdeg * bl_target + kon_sst2 * bl_target * bl_sst2 - koff_sst2 * bl_complex_sst2
    # Kin,sST2 = Kel,sST2 * (BL_sST2 + BL_IL-33/sST2); typical value 0.013832
    # nM/day. The subscripts on this relation are garbled in the published PDF
    # (it prints Kel,IL-33 * (BL_IL-33 + BL_IL-33/sST2)); the form used here is
    # the unique one that holds the untreated system exactly at baseline. See
    # the vignette "Assumptions and deviations" section.
    kin_sst2 <- kdeg_sst2 * (bl_sst2 + bl_complex_sst2)

    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp

    # -- 4. ODE system -------------------------------------------------------
    # Tozorakimab states are amounts (nmol); the four target-engagement states
    # are central-compartment concentrations (nM). Binding is confined to
    # plasma (Section 2.2), so the binding fluxes enter the tozorakimab
    # amount balance multiplied by vc. State expressions are written inline
    # inside d/dt() on purpose -- routing them through named intermediates can
    # silently zero a term in the nlmixr2 model-function form.
    d/dt(depot) <- -ka * depot
    d/dt(central) <- ka * depot - cl * (central / vc) -
      vc * (kon * (central / vc) * target - koff * complex) -
      q * (central / vc) + q * (peripheral1 / vp)
    d/dt(peripheral1) <- q * (central / vc) - q * (peripheral1 / vp)
    d/dt(target) <- kin_target - kdeg * target -
      kon * (central / vc) * target + koff * complex -
      kon_sst2 * target * sst2 + koff_sst2 * complex_sst2
    d/dt(sst2) <- kin_sst2 - kdeg_sst2 * sst2 -
      kon_sst2 * target * sst2 + koff_sst2 * complex_sst2
    d/dt(complex) <- kon * (central / vc) * target - koff * complex -
      kel * complex
    d/dt(complex_sst2) <- kon_sst2 * target * sst2 -
      koff_sst2 * complex_sst2 - kdeg_sst2 * complex_sst2

    # -- 5. Baselines and bioavailability -----------------------------------
    target(0) <- bl_target
    sst2(0) <- bl_sst2
    complex_sst2(0) <- bl_complex_sst2

    f(depot) <- exp(lfdepot + etalfdepot)

    # -- 6. Observations -----------------------------------------------------
    Cc <- central / vc
    Cc_complex <- complex
    Cc_complex_sst2 <- complex_sst2

    Cc ~ prop(propSd)
    Cc_complex ~ prop(propSd_complex)
    Cc_complex_sst2 ~ prop(propSd_Cc_complex_sst2)
  })
}
