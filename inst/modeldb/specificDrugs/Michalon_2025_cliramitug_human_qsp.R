Michalon_2025_cliramitug_human_qsp <- function() {
  description <- "QSP. Translational mechanistic PK/PD model of cardiac amyloid transthyretin (ATTR) depletion by cliramitug (NI006, ALXN2220) in a typical 70 kg patient with ATTR cardiomyopathy. Two-compartment serum PK allometrically scaled from rat, plus a cardiac-interstitium compartment holding antibody and ATTR; antibody-ATTR binding by the quasi-equilibrium approximation; ATTR formed at a constant rate, cleared spontaneously and by phagocytosis at a rate proportional to the fraction of ATTR bound; excess left ventricular wall thickness (LVWT) proportional to ATTR mass. Population-typical deterministic projection without between-subject variability or residual error."
  reference <- "Michalon A, Renaud L, Machacek M, Cortijo C, Udata C, Mercuri MF, Buller F, Hock C, Nitsch RM, Kahr PC, Grimm J. Prediction of Cardiac ATTR Depletion by NI006 (ALXN2220) Using Mechanistic PK/PD Modeling. Clin Pharmacol Ther. 2025;117(1):261. doi:10.1002/cpt.3455. Parameter values from Supplementary Table S1 (human PK also Table 1; disease progression also Table S2); ODEs and unit conversions from Supplementary Data S1 'Human PK/PD simulations' Mlxtran code."
  vignette <- "Michalon_2025_cliramitug"
  # heart_isf is the cardiac interstitial fluid (V3) holding total cliramitug
  # (free plus ATTR-bound); it is the only site of ATTR in the model and so is
  # paper-specific rather than the generic `isf` body-fluid compartment.
  paper_specific_compartments <- c("heart_isf")
  units <- list(time = "day", dosing = "mg", concentration = "ug/mL")

  covariateData <- list()

  compartmentData <- list(
    central = list(
      analyte = "cliramitug (NI006)",
      units = "mg",
      specimen = "serum",
      verified = TRUE
    ),
    peripheral1 = list(
      analyte = "cliramitug (NI006)",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    heart_isf = list(
      analyte = "cliramitug (NI006), total of free and ATTR-bound",
      units = "mg",
      specimen = "tissue",
      verified = TRUE
    ),
    total_target = list(
      analyte = "antibody-accessible cardiac amyloid transthyretin (ATTR), total of free and cliramitug-bound",
      units = "g",
      specimen = "tissue",
      verified = TRUE
    )
  )

  population <- list(
    species = "human",
    n_subjects = NA_integer_,
    n_studies = 0L,
    weight_median = "70 kg (standard individual; mg/kg doses are multiplied by 70)",
    disease_state = "Transthyretin amyloid cardiomyopathy (ATTR-CM), baseline LVWT 16.6 mm",
    dose_range = "Simulated 0.3, 1, 3, 10, 30 and 60 mg/kg IV every 4 weeks as a 2-hour infusion for 12 months (Figure 4); the phase Ib study NI006-101 (NCT04360434) tested the same doses",
    regions = "Not applicable (model-based first-in-human projection)",
    notes = "Forward projection, not a fit to patient data. Human PK scaled from the rat (Michalon_2025_cliramitug_rat) with exponents 1 on V1 and V2, 0.75 on Cl and 2/3 on Q (0.3 kg rat, 70 kg human). Cardiac interstitial volume, biodistribution coefficient and heart exchange clearance from Shah and Betts (2012, 2013). ATTR elimination rate at full occupancy from a 14-day ex vivo patient-myocardium phagocytosis assay; KD from surface plasmon resonance. Disease progression (ATTR formation and spontaneous elimination) fitted by least squares to published LVWT data from untreated ATTR-CM patients (Arvidsson 2015) and the patisiran APOLLO trial (Solomon 2018). Predictions were later compared without recalibration to the NI006-101 phase Ib imaging results (Figure 6)."
  )

  ini({
    # ---- Serum PK (Table S1 / Table 1 'Human value scaled from rat') ----
    lcl <- log(0.146); label("NI006 serum clearance, Cl (L/day)") # Table S1, 0.146 L/day (scaled rat PK parameter)
    lvc <- log(3.37); label("NI006 central volume, V1 (L)") # Table S1, 3.37 L (scaled rat PK parameter)
    lvp <- log(4.42); label("NI006 peripheral volume, V2 (L)") # Table S1, 4.42 L (scaled rat PK parameter)
    lq <- log(0.676); label("NI006 inter-compartmental exchange, Q (L/day)") # Table S1, 0.676 L/day (scaled rat PK parameter)

    # ---- Cardiac interstitium (literature IgG values) ----
    bc_heart <- fixed(0.102); label("Heart biodistribution coefficient, BCh (fraction of serum concentration)") # Table S1, NI006 heart BC 10.2% (Shah 2013)
    v_heart <- fixed(0.0488); label("Heart interstitial volume, Vh (L)") # Table S1, heart volume 0.0488 L (Shah 2012)
    q_heart <- fixed(0.372); label("Serum-heart interstitium exchange clearance, Qh (L/day)") # Table S1, NI006 heart exchange 0.372 L/day (Shah 2012)

    # ---- Binding and antibody-mediated ATTR elimination ----
    kd <- fixed(2.5); label("NI006-ATTR equilibrium dissociation constant, KD (nM)") # Table S1, KD 2.5 nM (SPR)
    lkelmax <- log(0.066); label("Maximum NI006-mediated ATTR elimination rate at 100% occupancy, kelmax (1/day)") # Table S1, 0.066 1/day (tissue amyloid depletion assay); 0.37 1/day xenograft value is the sensitivity case

    # ---- Disease progression (Table S2) and ATTR-LVWT link ----
    lkin <- log(0.43); label("LVWT growth rate from ATTR formation, kgrowth (mm/year)") # Table S1 / Table S2, growth rate 0.43 mm/year
    lkout <- log(0.053); label("Spontaneous ATTR (excess LVWT) elimination rate, kel (1/year)") # Table S1 / Table S2, 0.053 1/year
    lvwt_healthy <- fixed(8.5); label("LVWT in healthy individuals (mm)") # Table S2, baseline LVWT in healthy individuals 8.5 mm
    lvwt0 <- fixed(16.6); label("Initial LVWT of the simulated ATTR-CM patient (mm)") # Results 'Determination of a safe starting dose', baseline LVWT 16.6 mm (Solomon 2018, Clay 2006)
    attr0 <- fixed(1.36); label("Initial antibody-accessible cardiac ATTR mass (g)") # Table S1, initial ATTR mass 1.36 g (1% of 136 g excess LV mass)

    # ---- Molecular weights (Data S1) ----
    mw_ab <- fixed(150000); label("Molecular weight of NI006, IgG1 (g/mol)") # Data S1, MW of NI-301 = 150000 g/mol
    mw_ttr <- fixed(13750); label("Molecular weight of a TTR monomer (g/mol)") # Data S1, MW of TTR = 13750 g/mol
  })

  model({
    # Rates are reported per year; the model runs in days.
    days_per_year <- 365.25

    cl <- exp(lcl)
    vc <- exp(lvc)
    vp <- exp(lvp)
    q <- exp(lq)
    kelmax <- exp(lkelmax)
    kin <- exp(lkin) / days_per_year # mm/day
    kout <- exp(lkout) / days_per_year # 1/day

    # LVWT above the healthy value is proportional to ATTR mass (main-text
    # equation: LVWT = 8.5 + ATTR mass * (16.6 - 8.5) / 1.36). g_per_mm is the
    # inverse slope, Table S1 'k = LVWT per ATTR mass' ~ 6 mm/g.
    g_per_mm <- attr0 / (lvwt0 - lvwt_healthy)
    # ATTR formation rate in g/day (Data S1 TTRkprod_mass, expressed as mass).
    kin_attr <- kin * g_per_mm

    # Data S1 micro-constants
    kel <- cl / vc
    k12 <- q / vc
    k21 <- q / vp
    k1h <- bc_heart * q_heart / vc
    kh1 <- q_heart / v_heart

    # Quasi-equilibrium binding in the heart interstitium, in nM. mg/L times
    # 1e6/MW gives nM for the antibody; g/L times 1e9/MW for TTR monomers.
    ab_heart_nm <- heart_isf / v_heart * 1e6 / mw_ab
    attr_nm <- total_target / v_heart * 1e9 / mw_ttr
    free_heart_nm <- 0.5 * (ab_heart_nm - attr_nm - kd + sqrt((ab_heart_nm - attr_nm - kd)^2 + 4 * kd * ab_heart_nm))
    bound_heart_nm <- ab_heart_nm - free_heart_nm
    fbound <- bound_heart_nm / attr_nm

    # Free and ATTR-bound antibody amounts in the heart, back in mg
    free_heart <- free_heart_nm * v_heart * mw_ab / 1e6
    bound_heart <- bound_heart_nm * v_heart * mw_ab / 1e6

    # ATTR (and the antibody bound to it) is removed spontaneously at kout and
    # by phagocytosis at kelmax scaled by the fraction of ATTR bound.
    kelim_attr <- kout + kelmax * fbound

    d/dt(central) <- -kel * central - k12 * central + k21 * peripheral1 - k1h * central + kh1 * free_heart
    d/dt(peripheral1) <- k12 * central - k21 * peripheral1
    d/dt(heart_isf) <- k1h * central - kh1 * free_heart - kelim_attr * bound_heart
    d/dt(total_target) <- kin_attr - kelim_attr * total_target

    total_target(0) <- attr0

    # Serum concentration (mg/L = ug/mL), ATTR occupancy (%) and LVWT (mm)
    Cc <- central / vc
    Cheart <- heart_isf / v_heart
    occupancy <- 100 * fbound
    lvwt <- total_target / g_per_mm + lvwt_healthy

    # The paper simulates population-typical predictions only: no
    # between-subject variability and no residual error are attached.
  })
}
